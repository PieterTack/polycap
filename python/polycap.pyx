# Copyright (C) 2018 Pieter Tack, Tom Schoonjans and Laszlo Vincze
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cdef extern from "config.h":
    char *version "VERSION"

cdef extern from "polycap.h":
    void polycap_free(void *)

from error cimport *
from rng cimport *
from profile cimport *
from description cimport *
from transmission_efficiencies cimport *
from photon cimport *
from source cimport *
from libc.string cimport memcpy
from libc.stdlib cimport free
from cpython cimport Py_DECREF
from collections import namedtuple
cimport numpy as np
import numpy as np
import sys
import os
from pathlib import Path
import threading

__version__ = version.decode("utf-8")

np.import_array()

import urllib
import http.client
import uuid
from uuid import UUID
import platform

import logging
logger = logging.getLogger(__name__)

def __valid_uuid(_uuid):
    try:
        a = UUID(_uuid)
    except ValueError:
        return False
    return True

threading.Thread(target=__send_google_analytics_launch_event).start()

cdef extern from "Python.h":
    ctypedef void PyObject
    PyObject* PyErr_Occurred()
    void PyErr_SetString(object type, const char *message)

cdef extern from "xraylib.h" nogil:
    ctypedef enum xrl_error_code:
        XRL_ERROR_MEMORY
        XRL_ERROR_INVALID_ARGUMENT
        XRL_ERROR_IO
        XRL_ERROR_TYPE
        XRL_ERROR_UNSUPPORTED
        XRL_ERROR_RUNTIME

    ctypedef struct xrl_error:
        xrl_error_code code
        char *message

    int SymbolToAtomicNumber(const char *symbol, xrl_error **error)

    void xrlFree(void *)

    void xrl_error_free(xrl_error *)

    cdef struct compoundData:
        int nElements
        double nAtomsAll
        int *Elements
        double *massFractions
        double *nAtoms
        double molarMass

    void FreeCompoundData(compoundData *cd)

    compoundData* CompoundParser(const char compoundString[], xrl_error **error)

xrl_error_map = {
    XRL_ERROR_MEMORY: MemoryError,
    XRL_ERROR_INVALID_ARGUMENT: ValueError,
    XRL_ERROR_IO: IOError,
    XRL_ERROR_TYPE: TypeError,
    XRL_ERROR_UNSUPPORTED: NotImplementedError,
    XRL_ERROR_RUNTIME: RuntimeError,
}

# this is inspired by h5py...
cdef void xrl_set_exception(xrl_error *error) except *:
    if error == NULL:
        return
    eclass = xrl_error_map.get(error.code, RuntimeError) 
    PyErr_SetString(eclass, error.message)
    xrl_error_free(error)

polycap_error_map = {
    POLYCAP_ERROR_MEMORY: MemoryError,
    POLYCAP_ERROR_INVALID_ARGUMENT: ValueError,
    POLYCAP_ERROR_IO: IOError,
    POLYCAP_ERROR_OPENMP: IOError,
    POLYCAP_ERROR_TYPE: TypeError,
    POLYCAP_ERROR_UNSUPPORTED: NotImplementedError,
    POLYCAP_ERROR_RUNTIME: RuntimeError,
}

# this is inspired by h5py...
cdef void polycap_set_exception(polycap_error *error) except *:
    if error == NULL:
        return
    eclass = polycap_error_map.get(error.code, RuntimeError) 
    PyErr_SetString(eclass, error.message)
    polycap_error_free(error)

'''Class containing information about a polycapillary profile shape
'''
cdef class Profile:
    '''Codes to indicate the type of polycapillary external shape
    
    In each case, the single capillary shape is assumed conical around the central capillary axis
    -Conical shape: the external profile shape is a straight line between the polycapillary entrance and exit radii.
    -Paraboloidal shape: a third degree polynomial is fit through the polycapillary entrance and exit radii, as well as the linear extrapolation on each side towards the focal distances.
    -Ellipsoidal shape: an ellipse is described where the polycapillary side with largest radius has a horizontal tangent, whereas the tangent at the shortest radius side is directed towards the corresponding polycapillary focal distance.
    '''
    CONICAL = POLYCAP_PROFILE_CONICAL
    PARABOLOIDAL = POLYCAP_PROFILE_PARABOLOIDAL
    ELLIPSOIDAL = POLYCAP_PROFILE_ELLIPSOIDAL

    cdef polycap_profile *_profile

    def __cinit__(self,
	polycap_profile_type type,
	double length,
	double rad_ext_upstream,
	double rad_ext_downstream,
	double rad_int_upstream,
	double rad_int_downstream,
	double focal_dist_upstream,
	double focal_dist_downstream,
        bool ignore=False):
        '''Create a new profile for a given profile type with supplied polycapillary properties
        :param type: an integer or type that indicates the profile type
        :param length: the polycapillary length, as measured along the central axis [cm]
        :type length: double
        :param rad_ext_upstream: external upstream polycapillary radius (photons stream from upstream to downstream) [cm]
        :type rad_ext_upstream: double
        :param rad_ext_downstream: external downstream polycapillary radius (photons stream from upstream to downstream) [cm]
        :type rad_ext_downstream: double
        :param rad_int_upstream: internal upstream capillary radius (photons stream from upstream to downstream) [cm]
        :type rad_int_upstream: double
        :param rad_int_downstream: internal downstream capillary radius (photons stream from upstream to downstream) [cm]
        :type rad_int_downstream: double
        :param focal_dist_upstream: focal distance upstream of the polycapillary optic (photons stream from upstream to downstream) [cm]
        :type focal_dist_upstream: double
        :param focal_dist_downstream: focal distance downstream of the polycapillary optic (photons stream from upstream to downstream) [cm]
        :type focal_dist_downstream: double
        :param ignore: if set to True, a \c NULL Profile will be generated
        :type ignore: bool
        :return: a new :ref:``Profile``, or \c NULL if an error occurred
        '''
        if ignore is True:
            self._profile = NULL
            return

        cdef polycap_error *error = NULL
        self._profile = polycap_profile_new(
	    type,
	    length,
	    rad_ext_upstream,
	    rad_ext_downstream,
	    rad_int_upstream,
	    rad_int_downstream,
	    focal_dist_upstream,
	    focal_dist_downstream,
            &error)
        polycap_set_exception(error)

    def __dealloc__(self):
        '''Free the :ref:``Profile`` class and its associated data'''
        if self._profile is not NULL:
            polycap_profile_free(self._profile)

    @classmethod
    def new_from_arrays(cls,
            object ext not None,
            object cap not None,
            object z not None):
        '''Set profile shape using user defined arrays'''
        cdef polycap_error *error = NULL

        ext = np.asarray(ext, dtype=np.double)
        ext = np.atleast_1d(ext)
        cap = np.asarray(cap, dtype=np.double)
        cap = np.atleast_1d(cap)
        z = np.asarray(z, dtype=np.double)
        z = np.atleast_1d(z)

        # check array length consistency
        if(ext.size != z.size or ext.size != cap.size):
            raise ValueError("arrays must be of identical size.")

        cdef polycap_profile *profile = polycap_profile_new_from_arrays(z.size-1, <double*> np.PyArray_DATA(ext), <double*> np.PyArray_DATA(cap), <double*> np.PyArray_DATA(z), &error)
        polycap_set_exception(error)

        rv = Profile(Profile.CONICAL, 0, 0, 0, 0, 0, 0, 0, ignore=True)
        rv._profile = profile

        return rv

    @property
    def get_ext(self):
        '''Retrieve exterior profile from a :ref:``Profile`` class'''
        cdef size_t nid = 0
        cdef double *ext = NULL
        cdef polycap_error *error = NULL
        
        polycap_profile_get_ext(self._profile, &nid, &ext, &error)
        polycap_set_exception(error)
        cdef np.npy_intp dims[1]
        dims[0] = nid+1 

        rv = np.PyArray_EMPTY(1, dims, np.NPY_DOUBLE, False)
        memcpy(np.PyArray_DATA(rv), ext, sizeof(double) * nid+1)
        polycap_free(ext)

        return rv

    @property
    def get_cap(self):
        '''Retrieve capillary profile from a :ref:``Profile`` class'''
        cdef size_t nid = 0
        cdef double *cap = NULL
        cdef polycap_error *error = NULL
        
        polycap_profile_get_cap(self._profile, &nid, &cap, &error)
        polycap_set_exception(error)
        cdef np.npy_intp dims[1]
        dims[0] = nid+1 

        rv = np.PyArray_EMPTY(1, dims, np.NPY_DOUBLE, False)
        memcpy(np.PyArray_DATA(rv), cap, sizeof(double) * nid+1)
        polycap_free(cap)

        return rv

    @property
    def get_z(self):
        '''Retrieve length profile from a :ref:``Profile`` class'''
        cdef size_t nid = 0
        cdef double *z = NULL
        cdef polycap_error *error = NULL
        
        polycap_profile_get_z(self._profile, &nid, &z, &error)
        polycap_set_exception(error)
        cdef np.npy_intp dims[1]
        dims[0] = nid+1 

        rv = np.PyArray_EMPTY(1, dims, np.NPY_DOUBLE, False)
        memcpy(np.PyArray_DATA(rv), z, sizeof(double) * nid+1)
        polycap_free(z)

        return rv

'''Class containing a random number generator

The :ref:``Rng`` class is  mapped to either gsl_rng or easy_rng.``.
'''
cdef class Rng:
    cdef polycap_rng *_rng
    def __cinit__(self, seed=None):
        '''get a new rng with seed provided by caller
        :param seed: a seed provided by the caller, if not provided one is generated for the user
        :type seed: unsigned long int
        :return: a new :ref:``Rng``
        '''
        if seed is None:
            self._rng = polycap_rng_new()
            return
        seed = <unsigned long int?> seed
        self._rng = polycap_rng_new_with_seed(seed)

    def __dealloc__(self):
        '''free a ``Rng`` class'''
        if self._rng is not NULL:
            polycap_rng_free(self._rng)

def ensure_int(x):
    cdef xrl_error *error = NULL
    if isinstance(x, str):
        Z = SymbolToAtomicNumber(x.encode(), &error)
        xrl_set_exception(error)
    else:
        Z = int(x)
    return Z

'''Class containing information about a polycapillary description such as shape and composition
'''
cdef class Description:
    cdef polycap_description *_description
    cdef object _profile_py

    def __cinit__(self, 
        Profile profile not None,
        double sig_rough,
        int64_t n_cap,
        object composition not None,
        double density):
        '''Creates a new polycap_description by providing all its properties.
        :param profile: :ref:``Profile`` containing outer polycapillary and single capillary shape coordinates
        :type profile: polycap_profile
        :param sig_rough: Surface rougness of the capillaries [Angstrom]
        :type sig_rough: double
        :param n_cap: The amount of capillaries in the hexagonally packed polycapillary optic
        :type n_cap: int64_t
        :param composition: capillary material composition XRayLib dictionary or string
        :type composition: object
        :param density: Density of the capillary matrix [g/cm<sup>3</sup>]
        :type density: double
        :return: a new :ref:``Description``, or \c NULL if an error occurred	
        '''

        cdef np.npy_intp dims[1]
        cdef xrl_error *error_xrl = NULL

        if isinstance(composition, dict):
            if len(composition) == 0:
                raise ValueError("composition cannot be empty")
            # convert dict to numpy arrays
            iz = list(map(ensure_int, composition.keys()))
            iz = np.asarray(iz, dtype=np.int32)
            wi = list(map(lambda x: float(x), composition.values()))
            wi = np.asarray(wi, dtype=np.double)
            comp_len = len(composition)
        elif isinstance(composition, str):
            # try xraylib's compoundparser
            cd = CompoundParser(composition.encode(), &error_xrl)
            xrl_set_exception(error_xrl)
            if cd == NULL:
                raise ValueError("composition is not a valid chemical formula")
            dims[0] = cd.nElements
            iz = np.PyArray_EMPTY(1, dims, np.NPY_INT32, False)
            memcpy(np.PyArray_DATA(iz), cd.Elements, sizeof(int) * cd.nElements)
            wi = np.PyArray_EMPTY(1, dims, np.NPY_DOUBLE, False)
            memcpy(np.PyArray_DATA(wi), cd.massFractions, sizeof(double) * cd.nElements)
            comp_len = cd.nElements
            FreeCompoundData(cd)
        else:
            raise TypeError("composition must be a dictionary or a string")

        cdef polycap_error *error = NULL
        self._description = polycap_description_new(
            profile._profile,
            sig_rough,
            n_cap,
            comp_len,
            <int*> np.PyArray_DATA(iz),
            <double*> np.PyArray_DATA(wi),
            density,
            &error)
        polycap_set_exception(error)

#        self._profile_py = Profile()
#        self._profile_py._profile = polycap_description_get_profile(self.description)

    def __dealloc__(self):
        '''free a :ref:``Description`` class and associated data'''
        if self._description is not NULL:
            polycap_description_free(self._description)
    

'''Class containing all output information such as simulated photon coordinates, direction, energies, weights, ...
'''
cdef class TransmissionEfficiencies:
    cdef polycap_transmission_efficiencies *_trans_eff
    cdef object _energies_np
    cdef object _efficiencies_np
    cdef polycap_leak **_int_leaks
    cdef int64_t _n_int_leaks
    cdef polycap_leak **_ext_leaks
    cdef int64_t _n_ext_leaks
    cdef int64_t _n_start
    cdef polycap_vector3 *_start_coords
    cdef polycap_vector3 *_start_direction
    cdef polycap_vector3 *_start_elecv
    cdef polycap_vector3 *_src_start_coords
    cdef int64_t _n_exit
    cdef size_t _n_energies
    cdef polycap_vector3 *_exit_coords
    cdef polycap_vector3 *_exit_direction
    cdef polycap_vector3 *_exit_elecv
    cdef double **_exit_weights
    cdef int64_t *_n_refl
    cdef double *_d_travel

    def __cinit__(self):
        self._trans_eff = NULL
        self._energies_np = None
        self._efficiencies_np = None
        self._ext_leaks = NULL
        self._n_ext_leaks = 0
        self._int_leaks = NULL
        self._n_int_leaks = 0
        self._n_start = 0
        self._start_coords = NULL
        self._start_direction = NULL
        self._start_elecv = NULL
        self._n_exit = 0
        self._exit_coords = NULL
        self._exit_direction = NULL
        self._exit_elecv = NULL
        self._src_start_coords = NULL
        self._exit_weights = NULL
        self._n_refl = NULL
        self._d_travel = NULL

    def __dealloc__(self):
        '''free a :ref:``TransmissionEfficiencies`` class and all associated data'''
        if self._trans_eff is not NULL:
            polycap_transmission_efficiencies_free(self._trans_eff)

        if self._n_ext_leaks > 0:
            for i in range(self._n_ext_leaks):
                polycap_leak_free(self._ext_leaks[i])
            polycap_free(self._ext_leaks)
        
        if self._n_int_leaks > 0:
            for i in range(self._n_int_leaks):
                polycap_leak_free(self._int_leaks[i])
            polycap_free(self._int_leaks)
            
        if self._start_coords != NULL:
            polycap_free(self._start_coords)
        if self._start_direction != NULL:
            polycap_free(self._start_direction)
        if self._start_elecv != NULL:
            polycap_free(self._start_elecv)
        if self._src_start_coords != NULL:
            polycap_free(self._src_start_coords)

        if self._exit_weights != NULL:
            for i in range(self._n_exit):
                polycap_free(self._exit_weights[i])
            polycap_free(self._exit_weights)
        if self._exit_coords != NULL:
            polycap_free(self._exit_coords)
        if self._exit_direction != NULL:
            polycap_free(self._exit_direction)
        if self._exit_elecv != NULL:
            polycap_free(self._exit_elecv)
        if self._n_refl != NULL:
            polycap_free(self._n_refl)
        if self._d_travel != NULL:
            polycap_free(self._d_travel)

        #if self._energies_np is not None:
        #    Py_DECREF(self._energies_np)
        #if self._efficiencies_np is not None:
        #    Py_DECREF(self._efficiencies_np)

    def write_hdf5(self, str filename not None):
        '''Write :ref:``TransmissionEfficiencies`` data to a hdf5 file
        :param filename: a hdf5 file new, not None
	:type filename: str
        :return: true or false, or \c NULL if an error occurred
        '''
        cdef polycap_error *error = NULL
        polycap_transmission_efficiencies_write_hdf5(self._trans_eff, filename.encode(), &error)
        polycap_set_exception(error)

    @property
    def data(self):
        '''Extract data from a polycap_transmission_efficiencies struct. returned arrays should be freed by the user with polycap_free() or free().
        return : tuple of (self._energies_np, self._efficiencies_np)
        '''
        return (self._energies_np, self._efficiencies_np)

    @property
    def extleak_data(self):
        '''Retrieve exterior :ref:``Leak`` class array from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._n_ext_leaks == 0:
            polycap_transmission_efficiencies_get_extleak_data(self._trans_eff, &self._ext_leaks, &self._n_ext_leaks, &error)
            polycap_set_exception(error)

        if self._n_ext_leaks > 0:
            for i in range(self._n_ext_leaks):
                leak = Leak.create(self._ext_leaks[i])
                yield leak
        else:
            return None

    @property
    def intleak_data(self):
        '''Retrieve interior :ref:``Leak`` class array from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._n_int_leaks == 0:
            polycap_transmission_efficiencies_get_intleak_data(self._trans_eff, &self._int_leaks, &self._n_int_leaks, &error)
            polycap_set_exception(error)

        if self._n_int_leaks > 0:
            for i in range(self._n_int_leaks):
                leak = Leak.create(self._int_leaks[i])
                yield leak
        else:
            return None

    @property
    def exit_coords(self):
        '''Retrieve photon exit coordinates vector tuple from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._exit_coords == NULL:
            polycap_transmission_efficiencies_get_exit_data(self._trans_eff, &self._n_exit, &self._exit_coords, &self._exit_direction, &self._exit_elecv, &self._n_refl, &self._d_travel, &self._n_energies, &self._exit_weights, &error)
            polycap_set_exception(error)

        if self._exit_coords != NULL:
            for i in range(self._n_exit):
                yield vector2tuple(self._exit_coords[i])
        else:
            return None

    @property
    def exit_direction(self):
        '''Retrieve photon exit direction vector tuple from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._exit_direction == NULL:
            polycap_transmission_efficiencies_get_exit_data(self._trans_eff, &self._n_exit, &self._exit_coords, &self._exit_direction, &self._exit_elecv, &self._n_refl, &self._d_travel, &self._n_energies, &self._exit_weights, &error)
            polycap_set_exception(error)

        if self._exit_direction != NULL:
            for i in range(self._n_exit):
                yield vector2tuple(self._exit_direction[i])
        else:
            return None

    @property
    def exit_elecv(self):
        '''Retrieve photon exit electric vector vector tuple from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._exit_elecv == NULL:
            polycap_transmission_efficiencies_get_exit_data(self._trans_eff, &self._n_exit, &self._exit_coords, &self._exit_direction, &self._exit_elecv, &self._n_refl, &self._d_travel, &self._n_energies, &self._exit_weights, &error)
            polycap_set_exception(error)

        if self._exit_elecv != NULL:
            for i in range(self._n_exit):
                yield vector2tuple(self._exit_elecv[i])
        else:
            return None

    @property
    def n_refl(self):
        '''Retrieve photon amount of internal reflections from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._n_refl == NULL:
            polycap_transmission_efficiencies_get_exit_data(self._trans_eff, &self._n_exit, &self._exit_coords, &self._exit_direction, &self._exit_elecv, &self._n_refl, &self._d_travel, &self._n_energies, &self._exit_weights, &error)
            polycap_set_exception(error)

        if self._n_refl != NULL:
            for i in range(self._n_exit):
                yield self._n_refl[i]
        else:
            return None

    @property
    def exit_weights(self):
        '''Retrieve photon exit weight from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL
        cdef size_t n_energies = 0

        if self._exit_weights == NULL:
            polycap_transmission_efficiencies_get_exit_data(self._trans_eff, &self._n_exit, &self._exit_coords, &self._exit_direction, &self._exit_elecv, &self._n_refl, &self._d_travel, &self._n_energies, &self._exit_weights, &error)
            polycap_set_exception(error)

        cdef np.npy_intp dims[1]
        if self._exit_weights != NULL:
            dims[0] = self._n_energies
            for i in range(self._n_exit):

                weights_np = np.PyArray_EMPTY(1, dims, np.NPY_DOUBLE, False)
                # make read-only
                weights_np.flags.writeable = False
                memcpy(np.PyArray_DATA(weights_np), self._exit_weights[i], sizeof(double) * self._n_energies)
                yield weights_np
        else:
            return None

    @property
    def d_travel(self):
        '''Retrieve exit photon travel distance within optic from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._d_travel == NULL:
            polycap_transmission_efficiencies_get_exit_data(self._trans_eff, &self._n_exit, &self._exit_coords, &self._exit_direction, &self._exit_elecv, &self._n_refl, &self._d_travel, &self._n_energies, &self._exit_weights, &error)
            polycap_set_exception(error)

        if self._d_travel != NULL:
            for i in range(self._n_exit):
                yield self._d_travel[i]
        else:
            return None

    @property
    def start_coords(self):
        '''Retrieve photon start coordinates vector tuple from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._start_coords == NULL:
            polycap_transmission_efficiencies_get_start_data(self._trans_eff, &self._n_start, &self._n_exit, &self._start_coords, &self._start_direction, &self._start_elecv, &self._src_start_coords, &error)
            polycap_set_exception(error)

        if self._start_coords != NULL:
            for i in range(self._n_exit):
                yield vector2tuple(self._start_coords[i])
        else:
            return None

    @property
    def start_direction(self):
        '''Retrieve photon start direction vector tuple from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._start_direction == NULL:
            polycap_transmission_efficiencies_get_start_data(self._trans_eff, &self._n_start, &self._n_exit, &self._start_coords, &self._start_direction, &self._start_elecv, &self._src_start_coords, &error)
            polycap_set_exception(error)

        if self._start_direction != NULL:
            for i in range(self._n_exit):
                yield vector2tuple(self._start_direction[i])
        else:
            return None

    @property
    def start_elecv(self):
        '''Retrieve photon start electric vector vector tuple from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._start_elecv == NULL:
            polycap_transmission_efficiencies_get_start_data(self._trans_eff, &self._n_start, &self._n_exit, &self._start_coords, &self._start_direction, &self._start_elecv, &self._src_start_coords, &error)
            polycap_set_exception(error)

        if self._start_elecv != NULL:
            for i in range(self._n_exit):
                yield vector2tuple(self._start_elecv[i])
        else:
            return None

    @property
    def src_start_coords(self):
        '''Retrieve source start coordinates vector tuple from a :ref:``TransmissionEfficiencies`` class '''
        if self._trans_eff is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._src_start_coords == NULL:
            polycap_transmission_efficiencies_get_start_data(self._trans_eff, &self._n_start, &self._n_exit, &self._start_coords, &self._start_direction, &self._start_elecv, &self._src_start_coords, &error)
            polycap_set_exception(error)

        if self._src_start_coords != NULL:
            for i in range(self._n_exit):
                yield vector2tuple(self._src_start_coords[i])
        else:
            return None

    # factory method -> these objects cannot be newed, as they are produced via polycap_source_get_transmission_efficiencies
    @staticmethod
    cdef create(polycap_transmission_efficiencies *trans_eff):
        if trans_eff is NULL:
            return None
        rv = TransmissionEfficiencies()
        rv._trans_eff = trans_eff

        cdef size_t n_energies = 0
        cdef double *energies_arr = NULL
        cdef double *efficiencies_arr = NULL
        cdef polycap_error *error = NULL

        polycap_transmission_efficiencies_get_data(rv._trans_eff, &n_energies, &energies_arr, &efficiencies_arr, &error)
        polycap_set_exception(error)

        cdef np.npy_intp dims[1]
        dims[0] = n_energies

        rv._energies_np = np.PyArray_EMPTY(1, dims, np.NPY_DOUBLE, False)
        # make read-only
        #np.PyArray_CLEARFLAGS(rv._energies_np, np.NPY_ARRAY_WRITEABLE) # needs verification
        rv._energies_np.flags.writeable = False
        memcpy(np.PyArray_DATA(rv._energies_np), energies_arr, sizeof(double) * n_energies)
        polycap_free(energies_arr)

        rv._efficiencies_np= np.PyArray_EMPTY(1, dims, np.NPY_DOUBLE, False)
        # make read-only
        #np.PyArray_CLEARFLAGS(rv._efficiencies_np, np.NPY_ARRAY_WRITEABLE) # needs verification
        rv._efficiencies_np.flags.writeable = False
        memcpy(np.PyArray_DATA(rv._efficiencies_np), efficiencies_arr, sizeof(double) * n_energies)
        polycap_free(efficiencies_arr)

        return rv

cdef polycap_vector3 np2vector(np.ndarray[double, ndim=1] arr):
    '''Class describing a 3 dimensional vector where x and y are horizontal and vertical directions compared to the polycapillary respectively, and z is the direction along the polycapillary length
    :param arr: 1D numpy array [x,y,z]
    :type arr: double
    '''
    cdef polycap_vector3 rv
    rv.x = arr[0]
    rv.y = arr[1]
    rv.z = arr[2]
    return rv

VectorTuple = namedtuple('VectorTuple','x y z')

cdef vector2tuple(polycap_vector3 vec):
    return VectorTuple(vec.x, vec.y, vec.z)

'''Class containing information about the simulated leak events such as position and direction, energy and transmission weights.
'''
cdef class Leak:
    cdef object _coords
    cdef object _direction
    cdef object _elecv
    cdef object _weight
    cdef int64_t _n_refl
    cdef size_t _n_energies

    def __cinit__(self):
        self._coords = None
        self._direction = None
        self._elecv = None
        self._weight = None
        self._n_refl = 0
        self._n_energies = 0

    @staticmethod
    cdef create(polycap_leak *leak):
        if leak == NULL:
            return None
        rv = Leak()

        rv._coords = vector2tuple(leak.coords)
        rv._direction = vector2tuple(leak.direction)
        rv._elecv = vector2tuple(leak.elecv)

        cdef np.npy_intp dims[1]
        dims[0] = leak.n_energies
        rv._weight = np.PyArray_EMPTY(1, dims, np.NPY_DOUBLE, False)
        rv._weight.flags.writeable = False

        memcpy(np.PyArray_DATA(rv._weight), leak.weight, sizeof(double) * leak.n_energies)
        rv._n_refl = leak.n_refl
        rv._n_energies = leak.n_energies

        return rv
        # DO NOT FREE leak!


    @property
    def coords(self):
        '''Extract coords data from a polycap_leak struct. returned arrays should be freed by the user with polycap_free() or free().
        return : tuple of self._coords
        '''
        return self._coords

    @property
    def direction(self):
        '''Extract direction data from a polycap_leak struct. returned arrays should be freed by the user with polycap_free() or free().
        return : tuple of self._direction
        '''
        return self._direction

    @property
    def elecv(self):
        '''Extract electric vector data from a polycap_leak struct. returned arrays should be freed by the user with polycap_free() or free().
        return : tuple of self._elecv
        '''
        return self._elecv

    @property
    def weight(self):
        '''Extract weight data from a polycap_leak struct. returned arrays should be freed by the user with polycap_free() or free().
        return : tuple of self._weight
        '''
        return self._weight

    @property
    def n_refl(self):
        '''Extract amount of reflection data from a polycap_leak struct. returned arrays should be freed by the user with polycap_free() or free().
        return : tuple of self._n_refl
        '''
        return self._n_refl

'''Class containing information about the simulated photon such as position and direction, energy and transmission weights.
'''
cdef class Photon:
    cdef polycap_photon *_photon
    cdef polycap_leak **_ext_leaks
    cdef int64_t _n_ext_leaks
    cdef polycap_leak **_int_leaks
    cdef int64_t _n_int_leaks

    def __cinit__(self, 
        Description description,
        object start_coords,
        object start_direction,
        object start_electric_vector,
        bool ignore=False
        ):
        '''Creates a new polycap_photon with its initial position, direction and electric field vector.
        :param description: a :ref:``Description`` class
        :type description: Description
        :param start_coords: photon start coordinates array [xyz]
        :type start_coords: double array
        :param start_direction: photon start direction array [xyz]
        :type start_direction: double array
        :param start_electric_vector: photon start electric field vector array [xyz]
        :type start_electric_vector: double array
        :param ignore: if set to True, a \c NULL Photon will be generated
        :type ignore: bool
        :return: a new :ref:``Photon``, or \c NULL if an error occurred
        '''
        cdef polycap_vector3 start_coords_pc
        cdef polycap_vector3 start_direction_pc
        cdef polycap_vector3 start_electric_vector_pc
        cdef polycap_error *error = NULL

        self._ext_leaks = NULL
        self._n_ext_leaks = 0
        self._int_leaks = NULL
        self._n_int_leaks = 0

        if ignore is True:
            self._photon = NULL
        else:
            if description is None:
                raise ValueError("description cannot be None")

            start_coords = np.asarray(start_coords, dtype=np.double)
            if len(start_coords.shape) != 1 or start_coords.size != 3:
                raise ValueError("start_coords must have exactly three elements")

            start_direction = np.asarray(start_direction, dtype=np.double)
            if len(start_direction.shape) != 1 or start_direction.size != 3:
                raise ValueError("start_direction must have exactly three elements")

            start_electric_vector = np.asarray(start_electric_vector, dtype=np.double)
            if len(start_electric_vector.shape) != 1 or start_electric_vector.size != 3:
                raise ValueError("start_electric_vector must have exactly three elements")

            start_coords_pc = np2vector(start_coords)
            start_direction_pc = np2vector(start_direction)
            start_electric_vector_pc = np2vector(start_electric_vector)

            self._photon = polycap_photon_new(
                description._description,
                start_coords_pc,
                start_direction_pc,
                start_electric_vector_pc,
                &error)
            polycap_set_exception(error)

    def __dealloc__(self):
        '''Free a :ref:``Photon`` class and all associated data'''
        if self._photon is not NULL:
            polycap_photon_free(self._photon)

        if self._n_ext_leaks > 0:
            for i in range(self._n_ext_leaks):
                polycap_leak_free(self._ext_leaks[i])
            polycap_free(self._ext_leaks)
        
        if self._n_int_leaks > 0:
            for i in range(self._n_int_leaks):
                polycap_leak_free(self._int_leaks[i])
            polycap_free(self._int_leaks)

    def launch(self,
        object energies not None,
        bool leak_calc = False):
        '''Simulate a single photon trajectory for a given :ref:``Description``. For each single photon the transmission efficiency for all energies is calculated
        Weights memory is allocated by this function and should be freed by the user.

        :param energies: an array containing the discrete energies for which the transmission efficiency will be calculated [keV]
        :type energies: double array
        :param leak_calc: True: perform leak calculation; False: do not perform leak calculation
        :type leak_calc: bool
        :return: weights array that will contain the transmission efficiency values as a function of photon energy
        '''

        energies = np.asarray(energies, dtype=np.double)
        energies = np.atleast_1d(energies)
        if len(energies.shape) != 1:
            raise ValueError("energies must be a 1D array")

        cdef polycap_error *error = NULL
        cdef double *weights = NULL
           
        rv = polycap_photon_launch(self._photon, energies.size, <double*> np.PyArray_DATA(energies), &weights, leak_calc, &error)
        polycap_set_exception(error)
        if rv == 2:
            return None

        # copy weights to numpy array, free and return
        cdef np.npy_intp dims[1]
        dims[0] = energies.size 
        weights_np = np.PyArray_EMPTY(1, dims, np.NPY_DOUBLE, False)
        memcpy(np.PyArray_DATA(weights_np), weights, sizeof(double) * energies.size)
        polycap_free(weights)

        return weights_np

    @property
    def extleak_data(self):
        '''Retrieve exterior :ref:``Leak`` class array from a :ref:``Photon`` class '''
        if self._photon is NULL:
            return None

        cdef polycap_error *error = NULL

        if self._n_ext_leaks == 0:
            polycap_photon_get_extleak_data(self._photon, &self._ext_leaks, &self._n_ext_leaks, &error)
            polycap_set_exception(error)

        if self._n_ext_leaks > 0:
            for i in range(self._n_ext_leaks):
                leak = Leak.create(self._ext_leaks[i])
                yield leak
        else:
            return None

    @property
    def intleak_data(self):
        '''Retrieve interior :ref:``Leak`` class array from a :ref:``Photon`` class '''
        if self._photon is NULL:
            return None

        cdef polycap_error *error = NULL

        #logger.debug('Before calling C')
        if self._n_int_leaks == 0:
            polycap_photon_get_intleak_data(self._photon, &self._int_leaks, &self._n_int_leaks, &error)
            polycap_set_exception(error)
        #logger.debug('After calling C') 

        rv = list()

        if self._n_int_leaks > 0:
            for i in range(self._n_int_leaks):
                leak = Leak.create(self._int_leaks[i])
                yield leak
        else:
            return None

    @property
    def get_exit_coords(self):
        '''Retrieve exit coordinates from a :ref:``Photon`` class'''
        return self.exit_coords

    @property
    def exit_coords(self):
        '''Retrieve exit coordinates from a :ref:``Photon`` class'''
        return vector2tuple(polycap_photon_get_exit_coords(self._photon))

    @property
    def start_coords(self):
        '''Retrieve start coordinates from a :ref:``Photon`` class'''
        return vector2tuple(polycap_photon_get_start_coords(self._photon))

    @property
    def get_exit_direction(self):
        '''Retrieve exit direction from a :ref:``Photon`` class'''
        return self.exit_direction

    @property
    def exit_direction(self):
        '''Retrieve exit direction from a :ref:``Photon`` class'''
        return vector2tuple(polycap_photon_get_exit_direction(self._photon))

    @property
    def start_direction(self):
        '''Retrieve start direction from a :ref:``Photon`` class'''
        return vector2tuple(polycap_photon_get_start_direction(self._photon))

    @property
    def get_exit_electric_vector(self):
        '''Retrieve exit electric field vector from a :ref:``Photon`` class'''
        return self.exit_electric_vector

    @property
    def exit_electric_vector(self):
        '''Retrieve exit electric field vector from a :ref:``Photon`` class'''
        return vector2tuple(polycap_photon_get_exit_electric_vector(self._photon))

    @property
    def start_electric_vector(self):
        '''Retrieve start electric field vector from a :ref:``Photon`` class'''
        return vector2tuple(polycap_photon_get_start_electric_vector(self._photon))

    @property
    def i_refl(self):
        '''Retrieve i_refl from a :ref:``Photon`` class'''
        return np.double(polycap_photon_get_irefl(self._photon))

    @property
    def d_travel(self):
        '''Retrieve d_travel from a :ref:``Photon`` class'''
        return polycap_photon_get_dtravel(self._photon)

'''Class containing information on the source from which photons can be (randomly) selected
'''
cdef class Source:
    cdef polycap_source *_source

    def __cinit__(self, 
        Description description not None,
        double d_source,
        double src_x,
        double src_y,
        double src_sigx,
        double src_sigy,
        double src_shiftx,
        double src_shifty,
	double hor_pol,
	object energies not None):
        '''Creates a new Source class by providing all its properties
        :param description: a :ref:``Description`` class
        :type description: Description
        :param d_source: the distance between the source and polycapillary optic entrance window along the central axis [cm]
        :type d_source: double
        :param src_x: the source radius along the X (horizontal) direction [cm]
        :type src_x: double
        :param src_y: the source radius along the y (vertical) direction [cm]
        :type src_y: double
        :param src_sigx: the maximal divergence of photons along the X (horizontal) direction [rad]. Negative values in src_sigx or src_sigy represent homogeneous polycapillary optic illumination.
        :type src_sigx: double
        :param src_sigy: the maximal divergence of photons along the Y (vertical) direction [rad]. Negative values in src_sigx or src_sigy represent homogeneous polycapillary optic illumination.
        :type src_sigy: double
        :param src_shiftx: lateral shift of the source centre along the X (horizontal) direction with respect to the polycapillary optic central axis [cm]
        :type src_shiftx: double
        :param src_shifty: lateral shift of the source centre along the Y (vertical) direction with respect to the polycapillary optic central axis [cm]
        :type src_shifty: double
        :param hor_pol: the polarisation factor of the simulated source (-1 <= hor_pol <= 1)
        :type hor_pol: double
        :param energies: an array containing the discrete energies of which the source will emit photons
        :type energies: double array
        :return: a new :ref:``Source``, or \c NULL if an error occurred
        '''

        energies = np.asarray(energies, dtype=np.double)
        energies = np.atleast_1d(energies)
        if len(energies.shape) != 1:
            raise ValueError("energies must be a 1D array")


        cdef polycap_error *error = NULL
        self._source = polycap_source_new(
            description._description,
            d_source,
            src_x,
            src_y,
            src_sigx,
            src_sigy,
            src_shiftx,
            src_shifty,
            hor_pol,
            energies.size,
            <double*> np.PyArray_DATA(energies),
            &error)
        polycap_set_exception(error)

    def __dealloc__(self):
        '''free a :ref:``Source`` class and associated data'''
        if self._source is not NULL:
            polycap_source_free(self._source)

    def get_photon(self,
        Rng rng not None):
        '''Create a new random polycap_photon based on :ref:``Source``
        In the event of an error, \c NULL is returned and \c error is set appropriately.
        :param rng : a :ref:``Rng`` class, not None
        :type rng: Rng
    	:return : a new polycap_photon, or \c NULL if an error occurred
        '''

        cdef polycap_error *error = NULL
        cdef polycap_photon *photon = polycap_source_get_photon(
            self._source,
            rng._rng,
            &error)
        polycap_set_exception(error)

        rv = Photon(None, None, None, None, ignore=True)
        rv._photon = photon

        return rv

    def get_transmission_efficiencies(self,
        int max_threads,
        int n_photons,
        bool leak_calc = False):
        '''Obtain the transmission efficiencies for a given array of energies, and a full polycap_description.
        :param max_threads: the amount of threads to use. Set to -1 to use the maximum available amount of threads.
        :type max_threads: int
        :param n_photons: the amount of photons to simulate that reach the polycapillary end
        :type n_photons: int
        :param leak_calc: True: perform leak calculation; False: do not perform leak calculation
        :type leak_calc: bool
        :return: a new :ref:``TransmissionEfficiencies`` class, or \c NULL if an error occurred
        '''

        cdef polycap_error *error = NULL
        cdef polycap_transmission_efficiencies *transmission_efficiencies = polycap_source_get_transmission_efficiencies(
            self._source,
            max_threads,
            n_photons,
            leak_calc, #leak_calc option
            NULL, # polycap_progress_monitor
            &error)
        polycap_set_exception(error)

        return TransmissionEfficiencies.create(transmission_efficiencies)


