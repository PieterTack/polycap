from __future__ import print_function
from error cimport *
from rng cimport *
from profile cimport *
from description cimport *
from transmission_efficiencies cimport *
from photon cimport *
from libc.string cimport strdup, memcpy
from libc.stdlib cimport free
from cpython cimport Py_DECREF
cimport numpy as np
import numpy as np
import sys

np.import_array()

cdef extern from "Python.h":
    ctypedef void PyObject
    PyObject* PyErr_Occurred()
    void PyErr_SetString(object type, const char *message)

cdef extern from "polycap.h":
    void polycap_free(void *)

cdef extern from "xraylib.h":
    int SymbolToAtomicNumber(const char *symbol)
    void xrlFree(void *)

error_map = {
    POLYCAP_ERROR_MEMORY: MemoryError,
    POLYCAP_ERROR_INVALID_ARGUMENT: ValueError,
    POLYCAP_ERROR_IO: IOError,
    POLYCAP_ERROR_OPENMP: IOError,
    POLYCAP_ERROR_TYPE: TypeError,
    POLYCAP_ERROR_UNSUPPORTED: NotImplementedError,
    POLYCAP_ERROR_RUNTIME: RuntimeError,
}

# this is inspired by h5py...
cdef void set_exception(polycap_error *error) except *:
    if error == NULL:
        return
    eclass = error_map.get(error.code, RuntimeError) 
    PyErr_SetString(eclass, error.message)
    polycap_error_free(error)

cdef class Profile:

    CONICAL = POLYCAP_PROFILE_CONICAL
    PARABOLOIDAL = POLYCAP_PROFILE_PARABOLOIDAL
    ELLIPSOIDAL = POLYCAP_PROFILE_ELLIPSOIDAL

    cdef polycap_profile *profile

    def __cinit__(self,
	polycap_profile_type type,
	double length,
	double rad_ext_upstream,
	double rad_ext_downstream,
	double rad_int_upstream,
	double rad_int_downstream,
	double focal_dist_upstream,
	double focal_dist_downstream):
        cdef polycap_error *error = NULL
        self.profile = polycap_profile_new(
	    type,
	    length,
	    rad_ext_upstream,
	    rad_ext_downstream,
	    rad_int_upstream,
	    rad_int_downstream,
	    focal_dist_upstream,
	    focal_dist_downstream,
            &error)
        set_exception(error)


    def __dealloc__(self):
        if self.profile is not NULL:
            polycap_profile_free(self.profile)

cdef class Rng:
    cdef polycap_rng *rng
    def __cinit__(self, seed=None):
        if seed is None:
            self.rng = polycap_rng_new()
            return
        seed = <unsigned long int?> seed
        self.rng = polycap_rng_new_with_seed(seed)

    def __dealloc__(self):
        if self.rng is not NULL:
            polycap_rng_free(self.rng)

def ensure_int(x):
    if isinstance(x, str):
        Z = SymbolToAtomicNumber(x.encode())
        if Z == 0:
            raise ValueError("Invalid chemical element {} detected".format(x))
    else:
        Z = int(x)
    return Z

cdef class Description:
    cdef polycap_description *description
    cdef object profile_py

    def __cinit__(self, 
        double sig_rough,
        double sig_wave,
        double corr_length,
        int64_t n_cap,
        dict composition,
        double density,
        Profile profile):
        if profile is None:
            raise ValueError("profile cannot be None")
        if len(composition) == 0:
            raise ValueError("composition cannot be empty")
        # convert dict to numpy arrays

        iz = list(map(ensure_int, composition.keys()))
        iz = np.asarray(iz, dtype=np.int32)
        wi = list(map(lambda x: float(x), composition.values()))
        wi = np.asarray(wi, dtype=np.double)

        cdef polycap_error *error = NULL
        self.description = polycap_description_new(
            sig_rough,
            sig_wave,
            corr_length,
            n_cap,
            len(composition),
            <int*> np.PyArray_DATA(iz),
            <double*> np.PyArray_DATA(wi),
            density,
            profile.profile,
            &error)
        set_exception(error)

#        self.profile_py = Profile()
#        self.profile_py.profile = polycap_description_get_profile(self.description)

    def __dealloc__(self):
        if self.description is not NULL:
            polycap_description_free(self.description)
    

cdef class TransmissionEfficiencies:
    cdef polycap_transmission_efficiencies *trans_eff
    cdef object energies_np
    cdef object efficiencies_np

    def __cinit__(self):
        self.trans_eff = NULL
        self.energies_np = None
        self.efficiencies_np = None

    def __dealloc__(self):
        if self.trans_eff is not NULL:
            polycap_transmission_efficiencies_free(self.trans_eff)
        if self.energies_np is not None:
            Py_DECREF(self.energies_np)
        if self.efficiencies_np is not None:
            Py_DECREF(self.efficiencies_np)

    def write_hdf5(self, str filename):
        cdef polycap_error *error = NULL
        polycap_transmission_efficiencies_write_hdf5(self.trans_eff, filename, &error)
        set_exception(error)

    @property
    def data(self):
        return (self.energies_np, self.efficiencies_np)


    # factory method -> these objects cannot be newed, as they are produced via polycap_description_get_transmission_efficiencies
    @staticmethod
    cdef create(polycap_transmission_efficiencies *trans_eff):
        if trans_eff is NULL:
            return None
        rv = TransmissionEfficiencies()
        rv.trans_eff = trans_eff

        cdef size_t n_energies = 0
        cdef double *energies_arr = NULL
        cdef double *efficiencies_arr = NULL
        cdef polycap_error *error = NULL

        polycap_transmission_efficiencies_get_data(rv.trans_eff, &n_energies, &energies_arr, &efficiencies_arr, &error)
        set_exception(error)

        cdef np.npy_intp dims[1]
        dims[0] = n_energies

        rv.energies_np = np.PyArray_EMPTY(1, dims, np.NPY_DOUBLE, False)
        # make read-only
        np.PyArray_CLEARFLAG(rv.energies_np, np.NPY_ARRAY_WRITEABLE) # needs verification
        memcpy(np.PyArray_DATA(rv.energies_np), energies_arr, sizeof(double) * n_energies)
        polycap_free(energies_arr)

        rv.efficiencies_np= np.PyArray_EMPTY(1, dims, np.NPY_DOUBLE, False)
        # make read-only
        np.PyArray_CLEARFLAG(rv.efficiencies_np, np.NPY_ARRAY_WRITEABLE) # needs verification
        memcpy(np.PyArray_DATA(rv.efficiencies_np), efficiencies_arr, sizeof(double) * n_energies)
        polycap_free(efficiencies_arr)
        
        return rv

cdef polycap_vector3 np2vector(np.ndarray[double, ndim=1] arr):
    cdef polycap_vector3 rv
    rv.x = arr[0]
    rv.y = arr[1]
    rv.z = arr[2]
    return rv

cdef tuple vector2tuple(polycap_vector3 vec):
    return (vec.x, vec.y, vec.z)

cdef class Photon:
    cdef polycap_photon *photon

    def __cinit__(self, 
        Description description,
        Rng rng,
        object start_coords,
        object start_direction,
        object start_electric_vector,
        object energies,
        bool ignore=False
        ):
        cdef polycap_vector3 start_coords_pc
        cdef polycap_vector3 start_direction_pc
        cdef polycap_vector3 start_electric_vector_pc
        cdef polycap_error *error = NULL

        if ignore is True:
            self.photon = NULL
        else:
            if description is None:
                raise ValueError("description cannot be None")

            if rng is None:
                raise ValueError("rng cannot be None")

            start_coords = np.asarray(start_coords, dtype=np.double)
            if len(start_coords.shape) != 1 or start_coords.size != 3:
                raise ValueError("start_coords must have exactly three elements")

            start_direction = np.asarray(start_direction, dtype=np.double)
            if len(start_direction.shape) != 1 or start_direction.size != 3:
                raise ValueError("start_direction must have exactly three elements")

            start_electric_vector = np.asarray(start_electric_vector, dtype=np.double)
            if len(start_electric_vector.shape) != 1 or start_electric_vector.size != 3:
                raise ValueError("start_electric_vector must have exactly three elements")

            energies = np.asarray(energies, dtype=np.double)
            energies = np.atleast_1d(energies)
            if len(energies.shape) != 1:
                raise ValueError("energies must be a 1D array")
           
            start_coords_pc = np2vector(start_coords)
            start_direction_pc = np2vector(start_direction)
            start_electric_vector_pc = np2vector(start_electric_vector)

            self.photon = polycap_photon_new(
                description.description,
                rng.rng,
                start_coords_pc,
                start_direction_pc,
                start_electric_vector_pc,
                energies.size,
                <double*> np.PyArray_DATA(energies),
                &error)
            set_exception(error)

    def __dealloc__(self):
        if self.photon is not NULL:
            polycap_photon_free(self.photon)

    def launch(self):
        cdef polycap_error *error = NULL
        rv = polycap_photon_launch(self.photon, &error)
        set_exception(error)
        return rv

    def get_exit_coords(self):
        return vector2tuple(polycap_photon_get_exit_coords(self.photon))

    def get_exit_direction(self):
        return vector2tuple(polycap_photon_get_exit_direction(self.photon))

    def get_exit_electric_vector(self):
        return vector2tuple(polycap_photon_get_exit_electric_vector(self.photon))
