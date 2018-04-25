from error cimport *
from rng cimport *
from profile cimport *
from transmission_efficiencies cimport *
from libc.string cimport strdup, memcpy
from libc.stdlib cimport free
from cpython cimport Py_DECREF
cimport numpy as np
import numpy as np

np.import_array()

cdef extern from "Python.h":
    ctypedef void PyObject
    PyObject* PyErr_Occurred()
    void PyErr_SetString(object type, const char *message)

cdef extern from "polycap.h":
    void polycap_free(void *)
    
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
        if trans_eff == NULL:
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

