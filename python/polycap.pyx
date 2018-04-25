from error cimport *
from rng cimport *
from profile cimport *
from transmission_efficiencies cimport *
from libc.string cimport strdup, memcpy
from libc.stdlib cimport free
cimport numpy as np
import numpy as np

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
    POLYCAP_ERROR_OPENMP: IOError
}
# this is inspired by h5py...
cdef void set_exception(polycap_error *error) except *:
    if error == NULL:
        return
    eclass = error_map.get(error.code, RuntimeError) 
    #raise eclass(error.message)
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

    def __cinit__(self):
        trans_eff = NULL

    def __dealloc__(self):
        if self.trans_eff is not NULL:
            polycap_transmission_efficiencies_free(self.trans_eff)

    def write_hdf5(self, str filename):
        cdef polycap_error *error = NULL # FIXME in API!!!
        polycap_transmission_efficiencies_write_hdf5(self.trans_eff, filename)

    @property
    def data(self):
        cdef size_t n_energies = 0
        cdef double *energies_arr = NULL
        cdef double *efficiencies_arr = NULL
        cdef polycap_error *error = NULL
        
        polycap_transmission_efficiencies_get_data(self.trans_eff, &n_energies, &energies_arr, &efficiencies_arr, &error)
        set_exception(error)

        cdef np.ndarray[np.double_t, ndim=1] energies_np = np.empty((n_energies,), dtype=np.double)
        temp = <double*>energies_np.data
        memcpy(temp, energies_arr, sizeof(double) * n_energies)
        polycap_free(energies_arr)

        cdef np.ndarray[np.double_t, ndim=1] efficiencies_np = np.empty((n_energies,), dtype=np.double)
        temp = <double*>efficiencies_np.data
        memcpy(temp, efficiencies_arr, sizeof(double) * n_energies)
        polycap_free(efficiencies_arr)

        return (energies_np, efficiencies_np)

    # factory method -> these objects cannot be newed, as they are produced via polycap_description_get_transmission_efficiencies
    @staticmethod
    cdef create(polycap_transmission_efficiencies *trans_eff):
        if trans_eff == NULL:
            return None
        rv = TransmissionEfficiencies()
        rv.trans_eff = trans_eff
        
        return rv

