from cerror cimport *
from crng cimport *
from cprofile cimport *
from libc.string cimport strdup
from libc.stdlib cimport free

cdef extern from "Python.h":
    ctypedef void PyObject
    PyObject* PyErr_Occurred()
    void PyErr_SetString(object type, const char *message)


#cdef class PolycapException(Exception):
#    cdef cerror.polycap_error_code code
#    cdef char *message
#
#    def __cinit__(self):
#        self.message = NULL
#
#    cdef _set(self, cerror.polycap_error *error):
#        self.code = error.code
#        self.message = strdup(error.message)
#
#    def __dealloc__(self):
#        if self.message is not NULL:
#            free(self.message)
#
#    def __str__(self):
#        return self.message
#
#    # factory
#    @staticmethod
#    cdef create(cerror.polycap_error *error):
#        if error == NULL:
#            return None
#        rv = PolycapException()
#        rv._set(error)
#        return rv


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
