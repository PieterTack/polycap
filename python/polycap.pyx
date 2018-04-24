cimport cerror
cimport crng
from libc.string cimport strdup
from libc.stdlib cimport free

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

cdef class Rng:
    cdef crng.polycap_rng *rng
    def __cinit__(self, seed=None):
        if seed is None:
            self.rng = crng.polycap_rng_new()
            return
        seed = <unsigned long int?> seed
        self.rng = crng.polycap_rng_new_with_seed(seed)
    def __dealloc__(self):
        if self.rng is not NULL:
            crng.polycap_rng_free(self.rng)
