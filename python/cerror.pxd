cdef extern from "polycap-error.h" nogil:
    cdef enum polycap_error_code:
        POLYCAP_ERROR_MEMORY             # set in case of a memory allocation problem 
        POLYCAP_ERROR_INVALID_ARGUMENT   # set in case an invalid argument gets passed to a routine
        POLYCAP_ERROR_IO                 # set in case an error involving input/output occurred
        POLYCAP_ERROR_OPENMP             # set in case an error involving OpenMP occurred

    ctypedef struct polycap_error:
        polycap_error_code code
        char *message

    void polycap_free(void *)    
