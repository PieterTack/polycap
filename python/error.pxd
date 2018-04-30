cdef extern from "polycap-error.h" nogil:
    cdef enum polycap_error_code:
        POLYCAP_ERROR_MEMORY             # set in case of a memory allocation problem 
        POLYCAP_ERROR_INVALID_ARGUMENT   # set in case an invalid argument gets passed to a routine
        POLYCAP_ERROR_IO                 # set in case an error involving input/output occurred
        POLYCAP_ERROR_OPENMP             # set in case an error involving OpenMP occurred
        POLYCAP_ERROR_TYPE               # set in case an error involving type conversion occurred (HDF5 related) */
        POLYCAP_ERROR_UNSUPPORTED        # set in case an unsupported feature has been requested */
        POLYCAP_ERROR_RUNTIME            # set in case an unexpected runtime error occurred */

    ctypedef struct polycap_error:
        polycap_error_code code
        char *message

    void polycap_error_free(polycap_error *error)
