from error cimport polycap_error
from profile cimport polycap_profile
from libc.stdint cimport int64_t

cdef extern from "polycap-description.h" nogil:
    ctypedef struct polycap_description

    polycap_description* polycap_description_new(
        polycap_profile *profile,
        double sig_rough,
        double sig_wave,
        double corr_length,
        int64_t n_cap,
        unsigned int nelem,
        int iz[],
        double wi[],
        double density,
	polycap_error **error)

    const polycap_profile* polycap_description_get_profile(polycap_description *description)

    void polycap_description_free(polycap_description *description)
