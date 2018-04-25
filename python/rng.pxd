cdef extern from "polycap-rng.h" nogil:
    ctypedef struct polycap_rng

    polycap_rng* polycap_rng_new()
    polycap_rng* polycap_rng_new_with_seed(unsigned long int seed)

    void polycap_rng_free(polycap_rng *rng)
