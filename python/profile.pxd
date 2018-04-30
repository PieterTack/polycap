from error cimport polycap_error

cdef extern from "polycap-profile.h" nogil:
    ctypedef struct polycap_profile

    ctypedef enum polycap_profile_type:
        POLYCAP_PROFILE_CONICAL
        POLYCAP_PROFILE_PARABOLOIDAL
        POLYCAP_PROFILE_ELLIPSOIDAL


    polycap_profile* polycap_profile_new(
	polycap_profile_type type,
	double length,
	double rad_ext_upstream,
	double rad_ext_downstream,
	double rad_int_upstream,
	double rad_int_downstream,
	double focal_dist_upstream,
	double focal_dist_downstream,
	polycap_error **error)

    void polycap_profile_free(polycap_profile *profile)
