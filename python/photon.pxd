from error cimport polycap_error
from description cimport polycap_description
from rng cimport polycap_rng

cdef extern from "polycap-photon.h" nogil:

    ctypedef struct polycap_vector3:
        double x
        double y
        double z

    ctypedef struct polycap_photon

    polycap_photon* polycap_photon_new(
        polycap_description *description,
        polycap_rng *rng,
        polycap_vector3 start_coords,
        polycap_vector3 start_direction,
        polycap_vector3 start_electric_vector,
        size_t n_energies,
        double *energies,
        polycap_error **error)

    int polycap_photon_launch(polycap_photon *photon, polycap_error **error)

    polycap_vector3 polycap_photon_get_exit_coords(polycap_photon *photon)

    polycap_vector3 polycap_photon_get_exit_direction(polycap_photon *photon)

    polycap_vector3 polycap_photon_get_exit_electric_vector(polycap_photon *photon)

    void polycap_photon_free(polycap_photon *photon)
