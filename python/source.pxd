from error cimport polycap_error
from photon cimport polycap_photon
from description cimport polycap_description
from rng cimport polycap_rng
from transmission_efficiencies cimport polycap_transmission_efficiencies

cdef extern from "polycap-source.h" nogil:
    ctypedef struct polycap_source

    polycap_source* polycap_source_new(
        polycap_description *description,
        double d_source,
        double src_x,
        double src_y,
        double src_sigx,
        double src_sigy,
        double src_shiftx,
        double src_shifty,
        polycap_error **error)

    void polycap_source_free(polycap_source *source)

    polycap_photon* polycap_source_get_photon(
        polycap_source *source,
        polycap_rng *rng,
        size_t n_energies,
        double *energies,
        polycap_error **error)

    polycap_source* polycap_source_new_from_file(const char *filename, polycap_error **error)

    polycap_transmission_efficiencies* polycap_source_get_transmission_efficiencies(
        polycap_source *source,
        int max_threads,
        size_t n_energies,
        double *energies,
        int n_photons,
        polycap_error **error)

    const polycap_description* polycap_source_get_description(polycap_source *source)

