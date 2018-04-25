from error cimport polycap_error

cdef extern from "stdbool.h" nogil:
    ctypedef bint bool

cdef extern from "polycap-transmission-efficiencies.h" nogil:
    ctypedef struct polycap_transmission_efficiencies

    
    void polycap_transmission_efficiencies_free(polycap_transmission_efficiencies *efficiencies)

    void polycap_transmission_efficiencies_write_hdf5(polycap_transmission_efficiencies *efficiencies, const char *filename)

    bool polycap_transmission_efficiencies_get_data(polycap_transmission_efficiencies *efficiencies, size_t *n_energies, double **energies_arr, double **efficiencies_arr, polycap_error **error)
