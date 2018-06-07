
#ifndef POLYCAP_TRANSEFF_H
#define POLYCAP_TRANSEFF_H

#include "polycap-error.h"
#include "polycap-source.h"
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _polycap_transmission_efficiencies;

typedef struct _polycap_transmission_efficiencies   polycap_transmission_efficiencies;

// free a polycap_transmission_efficiencies struct
void polycap_transmission_efficiencies_free(polycap_transmission_efficiencies *efficiencies);

// write polycap_transmission_efficiencies data to hdf5 file
bool polycap_transmission_efficiencies_write_hdf5(polycap_source *source, polycap_transmission_efficiencies *efficiencies, const char *filename, polycap_error **error);

// extract data from struct. returned arrays should be freed with polycap_free
bool polycap_transmission_efficiencies_get_data(polycap_transmission_efficiencies *efficiencies, size_t *n_energies, double **energies_arr, double **efficiencies_arr, polycap_error **error);

#ifdef __cplusplus
}
#endif

#endif
