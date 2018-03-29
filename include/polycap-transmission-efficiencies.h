
#ifndef POLYCAP_TRANSEFF_H
#define POLYCAP_TRANSEFF_H

#include "polycap-error.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _polycap_transmission_efficiencies;
typedef struct _polycap_transmission_efficiencies   polycap_transmission_efficiencies;

// free a polycap_transmission_efficiencies struct
void polycap_transmission_efficiencies_free(polycap_transmission_efficiencies *efficiencies);

// write polycap_transmission_efficiencies data to hdf5 file
void polycap_transmission_efficiencies_write_hdf5(polycap_transmission_efficiencies *efficiencies, const char *filename);

#ifdef __cplusplus
}
#endif

#endif
