#ifndef POLYCAP_SOURCE_H
#define POLYCAP_SOURCE_H

#include "polycap-error.h"
#include "polycap-photon.h"
#include "polycap-description.h"
#include "polycap-rng.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _polycap_source;
typedef struct _polycap_source                      polycap_source;

// get a new polycap_source by providing all its properties
polycap_source* polycap_source_new(
	double d_source,
	double src_x,
	double src_y,
	double src_sigx,
	double src_sigy,
	double src_shiftx,
	double src_shifty,
	polycap_error **error);
//
// free a polycap_source struct
void polycap_source_free(polycap_source *source);

// construct a new random polycap_photon 
polycap_photon* polycap_source_get_photon(polycap_source *source, polycap_description *description, polycap_rng *rng, size_t n_energies, double *energies, polycap_error **error);

#ifdef __cplusplus
}
#endif

#endif
