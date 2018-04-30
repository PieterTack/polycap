#ifndef POLYCAP_DESCRIPTION_H
#define POLYCAP_DESCRIPTION_H

#include "polycap-error.h"
#include "polycap-profile.h"
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _polycap_description;
typedef struct _polycap_description                 polycap_description;

// get a new polycap_description by providing all its properties... perhaps a simpler variant of this function could be defined that would only set the most important parameters and use defaults for the others??
polycap_description* polycap_description_new(
	double sig_rough,
	double sig_wave,
	double corr_length,
	int64_t n_cap,
	unsigned int nelem,
	int iz[],
	double wi[],
	double density,
	polycap_profile *profile,
	polycap_error **error);

// get the polycap_profile from a polycap_description
const polycap_profile* polycap_description_get_profile(polycap_description *description);

// free a polycap_description struct
void polycap_description_free(polycap_description *description);

#ifdef __cplusplus
}
#endif

#endif
