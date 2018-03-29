#ifndef POLYCAP_RNG_H
#define POLYCAP_RNG_H

#include "polycap-error.h"
#include "polycap-description.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _polycap_rng; // our rng struct, which will be mapped to either gsl_rng or easy_rng
typedef struct _polycap_rng                         polycap_rng;

// get a new rng
polycap_rng* polycap_rng_new(unsigned long int seed);

// free the rng
void polycap_rng_free(polycap_rng *rng);

#ifdef __cplusplus
}
#endif

#endif
