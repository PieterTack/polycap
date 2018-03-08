
#include "polycap-private.h"
#include <stdlib.h>

/* public */

// get a new rng
polycap_rng* polycap_rng_new(unsigned long int seed) {
	polycap_rng *rng = malloc(sizeof(polycap_rng));
	rng->_rng = _polycap_rng_alloc(polycap_rng_mt19937);
	_polycap_rng_set(rng, seed);
	return rng;
}

// free the rng
void polycap_rng_free(polycap_rng *rng) {
	_polycap_rng_free(rng);
	free(rng);
}

/* private */
polycap_rng * polycap_rng_alloc(const polycap_rng_type * T) {
	polycap_rng *rng = malloc(sizeof(polycap_rng));
	rng->_rng = _polycap_rng_alloc(T);
	return rng;
}

void polycap_rng_set(const polycap_rng * r, unsigned long int s) {
	_polycap_rng_set(r, s);
}
double polycap_rng_uniform(const polycap_rng * r) {
	return _polycap_rng_uniform(r);
}

