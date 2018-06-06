#include "polycap-private.h"
#include <sys/time.h>
#include <stdlib.h>

/* public */

// get a new rng with seed provided by caller
polycap_rng* polycap_rng_new_with_seed(unsigned long int seed) {
	polycap_rng *rng = calloc(1, sizeof(polycap_rng));
	rng->_rng = _polycap_rng_alloc(polycap_rng_mt19937);
	_polycap_rng_set(rng, seed);
	return rng;
}

//get a new rng with seed from /dev/urandom or rand_s
polycap_rng* polycap_rng_new() {

#ifdef _WIN32
	unsigned int seed;
	if (rand_s(&seed) != 0)
#else
	unsigned long int seed;
	FILE *random_device = fopen("/dev/urandom","r");
	if(random_device != NULL) {
		fread(&seed, sizeof(unsigned long int), 1, random_device);
		fclose(random_device);
	}
	else
#endif
	{
		struct timeval tv;
		gettimeofday(&tv, NULL);
		seed = tv.tv_sec % tv.tv_usec;
	}

	return polycap_rng_new_with_seed(seed);
}


// free the rng
void polycap_rng_free(polycap_rng *rng) {
	if (rng == NULL)
		return;
	_polycap_rng_free(rng);
	free(rng);
}

/* private */
HIDDEN polycap_rng * polycap_rng_alloc(const polycap_rng_type * T) {
	polycap_rng *rng = calloc(1, sizeof(polycap_rng));
	rng->_rng = _polycap_rng_alloc(T);
	return rng;
}

HIDDEN void polycap_rng_set(const polycap_rng * r, unsigned long int s) {
	_polycap_rng_set(r, s);
}

HIDDEN double polycap_rng_uniform(const polycap_rng * r) {
	return _polycap_rng_uniform(r);
}

