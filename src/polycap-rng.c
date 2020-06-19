/*
 * Copyright (C) 2018 Pieter Tack, Tom Schoonjans and Laszlo Vincze
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include "config.h"
#ifdef HAVE_GETTIMEOFDAY
#include <sys/time.h>
#elif defined(HAVE__FTIME)
#include <sys/types.h>
#include <sys/timeb.h>
#endif
#include <stdlib.h>
#include "polycap-private.h"

/* public */

// get a new rng with seed provided by caller
polycap_rng* polycap_rng_new_with_seed(unsigned long int seed) {
	polycap_rng *rng = calloc(1, sizeof(polycap_rng));
	rng->_rng = _polycap_rng_alloc(polycap_rng_mt19937);
	_polycap_rng_set(rng, seed);
	return rng;
}

#if !defined(HAVE_GETTIMEOFDAY) && defined(HAVE__FTIME)
static int gettimeofday (struct timeval *tv, void *tz) {
	struct _timeb timebuf;
	_ftime (&timebuf);
	tv->tv_sec = timebuf.time;
	tv->tv_usec = timebuf.millitm * 1000;

	return 0;
}
#endif

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
polycap_rng * polycap_rng_alloc(const polycap_rng_type * T) {
	polycap_rng *rng = calloc(1, sizeof(polycap_rng));
	rng->_rng = _polycap_rng_alloc(T);
	return rng;
}

void polycap_rng_set(const polycap_rng * r, unsigned long int s) {
	_polycap_rng_set(r, s);
}

double polycap_rng_uniform(const polycap_rng * r) {
	return _polycap_rng_uniform(r);
}

