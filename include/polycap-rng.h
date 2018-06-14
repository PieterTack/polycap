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

#ifndef POLYCAP_RNG_H
#define POLYCAP_RNG_H

#include "polycap-error.h"
#include "polycap-description.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _polycap_rng; // our rng struct, which will be mapped to either gsl_rng or easy_rng
typedef struct _polycap_rng                         polycap_rng;

//get a new rng with seed from /dev/urandom or rand_s
polycap_rng* polycap_rng_new(void);

// get a new rng with seed provided by caller
polycap_rng* polycap_rng_new_with_seed(unsigned long int seed);

// free the rng
void polycap_rng_free(polycap_rng *rng);

#ifdef __cplusplus
}
#endif

#endif
