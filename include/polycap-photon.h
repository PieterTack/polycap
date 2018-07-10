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

#ifndef POLYCAP_PHOTON_H
#define POLYCAP_PHOTON_H

#include "polycap-error.h"
#include "polycap-description.h"
#include "polycap-rng.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	double x;
	double y;
	double z;
} polycap_vector3;


struct _polycap_leaks;
typedef struct _polycap_leaks                       polycap_leaks;	

struct _polycap_photon;
typedef struct _polycap_photon                      polycap_photon;

// construct a new polycap_photon with its initial position, direction, electric field vector and energy
polycap_photon* polycap_photon_new(
	polycap_description *description,
	polycap_rng *rng,
	polycap_vector3 start_coords,
	polycap_vector3 start_direction,
	polycap_vector3 start_electric_vector,
	polycap_error **error); //give full energy range as for each photon a full spectrum transmission is simulated

// simulate a single photon for a given polycap_description
int polycap_photon_launch(
	polycap_photon *photon, 
	size_t n_energies,
	double *energies,
	double **weights,
	polycap_error **error);

// get exit coordinates
polycap_vector3 polycap_photon_get_exit_coords(polycap_photon *photon);

// get exit direction
polycap_vector3 polycap_photon_get_exit_direction(polycap_photon *photon);

// get exit electric vector
polycap_vector3 polycap_photon_get_exit_electric_vector(polycap_photon *photon);

// free a polycap_photon
void polycap_photon_free(polycap_photon *photon);

#ifdef __cplusplus
}
#endif

#endif
