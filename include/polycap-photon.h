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

/** \file polycap-photon.h
 * \brief API for dealing with polycap_photon structures
 *
 * This header contains all functions and definitions that are necessary to create, manipulate and free polycap_photon structures that are used by polycap
 */

#ifndef POLYCAP_PHOTON_H
#define POLYCAP_PHOTON_H

#include "polycap-error.h"
#include "polycap-description.h"
#include "polycap-rng.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Struct describing a 3 dimensional vector where x and y are horizontal and vertical directions compared to the polycapillary respectively, and z is the direction along the polycapillary length
 */
typedef struct {
	double x; ///< horizontal X-axis vector contribution
	double y; ///< vertical Y-axis vector contribution
	double z; ///< Z-axis (perpendicular to X and Y) vector contribution
} polycap_vector3;


struct _polycap_leaks;
typedef struct _polycap_leaks                       polycap_leaks;	

struct _polycap_photon;
/** Struct containing information about the simulated photon such as position and direction, energy and transmission weights.
 *
 * When this struct is no longer required, it is the user's responsability to free the memory using polycap_photon_free().
 */
typedef struct _polycap_photon                      polycap_photon;

/** Creates a new polycap_photon with its initial position, direction and electric field vector.
 *
 * \param description a polycap_description
 * \param rng a random number generator pointer
 * \param start_coords photon start coordinates
 * \param start_direction photon start direction
 * \param start_electric_vector photon start electric field vector
 * \param error Struct containing information about an error
 * \returns a new polycap_photon
 */
polycap_photon* polycap_photon_new(
	polycap_description *description,
	polycap_rng *rng,
	polycap_vector3 start_coords,
	polycap_vector3 start_direction,
	polycap_vector3 start_electric_vector,
	polycap_error **error);

/** Simulate a single photon trajectory for a given polycap_description. For each single photon the transmission efficiency for all energies is calculated
 *
 * Weights memory is allocated by this function and should be freed by the user.
 *
 * \param photon a polycap_photon
 * \param n_energies the amount of discrete energies for which the transmission efficiency will be calculated
 * \param energies an array containing the discrete energies for which the transmission efficiency will be calculated [keV]
 * \param weights an array that will contain the transmission efficiency values
 * \returns an int: 0 if photon was absorbed by the polycapillary, 1 if photon reached the end of the polycapillary, -1 on error
 */
int polycap_photon_launch(
	polycap_photon *photon, 
	size_t n_energies,
	double *energies,
	double **weights,
	bool leak_calc,
	polycap_error **error);

/** Retrieve exit coordinates from a polycap_photon
 * 
 * \param photon a polycap_photon
 * \returns exit coordinates
 */
polycap_vector3 polycap_photon_get_exit_coords(polycap_photon *photon);

/** Retrieve exit direction vector from a polycap_photon
 * 
 * \param photon a polycap_photon
 * \returns exit direction vector
 */
polycap_vector3 polycap_photon_get_exit_direction(polycap_photon *photon);

/** Retrieve exit electric field vector from a polycap_photon
 * 
 * \param photon a polycap_photon
 * \returns exit electric field vector
 */
polycap_vector3 polycap_photon_get_exit_electric_vector(polycap_photon *photon);

/** Free a polycap_photon
 * 
 * \param photon a polycap_photon
 */
void polycap_photon_free(polycap_photon *photon);

#ifdef __cplusplus
}
#endif

#endif
