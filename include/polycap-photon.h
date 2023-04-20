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


struct _polycap_photon;
/** Struct containing information about the simulated photon such as position and direction, energy and transmission weights.
 *
 * When this struct is no longer required, it is the user's responsability to free the memory using polycap_photon_free().
 */
typedef struct _polycap_photon                      polycap_photon;

/** Struct containing information about the simulated photon leak events such as position and direction, energy and transmission weights.
 *
 * When this struct is no longer required, it is the user's responsability to free the memory using polycap_leak_free().
 */
typedef struct {
	polycap_vector3 coords;
  	polycap_vector3 direction;
  	polycap_vector3 elecv;
  	size_t n_energies;
  	double *weight;
  	int64_t n_refl;
} polycap_leak;


/** Creates a new polycap_photon with its initial position, direction and electric field vector.
 *
 * \param description a polycap_description
 * \param start_coords photon start coordinates
 * \param start_direction photon start direction
 * \param start_electric_vector photon start electric field vector
 * \param error a pointer to a \c NULL polycap_error, or \c NULL
 * \returns a new polycap_photon, or \c NULL if an error occurred
 */
POLYCAP_EXTERN
polycap_photon* polycap_photon_new(
	polycap_description *description,
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
 * \param leak_calc True: perform leak calculation; False: do not perform leak calculation
 * \param error a pointer to a \c NULL polycap_error, or \c NULL
 * \returns an int: 0 if photon was absorbed by the polycapillary, 1 if photon reached the end of the polycapillary, 2 if photon hits capillary wall on entrace, -2 if photon is not propagating towards optic entrance at start, -1 on error
 */
POLYCAP_EXTERN
int polycap_photon_launch(
	polycap_photon *photon, 
	size_t n_energies,
	double *energies,
	double **weights,
	bool leak_calc,
	polycap_error **error);

/** Retrieve start coordinates from a polycap_photon
 * 
 * \param photon a polycap_photon
 * \returns start coordinates
 */
POLYCAP_EXTERN
polycap_vector3 polycap_photon_get_start_coords(polycap_photon *photon);

/** Retrieve start direction vector from a polycap_photon
 * 
 * \param photon a polycap_photon
 * \returns start direction vector
 */
POLYCAP_EXTERN
polycap_vector3 polycap_photon_get_start_direction(polycap_photon *photon);

/** Retrieve start electric field vector from a polycap_photon
 * 
 * \param photon a polycap_photon
 * \returns start electric field vector
 */
POLYCAP_EXTERN
polycap_vector3 polycap_photon_get_start_electric_vector(polycap_photon *photon);

/** Retrieve exit coordinates from a polycap_photon
 * 
 * \param photon a polycap_photon
 * \returns exit coordinates
 */
POLYCAP_EXTERN
polycap_vector3 polycap_photon_get_exit_coords(polycap_photon *photon);

/** Retrieve exit direction vector from a polycap_photon
 * 
 * \param photon a polycap_photon
 * \returns exit direction vector
 */
POLYCAP_EXTERN
polycap_vector3 polycap_photon_get_exit_direction(polycap_photon *photon);

/** Retrieve exit electric field vector from a polycap_photon
 * 
 * \param photon a polycap_photon
 * \returns exit electric field vector
 *
 */
POLYCAP_EXTERN
polycap_vector3 polycap_photon_get_exit_electric_vector(polycap_photon *photon);

/** Retrieve d_travel [cm] from a polycap_photon
 *
 * \param photon a polycap_photon
 * \returns photon distance traveled within the optic [cm]
 *
 */
POLYCAP_EXTERN
double polycap_photon_get_dtravel(polycap_photon *photon);

/** Retrieve i_refl from a polycap_photon
 *
 * \param photon a polycap_photon
 * \returns amount of reflection events along photon path within the optic [int64_t]
 *
 */
POLYCAP_EXTERN
int64_t polycap_photon_get_irefl(polycap_photon *photon);

/** Retrieve extleak events from a polycap_photon
 *
 * \param photon a polycap photon
 * \param leaks a polycap_leak structure array to contain the returned extleak data
 * \param n_leaks a int64_t pointer that will contain the amount of returned extleaks
 * \param error a polycap_error
 * \returns true if succesful, false on error
 *
 */
POLYCAP_EXTERN
bool polycap_photon_get_extleak_data(polycap_photon *photon, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error);

/** Retrieve intleak events from a polycap_photon
 *
 * \param photon a polycap photon
 * \param leaks a polycap_leak structure array to contain the returned intleak data
 * \param n_leaks a int64_t pointer that will contain the amount of returned intleaks
 * \param error a polycap_error
 * \returns true if succesful, false on error
 *
 */
POLYCAP_EXTERN
bool polycap_photon_get_intleak_data(polycap_photon *photon, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error);

/** Free a polycap_photon
 * 
 * \param photon a polycap_photon
 */
POLYCAP_EXTERN
void polycap_photon_free(polycap_photon *photon);

/** Free a polycap_leak
 * 
 * \param leak a polycap_leak
 */
POLYCAP_EXTERN
void polycap_leak_free(polycap_leak *leak);

#ifdef __cplusplus
}
#endif

#endif
