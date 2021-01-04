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

/** \file polycap-transmission-efficiencies.h
 * \brief API for dealing with polycap_transmisson_efficiencies as used by polycap
 *
 * This header contains all functions and definitions that are necessary to create, manipulate and free polycap_transmission_efficiencies structures that are used by polycap.
 */

#ifndef POLYCAP_TRANSEFF_H
#define POLYCAP_TRANSEFF_H

#include "polycap-error.h"
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _polycap_transmission_efficiencies;
/** Struct containing all output information such as simulated photon coordinates, direction, energies, weights, ...
 *
 * When this struct is no longer required, it is the user's responsability to free the memory using polycap_transmission_efficiencies_free().
 */
typedef struct _polycap_transmission_efficiencies   polycap_transmission_efficiencies;

/** free a polycap_transmission_efficiencies struct
 *
 * \param efficiencies a polycap_transmission_efficiencies
 */
POLYCAP_EXTERN
void polycap_transmission_efficiencies_free(polycap_transmission_efficiencies *efficiencies);

/** Write polycap_transmission_efficiencies data to a hdf5 file
 *
 * \param efficiencies a polycap_transmission_efficiencies struct
 * \param filename a hdf5 file new
 * \param error a pointer to a \c NULL polycap_error, or \c NULL
 * \returns true or false, or \c NULL if an error occurred
 */
POLYCAP_EXTERN
bool polycap_transmission_efficiencies_write_hdf5(polycap_transmission_efficiencies *efficiencies, const char *filename, polycap_error **error);

/** Extract data from a polycap_transmission_efficiencies struct. returned arrays should be freed by the user with polycap_free() or free().
 *
 * \param efficiencies a polycap_transmission_efficiencies struct
 * \param n_energies a variable to contain the amount of simulated photon energies
 * \param energies_arr a variable to contain the simulated photon energies [keV]
 * \param efficiencies_arr a variable to contain the transmission efficiencies array
 * \param error a pointer to a \c NULL polycap_error, or \c NULL
 * \returns true or false
 */
POLYCAP_EXTERN
bool polycap_transmission_efficiencies_get_data(polycap_transmission_efficiencies *efficiencies, size_t *n_energies, double **energies_arr, double **efficiencies_arr, polycap_error **error);

/** Extract extleak data from a polycap_transmission_efficiencies struct.
 *
 * \param efficiencies a polycap_transmission_efficiencies struct
 * \param leaks a polycap_leak structure array to contain the returned intleak data
 * \param n_leaks a int64_t pointer that will contain the amount of returned intleaks
 * \param error a polycap_error
 * \returns true or false
 */
POLYCAP_EXTERN
bool polycap_transmission_efficiencies_get_extleak_data(polycap_transmission_efficiencies *efficiencies, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error);

/** Extract intleak data from a polycap_transmission_efficiencies struct.
 *
 * \param efficiencies a polycap_transmission_efficiencies struct
 * \param leaks a polycap_leak structure array to contain the returned intleak data
 * \param n_leaks a int64_t pointer that will contain the amount of returned intleaks
 * \param error a polycap_error
 * \returns true or false
 */
POLYCAP_EXTERN
bool polycap_transmission_efficiencies_get_intleak_data(polycap_transmission_efficiencies *efficiencies, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error);

/** Extract photon start data from a polycap_transmission_efficiencies struct.
 *
 * \param efficiencies a polycap_transmission_efficiencies struct
 * \param n_start a int64_t pointer that will contain the amount of returned start events
 * \param start_coords a polycap_vector3 structure array to contain the event start coordinates
 * \param start_direction a polycap_vector3 structure array to contain the event start direction
 * \param start_elecv a polycap_vector3 structure array to contain the event start electric vector
 * \param src_start_coords a polycap_vector3 structure array to contain the event source start coordinates. Z value is relative to optic entrance window at Z=0.
 * \param error a polycap_error
 * \returns true or false
 */
POLYCAP_EXTERN
bool polycap_transmission_efficiencies_get_start_data(polycap_transmission_efficiencies *efficiencies, int64_t *n_start, polycap_vector3 **start_coords, polycap_vector3 **start_direction, polycap_vector3 **start_elecv, polycap_vector3 **src_start_coords, polycap_error **error);

/** Extract photon exit data from a polycap_transmission_efficiencies struct.
 *
 * \param efficiencies a polycap_transmission_efficiencies struct
 * \param n_exit a int64_t pointer that will contain the amount of returned start events
 * \param exit_coords a polycap_vector3 structure array to contain the event start coordinates
 * \param exit_direction a polycap_vector3 structure array to contain the event start direction
 * \param exit_elecv a polycap_vector3 structure array to contain the event start electric vector
 * \param n_refl a int64_t array to contain the amount of internal reflections
 * \param d_travel a double array to contain the traveled distance within the optic [cm]
 * \param n_energies a size_t to contain the amount of simulated energies
 * \param exit_weights a double array to contain the exit efficiency weights for each event as a function of photon energy
 * \param error a polycap_error
 * \returns true or false
 */
POLYCAP_EXTERN
bool polycap_transmission_efficiencies_get_exit_data(polycap_transmission_efficiencies *efficiencies, int64_t *n_exit, polycap_vector3 **exit_coords, polycap_vector3 **exit_direction, polycap_vector3 **exit_elecv, int64_t **n_refl, double **d_travel, size_t *n_energies, double ***exit_weights, polycap_error **error);

#ifdef __cplusplus
}
#endif

#endif
