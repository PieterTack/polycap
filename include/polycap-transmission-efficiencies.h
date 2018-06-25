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
 *
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

typedef struct _polycap_transmission_efficiencies   polycap_transmission_efficiencies;

// free a polycap_transmission_efficiencies struct
void polycap_transmission_efficiencies_free(polycap_transmission_efficiencies *efficiencies);

// write polycap_transmission_efficiencies data to hdf5 file
bool polycap_transmission_efficiencies_write_hdf5(polycap_transmission_efficiencies *efficiencies, const char *filename, polycap_error **error);

// extract data from struct. returned arrays should be freed with polycap_free
bool polycap_transmission_efficiencies_get_data(polycap_transmission_efficiencies *efficiencies, size_t *n_energies, double **energies_arr, double **efficiencies_arr, polycap_error **error);

#ifdef __cplusplus
}
#endif

#endif
