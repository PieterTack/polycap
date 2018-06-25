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

/** \file polycap-description.h
 * \brief API for dealing with polycap description structures
 *
 * This header contains all functions and definitions that are necessary to create, manipulate and free description structures that are produced by polycap.
 * 
 */


#ifndef POLYCAP_DESCRIPTION_H
#define POLYCAP_DESCRIPTION_H

#include "polycap-error.h"
#include "polycap-profile.h"
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct _polycap_description;
/** Struct containing information about a polycapillary description such as shape and composition
 *
 * When this struct is no longer required, it is the user's responsability to free the memory using polycap_description_free().
 */
typedef struct _polycap_description                 polycap_description;

//TODO: figure out what sig_wave and corr_length actually represent and implement in code, or just get rid of them as parameters...
/** Creates a new polycap_description by providing all its properties.
 *
 * \param profile polycap_profile containing outer polycapillary and single capillary shape coordinates
 * \param sig_rough Surface rougness of the capillaries [$\AA$]
 * \param sig_wave Waviness of the capillaries; currently a dummy variable
 * \param corr_length Correlation length of the capillaries; currently a dummy variable
 * \param n_cap The amount of capillaries in the hexagonally packed polycapillary optic
 * \param nelem The amount of chemical elements present in the capillary matrix
 * \param iz Array of atomic numbers of the elements present in the capillary matrix
 * \param wi Array of weight percentages of the elements present in the capillary matrix
 * \param density Density of the capillary matrix [g/cm$^3$]
 * \param error Struct containing information about an error
 * \returns a new polycap_description
 */
polycap_description* polycap_description_new(
	polycap_profile *profile,
	double sig_rough,
	int64_t n_cap,
	unsigned int nelem,
	int iz[],
	double wi[],
	double density,
	polycap_error **error);

/** Extract the polycap_profile from a polycap_description
 * \param description polycap_description to extract a polycap_profile from
 * \returns a new polycap_profile
 */
const polycap_profile* polycap_description_get_profile(polycap_description *description);

/** free a polycap_description struct
 * \param description polycap_description to free
 */
void polycap_description_free(polycap_description *description);

#ifdef __cplusplus
}
#endif

#endif
