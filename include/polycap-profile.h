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

#ifndef POLYCAP_PROFILE_H
#define POLYCAP_PROFILE_H

#include "polycap-error.h"

#ifdef __cplusplus
extern "C" {
#endif

enum _polycap_profile_type {
	POLYCAP_PROFILE_CONICAL,
	POLYCAP_PROFILE_PARABOLOIDAL,
	POLYCAP_PROFILE_ELLIPSOIDAL,
};

struct _polycap_profile;

typedef enum   _polycap_profile_type                polycap_profile_type;
typedef struct _polycap_profile                     polycap_profile;

// get a new profile for a given type with properties
polycap_profile* polycap_profile_new(
	polycap_profile_type type,
	double length,
	double rad_ext_upstream,
	double rad_ext_downstream,
	double rad_int_upstream,
	double rad_int_downstream,
	double focal_dist_upstream,
	double focal_dist_downstream,
	polycap_error **error);

// get a new profile from Laszlo's ASCII files
// perhaps it would be a good idea to define a new filetype that would combine these three into a single file? Best to use something like XML for convenience... this could then be polycap_profile_new_from_xml
polycap_profile* polycap_profile_new_from_file(
	const char *single_cap_profile_file,
	const char *central_axis_file,
	const char *external_shape_file,
	polycap_error **error);

// free the polycap_profile structure and its associated data
void polycap_profile_free(polycap_profile *profile);

#ifdef __cplusplus
}
#endif

#endif
