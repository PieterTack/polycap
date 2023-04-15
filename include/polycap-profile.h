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

/** \file polycap-profile.h
 * \brief API for dealing with polycap_profile
 *
 * This header contains all functions and definitions that are necessary to create, manipulate and free a polycap_profile structure that is used by polycap.
 */

#ifndef POLYCAP_PROFILE_H
#define POLYCAP_PROFILE_H

#include "polycap-error.h"
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Codes to indicate the type of polycapillary external shape
 *
 * In each case, the single capillary shape is assumed conical around the central capillary axis
 *
 * Conical shape: the external profile shape is a straight line between the polycapillary entrance and exit radii.
 *
 * Paraboloidal shape: a third degree polynomial is fit through the polycapillary entrance and exit radii, as well as the linear extrapolation on each side towards the focal distances.
 *
 * Ellipsoidal shape: an ellipse is described where the polycapillary side with largest radius has a horizontal tangent, whereas the tangent at the shortest radius side is directed towards the corresponding polycapillary focal distance.
 *
 */
typedef enum {
	POLYCAP_PROFILE_CONICAL, ///< Conical external shape
	POLYCAP_PROFILE_PARABOLOIDAL, ///< Paraboloidal external shape
	POLYCAP_PROFILE_ELLIPSOIDAL, ///< Ellipsoidal external shape
} polycap_profile_type;

struct _polycap_profile;

/** Struct containing information about a polycapillary profile shape
 *  *
 *   * When this struct is no longer required, it is the user's responsability to free the memory using polycap_profile_free().
 *    */
typedef struct _polycap_profile                     polycap_profile;

/** Create a new profile for a given profile type with supplied polycapillary properties
 *
 * \param type an integer or type that indicates the profile type
 * \param length the polycapillary length, as measured along the central axis [cm]
 * \param rad_ext_upstream external upstream polycapillary radius (photons stream from upstream to downstream) [cm]
 * \param rad_ext_downstream external downstream polycapillary radius (photons stream from upstream to downstream) [cm]. This radius represents the radius of a circle circumscribing the hexagonal optic area.
 * \param rad_int_upstream internal upstream capillary radius (photons stream from upstream to downstream) [cm]. This radius represents the radius of a circle circumscribing the hexagonal optic area.
 * \param rad_int_downstream internal downstream capillary radius (photons stream from upstream to downstream) [cm]
 * \param focal_dist_upstream focal distance upstream of the polycapillary optic (photons stream from upstream to downstream) [cm]
 * \param focal_dist_downstream focal distance downstream of the polycapillary optic (photons stream from upstream to downstream) [cm]
 * \param error a pointer to a \c NULL polycap_error, or \c NULL
 * \returns a new polycap_profile, or \c NULL if an error occurred
 */
POLYCAP_EXTERN
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

/** Create a new profile given ASCII files correponding to the old polycap program format.
 *
 * Each file contains a 2 column data set, preceded by the amount of rows in the data set. The first column contains the Z-coordinate, running from 0 to polycapillary length, and the second column contains the corresponding radius. The Z-coordinates should be the same over all files.
 * In case of the central axis file a 3 column data set is expected, where the second and third column represent the X and Y coordinates of the polycapillary axis respectively.
 *
 * \param single_cap_profile_file filename of an ASCII file containing the single capillary radii [cm]. Default extension is *.prf
 * \param central_axis_file filename of an ASCII file containing the central polycapillary axis coordinates [cm]. Default extension is *.axs.
 * \param external_shape_file filename of and ASCII file containing the external polycapillary radii [cm]. Default extension is *.ext
 * \param error a pointer to a \c NULL polycap_error, or \c NULL
 * \returns a new polycap_profile, or \c NULL if an error occurred
 */
POLYCAP_EXTERN
polycap_profile* polycap_profile_new_from_file(
	const char *single_cap_profile_file,
	const char *central_axis_file,
	const char *external_shape_file,
	polycap_error **error);

/** Checks a profile for inconsistencies between inner capillary coordinates and the external radius.
 *
 * \param profile polycap_profile containing outer polycapillary and single capillary shape coordinates
 * \param n_cap amount of capillaries in the X-ray optic
 * \param error a pointer to a \c NULL polycap_error, or \c NULL
 * \returns an integer: 1 on success (valid profile), 0 on fail, -1 on error
 */
int polycap_profile_validate(polycap_profile *profile, int64_t n_cap, polycap_error **error);


/** Allows the user to set a profile shape given the appropriate shape parameter arrays
 *
 * \param nid amount of elements in the provided array. It is suggested to have at least 1000 elements.
 * \param ext external polycapillary shape profile radii [cm].
 * \param cap central capillary shape profile radii [cm].
 * \param z Z-coordinates running from 0 to polycapillary length [cm].
 * \param error a pointer to a \c NULL polycap_error, or \c NULL
 * \returns a polycap_profile
 */ 
POLYCAP_EXTERN
polycap_profile *polycap_profile_new_from_arrays(int nid, double *ext, double *cap, double *z, polycap_error **error);

/** Retrieve external polycapillary radius profile from a polycap_profile
 * 
 * \param profile a polycap_profile
 * \param nid amount of elements in the array.
 * \param ext external polycapillary shape profile radii [cm].
 * \param error a pointer to a \c NULL polycap_error, or \c NULL
 * \returns true on success, otherwise fail.
 */
POLYCAP_EXTERN
bool polycap_profile_get_ext(polycap_profile *profile, size_t *nid, double **ext, polycap_error **error);

/** Retrieve central capillary radius profile from a polycap_profile
 * 
 * \param profile a polycap_profile
 * \param nid amount of elements in the array.
 * \param cap central capillary shape profile radii [cm].
 * \param error a pointer to a \c NULL polycap_error, or \c NULL
 * \returns true on success, otherwise fail.
 */
POLYCAP_EXTERN
bool polycap_profile_get_cap(polycap_profile *profile, size_t *nid, double **cap, polycap_error **error);

/** Retrieve z array from a polycap_profile
 * 
 * \param profile a polycap_profile
 * \param nid amount of elements in the array.
 * \param z Z-coordinates running from 0 to polycapillary length [cm] 
 * \param error a pointer to a \c NULL polycap_error, or \c NULL
 * \returns true on success, otherwise fail.
 */
POLYCAP_EXTERN
bool polycap_profile_get_z(polycap_profile *profile, size_t *nid, double **z, polycap_error **error);

/** Free the polycap_profile structure and its associated data
 *
 * \param profile a polycap_profile
 */
POLYCAP_EXTERN
void polycap_profile_free(polycap_profile *profile);

#ifdef __cplusplus
}
#endif

#endif
