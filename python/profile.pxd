# Copyright (C) 2018 Pieter Tack, Tom Schoonjans and Laszlo Vincze
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from error cimport polycap_error

cdef extern from "polycap-profile.h" nogil:
    '''Struct containing information about a polycapillary profile shape
    When this struct is no longer required, it is the user's responsability to free the memory using polycap_profile_free().
    '''
    ctypedef struct polycap_profile

    '''Codes to indicate the type of polycapillary external shape
    In each case, the single capillary shape is assumed conical around the central capillary axis
    -Conical shape: the external profile shape is a straight line between the polycapillary entrance and exit radii.
    -Paraboloidal shape: a third degree polynomial is fit through the polycapillary entrance and exit radii, as well as the linear extrapolation on each side towards the focal distances.
    -Ellipsoidal shape: an ellipse is described where the polycapillary side with largest radius has a horizontal tangent, whereas the tangent at the shortest radius side is directed towards the corresponding polycapillary focal distance.
    '''
    ctypedef enum polycap_profile_type:
        POLYCAP_PROFILE_CONICAL
        POLYCAP_PROFILE_PARABOLOIDAL
        POLYCAP_PROFILE_ELLIPSOIDAL

    '''Create a new profile for a given profile type with supplied polycapillary properties
    @param type : an integer or type that indicates the profile type
    @param length : the polycapillary length, as measured along the central axis [cm]
    @param rad_ext_upstream : external upstream polycapillary radius (photons stream from upstream to downstream) [cm]
    @param rad_ext_downstream : external downstream polycapillary radius (photons stream from upstream to downstream) [cm]
    @param rad_int_upstream : internal upstream capillary radius (photons stream from upstream to downstream) [cm]
    @param rad_int_downstream : internal downstream capillary radius (photons stream from upstream to downstream) [cm]
    @param focal_dist_upstream : focal distance upstream of the polycapillary optic (photons stream from upstream to downstream) [cm]
    @param focal_dist_downstream : focal distance downstream of the polycapillary optic (photons stream from upstream to downstream) [cm]
    @param error : a pointer to a \c NULL polycap_error, or \c NULL
    @return : a new polycap_profile, or \c NULL if an error occurred
    '''
    polycap_profile* polycap_profile_new(
	polycap_profile_type type,
	double length,
	double rad_ext_upstream,
	double rad_ext_downstream,
	double rad_int_upstream,
	double rad_int_downstream,
	double focal_dist_upstream,
	double focal_dist_downstream,
	polycap_error **error)

    '''Free the polycap_profile structure and its associated data
    @param profile : a polycap_profile
    '''
    void polycap_profile_free(polycap_profile *profile)
