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
    ctypedef struct polycap_profile

    ctypedef enum polycap_profile_type:
        POLYCAP_PROFILE_CONICAL
        POLYCAP_PROFILE_PARABOLOIDAL
        POLYCAP_PROFILE_ELLIPSOIDAL

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

    void polycap_profile_free(polycap_profile *profile)
