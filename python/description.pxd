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
from profile cimport polycap_profile
from libc.stdint cimport int64_t

cdef extern from "polycap-description.h" nogil:
    ctypedef struct polycap_description

    polycap_description* polycap_description_new(
        polycap_profile *profile,
        double sig_rough,
        int64_t n_cap,
        unsigned int nelem,
        int iz[],
        double wi[],
        double density,
	polycap_error **error)

    const polycap_profile* polycap_description_get_profile(polycap_description *description)

    void polycap_description_free(polycap_description *description)
