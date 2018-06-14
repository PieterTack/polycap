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
from description cimport polycap_description
from rng cimport polycap_rng

cdef extern from "polycap-photon.h" nogil:

    ctypedef struct polycap_vector3:
        double x
        double y
        double z

    ctypedef struct polycap_photon

    polycap_photon* polycap_photon_new(
        polycap_description *description,
        polycap_rng *rng,
        polycap_vector3 start_coords,
        polycap_vector3 start_direction,
        polycap_vector3 start_electric_vector,
        size_t n_energies,
        double *energies,
        polycap_error **error)

    int polycap_photon_launch(polycap_photon *photon, polycap_error **error)

    polycap_vector3 polycap_photon_get_exit_coords(polycap_photon *photon)

    polycap_vector3 polycap_photon_get_exit_direction(polycap_photon *photon)

    polycap_vector3 polycap_photon_get_exit_electric_vector(polycap_photon *photon)

    void polycap_photon_free(polycap_photon *photon)
