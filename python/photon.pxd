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
from libc.stdint cimport int64_t

cdef extern from "stdbool.h" nogil:
    ctypedef bint bool

cdef extern from "polycap-photon.h" nogil:

    ctypedef struct polycap_vector3:
        double x
        double y
        double z

    ctypedef struct polycap_photon

    ctypedef struct polycap_leak:
        polycap_vector3 coords
        polycap_vector3 direction
        polycap_vector3 elecv
        size_t n_energies;
        double *weight
        int64_t n_refl

    polycap_photon* polycap_photon_new(
        polycap_description *description,
        polycap_vector3 start_coords,
        polycap_vector3 start_direction,
        polycap_vector3 start_electric_vector,
        polycap_error **error)

    int polycap_photon_launch(polycap_photon *photon, size_t n_energies, double *energies, double **weights, bint leak_calc, polycap_error **error)

    polycap_vector3 polycap_photon_get_start_coords(polycap_photon *photon)

    polycap_vector3 polycap_photon_get_start_direction(polycap_photon *photon)

    polycap_vector3 polycap_photon_get_start_electric_vector(polycap_photon *photon)

    polycap_vector3 polycap_photon_get_exit_coords(polycap_photon *photon)

    polycap_vector3 polycap_photon_get_exit_direction(polycap_photon *photon)

    polycap_vector3 polycap_photon_get_exit_electric_vector(polycap_photon *photon)

    double polycap_photon_get_dtravel(polycap_photon *photon)

    int64_t polycap_photon_get_irefl(polycap_photon *photon)

    bool polycap_photon_get_extleak_data(polycap_photon *photon, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error)

    bool polycap_photon_get_intleak_data(polycap_photon *photon, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error)

    void polycap_photon_free(polycap_photon *photon)

    void polycap_leak_free(polycap_leak *leak)
