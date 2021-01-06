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
from libc.stdint cimport int64_t

cdef extern from "stdbool.h" nogil:
    ctypedef bint bool


from photon cimport polycap_vector3, polycap_leak
#cdef extern from "polycap-photon.h" nogil:
#    ctypedef struct polycap_vector3:
#        double x
#        double y
#        double z
# 
#    ctypedef struct polycap_leak:
#        polycap_vector3 coords
#        polycap_vector3 direction
#        polycap_vector3 elecv
#        size_t n_energies

cdef extern from "polycap-transmission-efficiencies.h" nogil:
    ctypedef struct polycap_transmission_efficiencies

    void polycap_transmission_efficiencies_free(polycap_transmission_efficiencies *efficiencies)

    bool polycap_transmission_efficiencies_write_hdf5(polycap_transmission_efficiencies *efficiencies, const char *filename, polycap_error **error)

    bool polycap_transmission_efficiencies_get_data(polycap_transmission_efficiencies *efficiencies, size_t *n_energies, double **energies_arr, double **efficiencies_arr, polycap_error **error)

    bool polycap_transmission_efficiencies_get_extleak_data(polycap_transmission_efficiencies *efficiencies, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error)

    bool polycap_transmission_efficiencies_get_intleak_data(polycap_transmission_efficiencies *efficiencies, polycap_leak ***leaks, int64_t *n_leaks, polycap_error **error)

    bool polycap_transmission_efficiencies_get_start_data(polycap_transmission_efficiencies *efficiencies, int64_t *n_start, int64_t *n_exit, polycap_vector3 **start_coords, polycap_vector3 **start_direction, polycap_vector3 **start_elecv, polycap_vector3 **src_start_coords, polycap_error **error)

    bool polycap_transmission_efficiencies_get_exit_data(polycap_transmission_efficiencies *efficiencies, int64_t *n_exit, polycap_vector3 **exit_coords, polycap_vector3 **exit_direction, polycap_vector3 **exit_elecv, int64_t **n_refl, double **d_travel, size_t *n_energies, double *** exit_weights, polycap_error **error)
