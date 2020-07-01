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

cdef extern from "stdbool.h" nogil:
    ctypedef bint bool

cdef extern from "polycap-transmission-efficiencies.h" nogil:
    ctypedef struct polycap_transmission_efficiencies

    void polycap_transmission_efficiencies_free(polycap_transmission_efficiencies *efficiencies)

    bool polycap_transmission_efficiencies_write_hdf5(polycap_transmission_efficiencies *efficiencies, const char *filename, polycap_error **error)

    bool polycap_transmission_efficiencies_get_data(polycap_transmission_efficiencies *efficiencies, size_t *n_energies, double **energies_arr, double **efficiencies_arr, polycap_error **error)
