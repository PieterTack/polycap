# /*
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
# */

from error cimport polycap_error
from photon cimport polycap_photon
from description cimport polycap_description
from rng cimport polycap_rng
from transmission_efficiencies cimport polycap_transmission_efficiencies
from progress_monitor cimport polycap_progress_monitor

cdef extern from "polycap-source.h" nogil:
    ctypedef struct polycap_source

    polycap_source* polycap_source_new(
        polycap_description *description,
        double d_source,
        double src_x,
        double src_y,
        double src_sigx,
        double src_sigy,
        double src_shiftx,
        double src_shifty,
        polycap_error **error)

    void polycap_source_free(polycap_source *source)

    polycap_photon* polycap_source_get_photon(
        polycap_source *source,
        polycap_rng *rng,
        size_t n_energies,
        double *energies,
        polycap_error **error)

    polycap_source* polycap_source_new_from_file(const char *filename, polycap_error **error)

    polycap_transmission_efficiencies* polycap_source_get_transmission_efficiencies(
        polycap_source *source,
        int max_threads,
        size_t n_energies,
        double *energies,
        int n_photons,
        polycap_progress_monitor *progress_monitor,
        polycap_error **error)

    const polycap_description* polycap_source_get_description(polycap_source *source)

