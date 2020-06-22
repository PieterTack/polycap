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
from photon cimport polycap_photon
from description cimport polycap_description
from rng cimport polycap_rng
from transmission_efficiencies cimport polycap_transmission_efficiencies
from progress_monitor cimport polycap_progress_monitor

cdef extern from "polycap-source.h" nogil:
    """Struct containing information on the source from which photons can be (randomly) selected
    When this struct is no longer required, it is the user's responsability to free the memory using polycap_source_free().
    """
    ctypedef struct polycap_source

    """Creates a new polycap_source by providing all its properties
    @param description : a polycap_description
    @param d_source : the distance between the source and polycapillary optic entrance window along the central axis [cm]
    @param src_x : the source radius along the X (horizontal) direction [cm]
    @param src_y : the source radius along the y (vertical) direction [cm]
    @param src_sigx : the maximal divergence of photons along the X (horizontal) direction [rad]. Negative values in src_sigx or src_sigy represent homogeneous polycapillary optic illumination.
    @param src_sigy : the maximal divergence of photons along the Y (vertical) direction [rad]. Negative values in src_sigx or src_sigy represent homogeneous polycapillary optic illumination.
    @param src_shiftx : lateral shift of the source centre along the X (horizontal) direction with respect to the polycapillary optic central axis [cm]
    @param src_shifty : lateral shift of the source centre along the Y (vertical) direction with respect to the polycapillary optic central axis [cm]
    @param hor_pol : the polarisation factor of the simulated source (-1 <= hor_pol <= 1)
    @param n_energies : the amount of discrete energies of which the source will emit photons
    @param energies : an array containing the discrete energies of which the source will emit photons
    @param error : a pointer to a \c NULL polycap_error, or \c NULL
    @return : a new polycap_source, or \c NULL if an error occurred
    """
    polycap_source* polycap_source_new(
        polycap_description *description,
        double d_source,
        double src_x,
        double src_y,
        double src_sigx,
        double src_sigy,
        double src_shiftx,
        double src_shifty,
	double hor_pol,
	size_t n_energies,
	double *energies,
        polycap_error **error)

    """free a polycap_source struct
    @param source : a polycap_source instance
    """
    void polycap_source_free(polycap_source *source)

    """Create a new random polycap_photon based on polycap_source
    In the event of an error, \c NULL is returned and \c error is set appropriately.
    @param source : a polycap_source
    @param rng : a polycap_rng
    @param error : a pointer to a \c NULL polycap_error, or \c NULL
    @return : a new polycap_photon, or \c NULL if an error occurred
    """
    polycap_photon* polycap_source_get_photon(
        polycap_source *source,
        polycap_rng *rng,
        polycap_error **error)

    """Load a polycap_description from given ASCII *.inp input file correponding to the old polycap program format.
    @param filename : directory path to an ASCII input file. Default extension *.inp.
    @param error : a pointer to a \c NULL polycap_error, or \c NULL
    @return : a new polycap_source, or \c NULL if an error occurred
    """
    polycap_source* polycap_source_new_from_file(const char *filename, polycap_error **error)

    """Obtain the transmission efficiencies for a given array of energies, and a full polycap_description.
    Efficiencies are allocated by this function, and need to be freed with polycap_transmission_efficiencies_free().
    @param source : a polycap_source
    @param max_threads : the amount of threads to use. Set to -1 to use the maximum available amount of threads.
    @param n_photons : the amount of photons to simulate that reach the polycapillary end
    @param progress_monitor : a polycap_progress_monitor
    @param leak_calc : True: perform leak calculation; False: do not perform leak calculation
    @param error : a pointer to a \c NULL polycap_error, or \c NULL
    @return : a new polycap_transmission_efficiencies, or \c NULL if an error occurred
    """
    polycap_transmission_efficiencies* polycap_source_get_transmission_efficiencies(
        polycap_source *source,
        int max_threads,
        int n_photons,
	bint leak_calc,
        polycap_progress_monitor *progress_monitor,
        polycap_error **error)

    """Create new polycap_description from a polycap_source
    @param source : a polycap_source
    @return : a new polycap_description
    """
    const polycap_description* polycap_source_get_description(polycap_source *source)

