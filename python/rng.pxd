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

cdef extern from "polycap-rng.h" nogil:
    """Struct containing a rng
    The polycap_rng struct is  mapped to either gsl_rng or easy_rng. When this struct is no longer required, it is the user's responsability to free the memory using polycap_rng_free().
    """
    ctypedef struct polycap_rng

    """Create a new rng with seed from `/dev/urandom` or `rand_s`
    @return : a new polycap_rng
    """
    polycap_rng* polycap_rng_new()

    """get a new rng with seed provided by caller
    @param seed : a seed provided by the caller
    @return : a new polycap_rng
    """
    polycap_rng* polycap_rng_new_with_seed(unsigned long int seed)

    """free a polycap_rng structure
    @param rng : a polycap_rng
    """
    void polycap_rng_free(polycap_rng *rng)
