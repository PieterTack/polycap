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

cdef extern from "polycap-error.h" nogil:
    cdef enum polycap_error_code:
        POLYCAP_ERROR_MEMORY             # set in case of a memory allocation problem 
        POLYCAP_ERROR_INVALID_ARGUMENT   # set in case an invalid argument gets passed to a routine
        POLYCAP_ERROR_IO                 # set in case an error involving input/output occurred
        POLYCAP_ERROR_OPENMP             # set in case an error involving OpenMP occurred
        POLYCAP_ERROR_TYPE               # set in case an error involving type conversion occurred (HDF5 related) */
        POLYCAP_ERROR_UNSUPPORTED        # set in case an unsupported feature has been requested */
        POLYCAP_ERROR_RUNTIME            # set in case an unexpected runtime error occurred */

    ctypedef struct polycap_error:
        polycap_error_code code
        char *message

    void polycap_error_free(polycap_error *error)
