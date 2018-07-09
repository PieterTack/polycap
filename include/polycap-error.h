/*
 * Copyright (C) 2018 Pieter Tack, Tom Schoonjans and Laszlo Vincze
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */

/** \file polycap-error.h
 *  \brief API for dealing with errors in polycap
 *
 *  This header contains all functions and definitions that are necessary to create, manipulate and free the error objects that are produced by polycap.
 *  In the vast majority of cases, the user will not have to instantiate errors, but instead will obtain an error object by passing a reference to a `polycap_error` pointer that was initialized to \c NULL, to a function that supports it. Consider the following example that will demonstrate this feature.
 *
 *  \code{.c}
 *  polycap_profile *profile;
 *  polycap_error *error = NULL;
 *
 *  profile = polycap_profile_new(POLYCAP_PROFILE_CONICAL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &error);
 *
 *  if (err != NULL) {
 *      // in this case profile will be NULL!
 *      fprintf(stderr, "Some error occurred: %s", error->message);
 *      // deal with error...
 *      ...
 *      // free error
 *      polycap_clear_error(&error);
 *  } 
 *  \endcode
 *
 *  \note the polycap_error API is strongly based on [GError](https://developer.gnome.org/glib/stable/glib-Error-Reporting.html), as implemented in Glib. The user is recommended to have a look at the GError documentation for more information.
 *  \note there is no equivalent for this header in the Python API since whenever a function produces an error, the corresponding \c polycap_error instance will be converted into an appropriate exception, which will subsequently be raised.
 */

#ifndef POLYCAP_ERROR_H
#define POLYCAP_ERROR_H

#include <stdbool.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 *  This file is mostly copy-pasted from GLib's polycap_error methods...
 */ 
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ > 4)
#define GNUC_PRINTF( format_idx, arg_idx )    \
  __attribute__((__format__ (__printf__, format_idx, arg_idx)))
#else /* !__GNUC__ */
#define GNUC_PRINTF( format_idx, arg_idx )
#endif /* !__GNUC__ */
#endif

/** Codes to indicate the type of the error
 *
 */
enum polycap_error_code {
	POLYCAP_ERROR_MEMORY, ///< set in case of a memory allocation problem
	POLYCAP_ERROR_INVALID_ARGUMENT, ///< set in case an invalid argument gets passed to a routine
	POLYCAP_ERROR_IO, ///< set in case an error involving input/output occurred
	POLYCAP_ERROR_OPENMP, ///< set in case an error involving OpenMP occurred
	POLYCAP_ERROR_TYPE, ///< set in case an error involving type conversion occurred (HDF5 related)
	POLYCAP_ERROR_UNSUPPORTED, ///< set in case an unsupported feature has been requested
	POLYCAP_ERROR_RUNTIME, ///< set in case an unexpected runtime error occurred
};


/** Struct containing information about an error.
 * 
 * Typically the user will not have to deal with allocating and populating these structs,
 * as task will be accomplished by the polycap API. However, when such a struct is no longer required,
 * it is the user's responsability to free the memory using either polycap_error_free() or polycap_clear_error().
 */
typedef struct {
  enum polycap_error_code code; ///< this maps to an integer that will indicate the kind of error was encountered. 
  char *message; ///< a detailed error message
} polycap_error;

/** Creates a new polycap_error with the given code , and a message formatted with format 
 *
 * \param code error code
 * \param format printf()-style format for error message
 * \param ... parameters for message format
 * \returns a new polycap_error
 */
polycap_error* polycap_error_new(enum polycap_error_code code, const char *format, ...) GNUC_PRINTF (2, 3);

/** Creates a new polycap_error with the given code and message
 *
 * \param code error code
 * \param message an error message
 * \returns a new polycap_error
 */
polycap_error* polycap_error_new_literal(enum polycap_error_code code, const char *message);

/** Creates a new polycap_error with a va_list
 *
 * \param code error code
 * \param format printf()-style format for error message
 * \param args a va_list of arguments
 * \returns a new polycap_error, or \c NULL if an error occurred
 */
polycap_error* polycap_error_new_valist(enum polycap_error_code code, const char *format, va_list args) GNUC_PRINTF(2, 0);

/** Frees a polycap_error
 *
 * \param error a polycap_error
 */
void polycap_error_free(polycap_error *error);

/** Copies a polycap_error
 *
 * \param error a polycap_error
 * \returns a new polycap_error, or \c NULL if a \c NULL polycap_error was supplied
 */
polycap_error* polycap_error_copy(const polycap_error *error);

/** Matches a polycap_error_code to a polycap_error
 *
 * \param error a polycap_error
 * \param code error code
 * \returns true or false
 */
bool polycap_error_matches(const polycap_error *error, enum polycap_error_code code);

/** Sets a polycap_error
 *
 * \param err a pointer to a polycap_error
 * \param code error code
 * \param format printf()-style format for error message
 * \param ... parameters for message format
 */
void polycap_set_error(polycap_error **err, enum polycap_error_code code , const char *format, ...) GNUC_PRINTF (3, 4);

/** Sets an error message to a polycap_error
 *
 * \param err a pointer to a polycap_error
 * \param code error code
 * \param message an error message
 */
void polycap_set_error_literal(polycap_error **err, enum polycap_error_code code, const char *message);

/** Propagate an error
 *
 * \param dest a pointer to a polycap_error
 * \param src a polycap_error
 */
void polycap_propagate_error(polycap_error **dest, polycap_error *src);

/** Clears a polycap_error, and sets the pointer to \c NULL
 *
 * \param err a pointer to a polycap_error
 */
void polycap_clear_error(polycap_error **err);

#ifdef __cplusplus
}
#endif

#endif

