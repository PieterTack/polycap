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

#if __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ > 4)
#define GNUC_PRINTF( format_idx, arg_idx )    \
  __attribute__((__format__ (__printf__, format_idx, arg_idx)))
#else /* !__GNUC__ */
#define GNUC_PRINTF( format_idx, arg_idx )
#endif /* !__GNUC__ */

enum polycap_error_code {
	POLYCAP_ERROR_MEMORY, /* set in case of a memory allocation problem */
	POLYCAP_ERROR_INVALID_ARGUMENT, /* set in case an invalid argument gets passed to a routine */
	POLYCAP_ERROR_IO, /* set in case an error involving input/output occurred */
	POLYCAP_ERROR_OPENMP, /* set in case an error involving OpenMP occurred */
	POLYCAP_ERROR_TYPE, /* set in case an error involving type conversion occurred (HDF5 related) */
	POLYCAP_ERROR_UNSUPPORTED, /* set in case an unsupported feature has been requested */
	POLYCAP_ERROR_RUNTIME, /* set in case an unexpected runtime error occurred */
};


/**
 * polycap_error:
 * @code: error code, e.g. %POLYCAP_ERROR_MEMORY
 * @message: human-readable informative error message
 *
 * The `polycap_error` structure contains information about
 * an error that has occurred.
 */
typedef struct _polycap_error polycap_error;

struct _polycap_error
{
  enum polycap_error_code code;
  char *message;
};

polycap_error* polycap_error_new(enum polycap_error_code code, const char *format, ...) GNUC_PRINTF (2, 3);

polycap_error* polycap_error_new_literal(enum polycap_error_code code, const char *message);

polycap_error* polycap_error_new_valist(enum polycap_error_code code, const char *format, va_list args) GNUC_PRINTF(2, 0);

void polycap_error_free(polycap_error *error);

polycap_error* polycap_error_copy(const polycap_error *error);

bool polycap_error_matches(const polycap_error *error, enum polycap_error_code code);

void polycap_set_error(polycap_error **err, enum polycap_error_code code , const char *format, ...) GNUC_PRINTF (3, 4);

void polycap_set_error_literal(polycap_error **err, enum polycap_error_code code, const char *message);

void polycap_propagate_error(polycap_error **dest, polycap_error *src);

void polycap_clear_error(polycap_error **err);

#ifdef __cplusplus
}
#endif

#endif
