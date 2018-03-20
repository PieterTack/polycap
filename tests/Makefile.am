AM_CPPFLAGS = -I${top_srcdir}/src -I$(top_srcdir)/include

check_PROGRAMS = error
TESTS = $(check_PROGRAMS)

error_SOURCES = error.c
error_LDADD = ../src/libpolycap.la
error_CFLAGS = @OPENMP_CFLAGS@
error_LDFLAGS = @OPENMP_CFLAGS@
