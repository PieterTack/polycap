AM_CPPFLAGS = -I${top_srcdir}/src -I$(top_srcdir)/include

check_PROGRAMS = error profile
TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = OMP_CANCELLATION=true

error_SOURCES = error.c
error_LDADD = ../src/libpolycap.la
error_CFLAGS = @OPENMP_CFLAGS@
error_LDFLAGS = @OPENMP_CFLAGS@

profile_SOURCES = profile.c
profile_LDADD = ../src/libpolycap.la
profile_CFLAGS = @OPENMP_CFLAGS@
profile_LDFLAGS = @OPENMP_CFLAGS@
