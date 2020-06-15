
NULL=

lib_LTLIBRARIES = libpolycap.la
libpolycap_la_SOURCES = \
	polycap-profile.c \
	polycap-description.c \
	polycap-source.c \
	polycap-photon.c \
	polycap-capil.c \
	polycap-transmission-efficiencies.c \
	polycap-private.h \
	polycap-rng.c \
	polycap-error.c \
	polycap-aux.c \
	polycap-aux.h \
	$(NULL)
libpolycap_la_CPPFLAGS = @easyRNG_CFLAGS@ @gsl_CFLAGS@ @xraylib_CFLAGS@ -I$(top_srcdir)/include @HDF5_CFLAGS@
libpolycap_la_LIBADD = @easyRNG_LIBS@ @gsl_LIBS@ @xraylib_LIBS@ @HDF5_LIBS@ -lm
libpolycap_la_CFLAGS = @OPENMP_CFLAGS@ -Wno-error=attributes $(HIDDEN_VISIBILITY_CFLAGS)
libpolycap_la_LDFLAGS = @OPENMP_CFLAGS@ @LDFLAGS_LIBPOLYCAP@

check_LTLIBRARIES = libpolycap-check.la
libpolycap_check_la_SOURCES = $(libpolycap_la_SOURCES)
libpolycap_check_la_CPPFLAGS = $(libpolycap_la_CPPFLAGS) -DTEST_BUILD
libpolycap_check_la_LIBADD = $(libpolycap_la_LIBADD)
libpolycap_check_la_CFLAGS = $(libpolycap_la_CFLAGS)
libpolycap_check_la_LDFLAGS = $(libpolycap_la_LDFLAGS)

bin_PROGRAMS = polycap
polycap_SOURCES = main.c 
polycap_CFLAGS = @OPENMP_CFLAGS@ -Wno-error=attributes
polycap_CPPFLAGS = @easyRNG_CFLAGS@ @gsl_CFLAGS@ -I$(srcdir) -I$(top_srcdir)/include
polycap_LDADD = libpolycap.la
polycap_LDFLAGS = @OPENMP_CFLAGS@


