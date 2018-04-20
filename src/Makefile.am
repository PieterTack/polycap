
lib_LTLIBRARIES = libpolycap.la
libpolycap_la_SOURCES = polycap-profile.c polycap-description.c polycap-source.c polycap-photon.c polycap-capil.c polycap-transmission-efficiencies.c polycap-private.h polycap-rng.c polycap-error.c
libpolycap_la_CPPFLAGS = @easyRNG_CFLAGS@ @gsl_CFLAGS@ @xraylib_CFLAGS@ -I$(top_srcdir)/include @H5CC_CFLAGS@
libpolycap_la_LIBADD = @easyRNG_LIBS@ @gsl_LIBS@ @xraylib_LIBS@ @H5CC_LIBS@ -lm
libpolycap_la_CFLAGS = @OPENMP_CFLAGS@
libpolycap_la_LDFLAGS = @OPENMP_CFLAGS@

bin_PROGRAMS = polycap
polycap_SOURCES = main.c 
polycap_CFLAGS = @OPENMP_CFLAGS@
polycap_CPPFLAGS = @easyRNG_CFLAGS@ @gsl_CFLAGS@ -I$(srcdir) -I$(top_srcdir)/include
polycap_LDADD = libpolycap.la
polycap_LDFLAGS = @OPENMP_CFLAGS@


