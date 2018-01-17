bin_PROGRAMS = polycap
polycap_SOURCES = main.c polycap.c polycap.h polycap-private.h
polycap_LDADD = @easyRNG_LIBS@ @gsl_LIBS@ @xraylib_LIBS@ -lm
polycap_CPPFLAGS = @easyRNG_CFLAGS@ @gsl_CFLAGS@ @xraylib_CFLAGS@ -I$(srcdir)
polycap_CFLAGS = @OPENMP_CFLAGS@
polycap_LDFLAGS = @OPENMP_CFLAGS@

