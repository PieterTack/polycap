
lib_LTLIBRARIES = libpolycap.la
libpolycap_la_SOURCES = polycap.c polycap-profile.c polycap-private.h polycap-old.h polycap-rng.c
libpolycap_la_CPPFLAGS = @easyRNG_CFLAGS@ @gsl_CFLAGS@ @xraylib_CFLAGS@ -I$(top_srcdir)/include
libpolycap_la_LIBADD = @easyRNG_LIBS@ @gsl_LIBS@ @xraylib_LIBS@ -lm 
libpolycap_la_CFLAGS = @OPENMP_CFLAGS@
libpolycap_la_LDFLAGS = @OPENMP_CFLAGS@

bin_PROGRAMS = polycap
polycap_SOURCES = main.c 
polycap_CPPFLAGS = @easyRNG_CFLAGS@ @gsl_CFLAGS@ -I$(srcdir) -I$(top_srcdir)/include
polycap_LDADD = libpolycap.la


