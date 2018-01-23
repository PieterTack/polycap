polycapincludedir = $(includedir)/polycap
polycapinclude_HEADERS = polycap.h


lib_LTLIBRARIES = libpolycap.la
libpolycap_la_SOURCES = polycap.c polycap-profile.c polycap-private.h $(polycapinclude_HEADERS)
libpolycap_la_CPPFLAGS = @easyRNG_CFLAGS@ @gsl_CFLAGS@ @xraylib_CFLAGS@ -I$(srcdir)
libpolycap_la_LIBADD = @easyRNG_LIBS@ @gsl_LIBS@ @xraylib_LIBS@ -lm 
libpolycap_la_CFLAGS = @OPENMP_CFLAGS@
libpolycap_la_LDFLAGS = @OPENMP_CFLAGS@

bin_PROGRAMS = polycap
polycap_SOURCES = main.c 
polycap_CPPFLAGS = @easyRNG_CFLAGS@ @gsl_CFLAGS@
polycap_LDADD = libpolycap.la


