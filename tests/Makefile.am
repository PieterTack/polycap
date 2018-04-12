AM_CPPFLAGS = -I${top_srcdir}/src -I$(top_srcdir)/include -DEXAMPLE_DIR=\"$(top_srcdir)/example/\"

check_PROGRAMS = error profile description capil photon
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

description_SOURCES = description.c
description_LDADD = ../src/libpolycap.la
description_CFLAGS = @OPENMP_CFLAGS@ @easyRNG_CFLAGS@
description_LDFLAGS = @OPENMP_CFLAGS@

capil_SOURCES = capil.c
capil_LDADD = ../src/libpolycap.la
capil_CFLAGS = @OPENMP_CFLAGS@ @easyRNG_CFLAGS@
capil_LDFLAGS = @OPENMP_CFLAGS@

photon_SOURCES = photon.c
photon_LDADD = ../src/libpolycap.la
photon_CFLAGS = @OPENMP_CFLAGS@ @easyRNG_CFLAGS@
photon_LDFLAGS = @OPENMP_CFLAGS@
