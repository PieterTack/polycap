AM_CPPFLAGS = -I${top_srcdir}/src -I$(top_srcdir)/include -DEXAMPLE_DIR=\"$(top_srcdir)/example/\"

check_PROGRAMS = error profile description capil photon source
check_SCRIPTS =
if ENABLE_PYTHON
check_SCRIPTS += python.sh
endif

TESTS = $(check_PROGRAMS) $(check_SCRIPTS)

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

source_SOURCES = source.c
source_LDADD = ../src/libpolycap.la
source_CFLAGS = @OPENMP_CFLAGS@ @easyRNG_CFLAGS@
source_LDFLAGS = @OPENMP_CFLAGS@

python.sh: ../python/polycap.la python.py
	@echo "LD_LIBRARY_PATH=\"${top_builddir}/src/.libs\" DYLD_LIBRARY_PATH=\"${top_builddir}/src/.libs\" PYTHONPATH=\"${top_builddir}/python:${top_builddir}/python/.libs\" $(PYTHON) ${top_srcdir}/tests/python.py" > python.sh
	@chmod +x python.sh

EXTRA_DIST = python.py

clean-local:
	rm -rf python.sh
