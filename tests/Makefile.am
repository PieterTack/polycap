AM_CPPFLAGS = -I${top_srcdir}/src -I$(top_srcdir)/include -DEXAMPLE_DIR=\"$(top_srcdir)/example/\" -DTEST_BUILD @easyRNG_CFLAGS@ @gsl_CFLAGS@ @xraylib_CFLAGS@

check_PROGRAMS = version error profile description capil photon source leaks
check_SCRIPTS =
if ENABLE_PYTHON
check_SCRIPTS += python.sh
endif

TESTS = $(check_PROGRAMS) $(check_SCRIPTS)

version_SOURCES = version.c
version_LDADD = ../src/libpolycap-check.la
version_CFLAGS = @OPENMP_CFLAGS@
version_LDFLAGS = @OPENMP_CFLAGS@

error_SOURCES = error.c
error_LDADD = ../src/libpolycap-check.la
error_CFLAGS = @OPENMP_CFLAGS@
error_LDFLAGS = @OPENMP_CFLAGS@

profile_SOURCES = profile.c
profile_LDADD = ../src/libpolycap-check.la
profile_CFLAGS = @OPENMP_CFLAGS@
profile_LDFLAGS = @OPENMP_CFLAGS@

description_SOURCES = description.c
description_LDADD = ../src/libpolycap-check.la
description_CFLAGS = @OPENMP_CFLAGS@ @easyRNG_CFLAGS@
description_LDFLAGS = @OPENMP_CFLAGS@

capil_SOURCES = capil.c
capil_LDADD = ../src/libpolycap-check.la
capil_CFLAGS = @OPENMP_CFLAGS@ @easyRNG_CFLAGS@
capil_LDFLAGS = @OPENMP_CFLAGS@

photon_SOURCES = photon.c
photon_LDADD = ../src/libpolycap-check.la
photon_CFLAGS = @OPENMP_CFLAGS@ @easyRNG_CFLAGS@
photon_LDFLAGS = @OPENMP_CFLAGS@

source_SOURCES = source.c
source_LDADD = ../src/libpolycap-check.la
source_CFLAGS = @OPENMP_CFLAGS@ @easyRNG_CFLAGS@
source_LDFLAGS = @OPENMP_CFLAGS@

leaks_SOURCES = leaks.c
leaks_LDADD = ../src/libpolycap-check.la
leaks_CFLAGS = @OPENMP_CFLAGS@ @easyRNG_CFLAGS@
leaks_LDFLAGS = @OPENMP_CFLAGS@

python.sh: ../python/polycap.la python.py
	@echo "PATH=\"../src/.libs:$$PATH\" LD_LIBRARY_PATH=\"../src/.libs\" DYLD_LIBRARY_PATH=\"../src/.libs\" PYTHONPATH=\"../python/.libs\" $(PYTHON) ${top_srcdir}/tests/python.py" > python.sh
	@chmod +x python.sh

EXTRA_DIST = python.py

clean-local:
	rm -rf python.sh
