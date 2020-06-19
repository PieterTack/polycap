if ENABLE_PYTHON
pyexec_LTLIBRARIES = polycap.la 
polycap_la_CFLAGS = $(PYTHON_CFLAGS) $(CYTHON_ERROR_CFLAGS)
polycap_la_CPPFLAGS = @xraylib_CFLAGS@ -I$(top_srcdir)/include $(PYTHON_CPPFLAGS) $(NUMPY_HEADERS)
polycap_la_LIBADD = ../src/libpolycap.la
nodist_polycap_la_SOURCES = polycap-python.c
polycap_la_LDFLAGS = $(LDFLAGS_PYTHON) -avoid-version -module -shared -export-dynamic -shrext $(PYTHON_EXT)
if OS_WINDOWS
polycap_la_LDFLAGS += $(PYTHON_LDFLAGS)
endif
endif

NULL=

cython_sources = \
	error.pxd \
	rng.pxd \
	profile.pxd \
	transmission_efficiencies.pxd \
	description.pxd \
	photon.pxd \
	source.pxd \
	progress_monitor.pxd \
	polycap.pyx \
	$(NULL)

polycap-python.c: $(cython_sources)
	$(CYTHON) -X language_level=3,boundscheck=False,wraparound=False,cdivision=True -I${top_srcdir}/include -o ${builddir}/polycap-python.c ${srcdir}/polycap.pyx

clean-local:
	rm -rf polycap-python.c

distclean-local:
	rm -rf __pycache__

EXTRA_DIST = $(cython_sources) meson.build
