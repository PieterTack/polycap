if ENABLE_PYTHON
pyexec_LTLIBRARIES = polycap.la 
polycap_la_CFLAGS = $(PYTHON_CFLAGS) -I$(top_srcdir)/include $(PYTHON_CPPFLAGS) $(NUMPY_HEADERS)
polycap_la_LIBADD = ../src/libpolycap.la
nodist_polycap_la_SOURCES = polycap-python.c
polycap_la_LDFLAGS = $(LDFLAGS_PYTHON) -avoid-version -module -shared -export-dynamic -shrext $(PYTHON_EXT)
#if OS_WINDOWS
#polycap_la_LDFLAGS += $(PYTHON_LDFLAGS)
#endif
endif

cython_sources = cerror.pxd crng.pxd polycap.pyx

polycap-python.c: $(cython_sources)
	$(CYTHON) -X boundscheck=False,wraparound=False,cdivision=True -I${top_srcdir}/include -o ${builddir}/polycap-python.c ${srcdir}/polycap.pyx

clean-local:
	rm -rf polycap-python.c

distclean-local:
	rm -rf __pycache__

EXTRA_DIST = $(cython_sources)
