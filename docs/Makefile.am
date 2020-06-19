if HAVE_DOXYGEN

NULL =
DEPS = \
       $(top_srcdir)/include/polycap.h \
       $(top_srcdir)/include/polycap-description.h \
       $(top_srcdir)/include/polycap-error.h \
       $(top_srcdir)/include/polycap-profile.h \
       $(top_srcdir)/include/polycap-photon.h \
       $(top_srcdir)/include/polycap-progress-monitor.h \
       $(top_srcdir)/include/polycap-rng.h \
       $(top_srcdir)/include/polycap-source.h \
       $(top_srcdir)/include/polycap-transmission-efficiencies.h \
       $(EXTRA_DIST)

doxyfile.stamp: $(DEPS)
	$(DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

polycap_docdir = $(datadir)/doc/polycap/html


install-data-local:
	$(MKDIR_P) $(DESTDIR)$(polycap_docdir)
	$(MKDIR_P) $(DESTDIR)$(polycap_docdir)/search
	$(INSTALL_DATA) $(builddir)/html/*.css $(DESTDIR)$(polycap_docdir)
	$(INSTALL_DATA) $(builddir)/html/*.html $(DESTDIR)$(polycap_docdir)
	$(INSTALL_DATA) $(builddir)/html/*.js $(DESTDIR)$(polycap_docdir)
	$(INSTALL_DATA) $(builddir)/html/*.png $(DESTDIR)$(polycap_docdir)
	$(INSTALL_DATA) $(builddir)/html/search/*.css $(DESTDIR)$(polycap_docdir)/search
	$(INSTALL_DATA) $(builddir)/html/search/*.html $(DESTDIR)$(polycap_docdir)/search
	$(INSTALL_DATA) $(builddir)/html/search/*.js $(DESTDIR)$(polycap_docdir)/search
	$(INSTALL_DATA) $(builddir)/html/search/*.png $(DESTDIR)$(polycap_docdir)/search


CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp
clean-local:
	rm -rf html
uninstall-local:
	rm -rf $(DESTDIR)$(polycap_docdir)

endif

EXTRA_DIST = DoxygenLayout.xml extra_pages.dox meson.build

