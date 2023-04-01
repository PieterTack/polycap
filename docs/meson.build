doxygen_cdata = configuration_data()
doxygen_cdata.set('PACKAGE_VERSION', meson.project_version())
doxygen_cdata.set('top_srcdir', meson.project_source_root())
doxygen_cdata.set('srcdir', meson.current_source_dir())
doxygen_cdata.set('top_builddir', meson.project_build_root())
doxygen_cdata.set('builddir', meson.current_build_dir())

doxyfile = configure_file(input: 'Doxyfile.in',
                          output: 'Doxyfile',
                          configuration: doxygen_cdata)

datadir = join_paths(get_option('datadir'), 'doc', 'polycap')

html_target = custom_target('docs', input: [doxyfile, 'DoxygenLayout.xml', 'extra_pages.dox'] + libpolycap_sources + libpolycap_headers, output: 'html', command: [doxygen, doxyfile], install: true, install_dir: datadir)
