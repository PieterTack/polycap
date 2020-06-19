libpolycap_sources = files(
  'polycap-profile.c',
  'polycap-description.c',
  'polycap-source.c',
  'polycap-photon.c',
  'polycap-capil.c',
  'polycap-transmission-efficiencies.c',
  'polycap-private.h',
  'polycap-rng.c',
  'polycap-error.c',
  'polycap-aux.c',
  'polycap-aux.h',
)

core_c_args = [
  '-DHAVE_CONFIG_H',
  '-D_GNU_SOURCE',
]

libpolycap_error_flags = [
  '-Werror=implicit-function-declaration',
  '-Werror=deprecated',
  '-Werror=deprecated-declarations',
  '-Werror=implicit',
  '-Werror=unused-function',
  '-Werror=parentheses',
  '-Werror=unused-result',
]

libpolycap_error_flags = cc.get_supported_arguments(libpolycap_error_flags)


extra_include_dirs = include_directories('..', '.', '../include',)

polycap_lib = shared_library(
  'polycap',
  libpolycap_sources,
  version: polycap_libversion,
  darwin_versions: darwin_versions,
  dependencies: polycap_build_dep,
  install: true,
  c_args: core_c_args + libpolycap_error_flags,
  gnu_symbol_visibility: 'hidden',
  include_directories: extra_include_dirs,
  )

polycap_lib_dep = declare_dependency(
  link_with: polycap_lib,
  include_directories: extra_include_dirs,
  dependencies: polycap_build_dep,
  )

polycap_check_lib = shared_library(
  'polycap-check',
  libpolycap_sources,
  dependencies: polycap_build_dep,
  install: false,
  c_args: core_c_args + ['-DTEST_BUILD'] + libpolycap_error_flags,
  include_directories: extra_include_dirs,
  )

polycap_check_lib_dep = declare_dependency(
  link_with: polycap_check_lib,
  include_directories: extra_include_dirs,
  dependencies: polycap_build_dep,
  )

executable(
  'polycap',
  files('main.c'),
  dependencies: polycap_lib_dep,
  install: true,
  c_args: core_c_args + libpolycap_error_flags,
  )
