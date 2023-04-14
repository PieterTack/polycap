project('polycap', 'c',
  version: '1.2',
  license: 'GPL-3.0',
  meson_version: '>= 0.60.0',
  default_options: ['b_asneeded=false', # circumvent GSL trouble with cblas
    'c_std=c99'],
)

host_system = host_machine.system()

polycap_version_array = meson.project_version().split('.')
polycap_major_version = polycap_version_array[0].to_int()
polycap_minor_version = polycap_version_array[1].to_int()

# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
# The relation between libtool's current:revison:age interface versioning
# and the .so filename, .so.x.y.z, is
# x = current - age
# y = age
# z = revision
# If libtool_soversion is updated as described in libtool's documentation,
# x.y.z will usually *not* be equal to meson.project_version().
libtool_soversion = [3, 0, 1]
polycap_libversion = '@0@.@1@.@2@'.format(
  libtool_soversion[0] - libtool_soversion[2],
  libtool_soversion[2],
  libtool_soversion[1])

lib_current = libtool_soversion[0] + 1
lib_revision =  libtool_soversion[1]
darwin_versions = [lib_current, '@0@.@1@'.format(lib_current , lib_revision)]

# Use these instead of meson.source_root() and meson.build_root() in subdirectories.
# source_root() and build_root() are not useful, if this is a subproject.
project_source_root = meson.current_source_dir()
project_build_root = meson.current_build_dir()

cc = meson.get_compiler('c')

build_documentation_opt = get_option('build-documentation')
build_python_opt = get_option('build-python-bindings')

# Configuration data
config_h_data = configuration_data()
config_h_data.set_quoted('PACKAGE_TARNAME', meson.project_name())
config_h_data.set_quoted('PACKAGE_VERSION', meson.project_version())
config_h_data.set_quoted('VERSION', meson.project_version())

funcs = [
  'strdup',
  '_strdup',
  '_ftime',
  'gettimeofday',
  '_unlink',
  'unlink',
]

foreach f : funcs
  if cc.has_function(f)
    define = 'HAVE_' + f.underscorify().to_upper()
    config_h_data.set(define, 1)
  endif
endforeach

special_funcs = []

if host_system == 'windows'
  special_funcs += ['_vscprintf', '_scprintf', '_vsnprintf', '_snprintf']
else
  special_funcs += ['vasprintf', 'asprintf']
endif

legacy_stdio_dep =  cc.find_library('legacy_stdio_definitions', required: false)

foreach f : special_funcs
  if not cc.has_function(f, args : '-D_GNU_SOURCE', dependencies: legacy_stdio_dep)
    error('Required function ' + f + ' not found on this platform')
  endif
endforeach

if cc.get_id() != 'msvc' and cc.get_id() != 'clang-cl'
  config_h_data.set('HAVE_PROPER_COMPLEX_H', true)
endif

# Detect and set symbol visibility
if get_option('default_library') != 'static'
  if host_system == 'windows' or host_system == 'cygwin'
    config_h_data.set('WIN32', true)
    config_h_data.set('GSL_DLL', true)
    config_h_data.set('ZLIB_DLL', true)
    config_h_data.set('POLYCAP_DLL_EXPORT', true)
    if cc.get_id() == 'msvc' or cc.get_id() == 'clang-cl'
      config_h_data.set('POLYCAP_EXTERN', '__declspec(dllexport) extern')
    elif cc.has_argument('-fvisibility=hidden')
      config_h_data.set('POLYCAP_EXTERN', '__attribute__((visibility("default"))) __declspec(dllexport) extern')
    endif
  elif cc.has_argument('-fvisibility=hidden')
    config_h_data.set('POLYCAP_EXTERN', '__attribute__((visibility("default"))) extern')
  endif
endif

# Dependencies
polycap_pkg_config_requires = []
polycap_pkg_config_requires_private = []


xraylib_dep = dependency('libxrl', version: '>=4.0.0')
gsl_dep = dependency('gsl')
hdf5_dep = dependency('hdf5', language: 'c')
openmp_dep = dependency('openmp', language: 'c')
m_dep = cc.find_library('m', required : false)
thread_dep = dependency('threads')
polycap_build_dep = [xraylib_dep, gsl_dep, hdf5_dep, openmp_dep, m_dep, thread_dep]
if host_system == 'windows'
  zlib_dep = dependency('zlib')
  polycap_build_dep += zlib_dep # this is to avoild linking trouble with Conda, which is currently having trouble with processing paths in pkg-config files
endif
polycap_pkg_config_requires_private = [xraylib_dep, gsl_dep,]

easyRNG_dep = dependency('easyRNG', version: '>= 1.1', required: false)
if easyRNG_dep.found()
  config_h_data.set('HAVE_EASYRNG', true)
  polycap_build_dep += easyRNG_dep
  polycap_pkg_config_requires_private += easyRNG_dep
endif

configure_file(output : 'config.h', configuration : config_h_data)

subdir('src')
subdir('include')

if build_documentation_opt
    doxygen = find_program('doxygen')
    perl = find_program('perl')
    subdir('docs')
endif

pkgconfig = import('pkgconfig')

pkgconfig.generate(
  polycap_lib,
  filebase: 'libpolycap',
  name: 'polycap',
  description: 'A library to calculate X-ray transmission efficiency profiles for polycapillaries',
  version: meson.project_version(),
  libraries_private: [hdf5_dep, openmp_dep, m_dep],
  requires: polycap_pkg_config_requires,
  requires_private: polycap_pkg_config_requires_private,
  subdirs: 'polycap'
)


if build_python_opt
  pymod = import('python')
  python = pymod.find_installation(get_option('python'))

  if python.language_version().version_compare('< 3.6')
    error('Requires Python >= 3.6')
  endif

  python_dep = python.dependency()
  add_languages('cython', 'c', native: false)
  subdir('python')
endif

subdir('tests')
