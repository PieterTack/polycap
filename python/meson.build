numpy_include_script = '''
import numpy
print(numpy.get_include())
'''

rv = run_command(get_option('python'), '-c', numpy_include_script, check: false)
if rv.returncode() != 0
  error('Could not retrieve numpy include location. Please check that numpy has been installed.')
endif
numpy_header_location = rv.stdout().strip()

python_error_flags = [
  '-Wno-error=cpp',
  '-Wno-error=attributes',
  '-Wno-error=deprecated-declarations',
  '-Wno-error=unreachable-code',
  '-Wno-error=ignored-optimization-argument',
]

python_error_flags = cc.get_supported_arguments(python_error_flags)

polycap_ext = python.extension_module(
  'polycap',
  'polycap.pyx',
  dependencies : [polycap_lib_dep, python_dep],
  include_directories: extra_include_dirs,
  install: true,
  install_dir: python.get_install_dir(pure: false),
  c_args: core_c_args + ['-I' + numpy_header_location] + libpolycap_error_flags + python_error_flags,
  cython_args: ['-X', 'language_level=3,boundscheck=False,wraparound=False,cdivision=True'],
)

pydir = meson.current_build_dir()
