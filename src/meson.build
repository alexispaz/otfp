cc = meson.get_compiler('c')
m = cc.find_library('m', required: false)

tcl = dependency('tcl')

swig = find_program('swig')
swig_gen = generator(
  swig,
  output: ['@BASENAME@_wrap.c'],
  arguments: ['-tcl8', '-o', './@OUTPUT@','@INPUT@'],
  depfile: '@BASENAME@.d',
)

# Check with `pkg-config --libs lapack` if lapack is found.
# Check in PKG_CONFIG_PATH if proper .pc files are avaible.
# Choose library with `meson configure -Dblas=openblas -Dlapack=openblas`
blas_name = get_option('blas')
lapacke_name = get_option('lapacke')
lapack_name = get_option('lapack')
  
blas = dependency(blas_name)
lapack = dependency(lapack_name)

# LAPACKE dependency only if provieder need it
if lapacke_name == 'mkl_lapacke'
    add_project_arguments('-DUSE_MKL_LAPACKE', language: 'c')
    # MKL ships LAPACKE inside the libraries its mkl-*.pc files already link,
    # and no separate mkl_lapacke pkg-config or cmake package (MKL 2026).
    lapacke = blas
else
    lapacke = dependency(lapacke_name, method: 'pkg-config', required: true)
endif

install_headers('calcforces.tcl')
install_headers('cfacv.tcl')
install_headers('chapeau.tcl')
install_headers('cvs.tcl')
install_headers('data.tcl')
install_headers('restraints.tcl')
install_headers('mock.tcl')
install_headers('otfp_engine.h')

# CFACV Shared Library
 
cfacv_dep = [lapacke, lapack, blas, m, tcl]

chapeau_headers = files(
  'chapeau.h',
  'chapeau_internal.h',
)

chapeau_dense_src = files(
  'chapeau.c',
  'chapeau_io.c',
  'chapeau_solver.c',
  'chapeau_eval.c',
  'chapeau_accumulate.c',
  'chapeau_grid.c',
)
 
cfacv_src = files(
  'measurements.c',
  'measurements.h',
  'centers.c',
  'centers.h',
  'cfacv.c',
  'cfacv_dataspace.c',
  'cfacv_force.c',
  'cfacv_restart.c',
  'cfacv_restraint.c',
  'cfacv.h',
  'cfacv_tcl.h',
  'cfacv_internal.h',
  'otfp_engine.c',
  'otfp_engine.h',
  'cvs.c',
  'cvs.h',
  'cvs_tcl.h',
  'cvs_internal.h',
  )
 
cfacv_src += chapeau_headers

cfacv_src += chapeau_dense_src
                 
cfacv_src += swig_gen.process('cfacv.i')

cfacv_lib = shared_library(
  'cfacv',
  cfacv_src,
  dependencies: cfacv_dep,
  c_args : '-DVERSION="@0@"'.format(meson.project_version()),
  install: true,
  install_dir: 'lib/',
)

lifecycle_test = executable(
  'test-lifecycle',
  files('../tests/lifecycle.c'),
  include_directories: include_directories('.'),
  link_with: cfacv_lib,
  dependencies: cfacv_dep,
)

test('lifecycle', lifecycle_test)

embedding_test = executable(
  'test-embedding',
  files('../tests/embedding.c'),
  include_directories: include_directories('.'),
  link_with: cfacv_lib,
  dependencies: cfacv_dep,
)

test('embedding', embedding_test)

restart_test = executable(
  'test-restart',
  files('../tests/restart.c'),
  include_directories: include_directories('.'),
  link_with: cfacv_lib,
  dependencies: cfacv_dep,
)

test('restart', restart_test)

if fortran_enabled
  embedding_fortran_test = executable(
    'test-embedding-fortran',
    files('../tests/embedding_fortran.f90'),
    include_directories: include_directories('.'),
    link_with: cfacv_lib,
    dependencies: cfacv_dep,
  )

  test('embedding-fortran', embedding_fortran_test)
endif



cvs_src = files(
  'measurements.c',
  'measurements.h',
  'cvs.c',
  'cvs.h',
  'cvs_tcl.h',
  'cvs_internal.h',
)

cvs_src += chapeau_dense_src
cvs_src += chapeau_headers

cvs_src += swig_gen.process('cvs.i')

shared_library(
  'cvs',
  cvs_src,
  dependencies: [lapacke, lapack, blas, m, tcl],
  c_args : '-DVERSION="@0@"'.format(meson.project_version()),
  install: true,
  install_dir: 'lib/',
)  



catbinsp_src = files('catbinsp.c')
 
executable('catbinsp', 
  catbinsp_src, 
  dependencies: [m],
  install: true,
)
 


chapeau_tool_dep = [lapacke, lapack, blas, m]
chapeau_tools = {
  'chapsim': 'simetrize.c',
  'chapadd': 'addchapeaus.c',
  'chapcrop': 'cropchapeaus.c',
}

foreach tool_name, tool_source : chapeau_tools
  executable(
    tool_name,
    files(tool_source) + chapeau_dense_src + chapeau_headers,
    dependencies: chapeau_tool_dep,
    install: true,
  )
endforeach
