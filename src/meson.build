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
dss_name = get_option('dss')
  
blas = dependency(blas_name)
lapack = dependency(lapack_name)
dss = dependency(dss_name, required : get_option('sparse'))

# LAPACKE dependency only if provieder need it
if lapacke_name == 'mkl_lapacke'
    add_project_arguments('-DUSE_MKL', language: 'c')
    lapacke = dependency('mkl_lapacke', required: true)
else
    lapacke = dependency(lapacke_name, method: 'pkg-config', required: true)
endif

install_headers('calcforces.tcl')
install_headers('cfacv.tcl')
install_headers('chapeau.tcl')
install_headers('cvs.tcl')
install_headers('data.tcl')
install_headers('restraints.tcl')

# CFACV Shared Library
 
cfacv_dep = [lapacke, lapack, blas, m, tcl]
 
cfacv_src = files(
  'measurements.c',
  'measurements.h',
  'chapeau.h',
  'centers.c',
  'centers.h',
  'cfacv.c',
  'cfacv.h',
  'cfacv_obj.h',
  'cvs.c',
  'cvs.h',
  'cvs_obj.h',
  )
 
# Using sparse matrix requires MKL DSS
cfacv_src += files('chapeau_obj.h')
if get_option('sparse').enabled()
  
  if dss.found()
    cfacv_dep += dss
    # message('Using MKL ', mkl_name)
    add_project_arguments('-DUSE_MKL', language: 'c')
  else
    error('Intel MKL Direct Sparse Solver is required.')
  endif

  cfacv_src += files('chapeau_sparse.c')
  error('Direct Sparse Solver implementation is under construction.')

else
  cfacv_src += files('chapeau.c')
endif
                 
cfacv_src += swig_gen.process('cfacv.i')

shared_library(
  'cfacv',
  cfacv_src,
  dependencies: cfacv_dep,
  c_args : '-DVERSION="@0@"'.format(meson.project_version()),
  install: true,
  install_dir: 'lib/',
)



cvs_src = files(
  'measurements.c',
  'measurements.h',
  'cvs.c',
  'cvs.h',
  'cvs_obj.h',
  'chapeau.c',
  'chapeau.h',
  'chapeau_obj.h',
)

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
 


chapsim_src = files(
  'simetrize.c',
  'chapeau.c',
  'chapeau.h',
  'chapeau_obj.h',
)


executable('chapsim', 
  chapsim_src, 
  dependencies: [lapacke, lapack, blas, m],
  install: true,
)
 


chapconv_src = files(
  'convchapeau.c',
  'chapeau.c',
  'chapeau.h',
  'chapeau_obj.h',
  'chapeau_old.h',
  ) 

executable('chapconv', 
  chapconv_src, 
  dependencies: [lapacke, lapack, blas, m],
  install: true,
  )



chapadd_src = files(
  'addchapeaus.c',
  'chapeau.c',
  'chapeau.h',
  'chapeau_obj.h',
  )

executable('chapadd', 
  chapadd_src, 
  dependencies: [lapacke, lapack, blas, m],
  install: true,
  )



chapcrop_src = files(
  'cropchapeaus.c',
  'chapeau.c',
  'chapeau.h',
  'chapeau_obj.h',
  )

executable('chapcrop', 
  chapcrop_src, 
  dependencies: [lapacke, lapack, blas, m],
  install: true,
  )



