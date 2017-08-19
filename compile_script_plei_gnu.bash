#! /bin/bash
# Need to build SUNDIALS with -fPIC and set(CMAKE_MACOSX_RPATH 1)
export CANTERA_DIR=/nobackup/cma3/local/cantera/2.4.0
export SUN_INCLUDE=/nobackup/cma3/local/sundials/2.7.0/include
export SUN_LIB=/nobackup/cma3/local/sundials/2.7.0/lib
export BOOST_DIR=/nasa/pkgsrc/sles12/2016Q4/views/boost/1.62/include
export EIGEN_INCLUDE=/nobackup/cma3/local/eigen/3.2.9/eigen/

scons -j20 build prefix=$CANTERA_DIR \
  CXX=g++ CC=gcc FORTRAN=gfortran python_package=full \
  optimize_flags='-O3 -march=ivybridge' \
  env_vars='all' \
  sundials_include=$SUN_INCLUDE sundials_libdir=$SUN_LIB \
  boost_inc_dir=$BOOST_DIR f90_interface=y system_eigen=y extra_inc_dirs=$EIGEN_INCLUDE
if [ -z "SCONS_TEST" ]; then
    scons -j20 test
fi
scons -j20 install
