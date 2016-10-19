#! /bin/bash
# Need to build SUNDIALS with -fPIC and set(CMAKE_MACOSX_RPATH 1)
export CANTERA_DIR=/Users/wuhao/Codes/Cantera_Chemistry/Cantera/2.3.0
export SUN_INCLUDE=/Users/wuhao/Codes/tools/sundials/2.6.2/include
export SUN_LIB=/Users/wuhao/Codes/tools/sundials/2.6.2/lib
export BOOST_DIR=/usr/local/Cellar/boost/1.59.0/include
export EIGEN_INCLUDE=/usr/local/Cellar/eigen/3.2.6/include/eigen3

scons build prefix=$CANTERA_DIR CXX=clang++ CC=clang FORTRAN=gfortran-5 python_package=full \
    sundials_include=$SUN_INCLUDE sundials_libdir=$SUN_LIB matlab_path=/Applications/MATLAB_R2014b.app/ \
    boost_inc_dir=$BOOST_DIR f90_interface=y eigen_include=$EIGEN_INCLUDE
# scons test
scons install