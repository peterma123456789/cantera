#! /bin/bash
# Need to build SUNDIALS with -fPIC and set(CMAKE_MACOSX_RPATH 1)
export CANTERA_DIR=/Users/peterma/local/cantera/2.2.0
export SUN_INCLUDE=/usr/local/Cellar/sundials/2.7.0_1/include
export SUN_LIB=/usr/local/Cellar/sundials/2.7.0_1/lib
export BOOST_DIR=/usr/local/Cellar/boost/1.63.0/include

scons -j4 build prefix=$CANTERA_DIR CXX=clang++ CC=clang python_package=full \
    sundials_include=$SUN_INCLUDE sundials_libdir=$SUN_LIB \
    boost_inc_dir=$BOOST_DIR f90_interface=n \
    cxx_flags=-std=c++11
# scons test
scons -j4 install
