#! /bin/bash
# Need to build SUNDIALS with -fPIC and set(CMAKE_MACOSX_RPATH 1)

#export CANTERA_DIR=/Users/wuhao/Codes/Cantera_Chemistry/Cantera/merge_test
#export SUN_INCLUDE=/Users/wuhao/Codes/tools/sundials/2.6.2/include
#export SUN_LIB=/Users/wuhao/Codes/tools/sundials/2.6.2/lib
#export BOOST_DIR=/usr/local/Cellar/boost/1.62.0/include
#export EIGEN_INCLUDE=/usr/local/Cellar/eigen/3.2.8/include/eigen3
#
#scons -j5 build prefix=$CANTERA_DIR CXX=clang++ CC=clang python_package=full \
#    sundials_include=$SUN_INCLUDE sundials_libdir=$SUN_LIB matlab_path=/Applications/MATLAB_R2014b.app/ \
#    boost_inc_dir=$BOOST_DIR f90_interface=n system_eigen=y extra_inc_dirs=$EIGEN_INCLUDE

export CANTERA_DIR=/Users/peterma/local/cantera/2.4.0
export SUN_INCLUDE=/usr/local/Cellar/sundials/2.7.0_1/include
export SUN_LIB=/usr/local/Cellar/sundials/2.7.0_1/lib
export BOOST_DIR=/usr/local/Cellar/boost/1.63.0/include
export EIGEN_INCLUDE=/usr/local/Cellar/eigen/3.3.1/include/eigen3

scons -j4 build prefix=$CANTERA_DIR CXX=clang++ CC=clang python_package=full \
    sundials_include=$SUN_INCLUDE sundials_libdir=$SUN_LIB \
    boost_inc_dir=$BOOST_DIR f90_interface=n system_eigen=y extra_inc_dirs=$EIGEN_INCLUDE

# scons test
scons -j4 install
