#!/bin/bash -x

# on Mac OS, default compiler if problematic
#compiler=g++-4.9
#compiler=clang++
compiler=g++

where_make="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## first of all, install libLBFGS to a local directory
#${where_make}/main_code/liblbfgs-1.10/configure --prefix=${where_make}/main_code/lbfgs --srcdir=${where_make}/main_code/liblbfgs-1.10
#make clean
#make install

## remove temporary files created by installation of liblbfgs
#rm config.h config.log config.status libtool Makefile stamp-h1
#rm -r lib sample

where_lbfgs=${where_make}/main_code/lbfgs
where_eigen=${where_make}/main_code/eigen322

## theoretically, don't need lbfgs for infer_{MRF/CRF} executables, but because of the way the code is structured, included anyway

## previously used -flto flag
${compiler} -O3 -Ofast -funroll-loops -std=c++11 main_code/infer_MRF.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o infer_MRF
#${compiler} -O3 -Ofast -funroll-loops -std=c++11 main_code/infer_CRF.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o infer_CRF

#export OMPI_CXX=${compiler}
#mpic++ -O3 -Ofast -funroll-loops -std=c++11 main_code/learn_mpi.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o learn_mpi