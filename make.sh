#!/bin/bash

# you might need to set these variables
use_openmp=true
mpi_compiler=mpic++

compiler=g++-4.9
#compiler=clang++
#compiler=g++


where_make="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# make sure compiler can be found
if hash ${compiler} 2>/dev/null; then
    echo using specified ${compiler} compiler
else
    echo ERROR! cannot find compiler: ${compiler}.  Edit make.sh to specify your preferred compiler.
    exit 1
fi

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

${compiler} -O3 -Ofast -funroll-loops -std=c++11 main_code/infer_MRF.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o infer_MRF
${compiler} -O3 -Ofast -funroll-loops -std=c++11 main_code/infer_CRF.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o infer_CRF

if $use_openmp; then
    ${compiler} -O3 -Ofast -funroll-loops -std=c++11 main_code/learn_CRF.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o learn_CRF
else
    ${compiler} -O3 -Ofast -funroll-loops -fopenmp -std=c++11 main_code/learn_CRF.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o learn_CRF
fi

if hash ${mpi_compiler} 2>/dev/null; then
    echo using specified ${mpi_compiler} compiler for MPI
else
    echo Cannot find MPI compiler: ${mpi_compiler}.  Thus, will not produce learn_CRF_mp executable.
    echo Install openMPI and edit make.sh to specify mpi_compiler if you want to use MPI.
    echo You can just use learn_CRF instead if you prefer.
    exit 1
fi

export OMPI_CXX=${compiler}
mpic++ -O3 -Ofast -funroll-loops -std=c++11 main_code/learn_CRF_mpi.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o learn_mpi

