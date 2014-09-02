#!/bin/bash

# you might need to set these variables for your compiler
mpi_compiler=mpic++

#compiler=g++-4.9
#compiler=clang++ # works fine on MacOS, but has no openMP support
#compiler=gcc # works fine, but no openMP support with deault gcc on MacOS
compiler=g++

where_make="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## first of all, install libLBFGS to a local directory
echo compiling libLBFGS
rm -r ${where_make}/main_code/lbfgs >/dev/null 2>&1
${where_make}/main_code/liblbfgs-1.10/configure --prefix=${where_make}/main_code/lbfgs --srcdir=${where_make}/main_code/liblbfgs-1.10 >/dev/null 2>&1
make clean >/dev/null 2>&1
make install >/dev/null 2>&1

# check that libLBFGS installed
if ls ${where_make}/main_code/lbfgs >/dev/null 2>&1; then
    echo [libLBFGS compiled successfully];
else
    echo [libLBFGS compile failed!];
    exit 1
fi; 

# remove temporary files created by installation of liblbfgs
rm config.h config.log config.status libtool Makefile stamp-h1
rm -r lib sample

# make sure compiler can be found
if command -v ${compiler} >/dev/null 2>&1; then
    echo using specified '"'${compiler}'"' compiler
else
    echo ERROR! cannot find compiler: '"'${compiler}'"'.  Edit make.sh to specify your preferred compiler.
    exit 1
fi

where_lbfgs=${where_make}/main_code/lbfgs
where_eigen=${where_make}/main_code/eigen322

# compile a dummy program to see if compiler will compile with openMP flags.  If this doesn't compile openMP will be automatically disabled
rm ${where_make}/test_openmp 2> /dev/null
${compiler} -fopenmp ${where_make}/main_code/test_openmp.cpp -o ${where_make}/test_openmp 2>/dev/null
if command -v ${where_make}/test_openmp >/dev/null 2>&1; then
    echo openMP support detected
    use_openmp=true
    rm ${where_make}/test_openmp
else
    echo openMP not detected
    use_openmp=false
fi

## theoretically, don't need lbfgs for infer_{MRF/CRF} executables, but because of the way the code is structured, included anyway

${compiler} -O3 -Ofast -funroll-loops -std=c++11 main_code/infer_MRF.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o infer_MRF
if command -v ${where_make}/infer_MRF >/dev/null 2>&1; then
    echo "[infer_MRF compiled.]"
else
    echo "[infer_MRF compile failed!]"
    exit 1
fi

${compiler} -O3 -Ofast -funroll-loops -std=c++11 main_code/infer_CRF.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o infer_CRF
if command -v ${where_make}/infer_MRF >/dev/null 2>&1; then
    echo "[infer_CRF compiled.]"
else
    echo "[infer_CRF compile failed!]"
    exit 1
fi

if $use_openmp; then
    ${compiler} -O3 -Ofast -funroll-loops -std=c++11 main_code/learn_CRF.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o learn_CRF -fopenmp
else
    ${compiler} -O3 -Ofast -funroll-loops -std=c++11 main_code/learn_CRF.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o learn_CRF
fi

if ! command -v ${where_make}/learn_CRF >/dev/null 2>&1; then
    echo "[learn_CRF compile failed!]"
    exit 1
fi

if ${use_openmp}; then
    echo "[learn_CRF compiled with openMP multithreading.]"
else
    echo "[learn_CRF compiled without openMP multithreading.]"
fi

if command -v ${mpi_compiler} >/dev/null 2>&1; then
    echo using specified '"'${mpi_compiler}'"' compiler for MPI
else
    echo Cannot find MPI compiler: '"'${mpi_compiler}'"'.  Thus, will not produce learn_CRF_mp executable.
    echo Install openMPI and edit make.sh to specify mpi_compiler if you want to use MPI.
    echo You can just use learn_CRF instead if you prefer.
    exit 1
fi

export OMPI_CXX=${compiler}
mpic++ -O3 -Ofast -funroll-loops -std=c++11 main_code/learn_CRF_mpi.cpp -llbfgs -I${where_lbfgs}/include -L${where_lbfgs}/lib -I${where_eigen} -o learn_CRF_mpi

echo "[learn_CRF_mpi compiled.]"
