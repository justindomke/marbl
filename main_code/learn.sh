#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -p -100

export LD_LIBRARY_PATH=/share/apps/gcc/current/lib:/share/apps/gcc/current/lib64/:/share/apps/openmpi/current/gcc/lib
export PATH=/share/apps/gcc/current/bin:/share/apps/openmpi/current/gcc/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~
mpirun -np $NSLOTS ./learn_mpi -d data/*.txt -m models/*.txt -i 5 -a lbfgs -r .001 -w W.txt