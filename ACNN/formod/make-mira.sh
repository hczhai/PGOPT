#!/bin/bash

BLASLIB=-L/soft/libraries/alcf/current/gcc/BLAS/lib
LAPACKLIB=-L/soft/libraries/alcf/current/gcc/LAPACK/lib
for ffile in `ls *.f90`; do
    f2py -c -m ${ffile%.*} ${ffile} -llapack -lblas $BLASLIB $LAPACKLIB
done
