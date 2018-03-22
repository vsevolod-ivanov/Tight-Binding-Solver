#!/bin/bash

rm tight_binding.x
gcc -o TB_solve.o -c TB_solve.c
gcc -o mat_print.o -c mat_print.c
gcc -o ham_utils.o -c ham_utils.c
gcc -o hamiltonian.o -c hamiltonian.c
gcc -o tight_binding.x TB_solve.o mat_print.o hamiltonian.o ham_utils.o -L/usr/local/gfortran/lib/ -lgfortran -llapack -llapacke -lrefblas -lcblas
rm *.o
