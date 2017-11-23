#!/bin/bash

rm tight_binding.x
gcc -o TB_solve.o -c TB_solve.c
gcc -o mat_print.o -c mat_print.c
gcc -o tight_binding.x TB_solve.o mat_print.o  -llapack -llapacke -lrefblas -lcblas
rm *.o