#include <lapacke.h>
#include <cblas.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hamiltonian.h"

lapack_complex_double graphene_H(lapack_int n, double* spin, double* vec, lapack_complex_double* coef) {

	lapack_complex_double param;

	switch(n){
		case 0:
			param = coef[n];
		case 1:
			param = coef[n];
		case 2: 
			cblas_zdotu_sub(3, spin, 1, vec, 1, &param);
			param = param * I * coef[n];
	}
	return param;
}