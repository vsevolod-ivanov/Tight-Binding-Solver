#include <lapacke.h>
#include <cblas.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hamiltonian.h"

/* Computes the coefficient of the n-th NN term
 *
 */
lapack_complex_double graphene_H(lapack_int n, double* spin, double* vec, lapack_complex_double* coef){
	lapack_complex_double param = {0.0,0.0};
	switch(n){
		case 0 :
			param = coef[n];
			break;
		case 1 :
			param = coef[n];
			break;
		case 2 : 
			cblas_zdotu_sub(3, spin, 1, vec, 1, &param);
			param = param * I * coef[n];
			break;
		default :
			printf("Invalid NN depth!\n");
	}
	return param;
}
