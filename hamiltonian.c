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
	lapack_complex_double param;// = {{0.0,0.0}};
	//printf("\nFINDING COEFFICIENT\n\n");
	switch(n){
		case 0 :
			param = coef[n];
			break;
		case 1 :
			param = coef[n];
			break;
		case 2 : 
			//cblas_zdotu_sub(3, spin, 1, vec, 1, &param);
			param = spin[0]*vec[0]+spin[1]*vec[1]+spin[2]*vec[2];
			//printf( "Vec: (%6.4f, %6.4f, %6.4f)\n", vec[0], vec[1], vec[2] );
			//printf( "Spin: (%6.4f, %6.4f, %6.4f)\n", spin[0], spin[1], spin[2] );
			//printf( "Param: (%6.4f, %6.4f)\n", creal(param), cimag(param));
			param = param * I * coef[n];
			break;
		default :
			printf("Invalid NN depth!\n");
	}
	return param;
}
