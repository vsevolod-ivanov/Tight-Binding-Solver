/*
 *
 *
*/

#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
#include <lapacke.h> 
#include <cblas.h>
#include <math.h>

#include "hamiltonian.h"
//#include "hamiltonian.h"
#include "ham_utils.h"
#include "mat_print.h" 

/* Main program */
int main(int argc, char * * argv) {

	/* Locals */
	lapack_int n, lda, info, n_cell, 
		k_points, nndepth;
	int i, j, k, l, m; //loop indices
	double delta_k;
	double tedge;		//hopping across the edge (0.00 - 1.00)
	
	lapack_complex_double phi1, phi2, phi3; //hopping phases
	lapack_complex_double kx,ky;
	lapack_complex_double t2, spin;
	double cell_dims[3];
	double lat_vectors[9];
	double * atom_vectors;
	double * kpath;
	double * nn;
	double * kpoints;
	int * nncounts;
	int nntotal;
	char structfile[30];

	n_cell = 1;						//Number of unit cells in supercell
	k_points = 10;					//Number of k-points to plot
	t2 = 0.01;
	tedge = 1;
	spin = 1.0;

	/* Read from file */
	//read_coords("diamond.in", cell_dims, lat_vectors, &atom_vectors, &n, &nndepth);
	read_coords("cubic.struct", cell_dims, lat_vectors, &atom_vectors, &kpoints, &n, &nndepth);

	//print_dvector("Cell Dimensions" , 3, cell_dims);
	//print_rmatrix("Lattice Vectors" , 3, 3, lat_vectors, 3);
	//printf("Test");
	//printf("test: vec[%i,%i] = %8.5f\n", 1, 1 , atom_vectors[1*3+1]);
	//printf("%8.5f" , atom_vectors[3]);
	//print_rmatrix("Atom Locations" , 2, 3, atom_vectors, 3);

	//printf("number of atoms: %i\n", n);
	

	nncounts = (int * ) malloc((nndepth+1) * sizeof(int ));
    if (nncounts == NULL) {
    	printf("Memory allocation error when allocating NN counter array.\n");
    	exit(0);
	}
	memset(nncounts, 0, (nndepth+1) * sizeof(int ));
	find_neighbors(n, lat_vectors, atom_vectors, nndepth, &nn, nncounts);

	//printf("NN depth: %i\n", nndepth );
	
	//printf("NNcounts: (%i, %i )\n", nncounts[0], nncounts[1]);



	//Find total number of NN
	nntotal = 0;
	for( i = 0; i < nndepth+1; i++ ) {
		nntotal+=nncounts[i];
	}
	//printf("Total NN: %i\n", nntotal);
	//print_rmatrix("Neighbors", nntotal, 5, nn, 5);

	//Allocate memory for coefficients (to be used in Hamlitonian later) x2 because spin.
	lapack_complex_double* coeftable;
	coeftable = (lapack_complex_double * ) malloc(nntotal * 2 * sizeof(lapack_complex_double));
    if (coeftable == NULL) {
    	printf("Memory allocation error when allocating Hamiltonian coefficient table.\n");
    	exit(0);
	}
	memset(coeftable, 0, nntotal * 2 * sizeof(lapack_complex_double));

	build_C( nndepth, nntotal, nn, nncounts, coeftable, graphene_H);

	//Allocate Hamiltonian memory
	//x4 because spin in each dimension
	lapack_complex_double * H_tb;
	H_tb = (lapack_complex_double * ) malloc(n * n * 4 * sizeof(lapack_complex_double ));
    if (H_tb == NULL) {
    	printf("Memory allocation error when allocating H_tb.\n");
    	exit(0);
	}
	memset(H_tb, 0, n * n * 4 * sizeof(lapack_complex_double ));
	
	//*******************
	// Allocate temp H arrays, initializing to 0
	//*******************

	lapack_complex_double * evc;	//Eigenvectors
	double * eec;	//Eigenenergies
	lapack_complex_double * dHdk;	//Derivatives
	lapack_complex_double * tempP; //Temp array to store inner product
	double * bC; //band berry curvature 
	evc = (lapack_complex_double * ) malloc(n * n * 4 * sizeof(lapack_complex_double ));
    if (evc == NULL) {
    	printf("Memory allocation error when allocating H_tb.\n");
    	exit(0);
	}
	eec = (double * ) malloc(n * 2 * sizeof(double ));
    if (eec == NULL) {
    	printf("Memory allocation error when allocating EigenEnergy Array.\n");
    	exit(0);
	}
	dHdk = (lapack_complex_double * ) malloc(n * n * 4 * 3 * sizeof(lapack_complex_double ));
    if (dHdk == NULL) {
    	printf("Memory allocation error when allocating dHdk.\n");
    	exit(0);
	}
	tempP = (lapack_complex_double * ) malloc(n * n * 4 * sizeof(lapack_complex_double ));
    if (tempP == NULL) {
    	printf("Memory allocation error when allocating temporary matrix product storage.\n");
    	exit(0);
	}
	bC = (double * ) malloc(3 * n * 2 * sizeof(double ));
    if (bC == NULL) {
    	printf("Memory allocation error when allocating band Berry Curvature array.\n");
    	exit(0);
	}
	memset(evc, 0, n * n * 4 * sizeof(lapack_complex_double ));
	memset(eec, 0, n * 2 * sizeof(double ));
	memset(dHdk, 0, n * n * 4 * 3 * sizeof(lapack_complex_double ));
	memset(tempP, 0, n * n * 4 * sizeof(lapack_complex_double ));
	memset(eec, 0, 3 * n * 2 * sizeof(double ));

	//*******************
	// Compute step sizes for k-vectors
	//*******************
	lapack_int nkp = 10; //k point density 
	double kstep[3];
	kstep[0] = 0.25*M_PI/sqrt( pow(lat_vectors[0],2) + pow(lat_vectors[1],2) + pow(lat_vectors[2],2) );
	kstep[1] = 2.0*M_PI/sqrt( pow(lat_vectors[3],2) + pow(lat_vectors[4],2) + pow(lat_vectors[5],2) );
	kstep[2] = 2.0*M_PI/sqrt( pow(lat_vectors[6],2) + pow(lat_vectors[7],2) + pow(lat_vectors[8],2) );
	//divide into nkp steps
	kstep[0] /= ((double) nkp);
	kstep[1] /= ((double) nkp);
	kstep[2] /= ((double) nkp);

	double kvc[3]; //Stores the current k-vector
	double berryC[3]; //Stores the current berry curvature
	lapack_complex_double alpha = {1.0,0.0}; //Matrix multiplication is unscaled
	lapack_complex_double beta = {0.0,0.0};	//Don't copy existing matrix
	lapack_int hamSize = n * n * 4; //number of cells in Hamiltonian

	lapack_complex_double tempBC = 0.0; //temp variable for berry curvature
	
	for (i = -nkp-1; i < nkp+1; i++) {
	//for (i = 10; i < 11; i++) {
		for (j = -nkp-1; j < nkp+1; j++) {
		//for (j = 0; j < 1; j++) {
			//for (k = -nkp-1; k < nkp+1; k++) {
			for (k = 0; k < 1; k++) {
				kvc[0] = (double) i * kstep[0];
				kvc[1] = (double) j * kstep[1];
				kvc[2] = (double) k * kstep[2];

				//printf("\n*****\nk-vec: (%6.2f,%6.2f,%6.2f):\n", kvc[0],kvc[1],kvc[2]);
				
				memset(evc, 0, n * n * 4 * sizeof(lapack_complex_double ));
				memset(dHdk, 0, n * n * 4 * 3 * sizeof(lapack_complex_double ));
				memset(tempP, 0, n * n * 4 * sizeof(lapack_complex_double ));
				memset(eec, 0, n*2*sizeof(double));
				//build_H(n, kvc, nn, nntotal, evc, coeftable );
				build_testH(n, kvc, nn, nntotal, evc, coeftable );
				build_testdHdk(n, kvc, nn, nntotal, dHdk, coeftable );
				//build_dHdk(n, kvc, nn, nntotal, dHdk, coeftable );

				//print_matrix("Hamiltonian: ", 2*n, 2*n, evc, 2*n);
				// print_matrix("dHdkx: ", 2*n, 2*n, dHdk, 2*n);
				// print_matrix("dHdky: ", 2*n, 2*n, &dHdk[hamSize], 2*n);
				// print_matrix("dHdkz: ", 2*n, 2*n, &dHdk[2*hamSize], 2*n);

				//Solve H to get eigenvalues
				info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'L', n * 2, evc, n * 2, eec );
				if( info > 0 ) {
					printf( "The algorithm failed to compute eigenvalues.\n" );
					exit( 1 );
				}
				//print_matrix("Eigenvectors: ", 2*n, 2*n, evc, 2*n);
				//Multiply eigenvalues by each dHdk (might be able to do in one step?)
				cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, 2*n, 2*n, 2*n, &alpha, 
					dHdk, 2*n, evc, 2*n, &beta, tempP, 2*n);
				//print_matrix("first mult dHdkx: ", 2*n, 2*n, tempP, 2*n);
				//get inner product by multiplying by complex conjugate of eigenvalues
				//put result back into dHdk
				cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, 2*n, 2*n, 2*n, &alpha, 
					evc, 2*n, tempP, 2*n, &beta, dHdk, 2*n);
				//repeat for Y (using same tempP array)
				cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, 2*n, 2*n, 2*n, &alpha, 
					&dHdk[hamSize], 2*n, evc, 2*n, &beta, tempP, 2*n);
				cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, 2*n, 2*n, 2*n, &alpha, 
					evc, 2*n, tempP, 2*n, &beta, &dHdk[hamSize], 2*n);
				//Repeat for Z
				cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, 2*n, 2*n, 2*n, &alpha, 
					&dHdk[2*hamSize], 2*n, evc, 2*n, &beta, tempP, 2*n);
				cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, 2*n, 2*n, 2*n, &alpha, 
					evc, 2*n, tempP, 2*n, &beta, &dHdk[2*hamSize], 2*n);
				
				//print_dvector("EigenEnergies: ", 2*n, eec);
				//printf("k-vec: (%6.2f,%6.2f,%6.2f):\n", kvc[0],kvc[1],kvc[2]);
				// print_matrix("Mult-X result: ", 2*n, 2*n, dHdk, 2*n);
				// print_matrix("Mult-Y result: ", 2*n, 2*n, &dHdk[hamSize], 2*n);
				// print_matrix("Mult-Z result: ", 2*n, 2*n, &dHdk[2*hamSize], 2*n);
				
				memset(berryC, 0, 3*sizeof(double));
				memset(bC, 0, 3*n*2*sizeof(double));
				
				for (l = 0; l < n*2; l++) {
					for(m = 0; m < n*2; m++) {
						//Make sure that n'=/=n , i.e l=/=m
						if( l != m ){
							//Also helps make sure that there is no division by 0
							if( fabs(eec[l]-eec[m]) > 0.001 ){
								//Note the l/m switch in each term TODO
								//x => y * z*
								tempBC = cimag(dHdk[hamSize+l*n*2 + m]*conj(dHdk[2*hamSize+l*n*2 + m]));
								bC[l*3+0] += tempBC/pow((eec[l]-eec[m]),2);
								//y => z * x*
								tempBC = cimag(dHdk[2*hamSize+l*n*2 + m]*conj(dHdk[l*n*2 + m]));
								bC[l*3+1] += tempBC/pow((eec[l]-eec[m]),2);
								//z => x * y*
								tempBC = cimag(dHdk[l*n*2 + m]*conj(dHdk[hamSize+l*n*2 + m]));
								bC[l*3+2] += tempBC/pow((eec[l]-eec[m]),2);
							}
						} 
					}
				}
				
				//printf("curvature: (%7.6f,%7.6f,%7.6f):\n", berryC[0],berryC[1],berryC[2]);	

				//printf("%6.4f\t%6.4f\t%6.4f\n", kvc[0],kvc[1],kvc[2]);
				//printf("%6.4f\t%6.4f\t%6.4f\n", berryC[0],berryC[1],berryC[2]);	
				printf("%6.4f\t%6.4f\t%6.4f\n", bC[0],bC[1],bC[2]);
				//printf("%6.4f\t%6.4f\t%6.4f\n", bC[3],bC[4],bC[5]);	
				//printf("%6.4f\t%6.4f\n", eec[0], eec[1]);	
				
			}
		}
	}

	/* Arguments */
    for (i = 1; i < argc; i++) {
		
		//Get number of cells in supercell
        if (strcmp(argv[i], "-n") == 0) {
			n_cell = atoi(argv[i + 1]);
        	i++;
		}
		if (strcmp(argv[i], "-k") == 0) {
        	k_points = atoi(argv[i + 1]);
        	i++;
		}
		//percentage of hopping intra-strip
		if (strcmp(argv[i], "-edge") == 0) {
        	tedge = atof(argv[i + 1]);
        	i++;
		}
		//relative strength of t2/t (nnn hopping)
		if (strcmp(argv[i], "-t2") == 0) {
        	t2 = atof(argv[i + 1]);
        	i++;
		}
		//spin (+1 or -1)
		if (strcmp(argv[i], "-s") == 0) {
        	spin = atof(argv[i + 1]);
        	i++;
		}
		//file
		if (strcmp(argv[i], "-f") == 0) {
			memset(structfile, '\0', sizeof(structfile));
			strcpy(structfile, argv[i + 1]);
			//printf("Test: %s\n", structfile);
        	i++;
		}
	}

	/* Default Values */
	n = 2*n_cell;						//Total number of sites (also H size)
	delta_k = (double) 2.0*M_PI/((float)k_points);	//k-vector increment
	delta_k = delta_k/sqrt(3);
	t2 = t2 * I;
    // Local arrays
	

	
	exit(0);
	
}