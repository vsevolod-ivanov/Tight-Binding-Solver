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
#include "mat_print.h"
#include "hamiltonian.h"
#include "hamiltonian.h"
 
/* Main program */
int main(int argc, char * * argv) {

	/* Locals */
	lapack_int n, lda, info, n_cell, 
		k_points, nndepth;
	int i, j ,l;
	double delta_k;
	double tedge;		//hopping across the edge (0.00 - 1.00)
	
	lapack_complex_double phi1, phi2, phi3; //hopping phases
	lapack_complex_double kx,ky;
	lapack_complex_double t2, spin;
	double cell_dims[3];
	double lat_vectors[9];
	double * atom_vectors;
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
	read_coords("graphene.struct", cell_dims, lat_vectors, &atom_vectors, &n, &nndepth);

	print_dvector("Cell Dimensions" , 3, cell_dims);
	print_rmatrix("Lattice Vectors" , 3, 3, lat_vectors, 3);
	//printf("Test");
	//printf("test: vec[%i,%i] = %8.5f\n", 1, 1 , atom_vectors[1*3+1]);
	//printf("%8.5f" , atom_vectors[3]);
	print_rmatrix("Atom Locations" , 2, 3, atom_vectors, 3);

	printf("number of atoms: %i\n", n);
	printf("number of nn: %i\n", nndepth);

	nncounts = (int * ) malloc(nndepth * sizeof(int ));
    if (nncounts == NULL) {
    	printf("Memory allocation error when allocating NN counter array.\n");
    	exit(0);
	}
	memset(nncounts, 0, nndepth * sizeof(int ));
	find_neighbors(n, lat_vectors, atom_vectors, nndepth, &nn, nncounts);

	//printf("Number of NN: %i\n", )
	//print_rmatrix("Neighbors", 34, 5, nn, 5);
	printf("NNcounts: (%i, %i, %i )\n", nncounts[0], nncounts[1], nncounts[2]);

	//char * test;
	read_H("diamond.in");
	

	double ktest[3] = {1.2, 0.3, 0.0};	
	lapack_complex_double* Htest;


	//Find total number of NN
	nntotal = 0;
	for( i = 0; i < nndepth+1; i++ ) {
		nntotal+=nncounts[i];
	}

	//Allocate memory for coefficients (to be used in Hamlitonian later) x2 because spin.
	lapack_complex_double* coeftable;
	coeftable = (lapack_complex_double * ) malloc(nntotal * 2 * sizeof(lapack_complex_double));
    if (coeftable == NULL) {
    	printf("Memory allocation error when allocating Hamiltonian coefficient table.\n");
    	exit(0);
	}
	memset(coeftable, 0, nntotal * 2 * sizeof(lapack_complex_double));

	build_C( nndepth, nntotal, nn, nncounts, coeftable, graphene_H);


	printf("N: %i\n\n", n);
	lapack_complex_double * H_tb;
	//Allocate Hamiltonian memory
	//x4 because spin in each dimension
	H_tb = (lapack_complex_double * ) malloc(n * n * 4 * sizeof(lapack_complex_double ));
    if (H_tb == NULL) {
    	printf("Memory allocation error when allocating H_tb.\n");
    	exit(0);
	}

	memset(H_tb, 0, n * n * 4 * sizeof(lapack_complex_double ));
	
	build_H(n, ktest, nn, nntotal, H_tb, coeftable );
	print_matrix("Tight-Binding Matrix (before): ", 2*n, 2*n, H_tb, 2*n);
	

	double eigs[4] = {0.0, 0.0, 0.0, 0.0};
	//test diagonalization:
	info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'L', 2*n, H_tb, 2*n, eigs );
	
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	print_matrix("Tight-Binding Matrix (after): ", 2*n, 2*n, H_tb, 2*n);
	print_dvector("Eig: ", 2*n, eigs);

	//generate some k points to loop over:

	lapack_int numkpoints = 20; //temp value to store k point density 
	for (i = 1; i < argc; i++) {

	}



	

	//printf("nn[0]: %6.2f", nn[0]);
	//print_dvector("test", 5, nn[0]);
	//print_dvector("test", 5, nn[5]);


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

	//printf("n_cell: %i, kpoints: %i\n\n", n_cell, k_points);
	//printf("n: %i \n", n);
	//printf("0\n");
	/* Default Values */
	n = 2*n_cell;						//Total number of sites (also H size)
	//printf("1\n");
	delta_k = (double) 2.0*M_PI/((float)k_points);	//k-vector increment
	delta_k = delta_k/sqrt(3);
	t2 = t2 * I;
	//printf("2\n");
	//delta_k = 2.0*M_PI/((float)k_points);	//k-vector increment
	//printf("3\n");
	//printf("4\n");
	//k_points++;

	
	//printf("n: %i \n", n);

	//printf("k-vector increment: %6.2f\n",delta_k);

    // Local arrays
	
	//printf("test0");

	double * w; //[n]; //eigenvalue matrix for each k point
	double * bands; //bands[n*k_points]; //band matrix

	//printf("test0");

	
	lapack_complex_double k[2];

	//printf("test01");
	
	//nn hopping vectors
	lapack_complex_double 
		d1[2] = {{-1.0,0.0}, {0.0,0.0}}, 
		d2[2] = {{0.5,0.0}, {0.5*sqrt(3.0),0.0}},
		d3[2] = {{0.5,0.0}, {-0.5*sqrt(3.0),0.0}};

	//nnn hopping vectors
	lapack_complex_double 
		a1[2] = {{0.0,0.0}, {sqrt(3),0.0}}, 
		a2[2] = {{1.5,0.0}, {0.5*sqrt(3.0),0.0}},
		a3[2] = {{1.5,0.0}, {-0.5*sqrt(3.0),0.0}};

	//printf("test001\n\n");	
	
	//printf("test1");
	// Initialization 
	lda = n;
	/*
    H_tb = (lapack_complex_double * ) malloc(n * n * sizeof(lapack_complex_double ));
    if (H_tb == NULL) {
    	printf("Memory allocation error when allocating H_tb.\n");
    	exit(0);
	}*/
	w = (double * ) malloc(n * sizeof(double ));
    if (w == NULL) {
    	printf("Memory allocation error when allocating temp eigenvalue array.\n");
    	exit(0);
	}
	bands = (double * ) malloc(n * k_points * sizeof(double ));
    if (bands == NULL) {
    	printf("Memory allocation error when allocating band array.\n");
    	exit(0);
	}

	
	//printf("test2\n\n");
	//Set all to 0 for convenience
	memset(w, 0, n * sizeof(double ));
	memset(H_tb, 0, n * n * sizeof(lapack_complex_double ));
	memset(bands, 0, n * k_points * sizeof(double ));

	
	
	//START CALCULATION 

	
	kx = M_PI*0.0; // Want to be on the K, K' line

	lapack_complex_double p1_nn,p2_nn,p3_nn;		//nn phases
	lapack_complex_double p1_nnn,p2_nnn,p3_nnn;		//nnn phases
	int row_num, col_num;				//location in H_tb
	
	//k_points
	for (i = 0; i<1; i++) {

		ky = delta_k*i;
		k[0]  = kx; k[1] =ky;

		k[0] = 1.2;
		k[1] = 0.3;

		//printf("%6.2f\n\n",kx);
		//printf( "k: (%6.2f, %6.2f)\n\n", creal(k[0]), creal(k[1]));
		cblas_zdotu_sub(2, k, 1, d1, 1, &phi1);
		cblas_zdotu_sub(2, k, 1, d2, 1, &phi2);
		cblas_zdotu_sub(2, k, 1, d3, 1, &phi3);

		p1_nn = cexp(I*phi1);
		p2_nn = cexp(I*phi2);
		p3_nn = cexp(I*phi3);

		printf( "k: (%6.2f, %6.2f)\n\n", creal(k[0]), creal(k[1]));
		printf( "a1: (%6.2f, %6.2f)\n\n", creal(a1[0]), creal(a1[1]));
		cblas_zdotu_sub(2, k, 1, a1, 1, &phi1);
		cblas_zdotu_sub(2, k, 1, a2, 1, &phi2);
		cblas_zdotu_sub(2, k, 1, a3, 1, &phi3);

		printf( "phi1: (%6.2f, %6.2f)\n\n", creal(phi1), cimag(phi1));
		p1_nnn = cexp(I*phi1);
		//printf( "p1: (%6.2f, %6.2f)\n\n", creal(p1_nnn), cimag(p1_nnn));
		//printf( "p1*: (%6.2f, %6.2f)\n\n", creal(conj(p1_nnn)), cimag(conj(p1_nnn)));
		p2_nnn = cexp(I*phi2);
		p3_nnn = cexp(I*phi3);

		//printf( "Phase 1: (%6.2f, %6.2f)\n\n", creal(p1_nnn), cimag(p1_nnn));
		//printf( "Phase 2: (%6.2f, %6.2f)\n\n", creal(p2_nnn), cimag(p2_nnn));
		//printf( "Phase 3: (%6.2f, %6.2f)\n\n", creal(p3_nnn), cimag(p3_nnn));


		/*//Print Phases
		printf("\n\n");
		printf( "Phi's: (%6.2f, %6.2f), (%6.2f, %6.2f), (%6.2f, %6.2f)\n\n",
			creal(phi1), cimag(phi1), creal(phi2), cimag(phi2), creal(phi3), cimag(phi3) );

		printf(  "Phase's: (%6.2f, %6.2f), (%6.2f, %6.2f), (%6.2f, %6.2f)\n\n",
			creal(p1_nn), cimag(p1_nn), creal(p2_nn), cimag(p2_nn), creal	p3_nn), cimag	p3_nn) );

		printf(  "Phase's: (%6.2f, %6.2f), (%6.2f, %6.2f), (%6.2f, %6.2f)\n\n",
			creal(conj(p1_nn)), cimag(conj(p1_nn)), creal(p2_nn), cimag(p2_nn), creal	p3_nn), cimag	p3_nn) );
			
		*/
		//Generate H_tb
		memset(w, 0, n * sizeof(double ));
		memset(H_tb, 0, n * n * sizeof(lapack_complex_double ));
		//print_matrix("Tight-Binding Matrix (before): ", n, n, H_tb, lda);		

		//Add nn terms
		for (j = 0; j < n; j=j+2) {
			//first off diagonal term
			row_num = (j+1)%n;
			col_num = (j)%n;
			H_tb[row_num*n+col_num] += conj(p2_nn) + conj(p3_nn);
			//print_matrix("TB1 ", n, n, H_tb, lda);
			col_num = (j+2)%n;
			H_tb[row_num*n+col_num] += conj(p1_nn);
			//print_matrix("TB2 ", n, n, H_tb, lda);
			row_num = (j+2);
			if((row_num+1)<n){
				//printf("n: %i",n);
				//printf("rownum: %i" , row_num);
				col_num = (j+1)%n;
				H_tb[row_num*n+col_num] += conj(p1_nn);
				//printf("before crash");
				//print_matrix("TB3 ", n, n, H_tb, lda);
				//printf("after crash?");
				//col_num = (j+3)%n;
				//H_tb[row_num*n+col_num] += p1_nn;
			}
		}
		//Remove nn coupling between strips
		H_tb[(n-1)*n] *= tedge;

		//Add nnn terms
		for (j = 0; j < n; j=j+2) {
			//A sites

			//diag
			row_num = j;
			col_num = j;
			H_tb[row_num*n+col_num] += spin * t2 *( -p1_nnn + conj(p1_nnn) );
			//printf("%6.2f\n", creal(2 * spin * t2 ));
			//printf("%6.2f\n", creal( -p1_nnn + conj(p1_nnn) ));
			//printf("Diag: (%6.2f, %6.2f)\n", creal(H_tb[row_num*n+col_num]), cimag(H_tb[row_num*n+col_num]));
			
			//off-diag
			col_num = (j+2)%n;
			H_tb[row_num*n+col_num] += spin * t2 *( p2_nnn - p3_nnn);
			col_num = (j+4)%n;
			H_tb[row_num*n+col_num] += spin * t2 *( -conj(p2_nnn) + conj(p3_nnn) );
			
			//B sites
			row_num = (j+1);
			col_num = (j+1)%n;
			H_tb[row_num*n+col_num] += -2 * spin * t2 *( -p1_nnn + conj(p1_nnn) );
			//print_matrix("TB1 ", n, n, H_tb, lda);
			col_num = (j+3)%n;
			H_tb[row_num*n+col_num] += -spin * t2 *( p2_nnn - p3_nnn);
			col_num = (j+5)%n;
			H_tb[row_num*n+col_num] += -spin * t2 *( -conj(p2_nnn) + conj(p3_nnn) );
		}

		//Remove nn coupling between strips
		//1st coupling to 2nd to last row:n-2, col:0
		H_tb[(n-2)*n] *= tedge;
		//2nd coupling to last, row: n-1, col:1
		H_tb[(n-2)*n+1] *= tedge;



		// Display H_tb to check 
		print_matrix("Tight-Binding Matrix (after init): ", n, n, H_tb, lda);

		info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'L', n, H_tb, lda, w );
		
		if( info > 0 ) {
			printf( "The algorithm failed to compute eigenvalues.\n" );
			exit( 1 );
		}

		//Save Eigenvalues

		for(j = 0; j<n; j++)
		{
			bands[i*n+j] = w[j];
		}

		//Eigenvalues
		print_rmatrix( "Eigenvalues", 1, n, w, 1 );
	
		//Eigenvectors
		print_matrix("Eigenvectors, columnwise: ", n, n, H_tb, lda);

		//printf("---------finished k-point %i------------\n\n",i+1 );
	}
	
	//print_rmatrix( "", k_points, n, bands, n );

	//END CALCULATION 
	

	//printf("Done!\n\n");
	
	exit(0);
	
}