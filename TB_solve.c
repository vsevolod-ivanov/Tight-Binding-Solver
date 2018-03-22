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
#include <time.h>

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
	read_coords("TaAs.struct", cell_dims, lat_vectors, &atom_vectors, &kpoints, &n, &nndepth);

	//print_dvector("Cell Dimensions" , 3, cell_dims);
	//print_rmatrix("Lattice Vectors" , 3, 3, lat_vectors, 3);
	//printf("Test");
	//printf("test: vec[%i,%i] = %8.5f\n", 1, 1 , atom_vectors[1*3+1]);
	//printf("%8.5f" , atom_vectors[3]);
	//print_rmatrix("Atom Locations" , 32, 3, atom_vectors, 3);

	//printf("number of atoms: %i\n", n);

	int nmap;
	int * Hmap;
	int mapS[3] = {0,0,0}; //dimensions/shifts for the mapping
	read_Map( "TaAs_ham.plb", &nmap, mapS, &Hmap);
	//printf("Nmap: %i\n", nmap);

	//check that we gottem
	//printf("Mapshifts: (%i,%i,%i)\n", mapS[0], mapS[1], mapS[2]);
	int dM[3] = {2*mapS[0]+1, 2*mapS[1]+1, 2*mapS[2]+1};

	//printf("Max indexval: %i\n", dM[0]*dM[1]*dM[2]);

	lapack_complex_double* coefs;
	read_Coef("TaAs_ham_Hr_real.dat", "TaAs_ham_Hr_imag.dat", &coefs, n, nmap);
	
	//print a test value of Hmap
	//printf("Number of coef values: %i\n", nmap);
	//printf("index: %i\n", (-2+mapS[0])*dM[1]*dM[2]+(-4+mapS[1])*dM[2]+(0+mapS[2]));
	//printf("Test Entry (-2,-4,0): %i\n", Hmap[(-2+mapS[0])*dM[1]*dM[2]+(-4+mapS[1])*dM[2]+(0+mapS[2])] );
	

	nncounts = (int * ) malloc((nndepth+1) * sizeof(int ));
    if (nncounts == NULL) {
    	printf("Memory allocation error when allocating NN counter array.\n");
    	exit(0);
	}
	int* nncells;
	double* nnvecs;
	memset(nncounts, 0, (nndepth+1) * sizeof(int ));
	//printf("Test\n");
	//printf("Mapshifts: (%i,%i,%i)\n", mapS[0], mapS[1], mapS[2]);
	find_neighbors_cutoff(n, lat_vectors, atom_vectors, Hmap, mapS, 15.0, &nncells, &nnvecs, &nntotal );
	//find_neighbors(n, lat_vectors, atom_vectors, nndepth, &nn, nncounts);

	//printf("Test1\n");
	//printf("NN count: %i\n", nntotal );

	//Allocate Hamiltonian memory
	//x4 because spin in each dimension
	lapack_complex_double * H_tb;
	H_tb = (lapack_complex_double * ) malloc(n * n * sizeof(lapack_complex_double ));
    if (H_tb == NULL) {
    	printf("Memory allocation error when allocating H_tb.\n");
    	exit(0);
	}
	memset(H_tb, 0, n * n * sizeof(lapack_complex_double ));

	double * eig;	//Eigenenergies
	eig = (double * ) malloc(n * sizeof(double ));
    if (eig == NULL) {
    	printf("Memory allocation error when allocating EigenEnergy Array.\n");
    	exit(0);
	}
	memset(eig, 0, n * sizeof(double ));
	
	/*
	//quick band structure calc
	printf("Tes2\n");
	int numberk = 250;
	double ksx = 1.88799848;
	double ksy1 = 0.15173101;
	double ksy2 = 1.88189187;

	//[ 1.88799848 -0.          0.        ]
	//[ 0.15173101  1.88189187 -0.        ]
	//[ 0.30346257  0.27998014  0.98762127]

	
	ksx /= ((double)numberk);
	ksy1 /= ((double)numberk);
	ksy2 /= ((double)numberk);
	double testk[3] = {0.0, 0.0, 0.0};
	//for (i = 0; i < numberk*3; i++) {
	for (i = 0; i < 1; i++) {
		testk[0] = ((double)i)*(ksx+ksy1);
		testk[1] = ((double)i)*(ksy2);
		printf("Tes4\n");
		build_TaAsH(n, testk, nncells, nnvecs, nntotal, H_tb, coefs, Hmap, mapS );
		printf("Tes5\n");
		info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'L', n, H_tb, n, eig );
		if( info > 0 ) {
			printf( "The algorithm failed to compute eigenvalues.\n" );
			exit( 1 );
		}
		
		//printf("%9.6f\t%9.6f\t%9.6f\n",testk[0],testk[1],testk[2]);
		for(j=0; j<n; j++){
			printf("%9.6f\t", eig[j]);
		}
		printf("\n");

	}

	
	exit(0);
	*/
	
	
	
	//printf("NNcounts: (%i, %i )\n", nncounts[0], nncounts[1]);



	//Find total number of NN
	//nntotal = 0;
	//for( i = 0; i < nndepth+1; i++ ) {
	//	nntotal+=nncounts[i];
	//}
	//printf("Total NN: %i\n", nntotal);
	//print_rmatrix("Neighbors", nntotal, 5, nn, 5);

	/*
	//Allocate memory for coefficients (to be used in Hamlitonian later) x2 because spin.
	lapack_complex_double* coeftable;
	coeftable = (lapack_complex_double * ) malloc(nntotal * 2 * sizeof(lapack_complex_double));
    if (coeftable == NULL) {
    	printf("Memory allocation error when allocating Hamiltonian coefficient table.\n");
    	exit(0);
	}
	memset(coeftable, 0, nntotal * 2 * sizeof(lapack_complex_double));

	build_C( nndepth, nntotal, nn, nncounts, coeftable, graphene_H);


	*/
	
	
	
	//*******************
	// Allocate temp H arrays, initializing to 0
	//*******************

	lapack_complex_double * evc;	//Eigenvectors
	double * eec;	//Eigenenergies
	lapack_complex_double * dHdk;	//Derivatives
	lapack_complex_double * tempP; //Temp array to store inner product
	double * bC; //band berry curvature 
	evc = (lapack_complex_double * ) malloc(n * n  * sizeof(lapack_complex_double ));
    if (evc == NULL) {
    	printf("Memory allocation error when allocating H_tb.\n");
    	exit(0);
	}
	eec = (double * ) malloc(n  * sizeof(double ));
    if (eec == NULL) {
    	printf("Memory allocation error when allocating EigenEnergy Array.\n");
    	exit(0);
	}
	dHdk = (lapack_complex_double * ) malloc(n * n  * 3 * sizeof(lapack_complex_double ));
    if (dHdk == NULL) {
    	printf("Memory allocation error when allocating dHdk.\n");
    	exit(0);
	}
	tempP = (lapack_complex_double * ) malloc(n * n  * sizeof(lapack_complex_double ));
    if (tempP == NULL) {
    	printf("Memory allocation error when allocating temporary matrix product storage.\n");
    	exit(0);
	}
	bC = (double * ) malloc(3 * n  * sizeof(double ));
    if (bC == NULL) {
    	printf("Memory allocation error when allocating band Berry Curvature array.\n");
    	exit(0);
	}
	memset(bC, 0, 3 * n * sizeof(double ));
	memset(evc, 0, n * n * sizeof(lapack_complex_double ));
	memset(eec, 0, n * sizeof(double ));
	memset(dHdk, 0, n * n * 3 * sizeof(lapack_complex_double ));
	memset(tempP, 0, n * n  * sizeof(lapack_complex_double ));
	memset(eec, 0, 3 * n * sizeof(double ));

	//*******************
	// Compute step sizes for k-vectors
	//*******************

	//Tantalum Arsenide brillouin zone vectors
	// double bZvec[9] = { 1.03839916, -0.08372285, -0.29533053,
	// 					0.0        , 1.041769,   -0.29533053,
	// 					0.16690441, 0.15398908, 0.5431917  };


	double bZvec[9] = { 1.13279909, -0.09133402, -0.32217876,
						0.0        , 1.13647528, -0.32217876,
						0.18207754, 0.16798808, 0.59257276  };

	lapack_int nkp = 10; //k point density 

	//divide into nkp steps
	for(i = 0; i < 9; i++){
		bZvec[i] /= ((double) nkp);
	}

	double kstep[3];
	kstep[0] = sqrt( pow(bZvec[0],2) + pow(bZvec[1],2) + pow(bZvec[2],2) );
	kstep[1] = sqrt( pow(bZvec[3],2) + pow(bZvec[4],2) + pow(bZvec[5],2) );
	kstep[2] = sqrt( pow(bZvec[6],2) + pow(bZvec[7],2) + pow(bZvec[8],2) );
	// kstep[0] = 0.5*M_PI/sqrt( pow(lat_vectors[0],2) + pow(lat_vectors[1],2) + pow(lat_vectors[2],2) );
	// kstep[1] = 0.5*M_PI/sqrt( pow(lat_vectors[3],2) + pow(lat_vectors[4],2) + pow(lat_vectors[5],2) );
	// kstep[2] = 0.1*M_PI/sqrt( pow(lat_vectors[6],2) + pow(lat_vectors[7],2) + pow(lat_vectors[8],2) );

	//printf("x: %6.2f", sqrt( pow(lat_vectors[0],2) + pow(lat_vectors[1],2) + pow(lat_vectors[2],2) ));
	//printf("y: %6.2f", sqrt( pow(lat_vectors[3],2) + pow(lat_vectors[4],2) + pow(lat_vectors[5],2) ));
	//printf("z: %6.2f\n",sqrt( pow(lat_vectors[6],2) + pow(lat_vectors[7],2) + pow(lat_vectors[8],2) ));
	//divide into nkp steps
	// kstep[0] /= ((double) nkp);
	// kstep[1] /= ((double) nkp);
	// kstep[2] /= ((double) nkp);

	//printf("kvol: %9.6f\n", 1.0/(kstep[0]*kstep[1]*kstep[2]));
	//printf("kz step: %9.6f\n", kstep[2]);

	double kvc[3]; //Stores the current k-vector
	double berryC[3]; //Stores the current berry curvature
	lapack_complex_double alpha = {1.0,0.0}; //Matrix multiplication is unscaled
	lapack_complex_double beta = {0.0,0.0};	//Don't copy existing matrix
	lapack_int hamSize = n * n; //number of cells in Hamiltonian

	double tempBC = 0.0; //temp variable for berry curvature
	
	lapack_int numk = (2*nkp+2)*(2*nkp+2)*(2*nkp+2);

	//store all k points
	double * berryk;
	berryk = (double * ) malloc(3 * numk  * sizeof(double ));
    if (berryk == NULL) {
    	printf("Memory allocation error when allocating band Berry Curvature array.\n");
    	exit(0);
	}
	memset(berryk, 0, 3 * numk * sizeof(double ));
	//store all values
	double * berryv;
	berryv = (double * ) malloc(3 * n * numk  * sizeof(double ));
    if (berryv == NULL) {
    	printf("Memory allocation error when allocating band Berry Curvature array.\n");
    	exit(0);
	}
	memset(berryv, 0, 3 * n * numk * sizeof(double ));

	lapack_complex_double * berryevc;	//Eigenvectors at each k-point
	berryevc = (lapack_complex_double * ) malloc(numk * n * n  * sizeof(lapack_complex_double ));
    if (berryevc == NULL) {
    	printf("Memory allocation error when allocating berry eigenvector array.\n");
    	exit(0);
	}
	memset(berryevc, 0, numk * n * n * sizeof(lapack_complex_double ));



	//exit(0);
	
	
	// //do quick weyl calc
	// char filename[256];
	// readBerry("BC.out", berryv);
	// for (i = 0; i < n; i++) {
	// 	find_weyl3(nkp, berryv, kstep, bZvec, i, n);
	// 	sprintf(filename,"Berry/BC%02d.out",i);
	// 	//printf("%s\n",filename);
	// 	printBerryBand(filename, berryv, i,n, numk);
		

	// }
	// exit(0);

	clock_t allstart, end, weylend;
	double cpu_time_total;
	double avg_point_time;
	int time_h, time_m, time_s, time_r;
	
	allstart = clock();
	
	lapack_int berryi = 0;
	for (i = -nkp-1; i < nkp+1; i++) {
	//for (i = 220; i < 250; i++) {
		for (j = -nkp-1; j < nkp+1; j++) {
		//for (j = -20; j < 10; j++) {
			for (k = -nkp-1; k < nkp+1; k++) {
			//for (k = 0; k < 1; k++) {
			//for (k = nkp; k < nkp+1; k++) {

				//kvc[0] = ((double) i +0.5) * bZvec[0] + ((double) j+0.5) * bZvec[3] + ((double) k+0.5) * bZvec[6];
				//kvc[1] = ((double) i +0.5) * bZvec[1] + ((double) j+0.5) * bZvec[4] + ((double) k+0.5) * bZvec[7];
				//kvc[2] = ((double) i +0.5) * bZvec[2] + ((double) j+0.5) * bZvec[5] + ((double) k+0.5) * bZvec[8];

				// kvc[0] = ((double) i +0.0) * bZvec[0] + ((double) j+0.0) * bZvec[3] + ((double) k+0.0) * bZvec[6];
				// kvc[1] = ((double) i +0.0) * bZvec[1] + ((double) j+0.0) * bZvec[4] + ((double) k+0.0) * bZvec[7];
				//kvc[2] = ((double) i +0.0) * bZvec[2] + ((double) j+0.0) * bZvec[5] + ((double) k+0.0) * bZvec[8];

				kvc[0] = ((double) i) * bZvec[0] + ((double) j) * bZvec[3] + ((double) k) * bZvec[6] + 0.01;
				kvc[1] = ((double) i) * bZvec[1] + ((double) j) * bZvec[4] + ((double) k) * bZvec[7] + 0.02;
				kvc[2] = ((double) i) * bZvec[2] + ((double) j) * bZvec[5] + ((double) k) * bZvec[8] + 0.01;
				

				// kvc[0] = (double) i * kstep[0] + 0.005; //Shifted so no sym lines
				// kvc[1] = (double) j * kstep[1] + 0.005;
				// kvc[2] = (double) k * kstep[2] + 0.005;

				berryk[3*berryi+0] = kvc[0];
				berryk[3*berryi+1] = kvc[1];
				berryk[3*berryi+2] = kvc[2];
				
				//printf("\n*****\nk-vec: (%6.2f,%6.2f,%6.2f):\n", kvc[0],kvc[1],kvc[2]);
				
				memset(evc, 0, n * n * sizeof(lapack_complex_double ));
				build_TaAsH(n, kvc, nncells, nnvecs, nntotal, evc, coefs);
				
				memset(dHdk, 0, n * n * 3 * sizeof(lapack_complex_double ));
				build_TaAsdHdk(n, kvc, nncells, nnvecs, nntotal, dHdk, coefs);

				
				
				
			

				//build_testH(n, kvc, nn, nntotal, evc, coeftable );
				//build_testdHdk(n, kvc, nn, nntotal, dHdk, coeftable );
				//build_dHdk(n, kvc, nn, nntotal, dHdk, coeftable );

				//print_matrix("Hamiltonian: ", 2*n, 2*n, evc, 2*n);
				//print_matrix("dHdkx: ", n, n, dHdk, n);
				 //print_matrix("dHdky: ", n, n, &dHdk[hamSize], n);
				 //print_matrix("dHdkz: ", n, n, &dHdk[2*hamSize], n);

				memset(tempP, 0, n * n * sizeof(lapack_complex_double ));
				memset(eec, 0, n*sizeof(double));
				//Solve H to get eigenvalues
				info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'L', n, evc, n, eec );
				if( info > 0 ) {
					printf( "The algorithm failed to compute eigenvalues.\n" );
					exit( 1 );
				}
				//Save eigenvalues

				for (l = 0; l < n; l++) {
					for(m = 0; m < n; m++) {
						berryevc[n*n*berryi + l*n + m] = evc[l*n + m];
						//berryi gets incremented later
					}
				}


				//print_matrix("Eigenvectors: ", n, n, evc, n);
				//Multiply eigenvalues by each dHdk (might be able to do in one step?)
				cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, n, n, n, &alpha, 
					dHdk, n, evc, n, &beta, tempP, n);
				//print_matrix("first mult dHdkx: ", n, n, tempP, n);
				//get inner product by multiplying by complex conjugate of eigenvalues
				//put result back into dHdk
				cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, n, n, n, &alpha, 
					evc, n, tempP, n, &beta, dHdk, n);
				//repeat for Y (using same tempP array)
				cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, n, n, n, &alpha, 
					&dHdk[hamSize], n, evc, n, &beta, tempP, n);
				cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, n, n, n, &alpha, 
					evc, n, tempP, n, &beta, &dHdk[hamSize], n);
				//Repeat for Z
				cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, n, n, n, &alpha, 
					&dHdk[2*hamSize], n, evc, n, &beta, tempP, n);
				cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, n, n, n, &alpha, 
					evc, n, tempP, n, &beta, &dHdk[2*hamSize], n);
				
				//print_dvector("EigenEnergies: ", 2*n, eec);
				//printf("k-vec: (%6.2f,%6.2f,%6.2f):\n", kvc[0],kvc[1],kvc[2]);
				//print_matrix("Mult-X result: ", n, n, dHdk, n);
				// print_matrix("Mult-Y result: ", 2*n, 2*n, &dHdk[hamSize], 2*n);
				// print_matrix("Mult-Z result: ", 2*n, 2*n, &dHdk[2*hamSize], 2*n);

				//memset(berryC, 0, 3*sizeof(double));
				memset(bC, 0, 3*n*sizeof(double));
				
				for (l = 0; l < n; l++) {
					for(m = 0; m < n; m++) {
						//Make sure that n'=/=n , i.e l=/=m
						if( l != m ){
							//Also helps make sure that there is no division by 0
							//if( fabs(eec[l]-eec[m]) > 0.000001 ){
								//add small imaginary part to get rid of singularities
								//x => y * z*
								tempBC = cimag(dHdk[hamSize+l*n + m]*conj(dHdk[2*hamSize+l*n + m]));
								bC[l*3+0] += creal( tempBC/(pow((eec[l]-eec[m]),2) + 0.001*I) );
								
								//y => z * x*
								tempBC = cimag(dHdk[2*hamSize+l*n + m]*conj(dHdk[l*n + m]));
								bC[l*3+1] += creal( tempBC/(pow((eec[l]-eec[m]),2) + 0.001*I) );
								
								//z => x * y*
								tempBC = cimag(dHdk[l*n + m]*conj(dHdk[hamSize+l*n + m]));
								bC[l*3+2] += creal( tempBC/(pow((eec[l]-eec[m]),2) + 0.001*I) );
							//}
						} 
					}
				}

				
				// for (l = 0; l < n/2; l++) {
				// 	berryC[0] += bC[l*3+0];
				// 	berryC[1] += bC[l*3+1];
				// 	berryC[2] += bC[l*3+2];
				// }
				

				for (l = 0; l < n; l++) {
					berryv[3*n*berryi+3*l+0] = bC[l*3+0];
					berryv[3*n*berryi+3*l+1] = bC[l*3+1];
					berryv[3*n*berryi+3*l+2] = bC[l*3+2];
				}
				berryi ++;
				
				if(berryi%30 == 0){
					printf("Finished %i of %i. ", berryi, numk);
					weylend = clock();
					
					avg_point_time = ((double) (weylend-allstart)) / CLOCKS_PER_SEC;
					avg_point_time = (avg_point_time*((double)numk) )/((double) berryi);
					time_h = ((int) avg_point_time) / 3600;
					time_r = ((int) avg_point_time) % 3600;
					time_m = time_r / 60;
					time_s = time_r % 60;

					printf("Est. time: %ih %im %is \n", time_h, time_m, time_s);

				}



				
			}
		}
	}

	end = clock();
	cpu_time_total = ((double) (end-allstart)) / CLOCKS_PER_SEC;

	printf("FINISHED: Computed %i points in %6.2f seconds.\n", numk, cpu_time_total);

	find_weyl4(nkp, berryevc, n, 16, berryv);
	printf("TEST4 %i\n", berryi);

	printBerry("kp.out", berryk, berryi);
	printBerry("BC.out", berryv, berryi*n);



	//find_weyl3(nkp, berryv, kstep);

	exit(0);

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