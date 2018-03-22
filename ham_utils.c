#include <lapacke.h>
#include <cblas.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ham_utils.h"

void find_neighbors( lapack_int n, double* lat_vecs, double* atom_vecs,
	lapack_int nndepth, double ** nn, int* nn_counts ){

	lapack_int i ,j, k, l, m;

	double a[3] = {0.0, 0.0, 0.0}, b[3];
	double len;

	lapack_int rownum = 0;

	//Temporary neighbor array to store all possible NN.
	double * tempnn;

	//Want to search all adjacent unit cells nndepth deep
	lapack_int dim = (2*nndepth + 1); 
	lapack_int vol = dim*dim*dim;

	// printf("dim: %i\n", dim);
	// printf("vol: %i\n", vol);
	// printf("n: %i\n", n);
	// printf("nndepth: %i\n", nndepth);

	/*	Need 5 pieces of data in our array: 
	 *	1. Start atom;	2. End atom;	3-5. x, y, z hopping vector
	 *	
	 *	- Store atom numbers as double to avoid multiple arrays.
	 *	- The number of such entries is n * n * vol (all n*n pairs spanning vol cells).
	*/

	//	Initialize NN array.
	tempnn =  malloc( n * n * vol * 5 * sizeof(double ));
	if (tempnn == NULL) {
		printf("Memory allocation error when allocating atom vector array.\n");
		exit(0);
	}
	memset(tempnn, 0, n * n * vol * 5 * sizeof(double ));
	
	//loop through each atom in unit cell (or supercell!)
	for( m = 0; m < n; m++ ) {
		for( l = 0; l < n; l++ ) {
			//Generate adjacent unit cells (nndepth unit cells deep)
			for( i = -nndepth; i < nndepth+1; i++ ) {
				for( j = -nndepth; j < nndepth+1; j++ ) {
					for( k = -nndepth; k < nndepth+1; k++ ) {
						
						/*	Find the array index into which to write the nn pair, 
						   	and the vector connecting them. */
						rownum = m*n*vol + l*vol + (i+nndepth)*dim*dim+(j+nndepth)*dim+(k+nndepth);

						//printf("rownum: %i\n", rownum);
						
						tempnn[rownum*5]   = (double) (m+1); 		//Start atom
						tempnn[rownum*5+1] = (double) (l+1); 		//End   atom

						/*	{0,3,6} contain the x vals; {1,4,7} the y vals; {2,5,8} the z vals.
						 *	multiply by {i,j,k} to find the vector connecting unit cells,
						 *  then add the difference in atom locations (within unit cell)
						 *	to find the total vector.
						 */ 
						tempnn[rownum*5+2] = i*lat_vecs[0]+j*lat_vecs[3]+k*lat_vecs[6]+atom_vecs[l*3+0]-atom_vecs[m*3+0];
						tempnn[rownum*5+3] = i*lat_vecs[1]+j*lat_vecs[4]+k*lat_vecs[7]+atom_vecs[l*3+1]-atom_vecs[m*3+1];
						tempnn[rownum*5+4] = i*lat_vecs[2]+j*lat_vecs[5]+k*lat_vecs[8]+atom_vecs[l*3+2]-atom_vecs[m*3+2];
					}
				}
			}
		}
	}

	//Sort array (see compare_nn for the ordering)
	qsort((void*)tempnn, n*n*vol, 5*sizeof(double), compare_nn);


	//print_rmatrix("Neighbors", 34, 5, tempnn, 5);
	double nnvecL = 0.0;
	double nnvecLtemp = 0.0;
	//int nntotal = 0;
	j = -1;
	i = 0;

	while( 1 ){
		j++;
		//printf("Finding length of (%6.4f,%6.4f,%6.4f)\n", tempnn[5*j+2],tempnn[5*j+3],tempnn[5*j+4]);
		nnvecLtemp = pow(tempnn[5*j+2],2) + pow(tempnn[5*j+3],2) + pow(tempnn[5*j+4],2);
		//printf("Length  = %6.4f \n", nnvecLtemp);

		if ((nnvecLtemp - nnvecL) > 0.01 ){
			nnvecL = nnvecLtemp;
			i++;
			if ( i > nndepth)
			{
				break;
			} 
		}
		nn_counts[i] += 1;
		//printf("Incremented nn_counts[%i] to %i \n", i, nn_counts[i]);
		//nntotal ++;
	}
	//printf("NNcounts: (%i, %i, %i )\n", nn_counts[0], nn_counts[1], nn_counts[2]);
	//printf("NNtotal: %i \n", j);

	//Allocate trimmed NN array
	(*nn) =  malloc( j * 5 * sizeof(double ));
	if ((*nn) == NULL) {
		printf("Memory allocation error when allocating atom vector array.\n");
		exit(0);
	}
	memset((*nn), 0, j * 5 * sizeof(double ));
	for( i= 0; i < j * 5 ; i++ ) {
		(*nn)[i] = tempnn[i];
	}

	double testvec[3];
	v_perp((*nn), nn_counts, 1, testvec);
	//printf("\n\nPerpTest = (%6.2f, %6.2f, %6.2f)\n", testvec[0], testvec[1], testvec[2]);

}

void find_neighbors_cutoff( lapack_int n, double* lat_vecs , double* atom_vecs , int* H_map, lapack_int* map_S 
    ,double cutoff, lapack_int ** nn, double ** nn_vecs, lapack_int* nn_total ){

	lapack_int i ,j, k, l, m;

	//Temporary neighbor array to store all possible NN.
	lapack_int * tempnn;
	double * tempnn_vecs;

	//find how many unit cells to go in each direction
	//need to figure out a way to do this more intelligently
	lapack_int x_cutoff = 8; //ceil(cutoff/lat_vecs[0])
	lapack_int y_cutoff = 8;
	lapack_int z_cutoff = 6;

	lapack_int dimx = 2*x_cutoff+1;
	lapack_int dimy = 2*y_cutoff+1;
	lapack_int dimz = 2*z_cutoff+1;
	lapack_int vol = dimx*dimy*dimz;

	/*	Need 5 pieces of data in our array: 
	 *	1. Start atom;	2. End atom;	3-5. unit cell of end atom, N_x, N_y, N_z.
	 *	
	 *	- The number of such entries is n * n * vol (all n*n pairs spanning vol cells).
	*/

	//	Initialize NN array.
	tempnn =  malloc( n * n * vol * 3 * sizeof(lapack_int ));
	if (tempnn == NULL) {
		printf("Memory allocation error when allocating temp nn array.\n");
		exit(0);
	}
	tempnn_vecs = malloc( n * n * vol * 3 * sizeof(double));
	if (tempnn_vecs == NULL) {
		printf("Memory allocation error when allocating temp nn hopping vector array.\n");
		exit(0);
	}
	memset(tempnn, 0, n * n * vol * 3 * sizeof(lapack_int ));
	memset(tempnn_vecs, 0, n * n * vol * 3 * sizeof(double));

	double hop_V[3]={0.0,0.0,0.0}; 	//Hopping vector from start to end atom
	double hop_L=0.0;				//Hopping vector length

	int dM[3] = {2*map_S[0]+1, 2*map_S[1]+1, 2*map_S[2]+1}; //lookup table dimensions
	int indexval = 0;
	

	lapack_int rownum = 0;
	//loop through each atom in unit cell (or supercell!)
	for( m = 0; m < n; m++ ) {
		for( l = 0; l < n; l++ ) {
			//Generate adjacent unit cells (nndepth unit cells deep)
			for( i = -x_cutoff ; i < x_cutoff+1; i++ ) {
				for( j = -y_cutoff ; j < y_cutoff +1; j++ ) {
					for( k = -z_cutoff ; k < z_cutoff +1; k++ ) {
						
						/*	{0,3,6} contain the x vals; {1,4,7} the y vals; 
						 *	{2,5,8} the z vals. Multiply by {i,j,k} to find
						 *	the vector connecting unit cells, then add the 
						 * 	difference in atom locations (within unit cell)
						 *	to find the total vector.
						 */ 

						hop_V[0] = 	((double) i ) * lat_vecs[0]
									+((double) j ) * lat_vecs[3]
									+((double) k ) * lat_vecs[6]
									+atom_vecs[l*3+0]-atom_vecs[m*3+0];
						hop_V[1] = 	((double) i ) * lat_vecs[1]
									+((double) j ) * lat_vecs[4]
									+((double) k ) * lat_vecs[7]
									+atom_vecs[l*3+1]-atom_vecs[m*3+1];
						hop_V[2] = 	((double) i ) * lat_vecs[2]
									+((double) j ) * lat_vecs[5]
									+((double) k ) * lat_vecs[8]
									+atom_vecs[l*3+2]-atom_vecs[m*3+2];
						
						//Find length
						hop_L = sqrt(pow(hop_V[0],2) + pow(hop_V[1],2)
								+ pow(hop_V[2],2));

						if(hop_L < cutoff){
							//record the current pair
							tempnn[rownum*3]   = (m); 		//Start atom
							tempnn[rownum*3+1] = (l); 		//End   atom
							//find index for coefficient lookup table
							indexval = 	(i+map_S[0])*dM[1]*dM[2]+
										(j+map_S[1])*dM[2]+
										(k+map_S[2]);
							
							indexval = (n*H_map[indexval]+m)*n+l;
							tempnn[rownum*3+2] = indexval; //save in nn array for convenience

							tempnn_vecs[rownum*3+0] = hop_V[0];
							tempnn_vecs[rownum*3+1] = hop_V[1];
							tempnn_vecs[rownum*3+2] = hop_V[2];
							
							//increment rownum
							rownum += 1;
						}
					
					}
				}
			}
		}
	}

	//printf("rownum: %i\n", rownum);
	(*nn_total) = rownum;

	//Allocate trimmed NN array
	(*nn) =  malloc( rownum * 3 * sizeof(lapack_int ));
	if ((*nn) == NULL) {
		printf("Memory allocation error when allocating nn array.\n");
		exit(0);
	}
	memset((*nn), 0, rownum * 3 * sizeof(lapack_int ));
	//Allocate trimmed NN hopping vector array
	(*nn_vecs) =  malloc( rownum * 3 * sizeof(double ));
	if ((*nn_vecs) == NULL) {
		printf("Memory allocation error when allocating nn hopping vector array.\n");
		exit(0);
	}
	memset((*nn_vecs), 0, rownum * 3 * sizeof(double ));

	for( i= 0; i < rownum * 3 ; i++ ) {
		(*nn_vecs)[i] = tempnn_vecs[i];
		(*nn)[i] = tempnn[i];
	}
}

void find_weyl( lapack_int nkp, double* bC, double*  kdif){
	
	lapack_int s = 2*nkp+2; //array size
	lapack_int i, j, k;
	lapack_int x, y, z;

	double kvc[3] = {0.0,0.0,0.0};

	double tempbC[3] = {0.0,0.0,0.0};

	double bCsum = 0.0;

	lapack_int weylcount = 0;

	//double kdif[3] = {0.0,0.0,0.0};
	for (i = -nkp+1; i < nkp-1; i++) {
		for (j = -nkp+1; j < nkp-1; j++) {
			for (k = -nkp+1; k < nkp-1; k++) {
				x = i + nkp + 1;
				y = j + nkp + 1;
				z = k + nkp + 1;

				tempbC[0] =  1.0*bC[3*((x-2)*s*s+(y+0)*s+(z+0))+0]
						    -8.0*bC[3*((x-1)*s*s+(y+0)*s+(z+0))+0] 
							+8.0*bC[3*((x+1)*s*s+(y+0)*s+(z+0))+0]
							-1.0*bC[3*((x+2)*s*s+(y+0)*s+(z+0))+0];

				tempbC[1] =  1.0*bC[3*((x+0)*s*s+(y-2)*s+(z+0))+1]
						    -8.0*bC[3*((x+0)*s*s+(y-1)*s+(z+0))+1] 
							+8.0*bC[3*((x+0)*s*s+(y+1)*s+(z+0))+1]
							-1.0*bC[3*((x+0)*s*s+(y+2)*s+(z+0))+1];

				tempbC[2] =  1.0*bC[3*((x+0)*s*s+(y+0)*s+(z-2))+2]
						    -8.0*bC[3*((x+0)*s*s+(y+0)*s+(z-1))+2] 
							+8.0*bC[3*((x+0)*s*s+(y+0)*s+(z+1))+2]
							-1.0*bC[3*((x+0)*s*s+(y+0)*s+(z+2))+2];

				tempbC[0]/= (12.0*kdif[0]);
				tempbC[1]/= (12.0*kdif[1]);
				tempbC[2]/= (12.0*kdif[2]);

				bCsum = (tempbC[0] + tempbC[1] + tempbC[2]);//*kdif[0]*kdif[1]*kdif[2];

				kvc[0] = (double) i * kdif[0] + 0.025; //Shifted so no sym lines
				kvc[1] = (double) j * kdif[1] + 0.025;
				kvc[2] = (double) k * kdif[2] + 0.025;
				if( fabs(bCsum) >= 1.0){
				//if( (fabs(bCsum) >= 1.0) && (fabs(kvc[2])<0.1) ){
					printf("Weyl%i: (%6.2f,%6.2f,%6.2f): %6.2f\n", weylcount, kvc[0],kvc[1],kvc[2], bCsum);
					//printf("%9.6f\t%9.6f\t%9.6f\t%9.6f\n",kvc[0],kvc[1],kvc[2], bCsum);
					weylcount ++;
				}
				
			}
		}
	}

}

void find_weyl2( lapack_int nkp, double* bC, double*  kdif){
	
	lapack_int s = 2*nkp+2; //array size
	lapack_int i, j, k;
	lapack_int x, y, z;

	double kvc[3] = {0.0,0.0,0.0};

	double tempbC[3] = {0.0,0.0,0.0};

	double bCsum = 0.0;

	lapack_int weylcount = 0;

	lapack_int m,n;

	//double kdif[3] = {0.0,0.0,0.0};
	for (i = -nkp+1; i < nkp-1; i++) {
		for (j = -nkp+1; j < nkp-1; j++) {
			for (k = -nkp+1; k < nkp-1; k++) {
				x = i + nkp + 1;
				y = j + nkp + 1;
				z = k + nkp + 1;

				//top and bottom face
				for (m = -2; m < 3; m++) {
					for (n = -2; n < 3; n++) {
						bCsum+=bC[3*((x+m)*s*s+(y+n)*s+(z+2))+2]*kdif[0]*kdif[1];
						bCsum-=bC[3*((x+m)*s*s+(y+n)*s+(z-2))+2]*kdif[0]*kdif[1];
					}
				}

				//xz
				for (m = -2; m < 3; m++) {
					for (n = -2; n < 3; n++) {
						bCsum+=bC[3*((x+m)*s*s+(y+2)*s+(z+n))+1]*kdif[0]*kdif[2];
						bCsum-=bC[3*((x+m)*s*s+(y-2)*s+(z+n))+1]*kdif[0]*kdif[2];
					}
				}

				//yz
				for (m = -2; m < 3; m++) {
					for (n = -2; n < 3; n++) {
						bCsum+=bC[3*((x+2)*s*s+(y+m)*s+(z+n))+0]*kdif[1]*kdif[2];
						bCsum-=bC[3*((x-2)*s*s+(y+m)*s+(z+n))+0]*kdif[1]*kdif[2];
					}
				}

				kvc[0] = (double) i * kdif[0] + 0.03; //Shifted so no sym lines
				kvc[1] = (double) j * kdif[1] + 0.01;
				kvc[2] = (double) k * kdif[2] + 0.02;

				
				if( fabs(bCsum) >= 1.0){
				//if( (fabs(bCsum) >= 1.0) && (fabs(kvc[2])<0.1) ){
					printf("Weyl%i: (%6.2f,%6.2f,%6.2f): %6.2f\n", weylcount, kvc[0],kvc[1],kvc[2], bCsum);
					//printf("%9.6f\t%9.6f\t%9.6f\n",kvc[0],kvc[1],kvc[2]);
					weylcount ++;
				}

				//if weylcount = 85
				
			}
		}
	}

}

void find_weyl3( lapack_int nkp, double* bC, double*  kdif, double* bZvec, lapack_int bN, lapack_int btotal){
	
	// bN is band number
	// btotal is total number of bands
	lapack_int s = 2*nkp+2; //array size
	lapack_int i, j, k;
	lapack_int x, y, z;

	double kvc[3] = {0.0,0.0,0.0};

	double tempbC[3] = {0.0,0.0,0.0};

	double bCsum = 0.0;

	lapack_int weylcount = 0;

	double* divBC;
	divBC = (double * ) malloc( 8 * nkp * nkp * nkp * sizeof(double ));
    if (divBC == NULL) {
    	printf("Memory allocation error when allocating temp div Berry Curvature array.\n");
    	exit(0);
	}
	memset(divBC, 0, 8 * nkp * nkp * nkp * sizeof(double ));

	//double kdif[3] = {0.0,0.0,0.0};
	for (i = -nkp; i < nkp; i++) {
		for (j = -nkp; j < nkp; j++) {
			for (k = -nkp; k < nkp; k++) {
				x = i + nkp + 1;
				y = j + nkp + 1;
				z = k + nkp + 1;

				tempbC[0] = -1.0*bC[3*bN + 3*btotal*((x-1)*s*s+(y+0)*s+(z+0))+0] 
							+1.0*bC[3*bN + 3*btotal*((x+1)*s*s+(y+0)*s+(z+0))+0];

				tempbC[1] = -1.0*bC[3*bN + 3*btotal*((x+0)*s*s+(y-1)*s+(z+0))+1] 
							+1.0*bC[3*bN + 3*btotal*((x+0)*s*s+(y+1)*s+(z+0))+1];
							
				tempbC[2] = -1.0*bC[3*bN + 3*btotal*((x+0)*s*s+(y+0)*s+(z-1))+2] 
							+1.0*bC[3*bN + 3*btotal*((x+0)*s*s+(y+0)*s+(z+1))+2];

				tempbC[0]/= (2.0*kdif[0]);
				tempbC[1]/= (2.0*kdif[1]);
				tempbC[2]/= (2.0*kdif[2]);

				divBC[weylcount] = (tempbC[0] + tempbC[1] + tempbC[2])*kdif[0]*kdif[1]*kdif[2];

				weylcount ++;
				/*
				kvc[0] = (double) i * kdif[0] + 0.005; //Shifted so no sym lines
				kvc[1] = (double) j * kdif[1] + 0.005;
				kvc[2] = (double) k * kdif[2] + 0.005;
				if( fabs(bCsum) >= 1.0){
				//if( (fabs(bCsum) >= 1.0) && (fabs(kvc[2])<0.1) ){
					printf("Weyl%i: (%6.2f,%6.2f,%6.2f): %6.2f\n", weylcount, kvc[0],kvc[1],kvc[2], bCsum);
					//printf("%9.6f\t%9.6f\t%9.6f\t%9.6f\n",kvc[0],kvc[1],kvc[2], );
					
				}*/
				
			}
		}
	}

	//TaAs vectors for finding kpoints:
	/*double bZvec[9] = { 0.05191996, -0.00418614, -0.01476653,
						0.0       ,  0.05208845, -0.01476653,
						0.00834522, 0.00769945, 0.02715958 };*/

	//sum over 3x3x3 cubes
	weylcount = 0;
	s = 2*nkp;
	lapack_int m,n,o;
	for (i = -nkp+1; i < nkp-1; i++) {
		for (j = -nkp+1; j < nkp-1; j++) {
			for (k = -nkp+1; k < nkp-1; k++) {
				x = i + nkp;
				y = j + nkp;
				z = k + nkp;

				bCsum = 0.0;
				for (m = -1; m < 2; m++) {
					for (n = -1; n < 2; n++) {
						for (o = -1; o < 2; o++) {
							bCsum += divBC[(x+m)*s*s+(y+n)*s+(z+o)];
						}
					}
				}


				kvc[0] = ((double) i) * bZvec[0] + ((double) j) * bZvec[3] + ((double) k) * bZvec[6] + 0.01;
				kvc[1] = ((double) i) * bZvec[1] + ((double) j) * bZvec[4] + ((double) k) * bZvec[7] + 0.02;
				kvc[2] = ((double) i) * bZvec[2] + ((double) j) * bZvec[5] + ((double) k) * bZvec[8] + 0.01;
				
				
				//kvc[0] = ((double) i +0.5) * bZvec[0] + ((double) j+0.5) * bZvec[3] + ((double) k+0.5) * bZvec[6];
				//kvc[1] = ((double) i +0.5) * bZvec[1] + ((double) j+0.5) * bZvec[4] + ((double) k+0.5) * bZvec[7];
				//kvc[2] = ((double) i +0.5) * bZvec[2] + ((double) j+0.5) * bZvec[5] + ((double) k+0.5) * bZvec[8];
				
				// if (fabs(kvc[2]) > 0.8){
				// 	printf("i,j,k: (%i, %i, %i)\n", i, j, k);
				// }

				// kvc[0] = (double) i * kdif[0] + 0.005; //Shifted so no sym lines
				// kvc[1] = (double) j * kdif[1] + 0.005;
				// kvc[2] = (double) k * kdif[2] + 0.005;
				if( fabs(bCsum) >= 0.05){
				//if( (fabs(bCsum) >= 1.0) && (fabs(kvc[2])<0.1) ){
					//printf("Band%i, weyl#%i: (%6.2f,%6.2f,%6.2f): %9.6f\n", bN, weylcount, kvc[0],kvc[1],kvc[2], bCsum);
					printf("%9.6f\t%9.6f\t%9.6f\t%9.6f\n",kvc[0],kvc[1],kvc[2],bCsum);
					weylcount ++;
				}
			}
		}
	}

}

void find_weyl4( lapack_int nkp, lapack_complex_double * evc, lapack_int n , lapack_int bN, double* bC){
	// bN is band number
	// btotal is total number of bands
	lapack_int s = 2*nkp+2; //array size
	lapack_int i, j, k , l;
	lapack_int x, y, z;
	lapack_int kindexi; //array location corresponding to starting k-point
	lapack_int kindexf; //array location corresponding to ending k-point
	lapack_int BCcount = 0;

	//Multiplication consts
	lapack_complex_double alpha = {1.0,0.0}; //Matrix multiplication is unscaled
	lapack_complex_double beta = {0.0,0.0};	//Don't copy existing matrix

	//Compute link variables (3 per k point, per band)
	lapack_complex_double * tempP; //Temp array to store inner product
	tempP = (lapack_complex_double * ) malloc(n * n  * sizeof(lapack_complex_double ));
    if (tempP == NULL) {
    	printf("Memory allocation error when allocating temporary matrix product storage.\n");
    	exit(0);
	}
	memset(tempP, 0, n * n * sizeof(lapack_complex_double ));

	lapack_complex_double * berryc; //Temp array to store inner product
	berryc = (lapack_complex_double * ) malloc(3*n*8*nkp*nkp*nkp* sizeof(lapack_complex_double ));
    if (berryc == NULL) {
    	printf("Memory allocation error when allocating temporary matrix product storage.\n");
    	exit(0);
	}
	memset(berryc, 0, 3*n*8*nkp*nkp*nkp* sizeof(lapack_complex_double ));
	printf("\nTest0\n");

	for (i = -nkp; i < nkp; i++) {
		for (j = -nkp; j < nkp; j++) {
			for (k = -nkp; k < nkp; k++) {

				x = i + nkp + 1;
				y = j + nkp + 1;
				z = k + nkp + 1;

				//x direction first
				kindexi = (x)*s*s+(y)*s+(z);
				kindexf = (x+1)*s*s+(y)*s+(z);
				

				cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, n, n, n, &alpha, 
								&evc[kindexi*n*n], n, &evc[kindexf*n*n], n, &beta, tempP, n);
				for (l = 0; l < n; l++) {
					berryc[3*BCcount*n+l]=tempP[l*n+l]/cabs(tempP[l*n+l]);
				}
				//y dir
				kindexf = (x)*s*s+(y+1)*s+(z);
				cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, n, n, n, &alpha, 
								&evc[kindexi*n*n], n, &evc[kindexf*n*n], n, &beta, tempP, n);
				for (l = 0; l < n; l++) {
					berryc[3*BCcount*n+n+l]=tempP[l*n+l]/cabs(tempP[l*n+l]);
				}
				//z dir
				kindexf = (x)*s*s+(y)*s+(z+1);
				cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, n, n, n, &alpha, 
								&evc[kindexi*n*n], n, &evc[kindexf*n*n], n, &beta, tempP, n);
				for (l = 0; l < n; l++) {
					berryc[3*BCcount*n+2*n+l]=tempP[l*n+l]/cabs(tempP[l*n+l]);
					// printf("link var: (%9.6f,%9.6f)\t(%9.6f,%9.6f)\t(%9.6f,%9.6f)\n",
					// 	creal(berryc[3*BCcount*n+l]), cimag(berryc[3*BCcount*n+l]), 
					// 	creal(berryc[3*BCcount*n+n+l]), cimag(berryc[3*BCcount*n+n+l]), 
					// 	creal(berryc[3*BCcount*n+2*n+l]), cimag(berryc[3*BCcount*n+2*n+l]) );
				}
				BCcount++;
			}
		}
	}

	printf("Test1\n");

	lapack_complex_double tempBC[3] = { {0.0,0.0},{0.0,0.0},{0.0,0.0} };
	s = 2*nkp; // redefine size of matrix

	memset(bC, 0, 3*n*(2*nkp+2)*(2*nkp+2)*(2*nkp+2)* sizeof(double ));

	printf("Test2\n");

	for (i = -nkp; i < nkp-1; i++) {
		for (j = -nkp; j < nkp-1; j++) {
			for (k = -nkp; k < nkp-1; k++) {
				for (l = 0; l < n; l++) {
					x = i + nkp;
					y = j + nkp;
					z = k + nkp;
					//printf("kx: (%i/%i), ky: (%i/%i), kz: (%i/%i)\n", x, 2*nkp, y, 2*nkp, z, 2*nkp);

					kindexi = (x)*s*s+(y)*s+(z); //this is now index for berryc
					kindexf = (x)*s*s+(y+1)*s+(z);

					//printf("kindexi: %i\t, kindexf: %i\n", kindexi, kindexf);
					
					tempBC[0] = berryc[3*kindexi*n+n+l]*berryc[3*kindexf*n+2*n+l]/berryc[3*kindexi*n+2*n+l];
					kindexf = (x)*s*s+(y)*s+(z+1);
					//printf("kindexi: %i\t, kindexf: %i\n", kindexi, kindexf);
					tempBC[0] = tempBC[0]/berryc[3*kindexf*n+n+l];

					tempBC[1] = berryc[3*kindexi*n+2*n+l]*berryc[3*kindexf*n+l]/berryc[3*kindexi*n+l];
					kindexf = (x+1)*s*s+(y)*s+(z);
					//printf("kindexi: %i\t, kindexf: %i\n", kindexi, kindexf);
					tempBC[1] = tempBC[1]/berryc[3*kindexf*n+2*n+l];

					tempBC[2] = berryc[3*kindexi*n+l]*berryc[3*kindexf*n+n+l]/berryc[3*kindexi*n+n+l];
					kindexf = (x)*s*s+(y+1)*s+(z);
					//printf("kx: (%i/%i), ky: (%i/%i), kz: (%i/%i)\n", x, 2*nkp, y, 2*nkp, z, 2*nkp);
					//printf("kindexi: %i\t, kindexf: %i\n", kindexi, kindexf);
					tempBC[2] = tempBC[2]/berryc[3*kindexf*n+l];

					//tempBC[0] = clog(tempBC[0]);
					//tempBC[1] = clog(tempBC[1]);
					//tempBC[2] = clog(tempBC[2]);

					bC[3*n*kindexi+3*l+0] = cimag(clog(tempBC[0]))/M_PI;
					bC[3*n*kindexi+3*l+1] = cimag(clog(tempBC[1]))/M_PI;
					bC[3*n*kindexi+3*l+2] = cimag(clog(tempBC[2]))/M_PI;

					printf("BC: (%9.6f)\t(%9.6f)\t(%9.6f)\n",
							bC[3*n*kindexi+3*l+0], 
							bC[3*n*kindexi+3*l+1], 
							bC[3*n*kindexi+3*l+2] );
				
					// printf("BC: (%9.6f,%9.6f)\t(%9.6f,%9.6f)\t(%9.6f,%9.6f)\n",
					// 		creal(tempBC[0]), cimag(tempBC[0]), 
					// 		creal(tempBC[1]), cimag(tempBC[1]), 
					// 		creal(tempBC[2]), cimag(tempBC[2]) );

					// //M_PI
				}
			}
		}
	}


	// printf("BC: (%9.6f,%9.6f)\t(%9.6f,%9.6f)\t(%9.6f,%9.6f)\n",
	// 	creal(berryc[0]), cimag(berryc[0]), 
	// 	creal(berryc[n]), cimag(berryc[n]), 
	// 	creal(berryc[2*n]), cimag(berryc[2*n]));



	
}

/* Note: k must have 3 elements!
*/
void build_C(lapack_int nn_depth, int nn_total, double* nn, 
	int* nn_counts, lapack_complex_double * coef_table,
	lapack_complex_double (*ham)(lapack_int, double*, double*, lapack_complex_double*)){
	//Compute lookup table for coefficients
	lapack_int i, j;
	//lapack_int nn_total = 0;

	//for( i = 0; i < nn_depth+1; i++ ) {
	//	nn_total += nn_counts[i];
	//}

	double  up[3] 	= {0.0, 0.0, 1.0},
			down[3] = {0.0, 0.0,-1.0};
	/*
	lapack_complex_double test23 = {0.0, 0.0};
	
	double spintest[3] = {0.0, 0.0, 1.0};
	double vectest[3] = {0.0, 0.5, 0.866};
	printf( "testpre: (%6.2f, %6.2f)\n\n", creal(test23), cimag(test23));
	test23 = ham(2, spintest, vectest, coeff);
	printf( "testpost: (%6.4f, %6.4f)\n\n", creal(test23), cimag(test23));
	*/

	//TODO!!! Remove hardcoded coefficients!
	lapack_complex_double coeff[3] = {{0.0,0.0}, {1.0,0.0},{0.01,0.0}};
	lapack_int curdepth = 0;
	lapack_int nntemp = nn_counts[curdepth];
	//printf("TotalNN: %i\n", nn_total);
	double vperp[3];
	//double * nnpointer;
	//nn_total = nn_counts[curdepth];
	for( i = 0; i < nn_total; i++ ) {
		if (i >= nntemp){
			curdepth ++;
			nntemp += nn_counts[curdepth];
		}
		//nnpointer = &nn[5*i+2];
		//printf("Curdepth: %i\n", curdepth);
		//printf("Input vector: (%6.4f, %6.4f, %6.4f)\n", nnpointer[0], nnpointer[1], nnpointer[2]);
		v_perp(nn, nn_counts, i, vperp);
		//printf("Input vector: (%6.4f, %6.4f, %6.4f)\n", nnpointer[0], nnpointer[1], nnpointer[2]);
		//printf("Vperp: (%6.4f, %6.4f, %6.4f)\n", vperp[0], vperp[1], vperp[2]);
		coef_table[2*i] 	= ham(curdepth, up, vperp, coeff);
		coef_table[2*i+1] 	= ham(curdepth, down, vperp, coeff);

		//cblas_zdotu_sub(3, up, 1, &nn[5*i+2], 1, &test22);
		//printf("**********\nDone!!!\n");
		//printf( "Coeff%i: (%6.4f, %6.4f)\n", i+1, creal(coef_table[2*i]), cimag(coef_table[2*i]));
		
	}

	//Generate Hamiltonian
	//First do all spin up, then do all spin down	
}

void build_H(lapack_int H_size, double* k, double* nn, int nn_total, 
	lapack_complex_double* H_tb, lapack_complex_double* coef_table ){
	
	//int nn_depth_test = sizeof(*nn)/sizeof(double*);

	//printf("Size of nn: %i\n\n", H_size);
	//printf("\n***********\n");
	//printf("Other k: (%6.4f, %6.4f, %6.4f)\n", k[0], k[1], k[3]);
	lapack_int i, j;
	lapack_int x, y;
	lapack_complex_double phase = 0.0;
	//Assume H is initialized to all 0's
	for( i = 0; i < nn_total; i++ ) {
		//printf("\nx,y Coords: (%i, %i)\n", x , y);
		x = (int) nn[5*i+0];
		y = (int) nn[5*i+1];
		x--;
		y--;
		phase = k[0]*nn[5*i+2] + k[1]*nn[5*i+3] + k[2]*nn[5*i+4];
		//printf( "Angle %i: (%6.4f, %6.4f)\n", i+1, creal(phase), cimag(phase));
		phase = cexp(- I * phase);
		//printf( "Phase %i: (%6.4f, %6.4f)\n", i+1, creal(phase), cimag(phase));
		//printf( "Coeff %i: (%6.4f, %6.4f)\n", i+1, creal(coef_table[2*i]), cimag(coef_table[2*i]));

		
		//Spin-up
		H_tb[x*2*H_size + y] += coef_table[2*i]*phase;
		//printf("Coords (up): (%i, %i)\n", x*2*H_size , y);
		//printf("Coefficient: (%6.4f, %6.4f) \n", creal(coef_table[2*i]), cimag(coef_table[2*i]));
		//Spin-down
		H_tb[(x+H_size)*2*H_size + H_size + y] += coef_table[2*i+1]*phase;
	}


}

void build_TaAsH(lapack_int H_size, double* k, int* nn, double* nn_hops, int nn_total, 
	lapack_complex_double* H_tb, lapack_complex_double* c_table){
	
	lapack_int i;
	lapack_complex_double phase;// coef;
	//Set Hamiltonian to all 0's
	//memset(H_tb, 0, H_size * H_size * sizeof(lapack_complex_double ));
	
	//printf("nntotal: %i\n", nn_total);
	for( i = 0; i < nn_total; i++ ) {
		//For each nn pair
		//Compute phase:
		//phase = k[0]*nn_hops[3*i] + k[1]*nn_hops[3*i+1] + k[2]*nn_hops[3*i+2];
		//phase = cexp(- I * phase);

		phase = cexp(- I * (k[0]*nn_hops[3*i] + k[1]*nn_hops[3*i+1] + k[2]*nn_hops[3*i+2]) );
		//lookup coefficient
		//coef = c_table[ nn[3*i+2] ];

		//H_tb[nn[3*i]*H_size + nn[3*i+1]] += coef*phase;

		H_tb[nn[3*i]*H_size + nn[3*i+1]] += c_table[ nn[3*i+2] ]*phase;
		//printf("index: %i\n", nn[3*i]*H_size + nn[3*i+1]);
	}

}

void build_TaAsdHdk(lapack_int H_size, double* k, int* nn, double* nn_hops, int nn_total, 
	lapack_complex_double* dH_tb, lapack_complex_double* c_table){
	

	lapack_int i, index;
	lapack_complex_double coef;
	lapack_int H_shift = H_size * H_size; //Size of one Hamiltonian.
	
	//Set dHdk to all 0's
	memset(dH_tb, 0, H_size * H_size * 3 * sizeof(lapack_complex_double ));

	for( i = 0; i < nn_total; i++ ) {
		//For each nn pair
		//Compute phase and multiply by -I*coefficient
		coef = -I*c_table[ nn[3*i+2] ] * 
				cexp(- I * (k[0]*nn_hops[3*i] + 
						k[1]*nn_hops[3*i+1] + 
						k[2]*nn_hops[3*i+2]) );

		index = nn[3*i]*H_size + nn[3*i+1];
		
		//store x,y, and z derivatives
		dH_tb[index] += nn_hops[3*i]*coef;
		dH_tb[H_shift+index] += nn_hops[3*i+1]*coef;
		dH_tb[2*H_shift+index] += nn_hops[3*i+2]*coef;
	}

}

/* Generates Hamiltonian derivatives, dH/dk_x, dH/dk_y, dH/dk_z.
 * These are stored sequentially in the dH_tb array. 2*H_size refers
 * to the size of one dimension of dH/dk (2x for spin up/down).
 */
void build_dHdk(lapack_int H_size, double* k, double* nn, int nn_total, 
	lapack_complex_double* dH_tb, lapack_complex_double* coef_table ){
	
	lapack_int i, j;
	lapack_int x, y;
	lapack_complex_double phase = 0.0;
	lapack_int H_shift = H_size * H_size * 4; //Size of one Hamiltonian.
	//Initialize dH/dk to all 0's
	memset(dH_tb, 0, H_shift * 3 * sizeof(lapack_complex_double ));
	for( i = 0; i < nn_total; i++ ) {
		x = -1 + (int) nn[5*i+0];
		y = -1 + (int) nn[5*i+1];
		phase = k[0]*nn[5*i+2] + k[1]*nn[5*i+3] + k[2]*nn[5*i+4];
		phase = cexp(- I * phase);

		//Spin-up
		dH_tb[x*2*H_size + y] += -I*nn[5*i+2]*coef_table[2*i]*phase;
		dH_tb[H_shift+x*2*H_size + y] += -I*nn[5*i+3]*coef_table[2*i]*phase;
		dH_tb[2*H_shift+x*2*H_size + y] += -I*nn[5*i+4]*coef_table[2*i]*phase;

		//Spin-down
		dH_tb[(x+H_size)*2*H_size + H_size + y] += -I*nn[5*i+2]*coef_table[2*i+1]*phase;
		dH_tb[H_shift+(x+H_size)*2*H_size + H_size + y] += -I*nn[5*i+3]*coef_table[2*i+1]*phase;
		dH_tb[2*H_shift+(x+H_size)*2*H_size + H_size + y] += -I*nn[5*i+4]*coef_table[2*i+1]*phase;
	}


}

void build_testH(lapack_int H_size, double* k, double* nn, int nn_total, 
	lapack_complex_double* H_tb, lapack_complex_double* coef_table ){
	
	lapack_int i, j;
	lapack_int x, y;
	lapack_complex_double phase = 0.0;

	//***** coefs *****
	lapack_complex_double tx = -0.05;
	lapack_complex_double ty = 0.05;
	lapack_complex_double tz = 0.05;
	lapack_complex_double kw = 0.7854;
	lapack_complex_double mass = 0.1;

	//*** end coefs ***

	//Initialize dH/dk to all 0's
	lapack_int H_shift = H_size * H_size * 4; //Size of one Hamiltonian.
	memset(H_tb, 0, H_shift * sizeof(lapack_complex_double ));
	H_tb[0] += 2*tx*(cos(k[0]) - cos(kw));
	H_tb[3] += 2*tx*(cos(k[0]) - cos(kw));

	H_tb[1] += mass*(2 - cos(k[1]) - cos(k[2]));
	H_tb[2] += mass*(2 - cos(k[1]) - cos(k[2]));

	H_tb[1] += -I*(2*ty*sin(k[1]));
	H_tb[2] +=  I*(2*ty*sin(k[1]));

	H_tb[0] +=  (2*tz*sin(k[2])); //upper left
	H_tb[3] += -(2*tz*sin(k[2]));
	/*
	for( i = 1; i < 4; i++ ) {
	//for( i = 1; i < nn_total; i++ ) {
		//printf("\nx,y Coords: (%i, %i)\n", x , y);
		x = (int) nn[5*i+0];
		y = (int) nn[5*i+1];
		x--;
		y--;
		//Do by terms 
		
		//Const term
		H_tb[0] += 2*tx*(cos(k[0]*nn[5*i+2]) - cos(kw*nn[5*i+2])); //upper left
		H_tb[3] += 2*tx*(cos(k[0]*nn[5*i+2]) - cos(kw*nn[5*i+2])); //lower right
		// printf("Ham0: (%6.2f, %6.2f)\n", creal(H_tb[0]),cimag(H_tb[0]));
		// printf("Ham3: (%6.2f, %6.2f)\n", creal(H_tb[3]),cimag(H_tb[3]));
		//Sigma x
		H_tb[1] += mass*(2 - cos(k[1]*nn[5*i+3]) - cos(k[2]*nn[5*i+4])); //upper right
		H_tb[2] += mass*(2 - cos(k[1]*nn[5*i+3]) - cos(k[2]*nn[5*i+4])); //lower left
		// printf("Ham1: (%6.2f, %6.2f)\n", creal(H_tb[1]),cimag(H_tb[1]));
		// printf("Ham2: (%6.2f, %6.2f)\n", creal(H_tb[2]),cimag(H_tb[2]));
		//Sigma y
		H_tb[1] += -I*(2*ty*sin(k[1]*nn[5*i+3])); //upper right
		H_tb[2] +=  I*(2*ty*sin(k[1]*nn[5*i+3])); //lower left
		// printf("k: (%6.2f, %6.2f)\n", creal(k[1]),cimag(k[1]));
		// printf("nn: (%6.2f, %6.2f)\n", creal(nn[5*i+3]),cimag(nn[5*i+3]));
		// printf("sin(ky*ay): (%6.2f, %6.2f)\n", creal(sin(k[1]*nn[5*i+3])),cimag(sin(k[1]*nn[5*i+3])));
		// printf("Ham1: (%6.2f, %6.2f)\n", creal(H_tb[1]),cimag(H_tb[1]));
		// printf("Ham2: (%6.2f, %6.2f)\n", creal(H_tb[2]),cimag(H_tb[2]));
		// //Sigma z
		H_tb[0] +=  (2*tz*sin(k[2]*nn[5*i+4])); //upper left
		H_tb[3] += -(2*tz*sin(k[2]*nn[5*i+4])); //lower right
		// printf("Ham0: (%6.2f, %6.2f)\n", creal(H_tb[0]),cimag(H_tb[0]));
		// printf("Ham3: (%6.2f, %6.2f)\n", creal(H_tb[3]),cimag(H_tb[3]));

	}*/
}


void build_testdHdk(lapack_int H_size, double* k, double* nn, int nn_total, 
	lapack_complex_double* dH_tb, lapack_complex_double* coef_table ){
	
	lapack_int i, j;
	lapack_int x, y;
	lapack_complex_double phase = 0.0;

	//***** coefs *****
	lapack_complex_double tx = -0.05;
	lapack_complex_double ty = 0.05;
	lapack_complex_double tz = 0.05;
	lapack_complex_double kw = 0.7854;
	lapack_complex_double mass = 0.1;

	//*** end coefs ***

	//Initialize dH/dk to all 0's
	lapack_int H_shift = H_size * H_size * 4; //Size of one Hamiltonian.
	memset(dH_tb, 0, H_shift * 3 * sizeof(lapack_complex_double ));

	dH_tb[0] += 2*tx*(-sin(k[0]) ); //upper left
	dH_tb[3] += 2*tx*(-sin(k[0]) );	

	dH_tb[H_shift+1] += mass*(sin(k[1]) );
	dH_tb[H_shift+2] += mass*(sin(k[1]) );
	
	dH_tb[H_shift+1] += -I*(2*ty*cos(k[1]));
	dH_tb[H_shift+2] +=  I*(2*ty*cos(k[1]));
	
	dH_tb[2*H_shift+1] += mass*(sin(k[2]));
	dH_tb[2*H_shift+2] += mass*(sin(k[2]));
	
	dH_tb[2*H_shift+0] +=  (2*tz*cos(k[2]));
	dH_tb[2*H_shift+3] += -(2*tz*cos(k[2]));


	/*for( i = 1; i < 4; i++ ) {
	//for( i = 1; i < nn_total; i++ ) {
		//printf("\nx,y Coords: (%i, %i)\n", x , y);
		x = (int) nn[5*i+0];
		y = (int) nn[5*i+1];
		x--;
		y--;
		//Do by terms 
		
		//Dx Hamiltonian
		dH_tb[0] += 2*tx*(-nn[5*i+2]*sin(k[0]*nn[5*i+2]) ); //upper left
		dH_tb[3] += 2*tx*(-nn[5*i+2]*sin(k[0]*nn[5*i+2]) ); //lower right
		if (i == 0){
			printf("\n***dx terms***\n");
			printf("k[0] and a_x: %6.2f, %6.2f\n", k[0], nn[5*i+2]);
			printf("Ham val: (%6.2f, %6.2f)\n", creal(dH_tb[0]),cimag(dH_tb[0]));
		}
		//printf("Ham val: (%6.2f, %6.2f)\n", creal(dH_tb[0]),cimag(dH_tb[0]));
		//All other terms are 0

		//Dy Hamiltonian
		dH_tb[H_shift+0] += 0; //upper left
		dH_tb[H_shift+3] += 0; //lower right
		//Sigma x
		dH_tb[H_shift+1] += mass*(nn[5*i+3]*sin(k[1]*nn[5*i+3]) ); //upper right
		dH_tb[H_shift+2] += mass*(nn[5*i+3]*sin(k[1]*nn[5*i+3]) ); //lower left
		//Sigma y
		dH_tb[H_shift+1] += -I*(2*ty*nn[5*i+3]*cos(k[1]*nn[5*i+3])); //upper right
		dH_tb[H_shift+2] +=  I*(2*ty*nn[5*i+3]*cos(k[1]*nn[5*i+3])); //lower left
		//Sigma z
		dH_tb[H_shift+0] += 0; //upper left
		dH_tb[H_shift+3] += 0; //lower right

		//Dz Hamiltonian
		dH_tb[2*H_shift+0] += 0; //upper left
		dH_tb[2*H_shift+3] += 0; //lower right
		//Sigma x
		dH_tb[2*H_shift+1] += mass*(nn[5*i+4]*sin(k[2]*nn[5*i+4])); //upper right
		dH_tb[2*H_shift+2] += mass*(nn[5*i+4]*sin(k[2]*nn[5*i+4])); //lower left
		//Sigma y
		dH_tb[2*H_shift+1] += 0; //upper right
		dH_tb[2*H_shift+2] += 0; //lower left
		//Sigma z
		dH_tb[2*H_shift+0] +=  (2*tz*nn[5*i+4]*cos(k[2]*nn[5*i+4])); //upper left
		dH_tb[2*H_shift+3] += -(2*tz*nn[5*i+4]*cos(k[2]*nn[5*i+4])); //lower right

	}*/


}

void build_testH2(lapack_int H_size, double* k, double* nn, int nn_total, 
	lapack_complex_double* H_tb, lapack_complex_double* coef_table ){
	
	lapack_int i, j;
	lapack_int x, y;
	lapack_complex_double phase = 0.0;

	//***** coefs *****
	lapack_complex_double tx = -0.05;
	lapack_complex_double ty = 0.05;
	lapack_complex_double tz = 0.05;
	lapack_complex_double kw = 0.7854;
	lapack_complex_double mass = 0.1;
	//*** end coefs ***

	//Initialize dH/dk to all 0's
	lapack_int H_shift = H_size * H_size * 4; //Size of one Hamiltonian.
	memset(H_tb, 0, H_shift * sizeof(lapack_complex_double ));
	
	H_tb[1] += 2*tx*(cos(k[0]) - cos(kw)) + mass*(2 - cos(k[0]) - cos(k[1]));
	H_tb[2] += 2*tx*(cos(k[0]) - cos(kw)) + mass*(2 - cos(k[0]) - cos(k[1]));

	H_tb[1] += -I*(2*ty*sin(k[1]));
	H_tb[2] +=  I*(2*ty*sin(k[1]));

	H_tb[0] +=  (2*tz*sin(k[2]));
	H_tb[3] += -(2*tz*sin(k[2]));
}

void build_test2dHdk(lapack_int H_size, double* k, double* nn, int nn_total, 
	lapack_complex_double* dH_tb, lapack_complex_double* coef_table ){
	
	lapack_int i, j;
	lapack_int x, y;
	lapack_complex_double phase = 0.0;

	//***** coefs *****
	lapack_complex_double tx = -0.05;
	lapack_complex_double ty = 0.05;
	lapack_complex_double tz = 0.05;
	lapack_complex_double kw = 0.7854;
	lapack_complex_double mass = 0.1;

	//*** end coefs ***

	//Initialize dH/dk to all 0's
	lapack_int H_shift = H_size * H_size * 4; //Size of one Hamiltonian.
	memset(dH_tb, 0, H_shift * 3 * sizeof(lapack_complex_double ));
	//X
	dH_tb[1] += -2*tx*sin(k[0]) + mass*sin(k[0]) ;
	dH_tb[2] += -2*tx*sin(k[0]) + mass*sin(k[0]) ;
	//Y
	dH_tb[H_shift+1] += mass*sin(k[1]);
	dH_tb[H_shift+2] += mass*sin(k[1]);
	dH_tb[H_shift+1] += -I*(2*ty*cos(k[1]));
	dH_tb[H_shift+2] +=  I*(2*ty*cos(k[1]));
	//Z
	dH_tb[2*H_shift+0] +=  (2*tz*cos(k[2]));
	dH_tb[2*H_shift+3] += -(2*tz*cos(k[2]));
}

//void compute_BerryPhase()

int compare_nn( const void* vec1, const void* vec2){
	double* v1 = (double*)vec1;
	double* v2 = (double*)vec2;
	double d1, d2, xy_a1, xy_a2;
	double proj1, proj2;
	d1 = pow(v1[2],2) + pow(v1[3],2) + pow(v1[4],2);
	d2 = pow(v2[2],2) + pow(v2[3],2) + pow(v2[4],2);

	//Project onto xy-plane, and dot with (1,0)
	proj1 = v1[2]/sqrt( pow(v1[2],2) + pow(v1[3],2) );
	proj2 = v2[2]/sqrt( pow(v2[2],2) + pow(v2[3],2) );
	
	//compute angle (in rad) counterclockwise from (1,0)
	xy_a1 = acos(proj1);
	xy_a2 = acos(proj2);
	if(v1[3] < 0.0){xy_a1 = 2.0*M_PI - xy_a1;}
	if(v2[3] < 0.0){xy_a2 = 2.0*M_PI - xy_a2;}

	//Sort by hopping vector length first
	if (fabs(d1-d2)<0.01) {
		//Sort by starting atom 
		if(fabs(v1[0]-v2[0]) < 0.01){
			//Sort by end atom
			if(fabs(v1[1]-v2[1]) < 0.01){
				//Sort by xy-plane angle
				if(fabs(xy_a1-xy_a2) < 0.01){
					//Sort by z-value
					if(fabs(v1[5]-v2[5]) < 0.01){
						return 0;
					}
					else if (v1[5]<v2[5]){
						return -1;
					}
					else if (v1[5]>v2[5]){
						return 1;
					}
				}
				else if (xy_a1<xy_a2){
					return -1;
				}
				else if (xy_a1>xy_a2){
					return 1;
				}
			}
			else if (v1[1]<v2[1]){
				return -1;
			}
			else if (v1[1]>v2[1]){
				return 1;
			}
		}
		else if (v1[0]<v2[0]){
			return -1;
		}
		else if (v1[0]>v2[0]){
			return 1;
		}
		
	}
	else if (d1 < d2){
		return -1;
	}
	else if (d1 > d2){
		return 1;
	}

	return 0;
}

/* This finds the two NN hoppings that make up a given NNN hop.
 * It then takes their cross product and returns the resulting 
 * vector normal to the NNN hop.
 * NOTE THAT THIS ONLY WORKS FOR
 */
void v_perp(double* nn_array, int* nn_counts, int nn_index, double* vec){
	vec[0] = 0.0, vec[1] = 0.0, vec[2] = 0.0;
	double* v1;
	double* v2;
	double mag = 1.0;
	lapack_int i, j;
	if(nn_index < (nn_counts[0]+ nn_counts[1])){
		//Suppress output, since we return (0,0,0) anyways
		//printf("INVALID HOPPING VECTOR! (Not a NNN hop)\n");
	}
	//printf("SumShouldBe = (%6.2f, %6.2f, %6.2f)\n", nn_array[5*nn_index+2], nn_array[5*nn_index+3], nn_array[5*nn_index+4]);
	for(i = nn_counts[0]; i < nn_counts[0]+nn_counts[1]; i++){
		//printf("Test");
		v1 = &nn_array[5*i];
		for(j = nn_counts[0]; j < nn_counts[0]+nn_counts[1]; j++){
			v2 = &nn_array[5*j];
			//printf("Vector1 = (%6.2f, %6.2f, %6.2f)   ", v1[0], v1[1], v1[2]);
			//printf("Vector2 = (%6.2f, %6.2f, %6.2f)\n", v2[0], v2[1], v2[2]);
			if(	(((int) nn_array[5*nn_index]) == ((int) v1[0]) ) &&
				(((int) nn_array[5*nn_index+1]) == ((int) v2[1])) &&
				(fabs(nn_array[5*nn_index+2]-(v1[2]+v2[2]))<0.01) &&
				(fabs(nn_array[5*nn_index+3]-(v1[3]+v2[3]))<0.01) &&
				(fabs(nn_array[5*nn_index+4]-(v1[4]+v2[4]))<0.01)){
				/*printf("Nums:: (%i > %i) =?(%i > %i) + (%i > %i)\n", (int) nn_array[5*nn_index], 
					(int) nn_array[5*nn_index+1], (int) v1[0], (int) v1[1], (int) v2[0], (int) v2[1]);
				printf("Vector1 = (%6.2f, %6.2f, %6.2f)   ", v1[2], v1[3], v1[4]);
				printf("Vector2 = (%6.2f, %6.2f, %6.2f)\n", v2[2], v2[3], v2[4]);
				printf("VectorTotal = (%6.2f, %6.2f, %6.2f)\n", nn_array[5*nn_index+2], nn_array[5*nn_index+3], nn_array[5*nn_index+4]);
				*/
				//If match, compute cross product
				vec[0] = (v1[3]*v2[4] - v1[4]*v2[3]);
				vec[1] = (v1[4]*v2[2] - v1[2]*v2[4]);
				vec[2] = (v1[2]*v2[3] - v1[3]*v2[2]);
				//Divide by 1st magnitude
				mag = sqrt( pow(v1[2],2)+pow(v1[3],2)+pow(v1[4],2) );
				vec[0] /= mag;
				vec[1] /= mag;
				vec[2] /= mag;
				//Divide by 2nd magnitude
				mag = sqrt( pow(v2[2],2)+pow(v2[3],2)+pow(v2[4],2) );
				vec[0] /= mag;
				vec[1] /= mag;
				vec[2] /= mag;
				//Temp renormalization
				//vec[0] /= 0.866;
				//vec[1] /= 0.866;
				//vec[2] /= 0.866;
				//break out of loop (This seems sloppy :/)
				j = nn_counts[0]+nn_counts[1];
				i = nn_counts[0]+nn_counts[1];
			}
		}
	}
	//printf("Perp Vector: (%6.4f,%6.4f,%6.4f)\n", vec[0], vec[1], vec[2]);
	//return &vec[0];
}
