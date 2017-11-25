#include <lapacke.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mat_print.h"

/* Auxiliary routine: printing a matrix */
void print_matrix_rowmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );

        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", mat[i*ldm+j] );
                printf( "\n" );
        }
}


/* Auxiliary routine: printing a matrix */
void print_matrix_colmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );

        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", mat[i+j*ldm] );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a vector of integers */
void print_vector( char* desc, lapack_int n, lapack_int* vec ) {
        lapack_int j;
        printf( "\n %s\n", desc );
        for( j = 0; j < n; j++ ) printf( " %6i", vec[j] );
        printf( "\n" );
}

/* Auxiliary routine: printing a vector of doubles */
void print_dvector( char* desc, lapack_int n, double* vec ) {
	lapack_int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %6.2f", vec[j] );
	printf( "\n" );
}

/* Auxiliary routine: printing a matrix*/
void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda ) {
    lapack_int i, j;
    printf( "\n %s\n", desc );
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ )
		printf( " (%6.2f,%6.2f)", creal(a[i*lda+j]), cimag(a[i*lda+j]) );
    printf( "\n" );
    }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda ) {
	lapack_int i, j;
    printf( "\n %s\n", desc );
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
        printf( "\n" );
    }
}

void read_coords( char* filename, double* cell_vecs, 
	double* lat_vecs , double** atom_vecs, lapack_int* numatoms, lapack_int* nndepth){
	
	
	lapack_int i ,j , n;
	//(*numatoms) = 1;

	//printf("TEST: %i\n\n",  (*numatoms));
	//lapack_int numatoms;
	
	FILE *fp;
	fp = fopen(filename, "r");
	char word[1024];
	while( fscanf(fp, "%1023s", word) != EOF )
	{
		//printf("line-> %s\n", word);
		//Get number of cells in supercell
		if (strcmp(word, "[CELL]") == 0) {
			printf("Cell!\n");
			for( i = 0; i < 3; i++ ) {
				fscanf(fp, "%1023s", word);
				cell_vecs[i] = atof(word);
				//printf("cell[%i] %s\n", i, word);
			}
		}
		if (strcmp(word, "Natoms") == 0) {
			fscanf(fp, "%1023s", word);
			(*numatoms) = atoi(word);
			//printf("number of atoms: %i\n", *numatoms);
			n = (*numatoms);
		}
		if (strcmp(word, "NN_depth") == 0) {
			fscanf(fp, "%1023s", word);
			(*nndepth) = atoi(word);
		}
		if (strcmp(word, "[VECTORS]") == 0) {
			//printf("Cell!\n");
			for( i = 0; i < 3; i++ ) {
				for ( j = 0; j < 3; j++ ){
					fscanf(fp, "%1023s", word);
					lat_vecs[i*3+j] = atof(word);
				}
				
				//printf("cell[%i] %s\n", i, word);
			}
		}
		if (strcmp(word, "[ATOMS]") == 0) {
			//printf("Cell!\n");
			//(double *)
			*atom_vecs =  malloc( n * 3 * sizeof(double ));
			if (*atom_vecs == NULL) {
				printf("Memory allocation error when allocating atom vector array.\n");
				exit(0);
			}
			memset((*atom_vecs), 0, n * 3 * sizeof(double ));
			//printf("Test0: %8.5f\n", (*atom_vecs)[0]);
			//(*atom_vecs)[1] = 0.1;
			//printf("Test1: %8.5f\n", (*atom_vecs)[1]);
			//printf("Passed!, %i, %i\n", n , i);
			for( i = 0; i < n; i++ ) {
				fscanf(fp, "%1023s", word);
				//printf("test: %1023s", word);
				for ( j = 0; j < 3; j++ ){
					fscanf(fp, "%1023s", word);
					(*atom_vecs)[i*3+j] = atof(word);
					//printf("data[%i,%i] = %8.5f\n", i, j , (*atom_vecs)[i*3+j]);
				}
				//printf("cell[%i] %s\n", i, word);
			}

			//printf("test: vec[%i,%i] = %8.5f\n", 1, 1 , atom_vecs[1*3+1]);
		}
	}

	fclose(fp);
}

void read_H( char* filename){//, double* nn_counts, double** H_tb){

	lapack_int i ,j , n;

	printf("\n\n\n");
	//Need local arrays to store Hamiltonian params
	lapack_int nparams, h_length, nndepth;
	lapack_complex_double* param_vals;
	double temp1,temp2;
	
	lapack_int NUM_OF_LETTERS = 10;

	FILE *fp;
	fp = fopen(filename, "r");

	char** param_names;
	//char** ham;
	char (*ham)[10];
	char word[1024];

	/*
		Start by allocating and reading in 
		- the Hamiltonian
		- the parameters
		- the values of the parameters
	*/
	while( fscanf(fp, "%1023s", word) != EOF )
	{
		if (strcmp(word, "N_params") == 0) {
			fscanf(fp, "%1023s", word);
			nparams = atoi(word);
			//Allocate char array for parameters
			param_names = (char**)malloc(nparams * sizeof(char*));
		}
		if (strcmp(word, "H_length") == 0) {
			fscanf(fp, "%1023s", word);
			h_length = atoi(word);
			//ham = (char**)malloc(h_length * sizeof(char*));
			ham = malloc(h_length * sizeof(char[5]));
		}
		if (strcmp(word, "NN_depth") == 0) {
			fscanf(fp, "%1023s", word);
			nndepth = atoi(word);
		}
		//Read in parameter names and values for Hamiltonian
		if (strcmp(word, "[HPARAMS]") == 0) {
			
			param_vals = (lapack_complex_double * ) malloc( nparams * sizeof(lapack_complex_double) );
			if (param_vals == NULL) {
				printf("Memory allocation error when allocating Hamiltonian parameter array.\n");
				exit(0);
			}
			memset(param_vals, 0, nparams * sizeof(lapack_complex_double));

			for( i = 0; i < nparams; i++ ) {
				fscanf(fp, "%1023s", word);

				param_names[i] = (char*)malloc(sizeof(word));
				if (param_names[i] == NULL) {
					printf("Memory allocation error when allocating parameter name string for \"%s\".\n", word);
					exit(0);
				}

				strcpy(param_names[i], word);
				
				//Scan in real and imag part of the parameter.
				fscanf(fp, "%1023s", word);
				temp1 = atof(word); //read
				fscanf(fp, "%1023s", word);
				temp2 = atof(word); //imag

				param_vals[i] = temp1 + I*temp2;
				//printf( "Paramtest: (%6.2f,%6.2f)\n", creal( param_vals[i] ), cimag( param_vals[i] ) );
			}
		}
		if (strcmp(word, "[HAMILTONIAN]") == 0) {
			
			for( i = 0; i < h_length; i++ ) {
				fscanf(fp, "%1023s", word);
				
				/*ham[i] = (char*)malloc(sizeof(word));
				if (ham[i] == NULL) {
					printf("Memory allocation error when allocating string for \"%s\" in Hamiltonian.\n", word);
					exit(0);
				}*/
				strcpy(ham[i], word);	
			}
		}
	}
	fclose(fp);
	//printf("Hamiltonian: \n");
	//printf("%s\n", ham[0]);
	//printf("%s\n", ham[1]);
	for( j = 0; j < h_length; j++ ) {
		//printf("%s\n", ham[j]);
	}

	//printf("Test: %c", param_names[2][1]);
	//test = ham;

	//printf("Test: %s\n", test[0]);

	/*
		Constructing the Hamiltonian.
	*/

	//do onsite, then nn, then nnn ...
	
	for( i = 0; i < nndepth; i++){
		char stopstring[10]; 
		strcpy(stopstring, "Sigma_n");
		for (j = 0; j < i+1; j++)
		{
			strcat(stopstring, "n");
		}
		//printf("Stopstring: %s\n", stopstring);

		//for( j = 0; j < nn_counts[i]; i++){

		//}

	}

}

//void build_H()


void find_neighbors( lapack_int n, double* lat_vecs, double* atom_vecs,
	double nndepth, double ** nn, int* nn_counts ){

	lapack_int i ,j, k, l, m;

	double a[3] = {0.0, 0.0, 0.0}, b[3];
	double len;

	lapack_int rownum = 0;

	// TODO!! Make this into a passable parameter.
	//lapack_int nndepth = 2;

	//Temporary neighbor array to store all possible NN.
	double * tempnn;

	//Want to search all adjacent unit cells nndepth deep
	lapack_int dim = (2*nndepth + 1); 
	lapack_int vol = dim*dim*dim;

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
						rownum = m*n*vol + l*vol + (i+2)*dim*dim+(j+2)*dim+(k+2);
						
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

	double nnvecL = 0.0;
	double nnvecLtemp = 0.0;
	//int nntotal = 0;
	j = -1;
	i = 0;

	while( i < nndepth + 1 ){
		j++;
		nnvecLtemp = pow(tempnn[5*j+2],2) + pow(tempnn[5*j+3],2) + pow(tempnn[5*j+4],2);
		if ((nnvecLtemp - nnvecL) > 0.01 ){
			nnvecL = nnvecLtemp;
			i++;
		}
		nn_counts[i] += 1;
		//nntotal ++;
	}
	//printf("NNcounts: (%i, %i, %i )\n", nn_counts[0], nn_counts[1], nn_counts[2]);
	//printf("NNtotal: %i , %i\n", j);

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
}


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

double length_v( double* v1, double* v2){
	//Assumes that they are length 3!!!
	double t1, t2, t3;

	t1 = pow((v1[0] - v2[0]),2);
	t2 = pow((v1[1] - v2[1]),2);
	t3 = pow((v1[2] - v2[2]),2);

	return t1+t2+t3;
}