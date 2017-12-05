#include <lapacke.h>
#include <cblas.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mat_print.h"
//#include "hamiltonian.h"

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
		printf( " (%10.8f,%10.8f)", creal(a[i*lda+j]), cimag(a[i*lda+j]) );
    printf( "\n" );
    }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda ) {
	lapack_int i, j;
    printf( "\n %s\n", desc );
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ ) printf( " %6.4f", a[i*lda+j] );
        printf( "\n" );
    }
}

void read_coords( char* filename, double* cell_vecs, 
	double* lat_vecs , double** atom_vecs, double** k_path, lapack_int* numatoms, lapack_int* nndepth){
	
	
	lapack_int i ,j , n;
	lapack_int nkp;
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
			//printf("Cell!\n");
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
		if (strcmp(word, "[KPATH]") == 0) {
			fscanf(fp, "%1023s", word);
			nkp = atoi(word);
			*k_path =  malloc( nkp * 3 * sizeof(double ));

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

/* Useless :(
*/ 
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


double length_v( double* v1, double* v2){
	//Assumes that they are length 3!!!
	double t1, t2, t3;

	t1 = pow((v1[0] - v2[0]),2);
	t2 = pow((v1[1] - v2[1]),2);
	t3 = pow((v1[2] - v2[2]),2);

	return t1+t2+t3;
}