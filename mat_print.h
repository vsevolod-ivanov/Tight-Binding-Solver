#ifndef _IO_FUNCTIONS_
#define _IO_FUNCTIONS_


void print_matrix_rowmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm );
void print_matrix_colmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm );
void print_vector( char* desc, lapack_int n, lapack_int* vec );
void print_dvector( char* desc, lapack_int n, double* vec );
void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda );
void print_rmatrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda );
void read_coords( char* filename, double* cell_vecs, double* lat_vecs , 
    double** atom_vecs, double** k_path, lapack_int * numatoms, lapack_int * nndepth );

void read_Map( char* filename, int* n_map, int* map_S, int** H_map);

void read_Coef( char* filereal, char* fileimag, lapack_complex_double** coeftable,
    lapack_int n_states, lapack_int n_coefs);

void printBerry(char* filename, double* bV, int numV);
void printBerryBand(char* filename, double* bV, lapack_int bN,lapack_int btotal, lapack_int numV);
void readBerry(char* filename, double* bV);
//Computation

void read_H( char* filename);



double length_v( double* v1, double* v2);

#endif /* _IO_FUNCTIONS_ */