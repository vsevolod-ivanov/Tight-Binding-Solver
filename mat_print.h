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

//Computation

void read_H( char* filename);

double length_v( double* v1, double* v2);

#endif /* _IO_FUNCTIONS_ */