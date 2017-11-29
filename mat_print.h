#ifndef _LAPACKE_EXAMPLE_AUX_
#define _LAPACKE_EXAMPLE_AUX_


void print_matrix_rowmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm );
void print_matrix_colmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm );
void print_vector( char* desc, lapack_int n, lapack_int* vec );
void print_dvector( char* desc, lapack_int n, double* vec );
void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda );
void print_rmatrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda );
void read_coords( char* filename, double* cell_vecs, double* lat_vecs , 
    double** atom_vecs, lapack_int * numatoms, lapack_int * nndepth );

//Computation

void find_neighbors( lapack_int numatoms, double* lat_vecs , double* atom_vecs 
    ,lapack_int nndepth, double ** nn, int* nn_counts ); 
void read_H( char* filename);
int compare_nn( const void* vec1, const void* vec2);
double length_v( double* v1, double* v2);

void build_C(lapack_int nn_depth, int nn_total, double* nn, 
    int* nn_counts, lapack_complex_double* coef_table,
    lapack_complex_double (*ham)(lapack_int, double*, double*, lapack_complex_double*));
void build_H(lapack_int H_size, double* k, double* nn, int nn_total, 
    lapack_complex_double* H_tb, lapack_complex_double* coef_table );
void v_perp(double* nn_array, int* nn_counts, int nn_index, double* vec);

#endif /* _LAPACKE_EXAMPLE_AUX_*/