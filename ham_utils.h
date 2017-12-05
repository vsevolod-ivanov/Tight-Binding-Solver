#ifndef HAM_UTILS
#define HAM_UTILS


//lapack_complex_double graphene_H(lapack_int n, double* spin, double* vec, lapack_complex_double* coef);


void find_neighbors( lapack_int numatoms, double* lat_vecs , double* atom_vecs 
    ,lapack_int nndepth, double ** nn, int* nn_counts ); 

int compare_nn( const void* vec1, const void* vec2);

void build_C(lapack_int nn_depth, int nn_total, double* nn, 
    int* nn_counts, lapack_complex_double* coef_table,
    lapack_complex_double (*ham)(lapack_int, double*, double*, lapack_complex_double*));

void build_H(lapack_int H_size, double* k, double* nn, int nn_total, 
    lapack_complex_double* H_tb, lapack_complex_double* coef_table );

void build_dHdk(lapack_int H_size, double* k, double* nn, int nn_total, 
    lapack_complex_double* dH_tb, lapack_complex_double* coef_table );
        
    
void v_perp(double* nn_array, int* nn_counts, int nn_index, double* vec);

#endif /* HAM_UTILS */