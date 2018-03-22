#ifndef HAM_UTILS
#define HAM_UTILS


//lapack_complex_double graphene_H(lapack_int n, double* spin, double* vec, lapack_complex_double* coef);


void find_neighbors( lapack_int numatoms, double* lat_vecs , double* atom_vecs 
    ,lapack_int nndepth, double ** nn, int* nn_counts ); 

void find_neighbors_cutoff( lapack_int n, double* lat_vecs , double* atom_vecs, int* H_map, lapack_int* map_S  
    ,double cutoff, lapack_int ** nn, double ** nn_vecs, lapack_int* nn_total );

int compare_nn( const void* vec1, const void* vec2);

void build_C(lapack_int nn_depth, int nn_total, double* nn, 
    int* nn_counts, lapack_complex_double* coef_table,
    lapack_complex_double (*ham)(lapack_int, double*, double*, lapack_complex_double*));

void build_H(lapack_int H_size, double* k, double* nn, int nn_total, 
    lapack_complex_double* H_tb, lapack_complex_double* coef_table );

void build_testH(lapack_int H_size, double* k, double* nn, int nn_total, 
    lapack_complex_double* H_tb, lapack_complex_double* coef_table );

void build_testdHdk(lapack_int H_size, double* k, double* nn, int nn_total, 
    lapack_complex_double* H_tb, lapack_complex_double* coef_table );

void build_test2dHdk(lapack_int H_size, double* k, double* nn, int nn_total, 
    lapack_complex_double* dH_tb, lapack_complex_double* coef_table );

void build_testH2(lapack_int H_size, double* k, double* nn, int nn_total, 
    lapack_complex_double* H_tb, lapack_complex_double* coef_table );

void build_dHdk(lapack_int H_size, double* k, double* nn, int nn_total, 
    lapack_complex_double* dH_tb, lapack_complex_double* coef_table );
        
    
void v_perp(double* nn_array, int* nn_counts, int nn_index, double* vec);


void build_TaAsH(lapack_int H_size, double* k, int* nn, double* nn_hops, int nn_total, 
    lapack_complex_double* H_tb, lapack_complex_double* c_table );

void build_TaAsdHdk(lapack_int H_size, double* k, int* nn, double* nn_hops, int nn_total, 
    lapack_complex_double* dH_tb, lapack_complex_double* c_table );

void find_weyl( lapack_int nkp, double* bC, double*  kdif);
void find_weyl2( lapack_int nkp, double* bC, double*  kdif);
void find_weyl3( lapack_int nkp, double* bC, double*  kdif, double* bZvec, lapack_int bN, lapack_int btotal);

void find_weyl4( lapack_int nkp, lapack_complex_double * evc, lapack_int n, lapack_int bN, double* bC);

#endif /* HAM_UTILS */