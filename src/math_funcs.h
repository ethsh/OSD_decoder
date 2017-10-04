/// This file contains peripheral functions  the the OSD decoder defined in OSDDecoder.c
/// some functions are declared here but implemented in ecc.c. 
/// some of them are also implemented in math_funcs.c, but sould be implemented in ecc.c :)

#ifndef __MATH_FUNCS_
#define __MATH_FUNCS_

#include "ecc_c_headers.h"

typedef struct col_perms_s{
	int num_of_perms;
	int* a_indices;
	int* b_indices;
} col_perms_t;

void free_col_perms(col_perms_t *col_perms);
void init_col_perms_t(col_perms_t* col_perms, int len);

/// taken from osd_decoding_biawgn.c, with changes for our structs.
void OSDLeftSysMatrix_ethan(int_mtx_t* G, int_vec_t* perm1, col_perms_t* col_perms);


/// new funcs not implemented in ecc.h/c ///

// like init, but also malloc the gf2_vec_t*
void malloc_and_init_gf2_vec(gf2_vec_t** dst, int len);
void malloc_and_init_real_vec(real_vec_t** dst, int len);
void malloc_and_init_int_vec(int_vec_t** dst, int len);

// initialize a real_vec_t from a given real_vec_t (mallocing the data in dst)
void init_real_vec_from_rvec(real_vec_t* dst, real_vec_t* src);
// initialize an int_vec_t from a given int_vec_t  (mallocing the data in dst)
void init_int_vec_from_ivec(int_vec_t* dst, int_vec_t* src);
// initialize an gf2_vec_t from a given gf2_vec_t  (mallocing the data in dst)
void init_gf2_vec_from_gf2vec(gf2_vec_t* dst, gf2_vec_t* src);

// copied the content a src->data to dst->data, without reallocation dst->data
void copy_int_vec(int_vec_t* dst, int_vec_t* src);

// permutations applying
void apply_inv_perm_real(real_vec_t *lambda, int *ord);
void apply_inv_perm_int_mtx(int_mtx_t* mtx, int *ord);

/// switch two rows in a matrix
void switch_int_mtx_rows(int_mtx_t* mtx, int row_a, int row_b);

// vector arithmetics
void gf2_vec_add(gf2_vec_t* src1, gf2_vec_t* src2, gf2_vec_t* dst);



/// new functions, pretty similar to create_sys_BCH, but with a systematic G to the left
int_mtx_t ** create_sys_BCH_ethan(const char * s, int n);
/// taken from ecc.c, trying to be more efficient
void encode_symbols_ethan(gf2_vec_t *info, int_mtx_t *G, gf2_vec_t *symb);



#endif
