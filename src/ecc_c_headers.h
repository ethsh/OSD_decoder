#ifndef __ECC_C_HEADERS_
#define __ECC_C_HEADERS_

#include "ecc.h"

/// function that are implemented in ecc.c but not declared in ecc.h
int_mtx_t ** create_sys_BCH(const char * s, int n);
void encode_symbols(gf2_vec_t *info, int_mtx_t* G, gf2_vec_t * symb);
void init_measured_data(gf2_vec_t *symb, real_vec_t *data, double EBN0_db, double R);
void apply_perm_int_mtx(int_mtx_t* mtx, int *ord);
void apply_perm_real(real_vec_t *lambda, int *ord);
void apply_inv_perm_int(gf2_vec_t *dec_symb, int_vec_t* PHI);
void dev_poly(int_mtx_t *dividend, int_mtx_t *genpoly, int_mtx_t *quotient, int_mtx_t *reminder);
void init_symbol_data(gf2_vec_t *symb);


#endif