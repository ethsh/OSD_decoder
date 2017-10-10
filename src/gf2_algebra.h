#ifndef __GF2_ALGEBRA_
#define __GF2_ALGEBRA_

#include <memory.h>
#include "math_funcs.h"

#define WORD_SIZE (8)

static unsigned char masks[] = { 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1 };

typedef struct gf2_bit_vec_s {
	int n;
	int n_char;
	char *data;

	int n_int;
	int n_remain_start;
} gf2_bit_vec_t;

typedef struct gf2_bit_mtx_s {
	int
	n,		// number of cols
	m;		// number of rows
	//char *data;	// the matrix itself
	gf2_bit_vec_t **data;	// the matrix itself - every gf2_bit_vec_t is a row in the mtx
} gf2_bit_mtx_t;

// because aussian elimination will be faster this way
typedef struct gf2_byte_mtx_s {
	int
	n,		// number of cols
	m;		// number of rows
	char *data;	// the matrix itself
} gf2_byte_mtx_t;


void malloc_and_init_gf2_bit_vec(gf2_bit_vec_t** dst, int len);
void init_gf2_bit_vec(gf2_bit_vec_t* vec, int len);
void copy_gf2_bit_vec(gf2_bit_vec_t* src, gf2_bit_vec_t* dst);
void copy_gf2_bit_vec_range(gf2_bit_vec_t* src, gf2_bit_vec_t* dst, int start_index, int end_index);
void free_gf2_bit_vec(gf2_bit_vec_t* vec);
void add_gf2_bit_vec(gf2_bit_vec_t* src1, gf2_bit_vec_t* src2, gf2_bit_vec_t* dst);
void add_gf2_bit_vec_in_place(gf2_bit_vec_t* src1_dst, gf2_bit_vec_t* src2);
__inline char gf2_bit_vec_get_element(gf2_bit_vec_t* vec, int index);
__inline void gf2_bit_vec_set_element(gf2_bit_vec_t* vec, int index, char val);
void apply_inv_perm_gf2_bit_vec(gf2_bit_vec_t *dec_symb, int_vec_t* PHI);
void gf2_bit_vec_2_gf2_vec(gf2_bit_vec_t* src, gf2_vec_t* dst);
void print_gf2_bit_vec(gf2_bit_vec_t* vec);




void malloc_and_init_gf2_bit_mtx(gf2_bit_mtx_t **mtx, int m, int n);
void init_gf2_bit_mtx(gf2_bit_mtx_t *mtx, int m, int n);
void free_gf2_bit_mtx(gf2_bit_mtx_t *x);
__inline char gf2_bit_mtx_get_element(gf2_bit_mtx_t *mtx, int line_idx, int col_idx);
__inline void gf2_bit_mtx_set_element(gf2_bit_mtx_t *mtx, int line_idx, int col_idx, char in_data);
void copy_gf2_bit_mtx(gf2_bit_mtx_t* src, gf2_bit_mtx_t* dst);
void int_mtx_2_gf2_bit_mtx(int_mtx_t* src, gf2_bit_mtx_t** dst);
void switch_gf2_bit_mtx_rows(gf2_bit_mtx_t* mtx, int row_a, int row_b);
void encode_gf2_bit_symbols(gf2_bit_vec_t *info, gf2_bit_mtx_t *G, gf2_bit_vec_t *symb);
void gf2_bit_mtx_gaussian_elimination(gf2_bit_mtx_t* G, int_vec_t* perm1, col_perms_t* col_perms);
void apply_inv_perm_gf2_bit_mtx(gf2_bit_mtx_t* mtx, int *ord);
void print_gf2_bit_mtx(gf2_bit_mtx_t* mtx);



void malloc_and_init_gf2_byte_mtx(gf2_byte_mtx_t **mtx, int m, int n);
void init_gf2_byte_mtx(gf2_byte_mtx_t *mtx, int m, int n);
void free_gf2_byte_mtx(gf2_byte_mtx_t *x);
__inline char gf2_byte_mtx_get_element(gf2_byte_mtx_t *mtx, int line_idx, int col_idx);
__inline void gf2_byte_mtx_set_element(gf2_byte_mtx_t *mtx, int line_idx, int col_idx, char in_data);
void copy_gf2_byte_mtx(gf2_byte_mtx_t* src, gf2_byte_mtx_t* dst);
void int_mtx_2_gf2_byte_mtx(int_mtx_t* src, gf2_byte_mtx_t** dst);
void switch_gf2_byte_mtx_rows(gf2_byte_mtx_t* mtx, int row_a, int row_b);
void gf2_byte_mtx_gaussian_elimination(gf2_byte_mtx_t* G, int_vec_t* perm1, col_perms_t* col_perms);
void apply_inv_perm_gf2_byte_mtx(gf2_byte_mtx_t* mtx, int *ord);
void print_gf2_byte_mtx(gf2_byte_mtx_t* mtx);

void convert_gf2_byte_mtx_2_gf2_bit_mtx(gf2_byte_mtx_t* src, gf2_bit_mtx_t* dst);



__inline char gf2_bit_vec_get_element(gf2_bit_vec_t* vec, int index) {
	int char_index = index / WORD_SIZE;
	int bit_index = WORD_SIZE - 1 - (index % WORD_SIZE);
	//return (vec->data[char_index] >> (bit_index)) & 1;
	return (vec->data[char_index] & masks[bit_index]) != 0 ? 1 : 0;
}

__inline void gf2_bit_vec_set_element(gf2_bit_vec_t* vec, int index, char val) {
	int char_index = index / WORD_SIZE;
	int bit_index = WORD_SIZE - 1 - (index % WORD_SIZE);
	if (val == 0) {
		vec->data[char_index] &= ~masks[bit_index];
	}
	else {
		vec->data[char_index] |= masks[bit_index];
	}
}

__inline char gf2_bit_mtx_get_element(gf2_bit_mtx_t *mtx, int line_idx, int col_idx) {
	ecc_assert((line_idx >= 0) && (line_idx<mtx->m), "Line exceeds allowed (%d vs %d)\n", mtx->m, line_idx);
	ecc_assert((col_idx >= 0) && (col_idx <mtx->n), "Col exceeds allowed (%d vs %d)\n", mtx->n, col_idx);
	return gf2_bit_vec_get_element(mtx->data[line_idx], col_idx);
	//return mtx->data[line_idx * mtx->n + col_idx];
}
__inline void gf2_bit_mtx_set_element(gf2_bit_mtx_t *mtx, int line_idx, int col_idx, char in_data) {
	ecc_assert((line_idx >= 0) && (line_idx <= mtx->m), "Line exceeds allowed (%d vs %d)\n", mtx->m, line_idx);
	ecc_assert((col_idx >= 0) && (col_idx <= mtx->n), "Col exceeds allowed (%d vs %d)\n", mtx->n, col_idx);
	//mtx->data[line_idx * mtx->n + col_idx] = in_data;
	gf2_bit_vec_set_element(mtx->data[line_idx], col_idx, in_data);
}



__inline char gf2_byte_mtx_get_element(gf2_byte_mtx_t *mtx, int line_idx, int col_idx) {
	ecc_assert((line_idx >= 0) && (line_idx <= mtx->m), "Line exceeds allowed (%d vs %d)\n", mtx->m, line_idx);
	ecc_assert((col_idx >= 0) && (col_idx <= mtx->n), "Col exceeds allowed (%d vs %d)\n", mtx->n, col_idx);
	return mtx->data[line_idx * mtx->n + col_idx];
}

__inline void gf2_byte_mtx_set_element(gf2_byte_mtx_t *mtx, int line_idx, int col_idx, char in_data) {
	ecc_assert((line_idx >= 0) && (line_idx <= mtx->m), "Line exceeds allowed (%d vs %d)\n", mtx->m, line_idx);
	ecc_assert((col_idx >= 0) && (col_idx <= mtx->n), "Col exceeds allowed (%d vs %d)\n", mtx->n, col_idx);
	mtx->data[line_idx * mtx->n + col_idx] = in_data;
}


#endif