#include "gf2_algebra.h"

/// in bits arrays - the first bit is the MSB

/************* gf2_bit_vec_t functions *****************/

void malloc_and_init_gf2_bit_vec(gf2_bit_vec_t** dst, int len) {
	*dst = (gf2_bit_vec_t*)ecc_malloc(sizeof(gf2_bit_vec_t));
	if (len > 0) {
		init_gf2_bit_vec(*dst, len);
	}

}
void init_gf2_bit_vec(gf2_bit_vec_t* vec, int len) {
	vec->n = len;
	vec->n_char = len / WORD_SIZE + (len / WORD_SIZE == 0 ? 0 : 1);
	vec->data = (char*)ecc_malloc(vec->n_char * sizeof(char));
}

void copy_gf2_bit_vec(gf2_bit_vec_t* src, gf2_bit_vec_t* dst) {
	ecc_assert(src->n == dst->n, "in %s, trying to copy gf2 vecs in different sizes", __FUNCTION__);
	memcpy(dst->data, src->data, src->n_char * sizeof(char));
}

void copy_gf2_bit_vec_range(gf2_bit_vec_t* src, gf2_bit_vec_t* dst, int start_index, int end_index) {
	ecc_assert(src->n > end_index  && dst->n > end_index, "in %s, trying to copy gf2 vecs where end_index is bigger then lengths", __FUNCTION__);
	ecc_assert(start_index >= 0, "in %s, trying to copy gf2 vecs where start_index is negative", __FUNCTION__);
	for (size_t i = start_index; i <= end_index; i++)
	{
		gf2_bit_vec_set_element(dst, i, gf2_bit_vec_get_element(src, i));
	}
}

void free_gf2_bit_vec(gf2_bit_vec_t* vec) {
	if (vec)
	{
		if (vec->data)
			free(vec->data);
		free(vec);
	}
}


void add_gf2_bit_vec(gf2_bit_vec_t* src1, gf2_bit_vec_t* src2, gf2_bit_vec_t* dst) {
	ecc_assert(src1->n == src2->n, "in %s, src1 and 2 vecs in different sizes", __FUNCTION__);
	ecc_assert(src1->n == dst->n, "in %s, src1 and dst vecs in different sizes", __FUNCTION__);
	/*for (size_t i = 0; i < src1->n_char; i++)
	{
		dst->data[i] = src1->data[i] ^ src2->data[i];
	}*/
	size_t i;
	int n_char = src1->n_char;
	int n_int = n_char / sizeof(int);
	int n_remain_start = n_int * sizeof(int);
	int *dst_data_int = (int*)dst->data, *src1_data_int = (int*)src1->data, *src2_data_int = (int*)src2->data;
	for (i = 0; i < n_int; i++)
	{
		dst_data_int[i] = src1_data_int[i] ^ src2_data_int[i];
	}
	for (i = n_remain_start; i < n_char; i++)
	{
		dst->data[i] = src1->data[i] ^ src2->data[i];
	}
}

void add_gf2_bit_vec_in_place(gf2_bit_vec_t* src1_dst, gf2_bit_vec_t* src2) {
	add_gf2_bit_vec(src1_dst, src2, src1_dst);
}



void apply_inv_perm_gf2_bit_vec(gf2_bit_vec_t *dec_symb, int_vec_t* PHI) {
	gf2_bit_vec_t *tmp;
	int i;

	ecc_assert(dec_symb->n == PHI->n, "Length doesn't match %d vs. %d\n", dec_symb->n, PHI->n);

	malloc_and_init_gf2_bit_vec(&tmp, dec_symb->n);

	for (i = 0; i<dec_symb->n; i++) {
		gf2_bit_vec_set_element(tmp, PHI->data[i], gf2_bit_vec_get_element(dec_symb, i));
	}

	for (i = 0; i<dec_symb->n; i++)
		gf2_bit_vec_set_element(dec_symb, i, gf2_bit_vec_get_element(tmp, i));

	free_gf2_bit_vec(tmp);
}

void gf2_bit_vec_2_gf2_vec(gf2_bit_vec_t* src, gf2_vec_t* dst) {
	int i = 0;
	ecc_assert(src->n == dst->n, "Mismatched vectors %d vs, %d\n", src->n, dst->n);
	for (i = 0; i < dst->n; i++) {
		dst->data[i] = gf2_bit_vec_get_element(src, i);//src->data[i];
	}
}

void print_gf2_bit_vec(gf2_bit_vec_t* vec) {
	int i;

	for (i = 0; i < vec->n; i++) {
		printf("%d", gf2_bit_vec_get_element(vec, i));
	}
	printf("\n");

}



/************* gf2_bit_mtx_t functions *****************/

void malloc_and_init_gf2_bit_mtx(gf2_bit_mtx_t **mtx, int m, int n) {
	*mtx = (gf2_bit_mtx_t*)ecc_malloc(sizeof(gf2_bit_mtx_t));
	init_gf2_bit_mtx(*mtx, m, n);
}

void init_gf2_bit_mtx(gf2_bit_mtx_t *mtx, int m, int n) {
	mtx->n = n;
	mtx->m = m;
	mtx->data = (gf2_bit_vec_t**)ecc_malloc(sizeof(gf2_bit_vec_t*) * m);
	for (size_t i = 0; i < m; i++)
	{
		malloc_and_init_gf2_bit_vec(&(mtx->data[i]), n);
	}
}

void free_gf2_bit_mtx(gf2_bit_mtx_t *x) {
	if (x)
	{
		if (x->data) {
			for (size_t i = 0; i < x->m; i++)
			{
				free_gf2_bit_vec(x->data[i]);
			}
			free(x->data);
		}
		free(x);
	}
}

void copy_gf2_bit_mtx(gf2_bit_mtx_t* src, gf2_bit_mtx_t* dst) {
	ecc_assert(src->n == dst->n, "copy gf2 byte matrix - n not equal (%d vs. %d)\n", src->n ,dst->n);
	ecc_assert(src->m == dst->m, "copy gf2 byte matrix - m not equal  (%d vs. %d)\n", src->m ,dst->m);
	for (size_t i = 0; i < src->m; i++)
	{
		copy_gf2_bit_vec(src->data[i], dst->data[i]);
	}
}

void int_mtx_2_gf2_bit_mtx(int_mtx_t* src, gf2_bit_mtx_t** dst) {
	int n = src->n, m = src->m;
	malloc_and_init_gf2_bit_mtx(dst, m, n);
	gf2_bit_mtx_t* dst_mtx = *dst;
	for (size_t i = 0; i < n * m; i++)
	{
		gf2_bit_mtx_set_element(dst_mtx, i / src->n, i % src->n, src->data[i]);
	}
}

void switch_gf2_bit_mtx_rows(gf2_bit_mtx_t* mtx, int row_a, int row_b) {
	int row_len = mtx->n;
	gf2_bit_vec_t * tmp = mtx->data[row_a];
	mtx->data[row_a] = mtx->data[row_b];
	mtx->data[row_b] = tmp;
}

void encode_gf2_bit_symbols(gf2_bit_vec_t *info, gf2_bit_mtx_t *G, gf2_bit_vec_t *symb) {
	int i, j, sum;

	int n = G->n, m = G->m;

	// Checking all dimensions agree
	ecc_assert(G->n == symb->n, "symbol length doesn't match generator matrix (%d vs. %d)\n", G->n, symb->n);
	ecc_assert(G->m == info->n, "info   length doesn't match generator matrix (%d vs. %d)\n", G->m, info->n);

	//char* G_data = G->data;

	memset(symb->data, 0, symb->n_char * sizeof(char));

	for (j = 0; j < m; j++)	{
		if (1 == gf2_bit_vec_get_element(info, j)) {
			add_gf2_bit_vec_in_place(symb, G->data[j]);
		}
	}
}

void gf2_bit_mtx_gaussian_elimination(gf2_bit_mtx_t* G, int_vec_t* perm1, col_perms_t* col_perms) {
	int i, r, j;
	int temp_val = 0;
	int countperm = 0;
	int M = G->m, N = G->n;

	col_perms_t* col_perms_temp = (col_perms_t*)ecc_malloc(sizeof(col_perms_t));
	init_col_perms_t(col_perms_temp, N);

	//printf("starting OSDLeftSysMatrix. M and N are (%d, %d)\n", M, N);

	for (j = 0; j<M; j++)
	{
		/**** chercher le premier 1, dans la colonne j, situe sur la ligne i>=j ****/
		/**** et permuter les lignes j et i ***********/
		for (i = j; i < M; i++)
		if (gf2_bit_mtx_get_element(G, i, perm1->data[j]))
			break;

		if (i < M) /** on a trouve un 1 sur la colonne j ***/
		{
			switch_gf2_bit_mtx_rows(G, i, j);
		}
		else
			/***** pas de 1 sur la colonne j => permutation de colonnes **/
		{
			for (r = j + 1; r<N; r++)
			for (i = j; i<M; i++)
			if (gf2_bit_mtx_get_element(G, i, perm1->data[r]))
				goto exit_loop;

		exit_loop:

			if ((i >= M) && (r >= N))
			{
				fprintf(stderr, "OSDLeftSysMatrix(): All bits are unset !\n");
				fprintf(stderr, "OSDLeftSysMatrix(): The matrix rank has been modified.\n");
				M = j;
				goto Sortie;
			}

			switch_gf2_bit_mtx_rows(G, i, j);

			/*** permuter les colonnes r et j ****/
			for (i = 0; i<M; i++)
			{
				temp_val = gf2_bit_mtx_get_element(G, i, perm1->data[j]);
				gf2_bit_mtx_set_element(G, i, perm1->data[j], gf2_bit_mtx_get_element(G, i, perm1->data[r]));
				gf2_bit_mtx_set_element(G, i, perm1->data[r], temp_val);
			}
			/*** sauvegarder la permutation ****/
			col_perms_temp->a_indices[countperm] = j;
			col_perms_temp->b_indices[countperm] = r;
			countperm++;

		}/* end of pas de 1 sur colonne j */

		/****** mettre des 0 avant et apres en ajoutant la ligne j ***/
		for (i = 0; i<M; i++)
		if (gf2_bit_mtx_get_element(G, i, perm1->data[j]))
		if (i != j)
		for (r = 0; r < N; r++) {
			temp_val = gf2_bit_mtx_get_element(G, i, perm1->data[r]) ^ gf2_bit_mtx_get_element(G, j, perm1->data[r]);
			gf2_bit_mtx_set_element(G, i, perm1->data[r], temp_val);
		}

	}/* end of j loop */

Sortie:
	if (col_perms->a_indices)
		free(col_perms->a_indices);
	if (col_perms->b_indices)
		free(col_perms->b_indices);
	// printf("countperm is %d\n", countperm);
	init_col_perms_t(col_perms, countperm);
	memcpy(col_perms->a_indices, col_perms_temp->a_indices, sizeof(int)* col_perms->num_of_perms);
	memcpy(col_perms->b_indices, col_perms_temp->b_indices, sizeof(int)* col_perms->num_of_perms);
	free_col_perms(col_perms_temp);
}

void apply_inv_perm_gf2_bit_mtx(gf2_bit_mtx_t* mtx,  int *ord)
{
	gf2_bit_mtx_t *tmp;
	int i, j;

	malloc_and_init_gf2_bit_mtx(&tmp, mtx->m, mtx->n);

	for (i = 0; i<mtx->m; i++)
	for (j = 0; j<mtx->n; j++)
		gf2_bit_mtx_set_element(tmp, i, j, gf2_bit_mtx_get_element(mtx, i, ord[j]));

	copy_gf2_bit_mtx(tmp, mtx);

	free_gf2_bit_mtx(tmp);
}

void print_gf2_bit_mtx(gf2_bit_mtx_t* mtx) {
	int i, j;

	for (i = 0; i<mtx->m; i++)
	{
		printf("| ");
		for (j = 0; j < mtx->n; j++) {
			printf("%d", gf2_bit_mtx_get_element(mtx, i, j));
		}
		printf(" |\n");
	}

}



/************* gf2_byte_mtx_t functions *****************/

void malloc_and_init_gf2_byte_mtx(gf2_byte_mtx_t **mtx, int m, int n) {
	*mtx = (gf2_byte_mtx_t*)ecc_malloc(sizeof(gf2_byte_mtx_t));
	init_gf2_byte_mtx(*mtx, m, n);
}
void init_gf2_byte_mtx(gf2_byte_mtx_t *mtx, int m, int n) {
	mtx->n = n;
	mtx->m = m;
	mtx->data = (char*)ecc_malloc(sizeof(char) * n * m);
}
void free_gf2_byte_mtx(gf2_byte_mtx_t *x) {
	if (x)
	{
		if (x->data)
			free(x->data);
		free(x);
	}
}

void copy_gf2_byte_mtx(gf2_byte_mtx_t* src, gf2_byte_mtx_t* dst) {
	ecc_assert((dst->m == src->m), "Inconsistency with number of lines %d vs. %d\n", dst->m, src->m);
	ecc_assert((dst->n == src->n), "Inconsistency with number of lines %d vs. %d\n", dst->n, src->n);
	memcpy(dst->data, src->data, dst->m * dst->n * sizeof(char));
}

void int_mtx_2_gf2_byte_mtx(int_mtx_t* src, gf2_byte_mtx_t** dst) {
	int n = src->n, m = src->m;
	malloc_and_init_gf2_byte_mtx(dst, m, n);
	gf2_byte_mtx_t* dst_mtx = *dst;
	for (size_t i = 0; i < n * m; i++)
	{
		gf2_byte_mtx_set_element(dst_mtx, i / src->n, i % src->n, src->data[i]);
	}
}
void switch_gf2_byte_mtx_rows(gf2_byte_mtx_t* mtx, int row_a, int row_b) {
	int row_len = mtx->n;
	char* tmp_row = (char*)ecc_malloc(row_len * sizeof(char));
	memcpy(tmp_row, &(mtx->data[row_a * row_len]), row_len * sizeof(char));
	memcpy(&(mtx->data[row_a * row_len]), &(mtx->data[row_b * row_len]), row_len * sizeof(char));
	memcpy(&(mtx->data[row_b * row_len]), tmp_row, row_len * sizeof(char));
	free(tmp_row);
}

void gf2_byte_mtx_gaussian_elimination(gf2_byte_mtx_t* G, int_vec_t* perm1, col_perms_t* col_perms) {
	int i, r, j;
	int temp_val = 0;
	int countperm = 0;
	int M = G->m, N = G->n;

	col_perms_t* col_perms_temp = (col_perms_t*)ecc_malloc(sizeof(col_perms_t));
	init_col_perms_t(col_perms_temp, N);

	//printf("starting OSDLeftSysMatrix. M and N are (%d, %d)\n", M, N);

	for (j = 0; j<M; j++)
	{
		/**** chercher le premier 1, dans la colonne j, situe sur la ligne i>=j ****/
		/**** et permuter les lignes j et i ***********/
		for (i = j; i < M; i++)
		if (gf2_byte_mtx_get_element(G, i, perm1->data[j]))
			break;

		if (i < M) /** on a trouve un 1 sur la colonne j ***/
		{
			switch_gf2_byte_mtx_rows(G, i, j);
		}
		else
			/***** pas de 1 sur la colonne j => permutation de colonnes **/
		{
			for (r = j + 1; r<N; r++)
			for (i = j; i<M; i++)
			if (gf2_byte_mtx_get_element(G, i, perm1->data[r]))
				goto exit_loop;

		exit_loop:

			if ((i >= M) && (r >= N))
			{
				fprintf(stderr, "OSDLeftSysMatrix(): All bits are unset !\n");
				fprintf(stderr, "OSDLeftSysMatrix(): The matrix rank has been modified.\n");
				M = j;
				goto Sortie;
			}

			switch_gf2_byte_mtx_rows(G, i, j);

			/*** permuter les colonnes r et j ****/
			for (i = 0; i<M; i++)
			{
				temp_val = gf2_byte_mtx_get_element(G, i, perm1->data[j]);
				gf2_byte_mtx_set_element(G, i, perm1->data[j], gf2_byte_mtx_get_element(G, i, perm1->data[r]));
				gf2_byte_mtx_set_element(G, i, perm1->data[r], temp_val);
			}
			/*** sauvegarder la permutation ****/
			col_perms_temp->a_indices[countperm] = j;
			col_perms_temp->b_indices[countperm] = r;
			countperm++;

		}/* end of pas de 1 sur colonne j */

		/****** mettre des 0 avant et apres en ajoutant la ligne j ***/
		for (i = 0; i<M; i++)
		if (gf2_byte_mtx_get_element(G, i, perm1->data[j]))
		if (i != j)
		for (r = 0; r < N; r++) {
			temp_val = gf2_byte_mtx_get_element(G, i, perm1->data[r]) ^ gf2_byte_mtx_get_element(G, j, perm1->data[r]);
			gf2_byte_mtx_set_element(G, i, perm1->data[r], temp_val);
		}

	}/* end of j loop */

Sortie:
	if (col_perms->a_indices)
		free(col_perms->a_indices);
	if (col_perms->b_indices)
		free(col_perms->b_indices);
	// printf("countperm is %d\n", countperm);
	init_col_perms_t(col_perms, countperm);
	memcpy(col_perms->a_indices, col_perms_temp->a_indices, sizeof(int)* col_perms->num_of_perms);
	memcpy(col_perms->b_indices, col_perms_temp->b_indices, sizeof(int)* col_perms->num_of_perms);
	free_col_perms(col_perms_temp);
}

void apply_inv_perm_gf2_byte_mtx(gf2_byte_mtx_t* mtx, int *ord) {
	gf2_byte_mtx_t *tmp;
	int i, j;

	malloc_and_init_gf2_byte_mtx(&tmp, mtx->m, mtx->n);

	for (i = 0; i<mtx->m; i++)
	for (j = 0; j<mtx->n; j++)
		gf2_byte_mtx_set_element(tmp, i, j, gf2_byte_mtx_get_element(mtx, i, ord[j]));

	copy_gf2_byte_mtx(tmp, mtx);

	free_gf2_byte_mtx(tmp);
}

void print_gf2_byte_mtx(gf2_byte_mtx_t* mtx) {
	int i, j;

	for (i = 0; i<mtx->m; i++)
	{
		printf("| ");
		for (j = 0; j < mtx->n; j++) {
			printf("%d", gf2_byte_mtx_get_element(mtx, i, j));
		}
		printf(" |\n");
	}
}


void convert_gf2_byte_mtx_2_gf2_bit_mtx(gf2_byte_mtx_t* src, gf2_bit_mtx_t* dst) {
	int n = src->n, m = src->m;
	for (size_t i = 0; i < n * m; i++)
	{
		gf2_bit_mtx_set_element(dst, i / src->n, i % src->n, src->data[i]);
	}
}