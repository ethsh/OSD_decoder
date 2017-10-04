#include <memory.h>
#include "math_funcs.h"


void free_col_perms(col_perms_t *col_perms) {
	if (col_perms) {
		if (col_perms->a_indices)
			free(col_perms->a_indices);
		if (col_perms->b_indices)
			free(col_perms->b_indices);
		free(col_perms);
	}
}

void init_col_perms_t(col_perms_t* col_perms, int len) {
	if (0 == len) {
		col_perms->num_of_perms = 0;
		col_perms->a_indices = NULL;
		col_perms->b_indices = NULL;
		return;
	}
	col_perms->num_of_perms = len;
	col_perms->a_indices = (int*)ecc_malloc(len * sizeof(int));
	col_perms->b_indices = (int*)ecc_malloc(len * sizeof(int));
}


void malloc_and_init_gf2_vec(gf2_vec_t** dst, int len) {
	*dst = (gf2_vec_t*)ecc_malloc(sizeof(gf2_vec_t));
	if (len > 0) {
		init_gf2_vec(*dst, len);
	}
}

void malloc_and_init_real_vec(real_vec_t** dst, int len) {
	*dst = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	if (len > 0) {
		init_real_vec(*dst, len);
	}
}

void malloc_and_init_int_vec(int_vec_t** dst, int len) {
	*dst = (int_vec_t*)ecc_malloc(sizeof(int_vec_t));
	if (len > 0) {
		init_int_vec(*dst, len);
	}
}

void init_real_vec_from_rvec(real_vec_t* dst, real_vec_t* src) {
	dst->n = src->n;
	dst->data = (double*)ecc_malloc(dst->n * sizeof(double));
	memcpy(dst->data, src->data, dst->n * sizeof(double));
}

void init_int_vec_from_ivec(int_vec_t* dst, int_vec_t* src) {
	dst->n = src->n;
	dst->data = (int*)ecc_malloc(dst->n * sizeof(int));
	memcpy(dst->data, src->data, dst->n * sizeof(int));
}

void init_gf2_vec_from_gf2vec(gf2_vec_t* dst, gf2_vec_t* src) {
	dst->n = src->n;
	dst->data = (int*)ecc_malloc(dst->n * sizeof(int));
	memcpy(dst->data, src->data, dst->n * sizeof(int));
}


void copy_int_vec(int_vec_t* dst, int_vec_t* src) {
	memcpy(dst->data, src->data, dst->n * sizeof(int));
}


void apply_inv_perm_real(real_vec_t *lambda, int *ord)
{
	int i;
	double *tmp;

	tmp = (double*)ecc_malloc(sizeof(double)*lambda->n);
  
	for (i=0; i<lambda->n; i++)
		tmp[i] = lambda->data[ord[i]];

	for (i=0; i<lambda->n; i++)
		lambda->data[i] = tmp[i];

	free(tmp);
}

void apply_inv_perm_int_mtx(int_mtx_t* mtx, int *ord)
{
	int_mtx_t *tmp;
	int i, j;
  
	tmp = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  
	init_int_mtx(tmp, mtx->m, mtx->n);
  
	for (i=0; i<mtx->m; i++)
		for (j=0; j<mtx->n; j++)
			put_int_element(tmp, i, j, get_int_element(mtx, i, ord[j]));

	copy_int_mtx(mtx, tmp);

	free(tmp->data);
	free(tmp);
}

void switch_int_mtx_rows(int_mtx_t* mtx, int row_a, int row_b) {
	int row_len = mtx->n;
	int* temp_row = (int*)ecc_malloc(mtx->n * sizeof(int));
	memcpy(temp_row, mtx->data + (row_a * row_len), row_len * sizeof(int));
	memcpy(mtx->data + (row_a * row_len), mtx->data + (row_b * row_len), row_len * sizeof(int));
	memcpy(mtx->data + (row_b * row_len), temp_row, row_len * sizeof(int));
	free(temp_row);
}


void gf2_vec_add(gf2_vec_t* src1, gf2_vec_t* src2, gf2_vec_t* dst) {
	ecc_assert(src1->n == src2->n, "gf2_vec_add - Length doesn't match %d vs. %d\n", src1->n, src2->n);
	ecc_assert(src1->n == dst->n, "gf2_vec_add - Length doesn't match %d vs. %d\n", src1->n, dst->n);
	for (size_t i = 0; i < src1->n; i++)
	{
		dst->data[i] = src1->data[i] ^ src2->data[i];
	}
}


void OSDLeftSysMatrix_ethan(int_mtx_t* G, int_vec_t* perm1, col_perms_t* col_perms) {

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
			if (get_int_element(G, i, perm1->data[j])) 
				break;

		if (i < M) /** on a trouve un 1 sur la colonne j ***/
		{
			switch_int_mtx_rows(G, i, j);
		}
		else
			/***** pas de 1 sur la colonne j => permutation de colonnes **/
		{
			for (r = j + 1; r<N; r++)
				for (i = j; i<M; i++)
					if (get_int_element(G, i, perm1->data[r]))
						goto exit_loop;

		exit_loop:

			if ((i >= M) && (r >= N))
			{
				fprintf(stderr, "OSDLeftSysMatrix(): All bits are unset !\n");
				fprintf(stderr, "OSDLeftSysMatrix(): The matrix rank has been modified.\n");
				M = j;
				goto Sortie;
			}

			switch_int_mtx_rows(G, i, j);

			/*** permuter les colonnes r et j ****/
			for (i = 0; i<M; i++)
			{
				temp_val = get_int_element(G, i, perm1->data[j]);
				put_int_element(G, i, perm1->data[j], get_int_element(G, i, perm1->data[r]));
				put_int_element(G, i, perm1->data[r], temp_val);
			}
			/*** sauvegarder la permutation ****/
			col_perms_temp->a_indices[countperm] = j;
			col_perms_temp->b_indices[countperm] = r;
			countperm++;

		}/* end of pas de 1 sur colonne j */

		/****** mettre des 0 avant et apres en ajoutant la ligne j ***/
		for (i = 0; i<M; i++)
			if (get_int_element(G, i, perm1->data[j]))
				if (i != j)
				for (r = 0; r < N; r++) {
					temp_val = get_int_element(G, i, perm1->data[r]) ^ get_int_element(G, j, perm1->data[r]);
					put_int_element(G, i, perm1->data[r], temp_val);
				}

	}/* end of j loop */

Sortie:
	if (col_perms->a_indices)
		free(col_perms->a_indices);
	if (col_perms->b_indices)
		free(col_perms->b_indices);
	// printf("countperm is %d\n", countperm);
	init_col_perms_t(col_perms, countperm);
	memcpy(col_perms->a_indices, col_perms_temp->a_indices, sizeof(int) * col_perms->num_of_perms);
	memcpy(col_perms->b_indices, col_perms_temp->b_indices, sizeof(int)* col_perms->num_of_perms);
	free_col_perms(col_perms_temp);
}

// This functions creates systematic generator and parity check matrices for a BCH code with length n and generator polynomial represented by 's' - but now, I is in the left in G.
int_mtx_t** create_sys_BCH_ethan(const char* s, int n)
{
	// Asuming s represents the generating polynom as in Lin & Costello appendix
	int n_minus_k, len, k, i, j;
	int_mtx_t **p, *g, *h, *b, *genpoly, *quotient, *reminder, *dividend;
	char q[2];
	// 'p' is a pointer to both generator and parity check matrices
	p = (int_mtx_t**)ecc_malloc(sizeof(int_mtx_t*)* 2);
	g = p[0] = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	h = p[1] = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));

	genpoly = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	quotient = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	reminder = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	dividend = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	b = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));


	len = strlen(s);

	q[0] = s[0]; q[1] = 0;
	// Computing n-k according to the generator polynomial's degree
	n_minus_k = ((atoi(q) & (1 << 2)) ? 3 : ((atoi(q) & (1 << 1)) ? 2 : 1)) + (len - 1) * 3 - 1;
	k = n - n_minus_k;

	// Creating assistive matrices
	init_int_mtx(genpoly, 1, n_minus_k + 1);
	init_int_mtx(quotient, 1, n);
	init_int_mtx(reminder, 1, n);
	init_int_mtx(dividend, 1, n + 1);
	init_int_mtx(b, k, n - k);
	init_int_mtx(g, k, n);
	init_int_mtx(h, n - k, n);

	// forming gen poly (highest degree on the right)
	for (i = 0; i<n_minus_k + 1; i++) {
		//q[0] = s[len -1 - i/3];
		//put_int_element(genpoly, 0, i, ((atoi(q) & (1<< (i%3))) ? 1 : 0));

		q[0] = s[len - 1 - i / 3];
		put_int_element(genpoly, 0, n_minus_k - i, ((atoi(q) & (1 << (i % 3))) ? 1 : 0));

	}

	//print_int_mtx(genpoly);
	//printf("n = %d, k = %d\n", n, k);

	// In the following loop, we generate a matrix contains all remainders of the division of all x^w polynomials, n-k+1<=w<=n-1 with the generator polynomial
	// This will be used for systematic encoding of BCH codewords

	// creating b
	for (i = 0; i<k; i++) {
		dividend->n = n - k + i + 1;
		// Initialize the row to zeros
		for (j = 0; j<n; j++)
			dividend->data[j] = 0;

		// Assign '1' on the MSB to create x^(n-k+1+i)=x^w
		put_int_element(dividend, 0, n - k + i, 1);
		//print_int_mtx(dividend);
		// Calculate the division of x^w with the generator polynomial
		dev_poly(dividend, genpoly, quotient, reminder);

		//print_int_mtx(quotient);
		//print_int_mtx(reminder);
		//getchar();

		for (j = 0; j<n - k; j++)
			put_int_element(b, i, j, get_int_element(reminder, 0, j));

	}
	//print_int_mtx(b);	
	//getchar();

	// Assign the systematic part of the generator matrix (I)
	for (i = 0; i<k; i++)
		put_int_element(g, i, i, 1);
	// Assign the remainder part of the generator matrix to each row (P)
	for (i = 0; i<k; i++)
	for (j = 0; j<n - k; j++)
		put_int_element(g, i, j + k, get_int_element(b, i, j));
	//print_int_mtx(g);
	//getchar();

	// Assign P' to the parity check matrix
	for (i = 0; i<k; i++)
	for (j = 0; j<n - k; j++)
		put_int_element(h, j, i, get_int_element(b, i, j));
	// Assign I to the parity check matrix
	for (i = 0; i<n - k; i++)
		put_int_element(h, i, i + k, 1);

	//print_int_mtx(h);
	//getchar();


	free_int_mtx(genpoly);
	free_int_mtx(quotient);
	free_int_mtx(reminder);
	free_int_mtx(dividend);
	free_int_mtx(b);
	return p;
}


// This function encodes given information bits to a codewords using generator matrix G
void encode_symbols_ethan(gf2_vec_t *info, int_mtx_t *G, gf2_vec_t *symb)
{
	int i, j, sum;

	int n = G->n, m = G->m;
	int* G_data = G->data;

	// Checking all dimensions agree
	ecc_assert(G->n == symb->n, "symbol length doesn't match generator matrix (%d vs. %d)\n", G->n, symb->n);
	ecc_assert(G->m == info->n, "info   length doesn't match generator matrix (%d vs. %d)\n", G->m, info->n);

	memset(symb->data, 0, symb->n * sizeof(int));

	// efficient loop
	for (j = 0; j < m; j++)
	{
		if (1 == info->data[j]) {
			for (i = 0; i < n; i++)
			{
				symb->data[i] += get_int_element(G, j, i);
			}
			symb->data[i] &= 1;
		}
	}

	//// For every code bit
	//for (i = 0; i < n; i++) 	{
	//	sum = 0;
	//	// Standard encoding using matrix multiplication
	//	for (j = 0; j < m; j++)
	//		//sum += info->data[j] * get_int_element(G, j, i);
	//		sum += info->data[j] * G_data[j * n + i];

	//	symb->data[i] = sum % 2; // Result is over GF(2)
	//}
}
