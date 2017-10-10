#include <string.h>
#include <math.h>
#include <time.h>
#include "OSDDecoder.h"

typedef struct error_patterns_s{
	//gf2_bit_vec_t** all_error_pattern;
	int_vec_t** all_error_pattern; // error pattern is the locations of errors.
	int num_of_errors;
	int order;
	int k;
} error_patterns_t;

typedef struct osd_decoder_s{
	int n, k;
	
	gf2_byte_mtx_t *G, *G_perm;
	gf2_bit_mtx_t *G_perm_bit;

	col_perms_t* col_perms;
	int_vec_t* perm1;
	int_vec_t* total_perm;

	real_vec_t* curr_soft_symb;

	real_vec_t* curr_soft_symb_perm;
	real_vec_t* curr_confidence_perm;

	int max_pre_calced_order;
	error_patterns_t** pre_calced_errors;

	double curr_best_distance;
	gf2_bit_vec_t* curr_best_error_pattern;

	gf2_bit_vec_t* curr_hard_syms_codeword_perm;
	gf2_bit_vec_t* decoded_codeword;
} osd_decoder_t;


#ifdef _WIN32
real_vec_t* curr_conf_values_for_compare = NULL;

int cmp_decreasing_order(const void *arg1, const void *arg2)
{
	extern real_vec_t* curr_conf_values_for_compare;
	int* int_index_x = (int*)arg1;
	int* int_index_y = (int*)arg2;
	double* float_values = (double*)curr_conf_values_for_compare->data;

	double value_x = fabs(float_values[*int_index_x]);
	double value_y = fabs(float_values[*int_index_y]);

	if (value_x > value_y) return (-1);
	if (value_x < value_y) return (1);
	return 0;
}

#else
/// comparison func for qsort_r (in order to sort the indices by the soft symbols)
/// as taken from osd_decoding_biawgn.c
int cmp_decreasing_order(const void *index_x, const void *index_y, void *values) {
	int* int_index_x = (int*)index_x;
	int* int_index_y = (int*)index_y;
	double* float_values = (double*)values;

	double value_x = fabs(float_values[*int_index_x]);
	double value_y = fabs(float_values[*int_index_y]);

	if (value_x > value_y) return (-1);
	if (value_x < value_y) return (1);
	return 0;
}
#endif

void init_error_pattern(error_patterns_t* error_patterns, int order, int k);
void free_error_pattern(error_patterns_t* error_patterns);

void calc_confidence_permutaion(osd_decoder_t* osd_decoder);

void gaussian_eliminate(osd_decoder_t* osd_decoder);
void create_combined_perm(int_vec_t* perm1, col_perms_t* col_perms, int_vec_t* total_perm);

void apply_permutations(osd_decoder_t* osd_decoder);

void osd_perform_decoding(osd_decoder_t* osd_decoder, int full_osd_order);
void osd_perform_decoding_for_single_order(osd_decoder_t* osd_decoder, int order, gf2_bit_vec_t *initial_coded_word, gf2_bit_vec_t* out_best_error_pattern, double* out_best_error_distance);

int get_num_of_error_patterns(int osd_order, int k);


void create_all_error_patterns_recursive(int_vec_t** all_error_pattern, int osd_order, int k);
void recursive_create_error_pattern(int_vec_t** all_error_pattern, int* total_errors_allocated, int_vec_t* curr_indices, int indices_to_choose, int initial_index, int max_index);

void osd_retrieve_original_word(osd_decoder_t* osd_decoder, gf2_vec_t* decoded_word);



void init_osd_decoder(osd_decoder_t** osd_decoder, int_mtx_t* G, int pre_calc_error_order) {
	*osd_decoder = (osd_decoder_t*)ecc_malloc(sizeof(osd_decoder_t));

	osd_decoder_t* osd_decoder_p = *osd_decoder;

	osd_decoder_p->n = G->n;
	osd_decoder_p->k = G->m;

	int_mtx_2_gf2_byte_mtx(G, &(osd_decoder_p->G));
	malloc_and_init_gf2_byte_mtx(&(osd_decoder_p->G_perm), osd_decoder_p->k, osd_decoder_p->n);
	malloc_and_init_gf2_bit_mtx(&(osd_decoder_p->G_perm_bit), osd_decoder_p->k, osd_decoder_p->n);

	osd_decoder_p->col_perms = (col_perms_t*)ecc_malloc(sizeof(col_perms_t));
	init_col_perms_t(osd_decoder_p->col_perms, 0);

	malloc_and_init_int_vec(&(osd_decoder_p->perm1), osd_decoder_p->n);
	malloc_and_init_int_vec(&(osd_decoder_p->total_perm), osd_decoder_p->n);

	malloc_and_init_real_vec(&(osd_decoder_p->curr_soft_symb), osd_decoder_p->n);
	malloc_and_init_real_vec(&(osd_decoder_p->curr_soft_symb_perm), osd_decoder_p->n);
	malloc_and_init_real_vec(&(osd_decoder_p->curr_confidence_perm), osd_decoder_p->n);

	osd_decoder_p->curr_best_distance = 1e10;
	malloc_and_init_gf2_bit_vec(&(osd_decoder_p->curr_best_error_pattern), osd_decoder_p->k);
	malloc_and_init_gf2_bit_vec(&(osd_decoder_p->decoded_codeword), osd_decoder_p->n);
	malloc_and_init_gf2_bit_vec(&(osd_decoder_p->curr_hard_syms_codeword_perm), osd_decoder_p->n);

	osd_decoder_p->max_pre_calced_order = pre_calc_error_order;
	if (pre_calc_error_order >= 0) {
		osd_decoder_p->pre_calced_errors = (error_patterns_t**)ecc_malloc((pre_calc_error_order + 1) * sizeof(error_patterns_t*));
		for (size_t i = 0; i <= osd_decoder_p->max_pre_calced_order; i++)
		{
			osd_decoder_p->pre_calced_errors[i] = (error_patterns_t*)ecc_malloc(sizeof(error_patterns_t));
			init_error_pattern(osd_decoder_p->pre_calced_errors[i], i, osd_decoder_p->k);
		}
	}
}

void free_osd_decoder(osd_decoder_t* osd_decoder) {
	free_gf2_byte_mtx(osd_decoder->G);
	free_gf2_byte_mtx(osd_decoder->G_perm);
	free_gf2_bit_mtx(osd_decoder->G_perm_bit);

	free_col_perms(osd_decoder->col_perms);

	free_int_vec(osd_decoder->perm1);
	free_int_vec(osd_decoder->total_perm);

	free_real_vec(osd_decoder->curr_soft_symb);
	free_real_vec(osd_decoder->curr_soft_symb_perm);
	free_real_vec(osd_decoder->curr_confidence_perm);

	free_gf2_bit_vec(osd_decoder->curr_best_error_pattern);
	free_gf2_bit_vec(osd_decoder->decoded_codeword);
	free_gf2_bit_vec(osd_decoder->curr_hard_syms_codeword_perm);

	if (osd_decoder->max_pre_calced_order >= 0) {
		for (size_t i = 0; i <= osd_decoder->max_pre_calced_order; i++)
		{
			free_error_pattern(osd_decoder->pre_calced_errors[i]);
		}
		free(osd_decoder->pre_calced_errors);
	}

	free(osd_decoder);

#ifdef _WIN32
	extern real_vec_t* curr_conf_values_for_compare;
	if (NULL != curr_conf_values_for_compare) {
		free_real_vec(curr_conf_values_for_compare);
	}
#endif

}

void OSD_decode(osd_decoder_t* osd_decoder, real_vec_t* soft_symb, int osd_order, gf2_vec_t* decoded_word) {
	copy_rvec(soft_symb, osd_decoder->curr_soft_symb);
	// Sort rec_codeword
	calc_confidence_permutaion(osd_decoder);

	// gaussian eliminate G
	gaussian_eliminate(osd_decoder);

	// apply permutations
	apply_permutations(osd_decoder);

	// decode
	osd_perform_decoding(osd_decoder, osd_order);

	// retrieve the original codeword
	osd_retrieve_original_word(osd_decoder, decoded_word);
}



void init_error_pattern(error_patterns_t* error_patterns, int order, int k) {
	error_patterns->k = k;
	error_patterns->order = order;
	error_patterns->num_of_errors = get_num_of_error_patterns(error_patterns->order, error_patterns->k);
	//error_patterns->all_error_pattern = (gf2_bit_vec_t**)ecc_malloc(error_patterns->num_of_errors * sizeof(gf2_bit_vec_t*));
	error_patterns->all_error_pattern = (int_vec_t**)ecc_malloc(error_patterns->num_of_errors * sizeof(int_vec_t*));

	create_all_error_patterns_recursive(error_patterns->all_error_pattern, error_patterns->order, error_patterns->k);
}

void free_error_pattern(error_patterns_t* error_patterns) {
	for (size_t i = 0; i < error_patterns->num_of_errors; i++)
	{
		//free_gf2_bit_vec(error_patterns->all_error_pattern[i]);
		free_int_vec(error_patterns->all_error_pattern[i]);
	}
	free(error_patterns->all_error_pattern);
	free(error_patterns);
}

int get_num_of_error_patterns(int osd_order, int k) {
	int num_of_patterns = 1;
	for (size_t i = 0; i < osd_order; i++)
	{
		num_of_patterns *= (k - i);
	}
	for (size_t i = 1; i <= osd_order; i++)
	{
		num_of_patterns /= i;
	}
	return num_of_patterns;
}


void create_all_error_patterns_recursive(int_vec_t** all_error_pattern, int osd_order, int k) {
	int total_errors_allocated = 0;
	int_vec_t* curr_indices = (int_vec_t*)ecc_malloc(sizeof(int_vec_t));
	if (osd_order > 0)
		init_int_vec(curr_indices, osd_order);
	else {
		curr_indices->n = 0;
	}
	recursive_create_error_pattern(all_error_pattern, &total_errors_allocated, curr_indices, osd_order, 0, k);
	free_int_vec(curr_indices);
}

void recursive_create_error_pattern(int_vec_t** all_error_pattern, int* total_errors_allocated, int_vec_t* curr_indices, int indices_to_choose, int initial_index, int max_index) {
	if (0 == indices_to_choose){
		// create error pattern
		int curr_error_index = *total_errors_allocated;
		//malloc_and_init_gf2_bit_vec(&(all_error_pattern[curr_error_index]), max_index);
		all_error_pattern[curr_error_index] = (int_vec_t*)ecc_malloc(sizeof(int_vec_t));
		if (curr_indices->n > 0) {
			init_int_vec(all_error_pattern[curr_error_index], curr_indices->n);
			copy_int_vec(all_error_pattern[curr_error_index], curr_indices);
		}
		else {
			all_error_pattern[curr_error_index]->n = 0;
		}
		
		/*for (size_t i = 0; i < curr_indices->n; i++)
		{
			gf2_bit_vec_set_element(all_error_pattern[curr_error_index], curr_indices->data[i], 1);
		}*/
		
		*total_errors_allocated = curr_error_index + 1;
		return;
	}
	for (size_t i = initial_index; i < max_index - indices_to_choose + 1; i++)
	{
		curr_indices->data[indices_to_choose - 1] = i;
		recursive_create_error_pattern(all_error_pattern, total_errors_allocated, curr_indices, indices_to_choose - 1, i + 1, max_index);
	}
}

void calc_confidence_permutaion(osd_decoder_t* osd_decoder) {
	int i = 0;
	for (i = 0; i < osd_decoder->n; i++) {
		osd_decoder->perm1->data[i] = i;
	}

#ifdef _WIN32
	extern real_vec_t* curr_conf_values_for_compare;
	if (NULL == curr_conf_values_for_compare) {
		curr_conf_values_for_compare = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
		init_real_vec(curr_conf_values_for_compare, osd_decoder->n);
	}

	copy_rvec(osd_decoder->curr_soft_symb, curr_conf_values_for_compare);

	qsort((void *)osd_decoder->perm1->data, (size_t)osd_decoder->n, sizeof(int), cmp_decreasing_order);
#else
	qsort_r(osd_decoder->perm1->data, osd_decoder->n, sizeof(int), cmp_decreasing_order, osd_decoder->curr_soft_symb->data);
#endif
}


void gaussian_eliminate(osd_decoder_t* osd_decoder) {
	extern clock_t time_in_gussian_elimination;
	clock_t temp_time_in_gussian_elimination;
	temp_time_in_gussian_elimination = clock();

	copy_gf2_byte_mtx(osd_decoder->G, osd_decoder->G_perm);
	gf2_byte_mtx_gaussian_elimination(osd_decoder->G_perm, osd_decoder->perm1, osd_decoder->col_perms);

	time_in_gussian_elimination += (clock() - temp_time_in_gussian_elimination);
}

void create_combined_perm(int_vec_t* perm1, col_perms_t* col_perms, int_vec_t* total_perm) {
	int i = 0, temp_index = 0, index_a = 0, index_b = 0;
	copy_int_vec(total_perm, perm1);
	for (i = 0; i < col_perms->num_of_perms; i++) {
		temp_index = total_perm->data[col_perms->a_indices[i]];
		total_perm->data[col_perms->a_indices[i]] = total_perm->data[col_perms->b_indices[i]];
		total_perm->data[col_perms->b_indices[i]] = temp_index;
		/* index_a = total_perm->data[col_perms->a_indices[i]];
		index_b = total_perm->data[col_perms->b_indices[i]];
		
		temp_index = total_perm->data[index_a];
		total_perm->data[index_a] = total_perm->data[index_b];
		total_perm->data[index_b] = temp_index; */
	}
}


void apply_permutations(osd_decoder_t* osd_decoder) {
	create_combined_perm(osd_decoder->perm1, osd_decoder->col_perms, osd_decoder->total_perm);
	copy_rvec(osd_decoder->curr_soft_symb, osd_decoder->curr_soft_symb_perm);
	apply_inv_perm_real(osd_decoder->curr_soft_symb_perm, osd_decoder->total_perm->data);
	//apply_inv_perm_gf2_bit_mtx(osd_decoder->G_perm, osd_decoder->perm1->data);
	apply_inv_perm_gf2_byte_mtx(osd_decoder->G_perm, osd_decoder->perm1->data);
	convert_gf2_byte_mtx_2_gf2_bit_mtx(osd_decoder->G_perm, osd_decoder->G_perm_bit);
	for (size_t i = 0; i < osd_decoder->n; i++)
	{
		gf2_bit_vec_set_element(osd_decoder->curr_hard_syms_codeword_perm, i, 1 ? osd_decoder->curr_soft_symb_perm->data[i] > 0 : 0);
	}
}


void osd_perform_decoding(osd_decoder_t* osd_decoder, int full_osd_order) {
	double temp_error_distance;
	osd_decoder->curr_best_distance = 1e10;

	gf2_bit_vec_t *initial_coded_word = NULL, *temp_mrb = NULL;
	malloc_and_init_gf2_bit_vec(&initial_coded_word, osd_decoder->n);
	malloc_and_init_gf2_bit_vec(&temp_mrb, osd_decoder->k);
	copy_gf2_bit_vec_range(osd_decoder->curr_hard_syms_codeword_perm, temp_mrb, 0, osd_decoder->k - 1);
	encode_gf2_bit_symbols(temp_mrb, osd_decoder->G_perm_bit, initial_coded_word);


	gf2_bit_vec_t* temp_best_error_pattern = NULL;
	malloc_and_init_gf2_bit_vec(&temp_best_error_pattern, osd_decoder->k);
	for (size_t i = 0; i <= full_osd_order; i++) {
		osd_perform_decoding_for_single_order(osd_decoder, i, initial_coded_word, temp_best_error_pattern, &temp_error_distance);
		if (temp_error_distance < osd_decoder->curr_best_distance) {
			osd_decoder->curr_best_distance = temp_error_distance;
			copy_gf2_bit_vec(temp_best_error_pattern, osd_decoder->curr_best_error_pattern);
		}
		if (0 == temp_error_distance) {
			// found the best error pattern :)
			break;
		}
	}

	free_gf2_bit_vec(temp_best_error_pattern);
	free_gf2_bit_vec(initial_coded_word);
	free_gf2_bit_vec(temp_mrb);
}

void osd_perform_decoding_for_single_order(osd_decoder_t* osd_decoder, int order, gf2_bit_vec_t *initial_coded_word, gf2_bit_vec_t* out_best_error_pattern, double* out_best_error_distance) {
	double temp_error_distance = 0;
	*out_best_error_distance = 1e10;
	int n = osd_decoder->n, k = osd_decoder->k;

	// TEST
	double* tmp_abs_soft_ptr = osd_decoder->curr_soft_symb_perm->data;
	for (size_t i = 0; i < osd_decoder->curr_soft_symb_perm->n; i++)
	{
		tmp_abs_soft_ptr[i] = fabs(tmp_abs_soft_ptr[i]);
	}
	double curr_best_error = 1e10;

	int_vec_t *temp_error_pattern = NULL;
	gf2_bit_vec_t *temp_full_symbols = NULL, *error_locations_vec = NULL;
	malloc_and_init_gf2_bit_vec(&temp_full_symbols, osd_decoder->n);
	malloc_and_init_gf2_bit_vec(&error_locations_vec, osd_decoder->n);

	error_patterns_t* error_pattern = osd_decoder->pre_calced_errors[order];

	for (size_t i = 0; i < error_pattern->num_of_errors; i++)
	{
		copy_gf2_bit_vec(initial_coded_word, temp_full_symbols);
		temp_error_pattern = error_pattern->all_error_pattern[i];

		for (size_t j = 0; j < temp_error_pattern->n; j++)
		{
			add_gf2_bit_vec_in_place(temp_full_symbols, osd_decoder->G_perm_bit->data[temp_error_pattern->data[j]]);
		}

		add_gf2_bit_vec(temp_full_symbols, osd_decoder->curr_hard_syms_codeword_perm, error_locations_vec);

		for (size_t j = 0; j < temp_error_pattern->n; j++)
		{
			temp_error_distance += tmp_abs_soft_ptr[temp_error_pattern->data[j]];
		}

		curr_best_error = *out_best_error_distance;

		for (size_t l = k; l < n; l++)
		{
			
			if (0 != gf2_bit_vec_get_element(error_locations_vec, l)) {
				temp_error_distance += tmp_abs_soft_ptr[l];
				if (temp_error_distance > curr_best_error) {
					break;
				}
			}
		}

		if (temp_error_distance < *out_best_error_distance) {
			*out_best_error_distance = temp_error_distance;
			memset(out_best_error_pattern->data, 0, out_best_error_pattern->n_char * sizeof(char));
			for (size_t i = 0; i < temp_error_pattern->n; i++)
			{
				gf2_bit_vec_set_element(out_best_error_pattern, temp_error_pattern->data[i], 1);
			}
			//copy_gf2_bit_vec(temp_error_pattern, out_best_error_pattern);
		}
		if (0 == temp_error_distance) {
			// found the best error pattern :)
			break;
		}
		temp_error_distance = 0;
		//add_gf2_bit_vec_in_place(temp_info_symbols, temp_error_pattern);
	}

	free_gf2_bit_vec(error_locations_vec);
	free_gf2_bit_vec(temp_full_symbols);
}


void osd_retrieve_original_word(osd_decoder_t* osd_decoder, gf2_vec_t* decoded_word) {
	gf2_bit_vec_t* temp_info_symbols = NULL;
	malloc_and_init_gf2_bit_vec(&temp_info_symbols, osd_decoder->k);

	copy_gf2_bit_vec_range(osd_decoder->curr_hard_syms_codeword_perm, temp_info_symbols, 0, osd_decoder->k - 1);

	add_gf2_bit_vec_in_place(temp_info_symbols, osd_decoder->curr_best_error_pattern);

	encode_gf2_bit_symbols(temp_info_symbols, osd_decoder->G_perm_bit, osd_decoder->decoded_codeword);

	apply_inv_perm_gf2_bit_vec(osd_decoder->decoded_codeword, osd_decoder->total_perm);

	free_gf2_bit_vec(temp_info_symbols);

	gf2_bit_vec_2_gf2_vec(osd_decoder->decoded_codeword, decoded_word);
}
