#ifdef _WIN32
#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h>  
#endif

#include<stdio.h>
#include<stdlib.h>
#include"OSDDecoder.h"
#include "time.h"

typedef struct simulation_params_s {
	double start_EBN0_db;
	double end_EBN0_db;
	double EBN0_delta;
	int N;
	int max_errors;
	int_mtx_t *G;
	void (*decode_func)(real_vec_t *channel_word, gf2_vec_t *decoded_word);
}simulation_params_t;

static void dump_double_mtx_text(double_mtx_t *mtx, char *fname)
{
	FILE *fp;
	int type =2;
	int i, j;

	if (!(fp = fopen(fname, "w")))
		ecc_assert(0, "couldnt open %s for dumping", fname);
	
	for (i = 0; i < mtx->n; i++) {
		fprintf(fp, "|\t");
		for (j = 0; j < mtx->m; j++) {
			fprintf(fp, "%lf\t", mtx->data[i * mtx->m + j]);
		}
		fprintf(fp, "|\n");
	}

	fclose(fp);
}


osd_decoder_t* osd_decoder = NULL;
int osd_order = -1;
clock_t time_in_gussian_elimination = 0;

void osd_wrap_decode_func(real_vec_t *channel_word, gf2_vec_t *decoded_word) {
	extern osd_decoder_t* osd_decoder;
	extern int osd_order;
	OSD_decode(osd_decoder, channel_word, osd_order, decoded_word);
}

double_mtx_t * calc_ber(simulation_params_t *simulation_params);
char* get_poly_for_n_k(int n, int k);

int main(int argc, char ** argv) {
	if (argc < 9) {
		printf("Parameters are: n, k, osd_order, start_SNR, end_SNR, delta_SNR, max_iterations, max_errors");
		return -1;
	}
	simulation_params_t *simulation_params = NULL;
	int_mtx_t **pp = NULL;
	int_mtx_t *G = NULL, *H = NULL;
	int n, k;
	sscanf(argv[1], "%d", &n);
	sscanf(argv[2], "%d", &k);
	char *poly = get_poly_for_n_k(n, k);;
	if (poly == NULL) {
		printf("BCH(%d, %d) not supported. SORRY!", n, k);
		return -1;
	}

	double_mtx_t* ber = NULL;

	extern int osd_order;
	extern osd_decoder_t* osd_decoder;
	sscanf(argv[3], "%d", &osd_order);

	pp = create_sys_BCH_ethan(poly, n);
	G = pp[0];
	H = pp[1];
	k = G->m;

	init_osd_decoder(&osd_decoder, G, osd_order);

	simulation_params = (simulation_params_t*)ecc_malloc(sizeof(simulation_params_t));
	sscanf(argv[4], "%lf", &(simulation_params->start_EBN0_db));
	sscanf(argv[5], "%lf", &(simulation_params->end_EBN0_db));
	sscanf(argv[6], "%lf", &(simulation_params->EBN0_delta));
	sscanf(argv[7], "%d", &(simulation_params->N));
	sscanf(argv[8], "%d", &(simulation_params->max_errors));

	simulation_params->G = G;
	simulation_params->decode_func = osd_wrap_decode_func;

	ber = calc_ber(simulation_params);

	free_double_mtx(ber);
	free(simulation_params);
	free_int_mtx(G);
	free_int_mtx(H);
	free(pp);
	free_osd_decoder(osd_decoder);
	
#ifdef _WIN32
	_CrtDumpMemoryLeaks();
#endif
	getchar();
}

char* get_poly_for_n_k(int n, int k) {
	char *poly = NULL;
	if (n == 63) {
		if (k == 36) {
			poly = "1033500423";
		}
		else if (k == 45) {
			poly = "1701317";
		}
		else {
			return NULL;
		}
	}
	else if (n == 127) {
		if (k == 64){
			poly = "1206534025570773100045";
		}
		else if (k == 99) {
			poly = "3447023271";
		}
		else if (k == 85) {
			poly = "130704476322273";
		}
		else if (k == 92) {
			poly = "624730022327";
		}
		else {
			return NULL;
		}
	}
	else {
		return NULL;
	}
	return poly;
}

double_mtx_t * calc_ber(simulation_params_t *simulation_params)
{
	extern clock_t time_in_gussian_elimination;
	int_mtx_t *G;
	int N;
	double start_EBN0_db, end_EBN0_db, EBN0_delta;
	double n_errors, n_errors_uc;
	double word_errors;
	int word_err_flag;
	int i, j, iter, length, dec_uc;
	int erf_n_errors_uc, erf_n_errors;
	int max_word_errors;

	clock_t t_clk, sum_clk;

	double_mtx_t *ber;
	gf2_vec_t *info, *symb, *dec_symb;
	real_vec_t *data, *data_uc;

	N = simulation_params->N;
	end_EBN0_db = simulation_params->end_EBN0_db;
	start_EBN0_db = simulation_params->start_EBN0_db;
	EBN0_delta = simulation_params->EBN0_delta;
	G = simulation_params->G;
	max_word_errors = simulation_params->max_errors;

	length = ((end_EBN0_db - start_EBN0_db) / EBN0_delta) + 1;
	ber = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
	// Note there are 3 values allocated for each SNR value: The SNR value itself, the BER using the given code, and the BER when transmitting uncoded data
	init_double_mtx(ber, 6, length);

	malloc_and_init_gf2_vec(&info, G->m); // Vector that includes information bits
	malloc_and_init_gf2_vec(&dec_symb, G->n); // Vector that includes decoded bits (after channel transmission and decoding)
	malloc_and_init_gf2_vec(&symb, G->n); // Vector that includes encoded bits
	malloc_and_init_real_vec(&data, G->n); // Vector that includes real values - the measurements received after transmission over AWGN channel
	malloc_and_init_real_vec(&data_uc, G->n);
	
	// The following loop runs on all SNR values
	for (i = 0; i<length; i++) {
		printf("TG, Calculating SNR = %g, n=%d, k=%d\n", start_EBN0_db + EBN0_delta*i, G->n, G->m);

		put_double_element(ber, 0, i, start_EBN0_db + EBN0_delta*i);
		srand( (unsigned)time( NULL ) );
		
		// Initialize error counters
		n_errors = 0;
		n_errors_uc = 0;
		word_errors = 0;
		erf_n_errors_uc = erf_n_errors = 0;
		sum_clk = 0;
		time_in_gussian_elimination = 0;

		// Simulate N codewords
		for (iter = 0; iter<N && word_errors < max_word_errors; iter++) {
			word_err_flag = 0;

			// Print status every 100 codewords
			if (0 == (iter % 100)) printf("TG, SNR = %g, iter = %d out of %d, n=%d, k=%d, time=%lf\n", start_EBN0_db + EBN0_delta*i, iter, N, G->n, G->m, sum_clk / (double)CLOCKS_PER_SEC);
			// Generating information bits
			init_symbol_data(info);
			// Encoding information bits to a codeword
			encode_symbols(info, G, symb);
			// Transmit over AWGN channel   
			init_measured_data(symb, data, start_EBN0_db + EBN0_delta*i, ((double)(G->m) / (double)G->n)); // Transmission of the codeword
			init_measured_data(info, data_uc, start_EBN0_db + EBN0_delta*i, 1.0); // Transmission of only the information bits

			t_clk = clock();
			simulation_params->decode_func(data, dec_symb);
			sum_clk += (clock() - t_clk);

			// Now counting the number of error bits
			// coded data
			for (j = 0; j<G->n; j++) {
				if (symb->data[j] != dec_symb->data[j]) {
					n_errors++;
					if (!word_err_flag) {
						word_errors++;
					}
					word_err_flag = 1;
				}
			}
			// uncoded data
			for (j = 0; j<G->m; j++) {
				dec_uc = (data_uc->data[j] > 0) ? 1 : 0;
				n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
				erf_n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
			}
		}
		// Storing BER for the coded and uncoded transmissions
		put_double_element(ber, 1, i, n_errors / (double)(iter * G->n));
		put_double_element(ber, 2, i, word_errors / (double)(iter));
		put_double_element(ber, 3, i, sum_clk / CLOCKS_PER_SEC / (double)(iter));
		put_double_element(ber, 4, i, time_in_gussian_elimination / CLOCKS_PER_SEC / (double)(iter));
		put_double_element(ber, 5, i, n_errors_uc / (double)(iter * G->m));
	}

	free_gf2_vec(info);
	free_gf2_vec(symb);
	free_gf2_vec(dec_symb);
	free_real_vec(data);
	free_real_vec(data_uc);

	// Printing the BER matrix
	print_double_mtx(ber);
	dump_double_mtx_text(ber, "BER.txt");
	return ber;
}
