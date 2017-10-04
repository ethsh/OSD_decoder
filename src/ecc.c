#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <mcheck.h> 
#include <fstream> // For reading weight matrices
#include <cstdlib> // For reading weight matrices
#include <iostream> // For reading weight matrices
#include "ecc.h"
#include "scc.h"
#include "ldpc_class.h"
#include "timer.h"

#define sqr(x) ((x)*(x))

#define TEST_SC		0

#define rdtscll(val) do {                       \
     unsigned int a,d; \
     asm volatile("rdtsc" : "=a" (a), "=d" (d)); \
     (val) = ((unsigned long long)a) | (((unsigned long long)d)<<32); \
} while(0)

// #define __Eff_mRRD__
#define __use_soft_TG__


typedef struct ecc_s {
	int_mtx_t* (*create_permutation)();
	void* (*perm_get)(void* S, int pos);
	void  (*perm_operate)(void* S, int i, int j);
	int I1, I2, I3; 
	int_mtx_t *G;
	int_mtx_t *H;
	int t;
	int H_tx_present;
	int_mtx_t *H_tx;
	int N;
	int TG_iter;
	double start_EBN0_db;
	double end_EBN0_db;
	double EBN0_delta;
	double tg_alpha;
	double rrd_alpha;
	char ml_fname[MAX_STR];
	char tg_fname[MAX_STR];
	char rrd_fname[MAX_STR];
	char mbbp_fname[MAX_STR];
	char test_fname[MAX_STR];
	char jiang_fname[MAX_STR];
	char dbd_fname[MAX_STR];
	char alpha_fname[MAX_STR];
	double alpha_success;
	int use_genie;
	int num_edges_in_H;
	double_mtx_t *hidden_layers_weights;
	double_mtx_t *out_weights;
	int_mtx_t *accumulated_rows_sum_matrix;
	int_mtx_t *accumulated_columns_sum_matrix;

	int 	mbbp_i, 
		mbbp_l, 
		mbbp_c_rest;  // restriction on last column; used in golay case
	int_mtx_t *cgg;	


	int_mtx_t *cyclic_permutation;

}ecc_t;

typedef struct erf_s {
	double eb_n0;
	int nr_iterations, g_n, g_m, uncoded_errors, coded_errors, uncoded_frame_erros, coded_frame_errors;
	double uncoded_error_rate, coded_error_rate;
	double average_iter, average_iter_sqr, average_time, average_tg_time;
} erf_t;


#define ERF_QUNT 100
#define MAX_ERFS 20

void update_erf_file(char *ifname, double eb_n0, int quant, int g_n, int g_m, int errors_uc, int errors_c){
	FILE *fp;
	erf_t erf_array[MAX_ERFS], erf_array_new[MAX_ERFS];
	int pos;
	int cnt, found_flag, nr_elements;
	double tmp1, tmp2;
	char str[200], fname[200], *p;
	
	// Check for expected file name
	strcpy(fname, ifname);
	p = strstr(fname, ".dat");
	ecc_assert(p, "Fname %s doesn't contain suffix .dat\n", ifname);
	strcpy(p, ".erf");
	
	fp = fopen(fname, "rt");
	if (fp) {
		// file exist
		
		// analyze the file
		// First we check how many lines the file contains
		// Every line represents a different SNR
		cnt =0;
		while ((cnt<MAX_ERFS) && (!feof(fp))) {
			
			fscanf(fp, "%lg, %d, %d, %d, %d, %d, %lg, %lg\n", 
			        &erf_array[cnt].eb_n0, 
				&erf_array[cnt].nr_iterations, 
				&erf_array[cnt].g_n, 
				&erf_array[cnt].g_m, 
				&erf_array[cnt].uncoded_errors,
				&erf_array[cnt].coded_errors,
				&tmp1, &tmp2);
			cnt++;

		}
		nr_elements = cnt;
		
		ecc_assert(feof(fp), "File %s is to big; suspicious\n", fname);
		// Now we will check if the current SNR is already logged to the erf file
		cnt = 0; 
		found_flag = 0;
		while ((!found_flag) && (cnt<nr_elements)) {
			if (erf_array[cnt].eb_n0 == eb_n0) {
				found_flag = 1;
				// Then we have found a line in the file that belongs to the same SNR that is simulated now
				// First we check for consistency of the code's parameter (to avoid weird bugs)
				ecc_assert(erf_array[cnt].g_n == g_n, "Inconsistent n (%d vs. %d)\n", erf_array[cnt].g_n, g_n);
				ecc_assert(erf_array[cnt].g_m == g_m, "Inconsistent n (%d vs. %d)\n", erf_array[cnt].g_m, g_m);
				// Now we update the previously logged performance
				erf_array[cnt].nr_iterations += quant;
				erf_array[cnt].uncoded_errors += errors_uc;
				erf_array[cnt].coded_errors += errors_c;
				// And save the updated logs to a temporary array
				memcpy(&erf_array_new[0], &erf_array[0], sizeof(erf_t)*nr_elements);
			} else 
				cnt ++;
		}
		if (!found_flag) {
			// If the current SNR hasn't been logged yet, check where should it be logged and make room for current SNR's info
			pos = -1;
			if (eb_n0 < erf_array[0].eb_n0) {
				memcpy(&erf_array_new[1], &erf_array[0], sizeof(erf_t)*nr_elements);
				pos = 0;
			} else if (eb_n0 > erf_array[nr_elements-1].eb_n0) {
				memcpy(&erf_array_new[0], &erf_array[0], sizeof(erf_t)*nr_elements);
				pos = nr_elements;
			} else {
				for (cnt=0; cnt<nr_elements; cnt++) {
					if ((eb_n0 > erf_array[cnt].eb_n0) && (eb_n0 < erf_array[cnt+1].eb_n0)) {
						memcpy(&erf_array_new[0], &erf_array[0], sizeof(erf_t)*(cnt+1));
						memcpy(&erf_array_new[cnt+2], &erf_array[cnt+1], sizeof(erf_t)*(nr_elements - (cnt+1)));
						pos = cnt+1;
					}
				}
			}
			ecc_assert(pos != -1, "Error: Bug in design\n");
			// There is room now for current SNR statistics in the temporary array, update it
			nr_elements++;			
			erf_array_new[pos].eb_n0 = eb_n0;
			erf_array_new[pos].nr_iterations = quant;
			erf_array_new[pos].g_n = g_n;
			erf_array_new[pos].g_m = g_m; 
			erf_array_new[pos].uncoded_errors = errors_uc;
			erf_array_new[pos].coded_errors = errors_c;			
		} 
		fclose(fp);
		// saving back up for the file
		sprintf(str, "mv %s %s~", fname, fname);
		system(str);
		// Now save the updated data to the statistics logs file 
		fp = fopen(fname, "wt");
		for (cnt=0; cnt<nr_elements; cnt++)
			fprintf(fp, "%lg, %d, %d, %d, %d, %d, %lg, %lg\n", 
				erf_array_new[cnt].eb_n0, 
				erf_array_new[cnt].nr_iterations, 
				erf_array_new[cnt].g_n , 
				erf_array_new[cnt].g_m, 
				erf_array_new[cnt].uncoded_errors, 
				erf_array_new[cnt].coded_errors, 
				((double)erf_array_new[cnt].uncoded_errors)/ ((double)(erf_array_new[cnt].nr_iterations*erf_array_new[cnt].g_m)), 
				((double)erf_array_new[cnt].coded_errors) / ((double)(erf_array_new[cnt].nr_iterations*erf_array_new[cnt].g_n)));
		fclose(fp);
		//printf("Check file\n");
		//getchar();
		
	} else { 
		// file doesn't exist
		fp = fopen(fname, "w");
		fprintf(fp, "%lg, %d, %d, %d, %d, %d, %lg, %lg\n", 
			eb_n0, quant, g_n, g_m, errors_uc, errors_c, ((double)errors_uc)/ ((double)(quant*g_m)), ((double)errors_c) / ((double)(quant*g_n)));
		fclose(fp);
	}


}


void update_erf_file_with_iter(char *ifname, double eb_n0, int quant, int g_n, int g_m, int errors_uc, int errors_c, int iter_quant, int sum_iter, int sum_iter_sqr){
	FILE *fp;
	erf_t erf_array[MAX_ERFS], erf_array_new[MAX_ERFS];
	int pos;
	int cnt, found_flag, nr_elements;
	double tmp1, tmp2;
	// double average_iter, average_iter_sqr;
	char str[200], fname[200], *p;
	
	strcpy(fname, ifname);
	p = strstr(fname, ".dat");
	ecc_assert(p, "Fname %s doesn't contain suffix .dat\n", ifname);
	strcpy(p, ".erf");

	fp = fopen(fname, "rt");
	if (fp) {
		// file exist
		
		// analyze the file (reading its content)
		cnt =0;
		while ((cnt<MAX_ERFS) && (!feof(fp))) {
			
			fscanf(fp, "%lg, %d, %d, %d, %d, %d, %lg, %lg, %lg, %lg\n", 
			       &erf_array[cnt].eb_n0, 
			       &erf_array[cnt].nr_iterations, 
			       &erf_array[cnt].g_n, 
			       &erf_array[cnt].g_m, 
			       &erf_array[cnt].uncoded_errors,
			       &erf_array[cnt].coded_errors,
			       &tmp1, &tmp2, 
			       &erf_array[cnt].average_iter,
			       &erf_array[cnt].average_iter_sqr
			       );
			cnt++;

		}
		nr_elements = cnt;
		
		ecc_assert(feof(fp), "File %s is to big; suspicious\n", fname);
		cnt = 0; 
		found_flag = 0;
		while ((!found_flag) && (cnt<nr_elements)) {
			if (erf_array[cnt].eb_n0 == eb_n0) {
				found_flag = 1;
				//check for consistency
				ecc_assert(erf_array[cnt].g_n == g_n, "Inconsistent n (%d vs. %d)\n", erf_array[cnt].g_n, g_n);
				ecc_assert(erf_array[cnt].g_m == g_m, "Inconsistent n (%d vs. %d)\n", erf_array[cnt].g_m, g_m);
					
				erf_array[cnt].average_iter = (erf_array[cnt].average_iter * ((double)erf_array[cnt].nr_iterations) + ((double)sum_iter))/((double)(erf_array[cnt].nr_iterations + iter_quant));
				erf_array[cnt].average_iter_sqr = (erf_array[cnt].average_iter_sqr *((double) erf_array[cnt].nr_iterations) + ((double)sum_iter_sqr))/((double)(erf_array[cnt].nr_iterations + iter_quant));

				erf_array[cnt].nr_iterations += quant;
				erf_array[cnt].uncoded_errors += errors_uc;
				erf_array[cnt].coded_errors += errors_c;
				
				memcpy(&erf_array_new[0], &erf_array[0], sizeof(erf_t)*nr_elements);
			} else 
				cnt ++;
		}
		if (!found_flag) {
			pos = -1;
			if (eb_n0 < erf_array[0].eb_n0) {
				memcpy(&erf_array_new[1], &erf_array[0], sizeof(erf_t)*nr_elements);
				pos = 0;
			} else if (eb_n0 > erf_array[nr_elements-1].eb_n0) {
				memcpy(&erf_array_new[0], &erf_array[0], sizeof(erf_t)*nr_elements);
				pos = nr_elements;
			} else {
				for (cnt=0; cnt<nr_elements; cnt++) {
					if ((eb_n0 > erf_array[cnt].eb_n0) && (eb_n0 < erf_array[cnt+1].eb_n0)) {
						memcpy(&erf_array_new[0], &erf_array[0], sizeof(erf_t)*(cnt+1));
						memcpy(&erf_array_new[cnt+2], &erf_array[cnt+1], sizeof(erf_t)*(nr_elements - (cnt+1)));
						pos = cnt+1;
					}
				}
			}
			ecc_assert(pos != -1, "Error: Bug in design\n");
			nr_elements++;			
			erf_array_new[pos].eb_n0 = eb_n0;
			erf_array_new[pos].nr_iterations = quant;
			erf_array_new[pos].g_n = g_n;
			erf_array_new[pos].g_m = g_m; 
			erf_array_new[pos].uncoded_errors = errors_uc;
			erf_array_new[pos].coded_errors = errors_c;
			erf_array_new[pos].average_iter = ((double)sum_iter) / ((double)iter_quant);
			erf_array_new[pos].average_iter_sqr = ((double)sum_iter_sqr) / ((double)iter_quant);
		} 
		fclose(fp);
		// saving back up for the file
		sprintf(str, "mv %s %s~", fname, fname);
		system(str);
		fp = fopen(fname, "wt");
		for (cnt=0; cnt<nr_elements; cnt++)
			fprintf(fp, "%lg, %d, %d, %d, %d, %d, %lg, %lg, %lg, %lg\n", 
				erf_array_new[cnt].eb_n0, 
				erf_array_new[cnt].nr_iterations, 
				erf_array_new[cnt].g_n , 
				erf_array_new[cnt].g_m, 
				erf_array_new[cnt].uncoded_errors, 
				erf_array_new[cnt].coded_errors, 
				((double)erf_array_new[cnt].uncoded_errors)/ ((double)(erf_array_new[cnt].nr_iterations*erf_array_new[cnt].g_m)), 
				((double)erf_array_new[cnt].coded_errors) / ((double)(erf_array_new[cnt].nr_iterations*erf_array_new[cnt].g_n)), 
				erf_array_new[cnt].average_iter, 
				erf_array_new[cnt].average_iter_sqr);
		fclose(fp);
		//printf("Check file\n");
		//getchar();
		
	} else { 
		// file doesn't exist
		fp = fopen(fname, "w");
		fprintf(fp, "%lg, %d, %d, %d, %d, %d, %lg, %lg, %lg, %lg\n", 
			eb_n0, quant, g_n, g_m, errors_uc, errors_c, ((double)errors_uc)/ ((double)(quant*g_m)), ((double)errors_c) / ((double)(quant*g_n)), ((double)sum_iter) / ((double)quant), ((double)sum_iter_sqr) / ((double)iter_quant));
		fclose(fp);
	}


}


void update_erf_file_with_iter_and_fer(char *ifname, double eb_n0, int quant, int g_n, int g_m, int errors_uc, int errors_c, int iter_quant, int sum_iter, int sum_iter_sqr, int frame_errors_uc, int frame_errors_c, double sum_time, double sum_tg_time){
	FILE *fp;
	erf_t erf_array[MAX_ERFS], erf_array_new[MAX_ERFS];
	int pos;
	int cnt, found_flag, nr_elements;
	double tmp1, tmp2;
	// double average_iter, average_iter_sqr;
	char str[200], fname[200], *p;
	
	strcpy(fname, ifname);
	p = strstr(fname, ".dat");
	ecc_assert(p, "Fname %s doesn't contain suffix .dat\n", ifname);
	strcpy(p, ".erf");

	// printf("Starting to update the output file, data is: sum_time=%f and sum_tg_time=%f\n", sum_time, sum_tg_time);
	
	fp = fopen(fname, "rt");
	if (fp) {
		// file exist
		
		// analyze the file (reading its content)
		cnt =0;
		while ((cnt<MAX_ERFS) && (!feof(fp))) {
			
			fscanf(fp, "%lg, %d, %d, %d, %d, %d, %lg, %lg, %lg, %lg, %d, %d, %lg, %lg\n",
			       &erf_array[cnt].eb_n0, 
			       &erf_array[cnt].nr_iterations, // Number of frames transmitted before
			       &erf_array[cnt].g_n, 
			       &erf_array[cnt].g_m, 
			       &erf_array[cnt].uncoded_errors,
			       &erf_array[cnt].coded_errors,
			       &tmp1, // BER of uncoded transmission calculated so far
			       &tmp2, // BER of coded transmission calculated so far
			       &erf_array[cnt].average_iter, // BP decoding iterations
			       &erf_array[cnt].average_iter_sqr,
			       &erf_array[cnt].uncoded_frame_erros,
				   &erf_array[cnt].coded_frame_errors,
				   &erf_array[cnt].average_time,
				   &erf_array[cnt].average_tg_time
			       );
			cnt++;

		}
		nr_elements = cnt;
		
		ecc_assert(feof(fp), "File %s is to big; suspicious\n", fname);
		cnt = 0; 
		found_flag = 0;
		while ((!found_flag) && (cnt<nr_elements)) { // Searching for the SNR data from the file that matches the SNR we currently received data about
			if (erf_array[cnt].eb_n0 == eb_n0) {
				found_flag = 1;
				//check for consistency
				ecc_assert(erf_array[cnt].g_n == g_n, "Inconsistent n (%d vs. %d)\n", erf_array[cnt].g_n, g_n);
				ecc_assert(erf_array[cnt].g_m == g_m, "Inconsistent n (%d vs. %d)\n", erf_array[cnt].g_m, g_m);
					
				erf_array[cnt].average_time = (erf_array[cnt].average_time * ((double) erf_array[cnt].nr_iterations) + ((double) sum_time)) / ((double) (erf_array[cnt].nr_iterations + iter_quant));
				erf_array[cnt].average_tg_time = (erf_array[cnt].average_tg_time * erf_array[cnt].average_iter * ((double) erf_array[cnt].nr_iterations) + ((double) sum_tg_time)) / (((double) sum_iter) + erf_array[cnt].average_iter * ((double) erf_array[cnt].nr_iterations));
				erf_array[cnt].average_iter = (erf_array[cnt].average_iter * ((double)erf_array[cnt].nr_iterations) + ((double)sum_iter))/((double)(erf_array[cnt].nr_iterations + iter_quant));
				erf_array[cnt].average_iter_sqr = (erf_array[cnt].average_iter_sqr *((double) erf_array[cnt].nr_iterations) + ((double)sum_iter_sqr))/((double)(erf_array[cnt].nr_iterations + iter_quant));
				erf_array[cnt].nr_iterations += quant;
				erf_array[cnt].uncoded_errors += errors_uc;
				erf_array[cnt].coded_errors += errors_c;
				erf_array[cnt].uncoded_frame_erros += frame_errors_uc;
				erf_array[cnt].coded_frame_errors += frame_errors_c;
				
				memcpy(&erf_array_new[0], &erf_array[0], sizeof(erf_t)*nr_elements);
			} else 
				cnt ++;
		}
		if (!found_flag) {
			pos = -1;
			if (eb_n0 < erf_array[0].eb_n0) {
				memcpy(&erf_array_new[1], &erf_array[0], sizeof(erf_t)*nr_elements);
				pos = 0;
			} else if (eb_n0 > erf_array[nr_elements-1].eb_n0) {
				memcpy(&erf_array_new[0], &erf_array[0], sizeof(erf_t)*nr_elements);
				pos = nr_elements;
			} else {
				for (cnt=0; cnt<nr_elements; cnt++) {
					if ((eb_n0 > erf_array[cnt].eb_n0) && (eb_n0 < erf_array[cnt+1].eb_n0)) {
						memcpy(&erf_array_new[0], &erf_array[0], sizeof(erf_t)*(cnt+1));
						memcpy(&erf_array_new[cnt+2], &erf_array[cnt+1], sizeof(erf_t)*(nr_elements - (cnt+1)));
						pos = cnt+1;
					}
				}
			}
			ecc_assert(pos != -1, "Error: Bug in design\n");
			nr_elements++;			
			erf_array_new[pos].eb_n0 = eb_n0;
			erf_array_new[pos].nr_iterations = quant;
			erf_array_new[pos].g_n = g_n;
			erf_array_new[pos].g_m = g_m; 
			erf_array_new[pos].uncoded_errors = errors_uc;
			erf_array_new[pos].coded_errors = errors_c;
			erf_array_new[cnt].uncoded_frame_erros = frame_errors_uc;
			erf_array_new[cnt].coded_frame_errors = frame_errors_c;
			erf_array_new[pos].average_iter = ((double)sum_iter) / ((double)iter_quant);
			erf_array_new[pos].average_time = ((double) sum_time) / ((double) iter_quant);
			erf_array_new[pos].average_tg_time = ((double) sum_tg_time) / ((double) sum_iter);
			erf_array_new[pos].average_iter_sqr = ((double)sum_iter_sqr) / ((double)iter_quant);
		} 
		fclose(fp);
		// saving back up for the file
		sprintf(str, "mv %s %s~", fname, fname);
		system(str);
		fp = fopen(fname, "wt");
		for (cnt=0; cnt<nr_elements; cnt++)
			fprintf(fp, "%lg, %d, %d, %d, %d, %d, %lg, %lg, %lg, %lg, %d, %d, %lg, %lg\n",
				erf_array_new[cnt].eb_n0, 
				erf_array_new[cnt].nr_iterations, 
				erf_array_new[cnt].g_n , 
				erf_array_new[cnt].g_m, 
				erf_array_new[cnt].uncoded_errors, 
				erf_array_new[cnt].coded_errors, 
				((double)erf_array_new[cnt].uncoded_errors)/ ((double)(erf_array_new[cnt].nr_iterations*erf_array_new[cnt].g_m)), 
				((double)erf_array_new[cnt].coded_errors) / ((double)(erf_array_new[cnt].nr_iterations*erf_array_new[cnt].g_n)), 
				erf_array_new[cnt].average_iter, 
				erf_array_new[cnt].average_iter_sqr, 
				erf_array_new[cnt].uncoded_frame_erros, 
				erf_array_new[cnt].coded_frame_errors,
				erf_array_new[cnt].average_time,
				erf_array_new[cnt].average_tg_time);
		fclose(fp);
		//printf("Check file\n");
		//getchar();
		
	} else { 
		// file doesn't exist
		fp = fopen(fname, "w");
		fprintf(fp, "%lg, %d, %d, %d, %d, %d, %lg, %lg, %lg, %lg, %d, %d, %lg, %lg\n",
			eb_n0, quant, g_n, g_m, errors_uc, errors_c, ((double) errors_uc) / ((double) (quant*g_m)), ((double) errors_c) / ((double) (quant*g_n)), ((double) sum_iter) / ((double) quant), ((double) sum_iter_sqr) / ((double) iter_quant), frame_errors_uc, frame_errors_c, ((double) sum_time) / ((double) quant), ((double) sum_tg_time) / ((double) sum_iter));
		fclose(fp);
	}


}

// The following 2 functions generate standard normally - distributed random numbers (with expectation 0 and variance 1).
// randn() generates 1 random number and randn2() generates 2 independetly random number, both normally distributed.
// The random numbers are generated using the Box - Muller algorithm for generating normally distributed random numbers from 2 independent, uniformly distributed numbers.
static inline double randn()
{
	double x1, x2, w;
	double y1;
	do {
		x1 = 2.0 * ((double)rand())/(RAND_MAX+1.0) - 1.0;
		x2 = 2.0 * ((double)rand())/(RAND_MAX+1.0) - 1.0;


		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );
  
	w = sqrt( (-2.0 * log( w ) ) / w );
	y1 = x1 * w;
	return y1;
}


static inline void randn2(double *y1, double *y2)
{
	double x1, x2, w;
  
	do {
		x1 = 2.0 * ((double)rand())/(RAND_MAX+1.0) - 1.0;
		x2 = 2.0 * ((double)rand())/(RAND_MAX+1.0) - 1.0;


		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );
  
	w = sqrt( (-2.0 * log( w ) ) / w );
	*y1 = x1 * w;
	*y2 = x2 * w;
}

int_mtx_t* gaussian_elimination(int_mtx_t *q, int specify_cols, int_mtx_t *cols){
	
	int_mtx_t *p_mod; 
	int_mtx_t *internal_cols;
	int i, j, k, found_flag, data;
	
	p_mod = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx_from_imtx(p_mod, q);

	
	if (specify_cols) {
		ecc_assert(cols->m == 1, "ERROR: cols M dimenison must be 1 (cols->m = %d)\n", cols->m);
		ecc_assert(cols->n <= q->m, "ERROR: cols N dimension cannot be larger than p->n (%d vs. %d)\n", cols->n, q->m);
		internal_cols = cols;
	} else {
		internal_cols = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
		init_int_mtx(internal_cols, 1, q->m);
		for (i=0; i<internal_cols->n; i++)
			put_int_element(internal_cols, 0, i, i);
			
	}
	
	for (i=0; i< internal_cols->n; i++) {
		// pivoting
		

		// if the pivot element is zero	
		if (0 == get_int_element(p_mod, i, get_int_element(internal_cols, 0, i))) {

			// find the first non-zero element
			found_flag = 0;
			j=i+1;
			while ((!found_flag) && (j<p_mod->m)) {
				if (get_int_element(p_mod,j, get_int_element(internal_cols, 0, i))) {
					found_flag = 1;
					for (k=0; k<p_mod->n; k++) {
						data = get_int_element(p_mod, i, k) + get_int_element(p_mod, j, k);
						put_int_element(p_mod, i, k, data%2);
					}
				} else 
					j++;
			}
			ecc_assert(found_flag == 1, "Unable to locate pivot supplement element\n");	
		}


		// now we have the pivot element on the position
		for (j=0; j<p_mod->m; j++) {
			if ((j!=i) && (get_int_element(p_mod, j, get_int_element(internal_cols, 0, i)))) {
				for (k=0; k<p_mod->n; k++) {
					data = get_int_element(p_mod, i, k) + get_int_element(p_mod, j, k);
					put_int_element(p_mod, j, k, data%2);
					}
			}
		}
	}

	free_int_mtx(internal_cols);
	return p_mod;
}


int_mtx_t* gaussian_elimination_new(int_mtx_t *q, int specify_cols, int_mtx_t *cols){
	
	int_mtx_t *p_mod; 
	int_mtx_t *internal_cols;
	int i, j, k, found_flag, data, pivot_line;
	
	p_mod = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx_from_imtx(p_mod, q);

	
	if (specify_cols) {
		ecc_assert(cols->m == 1, "ERROR: cols M dimenison must be 1 (cols->m = %d)\n", cols->m);
		ecc_assert(cols->n <= q->m, "ERROR: cols N dimension cannot be larger than p->n (%d vs. %d)\n", cols->n, q->m);
		internal_cols = cols;
	} else {
		internal_cols = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
		init_int_mtx(internal_cols, 1, q->m);
		for (i=0; i<internal_cols->n; i++)
			put_int_element(internal_cols, 0, i, i);
			
	}
	
	pivot_line = 0;
	for (i=0; i< internal_cols->n; i++) {
		// pivoting
		

		
		
		// if the pivot element is zero
		found_flag = 1;
		if (0 == get_int_element(p_mod, pivot_line, get_int_element(internal_cols, 0, i))) {

			// find the first non-zero element
			found_flag = 0;
			j=pivot_line+1;
			while ((!found_flag) && (j<p_mod->m)) {
				if (get_int_element(p_mod,j, get_int_element(internal_cols, 0, i))) {
					found_flag = 1;
					for (k=0; k<p_mod->n; k++) {
						data = get_int_element(p_mod, pivot_line, k) + get_int_element(p_mod, j, k);
						put_int_element(p_mod, pivot_line, k, data%2);
					}
				} else 
					j++;
			}
			//ecc_assert(found_flag == 1, "Unable to locate pivot supplement element\n");	
		}


		// now we have the pivot element on the position
		if (found_flag) {
			for (j=0; j<p_mod->m; j++) {
				if ((j!=pivot_line) && (get_int_element(p_mod, j, get_int_element(internal_cols, 0, i)))) {
					for (k=0; k<p_mod->n; k++) {
						data = get_int_element(p_mod, pivot_line, k) + get_int_element(p_mod, j, k);
						put_int_element(p_mod, j, k, data%2);
					}
				}
			}
			pivot_line ++;
		}
//		printf("i = %d, col = %d, found_flag = %d, pivot_line = %d\n", i, get_int_element(internal_cols, 0, i), found_flag, pivot_line);
//		print_int_mtx(p_mod);
	}

	//free_int_mtx(internal_cols);
	return p_mod;
}


// This function generates a random bits vector (of GF(2) symbols)
void init_symbol_data(gf2_vec_t *symb)
{
	int i;
	
	for (i=0; i<symb->n; i++) // For every symbol
		// Randomize a value, and divide it with RAND_MAX (the result will always be between 0 and 1).
		// If the result is less than 0.5, the randomized bit (symbol) is 0. If it is 0.5 <= value <= 1, then the randomized bit (symbol) is 1.
		symb->data[i] = (int)round(((double)rand())/((double)RAND_MAX));
}

// This function encodes given information bits to a codewords using generator matrix G
void encode_symbols(gf2_vec_t *info, int_mtx_t *G, gf2_vec_t *symb)
{
	int i, j, sum;
 
 	// Checking all dimensions agree
	ecc_assert(G->n == symb->n, "symbol length doesn't match generator matrix (%d vs. %d)\n", G->n, symb->n);
	ecc_assert(G->m == info->n, "info   length doesn't match generator matrix (%d vs. %d)\n", G->m, info->n);

	// For every code bit
	for (i=0; i<G->n; i++) 	{
			sum = 0;
			// Standard encoding using matrix multiplication
			for (j=0; j<G->m; j++)
				sum += info->data[j] * get_int_element(G, j, i);
			symb->data[i] = sum %2; // Result is over GF(2)
		}    
}


// This function adds Gaussian noise to a vector of BPSK modulated bits
// The function received the codeword (symb), the desired SNR in dB (EBN0_dB) and the code's rate (R)
// The function returns a vector of received samples after transmission over AWGN channel (data)
void init_measured_data(gf2_vec_t *symb, real_vec_t *data, double EBN0_dB, double R)
{
	int i;
	double ESN0;

	// Convert the SNR from dB to Linear scale
	ESN0 = R*pow(10, EBN0_dB/10);
 
 	// For each bit (symbol), modulate it to BPSK where '0' --> (-1) and '1' --> (1)
 	// Then add the Gaussian noise, with variance of 1/(2*R*EBN0). This means that EBN0 is the SNR of Information bits (and not coded bits).
	for (i=0; i<data->n; i++)
		data->data[i] = (2.0*(double)symb->data[i] - 1) + (1.0/sqrt(2*ESN0))*randn();
}

// This function performs BP decoding of given AWGN channel output LLRs, which compose a codeword belongs to a code represented by H.
// The channel LLRs are given at lambda. This decoder performs n_iter BP iterations per codeword and uses alpha as damping coefficient.
// The decoded symbols are given at dec_symb and mu stores the decoded symbols' LLRs (output LLRs).
static void tanner_graph_decoder(int_mtx_t *H, real_vec_t *lambda, gf2_vec_t *dec_symb, real_vec_t *mu, int n_iter, double alpha, double *tg_time)
{
	int i, check_cnt, var_cnt, j, nd_cnt, z, w; 
	real_vec_t *sum_E_s;
	double_mtx_t *s_E, *s_E_prev, *E_s, *gamma;
	double sum, prod, tmp;
	Timer t;
	// Initialize all variables
	s_E		= (double_mtx_t*)ecc_malloc(sizeof(double_mtx_t));
	s_E_prev	= (double_mtx_t*)ecc_malloc(sizeof(double_mtx_t));
	E_s		= (double_mtx_t*)ecc_malloc(sizeof(double_mtx_t));
	gamma		= (double_mtx_t*)ecc_malloc(sizeof(double_mtx_t));
	sum_E_s	= (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
  
	init_double_mtx(s_E,      H->m, H->n);
	init_double_mtx(s_E_prev, H->m, H->n);
	init_double_mtx(E_s,      H->m, H->n);
	init_double_mtx(gamma,    H->m, H->n);
	init_real_vec(sum_E_s, H->n);

	t.start();
 	// Initialize matrix of LLRs with the channel LRRs values
 	// We assign LLR value at position (i,j) only if H(i,j) == 1. Otherwise we assign 0
	for (i=0; i<H->m; i++)
		for (j=0; j<H->n; j++)
			put_double_element(gamma, i, j, lambda->data[j]*get_int_element(H, i, j));
  
	//  print_int_mtx(H);
	//  printf("\n");
	 // print_double_mtx(gamma);
	 // print_double_mtx(E_s);
	 // print_double_mtx(s_E_prev);
	//  getchar();

	// Iterating as the number of iterations given
	// s_E - messages from variable nodes to check nodes
	// E_s - messages from check nodes to variable nodes
	for (i=0; i<n_iter; i++) {
		//printf("Calculating Variable nodes to check nodes\n");
	
		// passing message from variable node to check node
		for (check_cnt = 0; check_cnt < H->m; check_cnt ++) { // For every check node
			//printf("-=-=-=-=- check_cnt = %d\n", check_cnt);
			for (var_cnt = 0; var_cnt < H->n; var_cnt ++) { // For every variable node
				sum = 0;
				// if variable var_cnt connected to check node check_cnt
				if (get_int_element(H, check_cnt, var_cnt)) {		  
					// Than we calculate the message to be sent from variable node var_cnt to check_cnt
					// The message is simply a sum of all message arrived to the variable node from all other check nodes
					// Note: the code assumes E_s is initialized to 0. This may not always be true.
					for (j=0; j<H->m; j++)
						if (j!=check_cnt)  {
							sum += get_double_element(E_s, j, var_cnt);
							//printf("summing var_cnt = %d, check_cnt = %d, sum = %g\n", var_cnt, check_cnt, sum);
						}
					sum += get_double_element(gamma, check_cnt, var_cnt);
					//printf("summing the LLR = %g\n", sum);		  
				}
				//printf("putting element [%d %d] = %g\n", check_cnt, var_cnt, sum);
				//getchar();
				put_double_element(s_E, check_cnt, var_cnt, sum);	      
			}
		}
		
		// Damping extrinsic LLRs (note there is some difference with the common case - the intrinsic LLR is multiplied by (1-alpha) and isn't taken as is)
		// Note that s_E contains messages from variable nodes to check nodes
		// According to Ilan, this is done to avoid numeric instability.
		// handling alpha
		// s_E_mtx = alpha.*s_E_mtx  + (1-alpha).*s_E_mtx_prev;
		if (alpha < 1.0) {
			for (w=0; w<s_E->m; w++)
				for (z=0; z<s_E->n; z++) {
					if ((finite(get_double_element(s_E, w, z))) && (finite(get_double_element(s_E_prev, w, z))))
						put_double_element(s_E, w, z, alpha*get_double_element(s_E, w, z) + (1 - alpha)*get_double_element(s_E_prev, w, z));
					else
						ecc_assert(0, "Overflow issue\n");
					put_double_element(s_E_prev, w, z, get_double_element(s_E, w, z));
				}
		}
      

		//printf("printing s_E matrix\n");
		//print_double_mtx(s_E);
		//printf("\n");
		//getchar();


		//      printf("Calculating check nodes to Variable nodes\n");


		for (var_cnt = 0; var_cnt < H->n; var_cnt ++) {
			//printf("-=-=-=-=-=- var_cnt = %d\n", var_cnt);
			for (check_cnt = 0; check_cnt < H->m; check_cnt++) {
				prod = 1.0; nd_cnt = 0;
				if (get_int_element(H, check_cnt, var_cnt)) {
					// Than we calculate the message to be sent from check node check_cnt to variable node var_cnt
					// The message is a product of messages arrived to the check node from all other variable nodes connected to it, after applying tanh() on these messages
					for (j=0; j<H->n; j++) {
						if ((0 != get_double_element(s_E, check_cnt, j)) && (j!=var_cnt)) {
							tmp = get_double_element(s_E, check_cnt, j);
							//prod *= tanh(get_double_element(s_E, check_cnt, j)/2);
							prod *= tanh(clamp(get_double_element(s_E, check_cnt, j)/2, -15.0, 15.0)); // Note that tanh(0)=0 hence no need to treat zero believes in a different manner
							nd_cnt ++;
							//printf("multiplying var_cnt = %d, check_cnt = %d, prod = %g, tmp=%g\n", var_cnt, check_cnt, prod, tmp);
						}
					}
					if (0 == nd_cnt) {
						// If var_cnt is the only variable that is conencted to the check node, the message that will be sent is 0
						put_double_element(E_s, check_cnt, var_cnt, 0);	 
						//printf("putting element [%d %d] zero, nd_cnt = 0\n", check_cnt, var_cnt);
					}
					else if (1 == nd_cnt) {
						// If there is only 1 variable connected to check_cnt node besides var_cnt, the message that will be sent is the message received from this other variable node
						put_double_element(E_s, check_cnt, var_cnt, tmp);
						//printf("putting element [%d %d] %g, nd_cnt = 1\n", check_cnt, var_cnt, tmp);
					} else {
						// If there are multiple variables connected to check_cnt node besides var_cnt, the message that will be sent will be calculated in the standard way
						put_double_element(E_s, check_cnt, var_cnt, 2*atanh(prod));	 // no clamp	  
		      
						//printf("putting element [%d %d] %g, nd_cnt = %d\n", check_cnt, var_cnt, 2*atanh(prod), nd_cnt);
					}
				} else {
					put_double_element(E_s, check_cnt, var_cnt, 0);	      
					//printf("putting element [%d %d] zero\n", check_cnt, var_cnt);
				}
			}
		}

		//printf("printing E_s matrix\n");
		//print_double_mtx(E_s);
		//printf("\n");
		//getchar();


	}
  
	//  printf("printing E_s matrix\n");
	//  print_double_mtx(E_s);
	//  printf("\n");

	// Summing all messages received from the check nodes, per variable
	sum_double_mtx(E_s, sum_E_s);

	//  printf("printing sum_E_s vec\n");
	//  print_rvec(sum_E_s);
	//  printf("\n");
	//
	// Adding the channel LLRs to the previous summary
	// This produces the bit estimates for each bit (summary of messages from all nodes and channel LLRs)
	sum_rvec(sum_E_s, lambda, mu);

	//  printf("printing lambda vec\n");
	//  print_rvec(lambda);
	//  printf("\n");
	//
	//  printf("printing mu vec\n");
	//  print_rvec(mu);
	//  printf("\n");
	//
	for (i=0; i<H->n; i++)
		dec_symb->data[i] = ((mu->data[i] >= 0) ? 0 : 1);  


	t.stop();
	*tg_time = t.getElapsedTimeInMicroSec();
	free_double_mtx(s_E);
	free_double_mtx(s_E_prev);
	free_double_mtx(E_s);
	free_double_mtx(gamma);
	free_real_vec(sum_E_s);
	
}

// More efficient implemented version of BP from Illia
static void reduction_2_Eff_TG_Decoder_a(Decoder *dec, real_vec_t *lambda, gf2_vec_t *dec_symb, real_vec_t *mu, int n_iter, double alpha, double *tg_time)
{
	// create input to decoder => llr (Channel observations)
	for (int i = 0; i < dec->mBlocklength; i++)
	{
		dec->_LogLikelihoodRatio[i] = lambda->data[i];
	}

	// run decoder
	dec->Eff_TG_Decoder_a();

	// update output of decoder
	for (int i = 0; i < dec->mBlocklength; i++)
	{
		mu->data[i] = dec->OutputFromDecoder[i];
		dec_symb->data[i] = ((mu->data[i] >= 0) ? 0 : 1);
	}

	// update decoding time variable
	*tg_time = dec->mTGDec_a_DecodeTime;

}

// This function performs weighted (soft) BP decoding of given AWGN channel output LLRs, which compose a codeword belongs to a code represented by H.
// The channel LLRs are given at lambda. This decoder performs n_iter BP iterations per codeword and uses alpha as damping coefficient.
// The decoded symbols are given at dec_symb and mu stores the decoded symbols' LLRs (output LLRs).
static void soft_tanner_graph_decoder(int_mtx_t *H, real_vec_t *lambda, gf2_vec_t *dec_symb, real_vec_t *mu, int n_iter, double alpha, double *tg_time, double_mtx_t *W_check2var, double_mtx_t *W_output, int_mtx_t *accumulated_rows_sum_matrix, int_mtx_t *accumulated_columns_sum_matrix)
{
	int i, check_cnt, var_cnt, j, nd_cnt, z, w; 
	real_vec_t *sum_E_s;
	double_mtx_t *s_E, *s_E_prev, *E_s, *gamma;
	double sum, prod, tmp;
	Timer t;
	int k, l, weights_column_index, weight_index; // additions for using weights
	// int weights_column_index1, weight_index1;
	double weight, weighted_llr_from_check_nodes; // additions for using weights
	real_vec_t *weighted_sum_E_s;
	// Initialize all variables
	s_E		= (double_mtx_t*)ecc_malloc(sizeof(double_mtx_t));
	s_E_prev	= (double_mtx_t*)ecc_malloc(sizeof(double_mtx_t));
	E_s		= (double_mtx_t*)ecc_malloc(sizeof(double_mtx_t));
	gamma		= (double_mtx_t*)ecc_malloc(sizeof(double_mtx_t));
	sum_E_s	= (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	weighted_sum_E_s	= (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
  
	init_double_mtx(s_E,      H->m, H->n);
	init_double_mtx(s_E_prev, H->m, H->n);
	init_double_mtx(E_s,      H->m, H->n);
	init_double_mtx(gamma,    H->m, H->n);
	init_real_vec(sum_E_s, H->n);
	init_real_vec(weighted_sum_E_s, H->n);

	t.start();
 	// Initialize matrix of LLRs with the channel LRRs values
 	// We assign LLR value at position (i,j) only if H(i,j) == 1. Otherwise we assign 0
	for (i=0; i<H->m; i++)
		for (j=0; j<H->n; j++)
			put_double_element(gamma, i, j, lambda->data[j]*get_int_element(H, i, j));
  
	//  print_int_mtx(H);
	//  printf("\n");
	 // print_double_mtx(gamma);
	 // print_double_mtx(E_s);
	 // print_double_mtx(s_E_prev);
	//  getchar();

	// Iterating as the number of iterations given
	// s_E - messages from variable nodes to check nodes
	// E_s - messages from check nodes to variable nodes
	for (i=0; i<n_iter; i++) {
		//printf("Calculating Variable nodes to check nodes\n");
	
		// passing message from variable node to check node
		for (check_cnt = 0; check_cnt < H->m; check_cnt ++) { // For every check node
			//printf("-=-=-=-=- check_cnt = %d\n", check_cnt);
			for (var_cnt = 0; var_cnt < H->n; var_cnt ++) { // For every variable node
				sum = 0;
				// if variable var_cnt connected to check node check_cnt
				if (get_int_element(H, check_cnt, var_cnt)) {		  
					// Than we calculate the message to be sent from variable node var_cnt to check_cnt
					// The message is a weighted sum of all message arrived to the variable node from all other check nodes
					// Note: the code assumes E_s is initialized to 0. This may not always be true.
					
					// Calcuation of column index in W_check2var matrix
					weights_column_index = get_int_element(accumulated_columns_sum_matrix, check_cnt, var_cnt);
					// weights_column_index = 0;
					// for (k = 0; k < var_cnt; k++) {
					// 	for (l = 0; l < H->m; l++) {
					// 		weights_column_index += get_int_element(H, l, k);
					// 	}
					// }
					// for (l = 0; l < check_cnt; l++) {
					// 	weights_column_index += get_int_element(H, l, var_cnt);
					// }
					// ecc_assert((weights_column_index == weights_column_index1), "Weight column indices are different for check_cnt=%d and var_cnt=%d: %d vs %d\n", check_cnt, var_cnt, weights_column_index, weights_column_index1);
					for (j=0; j<H->m; j++) {
						if (get_int_element(H, j, var_cnt) && (j != check_cnt)) {
							// Then we have a message from check node j which we should sum over
							// First we calculate the index of this message's weight in the weight column
							weight_index = get_int_element(accumulated_rows_sum_matrix, j, var_cnt);
							// weight_index = 0;
							// for (l = 0; l < j; l++) {
							// 	for (k = 0; k < H->n; k++) {
							// 		weight_index += get_int_element(H, l, k);
							// 	}
							// }
							// for (k = 0; k < var_cnt; k++) {
							// 	weight_index += get_int_element(H, j, k);
							// }
							// ecc_assert((weight_index == weight_index1), "Weight indices are different for check_cnt=%d and var_cnt=%d: %d vs %d\n", j, var_cnt, weight_index, weight_index1);
							// Now get the message weight
							weight = get_double_element(W_check2var, weight_index, weights_column_index);
							sum += weight*get_double_element(E_s, j, var_cnt);
							//printf("summing var_cnt = %d, check_cnt = %d, sum = %g\n", var_cnt, check_cnt, sum);
						}
					}
					sum += get_double_element(gamma, check_cnt, var_cnt);
					//printf("summing the LLR = %g\n", sum);		  
				}
				//printf("putting element [%d %d] = %g\n", check_cnt, var_cnt, sum);
				//getchar();
				put_double_element(s_E, check_cnt, var_cnt, sum);	      
			}
		}
		
		// Damping extrinsic LLRs (note there is some difference with the common case - the intrinsic LLR is multiplied by (1-alpha) and isn't taken as is)
		// Note that s_E contains messages from variable nodes to check nodes
		// According to Ilan, this is done to avoid numeric instability.
		// handling alpha
		// s_E_mtx = alpha.*s_E_mtx  + (1-alpha).*s_E_mtx_prev;
		if (alpha < 1.0) {
			for (w=0; w<s_E->m; w++)
				for (z=0; z<s_E->n; z++) {
					if ((finite(get_double_element(s_E, w, z))) && (finite(get_double_element(s_E_prev, w, z))))
						put_double_element(s_E, w, z, alpha*get_double_element(s_E, w, z) + (1 - alpha)*get_double_element(s_E_prev, w, z));
					else
						ecc_assert(0, "Overflow issue\n");
					put_double_element(s_E_prev, w, z, get_double_element(s_E, w, z));
				}
		}
      

		//printf("printing s_E matrix\n");
		//print_double_mtx(s_E);
		//printf("\n");
		//getchar();


		//      printf("Calculating check nodes to Variable nodes\n");


		for (var_cnt = 0; var_cnt < H->n; var_cnt ++) {
			//printf("-=-=-=-=-=- var_cnt = %d\n", var_cnt);
			for (check_cnt = 0; check_cnt < H->m; check_cnt++) {
				prod = 1.0; nd_cnt = 0;
				if (get_int_element(H, check_cnt, var_cnt)) {
					// Than we calculate the message to be sent from check node check_cnt to variable node var_cnt
					// The message is a product of messages arrived to the check node from all other variable nodes connected to it, after applying tanh() on these messages
					for (j=0; j<H->n; j++) {
						if ((0 != get_double_element(s_E, check_cnt, j)) && (j!=var_cnt)) {
							tmp = get_double_element(s_E, check_cnt, j);
							//prod *= tanh(get_double_element(s_E, check_cnt, j)/2);
							prod *= tanh(clamp(get_double_element(s_E, check_cnt, j)/2, -15.0, 15.0)); // Note that tanh(0)=0 hence no need to treat zero believes in a different manner
							nd_cnt ++;
							//printf("multiplying var_cnt = %d, check_cnt = %d, prod = %g, tmp=%g\n", var_cnt, check_cnt, prod, tmp);
						}
					}
					if (0 == nd_cnt) {
						// If var_cnt is the only variable that is conencted to the check node, the message that will be sent is 0
						put_double_element(E_s, check_cnt, var_cnt, 0);	 
						//printf("putting element [%d %d] zero, nd_cnt = 0\n", check_cnt, var_cnt);
					}
					else if (1 == nd_cnt) {
						// If there is only 1 variable connected to check_cnt node besides var_cnt, the message that will be sent is the message received from this other variable node
						put_double_element(E_s, check_cnt, var_cnt, tmp);
						//printf("putting element [%d %d] %g, nd_cnt = 1\n", check_cnt, var_cnt, tmp);
					} else {
						// If there are multiple variables connected to check_cnt node besides var_cnt, the message that will be sent will be calculated in the standard way
						put_double_element(E_s, check_cnt, var_cnt, 2*atanh(prod));	 // no clamp	  
		      
						//printf("putting element [%d %d] %g, nd_cnt = %d\n", check_cnt, var_cnt, 2*atanh(prod), nd_cnt);
					}
				} else {
					put_double_element(E_s, check_cnt, var_cnt, 0);	      
					//printf("putting element [%d %d] zero\n", check_cnt, var_cnt);
				}
			}
		}

		//printf("printing E_s matrix\n");
		//print_double_mtx(E_s);
		//printf("\n");
		//getchar();


	}
  
	//  printf("printing E_s matrix\n");
	//  print_double_mtx(E_s);
	//  printf("\n");

	// Summing all messages received from the check nodes, per variable
	// sum_double_mtx(E_s, sum_E_s);
	// Now this will be a weighted sum
	for (var_cnt = 0; var_cnt < H->n; ++var_cnt) { // For every variable
		weighted_llr_from_check_nodes = 0;
		for (check_cnt = 0; check_cnt < H->m; ++check_cnt) {
			if (get_int_element(H, check_cnt, var_cnt)) {
				// Then we should perform summary on the message from check node check_cnt
				// First we will find its row index in the weights matrix (the column index is var_cnt)
				weight_index = get_int_element(accumulated_rows_sum_matrix, check_cnt, var_cnt);
				// weight_index = 0;
				// for (l = 0; l < check_cnt; l++) {
				// 	for (k = 0; k < H->n; k++) {
				// 		weight_index += get_int_element(H, l, k);
				// 	}
				// }
				// for (k = 0; k < var_cnt; k++) {
				// 	weight_index += get_int_element(H, check_cnt, k);
				// }
				// ecc_assert((weight_index == weight_index1), "Weight indices are different for output weights for check_cnt=%d and var_cnt=%d: %d vs %d\n", check_cnt, var_cnt, weight_index, weight_index1);
				// Second get the weight
				weight = get_double_element(W_output, weight_index, var_cnt);
				// Third add to the weighted sum
				weighted_llr_from_check_nodes += weight*get_double_element(E_s, check_cnt, var_cnt);
			}
		}
		sum_E_s->data[var_cnt] = weighted_llr_from_check_nodes;
	}

	//  printf("printing sum_E_s vec\n");
	//  print_rvec(sum_E_s);
	//  printf("\n");
	//
	// Adding the channel LLRs to the previous summary
	// This produces the bit estimates for each bit (summary of messages from all nodes and channel LLRs)
	sum_rvec(sum_E_s, lambda, mu);

	//  printf("printing lambda vec\n");
	//  print_rvec(lambda);
	//  printf("\n");
	//
	//  printf("printing mu vec\n");
	//  print_rvec(mu);
	//  printf("\n");
	//
	for (i=0; i<H->n; i++)
		dec_symb->data[i] = ((mu->data[i] >= 0) ? 0 : 1);  


	t.stop();
	*tg_time = t.getElapsedTimeInMicroSec();
	free_double_mtx(s_E);
	free_double_mtx(s_E_prev);
	free_double_mtx(E_s);
	free_double_mtx(gamma);
	free_real_vec(sum_E_s);
	
}

int calc_num_edges_in_H(ecc_t* ecc_p) {
	int_mtx_t *H = ecc_p->H;
	int num_edges = 0;

	for (int i = 0; i < H->m; ++i) {
		for (int j = 0; j < H->n; ++j) {
			num_edges += get_int_element(H, i, j);
		}
	}

	return num_edges;
}

void create_weight_matrices_from_files(ecc_t* ecc_p, char* hidden_layers_weights_filename, char* out_weights_filename) {
	std::ifstream hidden_weights_file(hidden_layers_weights_filename, std::ios::in);
	std::ifstream out_weights_file(out_weights_filename, std::ios::in);
    std::ios_base::iostate file_state;
    int weights_hidden_num_of_rows = ecc_p->num_edges_in_H;
    int weights_hidden_num_of_columns = ecc_p->num_edges_in_H;
    int out_weights_num_of_rows = ecc_p->num_edges_in_H;
    int out_weights_num_of_columns = ecc_p->H->n;

    // Check to see that the files were opened correctly
    if ((!hidden_weights_file.is_open()) | (!out_weights_file.is_open())) {
        std::cerr << "There was a problem opening the weights files! Weights are now unknown.\n";
        return; // exit
    }

    double weight = 0.0;
    // Store values from the text file into weights array
    for (int i = 0; i < weights_hidden_num_of_rows; ++i) {
        for (int j = 0; j < weights_hidden_num_of_columns; ++j) {
            if (hidden_weights_file >> weight) {
                put_double_element(ecc_p->hidden_layers_weights, i, j, weight);
            } else {
                file_state = hidden_weights_file.rdstate();
                if (file_state & std::ifstream::eofbit) {
                    std::cerr << "Hidden weights file error: reached EOF before expected, on row " << i << " column " << j << "\n";
                } else if (file_state & std::ifstream::failbit) {
                    std::cerr << "Hidden weights file error on extraction or imterpretation of characters from the file, on row " << i << " column" << j << "\n";
                } else if (file_state & std::ifstream::badbit) {
                    std::cerr << "Hidden weights file error on stream while attempting to read row " << i << " column " << j << "\n";
                } else {
                    std::cerr << "Unknown hidden weights file error error occured on row " << i << " column " << j << "\n";
                    std::cerr << file_state << std::endl;
                }
            }
        }
    }

    for (int i = 0; i < out_weights_num_of_rows; ++i) {
        for (int j = 0; j < out_weights_num_of_columns; ++j) {
            if (out_weights_file >> weight) {
                put_double_element(ecc_p->out_weights, i, j, weight);
            } else {
                file_state = out_weights_file.rdstate();
                if (file_state & std::ifstream::eofbit) {
                    std::cerr << "Out weights file error: reached EOF before expected, on row " << i << " column " << j << "\n";
                } else if (file_state & std::ifstream::failbit) {
                    std::cerr << "Out weights file error on extraction or imterpretation of characters from the file, on row " << i << " column" << j << "\n";
                } else if (file_state & std::ifstream::badbit) {
                    std::cerr << "Out weights file error on stream while attempting to read row " << i << " column " << j << "\n";
                } else {
                    std::cerr << "Unknown out weights file error error occured on row " << i << " column " << j << "\n";
                    std::cerr << file_state << std::endl;
                }
            }
        }
    }

}

void create_accumulated_sum_matrices_from_parity_check_matrix(int_mtx_t* H, int_mtx_t* accumulated_rows_sum_matrix, int_mtx_t* accumulated_columns_sum_matrix) {
	int i, j;
	int sum;

	// First we calculate accumulated sum according to rows
	sum = 0;
	for (i = 0; i < H->m; ++i) {
		for (j = 0; j < H->n; ++j) {
			// We start with putting sum to the matrix since the first usable index (which will be when the first '1' will show up in the parity check matrix) should be 0
			put_int_element(accumulated_rows_sum_matrix, i, j, sum);
			sum += get_int_element(H, i, j);
		}
	}

	// Second we calculate accumulated sum according to columns
	sum = 0;
	for (j = 0; j < H->n; ++j) {
		for (i = 0; i < H->m; ++i) {
			// We start with putting sum to the matrix since the first usable index (which will be when the first '1' will show up in the parity check matrix) should be 0
			put_int_element(accumulated_columns_sum_matrix, i, j, sum);
			sum += get_int_element(H, i, j);
		}
	}
	
}

int extended_hamming_844_I_tx[] = {1,0,0,0,
				   0,1,0,0,
				   0,0,1,0,
				   0,0,0,1};


int extended_hamming_844_G[] = {1,1,1,1,0,0,0,0,
				1,1,0,0,1,1,0,0,
				1,0,1,0,1,0,1,0,
				0,1,1,0,1,0,0,1};


int extended_hamming_844_H1[] = {1,1,1,1,1,1,1,1,
				 0,1,0,1,0,1,0,1,
				 0,0,1,1,0,0,1,1,
				 0,0,0,0,1,1,1,1};

int extended_hamming_844_H1_tx[] = {1,0,0,0,
				    0,0,0,1,
				    0,0,1,0,
				    0,1,0,0};


int extended_hamming_844_H2[] = {1,1,1,1,0,0,0,0,
				 0,0,1,1,1,1,0,0,
				 0,0,0,0,1,1,1,1,
				 0,1,1,0,0,1,1,0};



int extended_hamming_844_H2_tx[] = {1,1,0,0,
				    0,1,1,0,
				    0,1,0,0,
				    0,0,1,1};

int sys_extended_golay_Hg[] = {
	1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,    1,  0,  1,  0,  1,  1,  1,  0,  0,  0,  1,  1,
	0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,    1,  1,  1,  1,  1,  0,  0,  1,  0,  0,  1,  0,
	0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,    1,  1,  0,  1,  0,  0,  1,  0,  1,  0,  1,  1,
	0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,    1,  1,  0,  0,  0,  1,  1,  1,  0,  1,  1,  0,
	0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,    1,  1,  0,  0,  1,  1,  0,  1,  1,  0,  0,  1,
	0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,    0,  1,  1,  0,  0,  1,  1,  0,  1,  1,  0,  1,
	0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,    0,  0,  1,  1,  0,  0,  1,  1,  0,  1,  1,  1,
	0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,    1,  0,  1,  1,  0,  1,  1,  1,  1,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,    0,  1,  0,  1,  1,  0,  1,  1,  1,  1,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,    0,  0,  1,  0,  1,  1,  0,  1,  1,  1,  1,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,    1,  0,  1,  1,  1,  0,  0,  0,  1,  1,  0,  1,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,    0,  1,  0,  1,  1,  1,  0,  0,  0,  1,  1,  1};


int extended_golay_Hg[] = { 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1,
			    1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1,
			    0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
			    0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1,
			    1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
			    0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1,
			    1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1,
			    0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1,
			    0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1,
			    0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1,
			    0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1,
			    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};


int extended_golay_Hg_tag[] = { 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1,
				0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1,
				1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
				0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0,
				1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0,
				1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1,
				0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0,
				0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
				0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1,
				0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0,
				1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0};


// MacWilliams and Sloane version (pp 65)
int extended_golay_MS_Hg[] = { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 	1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0,    
			       1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 	0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1,    
			       1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 	1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0,    
			       1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 	0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0,    
			       1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 	0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0,    
			       1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 	0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1,    
			       1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 	1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1,    
			       1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 	1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1,   
			       1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 	1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0,    
			       1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 	0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1,    
			       1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 	1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1,    
			         									          
			       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1   } ;


int golay_perm[] = {2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,  1, 24, //perm 1
		    1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 23,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, //perm 2
		    24, 23, 12, 16, 18, 10, 20, 14, 21,  6, 17,  3, 22,  8, 19,  4, 11,  5, 15,  7,  9, 13,  2,  1,//perm 3
		    24, 18, 23, 16, 20,  9, 15, 13,  8, 21, 10,  2, 11, 22,  3,  4, 6,  7, 12, 19, 14, 17,  5, 1}; //perm 4


int g_63_39_m[] = {
  1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0
, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0
, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0
, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0
, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0
, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0
, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0
, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0
, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0
, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
};


int h_63_39_m[] = {

  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1
, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1
, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0
, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1
, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1
, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0
, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1
, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1
, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1
, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1
};


void ml_decoder(int_mtx_t *C, real_vec_t *data, gf2_vec_t *dec_symb)
{

	int i, j, min_pos;
	// finiding the minimum ML
	double min_dist, curr_dist;


	min_dist = 1e9; 
	min_pos = -1;
	for (i=0; i<C->m; i++) {
		curr_dist = 0;
		for (j=0; j<C->n; j++)
			curr_dist += sqr(2*((double)get_int_element(C, i, j)) - 1 - data->data[j]);

		if (curr_dist < min_dist) {
			min_dist = curr_dist;
			min_pos  = i;
		}
	}
  
	for (i=0; i<C->n; i++)
		dec_symb->data[i] = get_int_element(C, min_pos, i);

}

void calc_all_codewords(int_mtx_t *G, int_mtx_t *C)
{
	int N = (int)pow(2, G->m), i, j;
	gf2_vec_t *info, *symb;
  
	info     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	symb     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
  
	init_gf2_vec(info, G->m);
	init_gf2_vec(symb , G->n);
	init_int_mtx(C,      N, G->n);  // C will hold both code word and code info

	for (i=0; i<N; i++) {
		for (j=0; j<G->m; j++)	
			info->data[G->m - 1 - j] = (i & (1<<j)) ? 1: 0;

		encode_symbols(info, G, symb);
		for (j=0; j<G->n; j++)
			put_int_element(C, i, j, symb->data[j]);      
	}

}

double_mtx_t* calc_ML_ber(ecc_t *ecc_p)
{
	// double curr_EBN0, n_errors, n_errors_uc;
	double n_errors, n_errors_uc;
	// int i, j, iter, length, dec_uc, w;
	int i, j, iter, length, dec_uc;
	double_mtx_t *ber;
	gf2_vec_t *info, *symb, *dec_symb;
	real_vec_t *data, *data_uc;
	int_mtx_t *C;
	int erf_n_errors_uc, erf_n_errors;

 
	info     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	symb     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	dec_symb = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	data     = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	data_uc  = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	C        = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
 
	int_mtx_t *G		= ecc_p->G;
	int N			= ecc_p->N;
	double start_EBN0_db	= ecc_p->start_EBN0_db;
	double end_EBN0_db	= ecc_p->end_EBN0_db;
	double EBN0_delta	= ecc_p->EBN0_delta;

  
	length = floor((end_EBN0_db - start_EBN0_db)/EBN0_delta) + 1;
	ber = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
	init_double_mtx(ber, 3, length);
  
	init_gf2_vec(info, G->m);
	init_gf2_vec(dec_symb, G->n);
	init_gf2_vec(symb , G->n);
	init_real_vec(data , G->n);
	init_real_vec(data_uc , G->n);

	printf("Preparing all codewords\n");
	calc_all_codewords(G, C);
	printf("Finished preparing all codewordss\n");

	for (i=0; i<length; i++) {
		printf("ML, Calculating BER = %g, n=%d, k=%d\n", start_EBN0_db + EBN0_delta*i, G->n, G->m);
		put_double_element(ber, 0, i, start_EBN0_db + EBN0_delta*i);
		//srand( (unsigned)time( NULL ) );

		erf_n_errors_uc = erf_n_errors = n_errors_uc = 	n_errors = 0; 
		for (iter=0; iter<N; iter++) {
			if (0 == (iter%100)) printf("ML, BER = %g, iter = %d out of %d, n=%d, k=%d\n", start_EBN0_db + EBN0_delta*i, iter, N, G->n, G->m);
			init_symbol_data(info);
			encode_symbols(info, G, symb); 	         
			init_measured_data(symb, data, start_EBN0_db + EBN0_delta*i, ((double)(G->m)/(double)G->n));
			init_measured_data(info, data_uc, start_EBN0_db + EBN0_delta*i, 1.0);

			ml_decoder(C, data, dec_symb);

			// coded data
			for (j=0; j<G->n; j++) {
				n_errors += (symb->data[j] != dec_symb->data[j] ? 1 : 0);
				erf_n_errors += (symb->data[j] != dec_symb->data[j] ? 1 : 0);
			}
	  
			// uncoded data
			for (j=0; j<G->m; j++) {
				dec_uc = (data_uc->data[j] > 0) ? 1 : 0;
				n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
				erf_n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
			}

			if (((iter+1)%ERF_QUNT == 0)) { 
				update_erf_file(ecc_p->ml_fname,
						start_EBN0_db + EBN0_delta*i,  
						ERF_QUNT, 
						G->n, G->m, 
						erf_n_errors_uc,  
						erf_n_errors 
						);			
				erf_n_errors = erf_n_errors_uc = 0;
			}


		}
		put_double_element(ber, 1, i, n_errors/ (double)(N*G->n));
		put_double_element(ber, 2, i, n_errors_uc/ (double)(N*G->m));
	}
  
	dump_double_mtx(ber, ecc_p->ml_fname);
	print_double_mtx(ber);
	return ber;

}

// The following function calculates the error performance (BER) of standard BP decoding applied on a code transmitted over AWGN channel.
// The '_uc' suffix stands for uncoded transmission.
double_mtx_t * calc_TG_ber(ecc_t *ecc_p)
{
	int_mtx_t *G, *H;
	int N;
	double start_EBN0_db, end_EBN0_db, EBN0_delta, alpha;
	int TG_iter;
	int erf_n_errors_uc, erf_n_errors;
	//double curr_EBN0, n_errors, n_errors_uc;
	double n_errors, n_errors_uc;
	int i, j, iter, length, dec_uc, w;
	double_mtx_t *ber;
	gf2_vec_t *info, *symb, *dec_symb;
	real_vec_t *data, *data_uc;  
	real_vec_t *lambda, *mu;
	double tmp_tg_time;

	N		= ecc_p->N;
	end_EBN0_db	= ecc_p->end_EBN0_db;
	start_EBN0_db = ecc_p->start_EBN0_db;
	EBN0_delta	= ecc_p->EBN0_delta;
	alpha		= ecc_p->tg_alpha;
	G    		= ecc_p->G;
	H    		= ecc_p->H;
	TG_iter  	= ecc_p->TG_iter;


	lambda	= (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	mu		= (real_vec_t*)ecc_malloc(sizeof(real_vec_t)); 
	init_real_vec(lambda,  G->n);
	init_real_vec(mu, G->n);

	info     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	symb     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	dec_symb = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	data     = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	data_uc  = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));

  
	length = floor((end_EBN0_db - start_EBN0_db)/EBN0_delta) + 1;
	ber = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
	// Note there are 3 values allocated for each SNR value: The SNR value itself, the BER using the given code, and the BER when transmitting uncoded data
	init_double_mtx(ber, 3, length);
  
	init_gf2_vec(info, G->m); // Vector that includes information bits
	init_gf2_vec(dec_symb, G->n); // Vector that includes decoded bits (after channel transmission and decoding)
	init_gf2_vec(symb , G->n); // Vector that includes encoded bits
	init_real_vec(data , G->n); // Vector that includes real values - the measurements received after transmission over AWGN channel
	init_real_vec(data_uc , G->n);

	// The following loop runs on all SNR values
	for (i=0; i<length; i++) {
		printf("TG, Calculating SNR = %g, n=%d, k=%d\n", start_EBN0_db + EBN0_delta*i, G->n, G->m);
		put_double_element(ber, 0, i, start_EBN0_db + EBN0_delta*i);
		//srand( (unsigned)time( NULL ) );
		// Initialize error counters
		n_errors = 0; 
		n_errors_uc = 0;
		erf_n_errors_uc = erf_n_errors = 0;
		// Simulate N codewords
		for (iter=0; iter<N; iter++) {
			// Print status every 100 codewords
			if (0 == (iter%100)) printf("TG, SNR = %g, iter = %d out of %d, n=%d, k=%d\n", start_EBN0_db + EBN0_delta*i, iter, N, G->n, G->m);
			// Generating information bits
			init_symbol_data(info);
			// Encoding information bits to a codeword
			encode_symbols(info, G, symb); 	      
			// Transmit over AWGN channel   
			init_measured_data(symb, data, start_EBN0_db + EBN0_delta*i, ((double)(G->m)/(double)G->n)); // Transmission of the codeword
			init_measured_data(info, data_uc, start_EBN0_db + EBN0_delta*i, 1.0); // Transmission of only the information bits

			// for (w=0; w<lambda->n; w++)  lambda->data[w] = -4*data->data[w]; // Initialize channel LLR values into 'lambda' vector
			// No is 1/R*EbNo (when EbNo is in linear scale)
			for (w=0; w<lambda->n; w++)  lambda->data[w] = -4*((double)(G->m)/(double)G->n)*pow(10, (start_EBN0_db + EBN0_delta*i)/10)*data->data[w]; // Initialize channel LLR values into 'lambda' vector
			// tanner_graph_decoder(H, lambda, dec_symb, mu, TG_iter, alpha); // Applying BP decoding
#ifdef __use_soft_TG__
			soft_tanner_graph_decoder(H, lambda, dec_symb, mu, TG_iter, alpha, &tmp_tg_time, ecc_p->hidden_layers_weights, ecc_p->out_weights, ecc_p->accumulated_rows_sum_matrix, ecc_p->accumulated_columns_sum_matrix);
#else
			tanner_graph_decoder(H, lambda, dec_symb, mu, TG_iter, alpha, &tmp_tg_time); // Applying BP decoding
#endif

			// Now counting the number of error bits
			// coded data
			for (j=0; j<G->n; j++) {
				n_errors += (symb->data[j] != dec_symb->data[j] ? 1 : 0); // Note that we count errors over ALL code bits (and not only over information bits)
				erf_n_errors += (symb->data[j] != dec_symb->data[j] ? 1 : 0);
			}
	  
			// uncoded data
			for (j=0; j<G->m; j++) {
				dec_uc = (data_uc->data[j] > 0) ? 1 : 0;
				n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
				erf_n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
			}

			// Log the performance simulated so far to a file
			if (((iter+1)%ERF_QUNT == 0)) { 
				update_erf_file(ecc_p->tg_fname,
							  start_EBN0_db + EBN0_delta*i,  
							  ERF_QUNT, 
							  G->n, G->m, 
							  erf_n_errors_uc,  
							  erf_n_errors);			
				erf_n_errors = erf_n_errors_uc = 0;
			}
		}
		// Storing BER for the coded and uncoded transmissions
		put_double_element(ber, 1, i, n_errors/ (double)(N*G->n));
		put_double_element(ber, 2, i, n_errors_uc/ (double)(N*G->m));
	}
  	
  	// Saving the BER matrix to a file
	dump_double_mtx(ber, ecc_p->tg_fname);
	// Printing the BER matrix
	print_double_mtx(ber);
	return ber;
}

// Algoritm 3 operates on a list S of length N of genererating set; extracting element at specified 
// position from the list is done by pt_get; operating the operator on two elements is sone with 
// pt operate
void* Algorithm_3(void *S, int N, void* (*pt_get)(void*, int), void (*pt_operate)(void*, int, int))
{
	int i, j;

	// pair of random integers
	i = rand() % N; 
	do { j = rand() % N; } while (j ==i);
  
	// Martzi: debug prints
	//printf("i and j for permutation are -----\n");
	//printf("i is %d and j is %d\n", i, j);
	pt_operate(S, i, j);

	return pt_get(S, i); 
}

//Create the sutomorphism group of RM(1, m)
void* create_generating_set_RM1(int m, int N)
{
	int_mtx_t *mtx;
	int k, i, j, w, pos, rand_pos;

	k = m*m;
  
	ecc_assert(2*k < N, "N is not large enough\n");

	mtx = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t)*N);

	pos = 0;


	// 1. Creating the first k elements
  
	// 1.1 Createing the first k^2 - k
	for (i=0; i<m; i++)
		for (j=0; j<m; j++) {
			if (i!=j) {
				init_int_mtx(mtx+pos, m, m+1);
				for (w=0; w<m; w++) put_int_element(mtx+pos, w, w, 1);
				put_int_element(mtx+pos, i, j, 1);
				pos++;
			}

		}


	// 1.2 Createing the next k elements
	for (i=0; i<m; i++) {
		init_int_mtx(mtx+pos, m, m+1);
		for (w=0; w<m; w++) put_int_element(mtx+pos, w, w, 1);
		put_int_element(mtx+pos, i, m, 1);
		pos++;
	}


	// 2 Copying the first k as long as there enough room for the whole copy
	while (pos+k<N) {
		for (w=0; w<k; w++) { 
			init_int_mtx(mtx+pos, m, m+1);
			for (i=0; i<m; i++)
				for (j=0; j<m+1; j++)
					put_int_element(mtx+pos, i, j, get_int_element(mtx+pos-k, i, j)); 
			pos++;
		}
	}

	// 3 randomaize the rest of the elements
	while (pos<N) {
		init_int_mtx(mtx+pos, m, m+1);
		rand_pos = rand()%pos;
		for (i=0; i<m; i++)
			for (j=0; j<m+1; j++)
				put_int_element(mtx+pos, i, j, get_int_element(mtx+rand_pos, i, j)); 
		pos++;
	  
	}


	return (void*)mtx; 
}

void* perm_get(void* S, int pos)
{
	int_mtx_t *mtx = (int_mtx_t *)S; 
	return (void*)(mtx + pos);
}

void perm_operate(void* S, int i, int j)
{
	int_mtx_t *mtx = (int_mtx_t *)S; 
	int_mtx_t *ip, *jp, *tmp;
	int n, w;
  
	ecc_assert(mtx->m == 1, "error woth permutation matrices size: %d vs. required 1\n", mtx->m);

	n  = mtx->n;
	ip = mtx + i;
	jp = mtx + j;
	// Temporary variable to hold the permutation calculated
	tmp = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx(tmp, 1, n);

	//  printf("ip: ");
	//  print_int_mtx(ip);
	//  printf("\njp: ");
	//  print_int_mtx(jp);
 
	// Calculating the permutation according to given i and j
	for (w=0; w<n; w++)
		put_int_element(tmp, 0, get_int_element(jp, 0, w), get_int_element(ip, 0, w));

	//  printf("\ntmp: ");
	//  print_int_mtx(tmp);
	//  getchar();

	// Placing the permutation calculated back to S
	for (w=0; w<n; w++)
		put_int_element(ip, 0, w, get_int_element(tmp, 0, w));
  
	free_int_mtx(tmp);
}



  
void RM1_operate(void* S, int i, int j)
{
	int_mtx_t *mtx = (int_mtx_t *)S; 
	int_mtx_t *ip, *jp, *tmp;
	int w, v, z, m, sum;
  
	m  = mtx->m;
	ip = mtx + i;
	jp = mtx + j;
	tmp = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx(tmp, m, m+1);

	// 1. opertating [Ax|bx] on [Ay|by] on [Ax*Ay| Ax*by + by]

	// 1.1 multiplying AxAy
	for (w=0; w<m; w++)
		for (v=0; v<m; v++) {
			sum = 0; 
			for (z=0; z<m; z++)
				sum += get_int_element(jp, w, z) * get_int_element(ip, z, v);
			put_int_element(tmp, w, v, sum%2);
		}

	// 1.2 multiplying Ax*by
	for (w=0; w<m; w++) {
		sum = 0; 
		for (z=0; z<m; z++)
			sum += get_int_element(jp, w, z) * get_int_element(ip, z, m);
		sum += get_int_element(jp, w, m);	
		put_int_element(tmp, w, m, sum%2);
	}
  
	// 2 copying the result back to i
	for (w=0; w<m; w++)
		for (v=0; v<m+1; v++)
			put_int_element(ip, w, v,get_int_element(tmp, w, v));

	free_int_mtx(tmp);
}


int_mtx_t* create_generator_RM1(int m)
{
  
	int_mtx_t* p;
	int i, j, n;
  
	n = (int)pow(2,m);

	p = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx(p, m+1, n);

	for (i=0; i<n; i++) {
		// handling the first line
		put_int_element(p, 0, i, 1);
      
		// handling the other lines
		for (j=0; j<m; j++)
			put_int_element(p, m-j, i, ((i&(1<<j)) ? 1 : 0));
	}
	return p;
}

void RM1_get_permutation(int_mtx_t* orig, int_mtx_t* orig_tx, int_mtx_t* permuted, int *permutation)
{
	int n, m, i, j, z, sum, identical;
	int_mtx_t *tmp, *Hi, *Hi_mod;

	n = orig->n; 
	m = orig->m;
	tmp    = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	Hi     = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	Hi_mod = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));

	init_int_mtx(tmp, m, n);
	init_int_mtx(Hi,     m, n);
	init_int_mtx(Hi_mod, m, n);

	//  printf("orig:\n");
	//  print_int_mtx(orig);
	//  print_int_mtx(permuted);
	//  getchar();

	for (i=0; i<n; i++)
		put_int_element(tmp, 0, i, get_int_element(orig, 0, i));

	for (i=0; i<m-1; i++)
		for (j=0; j<n; j++) {
			sum = 0;
			for (z=0; z<m-1; z++)
				sum += get_int_element(permuted, i, z) * get_int_element(orig, z+1, j);
	
			sum += get_int_element(permuted, i, m-1);
			put_int_element(tmp, i+1, j, sum%2);
		}

	// performing the transform to the original matrix
	multiply_int_mtx(orig_tx, orig, Hi);
	multiply_int_mtx(orig_tx, tmp,  Hi_mod);

	//  printf("Hi\n");
	//  print_int_mtx(Hi);
	//  printf("Hi_mod\n");
	//  print_int_mtx(Hi_mod);
	//  getchar();

  
	// Creating the auto-morph matrix
	for (i=0; i<n; i++) {
		z = 0;
		do {
			identical = 1;
			for (j=0; j<m; j++)
				if  (get_int_element(Hi_mod, j, z) != get_int_element(Hi, j, i))
					identical =0;

		} while ((!identical)  && (++z<n));
		if (identical) 
			permutation[i] = z;
		else
			ecc_assert(0, "Couldn't find matching columns\n");
	}  

	//  printf("permutation:");
	//  for(i=0; i<n; i++)
	//    printf("%d ", permutation[i]);
	//  getchar();
	//  getchar();


	free_int_mtx(tmp);
	free_int_mtx(Hi);
	free_int_mtx(Hi_mod);
}

int_mtx_t *create_RM1_permutation(int m)
{
	// int_mtx_t *S, *tmp; 
	int_mtx_t *S; 
	int i, N;

	N=m*m*2+2;

	S = (int_mtx_t *) create_generating_set_RM1(m, N);
  
	for (i=0; i<60; i++) 
		Algorithm_3(S, N, perm_get, RM1_operate);
  
	return S;
}

int_mtx_t *create_RM1_3_permutation()
{
	return create_RM1_permutation(3);
}

int_mtx_t *create_permutation(int n, int nr_perm, int *perm, int perm_list_length)
{
	int_mtx_t *S;
	int i, j, lim, N, rand_pos, pos;
  

	ecc_assert(max(2*nr_perm+1, 10) < perm_list_length, "Error with perm_list_length, %d vs. %d\n", max(2*nr_perm+1, 10), perm_list_length);
	N = perm_list_length;
	S = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)*N); // This is the array of permutations (the Automorphism group)

	//reducing 1 from the permutation
	// We reduce 1 to be consistent with C array addressing (where the first element has index 0 and not 1 as in the permutations given)
	for (i=0; i<nr_perm*n; i++)
		perm[i]--;

	// Assigning (with repetitions) the nr_perm given permutations to an array of permutations S
	// handling the first 2*k (or more)
	lim = N / nr_perm;
	for (i=0; i<lim; i++)
		for (j=0; j<nr_perm; j++) {
			init_int_mtx_from_array(S + i*nr_perm + j, 1, n, perm + j*n);
		}

	// Since the number of permutations to be used N is not necessarily a multiplication of nr_perm, this section assigns the rest of (N modulo nr_perm) permutations required to S
	//handling the rest
	pos = lim*nr_perm;
	while (pos<N) {
		rand_pos = rand() % nr_perm;
		init_int_mtx_from_array(S + pos, 1, n, perm + rand_pos*n);
		pos++;	  
	}

	// Performing initialization of 60 Product - Replacement iterations to S
	for (i=0; i<60; i++) 
		Algorithm_3(S, N, perm_get, perm_operate);
  
	return S;

}

int_mtx_t *create_ext_golay_permutation()
{
	return create_permutation(24, 4, golay_perm, N_perm_list);
}

int *gt_perm;

int *create_BCH_permutation(int m)
{
	int n = (1<<m) - 1;
	int *perm = (int*)ecc_malloc(sizeof(int)*m*n);
	int i, j;
	for (j=0; j<m; j++) {
		for (i=0; i<n; i++) {
			perm[j*n + i] = ((1 + (1<<j)*i) % n) + 1;
			printf("%d ", perm[j*n + i]);
		}
		printf("\n");
	}
	gt_perm = perm;
	return perm;
}

int_mtx_t *create_BCH_31_permutation()
{
	int *perm  = create_BCH_permutation(5);
	int_mtx_t *p = create_permutation(31, 5, perm, N_perm_list);
	//free(perm);
	return p;
}

int_mtx_t *create_BCH_63_permutation()
{
	int *perm  = create_BCH_permutation(6);
	int_mtx_t *p = create_permutation(63, 6, perm, N_perm_list);
	free(perm);
	return p;
}

int_mtx_t *create_BCH_127_permutation()
{
	int *perm  = create_BCH_permutation(7);
	int_mtx_t *p = create_permutation(127, 7, perm, N_perm_list*10);
	free(perm);
	return p;
}



int is_zeros(gf2_vec_t* in_vec)
{
	int i, is_zero=1;
  
	for (i=0; i<in_vec->n; i++)
		if (0!=in_vec->data[i]) 
			is_zero = 0;

	free(in_vec->data);
	free(in_vec);
	return is_zero;

}


int is_finite(real_vec_t* in_vec)
{
	int i, is_finite=1;
  
	for (i=0; i<in_vec->n; i++)
		if ((!finite(in_vec->data[i])) || (my_isnan(in_vec->data[i])))
			is_finite = 0;
  
	return is_finite;
}

gf2_vec_t* gf2_mult(int_mtx_t *H, gf2_vec_t *dec_symb)
{
	gf2_vec_t* syndrome;
	int i, j, sum;

	ecc_assert(dec_symb->n == H->n, "Length doesn't match %d vs. %d\n", dec_symb->n, H->n);

	syndrome = (gf2_vec_t*)ecc_malloc(sizeof(gf2_vec_t));
	init_gf2_vec(syndrome, H->m);

	for (i=0; i<H->m; i++) {
		sum = 0;
		for (j=0; j<H->n; j++)
			sum += get_int_element(H, i, j) * dec_symb->data[j];
		syndrome->data[i] = sum % 2;
	}
	return syndrome;
}


void apply_inv_perm_int(gf2_vec_t *dec_symb, int_vec_t* PHI)
{
	gf2_vec_t *tmp;
	int i;

	ecc_assert(dec_symb->n == PHI->n, "Length doesn't match %d vs. %d\n", dec_symb->n, PHI->n);

	tmp = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	init_gf2_vec(tmp, dec_symb->n);

	for (i=0; i<dec_symb->n; i++) {
		tmp->data[PHI->data[i]] = dec_symb->data[i];
	}

	for (i=0; i<dec_symb->n; i++)
		dec_symb->data[i] = tmp->data[i];

	free_gf2_vec(tmp);
}

void apply_perm_real(real_vec_t *lambda, int *ord)
{
	int i;
	double *tmp;

	tmp = (double*)ecc_malloc(sizeof(double)*lambda->n);
  
	for (i=0; i<lambda->n; i++)
		tmp[ord[i]] = lambda->data[i];

	for (i=0; i<lambda->n; i++)
		lambda->data[i] = tmp[i];

	free(tmp);
}

void apply_perm_real_int_mtx(real_vec_t *lambda, int_mtx_t *perm)
{
	int i;
	double *tmp;

	tmp = (double*)ecc_malloc(sizeof(double)*lambda->n);
  
	for (i=0; i<lambda->n; i++)
		tmp[get_int_element(perm, 0, i)] = lambda->data[i];

	for (i=0; i<lambda->n; i++)
		lambda->data[i] = tmp[i];

	free(tmp);
}

void apply_perm_int_vec(int_vec_t* PHI, int *ord)
{
	int i;
	int *tmp;

	tmp = (int*)ecc_malloc(sizeof(int)*PHI->n);
  
	for (i=0; i<PHI->n; i++)
		tmp[ord[i]] = PHI->data[i];

	for (i=0; i<PHI->n; i++)
		PHI->data[i] = tmp[i];

	free(tmp);
}
void apply_perm_int_vec_int_mtx(int_vec_t* PHI, int_mtx_t *perm)
{
	int i;
	int *tmp;

	tmp = (int*)ecc_malloc(sizeof(int)*PHI->n);
  
	for (i=0; i<PHI->n; i++)
		tmp[get_int_element(perm, 0, i)] = PHI->data[i];

	for (i=0; i<PHI->n; i++)
		PHI->data[i] = tmp[i];

	free(tmp);
}

void apply_perm_int_mtx(int_mtx_t* mtx, int *ord)
{
	int_mtx_t *tmp;
	int i, j;
  

	tmp = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  
	init_int_mtx(tmp, mtx->m, mtx->n);
  
	for (i=0; i<mtx->m; i++)
		for (j=0; j<mtx->n; j++)
			put_int_element(tmp, i, ord[j], get_int_element(mtx, i, j));

	copy_int_mtx(mtx, tmp);

	free(tmp->data);
	free(tmp);

}




static void RRD_decoder(ecc_t *ecc_p, real_vec_t *lambda, gf2_vec_t *dec_symb, real_vec_t *mu,  int_mtx_t *per_list, int *total_nt_iter)
{
	// int i1, i2, i3; 
	int i2, i3; 
	double alpha; 
	int_vec_t *PHI; 
	// int *ord, i, j; 
	int *ord, i; 
	real_vec_t *s;
	int_mtx_t *auto_morph_element, *H_mod;
	gf2_vec_t* prev_dec_symb = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));

	init_gf2_vec(prev_dec_symb, dec_symb->n);
  
 
	int_mtx_t *H		= ecc_p->H;
	int_mtx_t *H_tx	= ecc_p->H_tx;
	double alpha0		= ecc_p->rrd_alpha;
	int I1		= ecc_p->I1;
	int I2		= ecc_p->I2;
	int I3		= ecc_p->I3;
	double tmp_tg_time;

	// Creating modified H
	// Martzi: Patch for RM code
	H_mod = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx(H_mod, H->m, H->n);
	if (ecc_p->H_tx_present)
		multiply_int_mtx(H_tx, H, H_mod);
	else
		copy_int_mtx(H_mod, H);

	// holding the permuted inputs
	s = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));	
	init_real_vec(s, H->n);

	// Place holders for the permutation
	ord = (int*)ecc_malloc(sizeof(int)*H->n);
	PHI = (int_vec_t*)ecc_malloc(sizeof(int_vec_t)*H->n); // Martzi: should we multiply by H->n? it looks unnecessary
	init_int_vec(PHI, H->n);
	for (i=0; i<H->n; i++) PHI->data[i] = i;
   
	// Starting with RRD decoding
	// Outer loop is for damping coefficient 'alpha' values
	alpha = alpha0;
	for (i3=0; i3<I3; i3++) 	{
		//printf("-=-=-=- i3 = %d\n", i3);
      
		for (i=0; i<H->n; i++) 	{
			PHI->data[i] = i;
			s->data[i] = lambda->data[i];
		}

		//printf("alpha: %g\n", alpha);
		//printf("lambda: "); print_rvec(lambda);

		for (i2=0; i2<I2; i2++)	{
			//printf("-=-=-=- i2 = %d\n", i2);
	  

			//printf("s [pre]: "); print_rvec(s);


			copy_gf2_vec(dec_symb, prev_dec_symb);
			// tanner_graph_decoder(H_mod, s, dec_symb, mu, I1, alpha);
			tanner_graph_decoder(H_mod, s, dec_symb, mu, I1, alpha, &tmp_tg_time);
			copy_rvec(mu, s);
	  

			//printf("s [post]: "); print_rvec(s);
			//printf("dec_symbols: "); print_gf2vec(dec_symb);
			//getchar();
	  
			if (!is_finite(s)) ecc_assert(0,  "Overflow issue\n");
			// Check if BP result is a valid codeword
			if (is_zeros(gf2_mult(H_mod, dec_symb))) {
				//printf("Is zero\n");
				//printf("PHI:");
				//for (i=0; i<H->n; i++) 
				//  printf("%d ", PHI->data[i]);
				//printf("\n");
				//getchar();

				// If so, apply inverse of all permutations and log the successful alpha value
				apply_inv_perm_int(dec_symb, PHI);

				ecc_p->alpha_success = alpha;

				//for (i=0; i<H->n; i++) 
				//  printf("%d ", dec_symb->data[i]);
				//printf("\n");
				//getchar();
	      
				free (ord);
				free_int_vec(PHI);
				free_real_vec(s);
				free_int_mtx(H_mod);
				free_gf2_vec(prev_dec_symb);
	      
	      		// For perfornance record, save the total number of BP iterations applied for current codeword
				*total_nt_iter = I1*(i2+1)*(i3+1);
				return;
			}
	 	  
 

			if (i2+1 < I2) 	{
				// Apply Product - Replacement algorithm to produce permutations from the automorphism group of the code
				// This also updates the automorphism group of the code
				//auto_morph_element = Algorithm_3(per_list, N_perm_list, perm_get, RM1_operate);
				auto_morph_element = (int_mtx_t *) Algorithm_3(per_list, N_perm_list, ecc_p->perm_get, ecc_p->perm_operate);
	      
				// special handling for RM1(3)
				if (ecc_p->H_tx_present) {
					RM1_get_permutation(H, H_tx, auto_morph_element, ord);	    
		  
					//printf("ord:");
					//for (i=0; i<H->n; i++) 
					//  printf("%d ", ord[i]);
					//printf("\n");
		  
					apply_perm_real(s, ord);
		  
					//printf("PHI (before):");
					//for (i=0; i<H->n; i++) 
					//  printf("%d ", PHI->data[i]);
					//printf("\n");
		  
					apply_perm_int_vec(PHI, ord);
		  
					//printf("PHI:");
					//for (i=0; i<H->n; i++) 
					//  printf("%d ", PHI->data[i]);
					//printf("\n");
					//getchar();
				} else	{
					apply_perm_real_int_mtx(s, auto_morph_element);
					apply_perm_int_vec_int_mtx(PHI, auto_morph_element);
				}

			}

		}

		alpha = alpha0 + (1-alpha0)*(((double)i3) + 1.0) / (((double)I3) - 1.0);      

	}
	ecc_assert(is_finite(s), "Overflow issue\n");
	if (!is_finite(s))
		copy_gf2_vec(prev_dec_symb, dec_symb);

	apply_inv_perm_int(dec_symb, PHI);

	*total_nt_iter = I1*(i2+1)*(i3+1);



	free (ord);
	free_int_vec(PHI);
	free_int_mtx(H_mod);
	free_real_vec(s);
	free_gf2_vec(prev_dec_symb);

}

static void create_cgg(int_mtx_t *cgg, int line, int_mtx_t *H, int c_rest)
{
	int i, j;
	
	if (c_rest) {
		// special case of extended golay
		for (i=0; i<cgg->n; i++) {		       
			for (j=0; j<cgg->n; j++)
				if (i<cgg->n-1) {
					if (j< cgg->n - c_rest)
						put_int_element(H, i, j, get_int_element(cgg, line, (j+i)%(cgg->n - c_rest)));		
					else						
						put_int_element(H, i, j, 0);		
				} else
					put_int_element(H, i, j, 1);
		}
	} else {
		for (i=0; i<cgg->n; i++) {
			for (j=0; j<cgg->n; j++)
				put_int_element(H, i, j, get_int_element(cgg, line, (j+i)%cgg->n));		
		}
	}
}


static void TEST_decoder(ecc_t *ecc_p, real_vec_t *lambda, gf2_vec_t *dec_symb, real_vec_t *mu, int_mtx_t *per_list, int *total_nt_iter, double *nr_time, double *sum_tg_time)
{
	
	Decoder Eff_BPdec("L_Ilan_sparse_BCH_63_45.txt", (ecc_p->G)->n - (ecc_p->G)->m, (ecc_p->G)->n);
	Eff_BPdec.SetDecoder("BP");
	Eff_BPdec.SetParameters(ecc_p->I1, ecc_p->tg_alpha);
	Eff_BPdec.SetBPSaturation(false, 30.0);
	// Eff_BPdec._LogLikelihoodRatio[8] = 0;
	double *mu1 = new double[lambda->n];
	double *mu2 = new double[lambda->n];

	int i, j, l, *idx, is_zero, zero_not_found, min_pos, sum;
	// int i1, i2;
	int i2;
	double *distance, min_distance;
	int_mtx_t  *H, *H_tx, *R, *S;
	int_vec_t *PHI; 
	int *ord; 
	int_mtx_t *auto_morph_element, *H_mod;;
	// int ex_golay_legal_lines[] = {0, 2, 4, 6, 7, 8, 10, 11, 12, 13, 15, 17, 18, 19, 20, 21, 25, 26, 29, 30, 31, 32};
	real_vec_t *s;
	Timer t;
	double tmp_tg_time;
	gf2_vec_t* prev_dec_symb = (gf2_vec_t *) ecc_malloc(sizeof(gf2_vec_t));

	int I1		= ecc_p->I1;
	int I2		= ecc_p->I2;

	H	= ecc_p->H;
	H_tx	= ecc_p->H_tx;



	// Creating modified H
	H_mod = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx(H_mod, H->m, H->n);
	if (ecc_p->H_tx_present)
		multiply_int_mtx(H_tx, H, H_mod);
	else
		copy_int_mtx(H_mod, H);

	init_gf2_vec(prev_dec_symb, dec_symb->n);

	ecc_assert(ecc_p->mbbp_l != 0, "ERROR: ecc_p->mbbp_l = %d\n", ecc_p->mbbp_l);
	ecc_assert(ecc_p->mbbp_i != 0, "ERROR: ecc_p->mbbp_l = %d\n", ecc_p->mbbp_i);
	

	R = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t)); // Holds the output of each decoder
	S = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t)); // For each decoded word, holds the parity check results for each parity check of the code
	idx = (int*)ecc_malloc(sizeof(int)*ecc_p->mbbp_l);
	memset(idx, 0, sizeof(int)*ecc_p->mbbp_l);
	distance = (double*)ecc_malloc(sizeof(double)*ecc_p->mbbp_l);
	

	init_int_mtx(R, ecc_p->mbbp_l, ecc_p->H->n);
	init_int_mtx(S, ecc_p->mbbp_l, ecc_p->H->n); // Martzi: The second dimension of S should be H->m and not H->n (as there are H->m parity checks in H)

	// holding the permuted inputs
	s = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));	
	init_real_vec(s, H->n);


	// Place holders for the permutation
	ord = (int*)ecc_malloc(sizeof(int)*H->n);
	PHI = (int_vec_t*)ecc_malloc(sizeof(int_vec_t)*H->n); // Martzi: the multiplication with H->n looks unneccesary
	init_int_vec(PHI, H->n);

	*total_nt_iter = 0;
	*sum_tg_time = 0;
	t.start();
	// For every decoder employed
	for (l = 0; l<ecc_p->mbbp_l; l++) {

		
		// initial permutation
		for (i=0; i<H->n; i++) 	{
			PHI->data[i] = i;
			s->data[i] = lambda->data[i];
		}
		
		// Martzi: debug prints
		//printf("LLRs in are -----\n");
		//print_rvec(s);

		
		// Outer loop
		for (i2=0; i2<I2; i2++)	{
			copy_gf2_vec(dec_symb, prev_dec_symb);
			tmp_tg_time = 0;
#ifndef __Eff_mRRD__
	#ifdef __use_soft_TG__
			soft_tanner_graph_decoder(H_mod, s, dec_symb, mu, I1, (H->n > 63 ? 1 : 0.7), &tmp_tg_time, ecc_p->hidden_layers_weights, ecc_p->out_weights, ecc_p->accumulated_rows_sum_matrix, ecc_p->accumulated_columns_sum_matrix);
	#else
			tanner_graph_decoder(H_mod, s, dec_symb, mu, I1, (H->n > 63 ? 1 : 0.7),&tmp_tg_time);
	#endif
#endif		
			// Martzi: debug prints	
			//printf("LLRs after BP are -----\n");
			//print_rvec(mu);
#ifdef __Eff_mRRD__
			reduction_2_Eff_TG_Decoder_a(&Eff_BPdec, s, dec_symb, mu, I1, (H->n > 63 ? 1 : 0.7), &tmp_tg_time);
#endif
			*sum_tg_time += tmp_tg_time;
			// printf("The time taken by %d BP iterations in usec is %f\n", I1, tmp_tg_time);
			copy_rvec(mu, s);
			if (!is_finite(s)) ecc_assert(0,  "Overflow issue\n");

			if (is_zeros(gf2_mult(H_mod, dec_symb))) {
				// Martzi: debug prints
				//printf("Current LLRs represent a valid codeword\n");
				break;

			}
	

			if (i2+1 < I2) 	{
				auto_morph_element = (int_mtx_t *)Algorithm_3(per_list, N_perm_list, ecc_p->perm_get, ecc_p->perm_operate);
				// Martzi: debug prints
				//printf("permutation list -----\n");
				//for (i=0; i<N_perm_list; i++)
					//print_int_mtx(per_list+ i);
	      
				if (ecc_p->H_tx_present) {
					RM1_get_permutation(H, H_tx, auto_morph_element, ord);	    
		  
					apply_perm_real(s, ord);
					apply_perm_int_vec(PHI, ord);
		  
				} else	{
					apply_perm_real_int_mtx(s, auto_morph_element);
					apply_perm_int_vec_int_mtx(PHI, auto_morph_element);
					// Martzi: debug prints
					//printf("LLRs after permutation are -----\n");
					//print_rvec(s);
					//printf("Cummulative permutation is -----\n");
					//print_ivec(PHI);
				}
			}

		}

		// Count number of iterations executed for the current codeword
		*total_nt_iter += I1*(i2+1);
		// Martzi: debug prints
		//printf("Total number of BP iters performed for current codeword is %d\n", *total_nt_iter);

		ecc_assert(is_finite(s), "Overflow issue\n");
		if (!is_finite(s)) {
			// Martzi: debug prints
			//printf("There is an overflow issue, using previously decoded symbols\n");
			copy_gf2_vec(prev_dec_symb, dec_symb);
		}
		
		apply_inv_perm_int(dec_symb, PHI);
		

		for (i=0; i<dec_symb->n; i++) 
			put_int_element(R, l, i, dec_symb->data[i]);

		for (i=0; i<H->m; i++) {
			sum = 0; 
			for (j=0; j<H->n; j++)
				sum += dec_symb->data[j] * get_int_element(H, i, j); // c*H_L'
			put_int_element(S, l, i, sum%2);
		}

		// creating S according to the originl paper

	}
	
	
	zero_not_found = 1;
	for (i=0; i<S->m; i++) {
		// For every decoder
		is_zero = 0;
		for (j=0; j<S->n; j++)
			// For every parity check result of the current decoder's output
			if (get_int_element(S, i, j)) {
				// Then the ouput is not a valid codeword
				is_zero = 1;
			}
		if (!is_zero) { 
			// Then decoder i has a valid codeword
			idx[i] = 1;
			zero_not_found = 0;
		}
	}
	
	if (zero_not_found) 
		// None of the decoders has output a valid codeword. We will consider all output codewords
		for (i=0; i<ecc_p->mbbp_l; i++) {
			idx[i] = 1;
			//printf("%d ", idx[i]);
		}

		
	// finiding the minimum distance (assuming lambda is -4*y)
	min_distance = 1e99;
	for (i=0; i<S->m; i++) {
		if (idx[i]) {
			distance[i] = 0;
			for (j=0; j<S->n; j++)
				distance[i] += sqr((-lambda->data[j]/4) - (2*(double)get_int_element(R, i, j) - 1));
			if (distance[i]<min_distance) {
				min_distance = distance[i];
				min_pos = i;
			}
			//printf("%g ", distance[i]);
		}		
	}
	//printf("\n min_distance = %g, min_pos = %d\n", min_distance , min_pos);
	
	for (i=0; i<dec_symb->n; i++)
		dec_symb->data[i] = get_int_element(R, min_pos, i);
	
	t.stop();
	*nr_time = t.getElapsedTimeInMicroSec();
	// printf("The time taken by mRRD for a single codeword in usec is %f\n", t.getElapsedTimeInMicroSec());
	// Martzi: debug prints
	//printf("Decoded word is -----\n");
	//print_gf2vec(dec_symb);
	
	//	getchar();
	
	free_int_mtx(R);
	free_int_mtx(S);
	free(idx);
	free(distance);

	free_gf2_vec(prev_dec_symb);
	free_int_mtx(H_mod);
	free_real_vec(s);
	free(ord);
	free_int_vec(PHI);
	delete[] mu1;
	delete[] mu2;
}


void hd_measured_data(real_vec_t* ch_symb, gf2_vec_t *hd_symb) 
{
	int i;
	for (i=0; i<ch_symb->n; i++) { 
		hd_symb->data[i] = (ch_symb->data[i] >=0 ? 1 : 0);
	}
	
}

double_mtx_t * calc_TEST_ber(ecc_t *ecc_p)
{
	// double curr_EBN0, n_errors, n_errors_uc;
	double n_errors, n_errors_uc;
	int erf_n_errors_uc, erf_n_errors, erf_sum_iter, erf_sum_iter_sqr;
	int erf_frame_errors_uc, erf_frame_errors;
	int i, j, iter, length, dec_uc, w;
	// double_mtx_t *ber, *alpha_history;
	double_mtx_t *ber;
	double_mtx_t *fer;
	gf2_vec_t *info, *symb, *dec_symb, *hd_dec_symb;
	real_vec_t *data, *data_uc;
	// int_mtx_t *S; 
	real_vec_t *lambda, *mu;
	unsigned  seed;
	int cum_iter, erf_cum_iter;
  
	double end_EBN0_db, start_EBN0_db, EBN0_delta;
	int N;	
	int_mtx_t *per_list; 
	// int_mtx_t *G, *H, *H_tx ;
	int_mtx_t *G, *H ;
	int total_nr_iter, sum_iter;
	//clock_t tS, tE;
	double erf_sum_time, sum_time, nr_time;
	double erf_sum_tg_time, sum_tg_time, tg_time;
	
	lambda		= (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	mu		= (real_vec_t*)ecc_malloc(sizeof(real_vec_t)); 
	info		= (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	symb		= (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	dec_symb	= (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	hd_dec_symb	= (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	data		= (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	data_uc		= (real_vec_t*)ecc_malloc(sizeof(real_vec_t));

	N		= ecc_p->N;
	end_EBN0_db	= ecc_p->end_EBN0_db;
	start_EBN0_db = ecc_p->start_EBN0_db;
	EBN0_delta	= ecc_p->EBN0_delta;
	G    		= ecc_p->G;
	H    		= ecc_p->H;


	init_real_vec(lambda,  G->n);
	init_real_vec(mu, G->n);
  
	length = floor((end_EBN0_db - start_EBN0_db)/EBN0_delta) + 1;
	ber = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
	init_double_mtx(ber, 3, length);

	fer = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
	init_double_mtx(fer, 3, length);


	//alpha_history = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
	//init_double_mtx(alpha_history, length, N);

	seed = (unsigned)time( NULL );
	//printf("seed = %ud\n", seed);
	//srand( seed);

	init_gf2_vec(info, G->m); 
	init_gf2_vec(dec_symb, G->n);
	init_gf2_vec(hd_dec_symb, G->n);
	init_gf2_vec(symb , G->n);
	init_real_vec(data , G->n);
	init_real_vec(data_uc , G->n);

	printf("G->n = %d, G->m = %d\n", G->n , G->m );
	
	per_list = ecc_p->create_permutation();
	
	// Martzi: debug prints
	//printf("permutation list -----\n");
	//for (i=0; i<N_perm_list; i++)
		//print_int_mtx(per_list+ i);
	//getchar();

	printf("G->n = %d, G->m = %d\n", G->n , G->m );

	for (i=0; i<length; i++) {
		sum_iter = 0; // Total BP iterations performed for all frames transmitted with current SNR
		printf("mRRD, Calculating for SNR = %g, n=%d, k=%d\n", start_EBN0_db + EBN0_delta*i, G->n, G->m);
		put_double_element(ber, 0, i, start_EBN0_db + EBN0_delta*i);
		put_double_element(fer, 0, i, start_EBN0_db + EBN0_delta*i);
           
		n_errors = 0; // Counts the total number of error bits over all the transmitted frames (when no genie is used), for coded transmission
		n_errors_uc = 0; // Counts the total number of error bits over all the transmitted frames, for uncoded transmission
		erf_n_errors_uc = erf_n_errors = 0; // Counts the number of error bits (for both coded and uncoded transmission) from all the frames transmitted since the last write to the ERF file
		erf_frame_errors_uc = erf_frame_errors = 0; // Counts the number of error frames (for both coded and uncoded transmission) from all the frames transmitted since the last write to the ERF file
		erf_sum_iter = erf_sum_iter_sqr = 0; // Counts the number of BP iterations (and the square of this number) performed on all the frames transmitted since the last write to the ERF file
		cum_iter = erf_cum_iter = 0; // Counts number of frames transmitted since the last write to the ERF file
		erf_sum_time = sum_time = 0;
		erf_sum_tg_time = sum_tg_time = 0;
		for (iter=0; iter<N; iter++) {
			if (0 == (iter%100)) printf("mRRD, SNR = %g, iter = %d out of %d, n=%d, k=%d (ave iter %d over %d iterations)\n", start_EBN0_db + EBN0_delta*i, iter, N, G->n, G->m, sum_iter/(cum_iter+1),(cum_iter+1));
			
			init_symbol_data(info);
			// Martzi: debug prints
			//printf("message is -----\n");
			//print_gf2vec(info);

			encode_symbols(info, G, symb); 
			// Martzi: debug prints
			//printf("codeword is -----\n");
			//print_gf2vec(symb);
		
			init_measured_data(symb, data, start_EBN0_db + EBN0_delta*i, ((double)(G->m)/(double)G->n));
			// Martzi: debug prints
			//printf("rx waveform is -----\n");
			//print_rvec(data);
			init_measured_data(info, data_uc, start_EBN0_db + EBN0_delta*i, 1.0);
			hd_measured_data(data, hd_dec_symb);
			if ((!ecc_p->use_genie) || ((ecc_p->use_genie) && (ecc_p->t < diff_gf2_vec(hd_dec_symb, symb)))) { 

				for (w=0; w<lambda->n; w++)  lambda->data[w] = -4*data->data[w];
		
				//			printf("iter = %d,seed = %u\n", iter, seed); 
				TEST_decoder(ecc_p, lambda, dec_symb, mu, per_list, &total_nr_iter, &nr_time, &tg_time);

				sum_iter += total_nr_iter;
				erf_sum_iter += total_nr_iter;
				erf_sum_iter_sqr += total_nr_iter*total_nr_iter;
				erf_sum_time += nr_time;
				erf_sum_tg_time += tg_time;
				// coded data
				
				for (j=0; j<G->n; j++) {					
					n_errors += (symb->data[j] != dec_symb->data[j] ? 1 : 0);
					erf_n_errors += (symb->data[j] != dec_symb->data[j] ? 1 : 0);					
				}
				// Martzi: debug prints
				//printf("Number of bit errors in all simulated codewords is %f\n", n_errors);
				
				for (j=0; j<G->n; j++) {
					if (symb->data[j] != dec_symb->data[j]) { 
						erf_frame_errors ++;
						break;
					}
				}
				// Martzi: debug prints
				//printf("Number of frame errors in all simulated codewords is %d\n", erf_frame_errors);
				// printf("iter = %d, erf_frame_errors = %d\n", iter, erf_frame_errors);

				cum_iter ++; 
				erf_cum_iter ++; 
			
//				printf("iter = %d, total_nr_iter = %d, sum_iter = %d, cum_iter = %d\n", 
//				       iter , total_nr_iter, sum_iter, cum_iter);
//				getchar();
			} 

			// uncoded data
			for (j=0; j<G->m; j++) {
				dec_uc = (data_uc->data[j] > 0) ? 1 : 0;
				n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
				erf_n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
			}

			
			for (j=0; j<G->m; j++) {
				dec_uc = (data_uc->data[j] > 0) ? 1 : 0;
				if (info->data[j] != dec_uc) {
					erf_frame_errors_uc ++;
					break;
				}					
			}
			
			if ((((iter+1)%ERF_QUNT == 0)) && (erf_cum_iter>0)) { 
				update_erf_file_with_iter_and_fer(ecc_p->test_fname,
								  start_EBN0_db + EBN0_delta*i,  
								  ERF_QUNT, 
								  G->n, G->m, 
								  erf_n_errors_uc,  
								  erf_n_errors ,
								  erf_cum_iter, // Frames since last write
								  erf_sum_iter, // BP iters since last write
								  erf_sum_iter_sqr,
								  erf_frame_errors_uc,
					erf_frame_errors,
					erf_sum_time,
					erf_sum_tg_time
								  );			
				erf_n_errors = erf_n_errors_uc = 0;
				erf_sum_iter = erf_sum_iter_sqr = 0;
				erf_cum_iter = 0;
				erf_frame_errors_uc = erf_frame_errors = 0;	  
				erf_sum_time = 0;
				erf_sum_tg_time = 0;
			}

		}
		// Martzi: There is a bug here. If we are using genie, then the number of frames transmitted is not neccesarily N and could be less
		// Hence the denominator should be (num_of_frames_taken_into_consideration * G->n)
		// This is true for element 1 of ber (of the coded transmittion) and not to element 2
		put_double_element(ber, 1, i, n_errors/ (double)(N*G->n));
		put_double_element(ber, 2, i, n_errors_uc/ (double)(N*G->m));

	}
 
	dump_double_mtx(ber,ecc_p->test_fname);

	print_double_mtx(ber);
	
	free_real_vec(lambda);
	free_real_vec(mu);
	free_gf2_vec(info);
	free_gf2_vec(symb);
	free_gf2_vec(dec_symb);
	free_gf2_vec(hd_dec_symb);
	free_real_vec(data);
	free_real_vec(data_uc);

	return ber;

}




static void MBBP_decoder(ecc_t *ecc_p, real_vec_t *lambda, gf2_vec_t *dec_symb, real_vec_t *mu)
{
	int i, j, l, *idx, is_zero, zero_not_found, min_pos, sum;
	double *distance, min_distance;
	// int_mtx_t *H, *R, *S, *Ht;
	int_mtx_t *H, *R, *S;
	double tmp_tg_time;

	int ex_golay_legal_lines[] = {0, 2, 4, 6, 7, 8, 10, 11, 12, 13, 15, 17, 18, 19, 20, 21, 25, 26, 29, 30, 31, 32};

	ecc_assert(ecc_p->mbbp_l <= ecc_p->cgg->m, "ERROR: l > m (%d vs. %d)\n", ecc_p->mbbp_l, ecc_p->cgg->m);
	ecc_assert(ecc_p->mbbp_l != 0, "ERROR: ecc_p->mbbp_l = %d\n", ecc_p->mbbp_l);
	ecc_assert(ecc_p->mbbp_i != 0, "ERROR: ecc_p->mbbp_l = %d\n", ecc_p->mbbp_i);
	
	H = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
	R = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
	S = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
	idx = (int*)ecc_malloc(sizeof(int)*ecc_p->mbbp_l);
	memset(idx, 0, sizeof(int)*ecc_p->mbbp_l);
	distance = (double*)ecc_malloc(sizeof(double)*ecc_p->mbbp_l);

	init_int_mtx(H, ecc_p->cgg->n, ecc_p->cgg->n);
	init_int_mtx(R, ecc_p->mbbp_l, ecc_p->cgg->n);
	init_int_mtx(S, ecc_p->mbbp_l, ecc_p->H->n);


	for (l=0; l<ecc_p->mbbp_l; l++) {

		// create H matrix
		if (ecc_p->mbbp_c_rest)
			// patch for golay
			create_cgg(ecc_p->cgg, ex_golay_legal_lines[l], H, ecc_p->mbbp_c_rest);
		else
			create_cgg(ecc_p->cgg, l, H, ecc_p->mbbp_c_rest);
			
		// tanner_graph_decoder(H, lambda, dec_symb, mu, ecc_p->mbbp_i, 1.0);
		tanner_graph_decoder(H, lambda, dec_symb, mu, ecc_p->mbbp_i, 1.0, &tmp_tg_time);
		for (i=0; i<dec_symb->n; i++) {
			put_int_element(R, l, i, dec_symb->data[i]);
			sum = 0; 
			for (j=0; j<H->n; j++)
				sum += dec_symb->data[j] * get_int_element(H, i, j); // c*H_L'
			put_int_element(S, l, i, sum%2);
		}

		
	}
	
	
	// RECHECK: inovation? works fine for BCH
	//Ht = transpose_int_mtx(ecc_p->H); 
	//multiply_int_mtx(R, Ht, S);

//	printf("R:\n");
//	print_int_mtx(R);
//	printf("S:\n");
//	print_int_mtx(S);
//	getchar();

	zero_not_found = 1;
	for (i=0; i<S->m; i++) {
		is_zero = 0;
		for (j=0; j<S->n; j++)
			if (get_int_element(S, i, j)) {
				is_zero = 1;
				zero_not_found = 0;
			}
		if (!is_zero)
			idx[i] = 1;
	}
	
	if (zero_not_found) 
		for (i=0; i<ecc_p->mbbp_l; i++) {
			idx[i] = 1;
			//printf("%d ", idx[i]);
		}

		
	// finiding the minimum distance (assuming lambda is -4*y)
	min_distance = 1e99;
	for (i=0; i<S->m; i++) {
		if (idx[i]) {
			distance[i] = 0;
			for (j=0; j<S->n; j++)
				distance[i] += sqr((-lambda->data[j]/4) - (2*(double)get_int_element(R, i, j) - 1));
			if (distance[i]<min_distance) {
				min_distance = distance[i];
				min_pos = i;
			}
			//printf("%g ", distance[i]);
		}		
	}
	//printf("\n min_distance = %g, min_pos = %d\n", min_distance , min_pos);
	
	for (i=0; i<dec_symb->n; i++)
		dec_symb->data[i] = get_int_element(R, min_pos, i);
	
	
	//	getchar();
	
	free_int_mtx(H);
	free_int_mtx(R);
	free_int_mtx(S);
	free(idx);
	free(distance);
}		


double_mtx_t * calc_MBBP_ber(ecc_t *ecc_p)
{
	// double curr_EBN0, n_errors, n_errors_uc;	
	double n_errors, n_errors_uc;	
	int erf_n_errors_uc, erf_n_errors;
	int i, j, iter, length, dec_uc, w;
	// double_mtx_t *ber, *alpha_history;
	double_mtx_t *ber;
	gf2_vec_t *info, *symb, *dec_symb;
	real_vec_t *data, *data_uc;
	// int_mtx_t *S; 
	real_vec_t *lambda, *mu;
	unsigned  seed;
  
	double end_EBN0_db, start_EBN0_db, EBN0_delta;
	int N;	

	// int_mtx_t *G, *H, *H_tx ;
	int_mtx_t *G, *H;
	
	lambda   = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	mu	   = (real_vec_t*)ecc_malloc(sizeof(real_vec_t)); 
	info     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	symb     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	dec_symb = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	data     = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	data_uc  = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));

	N		= ecc_p->N;
	end_EBN0_db	= ecc_p->end_EBN0_db;
	start_EBN0_db = ecc_p->start_EBN0_db;
	EBN0_delta	= ecc_p->EBN0_delta;
	G    		= ecc_p->G;
	H    		= ecc_p->H;

	ecc_assert(ecc_p->cgg, "ERROR: cgg wasn't defined\n");


	init_real_vec(lambda,  G->n);
	init_real_vec(mu, G->n);
  
	length = floor((end_EBN0_db - start_EBN0_db)/EBN0_delta) + 1;
	ber = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
	init_double_mtx(ber, 3, length);
	//alpha_history = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
	//init_double_mtx(alpha_history, length, N);

	seed = (unsigned)time( NULL );
	//printf("seed = %ud\n", seed);
	//srand( seed);

	init_gf2_vec(info, G->m); 
	init_gf2_vec(dec_symb, G->n);
	init_gf2_vec(symb , G->n);
	init_real_vec(data , G->n);
	init_real_vec(data_uc , G->n);
	
	for (i=0; i<length; i++) {
		printf("MBBP, Calculating MBBP = %g, n=%d, k=%d\n", start_EBN0_db + EBN0_delta*i, G->n, G->m);
		put_double_element(ber, 0, i, start_EBN0_db + EBN0_delta*i);
           
		n_errors = 0; 
		n_errors_uc = 0;
		erf_n_errors_uc = erf_n_errors = 0;

		for (iter=0; iter<N; iter++) {
			if (0 == (iter%100)) printf("MBBP, BER = %g, iter = %d out of %d, n=%d, k=%d\n", start_EBN0_db + EBN0_delta*i, iter, N, G->n, G->m);
			init_symbol_data(info);
			encode_symbols(info, G, symb); 
	  

			init_measured_data(symb, data, start_EBN0_db + EBN0_delta*i, ((double)(G->m)/(double)G->n));
			init_measured_data(info, data_uc, start_EBN0_db + EBN0_delta*i, 1.0);

			for (w=0; w<lambda->n; w++)  lambda->data[w] = -4*data->data[w];
			
			//printf("iter = %d,seed = %u\n", iter, seed); 
			MBBP_decoder(ecc_p, lambda, dec_symb, mu);
//
//			//print_gf2vec(dec_symb);
//			//getchar();
//		       
//	  
//
			// coded data
			for (j=0; j<G->n; j++) {
				n_errors += (symb->data[j] != dec_symb->data[j] ? 1 : 0);
				erf_n_errors += (symb->data[j] != dec_symb->data[j] ? 1 : 0);
			}
	  
			// uncoded data
			for (j=0; j<G->m; j++) {
				dec_uc = (data_uc->data[j] > 0) ? 1 : 0;
				n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
				erf_n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
			}
			
			if (((iter+1)%ERF_QUNT == 0)) { 
				update_erf_file(ecc_p->mbbp_fname,
						start_EBN0_db + EBN0_delta*i,  
						ERF_QUNT, 
						G->n, G->m, 
						erf_n_errors_uc,  
						erf_n_errors 
						);			
				erf_n_errors = erf_n_errors_uc = 0;
			}



		}
		put_double_element(ber, 1, i, n_errors/ (double)(N*G->n));
		put_double_element(ber, 2, i, n_errors_uc/ (double)(N*G->m));
	}
 
	dump_double_mtx(ber,			ecc_p->mbbp_fname);

	print_double_mtx(ber);

	free_real_vec(lambda);
	free_real_vec(mu);
	free_gf2_vec(info);
	free_gf2_vec(symb);
	free_gf2_vec(dec_symb);
	free_real_vec(data);
	free_real_vec(data_uc);

	return ber;

}


double_mtx_t * calc_RRD_ber(ecc_t *ecc_p)
{
	// double curr_EBN0, n_errors, n_errors_uc;
	double n_errors, n_errors_uc;
	int i, j, iter, length, dec_uc, w;
	double_mtx_t *ber, *alpha_history;
	gf2_vec_t *info, *symb, *dec_symb;
	real_vec_t *data, *data_uc;
	int_mtx_t *S; 
	real_vec_t *lambda, *mu;
	unsigned  seed;
  
	double end_EBN0_db, start_EBN0_db, EBN0_delta;
	int I1, I2, I3, N;
	int erf_n_errors_uc, erf_n_errors, erf_sum_iter, erf_sum_iter_sqr;

	int_mtx_t *G, *H, *H_tx ;
	int total_nr_iter, sum_iter;

	lambda   = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	mu	   = (real_vec_t*)ecc_malloc(sizeof(real_vec_t)); 
	info     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	symb     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	dec_symb = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	data     = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	data_uc  = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));

	I1		= ecc_p->I1;
	I2		= ecc_p->I2;
	I3		= ecc_p->I3;
	N		= ecc_p->N;
	end_EBN0_db	= ecc_p->end_EBN0_db;
	start_EBN0_db = ecc_p->start_EBN0_db;
	EBN0_delta	= ecc_p->EBN0_delta;
	G    		= ecc_p->G;
	H    		= ecc_p->H;
	H_tx 		= ecc_p->H_tx;


	init_real_vec(lambda,  G->n);
	init_real_vec(mu, G->n);
  
	length = floor((end_EBN0_db - start_EBN0_db)/EBN0_delta) + 1;
	ber = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
	init_double_mtx(ber, 3, length);
	alpha_history = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
	init_double_mtx(alpha_history, length, N);

	S = ecc_p->create_permutation();

#if 0
	int_mtx_t *auto_morph_element;
	int_mtx_t * H_t = transpose_int_mtx(ecc_p->H);
	int_mtx_t* dest = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx(dest, 21, 10);
	int_mtx_t* f = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx_from_array(f, 1, 31, gt_perm + 31*4);

	while (1) {
		printf("In here 1\n");
		auto_morph_element = Algorithm_3(S, N_perm_list, ecc_p->perm_get, ecc_p->perm_operate);
		print_int_mtx(auto_morph_element); getchar();
		//apply_perm_int_mtx(ecc_p->G, f->data);
		apply_perm_int_mtx(ecc_p->G, auto_morph_element->data);
 
		multiply_int_mtx(ecc_p->G, H_t, dest);
      
		print_int_mtx(dest); getchar();

	}
#endif

	//seed = 1173547539;
	seed = (unsigned)time( NULL );
	//printf("seed = %ud\n", seed);
	//srand( seed);

	init_gf2_vec(info, G->m); 
	init_gf2_vec(dec_symb, G->n);
	init_gf2_vec(symb , G->n);
	init_real_vec(data , G->n);
	init_real_vec(data_uc , G->n);

  
  	// For every SNR value
	for (i=0; i<length; i++) {
		sum_iter = 0;
		printf("RRD, Calculating SNR = %g, n=%d, k=%d \n", start_EBN0_db + EBN0_delta*i, G->n, G->m);
		put_double_element(ber, 0, i, start_EBN0_db + EBN0_delta*i);
           
		erf_n_errors_uc = erf_n_errors = n_errors_uc = 	n_errors = 0;
		erf_sum_iter = erf_sum_iter_sqr = 0;
		// For every transmitted codeword
		for (iter=0; iter<N; iter++) {
			// Print status if it is time for it (100 codewords were sent)
			if (0 == (iter%100)) printf("RRD, SNR = %g, iter = %d out of %d, n=%d, k=%d (ave iter %d)\n", start_EBN0_db + EBN0_delta*i, iter, N, G->n, G->m,  sum_iter/(iter+1));
			init_symbol_data(info); // Create information bits
			encode_symbols(info, G, symb); // Encoding information bits to a codeword
	  
			//print_gf2vec(info);
			//print_gf2vec(symb);

			// Transmission over AWGN channel (for coded and uncoded bits)
			init_measured_data(symb, data, start_EBN0_db + EBN0_delta*i, ((double)(G->m)/(double)G->n));
			init_measured_data(info, data_uc, start_EBN0_db + EBN0_delta*i, 1.0);

			// Initialize channel LLRs
			for (w=0; w<lambda->n; w++)  lambda->data[w] = -4*data->data[w];

			//	  printf("Original symbols"); print_gf2vec(symb);

			RRD_decoder(ecc_p, lambda, dec_symb, mu, S, &total_nr_iter);
			sum_iter += total_nr_iter;
			erf_sum_iter += total_nr_iter;
			erf_sum_iter_sqr += total_nr_iter*total_nr_iter;

			put_double_element(alpha_history, i, iter, ecc_p->alpha_success);

			//print_gf2vec(dec_symb);
			//getchar();
		       
	  

			// coded data
			for (j=0; j<G->n; j++) {
				n_errors += (symb->data[j] != dec_symb->data[j] ? 1 : 0);
				erf_n_errors += (symb->data[j] != dec_symb->data[j] ? 1 : 0);
			}
	  
			// uncoded data
			for (j=0; j<G->m; j++) {
				dec_uc = (data_uc->data[j] > 0) ? 1 : 0;
				n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
				erf_n_errors_uc += (info->data[j] != dec_uc ? 1 : 0);
			}
			
			if (((iter+1)%ERF_QUNT == 0)) { 
				update_erf_file_with_iter(ecc_p->rrd_fname,
							  start_EBN0_db + EBN0_delta*i,  
							  ERF_QUNT, 
							  G->n, G->m, 
							  erf_n_errors_uc,  
							  erf_n_errors ,
							  ERF_QUNT, 
							  erf_sum_iter,
							  erf_sum_iter_sqr
							  );			
				erf_n_errors = erf_n_errors_uc = 0;
				erf_sum_iter = erf_sum_iter_sqr = 0;
			}

		}
		put_double_element(ber, 1, i, n_errors/ (double)(N*G->n));
		put_double_element(ber, 2, i, n_errors_uc/ (double)(N*G->m));
	}
 
	dump_double_mtx(ber,			ecc_p->rrd_fname);
	dump_double_mtx(alpha_history,	ecc_p->alpha_fname);

	print_double_mtx(ber);

	free_real_vec(lambda);
	free_real_vec(mu);
	free_gf2_vec(info);
	free_gf2_vec(symb);
	free_gf2_vec(dec_symb);
	free_real_vec(data);
	free_real_vec(data_uc);
	
	return ber;

}


typedef double	KEY_T;

#define GT(x, y) ((x) > (y))

#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp

// Snipppet 1 (from ecc.c.bak)


// This function gets as input generating polynon, and retruns a matrix containg all extension field elements
int_mtx_t* create_GF2m(int generating_polynom)
{
	int_mtx_t *p;
	int i, j, m, *shift_reg_o, *shift_reg_n; 

	p = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));

	// finding the msb of the polynom
	m = 31; 
	while (!((1<<m) & generating_polynom)) m--;

	shift_reg_o = (int*)ecc_malloc(sizeof(int)*m);
	shift_reg_n = (int*)ecc_malloc(sizeof(int)*m);
	memset(shift_reg_o, 0, sizeof(int)*m);
	memset(shift_reg_n, 0, sizeof(int)*m);

	init_int_mtx(p, (1<<m), m);

	// initial condition for the shift reg
	j=1; 
	shift_reg_o[0] = 1;
	put_int_element(p, j, 0, 1);
	j++;
  
	// rotation machine
	do {

		shift_reg_n[0] = shift_reg_o[m-1];
		for (i=1; i<m; i++)
			shift_reg_n[i] = shift_reg_o[i-1] ^ (shift_reg_o[m-1]*((generating_polynom & (1<<i)) ? 1 : 0));
    
		for (i=0; i<m; i++) {
			put_int_element(p, j, i, shift_reg_n[i]);
			shift_reg_o[i] = shift_reg_n[i];
		}

		j++;
	} while (j< (1<<m));
 
	return p; 
}

// Generating the standard BCH parity check matrix
// Inputs are length of the code n, error correction capability t, and a matrix containing all elements of GF(n+1)
int_mtx_t* create_BCH_H(int n, int t, int_mtx_t *field_alpha)
{
	int_mtx_t *p = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	int i, j, m, idx;

	init_int_mtx(p, t*field_alpha->n, n);

	for (j=0; j<t; j++)
		for (i=0; i<n; i++) {
			// power of alpha
			idx = ((1+2*j)*i);

			idx = idx%n;

			for (m=0; m<field_alpha->n; m++) {
				put_int_element(p, j*field_alpha->n + m, i, get_int_element(field_alpha, idx+1, m));
			}
		}

	return p;
}

// Generating BCH generator matrix, by assigning k shifts of the generator polynomial given
// Inputs are the generator polynomial 's' and length of the code 'n'
int_mtx_t* create_BCH_G(const char *s, int n)
{
	// Asuming s represents the generating polynom as in Lin & Costello appendix
	int n_minus_k, len, k, i, j;
	int_mtx_t *p;
	char q[2];

	p = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));

	len = strlen(s);
  
	q[0] = s[0]; q[1] = 0;
	n_minus_k = ((atoi(q) & (1<<2)) ? 3 :  ((atoi(q) & (1<<1)) ? 2 : 1)) + (len-1)*3 - 1;
	k = n - n_minus_k;

	init_int_mtx(p, k, n);
  
	// forming the first line
	for (i=0; i<n_minus_k+1; i++) {
		q[0] = s[len -1 - i/3];
		put_int_element(p, 0, i, ((atoi(q) & (1<< (i%3))) ? 1 : 0));
	}
  
	for (j=1; j<k; j++)
		for (i=0; i<n_minus_k+1; i++)
			put_int_element(p, j, j+i, get_int_element(p, 0, i));

	return p;
}

int find_power(int_mtx_t *p) 
{
	int i;
	//printf("In find power\n");
	//print_int_mtx(p);
	for(i=p->n-1; i>=0; i--) {
		if (1 == get_int_element(p, 0, i))
			return i;
	}

	return -1; // Martzi: In case no power has been found
}

void dev_poly(int_mtx_t *dividend, int_mtx_t *genpoly, int_mtx_t *quotient, int_mtx_t *reminder)
{       	
	int div_power, poly_power, diff_power, i;

	for (i=0; i<quotient->n; i++) 
		quotient->data[i] = 0;

	for (i=0; i<reminder->n; i++) 
		reminder->data[i] = 0;

	poly_power = find_power(genpoly);
	div_power = find_power(dividend);


	for (i=0; i<=div_power; i++)
		reminder->data[i] = dividend->data[i];

	while (div_power >= poly_power) {		
		diff_power = div_power - poly_power;
		quotient->data[diff_power] = 1;

		//printf("quotient (B):");
		//print_int_mtx(quotient);
		//
		//printf("genpoly     :");
		//print_int_mtx(genpoly);
		
		
		for (i=0; i<=poly_power; i++)
			reminder->data[i+diff_power] = (reminder->data[i+diff_power] + genpoly->data[i])%2;
		div_power = find_power(reminder);

		//printf("quotient (A):");
		//print_int_mtx(quotient);
		//
		//printf("poly_power = %d, div_power = %d\n", poly_power, div_power);
		//getchar();

	}
}


// This functions creates systematic generator and parity check matrices for a BCH code with length n and generator polynomial represented by 's'
int_mtx_t** create_sys_BCH(const char* s, int n)
{
	// Asuming s represents the generating polynom as in Lin & Costello appendix
	int n_minus_k, len, k, i, j;
	int_mtx_t **p, *g, *h, *b, *genpoly, *quotient, *reminder, *dividend;
	char q[2];
	// 'p' is a pointer to both generator and parity check matrices
	p = (int_mtx_t**)ecc_malloc(sizeof(int_mtx_t*)*2);
	g = p[0] = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	h = p[1] = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));

	genpoly  = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	quotient = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	reminder = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	dividend = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	b        = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));


	len = strlen(s);
  
	q[0] = s[0]; q[1] = 0;
	// Computing n-k according to the generator polynomial's degree
	n_minus_k = ((atoi(q) & (1<<2)) ? 3 :  ((atoi(q) & (1<<1)) ? 2 : 1)) + (len-1)*3 - 1;
	k = n - n_minus_k;

	// Creating assistive matrices
	init_int_mtx(genpoly, 1, n_minus_k+1);	
	init_int_mtx(quotient, 1, n);	
	init_int_mtx(reminder, 1, n);	
	init_int_mtx(dividend, 1, n+1);	
	init_int_mtx(b, k, n-k);
	init_int_mtx(g, k, n);
	init_int_mtx(h, n-k, n);	

	// forming gen poly (highest degree on the right)
	for (i=0; i<n_minus_k+1; i++) {
		//q[0] = s[len -1 - i/3];
		//put_int_element(genpoly, 0, i, ((atoi(q) & (1<< (i%3))) ? 1 : 0));

		q[0] = s[len -1 - i/3];
		put_int_element(genpoly, 0, n_minus_k - i, ((atoi(q) & (1<< (i%3))) ? 1 : 0));

	}
	
	//print_int_mtx(genpoly);
	//printf("n = %d, k = %d\n", n, k);

	// In the following loop, we generate a matrix contains all remainders of the division of all x^w polynomials, n-k+1<=w<=n-1 with the generator polynomial
	// This will be used for systematic encoding of BCH codewords

	// creating b
	for (i=0; i<k; i++) {
		dividend->n = n-k+i +1;
		// Initialize the row to zeros
		for (j=0; j<n; j++) 
			dividend->data[j] = 0;
		
		// Assign '1' on the MSB to create x^(n-k+1+i)=x^w
		put_int_element(dividend, 0, n-k+i, 1);
		//print_int_mtx(dividend);
		// Calculate the division of x^w with the generator polynomial
		dev_poly(dividend, genpoly, quotient, reminder);

		//print_int_mtx(quotient);
		//print_int_mtx(reminder);
		//getchar();

		for (j=0; j<n-k; j++) 
			put_int_element(b, i, j, get_int_element(reminder, 0, j));
	
	}
	//print_int_mtx(b);	
	//getchar();
	
	// Assign the systematic part of the generator matrix (I)
	for (i=0; i<k; i++) 
		put_int_element(g, i, i + n-k, 1);
	// Assign the remainder part of the generator matrix to each row (P)
	for (i=0; i<k; i++) 
		for (j=0; j<n-k; j++) 
			put_int_element(g, i, j, get_int_element(b, i, j));
	//print_int_mtx(g);
	//getchar();

	// Assign P' to the parity check matrix
	for (i=0; i<k; i++)  
		for (j=0; j<n-k; j++) 
			put_int_element(h, j, i+n-k, get_int_element(b, i, j));
	// Assign I to the parity check matrix
	for (i=0; i<n-k; i++) 
		put_int_element(h, i, i, 1);
	
	//print_int_mtx(h);
	//getchar();


	free_int_mtx(genpoly );
	free_int_mtx(quotient);
	free_int_mtx(reminder);
	free_int_mtx(dividend);
	free_int_mtx(b       );
	return p;
}




void stand_alone_tg_decoder()
{
	int_mtx_t H1, H2, G;
	gf2_vec_t *info, *symb, *dec_symb;
	real_vec_t *data, *lambda, *mu;
	double EBN0_dB; 
	int w;
	double tmp_tg_time;

	//srand( (unsigned)time( NULL ) );
 

	EBN0_dB = 20.0;

	info     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	symb     = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	dec_symb = (gf2_vec_t *)ecc_malloc(sizeof(gf2_vec_t));
	data     = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	lambda   = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));
	mu       = (real_vec_t*)ecc_malloc(sizeof(real_vec_t));

	init_real_vec(lambda,  G.n);
	init_real_vec(mu, G.n);
	init_int_mtx_from_array(&G,  4, 8, extended_hamming_844_G); 
	init_int_mtx_from_array(&H1, 4, 8, extended_hamming_844_H1); 
	init_int_mtx_from_array(&H2, 4, 8, extended_hamming_844_H2); 
  

	init_gf2_vec(info, G.m); 
	init_gf2_vec(dec_symb, G.n);
	init_gf2_vec(symb , G.n);
	init_real_vec(data , G.n);


	init_symbol_data(info);
	encode_symbols(info, &G, symb); 
	init_measured_data(symb, data, EBN0_dB, ((double)(G.m))/((double)(G.n)));

	for (w=0; w<lambda->n; w++)  lambda->data[w] = -4*data->data[w];
	// tanner_graph_decoder(&H2, lambda, dec_symb, mu, 50, 1.0);
	tanner_graph_decoder(&H2, lambda, dec_symb, mu, 50, 1.0, &tmp_tg_time);

	//  tanner_graph_decoder(&H2, data, dec_symb, 50, 1.0);

	//  print_gf2vec(symb);
	//  print_gf2vec(dec_symb);
}



void permutation_testing()
{
	// int_mtx_t *S, *G, *tmp, H1_tx, H2_tx, H1, H2; 
	int_mtx_t *S, *G, H1_tx, H2_tx, H1, H2; 
	int i, m,*ord;
  
	init_int_mtx_from_array(&H1_tx, 4, 4, extended_hamming_844_H1_tx); 
	init_int_mtx_from_array(&H2_tx, 4, 4, extended_hamming_844_H2_tx); 
  
	init_int_mtx(&H1, 4, 8);
	init_int_mtx(&H2, 4, 8);
 

	m =3;

	G = create_generator_RM1(m);

	multiply_int_mtx(&H1_tx, G, &H1);
	multiply_int_mtx(&H2_tx, G, &H2);

	S = (int_mtx_t *) create_generating_set_RM1(m, N_perm_list);
  
	ord = (int*)ecc_malloc(sizeof(int)*G->n);

	for (i=0; i<60; i++) 
		Algorithm_3(S, N_perm_list, perm_get, RM1_operate);

	Algorithm_3(S, N_perm_list, perm_get, RM1_operate);
	RM1_get_permutation(G, eye(4), S, ord);
}


int_mtx_t * Algorithm_2(int_mtx_t *Hin)
{
	int_mtx_t *Htag;
	int r1_star, r2_star, g_star, Ng_star, Ng2_star, r1, r2, cnt, g, Ng, Ng2;
	scc_t *scc;
	int* tmp_row;
  
	scc = (scc_t*)ecc_malloc(sizeof(scc_t));
	Htag = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	tmp_row = (int*)ecc_malloc(sizeof(int)*Hin->n);

	init_int_mtx_from_imtx(Htag, Hin);
	r1_star = -1;
	r2_star = -1;
	scc_init(scc, Hin);
	scc_count(scc);
	g_star	= scc->g;
	Ng_star	= scc->Ng;
	Ng2_star	= scc->Ng2;
	scc_free(scc);

	//printf("g_star = %d, Ng_star = %d, Ng2_star = %d\n", g_star, Ng_star, Ng2_star);
	//getchar();

	do {
		printf("g_star = %d, Ng_star = %d, Ng2_star = %d\n", g_star, Ng_star, Ng2_star);
		if (r1_star != r2_star) {
			for (cnt = 0; cnt<Hin->n; cnt++)
				put_int_element(Htag, r2_star, cnt, (get_int_element(Htag, r2_star, cnt) + get_int_element(Htag, r1_star, cnt))%2);
		}
		r1_star = -1;
		r2_star = -1;
		for (r1 = 0; r1 < Hin->m; r1++)
			for (r2 = 0; r2 < Hin->m; r2++)
				if (r1 != r2) {
					for (cnt = 0; cnt<Hin->n; cnt++) {
						tmp_row[cnt] = get_int_element(Htag, r2, cnt);
						put_int_element(Htag, r2, cnt, (get_int_element(Htag, r2, cnt) + get_int_element(Htag, r1, cnt))%2);
					}
					scc_init(scc, Htag);
					scc_count(scc); 
					g	= scc->g;
					Ng	= scc->Ng;
					Ng2	= scc->Ng2;

					//printf("r1 = %d, r2 = %d, Ng = %d, Ng = %d, Ng2 = %d\n", r1, r2, g, Ng, Ng2);
					//getchar();

					if (g > g_star)	{
						g_star	= g;
						r1_star = r1;
						r2_star = r2;
						Ng_star = Ng;
						Ng2_star = Ng2;
					}
					else if ((g == g_star) && (Ng < Ng_star)) {
						r1_star = r1;
						r2_star = r2;
						Ng_star = Ng;
						Ng2_star = Ng2;
					}
					else if ((g == g_star) && (Ng == Ng_star)) {
						if (Ng2 < Ng2_star) {
							r1_star = r1;
							r2_star = r2;
							Ng2_star = Ng2;
						}
					}
					scc_free(scc);
					for (cnt = 0; cnt<Hin->n; cnt++)
						put_int_element(Htag, r2, cnt,  tmp_row[cnt]);
				}

		printf("r1_star = %d, r2_star = %d\n", r1_star, r2_star);
	} while ((-1 != r1_star) || (-1 != r2_star));

	scc_init(scc, Htag);
	scc_count(scc); 
	g	= scc->g;
	Ng	= scc->Ng;
	Ng2	= scc->Ng2;

	printf("final values: g = %d, Ng = %d, Ng2 = %d\n", g, Ng, Ng2);

	return Htag;
}

int_mtx_t * inv_Algorithm_2(int_mtx_t *Hin)
{
	int_mtx_t *Htag;
	int r1_star, r2_star, g_star, Ng_star, Ng2_star, r1, r2, cnt, g, Ng, Ng2;
	scc_t *scc;
	int* tmp_row;
  
	scc = (scc_t*)ecc_malloc(sizeof(scc_t));
	Htag = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	tmp_row = (int*)ecc_malloc(sizeof(int)*Hin->n);

	init_int_mtx_from_imtx(Htag, Hin);
	r1_star = -1;
	r2_star = -1;
	scc_init(scc, Hin);
	scc_count(scc);
	g_star	= scc->g;
	Ng_star	= scc->Ng;
	Ng2_star	= scc->Ng2;
	scc_free(scc);

	//printf("g_star = %d, Ng_star = %d, Ng2_star = %d\n", g_star, Ng_star, Ng2_star);
	//getchar();

	do {
		printf("g_star = %d, Ng_star = %d, Ng2_star = %d\n", g_star, Ng_star, Ng2_star);
		if (r1_star != r2_star) {
			for (cnt = 0; cnt<Hin->n; cnt++)
				put_int_element(Htag, r2_star, cnt, (get_int_element(Htag, r2_star, cnt) + get_int_element(Htag, r1_star, cnt))%2);
		}
		r1_star = -1;
		r2_star = -1;
		for (r1 = 0; r1 < Hin->m; r1++)
			for (r2 = 0; r2 < Hin->m; r2++)
				if (r1 != r2) {
					for (cnt = 0; cnt<Hin->n; cnt++) {
						tmp_row[cnt] = get_int_element(Htag, r2, cnt);
						put_int_element(Htag, r2, cnt, (get_int_element(Htag, r2, cnt) + get_int_element(Htag, r1, cnt))%2);
					}
					scc_init(scc, Htag);
					scc_count(scc); 
					g	= scc->g;
					Ng	= scc->Ng;
					Ng2	= scc->Ng2;

					//printf("r1 = %d, r2 = %d, Ng = %d, Ng = %d, Ng2 = %d\n", r1, r2, g, Ng, Ng2);
					//getchar();

					if (g < g_star)	{
						g_star	= g;
						r1_star = r1;
						r2_star = r2;
						Ng_star = Ng;
						Ng2_star = Ng2;
					}
					else if ((g == g_star) && (Ng > Ng_star)) {
						r1_star = r1;
						r2_star = r2;
						Ng_star = Ng;
						Ng2_star = Ng2;
					}
					else if ((g == g_star) && (Ng == Ng_star)) {
						if (Ng2 < Ng2_star) {
							r1_star = r1;
							r2_star = r2;
							Ng2_star = Ng2;
						}
					}
					scc_free(scc);
					for (cnt = 0; cnt<Hin->n; cnt++)
						put_int_element(Htag, r2, cnt,  tmp_row[cnt]);
				}

		printf("r1_star = %d, r2_star = %d\n", r1_star, r2_star);
	} while ((-1 != r1_star) || (-1 != r2_star));

	scc_init(scc, Htag);
	scc_count(scc); 
	g	= scc->g;
	Ng	= scc->Ng;
	Ng2	= scc->Ng2;

	printf("final values: g = %d, Ng = %d, Ng2 = %d\n", g, Ng, Ng2);

	return Htag;
}




void calc_weights(int_mtx_t *H, int M, char *fname)
{
	int_mtx_t *H_M, *H_M_T, *cw, *dest; 
	int i, j, k, tmp;
	i64 cnt, nr_pseudo_codewords;
	int is_zero;
	double *pseudo_cw, weight, sum, num, den;
	double *w1, *w2;  
	i64 quant = 10, w_cnt, w_size;
	FILE *fp;
	int type =2;
  

	w_size = quant;
	w_cnt = 0;
	w1 = (double*)ecc_malloc(sizeof(double)*w_size);

	H_M = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx(H_M, H->m*M, H->n*M);

	cw = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx(cw, 1, H->n*M);

	dest = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx(dest, 1, H->m*M);

	pseudo_cw = (double*)ecc_malloc(sizeof(double)*H->n);


	// creating the M fold matrix. The PI function I choose is (i+j)%Mwhere i is the col index and j is the row index
	for (j=0; j<H->m; j++)
		for (i=0; i<H->n; i++)
			if (get_int_element(H, j, i))
				for (k=0; k<M; k++)
					put_int_element(H_M, j*M + (i+j+k)%M, i*M + k, 1);

	//  print_int_mtx(H_M);
	//  getchar();
  
	// calc all weights
	H_M_T = transpose_int_mtx(H_M);
	nr_pseudo_codewords =  (((i64)1)<<(H->n*M));
	for (cnt=1; cnt<nr_pseudo_codewords; cnt++) {
		if (0 == (cnt%10000))
			printf("%lld out of %lld, w_cnt = %lld\n", cnt, nr_pseudo_codewords, w_cnt);

		for (i=0; i<H->n*M; i++)
			put_int_element(cw, 0, i, (((((i64)1)<<i) & cnt) ? 1 : 0));

		multiply_int_mtx(cw, H_M_T, dest);
      
		is_zero = 1; 
		for (i=0; i<dest->n; i++)
			if (get_int_element(dest, 0, i))
				is_zero = 0;

		if (is_zero) {
			// calc_weight
			for (i=0; i<H->n; i++) {
				sum =0;
				for (k=0; k<M; k++)
					sum += get_int_element(cw, 0, i*M + k);
	      
				pseudo_cw[i] = sum / M;
			}
			num = den = 0;
			for (i=0; i<H->n; i++) {
				num +=  pseudo_cw[i];
				den +=  pseudo_cw[i]*pseudo_cw[i];
			}
			weight = num*num/den;
	  
			//	  print_int_mtx(cw);
			//	  printf("-=-=-=-=-=-=-=-\n");
			//	  printf("weight = %g\n", weight);
			//	  getchar();

			w1[w_cnt++] = weight;
			if (w_cnt  == w_size) {
	      
				w2 = (double*)ecc_malloc(sizeof(double)*(w_size+quant));
				memcpy(w2, w1, sizeof(double)*w_size);
				w_size += quant;
				free(w1);
				w1 = w2; 
			}
		}
	}

	free_int_mtx(H_M);


	// saving file
	if (!(fp = fopen(fname, "wb")))
		ecc_assert(0, "couldnt open %s for dumping", fname);

	fwrite(&type, sizeof(int), 1, fp);
	tmp = 1;
	fwrite(&tmp, sizeof(int), 1, fp);
	tmp = (int)w_cnt; 
	ecc_assert(tmp == w_cnt, "Overflow in size: %d, %lld\n", tmp, w_cnt);
	fwrite(&tmp, sizeof(int), 1, fp);
	fwrite(w1, sizeof(double), tmp, fp);
	fclose(fp);

	free(w1);
}


char str1[MAX_STR];
char str2[MAX_STR];
char str3[MAX_STR];

void choose(FILE* fp, int total, int start, int num, int t_num, char *tmp, int* indices){
	int  i=0, j;

	if (t_num == num) {
		tmp[num] = 0; 

		// build the inequality
		memset(str1, 0, sizeof(MAX_STR));
		memset(str2, 0, sizeof(MAX_STR));
		memset(str3, 0, sizeof(MAX_STR));
		for (j=0; j<total; j++) {
			sprintf(str3, "%d", j+1);
			//			if (j == tmp[j] - '1') {
			if (strstr(tmp,str3)) {
				sprintf(str3, "+x%d ", indices[j]);
				strcat(str1, str3);
			} else {
				sprintf(str3, "-x%d ", indices[j]);	
				strcat(str1, str3);
			}
		}
		sprintf(str3, "<= %d", t_num -1);	
		strcat(str1, str3);
		fprintf(fp, "%s\n", str1);
		return;
	} else {		
		
		for (i=start+1; i<=total; i++) {
			tmp[num] = i+'0';
			choose(fp, total, i, num+1, t_num, tmp, indices);	
		}
	}
}


void generate_ieq_file(int_mtx_t *H, char* fname){
	FILE *fp;
	int i, j, w;
	char tmp[10];
	int *indices, sum;

	ecc_assert(fp = fopen(fname, "wt"), "can't open file %s for output\n", fname);

	indices = (int*)malloc(sizeof(int)*H->n);


	fprintf(fp, "DIM = %d\n", H->n);
	fprintf(fp, "VALID\n");
	for (i=0; i<H->n; i++) fprintf(fp, "0 ");
	fprintf(fp, "\n\n");
	fprintf(fp, "ELIMINATION_ORDER\n");
	for (i=0; i<H->n; i++) fprintf(fp, "%d ", i+1);
	fprintf(fp, "\n\n");
	fprintf(fp, "INEQUALITIES_SECTION\n");
	for (i=0; i<H->n; i++) fprintf(fp, "x%d <= 1\n", i+1);
	for (i=0; i<H->n; i++) fprintf(fp, "x%d >= 0\n", i+1);

	for (i=0; i<H->m; i++){
		sum = 0;
		for (j=0; j<H->n; j++)
			//sum += x[i*n + j];
			sum += get_int_element(H, i, j);

		
		for (j=0, w=0; j<H->n; j++)
			//if (x[i*n + j]) 
			if (get_int_element(H, i, j))
				indices[w++] = j+1;
		
		for (j=1; j<=sum; j+=2) {
			choose(fp, sum, 0, 0, j, tmp, indices);
		}
	}
	
	fprintf(fp, "END\n");

	fclose(fp);
}



int_mtx_t *build_cgg(int_mtx_t *h, int c_rest)
{
	int_mtx_t *C, *out_mtx;
	int i, j, k, w, z, curr_d, dmin;
	int *idx, n_dmin_words, end_flag, pos, last_pos, diff_flag, new_n_dmin_words;
	int lead_one;

	C        = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));


	// building all code book
	printf("Build CGG: Preparing all codewords\n");
	calc_all_codewords(h, C);
	printf("Build CGG: Finished preparing all codewordss\n");

	idx = (int*)ecc_malloc(sizeof(int)*C->m);

	// finding minimum distance
	dmin = (1<<30); 
	//printf("C [%d, %d]\n", C->m, C->n);
	for (i=1; i<C->m; i++) {
		curr_d = 0;
		for (j=0; j<C->n; j++)
			curr_d += get_int_element(C, i, j);
		idx[i] = curr_d;
		if (curr_d < dmin)
			dmin = curr_d;
	}
	
	// elliminating non-min-distance word
	n_dmin_words = 0;
	for (i=0; i<C->m; i++) {
		if (idx[i] != dmin)
			idx[i] = 0;
		else n_dmin_words++;
	}
	//printf("C[%d, %d] dmin = %d, n_dmin_word = %d\n", C->m, C->n, dmin, n_dmin_words);
	//getchar();
	
	// eliminating cyclic shifts variations
	end_flag = 0;
	last_pos = 1;
	new_n_dmin_words = n_dmin_words;
	for (i=0; i<n_dmin_words && !end_flag; i++){
		// find the first word
		pos = -1;
		for (j=last_pos; j<C->m; j++)
			if (idx[j] != 0) {
				pos = j;
				//printf(" j = %d\n", j);
				break;
			}
		
		if (pos == -1)
			break;
		else
			last_pos = pos+1;

		//printf("pos = %d, j = %d\n", pos, j);
				
		diff_flag = 0; // Martzi: to avoid use of diff_flag uninitialized
		// scanning all words from pos to end
		for (k=pos+1; k<C->m; k++)
			if (idx[k] != 0) {	
				// check if this word is a cyclic shift
				for (w=0; w<(C->n - c_rest); w++) {
					diff_flag = 0;
					for (z=0; z<(C->n - c_rest) && !diff_flag; z++)						
						if (get_int_element(C, pos, z) != get_int_element(C, k, (z + w)%(C->n - c_rest)))
							diff_flag = 1;					
					if (!diff_flag) 
						break;
				}
				if (!diff_flag) {
					/* for debug
					printf("Word @ %d: ",pos); 
					for (w=0; w<C->n; w++)
						printf("%d ", get_int_element(C, pos, w));
					printf("\n");			
					printf("Word @ %d: ",k); 
					for (w=0; w<C->n; w++)
						printf("%d ", get_int_element(C, k, w));
					printf("\n");					
					*/
					new_n_dmin_words --;
					//printf("new_n_dmin_words = %d\n", new_n_dmin_words);
					//getchar();
					idx[k] = 0;
				}
			}
	}
	

	//printf("new_n_dmin_words = %d\n", new_n_dmin_words);
	//getchar();

	// creating matrix with all dmin_words
	out_mtx = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
	init_int_mtx(out_mtx, new_n_dmin_words, C->n);
	pos = 0;
	for (k=0; k<C->m; k++){
		if (idx[k]) {
			lead_one = 0;
			while (get_int_element(C, k, (lead_one + C->n-1)%C->n)) lead_one++;

			//printf("k = %d, idx[k]=%d\n", k, idx[k]);
			for (w=0; w<C->n; w++)
				put_int_element(out_mtx, pos, w, get_int_element(C, k, (w + lead_one)%C->n));	
			pos++;
		}
	}
	
	//printf("in_mtx\n");
	//print_int_mtx(h);
	//printf("out_mtx\n");
	//print_int_mtx(out_mtx);
	//getchar();
	return out_mtx;
}





enum{
	extended_hamming, 
	extended_golay, 
	BCH_15_11_3,
	BCH_31_21_5,
	BCH_31_16_7,
	BCH_31_11_11,
	BCH_63_24_7,
	BCH_63_39_9,
	BCH_63_45_7,
	BCH_63_36_11, 
	BCH_127_113_5,
	BCH_127_106_7,
	BCH_127_64_21,
	
};

enum{
	original_matrix = 0,
	modified_matirx = 1,
};




int main_old(int argc, char *argv[])
{




#if 0

	int tmp[] = {1, 1, 1, 0, 0, 1, 1, 1};

	int tmp2[] = {1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1};

	int_mtx_t H1, H2, H3, H4, *H5;
	init_int_mtx_from_array(&H2, 4, 8, extended_hamming_844_H2); 
	init_int_mtx_from_array(&H1, 4, 8, extended_hamming_844_H1); 
	init_int_mtx_from_array(&H3, 2, 4, tmp); 
	init_int_mtx_from_array(&H4, 3, 7, tmp2); 

	//H5 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
	//calc_all_codewords(&H1, H5);
	//print_int_mtx(H5);
	//getchar();
  

	//print_int_mtx(&H4);
	//calc_weights(&H4, 4, "w1_test.dat");

	//  print_int_mtx(&H2);
	//  generate_ieq_file(&H2, "test_hamm.ieq");
  
	//generate_ieq_file(H5, "test_hamm2.ieq");


	// Kelly matrices
	int ha[] = {1,1,0,0,1,1,0,
		    0,1,1,1,1,0,0,
		    0,0,0,1,1,1,1};
  
	int hb[] = {1, 1, 1, 0, 0, 1, 0, 
		    0, 1, 1, 1, 0, 0, 1, 
		    1, 0, 1, 1, 1, 0, 0, 
		    0, 1, 0, 1, 1, 1, 0, 
		    0, 0, 1, 0, 1, 1, 1, 
		    1, 0, 0, 1, 0, 1, 1, 
		    1, 1, 0, 0, 1, 0, 1};

	int hc[] = {1,0,0,1,0,1,1,
		    0,0,1,0,1,1,1,
		    1,1,1,0,0,0,1,
		    1,0,1,1,1,0,0};


	int ha_asi[] = {1, 1, 0, 0, 1, 0, 1,
			0, 1, 1, 1, 0, 0, 1,
			0, 0, 1, 0, 1, 1, 1};


	int hb_asi[] = {1, 1, 0, 0, 1, 0, 1,
			1, 1, 1, 0, 0, 1, 0,
			0, 1, 1, 1, 0, 0, 1,
			1, 0, 1, 1, 1, 0, 0,
			0, 1, 0, 1, 1, 1, 0,
			0, 0, 1, 0, 1, 1, 1,
			1, 0, 0, 1, 0, 1, 1};

	int hc_asi[] = {1, 1, 0, 0, 1, 0, 1,
			1, 1, 1, 0, 0, 1, 0,	  
			0, 1, 1, 1, 0, 0, 1,
			1, 0, 1, 1, 1, 0, 0};



	int h_1[] = {1, 1, 1, 1, 0, 0, 0, 0,
		     0, 0, 1, 1, 1, 1, 0, 0, 
		     0, 0, 0, 0, 1, 1, 1, 1, 
		     0, 1, 1, 0, 0, 1, 1, 0};


	int h_2[] = {0, 1, 0, 1, 0, 1, 0, 1, 
		     1, 1, 0, 0, 0, 0, 1, 1, 
		     1, 0, 1, 0, 1, 0, 1, 0,
		     0, 1, 0, 1, 1, 0, 1, 0};
  
  

	int_mtx_t Ha, Hb, Hc, Ha_asi, Hb_asi, Hc_asi, H_1, H_2;
	init_int_mtx_from_array(&Ha, 3, 7, ha);
	init_int_mtx_from_array(&Hb, 7, 7, hb);
	init_int_mtx_from_array(&Hc, 4, 7, hc);
	init_int_mtx_from_array(&Ha_asi, 3, 7, ha_asi);
	init_int_mtx_from_array(&Hb_asi, 7, 7, hb_asi);
	init_int_mtx_from_array(&Hc_asi, 4, 7, hc_asi);
	init_int_mtx_from_array(&H_1, 4, 8, h_1);
	init_int_mtx_from_array(&H_2, 4, 8, h_2);

	print_int_mtx(&Ha_asi);
	generate_ieq_file(&Ha_asi, "test_ha_asi.ieq");
	print_int_mtx(&Hb_asi);
	generate_ieq_file(&Hb_asi, "test_hb_asi.ieq");
	print_int_mtx(&Hc_asi);
	generate_ieq_file(&Hc_asi, "test_hc_asi.ieq");

	//  print_int_mtx(&Ha);
	//  generate_ieq_file(&Ha, "test_ha.ieq");
	//  print_int_mtx(&Hb);
	//  generate_ieq_file(&Hb, "test_hb.ieq");
	//  print_int_mtx(&Hc);
	//  generate_ieq_file(&Hc, "test_hc.ieq");


	//  print_int_mtx(&H_1);
	//  generate_ieq_file(&H_1, "test_h_1.ieq");
	//  print_int_mtx(&H_2);
	//  generate_ieq_file(&H_2, "test_h_2.ieq");

	exit(1);
  
#else  
  
	// char bch_st1[100], bch_st2[100], bch_st3[100], bch_st4[100];
	int configuration, mtx_configuration; 
	int_mtx_t H1, H2, H1_tx, H2_tx, *G, golay_Hg, golay_Hg_tag;
	ecc_t ecc;
	// Damping coefficient
	// double alpha_val[] = {0.08, 0.1238, 0.25, .5};

	// Creating pointers to all codes' G and H matrices
	int_mtx_t 
		// *p_128, 
		*h_bch_127_113_5, *g_bch_127_113_5,
		*h_bch_127_64_21, *g_bch_127_64_21, 
		*h_bch_127_106_7, *g_bch_127_106_7; 

	int_mtx_t 
		// *p_64, 
		*h_bch_63_39_9, *g_bch_63_39_9, 
		*h_bch_63_24_7, *g_bch_63_24_7, 
		*h_bch_63_45_7, *g_bch_63_45_7, 
		*h_bch_63_36_11, *g_bch_63_36_11;

	int_mtx_t 
		*p_32, 
		*h_bch_31_21_5, *g_bch_31_21_5, 
		*h_bch_31_16_7, *g_bch_31_16_7,
		*h_bch_31_11_11, *g_bch_31_11_11,
		*tmp_H;
	int_mtx_t 
		*p_16, 
		*h_bch_15_11_3, *g_bch_15_11_3;
	int_mtx_t **pp;

	char hidden_layers_weights_filename [MAX_STR];
	char out_weights_filename [MAX_STR];

	//srand( (unsigned)time( NULL ) );


	// Extented Hamming (8, 4, 4)
	G = create_generator_RM1(3);

	//init_int_mtx_from_array(&G,  4, 8, extended_hamming_844_G); 
	init_int_mtx_from_array(&H2, 4, 8, extended_hamming_844_H2); 
	init_int_mtx_from_array(&H1, 4, 8, extended_hamming_844_H1); 
	init_int_mtx_from_array(&H1_tx, 4, 4, extended_hamming_844_H1_tx); 
	init_int_mtx_from_array(&H2_tx, 4, 4, extended_hamming_844_H2_tx);   
	init_int_mtx_from_array(&golay_Hg, 12, 24, extended_golay_Hg); 
	init_int_mtx_from_array(&golay_Hg_tag, 12, 24, extended_golay_Hg_tag); 


	p_16		= create_GF2m(0x13); // Generating GF(2^4) from x^4 = x + 1
	// Creating generator and parity check matrices of the (15,11,3) BCH code
	h_bch_15_11_3	= create_BCH_H(15, 1, p_16);
	g_bch_15_11_3	= create_BCH_G("23", 15);  // base on lin & costello apendix representation

	//print_int_mtx(h_bch_15_11_3);
	//generate_ieq_file(h_bch_15_11_3, "test_bch_15_11.ieq");
	//exit(1);

		
//	p_128		= create_GF2m(0x89);
//	h_bch_127_113_5	= create_BCH_H(127, 2, p_128);
//	g_bch_127_113_5	= create_BCH_G("41567", 127);  // base on lin & costello apendix representation

	pp = create_sys_BCH("41567", 127);
	g_bch_127_113_5 = pp[0];
	h_bch_127_113_5 = pp[1];

	pp = create_sys_BCH("11554743", 127);
	g_bch_127_106_7 = pp[0];
	h_bch_127_106_7 = pp[1];
	
	pp = create_sys_BCH("1206534025570773100045", 127);
	g_bch_127_64_21 = pp[0];
	h_bch_127_64_21 = pp[1];

//	p_64		= create_GF2m(0x43);
//	h_bch_63_39_9	= create_BCH_H(63, 4, p_64);
//	g_bch_63_39_9	= create_BCH_G("166623567", 63);  // base on lin & costello apendix representation
//	h_bch_63_45_7	= create_BCH_H(63, 3, p_64);
//	g_bch_63_45_7	= create_BCH_G("1701317", 63);
//	h_bch_63_36_11= create_BCH_H(63, 5, p_64);
//	g_bch_63_36_11= create_BCH_G("1033500423", 63);


	pp = create_sys_BCH("17323260404441", 63); 
	g_bch_63_24_7 = pp[0];
	h_bch_63_24_7 = pp[1];


	pp = create_sys_BCH("1701317", 63);
	g_bch_63_45_7 = pp[0];
	h_bch_63_45_7 = pp[1];


	pp = create_sys_BCH("166623567", 63);
	g_bch_63_39_9 = pp[0];
	h_bch_63_39_9 = pp[1];
		

	pp = create_sys_BCH("1033500423", 63);
	g_bch_63_36_11 = pp[0];
	h_bch_63_36_11 = pp[1];

	p_32		= create_GF2m(0x25);
	h_bch_31_21_5	= create_BCH_H(31, 2, p_32);
	g_bch_31_21_5	= create_BCH_G("3551", 31);  // base on lin & costello apendix representatio


	h_bch_31_16_7	= create_BCH_H(31, 3, p_32);
	g_bch_31_16_7	= create_BCH_G("107657", 31);  // base on lin & costello apendix representation


	h_bch_31_11_11= create_BCH_H(31, 5, p_32);
	g_bch_31_11_11= create_BCH_G("5423325", 31);  // base on lin & costello apendix representation

	//configuration  = BCH_15_11_3;
	//configuration  = BCH_31_21_5;
	//configuration  = BCH_31_16_7;
	//configuration  = BCH_31_11_11;
	//configuration		= BCH_63_39_9;
	//configuration		= BCH_63_24_7;
	// configuration		= BCH_63_45_7;
	//configuration		= BCH_63_36_11;
	//configuration		= extended_golay;
	//configuration  = BCH_127_113_5;
	 configuration  = BCH_127_64_21;
	// configuration  = BCH_127_106_7;

	// Determining the parity chcek matrix to be used:
	// original matrix - the created systematic parity check matrix
	// modified matrix - the matrix to be used is the created systematic parity check matrix after applying Halford & Chugg cycle reduction algorithm
	mtx_configuration	= modified_matirx; 
	//mtx_configuration	= 0; 

	ecc.mbbp_c_rest  = 0;

	if (extended_hamming == configuration){
		// Hamming (8, 4, 4)
		ecc.create_permutation	= create_RM1_3_permutation;
		ecc.perm_get			= perm_get;
		ecc.perm_operate		= RM1_operate;
		ecc.G				= G;
		ecc.H				= G; // self dual
		ecc.H_tx_present		= 1;
		ecc.H_tx			= &H2_tx;
		ecc.N				= 100000;
		ecc.I1			= 1; 
		ecc.I2			= 30;
		ecc.I3			= 20;
		ecc.start_EBN0_db		= 4;
		ecc.end_EBN0_db		= 9;
		ecc.EBN0_delta		= 0.5;
		ecc.rrd_alpha			= 0.08;
		ecc.tg_alpha			= 0.75;
		ecc.TG_iter			= 50;
      
		strcpy(ecc.tg_fname,	"H2_tg_ber.dat");
		strcpy(ecc.rrd_fname,	"H2_rrd_ber.dat");
		strcpy(ecc.mbbp_fname,	"H2_mbbp_ber.dat");
		strcpy(ecc.alpha_fname,	"H2_alpha_hist.dat");

	}else if (extended_golay == configuration){
		// Golay (24, 12, 8)
		ecc.create_permutation	= create_ext_golay_permutation;
		ecc.perm_get		= perm_get;
		ecc.perm_operate		= perm_operate;
		//ecc.G				= &golay_Hg;
		//ecc.H				= &golay_Hg; // self dual
		ecc.G				= &golay_Hg_tag;
		ecc.H				= &golay_Hg_tag; // self dual


		ecc.H_tx_present		= 0;
		ecc.N			= 150000;
		//ecc.N			= 10000;
		//ecc.N			= 2000;
		ecc.I1			= 2; 
		ecc.I2			= 30;
		ecc.I3			= 20;
		ecc.start_EBN0_db	= 3.5;
		ecc.end_EBN0_db		= 7;

		ecc.EBN0_delta		= 0.5;
		ecc.rrd_alpha			= 0.08;
		ecc.tg_alpha			= 0.75;
		ecc.TG_iter			= 50;


		ecc.mbbp_i	= 70;
		ecc.mbbp_l	= 20;
		ecc.mbbp_c_rest  = 1;
		ecc.tg_alpha			= 1;

		ecc.cgg = build_cgg(&golay_Hg, 1);
      
      
		sprintf(ecc.tg_fname,		"Hg_exgolay_tg_ber.dat");
		sprintf(ecc.rrd_fname,		"Hg_exgolay_rrd_ber.dat");
		sprintf(ecc.mbbp_fname,		"Hg_exgolay_mbbp_ber_%d_%d.dat", ecc.mbbp_i, ecc.mbbp_l);
		sprintf(ecc.alpha_fname,	"Hg_exgolay_alpha_hist.dat");
		sprintf(ecc.test_fname,		"Hg_exgolay_test_%d_%d_ber.dat", ecc.mbbp_l, ecc.I1*ecc.I2);
		sprintf(ecc.dbd_fname,		"Hg_exgolay_dbd_%d_%d_ber.dat", ecc.mbbp_l, ecc.I1*ecc.I2);
		sprintf(ecc.jiang_fname,	"Hg_exgolay_jiang_%d_ber.dat", ecc.I1*ecc.I2);
		
		// patch for DBD
		
		int_mtx_t tmp_perm;
		//		  int_mtx_t sys_golay_Hg;
		int w;
		//		  init_int_mtx_from_array(&sys_golay_Hg, 12, 24, sys_extended_golay_Hg); 
		init_int_mtx_from_array(&tmp_perm, 4, 24, golay_perm);
		for (w=0; w<tmp_perm.m*tmp_perm.n; w++) {
			tmp_perm.data[w]--;
		}
		ecc.cyclic_permutation = &tmp_perm;
		//
		//		  ecc.G				= &sys_golay_Hg;
		//		  ecc.H				= &sys_golay_Hg; // self dual
		//
		ecc.num_edges_in_H = 0;
		
	} else if ((BCH_31_21_5 == configuration) || 
		   (BCH_31_16_7 == configuration) || 
		   (BCH_31_11_11== configuration)){
		ecc.create_permutation	= create_BCH_31_permutation;
		ecc.perm_get		= perm_get;
		ecc.perm_operate		= perm_operate;
		ecc.I1			= 2; 
		ecc.I2			= 30;
		ecc.I3			= 20;

		switch(configuration) {
		case BCH_31_21_5:		
			ecc.mbbp_i	= 60;
			ecc.mbbp_l	= 20;
			ecc.t = 2;
			ecc.G = g_bch_31_21_5; 
			tmp_H = h_bch_31_21_5; 
			sprintf(hidden_layers_weights_filename, "BCH_31_21_5_hidden_layers_weights.txt");
			sprintf(out_weights_filename, "BCH_31_21_5_out_weights.txt");
			sprintf(ecc.ml_fname,		"H_bch_%d_%d%sml_ber.dat",	31, 21, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");  
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%stg_ber.dat",	31, 21, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.rrd_fname,		"H_bch_%d_%d%srrd_ber.dat",	31, 21, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.mbbp_fname,		"H_bch_%d_%d%smbbp_ber.dat",	31, 21, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
#if !TEST_SC
			sprintf(ecc.test_fname,		"H_bch_%d_%d%stest_%d_%d_ber.dat",	31, 21, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
			sprintf(ecc.dbd_fname,		"H_bch_%d_%d%sdbd_%d_%d_ber.dat",	31, 21, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
#else
			sprintf(ecc.test_fname,		"H_bch_%d_%d%stest_SC_%d_%d_ber.dat",	31, 21, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
			sprintf(ecc.dbd_fname,		"H_bch_%d_%d%sdbd_SC_%d_%d_ber.dat",	31, 21, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
#endif
			sprintf(ecc.jiang_fname,	"H_bch_%d_%d%sjiang_%d_ber.dat",	31, 21, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.I1*ecc.I2);
			sprintf(ecc.alpha_fname,	"H_bch_%d_%d%salpha_hist.dat",	31, 21, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			break;
		case BCH_31_16_7:
			ecc.t = 3;
			ecc.mbbp_i	= 60;
			ecc.mbbp_l	= 10;
			ecc.use_genie	= 0;
			ecc.G = g_bch_31_16_7; 
			tmp_H = h_bch_31_16_7; 
			sprintf(ecc.ml_fname,		"H_bch_%d_%d%sml_ber.dat",	31, 16, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");  
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%stg_ber.dat",	31, 16, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.rrd_fname,		"H_bch_%d_%d%srrd_ber.dat",	31, 16, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.mbbp_fname,		"H_bch_%d_%d%smbbp_ber.dat",	31, 16, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.test_fname,		"H_bch_%d_%d%stest_%d_%d_ber.dat",	31, 16, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
			sprintf(ecc.dbd_fname,		"H_bch_%d_%d%sdbd_%d_%d_ber.dat",	31, 16, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
			sprintf(ecc.jiang_fname,	"H_bch_%d_%d%sjiang_%d_ber.dat",	31, 16, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.I1*ecc.I2);
			sprintf(ecc.alpha_fname,	"H_bch_%d_%d%salpha_hist.dat",	31, 16, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			break;
		case BCH_31_11_11:
			ecc.t = 5;
			ecc.G = g_bch_31_11_11;
			tmp_H = h_bch_31_11_11;  
			sprintf(ecc.ml_fname,		"H_bch_%d_%d%sml_ber.dat",	31, 11, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");  
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%stg_ber.dat",	31, 11, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.rrd_fname,		"H_bch_%d_%d%srrd_ber.dat",	31, 11, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.mbbp_fname,		"H_bch_%d_%d%smbbp_ber.dat",	31, 11, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.alpha_fname,		"H_bch_%d_%d%salpha_hist.dat",	31, 11, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			break;
		}
		
	
		printf("BCH_31_21_5 H:\n");
		print_int_mtx(tmp_H);

		if (mtx_configuration	== modified_matirx)
			tmp_H = Algorithm_2(tmp_H);
	
		printf("BCH_31_21_5 G:\n");
		print_int_mtx(ecc.G);

		printf("BCH_31_21_5 H:\n");
		print_int_mtx(tmp_H);
//		exit(1);

		//ecc.cgg = build_cgg(tmp_H, 0);

		ecc.H			= tmp_H;
		// Initialize all weights matrices and variables
		ecc.num_edges_in_H = calc_num_edges_in_H(&ecc);
		ecc.hidden_layers_weights = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
		init_double_mtx(ecc.hidden_layers_weights, ecc.num_edges_in_H, ecc.num_edges_in_H);
		ecc.out_weights = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
		init_double_mtx(ecc.out_weights, ecc.num_edges_in_H, ecc.H->n);

#ifdef __use_soft_TG__
		create_weight_matrices_from_files(&ecc, hidden_layers_weights_filename, out_weights_filename);
		// To match Eliya's weight matrices representation, we must flip the columns of H and G
		flip_columns_int_mtx(ecc.H);
		flip_columns_int_mtx(ecc.G);
#endif
		// Create accumulated sum matrices to lower complexity
		// Initialization
		ecc.accumulated_rows_sum_matrix = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
		init_int_mtx(ecc.accumulated_rows_sum_matrix, ecc.H->m, ecc.H->n);
		ecc.accumulated_columns_sum_matrix = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
		init_int_mtx(ecc.accumulated_columns_sum_matrix, ecc.H->m, ecc.H->n);
		// Calculation
		create_accumulated_sum_matrices_from_parity_check_matrix(ecc.H, ecc.accumulated_rows_sum_matrix, ecc.accumulated_columns_sum_matrix);

		ecc.H_tx_present		= 0;
		ecc.N			= 60000;
		ecc.start_EBN0_db		= 3;
		ecc.end_EBN0_db		= 5.5;
		ecc.EBN0_delta		= 0.5;
		ecc.rrd_alpha			= 0.08;
		ecc.tg_alpha			= 0.75;
		ecc.TG_iter		= 50;
		ecc.use_genie = 0; // Martzi: Initialization for this variable was missing

     
		// patch for DBD
		int *perm  = create_BCH_permutation(5);
		int w;
		int_mtx_t tmp_perm;
		init_int_mtx_from_array(&tmp_perm, 1, (1<<5)-1, perm);	
		for (w=0; w<tmp_perm.m*tmp_perm.n; w++) {
			tmp_perm.data[w]--;
		}
		ecc.cyclic_permutation = &tmp_perm;


	}else if ((BCH_63_39_9  == configuration) ||
		  (BCH_63_24_7  == configuration) ||
		  (BCH_63_45_7  == configuration) ||
		  (BCH_63_36_11 == configuration) ){
		ecc.create_permutation	= create_BCH_63_permutation;
		ecc.perm_get		= perm_get;
		ecc.perm_operate		= perm_operate;
      
		ecc.I1			= 2; 
		ecc.I2			= 30;
		ecc.I3			= 20;

		switch(configuration) {
		case BCH_63_45_7:
			ecc.mbbp_i	= 60;
			ecc.mbbp_l	= 1;
			ecc.G = g_bch_63_45_7; 
			tmp_H = h_bch_63_45_7;
			sprintf(hidden_layers_weights_filename, "BCH_63_45_7_hidden_layers_weights_for_reduced_matrix_learning_on_5_iters.txt");
			// sprintf(hidden_layers_weights_filename, "BCH_63_45_7_hidden_layers_weights_1s_for_reduced_matrix_learning_on_5_iters.txt");
			sprintf(out_weights_filename, "BCH_63_45_7_out_weights_for_reduced_matrix_learning_on_5_iters.txt");
			// sprintf(out_weights_filename, "BCH_63_45_7_out_weights_1s_for_reduced_matrix_learning_on_5_iters.txt");
			sprintf(ecc.ml_fname,		"H_bch_%d_%d%sml_ber.dat",	63, 45, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");  
#ifdef __use_soft_TG__
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%ssoft_tg_ber.dat",	63, 45, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
#else
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%stg_ber.dat",	63, 45, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
#endif		
			sprintf(ecc.rrd_fname,		"H_bch_%d_%d%srrd_ber.dat",	63, 45, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.alpha_fname,	"H_bch_%d_%d%salpha_hist.dat",	63, 45, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");

			sprintf(ecc.mbbp_fname,		"H_bch_%d_%d%smbbp_ber.dat",	63, 45, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
#ifdef __use_soft_TG__
			sprintf(ecc.test_fname,		"H_bch_%d_%d%ssoft_mrrd_%d_%d_ber.dat",	63, 45, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
#else
			sprintf(ecc.test_fname,		"H_bch_%d_%d%smrrd_%d_%d_ber.dat",	63, 45, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
#endif		
			sprintf(ecc.dbd_fname,		"H_bch_%d_%d%sdbd_%d_%d_ber.dat",	63, 45, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);


			break;

		case BCH_63_24_7:
			ecc.mbbp_i	= 60;
			ecc.mbbp_l	= 20;
			ecc.G = g_bch_63_24_7; 
			tmp_H = h_bch_63_24_7; 
			sprintf(ecc.ml_fname,		"H_bch_%d_%d%sml_ber.dat",	63, 24, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");  
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%stg_ber.dat",	63, 24, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.rrd_fname,		"H_bch_%d_%d%srrd_ber.dat",	63, 24, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.alpha_fname,	"H_bch_%d_%d%salpha_hist.dat",	63, 24, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");

			sprintf(ecc.mbbp_fname,		"H_bch_%d_%d%smbbp_ber.dat",	63, 24, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.test_fname,		"H_bch_%d_%d%stest_%d_%d_ber.dat",	63, 24, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
			sprintf(ecc.dbd_fname,		"H_bch_%d_%d%sdbd_%d_%d_ber.dat",	63, 24, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
			break;


		case BCH_63_36_11:
			ecc.mbbp_i	= 60;
			ecc.mbbp_l	= 1;
			ecc.G = g_bch_63_36_11; 
			tmp_H = h_bch_63_36_11;
			sprintf(hidden_layers_weights_filename, "BCH_63_36_11_hidden_layers_weights_for_reduced_matrix_learning_on_5_iters.txt");
			// sprintf(hidden_layers_weights_filename, "BCH_63_36_11_hidden_layers_weights_1s_for_reduced_matrix_learning_on_5_iters.txt");
			sprintf(out_weights_filename, "BCH_63_36_11_out_weights_for_reduced_matrix_learning_on_5_iters.txt");
			// sprintf(out_weights_filename, "BCH_63_36_11_out_weights_1s_for_reduced_matrix_learning_on_5_iters.txt");
			sprintf(ecc.ml_fname,		"H_bch_%d_%d%sml_ber.dat",	63, 36, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");  
#ifdef __use_soft_TG__
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%ssoft_tg_ber.dat",	63, 36, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
#else
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%stg_ber.dat",	63, 36, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
#endif
			sprintf(ecc.rrd_fname,		"H_bch_%d_%d%srrd_ber.dat",	63, 36, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.alpha_fname,		"H_bch_%d_%d%salpha_hist.dat",	63, 36, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");

			sprintf(ecc.mbbp_fname,		"H_bch_%d_%d%smbbp_ber.dat",	63, 36, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
#ifdef __use_soft_TG__
			sprintf(ecc.test_fname,		"H_bch_%d_%d%ssoft_mrrd_%d_%d_ber.dat",	63, 36, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
#else
			sprintf(ecc.test_fname,		"H_bch_%d_%d%smrrd_%d_%d_ber.dat",	63, 36, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
#endif
			sprintf(ecc.dbd_fname,		"H_bch_%d_%d%sdbd_%d_%d_ber.dat",	63, 36, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);

			break;
		case BCH_63_39_9:
			ecc.mbbp_i	= 60;
			ecc.mbbp_l	= 20;	
			ecc.G = g_bch_63_39_9;
			tmp_H = h_bch_63_39_9;  
			sprintf(ecc.ml_fname,		"H_bch_%d_%d%sml_ber.dat",	63, 39, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");  
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%stg_ber.dat",	63, 39, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.rrd_fname,		"H_bch_%d_%d%srrd_ber.dat",	63, 39, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.alpha_fname,		"H_bch_%d_%d%salpha_hist.dat",	63, 39, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");

			sprintf(ecc.mbbp_fname,		"H_bch_%d_%d%smbbp_ber.dat",	63, 39, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.test_fname,		"H_bch_%d_%d%stest_%d_%d_ber.dat",	63, 39, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
			sprintf(ecc.dbd_fname,		"H_bch_%d_%d%sdbd_%d_%d_ber.dat",	63, 39, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
			break;
		}

		if (mtx_configuration	== modified_matirx)
			tmp_H = Algorithm_2(tmp_H);
		// print_int_mtx(tmp_H);
		//print_int_mtx(ecc.G);
		// exit(1);
		ecc.H			= tmp_H;
		
		// Initialize weights matrices
		ecc.num_edges_in_H = calc_num_edges_in_H(&ecc);
		ecc.hidden_layers_weights = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
		init_double_mtx(ecc.hidden_layers_weights, ecc.num_edges_in_H, ecc.num_edges_in_H);
		ecc.out_weights = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
		init_double_mtx(ecc.out_weights, ecc.num_edges_in_H, ecc.H->n);

		// Creating the matrices, if required
#ifdef __use_soft_TG__
		std::cout << "Creating weights matrices" << std::endl;
		create_weight_matrices_from_files(&ecc, hidden_layers_weights_filename, out_weights_filename);
		// To match Eliya's weight matrices representation, we must flip the columns of H and G
		flip_columns_int_mtx(ecc.H);
		// print_int_mtx(ecc.H);
		flip_columns_int_mtx(ecc.G);
#endif
		// Create accumulated sum matrices to lower complexity
		// Initialization
		ecc.accumulated_rows_sum_matrix = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
		init_int_mtx(ecc.accumulated_rows_sum_matrix, ecc.H->m, ecc.H->n);
		ecc.accumulated_columns_sum_matrix = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
		init_int_mtx(ecc.accumulated_columns_sum_matrix, ecc.H->m, ecc.H->n);
		// Calculation
		create_accumulated_sum_matrices_from_parity_check_matrix(ecc.H, ecc.accumulated_rows_sum_matrix, ecc.accumulated_columns_sum_matrix);
		// print_int_mtx(ecc.accumulated_rows_sum_matrix);
		// print_int_mtx(ecc.accumulated_columns_sum_matrix);

		ecc.H_tx_present		= 0;
	    ecc.N				= 300000;
	    // ecc.N				= 60000;
		// ecc.N			= 6000;
		// ecc.N			= 4;
		ecc.start_EBN0_db		= 1;
		ecc.end_EBN0_db		= 6;
		ecc.EBN0_delta		= 1;
		ecc.rrd_alpha			= 0.08;
		// ecc.tg_alpha			= 0.75;
		ecc.tg_alpha			= 1;
		// ecc.TG_iter		= 100;
		ecc.TG_iter		= 5;
		ecc.use_genie = 0; // Martzi: Initialization for this variable was missing


		// patch for DBD
		int *perm  = create_BCH_permutation(6);
		int w;
		int_mtx_t tmp_perm;
		init_int_mtx_from_array(&tmp_perm, 1, (1<<6)-1, perm);	
		for (w=0; w<tmp_perm.m*tmp_perm.n; w++) {
			tmp_perm.data[w]--;
		}
		ecc.cyclic_permutation = &tmp_perm;
     
	} else if ((BCH_127_113_5 == configuration) ||
				(BCH_127_106_7 == configuration) ||
				(BCH_127_64_21 == configuration)) {

		ecc.create_permutation	= create_BCH_127_permutation;
		ecc.perm_get		= perm_get;
		ecc.perm_operate	= perm_operate;
      
		ecc.I1			= 2; 
		ecc.I2			= 30;
		ecc.I3			= 20;


		switch(configuration) {
		case BCH_127_113_5:
			ecc.mbbp_i	= 60;
			ecc.mbbp_l	= 15;
			ecc.G = g_bch_127_113_5; 
			tmp_H = h_bch_127_113_5; 
			sprintf(ecc.ml_fname,		"H_bch_%d_%d%sml_ber.dat",	127, 113, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");  
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%stg_ber.dat",	127, 113, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.rrd_fname,		"H_bch_%d_%d%srrd_ber.dat",	127, 113, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.alpha_fname,	"H_bch_%d_%d%salpha_hist.dat",	127, 113, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");

			sprintf(ecc.mbbp_fname,		"H_bch_%d_%d%smbbp_ber.dat",	127, 113, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.test_fname,		"H_bch_%d_%d%stest_%d_%d_ber.dat",	127, 113, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
			sprintf(ecc.dbd_fname,		"H_bch_%d_%d%sdbd_%d_%d_ber.dat",	127, 113, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);

			break;
		case BCH_127_64_21:
			ecc.mbbp_i	= 60;
			ecc.mbbp_l	= 3;
			ecc.G = g_bch_127_64_21; 
			tmp_H = h_bch_127_64_21; 
			//sprintf(hidden_layers_weights_filename, "BCH_127_64_21_hidden_layers_weights_for_reduced_matrix_learning_on_5_iters.txt");
			sprintf(hidden_layers_weights_filename, "BCH_127_64_21_hidden_layers_weights_1s_for_reduced_matrix_learning_on_5_iters.txt");
			//sprintf(out_weights_filename, "BCH_127_64_21_out_weights_for_reduced_matrix_learning_on_5_iters.txt");
			sprintf(out_weights_filename, "BCH_127_64_21_out_weights_1s_for_reduced_matrix_learning_on_5_iters.txt");
			sprintf(ecc.ml_fname,		"H_bch_%d_%d%sml_ber.dat",	127, 64, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");  
#ifdef __use_soft_TG__
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%ssoft_tg_ber.dat",	127, 64, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.test_fname,		"H_bch_%d_%d%ssoft_mrrd_%d_%d_ber.dat",	127, 64, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
#else
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%stg_ber.dat",	127, 64, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.test_fname,		"H_bch_%d_%d%smrrd_%d_%d_ber.dat",	127, 64, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
#endif
			sprintf(ecc.rrd_fname,		"H_bch_%d_%d%srrd_ber.dat",	127, 64, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.alpha_fname,	"H_bch_%d_%d%salpha_hist.dat",	127, 64, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");

			sprintf(ecc.mbbp_fname,		"H_bch_%d_%d%smbbp_ber.dat",	127, 64, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.dbd_fname,		"H_bch_%d_%d%sdbd_%d_%d_ber.dat",	127, 64, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);

			break;
		case BCH_127_106_7:
			ecc.mbbp_i	= 60;
			ecc.mbbp_l	= 1;
			ecc.G = g_bch_127_106_7; 
			tmp_H = h_bch_127_106_7; 
			sprintf(ecc.ml_fname,		"H_bch_%d_%d%sml_ber.dat",	127, 106, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");  
			sprintf(ecc.tg_fname,		"H_bch_%d_%d%stg_ber.dat",	127, 106, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.rrd_fname,		"H_bch_%d_%d%srrd_ber.dat",	127, 106, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.alpha_fname,	"H_bch_%d_%d%salpha_hist.dat",	127, 106, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");

			sprintf(ecc.mbbp_fname,		"H_bch_%d_%d%smbbp_ber.dat",	127, 106, (mtx_configuration	== modified_matirx) ? "_mod_" : "_");
			sprintf(ecc.test_fname,		"H_bch_%d_%d%stest_%d_%d_ber.dat",	127, 106, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);
			sprintf(ecc.dbd_fname,		"H_bch_%d_%d%sdbd_%d_%d_ber.dat",	127, 106, (mtx_configuration	== modified_matirx) ? "_mod_" : "_", ecc.mbbp_l, ecc.I1*ecc.I2);

			break;
		
		}

		// printf("BCH_127_64_21 G:\n");
		// print_int_mtx(ecc.G);
		// print_int_mtx(tmp_H);

		if (mtx_configuration	== modified_matirx)
			tmp_H = Algorithm_2(tmp_H);
		
		ecc.H			= tmp_H;
		// Initialize all weights matrices and variables
		ecc.num_edges_in_H = calc_num_edges_in_H(&ecc);
		ecc.hidden_layers_weights = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
		init_double_mtx(ecc.hidden_layers_weights, ecc.num_edges_in_H, ecc.num_edges_in_H);
		ecc.out_weights = (double_mtx_t *)ecc_malloc(sizeof(double_mtx_t));
		init_double_mtx(ecc.out_weights, ecc.num_edges_in_H, ecc.H->n);

// Creating the matrices, if required
#ifdef __use_soft_TG__
		create_weight_matrices_from_files(&ecc, hidden_layers_weights_filename, out_weights_filename);
		// To match Eliya's weight matrices representation, we must flip the columns of H and G
		flip_columns_int_mtx(ecc.H);
		flip_columns_int_mtx(ecc.G);
#endif
		// Create accumulated sum matrices to lower complexity
		// Initialization
		ecc.accumulated_rows_sum_matrix = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
		init_int_mtx(ecc.accumulated_rows_sum_matrix, ecc.H->m, ecc.H->n);
		ecc.accumulated_columns_sum_matrix = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
		init_int_mtx(ecc.accumulated_columns_sum_matrix, ecc.H->m, ecc.H->n);
		// Calculation
		create_accumulated_sum_matrices_from_parity_check_matrix(ecc.H, ecc.accumulated_rows_sum_matrix, ecc.accumulated_columns_sum_matrix);

		ecc.H_tx_present		= 0;
	        ecc.N			= 60000;
		//ecc.N			= 50;
		ecc.start_EBN0_db		= 3;
		ecc.end_EBN0_db		= 5;
		ecc.EBN0_delta		= 0.5;
		ecc.rrd_alpha			= 0.08;
		ecc.tg_alpha			= 1;
		// ecc.tg_alpha			= 0.75;
		ecc.TG_iter		= 5;
		ecc.use_genie = 0; // Martzi: Initialization for this variable was missing

		// patch for DBD
		int *perm  = create_BCH_permutation(7);
		int w;
		int_mtx_t tmp_perm;
		init_int_mtx_from_array(&tmp_perm, 1, (1<<7)-1, perm);	
		for (w=0; w<tmp_perm.m*tmp_perm.n; w++) {
			tmp_perm.data[w]--;
		}
		ecc.cyclic_permutation = &tmp_perm;

	} else if (BCH_15_11_3 == configuration) {
		
		ecc.G = g_bch_15_11_3; 
		tmp_H = h_bch_15_11_3; 	

		if (mtx_configuration	== modified_matirx)
			tmp_H = Algorithm_2(tmp_H);
	}

	srand( (unsigned)time( NULL ) );
	
	//mtrace(); // for detecting memory leaks


	//stand_alone_tg_decoder();
	//permutation_testing();
	//calc_TG_ber(&ecc);
	//calc_RRD_ber(&ecc);
	//calc_ML_ber(&ecc);
	//calc_MBBP_ber(&ecc);
	calc_TEST_ber(&ecc);
	//calc_DBD_ber(&ecc);
	//calc_Jiang_ber(&ecc);


	return 1;

#endif
}
