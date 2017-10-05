#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ecc.h"

#define sqr(x) ((x)*(x))

#define TEST_SC		0

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
static __inline double randn()
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


static __inline void randn2(double *y1, double *y2)
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
  

	ecc_assert(fmax(2*nr_perm+1, 10) < perm_list_length, "Error with perm_list_length, %d vs. %d\n", fmax(2*nr_perm+1, 10), perm_list_length);
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
		if (in_vec->data[i] != in_vec->data[i])
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


void hd_measured_data(real_vec_t* ch_symb, gf2_vec_t *hd_symb) 
{
	int i;
	for (i=0; i<ch_symb->n; i++) { 
		hd_symb->data[i] = (ch_symb->data[i] >=0 ? 1 : 0);
	}
	
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
