/***
Last Modification: 14 October 2016: setting psc threshold to -1 will print values
Last Modification: 23 Aug 2016
*
*
*
***/

/*****
 * Check if(WHD<minWHD)
 * Replace by if((WHD-minWHD<epsilon) && (fabs(WHD-minWHD)<epsilon))
 * But there something goes wrong
 * 
 * ****/


#define STEP_COUNT_VERBOSE 0
#define ENABLE_PRESELECT 0

#include "mathcom.h"
#include "mathcod.h"


static int LeftSystematicMatrix();
static int InitOSDLeftSysMatrix(int M, int N);
static int OSDLeftSysMatrix();
extern int runOSD();

static double noisevariance, noisedeviation;
static long long int CnkFunctionLL(int N, int K);


//---------------------------------------------------------------
void qsort_r();
/* qsort_r comparison function for use with indices */
static int cmp_decreasing_order(int* index_x,int* index_y,double* values)
{
	if(values[*index_x] > values[*index_y]) return(-1);
	if(values[*index_x] < values[*index_y]) return(1);
	return(0);
}
//---------------------------------------------------------------

//---------------------------------------------------------------
typedef struct
{
	int a, b;
} PermTable;
//---------------------------------------------------------------

//---------------------------------------------------------------
		
// Definition linked list
typedef struct q_ll_entry q_ll_entry;

struct q_ll_entry
{
	Pint q_vector;
	q_ll_entry * next;
	int ind_1, ind_2;
};

typedef struct
{
	q_ll_entry *first, *last;
} q_ll;
//---------------------------------------------------------------




//---------------------------------------------------------------

typedef struct
{
	double param_lambda;
	int psc_threshold;
	int flag_enable_prepr_2;
	int nr_q_1_tau_entries;
	int max_osd_order;
	int param_tau;
	int flag_enable_prob_skip;
	int param_t;
	int param_theta;
	long long int max_STEPS;
} OSD_param_struct;

typedef struct
{
	Pint decoded_info_word;
	Pint order_0_word;
	Pint diff_par_order_0_hard_dec;
	Pint diff_par_order_0_hard_dec_e_2_sub;
	Pint decoded_word;
	q_ll *q_1_tau_array;
	q_ll_entry *q_ll_entries;
	GF2Matrix QMat;
	Pdouble conf_dif;
	GF2Matrix RedundancyBoxes;
	Pdouble WHDBoxes;
	Pint I_array;
	Pint e_2_sub;
	Pint e_2;
} OSD_memory_struct;



//---------------------------------------------------------------



static double Rate;
static int n, k;

int main(int argc, char **argv)
{
	int i, j, u;


	int osd_order_full, count_full_OSD;

	Pint sort_indices;
	Pint codeword_sorted, codeword_sorted_2; // used to sort the codeword (no in place sorting)
	Pdouble y_sorted, y_sorted_2, confidence_sorted, confidence_sorted_2;
	GF2Matrix GsysLeft_sorted;

	int countperms, newk;
	PermTable* perms;
	int temp; // used to apply the permutations


	int osd_order, max_osd_order;
	Pint decoded_word, decoded_info_word, osd_error, BestCodeword, BestCodeword_sorted;
	double minWHD, WHD;
	
	Pint  e_2_sub, e_2, order_0_word, diff_par_order_0_hard_dec, diff_par_order_0_hard_dec_e_2_sub;
	

	char chain[200], infile[200], datfile[200], code_name[100];
	FILE *resfptr,*paramfptr;
	int iter, maxiter, max_word_error, PrintFlag, PartialFlag, PartialCount, seed1, seed2;
	int errcnt_bit, errcnt_word, errcnt_ML;
	double snrdb1, snrdb2, snrstep, snrdb, snr, timer;
	ObjectId prng_id1, prng_id2, prng_id3;

	Pint infobits, codeword;
	Pdouble y, confidence;


	GF2Matrix GsysLeft, GenMatLeft;
	
	int flag_enable_prepr_2;
	int max_number_q_ll_entries, nr_q_1_tau_entries;
	int index_q_1_tau_array,index_q_1_tau_array_bis;
	q_ll *q_1_tau_array;
	q_ll_entry *q_ll_entries, *ll_current_pointer;	


	int flag_enable_prepr_1_full, param_t_full, param_theta_full;
	int flag_enable_prepr_2_full, param_tau_full;
	
	
	
	int flag_enable_prob_skip;
	double param_lambda;
	double prob_skip_thresh;
	

	Pint I_array;

	
	GF2Matrix RedundancyBoxes;
	Pdouble WHDBoxes;
	
	GF2Matrix QMat; 
	Pdouble conf_dif;
	
	long long int max_STEPS;
	
	Pint random_bias;
	Pdouble confidence_biased;
	int bias_iter, max_bias_iter, ret_val;
	double param_bias_magn;
	long long int *STEPs_processed;
	
	Pint sort_indices_interm, sort_indices_sorted;
	Pdouble confidence_sorted_interm;
	
	
	int psc_threshold, psc_s, psc_hw;
	
	
	OSD_param_struct OSD_params;
	OSD_memory_struct OSD_memory;	

	printf("#### 2-B OSD Decoding with preselect, skipping and split of a linear binary [n,k] code on AWGN or fading channel with BPSK input ####\n\n");



	if(argc==1) // Manual entry of parameters
	{
		printf("Manual parameter entry.\n Command line options: \"parameter file\" seed1 seed2\n\n");

		/******************************* Code parameters ******************************************************/
		printf("\e[1;4mCode parameters:\e[0m\n");

		n=128;
		printf("Code length n [%d]: ", n);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%d", &n);
		if((n<8)||(n>512))
		{
			fprintf(stderr,"%s: n out of range.\n",argv[0]);
			return(-1);
		}

		k=64;
		printf("Code dimension k [%d]: ", k);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%d", &k);
		if((k<1)||(k>=n))
		{
			fprintf(stderr,"%s: k out of range.\n",argv[0]);
			return(-1);
		}
		
		Rate=k/(double)n;
		printf("Coding rate is %d/%d=%1.4f \n", k, n, Rate);


		strcpy(code_name,"BCH");
		printf("Give a name to this code [%s]: ", code_name);
		ReadString(chain);
		if(chain[0]) strcpy(code_name, chain);


		GenMatLeft=AllocGF2Mat(k, n);
		GsysLeft=AllocGF2Mat(k, n);
		GsysLeft_sorted=AllocGF2Mat(k,n);


		strcpy(infile,"G_bch_n=128_k=64_dHmin=22_t=10_ext.dat");
		printf("Our convention is kxn for the generator matrix with identity on the left.\n");
		printf("If the input matrix G is not systematic, it will be converted into systematic form.\n");
		printf("Information bits are those associated to the identity (the first k columns of G)\n");
		printf("Input file name for reading the %dx%d generator matrix G [%s]: ", k, n, infile);
		ReadString(chain);
		if(chain[0]) strcpy(infile, chain);
		if(ReadGF2Matrix(GenMatLeft, k, n, infile, SAME)==-1)
		{
			fprintf(stderr,"%s: ReadGF2Matrix() failed.\n",argv[0]);
			return(-1);
		}

		u=TRUE;
		for(i=0; i<k; i++)
			for(j=0; j<k; j++) if(GenMatLeft[i][j]!= (i==j)) {
					u=FALSE;
					printf("Not systematic!\n");
					i=k;
					break;
				}
		if(!u)// convert the input matrix into systematic form
		{
			printf("LeftSystematicMatrix() applied on GenMatLeft...\n");
			if(LeftSystematicMatrix(GenMatLeft, k, n, &newk, &perms, &countperms)==-1)
			{
				fprintf(stderr,"%s: LeftSystematicMatrix() failed.\n", argv[0]);
				return(-1);
			}
			printf("Number of permutations is %d (not saved and not used later)\n", countperms);
			if(newk != k)
			{
				fprintf(stderr,"%s: input matrix does not have rank k=%d!\n", argv[0], k);
				return(-1);
			}
		}

		/********************************** End code parameters ************************************************/




		/*********************************** OSD parameters ***************************************************/
		printf("\n\e[1;4mOSD parameters:\e[0m\n");
	

		osd_order_full=4;
		printf("\e[4mFull OSD\e[0m order [%d]: ", osd_order_full);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%d", &osd_order_full);
		if((osd_order_full<2)||(osd_order_full>k-1))
		{
			fprintf(stderr,"%s: Full OSD order out of range.\n",argv[0]);
			return(-1);
		}
		
		flag_enable_prepr_1_full=1;
		printf("Enable preprocessing rule 1? (y/n) [%c] : ",flag_enable_prepr_1_full?'y':'n');
		ReadString(chain);
		if(chain[0] == (flag_enable_prepr_1_full?'n':'y')) flag_enable_prepr_1_full = !flag_enable_prepr_1_full;
		
		if(flag_enable_prepr_1_full)
		{
			param_t_full=32;
			printf("First preprocessing rule, parameter t [%d]: ", param_t_full);
			ReadString(chain);
			if(chain[0]) sscanf(chain,"%d", &param_t_full);
			
			param_theta_full=14;
			printf("First preprocessing rule, parameter theta [%d]: ", param_theta_full);
			ReadString(chain);
			if(chain[0]) sscanf(chain,"%d", &param_theta_full);
		}
		else
			param_theta_full = param_t_full = 1;
		
		flag_enable_prepr_2_full=0;
		printf("Enable preprocessing rule 2? (y/n) [%c] : ",flag_enable_prepr_2_full?'y':'n');
		ReadString(chain);
		if(chain[0] == (flag_enable_prepr_2_full?'n':'y')) flag_enable_prepr_2_full = !flag_enable_prepr_2_full;
		
		if(flag_enable_prepr_2_full)
		{
			param_tau_full=19;
			printf("Second preprocessing rule, parameter tau [%d]: ", param_tau_full);
			ReadString(chain);
			if(chain[0]) sscanf(chain,"%d", &param_tau_full);
		}

				
		max_bias_iter = 20;
		printf("Number of bias iterations [%d]: ", max_bias_iter);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%d", &max_bias_iter);
		if(max_bias_iter<1)
		{
			fprintf(stderr,"%s: max_bias_iter out of range.\n",argv[0]);
			return(-1);
		}
		
		param_bias_magn = 0.2;
		printf("Bias magnitude [%1.2f]: ", param_bias_magn);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%lf", &param_bias_magn);


		psc_threshold=12;
		printf("Probabilistic sufficient condition threshold [%d]: ", psc_threshold);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%d", &psc_threshold);

		flag_enable_prob_skip=1;
		printf("Enable probabilistic skipping rule? (y/n) [%c] : ",flag_enable_prob_skip?'y':'n');
		ReadString(chain);
		if(chain[0] == (flag_enable_prob_skip?'n':'y')) flag_enable_prob_skip = !flag_enable_prob_skip;
		if(flag_enable_prob_skip)
		{
			param_lambda=2.5;
			printf("Lambda [%1.2f]: ", param_lambda);
			ReadString(chain);
			if(chain[0]) sscanf(chain,"%lf", &param_lambda);
		}
		
		
		
		
		/********************************************************************************************************************************************************************************/
		
		max_osd_order = osd_order_full;
		max_STEPS=0;
		for(osd_order=0; osd_order<=max_osd_order; osd_order++) max_STEPS += CnkFunctionLL(k, osd_order);
		printf("Max number of generated STEPs [%lld]: ", max_STEPS);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%lld", &max_STEPS);
		if(max_STEPS<1)
		{
			fprintf(stderr,"%s: Max number of generated STEPs out of range.\n",argv[0]);
			return(-1);
		}
		/********************************************************************************************************************************************************************************/
		
	
		
		/******************************** End OSD parameters **************************************************/

		/******************************* Simulation parameters ************************************************/
		printf("\n\e[1;4mSimulation parameters:\e[0m\n");

		maxiter=100000000;
		printf("Max number of codewords to be transmitted [%d]: ", maxiter);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%d", &maxiter);


		snrdb1=1.5;
		printf("Start Eb/N0 in dB [%1.2f]: ", snrdb1);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%lf", &snrdb1);

		snrdb2=snrdb1;
		printf("End Eb/N0 in dB [%1.2f]: ", snrdb2);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%lf", &snrdb2);

		snrstep=0.5;
		printf("SNR step in dB [%1.2f]: ", snrstep);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%lf", &snrstep);

		seed1=13;
		printf("Seed for random generation of information bits [%d]: ", seed1);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%d", &seed1);
		prng_id1=InitPRNGMT(seed1);

		seed2=5;
		printf("Seed for random generation of AWG noise [%d]: ", seed2);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%d", &seed2);
		prng_id2=InitPRNGMT(seed2);
		prng_id3=InitPRNGMT(seed2); // TODO: own seed

		max_word_error=200;
		printf("Maximum number of word errors [%d]: ",  max_word_error);
		ReadString(chain);
		if(chain[0]) sscanf(chain,"%d", &max_word_error);


		sprintf(datfile,"perf_%s_%d_%d_awgn_rand_bias_OSD_preselect_skip_order_%d_seeds_%d_%d.dat", code_name, n, k, osd_order_full, seed1, seed2);
		resfptr=fopen(datfile,"w");
		printf("Output will be written into file %s \n", datfile);



		PartialFlag=TRUE;
		printf("Write partial results in output file y/n? [%c]: ",PartialFlag?'y':'n');
		ReadString(chain);
		if(chain[0]==(PartialFlag?'n':'y')) PartialFlag = !PartialFlag;
		if(PartialFlag)
		{
			PartialCount=10;
			printf("Period for writing partial results [%d]: ", PartialCount);
			ReadString(chain);
			if(chain[0]) sscanf(chain,"%d", &PartialCount);
		}


		PrintFlag=TRUE;
		printf("Display output to screen y/n? [%c]: ",PrintFlag?'y':'n');
		ReadString(chain);
		if(chain[0]==(PrintFlag?'n':'y')) PrintFlag = !PrintFlag;

		/**************************** End Simulation parameters *******************************************/
	}
	else if(argc==4)
	{
		printf("Reading the parameter file\n\n");
		u = 0;
		paramfptr = fopen(argv[1],"rt");
		u+=fscanf(paramfptr,"n=%d\n",&n);
		u+=fscanf(paramfptr,"k=%d\n",&k);
		u+=fscanf(paramfptr,"name=\"%[^\"]\"\n",code_name);
		u+=fscanf(paramfptr,"G=\"%[^\"]\"\n",infile);

		
		u+=fscanf(paramfptr,"osd_order_full=%d\n",&osd_order_full);
		u+=fscanf(paramfptr,"flag_enable_prepr_1_full=%d\n",&flag_enable_prepr_1_full);
		u+=fscanf(paramfptr,"param_t_full=%d\n",&param_t_full);
		u+=fscanf(paramfptr,"param_theta_full=%d\n",&param_theta_full);
		u+=fscanf(paramfptr,"flag_enable_prepr_2_full=%d\n",&flag_enable_prepr_2_full);
		u+=fscanf(paramfptr,"param_tau_full=%d\n",&param_tau_full);
		u+=fscanf(paramfptr,"max_bias_iter=%d\n",&max_bias_iter);
		u+=fscanf(paramfptr,"param_bias_magn=%lf\n",&param_bias_magn);	
		u+=fscanf(paramfptr,"psc_threshold=%d\n",&psc_threshold);	
		
		u+=fscanf(paramfptr,"flag_enable_prob_skip=%d\n",&flag_enable_prob_skip);
		u+=fscanf(paramfptr,"param_lambda=%lf\n",&param_lambda);
		u+=fscanf(paramfptr,"max_STEPS=%lld\n",&max_STEPS);
		
		u+=fscanf(paramfptr,"maxiter=%d\n",&maxiter);
		
		u+=fscanf(paramfptr,"snrdb1=%lf\n",&snrdb1);
		u+=fscanf(paramfptr,"snrdb2=%lf\n",&snrdb2);
		u+=fscanf(paramfptr,"snrstep=%lf\n",&snrstep);
		u+=fscanf(paramfptr,"max_word_error=%d\n",&max_word_error);
		u+=fscanf(paramfptr,"PartialFlag=%d\n",&PartialFlag);
		u+=fscanf(paramfptr,"PartialCount=%d\n",&PartialCount);
		u+=fscanf(paramfptr,"PrintFlag=%d\n",&PrintFlag);
		fclose(paramfptr);

		if(u!= 24) {
			fprintf(stderr, "Error while reading the parameter file\n");
			exit(EXIT_FAILURE);
		}

		
		
		
		
		/******************************* Code parameters ******************************************************/
		printf("\e[1;4mCode parameters:\e[0m\n");

		printf("Code length n = %d\n", n);
		if((n<8)||(n>512))
		{
			fprintf(stderr,"%s: n out of range.\n",argv[0]);
			return(-1);
		}

		printf("Code dimension k = %d\n", k);
		if((k<1)||(k>=n))
		{
			fprintf(stderr,"%s: k out of range.\n",argv[0]);
			return(-1);
		}
		
		Rate=k/(double)n;
		printf("Coding rate is %d/%d=%1.4f \n", k, n, Rate);


		printf("Code name = %s\n", code_name);


		GenMatLeft=AllocGF2Mat(k, n);
		GsysLeft=AllocGF2Mat(k, n);
		GsysLeft_sorted=AllocGF2Mat(k,n);


		printf("Input file name for reading the generator matrix = %s\n", infile);
		if(ReadGF2Matrix(GenMatLeft, k, n, infile, SAME)==-1)
		{
			fprintf(stderr,"%s: ReadGF2Matrix() failed.\n",argv[0]);
			return(-1);
		}

		u=TRUE;
		for(i=0; i<k; i++)
			for(j=0; j<k; j++) if(GenMatLeft[i][j]!= (i==j)) {
					u=FALSE;
					printf("Not systematic!\n");
					i=k;
					break;
				}
		if(!u)// convert the input matrix into systematic form
		{
			printf("LeftSystematicMatrix() applied on GenMatLeft...\n");
			if(LeftSystematicMatrix(GenMatLeft, k, n, &newk, &perms, &countperms)==-1)
			{
				fprintf(stderr,"%s: LeftSystematicMatrix() failed.\n", argv[0]);
				return(-1);
			}
			printf("Number of permutations is %d (not saved and not used later)\n", countperms);
			if(newk != k)
			{
				fprintf(stderr,"%s: input matrix does not have rank k=%d!\n", argv[0], k);
				return(-1);
			}
		}

		/********************************** End code parameters ************************************************/




		printf("\e[4mFull OSD\e[0m order = %d\n", osd_order_full);
		if(osd_order_full<1)
		{
			fprintf(stderr,"%s: full OSD order out of range.\n",argv[0]);
			return(-1);
		}
		
		printf("Enable preprocessing rule 1 = %c\n",flag_enable_prepr_1_full?'y':'n');
		
		if(flag_enable_prepr_1_full)
		{
			printf("First preprocessing rule, parameter t = %d\n", param_t_full);
			
			printf("First preprocessing rule, parameter theta = %d\n", param_theta_full);
		}
		else
			param_theta_full = param_t_full = 1;
		
		printf("Enable preprocessing rule 2 = %c\n",flag_enable_prepr_2_full?'y':'n');
		
		if(flag_enable_prepr_2_full)
			printf("Second preprocessing rule, parameter tau  =%d\n", param_tau_full);

				
		printf("Max bias iterations = %d\n",max_bias_iter);

		printf("Bias magnitude = %1.2f\n", param_bias_magn);
		
		printf("Probabilistic sufficient condition threshold = %d\n", psc_threshold);
		
		printf("Enable probabilistic skipping rule = %c\n",flag_enable_prob_skip?'y':'n');
		if(flag_enable_prob_skip)
			printf("Lambda = %1.2f\n", param_lambda);


		max_osd_order = osd_order_full;
		if(max_STEPS==0) for(osd_order=0; osd_order<=max_osd_order; osd_order++) max_STEPS += CnkFunctionLL(k, osd_order);
		printf("Max number of generated STEPs = %lld\n", max_STEPS);
	
		
		/******************************** End OSD parameters **************************************************/

		/******************************* Simulation parameters ************************************************/
		printf("\n\e[1;4mSimulation parameters:\e[0m\n");

		printf("Max number of codewords to be transmitted = %d\n", maxiter);
		

		printf("Start Eb/N0 = %1.2f dB\n", snrdb1);

		printf("End Eb/N0 = %1.2f dB\n", snrdb2);

		printf("SNR step = %1.2f dB\n", snrstep);

		seed1=atoi(argv[2]);
		printf("Seed for random generation of information bits = %d\n", seed1);
		prng_id1=InitPRNGMT(seed1);

		seed2=atoi(argv[3]);
		printf("Seed for random generation of AWG noise = %d\n", seed2);
		prng_id2=InitPRNGMT(seed2);
		prng_id3=InitPRNGMT(seed2); // TODO: own seed



		printf("Maximum number of word errors = %d\n",  max_word_error);
		
		// Ethan adding the directory name to datfile
		char *last_slash_in_path = NULL;
		char params_file_dir[256] = {0};
		last_slash_in_path = strrchr(argv[1], '/');
		if (last_slash_in_path != NULL) {
			int dir_name_len = last_slash_in_path - argv[1] + 1;
			memcpy(params_file_dir, argv[1], dir_name_len);
			params_file_dir[dir_name_len] = 0;
		}

		sprintf(datfile,"%sperf_%s_%d_%d_awgn_%d_seeds_%d_%d.dat", params_file_dir, code_name, n, k, osd_order_full, seed1, seed2);
		
		while (fopen(datfile,"r") != NULL) {
			strcat(datfile, "1");
		}

		resfptr=fopen(datfile,"w");
		printf("Output will be written into file %s \n", datfile);
		
	}
	else
	{
		fprintf(stderr, "Wrong number of command line options!\nUse no command line options for manual entry or use \"parameter file\" seed1 seed2\n\n");
		exit(EXIT_FAILURE);
	}








	/**************************** Memory allocation *******************************************/
	



	infobits=(Pint) calloc(k, sizeof(int));
	codeword=(Pint) calloc(n, sizeof(int));
	codeword_sorted=(Pint) calloc(n, sizeof(int));
	codeword_sorted_2=(Pint) calloc(n, sizeof(int));
	y=(Pdouble) calloc(n, sizeof(double));
	y_sorted=(Pdouble) calloc(n, sizeof(double));
	y_sorted_2=(Pdouble) calloc(n, sizeof(double));
	
	
	confidence=(Pdouble) calloc(n, sizeof(double));
	confidence_sorted=(Pdouble) calloc(n, sizeof(double));
	confidence_sorted_2=(Pdouble) calloc(n, sizeof(double));

	sort_indices=(Pint) calloc(n, sizeof(int));



	osd_error=(Pint) calloc(k, sizeof(int));
	BestCodeword=(Pint) calloc(n, sizeof(int));
	BestCodeword_sorted=(Pint) calloc(n, sizeof(int));
	
	STEPs_processed=(long long int *) calloc(max_bias_iter, sizeof(long long int));
	random_bias=(Pint) calloc(n, sizeof(int));
	confidence_biased=(Pdouble) calloc(n, sizeof(double));


	e_2_sub=(Pint) calloc(n-k, sizeof(int));
	e_2=(Pint) calloc(n-k, sizeof(int));
	order_0_word=(Pint) calloc(n, sizeof(int));
	diff_par_order_0_hard_dec=(Pint) calloc(n-k, sizeof(int));
	diff_par_order_0_hard_dec_e_2_sub=(Pint) calloc(n-k, sizeof(int));

	decoded_word=(Pint) calloc(n, sizeof(int));
	decoded_info_word=(Pint) calloc(k, sizeof(int));
	



	sort_indices_interm = (Pint) calloc(k,sizeof(int));
	sort_indices_sorted = (Pint) calloc(k,sizeof(int));
	confidence_sorted_interm = (Pdouble) calloc(k,sizeof(double));
	
	
	
	if(flag_enable_prepr_2_full)
	{
		nr_q_1_tau_entries = 1<<param_tau_full;	 
		max_number_q_ll_entries = CnkFunction(k,2);
	
		q_1_tau_array = (q_ll*) calloc(nr_q_1_tau_entries,sizeof(q_ll));
		q_ll_entries = (q_ll_entry*) calloc(max_number_q_ll_entries, sizeof(q_ll_entry));
		
		for(i=0; i<max_number_q_ll_entries; i++)
			q_ll_entries[i].q_vector = (Pint) calloc(n-k,sizeof(int));
	}
		
		
	/********************************************************************************************************************************************************************************/
	I_array = (Pint) calloc(k, sizeof(int));
	
	
	RedundancyBoxes=AllocGF2Mat(k,n-k);
	WHDBoxes=(Pdouble) calloc(k,sizeof(double));
	
	QMat=AllocGF2Mat(k-1, n-k);
	
	conf_dif=(Pdouble) calloc(k-1,sizeof(double));
	/********************************************************************************************************************************************************************************/

	/************************* End of memory allocation **************************************/


	/****************************** Initialization *******************************************/
	// initialization
	InitOSDLeftSysMatrix(k, n);

	OSD_params.param_lambda = param_lambda;
	OSD_params.psc_threshold = psc_threshold;
	OSD_params.nr_q_1_tau_entries = nr_q_1_tau_entries;
	OSD_params.flag_enable_prob_skip = flag_enable_prob_skip;
	OSD_params.max_STEPS = max_STEPS;
	
	printf("STEP number is %d", OSD_params.max_STEPS);


	OSD_memory.decoded_info_word = decoded_info_word;
	OSD_memory.order_0_word = order_0_word;
	OSD_memory.diff_par_order_0_hard_dec = diff_par_order_0_hard_dec;
	OSD_memory.diff_par_order_0_hard_dec_e_2_sub = diff_par_order_0_hard_dec_e_2_sub;
	OSD_memory.decoded_word = decoded_word;
	OSD_memory.q_1_tau_array = q_1_tau_array;
	OSD_memory.q_ll_entries = q_ll_entries;
	OSD_memory.QMat = QMat;
	OSD_memory.conf_dif = conf_dif;
	OSD_memory.RedundancyBoxes = RedundancyBoxes;
	OSD_memory.WHDBoxes = WHDBoxes;
	OSD_memory.I_array = I_array;
	OSD_memory.e_2_sub = e_2_sub;
	OSD_memory.e_2 = e_2;





	/**************************** End of Initialization ****************************************/


	// Ethan adding parameters into results file
	paramfptr = fopen(argv[1],"rt");
	char params_line[200] = {0};
	size_t params_line_len = 200, nread = 0;
	while (fgets(params_line, params_line_len, paramfptr) != NULL) {
		fprintf(resfptr,"* ");
		fprintf(resfptr, params_line);
		fflush(resfptr);
	}
	fprintf(resfptr,"\n");
	fclose(paramfptr);
	// End of Ethan's changes


	if(PrintFlag)
		printf("Eb/N0  ||  BER         ||    WER        ||  ML lower bound || total time  ||  time per word: \n");



	for(snrdb=snrdb1; snrdb <= snrdb2; snrdb += snrstep)// SNR per bit of information
	{
		timer = clock();


		snr=exp10(0.1*snrdb);
		noisevariance=0.5/snr;
		noisevariance /= Rate;
		noisedeviation=sqrt(noisevariance);
		ResetPRNGMT(prng_id1);// random generation of information
		ResetPRNGMT(prng_id2);// random generation of noise

		errcnt_bit=errcnt_word=errcnt_ML=0;

		count_full_OSD = 0;

		for(iter=0; iter<maxiter; iter++)
		{
			// generating n random bits, iid Bernoulli(1/2)
			PRNGMT_Binary(k, infobits, prng_id1);

			// Do the encoding via GenMat
			for(j=0; j<n; j++)
				for(codeword[j]=i=0; i<k; i++) codeword[j] ^= (infobits[i]&GenMatLeft[i][j]);
			
			// back to the original GsysLeft
			for(i=0; i<k; i++)
				for(j=0; j<n; j++) GsysLeft[i][j]=GenMatLeft[i][j];


			for(i=0; i<n; i++)
			{
				y[i]=(2.0*codeword[i]-1.0)+noisedeviation*PRNGMT_Gaussian(prng_id2);
				confidence[i] = fabs(y[i]);
			}
			
			
			/************************************* Start of preselect ************************************/

			// get indices of ordered confidence values
			for(i=0; i<n; i++) sort_indices[i] = i;
			qsort_r(sort_indices, n, sizeof(int), cmp_decreasing_order, confidence);


			// Sorting

			if(OSDLeftSysMatrix(sort_indices, GsysLeft, k, n, &newk, &perms, &countperms)==-1) // Gaussian elimination
			{
				fprintf(stderr,"%s: OSDLeftSysMatrix() failed.\n",argv[0]);
				return(-1);
			}
			
			// Apply the new ordering to the generator matrix so it can be used in the OSD - Note: perms already performed on generator
			for(i=0; i<k; i++)
				for(j=0; j<n; j++)
					GsysLeft_sorted[i][j]=GsysLeft[i][sort_indices[j]];



			//~ // Update index of sorted_channel_pos according to perms
			for(i=0; i<countperms; i++)
			{
				temp=sort_indices[perms[i].a];
				sort_indices[perms[i].a]=sort_indices[perms[i].b];
				sort_indices[perms[i].b]=temp;
			}
			


			// Apply the new ordering to y and the codeword so they correspond to the new generator matrix
			for(i=0; i<n; i++) {
				y_sorted[i] = y[sort_indices[i]];
				confidence_sorted[i] = confidence[sort_indices[i]];
			}




			count_full_OSD++;
			
			OSD_params.max_osd_order = osd_order_full;
			
			
			OSD_params.param_t=param_t_full;
			OSD_params.param_theta=param_theta_full;
			
			if(OSD_params.flag_enable_prepr_2=flag_enable_prepr_2_full);
				OSD_params.param_tau=param_tau_full;
				
			minWHD = 1000000.0;	
				
			for(bias_iter = 0; bias_iter < max_bias_iter; bias_iter++) STEPs_processed[bias_iter] = 0;
			
			for(bias_iter = 0; bias_iter < max_bias_iter; bias_iter++)
			{
				
				ret_val = runOSD(y_sorted, confidence_sorted, GsysLeft_sorted, n, k, BestCodeword_sorted, &minWHD, &(STEPs_processed[bias_iter]), &OSD_params, &OSD_memory);
				
				if(ret_val)
				{
					// put BestCodeword back in original order
					for(i=0;i<n;i++) BestCodeword[sort_indices[i]] = BestCodeword_sorted[i];
					if(ret_val == 2) goto check_errors; // Probabilistic stopping criterion met
				}
				
				// Add random bias
				PRNGMT_Binary(n, random_bias, prng_id3);
				for(i=0; i<n; i++) confidence_biased[i] = confidence[i] + param_bias_magn*(2*random_bias[i]-1);
				
				
				// get indices of ordered confidence values
				for(i=0; i<n; i++) sort_indices[i] = i;
				qsort_r(sort_indices, n, sizeof(int), cmp_decreasing_order, confidence_biased);


				// back to the original GsysLeft
				for(i=0; i<k; i++)
					for(j=0; j<n; j++) GsysLeft[i][j]=GenMatLeft[i][j];

				// Sorting

				if(OSDLeftSysMatrix(sort_indices, GsysLeft, k, n, &newk, &perms, &countperms)==-1) // Gaussian elimination
				{
					fprintf(stderr,"%s: OSDLeftSysMatrix() failed.\n",argv[0]);
					return(-1);
				}
				
				// Apply the new ordering to the generator matrix so it can be used in the OSD - Note: perms already performed on generator
				for(i=0; i<k; i++)
					for(j=0; j<n; j++)
						GsysLeft_sorted[i][j]=GsysLeft[i][sort_indices[j]];
				for(i=0; i<k; i++)
					for(j=0; j<n; j++)
						GsysLeft[i][j]=GsysLeft_sorted[i][j];
				

				
				// Update index of sorted_channel_pos according to perms
				for(i=0; i<countperms; i++)
				{
					temp=sort_indices[perms[i].a];
					sort_indices[perms[i].a]=sort_indices[perms[i].b];
					sort_indices[perms[i].b]=temp;
				}
				
				
				
				
				// Order new basis so monotonic property still valid!!!!!!!!!!

				for(i=0;i<k;i++) 
				{
					sort_indices_interm[i]=i;
					confidence_sorted_interm[i]=confidence[sort_indices[i]];
				}
				qsort_r(sort_indices_interm, k, sizeof(int), cmp_decreasing_order, confidence_sorted_interm);
				
				
				
				
				
				for(i=0; i<k; i++)
				{
					
					for(j=0; j<k; j++) GsysLeft_sorted[i][j] = i==j;
					for(j=k; j<n; j++) GsysLeft_sorted[i][j] = GsysLeft[sort_indices_interm[i]][j]; // memcopy is faster ??
					sort_indices_sorted[i] = sort_indices[sort_indices_interm[i]];
				}
				
				for(i=0; i<k; i++) sort_indices[i] = sort_indices_sorted[i];




				// Apply the new ordering to y and the codeword so they correspond to the new generator matrix
				for(i=0; i<n; i++) {
					y_sorted[i] = y[sort_indices[i]];
					confidence_sorted[i] = confidence[sort_indices[i]];
				}


				

				
			}
						
			
									
			/*** measure error rate on information bits only ***/
			check_errors:
			
#if STEP_COUNT_VERBOSE
			for(bias_iter = 0; bias_iter < max_bias_iter; bias_iter++) printf("STEPs_processed_B%d = %lld\n",bias_iter,STEPs_processed[bias_iter]);
			printf("\n\n");
#endif
			
			for(i=j=0; j<k; j++) if(BestCodeword[j]!=codeword[j]) ++i;
				
			if(i)
			{

				errcnt_bit +=i;
				++errcnt_word;

				for(WHD=j=0; j<n; j++)
					WHD += (codeword[j] == (y[j] > 0)) ? 0:confidence[j];

				if(WHD>minWHD)
					++errcnt_ML; // ML detector would also fail
					
				
			}

			if(PartialFlag)
				if(((iter+1)%PartialCount)==0)
				{
					fprintf(resfptr,"* Eb/N0=%1.2f  Pew=%d/%d=%1.2e  ML=%d/%d=%1.2e  Elapsed_time=%1.2f s\n",
					        snrdb, errcnt_word, iter+1, errcnt_word/((double)(iter+1)),errcnt_ML,iter+1,errcnt_ML/((double)iter+1),(clock()-timer)/CLOCKS_PER_SEC);
					fflush(resfptr);
				}

			if(errcnt_word >= max_word_error)
			{
				++iter;
				break;
			}
			

		}/* end of iter loop */


		
		if(PartialFlag)
		{
			fprintf(resfptr,"* Final ==> Eb/N0=%1.2f  Pew=%d/%d=%1.2e  ML=%d/%d=%1.2e  Elapsed_time=%1.2f s\n",
					        snrdb, errcnt_word, iter+1, errcnt_word/((double)(iter+1)),errcnt_ML,iter+1,errcnt_ML/((double)iter+1),(clock()-timer)/CLOCKS_PER_SEC);
			fflush(resfptr);
		}

		fprintf(resfptr,"%1.2f  %1.2e  %1.2e %1.2e %d %d %d %d %1.2e\n",
		        snrdb, errcnt_bit/((double)k*iter), errcnt_word/((double)iter),errcnt_ML/((double)iter),errcnt_bit,errcnt_word,errcnt_ML,iter,(clock()-timer)/(CLOCKS_PER_SEC*(double)iter));
		fflush(resfptr);
		if(PrintFlag)
			printf(" %1.2f      %1.2e          %1.2e          %1.2e         %1.2f s        %1.2e s\n",
			       snrdb, errcnt_bit/((double)k*iter), errcnt_word/((double)iter), errcnt_ML/((double)iter),(clock()-timer)/CLOCKS_PER_SEC,(clock()-timer)/(CLOCKS_PER_SEC*(double)iter));

	}/* end of snrdb loop */

	fclose(resfptr);
	EndPRNGMT(prng_id1);
	EndPRNGMT(prng_id2);
	return(0);
}/* end of main() */


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Gaussian elimination applied on a M*N binary matrix (M <=N)*/
/* Identity matrix is built on the left side */
/*---------------------------------------------------------------------------*/
static int LeftSystematicMatrix(Mat, M, N, newM, perms, countperms)
GF2Matrix Mat;
int M, N, *newM, *countperms;
PermTable** perms;
{
	int i, r, j;
	char *temp;
	char *coltemp;
	int countperm=0;
	PermTable *table;


	if(Mat==NULL)
	{
		fprintf(stderr,"LeftSystematicMatrix(): null argument.\n");
		return(-1);
	}

	if((N>100000)||(N<2))
	{
		fprintf(stderr,"LeftSystematicMatrix(): N out of range.\n");
		return(-1);
	}

	if((M>N)||(M<2))
	{
		fprintf(stderr,"LeftSystematicMatrix(): M out of range.\n");
		return(-1);
	}


	coltemp=(char *) malloc(M);
	table=(PermTable*) calloc(N, sizeof(PermTable));
	if((table==NULL)||(coltemp==NULL))
	{
		fprintf(stderr,"LeftSystematicMatrix(): memory allocation error.\n");
		return(-1);
	}


	for(j=0; j<M; j++)
	{
		/**** chercher le premier 1, dans la colonne j, situe sur la ligne i>=j ****/
		/**** et permuter les lignes j et i ***********/
		for(i=j; i<M; i++)
			if(Mat[i][j]) break;

		if(i<M) /** on a trouve un 1 sur la colonne j ***/
		{
			temp=Mat[i];
			Mat[i]=Mat[j];
			Mat[j]=temp;
		}
		else
			/***** pas de 1 sur la colonne j => permutation de colonnes **/
		{
			for(r=j+1; r<N; r++)
				if(Mat[j][r]) break;

			if(r >= N)
			{
				fprintf(stderr,"LeftSystematicMatrix(): All bits are unset !\n");
				fprintf(stderr,"LeftSystematicMatrix(): The matrix rank has been modified.\n");
				M=j;
				goto Sortie;
			}
			/*** permuter les colonnes r et j ****/
			for(i=0; i<M; i++)
			{
				coltemp[i]=Mat[i][j];
				Mat[i][j]=Mat[i][r];
				Mat[i][r]=coltemp[i];
			}
			/*** sauvegarder la permutation ****/
			table[countperm].a=j;
			table[countperm].b=r;
			++countperm;
		}

		/****** mettre des 0 avant et apres en ajoutant la ligne j ***/
		for(i=0; i<M; i++)
			if(Mat[i][j])
				if(i != j)
					for(r=0; r<N; r++) Mat[i][r] ^= Mat[j][r];
	}

Sortie:
	free((void *)coltemp);
	*countperms=countperm;
	*perms=table;
	*newM=M;

	return(0);
}/* end of LeftSystematicMatrix() */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static char *coltemp;
static PermTable *table;

static int InitOSDLeftSysMatrix(int M, int N)
{

	coltemp=(char *) malloc(M);
	table=(PermTable*) calloc(N, sizeof(PermTable));

	if((table==NULL)||(coltemp==NULL))
	{
		fprintf(stderr,"OSDLeftSysMatrix(): memory allocation error.\n");
		return(-1);
	}

	return(0);
}/* end of InitOSDLeftSysMatrix() */

/*---------------------------------------------------------------------------*/
/* Gaussian elimination for OSD applied on a M*N binary matrix (M <=N)*/
/* Identity matrix is built on the left side */
/*---------------------------------------------------------------------------*/
// created: 18 May 2015
static int OSDLeftSysMatrix(ch_pos, Mat, M, N, newM, perms, countperms)
Pint ch_pos;// sorted channel positions, not modified by this function, use perms to synchronize
GF2Matrix Mat;
int M, N, *newM, *countperms;
PermTable** perms;// Notice: a unique array is allocated by this function despite multiple calls
{
	int i, r, j;
	char *temp;
	int countperm=0;

	if(Mat==NULL)
	{
		fprintf(stderr,"OSDLeftSysMatrix(): null argument.\n");
		return(-1);
	}

	if((N>100000)||(N<2))
	{
		fprintf(stderr,"OSDLeftSysMatrix(): N out of range.\n");
		return(-1);
	}

	if((M>N)||(M<2))
	{
		fprintf(stderr,"OSDLeftSysMatrix(): M out of range.\n");
		return(-1);
	}

	for(j=0; j<M; j++)
	{
		/**** chercher le premier 1, dans la colonne j, situe sur la ligne i>=j ****/
		/**** et permuter les lignes j et i ***********/
		for(i=j; i<M; i++)
			if(Mat[i][ch_pos[j]]) break;

		if(i<M) /** on a trouve un 1 sur la colonne j ***/
		{
			temp=Mat[i];
			Mat[i]=Mat[j];
			Mat[j]=temp;
		}
		else
			/***** pas de 1 sur la colonne j => permutation de colonnes **/
		{
			#if 0
			// Old version - slightly faster but confidence values not strictly sorted
			for(r=j+1; r<N; r++)
				if(Mat[j][ch_pos[r]]) break;

			if(r >= N)
			{
				fprintf(stderr,"OSDLeftSysMatrix(): All bits are unset !\n");
				fprintf(stderr,"OSDLeftSysMatrix(): The matrix rank has been modified.\n");
				M=j;
				goto Sortie;
			}
			/*** permuter les colonnes r et j ****/
			for(i=0; i<M; i++)
			{
				coltemp[i]=Mat[i][ch_pos[j]];
				Mat[i][ch_pos[j]]=Mat[i][ch_pos[r]];
				Mat[i][ch_pos[r]]=coltemp[i];
			}
			/*** sauvegarder la permutation ****/
			table[countperm].a=j;
			table[countperm].b=r;
			++countperm;
			#endif
			
			#if 1
			for(r=j+1; r<N; r++)
				for(i=j; i<M; i++)			
					if(Mat[i][ch_pos[r]]) goto exit_loop;
			
			exit_loop:

			if((i >= M) && (r >= N))
			{
				fprintf(stderr,"OSDLeftSysMatrix(): All bits are unset !\n");
				fprintf(stderr,"OSDLeftSysMatrix(): The matrix rank has been modified.\n");
				M=j;
				goto Sortie;
			}
			

			temp=Mat[i];
			Mat[i]=Mat[j];
			Mat[j]=temp;
			

			/*** permuter les colonnes r et j ****/
			for(i=0; i<M; i++)
			{
				coltemp[i]=Mat[i][ch_pos[j]];
				Mat[i][ch_pos[j]]=Mat[i][ch_pos[r]];
				Mat[i][ch_pos[r]]=coltemp[i];
			}
			/*** sauvegarder la permutation ****/
			table[countperm].a=j;
			table[countperm].b=r;
			++countperm;
			#endif


		}/* end of pas de 1 sur colonne j */

		/****** mettre des 0 avant et apres en ajoutant la ligne j ***/
		for(i=0; i<M; i++)
			if(Mat[i][ch_pos[j]])
				if(i != j)
					for(r=0; r<N; r++) Mat[i][ch_pos[r]] ^= Mat[j][ch_pos[r]];

	}/* end of j loop */

Sortie:
	//free((void *)coltemp);
	*countperms=countperm;
	*perms=table;
	*newM=M;

	return(0);
}/* end of OSDLeftSysMatrix() */

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
static long long int CnkFunctionLL(int N, int K)
{
   // This function gets the total number of unique combinations based upon N and K.
   // N is the total number of items.
   // K is the size of the group.
   // Total number of unique combinations = N! / ( K! (N - K)! ).
   // This function is less efficient, but is more likely to not overflow when N and K are large.
   // Taken from:  http://blog.plover.com/math/choose.html
   //
   long long r = 1;
   long long d;
   if (K > N) return 0;
   for (d = 1; d <= K; d++)
   {
      r *= N--;
      r /= d;
   }
   return r;
}



















