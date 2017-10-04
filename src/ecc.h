#ifndef __ECC_H_
#define __ECC_H_
#include <algorithm> 
#include <iostream> // For prints
using std::min;
using std::max;

#define likely(x) __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect((x),0)
#define barrier() __asm__ __volatile__("": : :"memory")

#if (!BERSIM)
#define ecc_assert(cond,msg,a...)                                                              \
        do {                                                                                    \
                if (unlikely(!(cond))) {                                \
                        printf("assertion failed at file %s, function %s, line %d - ",  \
                                    __FILE__, __FUNCTION__, __LINE__);          \
                        printf(msg, ##a);                                                       \
                        printf("\n");                                                           \
                        exit(1);                                                                        \
                }                                                                               \
        } while (0)
#else
#define ecc_assert(cond,msg,a...)   
#endif

//static __inline void *ecc_malloc_err(int size, char* file, int line)
// Never put non-static function definitions in a header file! either make the function
// static or move the definition to a .c file
static void *ecc_malloc_err(int size, const char* file, int line)
{
        void *p;
        ecc_assert(size != 0,"Size is equall 0, file %s, line %d\n", file, line);
        p = malloc(size);
        ecc_assert(p != NULL, "unable to allocate %d bytes, out of memory, file %s, line %d\n", size, file, line);
        memset(p,0,(int) size);

        return p;
}

#define ecc_malloc(size) ecc_malloc_err((size),__FILE__,__LINE__)

// #define max(x, y) (((x)>=(y)) ? (x) : (y))
// #define min(x, y) (((x) <(y)) ? (x) : (y))

#define clamp(x, min_val, max_val) (min(max((x), (min_val)), (max_val)))

#define N_perm_list	20
#define MAX_STR		256

typedef long long i64;

// Defining an (n-k) by n matrix (integer)
typedef struct int_mtx_s {
  int 
    n,		// number of cols (n)
    m;		// number of rows (n-k)
  int *data;	// the matrix itself	
} int_mtx_t;

// Defining an (n-k) by n matrix (double)
typedef struct double_mtx_s {
  int 
    n,		// number of cols (n)
    m;		// number of rows (n-k)
  double *data;	// the matrix itself	
} double_mtx_t;

static __inline int get_int_element(int_mtx_t *mtx, int line_idx, int col_idx)
{
  ecc_assert((line_idx>=0) && (line_idx<mtx->m), "Line exceeds allowed (%d vs %d)\n", mtx->m, line_idx);
  ecc_assert((col_idx >=0) && (col_idx <mtx->n), "Col exceeds allowed (%d vs %d)\n", mtx->n, col_idx);
  return mtx->data[line_idx * mtx->n + col_idx];
}

static __inline double get_double_element(double_mtx_t *mtx, int line_idx, int col_idx)
{
  ecc_assert((line_idx>=0) && (line_idx<=mtx->m), "Line exceeds allowed (%d vs %d)\n", mtx->m, line_idx);
  ecc_assert((col_idx >=0) && (col_idx <=mtx->n), "Col exceeds allowed (%d vs %d)\n", mtx->n, col_idx);
  return mtx->data[line_idx * mtx->n + col_idx];
}

static __inline void put_double_element(double_mtx_t *mtx, int line_idx, int col_idx, double in_data)
{
  ecc_assert((line_idx>=0) && (line_idx<=mtx->m), "Line exceeds allowed (%d vs %d)\n", mtx->m, line_idx);
  ecc_assert((col_idx >=0) && (col_idx <=mtx->n), "Col exceeds allowed (%d vs %d)\n", mtx->n, col_idx);
  mtx->data[line_idx * mtx->n + col_idx] = in_data;
}

static __inline void print_double_mtx(double_mtx_t *mtx)
{
  int i, j;

  for (i=0; i<mtx->m; i++)
    {
      printf("| ");
      for (j=0; j<mtx->n; j++)
	printf("%12.6g ", get_double_element(mtx, i, j));
      printf(" |\n");
    }
}


static __inline void print_int_mtx(int_mtx_t *mtx)
{
  int i, j;

  for (i=0; i<mtx->m; i++)
    {
      printf("| ");
      for (j=0; j<mtx->n; j++)
	printf("%d ", get_int_element(mtx, i, j));
      printf(" |\n");
    }
}


static __inline int sum_int_mtx(int_mtx_t *mtx)
{
	int i, j, s = 0;
	
	for (i=0; i<mtx->m; i++) {
		for (j=0; j<mtx->n; j++)
			s += get_int_element(mtx, i, j);
	}
	return s;
}





static inline void init_double_mtx(double_mtx_t *mtx, int m, int n)
{
  mtx->n = n;
  mtx->m = m;
  mtx->data = (double*)ecc_malloc(sizeof(double)*n*m);
}

static inline void init_int_mtx(int_mtx_t *mtx, int m, int n)
{
  mtx->n = n;
  mtx->m = m;
  mtx->data = (int*)ecc_malloc(sizeof(int)*n*m);
  memset(mtx->data, 0, sizeof(int)*n*m); 
}

static inline void init_int_mtx_from_imtx(int_mtx_t *mtx, int_mtx_t *imtx)
{
  mtx->n = imtx->n;
  mtx->m = imtx->m;
  mtx->data = (int*)ecc_malloc(sizeof(int)*mtx->n*mtx->m);
  memcpy(mtx->data, imtx->data, sizeof(int)*mtx->n*mtx->m);
}


static __inline void put_int_element(int_mtx_t *mtx, int line_idx, int col_idx, int in_data)
{
  ecc_assert((line_idx>=0) && (line_idx<mtx->m), "Line exceeds allowed (%d vs %d)\n", mtx->m, line_idx);
  ecc_assert((col_idx >=0) && (col_idx <mtx->n), "Line exceeds allowed (%d vs %d)\n", mtx->n, col_idx);
  mtx->data[line_idx * mtx->n + col_idx] = in_data;
}

static __inline void copy_int_mtx(int_mtx_t *dest, int_mtx_t *source)
{
  int i=0, j=0; 

  ecc_assert( (dest->m == source->m), "Inconsistency with number of lines %d vs. %d\n", dest->m ,source->m);
  ecc_assert( (dest->n == source->n), "Inconsistency with number of lines %d vs. %d\n", dest->n ,source->n);

  for (i=0; i<dest->m; i++)
    for (j=0; j<dest->n; j++)
      put_int_element(dest, i, j, get_int_element(source, i, j));
}

static __inline void copy_double_mtx(double_mtx_t *dest, double_mtx_t *source)
{
  int i, j;

  ecc_assert( (dest->m == source->m), "Inconsistency with number of lines %d vs. %d\n", dest->m ,source->m);
  ecc_assert( (dest->n == source->n), "Inconsistency with number of lines %d vs. %d\n", dest->n ,source->n);

  for (i=0; i<dest->m; i++)
    for (j=0; j<dest->n; j++)
      put_double_element(dest, i, j, get_double_element(source, i, j));
}





static void dump_int_mtx(int_mtx_t *mtx, char *fname)
{
  FILE *fp;
  int type =1;

  if (!(fp = fopen(fname, "wb")))
    ecc_assert(0, "couldnt open %s for dumping", fname);

  fwrite(&type, sizeof(int), 1, fp);
  fwrite(&mtx->m, sizeof(int), 1, fp);
  fwrite(&mtx->n, sizeof(int), 1, fp);
  fwrite(mtx->data, sizeof(int), mtx->m*mtx->n, fp);
  fclose(fp);
}

static void dump_double_mtx(double_mtx_t *mtx, char *fname)
{
  FILE *fp;
  int type =2;

  if (!(fp = fopen(fname, "wb")))
    ecc_assert(0, "couldnt open %s for dumping", fname);

  fwrite(&type, sizeof(int), 1, fp);
  fwrite(&mtx->m, sizeof(int), 1, fp);
  fwrite(&mtx->n, sizeof(int), 1, fp);
  fwrite(mtx->data, sizeof(double), mtx->m*mtx->n, fp);
  fclose(fp);
}

static __inline int_mtx_t * eye(int n)
{
  int_mtx_t* mtx = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  int i;

  init_int_mtx(mtx, n, n);

  for (i=0; i<n; i++) 
    put_int_element(mtx, i, i, 1);
  
  return mtx;
}

static __inline void multiply_int_mtx(int_mtx_t* src0, int_mtx_t* src1, int_mtx_t* dest)
{
  int i, j, z, sum;
  
  ecc_assert(src0->n == src1->m, "matrices size don't match (%d, %d) * (%d, %d)\n", 
	     src0->m, src0->n, src1->m, src1->n);
  ecc_assert((src0->m == dest->m) && (src1->n == dest->n), "src and dest don't match (%d, %d) * (%d, %d) != (%d, %d)\n", 
	     src0->m, src0->n, src1->m, src1->n, dest->m, dest->n);

  for (i=0; i<src0->m; i++)
    for (j=0; j<src1->n; j++)
      {
	sum = 0;
	for (z=0; z<src0->n; z++)
	  sum += get_int_element(src0, i, z) * get_int_element(src1, z, j);
	put_int_element(dest,     i, j, sum%2);
      }
}

static __inline void multiply_int_mtx_nm(int_mtx_t* src0, int_mtx_t* src1, int_mtx_t* dest)
{
  int i, j, z, sum;
  
  ecc_assert(src0->n == src1->m, "matrices size don't match (%d, %d) * (%d, %d)\n", 
	     src0->m, src0->n, src1->m, src1->n);
  ecc_assert((src0->m == dest->m) && (src1->n == dest->n), "src and dest don't match (%d, %d) * (%d, %d) != (%d, %d)\n", 
	     src0->m, src0->n, src1->m, src1->n, dest->m, dest->n);

  for (i=0; i<src0->m; i++)
    for (j=0; j<src1->n; j++)
      {
	sum = 0;
	for (z=0; z<src0->n; z++)
	  sum += get_int_element(src0, i, z) * get_int_element(src1, z, j);
	put_int_element(dest,     i, j, sum);
      }
}

static __inline void multiply_int_mtx_diag(int_mtx_t* src0, int_mtx_t* src1, int_mtx_t* dest)
{
  int i, j;
  
  multiply_int_mtx_nm(src0, src1, dest);
  ecc_assert(dest->m == dest->n, "Matrix must be sqyare %d, %d\n", dest->m, dest->n);

  for (i=0; i<src0->m; i++)
    for (j=0; j<src1->n; j++)
      if (i!=j) put_int_element(dest,     i, j, 0);
}

static __inline void multiply_int_mtx_zero(int_mtx_t* src0, int_mtx_t* src1, int_mtx_t* dest)
{
  int i;
  
  multiply_int_mtx_nm(src0, src1, dest);
  ecc_assert(dest->m == dest->n, "Matrix must be sqyare %d, %d\n", dest->m, dest->n);
  
  for (i=0; i<dest->m; i++)    
    put_int_element(dest,     i, i, 0);
}

static __inline void direct_product_int_mtx(int_mtx_t* src0, int_mtx_t* src1, int_mtx_t* dest)
{
  int i, j;
 
  ecc_assert((src0->m == src1->m) && (src0->m == dest->m), "error with m dimension %d, %d, %d\n", src0->m, src1->m, dest->m);
  ecc_assert((src0->n == src1->n) && (src0->n == dest->n), "error with n dimension %d, %d, %d\n", src0->n, src1->n, dest->n);
 
  for (i=0; i<dest->m; i++)    
    for (j=0; j<dest->n; j++)        
      put_int_element(dest,     i, j, get_int_element(src0, i, j)* get_int_element(src1, i, j));
}


static __inline int trace_int_mtx(int_mtx_t* src0)
{
  int i, sum = 0;
  ecc_assert(src0->n == src0->m, "matrix must be square: %d %d\n", src0->n, src0->m);
  
  for (i=0; i<src0->n; i++)
    sum += get_int_element(src0, i, i);
  return sum;
}

static __inline void diagonal_int_mtx(int_mtx_t* src0, double *p)
{
  int i;
  
  for (i=0; i<src0->n; i++)
    p[i] = get_int_element(src0, i, i);
}


static __inline void subtract_int_mtx(int_mtx_t* src0, int_mtx_t* src1, int_mtx_t* dest)
{
  int diff, i, j; 
  ecc_assert( (src0->m == src1->m) && (src0->m == dest->m), "matrices size don't match (%d, %d, %d)\n", src0->m, src1->m, dest->m);
  ecc_assert( (src0->n == src1->n) && (src0->n == dest->n), "matrices size don't match (%d, %d, %d)\n", src0->n, src1->n, dest->n);

  for (i=0; i<src0->m; i++)
    for (j=0; j<src1->n; j++)
      {
	diff = get_int_element(src0, i,j) - get_int_element(src1, i,j);
	put_int_element(dest, i, j, diff);
      }
}

static __inline int_mtx_t* transpose_int_mtx(int_mtx_t* src)
{
  int_mtx_t *p;
  int i, j;

  p = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  init_int_mtx(p, src->n, src->m);

  for (i=0; i<src->n; i++)
    for (j=0; j<src->m; j++)
      put_int_element(p, i, j, get_int_element(src, j, i)); 

  return p;

}

typedef struct gf2_vec_s {
  int n;
  int *data;	
} gf2_vec_t;

typedef struct real_vec_s {
  int n;
  double *data;	
} real_vec_t;


typedef struct int_vec_s {
  int n;
  int *data;	
} int_vec_t;

// Initializing the given matrix: assigning number of rows, columns, and allocating memory for the data
// Then assigning the data values given at '*data' to the created matrix
static __inline int init_int_mtx_from_array(int_mtx_t *H, int m, int n, int *data)
{
  H->n = n;
  H->m = m;
  H->data = (int*)ecc_malloc(sizeof(int)*m*n);
  memcpy(H->data, data, sizeof(int)*m*n);
  return 1;
}

static __inline int init_gf2_vec(gf2_vec_t *vec, int n)
{
  vec->n = n;
  vec->data = (int*)ecc_malloc(sizeof(int)*n);
  return 1;
}

static __inline int init_real_vec(real_vec_t *vec, int n)
{
  vec->n = n;
  vec->data = (double*)ecc_malloc(sizeof(double)*n);
  return 1;
}

static __inline int init_int_vec(int_vec_t *vec, int n)
{
  vec->n = n;
  vec->data = (int*)ecc_malloc(sizeof(int)*n);
  return 1;
}


static __inline void free_double_mtx(double_mtx_t *x)
{
  if (x)
    {
      if (x->data)
	free(x->data);
      free(x);
    }
}

static __inline void free_int_mtx(int_mtx_t *x)
{
  if (x)
    {
      if (x->data)
	free(x->data);
      free(x);
    }
}

static __inline void free_real_vec(real_vec_t *x)
{
  if (x)
    {
      if (x->data)
	free(x->data);
      free(x);
    }
}

static __inline void free_gf2_vec(gf2_vec_t *x)
{
  if (x)
    {
      if (x->data)
	free(x->data);
      free(x);
    }
}


static __inline void free_int_vec(int_vec_t *x)
{
  if (x)
    {
      if (x->data)
	free(x->data);
      free(x);
    }
}
 

static __inline void sum_double_mtx(double_mtx_t *E_s, real_vec_t *sum_E_s)
{
  int j, i; 
  double sum;
  
  for (i=0; i<E_s->n; i++)
    {
      sum = 0;
      for (j=0; j<E_s->m; j++)
	sum += get_double_element(E_s, j, i);
      sum_E_s->data[i] = sum;
    }
}

static __inline void copy_rvec(real_vec_t *src0, real_vec_t *dest)
{
  int i;
  ecc_assert(src0->n == dest->n, "Vecors must be in the same length\n");
    
  for (i=0; i<src0->n; i++)
    dest->data[i] = src0->data[i];
}


static __inline void sum_rvec(real_vec_t *src0, real_vec_t *src1, real_vec_t *dest)
{
  int i;
  ecc_assert(src0->n == src1->n, "Vecors must be in the same length\n");
  ecc_assert(src0->n == dest->n, "Vecors must be in the same length\n");

  for (i=0; i<src0->n; i++)
    dest->data[i] = src0->data[i] + src1->data[i];
}

static __inline void copy_gf2_vec(gf2_vec_t *src, gf2_vec_t *dest)
{
  int i=0;
  ecc_assert(src->n == dest->n, "Mismatched vectors %d vs, %d\n", src->n, dest->n);
  for (i=0; i<dest->n; i++)
    dest->data[i] = src->data[i];
}

static __inline int diff_gf2_vec(gf2_vec_t *src, gf2_vec_t *dest)
{
  int i=0, diff = 0;
  ecc_assert(src->n == dest->n, "Mismatched vectors %d vs, %d\n", src->n, dest->n);
  for (i=0; i<dest->n; i++) { 
	  if (dest->data[i] != src->data[i]) 
		  diff++;
  }
  return diff;
}


static __inline void print_rvec(real_vec_t *vec)
{
  int i;

  for (i=0; i<vec->n; i++)
    printf("%06.3g ", vec->data[i]);
  printf("\n");
}


static __inline void print_ivec(int_vec_t *vec)
{
  int i;

  for (i=0; i<vec->n; i++)
    printf("%d ", vec->data[i]);
  printf("\n");
}


static __inline void print_gf2vec(gf2_vec_t *vec)
{
  int i;

  for (i=0; i<vec->n; i++)
    printf("%d ", vec->data[i]);
  printf("\n");
}

static __inline void flip_columns_int_mtx(int_mtx_t *mtx)
{
  int i, j;
  int_mtx_t *temp_mtx;
  temp_mtx = (int_mtx_t *)ecc_malloc(sizeof(int_mtx_t));
  init_int_mtx(temp_mtx, mtx->m, mtx->n);
  copy_int_mtx(temp_mtx, mtx);

  for (i=0; i<mtx->m; i++)
    for (j=0; j<mtx->n; j++)
      put_int_element(mtx, i, j, get_int_element(temp_mtx, i, mtx->n - 1 - j));

  free_int_mtx(temp_mtx);
}

#endif



