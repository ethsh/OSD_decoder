#ifndef __SCC_H_
#define __SCC_H_


// Short Cycle Counter
typedef struct scc_s {
  int U, W;
  int g;			// girth
  int Ng;			// N_g
  int Ng2;			// N_{g+2}
  int Ng4;			// N_{g+4}

  // Vectors that store the number of cycles of length
  // g, g+2 and g+4 incident on each vertex in U.
  double* Ng_per_u;
  double* Ng2_per_u;
  double* Ng4_per_u;

  int_mtx_t *E, *ET;
  int_mtx_t *L_U_temp, *L_W_temp;

  // Matrices used to count short cycles (4,6,8).
  int_mtx_t *P_U_2, *P_W_2;				 // P_2^\mathcal{U,W}
  int_mtx_t *P_U_2_c2, *P_W_2_c2;		 // \binom{P_2^\mathcal{U,W}}{2}
  int_mtx_t *P_U_3, *P_W_3;				 // P_3^\mathcal{U,W}
  int_mtx_t *P_U_4, *P_W_4;				 // P_4^\mathcal{U,W}
  int_mtx_t *P_U_5, *P_W_5;				 // P_5^\mathcal{U,W}
  int_mtx_t *P_U_6, *P_W_6;				 // P_6^\mathcal{U,W}
  int_mtx_t *P_U_7, *P_W_7;				 // P_7^\mathcal{U,W}
  int_mtx_t *L_U_0_2_m1, *L_W_0_2_m1;    // \max{L_{(0,2)}^\mathcal{U,W}-1,0}
  int_mtx_t *L_U_0_2_m2, *L_W_0_2_m2;    // \max{L_{(0,2)}^\mathcal{U,W}-2,0}
  int_mtx_t *L_U_1_2, *L_W_1_2;			 // L_{(1,2)}^\mathcal{U,W}
  int_mtx_t *L_U_0_4, *L_W_0_4;			 // L_{(0,4)}^\mathcal{U,W}
  int_mtx_t *L_U_2_2, *L_W_2_2;			 // L_{(2,2)}^\mathcal{U,W}
  int_mtx_t *L_U_1_4, *L_W_1_4;			 // L_{(1,4)}^\mathcal{U,W}
  int_mtx_t *L_U_3_2, *L_W_3_2;			 // L_{(3,2)}^\mathcal{U,W}
  int_mtx_t *L_U_0_6, *L_W_0_6;			 // L_{(0,6)}^\mathcal{U,W}
  int_mtx_t *L_U_2_4, *L_W_2_4;			 // L_{(2,4)}^\mathcal{U,W}
  int_mtx_t *L_U_4_2, *L_W_4_2;			 // L_{(4,2)}^\mathcal{U,W}
  int_mtx_t *L_U_1_6, *L_W_1_6;			 // L_{(1,6)}^\mathcal{U,W}
  int_mtx_t *L_U_3_4;					 // L_{(3,4)}^\mathcal{U}
  int_mtx_t *L_U_5_2, *L_W_5_2;			 // L_{(5,2)}^\mathcal{U,W}
  int_mtx_t *L_U_0_8, *L_W_0_8;			 // L_{(0,8)}^\mathcal{U,W}
	

  // Long cycle matrices.
  
} scc_t;


void scc_init(scc_t *scc, int_mtx_t *Ein);
void scc_free(scc_t *scc);
void scc_count(scc_t *scc);

#endif
