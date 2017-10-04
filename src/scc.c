#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include "ecc.h"
#include "scc.h"





void scc_init(scc_t *scc, int_mtx_t *Ein)
{
  scc->U		= Ein->m;
  scc->W		= Ein->n;
  scc->Ng_per_u		= (double*)ecc_malloc(sizeof(double)*scc->U);
  scc->Ng2_per_u	= (double*)ecc_malloc(sizeof(double)*scc->U);
  scc->Ng4_per_u 	= (double*)ecc_malloc(sizeof(double)*scc->U);

  memset(scc->Ng_per_u,  0,scc->U*sizeof(double));
  memset(scc->Ng2_per_u ,0,scc->U*sizeof(double));
  memset(scc->Ng4_per_u ,0,scc->U*sizeof(double));


  scc->E  = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  init_int_mtx_from_imtx(scc->E, Ein);; 
  scc->ET = transpose_int_mtx(Ein);
  
  scc->L_U_temp = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  scc->L_W_temp = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  scc->P_U_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  scc->P_W_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  scc->P_U_2_c2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  scc->P_W_2_c2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  scc->P_U_3 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));
  scc->P_W_3 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->P_U_4 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->P_W_4 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->P_U_5 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->P_W_5 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->P_U_6 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->P_W_6 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->P_U_7 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->P_W_7 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_0_2_m1 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_0_2_m1 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));    
  scc->L_U_0_2_m2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_0_2_m2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));    
  scc->L_U_1_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_1_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_0_4 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_0_4 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_2_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_2_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_1_4 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_1_4 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_3_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_3_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_0_6 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_0_6 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_2_4 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_2_4 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_4_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_4_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_1_6 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_1_6 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_3_4 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_5_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_5_2 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		
  scc->L_U_0_8 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t)); 
  scc->L_W_0_8 = (int_mtx_t*)ecc_malloc(sizeof(int_mtx_t));		

}

void scc_free(scc_t *scc)
{
  free(scc->Ng_per_u);
  free(scc->Ng2_per_u);
  free(scc->Ng4_per_u);
  free_int_mtx(scc->E);
  free_int_mtx(scc->ET);
  free_int_mtx(scc->L_U_temp);
  free_int_mtx(scc->L_W_temp);
  free_int_mtx(scc->P_U_2);
  free_int_mtx(scc->P_W_2);
  free_int_mtx(scc->P_U_2_c2);
  free_int_mtx(scc->P_W_2_c2);
  free_int_mtx(scc->P_U_3);
  free_int_mtx(scc->P_W_3);
  free_int_mtx(scc->P_U_4);
  free_int_mtx(scc->P_W_4);
  free_int_mtx(scc->P_U_5);
  free_int_mtx(scc->P_W_5);
  free_int_mtx(scc->P_U_6);
  free_int_mtx(scc->P_W_6);
  free_int_mtx(scc->P_U_7);
  free_int_mtx(scc->P_W_7);
  free_int_mtx(scc->L_U_0_2_m1);
  free_int_mtx(scc->L_W_0_2_m1);
  free_int_mtx(scc->L_U_0_2_m2);
  free_int_mtx(scc->L_W_0_2_m2);
  free_int_mtx(scc->L_U_1_2);
  free_int_mtx(scc->L_W_1_2);
  free_int_mtx(scc->L_U_0_4);
  free_int_mtx(scc->L_W_0_4);
  free_int_mtx(scc->L_U_2_2);
  free_int_mtx(scc->L_W_2_2);
  free_int_mtx(scc->L_U_1_4);
  free_int_mtx(scc->L_W_1_4);
  free_int_mtx(scc->L_U_3_2);
  free_int_mtx(scc->L_W_3_2);
  free_int_mtx(scc->L_U_0_6);
  free_int_mtx(scc->L_W_0_6);
  free_int_mtx(scc->L_U_2_4);
  free_int_mtx(scc->L_W_2_4);
  free_int_mtx(scc->L_U_4_2);
  free_int_mtx(scc->L_W_4_2);
  free_int_mtx(scc->L_U_1_6);
  free_int_mtx(scc->L_W_1_6);
  free_int_mtx(scc->L_U_3_4);
  free_int_mtx(scc->L_U_5_2);
  free_int_mtx(scc->L_W_5_2);
  free_int_mtx(scc->L_U_0_8);
  free_int_mtx(scc->L_W_0_8);

}

void scc_process_P_U_2(scc_t *scc)
{
  int rr, oo, cc, pp, val;
  
  init_int_mtx(scc->P_U_2,	scc->E->m, scc->ET->n);
  init_int_mtx(scc->P_U_2_c2,	scc->E->m, scc->ET->n);
  init_int_mtx(scc->L_U_0_2_m1, scc->E->m, scc->ET->n);
  init_int_mtx(scc->L_U_0_2_m2, scc->E->m, scc->ET->n);

  //P_U_2 = E*E' - L_U_0_2
  multiply_int_mtx_nm(scc->E, scc->ET, scc->P_U_2);
  
  for( rr = 0, oo = 0; rr < scc->P_U_2->m; rr++, oo += scc->P_U_2->n )
    {
      for( cc = 0, pp = oo; cc < scc->P_U_2->n; cc++, pp++ )
	{
	  val = scc->P_U_2->data[pp];
	  
	  if( val > 1 ) 
	    scc->P_U_2_c2->data[pp] = val*(val-1)/2;
	  
	  if( rr == cc )
	    {
	      scc->P_U_2->data[pp] = 0; 
	      scc->P_U_2_c2->data[pp] = 0;
	      if( val > 1 )
		{
		  scc->L_U_0_2_m1->data[pp] = val-1;
		  if( val > 2 ) 
		    scc->L_U_0_2_m2->data[pp] = val-2;
		} 				
	    }
	}
    }
}

void scc_process_P_W_2(scc_t *scc)
{

  int rr, oo, cc, pp, val;
  
  init_int_mtx(scc->P_W_2,	scc->ET->m, scc->E->n);
  init_int_mtx(scc->P_W_2_c2,	scc->ET->m, scc->E->n);
  init_int_mtx(scc->L_W_0_2_m1, scc->ET->m, scc->E->n);
  init_int_mtx(scc->L_W_0_2_m2, scc->ET->m, scc->E->n);

  //P_W_2 = E'*E - L_W_0_2
  multiply_int_mtx_nm(scc->ET, scc->E, scc->P_W_2);

  for( rr = 0, oo = 0; rr < scc->P_W_2->m; rr++, oo += scc->P_W_2->n )
    {
      for( cc = 0, pp = oo; cc < scc->P_W_2->n; cc++, pp++ )
	{
	  val = scc->P_W_2->data[pp];
	  
	  if( val > 1 ) 
	    scc->P_W_2_c2->data[pp] = val*(val-1)/2;
	  
	  if( rr == cc )
	    {
	      scc->P_W_2->data[pp] = 0; 
	      scc->P_W_2_c2->data[pp] = 0;
	      if( val > 1 )
		{
		  scc->L_W_0_2_m1->data[pp] = val-1;
		  if( val > 2 ) 
		    scc->L_W_0_2_m2->data[pp] = val-2;
		} 				
	    }
	}
    }
}


int scc_count_four_cycles(scc_t *scc)
{	

  //// Compute P_U_2, P_W_2, L_U_0_2_m1, L_W_0_2_m1, 
  //// L_U_0_2_m2, L_W_0_2_m2, P_U_2_c2, and P_W_2_c2
  //// simulataneously for speed.
  //process_P_U_2();
  //process_P_W_2();
  scc_process_P_U_2(scc);
  scc_process_P_W_2(scc);
  
  //// Compute L_U_1_2 and L_W_1_2.
  //L_U_1_2_.matrix_mult(E_,L_W_0_2_m1_);
  //L_W_1_2_.matrix_mult(ET_,L_U_0_2_m1_);
  init_int_mtx(scc->L_U_1_2, scc->E->m, scc->L_W_0_2_m1->n);
  init_int_mtx(scc->L_W_1_2, scc->ET->m, scc->L_U_0_2_m1->n);
  multiply_int_mtx_nm(scc->E,  scc->L_W_0_2_m1, scc->L_U_1_2);
  multiply_int_mtx_nm(scc->ET, scc->L_U_0_2_m1, scc->L_W_1_2);

  //// Compute P_U_3, P_W_3.
  //P_U_3_.matrix_mult(P_U_2_,E_);
  //P_U_3_ -= L_U_1_2_;
  //P_W_3_.transpose(P_U_3_);
  init_int_mtx(scc->P_U_3, scc->P_U_2->m, scc->E->n);
  multiply_int_mtx_nm(scc->P_U_2, scc->E, scc->P_U_3);
  subtract_int_mtx(scc->P_U_3, scc->L_U_1_2, scc->P_U_3);
  scc->P_W_3 = transpose_int_mtx(scc->P_U_3);
	
  //// Compute L_U_0_4, L_W_0_4.
  //L_U_0_4_.mx_mult_diag(P_U_3_,ET_);
  //L_W_0_4_.mx_mult_diag(P_W_3_,E_);
  init_int_mtx(scc->L_U_0_4, scc->P_U_3->m, scc->ET->n);
  init_int_mtx(scc->L_W_0_4, scc->P_W_3->m, scc->E->n);
  multiply_int_mtx_diag(scc->P_U_3, scc->ET, scc->L_U_0_4);
  multiply_int_mtx_diag(scc->P_W_3, scc->E,  scc->L_W_0_4);
  

  //int temp = L_U_0_4_.int_trace();
  //if( temp == 0 ) return 0;
  //else
  //  {
  //    g_  = 4;
  //    Ng_ = temp/4;
  //    L_U_0_4_.diagonal(Ng_per_u_);
  //    return 1;
  //  }
  int temp = trace_int_mtx(scc->L_U_0_4);
  if( temp == 0 ) return 0;
  else
    {
      scc->g  = 4;
      scc->Ng = temp/4;
      //      printf("Ng = %d\n", scc->Ng);
      diagonal_int_mtx(scc->L_U_0_4, scc->Ng_per_u);
      return 1;
    }
}

void scc_count_six_eight_cycles(scc_t *scc)
{
  //// Compute L_U_2_2, L_W_2_2.
  //L_U_2_2_.mx_mult_zero(E_,L_W_1_2_);
  //L_W_2_2_.mx_mult_zero(ET_,L_U_1_2_);
  
  init_int_mtx(scc->L_U_2_2, scc->E->m,  scc->L_W_1_2->n);
  init_int_mtx(scc->L_W_2_2, scc->ET->m, scc->L_U_1_2->n);
  multiply_int_mtx_zero(scc->E,  scc->L_W_1_2, scc->L_U_2_2);
  multiply_int_mtx_zero(scc->ET, scc->L_U_1_2, scc->L_W_2_2);


  //// Comput P_U_4, P_W_4.
  //P_U_4_.matrix_mult(P_U_3_,ET_); 
  //P_U_4_ -= L_U_0_4_; 
  //P_U_4_ -= L_U_2_2_;
  init_int_mtx(scc->P_U_4, scc->P_U_3->m,  scc->ET->n);
  multiply_int_mtx_nm( scc->P_U_3, scc->ET, scc->P_U_4);
  subtract_int_mtx(scc->P_U_4, scc->L_U_0_4, scc->P_U_4);
  subtract_int_mtx(scc->P_U_4, scc->L_U_2_2, scc->P_U_4);

  //P_W_4_.matrix_mult(P_W_3_,E_);  
  //P_W_4_ -= L_W_0_4_; 
  //P_W_4_ -= L_W_2_2_;
  init_int_mtx(scc->P_W_4, scc->P_W_3->m, scc->E->n);
  multiply_int_mtx_nm(scc->P_W_3, scc->E, scc->P_W_4);
  subtract_int_mtx(scc->P_W_4, scc->L_W_0_4, scc->P_W_4);
  subtract_int_mtx(scc->P_W_4, scc->L_W_2_2, scc->P_W_4);

  //// Compute L_U_1_4, L_W_1_4.
  //L_U_1_4_.matrix_mult(E_,L_W_0_4_);  
  //L_U_temp_ = P_U_3_*E_;
  //L_U_1_4_ -= L_U_temp_;
  //L_U_1_4_ -= L_U_temp_; 
  init_int_mtx(scc->L_U_1_4, scc->E->m, scc->L_W_0_4->n);
  multiply_int_mtx_nm(scc->E, scc->L_W_0_4, scc->L_U_1_4);
  init_int_mtx(scc->L_U_temp, scc->E->m, scc->E->n);
  direct_product_int_mtx(scc->P_U_3, scc->E, scc->L_U_temp);
  subtract_int_mtx(scc->L_U_1_4, scc->L_U_temp, scc->L_U_1_4);
  subtract_int_mtx(scc->L_U_1_4, scc->L_U_temp, scc->L_U_1_4);
  
  //L_W_1_4_.matrix_mult(ET_,L_U_0_4_); 
  //L_W_temp_ = P_W_3_*ET_;
  //L_W_1_4_ -= L_W_temp_;
  //L_W_1_4_ -= L_W_temp_;
  init_int_mtx(scc->L_W_1_4, scc->ET->m, scc->L_U_0_4->n);
  multiply_int_mtx_nm(scc->ET, scc->L_U_0_4, scc->L_W_1_4);
  init_int_mtx(scc->L_W_temp, scc->ET->m, scc->ET->n);
  direct_product_int_mtx(scc->P_W_3, scc->ET, scc->L_W_temp);
  subtract_int_mtx(scc->L_W_1_4, scc->L_W_temp, scc->L_W_1_4);
  subtract_int_mtx(scc->L_W_1_4, scc->L_W_temp, scc->L_W_1_4);

  //// Compute L_U_3_2, L_W_3_2.
  //L_U_3_2_.matrix_mult(P_U_3_,L_W_0_2_m1_); 
  //L_U_3_2_ -= L_U_temp_;
  init_int_mtx(scc->L_U_3_2, scc->P_U_3->m, scc->L_W_0_2_m1->n);
  multiply_int_mtx_nm(scc->P_U_3, scc->L_W_0_2_m1, scc->L_U_3_2);
  subtract_int_mtx(scc->L_U_3_2, scc->L_U_temp, scc->L_U_3_2);
  
  //L_W_3_2_.matrix_mult(P_W_3_,L_U_0_2_m1_); 
  //L_W_3_2_ -= L_W_temp_;
  init_int_mtx(scc->L_W_3_2, scc->P_W_3->m, scc->L_U_0_2_m1->n);
  multiply_int_mtx_nm(scc->P_W_3, scc->L_U_0_2_m1, scc->L_W_3_2);
  subtract_int_mtx(scc->L_W_3_2, scc->L_W_temp, scc->L_W_3_2);
  
  //// Compute P_U_5, P_W_5.  
  //P_U_5_.matrix_mult(P_U_4_,E_);  
  //P_U_5_ -= L_U_1_4_; 
  //P_U_5_ -= L_U_3_2_; 
  init_int_mtx(scc->P_U_5, scc->P_U_4->m, scc->E->n);
  multiply_int_mtx_nm(scc->P_U_4, scc->E, scc->P_U_5);
  subtract_int_mtx(scc->P_U_5, scc->L_U_1_4, scc->P_U_5);
  subtract_int_mtx(scc->P_U_5, scc->L_U_3_2, scc->P_U_5);

  //P_W_5_.transpose(P_U_5_); 
  scc->P_W_5 = transpose_int_mtx(scc->P_U_5);
	
  //// Compute L_U_0_6, L_W_0_6.
  //L_U_0_6_.mx_mult_diag(P_U_5_,ET_);
  //L_W_0_6_.mx_mult_diag(P_W_5_,E_);
  init_int_mtx(scc->L_U_0_6, scc->P_U_5->m, scc->ET->n);
  multiply_int_mtx_diag(scc->P_U_5, scc->ET, scc->L_U_0_6);

  init_int_mtx(scc->L_W_0_6, scc->P_W_5->m, scc->E->n);
  multiply_int_mtx_diag(scc->P_W_5, scc->E, scc->L_W_0_6);

  //int temp = L_U_0_6_.int_trace();
  //if( g_ == 4 ) 
  //  {
  //    Ng2_ = temp/6; 
  //    L_U_0_6_.diagonal(Ng2_per_u_);
  //  }
  //else if( temp > 0 ) g_ = 6;

  int temp = trace_int_mtx(scc->L_U_0_6);
  if (scc->g == 4)
    {
      scc->Ng2 = temp/6;
      //printf("Ng2 = %d\n",scc->Ng2);
      diagonal_int_mtx(scc->L_U_0_6, scc->Ng2_per_u);
    }
  else if( temp > 0 ) scc->g = 6;

  // Meanwhile the rest is not translated. 

#if 0
	// Compute L_U_2_4, L_W_2_4.
	L_U_2_4_.mx_mult_zero(E_,L_W_1_4_);  
	if( g_ == 4 ) 
	{
		L_U_temp_ = P_U_2_.mx_choose_3(6.0);
		L_U_2_4_ -= L_U_temp_;
	}
	
	L_W_1_4_.delete_data();
	
	L_W_2_4_.mx_mult_zero(ET_,L_U_1_4_); 
	if( g_ == 4 ) 
	{
		L_W_temp_ = P_W_2_.mx_choose_3(6.0);
		L_W_2_4_ -= L_W_temp_;
	}

	// Compute L_U_4_2, L_W_4_2.
	L_U_4_2_.mx_mult_zero(E_,L_W_3_2_);  
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_2_2_); 
	L_U_4_2_ -= L_U_temp_; 
	if( g_ == 4 )
	{
		L_U_4_2_ += P_U_2_c2_;
		L_U_4_2_ += P_U_2_c2_;
	}
	
	L_U_2_2_.delete_data();
	
	L_W_4_2_.mx_mult_zero(ET_,L_U_3_2_); 
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_2_2_); 
	L_W_4_2_ -= L_W_temp_; 
	if( g_ == 4 )
	{
		L_W_4_2_ += P_W_2_c2_;
		L_W_4_2_ += P_W_2_c2_;
	}
	
	L_W_2_2_.delete_data();
	
	// Compute P_U_6, P_W_6.
	P_U_6_.matrix_mult(P_U_5_,ET_); 
	P_U_6_ -= L_U_0_6_; 
	P_U_6_ -= L_U_2_4_; 
	P_U_6_ -= L_U_4_2_;
	L_U_2_4_.delete_data();
	
	P_W_6_.matrix_mult(P_W_5_,E_); 
	P_W_6_ -= L_W_0_6_; 
	P_W_6_ -= L_W_2_4_; 
	P_W_6_ -= L_W_4_2_;
		
	// Compute L_U_1_6, L_W_1_6.
	L_U_1_6_.matrix_mult(E_,L_W_0_6_); 
	L_U_temp_  = P_U_5_;
	L_U_temp_ *= E_;
	L_U_temp_ *= 2.0;
	L_U_1_6_  -= L_U_temp_;
	if( g_ == 4 )
	{
		L_U_temp_  = P_U_3_.mx_choose_2(2.0);
		L_U_temp_ *= E_;
		L_U_1_6_  -= L_U_temp_;
		L_U_temp_.matrix_mult(P_U_2_c2_,E_);
		L_U_temp_ -= P_U_3_;
		L_U_temp_ *= E_;
		L_U_temp_ *= 2.0;
		L_U_1_6_  += L_U_temp_;
		L_U_temp_.matrix_mult(E_,P_W_2_c2_);
		L_U_temp_ -= P_U_3_;
		L_U_temp_ *= E_;
		L_U_temp_ *= 2.0;
		L_U_1_6_  += L_U_temp_;
	}

	L_W_1_6_.matrix_mult(ET_,L_U_0_6_); 
	L_W_temp_  = P_W_5_;
	L_W_temp_ *= ET_;
	L_W_temp_ *= 2.0;
	L_W_1_6_  -= L_W_temp_;
	if( g_ == 4 )
	{
		L_W_temp_  = P_W_3_.mx_choose_2(2.0);
		L_W_temp_ *= ET_;
		L_W_1_6_  -= L_W_temp_;
		L_W_temp_.matrix_mult(P_W_2_c2_,ET_);
		L_W_temp_ -= P_W_3_;
		L_W_temp_ *= ET_;
		L_W_temp_ *= 2.0;
		L_W_1_6_  += L_W_temp_;
		L_W_temp_.matrix_mult(ET_,P_U_2_c2_);
		L_W_temp_ -= P_W_3_;
		L_W_temp_ *= ET_;
		L_W_temp_ *= 2.0;
		L_W_1_6_  += L_W_temp_;
	}

	// Compute L_U_3_4.
	if( g_ == 4 )
	{
		L_U_3_4_.matrix_mult(E_,L_W_2_4_); 
		L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_1_4_); 
		L_U_3_4_  -= L_U_temp_;
		L_U_temp_  = P_U_3_.mx_choose_2(4.0);
		L_U_temp_ *= E_;
		L_U_3_4_  -= L_U_temp_;
		L_U_temp_.matrix_mult(P_U_2_c2_,E_);
		L_U_temp_ -= P_U_3_;
		L_U_temp_ *= E_;
		L_U_temp_ *= 4.0;
		L_U_3_4_  += L_U_temp_;
		L_U_temp_.matrix_mult(E_,P_W_2_c2_);
		L_U_temp_ -= P_U_3_;
		L_U_temp_ *= E_;
		L_U_temp_ *= 6.0;
		L_U_3_4_  += L_U_temp_;
	}
	
	else{ L_U_3_4_ = L_U_1_4_; } // Both 0 if g_ > 4.
	
	L_U_1_4_.delete_data();

	// Compute L_U_5_2 and L_W_5_2.
	L_U_5_2_.matrix_mult(E_,L_W_4_2_);
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_3_2_);
	L_U_5_2_ -= L_U_temp_;
	L_U_temp_ = P_U_5_*E_;
	L_U_5_2_ -= L_U_temp_;
	if( g_ == 4 )
	{
		L_U_temp_  = P_U_3_*E_;
		L_U_5_2_  += L_U_temp_;
		L_U_5_2_  += L_U_temp_; 
		L_U_temp_.matrix_mult(L_U_0_4_,L_U_1_2_);
		L_U_5_2_  -= L_U_temp_;
		L_U_temp_  = L_U_3_2_*E_;
		L_U_5_2_  += L_U_temp_;
		L_W_temp_.matrix_mult(P_U_3_,L_W_0_2_m2_);
		L_U_temp_  = L_W_temp_*E_;
		L_U_5_2_  += L_U_temp_;
		L_U_5_2_  += L_U_temp_;
		L_W_temp_.matrix_mult(P_U_2_c2_,E_);
		L_W_temp_ -= P_U_3_;
		L_U_temp_  = L_W_temp_*E_;
		L_U_5_2_  += L_U_temp_;
		L_U_5_2_  += L_U_temp_;
	}

	P_U_2_c2_.delete_data(); 

	L_W_5_2_.matrix_mult(ET_,L_U_4_2_);
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_3_2_);
	L_W_5_2_ -= L_W_temp_;
	L_W_temp_ = P_W_5_*ET_;
	L_W_5_2_ -= L_W_temp_;
	if( g_ == 4 )
	{
		L_W_temp_.matrix_mult(L_W_0_4_,L_W_1_2_);
		L_W_5_2_  -= L_W_temp_;
		L_W_temp_  = L_W_3_2_;
		L_U_temp_.matrix_mult(P_W_3_,L_U_0_2_m2_);
		L_W_temp_ += L_U_temp_;
		L_W_temp_ += L_U_temp_;
		L_U_temp_.matrix_mult(P_W_2_,ET_);
		L_W_temp_ += L_U_temp_;
		L_W_temp_ += L_U_temp_;
		L_W_5_2_  += L_W_temp_;	
	}	
	
	P_W_2_c2_.delete_data();
	
	// Compute P_U_7 and P_W_7.
	P_U_7_.matrix_mult(P_U_6_,E_); 
	P_U_7_ -= L_U_1_6_; 
	P_U_7_ -= L_U_3_4_; 
	P_U_7_ -= L_U_5_2_;	
	L_U_3_4_.delete_data();
	
	P_W_7_.transpose(P_U_7_);
	
	// Compute L_U_0_8 and L_W_0_8.
	L_U_0_8_.mx_mult_diag(P_U_7_,ET_);
	L_W_0_8_.mx_mult_diag(P_W_7_,E_);
	
	temp = L_U_0_8_.int_trace();
	if( g_ == 4 ) 
	{
		Ng4_ = temp/8;
		L_U_0_8_.diagonal(Ng4_per_u_);
	}
	
	else if( temp && g_ != 6 ) g_ = 8;	 
#endif
}

void scc_count(scc_t *scc)
{
  scc->g   = 1000000;
  scc->Ng  = 0;
  scc->Ng2 = 0;
  scc->Ng4 = 0;
	
  // Count 4 cycles first to determine girth.
  if( scc_count_four_cycles(scc))
    {
      // girth == 4
      scc_count_six_eight_cycles(scc);
    }
  
  else
    {
      ecc_assert(0, "full support of dirth 6 is not fully supported\n");
      // Count 6 and 8 cycles to determine the girth.
    }
}

int tH[] = {1,1,0,
	    1,1,1, 
	    0,1,1};
int_mtx_t H;

void test_scc()
{
  scc_t *scc;
  scc = (scc_t*)ecc_malloc(sizeof(scc_t));
  init_int_mtx_from_array(&H, 3, 3, tH);

  scc_init(scc, &H);
  scc_count(scc);

}

void test_scc2(int_mtx_t* H)
{
  scc_t *scc;
  scc = (scc_t*)ecc_malloc(sizeof(scc_t));
  scc_init(scc, H);
  scc_count(scc);

}
