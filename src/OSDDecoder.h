#ifndef __OSD_DECODER__
#define __OSD_DECODER__

#include <stdio.h>
#include <memory.h>
#include "ecc.h"
#include "math_funcs.h"
#include "gf2_algebra.h"


typedef struct osd_decoder_s osd_decoder_t;


// initialize an osd_decoder_t. allocated the memory for *osd_decoder, copied G into osd_decoder (by value).
// pre_calc_error_order is the max order this decoder will be able to calculate (errors are precalculated)
void init_osd_decoder(osd_decoder_t** osd_decoder, int_mtx_t* G, int pre_calc_error_order);

// decode the given word soft_symb. 'osd_order' must be smaller or equal to 'pre_calc_error_order' in init_osd_decoder.
// 'decoded_word' must be allocated to the right size by the user, and will contain the decoded word when returning.
void OSD_decode(osd_decoder_t* osd_decoder, real_vec_t* soft_symb, int osd_order, gf2_vec_t* decoded_word);

// free an osd_decoder_t struct
void free_osd_decoder(osd_decoder_t* osd_decoder);

#endif
