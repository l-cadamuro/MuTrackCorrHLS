#include "algo_parts.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

// void simple_algo_hw( ap_int<32> inA, ap_int<32> inB, ap_int<32> &outA ) {

//         ap_int<32> offset = 15;
//         outA = inA + inB + offset;

// }


void is_in_boundaries_th(ap_int<THETA_W_SIZE> val, ap_int<THETA_W_SIZE> b_low, ap_int<THETA_W_SIZE> b_high, ap_uint<1> &out_bool)
{
    is_in_boundaries<THETA_W_SIZE>(val, b_low, b_high, out_bool);   
}