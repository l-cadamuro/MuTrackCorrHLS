#ifndef ALGO_PARTS_H
#define ALGO_PARTS_H

#include "ap_int.h"

// the input data formats properties are defined here
// FIXME: move to a specific file to include here?
#define MU_W_SIZE 100
#define TRK_W_SIZE 100

#define THETA_W_SIZE 6
#define PHI_W_SIZE 6

// check if value is in the interval [b_low, b_high)
template <int N>
void is_in_boundaries(ap_int<N> val, ap_int<N> b_low, ap_int<N> b_high, ap_uint<1> &out_bool){
    if (b_low <= val && val < b_high)
        out_bool = 1;
    else
        out_bool = 0;
}

void is_in_boundaries_th(ap_int<THETA_W_SIZE> val, ap_int<THETA_W_SIZE> b_low, ap_int<THETA_W_SIZE> b_high, ap_uint<1> &out_bool);

// void simple_algo_hw( ap_int<32> inA, ap_int<32> inB, ap_int<32> &outA );

#endif
