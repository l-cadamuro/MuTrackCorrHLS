// this file contains a series of unmaintained functins, code, etc... that I use to test and learn HLS
// not to be used as part of the MuTk algorithm

#ifndef PLAYGROUND_H
#define PLAYGROUND_H

#define AP_INT_MAX_W 2048

#include "ap_fixed.h"
#include "ap_int.h"

// the input data formats properties are defined here
// FIXME: move to a specific file to include here?
#define MU_W_SIZE 40
#define TRK_W_SIZE 100
#define TKMU_W_SIZE 100

// size of input and output collections
#define N_MU 12
#define N_TRK 15
#define N_TKMU 12 // by default, the same as N_MU


void func_with_inner_val(ap_int<8> v_in, ap_int<8> &v_out);
void func_with_inner_val_present(ap_int<8> v_in, ap_int<8> &v_out, ap_uint<1> is_first);

void sum_array(ap_int<16> arr_in[50], ap_int<16> &sum_out);

void pipelined_transfer(ap_uint<16> in_info, ap_uint<16> &out_info);
void pipelined_transfer_corrdf(
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info,
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> &out_info);

// void pipelined_transfer_l(ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info,
//     ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
//     ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
//     ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2);

#endif