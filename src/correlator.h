#ifndef CORRELATOR_H
#define CORRELATOR_H

// #define AP_INT_MAX_W 2048
// #include "dataformats.h"
#include "algo_parts.h"
#include "ap_int.h"

// 2048 bit width in the UF setup -> keep 15 tracks x 100 bits + 12 EMTF x 40 bits = 1980 bits in input
// starting easy with no pipelining of the event information
// to make the pipeline, have the input word contain various tracks but the same EMTFs for several addresses
void correlator (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info[100], ap_uint<TKMU_W_SIZE*N_TKMU> out_info[100],
    ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
    ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
    ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2);

// no loop on bx inside
void correlator_one (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info,
    ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
    ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
    ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2);

// just copy N_MU tracks to tkmu output collection - same interface as correlator_one for easy debug
void passthrough (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info,
    ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
    ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
    ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2);

#endif
