#ifndef ALGO_PARTS_H
#define ALGO_PARTS_H

#include "ap_int.h"
#include "dataformats.h"

#define PT_LUT_SIZE 12

/*
// check if value is in the interval [b_low, b_high)
template <int N>
void is_in_boundaries(ap_int<N> val, ap_int<N> b_low, ap_int<N> b_high, ap_uint<1> &out_bool){
    if (b_low <= val && val < b_high)
        out_bool = 1;
    else
        out_bool = 0;
}
*/

template <int N>
bool is_in_boundaries(ap_int<N> val, ap_int<N> b_low, ap_int<N> b_high){
    if (b_low <= val && val < b_high)
        return true;
    else
        return false;
}

// void is_in_boundaries_th(ap_int<THETA_W_SIZE> val, ap_int<THETA_W_SIZE> b_low, ap_int<THETA_W_SIZE> b_high, ap_uint<1> &out_bool);
// void is_in_boundaries_ph(ap_int<PHI_W_SIZE>   val,  ap_int<PHI_W_SIZE>  b_low, ap_int<PHI_W_SIZE>   b_high, ap_uint<1> &out_bool);

bool is_in_boundaries_th(ap_int<THETA_W_SIZE> val, ap_int<THETA_W_SIZE> b_low, ap_int<THETA_W_SIZE> b_high);
bool is_in_boundaries_ph(ap_int<PHI_W_SIZE>   val,  ap_int<PHI_W_SIZE>  b_low, ap_int<PHI_W_SIZE>   b_high);


// void make_tkmu(ap_int<MU_W_SIZE> in_mu, ap_int<TRK_W_SIZE> in_tk, ap_int<TRK_W_SIZE> &out_tkmu, ap_uint<1> queue_refresh); // FIXME: does a better data type exist for the queue_refresh switch?

// void eta_to_theta(ap_int)

// void simple_algo_hw( ap_int<32> inA, ap_int<32> inB, ap_int<32> &outA );

// match a single EMTF to a list of tracks
// void build_one_tkmu (ap_uint<MU_W_SIZE> muons[N_MU], ap_uint<TRK_W_SIZE> tracks[N_TRK], ap_uint<TKMU_W_SIZE> tkmus[N_TKMU]);

void correlate_single(ap_uint<MU_W_SIZE> mu, ap_uint<TRK_W_SIZE> trk, ap_uint<TKMU_W_SIZE> &tkmu, 
    ap_uint<10> this_phi_low_bound, ap_uint<10> this_phi_high_bound, ap_uint<10> this_theta_low_bound, ap_uint<10> this_theta_high_bound,
    ap_uint<TRK_PT_W_SIZE> &curr_sort_par);

// main function that builds the tkmu from the entire list of muons and tracks
void build_tkmu(ap_uint<MU_W_SIZE> muons[N_MU], ap_uint<TRK_W_SIZE> tracks[N_TRK], ap_uint<TKMU_W_SIZE> tkmus[N_TKMU]);
void build_one_tkmu (
    ap_uint<MU_W_SIZE> muon, ap_uint<TRK_W_SIZE> track,
    ap_uint<TRK_PT_W_SIZE> in_sorting_par,
    ap_uint<10> phi_low_bound,
    ap_uint<10> phi_high_bound,
    ap_uint<10> theta_low_bound,
    ap_uint<10> theta_high_bound,
    ap_uint<TKMU_W_SIZE> in_tkmu,
    ap_uint<TKMU_W_SIZE> &tkmu, ap_uint<TRK_PT_W_SIZE> &out_sorting_par);


#endif
