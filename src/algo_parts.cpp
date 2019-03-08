#include "algo_parts.h"
#include "matching_LUTs.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

// void simple_algo_hw( ap_int<32> inA, ap_int<32> inB, ap_int<32> &outA ) {

//         ap_int<32> offset = 15;
//         outA = inA + inB + offset;

// }

/*
void is_in_boundaries_th(ap_int<THETA_W_SIZE> val, ap_int<THETA_W_SIZE> b_low, ap_int<THETA_W_SIZE> b_high, ap_uint<1> &out_bool)
{
    is_in_boundaries<THETA_W_SIZE>(val, b_low, b_high, out_bool);   
}

void is_in_boundaries_ph(ap_int<PHI_W_SIZE> val, ap_int<PHI_W_SIZE> b_low, ap_int<PHI_W_SIZE> b_high, ap_uint<1> &out_bool)
{
    is_in_boundaries<PHI_W_SIZE>(val, b_low, b_high, out_bool);   
}
*/

bool is_in_boundaries_th(ap_int<THETA_W_SIZE> val, ap_int<THETA_W_SIZE> b_low, ap_int<THETA_W_SIZE> b_high)
{
    return is_in_boundaries<THETA_W_SIZE>(val, b_low, b_high);   
}

bool is_in_boundaries_ph(ap_int<PHI_W_SIZE> val, ap_int<PHI_W_SIZE> b_low, ap_int<PHI_W_SIZE> b_high)
{
    return is_in_boundaries<PHI_W_SIZE>(val, b_low, b_high);   
}


/*
void make_tkmu(ap_int<MU_W_SIZE> in_mu, ap_int<TRK_W_SIZE> in_tk, ap_int<TRK_W_SIZE> &out_tkmu, ap_uint<1> queue_refresh)
{
    // --- decode the information from the inout tracks

    // FIXME: using an integer representation with some kind of LSB
    // is it better to use ap_fixed? How to handle the conversion from the track word in that case?

    // FIXME: replace the assignments below with a proper decoding of the inputs tracks
    ap_uint<PT_W_SIZE> trk_pt = 25;
    ap_uint<PT_W_SIZE> mu_pt  = 35;

    ap_uint<THETA_W_SIZE> trk_theta = 25;
    ap_uint<THETA_W_SIZE> mu_theta  = 35;

    ap_uint<PHI_W_SIZE> trk_phi = 25;
    ap_uint<PHI_W_SIZE> mu_phi  = 35;

    // --- use the track pT information to find all the necessary indices etc

    // FIXME
    ap_uint<PT_LUT_SIZE> pt_lut_addr = 128; 

    ap_uint<1> match_th;
    ap_uint<1> match_ph;

    // --- compute the deltas between these objects

    //ap_int<THETA_W_SIZE> val, ap_int<THETA_W_SIZE> b_low, ap_int<THETA_W_SIZE> b_high, ap_uint<1> &out_bool);
    is_in_boundaries_th()

    // --- 
}
*/


void correlate_single(ap_uint<MU_W_SIZE> mu, ap_uint<TRK_W_SIZE> trk, ap_uint<TKMU_W_SIZE> &tkmu, 
    ap_uint<10> this_phi_low_bound, ap_uint<10> this_phi_high_bound, ap_uint<10> this_theta_low_bound, ap_uint<10> this_theta_high_bound,
    ap_uint<TRK_PT_W_SIZE> &curr_sort_par)
{
    // ---- dtheta
    ap_int<TRK_THETA_W_SIZE> trk_theta = get_trk_theta(trk);
    ap_int<MU_THETA_W_SIZE> mu_theta = get_mu_theta(mu);

    // FIXMEL which precision? All the same? Convert to common one?
    ap_int<MU_THETA_W_SIZE> delta_theta = trk_theta - mu_theta;
    if (delta_theta < 0)
        delta_theta *= -1;


    // ---- dphi
    ap_uint<TRK_PHI_W_SIZE> trk_phi = get_trk_phi(trk);
    ap_uint<MU_PHI_W_SIZE> mu_phi = get_mu_phi(mu);

    // FIXMEL which precision? All the same? Convert to common one?
    // FIXME 2 : need to normalize in -pi, pi
    ap_int<MU_PHI_W_SIZE> delta_phi = trk_phi - mu_phi;
    if (delta_phi < 0)
        delta_phi *= -1;

    if (is_in_boundaries_th(delta_theta, this_theta_low_bound, this_theta_high_bound) &&
        is_in_boundaries_ph(delta_phi, this_phi_low_bound, this_phi_high_bound) )
    {
        ap_uint<TRK_PT_W_SIZE> sort_par = get_trk_pt(trk);
        if (sort_par > curr_sort_par)
        {
            curr_sort_par = sort_par;
            tkmu = trk; // FIXME: a proper reformatting of the values   
        }
    }
}


// NOTE: all the arrays in C are anyway passed by reference, so this function will modify their content
void build_tkmu(ap_uint<MU_W_SIZE> muons[N_MU], ap_uint<TRK_W_SIZE> tracks[N_TRK], ap_uint<TKMU_W_SIZE> tkmus[N_TKMU])
{
    // for (unsigned int itkmu = 0; itkmu < N_TKMU; ++itkmu)
    //     tkmus[itkmu] = itkmu;

    // this is a new event -> reset everything

    // reset the criteria to find the best track
    ap_uint<TRK_PT_W_SIZE> sorting_par[N_TKMU];
    for (size_t i = 0; i < N_TKMU; ++i){
        #pragma HLS unroll
        sorting_par[N_TKMU] = 0;
    }
    
    matching_LUTs mluts;

    // tracks incoming - for each one, pipeline the operations
    for (size_t itrk = 0; itrk < N_TRK; ++itrk)
    {
        #pragma HLS pipeline
        // find the address - FIXME : move address function together with the LUT?
        // -> not sure how to return the two array idxs
        // will using references cause another block to be created?

        // LUTs format : phi_low_bounds [16][512] -> 4 bits x 9 bits
        // FIXME: some saturation may be needed
        // FIXME: theta sign also needed (1 bit)
        ap_uint<4> addr_th;
        ap_uint<4> addr_pt;

        ap_int<TRK_THETA_W_SIZE> abs_trk_theta = get_trk_theta(tracks[itrk]);
        if (abs_trk_theta < 0)
            abs_trk_theta *= -1;

        // FIXME; these boundaries are just to write the code
        // need to change with the actual bin boundaries
        // 5 bits (+sign) -> 32 values and I use 10 bins
        if      (abs_trk_theta < 10) addr_th = 0;
        else if (abs_trk_theta < 12) addr_th = 1;
        else if (abs_trk_theta < 14) addr_th = 2;
        else if (abs_trk_theta < 16) addr_th = 3;
        else if (abs_trk_theta < 18) addr_th = 4;
        else if (abs_trk_theta < 20) addr_th = 5;
        else if (abs_trk_theta < 22) addr_th = 6;
        else if (abs_trk_theta < 24) addr_th = 7;
        else if (abs_trk_theta < 26) addr_th = 8;
        else                         addr_th = 9;

        addr_pt = get_trk_pt(tracks[itrk]) << 1; // since I express pt in 0.5 steps this is easy, just pt*2

        ap_uint<10> this_phi_low_bound    = mluts.phi_low_bounds    [addr_th][addr_pt];
        ap_uint<10> this_phi_high_bound   = mluts.phi_high_bounds   [addr_th][addr_pt];
        ap_uint<10> this_theta_low_bound  = mluts.theta_low_bounds  [addr_th][addr_pt];
        ap_uint<10> this_theta_high_bound = mluts.theta_high_bounds [addr_th][addr_pt];

        // make the correlation - fully parallel
        for (size_t imu = 0; imu < N_MU; ++imu)
        {
            #pragma HLS unroll
            correlate_single(muons[imu], tracks[itrk], tkmus[imu], 
            this_phi_low_bound, this_phi_high_bound, this_theta_low_bound, this_theta_high_bound,
            sorting_par[imu]);
        }
    }
}

void build_one_tkmu (
    ap_uint<MU_W_SIZE> muon, ap_uint<TRK_W_SIZE> track,
    ap_uint<TRK_PT_W_SIZE> in_sorting_par,
    ap_uint<10> phi_low_bound,
    ap_uint<10> phi_high_bound,
    ap_uint<10> theta_low_bound,
    ap_uint<10> theta_high_bound,
    ap_uint<TKMU_W_SIZE> in_tkmu,
    ap_uint<TKMU_W_SIZE> &tkmu, ap_uint<TRK_PT_W_SIZE> &out_sorting_par)
{
    // ---- dtheta
    ap_int<TRK_THETA_W_SIZE> trk_theta = get_trk_theta(track);
    ap_int<MU_THETA_W_SIZE> mu_theta = get_mu_theta(muon);

    // FIXMEL which precision? All the same? Convert to common one?
    ap_int<MU_THETA_W_SIZE> delta_theta = trk_theta - mu_theta;
    if (delta_theta < 0)
        delta_theta *= -1;

    // ---- dphi
    ap_uint<TRK_PHI_W_SIZE> trk_phi = get_trk_phi(track);
    ap_uint<MU_PHI_W_SIZE> mu_phi = get_mu_phi(muon);

    // FIXMEL which precision? All the same? Convert to common one?
    // FIXME 2 : need to normalize in -pi, pi
    ap_int<MU_PHI_W_SIZE> delta_phi = trk_phi - mu_phi;
    if (delta_phi < 0)
        delta_phi *= -1;

    if (is_in_boundaries_th(delta_theta, theta_low_bound, theta_high_bound) &&
        is_in_boundaries_ph(delta_phi, phi_low_bound, phi_high_bound) )
    {
        ap_uint<TRK_PT_W_SIZE> sort_par = get_trk_pt(track);
        if (sort_par > in_sorting_par)
        {
            out_sorting_par = sort_par;
            tkmu = track; // FIXME: here the correct assignment of qualities
        }
        else
        {
            out_sorting_par = in_sorting_par;
            tkmu = in_tkmu;
        }
    }   
}
