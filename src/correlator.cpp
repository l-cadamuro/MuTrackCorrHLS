#include "correlator.h"
#include "matching_LUTs.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

void correlator (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info[100], ap_uint<TKMU_W_SIZE*N_TKMU> out_info[100],
    ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
    ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
    ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2)

{
    for (unsigned int ibx = 0; ibx < 100; ++ibx)
    {
        // spy_trk1 = ibx;
        // spy_trk2 = ibx+1;
        // spy_mu1  = ibx;
        // spy_mu2  = ibx+1;

        // spy_tkmu1 = 100+ibx;
        // spy_tkmu2 = 100+ibx+1;

        // original function body
        spy_mu1   = in_info[ibx].range(MU_W_SIZE-1, 0);
        spy_mu2   = in_info[ibx].range(MU_W_SIZE*2-1, MU_W_SIZE);
        spy_trk1  = in_info[ibx].range(TRK_W_SIZE-1 + MU_W_SIZE*N_MU, MU_W_SIZE*N_MU);
        spy_trk2  = in_info[ibx].range(TRK_W_SIZE*2-1 + MU_W_SIZE*N_MU, TRK_W_SIZE + MU_W_SIZE*N_MU);

        // dispatch all the emtf info to the output - just a passthrough
        for (unsigned int iemtf = 0; iemtf < N_MU; ++iemtf)
        {
            // ap_uint<TKMU_W_SIZE> my_tkmu;
            // NOTE that I fill the output according the TKMU word size
            out_info[ibx].range((iemtf+1)*TKMU_W_SIZE-1, iemtf*TKMU_W_SIZE) = in_info[ibx].range((iemtf+1)*MU_W_SIZE-1, iemtf*MU_W_SIZE);
        }
        spy_tkmu1 = out_info[ibx].range(TKMU_W_SIZE-1, 0);
        spy_tkmu2 = out_info[ibx].range(TKMU_W_SIZE*2-1, TKMU_W_SIZE);


        // some more debugging
        /*
        spy_mu1   = in_info[ibx].range(MU_W_SIZE-1, 0);
        spy_mu2   = in_info[ibx].range(MU_W_SIZE*2-1, MU_W_SIZE);
        spy_trk1  = in_info[ibx].range(TRK_W_SIZE-1 + MU_W_SIZE*N_MU, MU_W_SIZE*N_MU);
        spy_trk2  = in_info[ibx].range(TRK_W_SIZE*2-1 + MU_W_SIZE*N_MU, TRK_W_SIZE + MU_W_SIZE*N_MU);


        // dispatch all the emtf info to the output - just a passthrough
        for (unsigned int iemtf = 0; iemtf < N_MU; ++iemtf)
        {
            // ap_uint<TKMU_W_SIZE> my_tkmu;
            // NOTE that I fill the output according the TKMU word size
            out_info[ibx].range((iemtf+1)*TKMU_W_SIZE-1, iemtf*TKMU_W_SIZE) = ibx;
        }
        spy_tkmu1 = out_info[ibx].range(TKMU_W_SIZE-1, 0);
        spy_tkmu2 = out_info[ibx].range(TKMU_W_SIZE*2-1, TKMU_W_SIZE);
        */
    
    }
}


void correlator_one (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info,
    ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
    ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
    ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2)
{
    // unpack my input DF into the various muons
    ap_uint<MU_W_SIZE> muons[N_MU]; // input to corr
    ap_uint<TRK_W_SIZE> tracks[N_TRK]; // input to corr
    ap_uint<TKMU_W_SIZE> tkmus[N_TKMU]; // output to corr

    for (size_t imu = 0; imu < N_MU; ++imu)
    {
        #pragma HLS unroll
        muons[imu] = in_info.range(MU_W_SIZE*(imu+1)-1, MU_W_SIZE*imu);
    }

    for (size_t itrk = 0; itrk < N_TRK; ++itrk)
    {
        #pragma HLS unroll
        tracks[itrk] = in_info.range(TRK_W_SIZE*(itrk+1)-1 + MU_W_SIZE*N_MU, TRK_W_SIZE*itrk + MU_W_SIZE*N_MU);
    }

    // attach spy signals
    spy_mu1 = muons[0];
    spy_mu2 = muons[1];
    spy_trk1 = tracks[0];
    spy_trk2 = tracks[1];

    // launch the correlator
    // build_tkmu(muons, tracks, tkmus);

    // ----------------------------

    // the matchign LUTs
    matching_LUTs mluts;

    // prepare the sorting parameter criteria
    ap_uint<TRK_PT_W_SIZE> sorting_par[N_TKMU];
    for (size_t i = 0; i < N_TKMU; ++i){
        #pragma HLS unroll
        sorting_par[N_TKMU] = 0;
    }

    // make a pipeline on the input tracks
    for (size_t itrk = 0; itrk < N_TRK; ++itrk)
    {
        #pragma HLS pipeline

        // build the matching LUT address
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
            build_one_tkmu(
                muons[imu], tracks[itrk],
                sorting_par[imu],
                this_phi_low_bound,
                this_phi_high_bound,
                this_theta_low_bound,
                this_theta_high_bound,
                tkmus[imu],
                tkmus[imu], sorting_par[imu]
            );
        }

        /*
        build_one_tkmu(
            muons[0], tracks[itrk],
            sorting_par[0],
            this_phi_low_bound,
            this_phi_high_bound,
            this_theta_low_bound,
            this_theta_high_bound,
            tkmus[0],
            tkmus[0], sorting_par[0]
        );
        */
    }

    // ----------------------------

    // attach spy signals
    spy_tkmu1 = tkmus[0];
    spy_tkmu2 = tkmus[1];

    // repack the tkmus inside the single output word
    for (size_t itkmu = 0; itkmu < N_MU; ++itkmu)
    {
        #pragma HLS unroll
        out_info.range(TKMU_W_SIZE*(itkmu+1)-1, TKMU_W_SIZE*itkmu) = tkmus[itkmu];
    }
}


// void correlator_one (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info,
//     ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
//     ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
//     ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2)
// {
//     // #pragma HLS latency min=2 max=2

//     // define here my big LUT for the matching
//     // 12 input bits (3 theta + 9 pt) , 11 output bits
//     // 13 input bits (4 theta + 9 pt) , 10 output bits    
//     ap_uint<10> LUT_th_up[8192];
//     ap_uint<10> LUT_th_do[8192];
//     ap_uint<10> LUT_phi_up[8192];
//     ap_uint<10> LUT_phi_do[8192];

//     // setup the LUTs - dummy for now
//     // NOTE: values will restart from 0 after 2048 because the word size is only 11 bits
//     for (unsigned int ilut = 0; ilut < 8192; ++ilut)
//     {
//         LUT_th_up[ilut]  = ilut;
//         LUT_th_do[ilut]  = ilut;
//         LUT_phi_up[ilut] = ilut;
//         LUT_phi_do[ilut] = ilut;
//         // printf("MAKING LUT idx = ilut --> %i %i %i %i\n", ilut, LUT_th_up[ilut].to_uint(), LUT_th_do[ilut].to_uint(), LUT_phi_up[ilut].to_uint(), LUT_phi_do[ilut].to_uint());
//     }

//     // spy_trk1 = ibx;
//     // spy_trk2 = ibx+1;
//     // spy_mu1  = ibx;
//     // spy_mu2  = ibx+1;

//     // spy_tkmu1 = 100+ibx;
//     // spy_tkmu2 = 100+ibx+1;

//     // original function body
//     spy_mu1   = in_info.range(MU_W_SIZE-1, 0);
//     spy_mu2   = in_info.range(MU_W_SIZE*2-1, MU_W_SIZE);
//     spy_trk1  = in_info.range(TRK_W_SIZE-1 + MU_W_SIZE*N_MU, MU_W_SIZE*N_MU);
//     spy_trk2  = in_info.range(TRK_W_SIZE*2-1 + MU_W_SIZE*N_MU, TRK_W_SIZE + MU_W_SIZE*N_MU);

//     // dispatch all the emtf info to the output - just a passthrough
//     for (unsigned int iemtf = 0; iemtf < N_MU; ++iemtf)
//     {
//         #pragma HLS unroll
//         // ap_uint<TKMU_W_SIZE> my_tkmu;
//         // NOTE that I fill the output according the TKMU word size
        
//         // out_info.range((iemtf+1)*TKMU_W_SIZE-1, iemtf*TKMU_W_SIZE) = in_info.range((iemtf+1)*MU_W_SIZE-1, iemtf*MU_W_SIZE) + 1;

//         // a simple sketch case for the usage of the LUTs
        
//         unsigned int iaddr = (iemtf < 6 ? iemtf : 4000 + iemtf); // not sure if +4000 is needed, I want to avoid that HLS trims my array
//         ap_uint<11>  val1   = LUT_th_up[iaddr];
//         ap_uint<11>  val2   = LUT_th_do[iaddr];
//         ap_uint<11>  val3   = LUT_phi_up[iaddr];
//         ap_uint<11>  val4   = LUT_phi_do[iaddr];

//         // printf("@@@@@@ addr %i for EMTF %i gives %i %i %i %i\n", iaddr, iemtf, val1.to_uint(), val2.to_uint(), val3.to_uint(), val4.to_uint());
//         out_info.range((iemtf+1)*TKMU_W_SIZE-1, iemtf*TKMU_W_SIZE) = in_info.range((iemtf+1)*MU_W_SIZE-1, iemtf*MU_W_SIZE) + val1 + val2 + val3 + val4;
//     }
//     // spy_tkmu1 = out_info.range(TKMU_W_SIZE-1, 0);
//     // spy_tkmu2 = out_info.range(TKMU_W_SIZE*2-1, TKMU_W_SIZE);

//     // will this work as a register and effecively insert a delay of 1 clk?
//     // answer: no it does not :-( ==> because order is wrong!
//     // static ap_uint<TKMU_W_SIZE> reg_tkmu_1;
//     // static ap_uint<TKMU_W_SIZE> reg_tkmu_2;
//     // reg_tkmu_1 = out_info.range(TKMU_W_SIZE-1, 0);
//     // reg_tkmu_2 = out_info.range(TKMU_W_SIZE*2-1, TKMU_W_SIZE);
//     // spy_tkmu1 = reg_tkmu_1;
//     // spy_tkmu2 = reg_tkmu_2;

//     // trying reverse order
//     // will this work as a register and effecively insert a delay of 1 clk?
//     // answer: this works!
//     static ap_uint<TKMU_W_SIZE> reg_tkmu_1;
//     static ap_uint<TKMU_W_SIZE> reg_tkmu_2;
//     spy_tkmu1  = reg_tkmu_1;
//     spy_tkmu2  = reg_tkmu_2;
//     reg_tkmu_1 = out_info.range(TKMU_W_SIZE-1, 0);
//     reg_tkmu_2 = out_info.range(TKMU_W_SIZE*2-1, TKMU_W_SIZE);


//     /*
//     // so let's try a pipeline
//     // seems to work. But NB: variables stay on the ILA for more clock cycles
//     // it's not a full pipeline, the assignments get resolved in multiple clock cycles
//     const int npipe = 3;
//     static ap_uint<TKMU_W_SIZE> pipe_tkmu_1 [npipe];
//     static ap_uint<TKMU_W_SIZE> pipe_tkmu_2 [npipe];

//     spy_tkmu1 = pipe_tkmu_1[npipe-1];
//     spy_tkmu2 = pipe_tkmu_2[npipe-1];
//     for (unsigned int ipipe = npipe-1; ipipe > 0; --ipipe)
//     {
//         pipe_tkmu_1[ipipe] = pipe_tkmu_1[ipipe-1];
//         pipe_tkmu_2[ipipe] = pipe_tkmu_2[ipipe-1];
//     }
//     pipe_tkmu_1[0] = out_info.range(TKMU_W_SIZE-1, 0);
//     pipe_tkmu_2[0] = out_info.range(TKMU_W_SIZE*2-1, TKMU_W_SIZE);
//     */

//     // some more debugging
//     /*
//     spy_mu1   = in_info.range(MU_W_SIZE-1, 0);
//     spy_mu2   = in_info.range(MU_W_SIZE*2-1, MU_W_SIZE);
//     spy_trk1  = in_info.range(TRK_W_SIZE-1 + MU_W_SIZE*N_MU, MU_W_SIZE*N_MU);
//     spy_trk2  = in_info.range(TRK_W_SIZE*2-1 + MU_W_SIZE*N_MU, TRK_W_SIZE + MU_W_SIZE*N_MU);


//     // dispatch all the emtf info to the output - just a passthrough
//     for (unsigned int iemtf = 0; iemtf < N_MU; ++iemtf)
//     {
//         // ap_uint<TKMU_W_SIZE> my_tkmu;
//         // NOTE that I fill the output according the TKMU word size
//         out_info.range((iemtf+1)*TKMU_W_SIZE-1, iemtf*TKMU_W_SIZE) = ibx;
//     }
//     spy_tkmu1 = out_info.range(TKMU_W_SIZE-1, 0);
//     spy_tkmu2 = out_info.range(TKMU_W_SIZE*2-1, TKMU_W_SIZE);
//     */
    
// }