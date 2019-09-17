#include "correlator.h"
#include "matching_LUTs.h"
#include "matcher.h"
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

    // output array must be clean
    for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
    {
        #pragma HLS unroll
        tkmus[itkmu] = 0x0; // FIXME: shoudl this be done at the testbench level instead?
    }

    #pragma HLS ARRAY_PARTITION variable=muons  complete
    #pragma HLS ARRAY_PARTITION variable=tracks complete
    #pragma HLS ARRAY_PARTITION variable=tkmus  complete

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
    #pragma HLS ARRAY_PARTITION variable=sorting_par complete

    // ap_uint<TRK_PT_W_SIZE> sorting_par_4debug[N_TKMU];
    // #pragma HLS ARRAY_PARTITION variable=sorting_par_4debug complete

    
    for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu){
        #pragma HLS unroll
        sorting_par[itkmu] = 0;
    }


    // make a pipeline on the input tracks
    for (size_t itrk = 0; itrk < N_TRK; ++itrk)
    {
        #pragma HLS pipeline

        // build the matching LUT address
        ap_uint<4> addr_th;
        ap_uint<9> addr_pt;

        ap_uint<TRK_THETA_W_SIZE> abs_trk_theta = get_trk_theta(tracks[itrk]);
        // ap_int<TRK_THETA_W_SIZE> abs_trk_theta = get_trk_theta(tracks[itrk]);
        // if (abs_trk_theta < 0)
        //     abs_trk_theta *= -1;

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

        addr_pt = get_trk_pt(tracks[itrk]); // since I express pt in 0.5 steps this is easy, the address is the pt

        ap_uint<10> this_phi_low_bound    = mluts.phi_low_bounds    [addr_th][addr_pt];
        ap_uint<10> this_phi_high_bound   = mluts.phi_high_bounds   [addr_th][addr_pt];
        ap_uint<10> this_theta_low_bound  = mluts.theta_low_bounds  [addr_th][addr_pt];
        ap_uint<10> this_theta_high_bound = mluts.theta_high_bounds [addr_th][addr_pt];

        // if (get_trk_pt(tracks[itrk]) > 20 && get_trk_pt(tracks[itrk]) < 100)
        // {
        //     printf("=== debugging the LUT lookup in the 20 - 100 hw pt gap ===\n");
        //     printf(" -- trk theta = %10u  trk pt = %10u \n", get_trk_theta(tracks[itrk]).to_uint64(), get_trk_pt(tracks[itrk]).to_uint64());
        //     printf(" -- addr theta = %10u addr pt = %10u \n", addr_th.to_uint64(), addr_pt.to_uint64());
        //     printf(" -- philow   = %10u, phiup   = %10u\n",
        //         this_phi_low_bound.to_uint64(),
        //         this_phi_high_bound.to_uint64());
        //     printf(" -- thetalow = %10u, thetaup = %10u\n",
        //         this_theta_low_bound.to_uint64(),
        //         this_theta_high_bound.to_uint64());
        // }


        // #ifndef __SYNTHESIS__
        // printf("=== debugging the main LUT lookup ===\n");
        // printf(" -- trk theta = %10u  trk pt = %10u \n", get_trk_theta(tracks[itrk]).to_uint64(), get_trk_pt(tracks[itrk]).to_uint64());
        // printf(" -- addr theta = %10u addr pt = %10u \n", addr_th.to_uint64(), addr_pt.to_uint64());
        // printf(" -- philow   = %10u, phiup   = %10u\n",
        //     this_phi_low_bound.to_uint64(),
        //     this_phi_high_bound.to_uint64());
        // printf(" -- thetalow = %10u, thetaup = %10u\n",
        //     this_theta_low_bound.to_uint64(),
        //     this_theta_high_bound.to_uint64());
        // printf(" >>> dump of the full LUT [0] [0-10]\n");
        // for (size_t il = 0; il < 10; ++il)
        //     printf("%i .. %10u\n", il, mluts.phi_low_bounds[0][il].to_uint64());
        // printf(" >>> dump of the full LUT [15] [0-10]\n");
        // for (size_t il = 0; il < 10; ++il)
        //     printf("%i .. %10u\n", il, mluts.phi_low_bounds[15][il].to_uint64());
        // #endif

        // // make the correlation - fully parallel
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

        // for (size_t imu = 0; imu < N_MU; ++imu)
        // {
        //     #pragma HLS unroll
        //     build_one_tkmu(
        //         muons[imu], tracks[itrk],
        //         sorting_par[imu],
        //         this_phi_low_bound,
        //         this_phi_high_bound,
        //         this_theta_low_bound,
        //         this_theta_high_bound,
        //         tkmus[imu],
        //         tkmus[imu], sorting_par_4debug[imu]
        //     );
        // }



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

        build_one_tkmu(
            muons[1], tracks[itrk],
            sorting_par[1],
            this_phi_low_bound,
            this_phi_high_bound,
            this_theta_low_bound,
            this_theta_high_bound,
            tkmus[1],
            tkmus[1], sorting_par[1]
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

    // // debug output -- this works
    // for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
    // {
    //     #pragma HLS unroll
    //     out_info.range(TKMU_W_SIZE*(itkmu+1)-1, TKMU_W_SIZE*itkmu) = itkmu;
    // }

}


void arr_correlator_one (
    ap_uint<TRK_W_SIZE> tracks [N_TRK], ap_uint<MU_W_SIZE> muons [N_MU], // inputs
    ap_uint<TKMU_W_SIZE> tkmus [N_TKMU]) // outputs
{
    #pragma HLS interface ap_ctrl_none port=return // not sure if needed at all?

    // // output array must be clean
    // for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
    // {
    //     #pragma HLS unroll
    //     tkmus[itkmu] = 0x0; // FIXME: shoudl this be done at the testbench level instead?
    // }

    // #pragma HLS ARRAY_PARTITION variable=muons  complete
    // #pragma HLS ARRAY_PARTITION variable=tracks complete
    // #pragma HLS ARRAY_PARTITION variable=tkmus  complete

    // launch the correlator
    // build_tkmu(muons, tracks, tkmus);

    // ----------------------------

    // the matchign LUTs
    matching_LUTs mluts;

    // prepare the sorting parameter criteria
    ap_uint<TRK_PT_W_SIZE> sorting_par[N_TKMU];
    #pragma HLS ARRAY_PARTITION variable=sorting_par complete
    
    for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu){
        #pragma HLS unroll
        sorting_par[itkmu] = 0;
    }

    // make a pipeline on the input tracks
    for (size_t itrk = 0; itrk < N_TRK; ++itrk)
    {
        #pragma HLS pipeline

        // build the matching LUT address
        ap_uint<4> addr_th;
        ap_uint<9> addr_pt;

        ap_uint<TRK_THETA_W_SIZE> abs_trk_theta = get_trk_theta(tracks[itrk]);
        // ap_int<TRK_THETA_W_SIZE> abs_trk_theta = get_trk_theta(tracks[itrk]);
        // if (abs_trk_theta < 0)
        //     abs_trk_theta *= -1;

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

        addr_pt = get_trk_pt(tracks[itrk]); // since I express pt in 0.5 steps this is easy, the address is the pt

        ap_uint<10> this_phi_low_bound    = mluts.phi_low_bounds    [addr_th][addr_pt];
        ap_uint<10> this_phi_high_bound   = mluts.phi_high_bounds   [addr_th][addr_pt];
        ap_uint<10> this_theta_low_bound  = mluts.theta_low_bounds  [addr_th][addr_pt];
        ap_uint<10> this_theta_high_bound = mluts.theta_high_bounds [addr_th][addr_pt];

        // // make the correlation - fully parallel
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
    }
}



// handles one track arriving from each sector
void arr_correlator_mult (
    ap_uint<TRK_W_SIZE> all_tracks [N_TRK_SECTORS][N_TRK], ap_uint<MU_W_SIZE> muons [N_MU], // inputs
    ap_uint<TKMU_W_SIZE> tkmus [N_TKMU])
{
 
    ap_uint<TKMU_W_SIZE> all_tkmus [N_TRK_SECTORS][N_TKMU]; // N_TRK_SECTOR copies of the output tkmu
    
    // init to zeros all the tkmus
    for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec){
        for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu){
            #pragma HLS unroll
            all_tkmus[isec][itkmu] = 0x0;
        }
    }

    // distribute the correlator over the sectors
    for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
    {
        #pragma HLS unroll
        arr_correlator_one (all_tracks[isec], muons, all_tkmus[isec]);
    }

    // sort the best tkmu out of all the sectors
    // NOTE: using here as sorting parameter the best pT, must be aligned with what is in the tkmu correlator_one definition

    for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
    {
        #pragma HLS unroll // parallel over all tkmu, since each position corresponds to a different standalone muon

        for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
        {
            #pragma HLS pipeline // pipelined over the sectors -> will take N_TRK_SECTORS clocks to do. NOTE: a better sorting could be implemented (e.g. bitonic)
            if ( get_tkmu_pt(all_tkmus[isec][itkmu]) > get_tkmu_pt(tkmus[itkmu]) )
                tkmus[itkmu] = all_tkmus[isec][itkmu];
        }
    }   
}


void top_arr_correlator_one (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info)
{
    // unpack my input DF into the various muons
    ap_uint<MU_W_SIZE> muons[N_MU]; // input to corr
    ap_uint<TRK_W_SIZE> tracks[N_TRK]; // input to corr
    ap_uint<TKMU_W_SIZE> tkmus[N_TKMU]; // output to corr

    // output array must be clean
    for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
    {
        #pragma HLS unroll
        tkmus[itkmu] = 0x0; // FIXME: shoudl this be done at the testbench level instead?
    }

    #pragma HLS ARRAY_PARTITION variable=muons  complete
    #pragma HLS ARRAY_PARTITION variable=tracks complete
    #pragma HLS ARRAY_PARTITION variable=tkmus  complete

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

    // run the correlator
    arr_correlator_one(tracks, muons, tkmus);

    // repack the tkmus inside the single output word
    for (size_t itkmu = 0; itkmu < N_MU; ++itkmu)
    {
        #pragma HLS unroll
        out_info.range(TKMU_W_SIZE*(itkmu+1)-1, TKMU_W_SIZE*itkmu) = tkmus[itkmu];
    }
}

// void top_arr_correlator_mult (ap_uint<TRK_W_SIZE*N_TRK*N_TRK_SECTORS + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info)
// {
//     // unpack my input DF into the various muons
//     ap_uint<MU_W_SIZE> muons[N_MU]; // input to corr
//     ap_uint<TRK_W_SIZE> all_tracks[N_TRK_SECTORS][N_TRK]; // input to corr
//     ap_uint<TKMU_W_SIZE> tkmus[N_TKMU]; // output to corr

//     // output array must be clean
//     for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
//     {
//         #pragma HLS unroll
//         tkmus[itkmu] = 0x0; // FIXME: shoudl this be done at the testbench level instead?
//     }

//     #pragma HLS ARRAY_PARTITION variable=muons  complete
//     #pragma HLS ARRAY_PARTITION variable=all_tracks complete
//     #pragma HLS ARRAY_PARTITION variable=tkmus  complete

//     for (size_t imu = 0; imu < N_MU; ++imu)
//     {
//         #pragma HLS unroll
//         muons[imu] = in_info.range(MU_W_SIZE*(imu+1)-1, MU_W_SIZE*imu);
//     }

//     for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
//         for (size_t itrk = 0; itrk < N_TRK; ++itrk)
//         {
//             #pragma HLS unroll
//             all_tracks[isec][itrk] = in_info.range(TRK_W_SIZE*(itrk+1)-1 + TRK_W_SIZE*N_TRK*isec + MU_W_SIZE*N_MU, TRK_W_SIZE*itrk + TRK_W_SIZE*N_TRK*isec + MU_W_SIZE*N_MU);
//         }

//     // run the correlator    
//     arr_correlator_mult(all_tracks, muons, tkmus);

//     // repack the tkmus inside the single output word
//     for (size_t itkmu = 0; itkmu < N_MU; ++itkmu)
//     {
//         #pragma HLS unroll
//         out_info.range(TKMU_W_SIZE*(itkmu+1)-1, TKMU_W_SIZE*itkmu) = tkmus[itkmu];
//     }
// }


/*
void correlator_mult (ap_uint<TRK_W_SIZE*N_TRK*N_TRK_SECTORS + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info,
    ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
    ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
    ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2)
{
    // unpack the big track word into single ones
    // unpack my input DF into the various muons
    ap_uint<MU_W_SIZE>   muons[N_MU];                  // input to corr
    ap_uint<TRK_W_SIZE>  tracks[N_TRK_SECTORS][N_TRK]; // input to corr
    ap_uint<TKMU_W_SIZE> tkmus[N_TRK_SECTORS][N_TKMU]; // output to corr - there will be one list for each sector

    #pragma HLS ARRAY_PARTITION variable=muons        complete
    #pragma HLS ARRAY_PARTITION variable=tracks dim=0 complete
    #pragma HLS ARRAY_PARTITION variable=tkmus  dim=0 complete

    for (size_t imu = 0; imu < N_MU; ++imu)
    {
        #pragma HLS unroll
        muons[imu] = in_info.range(MU_W_SIZE*(imu+1)-1, MU_W_SIZE*imu);
    }

    for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
    {
        for (size_t itrk = 0; itrk < N_TRK; ++itrk)
        {
            #pragma HLS unroll
            tracks[isec][itrk] = in_info.range(TRK_W_SIZE*(itrk+1)-1 + TRK_W_SIZE*N_TRK*isec + MU_W_SIZE*N_MU, TRK_W_SIZE*itrk + TRK_W_SIZE*N_TRK*isec + MU_W_SIZE*N_MU);
        }
    }

    // now repack all into separate arrays to be fed to the correlator
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info_local[N_TRK_SECTORS];
    ap_uint<TKMU_W_SIZE*N_TKMU> out_info_local[N_TRK_SECTORS];

    for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
    {
        #pragma HLS unroll

    }

    // // attach spy signals
    // spy_mu1 = muons[0];
    // spy_mu2 = muons[1];
    // spy_trk1 = tracks[0];
    // spy_trk2 = tracks[1];    
}
*/





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


// void arr_correlator_mult_matcher (
//     ap_uint<TRK_W_SIZE>  all_tracks [N_TRK_SECTORS], // inputs : one track incoming per sector
//     ap_uint<MU_W_SIZE>   muons [N_MU],   // inputs : all the muons ready at the same time
//     ap_uint<TKMU_W_SIZE> tkmus [N_TKMU], // output TkMus
//     ap_uint<1>           new_bx
// )

// {
//     // define as many matchers as sectors x muons 
//     matcher match_units[N_TKMU];

//     // define the maatchign LUTs, one per incoming track (for the lookup)
//     matching_LUTs mluts[N_TRK_SECTORS];


//     if (new_bx) // re-init everything
//     {
//         for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
//         {
//             #pragma HLS unroll
//             match_units[itkmu].init(muons[imu]); /// NOTE!! This implicitly requires N_TKMU = N_MU
//         }
//     }

//     for (size_t itrk = 0; itrk < N_TRK_SECTORS; ++itrk)
//     {
//         #pragma HLS unroll
//         // build the matching LUT address
//         ap_uint<4> addr_th;
//         ap_uint<9> addr_pt;

//         ap_uint<TRK_THETA_W_SIZE> abs_trk_theta = all_tracks[itrk].theta;
//         // ap_int<TRK_THETA_W_SIZE> abs_trk_theta = get_trk_theta(tracks[itrk]);
//         // if (abs_trk_theta < 0)
//         //     abs_trk_theta *= -1;

//         // FIXME; these boundaries are just to write the code
//         // need to change with the actual bin boundaries
//         // 5 bits (+sign) -> 32 values and I use 10 bins
        
//         if      (abs_trk_theta < 10) addr_th = 0;
//         else if (abs_trk_theta < 12) addr_th = 1;
//         else if (abs_trk_theta < 14) addr_th = 2;
//         else if (abs_trk_theta < 16) addr_th = 3;
//         else if (abs_trk_theta < 18) addr_th = 4;
//         else if (abs_trk_theta < 20) addr_th = 5;
//         else if (abs_trk_theta < 22) addr_th = 6;
//         else if (abs_trk_theta < 24) addr_th = 7;
//         else if (abs_trk_theta < 26) addr_th = 8;
//         else                         addr_th = 9;

//         addr_pt = all_tracks[itrk].pt; // since I express pt in 0.5 steps this is easy, the address is the pt

//         ap_uint<10> this_phi_low_bound    = mluts.phi_low_bounds    [addr_th][addr_pt];
//         ap_uint<10> this_phi_high_bound   = mluts.phi_high_bounds   [addr_th][addr_pt];
//         ap_uint<10> this_theta_low_bound  = mluts.theta_low_bounds  [addr_th][addr_pt];
//         ap_uint<10> this_theta_high_bound = mluts.theta_high_bounds [addr_th][addr_pt];

//         for (size_t imatcher = 0; imatcher < N_TKMU; ++imatcher)
//         {
//             #pragma HLS unroll
//             match_units[imatcher].process(
//                 all_tracks[itrk],
//                 this_phi_low_bound,
//                 this_phi_high_bound,
//                 this_theta_low_bound,
//                 this_theta_high_bound
//             );
//         }
//     }

//     // the current tkmu candidate gets copied to the output line
//     for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
//     {
//         #pragma HLS unroll
//         tkmus[itkmu] = match_units[itkmu].tkmu_; /// NOTE!! This implicitly requires N_TKMU = N_MU
//     }
// }

void correlator_stream (
    muon_t                in_muons  [N_MU],          // inputs : all the muons ready at the same time
    hls::stream<track_t>  in_tracks [N_TRK_SECTORS], // inputs : one track incoming per sector as a stream
    tkmu_t                out_tkmu  [N_TKMU],        // output TkMus
    ap_uint<1>            &corr_done)        
{
    // ----------------------------------------
    // -- definitions of blocks --

    #pragma HLS data_pack variable=in_muons
    #pragma HLS data_pack variable=in_tracks
    #pragma HLS data_pack variable=out_tkmu

    // define the matching units. They must be persistent outside the function -> static
    static matcher match_units[N_TKMU];
    // independent units -> partition the array
    #pragma HLS ARRAY_PARTITION variable=match_units  complete

    // one set of matching LUT per incoming track
    matching_LUTs mluts[N_TRK_SECTORS];
    #pragma HLS ARRAY_PARTITION variable=mluts  complete

    // and output buffer decoupled from the matcher - so that the data is persistent
    static tkmu_t     out_buffer [N_TKMU];
    static ap_uint<1> corr_done_buffer;
    #pragma HLS ARRAY_PARTITION variable=out_buffer  complete
    // connect output to the buffers
    for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
    {
        #pragma HLS unroll
        out_tkmu[itkmu] = out_buffer[itkmu];
    }
    corr_done = corr_done_buffer;

    // ----------------------------------------
    // -- track processing --
    // process tracks if the stream is full. FIXME: is it OK to check just the first stream if empty?
    // Expect all streams to be filled in the same way, even with empty candidates
    if (!(in_tracks[0]).empty())
    {
        // #ifndef __SYNTHESIS__
        // std::cout << " ... TRACK STREAM HAS DATA, PROCESSING ..." << std::endl;
        // #endif

        for (size_t itrk = 0; itrk < N_TRK_SECTORS; ++itrk)
        {
            #pragma HLS unroll
            
            track_t trk = in_tracks[itrk].read();

            // std::cout << " ....... sector " << itrk << " with track pt = " << trk.pt.to_uint64() << " theta = " << trk.theta.to_uint64() << " phi = " << trk.phi.to_uint64() << std::endl;

            // build the matching LUT address
            ap_uint<4> addr_th;
            ap_uint<9> addr_pt;

            ap_uint<TRK_THETA_W_SIZE> abs_trk_theta = trk.theta;
            // ap_int<TRK_THETA_W_SIZE> abs_trk_theta = get_trk_theta(tracks[itrk]);
            // if (abs_trk_theta < 0)
            //     abs_trk_theta *= -1;

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

            addr_pt = trk.pt; // since I express pt in 0.5 steps this is easy, the address is the pt

            ap_uint<10> this_phi_low_bound    = mluts[itrk].phi_low_bounds    [addr_th][addr_pt];
            ap_uint<10> this_phi_high_bound   = mluts[itrk].phi_high_bounds   [addr_th][addr_pt];
            ap_uint<10> this_theta_low_bound  = mluts[itrk].theta_low_bounds  [addr_th][addr_pt];
            ap_uint<10> this_theta_high_bound = mluts[itrk].theta_high_bounds [addr_th][addr_pt];

            for (size_t imatcher = 0; imatcher < N_TKMU; ++imatcher)
            {
                #pragma HLS unroll
                match_units[imatcher].process(
                    trk,
                    this_phi_low_bound,
                    this_phi_high_bound,
                    this_theta_low_bound,
                    this_theta_high_bound
                );
            }
        }

        corr_done_buffer = 0;
    }

    // ----------------------------------------
    // -- output and mu processing --
    
    // if stream is empty, flush out TkMu and reinit with new muons
    else // re-init everything
    {
        // #ifndef __SYNTHESIS__
        // std::cout << " ... TRACK STREAM IS EMPTY, MAKING CORRESPONDING OPS ..." << std::endl;
        // #endif

        // FIXME: implement sorting across various matchers here here
        // std::cout << "nel blocco mu" << std::endl;

        // copy the TkMu to the output buffers
        for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
        {
            #pragma HLS unroll
            out_buffer[itkmu] = match_units[itkmu].tkmu_;
            corr_done_buffer = 1;

        }

        // reinitialize everything
        for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
        {
            #pragma HLS unroll
            match_units[itkmu].init(in_muons[itkmu]); /// NOTE!! This implicitly requires N_TKMU = N_MU
        }
    }
}



void correlator_nostream (
    muon_t                in_muons  [N_MU],          // inputs : all the muons ready at the same time
    track_t               in_tracks [N_TRK_SECTORS], // inputs : one track incoming per sector as a stream
    ap_uint<1>            sending_trk,               // sending tracks -> start the correlator
    tkmu_t                out_tkmu  [N_TKMU],        // output TkMus
    ap_uint<1>            &corr_done)        
{
    // ----------------------------------------
    // -- definitions of blocks --

    #pragma HLS data_pack variable=in_muons
    #pragma HLS data_pack variable=in_tracks
    #pragma HLS data_pack variable=out_tkmu

    // define the matching units. They must be persistent outside the function -> static
    static matcher match_units[N_TRK_SECTORS][N_TKMU];
    // independent units -> partition the array
    #pragma HLS ARRAY_PARTITION variable=match_units  complete dim=0

    // one set of matching LUT per incoming track
    matching_LUTs mluts[N_TRK_SECTORS];
    #pragma HLS ARRAY_PARTITION variable=mluts  complete

    // and output buffer decoupled from the matcher - so that the data is persistent
    static tkmu_t     out_buffer [N_TKMU];
    static ap_uint<1> corr_done_buffer;
    #pragma HLS ARRAY_PARTITION variable=out_buffer  complete
    // connect output to the buffers
    for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
    {
        #pragma HLS unroll
        out_tkmu[itkmu] = out_buffer[itkmu];
    }
    corr_done = corr_done_buffer;

    // ----------------------------------------
    // -- track processing --
    // process tracks if the stream is full. FIXME: is it OK to check just the first stream if empty?
    // Expect all streams to be filled in the same way, even with empty candidates
    if (sending_trk)
    {
        // #ifndef __SYNTHESIS__
        // std::cout << " ... TRACK STREAM HAS DATA, PROCESSING ..." << std::endl;
        // #endif

        for (size_t itrk = 0; itrk < N_TRK_SECTORS; ++itrk)
        {
            #pragma HLS unroll
            
            // track_t trk = in_tracks[itrk].read();
            track_t trk = in_tracks[itrk];

            // std::cout << " ....... sector " << itrk << " with track pt = " << trk.pt.to_uint64() << " theta = " << trk.theta.to_uint64() << " phi = " << trk.phi.to_uint64() << std::endl;

            // build the matching LUT address
            ap_uint<4> addr_th;
            ap_uint<9> addr_pt;

            ap_uint<TRK_THETA_W_SIZE> abs_trk_theta = trk.theta;
            // ap_int<TRK_THETA_W_SIZE> abs_trk_theta = get_trk_theta(tracks[itrk]);
            // if (abs_trk_theta < 0)
            //     abs_trk_theta *= -1;

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

            addr_pt = trk.pt; // since I express pt in 0.5 steps this is easy, the address is the pt

            ap_uint<10> this_phi_low_bound    = mluts[itrk].phi_low_bounds    [addr_th][addr_pt];
            ap_uint<10> this_phi_high_bound   = mluts[itrk].phi_high_bounds   [addr_th][addr_pt];
            ap_uint<10> this_theta_low_bound  = mluts[itrk].theta_low_bounds  [addr_th][addr_pt];
            ap_uint<10> this_theta_high_bound = mluts[itrk].theta_high_bounds [addr_th][addr_pt];

            for (size_t imatcher = 0; imatcher < N_TKMU; ++imatcher)
            {
                #pragma HLS unroll
                match_units[itrk][imatcher].process(
                    trk,
                    this_phi_low_bound,
                    this_phi_high_bound,
                    this_theta_low_bound,
                    this_theta_high_bound
                );
            }
        }

        corr_done_buffer = 0;
    }

    // ----------------------------------------
    // -- output and mu processing --
    
    // if stream is empty, flush out TkMu and reinit with new muons
    else // re-init everything
    {
        // #ifndef __SYNTHESIS__
        // std::cout << " ... TRACK STREAM IS EMPTY, MAKING CORRESPONDING OPS ..." << std::endl;
        // #endif

        // FIXME: implement sorting across various matchers here here
        // std::cout << "nel blocco mu" << std::endl;

        // copy the TkMu to the output buffers
        for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
        {
            #pragma HLS unroll
            out_buffer[itkmu] = match_units[0][itkmu].tkmu_; // FIXME: sorting
            corr_done_buffer = 1;

        }

        // reinitialize everything
        for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
        {
            for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
            {
                #pragma HLS unroll
                match_units[isec][itkmu].init(in_muons[itkmu]); /// NOTE!! This implicitly requires N_TKMU = N_MU
            }
        }
    }
}




void passthrough (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info,
    ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
    ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
    ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2)
{
    // unpack my input DF into the various muons
    ap_uint<MU_W_SIZE> muons[N_MU]; // input to corr
    ap_uint<TRK_W_SIZE> tracks[N_TRK]; // input to corr
    ap_uint<TKMU_W_SIZE> tkmus[N_TKMU]; // output to corr

    // output array must be clean
    for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
    {
        #pragma HLS unroll
        tkmus[itkmu] = 0x0; // FIXME: shoudl this be done at the testbench level instead?
    }

    #pragma HLS ARRAY_PARTITION variable=muons  complete
    #pragma HLS ARRAY_PARTITION variable=tracks complete
    #pragma HLS ARRAY_PARTITION variable=tkmus  complete

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

    // copy the pt/eta/phi of the first N_TKMU tracks to the tkmus output
    // NOTE: need to have N_TKMU <= N_TRK, but I am not checking array bounds here

    // use a set of static registers to introduce one clk latency so that the logic gets clocked to have the same interface as the correlator
    static ap_uint<TKMU_W_SIZE> registers_1 [N_TKMU];
    static ap_uint<TKMU_W_SIZE> registers_2 [N_TKMU];

    for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
    {
        #pragma HLS unroll
        // ref_to_tkmu_pt(tkmus[itkmu])    = get_trk_pt(tracks[itkmu]); // FIXME: here I am assuming the same bit width
        // ref_to_tkmu_phi(tkmus[itkmu])   = get_trk_phi(tracks[itkmu]); // FIXME: here I am assuming the same bit width
        // ref_to_tkmu_theta(tkmus[itkmu]) = get_trk_theta(tracks[itkmu]); // FIXME: here I am assuming the same bit width

        ref_to_tkmu_pt(tkmus[itkmu])    = get_trk_pt(registers_2[itkmu]); // FIXME: here I am assuming the same bit width
        ref_to_tkmu_phi(tkmus[itkmu])   = get_trk_phi(registers_2[itkmu]); // FIXME: here I am assuming the same bit width
        ref_to_tkmu_theta(tkmus[itkmu]) = get_trk_theta(registers_2[itkmu]); // FIXME: here I am assuming the same bit width

        registers_2[itkmu] = registers_1[itkmu];
        registers_1[itkmu] = tracks[itkmu];


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


void BRAM_to_corr (ap_uint<MTF7_BRAM_SIZE> in_info, ap_uint<1> in_muons, ap_uint<MTF7_BRAM_SIZE> &out_info, ap_uint<1> &out_valid)
{
    // #pragma HLS dataflow

    // first word: the muon data
    muon_t muons[N_MU];

    // the destination array
    tkmu_t tkmus[N_TKMU];

    // hls::stream<track_t> tracks[N_TRK_SECTORS]  = {
    //     hls::stream<track_t>("SECTOR_1"),
    //     hls::stream<track_t>("SECTOR_2"),
    //     hls::stream<track_t>("SECTOR_3"),
    //     hls::stream<track_t>("SECTOR_4"),
    //     hls::stream<track_t>("SECTOR_5"),
    //     hls::stream<track_t>("SECTOR_6"),
    //     hls::stream<track_t>("SECTOR_7"),
    //     hls::stream<track_t>("SECTOR_8"),
    //     hls::stream<track_t>("SECTOR_9")
    // };
    
    static hls::stream<track_t> tracks[N_TRK_SECTORS];

    #pragma HLS ARRAY_PARTITION variable=muons   complete
    #pragma HLS ARRAY_PARTITION variable=tracks  complete
    #pragma HLS ARRAY_PARTITION variable=tkmus   complete
    // hls::stream<track_t> tra("SECTOR_1");

    // #pragma HLS stream variable=tracks[0]
    // #pragma HLS stream variable=tracks[1]
    // #pragma HLS stream variable=tracks[2]
    // #pragma HLS stream variable=tracks[3]
    // #pragma HLS stream variable=tracks[4]
    // #pragma HLS stream variable=tracks[5]
    // #pragma HLS stream variable=tracks[6]
    // #pragma HLS stream variable=tracks[7]
    // #pragma HLS stream variable=tracks[9]

    // static size_t info_counter = 0; // every 1 + N_TRK_PER_SECTOR it's a new event

    // if (info_counter == 0)
    if (in_muons == 1)
    {    
        // convert the word into a muon
        // ap_uint<MTF7_BRAM_SIZE> all_mu_word = in_bram[0 + in_offset];
        for (size_t imu = 0; imu < N_MU; ++imu)
        {
          #pragma HLS unroll

          ap_uint<MU_W_SIZE> mu_word = in_info.range((imu+1)*MU_W_SIZE-1, imu*MU_W_SIZE);
          muons[imu].pt         = get_mu_pt(mu_word);
          muons[imu].theta      = get_mu_theta(mu_word);
          muons[imu].theta_sign = get_mu_theta_sign(mu_word);
          muons[imu].phi        = get_mu_phi(mu_word);
        }

        // for (size_t imu = 0; imu < N_MU; ++imu)
        // {
        //   std::cout << " ... mu " << imu
        //             << " pt = " << muons[imu].pt.to_uint64()
        //             << " theta = " << muons[imu].theta.to_uint64()
        //             << " theta_sign = " << muons[imu].theta_sign.to_uint64()
        //             << " phi = " << muons[imu].phi.to_uint64() << std::endl;
        // }
    }

    
    else
    {
        // second and following word: these are the tracks
        // each memory address carries one incoming set of track for a clk cycle
        // NOTE: does not necessarily has to be one BRAM address / bx
        // at the testbench level, multiple addresses can be combined if more bits are needed

        // printf("NUMBER OF TRACK: %i\n", N_TRK);

        // first round: inject ths muons using empty track streasm
        // std::cout << ".... injecting muons" << std::endl;
        // correlator_stream(muons, tracks, tkmus);

        // std::cout << " ---------- track TMT number " << itrk << " ------------" << std::endl;

        // for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
        // {
        //     // PATCH: set properties of tracks by hand
        //     // to be read from the input pattern file
        //     track_t this_trk;
        //     this_trk.pt    = 100*itrk + isec;
        //     this_trk.theta = 100*itrk + isec + 1;
        //     this_trk.phi   = 100*itrk + isec + 2;
        //     this_trk.nstubs = 5;
        //     this_trk.chisq  = 30;
        //     tracks[isec].write(this_trk);
        // }

        // ap_uint<MTF7_BRAM_SIZE> all_trk_word = in_bram[1 + itrk + in_offset];

        for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
        {
          #pragma HLS unroll
          // #pragma HLS dataflow

          ap_uint<TRK_W_SIZE> trk_word = in_info.range((isec+1)*TRK_W_SIZE-1, isec*TRK_W_SIZE);
          track_t this_trk;
          this_trk.pt         = get_trk_pt(trk_word);
          this_trk.theta      = get_trk_theta(trk_word);
          this_trk.theta_sign = get_trk_theta_sign(trk_word);
          this_trk.phi        = get_trk_phi(trk_word);
          this_trk.charge     = get_trk_charge(trk_word);
          this_trk.chisq      = get_trk_chisq(trk_word);
          this_trk.nstubs     = get_trk_nstubs(trk_word);

          // std::cout << " ... tracks - sector = " << isec << " tmt = " << itrk
          //           << " pt = "         << this_trk.pt.to_uint64()
          //           << " theta = "      << this_trk.theta.to_uint64()
          //           << " theta_sign = " << this_trk.theta_sign.to_uint64()
          //           << " phi = "        << this_trk.phi.to_uint64()
          //           << " charge = "     << this_trk.charge.to_uint64()
          //           << " chisq = "      << this_trk.chisq.to_uint64()
          //           << " nstubs = "     << this_trk.nstubs.to_uint64()
          //           << std::endl;

          tracks[isec].write(this_trk);
        }

        // std::cout << ".... running correlator on a new set of tracks" << std::endl;
    }

    ap_uint<1> corr_done;
    out_valid = corr_done; // connect corr_done to out_valid
    correlator_stream(muons, tracks, tkmus, corr_done);
    // info_counter == 1 + N_TRK_PER_SEC -1 ? info_counter = 0 : ++info_counter;

    // std::cout << "CORR DONE " << corr_done.to_uint() << std::endl;
    // for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
    // {
    //   std::cout << " ... tkmu " << itkmu
    //             << " pt = " << tkmus[itkmu].pt.to_uint64()
    //             << " theta = " << tkmus[itkmu].theta.to_uint64()
    //             << " theta_sign = " << tkmus[itkmu].theta_sign.to_uint64()
    //             << " phi = " << tkmus[itkmu].phi.to_uint64() << std::endl;
    // }


    /// repack the output in out_bram and write it to file in output
    // ap_uint<MTF7_BRAM_SIZE> out_info;
    if (corr_done)
    // if (true)
    {
        // for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
        // {
        //   std::cout << " ... tkmu " << itkmu
        //             << " pt = " << tkmus[itkmu].pt.to_uint64()
        //             << " theta = " << tkmus[itkmu].theta.to_uint64()
        //             << " theta_sign = " << tkmus[itkmu].theta_sign.to_uint64()
        //             << " phi = " << tkmus[itkmu].phi.to_uint64() << std::endl;
        // }

        for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
        {
          #pragma HLS unroll
          out_info.range(TKMU_W_SIZE*(itkmu+1)-1, TKMU_W_SIZE*itkmu) = build_tkmu_word(tkmus[itkmu]);
        }
    }
    
    // out_bram[itest] = out_info;

}





void BRAM_to_corr_nostream (ap_uint<MTF7_BRAM_SIZE> in_info, ap_uint<1> in_muons, ap_uint<MTF7_BRAM_SIZE> &out_info, ap_uint<1> &out_valid)
{
    // #pragma HLS dataflow

    // printf("e poi ... %0x\n",in_info.range(31,0).to_uint64() );

    // std::cout << " ........ inside corr ||| "
    //         << " in_info_last : " << std::hex << in_info.range(63,0).to_uint64()
    //         << " is_muons : " << in_muons.to_uint()
    //         << " out_info_last : " << std::hex << out_info.range(63,0).to_uint64()
    //         << " out_valid : " << out_valid 
    //         << std::endl;


    // first word: the muon data
    muon_t muons[N_MU];

    // the destination array
    tkmu_t tkmus[N_TKMU];

    // hls::stream<track_t> tracks[N_TRK_SECTORS]  = {
    //     hls::stream<track_t>("SECTOR_1"),
    //     hls::stream<track_t>("SECTOR_2"),
    //     hls::stream<track_t>("SECTOR_3"),
    //     hls::stream<track_t>("SECTOR_4"),
    //     hls::stream<track_t>("SECTOR_5"),
    //     hls::stream<track_t>("SECTOR_6"),
    //     hls::stream<track_t>("SECTOR_7"),
    //     hls::stream<track_t>("SECTOR_8"),
    //     hls::stream<track_t>("SECTOR_9")
    // };
    
    track_t tracks[N_TRK_SECTORS];

    #pragma HLS ARRAY_PARTITION variable=muons   complete
    #pragma HLS ARRAY_PARTITION variable=tracks  complete
    #pragma HLS ARRAY_PARTITION variable=tkmus   complete
    // hls::stream<track_t> tra("SECTOR_1");

    // #pragma HLS stream variable=tracks[0]
    // #pragma HLS stream variable=tracks[1]
    // #pragma HLS stream variable=tracks[2]
    // #pragma HLS stream variable=tracks[3]
    // #pragma HLS stream variable=tracks[4]
    // #pragma HLS stream variable=tracks[5]
    // #pragma HLS stream variable=tracks[6]
    // #pragma HLS stream variable=tracks[7]
    // #pragma HLS stream variable=tracks[9]

    // static size_t info_counter = 0; // every 1 + N_TRK_PER_SECTOR it's a new event

    // if (info_counter == 0)
    if (in_muons == 1)
    {    
        // convert the word into a muon
        // ap_uint<MTF7_BRAM_SIZE> all_mu_word = in_bram[0 + in_offset];
        for (size_t imu = 0; imu < N_MU; ++imu)
        {
          #pragma HLS unroll

          ap_uint<MU_W_SIZE> mu_word = in_info.range((imu+1)*MU_W_SIZE-1, imu*MU_W_SIZE);
          muons[imu].pt         = get_mu_pt(mu_word);
          muons[imu].theta      = get_mu_theta(mu_word);
          muons[imu].theta_sign = get_mu_theta_sign(mu_word);
          muons[imu].phi        = get_mu_phi(mu_word);
        }

        // for (size_t imu = 0; imu < N_MU; ++imu)
        // {
        //   std::cout << " ... mu " << imu
        //             << " pt = " << muons[imu].pt.to_uint64()
        //             << " theta = " << muons[imu].theta.to_uint64()
        //             << " theta_sign = " << muons[imu].theta_sign.to_uint64()
        //             << " phi = " << muons[imu].phi.to_uint64() << std::endl;
        // }
    }

    
    else
    {
        // second and following word: these are the tracks
        // each memory address carries one incoming set of track for a clk cycle
        // NOTE: does not necessarily has to be one BRAM address / bx
        // at the testbench level, multiple addresses can be combined if more bits are needed

        // printf("NUMBER OF TRACK: %i\n", N_TRK);

        // first round: inject ths muons using empty track streasm
        // std::cout << ".... injecting muons" << std::endl;
        // correlator_stream(muons, tracks, tkmus);

        // std::cout << " ---------- track TMT number " << itrk << " ------------" << std::endl;

        // for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
        // {
        //     // PATCH: set properties of tracks by hand
        //     // to be read from the input pattern file
        //     track_t this_trk;
        //     this_trk.pt    = 100*itrk + isec;
        //     this_trk.theta = 100*itrk + isec + 1;
        //     this_trk.phi   = 100*itrk + isec + 2;
        //     this_trk.nstubs = 5;
        //     this_trk.chisq  = 30;
        //     tracks[isec].write(this_trk);
        // }

        // ap_uint<MTF7_BRAM_SIZE> all_trk_word = in_bram[1 + itrk + in_offset];

        for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
        {
          #pragma HLS unroll
          // #pragma HLS dataflow

          ap_uint<TRK_W_SIZE> trk_word = in_info.range((isec+1)*TRK_W_SIZE-1, isec*TRK_W_SIZE);
          track_t this_trk;
          this_trk.pt         = get_trk_pt(trk_word);
          this_trk.theta      = get_trk_theta(trk_word);
          this_trk.theta_sign = get_trk_theta_sign(trk_word);
          this_trk.phi        = get_trk_phi(trk_word);
          this_trk.charge     = get_trk_charge(trk_word);
          this_trk.chisq      = get_trk_chisq(trk_word);
          this_trk.nstubs     = get_trk_nstubs(trk_word);

          // std::cout << " ... tracks - sector = " << isec << " tmt = " << itrk
          //           << " pt = "         << this_trk.pt.to_uint64()
          //           << " theta = "      << this_trk.theta.to_uint64()
          //           << " theta_sign = " << this_trk.theta_sign.to_uint64()
          //           << " phi = "        << this_trk.phi.to_uint64()
          //           << " charge = "     << this_trk.charge.to_uint64()
          //           << " chisq = "      << this_trk.chisq.to_uint64()
          //           << " nstubs = "     << this_trk.nstubs.to_uint64()
          //           << std::endl;

          tracks[isec] = this_trk;
        }

        // std::cout << ".... running correlator on a new set of tracks" << std::endl;
    }

    // ap_uint<1> corr_done;
    // out_valid = corr_done; // connect corr_done to out_valid
    correlator_nostream(muons, tracks, !in_muons, tkmus, out_valid);
    // info_counter == 1 + N_TRK_PER_SEC -1 ? info_counter = 0 : ++info_counter;

    // std::cout << "CORR DONE " << corr_done.to_uint() << std::endl;
    // for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
    // {
    //   std::cout << " ... tkmu " << itkmu
    //             << " pt = " << tkmus[itkmu].pt.to_uint64()
    //             << " theta = " << tkmus[itkmu].theta.to_uint64()
    //             << " theta_sign = " << tkmus[itkmu].theta_sign.to_uint64()
    //             << " phi = " << tkmus[itkmu].phi.to_uint64() << std::endl;
    // }


    /// repack the output in out_bram and write it to file in output
    // ap_uint<MTF7_BRAM_SIZE> out_info;
    if (out_valid)
    // if (true)
    {
        // for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
        // {
        //   std::cout << " ... tkmu " << itkmu
        //             << " pt = " << tkmus[itkmu].pt.to_uint64()
        //             << " theta = " << tkmus[itkmu].theta.to_uint64()
        //             << " theta_sign = " << tkmus[itkmu].theta_sign.to_uint64()
        //             << " phi = " << tkmus[itkmu].phi.to_uint64() << std::endl;
        // }

        for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
        {
          #pragma HLS unroll
          out_info.range(TKMU_W_SIZE*(itkmu+1)-1, TKMU_W_SIZE*itkmu) = build_tkmu_word(tkmus[itkmu]);
        }
    }
    
    // out_bram[itest] = out_info;

}