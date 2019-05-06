#include <cstdio>
#include "src/correlator.h"
#include "ap_int.h"
#include "ap_fixed.h"

#define NTEST 20

void inject_simple_idx (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> *in_bram, unsigned int nbx)
{
    // create the patternr from A.M. : chunks of 32 bits that contain 0,1,2,3,4...
    for (unsigned int ibx = 0; ibx < nbx; ++ibx)
    {
        ap_uint<2048> bramword; // because in the cpp PCIe code I have an I/O of 2048
        int offset = ibx*64;
        for (uint ipart = 0; ipart < 64; ++ipart)
        {
            bramword.range((ipart+1)*32-1, ipart*32) = ipart+offset;
        }
        
        // now trim down to the  in_bram size
        in_bram[ibx] = bramword.range(TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU-1, 0);
    }
}

void inject_null (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> *in_bram, unsigned int nbx)
{
    for (unsigned int ibx = 0; ibx < nbx; ++ibx)
        in_bram[ibx] = 0; 
}

void inject_debugdata_correlator (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> *in_bram, unsigned int nbx)
{
    for (unsigned int ibx = 0; ibx < nbx; ++ibx)
    {
        ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> input_word;

        // define muons
        for (uint imu = 0; imu < N_MU; ++imu)
        {
            ap_uint<MU_W_SIZE> my_mu = 0xf00*ibx + imu;
            input_word.range((imu+1)*MU_W_SIZE-1, imu*MU_W_SIZE) = my_mu;
            ap_uint<MU_W_SIZE> xcheck = input_word.range((imu+1)*MU_W_SIZE-1, imu*MU_W_SIZE);
            // printf(">> Muon idx %i --> %i (vs %i)\n", imu, my_mu.to_int(), xcheck.to_int());
        }

        const int goffset = MU_W_SIZE*N_MU;
        // define tracks
        for (uint itrk = 0; itrk < N_TRK; ++itrk)
        {
            ap_uint<TRK_W_SIZE> my_trk = 0xf00*ibx + itrk;
            input_word.range((itrk+1)*TRK_W_SIZE-1 + goffset, itrk*TRK_W_SIZE + goffset) = my_trk;
            ap_uint<TRK_W_SIZE> xcheck = input_word.range((itrk+1)*TRK_W_SIZE-1 + goffset, itrk*TRK_W_SIZE + goffset);
            // printf(">> Track idx %i --> %i (vs %i)\n", itrk, my_trk.to_int(), xcheck.to_int());
        }

        in_bram[ibx] = input_word;
    }
}


int main()
{
    const unsigned int nbx = NTEST;

    // my input BRAM
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_bram[nbx];
    // my output BRAM
    ap_uint<TKMU_W_SIZE*N_TKMU>                out_bram[nbx];


    inject_debugdata_correlator(in_bram, nbx);
    // inject_simple_idx(in_bram, nbx);
    // inject_null(in_bram, nbx);
    


    // output BRAM values should all be 0 instead
    for (unsigned int ibx = 0; ibx < nbx; ++ibx)
        out_bram[ibx] = 0;


    // let's debug the content:
    ap_uint<TRK_W_SIZE> check_trk1;
    ap_uint<TRK_W_SIZE> check_trk2;
    ap_uint<MU_W_SIZE>  check_mu1 ;
    ap_uint<MU_W_SIZE>  check_mu2 ;
    int ibx;
    const int goffset = MU_W_SIZE*N_MU;

    // now test the algo
    ap_uint<TRK_W_SIZE>  spy_trk1;
    ap_uint<TRK_W_SIZE>  spy_trk2;

    ap_uint<MU_W_SIZE>   spy_mu1;
    ap_uint<MU_W_SIZE>   spy_mu2;

    ap_uint<TKMU_W_SIZE> spy_tkmu1;
    ap_uint<TKMU_W_SIZE> spy_tkmu2;

    for (uint itest = 0; itest < NTEST; ++itest)
    {
        passthrough(in_bram[itest], out_bram[itest],
            spy_trk1,  spy_trk2,
            spy_mu1,   spy_mu2,
            spy_tkmu1, spy_tkmu2);
        printf(">>> checking BX %2i/%3i :: mu1 = %13x, mu2 = %13x, trk1 = %13x, trk2 = %13x, tkmu1 = %13x, tkmu2 = %13x, (outw = %llx)\n",
        itest, NTEST,
        spy_mu1.to_int(),   spy_mu2.to_int(),
        spy_trk1.to_int(),  spy_trk2.to_int(),
        spy_tkmu1.to_int(), spy_tkmu2.to_int(),
        out_bram[itest].to_uint64() // will trim down to the 64 lsb, but OK for debugging
        );
    }

}