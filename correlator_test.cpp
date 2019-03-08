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

    // define the values in the input bram: for simplicity, assign each trk a number from 1 to ntrk
    // offset from event to event by f00 (3840)
    
    /*
    for (unsigned int ibx = 0; ibx < nbx; ++ibx)
    {
        ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> input_word;

        // define muons
        for (uint imu = 0; imu < N_MU; ++imu)
        {
            ap_uint<MU_W_SIZE> my_mu = 0xf00*ibx + imu;
            input_word.range((imu+1)*MU_W_SIZE-1, imu*MU_W_SIZE) = my_mu;
            ap_uint<MU_W_SIZE> xcheck = input_word.range((imu+1)*MU_W_SIZE-1, imu*MU_W_SIZE);
            printf(">> Muon idx %i --> %i (vs %i)\n", imu, my_mu.to_int(), xcheck.to_int());
        }

        const int goffset = MU_W_SIZE*N_MU;
        // define tracks
        for (uint itrk = 0; itrk < N_TRK; ++itrk)
        {
            ap_uint<TRK_W_SIZE> my_trk = 0xf00*ibx + itrk;
            input_word.range((itrk+1)*TRK_W_SIZE-1 + goffset, itrk*TRK_W_SIZE + goffset) = my_trk;
            ap_uint<TRK_W_SIZE> xcheck = input_word.range((itrk+1)*TRK_W_SIZE-1 + goffset, itrk*TRK_W_SIZE + goffset);
            printf(">> Track idx %i --> %i (vs %i)\n", itrk, my_trk.to_int(), xcheck.to_int());
        }

        in_bram[ibx] = input_word;
    }
    */

    inject_debugdata_correlator(in_bram, nbx);
    // inject_simple_idx(in_bram, nbx);
    // inject_null(in_bram, nbx);
    

    /*
    // simple one => just inject the number to compare with the c++ code
    for (unsigned int ibx = 0; ibx < nbx; ++ibx)
    {
        ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> input_word;
        // split the input word in 64 chunks of 32 bits and assign to each one its content - from 0 to 64
        // NOTE: the correlator actually takes a bucket of TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU = 1980 bits, which make 61 full words
        // so just fill for 61 of them
        for (uint ipart = 0; ipart < 61; ++ipart)
            // input_word.range((ipart+1)*32-1, ipart*32) = ipart; // 0, 1, 2, ..
            input_word.range((ipart+1)*32-1, ipart*32) = 0xff - (ipart%32); // 255, 254, 253 ...

        in_bram[ibx] = input_word;
        // NOTE: to_int will trim to the lower bits
        // printf("All'indice %i : uint = %i\n", ibx, input_word.range(63,32).to_uint());
        // printf("All'indice %i : uint64 = %04x\n", ibx, input_word.range(63,0).to_uint64());
        // printf("All'indice %i : uint64 = %02x%02x\n", ibx, input_word.range(63,32).to_int(), input_word.range(31,0).to_int() ); // this works, but is ugly
        // if (input_word.range(63,0).to_uint64() > 70000)
        //     printf("Ha! You are indeed larger than 255\n");

        // // make my print by hand, just to check if that works as expected
        // for (int ip = 63; ip >= 0; --ip)
        // {
        //     // printf("ip = %i\n", ip);
        //     bool pass = (input_word & (1 << ip));
        //     printf("%s", (pass ? "1" : "0") );
        // }
        // printf("DONE\n");

        printf("All'indice %i : uint64 = %16llx\n", ibx, input_word.range(63,0).to_uint64());
    }
    */
    

    /*
    // just check I did not overwrite anything by mistake
    for (unsigned int ibx = 0; ibx < nbx; ++ibx)
    {
        printf("\n--------------- This is a cross check for bx = %i ---------------- \n", ibx);

        ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> input_word = in_bram[ibx];

        for (uint imu = 0; imu < N_MU; ++imu)
        {
            ap_uint<MU_W_SIZE> xcheck = input_word.range((imu+1)*MU_W_SIZE-1, imu*MU_W_SIZE);
            printf(">> XC: Muon idx %i --> %i\n", imu, xcheck.to_int());
        }

        const int goffset = MU_W_SIZE*N_MU;
        // define tracks
        for (uint itrk = 0; itrk < N_TRK; ++itrk)
        {
            ap_uint<TRK_W_SIZE> xcheck = input_word.range((itrk+1)*TRK_W_SIZE-1 + goffset, itrk*TRK_W_SIZE + goffset);
            printf(">> XC: Track idx %i --> %i\n", itrk, xcheck.to_int());
        }
    }
    */

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

    /*
    // test bx 0
    ibx = 0;
    check_mu1   = in_bram[ibx].range(MU_W_SIZE-1, 0);
    check_mu2   = in_bram[ibx].range(MU_W_SIZE*2-1, MU_W_SIZE);
    check_trk1  = in_bram[ibx].range(TRK_W_SIZE-1 + MU_W_SIZE*N_MU, MU_W_SIZE*N_MU);
    check_trk2  = in_bram[ibx].range(TRK_W_SIZE*2-1 + MU_W_SIZE*N_MU, TRK_W_SIZE + MU_W_SIZE*N_MU);
    // check_mu1  = in_bram[ibx].range((0+1)*MU_W_SIZE-1, 0*MU_W_SIZE);
    // check_mu2  = in_bram[ibx].range((1+1)*MU_W_SIZE-1, 1*MU_W_SIZE);
    // check_trk1 = in_bram[ibx].range((0+1)*TRK_W_SIZE-1 + goffset, 0*TRK_W_SIZE + goffset);
    // check_trk2 = in_bram[ibx].range((1+1)*TRK_W_SIZE-1 + goffset, 1*TRK_W_SIZE + goffset);
    printf("BX %2i/%3i :: mu1 = %i, mu2 = %i, trk1 = %i, trk2 = %i\n",
        ibx, nbx,
        check_mu1.to_int(), check_mu2.to_int(), check_trk1.to_int(), check_trk2.to_int());

    // test bx 1
    ibx = 1;
    check_mu1   = in_bram[ibx].range(MU_W_SIZE-1, 0);
    check_mu2   = in_bram[ibx].range(MU_W_SIZE*2-1, MU_W_SIZE);
    check_trk1  = in_bram[ibx].range(TRK_W_SIZE-1 + MU_W_SIZE*N_MU, MU_W_SIZE*N_MU);
    check_trk2  = in_bram[ibx].range(TRK_W_SIZE*2-1 + MU_W_SIZE*N_MU, TRK_W_SIZE + MU_W_SIZE*N_MU);
    // check_mu1  = in_bram[ibx].range((0+1)*MU_W_SIZE-1, 0*MU_W_SIZE);
    // check_mu2  = in_bram[ibx].range((1+1)*MU_W_SIZE-1, 1*MU_W_SIZE);
    // check_trk1 = in_bram[ibx].range((0+1)*TRK_W_SIZE-1 + goffset, 0*TRK_W_SIZE + goffset);
    // check_trk2 = in_bram[ibx].range((1+1)*TRK_W_SIZE-1 + goffset, 1*TRK_W_SIZE + goffset);
    printf("BX %2i/%3i :: mu1 = %i, mu2 = %i, trk1 = %i, trk2 = %i\n",
        ibx, nbx,
        check_mu1.to_int(), check_mu2.to_int(), check_trk1.to_int(), check_trk2.to_int());

    // test bx 90
    ibx = 90;
    check_mu1   = in_bram[ibx].range(MU_W_SIZE-1, 0);
    check_mu2   = in_bram[ibx].range(MU_W_SIZE*2-1, MU_W_SIZE);
    check_trk1  = in_bram[ibx].range(TRK_W_SIZE-1 + MU_W_SIZE*N_MU, MU_W_SIZE*N_MU);
    check_trk2  = in_bram[ibx].range(TRK_W_SIZE*2-1 + MU_W_SIZE*N_MU, TRK_W_SIZE + MU_W_SIZE*N_MU);
    // check_mu1  = in_bram[ibx].range((0+1)*MU_W_SIZE-1, 0*MU_W_SIZE);
    // check_mu2  = in_bram[ibx].range((1+1)*MU_W_SIZE-1, 1*MU_W_SIZE);
    // check_trk1 = in_bram[ibx].range((0+1)*TRK_W_SIZE-1 + goffset, 0*TRK_W_SIZE + goffset);
    // check_trk2 = in_bram[ibx].range((1+1)*TRK_W_SIZE-1 + goffset, 1*TRK_W_SIZE + goffset);
    printf("BX %2i/%3i :: mu1 = %i, mu2 = %i, trk1 = %i, trk2 = %i\n",
        ibx, nbx,
        check_mu1.to_int(), check_mu2.to_int(), check_trk1.to_int(), check_trk2.to_int());
    */

    // now test the algo
    ap_uint<TRK_W_SIZE>  spy_trk1;
    ap_uint<TRK_W_SIZE>  spy_trk2;

    ap_uint<MU_W_SIZE>   spy_mu1;
    ap_uint<MU_W_SIZE>   spy_mu2;

    ap_uint<TKMU_W_SIZE> spy_tkmu1;
    ap_uint<TKMU_W_SIZE> spy_tkmu2;

    for (uint itest = 0; itest < NTEST; ++itest)
    {
        correlator_one(in_bram[itest], out_bram[itest],
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
        // printf(">>> checking BX %i/%i :: IN = %16llx,      OUT = %16llx\n",
        // itest, NTEST,
        // in_bram[itest].to_uint64(), out_bram[itest].to_uint64()
        // );

    }

}