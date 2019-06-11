#include <cstdio>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include "src/correlator.h"
#include "ap_int.h"
#include "ap_fixed.h"

#define NTEST 1000

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

// reads a line of 64 chunks of 32 bits from a file and converts them to the write buffer
// each chunk must be in hex format
// Assuming that the word (1 line in the txt file) is written as: [chunk63] [chunk62] [chunk61] ... [chunk1] [chunk0], then
//
// reverse = false: ==> the word is written to the buffer in the same way you can "read" it from the txt file
// buffer [0] = [chunk0]
// buffer [1] = [chunk1]
// ...
// buffer [63] = [chunk63]
//
// reverse = true: ==> the word chunks order is inverted w.r.t the original one in the txt file
// buffer [0] = [chunk63]
// buffer [1] = [chunk62]
// ...
// buffer [63] = [chunk0]
//
// NB!! the input txt file is read from left to right, so the incoming arrival order is [chunk63] > [chunk62] > ... > [chunk1] [chunk0]
void inject_load_from_file (
        ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> *in_bram, unsigned int nbx,
        std::string input_file_name, bool reverse=false)
{

  // first clean up all elements of in_bram
  for (size_t i = 0; i < nbx; ++i)
    in_bram[i] = 0;

  std::fstream input_file(input_file_name.c_str());

  // each line contains 64 buckets of 32 bits, encoded in hex format
  std::string line;
  int nlines = 0;
  while (std::getline(input_file, line))
  {
      std::istringstream iss(line);
      int nwords = 0;
      uint32_t wtemp;
      uint32_t words[64];
      ap_uint<2048> tmp_word; // patterns are for a 2048 word - fill one, then trim down to the real workd size

      while (iss >> std::hex >> wtemp)
      {
        // buffer size is limited to 64 -> if this is already exceeded, throw error message
        if (nwords >= 64){
          std::ostringstream strnlines;
          strnlines << nlines;
          throw std::runtime_error(std::string("Error in parsing the patter input file at line: ") + strnlines.str());
        }

        if (reverse)
          words[nwords] = wtemp;
        else
          words[63-nwords] = wtemp;
        
        ++nwords;
      }
      std::cout << "........ wpatt_load_from_file : from line " << nlines << " read " << nwords << " words" << std::endl;
      
      // now copy the read line to the ap input buffer
      // unsigned int offset = nlines*64;
      for (size_t i = 0; i < 64; ++i){
        tmp_word.range((i+1)*32-1, 32*i) = words[i];
        // in_bram[nlines]
        // write_buf[i + offset] = words[i];
      }

      // now trim the input word to the actual size
      in_bram[nlines] = tmp_word.range(TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU -1, 0);

      // and set the start and reset bits for this word
      //set_start_reset_bits_tobuffer (write_buf, nlines);

      // update the line count
      ++nlines;

      // this was the last element in the in_bram array
      if (nlines >= nbx)
        break;
  }

  std::cout << "... wpatt_load_from_file : read " << nlines << " input patterns" << std::endl;
}

void dump_output_to_file (ap_uint<TKMU_W_SIZE*N_TKMU> *out_bram, unsigned int nbx, std::string output_file_name)
{
    FILE * pFile = fopen (output_file_name.c_str(),"w");
    fprintf (pFile, "iEv itkmu pt theta thetasign phi charge\n");
    for (size_t ibx = 0; ibx < nbx; ++ibx)
    {
        for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
        {
            ap_uint<TKMU_W_SIZE> tkmu = out_bram[ibx].range(TKMU_W_SIZE*(itkmu+1)-1, TKMU_W_SIZE*itkmu);
            fprintf (pFile, "%5i %2i %10u %10u %1u %10u %1u\n",
                ibx, itkmu,
                get_tkmu_pt(tkmu).to_uint64(),
                get_tkmu_theta(tkmu).to_uint64(),
                get_tkmu_theta_sign(tkmu).to_uint64(),
                get_tkmu_phi(tkmu).to_uint64(),
                get_tkmu_charge(tkmu).to_uint64());
        }
    }
    fclose (pFile);
}



/////// ancillary functions to print to screen an ap::uint value
template <size_t N>
std::string ap_uint_to_string (ap_uint<N> val, bool pad_space = true)
{
    size_t left = N;
    size_t offset = 0;
    std::string result = "";
    std::string sep = pad_space ? " " : "";
    // the ones I cannot fit in uint64
    while (left > 64)
    {

        std::string out_string;
        std::stringstream ss;
        ss << std::hex << std::setfill('0') << std::setw(64/4) << val.range(63+offset, 0+offset).to_uint64();
        out_string = ss.str();
        result = out_string + sep + result;
        // std::cout << " ..... " << 63+offset << " " << 0+offset << " ,,, " << out_string << std::endl;
        left -= 64;
        offset += 64;
    }

    // the final chiunk
    std::string out_string;
    std::stringstream ss;
    size_t tot_s = N - offset;
    ss << std::hex << std::setfill('0') << std::setw(tot_s/4) << val.range(N-1, 0+offset).to_uint64();
    out_string = ss.str();
    result = out_string + sep + result;

    return result;
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

    // PU 0 I/O
    // std::string in_pattern_filename = "/home/lcadamur/CorrPCIeGit/MuTrackCorrMTF7Utils/patterns/test_patterns_PU0.txt";
    // std::string out_tkmu_filename   = "/home/lcadamur/MuTrackCorrHLS/MuTrackCorrHLS/corr_output_PU0.txt";

    // PU 200 I/O
    std::string in_pattern_filename = "/home/lcadamur/CorrPCIeGit/MuTrackCorrMTF7Utils/patterns/test_patterns.txt";
    std::string out_tkmu_filename   = "/home/lcadamur/MuTrackCorrHLS/MuTrackCorrHLS/corr_output.txt";


    // inject_debugdata_correlator(in_bram, nbx);
    // inject_simple_idx(in_bram, nbx);
    // inject_null(in_bram, nbx);
    inject_load_from_file (in_bram, nbx, in_pattern_filename.c_str());

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
        if (itest < 10)
        {
            printf(">>> checking BX %2i/%3i :: mu1 = %13x, mu2 = %13x, trk1 = %13x, trk2 = %13x, tkmu1 = %13x, tkmu2 = %13x, (outw = %llx)\n",
            itest, NTEST,
            spy_mu1.to_int(),   spy_mu2.to_int(),
            spy_trk1.to_int(),  spy_trk2.to_int(),
            spy_tkmu1.to_int(), spy_tkmu2.to_int(),
            out_bram[itest].to_uint64() // will trim down to the 64 lsb, but OK for debugging
            );
            printf("     .... mu1   .... pt = %10u, theta = (%1u) %10u, phi = %10u, charge = %1u, word = %s\n",
                get_mu_pt(spy_mu1).to_uint64(),
                get_mu_theta_sign(spy_mu1).to_uint(), 
                get_mu_theta(spy_mu1).to_uint64(),
                get_mu_phi(spy_mu1).to_uint64(),
                get_mu_charge(spy_mu1).to_uint64(),
                ap_uint_to_string<MU_W_SIZE>(spy_mu1).c_str());
            printf("     .... mu2   .... pt = %10u, theta = (%1u) %10u, phi = %10u, charge = %1u, word = %s\n",
                get_mu_pt(spy_mu2).to_uint64(),
                get_mu_theta_sign(spy_mu2).to_uint(),
                get_mu_theta(spy_mu2).to_uint64(),
                get_mu_phi(spy_mu2).to_uint64(),
                get_mu_charge(spy_mu2).to_uint64(),
                ap_uint_to_string<MU_W_SIZE>(spy_mu2).c_str());
            printf("     .... trk1  .... pt = %10u, theta = (%1u) %10u, phi = %10u, charge = %1u, chisq = %10u, nstubs = %10u, word = %s\n",
                get_trk_pt(spy_trk1).to_uint64(),
                get_trk_theta_sign(spy_trk1).to_uint(),
                get_trk_theta(spy_trk1).to_uint64(),
                get_trk_phi(spy_trk1).to_uint64(),
                get_trk_charge(spy_trk1).to_uint64(),
                get_trk_chisq(spy_trk1).to_uint64(),
                get_trk_nstubs(spy_trk1).to_uint64(),
                ap_uint_to_string<TRK_W_SIZE>(spy_trk1).c_str());
            printf("     .... trk2  .... pt = %10u, theta = (%1u) %10u, phi = %10u, charge = %1u, chisq = %10u, nstubs = %10u, word = %s\n",
                get_trk_pt(spy_trk2).to_uint64(),
                get_trk_theta_sign(spy_trk2).to_uint(),
                get_trk_theta(spy_trk2).to_uint64(),
                get_trk_phi(spy_trk2).to_uint64(),
                get_trk_charge(spy_trk2).to_uint64(),
                get_trk_chisq(spy_trk2).to_uint64(),
                get_trk_nstubs(spy_trk2).to_uint64(),
                ap_uint_to_string<TRK_W_SIZE>(spy_trk2).c_str());
            printf("     .... tkmu1 .... pt = %10u, theta = %10u, phi = %10u, word = %s\n",
                get_tkmu_pt(spy_tkmu1).to_uint64(),
                get_tkmu_theta(spy_tkmu1).to_uint64(),
                get_tkmu_phi(spy_tkmu1).to_uint64(),
                ap_uint_to_string<TKMU_W_SIZE>(spy_tkmu1).c_str());
            printf("     .... tkmu2 .... pt = %10u, theta = %10u, phi = %10u, word = %s\n",
                get_tkmu_pt(spy_tkmu2).to_uint64(),
                get_tkmu_theta(spy_tkmu2).to_uint64(),
                get_tkmu_phi(spy_tkmu2).to_uint64(),
                ap_uint_to_string<TKMU_W_SIZE>(spy_tkmu2).c_str());
        }
        else if (itest % 100 == 0)
            printf(">>> checking BX %2i/%3i :: mu1 = %13x, mu2 = %13x, trk1 = %13x, trk2 = %13x, tkmu1 = %13x, tkmu2 = %13x, (outw = %llx)\n",
            itest, NTEST,
            spy_mu1.to_int(),   spy_mu2.to_int(),
            spy_trk1.to_int(),  spy_trk2.to_int(),
            spy_tkmu1.to_int(), spy_tkmu2.to_int(),
            out_bram[itest].to_uint64() // will trim down to the 64 lsb, but OK for debugging
            );        // printf(">>> checking BX %i/%i :: IN = %16llx,      OUT = %16llx\n",
        // itest, NTEST,
        // in_bram[itest].to_uint64(), out_bram[itest].to_uint64()
        // );
    }

    printf("... dumping result to file\n");
    // note: normally file woudl go under "projCorrelator/solution1/csim/build/"
    // so give here an abs path
    dump_output_to_file (out_bram, NTEST, out_tkmu_filename.c_str());
    printf("... done\n");
}