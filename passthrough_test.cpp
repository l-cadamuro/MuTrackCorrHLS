#include <cstdio>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
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


int main()
{
    const unsigned int nbx = NTEST;

    // my input BRAM
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_bram[nbx];
    // my output BRAM
    ap_uint<TKMU_W_SIZE*N_TKMU>                out_bram[nbx];


    // PU 0 I/O
    // std::string in_pattern_filename = "/home/lcadamur/CorrPCIeGit/MuTrackCorrMTF7Utils/patterns/test_patterns_PU0.txt";
    // std::string out_tkmu_filename   = "/home/lcadamur/MuTrackCorrHLS/MuTrackCorrHLS/corr_output_PU0.txt";

    // PU 200 I/O
    std::string in_pattern_filename = "/home/lcadamur/CorrPCIeGit/MuTrackCorrMTF7Utils/patterns/test_patterns.txt";
    std::string out_tkmu_filename   = "/home/lcadamur/MuTrackCorrHLS/MuTrackCorrHLS/passthrough_output.txt";


    // inject_debugdata_correlator(in_bram, nbx);
    // inject_simple_idx(in_bram, nbx);
    // inject_null(in_bram, nbx);
    inject_load_from_file (in_bram, nbx, in_pattern_filename.c_str());
    

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

    printf("... dumping result to file\n");
    // note: normally file woudl go under "projCorrelator/solution1/csim/build/"
    // so give here an abs path
    dump_output_to_file (out_bram, NTEST, out_tkmu_filename.c_str());
    printf("... done\n");

}