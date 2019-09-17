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
// #define MTF7_BRAM_SIZE 2048

namespace patt_inj
{
    template <size_t MEMSIZE>
    void inject_dummy (ap_uint<MEMSIZE> *in_bram, unsigned int nbx)
    {
        for (unsigned int ibx = 0; ibx < nbx; ++ibx)
            in_bram[nbx] = ibx;
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

    template <size_t MEMSIZE>
    void inject_load_from_file (
            ap_uint<MEMSIZE> *in_bram, unsigned int nbx,
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

            // if (nlines == 0)
            //   printf("bram %i --> %x\n", (reverse ? nwords : 63-nwords) , wtemp);
            
            ++nwords;
          }
          // std::cout << "........ wpatt_load_from_file : from line " << nlines << " read " << nwords << " words" << std::endl;
          
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
}

template <size_t MEMSIZE> 
void dump_output_to_file (ap_uint<MEMSIZE> *out_bram, unsigned int nbx, std::string output_file_name)
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
    std::cout << "... starting testbench" << std::endl;

    // the I/O from the MTF7 BRAM
    const size_t lines_to_read = NTEST*(1+N_TRK_PER_SEC);  // one test needs 1+N_TRK lines (one event)
    ap_uint<MTF7_BRAM_SIZE> in_bram[lines_to_read];
    // ap_uint<MTF7_BRAM_SIZE> out_bram[lines_to_read];
    ap_uint<MTF7_BRAM_SIZE> out_bram[NTEST]; // but in output all TkMu fit in one line -> out bram size is NTEST

    std::cout << "... injecting patterns" << std::endl;

    // PU 0 I/O
    // std::string in_pattern_filename = "";
    // std::string out_tkmu_filename   = "";

    // PU 200 I/O
    std::string in_pattern_filename = "/home/lcadamur/CorrPCIeGit/MuTrackCorrMTF7Utils/patterns/mu_track_patterns_stream.txt";
    std::string out_tkmu_filename   = "/home/lcadamur/MuTrackCorrHLS/MuTrackCorrHLS/corr_output_stream.txt";


    // inject the pattern data
    // patt_inj::inject_dummy<MTF7_BRAM_SIZE> (in_bram, NTEST);
    printf("[INFO] To process %i events will need to read %i lines\n", NTEST, lines_to_read);
    patt_inj::inject_load_from_file<MTF7_BRAM_SIZE> (in_bram, lines_to_read, in_pattern_filename);

    // printf("BRAM AT 0 : %x\n", in_bram[0].range(31,0).to_uint()); // LSB word
    // printf("BRAM AT 1 : %x\n", in_bram[1].range(31,0).to_uint()); // LSB word
    // printf("BRAM AT 2 : %x\n", in_bram[2].range(31,0).to_uint()); // LSB word

    // // muons of evt 2
    // printf("BRAM AT 38 : %x\n", in_bram[38].range(31,0).to_uint()); // LSB word

    // printf("- Last element of input BRAM getting printed\n");
    // for (size_t ib = 0; ib < 100; ++ib)
    // {
    //   printf("BRAM AT %i : %08x\n", ib, in_bram[ib].range(31,0).to_uint()); // LSB word
    // }

    std::cout << "... starting test on " << NTEST << " cycles" << std::endl;

    // run the test
    for (uint itest = 0; itest < NTEST; ++itest)
    {
        std::cout << ".... Test number " << itest << std::endl;

        // convert the BRAM input into hls_stream data

        // // first word: the muon data
        // muon_t muons[N_MU];

        // // the destination array
        // tkmu_t tkmus[N_TKMU];

        // size_t in_offset = itest * (1 + N_TRK_PER_SEC); // every event receives 1 line for muons + N_TRK_PER_SEC lines for tracks

        // // printf("%i --> offset %i\n", itest, in_offset);
        
        // // convert the word into a muon
        // ap_uint<MTF7_BRAM_SIZE> all_mu_word = in_bram[0 + in_offset];
        // for (size_t imu = 0; imu < N_MU; ++imu)
        // {
        //   ap_uint<MU_W_SIZE> mu_word = all_mu_word.range((imu+1)*MU_W_SIZE-1, imu*MU_W_SIZE);
        //   muons[imu].pt         = get_mu_pt(mu_word);
        //   muons[imu].theta      = get_mu_theta(mu_word);
        //   muons[imu].theta_sign = get_mu_theta_sign(mu_word);
        //   muons[imu].phi        = get_mu_phi(mu_word);
        // }

        // // for (size_t imu = 0; imu < N_MU; ++imu)
        // // {
        // //   std::cout << " ... mu " << imu
        // //             << " pt = " << muons[imu].pt.to_uint64()
        // //             << " theta = " << muons[imu].theta.to_uint64()
        // //             << " theta_sign = " << muons[imu].theta_sign.to_uint64()
        // //             << " phi = " << muons[imu].phi.to_uint64() << std::endl;
        // // }

        // // second and following word: these are the tracks
        // // each memory address carries one incoming set of track for a clk cycle
        // // NOTE: does not necessarily has to be one BRAM address / bx
        // // at the testbench level, multiple addresses can be combined if more bits are needed

        // // printf("NUMBER OF TRACK: %i\n", N_TRK);
        // hls::stream<track_t> tracks[N_TRK_SECTORS]  = {
        //     {"SECTOR_1"},
        //     {"SECTOR_2"},
        //     {"SECTOR_3"},
        //     {"SECTOR_4"},
        //     {"SECTOR_5"},
        //     {"SECTOR_6"},
        //     {"SECTOR_7"},
        //     {"SECTOR_8"},
        //     {"SECTOR_9"}
        // };

        // // first round: inject ths muons using empty track streasm
        // std::cout << ".... injecting muons" << std::endl;
        // correlator_stream(muons, tracks, tkmus);

        // // second round: write data on the input track stream
        // for (size_t itrk = 0; itrk < N_TRK_PER_SEC; ++itrk)
        // {
        //     // std::cout << " ---------- track TMT number " << itrk << " ------------" << std::endl;

        //     // for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
        //     // {
        //     //     // PATCH: set properties of tracks by hand
        //     //     // to be read from the input pattern file
        //     //     track_t this_trk;
        //     //     this_trk.pt    = 100*itrk + isec;
        //     //     this_trk.theta = 100*itrk + isec + 1;
        //     //     this_trk.phi   = 100*itrk + isec + 2;
        //     //     this_trk.nstubs = 5;
        //     //     this_trk.chisq  = 30;
        //     //     tracks[isec].write(this_trk);
        //     // }

        //     ap_uint<MTF7_BRAM_SIZE> all_trk_word = in_bram[1 + itrk + in_offset];
        //     for (size_t isec = 0; isec < N_TRK_SECTORS; ++isec)
        //     {
        //       ap_uint<TRK_W_SIZE> trk_word = all_trk_word.range((isec+1)*TRK_W_SIZE-1, isec*TRK_W_SIZE);
        //       track_t this_trk;
        //       this_trk.pt        = get_trk_pt(trk_word);
        //       this_trk.theta     = get_trk_theta(trk_word);
        //       this_trk.theta_sign = get_trk_theta_sign(trk_word);
        //       this_trk.phi       = get_trk_phi(trk_word);
        //       this_trk.charge    = get_trk_charge(trk_word);
        //       this_trk.chisq     = get_trk_chisq(trk_word);
        //       this_trk.nstubs    = get_trk_nstubs(trk_word);

        //       // std::cout << " ... tracks - sector = " << isec << " tmt = " << itrk
        //       //           << " pt = "         << this_trk.pt.to_uint64()
        //       //           << " theta = "      << this_trk.theta.to_uint64()
        //       //           << " theta_sign = " << this_trk.theta_sign.to_uint64()
        //       //           << " phi = "        << this_trk.phi.to_uint64()
        //       //           << " charge = "     << this_trk.charge.to_uint64()
        //       //           << " chisq = "      << this_trk.chisq.to_uint64()
        //       //           << " nstubs = "     << this_trk.nstubs.to_uint64()
        //       //           << std::endl;

        //       tracks[isec].write(this_trk);
        //     }

        //     // std::cout << ".... running correlator on a new set of tracks" << std::endl;
        //     correlator_stream(muons, tracks, tkmus);
        // }

        // for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
        // {
        //   std::cout << " ... tkmu " << itkmu
        //             << " pt = " << tkmus[itkmu].pt.to_uint64()
        //             << " theta = " << tkmus[itkmu].theta.to_uint64()
        //             << " theta_sign = " << tkmus[itkmu].theta_sign.to_uint64()
        //             << " phi = " << tkmus[itkmu].phi.to_uint64() << std::endl;
        // }


        // /// repack the output in out_bram and write it to file in output
        // ap_uint<MTF7_BRAM_SIZE> out_info;
        // for (size_t itkmu = 0; itkmu < N_TKMU; ++itkmu)
        // {
        //   // #pragma HLS unroll
        //   out_info.range(TKMU_W_SIZE*(itkmu+1)-1, TKMU_W_SIZE*itkmu) = build_tkmu_word(tkmus[itkmu]);
        // }
        
        // out_bram[itest] = out_info;

        ap_uint<1> out_valid = 0;
        ap_uint<MTF7_BRAM_SIZE> in_info = 0x0;
        ap_uint<MTF7_BRAM_SIZE> out_info = 0x0;
        size_t in_offset = itest * (1 + N_TRK_PER_SEC); // every event receives 1 line for muons + N_TRK_PER_SEC lines for tracks
        for (size_t iline = 0; iline < 1 + N_TRK_PER_SEC; ++iline)
        {


          in_info = in_bram[iline + in_offset];
          ap_uint<1> is_muons = (iline == 0 ? 1 : 0);
          // BRAM_to_corr (in_info, is_muons, out_info, out_valid);

          // std::cout << " ------ line nr " << iline
          //           << " in_info_last : " << std::hex << in_info.range(63,0).to_uint64()
          //           << " is_muons : " << is_muons.to_uint()
          //           << " out_info_last : " << std::hex << out_info.range(63,0).to_uint64()
          //           << " out_valid : " << out_valid 
          //           << std::endl;
          // std::cout << "ciao " << iline << " " << (unsigned int) in_info.range(31,0) << std::endl;
          // printf("%0x\n", in_info.range(31,0).to_uint64() );

          BRAM_to_corr_nostream (in_info, is_muons, out_info, out_valid);

          // std::cout << " .... post exec ... ------ line nr " << iline
          //           << " in_info_last : " << std::hex << in_info.range(63,0).to_uint64()
          //           << " is_muons : " << is_muons.to_uint()
          //           << " out_info_last : " << std::hex << out_info.range(63,0).to_uint64()
          //           << " out_valid : " << out_valid 
          //           << std::endl;


          // std::cout << " -------------- out_valid : " << out_valid.to_uint() << " out_info_last : " << std::hex << out_info.range(63,0).to_uint64() << std::endl;

          if (out_valid) // out_valid is raised while transmitting the tracks from the next event because of latency
            out_bram[itest] = out_info;
        }

        // std::cout << "AFTER RUNNING THE CORRELATOR out_valid = " << out_valid.to_uint() << " last part of word : " << out_info.range(63,0).to_uint64() << std::endl;

        // if (out_valid)
        // // if (true)
        // {
        //   out_bram[itest] = out_info;
        // }
    }

    // for (size_t ib = 0; ib < NTEST; ++ib)
    //   std::cout << "cross check " << out_bram[ib].range(63,0).to_uint64() << std::endl;

    printf("... dumping result to file\n");
    // note: normally file woudl go under "projCorrelator/solution1/csim/build/"
    // so give here an abs path
    dump_output_to_file<MTF7_BRAM_SIZE> (out_bram, NTEST, out_tkmu_filename.c_str());
    printf("... done\n");

}
