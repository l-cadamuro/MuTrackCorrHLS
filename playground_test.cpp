#include <cstdio>
#include <iostream>
#include "src/playground.h"
#include "ap_int.h"
#include "ap_fixed.h"

#define NTEST 20
#define NTOSUM 10

int main()
{
    /*
    // ap_uint<16> inputs[NTEST];
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> inputs[NTEST];
    for (uint i = 0; i < NTEST; ++i)
        inputs[i] = i;

    // ap_uint<16> outputs[NTEST];
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> outputs[NTEST];
    for (uint i = 0; i < NTEST; ++i)
        outputs[i] = 0;

    for (uint itest = 0; itest < NTEST; ++itest)
    {
        pipelined_transfer_corrdf(inputs[itest], outputs[itest]);
        printf(">> Clock %i :: in -> %i, out %i\n", itest, inputs[itest].to_int(), outputs[itest].to_int());
    }
    */

    // for this test:
    // write 3 numbers in the input word
    // 2 1 0
    // 3 2 1
    // 4 3 2
    // ...
    /// and so on for NTOSUM lines

    ap_uint<2048> input_word;
    ap_uint<2048> output_word;
    ap_uint<32> offset = 0;

    for (uint itest = 0; itest < NTEST; ++itest)
    {

        std::cout << "====== TEST CYCLE ==== " << itest << std::endl;

        input_word = 0;
        output_word = 0;

        ap_uint<32> p1 = offset + 0;
        ap_uint<32> p2 = offset + 1;
        ap_uint<32> p3 = offset + 2;

        ap_uint<1> start = 0;
        ap_uint<1> run   = 0;
        ap_uint<1> stop  = 0;

        for (uint ipart = 0; ipart < NTOSUM; ++ipart)
        {
            start = (ipart == 0 ? 1 : 0);
            stop  = (ipart == NTOSUM-1 ? 1 : 0);
            run   = 1; 

            input_word.range(31,0) = p1;
            input_word.range(31+32,0+32) = p2;
            input_word.range(31+32+32,0+32+32) = p3;
            input_word.range(96,96) = start;
            input_word.range(97,97) = run;
            input_word.range(98,98) = stop;

            stream_wrapper (input_word, output_word);

            // check the result
            ap_uint<512> sum1;
            ap_uint<512> sum2;
            ap_uint<512> sum3;
            
            sum1 = output_word.range(511, 0);
            sum2 = output_word.range(511+512, 0+512);
            sum3 = output_word.range(511+512+512, 0+512+512);

            std::cout << ".. part " << ipart << " -- start = " << start.to_uint() << " , run = " << run.to_uint() << " , stop = " << stop.to_uint() << std::endl;
            std::cout << "....... p1 = " << p1.to_uint() << " -- p2 = " << p2.to_uint() << " , p3 = " << p3.to_uint() << std::endl;
            std::cout << "Sum 1 = " << sum1.to_uint64() << std::endl;
            std::cout << "Sum 2 = " << sum2.to_uint64() << std::endl;
            std::cout << "Sum 3 = " << sum3.to_uint64() << std::endl;

            p1 += 1;
            p2 += 1;
            p3 += 1;
        }

        offset += 10;

        std::cout << std::endl << std::endl << std::endl;
    }
}