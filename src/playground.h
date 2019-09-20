// this file contains a series of unmaintained functins, code, etc... that I use to test and learn HLS
// not to be used as part of the MuTk algorithm

#ifndef PLAYGROUND_H
#define PLAYGROUND_H

#define AP_INT_MAX_W 2048

#include "ap_fixed.h"
#include "ap_int.h"
#include "hls_stream.h"

// the input data formats properties are defined here
// FIXME: move to a specific file to include here?
#define MU_W_SIZE 40
#define TRK_W_SIZE 100
#define TKMU_W_SIZE 100

// size of input and output collections
#define N_MU 12
#define N_TRK 15
#define N_TKMU 12 // by default, the same as N_MU


void func_with_inner_val(ap_int<8> v_in, ap_int<8> &v_out);
void func_with_inner_val_present(ap_int<8> v_in, ap_int<8> &v_out, ap_uint<1> is_first);

void sum_array(ap_int<16> arr_in[50], ap_int<16> &sum_out);

void pipelined_transfer(ap_uint<16> in_info, ap_uint<16> &out_info);
void pipelined_transfer_corrdf(
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info,
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> &out_info);

// toy of the MTF7 worflow : unpack pattern, stream to correlator, repack
void stream_wrapper (ap_uint<2048> input, ap_uint<2048> &output); // my top function
void unpacker(ap_uint<2048> input_word, ap_uint<1> &start, ap_uint<1> &run, ap_uint<1> &end, ap_uint<32> &out_w1, ap_uint<32> &out_w2, ap_uint<32> &out_w3);
void worker(hls::stream<ap_uint<1> > &start, hls::stream<ap_uint<1> > &run, hls::stream<ap_uint<1> > &stop,
    hls::stream<ap_uint<32> > &link1, hls::stream<ap_uint<32> > &link2, hls::stream<ap_uint<32> > &link3,
    ap_uint<512> &sum1, ap_uint<512> &sum2, ap_uint<512> &sum3); // the worker just sums the streamed content for each link 
void repacker(ap_uint<512> in_w1, ap_uint<512> in_w2, ap_uint<512> in_w3, ap_uint<2048> &output_word);



template <unsigned int depth, unsigned int Tsize>
void pipeline_delay(ap_uint<Tsize> in, ap_uint<Tsize> &out)
{
    static ap_uint<Tsize> pipe[depth];
    out = pipe[0];
    for (size_t ireg = 0; ireg < depth-1; ++ireg)
    {
        #pragma HLS unroll
        pipe[ireg] = pipe[ireg+1];
    }
    pipe[depth-1] = in;
}

template <unsigned int depth, unsigned int Tsize>
void pipeline_delay_1(ap_uint<Tsize> in, ap_uint<Tsize> &out)
{
    static ap_uint<Tsize> pipe[depth];
    out = pipe[0];
    for (size_t ireg = 0; ireg < depth-1; ++ireg)
    {
        #pragma HLS unroll
        pipe[ireg] = pipe[ireg+1];
    }
    pipe[depth-1] = in;
}

template <unsigned int depth, unsigned int Tsize>
void pipeline_delay_2(ap_uint<Tsize> in, ap_uint<Tsize> &out)
{
    static ap_uint<Tsize> pipe[depth];
    out = pipe[0];
    for (size_t ireg = 0; ireg < depth-1; ++ireg)
    {
        #pragma HLS unroll
        pipe[ireg] = pipe[ireg+1];
    }
    pipe[depth-1] = in;
}

template <unsigned int depth, unsigned int Tsize>
void pipeline_delay_3(ap_uint<Tsize> in, ap_uint<Tsize> &out)
{
    static ap_uint<Tsize> pipe[depth];
    out = pipe[0];
    for (size_t ireg = 0; ireg < depth-1; ++ireg)
    {
        #pragma HLS unroll
        pipe[ireg] = pipe[ireg+1];
    }
    pipe[depth-1] = in;
}

template <unsigned int depth, unsigned int Tsize>
void pipeline_delay_4(ap_uint<Tsize> in, ap_uint<Tsize> &out)
{
    static ap_uint<Tsize> pipe[depth];
    out = pipe[0];
    for (size_t ireg = 0; ireg < depth-1; ++ireg)
    {
        #pragma HLS unroll
        pipe[ireg] = pipe[ireg+1];
    }
    pipe[depth-1] = in;
}

template <unsigned int depth, unsigned int Tsize>
void pipeline_delay_5(ap_uint<Tsize> in, ap_uint<Tsize> &out)
{
    static ap_uint<Tsize> pipe[depth];
    out = pipe[0];
    for (size_t ireg = 0; ireg < depth-1; ++ireg)
    {
        #pragma HLS unroll
        pipe[ireg] = pipe[ireg+1];
    }
    pipe[depth-1] = in;
}


// void pipelined_transfer_l(ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info,
//     ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
//     ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
//     ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2);

#endif