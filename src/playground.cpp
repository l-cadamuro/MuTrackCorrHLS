#include "playground.h"
#include "matching_LUTs.h"


// NOTE: if I only put the source code below in playground.h, csym works fine but csynth fails saying that "func_with_inner_val" is not in the top function
// I guess vivado cares abotu the compiled code in this src and cannot build the top function out of the interface declaration in the .h file

ap_int<8> my_ext_nr = 0; // putting this instruction in the .h file makes compilation fail, why?

void func_with_inner_val(ap_int<8> v_in, ap_int<8> &v_out)
{
    v_out = v_in + my_ext_nr;
    my_ext_nr += 1;
}


// the function below does not work properly, the variable defined inside gets re-initialised every time
void func_with_inner_val_present(ap_int<8> v_in, ap_int<8> &v_out, ap_uint<1> is_first)
{
    ap_int<8> my_int_nr;
    
    if (is_first)
        my_int_nr = 0;
    else
        my_int_nr += 1;
    
    v_out = v_in + my_int_nr;
}

// sums up all the elements in the array and returns the sum out
void sum_array(ap_int<16> arr_in[50], ap_int<16> &sum_out)
{
    // ap_int<16> sum = 0;    
    for (int i = 0; i < 50; ++i)
    {
        sum_out += arr_in[i];
    }
    // sum_out = sum;
}

void pipelined_transfer_corrdf(
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info,
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> &out_info)
    // ap_uint<TKMU_W_SIZE*N_TKMU> &out_info)
{
    // simply shift the full word across N register
    const int nreg = 4;
    static ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> registers [nreg];

    out_info = registers[0];
    for (size_t ireg = 0; ireg < nreg-1; ++ireg)
    {
        #pragma HLS unroll
        registers[ireg] = registers[ireg+1];
    }
    registers[nreg-1] = in_info;
    // registers[nreg-1] = in_info + 1;
    // registers[nreg-1] = in_info + phi_high_bounds[0][10]; ///0, 10

}

void pipelined_transfer(ap_uint<16> in_info, ap_uint<16> &out_info)
{
    /*
    const static ap_int<16> phi_high_bounds [3][3] = {
        {
            0x1,
            0x2,
            0x3,
        },
        {
            0x10,
            0x20,
            0x30,
        },
        {
            0x100,
            0x200,
            0x300,
        },
    };
    */

    // data here get "hold" for a long time
    /*
    const int nreg = 4;
    static ap_uint<16> registers [nreg]; // 8 steps -> 8 clk latency

    out_info = registers[0];
    for (size_t ireg = 0; ireg < nreg-1; ++ireg)
    {
        #pragma HLS pipeline
        registers[ireg] = registers[ireg+1];        
    }
    registers[nreg-1] = in_info;
    */

    /*
    // trying the rewind option and more, they do not work
    const int nreg = 4;
    static ap_uint<16> registers [nreg]; // 8 steps -> 8 clk latency

    out_info = registers[0];
    for (size_t ireg = 0; ireg < nreg-1; ++ireg)
    {
        // #pragma HLS pipeline II=1 enable_flush rewind
        #pragma HLS pipeline
        registers[ireg] = registers[ireg+1];        
    }
    registers[nreg-1] = in_info;
    */

    /*
    const int nreg = 4;
    static ap_uint<16> registers [nreg];

    for (size_t ireg = 0; ireg < nreg; ++ireg)
    {
        #pragma HLS pipeline II=1 enable_flush rewind

        if (ireg == 0){
            out_info = registers[ireg];
            registers[ireg] = registers[ireg+1];
        }
        else if (ireg < nreg-1){
            registers[ireg] = registers[ireg+1];
        }
        else{
            registers[ireg] = in_info;
        }
    }
    */


    /*
    // trying what I woudl do in an HDL - linking wires before the assignments
    // --> this gives me 3 clk latency intead of 4 as in the example [example_working] below
    // plus is not working since values are still hold in memory for a long time
    const int nreg = 4;
    static ap_uint<16> registers [nreg];


    out_info = registers[0];
    registers[nreg-1] = in_info;

    for (size_t ireg = 0; ireg < nreg-1; ++ireg)
    {
        #pragma HLS pipeline II=1 enable_flush rewind
        registers[ireg] = registers[ireg+1];        
    }
    */


    /*
    const int nreg = 4;
    static ap_uint<16> registers [nreg];
    #pragma HLS ARRAY_PARTITION variable=registers complete dim=0

    out_info = registers[0];
    for (size_t ireg = 0; ireg < nreg-1; ++ireg)
    {
        #pragma HLS pipeline
        registers[ireg] = registers[ireg+1];        
    }
    registers[nreg-1] = in_info;
    */

    // trying to unroll instead
    const int nreg = 4;
    static ap_uint<16> registers [nreg];

    out_info = registers[0];
    for (size_t ireg = 0; ireg < nreg-1; ++ireg)
    {
        #pragma HLS unroll
        registers[ireg] = registers[ireg+1];
    }
    // registers[nreg-1] = in_info + 1;
    static matching_LUTs mluts;
    registers[nreg-1] = in_info + mluts.phi_high_bounds[0][10]; ///0, 10


    
    /*
    // try to make the chain by hand - 4 registers
    // [example_working]
    static ap_uint<16> r0;
    static ap_uint<16> r1;
    static ap_uint<16> r2;
    static ap_uint<16> r3;

    out_info = r0;
    r0 = r1;
    r1 = r2;
    r2 = r3;
    r3 = in_info;
    */
}

void stream_wrapper (ap_uint<2048> input, ap_uint<2048> &output)
{
    ap_uint<1> start, run, stop;
    ap_uint<32> w1, w2, w3;
    unpacker(input, start, run, stop, w1, w2, w3);

    ap_uint<512> sum1;
    ap_uint<512> sum2;
    ap_uint<512> sum3;

    // the FIFOs    
    static hls::stream<ap_uint<1> > s_start ("s_start");
    static hls::stream<ap_uint<1> > s_run   ("s_run");
    static hls::stream<ap_uint<1> > s_stop  ("s_stop");
    static hls::stream<ap_uint<32> > s_w1   ("s_w1");
    static hls::stream<ap_uint<32> > s_w2   ("s_w2");
    static hls::stream<ap_uint<32> > s_w3   ("s_w3");

    #pragma HLS STREAM variable = s_start depth = 4
    #pragma HLS STREAM variable = s_run   depth = 4
    #pragma HLS STREAM variable = s_stop  depth = 4
    #pragma HLS STREAM variable = s_w1    depth = 4
    #pragma HLS STREAM variable = s_w2    depth = 4
    #pragma HLS STREAM variable = s_w3    depth = 4

    s_start.write(start);
    s_run.write(run);
    s_stop.write(stop);
    s_w1.write(w1);
    s_w2.write(w2);
    s_w3.write(w3);

    worker(s_start, s_run, s_stop, s_w1, s_w2, s_w3, sum1, sum2, sum3);

    repacker(sum1, sum2, sum3, output);
}

///////////////////////////////////
/////// stream test
///////////////////////////////////

void unpacker(ap_uint<2048> input_word, ap_uint<1> &start, ap_uint<1> &run, ap_uint<1> &end, ap_uint<32> &out_w1, ap_uint<32> &out_w2, ap_uint<32> &out_w3)
{
    #pragma HLS INLINE

    out_w1 = input_word.range(31,0);
    out_w2 = input_word.range(63,32);
    out_w3 = input_word.range(95,64);
    start  = input_word.range(96,96);
    run    = input_word.range(97,97);
    end    = input_word.range(98,98);
}

void worker(
    hls::stream<ap_uint<1> > &start, hls::stream<ap_uint<1> > &run, hls::stream<ap_uint<1> > &stop,
    hls::stream<ap_uint<32> > &link1, hls::stream<ap_uint<32> > &link2, hls::stream<ap_uint<32> > &link3,
    ap_uint<512> &sum1, ap_uint<512> &sum2, ap_uint<512> &sum3) // the worker just sums the streamed content for each link 
{
    // temporary internal sums
    static ap_uint<512> part_sum1;
    static ap_uint<512> part_sum2;
    static ap_uint<512> part_sum3;

    // output buffer registers
    static ap_uint<512> buf_sum1; 
    static ap_uint<512> buf_sum2; 
    static ap_uint<512> buf_sum3;

    // attach the output to the buffers
    sum1 = buf_sum1;
    sum2 = buf_sum2;
    sum3 = buf_sum3;

    // // -------------------- with delay
    // ap_uint<32> this_w1_in   = link1.read();    
    // ap_uint<32> this_w2_in   = link2.read();    
    // ap_uint<32> this_w3_in   = link3.read();

    // ap_uint<1> this_start_in = start.read();
    // ap_uint<1> this_run_in   = run.read();
    // ap_uint<1> this_stop_in  = stop.read();

    // ap_uint<1> this_start;
    // ap_uint<1> this_run;  
    // ap_uint<1> this_stop; 

    // ap_uint<32> this_w1;
    // ap_uint<32> this_w2;
    // ap_uint<32> this_w3;

    // // NB: if I reuse pipeline_delay, the static member is going to be reused everi time, generating a big confusion :-(
    // // my current shortcut is to make 6 copies of this function
    // pipeline_delay<5, 1> (this_start_in, this_start);
    // pipeline_delay_1<5, 1> (this_run_in,   this_run);
    // pipeline_delay_2<5, 1> (this_stop_in,  this_stop);

    // pipeline_delay_3<5, 32> (this_w1_in, this_w1);
    // pipeline_delay_4<5, 32> (this_w2_in, this_w2);
    // pipeline_delay_5<5, 32> (this_w3_in, this_w3);

    // -------------------- no delay
    ap_uint<1> this_start = start.read();
    ap_uint<1> this_run   = run.read();
    ap_uint<1> this_stop  = stop.read();

    ap_uint<32> this_w1  = link1.read();
    ap_uint<32> this_w2  = link2.read();
    ap_uint<32> this_w3  = link3.read();

    // #ifndef __SYNTHESIS__
    // std::cout
    //     << "w1 = " << this_w1.to_uint() << " "
    //     << "w2 = " << this_w2.to_uint() << " "
    //     << "w3 = " << this_w3.to_uint() << " "
    // << std::endl;
    // #endif

    // reset all
    if (this_start == 1)
    {
        part_sum1 = 0;
        part_sum2 = 0;
        part_sum3 = 0;
    }

    if (this_run == 1)
    {
        part_sum1 += this_w1;
        part_sum2 += this_w2;
        part_sum3 += this_w3;        
    }

    if (this_stop == 1)
    {
        // copy the values
        buf_sum1 = part_sum1;
        buf_sum2 = part_sum2;
        buf_sum3 = part_sum3;
    }
    // else
    // {
    //     buf_sum1 = buf_sum1; // hold ancient value
    //     buf_sum2 = buf_sum2; // hold ancient value
    //     buf_sum3 = buf_sum3; // hold ancient value
    // }

}

void repacker(ap_uint<512> in_w1, ap_uint<512> in_w2, ap_uint<512> in_w3, ap_uint<2048> &output_word)
{
    #pragma HLS INLINE
    
    output_word.range(511, 0) = in_w1;
    output_word.range(511+512, 0+512) = in_w2;
    output_word.range(511+512+512, 0+512+512) = in_w3;
}
