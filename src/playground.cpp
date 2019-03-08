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
    registers[nreg-1] = in_info + phi_high_bounds[0][10]; ///0, 10


    
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