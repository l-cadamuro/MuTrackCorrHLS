#include <cstdio>
#include "src/algo_parts.h"
#include "src/playground.h"
#include "ap_int.h"
#include "ap_fixed.h"

#define NTEST 5


/*
int main() {

    // ap_int<32> b_low;
    // ap_int<32> b_high;
    // ap_int<32> val;
    ap_int<THETA_W_SIZE> b_low;
    ap_int<THETA_W_SIZE> b_high;
    ap_int<THETA_W_SIZE> val;
    ap_uint<1> out_bool;

    ap_fixed<16,9> fp_nr;

    b_low  = 18;
    b_high = 20;
    val    = 17;

    fp_nr = 3.1415926536;

    for (int test = 1; test <= NTEST; ++test) {

        // b_low  = 18;
        // b_high = 28;
        // val    = 24;

        // is_in_boundaries<32>(val, b_low, b_high, out_bool);
        is_in_boundaries_th(val, b_low, b_high, out_bool);
        printf( ">>> %i ) %i <= %i < %i ? ... %i \n", test, int(b_low), int(val), int(b_high), int(out_bool) );

        float start = 0.0;
        float end   = 1.0;
        float val   = start + 1.*test*(end-start)/NTEST;
        // fp_nr = val;
        fp_nr = 2.9;

        printf("My fixed point number : %f -> %f, integer = %i, double = %f \n", float(val), float(fp_nr), fp_nr.to_int(), fp_nr.to_double());

        val = val+1;
    }

    return 0;
}
*/

// int main()
// {
//     ap_int<8> v_in = 0;
//     ap_int<8> v_out;

//     for (int test = 1; test <= NTEST; ++test) {
//         // func_with_inner_val(v_in, v_out);
//         func_with_inner_val_present(v_in, v_out, (test == 1 ? 1 : 0));
//         printf(">> %i : IN = %i, OUT = %i\n", test, int(v_in), int(v_out));
//     }
// }


// // --------------------------------------------------------

// // a model of the in/out of the tracks
// int main()
// {
//     ap_uint<MU_W_SIZE> muons [N_MU];
//     ap_uint<TRK_W_SIZE> tracks[N_TRK];
//     ap_uint<TKMU_W_SIZE> tkmus [N_TKMU];
    
//     ap_int<10> prova;
//     prova = 100;
//     printf("PROVA = %i\n", int(prova));

//     for (int test = 1; test <= NTEST; ++test) {

//         printf("... iteration %i\n", test);
//         for (int i = 0; i < N_TKMU; ++i)
//         {
//             tkmus[i] = 0;
//         }

//         build_tkmu(muons, tracks, tkmus);

//         for (int i = 0; i < N_TKMU; ++i)
//             printf("%i) -> %i\n", i, int(tkmus[i]));
//     }
// }


// --------------------------------------------------------

// study parallelism
int main()
{
    ap_int<16> arr_in[50];
    for (int i = 0; i < 50; ++i)
    {
        arr_in[i] = 2;
    }
    ap_int<16> sum_out = 0; //ffff
    printf("SUM IS (before) = %i\n", int(sum_out));

    sum_array(arr_in, sum_out);

    printf("SUM IS = %i\n", int(sum_out));
}
