#include <cstdio>
#include "src/algo_parts.h"
#include "ap_int.h"

#define NTEST 5


int main() {

    // ap_int<32> b_low;
    // ap_int<32> b_high;
    // ap_int<32> val;
    ap_int<THETA_W_SIZE> b_low;
    ap_int<THETA_W_SIZE> b_high;
    ap_int<THETA_W_SIZE> val;
    ap_uint<1> out_bool;

    b_low  = 18;
    b_high = 20;
    val    = 17;

    for (int test = 1; test <= NTEST; ++test) {

        // b_low  = 18;
        // b_high = 28;
        // val    = 24;

        // is_in_boundaries<32>(val, b_low, b_high, out_bool);
        is_in_boundaries_th(val, b_low, b_high, out_bool);
        printf( ">>> %i ) %i <= %i < %i ? ... %i \n", test, int(b_low), int(val), int(b_high), int(out_bool) );

        val = val+1;
    }

    return 0;
}