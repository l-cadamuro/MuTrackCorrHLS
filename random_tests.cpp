#include <cstdio>
#include "src/algo_parts.h"
#include "ap_int.h"

#define NTEST 1


int main() {

    // ap_int<32> b_low;
    // ap_int<32> b_high;
    // ap_int<32> val;
    ap_int<THETA_W_SIZE> b_low;
    ap_int<THETA_W_SIZE> b_high;
    ap_int<THETA_W_SIZE> val;
    ap_uint<1> out_bool;

    for (int test = 1; test <= NTEST; ++test) {

        b_low  = 55;
        b_high = 130;
        val    = 100;

        // is_in_boundaries<32>(val, b_low, b_high, out_bool);
        is_in_boundaries_th(val, b_low, b_high, out_bool);
        printf( "%i <= %i < %i ? ... %i \n", int(b_low), int(val), int(b_high), int(out_bool) );
    }

    return 0;
}