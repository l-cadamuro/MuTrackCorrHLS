#include "simple_algo.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

void simple_algo_hw( ap_int<32> inA, ap_int<32> inB, ap_int<32> &outA ) {

		ap_int<32> offset = 15;
		outA = inA + inB + offset;

}
