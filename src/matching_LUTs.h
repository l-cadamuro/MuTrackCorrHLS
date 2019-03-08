#ifndef MATCHING_LUTS_H
#define MATCHING_LUTS_H
#include "ap_int.h"

struct matching_LUTs {
   const ap_uint<10> phi_low_bounds [16][512] ;
   const ap_uint<10> phi_high_bounds [16][512] ;
   const ap_uint<10> theta_low_bounds [16][512] ;
   const ap_uint<10> theta_high_bounds [16][512] ;
};

#endif
