#ifndef DATAFORMATS_H
#define DATAFORMATS_H

#include "ap_int.h"

// the input data formats properties are defined here
// FIXME: move to a specific file to include here?
#define MU_W_SIZE 40
#define TRK_W_SIZE 100
#define TKMU_W_SIZE 100

// size of input and output collections
#define N_MU 12
#define N_TRK 15
#define N_TKMU 12 // by default, the same as N_MU

// FIXME: set to meaningful values

// the parts of the word - you have to ensure that they fit in the declared variable sizes
// assume padding HSB --> [. THETA .][. PHI .][. PT .] <-- LSB

// NOTE: pt    is unsigned (0+)
// NOTE: phi   is unsigned (0 - 2pi)
// NOTE: theta is signed   (+/-)

#define MU_THETA_W_SIZE 6 // includes sign bit
#define MU_PHI_W_SIZE 6
#define MU_PT_W_SIZE 15

#define TRK_THETA_W_SIZE 6 // includes sign bit
#define TRK_PHI_W_SIZE 6
#define TRK_PT_W_SIZE 15

#define TKMU_THETA_W_SIZE 6 // includes sign bit
#define TKMU_PHI_W_SIZE 6
#define TKMU_PT_W_SIZE 15

// FIXME: to be made uniform
#define THETA_W_SIZE 6
#define PHI_W_SIZE 6

// unpackers
ap_uint<MU_PT_W_SIZE> get_mu_pt(ap_uint<MU_W_SIZE> word);
ap_uint<MU_PHI_W_SIZE> get_mu_phi(ap_uint<MU_W_SIZE> word);
ap_int<MU_THETA_W_SIZE> get_mu_theta(ap_uint<MU_W_SIZE> word);
////////
ap_uint<TRK_PT_W_SIZE> get_trk_pt(ap_uint<TRK_W_SIZE> word);
ap_uint<TRK_PHI_W_SIZE> get_trk_phi(ap_uint<TRK_W_SIZE> word);
ap_int<TRK_THETA_W_SIZE> get_trk_theta(ap_uint<TRK_W_SIZE> word);
////////
ap_uint<TKMU_PT_W_SIZE> get_tkmu_pt(ap_uint<TKMU_W_SIZE> word);
ap_uint<TKMU_PHI_W_SIZE> get_tkmu_phi(ap_uint<TKMU_W_SIZE> word);
ap_int<TKMU_THETA_W_SIZE> get_tkmu_theta(ap_uint<TKMU_W_SIZE> word);

#endif