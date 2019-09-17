#ifndef DATAFORMATS_H
#define DATAFORMATS_H

#define AP_INT_MAX_W 2048
// #define AP_INT_MAX_W 16384 // just to fit for test the 9 sectors --> simulation becomes too slow
#include "ap_int.h"

// the input data formats properties are defined here
// FIXME: move to a specific file to include here?
#define MU_W_SIZE 40
// #define MU_W_SIZE 100 // Darin's suggestion
#define TRK_W_SIZE 100
#define TKMU_W_SIZE 100

// size of input and output collections
#define N_MU 12
#define N_TRK 15
// #define N_TRK 8 // with 100 bits per muons, reduce the number of tracks to fit in the 2048 bits
#define N_TKMU 12 // by default, the same as N_MU

#define N_TRK_SECTORS 9 // the nine tracker nonets
#define N_TRK_PER_SEC 18 // number of tracks from each sector

// FIXME: set to meaningful values

// the parts of the word - you have to ensure that they fit in the declared variable sizes
// for the track:
// assume padding HSB --> [. NSTUBS.][. CHI SQ.][. CHARGE.][. THETA SIGN.][. THETA .][. PHI .][. PT .] <-- LSB
// for the muon
// assume padding HSB -->                       [. CHARGE.][. THETA SIGN.][. THETA .][. PHI .][. PT .] <-- LSB

// NOTE: pt    is unsigned (0+)
// NOTE: phi   is unsigned (0 - 2pi)
// NOTE: theta is unsigned + sign bit  (+/-)

#define MU_THETA_W_SIZE 11 // does not include sign bit
#define MU_PHI_W_SIZE 12
#define MU_PT_W_SIZE 15

#define TRK_THETA_W_SIZE 15 // does not include sign bit
#define TRK_PHI_W_SIZE 15
#define TRK_PT_W_SIZE 15
#define TRK_CHISQ_W_SIZE 10
#define TRK_NSTUBS_W_SIZE 10

#define TKMU_THETA_W_SIZE 15 // includes sign bit
#define TKMU_PHI_W_SIZE 15
#define TKMU_PT_W_SIZE 15

// FIXME: to be made uniform
// this is used for the deltas of tkmu and mu
#define THETA_W_SIZE 11
#define PHI_W_SIZE 11
#define MW_LUT_ANGLE_W_SIZE 11 // number of bits used for the matching LUTs boundaries. FIXME: now only works if == theta_w_size!!

// for the interface with the BRAM
#define MTF7_BRAM_SIZE 2048

// -----------------------------------------------------------
// the structures defining input objects
struct muon_t {
    ap_uint<MU_PT_W_SIZE>    pt;
    ap_uint<MU_THETA_W_SIZE> theta;
    ap_uint<1>               theta_sign;
    ap_uint<MU_PHI_W_SIZE>   phi;
};
inline void clear(muon_t & c) {
    c.pt = 0; c.theta = 0; c.theta_sign = 0; c.phi = 0;}


struct track_t {
    ap_uint<TRK_PT_W_SIZE>     pt;
    ap_uint<TRK_THETA_W_SIZE>  theta;
    ap_uint<1>                 theta_sign;
    ap_uint<1>                 charge;
    ap_uint<TRK_PHI_W_SIZE>    phi;
    ap_uint<TRK_CHISQ_W_SIZE>  chisq;
    ap_uint<TRK_NSTUBS_W_SIZE> nstubs;
};
inline void clear(track_t & c) {
    c.pt = 0; c.theta = 0; c.theta_sign = 0; c.charge = 0; c.phi = 0; c.chisq = 0; c.nstubs = 0;}


struct tkmu_t {
    ap_uint<TKMU_PT_W_SIZE>     pt;
    ap_uint<TKMU_THETA_W_SIZE>  theta;
    ap_uint<1>                  theta_sign;
    ap_uint<1>                  charge;
    ap_uint<TKMU_PHI_W_SIZE>    phi;
};
inline void clear(tkmu_t & c) {
    c.pt = 0; c.theta = 0; c.theta_sign = 0; c.charge = 0; c.phi = 0;}


// -----------------------------------------------------------

// unpackers
ap_uint<MU_PT_W_SIZE>    get_mu_pt(ap_uint<MU_W_SIZE> word);
ap_uint<MU_PHI_W_SIZE>   get_mu_phi(ap_uint<MU_W_SIZE> word);
ap_uint<MU_THETA_W_SIZE> get_mu_theta(ap_uint<MU_W_SIZE> word);
ap_uint<1>               get_mu_theta_sign(ap_uint<MU_W_SIZE> word);
ap_uint<1>               get_mu_charge(ap_uint<MU_W_SIZE> word);
////////
ap_uint<TRK_PT_W_SIZE>     get_trk_pt(ap_uint<TRK_W_SIZE> word);
ap_uint<TRK_PHI_W_SIZE>    get_trk_phi(ap_uint<TRK_W_SIZE> word);
ap_uint<TRK_THETA_W_SIZE>  get_trk_theta(ap_uint<TRK_W_SIZE> word);
ap_uint<1>                 get_trk_theta_sign(ap_uint<TRK_W_SIZE> word);
ap_uint<1>                 get_trk_charge(ap_uint<TRK_W_SIZE> word);
ap_uint<TRK_CHISQ_W_SIZE>  get_trk_chisq(ap_uint<TRK_W_SIZE> word);
ap_uint<TRK_NSTUBS_W_SIZE> get_trk_nstubs(ap_uint<TRK_W_SIZE> word);
////////
ap_uint<TKMU_PT_W_SIZE>    get_tkmu_pt(ap_uint<TKMU_W_SIZE> word);
ap_uint<TKMU_PHI_W_SIZE>   get_tkmu_phi(ap_uint<TKMU_W_SIZE> word);
ap_uint<TKMU_THETA_W_SIZE> get_tkmu_theta(ap_uint<TKMU_W_SIZE> word);
ap_uint<1>                 get_tkmu_theta_sign(ap_uint<TKMU_W_SIZE> word);
ap_uint<1>                 get_tkmu_charge(ap_uint<TKMU_W_SIZE> word);

// packers
ap_uint<TKMU_W_SIZE>       build_tkmu_word(tkmu_t tkmu);

// ref getters (nota: ap_range_ref<full_width, is_signed> )
ap_range_ref<TKMU_W_SIZE, false> ref_to_tkmu_pt(ap_uint<TKMU_W_SIZE>& word);
ap_range_ref<TKMU_W_SIZE, false> ref_to_tkmu_phi(ap_uint<TKMU_W_SIZE>& word);
ap_range_ref<TKMU_W_SIZE, false> ref_to_tkmu_theta(ap_uint<TKMU_W_SIZE>& word);

#endif