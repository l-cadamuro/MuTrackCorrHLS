#include "dataformats.h"

// unpackers
ap_uint<MU_PT_W_SIZE> get_mu_pt(ap_uint<MU_W_SIZE> word){
    return word.range(MU_PT_W_SIZE-1, 0);
}
ap_uint<MU_PHI_W_SIZE> get_mu_phi(ap_uint<MU_W_SIZE> word){
    return word.range(MU_PT_W_SIZE + MU_PHI_W_SIZE -1, MU_PT_W_SIZE);
}
ap_int<MU_THETA_W_SIZE> get_mu_theta(ap_uint<MU_W_SIZE> word){
    return word.range(MU_PT_W_SIZE + MU_PHI_W_SIZE + MU_THETA_W_SIZE -1, MU_PT_W_SIZE + MU_PHI_W_SIZE);
}
////////
ap_uint<TRK_PT_W_SIZE> get_trk_pt(ap_uint<TRK_W_SIZE> word){
    return word.range(TRK_PT_W_SIZE-1, 0);
}
ap_uint<TRK_PHI_W_SIZE> get_trk_phi(ap_uint<TRK_W_SIZE> word){
    return word.range(TRK_PT_W_SIZE + TRK_PHI_W_SIZE -1, TRK_PT_W_SIZE);
}
ap_int<TRK_THETA_W_SIZE> get_trk_theta(ap_uint<TRK_W_SIZE> word){
    return word.range(TRK_PT_W_SIZE + TRK_PHI_W_SIZE + TRK_THETA_W_SIZE -1, TRK_PT_W_SIZE + TRK_PHI_W_SIZE);
}
////////
ap_uint<TKMU_PT_W_SIZE> get_tkmu_pt(ap_uint<TKMU_W_SIZE> word){
    return word.range(TKMU_PT_W_SIZE-1, 0);
}
ap_uint<TKMU_PHI_W_SIZE> get_tkmu_phi(ap_uint<TKMU_W_SIZE> word){
    return word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE -1, TKMU_PT_W_SIZE);
}
ap_int<TKMU_THETA_W_SIZE> get_tkmu_theta(ap_uint<TKMU_W_SIZE> word){
    return word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE -1, TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE);
}
