#include "dataformats.h"

// unpackers
ap_uint<MU_PT_W_SIZE> get_mu_pt(ap_uint<MU_W_SIZE> word){
    return word.range(MU_PT_W_SIZE-1, 0);
}
ap_uint<MU_PHI_W_SIZE> get_mu_phi(ap_uint<MU_W_SIZE> word){
    return word.range(MU_PT_W_SIZE + MU_PHI_W_SIZE -1, MU_PT_W_SIZE);
}
ap_uint<MU_THETA_W_SIZE> get_mu_theta(ap_uint<MU_W_SIZE> word){
    return word.range(MU_PT_W_SIZE + MU_PHI_W_SIZE + MU_THETA_W_SIZE -1, MU_PT_W_SIZE + MU_PHI_W_SIZE);
}
ap_uint<1> get_mu_theta_sign(ap_uint<MU_W_SIZE> word){
    return word.range(MU_PT_W_SIZE + MU_PHI_W_SIZE + MU_THETA_W_SIZE, MU_PT_W_SIZE + MU_PHI_W_SIZE + MU_THETA_W_SIZE);
}
ap_uint<1> get_mu_charge(ap_uint<MU_W_SIZE> word){
    return word.range(MU_PT_W_SIZE + MU_PHI_W_SIZE + MU_THETA_W_SIZE +1, MU_PT_W_SIZE + MU_PHI_W_SIZE + MU_THETA_W_SIZE +1);
}
////////
ap_uint<TRK_PT_W_SIZE> get_trk_pt(ap_uint<TRK_W_SIZE> word){
    return word.range(TRK_PT_W_SIZE-1, 0);
}
ap_uint<TRK_PHI_W_SIZE> get_trk_phi(ap_uint<TRK_W_SIZE> word){
    return word.range(TRK_PT_W_SIZE + TRK_PHI_W_SIZE -1, TRK_PT_W_SIZE);
}
ap_uint<TRK_THETA_W_SIZE> get_trk_theta(ap_uint<TRK_W_SIZE> word){
    return word.range(TRK_PT_W_SIZE + TRK_PHI_W_SIZE + TRK_THETA_W_SIZE -1, TRK_PT_W_SIZE + TRK_PHI_W_SIZE);
}
ap_uint<1> get_trk_theta_sign(ap_uint<TRK_W_SIZE> word){
    return word.range(TRK_PT_W_SIZE + TRK_PHI_W_SIZE + TRK_THETA_W_SIZE, TRK_PT_W_SIZE + TRK_PHI_W_SIZE + TRK_THETA_W_SIZE);
}
ap_uint<1> get_trk_charge(ap_uint<TRK_W_SIZE> word){
    return word.range(TRK_PT_W_SIZE + TRK_PHI_W_SIZE + TRK_THETA_W_SIZE +1, TRK_PT_W_SIZE + TRK_PHI_W_SIZE + TRK_THETA_W_SIZE +1);
}
ap_uint<TRK_CHISQ_W_SIZE> get_trk_chisq(ap_uint<TRK_W_SIZE> word){
    return word.range(TRK_PT_W_SIZE + TRK_PHI_W_SIZE + TRK_THETA_W_SIZE +1 +1 + TRK_CHISQ_W_SIZE -1, TRK_PT_W_SIZE + TRK_PHI_W_SIZE + TRK_THETA_W_SIZE +1 + 1);
}
ap_uint<TRK_NSTUBS_W_SIZE> get_trk_nstubs(ap_uint<TRK_W_SIZE> word){
    return word.range(TRK_PT_W_SIZE + TRK_PHI_W_SIZE + TRK_THETA_W_SIZE +1 +1 + TRK_CHISQ_W_SIZE + TRK_NSTUBS_W_SIZE -1, TRK_PT_W_SIZE + TRK_PHI_W_SIZE + TRK_THETA_W_SIZE +1 +1 + TRK_CHISQ_W_SIZE);
}

////////
ap_uint<TKMU_PT_W_SIZE> get_tkmu_pt(ap_uint<TKMU_W_SIZE> word){
    return word.range(TKMU_PT_W_SIZE-1, 0);
}
ap_uint<TKMU_PHI_W_SIZE> get_tkmu_phi(ap_uint<TKMU_W_SIZE> word){
    return word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE -1, TKMU_PT_W_SIZE);
}
ap_uint<TKMU_THETA_W_SIZE> get_tkmu_theta(ap_uint<TKMU_W_SIZE> word){
    return word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE -1, TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE);
}
ap_uint<1> get_tkmu_theta_sign(ap_uint<TKMU_W_SIZE> word){
    return word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE, TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE);
}
ap_uint<1> get_tkmu_charge(ap_uint<TKMU_W_SIZE> word){
    return word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE +1, TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE +1);
}

///////////////////////////
// packers

ap_uint<TKMU_W_SIZE> build_tkmu_word(tkmu_t tkmu)
{
    ap_uint<TKMU_W_SIZE> word;
    word.range(TKMU_PT_W_SIZE-1, 0) = tkmu.pt;
    word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE -1, TKMU_PT_W_SIZE) = tkmu.phi;
    word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE -1, TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE) = tkmu.theta;
    word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE, TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE) = tkmu.theta_sign;
    word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE +1, TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE +1) = tkmu.charge;
    return word;
}

///////////////////////////
// ref getters

ap_range_ref<TKMU_W_SIZE, false> ref_to_tkmu_pt(ap_uint<TKMU_W_SIZE>& word)
{
    return word.range(TKMU_PT_W_SIZE-1, 0);    
}
////////
ap_range_ref<TKMU_W_SIZE, false> ref_to_tkmu_phi(ap_uint<TKMU_W_SIZE>& word)
{
    return word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE -1, TKMU_PT_W_SIZE);
}
////////
ap_range_ref<TKMU_W_SIZE, false> ref_to_tkmu_theta(ap_uint<TKMU_W_SIZE>& word)
{
    return word.range(TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE + TKMU_THETA_W_SIZE -1, TKMU_PT_W_SIZE + TKMU_PHI_W_SIZE);
}
////////



