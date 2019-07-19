#include "matcher.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <int N>
inline bool is_within(ap_uint<N> val, ap_uint<N> b_low, ap_uint<N> b_high){
    // printf("xxxxxxx in boundaries debug, val = %10u, blow = %10u, bup = %10u\n", val.to_uint64(), b_low.to_uint64(), b_high.to_uint64());

    if (b_low <= val && val < b_high){
        // printf("   >>>> ritorno True\n");
        return true;
    }
    else{
        // printf("   >>>> ritorno False\n");
        return false;
    }
}

inline bool is_within_th(ap_uint<THETA_W_SIZE> val, ap_uint<THETA_W_SIZE> b_low, ap_uint<THETA_W_SIZE> b_high){
    return is_within<THETA_W_SIZE>(val, b_low, b_high);
}

inline bool is_within_ph(ap_uint<PHI_W_SIZE>   val,  ap_uint<PHI_W_SIZE>  b_low, ap_uint<PHI_W_SIZE>   b_high){
    return is_within<PHI_W_SIZE>(val, b_low, b_high);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void matcher::init(muon_t mu)
{
    mu_ = mu;        // attach the muon to this sorter
    sort_par_ = 0x0; // reset to default value of sorting parameter
    clear(tkmu_);     // init to an invalid tkmu candidate
}

void matcher::process(track_t trk,
            ap_uint<MW_LUT_ANGLE_W_SIZE> phi_low_bound,
            ap_uint<MW_LUT_ANGLE_W_SIZE> phi_high_bound,
            ap_uint<MW_LUT_ANGLE_W_SIZE> theta_low_bound,
            ap_uint<MW_LUT_ANGLE_W_SIZE> theta_high_bound)
{
    // // ---- dtheta
    // ap_uint<TRK_THETA_W_SIZE> trk_theta = get_trk_theta(track);
    // ap_uint<MU_THETA_W_SIZE>  mu_theta = get_mu_theta(muon);

    // // FIXME: which precision? All the same? Convert to common one?
    // // would be way better to move to ap_fixed to get this handled automatically
    ap_uint<MU_THETA_W_SIZE> delta_theta = (trk.theta > mu_.theta) ? (trk.theta - mu_.theta) : (mu_.theta - trk.theta);
    // if (delta_theta < 0)
    //     delta_theta *= -1;

    // ---- dphi
    // ap_uint<TRK_PHI_W_SIZE> trk_phi = get_trk_phi(track);
    // ap_uint<MU_PHI_W_SIZE> mu_phi = get_mu_phi(muon);

    // FIXMEL which precision? All the same? Convert to common one?
    ap_uint<MU_PHI_W_SIZE> delta_phi = (trk.phi > mu_.phi) ? (trk.phi - mu_.phi) : (mu_.phi - trk.phi);
    // currently phi is in units of 0.002, meaning that pi = 1571
    // so subtract this offset to normalize in case
    if (delta_phi > 1571)
        delta_phi -= 1571;

    // if (delta_phi < 0)
    //     delta_phi *= -1;

    bool valid_trk   = (trk.pt > 0);
    bool valid_mu    = (mu_.pt > 0);
    bool trk_nstubs  = (trk.nstubs > 3);
    bool trk_chi2    = (trk.chisq < 100);;
    bool same_endcap = (trk.theta_sign == mu_.theta_sign);
    bool theta_bound = is_within_th(delta_theta, theta_low_bound, theta_high_bound);
    bool phi_bound   = is_within_ph(delta_phi, phi_low_bound, phi_high_bound);
    bool dphi_sign   = (trk.charge > 0 ? (trk.phi > mu_.phi) : (mu_.phi > trk.phi) ); //// FIXME!!! does not account for diffs across the 0/2pi borders

    if (
            valid_trk   &&
            valid_mu    &&
            trk_nstubs  &&
            trk_chi2    &&
            same_endcap &&
            theta_bound &&
            phi_bound   &&
            dphi_sign   &&
            true
        )
    {
        ap_uint<TRK_PT_W_SIZE> sort_par = trk.pt;
        
        // #ifndef __SYNTHESIS__
        // printf("It's a match! Muon of pt = %10u with track of pt = %10u\n", get_mu_pt(muon).to_uint64(), sort_par.to_uint64());
        // #endif
        
        if (sort_par > sort_par_)
        {
            // #ifndef __SYNTHESIS__
            // printf("And it's better than the previous one! Muon of pt = %10u with track of pt = %10u\n", get_mu_pt(muon).to_uint64(), sort_par.to_uint64());
            // #endif
            sort_par_ = sort_par;
            // FIXME: here the correct assignment of qualities
            tkmu_.pt         = trk.pt;
            tkmu_.theta      = trk.theta;
            tkmu_.theta_sign = trk.theta_sign;
            tkmu_.charge     = trk.charge;
            tkmu_.phi        = trk.phi;
        }
        // else
        // {
        //     // #ifndef __SYNTHESIS__
        //     // printf("but NOT better than the previous one! Muon of pt = %10u with track of pt = %10u\n", get_mu_pt(muon).to_uint64(), sort_par.to_uint64());
        //     // #endif
        //     out_sorting_par = in_sorting_par;
        //     tkmu = in_tkmu;
        // }
    }

}
