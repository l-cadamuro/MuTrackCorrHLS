#include "dataformats.h"

class matcher {
    public:
        void init(muon_t mu);
        void process(track_t trk,
            ap_uint<MW_LUT_ANGLE_W_SIZE> phi_low_bound,
            ap_uint<MW_LUT_ANGLE_W_SIZE> phi_high_bound,
            ap_uint<MW_LUT_ANGLE_W_SIZE> theta_low_bound,
            ap_uint<MW_LUT_ANGLE_W_SIZE> theta_high_bound);

        // variables used inside the matcher
        ap_uint <TRK_PT_W_SIZE> sort_par_;
        muon_t mu_;
        tkmu_t tkmu_;
};
