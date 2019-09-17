#ifndef CORRELATOR_H
#define CORRELATOR_H

// #define AP_INT_MAX_W 2048
// #include "dataformats.h"
#include "algo_parts.h"
#include "ap_int.h"
#include "hls_stream.h"

// 2048 bit width in the UF setup -> keep 15 tracks x 100 bits + 12 EMTF x 40 bits = 1980 bits in input
// starting easy with no pipelining of the event information
// to make the pipeline, have the input word contain various tracks but the same EMTFs for several addresses
void correlator (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info[100], ap_uint<TKMU_W_SIZE*N_TKMU> out_info[100],
    ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
    ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
    ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2);

// no loop on bx inside
void correlator_one (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info,
    ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
    ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
    ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2);

// as above, but interfaced with arrays
void arr_correlator_one (
    ap_uint<TRK_W_SIZE> tracks [N_TRK], ap_uint<MU_W_SIZE> muons [N_MU], // inputs
    ap_uint<TKMU_W_SIZE> tkmus [N_TKMU]); // outputs

// handles one track arriving from each sector
void arr_correlator_mult (
    ap_uint<TRK_W_SIZE> all_tracks [N_TRK_SECTORS][N_TRK], ap_uint<MU_W_SIZE> muons [N_MU], // inputs
    ap_uint<TKMU_W_SIZE> tkmus [N_TKMU]); // outputs

//////////////// ----- versions based on the matcher class
// at each event, the input port contains the muons, the tracks, and a flag signallinf ig this is a new bx

// void arr_correlator_mult_matcher (
//     ap_uint<TRK_W_SIZE>  all_tracks [N_TRK_SECTORS], // inputs : one track incoming per sector
//     ap_uint<MU_W_SIZE>   muons [N_MU],   // inputs : all the muons ready at the same time
//     ap_uint<TKMU_W_SIZE> tkmus [N_TKMU], // output TkMus
//     ap_uint<1>           new_bx
// );

// OLD
// muons are exposed in the input registers, and tracks are streamed in 
// assume that new_bx == 1 -> only muons are read
// and when new_bx == 0 -> read track stream

// NEW: correlator has N_MU in input ports and then a stream
// if the stream is empty, will flush out the TkMu and reinitialise the muons
// otherwise it will simply process the stream

void correlator_stream (
    muon_t                in_muons  [N_MU],          // inputs : all the muons ready at the same time
    hls::stream<track_t>  in_tracks [N_TRK_SECTORS], // inputs : one track incoming per sector as a stream
    tkmu_t                out_tkmu  [N_TKMU],        // output TkMus
    ap_uint<1>            &corr_done);        
    // ap_uint<1>            new_bx);                   // the signal to declare a new incoming event

// the interface between the BRAM and the correlator
// the input is a full MTF7 BRAM line that contains either muons or tracks
// the output is the packed version of the tkmus
// in_muons = 1 -> in_info contains muons; in_muons = 0 -> in_info contains tracks
void BRAM_to_corr (ap_uint<MTF7_BRAM_SIZE> in_info, ap_uint<1> in_muons, ap_uint<MTF7_BRAM_SIZE> &out_info, ap_uint<1> &out_valid);

// // handles one track arriving from each sector
// // spies are attached to the first track sector
// void correlator_mult (ap_uint<TRK_W_SIZE*N_TRK*N_TRK_SECTORS + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info,
//     ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
//     ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
//     ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2);


void BRAM_to_corr_nostream (ap_uint<MTF7_BRAM_SIZE> in_info, ap_uint<1> in_muons, ap_uint<MTF7_BRAM_SIZE> &out_info, ap_uint<1> &out_valid);
void correlator_nostream (
    muon_t                in_muons  [N_MU],          // inputs : all the muons ready at the same time
    track_t               in_tracks [N_TRK_SECTORS], // inputs : one track incoming per sector as a stream
    ap_uint<1>            sending_trk,               // sending tracks -> start the correlator
    tkmu_t                out_tkmu  [N_TKMU],        // output TkMus
    ap_uint<1>            &corr_done);        


// just copy N_MU tracks to tkmu output collection - same interface as correlator_one for easy debug
void passthrough (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info,
    ap_uint<TRK_W_SIZE>  &spy_trk1,  ap_uint<TRK_W_SIZE>  &spy_trk2,
    ap_uint<MU_W_SIZE>   &spy_mu1,   ap_uint<MU_W_SIZE>   &spy_mu2,
    ap_uint<TKMU_W_SIZE> &spy_tkmu1, ap_uint<TKMU_W_SIZE> &spy_tkmu2);

// ---------------------------------------------------------------------------

// encapsulate packing/unpacking of my dataformat
void top_arr_correlator_one (ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info);

// encapsulate packing/unpacking of my dataformat
// void top_arr_correlator_mult (ap_uint<TRK_W_SIZE*N_TRK*N_TRK_SECTORS + MU_W_SIZE*N_MU> in_info, ap_uint<TKMU_W_SIZE*N_TKMU> &out_info);


#endif
