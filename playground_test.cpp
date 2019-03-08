#include <cstdio>
#include "src/playground.h"
#include "ap_int.h"
#include "ap_fixed.h"

#define NTEST 20

int main()
{
    // ap_uint<16> inputs[NTEST];
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> inputs[NTEST];
    for (uint i = 0; i < NTEST; ++i)
        inputs[i] = i;

    // ap_uint<16> outputs[NTEST];
    ap_uint<TRK_W_SIZE*N_TRK + MU_W_SIZE*N_MU> outputs[NTEST];
    for (uint i = 0; i < NTEST; ++i)
        outputs[i] = 0;

    for (uint itest = 0; itest < NTEST; ++itest)
    {
        pipelined_transfer_corrdf(inputs[itest], outputs[itest]);
        printf(">> Clock %i :: in -> %i, out %i\n", itest, inputs[itest].to_int(), outputs[itest].to_int());
    }
}