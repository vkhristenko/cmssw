#include "RecoLocalCalo/HcalRecAlgos/interface/MahiGPU.h"

#include <cuda/api_wrappers.h>

namespace hcal { namespace mahi {

// TODO: add/validate restrict (will increase #registers in use by the kernel)
template<int NSAMPLES>
__global__
void kernel_prep1d(
        uint16_t const* dataf01HE,
        uint16_t const* dataf5HB,
        uint32_t const* idsf01HE,
        uint32_t const* idsf5HB,
        uint32_t const nchannels) {
    // constants
    constexpr auto nsamples = NSAMPLES;

    // indices
    auto const tx = threadIdx.x + blockIdx.x*blockDim.x;
    int const linearCh = tx / nsamples;
    auto const sample = threadIdx.x % nsamples;
    auto const nchannels_per_block = blockDim.x / nsamples;

    // remove 
    if (linearCh >= nchannels) return;


}

void entryPoint(
        InputDataGPU const& inputGPU,
        ConditionsProducts const& conditions,
        ConfigParameters const& configParameters,
        cuda::stream_t<>& cudaStream) {
    // 
}

}}
