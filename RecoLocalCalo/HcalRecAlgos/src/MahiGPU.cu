#include "RecoLocalCalo/HcalRecAlgos/interface/MahiGPU.h"

#include <cuda/api_wrappers.h>

namespace hcal { namespace mahi {

void entryPoint(
        InputDataGPU const& inputGPU,
        ConditionsProducts const& conditions,
        ConfigParameters const& configParameters,
        cuda::stream_t<>& cudaStream) {
    // 
}

}}
