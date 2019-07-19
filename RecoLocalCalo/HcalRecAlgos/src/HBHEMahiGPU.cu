#include "RecoLocalCalo/HcalRecAlgos/interface/HBHEMahiGPU.h"

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
