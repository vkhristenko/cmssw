#ifndef RecoLocalCalo_HcalRecAlgos_interface_MahiGPU_h
#define RecoLocalCalo_HcalRecAlgos_interface_MahiGPU_h

#include "RecoLocalCalo/HcalRecAlgos/interface/DeclsForKernels.h"

#include <cuda/api_wrappers.h>

namespace hcal { namespace mahi {

void entryPoint(
        InputDataGPU const&, OutputDataGPU&,
        ConditionsProducts const&, ScratchDataGPU&,
        ConfigParameters const&, cuda::stream_t<>&);

}}

#endif // RecoLocalCalo_HcalRecAlgos_interface_MahiGPU_h
