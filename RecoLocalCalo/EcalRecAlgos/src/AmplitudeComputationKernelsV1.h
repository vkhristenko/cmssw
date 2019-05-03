#ifndef RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernelsV1
#define RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernelsV1

#include <iostream>
#include <limits>

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

#include "DataFormats/Math/interface/approx_exp.h"
#include "DataFormats/Math/interface/approx_log.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/DeclsForKernels.h"

#include "cuda.h"

//#define DEBUG

//#define ECAL_RECO_CUDA_DEBUG

namespace ecal { namespace multifit {

void minimization_procedure(
        EventInputDataCPU const& eventInputCPU, EventInputDataGPU& eventInputGPU,
        EventOutputDataGPU& eventOutputGPU, EventDataForScratchGPU& scratch,
        ConditionsProducts const& conditions,
        ConfigurationParameters const& configParameters,
        cuda::stream_t<>& cudaStream);

}}

#endif // RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernelsV2
