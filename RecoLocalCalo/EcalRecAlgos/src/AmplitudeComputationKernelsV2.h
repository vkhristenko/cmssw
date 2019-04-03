#ifndef RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernelsV2
#define RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernelsV2

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

namespace v2 {

void minimization_procedure(
        device_data& d_data, 
        host_data& h_data);

}

}}

#endif // RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernelsV2
