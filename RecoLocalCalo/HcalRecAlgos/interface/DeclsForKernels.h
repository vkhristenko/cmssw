#ifndef RecoLocalCalo_HcalRecAlgos_interface_DeclsForKernels_h
#define RecoLocalCalo_HcalRecAlgos_interface_DeclsForKernels_h

#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "CUDADataFormats/HcalDigi/interface/DigiCollection.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalRecoParamsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalGainWidthsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalGainsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalLUTCorrsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalPedestalWidthsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalPedestalsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalQIECodersGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalRespCorrsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalTimeCorrsGPU.h"

namespace hcal { namespace mahi {

struct ConditionsProducts {
    HcalGainWidthsGPU::Product const& gainWidths;
    HcalGainsGPU::Product const& gains;
    HcalLUTCorrsGPU::Product const& lutCorrs;
    HcalPedestalWidthsGPU::Product const& pedestalWidths;
    HcalPedestalsGPU::Product const& pedestals;
    HcalQIECodersGPU::Product const& qieCoders;
    HcalRecoParamsGPU::Product const& recoParams;
    HcalRespCorrsGPU::Product const& respCorrs;
    HcalTimeCorrsGPU::Product const& timeCorrs;
};

struct ConfigParameters {
    uint32_t maxChannels;
};

struct InputDataGPU {
    DigiCollection<Flavor01> const& f01HEDigis;
    DigiCollection<Flavor5> const& f5HBDigis;
};

}}

#endif // RecoLocalCalo_HcalRecAlgos_interface_DeclsForKernels_h
