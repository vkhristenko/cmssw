#ifndef RecoLocalCalo_HcalRecAlgos_interface_DeclsForKernels_h
#define RecoLocalCalo_HcalRecAlgos_interface_DeclsForKernels_h

#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "CUDADataFormats/HcalDigi/interface/DigiCollection.h"

#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalRecoParamsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalGainWidthsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalGainsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalLUTCorrsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalConvertedPedestalWidthsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalConvertedPedestalsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalQIECodersGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalRespCorrsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalTimeCorrsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalQIETypesGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSiPMParametersGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSiPMCharacteristicsGPU.h"

namespace hcal { namespace mahi {

struct ConditionsProducts {
    HcalGainWidthsGPU::Product const& gainWidths;
    HcalGainsGPU::Product const& gains;
    HcalLUTCorrsGPU::Product const& lutCorrs;
    HcalConvertedPedestalWidthsGPU::Product const& pedestalWidths;
    HcalConvertedPedestalsGPU::Product const& pedestals;
    HcalQIECodersGPU::Product const& qieCoders;
    HcalRecoParamsGPU::Product const& recoParams;
    HcalRespCorrsGPU::Product const& respCorrs;
    HcalTimeCorrsGPU::Product const& timeCorrs;
    HcalQIETypesGPU::Product const& qieTypes;
    HcalSiPMParametersGPU::Product const& sipmParameters;
    HcalSiPMCharacteristicsGPU::Product const& sipmCharacteristics;
    HcalTopology const* topology;
    HcalDDDRecConstants const* recConstants;
    uint32_t offsetForHashes;
};

struct ConfigParameters {
    uint32_t maxChannels;
    uint32_t kprep1dChannelsPerBlock;
    int sipmQTSShift;
    int sipmQNTStoSum;
    int firstSampleShift;
};

struct InputDataGPU {
    DigiCollection<Flavor01> const& f01HEDigis;
    DigiCollection<Flavor5> const& f5HBDigis;
};

}}

#endif // RecoLocalCalo_HcalRecAlgos_interface_DeclsForKernels_h
