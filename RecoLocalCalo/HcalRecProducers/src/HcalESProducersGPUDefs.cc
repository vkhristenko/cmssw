#include "HcalESProducerGPU.h"

#include "CondFormats/HcalObjects/interface/HcalRecoParams.h"
#include "CondFormats/HcalObjects/interface/HcalPedestals.h"
#include "CondFormats/HcalObjects/interface/HcalGains.h"
#include "CondFormats/HcalObjects/interface/HcalLUTCorrs.h"
#include "CondFormats/HcalObjects/interface/HcalRespCorrs.h"
#include "CondFormats/HcalObjects/interface/HcalTimeCorrs.h"
#include "CondFormats/HcalObjects/interface/HcalPedestalWidths.h"
#include "CondFormats/HcalObjects/interface/HcalGainWidths.h"
#include "CondFormats/HcalObjects/interface/HcalQIEData.h"
#include "CondFormats/HcalObjects/interface/HcalQIETypes.h"
#include "CondFormats/HcalObjects/interface/HcalSiPMParameters.h"
#include "CondFormats/HcalObjects/interface/HcalSiPMCharacteristics.h"

#include "CondFormats/DataRecord/interface/HcalRecoParamsRcd.h"
#include "CondFormats/DataRecord/interface/HcalPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/HcalGainsRcd.h"
#include "CondFormats/DataRecord/interface/HcalLUTCorrsRcd.h"
#include "CondFormats/DataRecord/interface/HcalRespCorrsRcd.h"
#include "CondFormats/DataRecord/interface/HcalTimeCorrsRcd.h"
#include "CondFormats/DataRecord/interface/HcalPedestalWidthsRcd.h"
#include "CondFormats/DataRecord/interface/HcalGainWidthsRcd.h"
#include "CondFormats/DataRecord/interface/HcalQIEDataRcd.h"
#include "CondFormats/DataRecord/interface/HcalQIETypesRcd.h"
#include "CondFormats/DataRecord/interface/HcalSiPMParametersRcd.h"
#include "CondFormats/DataRecord/interface/HcalSiPMCharacteristicsRcd.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalRecoParamsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalPedestalsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalGainsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalLUTCorrsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalRespCorrsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalTimeCorrsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalPedestalWidthsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalGainWidthsGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalQIECodersGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalQIETypesGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSiPMParametersGPU.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSiPMCharacteristicsGPU.h"

#include <iostream>

using HcalRecoParamsGPUESProducer = HcalESProducerGPU<HcalRecoParamsGPU,
                                                      HcalRecoParams,
                                                      HcalRecoParamsRcd>;
using HcalPedestalsGPUESProducer = HcalESProducerGPU<HcalPedestalsGPU,
                                                     HcalPedestals,
                                                     HcalPedestalsRcd>;

using HcalGainsGPUESProducer = HcalESProducerGPU<HcalGainsGPU,
                                                 HcalGains,
                                                 HcalGainsRcd>;

using HcalLUTCorrsGPUESProducer = HcalESProducerGPU<HcalLUTCorrsGPU,
                                                     HcalLUTCorrs,
                                                     HcalLUTCorrsRcd>;

using HcalRespCorrsGPUESProducer = HcalESProducerGPU<HcalRespCorrsGPU,
                                                     HcalRespCorrs,
                                                     HcalRespCorrsRcd>;

using HcalTimeCorrsGPUESProducer = HcalESProducerGPU<HcalTimeCorrsGPU,
                                                     HcalTimeCorrs,
                                                     HcalTimeCorrsRcd>;

using HcalPedestalWidthsGPUESProducer = HcalESProducerGPU<HcalPedestalWidthsGPU,
                                                     HcalPedestalWidths,
                                                     HcalPedestalWidthsRcd>;

using HcalGainWidthsGPUESProducer = HcalESProducerGPU<HcalGainWidthsGPU,
                                                     HcalGainWidths,
                                                     HcalGainWidthsRcd>;

using HcalQIECodersGPUESProducer = HcalESProducerGPU<HcalQIECodersGPU,
                                                     HcalQIEData,
                                                     HcalQIEDataRcd>;

using HcalQIETypesGPUESProducer = HcalESProducerGPU<
    HcalQIETypesGPU,
    HcalQIETypes,
    HcalQIETypesRcd>;

using HcalSiPMParametersGPUESProducer = HcalESProducerGPU<
    HcalSiPMParametersGPU,
    HcalSiPMParameters,
    HcalSiPMParametersRcd>;

using HcalSiPMCharacteristicsGPUESProducer = HcalESProducerGPU<
    HcalSiPMCharacteristicsGPU,
    HcalSiPMCharacteristics,
    HcalSiPMCharacteristicsRcd>;

DEFINE_FWK_EVENTSETUP_MODULE(HcalRecoParamsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(HcalPedestalsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(HcalGainsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(HcalLUTCorrsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(HcalRespCorrsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(HcalTimeCorrsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(HcalPedestalWidthsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(HcalGainWidthsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(HcalQIECodersGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(HcalQIETypesGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(HcalSiPMParametersGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(HcalSiPMCharacteristicsGPUESProducer);
