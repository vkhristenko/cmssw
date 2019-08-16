#include "EcalESProducerGPU.h"

// for uncalibrechit
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPulseShapesRcd.h"
#include "CondFormats/DataRecord/interface/EcalPulseCovariancesRcd.h"
#include "CondFormats/DataRecord/interface/EcalSamplesCorrelationRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeBiasCorrectionsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"
// for rechit
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"


// for uncalibrechit
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalPedestalsGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalGainRatiosGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalPulseShapesGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalPulseCovariancesGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSamplesCorrelationGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalTimeBiasCorrectionsGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalTimeCalibConstantsGPU.h"
// for rechit
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalADCToGeVConstantGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalIntercalibConstantsGPU.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalChannelStatusGPU.h"


#include <iostream>

using EcalPedestalsGPUESProducer = EcalESProducerGPU<EcalPedestalsGPU,
                                                     EcalPedestals,
                                                     EcalPedestalsRcd>;
using EcalGainRatiosGPUESProducer = EcalESProducerGPU<EcalGainRatiosGPU,
                                                      EcalGainRatios,
                                                      EcalGainRatiosRcd>;
using EcalPulseShapesGPUESProducer = EcalESProducerGPU<EcalPulseShapesGPU,
                                                       EcalPulseShapes,
                                                       EcalPulseShapesRcd>;
using EcalPulseCovariancesGPUESProducer = EcalESProducerGPU<EcalPulseCovariancesGPU,
                                                            EcalPulseCovariances,
                                                            EcalPulseCovariancesRcd>;
using EcalSamplesCorrelationGPUESProducer = EcalESProducerGPU<
    EcalSamplesCorrelationGPU,
    EcalSamplesCorrelation,
    EcalSamplesCorrelationRcd>;

using EcalTimeBiasCorrectionsGPUESProducer = EcalESProducerGPU<
    EcalTimeBiasCorrectionsGPU,
    EcalTimeBiasCorrections,
    EcalTimeBiasCorrectionsRcd>;

using EcalTimeCalibConstantsGPUESProducer = EcalESProducerGPU<
    EcalTimeCalibConstantsGPU,
    EcalTimeCalibConstants,
    EcalTimeCalibConstantsRcd>;

using EcalADCToGeVConstantGPUESProducer = EcalESProducerGPU<EcalADCToGeVConstantGPU,
    EcalADCToGeVConstant,
    EcalADCToGeVConstantRcd>;

using EcalIntercalibConstantsGPUESProducer = EcalESProducerGPU<EcalIntercalibConstantsGPU,
    EcalIntercalibConstants,
    EcalIntercalibConstantsRcd>;

using EcalChannelStatusGPUESProducer = EcalESProducerGPU<EcalChannelStatusGPU,
    EcalChannelStatus,
    EcalChannelStatusRcd>;
    
DEFINE_FWK_EVENTSETUP_MODULE(EcalPedestalsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalGainRatiosGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalPulseShapesGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalPulseCovariancesGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalSamplesCorrelationGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalTimeBiasCorrectionsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalTimeCalibConstantsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalADCToGeVConstantGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalIntercalibConstantsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalChannelStatusGPUESProducer);



