#include "EcalESProducerGPU.h"

#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPulseShapesRcd.h"
#include "CondFormats/DataRecord/interface/EcalPulseCovariancesRcd.h"
#include "CondFormats/DataRecord/interface/EcalSamplesCorrelationRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeBiasCorrectionsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"

#include "EcalPedestalsGPU.h"
#include "EcalGainRatiosGPU.h"
#include "EcalPulseShapesGPU.h"
#include "EcalPulseCovariancesGPU.h"
#include "EcalSamplesCorrelationGPU.h"
#include "EcalTimeBiasCorrectionsGPU.h"
#include "EcalTimeCalibConstantsGPU.h"

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

DEFINE_FWK_EVENTSETUP_MODULE(EcalPedestalsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalGainRatiosGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalPulseShapesGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalPulseCovariancesGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalSamplesCorrelationGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalTimeBiasCorrectionsGPUESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(EcalTimeCalibConstantsGPUESProducer);
