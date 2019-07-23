#ifndef RecoLocalCalo_HcalRecProducers_src_HcalCombinedRecordsGPU_h
#define RecoLocalCalo_HcalRecProducers_src_HcalCombinedRecordsGPU_h

#include "CondFormats/DataRecord/interface/HcalPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/HcalQIEDataRcd.h"
#include "CondFormats/DataRecord/interface/HcalQIETypesRcd.h"

#include "RecoLocalCalo/HcalRecProducers/src/HcalCombinedRecord.h"

using HcalConvertedPedestalsRcd = HcalCombinedRecord<
    HcalPedestalsRcd,
    HcalQIEDataRcd,
    HcalQIETypesRcd
>;

#endif // RecoLocalCalo_HcalRecProducers_src_HcalCombinedRecordsGPU_h
