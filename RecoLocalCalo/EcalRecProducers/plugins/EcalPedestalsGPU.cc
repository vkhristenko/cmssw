#include "EcalPedestalsGPU.h"

#include "FWCore/Utilities/interface/typelookup.h"

EcalPedestalsGPU::EcalPedestalsGPU(EcalPedestals const& pedestals) 
    : mean_x12(pedestals.size())
    , rms_x12(pedestals.size())
    , mean_x6(pedestals.size())
    , rms_x6(pedestals.size())
    , mean_x1(pedestals.size())
    , rms_x1(pedestals.size())
{
    // fill in eb
    auto const& barrelValues = pedestals.barrelItems();
    for (unsigned int i=0; i<barrelValues.size(); i++) {
        mean_x12[i] = barrelValues[i].mean_x12;
        rms_x12[i] = barrelValues[i].rms_x12;
        mean_x6[i] = barrelValues[i].mean_x6;
        rms_x6[i] = barrelValues[i].rms_x6;
        mean_x1[i] = barrelValues[i].mean_x1;
        rms_x1[i] = barrelValues[i].rms_x1;
    }
    
    // fill in ee
    auto const& endcapValues = pedestals.endcapItems();
    for (unsigned int i=0; i<endcapValues.size(); i++) {
        mean_x12[i] = endcapValues[i].mean_x12;
        rms_x12[i] = endcapValues[i].rms_x12;
        mean_x6[i] = endcapValues[i].mean_x6;
        rms_x6[i] = endcapValues[i].rms_x6;
        mean_x1[i] = endcapValues[i].mean_x1;
        rms_x1[i] = endcapValues[i].rms_x1;
    }
}


TYPELOOKUP_DATA_REG(EcalPedestalsGPU);
