#ifndef EventFilter_HcalRawToDigi_interface_DecodeGPU_h
#define EventFilter_HcalRawToDigi_interface_DecodeGPU_h

#include <cuda/api_wrappers.h>

#include "EventFilter/HcalRawToDigi/plugins/DeclsForKernels.h"

namespace hcal { namespace raw {

void entryPoint(
        InputDataCPU const&, InputDataGPU&, OutputDataGPU&,
        ScratchDataGPU&, OutputDataCPU&,
        ConditionsProducts const&, ConfigurationParameters const&,
        cuda::stream_t<> &cudaStream,
        uint32_t const, uint32_t const);

}}

#endif // EventFilter_HcalRawToDigi_interface_DecodeGPU_h
