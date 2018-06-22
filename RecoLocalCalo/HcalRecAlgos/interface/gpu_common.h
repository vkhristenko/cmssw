#ifndef RecoLocalCalo_HcalRecAlgos_interface_gpu_common_h
#define RecoLocalCalo_HcalRecAlgos_interface_gpu_common_h

#include <cuda.h>
#include <cuda_runtime.h>

#include <cassert>

namespace hcal { namespace cuda {
    void assert_if_error();
}}

#endif // RecoLocalCalo_HcalRecAlgos_interface_gpu_common_h
