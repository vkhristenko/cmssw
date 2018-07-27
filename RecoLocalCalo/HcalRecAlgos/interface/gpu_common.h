#ifndef RecoLocalCalo_HcalRecAlgos_interface_gpu_common_h
#define RecoLocalCalo_HcalRecAlgos_interface_gpu_common_h

#include <cuda.h>
#include <cuda_runtime.h>

#include <cassert>
#include <vector>

namespace hcal { namespace cuda {
    void assert_if_error();

    template<typename T>
    inline 
    void copy_vector_to_device(std::vector<T> const& src, T *dest) {
        cudaMemcpy(dest, src.data(), src.size() * sizeof(T), 
            cudaMemcpyHostToDevice);
    }

}}

#endif // RecoLocalCalo_HcalRecAlgos_interface_gpu_common_h
