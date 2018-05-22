#ifndef RecoLocalCalo_HcalRecProducers_interface_test_kernels_h
#define RecoLocalCalo_HcalRecProducers_interface_test_kernels_h

#include <cuda_runtime.h>
#include <functional>

namespace hcal {
namespace test {
    void vector_add(int*, int*, int*, int const, cudaStream_t);

    template<typename T, typename F>
    void init_host_vector(T* xs, int const N, F func) {
        for (auto i=0; i<N; i++) {
            xs[i] = func(i);
        }
    }
}
}

#endif
