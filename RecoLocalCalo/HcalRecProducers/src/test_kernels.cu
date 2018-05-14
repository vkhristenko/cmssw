#include "RecoLocalCalo/HcalRecProducers/interface/test_kernels.h"

namespace hcal {
namespace test {
    __global__ void cu_vector_add(int* a, int *b, int *c) {
        int id = blockIdx.x;
        c[id] = a[id] + b[id];
    }

    void vector_add(int* a, int* b, int* c, int const n, cudaStream_t stream) {
        cu_vector_add<<<n, 1, 0, stream>>>(a, b, c);
    }
}
}
