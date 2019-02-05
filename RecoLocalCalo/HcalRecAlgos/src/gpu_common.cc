#include "RecoLocalCalo/HcalRecAlgos/interface/gpu_common.h"

#include <iostream>
#include <vector>

namespace hcal { namespace cuda {
    void assert_if_error() {
        auto check = [](auto code) {
            if (code != cudaSuccess) {
                std::cout << cudaGetErrorString(code) << std::endl;
                assert(false);
            }
        };

        check(cudaGetLastError());
    }
}}
