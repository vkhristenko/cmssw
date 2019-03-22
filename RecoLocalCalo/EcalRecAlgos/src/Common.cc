#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

#include <cuda.h>
#include <cuda_runtime.h>

#include <iostream>
#include <cassert>

namespace ecal { namespace cuda {

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
