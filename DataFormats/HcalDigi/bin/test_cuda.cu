#include <iostream>

#include "cuda_runtime.h"
#include "cuda.h"

#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"

class Test {
public:
    int v;

    HOST DEVICE Test() {v=0;}
};

__global__ void test() {
    auto x = Test();
    auto f = HBHEDataFrame();
//    auto f = HBHEDataFrame();
}

int main(int argc, char** argv) {
    std::cout << "hello world" << std::endl;

#ifdef __CUDACC__
    std::cout << "cudacc is set " << std::endl;
#else
    std::cout << "cudacc is not set" << std::endl;
#endif

    int nDevices;
    cudaGetDeviceCount(&nDevices);
    std::cout << "nDevices = " << nDevices << std::endl;

    return 0;
}
