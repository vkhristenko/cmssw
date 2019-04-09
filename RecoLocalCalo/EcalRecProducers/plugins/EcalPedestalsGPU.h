#include "CondFormats/EcalObjects/interface/EcalPedestals.h"

#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"

class Pedestals {};

class EcalPedestalsGPU {
public:
    using Product = Pedestals;

    EcalPedestalsGPU(EcalPedestals const&) {}

private:
    std::vector<float, CUDAHostAllocator<float>> mean_x12;
    std::vector<float, CUDAHostAllocator<float>> rms_x12;
    std::vector<float, CUDAHostAllocator<float>> mean_x6;
    std::vector<float, CUDAHostAllocator<float>> rms_x6;
    std::vector<float, CUDAHostAllocator<float>> mean_x1;
    std::vector<float, CUDAHostAllocator<float>> rms_x1;
};
