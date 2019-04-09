#include "CondFormats/EcalObjects/interface/EcalPedestals.h"

#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"

class Pedestals {
    float *mean_x12, *mean_x6, *mean_x1;
    float *rms_x12, *rms_x6, *rms_x1;
};

class EcalPedestalsGPU {
public:
    using Product = Pedestals;

    EcalPedestalsGPU(EcalPedestals const&);

public:
    // in the future, we need to arrange so to avoid this copy on the host
    // store eb first then ee
    // maintain offset
    int offsetEB;
    std::vector<float, CUDAHostAllocator<float>> mean_x12;
    std::vector<float, CUDAHostAllocator<float>> rms_x12;
    std::vector<float, CUDAHostAllocator<float>> mean_x6;
    std::vector<float, CUDAHostAllocator<float>> rms_x6;
    std::vector<float, CUDAHostAllocator<float>> mean_x1;
    std::vector<float, CUDAHostAllocator<float>> rms_x1;
};
