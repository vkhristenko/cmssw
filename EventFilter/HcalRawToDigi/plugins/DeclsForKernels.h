#ifndef EventFilter_HcalRawToDigi_interface_DeclsForKernels_h
#define EventFilter_HcalRawToDigi_interface_DeclsForKernels_h

#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include <vector>

namespace hcal { namespace raw {

constexpr int32_t empty_event_size = 32;
constexpr uint32_t utca_nfeds_max = 50;
constexpr uint32_t nbytes_per_fed_max = 10 * 1024;

struct ConfigurationParameters {
    uint32_t maxChannels;
};

struct InputDataCPU {
    std::vector<unsigned char, CUDAHostAllocator<unsigned char>> data;
    std::vector<uint32_t, CUDAHostAllocator<uint32_t>> offsets;
    std::vector<int, CUDAHostAllocator<int>> feds;

    void allocate() {
        data.resize(utca_nfeds_max * sizeof(unsigned char) * nbytes_per_fed_max);
        offsets.resize(utca_nfeds_max, 0);
        feds.resize(utca_nfeds_max, 0);
    }
};

struct InputDataGPU {
    unsigned char *data = nullptr;
    uint32_t *offsets = nullptr;
    int *feds = nullptr;

    void allocate() {
        cudaCheck( cudaMalloc((void**)&data,
            sizeof(unsigned char) * nbytes_per_fed_max * utca_nfeds_max) );
        cudaCheck( cudaMalloc((void**)&offsets,
            sizeof(uint32_t) * utca_nfeds_max) );
        cudaCheck( cudaMalloc((void**)&feds,
            sizeof(int) * utca_nfeds_max) );
    }

    void deallocate() {
        if (data) {
            cudaCheck( cudaFree(data) );
            cudaCheck( cudaFree(offsets) );
            cudaCheck( cudaFree(feds) );
        }
    }
};

}}

#endif // EventFilter_HcalRawToDigi_interface_DeclsForKernels_h
