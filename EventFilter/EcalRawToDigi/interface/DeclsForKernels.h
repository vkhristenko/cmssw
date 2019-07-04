#ifndef EventFilter_EcalRawToDigi_interface_DeclsForKernels_h
#define EventFilter_EcalRawToDigi_interface_DeclsForKernels_h

#include <vector>

#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "EventFilter/EcalRawToDigi/interface/DCCRawDataDefinitions.h"

namespace ecal { namespace raw {

constexpr auto empty_event_size = EMPTYEVENTSIZE;
constexpr uint32_t nfeds_max = 54;
constexpr uint32_t nbytes_per_fed_max = 10 * 1024;

struct InputDataCPU {
    std::vector<unsigned char, CUDAHostAllocator<unsigned char>> data;
    std::vector<uint32_t, CUDAHostAllocator<uint32_t>> offsets;
    std::vector<int, CUDAHostAllocator<int>> feds;

    void allocate() {
        // 2KB per FED resize
        data.resize(nfeds_max * sizeof(unsigned char) * nbytes_per_fed_max);
        offsets.resize(nfeds_max, 0);
        feds.resize(nfeds_max, 0);
    }
};

struct InputDataGPU {
    unsigned char *data=nullptr;
    uint32_t *offsets=nullptr;
    int *feds=nullptr;

    void allocate() {
        cudaCheck( cudaMalloc((void**)&data, 
            sizeof(unsigned char) * nbytes_per_fed_max * nfeds_max) );
        cudaCheck( cudaMalloc((void**)&offsets,
            sizeof(uint32_t) * nfeds_max) );
        cudaCheck( cudaMalloc((void**)&feds,
            sizeof(int) * nfeds_max) );
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

#endif // EventFilter_EcalRawToDigi_interface_DeclsForKernels_h
