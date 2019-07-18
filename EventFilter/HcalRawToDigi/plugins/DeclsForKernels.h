#ifndef EventFilter_HcalRawToDigi_interface_DeclsForKernels_h
#define EventFilter_HcalRawToDigi_interface_DeclsForKernels_h

#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "EventFilter/HcalRawToDigi/plugins/ElectronicsMappingGPU.h"

#include "CUDADataFormats/HcalDigi/interface/DigiCollection.h"

#include <vector>

namespace hcal { namespace raw {

constexpr int32_t empty_event_size = 32;
constexpr uint32_t utca_nfeds_max = 50;
constexpr uint32_t nbytes_per_fed_max = 10 * 1024;

// each collection corresponds to a particular flavor with a certain number of 
// samples per digi
constexpr uint32_t numOutputCollections = 2;
constexpr uint8_t OutputF01HE = 0;
constexpr uint8_t OutputF5HB = 1;

struct ConfigurationParameters {
    uint32_t maxChannelsF01HE;
    uint32_t maxChannelsF5HB;
    uint32_t nsamplesF01HE;
    uint32_t nsamplesF5HB;
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

struct OutputDataCPU {
    std::vector<uint32_t, CUDAHostAllocator<uint32_t>> nchannels;

    void allocate() {
        nchannels.resize(numOutputCollections);
    }
};

struct ScratchDataGPU {
    // depends on tHE number of output collections
    // that is a statically known predefined number!!!
    uint32_t *pChannelsCounters = nullptr;

    void allocate(ConfigurationParameters const&) {
        cudaCheck( cudaMalloc((void**)&pChannelsCounters,
            sizeof(uint32_t) * numOutputCollections) );
    }

    void deallocate(ConfigurationParameters const&) {
        if (pChannelsCounters) {
            cudaCheck( cudaFree(pChannelsCounters) );
        }
    }
};

struct OutputDataGPU {
    // qie 11 HE
    uint16_t *digisF01HE = nullptr;
    uint32_t *idsF01HE = nullptr;

    // qie 8 HB
    uint16_t *digisF5HB = nullptr;
    uint32_t *idsF5HB = nullptr;

    void allocate(ConfigurationParameters const& config) {
        cudaCheck( cudaMalloc((void**)&digisF01HE,
            config.maxChannelsF01HE * sizeof(uint16_t) *
            config.nsamplesF01HE * Flavor01::WORDS_PER_SAMPLE) );
        cudaCheck( cudaMalloc((void**)&idsF01HE,
            sizeof(uint32_t) * config.maxChannelsF01HE) );
        
        cudaCheck( cudaMalloc((void**)&digisF5HB,
            config.maxChannelsF5HB * sizeof(uint16_t) *
            config.nsamplesF5HB * Flavor5::WORDS_PER_SAMPLE) );
        cudaCheck( cudaMalloc((void**)&idsF5HB,
            sizeof(uint32_t) * config.maxChannelsF5HB) );
    }

    void deallocate(ConfigurationParameters const& config) {
        if (digisF01HE) {
            cudaCheck( cudaFree(digisF01HE) );
            cudaCheck( cudaFree(idsF01HE) );

            cudaCheck( cudaFree(digisF5HB) );
            cudaCheck( cudaFree(idsF5HB) );
        }    
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

struct ConditionsProducts {
    ElectronicsMappingGPU::Product const& eMappingProduct;
};

}}

#endif // EventFilter_HcalRawToDigi_interface_DeclsForKernels_h
