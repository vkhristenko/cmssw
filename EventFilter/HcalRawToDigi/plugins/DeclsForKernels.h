#ifndef EventFilter_HcalRawToDigi_interface_DeclsForKernels_h
#define EventFilter_HcalRawToDigi_interface_DeclsForKernels_h

#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "EventFilter/HcalRawToDigi/plugins/ElectronicsMappingGPU.h"

#include <vector>

namespace hcal { namespace raw {

constexpr int32_t empty_event_size = 32;
constexpr uint32_t utca_nfeds_max = 50;
constexpr uint32_t nbytes_per_fed_max = 10 * 1024;

struct Flavor01 {
    static constexpr int WORDS_PER_SAMPLES = 1;
    static constexpr int HEADER_WORDS = 1;
};

struct Flavor2 {
    static constexpr int WORDS_PER_SAMPLES = 2;
    static constexpr int HEADER_WORDS = 1;
};

struct Flavor3 {
    static constexpr int WORDS_PER_SAMPLES = 1;
    static constexpr int HEADER_WORDS = 1;
};

struct Flavor4 {
    static constexpr int WORDS_PER_SAMPLES = 1;
    static constexpr int HEADER_WORDS = 1;
};

struct Flavor5 {
    static constexpr float WORDS_PER_SAMPLES = 0.5;
    static constexpr int HEADER_WORDS = 1;
};

struct ConfigurationParameters {
    uint32_t maxChannelsQIE11HE;
    uint32_t maxChannelsQIE8HB;
    uint32_t nsamples_qie10_hf;
    uint32_t nsamples_qie11_he;
    uint32_t nsamples_qie8_hb;
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

struct OutputDataGPU {
    // qie 11 he
    uint16_t *digisQIE11HE = nullptr;
    uint32_t *idsQIE11HE = nullptr;

    // qie 8 hb
    uint16_t *digisQIE8HB = nullptr;
    uint32_t *idsQIE8HB = nullptr;

    void allocate(ConfigurationParameters const& config) {
        cudaCheck( cudaMalloc((void**)&digisQIE11HE,
            config.maxChannelsQIE11HE * sizeof(uint16_t) *
            config.nsamples_qie11_he * Flavor01::WORDS_PER_SAMPLES) );
        cudaCheck( cudaMalloc((void**)&idsQIE11HE,
            sizeof(uint32_t) * config.maxChannelsQIE11HE) );
        
        cudaCheck( cudaMalloc((void**)&digisQIE8HB,
            config.maxChannelsQIE8HB * sizeof(uint16_t) *
            config.nsamples_qie8_hb * Flavor5::WORDS_PER_SAMPLES) );
        cudaCheck( cudaMalloc((void**)&idsQIE11HE,
            sizeof(uint32_t) * config.maxChannelsQIE11HE) );
    }

    void deallocate(ConfigurationParameters const& config) {
        if (digisQIE11HE) {
            cudaCheck( cudaFree(digisQIE11HE) );
            cudaCheck( cudaFree(idsQIE11HE) );

            cudaCheck( cudaFree(digisQIE8HB) );
            cudaCheck( cudaFree(idsQIE8HB) );
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
