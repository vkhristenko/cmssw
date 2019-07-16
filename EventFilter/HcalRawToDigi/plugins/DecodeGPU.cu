#include "EventFilter/HcalRawToDigi/plugins/DecodeGPU.h"

namespace hcal { namespace raw {

__global__
void kernel_rawdecode_test(
        unsigned char const* data,
        uint32_t const* offsets,
        int const* feds,
        uint32_t const nBytesTotal) {
    auto const ifed = blockIdx.x;
    auto const fed = feds[ifed];
    auto const offset = offsets[ifed];
    auto const size = ifed==gridDim.x-1 
        ? nBytesTotal - offset 
        : offsets[ifed+1] - offset;

#ifdef HCAL_RAWDECODE_GPUDEBUG
    printf("ifed = %d fed = %d offset = %u size = %u\n", ifed, fed, offset, size);
#endif

    // offset to the right raw buffer
    uint64_t const* buffer = reinterpret_cast<uint64_t const*>(data + offset);

    // 
    // fed header
    //
    auto const fed_header = buffer[0];
    uint32_t const fed_id = (fed_header >> 8) & 0xfff;
    uint32_t const bx = (fed_header >> 20) & 0xfff;
    uint32_t const lv1 = (fed_header >> 32) & 0xffffff;
    uint8_t const trigger_type = (fed_header >> 56) & 0xf;
    uint8_t const bid_fed_header = (fed_header >> 60) & 0xf;

#ifdef HCAL_RAWDECODE_GPUDEBUG
    printf("fed = %d fed_id = %u bx = %u lv1 = %u trigger_type = %u bid = %u\n",
        fed, fed_id, bx, lv1, trigger_type, bid_fed_header);
#endif
}

void entryPoint(
        InputDataCPU const& inputCPU, InputDataGPU& inputGPU,
        cuda::stream_t<> &cudaStream,
        uint32_t const nfedsWithData,
        uint32_t const nbytesTotal) {
    // transfer
    cudaCheck( cudaMemcpyAsync(inputGPU.data,
        inputCPU.data.data(),
        nbytesTotal * sizeof(unsigned char),
        cudaMemcpyHostToDevice,
        cudaStream.id()) );
    cudaCheck( cudaMemcpyAsync(inputGPU.offsets,
        inputCPU.offsets.data(),
        nfedsWithData * sizeof(uint32_t),
        cudaMemcpyHostToDevice,
        cudaStream.id()) );
    cudaCheck( cudaMemcpyAsync(inputGPU.feds,
        inputCPU.feds.data(),
        nfedsWithData * sizeof(int),
        cudaMemcpyHostToDevice,
        cudaStream.id()) );

    kernel_rawdecode_test<<<nfedsWithData, 1, 0, cudaStream.id()>>>(
        inputGPU.data,
        inputGPU.offsets,
        inputGPU.feds,
        nbytesTotal);
    cudaCheck( cudaGetLastError() );
}

}}
