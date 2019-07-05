#include "EventFilter/EcalRawToDigi/interface/UnpackGPU.h"

namespace ecal { namespace raw {

__global__
void kernel_unpack_test(
        unsigned char const* data,
        uint32_t const* offsets,
        int const* feds,
        uint32_t const nbytesTotal) {
    // indices
    auto const ifed = blockIdx.x;

    // offset in bytes
    auto const offset = offsets[ifed];
    // fed id
    auto const fed = feds[ifed];
    // size
    auto const size = ifed==gridDim.x-1 ? nbytesTotal - offset : offsets[ifed+1] - offset;

    // FIXME: debugging
    printf("ifed = %u fed = %d offset = %u size = %u\n", ifed, fed, offset, size);

    // offset to the right raw buffer
    uint64_t const* buffer = reinterpret_cast<uint64_t const*>(data + offset);

    // fed header
    auto const fed_header = buffer[0];
    uint32_t fed_id = (fed_header >> 8) & 0xfff;
    uint32_t bx = (fed_header >> 24) & 0xff;
    uint32_t lv1 = (fed_header >> 32) & 0xffffff;
    uint8_t const bid_fed_header = (fed_header >> 60) & 0xf;

    printf("fed = %d fed_id = %u bx = %u lv1 = %u  bid = 0x%u\n",
        fed, fed_id, bx, lv1, bid_fed_header);

    // dcc header
    auto const dcc_header = buffer[1];
    uint32_t event_length = dcc_header & 0xffffff;
    uint8_t dcc_errors = (dcc_header >> 24) & 0xff;
    uint32_t run_number = (dcc_header >> 32) & 0xffffff;
    uint8_t const word_dcc = (dcc_header >> 56) & 0x3f;
    printf("fed = %d size = %u event_length = %u dcc_errors = %hhu run_number = %u word_dcc = 0x%u\n",
        fed, size, 8*event_length, dcc_errors, run_number, word_dcc);
}

void entryPoint(
        InputDataCPU const& inputCPU, 
        InputDataGPU& inputGPU,
        ConditionsProducts const& conditions,
        cuda::stream_t<>& cudaStream,
        uint32_t const nfedsWithData,
        uint32_t const nbytesTotal) {
    std::cout << "nfeds with data = " << nfedsWithData
              << " nbytesTotal = " << nbytesTotal
              << " data bytes total = " << inputCPU.data.size() << std::endl;

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

    kernel_unpack_test<<<nfedsWithData,1, 0, cudaStream.id()>>>(
        inputGPU.data,
        inputGPU.offsets,
        inputGPU.feds,
        nbytesTotal
    );
    cudaCheck( cudaGetLastError() );
}

}}
