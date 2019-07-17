#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"

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

    // FIXME: for debugging
    if (ifed > 0)
        return;

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

    // amc 13 header
    auto const amc13word = buffer[1];
    uint8_t const namc = (amc13word >> 52) & 0xf;
    uint8_t const amc13version = (amc13word >> 60) & 0xf;
    uint32_t const amc13OrbitNumber = (amc13word >> 4) & 0xffffffffu;

#ifdef HCAL_RAWDECODE_GPUDEBUG
    printf("fed = %d namc = %u amc13version = %u amc13OrbitNumber = %u\n",
        fed, namc, amc13version, amc13OrbitNumber);
#endif

    uint32_t amcoffset = 0;
    for (uint8_t iamc=0u; iamc < namc; ++iamc) {
        auto const word = buffer[2 + iamc];
        uint16_t const amcid = word & 0xffff;
        int const slot = (word >> 16) & 0xf;
        int const amcBlockNumber = (word >> 20) & 0xff;
        int const amcSize = (word >> 32) & 0xffffff;

#ifdef HCAL_RAWDECODE_GPUDEBUG
        printf("fed = %d amcid = %u slot = %d amcBlockNumber = %d amcSize = %d\n",
            fed, amcid, slot, amcBlockNumber, amcSize);
#endif

        bool const amcmore = ((word >> 61) & 0x1) != 0;
        bool const amcSegmented = ((word >> 60) & 0x1) != 0;
        bool const amcLengthOk = ((word >> 62) & 0x1) != 0;
        bool const amcCROk = ((word >> 56) & 0x1) != 0;
        bool const amcDataPresent = ((word >> 58) & 0x1) != 0;
        bool const amcDataValid = ((word >> 56) & 0x1) != 0;
        bool const amcEnabled = ((word >> 59) & 0x1) != 0;

#ifdef HCAL_RAWDECODE_GPUDEBUG
        printf("fed = %d amcmore = %d amcSegmented = %d, amcLengthOk = %d amcCROk = %d\n>> amcDataPresent = %d amcDataValid = %d amcEnabled = %d\n", 
            fed, static_cast<int>(amcmore),
            static_cast<int>(amcSegmented), static_cast<int>(amcLengthOk), 
            static_cast<int>(amcCROk), static_cast<int>(amcDataPresent), 
            static_cast<int>(amcDataValid), static_cast<int>(amcEnabled));
#endif

        // get to the payload
        auto const* payload64 = buffer + 2 + namc + amcoffset;
        auto const* payload16 = reinterpret_cast<uint16_t const*>(payload64);
        amcoffset += amcSize;

        // uhtr header v1 1st 64 bits
        auto const payload64_w0 = payload64[0];
        uint32_t const data_length64 = payload64_w0 & 0xfffff;
        uint16_t bcn = (payload64_w0 >> 20) & 0xfff;
        uint32_t evn = (payload64_w0 >> 32) & 0xffffff;

#ifdef HCAL_RAWDECODE_GPUDEBUG
        printf("fed = %d data_length64 = %u bcn = %u evn = %u\n",
            fed, data_length64, bcn, evn);
#endif

        // uhtr header v1 2nd 64 bits
        auto const payload64_w1 = payload64[1];
        uint8_t const uhtrcrate = payload64_w1 & 0xff;
        uint8_t const uhtrslot = (payload64_w1 >> 8) & 0xf;
        uint8_t const presamples = (payload64_w1 >> 12) & 0xf;
        uint16_t const orbitN = (payload64_w1 >> 16) & 0xffff;
        uint8_t const firmFlavor = (payload64_w1 >> 32) & 0xff;
        uint8_t const eventType = (payload64_w1 >> 40) & 0xf;
        uint8_t const payloadFormat = (payload64_w1 >> 44) & 0xf;

#ifdef HCAL_RAWDECODE_GPUDEBUG
        printf("fed = %d crate = %u slot = %u presamples = %u\n>>> orbitN = %u firmFlavor = %u eventType = %u payloadFormat = %u\n",
            fed, uhtrcrate, uhtrslot, presamples, orbitN, firmFlavor, eventType, payloadFormat);
#endif

        // this should be filtering out uMNio...
        if (payloadFormat != 1)
            continue;

        // skip uhtr header words
        auto const channelDataSize = data_length64 - 2; // 2 uhtr header v1 words
        auto const* channelDataBuffer64Start = payload64 + 2; // 2 uhtr header v2 wds
        auto const* channelDataBuffer64End = channelDataBuffer64Start 
            + channelDataSize;
        auto const* ptr = reinterpret_cast<uint16_t const*>(channelDataBuffer64Start);
        auto const* end = reinterpret_cast<uint16_t const*>(channelDataBuffer64End); 

        // iterate thru the channel data
        uint8_t wordsPerChannel = 0;
        while (ptr != end) {
            uint8_t const fw_lastbit = (ptr[0] >> 15) & 0x1;
            if (fw_lastbit != 1) {
                wordsPerChannel++;
                ptr++;
                continue;
            } else {

#ifdef HCAL_RAWDECODE_GPUDEBUG
                printf(">>> wordsPerChannel = %u\n", wordsPerChannel);
#endif

                wordsPerChannel = 1;
            }
            uint8_t const flavor = (ptr[0] >> 12) & 0x7;
            uint8_t const channelid = ptr[0] & 0xff;

            switch (flavor) {
            case 0:
            case 1:
            case 3:
            {
                uint8_t fiber = (channelid >> 3) & 0x1f;
                uint8_t fchannel = channelid & 0x7;
                HcalElectronicsId eid{uhtrcrate, uhtrslot, fiber, fchannel, false};

#ifdef HCAL_RAWDECODE_GPUDEBUG
                printf("flavor = %u crate = %u slot = %u channelid = %u fiber = %u fchannel = %u\n",
                    flavor, uhtrcrate, uhtrslot, channelid, fiber, fchannel);
#endif

                break;
            }
            case 2:
            {
                uint8_t fiber = (channelid >> 3) & 0x1f;
                uint8_t fchannel = channelid & 0x7;
                HcalElectronicsId eid{uhtrcrate, uhtrslot, fiber, fchannel, false};

#ifdef HCAL_RAWDECODE_GPUDEBUG
                printf("flavor = %u crate = %u slot = %u channelid = %u fiber = %u fchannel = %u\n",
                    flavor, uhtrcrate, uhtrslot, channelid, fiber, fchannel);
#endif

                break;
            }
            case 4:
            {
                uint8_t link = (channelid >> 4) & 0xf;
                uint8_t tower = channelid & 0xf;
                HcalElectronicsId eid{uhtrcrate, uhtrslot, link, tower, true};

#ifdef HCAL_RAWDECODE_GPUDEBUG
                printf("flavor = %u crate = %u slot = %u channelid = %u link = %u tower = %u\n",
                    flavor, uhtrcrate, uhtrslot, channelid, link, tower);
#endif

                break;
            }
            case 5:
            case 7:
            {
                uint8_t const fiber = (channelid >> 2) & 0x3f;
                uint8_t const fchannel = channelid & 0x3;
                HcalElectronicsId eid{uhtrcrate, uhtrslot, fiber, fchannel, false};

#ifdef HCAL_RAWDECODE_GPUDEBUG
                printf("flavor = %u crate = %u slot = %u channelid = %u fiber = %u fchannel = %u\n",
                    flavor, uhtrcrate, uhtrslot, channelid, fiber, fchannel);
#endif

                break;
            }
            default:
#ifdef HCAL_RAWDECODE_GPUDEBUG
                printf("flavor = %u crate = %u slot = %u channelid = %u\n",
                    flavor, uhtrcrate, uhtrslot, channelid);
#else
                ;
#endif
            }

/*
#ifdef HCAL_RAWDECODE_GPUDEBUG
            printf("fed = %d crate = %u slot = %u flavor = %u fw_lastbit = %u channelid = %u\n", fed, uhtrcrate, uhtrslot, flavor, fw_lastbit, channelid);
#endif
*/

            ptr++;
        }
    }
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
