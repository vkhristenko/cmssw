#ifndef CUDADataFormats_HcalDigi_interface_DigiCollection_h
#define CUDADataFormats_HcalDigi_interface_DigiCollection_h

namespace hcal {

struct Flavor01 {
    using adc_type = uint8_t;
    using tdc_type = uint8_t;

    static constexpr int WORDS_PER_SAMPLE = 1;
    static constexpr int HEADER_WORDS = 1;

    static constexpr adc_type adc(uint16_t const* const sample_start) {
        return (*sample_start & 0xff);
    }

    static constexpr tdc_type tdc(uint16_t const* const sample_start) {
        return (*sample_start >> 8) & 0x3f;
    }
};

struct Flavor2 {
    static constexpr int WORDS_PER_SAMPLE = 2;
    static constexpr int HEADER_WORDS = 1;
};

struct Flavor3 {
    static constexpr int WORDS_PER_SAMPLE = 1;
    static constexpr int HEADER_WORDS = 1;
};

struct Flavor4 {
    static constexpr int WORDS_PER_SAMPLE = 1;
    static constexpr int HEADER_WORDS = 1;
};

struct Flavor5 {
    using adc_type = uint8_t;

    static constexpr float WORDS_PER_SAMPLE = 0.5;
    static constexpr int SAMPLES_PER_WORD = 2;
    static constexpr int HEADER_WORDS = 1;

    static constexpr adc_type adc(
            uint16_t const* const sample_start, uint8_t const shifter) {
        return ((*sample_start >> shifter*8) & 0x7f);
    }
};

template<typename Flavor>
constexpr
typename Flavor::adc_type adc_for_sample(
        uint16_t const* const dfstart, uint32_t const sample) {
    return Flavor::adc(
        dfstart + Flavor::HEADER_WORDS + sample * Flavor::WORDS_PER_SAMPLE);
}

template<typename Flavor>
constexpr 
typename Flavor::tdc_type tdc_for_sample(
        uint16_t const* const dfstart, uint32_t const sample) {
    return Flavor::tdc(
        dfstart + Flavor::HEADER_WORDS + sample*Flavor::WORDS_PER_SAMPLE);
}

template<>
constexpr
Flavor5::adc_type adc_for_sample<Flavor5>(
        uint16_t const* const dfstart, uint32_t const sample) {
    // avoid using WORDS_PER_SAMPLE and simply shift 
    return Flavor5::adc(
        dfstart + Flavor5::HEADER_WORDS + (sample >> 1),
        sample%2);
}

template<typename Flavor>
constexpr uint32_t compute_stride(uint32_t const nsamples) { 
    return static_cast<uint32_t>(nsamples * Flavor::WORDS_PER_SAMPLE) + 
        Flavor::HEADER_WORDS;
}

template<typename Flavor>
constexpr uint32_t compute_nsamples(uint32_t const nwords) {
    return (nwords - Flavor::HEADER_WORDS) / Flavor::WORDS_PER_SAMPLE;
}

template<>
constexpr uint32_t compute_nsamples<Flavor5>(uint32_t const nwords) {
    return (nwords - Flavor5::HEADER_WORDS) * Flavor5::SAMPLES_PER_WORD;
}

//
// this is basically a view 
// it does not own the actual memory -> does not reclaim
//
template<typename Flavor>
struct DigiCollection {
    DigiCollection() = default;
    DigiCollection(uint32_t *ids, uint16_t *data, uint32_t ndigis, uint32_t stride)
        : ids{ids}, data{data}, ndigis{ndigis}, stride{stride}
    {}
    DigiCollection(DigiCollection const&) = default;
    DigiCollection& operator=(DigiCollection const&) = default;

    DigiCollection(DigiCollection&&) = default;
    DigiCollection& operator=(DigiCollection&&) = default;

    // stride is statically known
    uint32_t *ids=nullptr;
    uint16_t *data=nullptr;
    uint32_t ndigis;
    uint32_t stride;
};

}

#endif // CUDADataFormats_HcalDigi_interface_DigiCollection_h
