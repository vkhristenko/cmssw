#ifndef CUDADataFormats_HcalDigi_interface_DigiCollection_h
#define CUDADataFormats_HcalDigi_interface_DigiCollection_h

namespace hcal {

struct Flavor01 {
    static constexpr int WORDS_PER_SAMPLE = 1;
    static constexpr int HEADER_WORDS = 1;
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
    static constexpr float WORDS_PER_SAMPLE = 0.5;
    static constexpr int HEADER_WORDS = 1;
};

//
// this is basically a view 
// it does not own the actual memory -> does not reclaim
//
template<typename Flavor>
struct DigiCollection {
    DigiCollection() = default;
    DigiCollection(uint32_t *ids, uint16_t *data, uint32_t ndigis, uint32_t nsamples)
        : ids{ids}, data{data}, ndigis{ndigis}, nsamples{nsamples}
    {}
    DigiCollection(DigiCollection const&) = default;
    DigiCollection& operator=(DigiCollection const&) = default;

    DigiCollection(DigiCollection&&) = default;
    DigiCollection& operator=(DigiCollection&&) = default;

    // stride is statically known
    uint32_t *ids=nullptr;
    uint16_t *data=nullptr;
    uint32_t ndigis;
    uint32_t nsamples;
};

}

#endif // CUDADataFormats_HcalDigi_interface_DigiCollection_h
