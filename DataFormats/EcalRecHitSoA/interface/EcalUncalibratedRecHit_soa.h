#ifndef DataFormats_EcalRecHitSoA_interface_EcalUncalibratedRecHit_soa_h
#define DataFormats_EcalRecHitSoA_interface_EcalUncalibratedRecHit_soa_h

#include <vector>

namespace ecal {

namespace Tag {

struct soa {};
struct aos {};

}

template<typename T, typename L = Tag::soa>
struct type_wrapper {
    using type = std::vector<T>;
};

template<typename T>
struct type_wrapper<T, Tag::aos> {
    using type = T;
};

template<typename L = Tag::soa>
struct UncalibratedRecHit {
    UncalibratedRecHit() = default;
    UncalibratedRecHit(const UncalibratedRecHit&) = default;
    UncalibratedRecHit& operator=(const UncalibratedRecHit&) = default;

    UncalibratedRecHit(UncalibratedRecHit&&) = default;
    UncalibratedRecHit& operator=(UncalibratedRecHit&&) = default;

    typename type_wrapper<float, L>::type amplitude;
    typename type_wrapper<float, L>::type chi2;
    typename type_wrapper<uint32_t, L>::type did;
};

using SoAUncalibratedRecHitCollection = UncalibratedRecHit<Tag::soa>;

}

#endif // RecoLocalCalo_EcalRecAlgos_interface_EcalUncalibratedRecHit_soa_h
