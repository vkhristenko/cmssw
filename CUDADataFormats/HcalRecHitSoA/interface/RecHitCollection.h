#ifndef CUDADataFormats_HcalRecHitCollectionSoA_interface_RecHitCollection_h
#define CUDADataFormats_HcalRecHitCollectionSoA_interface_RecHitCollection_h

#include <vector>

#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"

namespace hcal {

namespace Tag {

struct soa {};
struct ptr {};

}

template<typename T, typename L = Tag::soa>
struct type_wrapper {
    using type = std::vector<T, CUDAHostAllocator<T>>;
};

template<typename T>
struct type_wrapper<T, Tag::ptr> {
    using type = T*;
};

namespace Detail {

// empty base 
template<typename T>
struct Base {};

// add number of values for ptr case
template<>
struct Base<::hcal::Tag::ptr> {
    uint32_t size;
};

}

template<typename L = Tag::soa>
struct RecHitCollection : public Detail::Base<L> {
    RecHitCollection() = default;
    RecHitCollection(const RecHitCollection&) = default;
    RecHitCollection& operator=(const RecHitCollection&) = default;

    RecHitCollection(RecHitCollection&&) = default;
    RecHitCollection& operator=(RecHitCollection&&) = default;

    typename type_wrapper<float, L>::type energy;
    typename type_wrapper<float, L>::type time;
    typename type_wrapper<uint32_t, L>::type did;

    template<typename U = L>
    typename std::enable_if<std::is_same<U, Tag::soa>::value, void>::type 
    resize(size_t size) {
        energy.resize(size);
        time.resize(size);
        did.resize(size);
    }
};

}

#endif // RecoLocalCalo_HcalRecAlgos_interface_RecHitCollection_h
