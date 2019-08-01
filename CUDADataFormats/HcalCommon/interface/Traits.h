#ifndef CUDADataFormats_HcalCommon_interface_Traits_h
#define CUDADataFormats_HcalCommon_interface_Traits_h

#include <vector>

namespace hcal { namespace common {

enum class Tag : int8_t {
    Vec = 0,
    Ptr = 1,
};

template<Tag tag>
struct AddSize {};

template<>
struct AddSize<Tag::Ptr> {
    uint32_t size;
};

struct ViewStoragePolicy {
    static constexpr Tag TagValue = Tag::Ptr;

    template<typename T>
    struct StorageSelector {
        using type = T*;
    };
};

template<template<typename> typename Allocator = std::allocator>
struct VecStoragePolicy {
    static constexpr auto TagValue = Tag::Vec;

    template<typename T>
    struct StorageSelector {
        using type = std::vector<T, Allocator<T>>;
    };
};

}}

#endif // CUDADataFormats_HcalCommon_interface_Traits_h
