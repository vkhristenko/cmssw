#ifndef RecoLocalCalo_EcalRecAlgos_interface_Common_h
#define RecoLocalCalo_EcalRecAlgos_interface_Common_h

#include <cstdint>
#include <cmath>
#include <cassert>

#include <cuda.h>
#include <cuda_runtime.h>

// a workaround for std::abs not being a constexpr function
namespace ecal {

template<typename T>
constexpr T abs(T const& value) {
    return ::std::max(value, -value);
}

// temporary
namespace mgpa {

constexpr int adc(uint16_t sample) { return sample & 0xfff; }
constexpr int gainId(uint16_t sample) { return (sample>>12) & 0x3; }

}

}

#define AssertIfError \
    { \
        auto code = cudaGetLastError(); \
        if (code != cudaSuccess) { \
            std::cout << cudaGetErrorString(code) << std::endl; \
            assert(false); \
        } \
    }

#endif
