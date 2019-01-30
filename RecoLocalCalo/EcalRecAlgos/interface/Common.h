#ifndef RecoLocalCalo_EcalRecAlgos_interface_Common_h
#define RecoLocalCalo_EcalRecAlgos_interface_Common_h

#include <cmath>

// a workaround for std::abs not being a constexpr function
namespace ecal {

template<typename T>
constexpr T abs(T const& value) {
    return ::std::max(value, -value);
}

}

namespace ecal { namespace cuda {

void assert_if_error();

}}

#endif
