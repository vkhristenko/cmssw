#ifndef RecoLocalCalo_Common_interface_inplace_fnnls_h
#define RecoLocalCalo_Common_interface_inplace_fnnls_h

#include "RecoLocalCalo/Common/interface/data_types.h"

void
inplace_fnnls(const FixedMatrix& A,
              const FixedVector& b,
              FixedVector& x,
              const double eps = 1e-11,
              const unsigned int max_iterations = 1000);

#endif // RecoLocalCalo_Common_interface_inplace_fnnls_h
