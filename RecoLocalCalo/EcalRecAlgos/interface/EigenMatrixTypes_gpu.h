#ifndef RecoLocalCalo_EcalRecAlgos_EigenMatrixTypes_gpu_h
#define RecoLocalCalo_EcalRecAlgos_EigenMatrixTypes_gpu_h

#include <Eigen/Dense>
#include <array>

#include "CUDADataFormats/EcalRecHitSoA/interface/RecoTypes.h"

namespace ecal { namespace multifit {

constexpr int SampleVectorSize = 10;
constexpr uint32_t noiseCovIsZeroMask = 0xffffffff;
constexpr uint32_t idIsInvalidMask = 0xffffffff;

using data_type = ::ecal::reco::ComputationScalarType;

typedef Eigen::Matrix<data_type,SampleVectorSize,1> SampleVector;
typedef Eigen::Matrix<char, SampleVectorSize,1> SampleGainVector;
typedef Eigen::Matrix<data_type,SampleVectorSize,SampleVectorSize> SampleMatrix;

}}

#endif
