#ifndef RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationCommonKernels
#define RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationCommonKernels

#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes_gpu.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/DeclsForKernels.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

class EcalPulseShape;
            // this flag setting is applied to all of the cases
class EcalPulseCovariance;
class EcalUncalibratedRecHit;

namespace ecal { namespace multifit {

///
/// assume kernel launch configuration is 
/// (MAXSAMPLES * nchannels, blocks)
/// TODO: is there a point to split this kernel further to separate reductions
/// 
__global__
void kernel_prep_1d_and_initialize(EcalPulseShape const* shapes_in,
                    uint16_t const* digis_in,
                    uint32_t const* dids,
                    SampleVector* amplitudes,
                    SampleVector* amplitudesForMinimization,
                    SampleGainVector* gainsNoise,
                    float const* mean_x1,
                    float const* mean_x12,
                    float const* rms_x12,
                    float const* mean_x6,
                    float const* gain6Over1,
                    float const* gain12Over6,
                    bool* hasSwitchToGain6,
                    bool* hasSwitchToGain1,
                    bool* isSaturated,
                    ::ecal::reco::StorageScalarType* energies,
                    ::ecal::reco::StorageScalarType* chi2,
                    ::ecal::reco::StorageScalarType* pedestal,
                    uint32_t *flags,
                    uint32_t *v2rmapping,
                    char *npassive,
                    char *samplesMapping,
                    uint32_t offsetForHashes,
                    bool const gainSwitchUseMaxSampleEB,
                    bool const gainSwitchUseMaxSampleEE,
                    int const nchannels);

}}

#endif // RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationCommonKernels
