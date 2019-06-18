#ifndef RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernels
#define RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernels

#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes_gpu.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/DeclsForKernels.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

class EcalPulseShape;
class EcalPulseCovariance;
class EcalUncalibratedRecHit;

namespace ecal { namespace multifit {

void minimization_procedure(
        EventInputDataCPU const& eventInputCPU, EventInputDataGPU& eventInputGPU,
        EventOutputDataGPU& eventOutputGPU, EventDataForScratchGPU& scratch,
        ConditionsProducts const& conditions,
        ConfigurationParameters const& configParameters,
        cuda::stream_t<>& cudaStream,
        unsigned int const offsetForHashes);

///
/// TODO: trivial impl for now, there must be a way to improve
///
/// Conventions:
///   - amplitudes -> solution vector, what we are fitting for
///   - samples -> raw detector responses
///   - passive constraint - satisfied constraint
///   - active constraint - unsatisfied (yet) constraint
///
__global__
void kernel_minimize(SampleMatrix const* noisecov,
                     EcalPulseCovariance const* full_pulse_cov,
                     BXVectorType *bxs,
                     SampleVector const* samples,
                     SampleVector* amplitudes,
                     PulseMatrixType* pulse_matrix, 
                     ::ecal::reco::StorageScalarType* chi2s,
                     uint32_t const* dids,
                     uint32_t const* v2rmapping,
                     uint32_t const* noiseCovIsZero,
                     int nchannels,
                     int max_iterations,
                     unsigned int offsetForHashes);

}}

#endif // RecoLocalCalo_EcalRecAlgos_src_AmplitudeComputationKernelsV1
