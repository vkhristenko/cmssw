#include <iostream>
#include <limits>

#include "cuda.h"

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "DataFormats/Math/interface/approx_exp.h"
#include "DataFormats/Math/interface/approx_log.h"

#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "AmplitudeComputationKernels.h"
#include "AmplitudeComputationCommonKernels.h"
#include "KernelHelpers.h"

namespace ecal { namespace multifit { 
///
/// The following conventions apply to all of the kernels below
///     - virtual channel id - channel id for the current iteration
///     - real  channel id - channel id for the original array
/// 

/// launch ctx:
/// 10 x 10 x nchannels per block
// FIXME: add __restrict__ whenever is needed
__global__
void kernel_newiter_update_covariance_compute_cholesky(
        EcalPulseCovariance const* g_PulseCovariance,
        SampleVector const* g_amplitudes,
        SampleVector::Scalar* g_L,
        uint32_t const* v2ridmapping,
        uint32_t *dids,
        uint32_t *pChannelsCounter,
        SampleGainVector const* gainNoise,
        float const* rms_x12,
        float const* rms_x6,
        float const* rms_x1,
        float const* gain12Over6,
        float const* gain6Over1,
        double const* G12SamplesCorrelationEB,
        double const* G6SamplesCorrelationEB,
        double const* G1SamplesCorrelationEB,
        double const* G12SamplesCorrelationEE,
        double const* G6SamplesCorrelationEE,
        double const* G1SamplesCorrelationEE,
        bool const* hasSwitchToGain6,
        bool const* hasSwitchToGain1,
        bool const* isSaturated,
        uint32_t const offsetForHashes,
        uint32_t const iteration,
        uint32_t const nchannels) {
    // constants / type aliases
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    constexpr auto npulses = nsamples;
    constexpr auto template_samples = EcalPulseShape::TEMPLATESAMPLES;
    using DataType = SampleVector::Scalar;
    using ConstDataType = DataType const;
    constexpr auto nvaluesForL = MapSymM<DataType, nsamples>::total;
    constexpr bool simplifiedNoiseModelForGainSwitch = true; 
    constexpr float addPedestalUncertainty = 0.f;
    constexpr bool dynamicPedestal = false;

    // indices
    int const tx = threadIdx.x;
    int const ty = threadIdx.y;
    auto const vch_per_block = threadIdx.z;
    auto const gvch = vch_per_block + blockIdx.x*blockDim.z;
    if (gvch >= nchannels) return;

    // done for each new iteration...
    if (tx==0 && ty==0 && gvch==0) {
        *pChannelsCounter = 0;
    }

    auto const grch = v2ridmapping[gvch];

    // non-divergent branch!
    // only for the first iteration...
    if (iteration == 0) {
        // we need to filter out channels that either
        // - have values that are precomputed
        // - noise matrix is 0
        if (grch==idIsInvalidMask)
            return;
    }
    
    auto const did = DetId{dids[grch]};
    auto const isBarrel = did.subdetId() == EcalBarrel;
    auto const hashedId = isBarrel
        ? hashedIndexEB(did.rawId())
        : offsetForHashes + hashedIndexEE(did.rawId());
    bool tmp0 = hasSwitchToGain6[grch];
    bool tmp1 = hasSwitchToGain1[grch];
    auto const* G12SamplesCorrelation = isBarrel
        ? G12SamplesCorrelationEB
        : G12SamplesCorrelationEE;
    auto const* G6SamplesCorrelation = isBarrel
        ? G6SamplesCorrelationEB
        : G6SamplesCorrelationEE;
    auto const* G1SamplesCorrelation = isBarrel
        ? G1SamplesCorrelationEB
        : G1SamplesCorrelationEE;
    bool tmp2 = isSaturated[grch];
    bool hasGainSwitch = tmp0 || tmp1 || tmp2;
    auto const vidx = ecal::abs<int>(ty - tx);

    // non-divergent branch for all threads per cnannel
    float noiseValue = 0;
    if (hasGainSwitch) {
        // non-divergent branch - all threads per block
        if (simplifiedNoiseModelForGainSwitch) {
            int isample_max = 5; // according to cpu defs
            int gainidx = gainNoise[grch][isample_max];

            // non-divergent branches
            if (gainidx==0)
                noiseValue = rms_x12[hashedId]*rms_x12[hashedId]
                    * G12SamplesCorrelation[vidx];
            if (gainidx==1) 
                noiseValue = gain12Over6[hashedId]*gain12Over6[hashedId] 
                    * rms_x6[hashedId]*rms_x6[hashedId]
                    * G6SamplesCorrelation[vidx];
            if (gainidx==2)
                noiseValue = gain12Over6[hashedId]*gain12Over6[hashedId]
                    * gain6Over1[hashedId]*gain6Over1[hashedId] 
                    * rms_x1[hashedId]*rms_x1[hashedId]
                    * G1SamplesCorrelation[vidx];
            if (!dynamicPedestal && addPedestalUncertainty>0.f)
                noiseValue += addPedestalUncertainty*addPedestalUncertainty;
        } else {
            int gainidx=0;
            char mask = gainidx;
            int pedestal = gainNoise[grch][ty] == mask ? 1 : 0;
            noiseValue += /* gainratio is 1*/ rms_x12[hashedId]*rms_x12[hashedId]
                * pedestal* G12SamplesCorrelation[vidx];
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noiseValue += /* gainratio is 1 */
                    addPedestalUncertainty*addPedestalUncertainty*pedestal;
            }

            //
            gainidx=1;
            mask = gainidx;
            pedestal = gainNoise[grch][ty] == mask ? 1 : 0;
            noiseValue += gain12Over6[hashedId]*gain12Over6[hashedId]
                *rms_x6[hashedId]*rms_x6[hashedId]*pedestal
                * G6SamplesCorrelation[vidx];
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noiseValue += gain12Over6[hashedId]*gain12Over6[hashedId]
                    *addPedestalUncertainty*addPedestalUncertainty
                    *pedestal;
            }
            
            //
            gainidx=2;
            mask = gainidx;
            pedestal = gainNoise[grch][ty] == mask ? 1 : 0;
            float tmp = gain6Over1[hashedId] * gain12Over6[hashedId];
            noiseValue+= tmp*tmp * rms_x1[hashedId]*rms_x1[hashedId]
                *pedestal* G1SamplesCorrelation[vidx];
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noiseValue += tmp*tmp * addPedestalUncertainty*addPedestalUncertainty
                    * pedestal;
            }
        }
    } else {
        auto rms = rms_x12[hashedId];
        noiseValue = rms*rms * G12SamplesCorrelation[vidx];
        if (!dynamicPedestal && addPedestalUncertainty>0.f) {
            //----  add fully correlated component to noise 
            // covariance to inflate pedestal uncertainty
            noiseValue += addPedestalUncertainty*addPedestalUncertainty;
        }
    }
    
    // shared mem configuration
    extern __shared__ char smem[];
    // nchannels per block * 12 * 12 (TEMPLATESAMPLES)
    //DataType* __shrPulseCovariance = reinterpret_cast<DataType*>(smem);
    // nchannels per block * 10
    DataType* __shrCovarianceDiagValues = reinterpret_cast<DataType*>(smem);
    //    __shrPulseCovariance + 
    //    blockDim.z * template_samples * template_samples;
    // nchannels per block * 55
    DataType* __shrL = __shrCovarianceDiagValues + 
        blockDim.z * nsamples;
    // nchannels per block * 10 * 10
    DataType* __shrLSums = __shrL + 
        blockDim.z * nvaluesForL;

    // map global mem
    MapM<float const, template_samples, Eigen::RowMajor> PulseCovariance
    {reinterpret_cast<float const*>(g_PulseCovariance + hashedId)};
    MapV<ConstDataType> amplitudes{g_amplitudes[grch].data()};
    MapSymM<DataType, nsamples, Eigen::RowMajor> L
    {g_L + nvaluesForL * grch};

    // map allocated shared mem
    MapV<DataType> shrCovarianceDiagValues
    {__shrCovarianceDiagValues + vch_per_block*nsamples};
    MapSymM<DataType, nsamples, Eigen::RowMajor> shrL
    {__shrL + vch_per_block * nvaluesForL};
    MapSymM<DataType, nsamples> shrLSums
    {__shrLSums + vch_per_block * nvaluesForL};

    // preload some values
    //auto noiseValue = noise_value;
    if (ty>=tx)
        shrLSums(ty, tx) = 0;

    // compute the updated total covariance matrix
    #pragma unroll
    for (int ipulse=0; ipulse<npulses; ipulse++) {
        auto const amplitude = amplitudes(ipulse);
        if (amplitude == 0)
            continue;

        // note, we do not permute results anymore!
        // therefore straightforward mapping
        auto const bx = ipulse - 5;
        auto const first_sample_t = std::max(0, bx + 3);
        if (!(tx >= first_sample_t && ty >= first_sample_t))
            continue;
        
        // note: we are no longer using 19 x 19 pulse covariance matrix,
        // just what conditions provide directly
        auto const offset = -3 - bx;
        auto const amp_sq = amplitude * amplitude;
        auto const nsample_pulse = nsamples - first_sample_t;

        // update matrix eleement
        noiseValue += amp_sq * PulseCovariance(ty + offset, tx + offset);
    }

    // we need to store to shared mem diagonal values
    if (ty == tx)
        shrCovarianceDiagValues(tx) = noiseValue;
    __syncthreads();

    //
    // cholesky decomposition of 10 x 10
    // noiseValue is the (i, j) element of a matrix to decompose
    //

    // column 0
    // compute L(ty, 0) for ty >= 1
    // note, important we get m_j_j = m(tx) -> we need same value for 
    // all values of a given column
    auto const m_j_j = shrCovarianceDiagValues(tx);
    auto const m_i_j = noiseValue;
    if (tx == 0 && ty == 0)
        shrL(0, 0) = std::sqrt(m_j_j);
    else if (tx==0 && ty >=1)
        shrL(ty, 0) = noiseValue / std::sqrt(m_j_j);
    __syncthreads();

    // TODO verify that hte loop is unrolled
    #pragma unroll
    for (int column=1; column<nsamples; ++column) {
        if (tx==column && ty>=column) {
            // compute L(j, j) = sqrt(M(j, j) - Sum[k](L(j, k) * L(j, k)))
            auto const sumsq = shrLSums(column, column) + 
                shrL(column, column-1) * shrL(column, column-1);
            auto const l_j_j = std::sqrt(m_j_j - sumsq);
            if (ty == column)
                shrL(column, column) = l_j_j;
            else {
                auto const sumsq_i_j = shrLSums(ty, column) + 
                    shrL(ty, column-1) * shrL(column, column-1);
                auto const tmp = m_i_j - sumsq_i_j;
                auto const l_i_column = tmp / l_j_j;
                shrL(ty, column) = l_i_column;
            }
        }

        if (tx>column && ty>=tx)
            shrLSums(ty, tx) += shrL(ty, column-1) * shrL(tx, column-1);
        __syncthreads();
    }

    // store back to global
    if (ty>=tx)
        L(ty, tx) = shrL(ty, tx);
}

/// launch ctx:
/// 10 x 10 x nchannels per block
// FIXME: add __restrict__ whenever is needed
__global__
void kernel_solve_mm_mv_mults(
        EcalPulseShape const* g_pulse_shape,
        SampleVector::Scalar const* g_L,
        SampleVector const* g_s,
        SampleVector::Scalar* g_AtA,
        SampleVector::Scalar* g_Atb,
        uint32_t const* v2ridmapping,
        uint32_t *dids,
        uint32_t const offsetForHashes,
        uint32_t const iteration,
        uint32_t const nchannels) {
    // constants, typedefs
    using DataType = SampleVector::Scalar;
    using ConstDataType = DataType const;
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    constexpr auto nvaluesForL = MapSymM<DataType, nsamples>::total;
    constexpr auto nvaluesForAtA = nvaluesForL;
    constexpr auto nvaluesForPulseShape=EcalPulseShape::TEMPLATESAMPLES;

    // indices
    auto const tx = threadIdx.x;
    auto const ty = threadIdx.y;
    auto const lch = threadIdx.z;
    auto const gvch = lch + blockIdx.x * blockDim.z;
    if (gvch >= nchannels) return;

    // get teh real ch id
    auto const grch = v2ridmapping[gvch];
    
    // non-divergent branch!
    // only for the first iteration...
    if (iteration == 0) {
        // we need to filter out channels that either
        // - have values that are precomputed
        // - noise matrix is 0
        if (grch==idIsInvalidMask)
            return;
    }
    
    // get did and hashedId
    auto const did = DetId{dids[grch]};
    auto const isBarrel = did.subdetId() == EcalBarrel;
    auto const hashedId = isBarrel
        ? hashedIndexEB(did.rawId())
        : offsetForHashes + hashedIndexEE(did.rawId());

    // configure shared mem
    extern __shared__ char smem[];
    DataType* __shrL = reinterpret_cast<DataType*>(smem);
    float* __shrP = __shrL + nvaluesForL * blockDim.z;
    DataType* __shrs = __shrP + nvaluesForPulseShape * blockDim.z;
    DataType* __shrAtmp = __shrs + nsamples * blockDim.z;
    DataType* __shrbtmp = __shrAtmp + nsamples * nsamples * blockDim.z;

    // map global mem
    // NOTE: this guy needs to be in sync with what is used above/below
    MapSymM<ConstDataType, nsamples, Eigen::RowMajor> L
    {g_L + grch*nvaluesForL};
    MapMForPM<float> pulse_shape {(float*)&g_pulse_shape[hashedId]};
    MapV<ConstDataType> s{g_s[grch].data()};
    MapSymM<DataType, nsamples, Eigen::ColMajor> AtA
    {g_AtA + grch * nvaluesForAtA};
    MapV<DataType> Atb{g_Atb + grch * nsamples};

    // map shared mem
    MapSymM<DataType, nsamples, Eigen::RowMajor> shrL
    {__shrL + lch * nvaluesForL};
    MapMForPM<float> shrP
    {__shrP + lch*nvaluesForPulseShape};
    MapV<DataType> shrs{__shrs + lch*nsamples};
    MapM<DataType, nsamples, Eigen::RowMajor> shrAtmp
    {__shrAtmp + lch*nsamples*nsamples};
    MapV<DataType> shrbtmp{__shrbtmp + lch*nsamples};

    // move into shared mem
    //shrP(ty, tx) = P(ty, tx);
    if (ty>=tx)
        shrL(ty, tx) = L(ty, tx);
    if (ty == 0)
        shrs(tx) = s(tx);
    // store to shared memory. note we are using plain buffers
    if (ty == 1)
        shrP.data[tx] = pulse_shape.data[tx];
    if (ty == 2 && tx<2)
        shrP.data[10 + tx] = pulse_shape.data[10 + tx];
    __syncthreads();

    // run forward substitution
    // use 10 thraeds for matrix and 1 thread vector
    // TODO: would it be better with more threads + synching???
    if (ty==0) {
        ForwardSubstitutionUnrolled<DataType, nsamples>::compute(
            shrL, shrP, shrAtmp, tx);
    } else if (ty==1 && tx == 1) {
        ForwardSubstitutionUnrolled<DataType, nsamples>::compute(
            shrL, shrs, shrbtmp);
    }
    __syncthreads();

    // matrix matrix and matrix vector mult
    // shrAtmp is matrix A
    // shrbtmp is vector b
    DataType ata_i_j = 0;
    #pragma unroll
    for (int i=0; i<nsamples; ++i) {
        auto const A_i_ty = shrAtmp(i, ty);
        ata_i_j += shrAtmp(i, ty) * shrAtmp(i, tx);
    }

    // store back to global. note AtA is symmetric only L is stored
    if (ty>=tx)
        AtA(ty, tx) = ata_i_j;
    if (ty == 0) {
        DataType sum = 0;
        #pragma unroll
        for (int i=0; i<nsamples; i++)
            sum += shrAtmp(i, tx) * shrbtmp(i);

        // store back to global
        Atb(tx) = sum;
    }
}

// must have single channel per block
// cause __syncthreads will deadlock threads for other channel...
__global__
void kernel_minimize_fused(
        EcalPulseCovariance const* g_PulseCovariance,
        EcalPulseShape const* g_pulse_shape,
        SampleVector *g_x,
        SampleVector const* g_s,
        char const* g_mapping,
        char const* g_npassive,
        ::ecal::reco::StorageScalarType *energies,
        ::ecal::reco::StorageScalarType *g_chi2,
        uint32_t const* v2ridmapping,
        uint32_t *dids,
        SampleGainVector const* gainNoise,
        float const* rms_x12,
        float const* rms_x6,
        float const* rms_x1,
        float const* gain12Over6,
        float const* gain6Over1,
        double const* G12SamplesCorrelationEB,
        double const* G6SamplesCorrelationEB,
        double const* G1SamplesCorrelationEB,
        double const* G12SamplesCorrelationEE,
        double const* G6SamplesCorrelationEE,
        double const* G1SamplesCorrelationEE,
        bool const* hasSwitchToGain6,
        bool const* hasSwitchToGain1,
        bool const* isSaturated,
        uint32_t const offsetForHashes,
        uint32_t const startingIteration,
        uint32_t const nchannels) {
    // constants / type aliases
    using DataType = SampleVector::Scalar;
    using ConstDataType = DataType const;
    constexpr double eps = 1e-11;
    constexpr unsigned int max_iterations = 500;
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    constexpr auto nvaluesForL = MapSymM<DataType, nsamples>::total;
    constexpr auto nvaluesForAtA = nvaluesForL;
    constexpr auto nvaluesForPulseShape=EcalPulseShape::TEMPLATESAMPLES;
    constexpr bool simplifiedNoiseModelForGainSwitch = true;
    constexpr float addPedestalUncertainty = 0.f;
    constexpr bool dynamicPedestal = false;
    
    // indices
    auto const tx = threadIdx.x;
    auto const ty = threadIdx.y;
    auto const gvch = blockIdx.x;
    if (gvch >= nchannels) return;
    auto const grch = v2ridmapping[gvch];

    // non-divergent branch!
    // only for the zero iteration...
    if (startingIteration == 0) {
        // we need to filter out channels that either
        // - have values that are precomputed
        // - noise matrix is 0
        if (grch==idIsInvalidMask)
            return;
    }
   
    // conditions
    auto const did = DetId{dids[grch]};
    auto const isBarrel = did.subdetId() == EcalBarrel;
    auto const hashedId = isBarrel
        ? hashedIndexEB(did.rawId())
        : offsetForHashes + hashedIndexEE(did.rawId());
    bool tmp0 = hasSwitchToGain6[grch];
    bool tmp1 = hasSwitchToGain1[grch];
    auto const* G12SamplesCorrelation = isBarrel
        ? G12SamplesCorrelationEB
        : G12SamplesCorrelationEE;
    auto const* G6SamplesCorrelation = isBarrel
        ? G6SamplesCorrelationEB
        : G6SamplesCorrelationEE;
    auto const* G1SamplesCorrelation = isBarrel
        ? G1SamplesCorrelationEB
        : G1SamplesCorrelationEE;
    bool tmp2 = isSaturated[grch];
    bool hasGainSwitch = tmp0 || tmp1 || tmp2;
    auto const vidx = ecal::abs<int>(ty - tx);
   
    // compute Noise Covariance Matrix
    // non-divergent branch for all threads per cnannel
    float noiseValue = 0;
    if (hasGainSwitch) {
        // non-divergent branch - all threads per block
        if (simplifiedNoiseModelForGainSwitch) {
            int isample_max = 5; // according to cpu defs
            int gainidx = gainNoise[grch][isample_max];

            // non-divergent branches
            if (gainidx==0)
                noiseValue = rms_x12[hashedId]*rms_x12[hashedId]
                    * G12SamplesCorrelation[vidx];
            if (gainidx==1) 
                noiseValue = gain12Over6[hashedId]*gain12Over6[hashedId] 
                    * rms_x6[hashedId]*rms_x6[hashedId]
                    * G6SamplesCorrelation[vidx];
            if (gainidx==2)
                noiseValue = gain12Over6[hashedId]*gain12Over6[hashedId]
                    * gain6Over1[hashedId]*gain6Over1[hashedId] 
                    * rms_x1[hashedId]*rms_x1[hashedId]
                    * G1SamplesCorrelation[vidx];
            if (!dynamicPedestal && addPedestalUncertainty>0.f)
                noiseValue += addPedestalUncertainty*addPedestalUncertainty;
        } else {
            int gainidx=0;
            char mask = gainidx;
            int pedestal = gainNoise[grch][ty] == mask ? 1 : 0;
            noiseValue += /* gainratio is 1*/ rms_x12[hashedId]*rms_x12[hashedId]
                * pedestal* G12SamplesCorrelation[vidx];
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noiseValue += /* gainratio is 1 */
                    addPedestalUncertainty*addPedestalUncertainty*pedestal;
            }

            //
            gainidx=1;
            mask = gainidx;
            pedestal = gainNoise[grch][ty] == mask ? 1 : 0;
            noiseValue += gain12Over6[hashedId]*gain12Over6[hashedId]
                *rms_x6[hashedId]*rms_x6[hashedId]*pedestal
                * G6SamplesCorrelation[vidx];
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noiseValue += gain12Over6[hashedId]*gain12Over6[hashedId]
                    *addPedestalUncertainty*addPedestalUncertainty
                    *pedestal;
            }
            
            //
            gainidx=2;
            mask = gainidx;
            pedestal = gainNoise[grch][ty] == mask ? 1 : 0;
            float tmp = gain6Over1[hashedId] * gain12Over6[hashedId];
            noiseValue+= tmp*tmp * rms_x1[hashedId]*rms_x1[hashedId]
                *pedestal* G1SamplesCorrelation[vidx];
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noiseValue += tmp*tmp * addPedestalUncertainty*addPedestalUncertainty
                    * pedestal;
            }
        }
    } else {
        auto rms = rms_x12[hashedId];
        noiseValue = rms*rms * G12SamplesCorrelation[vidx];
        if (!dynamicPedestal && addPedestalUncertainty>0.f) {
            //----  add fully correlated component to noise 
            // covariance to inflate pedestal uncertainty
            noiseValue += addPedestalUncertainty*addPedestalUncertainty;
        }
    }

    // configure shared memory for what can not be reused between iterations
    __shared__ DataType __shrs[nsamples];
    __shared__ DataType __shrx[nsamples];
    __shared__ float __shrP[nvaluesForPulseShape]; 
    __shared__ DataType shr_current_chi2;
    __shared__ bool shrFinished;
    __shared__ char __shrmapping[nsamples];

    __shared__ DataType __shrL[nvaluesForL];
    __shared__ DataType __shrAtmp[nsamples*nsamples];
    __shared__ DataType __shrbtmp[nsamples];
    __shared__ DataType __shrAtA[nvaluesForAtA];
    __shared__ DataType __shrAtb[nsamples];

    // configure shared memory for what can be reused - scratch
    DataType *__shrLSums = __shrAtmp;
    DataType *__shrCovarianceDiagValues = __shrbtmp;
    DataType *__shrLtmp = __shrAtmp;
    DataType *__shrstmp = __shrbtmp;
    DataType *__shrtmp = __shrAtmp + nvaluesForL;
    DataType *__shrv = __shrAtmp;
   
    // map global mem
    MapM<float const, nvaluesForPulseShape, Eigen::RowMajor> PulseCovariance
    {reinterpret_cast<float const*>(g_PulseCovariance + hashedId)};
    MapMForPM<float> pulse_shape{(float*)&g_pulse_shape[hashedId]};
    MapV<ConstDataType> s{g_s[grch].data()};
    MapV<DataType> x{g_x[grch].data()};
    MapV<char const> mapping{g_mapping + grch*nsamples};
    
    // map allocated and to be reused shared mem
    MapV<DataType> shrCovarianceDiagValues{__shrCovarianceDiagValues};
    MapSymM<DataType, nsamples, Eigen::RowMajor> shrL{__shrL};
    MapSymM<DataType, nsamples> shrLSums{__shrLSums};
    MapV<DataType> shrx{__shrx};
    MapV<DataType> shrs{__shrs};
    MapV<char> shrmapping{__shrmapping};
    MapMForPM<float> shrP{__shrP};
    
    MapM<DataType, nsamples, Eigen::RowMajor> shrAtmp{__shrAtmp};
    MapV<DataType> shrbtmp{__shrbtmp};
    MapSymM<DataType, nsamples, Eigen::ColMajor> shrAtA{__shrAtA};
    MapV<DataType> shrAtb{__shrAtb};
        
    MapSymMWithCheck<ConstDataType, nsamples, Eigen::ColMajor> 
    shrAtAWithCheck{shrAtA.data};
    MapV<DataType> shrstmp{__shrstmp};
    MapV<DataType> shrv{__shrv};
    MapV<DataType> shrtmp{__shrtmp};
    MapSymM<DataType, nsamples, Eigen::RowMajor> shrLtmp{__shrLtmp};
    
    // initialization
    // pull from global or init to 0 depending on config
    if (ty==0 && tx==0) {
        shrFinished = false;
        shr_current_chi2 = startingIteration==0 ? 0 : g_chi2[grch];
    }

    // mov pulse shape from global to shared
    if (ty == 1)
        shrP.data[tx] = pulse_shape.data[tx];
    if (ty == 2 && tx<2)
        shrP.data[10 + tx] = pulse_shape.data[10 + tx];

    // initialize amplitudes
    // if we are running minimization fully from scratch, set 0s
    // if we are continueing, load/store from global mem to shared
    if (ty==4)
        shrx(tx) = startingIteration==0 ? 0 : x(tx) ;
    if (ty == 3)
        shrs(tx) = s(tx);
    if (ty==5) 
        shrmapping(tx) = startingIteration ? tx : mapping(tx);
    __syncthreads();

    // init number of passive values in the set
    auto nPassive = startingIteration==0 ? 0 : g_npassive[grch];

    // main outer loop
    auto iteration = startingIteration;
    for (; iteration<50; ++iteration) {
        // compute the updated total covariance matrix
        #pragma unroll
        for (int ipulse=0; ipulse<nsamples; ipulse++) {
            auto const amplitude = shrx(ipulse);
            if (amplitude == 0)
                continue;

            // note, we do not permute results anymore!
            // therefore straightforward mapping
            auto const bx = ipulse - 5;
            auto const first_sample_t = std::max(0, bx + 3);
            if (!(tx >= first_sample_t && ty >= first_sample_t))
                continue;
            
            // note: we are no longer using 19 x 19 pulse covariance matrix,
            // just what conditions provide directly
            auto const offset = -3 - bx;
            auto const amp_sq = amplitude * amplitude;
            auto const nsample_pulse = nsamples - first_sample_t;

            // update matrix eleement
            noiseValue += amp_sq * PulseCovariance(ty + offset, tx + offset);
        }

        // we need to store to shared mem diagonal values
        if (ty == tx)
            shrCovarianceDiagValues(tx) = noiseValue;
        if (ty>=tx)
            shrLSums(ty, tx) = 0;
        __syncthreads();

        //
        // cholesky decomposition of 10 x 10
        // noiseValue is the (i, j) element of a matrix to decompose
        //

        // column 0
        // compute L(ty, 0) for ty >= 1
        // note, important we get m_j_j = m(tx) -> we need same value for 
        // all values of a given column
        auto const m_j_j = shrCovarianceDiagValues(tx);
        auto const m_i_j = noiseValue;
        if (tx == 0 && ty == 0)
            shrL(0, 0) = std::sqrt(m_j_j);
        else if (tx==0 && ty >=1)
            shrL(ty, 0) = noiseValue / std::sqrt(m_j_j);
        __syncthreads();

        // TODO verify that hte loop is unrolled
        #pragma unroll
        for (int column=1; column<nsamples; ++column) {
            if (tx==column && ty>=column) {
                // compute L(j, j) = sqrt(M(j, j) - Sum[k](L(j, k) * L(j, k)))
                auto const sumsq = shrLSums(column, column) + 
                    shrL(column, column-1) * shrL(column, column-1);
                auto const l_j_j = std::sqrt(m_j_j - sumsq);
                if (ty == column)
                    shrL(column, column) = l_j_j;
                else {
                    auto const sumsq_i_j = shrLSums(ty, column) + 
                        shrL(ty, column-1) * shrL(column, column-1);
                    auto const tmp = m_i_j - sumsq_i_j;
                    auto const l_i_column = tmp / l_j_j;
                    shrL(ty, column) = l_i_column;
                }
            }

            if (tx>column && ty>=tx)
                shrLSums(ty, tx) += shrL(ty, column-1) * shrL(tx, column-1);
            __syncthreads();
        }

        // shrL contains main loop Cholesky Decomp
        // shrLSums can be reused from this point
        // shrCovarianceDiagValues can be reused from this point

        // run forward substitution
        // use 10 thraeds for matrix and 1 thread vector
        // TODO: would it be better with more threads + synching???
        if (ty==0) {
            ForwardSubstitutionUnrolled<DataType, nsamples>::compute(
                shrL, shrP, shrAtmp, tx);
        } else if (ty==1 && tx == 1) {
            ForwardSubstitutionUnrolled<DataType, nsamples>::compute(
                shrL, shrs, shrbtmp);
        }
        __syncthreads();

        // matrix matrix and matrix vector mult
        // shrAtmp is matrix A
        // shrbtmp is vector b
        DataType ata_i_j = 0;
        #pragma unroll
        for (int i=0; i<nsamples; ++i) {
            auto const A_i_ty = shrAtmp(i, ty);
            ata_i_j += shrAtmp(i, ty) * shrAtmp(i, tx);
        }

        // store back to global. note AtA is symmetric only L is stored
        if (ty>=tx)
            shrAtA(ty, tx) = ata_i_j;
        if (ty == 0) {
            DataType sum = 0;
            #pragma unroll
            for (int i=0; i<nsamples; i++)
                sum += shrAtmp(i, tx) * shrbtmp(i);

            // store back to global
            shrAtb(tx) = sum;
        }
        __syncthreads();

        //
        // perform fnnls
        // 

        // TODO: can we use more than 1 thread in here!?
        if (tx==0 && ty==0) {
            // used to break from the loop
            unsigned int w_max_idx_prev = 0;
            DataType w_max_prev = 0;
            constexpr double eps = 1e-11;
            auto eps_to_use = eps;

            // main loop
            for (int iter = 0; iter < max_iterations; ++iter) {
                if (iter>0 || nPassive==0) {
                    const auto nActive = nsamples - nPassive;

                    if(!nActive)
                      break;

                    //  
                    unsigned int w_max_idx = -1;
                    auto max_w {static_cast<DataType>(-1)};
                    for (unsigned int i=nsamples-nActive; i<nsamples; i++) {
                        auto sum_per_row{static_cast<DataType>(0)};
                        auto const real_i = shrmapping(i);
                        auto const atb = shrAtb(real_i);
                        #pragma unroll
                        for (unsigned int k=0; k<nsamples; ++k)
                            // note, we do not need to look up k in the mapping
                            // both AtA and x have swaps applied -> therefore dot product will 
                            // not change per row
                            sum_per_row += shrAtAWithCheck(real_i, k) * shrx(k);

                        // compute gradient value and check if it is greater than max
                        auto const wvalue = atb - sum_per_row;
                        if (max_w < wvalue) {
                            max_w = wvalue;
                            w_max_idx = i;
                        }
                    }

                    // check for convergence
                    if (max_w < eps_to_use || 
                        (w_max_idx==w_max_idx_prev && max_w==w_max_prev))
                      break;

                    // update
                    w_max_prev = max_w;
                    w_max_idx_prev = w_max_idx;

                    // swap samples to maintain passive vs active set
                    Eigen::numext::swap(
                            shrmapping(nPassive), 
                            shrmapping(w_max_idx));
                    ++nPassive;
                }

                // inner loop
                while (nPassive > 0) {
                  switch (nPassive) {
                  case 1:
                      FusedCholeskySolver<DataType, 1>::compute(
                        shrAtAWithCheck, shrAtb, shrstmp, shrmapping);
                      break;
                  case 2:
                      FusedCholeskySolver<DataType, 2>::compute(
                        shrAtAWithCheck, shrAtb, shrstmp, shrmapping);
                      break;
                  case 3:
                      FusedCholeskySolver<DataType, 3>::compute(
                        shrAtAWithCheck, shrAtb, shrstmp, shrmapping);
                      break;
                  case 4:
                      FusedCholeskyForwardSubstUnrolled<DataType, 4>::compute(
                        shrAtAWithCheck, shrAtb, shrLtmp, shrtmp, shrmapping);
                      BackwardSubstitutionUnrolled<DataType, 4>::compute(
                        shrLtmp, shrtmp, shrstmp);
                      break;
                  case 5:
                      FusedCholeskyForwardSubstUnrolled<DataType, 5>::compute(
                        shrAtAWithCheck, shrAtb, shrLtmp, shrtmp, shrmapping);
                      BackwardSubstitutionUnrolled<DataType, 5>::compute(
                        shrLtmp, shrtmp, shrstmp);
                      break;
                  case 6:
                      FusedCholeskyForwardSubstUnrolled<DataType, 6>::compute(
                        shrAtAWithCheck, shrAtb, shrLtmp, shrtmp, shrmapping);
                      BackwardSubstitutionUnrolled<DataType, 6>::compute(
                        shrLtmp, shrtmp, shrstmp);
                      break;
                  case 7:
                      FusedCholeskyForwardSubstUnrolled<DataType, 7>::compute(
                        shrAtAWithCheck, shrAtb, shrLtmp, shrtmp, shrmapping);
                      BackwardSubstitutionUnrolled<DataType, 7>::compute(
                        shrLtmp, shrtmp, shrstmp);
                      break;
                  case 8:
                      FusedCholeskyForwardSubstUnrolled<DataType, 8>::compute(
                        shrAtAWithCheck, shrAtb, shrLtmp, shrtmp, shrmapping);
                      BackwardSubstitutionUnrolled<DataType, 8>::compute(
                        shrLtmp, shrtmp, shrstmp);
                      break;
                  default:
                      FusedCholeskyForwardSubst<DataType>::compute(
                        shrAtAWithCheck, shrAtb, shrLtmp, shrtmp, shrmapping, nPassive);
                      BackwardSubstitution<DataType>::compute(
                        shrLtmp, shrtmp, shrstmp, nPassive);
                  }
                
                  bool hasNegative = false;
                  // TODO: how to go around instabilities in cholesky+solvesr?
                  // 1) there are situations where cholesky decomposition + solvers yield nans
                  //   nans are treated by checking for them and breaking out
                  // 2) there are situations, in particularly for the Cholesky decomp,
                  //   when computing L_i_i = std::sqrt(M_ii - Sum) is non-deterministic
                  //   and not stable. This is different from plain nans as there is a result
                  //   value that appears to be normal, but differes from cpu 
                  bool hasNans = false;
                  for (int ii=0; ii<nPassive; ++ii) {
                      auto const s_ii = shrstmp(ii);
                      hasNegative |= s_ii <= 0;
                      hasNans |= isnan(s_ii);
                  }
                  if (hasNans)
                      break;
                  if (!hasNegative) {
                      for (int i=0; i<nPassive; ++i) {
                          // note, s contains passive/active set layout
                          // and x contains unpermuted final values in their repective pos
                          auto const real_i = shrmapping(i);
                          shrx(real_i) = shrstmp(i);
                      }
                      break;
                  }

                  auto alpha = std::numeric_limits<DataType>::max();
                  char alpha_idx=0, real_alpha_idx=0;

                  for (auto i = 0; i < nPassive; ++i) {
                    auto const s_i = shrstmp(i);
                    if (s_i <= 0.) {
                      auto const real_i = shrmapping(i);
                      auto const x_i = shrx(real_i);
                      auto const ratio = x_i / (x_i - s_i);
                      if (ratio < alpha) {
                        alpha = ratio;
                        alpha_idx = i;
                        real_alpha_idx = real_i;
                      }
                    }
                  }

                  for (int ii=0; ii<nPassive; ++ii) {
                    auto const real_i = shrmapping(ii);
                    auto const x_ii = shrx(real_i);
                    shrx(real_i) += alpha * (shrstmp(ii) - x_ii);
                  }
                  shrx(real_alpha_idx) = 0;
                  --nPassive;

                  Eigen::numext::swap(
                        shrmapping(nPassive), 
                        shrmapping(alpha_idx));
                }

                // this is as it was in cpu version
                // note: iter was incremented first and then the check was done
                if ((iter+1) % 16 == 0)
                    eps_to_use *= 2;
            }
        }
        __syncthreads();

        //
        // compute chi2
        //
        if (ty==0) {
            DataType sum{0};
            auto const s_sample = shrs(tx);
            #pragma unroll
            for (int i=0; i<nsamples; ++i)
                sum += shrP(tx, i) * shrx(i);
            auto const v_sample = sum - s_sample;
            shrv(tx) = v_sample;
        }
        __syncthreads();

        // use single thread in here for now
        // FIXME: should we add?
        if (ty==0 && tx==0) {
            // f/b subst solver
            ForwardSubstitutionUnrolled<DataType, nsamples>::compute(
                shrL, shrv, shrtmp);

            // sum will be the chi2 
            DataType chi2_new{0};
            auto const chi2_old =  shr_current_chi2;
            #pragma unroll
            for (int i=0; i<nsamples; i++)
                chi2_new += shrtmp(i)*shrtmp(i);
            shr_current_chi2 = chi2_new;
            auto const chi2_diff = chi2_new - chi2_old;

            // if this channel finished - set the flag
            if (ecal::abs(chi2_diff) < 1e-3) {
                shrFinished=true;
            }
        }
        __syncthreads();

        // manage state for all threads per channel
        // note, non-divergent branch per block 
        // all threads per channel will follow or not
        if (shrFinished) {
            // store to global once finished minimization
            if (ty==0)
                x(tx) = shrx(tx);
            if (ty==1 && tx==0)
                g_chi2[grch] = shr_current_chi2;
            if (ty==1 && tx==1)
                energies[grch] = shrx(5); // TODO remove hardcoded values
            return;
        }
    }
}

/// launch ctx:
/// 1 thread per channel
// TODO: add __restrict__ 
// TODO: do we need shared mem in here?
__global__
void kernel_fnnls(
        SampleVector::Scalar const* g_AtA,
        SampleVector::Scalar const* g_Atb,
        SampleVector::Scalar *g_L,
        SampleVector *xs,
        char *g_mapping,
        char *g_npassive,
        uint32_t const* v2ridmapping,
        uint32_t const nchannels,
        uint32_t const ml_iteration) {
    // 
    using DataType = SampleVector::Scalar;
    using ConstDataType = DataType const;
    constexpr double eps = 1e-11;
    constexpr unsigned int max_iterations = 500;
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    constexpr auto nvaluesForL = MapSymM<DataType, nsamples>::total;
    constexpr auto nvaluesForAtA = nvaluesForL;

    // indices
    auto const vch = threadIdx.x + blockIdx.x * blockDim.x;
    if (vch >= nchannels) return;
    auto const rch = v2ridmapping[vch];

    // non-divergent branch!
    // only for the first iteration...
    if (ml_iteration == 0) {
        // we need to filter out channels that either
        // - have values that are precomputed
        // - noise matrix is 0
        if (rch==idIsInvalidMask)
            return;
    }

    // map global mem
    MapSymMWithCheck<ConstDataType, nsamples, Eigen::ColMajor> 
    AtA{g_AtA + nvaluesForAtA*rch};
    MapSymM<DataType, nsamples, Eigen::RowMajor> L{g_L + nvaluesForL*rch};
    MapV<ConstDataType> Atb{g_Atb + nsamples*rch};
    MapV<DataType> x{xs[rch].data()};
    MapV<char> mapping{g_mapping + rch*nsamples};

    // load number of elements in the passive set
    auto nPassive = g_npassive[rch];

    // used to break from the loop
    unsigned int w_max_idx_prev = 0;
    DataType w_max_prev = 0;
    double eps_to_use = eps;

    // main loop
    for (int iter = 0; iter < max_iterations; ++iter) {
        if (iter>0 || nPassive==0) {
            const auto nActive = nsamples - nPassive;

            if(!nActive)
              break;

            //  
            unsigned int w_max_idx = -1;
            auto max_w {static_cast<DataType>(-1)};
            for (unsigned int i=nsamples-nActive; i<nsamples; i++) {
                auto sum_per_row{static_cast<DataType>(0)};
                auto const real_i = mapping(i);
                auto const atb = Atb(real_i);
                #pragma unroll
                for (unsigned int k=0; k<nsamples; ++k)
                    // note, we do not need to look up k in the mapping
                    // both AtA and x have swaps applied -> therefore dot product will 
                    // not change per row
                    sum_per_row += AtA(real_i, k) * x(k);

                // compute gradient value and check if it is greater than max
                auto const wvalue = atb - sum_per_row;
                if (max_w < wvalue) {
                    max_w = wvalue;
                    w_max_idx = i;
                }
            }

            // check for convergence
            if (max_w < eps_to_use || 
                (w_max_idx==w_max_idx_prev && max_w==w_max_prev))
              break;

            // update
            w_max_prev = max_w;
            w_max_idx_prev = w_max_idx;

            // swap samples to maintain passive vs active set
            Eigen::numext::swap(
                    mapping(nPassive), 
                    mapping(w_max_idx));
            ++nPassive;
        }

        // inner loop
        DataType __s[nsamples], __tmp[nsamples];
        MapV<DataType> s{__s}, tmp{__tmp};
        while (nPassive > 0) {
          switch (nPassive) {
          case 1:
              FusedCholeskySolver<DataType, 1>::compute(AtA, Atb, s, mapping);
              break;
          case 2:
              FusedCholeskySolver<DataType, 2>::compute(AtA, Atb, s, mapping);
              break;
          case 3:
              FusedCholeskySolver<DataType, 3>::compute(AtA, Atb, s, mapping);
              break;
          case 4:
              FusedCholeskyForwardSubstUnrolled<DataType, 4>::compute(AtA, Atb, L, tmp, 
                mapping);
              BackwardSubstitutionUnrolled<DataType, 4>::compute(L, tmp, s);
              break;
          case 5:
              FusedCholeskyForwardSubstUnrolled<DataType, 5>::compute(AtA, Atb, L, tmp,
                mapping);
              BackwardSubstitutionUnrolled<DataType, 5>::compute(L, tmp, s);
              break;
          case 6:
              FusedCholeskyForwardSubstUnrolled<DataType, 6>::compute(AtA, Atb, L, tmp,
                mapping);
              BackwardSubstitutionUnrolled<DataType, 6>::compute(L, tmp, s);
              break;
          case 7:
              FusedCholeskyForwardSubstUnrolled<DataType, 7>::compute(AtA, Atb, L, tmp,
                mapping);
              BackwardSubstitutionUnrolled<DataType, 7>::compute(L, tmp, s);
              break;
          case 8:
              FusedCholeskyForwardSubstUnrolled<DataType, 8>::compute(AtA, Atb, L, tmp,
                mapping);
              BackwardSubstitutionUnrolled<DataType, 8>::compute(L, tmp, s);
              break;
          default:
              FusedCholeskyForwardSubst<DataType>::compute(AtA, Atb, L, tmp,
                mapping, nPassive);
              BackwardSubstitution<DataType>::compute(L, tmp, s, nPassive);
          }
        
          bool hasNegative = false;
          // TODO: how to go around instabilities in cholesky+solvesr?
          // 1) there are situations where cholesky decomposition + solvers yield nans
          //   nans are treated by checking for them and breaking out
          // 2) there are situations, in particularly for the Cholesky decomp,
          //   when computing L_i_i = std::sqrt(M_ii - Sum) is non-deterministic
          //   and not stable. This is different from plain nans as there is a result
          //   value that appears to be normal, but differes from cpu 
          bool hasNans = false;
          for (int ii=0; ii<nPassive; ++ii) {
              auto const s_ii = s(ii);
              hasNegative |= s_ii <= 0;
              hasNans |= isnan(s_ii);
          }
          if (hasNans)
              break;
          if (!hasNegative) {
              for (int i=0; i<nPassive; ++i) {
                  // note, s contains passive/active set layout
                  // and x contains unpermuted final values in their repective pos
                  auto const real_i = mapping(i);
                  x(real_i) = s(i);
              }
              break;
          }

          auto alpha = std::numeric_limits<DataType>::max();
          char alpha_idx=0, real_alpha_idx=0;

          for (auto i = 0; i < nPassive; ++i) {
            auto const s_i = s(i);
            if (s_i <= 0.) {
              auto const real_i = mapping(i);
              auto const x_i = x(real_i);
              auto const ratio = x_i / (x_i - s_i);
              if (ratio < alpha) {
                alpha = ratio;
                alpha_idx = i;
                real_alpha_idx = real_i;
              }
            }
          }

          for (int ii=0; ii<nPassive; ++ii) {
            auto const real_i = mapping(ii);
            auto const x_ii = x(real_i);
            x(real_i) += alpha * (s(ii) - x_ii);
          }
          x(real_alpha_idx) = 0;
          --nPassive;

          Eigen::numext::swap(
                mapping(nPassive), 
                mapping(alpha_idx));
        }

        // this is as it was in cpu version
        // note: iter was incremented first and then the check was done
        if ((iter+1) % 16 == 0)
            eps_to_use *= 2;
    }

    // store back to global
    g_npassive[rch] = nPassive;
}

// TODO: add __restrict__
// launch ctx: 10 threads per channel x nchannels per block
__global__
void kernel_compute_chi2(
        SampleVector::Scalar const* g_L,
        EcalPulseShape const* g_pulse_shape,
        SampleVector const* g_x,
        ::ecal::reco::StorageScalarType *energies,
        SampleVector const* g_s,
        ::ecal::reco::StorageScalarType* g_chi2,
        ::ecal::reco::StorageScalarType* g_chi2_prev,
        uint32_t const* input_v2ridmapping,
        uint32_t *output_v2ridmapping,
        uint32_t *dids,
        uint32_t const offsetForHashes,
        uint32_t *pChannelsLeft,
        uint32_t const nchannels,
        uint32_t const iteration) {
    // constants, typedefs
    using DataType = SampleVector::Scalar;
    using ConstDataType = DataType const;
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    constexpr auto nvaluesForL = MapSymM<DataType, nsamples>::total;
    constexpr auto nvaluesForPulseShape=EcalPulseShape::TEMPLATESAMPLES;

    // indices
    auto const gtx = threadIdx.x + blockDim.x * blockIdx.x;
    auto const gvch = gtx / nsamples;
    if (gvch >= nchannels) return;
   
    auto const grch = input_v2ridmapping[gvch];
    auto const lch = threadIdx.x  / nsamples;
    auto const sample = threadIdx.x % nsamples;

    // non-divergent branch!
    // only for the first iteration...
    if (iteration == 0) {
        // we need to filter out channels that either
        // - have values that are precomputed
        // - noise matrix is 0
        if (grch==idIsInvalidMask)
            return;
    }
    
    // get did and hashedId
    auto const did = DetId{dids[grch]};
    auto const isBarrel = did.subdetId() == EcalBarrel;
    auto const hashedId = isBarrel
        ? hashedIndexEB(did.rawId())
        : offsetForHashes + hashedIndexEE(did.rawId());

    // shared mem
    extern __shared__ char smem[];
    DataType* __shrL = reinterpret_cast<DataType*>(smem);
    float* __shrP = __shrL + nvaluesForL * blockDim.x / nsamples;
    DataType *__shrv = __shrP + nvaluesForPulseShape*blockDim.x/nsamples;
    DataType *__shrtmpx = __shrv + blockDim.x;

    // map global mem
    MapSymM<ConstDataType, nsamples, Eigen::RowMajor> L
    {g_L + nvaluesForL * grch};
    MapMForPM<float> pulse_shape{(float*)&g_pulse_shape[hashedId]};
    MapV<ConstDataType> x{g_x[grch].data()};
    MapV<ConstDataType> s{g_s[grch].data()};

    // map shared mem
    MapSymM<DataType, nsamples, Eigen::RowMajor> shrL
    {__shrL + nvaluesForL*lch};
    MapMForPM<float> shrP{__shrP + nvaluesForPulseShape*lch};
    MapV<DataType> shrv{__shrv + nsamples*lch};
    MapV<DataType> shrtmpx{__shrtmpx + nsamples*lch};

    // store to shared mem
    shrtmpx(sample) = x(sample);
    if (sample == 5)
        energies[grch] = shrtmpx(sample); // store to global soi amplitude
    // use plain buffers (global mem) for the pulse shape
    shrP.data[sample] = pulse_shape.data[sample];
    if (sample<2)
        shrP.data[sample + 10] = pulse_shape.data[sample + 10];
    #pragma unroll
    for (int i=0; i<nsamples; i++) {
        if (sample >= i)
            shrL(sample, i) = L(sample, i);
    }
    __syncthreads();

    // compute P*x - s
    DataType sum{0};
    auto const s_sample = s(sample);
    #pragma unroll
    for (int i=0; i<nsamples; ++i)
        sum += shrP(sample, i) * shrtmpx(i);
    auto const v_sample = sum - s_sample;
    shrv(sample) = v_sample;
    __syncthreads();


    // use single thread in here for now
    // FIXME: should we add?
    if (sample == 0) {
        // f/b subst solver
        ForwardSubstitutionUnrolled<DataType, nsamples>::compute(
            shrL, shrv, shrtmpx);

        // sum will be the chi2 
        DataType chi2_new{0};
        auto const chi2_old =  g_chi2[grch];
        #pragma unroll
        for (int i=0; i<nsamples; i++)
            chi2_new += shrtmpx(i)*shrtmpx(i);
        g_chi2[grch] = chi2_new;
        auto const chi2_diff = chi2_new - chi2_old;

        // state management logic
        // if this ch needs to continue, inc the counter atomically
        if (ecal::abs(chi2_diff) >= 1e-3) {
            // FIXME: remove hardcodded values
            if (iteration>=9) {
                auto const chi2_prev = g_chi2_prev[grch];
                // TODO debugging
                printf("chid = %d chi2_new = %f chi2_old = %f chi2_prev = %f\n",
                    grch, chi2_new, chi2_old, chi2_prev);
                for (int i=0; i<10;i ++)
                    printf("x(%d) = %f\n", i, x(i));
                // if there is a rotation of chi2 values stop here
                //if (ecal::abs(chi2_new - chi2_prev) / chi2_prev < 1e-3)
                //    return;
            }

            g_chi2_prev[grch] = chi2_old;
            auto const pos = atomicAdd(pChannelsLeft, 1);
            output_v2ridmapping[pos] = grch;
        }
    }
}

/*
__global__
void kernel_minimization_launcher() {
    // we have a single guy in here
    if (threadIdx.x >= 1) return;

    int iter = 0;
    auto channelsLeft = totalChannels;
    while (true) {
        if (iter >= max_iterations)
            break;

        // 
        kernel_update_covariance_compute_cholesky();

        // forward / backward substitution
        kernel_solve();

        // matrix matrix and matrix vector mults
        kernel_mm_mv_mults();

        // fnnls 
        kernel_fnnls();

        // compute chi2
        kernel_compute_chi2();
        cudaDeviceSynchronize();

        // check exit
        channelsLeft = *channelsCounter;
        if (channelsLeft == 0)
            return;
    }
}
*/

void minimization_procedure_fused(
        EventInputDataCPU const& eventInputCPU, EventInputDataGPU& eventInputGPU,
        EventOutputDataGPU& eventOutputGPU, EventDataForScratchGPU& scratch,
        ConditionsProducts const& conditions,
        ConfigurationParameters const& configParameters,
        cuda::stream_t<>& cudaStream,
        unsigned int offsetForHashes) {
    // constants
    unsigned int totalChannels = eventInputCPU.ebDigis.size() 
        + eventInputCPU.eeDigis.size();
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    using DataType = SampleVector::Scalar;

    // launch 
    dim3 const ts{nsamples, nsamples, 1};
    uint32_t const bs = totalChannels;
    kernel_minimize_fused<<<bs, ts, 0, cudaStream.id()>>>(
        conditions.pulseCovariances.values,
        conditions.pulseShapes.values,
        (SampleVector*)eventOutputGPU.amplitudesAll,
        scratch.samples,
        scratch.samplesMapping,
        scratch.npassive,
        eventOutputGPU.amplitude,
        eventOutputGPU.chi2,
        scratch.v2rmapping_1,
        eventInputGPU.ids,
        scratch.gainsNoise,
        conditions.pedestals.rms_x12,
        conditions.pedestals.rms_x6,
        conditions.pedestals.rms_x1,
        conditions.gainRatios.gain12Over6,
        conditions.gainRatios.gain6Over1,
        conditions.samplesCorrelation.EBG12SamplesCorrelation,
        conditions.samplesCorrelation.EBG6SamplesCorrelation,
        conditions.samplesCorrelation.EBG1SamplesCorrelation,
        conditions.samplesCorrelation.EEG12SamplesCorrelation,
        conditions.samplesCorrelation.EEG6SamplesCorrelation,
        conditions.samplesCorrelation.EEG1SamplesCorrelation,
        scratch.hasSwitchToGain6,
        scratch.hasSwitchToGain1,
        scratch.isSaturated,
        offsetForHashes,
        0,
        totalChannels
    );
    cudaCheck( cudaGetLastError() );
}

void minimization_procedure_splitted_host_launch(
        EventInputDataCPU const& eventInputCPU, EventInputDataGPU& eventInputGPU,
        EventOutputDataGPU& eventOutputGPU, EventDataForScratchGPU& scratch,
        ConditionsProducts const& conditions,
        ConfigurationParameters const& configParameters,
        cuda::stream_t<>& cudaStream,
        unsigned int offsetForHashes) {
    unsigned int totalChannels = eventInputCPU.ebDigis.size() 
        + eventInputCPU.eeDigis.size();
    uint32_t nchannels = totalChannels;
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    using DataType = SampleVector::Scalar;

    // FIXME: all the constants below need to be propagated properly
    // once results are valid
    for (int iteration=0; iteration<50; ++iteration) {
        std::cout << "iteration = " << iteration 
                  << "  nchannels = " << nchannels
                  << std::endl;

        //
        dim3 const ts1{nsamples, nsamples, nchannels>10 ? 10 : nchannels};
        uint32_t const bs1 = nchannels>10
            ? (nchannels + 10 - 1) / 10
            : 1;
        uint32_t const shrBytes1 = 10 * (sizeof(DataType)*(55 + 55 + 10));
        kernel_newiter_update_covariance_compute_cholesky<<<bs1, ts1, shrBytes1, cudaStream.id()>>>(
            conditions.pulseCovariances.values,
            (SampleVector const*)eventOutputGPU.amplitudesAll,
            scratch.decompMatrixMainLoop,
            iteration%2==0 ? scratch.v2rmapping_1 : scratch.v2rmapping_2,
            eventInputGPU.ids,
            scratch.pChannelsCounter,
            scratch.gainsNoise,
            conditions.pedestals.rms_x12,
            conditions.pedestals.rms_x6,
            conditions.pedestals.rms_x1,
            conditions.gainRatios.gain12Over6,
            conditions.gainRatios.gain6Over1,
            conditions.samplesCorrelation.EBG12SamplesCorrelation,
            conditions.samplesCorrelation.EBG6SamplesCorrelation,
            conditions.samplesCorrelation.EBG1SamplesCorrelation,
            conditions.samplesCorrelation.EEG12SamplesCorrelation,
            conditions.samplesCorrelation.EEG6SamplesCorrelation,
            conditions.samplesCorrelation.EEG1SamplesCorrelation,
            scratch.hasSwitchToGain6,
            scratch.hasSwitchToGain1,
            scratch.isSaturated,
            offsetForHashes,
            iteration,
            nchannels);
        cudaCheck( cudaGetLastError() );

        //
        auto const ts2 = ts1;
        auto const bs2 = bs1;
        uint32_t const shrBytes2 = 10 * (sizeof(DataType)*(55 + 10 + 100 + 10) + 
            sizeof(float)*EcalPulseShape::TEMPLATESAMPLES);
        kernel_solve_mm_mv_mults<<<bs2, ts2, shrBytes2, cudaStream.id()>>>(
            conditions.pulseShapes.values,
            scratch.decompMatrixMainLoop,
            scratch.samples,
            scratch.AtA,
            scratch.Atb,
            iteration%2==0 ? scratch.v2rmapping_1 : scratch.v2rmapping_2,
            eventInputGPU.ids,
            offsetForHashes,
            iteration,
            nchannels);
        cudaCheck( cudaGetLastError() );

        //
        // FIXME: rename kernelMinimizeThreads
        uint32_t ts3 = nchannels>configParameters.kernelMinimizeThreads[0]
            ? configParameters.kernelMinimizeThreads[0]
            : nchannels;
        uint32_t bs3 = nchannels>ts3
            ? (nchannels + ts3 - 1) / ts3
            : 1;
        kernel_fnnls<<<bs3, ts3, 0, cudaStream.id()>>>(
            scratch.AtA,
            scratch.Atb,
            scratch.decompMatrixFnnls,
            (SampleVector*)eventOutputGPU.amplitudesAll,
            scratch.samplesMapping,
            scratch.npassive,
            iteration%2==0 ? scratch.v2rmapping_1 : scratch.v2rmapping_2,
            nchannels,
            iteration);
        cudaCheck( cudaGetLastError() );

        //
        uint32_t const ts4 = nchannels>32 ? 10*32 : nchannels*10;
        uint32_t const bs4 = nchannels>ts4
            ? (nchannels*10 + ts4 - 1) / ts4
            : 1;
        uint32_t const shrBytes4 = 32 * (sizeof(DataType)*(55 + 10 + 10)
            + sizeof(float)*EcalPulseShape::TEMPLATESAMPLES);
        kernel_compute_chi2<<<bs4, ts4, shrBytes4, cudaStream.id()>>>(
            scratch.decompMatrixMainLoop,
            conditions.pulseShapes.values,
            (SampleVector const*)eventOutputGPU.amplitudesAll,
            eventOutputGPU.amplitude,
            scratch.samples,
            eventOutputGPU.chi2,
            scratch.chi2_prev,
            iteration%2==0 ? scratch.v2rmapping_1 : scratch.v2rmapping_2,
            iteration%2==0 ? scratch.v2rmapping_2 : scratch.v2rmapping_1,
            eventInputGPU.ids,
            offsetForHashes,
            scratch.pChannelsCounter,
            nchannels,
            iteration);
        cudaCheck( cudaGetLastError() );

        // 
        cudaCheck( cudaMemcpyAsync(&nchannels, scratch.pChannelsCounter, 
            sizeof(uint32_t), cudaMemcpyDeviceToHost, cudaStream.id()) );
        cudaCheck( cudaStreamSynchronize(cudaStream.id()) );
        
        // 
        if (nchannels == 0)
            return;
    }
}

}}
