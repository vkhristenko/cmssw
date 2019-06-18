#include <iostream>
#include <limits>

#include "cuda.h"

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "DataFormats/Math/interface/approx_exp.h"
#include "DataFormats/Math/interface/approx_log.h"

#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "inplace_fnnls.h"
#include "AmplitudeComputationKernels.h"
#include "AmplitudeComputationCommonKernels.h"
#include "KernelHelpers.h"

namespace ecal { namespace multifit { 

#define PRINT_MATRIX_10x10(M)\
            printf("%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f %f\n", \
                M(0, 0), M(1, 0), M(2, 0), M(3, 0), M(4, 0), \
                M(5, 0), M(6, 0), M(7, 0), M(8, 0), M(9, 0), \
                M(0, 1), M(1, 1), M(2, 1), M(3, 1), M(4, 1), \
                M(5, 1), M(6, 1), M(7, 1), M(8, 1), M(9, 1), \
                M(0, 2), M(1, 2), M(2, 2), M(3, 2), M(4, 2), \
                M(5, 2), M(6, 2), M(7, 2), M(8, 2), M(9, 2), \
                M(0, 3), M(1, 3), M(2, 3), M(3, 3), M(4, 3), \
                M(5, 3), M(6, 3), M(7, 3), M(8, 3), M(9, 3), \
                M(0, 4), M(1, 4), M(2, 4), M(3, 4), M(4, 4), \
                M(5, 4), M(6, 4), M(7, 4), M(8, 4), M(9, 4), \
                M(0, 5), M(1, 5), M(2, 5), M(3, 5), M(4, 5), \
                M(5, 5), M(6, 5), M(7, 5), M(8, 5), M(9, 5), \
                M(0, 6), M(1, 6), M(2, 6), M(3, 6), M(4, 6), \
                M(5, 6), M(6, 6), M(7, 6), M(8, 6), M(9, 6), \
                M(0, 7), M(1, 7), M(2, 7), M(3, 7), M(4, 7), \
                M(5, 7), M(6, 7), M(7, 7), M(8, 7), M(9, 7), \
                M(0, 8), M(1, 8), M(2, 8), M(3, 8), M(4, 8), \
                M(5, 8), M(6, 8), M(7, 8), M(8, 8), M(9, 8), \
                M(0, 9), M(1, 9), M(2, 9), M(3, 9), M(4, 9), \
                M(5, 9), M(6, 9), M(7, 9), M(8, 9), M(9, 9) \
            )

__device__
__forceinline__
bool update_covariance(SampleMatrix const& noisecov,
                       EcalPulseCovariance const& pulse_cov,
                       SampleMatrix& inverse_cov,
                       BXVectorType const& bxs,
                       SampleDecompLLT& covariance_decomposition,
                       SampleVector const& amplitudes) {
    constexpr int nsamples = SampleVector::RowsAtCompileTime;
    constexpr int npulses = BXVectorType::RowsAtCompileTime;

    inverse_cov = noisecov;

    for (unsigned int ipulse=0; ipulse<npulses; ipulse++) {
        if (amplitudes.coeff(ipulse) == 0) 
            continue;

        int bx = bxs.coeff(ipulse);
        int first_sample_t = std::max(0, bx+3);
        int offset = -3 - bx;

        auto const value = amplitudes.coeff(ipulse);
        auto const value_sq = value*value;

        unsigned int nsample_pulse = nsamples - first_sample_t;
        for (int i=first_sample_t; i<first_sample_t+nsample_pulse; ++i)
            for (int j=first_sample_t; j<first_sample_t+nsample_pulse; ++j)
                inverse_cov(i, j) += value_sq * pulse_cov.covval[i+offset][j+offset];

        /*
        inverse_cov.block(first_sample_t, first_sample_t, 
                          nsample_pulse, nsample_pulse)
            += value_sq * full_pulse_cov.block(first_sample_t + offset,
                                               first_sample_t + offset,
                                               nsample_pulse,
                                               nsample_pulse);
                                               */
    }

    return true;
}

///
/// The following conventions apply to all of the kernels below
///     - virtual channel id - channel id for the current iteration
///     - real  channel id - channel id for the original array
/// 

/// launch ctx
// FIXME: add __restrict__ whenever is needed
// TODO: add boundaries checking for threads
__global__
void kernel_update_covariance_compute_cholesky(
        EcalPulseCovariance const* g_PulseCovariance,
        SampleMatrix const* g_noiseCovariance,
        SampleVector const* g_amplitudes,
        SampleVector::Scalar* g_L,
        uint32_t *v2ridmapping,
        uint32_t *dids,
        uint32_t offsetForHashes) {
    // constants
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    constexpr auto npulses = nsamples;
    constexpr auto template_samples = EcalPulseShape::TEMPLATESAMPLES;
    using DataType = SampleVector::Scalar;
    using ConstDataType = DataType const;
    constexpr auto nvaluesForL = MapSymM<DataType, nsamples>::total;

    // indices
    auto const tx = threadIdx.x;
    auto const ty = threadIdx.y;
    auto const vch_per_block = threadIdx.z;
    auto const gvch = vch_per_block + blockIdx.x*blockDim.z;
    auto const grch = v2ridmapping[gvch];
    auto const did = DetId{dids[grch]};
    auto const isBarrel = did.subdetId() == EcalBarrel;
    auto const hashedId = isBarrel
        ? hashedIndexEB(did.rawId())
        : offsetForHashes + hashedIndexEE(did.rawId());

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
    // TODO: do we need 10 x 10 here??? oro 55 is enough???
    DataType* __shrLSums = __shrL + 
        blockDim.z * nvaluesForL;

    // map global mem
    MapM<ConstDataType, template_samples, Eigen::RowMajor> PulseCovariance
    {reinterpret_cast<DataType const*>(g_PulseCovariance + hashedId)};
    MapM<ConstDataType, nsamples> noiseCovariance
    {g_noiseCovariance[grch].data()};
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
    auto noiseValue = noiseCovariance(ty, tx);
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
        auto const offset = 3 - bx;
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

    // TODO: verify that the loop is unrolled
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

        if (tx>=column && ty>=tx)
            shrLSums(ty, tx) += shrL(ty, column-1) * shrL(tx, column-1);
        __syncthreads();
    }

    // store back to global
    L(ty, tx) = shrL(ty, tx);
}

// FIXME: add __restrict__ whenever is needed
// TODO: add boundaries checking for threads
__global__
void kernel_solve_mm_mv_mults(
        PulseMatrixType const* g_pulseMatrix,
        SampleVector::Scalar const* g_L,
        SampleVector const* g_s,
        SampleVector::Scalar* g_AtA,
        SampleVector::Scalar* g_Atb,
        uint32_t const* v2ridmapping) {
    // constants, typedefs
    using DataType = SampleVector::Scalar;
    using ConstDataType = DataType const;
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    constexpr auto nvaluesForL = MapSymM<DataType, nsamples>::total;

    // indices
    auto const tx = threadIdx.x;
    auto const ty = threadIdx.y;
    auto const lch = threadIdx.z;
    auto const gvch = lch + blockIdx.x * blockDim.z;
    auto const grch = v2ridmapping[gvch];
//    auto const did = DetId{dids[grch]};

    // configure shared mem
    extern __shared__ char smem[];
    DataType* __shrL = reinterpret_cast<DataType*>(smem);
    DataType* __shrP = __shrL + nvaluesForL * blockDim.z;
    DataType* __shrs = __shrP + nsamples * nsamples * blockDim.z;
    DataType* __shrAtmp = __shrs + nsamples * blockDim.z;
    DataType* __shrbtmp = __shrAtmp + nsamples * nsamples * blockDim.z;

    // map global mem
    // NOTE: this guy needs to be in sync with what is used above/below
    MapSymM<ConstDataType, nsamples, Eigen::RowMajor> L
    {g_L + grch*nvaluesForL};
    // FIXME: note, we do not need to have a matrix in global mem
    MapM<ConstDataType, nsamples> P{g_pulseMatrix[grch].data()};
    MapV<ConstDataType> s{g_s[grch].data()};
    // FIXME: we should use symmetric matrix in here
    MapM<DataType, nsamples, Eigen::RowMajor> AtA
    {g_AtA + grch * nsamples*nsamples};
    MapV<DataType> Atb{g_Atb + grch * nsamples};

    // map shared mem
    MapSymM<DataType, nsamples, Eigen::RowMajor> shrL
    {__shrL + lch * nvaluesForL};
    MapM<DataType, nsamples, Eigen::RowMajor> shrP
    {__shrP + lch*nsamples*nsamples};
    MapV<DataType> shrs{__shrs + lch*nsamples};
    MapM<DataType, nsamples, Eigen::RowMajor> shrAtmp
    {__shrAtmp + lch*nsamples*nsamples};
    MapV<DataType> shrbtmp{__shrbtmp + lch*nsamples};

    // move into shared mem
    shrP(ty, tx) = P(ty, tx);
    if (ty>=tx)
        shrL(ty, tx) = L(ty, tx);
    if (ty == 0)
        shrs(tx) = s(tx);
    __syncthreads();

    // run forward and backward substitution
    // use 10 thraeds for matrix and 1 thread vector
    // TODO: would it be better with more threads + synching???
    if (ty==0) {
        ForwardSubstitutionUnrolled<DataType, nsamples>::compute(
            shrL, shrP, shrAtmp, tx);
        // note, we are reusing shrP as the result space
        // for backward substitution, there are races!
        BackwardSubstitutionUnrolled<DataType, nsamples>::compute(
            shrL, shrAtmp, shrP, tx);
    } else if (ty==1 && tx == 1) {
        ForwardSubstitutionUnrolled<DataType, nsamples>::compute(
            shrL, shrs, shrbtmp);
        BackwardSubstitutionUnrolled<DataType, nsamples>::compute(
            shrL, shrbtmp, shrs);
    }
    __syncthreads();

    // matrix matrix and matrix vector mult
    // shrP is matrix A
    // shrs is vector b
    DataType ata_i_j = 0;
    #pragma unroll
    for (int i=0; i<nsamples; ++i) {
        auto const A_i_ty = shrP(i, ty);
        ata_i_j += shrP(i, ty) * shrP(i, tx);
    }

    // store back to global
    AtA(ty, tx) = ata_i_j;
    if (ty == 0) {
        DataType sum = 0;
        #pragma unroll
        for (int i=0; i<nsamples; i++)
            sum += shrP(i, tx) * shrs(i);

        // store back to global
        Atb(tx) = sum;
    }
}

// TODO: add __restrict__ 
__global__
void kernel_fnnls(
        SampleVector::Scalar const* g_AtA,
        SampleVector::Scalar const* g_Atb,
        SampleVector::Scalar *g_L,
        SampleVector *xs,
        char *g_mapping,
        char *g_npassive,
        uint32_t const* v2ridmapping,
        uint32_t const nchannels) {
    // 
    using DataType = SampleVector::Scalar;
    using ConstDataType = DataType const;
    constexpr double eps = 1e-11;
    constexpr unsigned int max_iterations = 1000;
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    constexpr auto nvaluesForL = MapSymM<DataType, nsamples>::total;

    // indices
    auto const vch = threadIdx.x + blockIdx.x * blockDim.x;
    auto const rch = v2ridmapping[vch];

    if (vch >= nchannels) return;

    // map global mem
    MapM<ConstDataType, nsamples, Eigen::RowMajor> AtA{g_AtA + nsamples*nsamples*rch};
    MapSymM<DataType, nsamples, Eigen::RowMajor> L{g_L + nvaluesForL*rch};
    MapV<ConstDataType> Atb{g_Atb + nsamples*rch};
    MapV<DataType> x{xs[rch].data()};
    MapV<char> mapping{g_mapping + rch*nsamples};

    // load number of elements in the passive set
    auto nPassive = g_npassive[rch];

    // main loop
    for (auto iter = 0; iter < max_iterations; ++iter) {
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
        if (max_w < eps)
          break;

        Eigen::numext::swap(
                mapping(nPassive), 
                mapping(w_max_idx));
        ++nPassive;

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
          for (int ii=0; ii<nPassive; ++ii) {
              hasNegative |= s(ii) <= 0;
          }
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
            if (s(i) <= 0.) {
              auto const real_i = mapping(i);
              auto const x_i = x(real_i);
              auto const ratio = x_i / (x_i - s(i));
              if (ratio < alpha) {
                alpha = ratio;
                alpha_idx = i;
                real_alpha_idx = real_i;
              }
            }
          }

          if (std::numeric_limits<DataType>::max() == alpha) {
            for (int i=0; i<nPassive; ++i) {
                auto const real_i = mapping(i);
                x(real_i) = s(i);
            }
            break;
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
    }

    // store back to global
    g_npassive[rch] = nPassive;
}

// TODO: add __restrict__
// launch ctx: 10 threads per channel
__global__
void kernel_compute_chi2(
        SampleVector::Scalar const* g_L,
        PulseMatrixType const* g_P,
        SampleVector const* g_x,
        SampleVector const* g_s,
        ::ecal::reco::StorageScalarType* g_chi2,
        uint32_t const* input_v2ridmapping,
        uint32_t *output_v2ridmapping,
        uint32_t *pChannelsLeft,
        uint32_t const nchannels) {
    // constants, typedefs
    using DataType = SampleVector::Scalar;
    using ConstDataType = DataType const;
    constexpr auto nsamples = SampleVector::RowsAtCompileTime;
    constexpr auto nvaluesForL = MapSymM<DataType, nsamples>::total;

    // indices
    auto const gtx = threadIdx.x + blockDim.x * blockIdx.x;
    auto const gvch = gtx / nsamples;
    auto const grch = input_v2ridmapping[gvch];
    auto const lch = threadIdx.x  / nsamples;
    auto const sample = threadIdx.x % nsamples;

    if (gvch >= nchannels) return;

    // shared mem
    extern __shared__ char smem[];
    DataType* __shrL = reinterpret_cast<DataType*>(smem);
    // FIXME: Pulse matrix should be a vector with proper indexing, not a matrix
    // nchannels per block (threads per block / 10) * values for L m
    DataType* __shrP = __shrL + nvaluesForL * blockDim.x / nsamples;
    DataType *__shrv = __shrP + nsamples*blockDim.x;
    DataType *__shrtmpx = __shrv + blockDim.x;

    // map global mem
    MapSymM<ConstDataType, nsamples, Eigen::RowMajor> L
    {g_L + nvaluesForL * grch};
    MapM<ConstDataType, nsamples> P
    {g_P[grch].data()};
    MapV<ConstDataType> x{g_x[grch].data()};
    MapV<ConstDataType> s{g_s[grch].data()};

    // map shared mem
    MapSymM<DataType, nsamples, Eigen::RowMajor> shrL
    {__shrL + nvaluesForL*lch};
    MapM<DataType, nsamples, Eigen::RowMajor> shrP
    {__shrP + nsamples*nsamples*lch};
    MapV<DataType> shrv{__shrv + nsamples*lch};
    MapV<DataType> shrtmpx{__shrtmpx + nsamples*lch};

    // mov to shared mem
    auto const x_sample = x(sample);
    auto const s_sample = s(sample);
    #pragma unroll
    for (int i=0; i<nsamples; ++i) {
        shrP(i, sample) = P(i, sample) * x_sample;
        if (sample>=i)
            shrL(sample, i) = L(sample, i);
    }
    __syncthreads();

    // prepare input for forward subst
    DataType sum{0};
    #pragma unroll
    for (int i=0; i<nsamples; ++i)
        sum += P(sample, i);
    auto const v_sample = sum - s_sample;
    shrv(sample) = v_sample;
    __syncthreads();

    // use single thread in here for now
    // FIXME: should we add?
    if (sample == 0) {
        // f/b subst solver
        ForwardSubstitutionUnrolled<DataType, nsamples>::compute(
            shrL, shrv, shrtmpx);
        BackwardSubstitutionUnrolled<DataType, nsamples>::compute(
            shrL, shrtmpx, shrv);

        // sum will be the chi2 
        DataType chi2_new{0};
        auto const chi2_old =  g_chi2[grch];
        #pragma unroll
        for (int i=0; i<nsamples; i++)
            chi2_new += shrv(i)*shrv(i);
        g_chi2[grch] = chi2_new;
        auto const chi2_diff = chi2_new - chi2_old;

        // state management logic
        // if this ch needs to continue, inc the counter atomically
        if (ecal::abs(chi2_diff) >= 1e-3) {
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

__device__
__forceinline__
SampleVector::Scalar compute_chi2(SampleDecompLLT& covariance_decomposition,
                   PulseMatrixType const& pulse_matrix,
                   SampleVector const& amplitudes,
                   SampleVector const& samples) {
    return covariance_decomposition.matrixL()
        .solve(pulse_matrix * amplitudes - samples)
        .squaredNorm();
}

///
/// launch ctx parameters are (nchannels / block, blocks)
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
                     EcalPulseCovariance const* pulse_cov,
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
                     unsigned int offsetForHashes) {
    int idx = threadIdx.x + blockDim.x*blockIdx.x;
    if (idx < nchannels) {
        // this channel had values precomputed
        if (v2rmapping[idx] == 0xffffffff)
            return;
        /*
        if (static_cast<MinimizationState>(acState[idx]) == 
            MinimizationState::Precomputed)
            return;
            */

        auto const did = DetId{dids[idx]};
        auto const isBarrel = did.subdetId() == EcalBarrel;
        auto const hashedId = isBarrel
            ? hashedIndexEB(did.rawId())
            : offsetForHashes + hashedIndexEE(did.rawId());


        // inits
        bool status = false;
        int iter = 0;
        int npassive = 0;
        
        // FIXME: 
        // need to identify this earrlier and set the state properly
        if (noiseCovIsZero[idx] == 0xffff)
            return;
        /*
        if (noisecov[idx].isZero(0))
            return;
            */

        // inits
        SampleDecompLLT covariance_decomposition;
        SampleMatrix inverse_cov;
        SampleVector::Scalar chi2 = 0, chi2_now = 0;

#ifdef ECAL_MULTIFIT_KERNEL_MINIMIZE_V1
//    PRINT_MATRIX_10x10(noisecov[idx]);
#endif

        // loop until ocnverge
        while (true) {
            if (iter >= max_iterations)
                break;

            status = update_covariance(
                noisecov[idx], 
                pulse_cov[hashedId],
                inverse_cov,
                bxs[idx],
                covariance_decomposition,
                amplitudes[idx]);

            // compute actual covariance decomposition
            covariance_decomposition.compute(inverse_cov);

            // prepare input matrices for fnnls
            SampleMatrix A = covariance_decomposition.matrixL()
                .solve(pulse_matrix[idx]);
            SampleVector b = covariance_decomposition.matrixL()
                .solve(samples[idx]);
            
            status = inplace_fnnls(
                A, b, amplitudes[idx],
                npassive, bxs[idx], pulse_matrix[idx]);
                
            chi2_now = compute_chi2(
                covariance_decomposition,
                pulse_matrix[idx],
                amplitudes[idx],
                samples[idx]);
            auto deltachi2 = chi2_now - chi2;


#ifdef ECAL_MULTIFIT_KERNEL_MINIMIZE_V1
            if (iter > 10) {
                printf("idx = %d iter = %d chi2 = %f chi2old = %f\n",
                    idx, iter, chi2_now, chi2);
                
                printf("noisecov(0, i): %f %f %f %f %f %f %f %f %f %f\n",
                    noisecov[idx](0, 0),
                    noisecov[idx](0, 1),
                    noisecov[idx](0, 2),
                    noisecov[idx](0, 3),
                    noisecov[idx](0, 4),
                    noisecov[idx](0, 5),
                    noisecov[idx](0, 6),
                    noisecov[idx](0, 7),
                    noisecov[idx](0, 8),
                    noisecov[idx](0, 9));

                printf("ampls: %f %f %f %f %f %f %f %f %f %f\n",
                    amplitudes[idx](0),
                    amplitudes[idx](1),
                    amplitudes[idx](2),
                    amplitudes[idx](3),
                    amplitudes[idx](4),
                    amplitudes[idx](5),
                    amplitudes[idx](6),
                    amplitudes[idx](7),
                    amplitudes[idx](8),
                    amplitudes[idx](9));
            }
#endif

            chi2 = chi2_now;

            if (ecal::abs(deltachi2) < 1e-3)
                break;

            //---- AM: TEST
            //---- it was 3 lines above, now here as in the CPU version
            ++iter;
            
        }

        // the rest will be set later
        chi2s[idx] = chi2;
    }
}

void minimization_procedure(
        EventInputDataCPU const& eventInputCPU, EventInputDataGPU& eventInputGPU,
        EventOutputDataGPU& eventOutputGPU, EventDataForScratchGPU& scratch,
        ConditionsProducts const& conditions,
        ConfigurationParameters const& configParameters,
        cuda::stream_t<>& cudaStream,
        unsigned int offsetForHashes) {
    unsigned int totalChannels = eventInputCPU.ebDigis.size() 
        + eventInputCPU.eeDigis.size();
//    unsigned int threads_min = conf.threads.x;
    // TODO: configure from python
    unsigned int threads_min = configParameters.kernelMinimizeThreads[0];
    unsigned int blocks_min = threads_min > totalChannels
        ? 1
        : (totalChannels + threads_min - 1) / threads_min;
    kernel_minimize<<<blocks_min, threads_min, 0, cudaStream.id()>>>(
        scratch.noisecov,
        conditions.pulseCovariances.values,
        scratch.activeBXs,
        scratch.samples,
        (SampleVector*)eventOutputGPU.amplitudesAll,
        scratch.pulse_matrix,
        eventOutputGPU.chi2,
        eventInputGPU.ids,
        scratch.v2rmapping_1,
        scratch.noiseCovIsZero,
        totalChannels,
        50,
        offsetForHashes);
    cudaCheck(cudaGetLastError());

    //
    // permute computed amplitudes
    // and assign the final uncalibared energy value
    //
    unsigned int threadsPermute = 32 * EcalDataFrame::MAXSAMPLES; // 32 * 10
    unsigned int blocksPermute = threadsPermute > 10 * totalChannels
        ? 1
        : (10 * totalChannels + threadsPermute - 1) / threadsPermute;
    int bytesPermute = threadsPermute * sizeof(SampleVector::Scalar);
    kernel_permute_results<<<blocksPermute, threadsPermute, 
                             bytesPermute, cudaStream.id()>>>(
        (SampleVector*)eventOutputGPU.amplitudesAll,
        scratch.activeBXs,
        eventOutputGPU.amplitude,
        scratch.v2rmapping_1,
        totalChannels);
    cudaCheck(cudaGetLastError());
}

}}
