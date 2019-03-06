#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgo_gpu_new.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalXtalGroupId.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovariances.h"
#include "CondFormats/EcalObjects/interface/EcalSampleMask.h"

#include <iostream>

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/inplace_fnnls.h"

#include "cuda.h"

//#define DEBUG

namespace ecal { namespace multifit { namespace v1 {

///
/// assume kernel launch configuration is 
/// (MAXSAMPLES * nchannels, blocks)
/// TODO: is there a point to split this kernel further to separate reductions
/// 
__global__
void kernel_prep_1d(EcalPulseShape const* shapes_in,
                    FullSampleVector* shapes_out, 
                    uint16_t const* digis_in,
                    SampleVector* amplitudes,
                    SampleGainVector* gainsNoise,
                    SampleGainVector* gainsPedestal,
                    float const* mean_x1,
                    float const* mean_x12,
                    float const* mean_x6,
                    float const* gain6Over1,
                    float const* gain12Over6,
                    bool* hasSwitchToGain6,
                    bool* hasSwitchToGain1,
                    bool* isSaturated,
                    float* energies,
                    bool gainSwitchUseMaxSample,
                    int nchannels) {
    constexpr bool dynamicPedestal = false;
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;
    int tx = threadIdx.x + blockIdx.x*blockDim.x;
    int nchannels_per_block = blockDim.x / nsamples;
    int total_threads = nchannels * nsamples;
    int ch = tx / nsamples;

    if (ch < nchannels) {
        // array of 10 x channels per block
        // TODO: any other way of doing simple reduction
        // assume bool is 1 byte, should be quite safe
        extern __shared__ char shared_mem[];
        bool* shr_hasSwitchToGain6 = reinterpret_cast<bool*>(
            shared_mem);
        bool* shr_hasSwitchToGain1 = shr_hasSwitchToGain6 + 
            nchannels_per_block*nsamples;
        bool* shr_hasSwitchToGain0 = shr_hasSwitchToGain1 + 
            nchannels_per_block*nsamples;
        bool* shr_isSaturated = shr_hasSwitchToGain0 + 
            nchannels_per_block*nsamples;
        char* shr_counts = reinterpret_cast<char*>(
            shr_isSaturated) + nchannels_per_block*nsamples;
        
        //
        // pulse shape template
        //
        int sample = threadIdx.x % nsamples;
        for (int isample=sample; isample<EcalPulseShape::TEMPLATESAMPLES; 
            isample+=nsamples)
            shapes_out[ch](isample + 7) = shapes_in[ch].pdfval[isample];

        //
        // amplitudes
        //
        int adc = ecalMGPA::adc(digis_in[tx]);
        int gainId = ecalMGPA::gainId(digis_in[tx]);
        float amplitude = 0.f;
        float pedestal = 0.f;
        float gainratio = 0.f;

        shr_hasSwitchToGain6[threadIdx.x] = gainId == EcalMgpaBitwiseGain6;
        shr_hasSwitchToGain1[threadIdx.x] = gainId == EcalMgpaBitwiseGain1;
        shr_hasSwitchToGain0[threadIdx.x] = gainId == EcalMgpaBitwiseGain0 ? 1 : 0;
        shr_counts[threadIdx.x] = 0;
        __syncthreads();
        
        // non-divergent branch (except for the last 4 threads)
        if (threadIdx.x<=blockDim.x-5) {
#pragma unroll
            for (int i=0; i<5; i++)
                shr_counts[threadIdx.x] += 
                    shr_hasSwitchToGain0[threadIdx.x+i];
        }
        shr_isSaturated[threadIdx.x] = shr_counts[threadIdx.x] == 5;

        //
        // unrolled reductions
        // TODO
        //
        if (sample < 5) {
            shr_hasSwitchToGain6[threadIdx.x] |= 
                shr_hasSwitchToGain6[threadIdx.x + 5];
            shr_hasSwitchToGain1[threadIdx.x] |= 
                shr_hasSwitchToGain1[threadIdx.x + 5];
        }
        __syncthreads();
        
        if (sample<2) {
            // note, both threads per channel take value [3] twice to avoid another if
            shr_hasSwitchToGain6[threadIdx.x] |= 
                shr_hasSwitchToGain6[threadIdx.x+2] | 
                shr_hasSwitchToGain6[threadIdx.x+3];
            shr_hasSwitchToGain1[threadIdx.x] |= 
                shr_hasSwitchToGain1[threadIdx.x+2] | 
                shr_hasSwitchToGain1[threadIdx.x+3];

            // sample < 2 -> first 2 threads of each channel will be used here
            // => 0 -> will compare 3 and 4
            // => 1 -> will compare 4 and 5
            shr_isSaturated[threadIdx.x+3] |= shr_isSaturated[threadIdx.x+4];
        }
        __syncthreads();

        if (sample==0) {
            shr_hasSwitchToGain6[threadIdx.x] |= shr_hasSwitchToGain6[threadIdx.x+1];
            shr_hasSwitchToGain1[threadIdx.x] |= shr_hasSwitchToGain1[threadIdx.x+1];

            hasSwitchToGain6[ch] = shr_hasSwitchToGain6[threadIdx.x];
            hasSwitchToGain1[ch] = shr_hasSwitchToGain1[threadIdx.x];

            shr_isSaturated[threadIdx.x+3] |= shr_isSaturated[threadIdx.x+4];
            isSaturated[ch] = shr_isSaturated[threadIdx.x+3];
        }

        // TODO: divergent branch
        if (gainId==0 || gainId==3) {
            pedestal = mean_x1[ch];
            gainratio = gain6Over1[ch] * gain12Over6[ch];
            gainsNoise[ch](sample) = 2;
            gainsPedestal[ch](sample) = dynamicPedestal ? 2 : -1;
        } else if (gainId==1) {
            pedestal = mean_x12[ch];
            gainratio = 1.;
            gainsNoise[ch](sample) = 0;
            gainsPedestal[ch](sample) = dynamicPedestal ? 0 : -1;
        } else if (gainId==2) {
            pedestal = mean_x6[ch];
            gainratio = gain12Over6[ch];
            gainsNoise[ch](sample)  = 1;
            gainsPedestal[ch](sample) = dynamicPedestal ? 1 : -1;
        }

        // TODO: compile time constant -> branch should be non-divergent
        if (dynamicPedestal)
            amplitude = static_cast<float>(adc) * gainratio;
        else
            amplitude = (static_cast<float>(adc) - pedestal) * gainratio;

        amplitudes[ch][sample] = amplitude;

#ifdef ECAL_RECO_CUDA_DEBUG
        printf("%d %d %d %d %f %f %f\n", tx, ch, sample, adc, amplitude,
            pedestal, gainratio);
        if (adc==0)
            printf("adc is zero\n");
#endif

        // TODO this needs to be fixed
        float max_amplitude = 1;
        float shape_value = 1;

        if (sample == 0) {
            // note, no syncing as the same thread will be accessing here
            bool hasGainSwitch = shr_hasSwitchToGain6[threadIdx.x]
                | shr_hasSwitchToGain1[threadIdx.x]
                | shr_isSaturated[threadIdx.x+3];
            // TODO: divergent branch (groups of 10 will enter)
            if (hasGainSwitch && gainSwitchUseMaxSample) {
                // thread for sample=0 will access the right guys
                energies[ch] = max_amplitude / shape_value;
            }
        }
    }
}

///
/// assume kernel launch configuration is 
/// ([MAXSAMPLES, MAXSAMPLES], nchannels)
///
__global__
void kernel_prep_2d(EcalPulseCovariance const* pulse_cov_in,
                    FullSampleMatrix* pulse_cov_out,
                    SampleGainVector const* gainNoise,
                    SampleMatrixD const* noisecorrs,
                    float const* rms_x12,
                    float const* rms_x6,
                    float const* rms_x1,
                    float const* gain12Over6,
                    float const* gain6Over1,
                    SampleMatrix* noisecov,
                    PulseMatrixType* pulse_matrix,
                    FullSampleVector const* pulse_shape,
                    BXVectorType const* bxs,
                    bool const* hasSwitchToGain6,
                    bool const* hasSwitchToGain1,
                    bool const* isSaturated) {
    int ch = blockIdx.x;
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;
    constexpr float addPedestalUncertainty = 0.f;
    constexpr bool dynamicPedestal = false;
    constexpr bool simplifiedNoiseModelForGainSwitch = true;
    constexpr int template_samples = EcalPulseShape::TEMPLATESAMPLES;

    for (int iy=ty, ix=tx; ix<=template_samples && iy<=template_samples; 
        ix+=nsamples, iy+=nsamples)
        pulse_cov_out[ch](iy+7, ix+7) = pulse_cov_in[ch].covval[iy][ix];
    
    bool tmp0 = hasSwitchToGain6[ch];
    bool tmp1 = hasSwitchToGain1[ch];
    bool tmp2 = isSaturated[ch];
    bool hasGainSwitch = tmp0 | tmp1 | tmp2;
    // non-divergent branch for all threads per block
    if (hasGainSwitch) {
        // TODO: did not include simplified noise model
        float noise_value = 0;

        // non-divergent branch - all threads per block
        // TODO: all of these constants indicate that 
        // that these parts could be splitted into completely different 
        // kernels and run one of them only depending on the config
        if (simplifiedNoiseModelForGainSwitch) {
            int isample_max = 5; // according to cpu defs
            int gainidx = gainNoise[ch][isample_max];

            // non-divergent branches
            if (gainidx==0)
                noise_value = rms_x12[ch]*rms_x12[ch]*noisecorrs[0](ty, tx);
            if (gainidx==1) 
                noise_value = gain12Over6[ch]*gain12Over6[ch] * rms_x6[ch]*rms_x6[ch]
                    *noisecorrs[1](ty, tx);
            if (gainidx==2)
                noise_value = gain6Over1[ch]*gain6Over1[ch] * rms_x1[ch]*rms_x1[ch]
                    *noisecorrs[2](ty, tx);
            if (!dynamicPedestal && addPedestalUncertainty>0.f)
                noise_value += addPedestalUncertainty*addPedestalUncertainty;
        } else {
            // 
            int gainidx=0;
            char mask = gainidx;
            int pedestal = gainNoise[ch][ty] == mask ? 1 : 0;
            noise_value += /* gainratio is 1*/ rms_x12[ch]*rms_x12[ch]
                *pedestal*noisecorrs[0](ty, tx);
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noise_value += /* gainratio is 1 */
                    addPedestalUncertainty*addPedestalUncertainty*pedestal;
            }

            //
            gainidx=1;
            mask = gainidx;
            pedestal = gainNoise[ch][ty] == mask ? 1 : 0;
            noise_value += gain12Over6[ch]*gain12Over6[ch]
                *rms_x6[ch]*rms_x6[ch]*pedestal*noisecorrs[1](ty, tx);
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noise_value += gain12Over6[ch]*gain12Over6[ch]
                    *addPedestalUncertainty*addPedestalUncertainty
                    *pedestal;
            }
            
            //
            gainidx=2;
            mask = gainidx;
            pedestal = gainNoise[ch][ty] == mask ? 1 : 0;
            float tmp = gain6Over1[ch] * gain12Over6[ch];
            noise_value += tmp*tmp * rms_x1[ch]*rms_x1[ch]
                *pedestal*noisecorrs[2](ty, tx);
            // non-divergent branch
            if (!dynamicPedestal && addPedestalUncertainty>0.f) {
                noise_value += tmp*tmp * addPedestalUncertainty*addPedestalUncertainty
                    * pedestal;
            }
        }

        noisecov[ch](ty, tx) = noise_value;
    } else {
        auto rms = rms_x12[ch];
        float noise_value = rms*rms * noisecorrs[0](ty, tx);
        if (!dynamicPedestal && addPedestalUncertainty>0.f)
            noise_value += addPedestalUncertainty*addPedestalUncertainty;
        noisecov[ch](ty, tx) = noise_value;
    }

    // pulse matrix
    int bx = (*bxs)(tx);
    int offset = 7 - 3 - bx;
    float value = pulse_shape[ch](offset + ty);
    pulse_matrix[ch](ty, tx) = value;
}

__device__
bool update_covariance(SampleMatrix const& noisecov,
                       FullSampleMatrix const& full_pulse_cov,
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
        int offset = 7 - 3 - bx;

        float value = amplitudes.coeff(ipulse);
        float value_sq = value*value;

        unsigned int nsample_pulse = nsamples - first_sample_t;
        inverse_cov.block(first_sample_t, first_sample_t, 
                          nsample_pulse, nsample_pulse)
            += value_sq * full_pulse_cov.block(first_sample_t + offset,
                                               first_sample_t + offset,
                                               nsample_pulse,
                                               nsample_pulse);
    }

    covariance_decomposition.compute(inverse_cov);
    return true;
}

__device__
float compute_chi2(SampleDecompLLT& covariance_decomposition,
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
                     FullSampleMatrix const* full_pulse_cov,
                     BXVectorType const* bxs,
                     SampleVector const* samples,
                     SampleVector* amplitudes,
                     float* energies,
                     PulseMatrixType* pulse_matrix, 
                     bool* statuses,
                     float* chi2s,
                     bool* isSaturated,
                     bool* hasSwitchToGain6,
                     bool* hasSwitchToGain1,
                     int nchannels,
                     int max_iterations, 
                     bool gainSwitchUseMaxSample) {
    int idx = threadIdx.x + blockDim.x*blockIdx.x;
    if (idx < nchannels) {
        bool hasGainSwitch = isSaturated[idx] 
            | hasSwitchToGain6[idx]
            | hasSwitchToGain1[idx];
        // TODO: gainSwitchUseMaxSimple depends on eb/ee
        // in principle can be splitted/removed into different kernels
        // for ee non-divergent branch
        if (hasGainSwitch && gainSwitchUseMaxSample)
            return;
        bool status = false;
        int iter = 0;
        SampleDecompLLT covariance_decomposition;
        SampleMatrix inverse_cov;
        int npassive = 0;
        amplitudes[idx] = SampleVector::Zero();
        float chi2 = 0;
        float chi2_now = 0;
        while (true) {
            if (iter >= max_iterations)
                break;

            // TODO
            status = update_covariance(
                noisecov[idx], 
                full_pulse_cov[idx],
                inverse_cov,
                *bxs,
                covariance_decomposition,
                amplitudes[idx]);
            if (!status) 
                break;

            // TODO
            SampleMatrix A = covariance_decomposition.matrixL()
                .solve(pulse_matrix[idx]);
            SampleVector b = covariance_decomposition.matrixL()
                .solve(samples[idx]);
            
            status = inplace_fnnls(
                A, b, amplitudes[idx],
                npassive);
                
            if (!status)
                break;

            // TODO
            chi2_now = compute_chi2(
                covariance_decomposition,
                pulse_matrix[idx],
                amplitudes[idx],
                samples[idx]);
            float deltachi2 = chi2_now - chi2;
            chi2 = chi2_now;
            ++iter;

            if (ecal::abs(deltachi2) < 1e-3)
                break;
        }

        float energy = amplitudes[idx](5);
#ifdef ECAL_RECO_CUDA_DEBUG
        printf("%d %d %f %f %f\n", idx, iter, energy, chi2, chi2_now);
        if (iter==1 && energy==0) {
            printf("%d sol amplitudes %f %f %f %f %f %f %f %f %f %f\n", 
                idx,
                amplitudes[idx](0),
                amplitudes[idx](1),
                amplitudes[idx](2),
                amplitudes[idx](3),
                amplitudes[idx](4),
                amplitudes[idx](5),
                amplitudes[idx](6),
                amplitudes[idx](7),
                amplitudes[idx](8),
                amplitudes[idx](9)
            );
            printf("%d samples %f %f %f %f %f %f %f %f %f %f\n", 
                idx,
                samples[idx](0),
                samples[idx](1),
                samples[idx](2),
                samples[idx](3),
                samples[idx](4),
                samples[idx](5),
                samples[idx](6),
                samples[idx](7),
                samples[idx](8),
                samples[idx](9)
            );
        }
#endif
        energies[idx] = energy; // according to bxs vector bxs[5] = 0
        statuses[idx] = status;
        chi2s[idx] = chi2;
    }
}

///
/// Build an Ecal RecHit.
/// TODO: Use SoA data structures on the host directly
/// the reason for removing this from minimize kernel is to isolate the minimize + 
/// again, building an aos rec hit involves strides... -> bad memory access pattern
///
__global__
void kernel_build_rechit(
    float const* energies,
    float const* chi2s,
    uint32_t* dids,
    EcalUncalibratedRecHit* rechits,
    int nchannels) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx < nchannels) {
        rechits[idx] = EcalUncalibratedRecHit{dids[idx], energies[idx],
            0, 0, chi2s[idx], 0};
    }
}

void scatter(host_data& h_data, device_data& d_data, conf_data const& conf) {

/*
void scatter(EcalDigiCollection const& digis,
             EcalUncalibratedRecHitCollection& rechits,
             std::vector<EcalPedestal> const& vpedestals,
             std::vector<EcalMGPAGainRatio> const& vgains,
             std::vector<EcalXtalGroupId> const& vxtals,
             std::vector<EcalPulseShape> const& vpulses,
             std::vector<EcalPulseCovariance> const& vcovariances,
             SampleMatrixGainArray const& noisecors,
             device_data &d_data) {
*/
    auto const& ids = h_data.digis->ids();
    auto const& digis_data = h_data.digis->data();
    using digis_type = std::vector<uint16_t>;
    using dids_type = std::vector<uint32_t>;
    bool barrel = 
        DetId{h_data.digis->begin()->id()}
            .subdetId() == EcalBarrel;
    bool gainSwitchUseMaxSample = barrel; // accodring to the cpu setup
    
    //
    // TODO: remove per event alloc/dealloc -> do once at the start
    //
    /*
    cudaMalloc((void**)&d_digis_data,
        digis_data.size() * sizeof(digis_type::value_type));
    cudaMalloc((void**)&d_ids,
        ids.size() * sizeof(dids_type::value_type));
    cudaMalloc((void**)&d_pedestals,
        vpedestals.size() * sizeof(EcalPedestal));
    cudaMalloc((void**)&d_gains, 
        vgains.size() * sizeof(EcalMGPAGainRatio));
    cudaMalloc((void**)&d_xtals,
        vxtals.size() * sizeof(EcalXtalGroupId));
    cudaMalloc((void**)&d_shapes,
        vpulses.size() * sizeof(EcalPulseShape));
    cudaMalloc((void**)&d_covariances,
        vcovariances.size() * sizeof(EcalPulseCovariance));
    cudaMalloc((void**)&d_rechits,
        rechits.size() * sizeof(EcalUncalibratedRecHit));
    cudaMalloc((void**)&d_noisecors,
        noisecors.size() * sizeof(SampleMatrix));
    ecal::cuda::assert_if_error();
    */

    // 
    // copy to the device
    // TODO: can conditions be copied only once when updated?
    //
    cudaMemcpy(d_data.digis_data, digis_data.data(),
        digis_data.size() * sizeof(digis_type::value_type),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.ids, ids.data(),
        ids.size() * sizeof(dids_type::value_type),
        cudaMemcpyHostToDevice);
 
//    cudaMemcpy(d_data.pedestals, h_data.pedestals->data(),
//        h_data.pedestals->size() * sizeof(EcalPedestal),
//        cudaMemcpyHostToDevice);

    cudaMemcpy(d_data.mean_x12, h_data.ped_data.mean_x12.data(),
        h_data.ped_data.mean_x12.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.rms_x12, h_data.ped_data.rms_x12.data(),
        h_data.ped_data.rms_x12.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.mean_x6, h_data.ped_data.mean_x6.data(),
        h_data.ped_data.mean_x6.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.rms_x6, h_data.ped_data.rms_x6.data(),
        h_data.ped_data.rms_x6.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.mean_x1, h_data.ped_data.mean_x1.data(),
        h_data.ped_data.mean_x1.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.rms_x1, h_data.ped_data.rms_x1.data(),
        h_data.ped_data.rms_x1.size() * sizeof(float),
        cudaMemcpyHostToDevice);

//    cudaMemcpy(d_data.gains, h_data.gains->data(),
//        h_data.gains->size() * sizeof(EcalMGPAGainRatio),
//        cudaMemcpyHostToDevice);

    cudaMemcpy(d_data.gain12Over6, h_data.gainratio_data.gain12Over6.data(),
        h_data.gainratio_data.gain12Over6.size() * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.gain6Over1, h_data.gainratio_data.gain6Over1.data(),
        h_data.gainratio_data.gain6Over1.size() * sizeof(float),
        cudaMemcpyHostToDevice);

    cudaMemcpy(d_data.xtals, h_data.xtals->data(),
        h_data.xtals->size() * sizeof(EcalXtalGroupId),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.pulses, h_data.pulse_shapes->data(),
        h_data.pulse_shapes->size() * sizeof(EcalPulseShape),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.covariances, h_data.pulse_covariances->data(),
        h_data.pulse_covariances->size() * sizeof(EcalPulseCovariance),
        cudaMemcpyHostToDevice);

    cudaMemcpy(d_data.noisecorrs, h_data.noisecorrs->data(),
        h_data.noisecorrs->size() * sizeof(SampleMatrixD),
        cudaMemcpyHostToDevice);

    cudaMemcpy(d_data.sample_mask, h_data.sample_mask,
        sizeof(EcalSampleMask),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.EBTimeCorrAmplitudeBins, 
        h_data.time_bias_corrections->EBTimeCorrAmplitudeBins.data(),
        sizeof(float) * h_data.time_bias_corrections->EBTimeCorrAmplitudeBins.size(),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.EBTimeCorrShiftBins, 
        h_data.time_bias_corrections->EBTimeCorrShiftBins.data(),
        sizeof(float) * h_data.time_bias_corrections->EBTimeCorrShiftBins.size(),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.EETimeCorrAmplitudeBins, 
        h_data.time_bias_corrections->EETimeCorrAmplitudeBins.data(),
        sizeof(float) * h_data.time_bias_corrections->EETimeCorrAmplitudeBins.size(),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.EETimeCorrShiftBins, 
        h_data.time_bias_corrections->EETimeCorrShiftBins.data(),
        sizeof(float) * h_data.time_bias_corrections->EETimeCorrShiftBins.size(),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.weights,
        h_data.weights->data(),
        sizeof(EMatrix) * h_data.weights->size(),
        cudaMemcpyHostToDevice);
    cudaMemcpy(d_data.bxs, h_data.bxs,
        sizeof(BXVectorType),
        cudaMemcpyHostToDevice);
    ecal::cuda::assert_if_error();

    int nthreads_per_block = conf.threads.x;
    int nblocks = (h_data.digis->size() + nthreads_per_block - 1) / nthreads_per_block;

    // 
    // 1d preparation kernel
    //
    unsigned int nchannels_per_block = 32;
    unsigned int threads_1d = 10 * nchannels_per_block;
    unsigned int blocks_1d = threads_1d > 10*h_data.digis->size() 
        ? 1 : (h_data.digis->size()*10 + threads_1d - 1) / threads_1d;
    int shared_bytes = nchannels_per_block * EcalDataFrame::MAXSAMPLES * (
        sizeof(bool) + sizeof(bool) + sizeof(bool) + sizeof(bool) + sizeof(char)
    );
    std::cout << "nchannels = " << h_data.digis->size() << std::endl;
    std::cout << "shared memory per block = " << shared_bytes << "B\n";
    kernel_prep_1d<<<blocks_1d, threads_1d, shared_bytes>>>(
        d_data.pulses, d_data.epulses,
        d_data.digis_data, d_data.samples,
        d_data.gainsNoise,
        d_data.gainsPedestal,
        d_data.mean_x1,
        d_data.mean_x12,
        d_data.mean_x6,
        d_data.gain6Over1,
        d_data.gain12Over6,
        d_data.hasSwitchToGain6,
        d_data.hasSwitchToGain1,
        d_data.isSaturated,
        d_data.energies,
        gainSwitchUseMaxSample,
        h_data.digis->size());
    cudaDeviceSynchronize();
    ecal::cuda::assert_if_error();

    std::cout << " after kernel prep 1d\n";

    //
    // 2d preparation kernel
    //
    int blocks_2d = h_data.digis->size();
    dim3 threads_2d{10, 10};
    kernel_prep_2d<<<blocks_2d, threads_2d>>>(
        d_data.covariances, d_data.pulse_covariances,
        d_data.gainsNoise,
        d_data.noisecorrs,
        d_data.rms_x12,
        d_data.rms_x6,
        d_data.rms_x1,
        d_data.gain12Over6,
        d_data.gain6Over1,
        d_data.noisecov,
        d_data.pulse_matrix,
        d_data.epulses,
        d_data.bxs,
        d_data.hasSwitchToGain6,
        d_data.hasSwitchToGain1,
        d_data.isSaturated);
    cudaDeviceSynchronize();
    ecal::cuda::assert_if_error();

    std::cout << "after kernel prep 2d\n";
//    kernel_minimize<<<>>>();

    unsigned int threads_min = conf.threads.x;
    unsigned int blocks_min = threads_min > h_data.digis->size()
        ? 1 : (h_data.digis->size() + threads_min - 1) / threads_min;
    kernel_minimize<<<blocks_min, threads_min>>>(
        d_data.noisecov,
        d_data.pulse_covariances,
        d_data.bxs,
        d_data.samples,
        d_data.amplitudes,
        d_data.energies,
        d_data.pulse_matrix,
        d_data.statuses,
        d_data.chi2,
        d_data.isSaturated,
        d_data.hasSwitchToGain6,
        d_data.hasSwitchToGain1,
        h_data.digis->size(),
        50,
        gainSwitchUseMaxSample);
    cudaDeviceSynchronize();
    ecal::cuda::assert_if_error();
/*
    kernel_build_rechit<<<blocks_min, threads_min>>>(
        d_data.energies,
        d_data.chi2,
        d_data.ids,
        d_data.rechits,
        h_data.digis->size());
    cudaDeviceSynchronize();
    ecal::cuda::assert_if_error();
    */

  /*  kernel_reconstruct<<<nblocks, nthreads_per_block>>>(
        d_data.digis_data,
        d_data.ids,*/
        /* d_rechits, */
/*        d_data.pedestals,
        d_data.gains,
        d_data.xtals,
        d_data.pulses,
        d_data.covariances,
        d_data.rechits,
        d_data.noisecors,
        d_data.sample_mask,
        d_data.EBTimeCorrAmplitudeBins, 
        h_data.time_bias_corrections->EBTimeCorrAmplitudeBins.size(),
        d_data.EBTimeCorrShiftBins, 
        h_data.time_bias_corrections->EBTimeCorrShiftBins.size(),
        d_data.EETimeCorrAmplitudeBins, 
        h_data.time_bias_corrections->EETimeCorrAmplitudeBins.size(),
        d_data.EETimeCorrShiftBins, 
        h_data.time_bias_corrections->EETimeCorrShiftBins.size(),
        d_data.weights,
        h_data.digis->size()
    );*/

    //
    // transfer the results back
    //
    cudaMemcpy(&(*h_data.rechits_soa.amplitude.begin()),
               d_data.energies,
               h_data.rechits_soa.amplitude.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*h_data.rechits_soa.chi2.begin()),
               d_data.chi2,
               h_data.rechits_soa.chi2.size() * sizeof(float),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&(*h_data.rechits_soa.did.begin()),
               d_data.ids,
               h_data.rechits_soa.did.size() * sizeof(uint32_t),
               cudaMemcpyDeviceToHost);

//    cudaMemcpy(&(*h_data.rechits->begin()), d_data.rechits,
//        h_data.rechits->size() * sizeof(EcalUncalibratedRecHit),
//        cudaMemcpyDeviceToHost);

    // 
    // free all the device ptrs
    // TODO: remove per event dealloc
    //
    /*
    cudaFree(d_digis_data);
    cudaFree(d_ids);
    cudaFree(d_pedestals);
    cudaFree(d_gains);
    cudaFree(d_xtals);
    cudaFree(d_shapes);
    cudaFree(d_covariances);
    cudaFree(d_rechits);
    cudaFree(d_noisecors);
    ecal::cuda::assert_if_error();
    */
}

}}}

/*
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"

EcalUncalibRecHitMultiFitAlgo::EcalUncalibRecHitMultiFitAlgo() : 
  _computeErrors(true),
  _doPrefit(false),
  _prefitMaxChiSq(1.0),
  _dynamicPedestals(false),
  _mitigateBadSamples(false),
  _selectiveBadSampleCriteria(false),
  _addPedestalUncertainty(0.),
  _simplifiedNoiseModelForGainSwitch(true),
  _gainSwitchUseMaxSample(false){
    
  _singlebx.resize(1);
  _singlebx << 0;
  
  _pulsefuncSingle.disableErrorCalculation();
  _pulsefuncSingle.setMaxIters(1);
  _pulsefuncSingle.setMaxIterWarnings(false);
    
}

/// compute rechits
EcalUncalibratedRecHit EcalUncalibRecHitMultiFitAlgo::makeRecHit(const EcalDataFrame& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const SampleMatrixGainArray &noisecors, const FullSampleVector &fullpulse, const FullSampleMatrix &fullpulsecov, const BXVector &activeBX) {

  uint32_t flags = 0;
  
  const unsigned int nsample = EcalDataFrame::MAXSAMPLES;
  
  double maxamplitude = -std::numeric_limits<double>::max();
  const unsigned int iSampleMax = 5;
  const unsigned int iFullPulseMax = 9;
  
  double pedval = 0.;
    
  SampleVector amplitudes;
  SampleGainVector gainsNoise;
  SampleGainVector gainsPedestal;
  SampleGainVector badSamples = SampleGainVector::Zero();
  bool hasSaturation = dataFrame.isSaturated();
  bool hasGainSwitch = hasSaturation || dataFrame.hasSwitchToGain6() || dataFrame.hasSwitchToGain1();
  
  //no dynamic pedestal in case of gain switch, since then the fit becomes too underconstrained
  bool dynamicPedestal = _dynamicPedestals && !hasGainSwitch;
  
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
        
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    
    double amplitude = 0.;
    int gainId = sample.gainId();
    
    double pedestal = 0.;
    double gainratio = 1.;
    
    if (gainId==0 || gainId==3) {
      pedestal = aped->mean_x1;
      gainratio = aGain->gain6Over1()*aGain->gain12Over6();
      gainsNoise[iSample] = 2;
      gainsPedestal[iSample] = dynamicPedestal ? 2 : -1;  //-1 for static pedestal
    }
    else if (gainId==1) {
      pedestal = aped->mean_x12;
      gainratio = 1.;
      gainsNoise[iSample] = 0;
      gainsPedestal[iSample] = dynamicPedestal ? 0 : -1; //-1 for static pedestal
    }
    else if (gainId==2) {
      pedestal = aped->mean_x6;
      gainratio = aGain->gain12Over6();
      gainsNoise[iSample] = 1;
      gainsPedestal[iSample] = dynamicPedestal ? 1 : -1; //-1 for static pedestals
    }

    if (dynamicPedestal) {
      amplitude = (double)(sample.adc())*gainratio;
    }
    else {
      amplitude = ((double)(sample.adc()) - pedestal) * gainratio;
    }
    
    if (gainId == 0) {
       edm::LogError("EcalUncalibRecHitMultiFitAlgo")<< "Saturation encountered.  Multifit is not intended to be used for saturated channels.";
      //saturation
      if (dynamicPedestal) {
        amplitude = 4095.*gainratio;
      }
      else {
        amplitude = (4095. - pedestal) * gainratio;
      }
    }
        
    amplitudes[iSample] = amplitude;
    
    if (iSample==iSampleMax) {
      maxamplitude = amplitude;
      pedval = pedestal;
    }
        
  }

  double amplitude, amperr, chisq;
  bool status = false;
    
  //special handling for gain switch, where sample before maximum is potentially affected by slew rate limitation
  //optionally apply a stricter criteria, assuming slew rate limit is only reached in case where maximum sample has gain switched but previous sample has not
  //option 1: use simple max-sample algorithm
  if (hasGainSwitch && _gainSwitchUseMaxSample) {
    double maxpulseamplitude = maxamplitude / fullpulse[iFullPulseMax];
    EcalUncalibratedRecHit rh( dataFrame.id(), maxpulseamplitude, pedval, 0., 0., flags );
    rh.setAmplitudeError(0.);
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      int bx = _pulsefunc.BXs().coeff(ipulse);
      if (bx!=0) {
        rh.setOutOfTimeAmplitude(bx+5, 0.0);
      }
    }
    return rh;
  }

  //option2: A floating negative single-sample offset is added to the fit
  //such that the affected sample is treated only as a lower limit for the true amplitude
  bool mitigateBadSample = _mitigateBadSamples && hasGainSwitch && iSampleMax>0;
  mitigateBadSample &= (!_selectiveBadSampleCriteria || (gainsNoise.coeff(iSampleMax-1)!=gainsNoise.coeff(iSampleMax)) );
  if (mitigateBadSample) {
    badSamples[iSampleMax-1] = 1;
  }
  
  //compute noise covariance matrix, which depends on the sample gains
  SampleMatrix noisecov;
  if (hasGainSwitch) {
    std::array<double,3> pedrmss = {{aped->rms_x12, aped->rms_x6, aped->rms_x1}};
    std::array<double,3> gainratios = {{ 1., aGain->gain12Over6(), aGain->gain6Over1()*aGain->gain12Over6()}};
    if (_simplifiedNoiseModelForGainSwitch) {
      int gainidxmax = gainsNoise[iSampleMax];
      noisecov = gainratios[gainidxmax]*gainratios[gainidxmax]*pedrmss[gainidxmax]*pedrmss[gainidxmax]*noisecors[gainidxmax];
      if (!dynamicPedestal && _addPedestalUncertainty>0.) {
        //add fully correlated component to noise covariance to inflate pedestal uncertainty
        noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
      }
    }
    else {
      noisecov = SampleMatrix::Zero();
      for (unsigned int gainidx=0; gainidx<noisecors.size(); ++gainidx) {
        SampleGainVector mask = gainidx*SampleGainVector::Ones();
        SampleVector pedestal = (gainsNoise.array()==mask.array()).cast<SampleVector::value_type>();
        if (pedestal.maxCoeff()>0.) {
          //select out relevant components of each correlation matrix, and assume no correlation between samples with
          //different gain
          noisecov += gainratios[gainidx]*gainratios[gainidx]*pedrmss[gainidx]*pedrmss[gainidx]*pedestal.asDiagonal()*noisecors[gainidx]*pedestal.asDiagonal();
          if (!dynamicPedestal && _addPedestalUncertainty>0.) {
            //add fully correlated component to noise covariance to inflate pedestal uncertainty
            noisecov += gainratios[gainidx]*gainratios[gainidx]*_addPedestalUncertainty*_addPedestalUncertainty*pedestal.asDiagonal()*SampleMatrix::Ones()*pedestal.asDiagonal();
          }
        }
      }
    }
  }
  else {
    noisecov = aped->rms_x12*aped->rms_x12*noisecors[0];
    if (!dynamicPedestal && _addPedestalUncertainty>0.) {
      //add fully correlated component to noise covariance to inflate pedestal uncertainty
      noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
    }
  }
  
  //optimized one-pulse fit for hlt
  bool usePrefit = false;
  if (_doPrefit) {
    status = _pulsefuncSingle.DoFit(amplitudes,noisecov,_singlebx,fullpulse,fullpulsecov,gainsPedestal,badSamples);
    amplitude = status ? _pulsefuncSingle.X()[0] : 0.;
    amperr = status ? _pulsefuncSingle.Errors()[0] : 0.;
    chisq = _pulsefuncSingle.ChiSq();
    
    if (chisq < _prefitMaxChiSq) {
      usePrefit = true;
    }
  }
  
  if (!usePrefit) {
  
    if(!_computeErrors) _pulsefunc.disableErrorCalculation();
    status = _pulsefunc.DoFit(amplitudes,noisecov,activeBX,fullpulse,fullpulsecov,gainsPedestal,badSamples);
    chisq = _pulsefunc.ChiSq();
    
    if (!status) {
      edm::LogWarning("EcalUncalibRecHitMultiFitAlgo::makeRecHit") << "Failed Fit" << std::endl;
    }

    unsigned int ipulseintime = 0;
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      if (_pulsefunc.BXs().coeff(ipulse)==0) {
        ipulseintime = ipulse;
        break;
      }
    }
    
    amplitude = status ? _pulsefunc.X()[ipulseintime] : 0.;
    amperr = status ? _pulsefunc.Errors()[ipulseintime] : 0.;
  
  }
  
  double jitter = 0.;
  
  EcalUncalibratedRecHit rh( dataFrame.id(), amplitude , pedval, jitter, chisq, flags );
  rh.setAmplitudeError(amperr);
  
  if (!usePrefit) {
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      int bx = _pulsefunc.BXs().coeff(ipulse);
      if (bx!=0 && std::abs(bx)<100) {
        rh.setOutOfTimeAmplitude(bx+5, status ? _pulsefunc.X().coeff(ipulse) : 0.);
      }
      else if (bx==(100+gainsPedestal[iSampleMax])) {
        rh.setPedestal(status ? _pulsefunc.X().coeff(ipulse) : 0.);
      }
    }
  }
  
  return rh;
}
*/
