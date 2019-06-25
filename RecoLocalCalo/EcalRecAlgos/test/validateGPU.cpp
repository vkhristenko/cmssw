#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "CUDADataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"
    
constexpr float eps_diff = 1e-3;

struct Histos {
    TH1D *hSOIAmplitudesEBGPU, *hSOIAmplitudesEEGPU;
    TH1D *hSOIAmplitudesEBCPU, *hSOIAmplitudesEECPU;
    TH1D *hChi2EBGPU, *hChi2EEGPU;
    TH1D *hChi2EBCPU, *hChi2EECPU;
    TH2D *hSOIAmplitudesEBGPUvsCPU, *hSOIAmplitudesEEGPUvsCPU;
    TH2D *hChi2EBGPUvsCPU, *hChi2EEGPUvsCPU;
    TH1D *hAmplitudesEB, *hAmplitudesEE;
    TH1D *hAmplitudesAbsDiffEB, *hAmplitudesAbsDiffEE;
    TH1D *hChi2AbsDiffEB, *hChi2AbsDiffEE;
    TH1I *hAmplBitsEB, *hAmplBitsEE;
    TH1I *hChi2BitsEB, *hChi2BitsEE;
};

struct Stats {
    int ndiscrs{0}, n1pdiscrs{0};
};

template<typename TCPU, typename TGPU>
void accumulate(
        TCPU const& cpu, TGPU const& gpu, 
        Histos& histos, Stats& stats, int ie, int det, int offset) {
    // make sure the size of collections is the same
    assert(gpu.amplitude.size() == cpu.size());
    assert(gpu.amplitude.size() == cpu.size());

    auto const nchannels = cpu.size();

    // iterate and accumulate
    for (uint32_t i=0; i<nchannels; ++i) {
        auto const soi_amp_gpu = gpu.amplitude[i];
        auto const soi_amp_cpu = cpu[i].amplitude();
        auto const chi2_gpu = gpu.chi2[i];
        auto const chi2_cpu = cpu[i].chi2();
        auto const abs_diff_amp = soi_amp_cpu == 0 
            ? 0
            : (soi_amp_gpu - soi_amp_cpu) / soi_amp_cpu;
        auto const abs_diff_chi2 = chi2_cpu == 0
            ? 0
            : (chi2_gpu - chi2_cpu) / chi2_cpu;

        // direct reinterpret from float to uint32_t 
        // -> breaks strict aliasing rules
        auto const* amp_bits_gpu = reinterpret_cast<uint32_t const*>(reinterpret_cast<char const*>(&soi_amp_gpu));
        auto const* amp_bits_cpu = reinterpret_cast<uint32_t const*>(reinterpret_cast<char const*>(&soi_amp_cpu));
        auto const* chi2_bits_gpu = reinterpret_cast<uint32_t const*>(reinterpret_cast<char const*>(&chi2_gpu));
        auto const* chi2_bits_cpu = reinterpret_cast<uint32_t const*>(reinterpret_cast<char const*>(&chi2_cpu));
        auto const amp_xor = 
            *amp_bits_gpu ^ 
            *amp_bits_cpu;
        auto const chi2_xor = 
            *chi2_bits_gpu ^ 
            *chi2_bits_cpu;

        // propogate bits diff
        if (det == 0) {
            for (int ibit=0; ibit<32; ibit++) {
                uint32_t mask = 1 << ibit;
                if ((amp_xor & mask) != 0x0)
                    histos.hAmplBitsEB->Fill(ibit);
                if ((chi2_xor & mask) != 0x0)
                    histos.hChi2BitsEB->Fill(ibit);
            }
        } else {
            for (int ibit=0; ibit<32; ibit++) {
                uint32_t mask = 1 << ibit;
                if ((amp_xor & mask) != 0x0)
                    histos.hAmplBitsEE->Fill(ibit);
                if ((chi2_xor & mask) != 0x0)
                    histos.hChi2BitsEE->Fill(ibit);
            }
        }

        if (det == 0) {
            histos.hSOIAmplitudesEBGPU->Fill(soi_amp_gpu);
            histos.hSOIAmplitudesEBCPU->Fill(soi_amp_cpu);
            histos.hSOIAmplitudesEBGPUvsCPU->Fill(soi_amp_cpu, soi_amp_gpu);
            histos.hChi2EBGPU->Fill(chi2_gpu);
            histos.hChi2EBCPU->Fill(chi2_cpu);
            histos.hChi2EBGPUvsCPU->Fill(chi2_cpu, chi2_gpu);
            histos.hAmplitudesAbsDiffEB->Fill(abs_diff_amp);
            histos.hChi2AbsDiffEB->Fill(abs_diff_chi2);
        } else {
            histos.hSOIAmplitudesEEGPU->Fill(soi_amp_gpu);
            histos.hSOIAmplitudesEECPU->Fill(soi_amp_cpu);
            histos.hSOIAmplitudesEEGPUvsCPU->Fill(soi_amp_cpu, soi_amp_gpu);
            histos.hChi2EEGPU->Fill(chi2_gpu);
            histos.hChi2EECPU->Fill(chi2_cpu);
            histos.hChi2EEGPUvsCPU->Fill(chi2_cpu, chi2_gpu);
            histos.hAmplitudesAbsDiffEE->Fill(abs_diff_amp);
            histos.hChi2AbsDiffEE->Fill(abs_diff_chi2);
        }

        if (std::abs(soi_amp_gpu - soi_amp_cpu) >= eps_diff ||
            std::isnan(soi_amp_gpu)) {
            printf("--- amp diff %s eventid = %d chid = %d ---\n", 
                det==0 ? "eb" : "ee", ie, offset+i);
            printf(">>> amp_gpu = %f amp_cpu %f chi2_gpu = %f chi2_cpu = %f\n",
                soi_amp_gpu, soi_amp_cpu, chi2_gpu, chi2_cpu);
            printf(">>> bit-wise diff: amp = 0x%08x chi2 = 0x%08x\n", amp_xor, chi2_xor);
            stats.ndiscrs++;
        }
        
        if (std::abs(chi2_gpu - chi2_cpu) >= eps_diff || std::isnan(chi2_gpu)) {
            printf("--- chi2 diff %s eventid = %d chid = %d ---\n",
                det==0 ? "eb" : "ee", ie, i);
            printf(">>> amp_gpu = %f amp_cpu %f chi2_gpu = %f chi2_cpu = %f\n",
                soi_amp_gpu, soi_amp_cpu, chi2_gpu, chi2_cpu);
            printf(">>> bit-wise diff: amp = 0x%08x chi2 = 0x%08x\n", amp_xor, chi2_xor);
            stats.ndiscrs++;
        }
        if (std::isnan(chi2_gpu) || std::isnan(soi_amp_gpu))
            printf("*** nan ***\n");
        
        if (std::abs(abs_diff_amp)*100 > 1 ||
            std::abs(abs_diff_chi2)*100 > 1)
            stats.n1pdiscrs++; 
    }
}

int main(int argc, char *argv[]) {
    if (argc<3) {
        std::cout << "run with: ./validateGPU <path to input file> <output file>\n";
        exit(0);
    }

    edm::Wrapper<ecal::UncalibratedRecHit<ecal::Tag::soa>> *wgpuEB=nullptr;
    edm::Wrapper<ecal::UncalibratedRecHit<ecal::Tag::soa>> *wgpuEE=nullptr;
    edm::Wrapper<EBUncalibratedRecHitCollection> *wcpuEB = nullptr;
    edm::Wrapper<EEUncalibratedRecHitCollection> *wcpuEE = nullptr;


    std::string fileName = argv[1];
    std::string outFileName = argv[2];

    // output
    TFile rfout{outFileName.c_str(), "recreate"};
    Histos histos;

    int nbins = 100; int last = 1000;
    histos.hSOIAmplitudesEBGPU = new TH1D("hSOIAmplitudesEBGPU", "hSOIAmplitudesEBGPU",
        nbins, 0, last);
    histos.hSOIAmplitudesEEGPU = new TH1D("hSOIAmplitudesEEGPU", "hSOIAmplitudesEEGPU",
        nbins, 0, last);
    histos.hSOIAmplitudesEBCPU = new TH1D("hSOIAmplitudesEBCPU", "hSOIAmplitudesEBCPU",
        nbins, 0, last);
    histos.hSOIAmplitudesEECPU = new TH1D("hSOIAmplitudesEECPU", "hSOIAmplitudesEECPU",
        nbins, 0, last);

    histos.hAmplBitsEB = new TH1I("AmplBitsEB", "AmplBitsEB", 32, 0, 32);
    histos.hAmplBitsEE = new TH1I("AmplBitsEE", "AmplBitsEE", 32, 0, 32);
    histos.hChi2BitsEB = new TH1I("Chi2BitsEB", "Chi2BitsEB", 32, 0, 32);
    histos.hChi2BitsEE = new TH1I("Chi2BitsEE", "Chi2BitsEE", 32, 0, 32);
    
    int nbins_absdiff = 100;
    histos.hAmplitudesAbsDiffEB = new TH1D("hAmplitudesAbsDiffEB", 
        "hAmplitudesAbsDiffEB",
        nbins_absdiff, -1, 1);
    histos.hAmplitudesAbsDiffEE = new TH1D("hAmplitudesAbsDiffEE", 
        "hAmplitudesAbsDiffEE",
        nbins_absdiff, -1, 1);
    histos.hChi2AbsDiffEB = new TH1D("hChi2AbsDiffEB", "hChi2AbsDiffEB",
        nbins_absdiff, -1, 1);
    histos.hChi2AbsDiffEE = new TH1D("hChi2AbsDiffEE", "hChi2AbsDiffEE",
        nbins_absdiff, -1, 1);
    
    int nbins_chi2 = 100; int last_chi2 = 100;
    histos.hChi2EBGPU = new TH1D("hChi2EBGPU", "hChi2EBGPU",
        nbins_chi2, 0, last_chi2);
    histos.hChi2EEGPU = new TH1D("hChi2EEGPU", "hChi2EEGPU",
        nbins_chi2, 0, last_chi2);
    histos.hChi2EBCPU = new TH1D("hChi2EBCPU", "hChi2EBCPU",
        nbins_chi2, 0, last_chi2);
    histos.hChi2EECPU = new TH1D("hChi2EECPU", "hChi2EECPU",
        nbins_chi2, 0, last_chi2);
    
    histos.hSOIAmplitudesEBGPUvsCPU = new TH2D("hSOIAmplitudesEBGPUvsCPU", 
        "hSOIAmplitudesEBGPUvsCPU", nbins, 0, last, nbins, 0, last);
    histos.hSOIAmplitudesEEGPUvsCPU = new TH2D("hSOIAmplitudesEEGPUvsCPU", 
        "hSOIAmplitudesEEGPUvsCPU", nbins, 0, last, nbins, 0, last);
    
    histos.hChi2EBGPUvsCPU = new TH2D("hChi2EBGPUvsCPU", 
        "hChi2EBGPUvsCPU", nbins_chi2, 0, last_chi2, nbins_chi2, 0, last_chi2);
    histos.hChi2EEGPUvsCPU = new TH2D("hChi2EEGPUvsCPU", 
        "hChi2EEGPUvsCPU", nbins_chi2, 0, last_chi2, nbins_chi2, 0, last_chi2);

    // input
    std::cout << "validating file " << fileName << std::endl;
    TFile rf{fileName.c_str()};
    TTree *rt = (TTree*)rf.Get("Events");
    rt->SetBranchAddress("ecalTagsoaecalUncalibratedRecHit_ecalUncalibRecHitProducerGPU_EcalUncalibRecHitsEB_RECO.", &wgpuEB);
    rt->SetBranchAddress("ecalTagsoaecalUncalibratedRecHit_ecalUncalibRecHitProducerGPU_EcalUncalibRecHitsEE_RECO.", &wgpuEE);
    rt->SetBranchAddress("EcalUncalibratedRecHitsSorted_ecalMultiFitUncalibRecHit_EcalUncalibRecHitsEB_RECO.", &wcpuEB);
    rt->SetBranchAddress("EcalUncalibratedRecHitsSorted_ecalMultiFitUncalibRecHit_EcalUncalibRecHitsEE_RECO.", &wcpuEE);

    Stats stats;
    int nchannelsTotal = 0;

    // accumulate
    auto const nentries = rt->GetEntries();
    std::cout << "#events to validate over: " << nentries << std::endl;
    for (int ie=0; ie<nentries; ++ie) {
        rt->GetEntry(ie);

        auto const neb = wcpuEB->bareProduct().size();
        auto const nee = wcpuEE->bareProduct().size();

        accumulate(wcpuEB->bareProduct(), wgpuEB->bareProduct(),
            histos, stats, ie, 0, 0);
        accumulate(wcpuEE->bareProduct(), wgpuEE->bareProduct(),
            histos, stats, ie, 1, neb);

        nchannelsTotal += neb + nee;
    }

    std::cout 
        << "-----------------------\n" 
        << "--- summary\n"
        << "--- events = " << nentries << std::endl
        << "--- nchannelsTotal = " << nchannelsTotal << std::endl
        << "--- num discrs (amp or chi2) = " << stats.ndiscrs << std::endl
        << "--- percentage of discrs (num discrs / nchannels) = " << static_cast<double>(stats.ndiscrs) / static_cast<double>(nchannelsTotal) << std::endl
        << "--- num of discrs with abs diff >= 1% = " << stats.n1pdiscrs << std::endl
        << "--- num of bit diffs for ampl = "
        << histos.hAmplBitsEB->GetEntries() + histos.hAmplBitsEE->GetEntries() << std::endl
        << "--- percentage of bits discrs for ampl = "
        << static_cast<double>(histos.hAmplBitsEB->GetEntries() + histos.hAmplBitsEE->GetEntries()) / static_cast<double>(32 * nchannelsTotal) << std::endl
        << "--- num of bit diffs for chi2 = "
        << histos.hChi2BitsEB->GetEntries() + histos.hChi2BitsEE->GetEntries() << std::endl
        << "--- percentage of bits discrs for chi2 = "
        << static_cast<double>(histos.hChi2BitsEB->GetEntries() + histos.hChi2BitsEE->GetEntries()) / static_cast<double>(32 * nchannelsTotal) << std::endl;

    rf.Close();
    rfout.Write();
    rfout.Close();
    return 0;
}
