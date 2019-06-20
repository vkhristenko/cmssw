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
    TH1D *hSOIAmplitudesEBGPU, *hSOIAmplitudesEEGPU;
    TH1D *hSOIAmplitudesEBCPU, *hSOIAmplitudesEECPU;
    TH2D *hSOIAmplitudesEBGPUvsCPU, *hSOIAmplitudesEEGPUvsCPU;
    TH1D *hAmplitudesEB, *hAmplitudesEE;

    int nbins = 100; int last = 1000;
    hSOIAmplitudesEBGPU = new TH1D("hSOIAmplitudesEBGPU", "hSOIAmplitudesEBGPU",
        nbins, 0, last);
    hSOIAmplitudesEEGPU = new TH1D("hSOIAmplitudesEEGPU", "hSOIAmplitudesEEGPU",
        nbins, 0, last);
    hSOIAmplitudesEBCPU = new TH1D("hSOIAmplitudesEBCPU", "hSOIAmplitudesEBCPU",
        nbins, 0, last);
    hSOIAmplitudesEECPU = new TH1D("hSOIAmplitudesEECPU", "hSOIAmplitudesEECPU",
        nbins, 0, last);
    
    hSOIAmplitudesEBGPUvsCPU = new TH2D("hSOIAmplitudesEBGPUvsCPU", 
        "hSOIAmplitudesEBGPUvsCPU", nbins, 0, last, nbins, 0, last);
    hSOIAmplitudesEEGPUvsCPU = new TH2D("hSOIAmplitudesEEGPUvsCPU", 
        "hSOIAmplitudesEEGPUvsCPU", nbins, 0, last, nbins, 0, last);

    // input
    std::cout << "validating file " << fileName << std::endl;
    TFile rf{fileName.c_str()};
    TTree *rt = (TTree*)rf.Get("Events");
    rt->SetBranchAddress("ecalTagsoaecalUncalibratedRecHit_ecalUncalibRecHitProducerGPU_EcalUncalibRecHitsEB_RECO.", &wgpuEB);
    rt->SetBranchAddress("ecalTagsoaecalUncalibratedRecHit_ecalUncalibRecHitProducerGPU_EcalUncalibRecHitsEE_RECO.", &wgpuEE);
    rt->SetBranchAddress("EcalUncalibratedRecHitsSorted_ecalMultiFitUncalibRecHit_EcalUncalibRecHitsEB_RECO.", &wcpuEB);
    rt->SetBranchAddress("EcalUncalibratedRecHitsSorted_ecalMultiFitUncalibRecHit_EcalUncalibRecHitsEE_RECO.", &wcpuEE);

    // accumulate
    auto const nentries = rt->GetEntries();
    std::cout << "#events to validate over: " << nentries << std::endl;
    for (int i=0; i<nentries; ++i) {
        rt->GetEntry(i);

        assert(wgpuEB->bareProduct().amplitude.size() == wcpuEB->bareProduct().size());
        assert(wgpuEE->bareProduct().amplitude.size() == wcpuEE->bareProduct().size());
        auto const neb = wcpuEB->bareProduct().size();
        auto const nee = wcpuEE->bareProduct().size();

        for (uint32_t i=0; i<neb; ++i) {
            auto const soi_amp_gpu = wgpuEB->bareProduct().amplitude[i];
            auto const soi_amp_cpu = wcpuEB->bareProduct()[i].amplitude();

            hSOIAmplitudesEBGPU->Fill(soi_amp_gpu);
            hSOIAmplitudesEBCPU->Fill(soi_amp_cpu);
            hSOIAmplitudesEBGPUvsCPU->Fill(soi_amp_cpu, soi_amp_gpu);
        }

        for (uint32_t i=0; i<nee; ++i) {
            auto const soi_amp_gpu = wgpuEE->bareProduct().amplitude[i];
            auto const soi_amp_cpu = wcpuEE->bareProduct()[i].amplitude();

            hSOIAmplitudesEEGPU->Fill(soi_amp_gpu);
            hSOIAmplitudesEECPU->Fill(soi_amp_cpu);

            hSOIAmplitudesEEGPUvsCPU->Fill(soi_amp_cpu, soi_amp_gpu);
        }
    }

    rf.Close();
    rfout.Write();
    rfout.Close();
    return 0;
}
