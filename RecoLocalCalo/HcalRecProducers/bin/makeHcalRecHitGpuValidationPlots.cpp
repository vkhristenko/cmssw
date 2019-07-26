#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TPaveStats.h>

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "CUDADataFormats/HcalRecHitSoA/interface/RecHitCollection.h"

#define CREATE_HIST_1D(varname, nbins, first, last) \
    auto varname = new TH1D(#varname, #varname, nbins, first, last)

#define CREATE_HIST_2D(varname, nbins, first, last) \
    auto varname = new TH2D(#varname, #varname, nbins, first, last, nbins, first, last)

int main(int argc, char *argv[]) {
    if (argc<3) {
        std::cout << "run with: ./<exe> <path to input file> <path to output file>\n";
        exit(0);
    }

    std::string inFileName{argv[1]};
    std::string outFileName{argv[2]};
    
    // branches to use
    edm::Wrapper<HBHERecHitCollection> *wcpu=nullptr;
    edm::Wrapper<hcal::RecHitCollection<hcal::Tag::soa>> *wgpu=nullptr;

    // prep output 
    TFile rfout{outFileName.c_str(), "recreate"};

    CREATE_HIST_1D(hEnergyM0GPU, 1000, 0, 100);
    CREATE_HIST_1D(hEnergyM0CPU, 1000, 0, 100);
    CREATE_HIST_2D(hEnergyM0GPUvsCPU, 1000, 0, 100);

    // prep input
    TFile rfin{inFileName.c_str()};
    TTree *rt = (TTree*)rfin.Get("Events");
    rt->SetBranchAddress("hcalTagsoahcalRecHitCollection_hcalCPURecHitsProducer_recHitsM0LabelOut_RECO.", &wgpu);
    rt->SetBranchAddress("HBHERecHitsSorted_hbheprereco__RECO.", &wcpu);

    // accumulate
    auto const nentries = rt->GetEntries();
    std::cout << ">>> nentries = " << nentries << std::endl;
    for (int ie=0; ie<nentries; ++ie) {
        rt->GetEntry(ie);

        auto const& gpuProduct = wgpu->bareProduct();
        auto const& cpuProduct = wcpu->bareProduct();

        auto const ncpu = cpuProduct.size();
        auto const ngpu = gpuProduct.energy.size();

        if (ngpu != ncpu) {
            std::cerr << "*** mismatch in number of rec hits for event "
                      << ie
                      << std::endl
                      << ">>> ngpu = " << ngpu << std::endl
                      << ">>> ncpu = " << ncpu
                      << std::endl;
        }

        for (uint32_t ich=0; ich<ncpu; ich++) {
            auto const& cpurh = cpuProduct[ich];
            auto const& did = cpurh.id();
            auto iter2idgpu = std::find(
                gpuProduct.did.begin(), gpuProduct.did.end(), did.rawId());

            if (iter2idgpu == gpuProduct.did.end()) {
                std::cerr << "missing " << did << std::endl;
                continue;
            }

            assert(*iter2idgpu == did.rawId());

            auto const ichgpu = iter2idgpu - gpuProduct.did.begin();
            auto const gpu_energy_m0 = gpuProduct.energy[ichgpu];
            auto const cpu_energy_m0 = cpurh.eraw();

            hEnergyM0GPU->Fill(gpu_energy_m0);
            hEnergyM0CPU->Fill(cpu_energy_m0);
            hEnergyM0GPUvsCPU->Fill(cpu_energy_m0, gpu_energy_m0);
        }
    }
        
    {
        TCanvas c{"plots", "plots", 4200, 6200};
        c.Divide(2, 1);
        c.cd(1);
        {
            gPad->SetLogy();
            hEnergyM0CPU->SetLineColor(kBlack);
            hEnergyM0CPU->SetLineWidth(1.);
            hEnergyM0CPU->Draw("");
            hEnergyM0GPU->SetLineColor(kBlue);
            hEnergyM0GPU->SetLineWidth(1.);
            hEnergyM0GPU->Draw("sames");
            gPad->Update();
            auto stats = (TPaveStats*)hEnergyM0GPU->FindObject("stats");
            auto y2 = stats->GetY2NDC();
            auto y1 = stats->GetY1NDC();
            stats->SetY2NDC(y1);
            stats->SetY1NDC(y1 - (y2-y1));
        }
        c.cd(2);
        hEnergyM0GPUvsCPU->Draw("colz");
        c.SaveAs("plots.pdf");
    }

    rfin.Close();
    rfout.Write();
    rfout.Close();
}
