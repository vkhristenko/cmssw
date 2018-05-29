#include <cuda_runtime.h>
#include <cuda.h>

#include <iostream>
#include <assert.h>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

__global__ void kernel_test_hcal_qiesample(HcalQIESample* sample, uint16_t value) {
    printf("kernel: testing hcal qie sampel\n");
    printf("%f %f %f\n", nominal_adc2fc[0], nominal_adc2fc[1], nominal_adc2fc[2]);

    HcalQIESample tmp{value};
    *sample = tmp;
}

__global__ void kernel_test_hcal_hbhedf(HBHEDataFrame *df) {
    printf("kernel: testing hcal hbhe dataframe\n");
    df->setSize(10);
    for (auto i=0; i<10; i++)
        df->setSample(i, HcalQIESample(100));
    df->setReadoutIds(HcalElectronicsId(100));
}

void test_hcal_qiesample() {
    HcalQIESample h_sample, h_test_sample0{100}, h_test_sample1;
    HcalQIESample *d_sample;

    cudaMalloc((void**)&d_sample, sizeof(HcalQIESample));
    cudaMemcpy(d_sample, &h_sample, sizeof(HcalQIESample), cudaMemcpyHostToDevice);
    kernel_test_hcal_qiesample<<<1,1>>>(d_sample, 100);
    cudaMemcpy(&h_sample, d_sample, sizeof(HcalQIESample), cudaMemcpyDeviceToHost);

    assert(h_sample() == h_test_sample0());
    assert(h_sample() != h_test_sample1());
}

template<typename TDF>
__global__ void kernel_test_hcal_digis(TDF *pdfs, uint32_t* out) {
    int id = threadIdx.x;
    uint32_t sum = 0;
    for (auto i=0; i<10; i++)
        sum += pdfs[id].sample(i).raw();
    out[id] = sum;
}

template<typename TDF>
void test_hcal_digis() {
    constexpr int n = 10;
    edm::SortedCollection<TDF> coll{n};
    TDF *d_dfs;
    uint32_t *d_out;
    uint32_t h_out[n], h_test_out[n];
    for (auto i=0; i<n; i++) {
        TDF &df = coll[i];
        df.setSize(10);
        h_test_out[i] = 0;
        uint32_t test = 0;
        for (auto j=0; j<10; j++) {
            df.setSample(j, HcalQIESample(100));
            h_test_out[i] += df.sample(j).raw();
            test += df.sample(j).raw();
        }
    }

    cudaMalloc((void**)&d_dfs, n * sizeof(TDF));
    cudaMalloc((void**)&d_out, n * sizeof(uint32_t));
    cudaMemcpy(d_dfs, &(*coll.begin()), n * sizeof(TDF), 
        cudaMemcpyHostToDevice);
    kernel_test_hcal_digis<<<1, n>>>(d_dfs, d_out);
    cudaMemcpy(&h_out, d_out, n * sizeof(uint32_t), cudaMemcpyDeviceToHost);

    std::cout << "collection size = " << coll.size() << std::endl;

    // comparison
    for (auto i=0; i<n; i++) {
        std::cout << h_out[i] << " == " << h_test_out[i] << std::endl;
        assert(h_out[i] == h_test_out[i]);
    }
}

void test_hcal_hbhedf() {
    HBHEDataFrame h_df, h_test_df;
    HBHEDataFrame *d_df;

    h_test_df.setSize(10);
    for (auto i=0; i<10; i++)
        h_test_df.setSample(i, HcalQIESample(100));
    h_test_df.setReadoutIds(HcalElectronicsId(100));

    cudaMalloc((void**)&d_df, sizeof(HBHEDataFrame));
    cudaMemcpy(d_df, &h_df, sizeof(HBHEDataFrame), cudaMemcpyHostToDevice);
    kernel_test_hcal_hbhedf<<<1,1>>>(d_df);
    cudaMemcpy(&h_df, d_df, sizeof(HBHEDataFrame), cudaMemcpyDeviceToHost);

    assert(h_df.size() == h_test_df.size());
    assert(h_df.elecId() == h_test_df.elecId());
    for (auto i=0; i<10; i++)
        assert(h_df[i].raw() == h_test_df[i].raw());
}

int main(int argc, char** argv) {
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    std::cout << "nDevices = " << nDevices << std::endl;

    if (nDevices > 0) {
        test_hcal_qiesample();
        test_hcal_hbhedf();
        test_hcal_digis<HBHEDataFrame>();
        test_hcal_digis<HFDataFrame>();
        test_hcal_digis<HODataFrame>();
    }
}
