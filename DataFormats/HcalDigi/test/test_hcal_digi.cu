#include <cuda_runtime.h>
#include <cuda.h>

#include <iostream>
#include <assert.h>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"

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
    test_hcal_qiesample();
    test_hcal_hbhedf();
}
