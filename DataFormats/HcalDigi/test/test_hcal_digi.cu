#include <cuda_runtime.h>
#include <cuda.h>

#include <iostream>
#include <assert.h>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"

__global__ void test_hcal_qiesample() {
    printf("hello\n");
    printf("%f %f %f\n", nominal_adc2fc[0], nominal_adc2fc[1], nominal_adc2fc[2]);
}

__global__ void test_hcal_hbhedf() {

}

int main(int argc, char** argv) {
    test_hcal_qiesample<<<1,1>>>();
    cudaDeviceSynchronize();
}
