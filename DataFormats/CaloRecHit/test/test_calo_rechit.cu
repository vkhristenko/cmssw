#include <cuda.h>
#include <cuda_runtime.h>

#include <iostream>
#include <cassert>

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"

/*
#ifdef __CUDA_ARCH__
__constant__ uint32_t calo_rechit_masks[] = { 
    0x00000000u,0x00000001u,0x00000003u,0x00000007u,0x0000000fu,0x0000001fu,
    0x0000003fu,0x0000007fu,0x000000ffu,0x000001ffu,0x000003ffu,0x000007ffu,
    0x00000fffu,0x00001fffu,0x00003fffu,0x00007fffu,0x0000ffffu,0x0001ffffu,
    0x0003ffffu,0x0007ffffu,0x000fffffu,0x001fffffu,0x003fffffu,0x007fffffu,
    0x00ffffffu,0x01ffffffu,0x03ffffffu,0x07ffffffu,0x0fffffffu,0x1fffffffu,
    0x3fffffffu,0x7fffffffu,0xffffffffu};
#endif
*/

__global__ void kernel_test_calo_rechit(CaloRecHit* other) {
    CaloRecHit rh{DetId(0), 10, 1, 0, 0};
    other->setEnergy(rh.energy());
    other->setTime(rh.time());
    other->setFlagField(10, 31, 1);
}

void test_calo_rechit() {
    auto check_error = [](auto code) {
        if (code != cudaSuccess) {
            std::cout << cudaGetErrorString(code) << std::endl;
            assert(false);
        }
    };

    CaloRecHit h_rh, h_rh_test{DetId(0), 10, 1, 0, 0};
    h_rh_test.setFlagField(10, 31, 1);
    CaloRecHit *d_rh;

    cudaMalloc((void**)&d_rh, sizeof(CaloRecHit));
    cudaMemcpy(d_rh, &h_rh, sizeof(CaloRecHit), cudaMemcpyHostToDevice);
    kernel_test_calo_rechit<<<1,1>>>(d_rh);
    cudaDeviceSynchronize();
    check_error(cudaGetLastError());
    cudaMemcpy(&h_rh, d_rh, sizeof(CaloRecHit), cudaMemcpyDeviceToHost);

    std::cout << h_rh << std::endl;
    std::cout << h_rh_test << std::endl;
    assert(h_rh.energy() == h_rh_test.energy());
    assert(h_rh.time() == h_rh_test.time());
    assert(h_rh.flags() == h_rh_test.flags());
    assert(h_rh.aux() == h_rh_test.aux());
    assert(h_rh.detid() == h_rh_test.detid());
}

int main(int argc, char** argv) {
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    std::cout << "nDevices = " << nDevices << std::endl;

    if (nDevices > 0)
        test_calo_rechit();

    std::cout << "all good!" << std::endl;
    return 0;
}
