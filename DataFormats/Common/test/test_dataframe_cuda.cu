#include <cuda.h>
#include <cuda_runtime.h>

#include "DataFormats/Common/interface/DataFrame.h"
#include "DataFormats/Common/interface/DataFrameContainer.h"

__global__ void kernel_test_df() {
    printf("kernel: test_df\n");
    edm::DataFrame df0;
    edm::DataFrame df1 {0, nullptr, 0};
}

__global__ void kernel_test_gen_df(unsigned int* ids, unsigned short datas) {

}

void test_df() {
    edm::DataFrame df0;
    edm::DataFrame df1{0, nullptr, 0};
    edm::DataFrameContainer dfc{10, 2, 10};
    for (auto i=0; i<10; i++) {
        auto tmp = dfc[i];
        for (auto j=0; j<10; j++)
            tmp[j] = 100;
    }

    kernel_test_df<<<1,1>>>();
    cudaDeviceSynchronize();
}

int main(int argc, char** argv) {
    test_df();
}
