#include "TestGPU/Dummy/interface/gpu_kernels.h"

#include <stdio.h>

namespace testgpu {

//
// Vector Addition Kernel
//
template<typename T>
__global__
void vectorAdd(T *a, T *b, T *c) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    for (size_t j=0; j<1000000; j++)
        c[i] = a[i] + b[i];
}

//
// Vector Addition Kernel Wrapper
//
template<typename T>
void wrapperVectorAdd(T* d_a, T* d_b, T* d_c, cudaStream_t stream, int N) {
    int threadsPerBlock {256};
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    vectorAdd<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_a, d_b, d_c);
}

//
// Macros to simplify the template instantiation
//
#define WRAPPERVECTORADD(TYPE) \
    template void wrapperVectorAdd<TYPE>(TYPE*, TYPE*, TYPE*, cudaStream_t, int)

//
// NOTE:
// -----
// We have to instantiate tempaltes explicitly given that kernels will be compiled 
// as separate compilation units and linked afterwards.
// 
//
WRAPPERVECTORADD(int);
WRAPPERVECTORADD(float);
WRAPPERVECTORADD(double);
WRAPPERVECTORADD(long);

//
// Standalone function that allocates/copies/launches/frees and prints the results
//
void launch_on_gpu() {
    int const NUM_VALUES = 10000;
    printf("start launch_on_gpu\n");
    int h_a[NUM_VALUES], h_b[NUM_VALUES], h_c[NUM_VALUES];
    for (auto i=0; i<NUM_VALUES; i++) {
        h_a[i] = i;
        h_b[i] = i*i;
    }

    int *d_a, *d_b, *d_c;
    cudaMalloc(&d_a, NUM_VALUES*sizeof(int));
    cudaMalloc(&d_b, NUM_VALUES*sizeof(int));
    cudaMalloc(&d_c, NUM_VALUES*sizeof(int));

    cudaMemcpy(d_a, h_a, NUM_VALUES*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b, NUM_VALUES*sizeof(int), cudaMemcpyHostToDevice);

    int threadsPerBlock {256};
    int blocksPerGrid = (NUM_VALUES + threadsPerBlock - 1) / threadsPerBlock;
    vectorAdd<<<blocksPerGrid, threadsPerBlock>>>(d_a, d_b, d_c);

    cudaMemcpy(h_c, d_c, NUM_VALUES*sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    for (auto i=0; i<10; i++) {
        printf("c[%d] = %d\n", i, h_c[i]);
    }

    printf("\n");
    printf("stop launch_on_gpu\n");
}

}
