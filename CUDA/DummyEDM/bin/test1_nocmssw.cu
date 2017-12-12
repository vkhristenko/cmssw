#include <stdio.h>
#include <future>
#include <thread>
#include <chrono>
#include <iostream>

#define N 1000000

__constant__ int factor = 0;

__global__ 
void vectorAdd(int *a, int *b, int *c) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    c[i] = factor*(a[i] + b[i]);
}

__global__
void matrixAdd(int **a,int **b, int**c) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    c[i][j] = a[i][j] + b[i][j];
}

#define PRINT(x) \
    std::cout << #x " = " << x << std::endl

void func(const char* ptr) {
    std::cout << "ptr = " << ptr << std::endl;
}

int main(int argc, char** argv) {
    auto startTime = std::chrono::high_resolution_clock::now();
    printf("Hello World\n");
    int h_a[N], h_b[N], h_c[N];
    for (auto i=0; i<N; i++) {
        h_a[i] = i;
        h_b[i] = i*i;
    }

    int lfactor = 10;
    cudaMemcpyToSymbol(factor, &lfactor, sizeof(int));
    printf("lfactory = %d\n", lfactor);

    // allocation
    int *d_a, *d_b, *d_c;
    cudaMalloc(&d_a, N*sizeof(int));
    cudaMalloc(&d_b, N*sizeof(int));
    cudaMalloc(&d_c, N*sizeof(int));

    std::chrono::duration<double, std::milli> timeToCopyH2D = 
        std::chrono::high_resolution_clock::now() - startTime;
    printf("Copy from Host to Device\n");
    // copy from host to device memory
    cudaMemcpy(d_a, h_a, N*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b, N*sizeof(int), cudaMemcpyHostToDevice);
    std::chrono::duration<double, std::milli> timeToAfterCopyH2D = 
        std::chrono::high_resolution_clock::now() - startTime;

    // vector addition
    int threadsPerBlock(256);
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock; 
    printf("Start the Kernel\n");
    std::chrono::duration<double, std::milli> timeToJustBeforeKernelH2D = 
        std::chrono::high_resolution_clock::now() - startTime;
    vectorAdd<<<blocksPerGrid, threadsPerBlock>>>(d_a, d_b, d_c);
    std::chrono::duration<double, std::milli> timeToJustAfterKernel = 
        std::chrono::high_resolution_clock::now() - startTime;
    printf("Finished the Kernel\n");

    // copy the result output
    printf("Copy from Device to Host\n");
    std::chrono::duration<double, std::milli> timeToBeforeCopyD2H = 
        std::chrono::high_resolution_clock::now() - startTime;
    cudaMemcpy(h_c, d_c, N*sizeof(int), cudaMemcpyDeviceToHost);
    std::chrono::duration<double, std::milli> timeToAfterCopyD2H = 
        std::chrono::high_resolution_clock::now() - startTime;
    printf("Finished copying...\n");

    PRINT(timeToCopyH2D.count());
    PRINT(timeToAfterCopyH2D.count());
    PRINT(timeToJustBeforeKernelH2D.count());
    PRINT(timeToJustAfterKernel.count());
    PRINT(timeToBeforeCopyD2H.count());
    PRINT(timeToAfterCopyD2H.count());

    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    for (auto i=0; i<10; i++) {
        printf("c[%d] = %d\n", i, h_c[i]);
    }
    printf("\n");
    
    // matrix addition
    /*
    dim3 threadsPerBlock(10, 10);
    dim3 numBlocks(N / threadsPerBlock.x, N / threadsPerBlock.y);
    maxtrixAdd<<<numBlocks, threadsPerBlock>>>(dma, dmb, dmc);
    */

    printf("Goodbye World\n");
}
