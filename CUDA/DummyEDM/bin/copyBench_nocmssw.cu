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

#define SIZE 10 * 1024 * 1024

float cuda_malloc_test(int size, bool up) {
    cudaEvent_t start, stop;
    int *a, *d_a;
    float elapsedTime;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // allocate host array
    a = new int[size];
    // allocate device array
    cudaMalloc((void**)&d_a, size* sizeof(int));

    // revcord the start of the transfer
    cudaEventRecord(start, 0);
    for (auto i=0; i<100; i++) {
        if (up)
            cudaMemcpy(d_a, a, size*sizeof(*d_a), cudaMemcpyHostToDevice);
        else
            cudaMemcpy(a, d_a, size*sizeof(*d_a), cudaMemcpyDeviceToHost);
    }

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    free(a);
    cudaFree(d_a);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    return elapsedTime;
}

float cuda_host_alloc_test(int size, bool up) {
    cudaEvent_t start, stop;
    int *a, *d_a;
    float elapsedTime;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaHostAlloc((void**)&a, size * sizeof(int), cudaHostAllocDefault);
    cudaMalloc((void**)&d_a, size * sizeof(int));

    cudaEventRecord(start, 0);
    for(auto i = 0; i<100; i++) {
        if (up)
            cudaMemcpy(d_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
        else
            cudaMemcpy(a, d_a, size * sizeof(int), cudaMemcpyDeviceToHost);
    }

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    cudaFreeHost(a);
    cudaFree(d_a);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    return elapsedTime;
}

int main(int argc, char** argv) {
    // start time
    auto startTime = std::chrono::high_resolution_clock::now();
    printf("Hello World\n");

    // get the number of devices
    int numDevices;
    cudaGetDeviceCount(&numDevices);
    PRINT(numDevices);

    cudaDeviceProp prop;
    for (auto i=0 ; i<numDevices; i++) {
        cudaGetDeviceProperties(&prop, i);
        PRINT(prop.name);
        PRINT(prop.totalGlobalMem);
        PRINT(prop.sharedMemPerBlock);
        PRINT(prop.regsPerBlock);
        PRINT(prop.warpSize);
        PRINT(prop.memPitch);
        PRINT(prop.maxThreadsPerBlock);
        PRINT(prop.maxThreadsDim[0]);
        PRINT(prop.maxThreadsDim[1]);
        PRINT(prop.maxThreadsDim[2]);
        PRINT(prop.maxGridSize[0]);
        PRINT(prop.maxGridSize[1]);
        PRINT(prop.maxGridSize[2]);
        PRINT(prop.totalConstMem);
        PRINT(prop.major);
        PRINT(prop.minor);
        PRINT(prop.clockRate);
        PRINT(prop.textureAlignment);
        PRINT(prop.deviceOverlap);
        PRINT(prop.multiProcessorCount);
        PRINT(prop.kernelExecTimeoutEnabled);
        PRINT(prop.integrated);
        PRINT(prop.canMapHostMemory);
        PRINT(prop.computeMode);
        PRINT(prop.maxTexture1D);
        PRINT(prop.maxTexture2D[0]);
        PRINT(prop.maxTexture2D[1]);
        PRINT(prop.maxTexture3D[0]);
        PRINT(prop.maxTexture3D[1]);
        PRINT(prop.maxTexture3D[2]);
//        PRINT(prop.maxTexture2DArray[0]);
//        PRINT(prop.maxTexture2DArray[1]);
//        PRINT(prop.maxTexture2DArray[2]);
        PRINT(prop.concurrentKernels);
    }

    float elapsedTime;
    float MB = (float)100 * SIZE  * sizeof(int)/1024/1024;

    // test with cudaMalloc host -> .device
    elapsedTime = cuda_malloc_test(SIZE, true);
    printf("Time using cudaMalloc: %3.1f ms\n", elapsedTime);
    printf("\tMB/s udring copy host -> device: %3.1f\n", MB/(elapsedTime/1000));

    // test with cudaMalloc device -> host
    elapsedTime = cuda_malloc_test(SIZE, false);
    printf("Time using cudaMalloc: %3.1f ms\n", elapsedTime);
    printf("\tMB/s during copy device -> host: %3.1f\n", MB/(elapsedTime/1000));

    // test with cudaHostAlloc host -> device
    elapsedTime = cuda_host_alloc_test(SIZE, true);
    printf("Time using cudaHostAlloc: %3.1f ms\n", elapsedTime);
    printf("\tMB/s udring copy host -> device: %3.1f\n", MB/(elapsedTime/1000));

    // test with cudaHostAlloc device -> host
    elapsedTime = cuda_host_alloc_test(SIZE, false);
    printf("Time using cudaHostAlloc: %3.1f ms\n", elapsedTime);
    printf("\tMB/s during copy device -> host: %3.1f\n", MB/(elapsedTime/1000));

    // stop time
    auto stopTime = std::chrono::high_resolution_clock::now();
    PRINT((stopTime - startTime).count());

    printf("Goodbye World\n");
}
