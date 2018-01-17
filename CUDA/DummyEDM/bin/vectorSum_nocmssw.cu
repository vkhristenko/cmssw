#include <stdio.h>
#include <future>
#include <thread>
#include <chrono>
#include <iostream>
#include <iterator>
#include <cstring>

#define N 1000000
#define SIZE 100

__constant__ int factor = 1;

//
// NOTE: while loop is for the case when number of elements in the array exceeds the 
// number of blocks possible total in device to be launched
//
__global__ 
void vectorAdd(int *a, int *b, int *c) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    while (i < N) {
        c[i] = factor*(a[i] + b[i]);
        i += blockDim.x * gridDim.x;
    }
}

__global__
void matrixAdd(int **a,int **b, int**c) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    c[i][j] = a[i][j] + b[i][j];
}

#define PRINT(x) \
    std::cout << #x " = " << x << std::endl

void printImage(char* ptr) {
    for (auto i=0; i<SIZE; i++) {
        char cpyPtr[SIZE+1];
        std::memcpy((void*)cpyPtr, (void*)(ptr + SIZE*i), SIZE);
        cpyPtr[SIZE] = '\0';
        printf("%s\n", cpyPtr);
    }

    printf("\n");
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

    int h_a[N], h_b[N], h_c[N];
    int *d_a, *d_b, *d_c;

    cudaMalloc(&d_a, N*sizeof(int));
    cudaMalloc(&d_b, N*sizeof(int));
    cudaMalloc(&d_c, N*sizeof(int));

    for (auto i=0; i<N; i++) {
        h_a[i ] = i;
        h_b[i] = i*i;
    }

    cudaMemcpy(d_a, h_a, N*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b, N*sizeof(int), cudaMemcpyHostToDevice);


    vectorAdd<<<(N+127)/128, 128>>>(d_a, d_b, d_c);

    cudaMemcpy(h_c, d_c, N*sizeof(int), cudaMemcpyDeviceToHost);

    for (auto i=0; i<N; i++) {
        if (i < 10)
            printf("c[%d] = %d\n", i, h_c[i]);
    }

    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    // stop time
    auto stopTime = std::chrono::high_resolution_clock::now();
    PRINT((stopTime - startTime).count());

    printf("Goodbye World\n");
}
