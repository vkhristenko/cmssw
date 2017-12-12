#include <stdio.h>
#include <future>
#include <thread>
#include <chrono>
#include <iostream>

#define N 1000000
#define NUM_THREADS_PER_BLOCK 256
#define NUM_BLOCKS_PER_GRID 1024
//#define NUM_BLOCKS_PER_GRID (N + NUM_THREADS_PER_BLOCK-1) / NUM_THREADS_PER_BLOCK;

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

__global__ 
void dotProduct(float *a, float *b, float *c) {
    // shared memory!
    __shared__ float cache[NUM_THREADS_PER_BLOCK];
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int cacheIndex = threadIdx.x;

    float runningSum = 0;
    while (tid < N) {
        runningSum += a[tid] * b[tid];
        tid += blockDim.x * gridDim.x;
    }

    // store the current running sum for the threads
    cache[cacheIndex] = runningSum;

    // sync all the threads before starting to cooperate
    __syncthreads();

    // reduction
    int i = blockDim.x/2; // number of threads per block
    while (i != 0) {
        if (cacheIndex < i)
            cache[cacheIndex] += cache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }

    if (cacheIndex == 0)
        c[blockIdx.x] = cache[0];
}

#define PRINT(x) \
    std::cout << #x " = " << x << std::endl

void func(const char* ptr) {
    std::cout << "ptr = " << ptr << std::endl;
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

    float h_a[N], h_b[N], h_c[NUM_BLOCKS_PER_GRID];
    float *d_a, *d_b, *d_c;

    cudaMalloc(&d_a, N*sizeof(float));
    cudaMalloc(&d_b, N*sizeof(float));
    cudaMalloc(&d_c, NUM_BLOCKS_PER_GRID*sizeof(float));

    for (auto i=0; i<N; i++) {
        h_a[i] = i;
        h_b[i ] = i*2;
    }

    cudaMemcpy(d_a, h_a, N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b, N*sizeof(float), cudaMemcpyHostToDevice);

    dotProduct<<<NUM_BLOCKS_PER_GRID, NUM_THREADS_PER_BLOCK>>>(d_a, d_b, d_c);

    cudaMemcpy(h_c, d_c, NUM_BLOCKS_PER_GRID*sizeof(float), cudaMemcpyDeviceToHost);
    float sum = 0;
    for (auto i=0; i<NUM_BLOCKS_PER_GRID; i++)
        sum += h_c[i];

    cudaFree(d_c);
    cudaFree(d_a);
    cudaFree(d_b);

#define sum_squares(x) (x*(x+1)*(2*x+1)/6)
    printf("Doues GPU version equal CPU version: %.6g = %.6g\n", sum, 2*sum_squares((float)(N-1)));

    // stop time
    auto stopTime = std::chrono::high_resolution_clock::now();
    PRINT((stopTime - startTime).count());

    printf("Goodbye World\n");
}
