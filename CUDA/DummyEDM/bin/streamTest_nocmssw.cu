#include <stdio.h>
#include <future>
#include <thread>
#include <chrono>
#include <iostream>

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

#define N 1024*1024
#define FULL_DATA_SIZE N*20

__global__ void kernel(int *a, int *b, int *c) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < N) {
        int idx1 = (idx+1) % 256;
        int idx2 = (idx+2) % 256;
        float as = (a[idx] + a[idx1] + a[idx2]) / 3.0f;
        float bs = (b[idx] + b[idx1] + b[idx2]) / 3.0f;
        c[idx] = (as + bs) / 2;
    }
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

    cudaEvent_t start, stop;
    float elapsedTime;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    cudaStream_t stream0, stream1;
    cudaStreamCreate(&stream0);
    cudaStreamCreate(&stream1);

    int *h_a, *h_b, *h_c;
    int *d_a0, *d_b0, *d_c0;
    int *d_a1, *d_b1, *d_c1;

    cudaMalloc((void**)&d_a0, N * sizeof(int));
    cudaMalloc((void**)&d_b0, N * sizeof(int));
    cudaMalloc((void**)&d_c0, N * sizeof(int));
    cudaMalloc((void**)&d_a1, N * sizeof(int));
    cudaMalloc((void**)&d_b1, N * sizeof(int));
    cudaMalloc((void**)&d_c1, N * sizeof(int));

    cudaHostAlloc((void**)&h_a, FULL_DATA_SIZE * sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_b, FULL_DATA_SIZE * sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**)&h_c, FULL_DATA_SIZE * sizeof(int), cudaHostAllocDefault);

    for (auto i =0; i<FULL_DATA_SIZE; i++) {
        h_a[i] = i;
        h_b[i] = i*i;
    }

    for (auto i=0; i<FULL_DATA_SIZE; i+=2*N) {
        // copy a for both streams
        cudaMemcpyAsync(d_a0, h_a + i, N * sizeof(int),
                        cudaMemcpyHostToDevice, stream0);
        cudaMemcpyAsync(d_a1, h_a + i + N,
                        N * sizeof(int), cudaMemcpyHostToDevice, stream1);

        // copy b for both streams
        cudaMemcpyAsync(d_b0, h_b + i, N * sizeof(int),
                        cudaMemcpyHostToDevice, stream0);
        cudaMemcpyAsync(d_a1, h_a + i + N,
                        N * sizeof(int), cudaMemcpyHostToDevice, stream1);

        // execute kernels for both streams
        kernel<<<N/256, 256, 0, stream0>>>(d_a0, d_b0, d_c0);
        kernel<<<N/256, 256, 0, stream1>>>(d_a1, d_b1, d_c1);

        // copy c back for both streams
        cudaMemcpyAsync(h_c + i, d_c0, N * sizeof(int),
                        cudaMemcpyDeviceToHost, stream0);
        cudaMemcpyAsync(h_c + i + N, d_c1, N * sizeof(int),
                        cudaMemcpyDeviceToHost, stream1);
    }

    // CPU to wait until GPU has finished
    cudaStreamSynchronize(stream0);
    cudaStreamSynchronize(stream1);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time taken: %3.1f ms\n", elapsedTime);
    for (auto i=0; i<10; i++)
        printf("c[%d] = %d\n", i, h_c[i]);

    cudaFreeHost(h_a);
    cudaFreeHost(h_b);
    cudaFreeHost(h_c);
    cudaFree(d_a0);
    cudaFree(d_b0);
    cudaFree(d_c0);
    cudaFree(d_a1);
    cudaFree(d_b1);
    cudaFree(d_c1);

    cudaStreamDestroy(stream0);
    cudaStreamDestroy(stream1);

    // stop time
    auto stopTime = std::chrono::high_resolution_clock::now();
    PRINT((stopTime - startTime).count());

    printf("Goodbye World\n");
}
