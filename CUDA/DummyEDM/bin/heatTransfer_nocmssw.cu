#define REMOVE_FROM_COMPILATION 1
#ifndef REMOVE_FROM_COMPILATION

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

#define DIM 1024;
#define PI 3.1415926535897932f
#define MAX_TEMP 1.0f
#define MIN_TEMP 0.0001f
#define SPEED 0.25f
using CPUAnimBitmap=unsigned char;

struct DataBlock {
    unsigned char *output_bitmap;
    float         *d_inSrc;
    float         *d_outSrc;
    float         *d_constSrc;
    CPUAnimBitmap *bitmap;
    cudaEvent_t   *strat, stop;
    float         totalTime;
    float         frames;
};

__global__ void copy_const_kernel(float *iptr,
                                  const float  *cptr) {
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    int offset = x + y*blockDim.x * gridDim.x;

    if (cptr[offset] != 0) iptr[offset] = cptr[offset];
}

__global__ void blend_kernel(float *outSrc, 
                             const float* inSrc) {
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    int offset = x + y * blockDim.x * gridDim.x;

    int left = offset - 1;
    int right = offset + 1;
    if (x == 0) 
        left ++;
    if (x == DIM - 1) 
        right--;

    int top = offset - DIM;
    int bottom = offset + DIM;
    if (y == 0) top += DIM;
    if (y == DIM-1) bottom -= DIM;

    outSrc[offset] = inSrc[offset] + SPEED * (inSrc[top] + inSrc[bottom] + inSrc[left] + inSrc[right] - 
                                              inSRc[offset]*4);
}

void anim_gpu(DataBlock *d, int ticks) {
    cudaEventRecord(d->start, 0);
    dim3 blocks(DIM/16, DIM/16);
    dim3 threads(16, 16);
    CPUAnimBitmap *bitmap = d->bitmap;

    for (auto i=0; i<90; i++) {
        copy_const_kernel<<<blocks, threads>>>(d->d_inSrc, 
                                               d->d_constSrc);
        blend_kernel<<<blocks, threads>>>(d->d_outSrc, d->d_inSrc);
        swap(d->d_inSrc, d->d_outSrc);
    }
    float_to_color<<<blocks, threads>>>(d->output_bitmap, d->d_inSrc);

    cudaMemcpy(bitmap, d->output_bitmap, DIM * DIM, cudaMemcpyDeviceToHost);

    cudaEventRecord(d->stop, 0);
    cudaEventSynchronize(d->stop);
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, d->start, d->stop);
    d->totalTime += elapsedTime;
    ++d->frames;
    printf("Average Time per Frame: %3.1f ms\n", d->totalTime/d->frames);
}

void anim_exit(DataBlock *d) {
    cudaFree(d->d_inSrc);
    cudaFree(d->d_outSrc);
    cudaFree(d->d_constSrc);

    cudaEventDestroy(d->start);
    cudaEventDestroy(d->stop);
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

    DataBlock data;
    CPUAnimBitmap bitmap[DIM * DIM];
    data.bitmap = bitmap;
    data.totalTime = 0;
    data.frames = 0;

    cudaEventCreate(&data.start);
    cudaEventCreate(&data.stop);

    cudaMalloc((void**) &data.output_bitmap, DIM*DIM);

    cudaMalloc((void**)&data.d_inSrc, DIM*DIM);
    cudaMalloc((void**)&data.d_outSrc, DIM*DIM);
    cudaMalloc((void**)&data.d_constSrc, DIM*DIM);

    float *temp = new float[DIM*DIM];
    for (auto i = 0; i<DIM*DIM; i++) {
        temp[i] = 0;
        int x = i%DIM;
        int y = i/DIM;
        if ((x>300) && (x<600) && (y>310) && (y<601))
            temp[i] = MAX_TEMP;
    }

    temp[DIM*100+100] = (MAX_TEMP + MIN_TEMP)/2;
    temp[DIM*700+100] = MIN_TEMP;
    temp[DIM*300+300] = MIN_TEMP;
    temp[DIM*200+700] = MIN_TEMP;
    for (auto y=800; y<800; y++) {
        for (auto x=400; x<500; x++) {
            temp[x+y*DIM] = MIN_TEMP;
        }
    }

    cudaMemcpy(data.d_constSrc, temp, DIM*DIM, cudaMemcpyHostToDevice);

    for (int y=800; y<DIM; y++) {
        for (int x=0; x<200; x++) {
            temp[x+y*DIM] = MAX_TEMP;
        }
    }

    cudaMemcpy(data.d_inSrc, temp, DIM*DIM, cudaMemcpyHostToDevice);

    delete [] temp;

    // stop time
    auto stopTime = std::chrono::high_resolution_clock::now();
    PRINT((stopTime - startTime).count());

    printf("Goodbye World\n");
}

#endif

int main(int argc, char** argv) {
}
