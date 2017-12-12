#include <stdio.h>
#include <future>
#include <thread>
#include <chrono>
#include <iostream>
#include <iterator>
#include <cstring>

#define N 1000000
#define SIZE 100

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

struct cuComplex {
    float r;
    float i;
    __host__ __device__ cuComplex(float a, float b) : r(a), i(b) {}

    __host__ __device__ float magnitude2(void) { return r*r + i*i; }
    __host__ __device__ cuComplex operator*(const cuComplex& a) {
        return cuComplex(r*a.r - i*a.i, i*a.r + r*a.i);
    }
    __host__ __device__ cuComplex operator+(const cuComplex& a) {
        return cuComplex(r+a.r, i+a.i);
    }
};

/*
struct d_cuComplex {
    float r;
    float i;
    __device__ d_cuComplex(float a, float b) : r(a), i(b) {}

    __device__ float magnitude2(void) { return r*r + i*i; }
    __device__ d_cuComplex operator*(const d_cuComplex& a) {
        return d_cuComplex(r*a.r - i*a.i, i*a.r + r*a.i);
    }
    __device__ d_cuComplex operator+(const d_cuComplex& a) {
        return d_cuComplex(r+a.r, i+a.i);
    }
};*/

/*
int julia(int x, int y) {
    const float scale = 1.5;
    float jx = scale * (float)(SIZE/2 - x)/(SIZE/2);
    float jy = scale * (float)(SIZE/2 - y)/(SIZE/2);

//    cuComplex c(-0.4, 0.6);
    cuComplex c(-0.8, 0.156);
    cuComplex a(jx, jy);
    int i = 0;
    for (i=0; i<200; i++) {
        a = a*a + c;
        if (a.magnitude2() > 1000)
            return 0;
    }

    return 1;
}
*/

__host__ __device__ int julia(int x, int y) {
    const float scale = 1.5;
    float jx = scale * (float)(SIZE/2 - x)/(SIZE/2);
    float jy = scale * (float)(SIZE/2 - y)/(SIZE/2);

    cuComplex c(-0.8, 0.156);
    cuComplex a(jx, jy);
    int i = 0;
    for (i=0; i<200; i++) {
        a = a*a + c;
        if (a.magnitude2() > 1000)
            return 0;
    }

    return 1;
}

void kernel(char* ptr) {
    for (int y=0; y<SIZE; y++) {
        for (int x=0; x<SIZE; x++) {
            int offset = x+y*SIZE;

            int juliaValue = julia(x,y);
            ptr[offset] = juliaValue == 1 ? 'x' : ' ';
        }
    }
}

void printImage(char* ptr) {
    for (auto i=0; i<SIZE; i++) {
        char cpyPtr[SIZE+1];
        std::memcpy((void*)cpyPtr, (void*)(ptr + SIZE*i), SIZE);
        cpyPtr[SIZE] = '\0';
        printf("%s\n", cpyPtr);
    }

    printf("\n");
}

__global__ void julia_kernel(char *image) {
    int x = blockIdx.x;
    int y = blockIdx.y;
    int offset = x+y * gridDim.x;

    int juliaValue = julia(x,y);
    image[offset] = juliaValue == 1 ? 'x' : ' ';
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

    char image[SIZE * SIZE];
    kernel(image);
    printImage(image);

    char h_image[SIZE * SIZE];
    char *d_image;
    cudaMalloc((void**)&d_image, SIZE*SIZE);
    dim3 grid(SIZE, SIZE);

    julia_kernel<<<grid, 1>>>(d_image);
    cudaMemcpy(h_image, d_image, SIZE * SIZE, cudaMemcpyDeviceToHost);
    printImage(h_image);

    // stop time
    auto stopTime = std::chrono::high_resolution_clock::now();
    PRINT((stopTime - startTime).count());

    printf("Goodbye World\n");
}
