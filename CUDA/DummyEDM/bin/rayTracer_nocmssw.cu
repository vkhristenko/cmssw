#include <stdio.h>
#include <future>
#include <thread>
#include <chrono>
#include <iostream>
#include <cstring>
#include <stdio.h>

#define N 1000000
#define SIZE 100

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                          \
     cudaError_t e=cudaGetLastError();                                 \
     if(e!=cudaSuccess) {                                              \
            printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
            exit(0); \
          }                                                                 \
}

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

#define INF 2e10f

struct Sphere {
    float r,b,g;
    float radius;
    float x,y,z;

    __device__ float hit(float ox, float oy, float *n) {
        float dx = ox - x;
        float dy = oy - y;
        printf("%f %f %f %f\n", x, y, dx, dy);
        printf("%f %f %f\n", dx*dx, dy*dy, radius*radius);
        if (dx*dx + dy*dy < radius*radius) {
            float dz = sqrtf(radius*radius - dx*dx - dy*dy);
            *n = dz / sqrtf(radius*radius);
            printf("n = %f\n", *n);
            return dz + z;
        }
        return -INF;
    }
};

#define rnd(x) (x*rand() / RAND_MAX)
#define SPHERES 20
__constant__ Sphere s[SPHERES];

__global__ void kernel(char *image) {
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    int offset = x+ y*blockDim.x*gridDim.x;

    float ox = (x - SIZE/2);
    float oy = (y - SIZE/2);

    float r=0, g=0, b=0;
    float maxz = -INF;
    for (int i=0; i<SPHERES; i++) {
        float n;
        float t = s[i].hit(ox, oy, &n);
        printf("t = %f\n", t);
        if (t > maxz) {
            float fscale = n;
            r = s[i].r * fscale;
            g = s[i].g * fscale;
            b = s[i].b * fscale;
            printf("%f\n %f\n", s[i].r, fscale);
        }
    }

    image[offset] = r>0 && g>0 && b>0 ? 'x' : ' '; 
}

void printImage(char *ptr) {
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
    cudaCheckError();
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

    // capture the start/stop time
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    cudaCheckError();

    char *d_image;
    cudaMalloc(&d_image, SIZE*SIZE*sizeof(char));
    cudaCheckError();

    Sphere *h_s = new Sphere[SPHERES];
    for (auto i=0; i<SPHERES; i++) {
        h_s[i].r = rnd(1.0f);
        h_s[i].g = rnd(1.0f);
        h_s[i].b = rnd(1.0f);

        h_s[i].x = rnd((float)SIZE) - SIZE/2;
        h_s[i].y = rnd((float)SIZE) - SIZE/2;
        h_s[i].z = rnd((float)SIZE) - SIZE/2;
        h_s[i].radius = rnd((float)SIZE/10) + 2;
    }

    cudaMemcpyToSymbol(s, h_s, sizeof(Sphere) * SPHERES);
    cudaCheckError();
//    cudaMemcpy(s, h_s, sizeof(Sphere) * SPHERES, 
//            cudaMemcpyHostToDevice);
    delete [] h_s;

    dim3 grids(SIZE/16, SIZE/16);
    dim3 threads(16,16);
    kernel<<<grids, threads>>>(d_image);
    cudaCheckError();

    char h_image[SIZE*SIZE];
    cudaMemcpy(h_image, d_image, SIZE*SIZE*sizeof(char),
            cudaMemcpyDeviceToHost);
    cudaCheckError();

    cudaFree(d_image);
    cudaCheckError();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);

    float elapsed;
    cudaEventElapsedTime(&elapsed, start, stop);
    printf("elapsed = %f\n", elapsed);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    printImage(h_image);

    // stop time
    auto stopTime = std::chrono::high_resolution_clock::now();
    PRINT((stopTime - startTime).count());

    printf("Goodbye World\n");
}
