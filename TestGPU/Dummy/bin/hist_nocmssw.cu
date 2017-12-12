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

#define NUM_BINS 256
#define SIZE 100 * 1024 * 1024
#define rnd(x) (x*rand() / RAND_MAX)
#define CPU 0
#undef CPU

unsigned char* generate_random_buffer(int const size) {
    unsigned char * buffer = new unsigned char[size];
    for (auto i=0; i<size; i++) {
        buffer[i] = rand() % 256;
    }

    return buffer;
}

__global__ void kernel(unsigned char *buffer, int const size, unsigned int *histo) {
    __shared__ unsigned int temp[NUM_BINS];
    temp[threadIdx.x] = 0;
    __syncthreads();

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    while (i < size) {
        atomicAdd(&temp[buffer[i]], 1);
        i += stride;
    }
    __syncthreads();

    atomicAdd(&(histo[threadIdx.x]), temp[threadIdx.x]);

    /*
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    while (i < size) {
        atomicAdd(&(histo[buffer[i]]), 1);
        i += stride;
    }
    */
}

int main(int argc, char** argv) {
    // start time
    auto startTime = std::chrono::high_resolution_clock::now();
    printf("Hello World\n");

    // stop time
    auto stopTime = std::chrono::high_resolution_clock::now();
    PRINT((stopTime - startTime).count());

    unsigned char *buffer = generate_random_buffer(SIZE);
    unsigned int histo[NUM_BINS];

#ifdef CPU
    for (auto i=0; i<NUM_BINS; i++)
        histo[i] = 0;

    for (auto i=0; i<SIZE; i++)
        histo[buffer[i]]++;

    long histoCount = 0;
    for (auto i=0; i<NUM_BINS; i++)
        histoCount += histo[i];
    printf("Histogram Sum: %ld\n", histoCount);
#else
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    int blocks = prop.multiProcessorCount;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    unsigned char *d_buffer;
    unsigned int *d_histo;
    cudaMalloc((void**)&d_buffer, SIZE);
    cudaMemcpy(d_buffer, buffer, SIZE, cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_histo, NUM_BINS * sizeof(unsigned int));
    cudaMemset(d_histo, 0, NUM_BINS * sizeof(unsigned int));

    // kernel
    kernel<<<blocks*2, 256>>>(d_buffer, SIZE, d_histo);

    cudaMemcpy(histo, d_histo, NUM_BINS * sizeof(unsigned int), cudaMemcpyDeviceToHost);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);

    printf("Time to generate %3.1f ms\n", elapsedTime);

    long histoCount = 0;
    for (auto i=0; i<NUM_BINS; i++)
        histoCount += histo[i];
    printf("Hisotgram Sum: %ld\n", histoCount);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    cudaFree(d_histo);
    cudaFree(d_buffer);
#endif

    free( buffer);

    printf("Goodbye World\n");
}
