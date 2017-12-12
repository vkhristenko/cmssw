#ifndef TestGPU_Dummy_interface_gpu_kernels_h
#define TestGPU_Dummy_interface_gpu_kernels_h

#include <cuda_runtime.h>

namespace testgpu {

//
// simple test 
//
void launch_on_gpu();

//
// a wrapper for the kernel for vector addition.
// NOTE:
// -----
// cudaStream_t is an alias for:
// typedef __device_builtin__ struct CUstream_st *cudaStream_t;
// Therefore no need to pass by reference!
//
// NOTE:
// ------
//
template<typename T>
void wrapperVectorAdd(T*, T*, T*, cudaStream_t, int N=1000);

}

#endif
