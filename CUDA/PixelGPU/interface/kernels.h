#ifndef CUDA_PixelGPU_interface_kernels_h
#define CUDA_PixelGPU_interface_kernels_h

#include "CUDA/DataFormats/interface/DigiFrame.h"
#include <cuda.h>
#include <cuda_runtime.h>

namespace testpixel {
    using Wor64 = unsigned long;
    using Word32 = unsigned int;
    using DataWord = Word32;
    using DataType = DigiFrame;
    
    //
    // kernel wrapper
    //
    void wrap_raw2digi_simple(cudaStream_t, DataWord*, DataType*, int);
}

#endif
