#include "CUDA/PixelGPU/interface/kernels.h"
#include "CUDA/PixelGPU/interface/constants.h"
#include <stdio.h>

namespace testpixel {
    //
    // kernel itsel
    //
    __global__ void raw2digi_simple(DataWord* words, DataType* digis, int size) {
        int i = threadIdx.x + blockIdx.x * blockDim.x;
        if (i < size) {
            DataWord word = words[i];

            // it's possible that the word is 0
            digis[i].m_link = (word >> SHIFT_LINK) & MASK_LINK;
            digis[i].m_roc = (word >> SHIFT_ROC) & MASK_ROC;
            digis[i].m_dcol = (word >> SHIFT_DCOL) & MASK_DCOL;
            digis[i].m_pixel = (word >> SHIFT_PIXEL) & MASK_PIXEL;
            digis[i].m_adc = (word) & MASK_ADC;
        }
    }

    //
    // kernel wrapper
    //
    void wrap_raw2digi_simple(cudaStream_t stream,
                              DataWord* word, DataType* digis, int size) {
       int threadsPerBlock {512};
       int blocksPerGrid = (size + threadsPerBlock - 1) / threadsPerBlock;
       raw2digi_simple<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(word, digis, size);
    }
}
