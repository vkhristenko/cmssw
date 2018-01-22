#ifndef CUDA_PixelGPU_interaface_constants_h
#define CUDA_PixelGPU_interface_constants_h

namespace testpixel {
    constexpr int BITS_LINK = 6;
    constexpr int BITS_ROC = 5;
    constexpr int BITS_DCOL = 5;
    constexpr int BITS_PIXEL = 8;
    constexpr int BITS_ADC = 8;

    constexpr int SHIFT_LINK = BITS_ROC + BITS_DCOL + BITS_PIXEL + BITS_ADC;
    constexpr int SHIFT_ROC = BITS_DCOL + BITS_PIXEL + BITS_ADC;
    constexpr int SHIFT_DCOL = BITS_PIXEL + BITS_ADC;
    constexpr int SHIFT_PIXEL = BITS_ADC;

    constexpr int MASK_LINK = 0x0000003F;
    constexpr int MASK_ROC = 0x0000001F;
    constexpr int MASK_DCOL = MASK_ROC;
    constexpr int MASK_PIXEL = 0x000000FF;
    constexpr int MASK_ADC = MASK_PIXEL;
}

#endif
