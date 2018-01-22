#ifndef CUDA_DataFormats_DigiFrame_h
#define CUDA_DataFormats_DigiFrame_h

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 

#include <vector>

namespace testpixel {

class DigiFrame {
public:
    CUDA_CALLABLE_MEMBER DigiFrame(int link, int roc, int dcol, int pixel, int adc) :
        m_link(link), m_roc(roc), m_dcol(dcol), m_pixel(pixel), m_adc(adc) 
    {}

    CUDA_CALLABLE_MEMBER DigiFrame() = default;

    // copy ctor
    CUDA_CALLABLE_MEMBER DigiFrame(DigiFrame const&) = default;
    // move ctor
    CUDA_CALLABLE_MEMBER DigiFrame(DigiFrame&& rhs) :
        m_link(std::move(rhs.m_link)),
        m_roc(std::move(rhs.m_roc)),
        m_dcol(std::move(rhs.m_dcol)),
        m_pixel(std::move(rhs.m_pixel)),
        m_adc(std::move(rhs.m_adc))
    {}

    int m_link;
    int m_roc;
    int m_dcol;
    int m_pixel;
    int m_adc;
};

typedef std::vector<testpixel::DigiFrame> DigiFrames;

} // end of namespace testgpu

#endif
