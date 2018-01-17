#ifndef CUDA_DummyEDM_Vector_h
#define CUDA_DummyEDM_Vector_h

#include <vector>

namespace testgpu {

template<typename T>
class Vector {
public:
    std::vector<T> m_values;
};

typedef testgpu::Vector<int> VectorOfInt;
typedef testgpu::Vector<float> VectorOfFloat;
typedef testgpu::Vector<double> VectorOfDouble;

}

#endif
