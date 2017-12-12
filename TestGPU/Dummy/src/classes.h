#include "TestGPU/Dummy/interface/Vector.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace {
    struct dictionary {
        testgpu::VectorOfInt vints;
        testgpu::VectorOfFloat vfloats;
        testgpu::VectorOfDouble vdoubles;
    };
}
