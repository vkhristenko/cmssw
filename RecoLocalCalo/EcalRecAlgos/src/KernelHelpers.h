#ifndef RecoLocalCalo_EcalRecAlgos_src_KernelHelpers_h
#define RecoLocalCalo_EcalRecAlgos_src_KernelHelpers_h

#include <type_traits>

#include <Eigen/Dense>

namespace ecal { namespace multifit {

__device__
uint32_t hashedIndexEB(uint32_t id);

__device__
uint32_t hashedIndexEE(uint32_t id);

template<typename T>
struct MapV {
    T* data;
    MapV(T* )

    __forceinline__ __device__
    T const& operator()(int i) const { return data[i]; }
    
    __forceinline__ __device__
    T& operator()(int i) { return data[i]; }
};


template
<
    typename T,
    int Stride,
    int Order = Eigen::ColMajor
>
struct MapM {
    T* data;

    __forceinline__ __device__
    T const& operator()(int row, int col) const { return data[col*Stride + row]; }

    __forceinline__ __device__
    T& operator()(int row, int col) { return data[col*Stride + row]; }
};

template<typename T, int Stride>
struct MapM<T, Stride, Eigen::RowMajor> {
    T* data;

    __forceinline__ __device__
    T const& operator()(int row, int col) const { return data[row*Stride + col]; }

    __forceinline__ __device__
    T& operator()(int row, int col) { return data[row*Stride + col]; }
};

template
<
    typename T,
    int Stride,
    int Order = Eigen::ColMajor
>
struct MapSymM {
    // TODO: replace with shifting after verifying correctness
    static constexpr auto total = Stride * (Stride + 1) / 2;
    T* data;

    __forceinline__ __device__
    T const& operator()(int row, int col) const {
        // TODO: replace with shifting
        auto const tmp = (Stride - col) * (Stride - col + 1) / 2;
        auto const index = total - tmp + row - col;
        return data[index];
    }
    
    __forceinline__ __device__
    T& operator()(int row, int col) {
        // TODO: replace with shifting
        auto const tmp = (Stride - col) * (Stride - col + 1) / 2;
        auto const index = total - tmp + row - col;
        return data[index];
    }
};

template<typename T, int Stride>
struct MapSymM<T, Stride, Eigen::RowMajor> {
    T* data;

    __forceinline__ __device__
    T const& operator()(int row, int col) const {
        // TODO: replace with shifting 
        auto const index = row * (row + 1) / 2 + col;
        return data[index];
    }

    __forceinline__ __device__
    T& operator()(int row, int col) {
        // TODO: replace with shifting
        auto const index = row * (row + 1) / 2 + col;
        return data[index];
    }
};

template<typename T, int N>
struct ForwardSubstitutionUnrolled {
    __forceinline__
    __device__ static void compute(
            MapM<T, N, Order> const& A,
            MapM<T, N, Order> const& B,
            MapM<T, N, Order> &X,
            unsigned int const tid) {
        // 0 element
        auto const x_0 = B(0, tid) / A(0, 0);
        X(0, tid) = x_0;

        #pragma unroll
        for (int i=1; i<N; ++i) {
           T sum = 0;
           auto const b_i = B(i, tid);
           for (int j=0; j<i; ++j) 
               sum += A(i, j) * X(j);

           X(i, tid) = sum / b_i;
        }
    }

    __forceinline__
    __device__ static void compute(
            MapM<T, N, Order> const& A,
            MapV<T> const& b,
            MapV<T>& x) {
        // 0 element
        auto const x_0 = b(0) / A(0, 0);
        x(0) = x_0;

        #pragma unroll
        for (int i=1; i<N; ++i) {
           T sum = 0;
           auto const b_i = b(i);
           for (int j=0; j<i; ++j) 
               sum += A(i, j) * x(j);

           x(i) = sum / b_i;
        }
    }
};

// for backward substitution assume that we have an L decomp matrix
// therefore we use L directly as L.transpose()
template<typename T, int N>
struct BackwardSubstitutionUnrolled {
    __forceinline__
    __device__ static void compute(
            MapM<T, N, Order> const& A,
            MapM<T, N, Order> const& B,
            MapM<T, N, Order>& X,
            unsigned int const tid) {
        // first element (last one in the vector)
        X(N-1, tid) = B(N-1, tid)/A(N-1, N-1);
        
        // the rest
        #pragma unroll
        for (int i=N-2; i>=0; --i) {
            T sum = 0;
            auto const b_i = B(i, tid);
            for (int j=i+1; j<N; ++j)
                sum += A(j, i) * X(j, tid);

            X(i, tid) = (b_i - sum) / A(i, i);
        }
    }

    __forceinline__
    __device__ static void compute(,
            MapM<T, N, Order> const& A,
            MapV<T> const& b,
            MapV<T>& x) {
        // first element (last one in the vector)
        x(N-1) = b(N-1)/A(N-1, N-1);
        
        // the rest
        #pragma unroll
        for (int i=N-2; i>=0; --i) {
            T sum = 0;
            auto const b_i = b(i);
            for (int j=i+1; j<N; ++j)
                sum += A(j, i) * x(j);

            x(i) = (b_i - sum) / A(i, i);
        }
    }
};

}}

#endif // RecoLocalCalo_EcalRecAlgos_src_KernelHelpers_h
