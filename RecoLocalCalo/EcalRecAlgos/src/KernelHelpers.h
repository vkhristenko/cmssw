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
    using type = T;
    using base_type = typename std::remove_const<type>::type;

    type* data;
    __forceinline__ __device__
    MapV(type* data) : data{data} {}

    __forceinline__ __device__
    base_type const& operator()(int i) const { return data[i]; }
    
    template<typename U = T>
    __forceinline__ __device__
    typename std::enable_if<std::is_same<base_type, U>::value, base_type>::type&
    operator()(int i) { return data[i]; }
};


template
<
    typename T,
    int Stride,
    int Order = Eigen::ColMajor
>
struct MapM {
    using type = T;
    using base_type = typename std::remove_const<type>::type;

    type* data;
    __forceinline__ __device__
    MapM(type* data) : data{data} {}

    __forceinline__ __device__
    base_type const& operator()(int row, int col) const { return data[col*Stride + row]; }

    template<typename U = T>
    __forceinline__ __device__
    typename std::enable_if<std::is_same<base_type, U>::value, base_type>::type&
    operator()(int row, int col) { return data[col*Stride + row]; }
};

template<typename T, int Stride>
struct MapM<T, Stride, Eigen::RowMajor> {
    using type = T;
    using base_type = typename std::remove_const<type>::type;

    type* data;
    __forceinline__ __device__
    MapM(type* data) : data{data} {}

    __forceinline__ __device__
    base_type const& operator()(int row, int col) const { return data[row*Stride + col]; }

    template<typename U = T>
    __forceinline__ __device__
    typename std::enable_if<std::is_same<base_type, U>::value, base_type>::type&
    operator()(int row, int col) { return data[row*Stride + col]; }
};

template
<
    typename T,
    int Stride,
    int Order = Eigen::ColMajor
>
struct MapSymM {
    using type = T;
    using base_type = typename std::remove_const<type>::type;

    // TODO: replace with shifting after verifying correctness
    static constexpr auto total = Stride * (Stride + 1) / 2;
    T* data;
    __forceinline__ __device__
    MapSymM(T* data) : data{data} {}

    __forceinline__ __device__
    T const& operator()(int row, int col) const {
        // TODO: replace with shifting
        auto const tmp = (Stride - col) * (Stride - col + 1) / 2;
        auto const index = total - tmp + row - col;
        return data[index];
    }
    
    template<typename U = T>
    __forceinline__ __device__
    typename std::enable_if<std::is_same<base_type, U>::value, base_type>::type&
    operator()(int row, int col) {
        // TODO: replace with shifting
        auto const tmp = (Stride - col) * (Stride - col + 1) / 2;
        auto const index = total - tmp + row - col;
        return data[index];
    }
};

template<typename T, int Stride>
struct MapSymM<T, Stride, Eigen::RowMajor> {
    using type = T;
    using base_type = typename std::remove_const<type>::type;

    static constexpr auto total = Stride * (Stride + 1) / 2;
    T* data;
    __forceinline__ __device__
    MapSymM(T* data) : data{data} {}

    __forceinline__ __device__
    T const& operator()(int row, int col) const {
        // TODO: replace with shifting 
        auto const index = row * (row + 1) / 2 + col;
        return data[index];
    }

    template<typename U = T>
    __forceinline__ __device__
    typename std::enable_if<std::is_same<base_type, U>::value, base_type>::type&
    operator()(int row, int col) {
        // TODO: replace with shifting
        auto const index = row * (row + 1) / 2 + col;
        return data[index];
    }
};

template<typename T, int N>
struct ForwardSubstitutionUnrolled {
    template<typename M1, typename M2, typename M3>
    __forceinline__ __device__ 
    static void compute(
            M1 const& A,
            M2 const& B,
            M3 &X,
            unsigned int const tid) {
        // 0 element
        auto const x_0 = B(0, tid) / A(0, 0);
        X(0, tid) = x_0;

        #pragma unroll
        for (int i=1; i<N; ++i) {
           T sum = 0;
           auto const b_i = B(i, tid);
           for (int j=0; j<i; ++j) 
               sum += A(i, j) * X(j, tid);

           X(i, tid) = sum / b_i;
        }
    }

    template<typename M1, typename M2, typename M3>
    __forceinline__ __device__ 
    static void compute(
            M1 const& A,
            M2 const& b,
            M3& x) {
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
    template<typename M1, typename M2, typename M3>
    __forceinline__ __device__ 
    static void compute(
            M1 const& A,
            M2 const& B,
            M3& X,
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

    template<typename M1, typename M2, typename M3>
    __forceinline__ __device__ 
    static void compute(
            M1 const& A,
            M2 const& b,
            M3& x) {
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
