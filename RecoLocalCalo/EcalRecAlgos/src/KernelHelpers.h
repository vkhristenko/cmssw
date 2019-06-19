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
    base_type const& operator()(int const i) const { return data[i]; }
    
    template<typename U = T>
    __forceinline__ __device__
    typename std::enable_if<std::is_same<base_type, U>::value, base_type>::type&
    operator()(int const i) { return data[i]; }
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
    base_type const& operator()(int const row, int const col) const { return data[col*Stride + row]; }

    template<typename U = T>
    __forceinline__ __device__
    typename std::enable_if<std::is_same<base_type, U>::value, base_type>::type&
    operator()(int const row, int const col) { return data[col*Stride + row]; }
};

template<typename T, int Stride>
struct MapM<T, Stride, Eigen::RowMajor> {
    using type = T;
    using base_type = typename std::remove_const<type>::type;

    type* data;
    __forceinline__ __device__
    MapM(type* data) : data{data} {}

    __forceinline__ __device__
    base_type const& operator()(int const row, int const col) const { return data[row*Stride + col]; }

    template<typename U = T>
    __forceinline__ __device__
    typename std::enable_if<std::is_same<base_type, U>::value, base_type>::type&
    operator()(int const row, int const col) { return data[row*Stride + col]; }
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
    T const& operator()(int const row, int const col) const {
        // TODO: replace with shifting
        auto const tmp = (Stride - col) * (Stride - col + 1) / 2;
        auto const index = total - tmp + row - col;
        return data[index];
    }
    
    template<typename U = T>
    __forceinline__ __device__
    typename std::enable_if<std::is_same<base_type, U>::value, base_type>::type&
    operator()(int const row, int const col) {
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
    T const& operator()(int const row, int const col) const {
        // TODO: replace with shifting 
        auto const index = row * (row + 1) / 2 + col;
        return data[index];
    }

    template<typename U = T>
    __forceinline__ __device__
    typename std::enable_if<std::is_same<base_type, U>::value, base_type>::type&
    operator()(int const row, int const col) {
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

           X(i, tid) = (b_i - sum) / A(i, i);
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

           x(i) = (b_i - sum) / A(i, i);
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

template<typename T>
struct BackwardSubstitution {
    template<typename M1, typename M2, typename M3>
    __forceinline__ __device__ 
    static void compute(
            M1 const& A,
            M2 const& B,
            M3& X,
            unsigned int const tid,
            int view) {
        // first element (last one in the vector)
        X(view-1, tid) = B(view-1, tid)/A(view-1, view-1);
        
        // the rest
        for (int i=view-2; i>=0; --i) {
            T sum = 0;
            auto const b_i = B(i, tid);
            for (int j=i+1; j<view; ++j)
                sum += A(j, i) * X(j, tid);

            X(i, tid) = (b_i - sum) / A(i, i);
        }
    }

    template<typename M1, typename M2, typename M3>
    __forceinline__ __device__ 
    static void compute(
            M1 const& A,
            M2 const& b,
            M3& x,
            int view) {
        // first element (last one in the vector)
        x(view-1) = b(view-1)/A(view-1, view-1);
        
        // the rest
        for (int i=view-2; i>=0; --i) {
            T sum = 0;
            auto const b_i = b(i);
            for (int j=i+1; j<view; ++j)
                sum += A(j, i) * x(j);

            x(i) = (b_i - sum) / A(i, i);
        }
    }
};

template<typename T, int N>
struct FusedCholeskySolver;

template<typename T>
struct FusedCholeskySolver<T, 1> {
    template<typename M1, typename M2, typename M3>
    __forceinline__ __device__ 
    static void compute(
            M1 const& M,
            M2 const& b,
            M3 &x,
            MapV<char> const& mapping) {
        auto const real_0 = mapping(0);
        auto const x_0 = b(real_0) / M(real_0, real_0);
        x(0) = x_0;
    }
};

template<typename T>
struct FusedCholeskySolver<T, 2> {
    template<typename M1, typename M2, typename M3>
    __forceinline__ __device__ 
    static void compute(
            M1 const& M,
            M2 const& b,
            M3 &x,
            MapV<char> const& mapping) {
        // element 0
        auto const real_0 = mapping(0);
        auto const real_1 = mapping(1);
        auto const l_0_0 = std::sqrt(M(real_0, real_0));
        auto const interm_0 = b(real_0) / l_0_0;

        // element 1
        auto const l_1_0 = M(real_1, real_0) / l_0_0;
        auto const l_1_1 = std::sqrt(M(real_1, real_1) - l_1_0*l_1_0);
        auto const interm_1 = (b(real_1) - interm_0 * l_1_0) / l_1_1;
        auto const x_1 = interm_1 / l_1_1;
        x(1) = x_1;
        auto const x_0 = (interm_0 - l_1_0 * x_1) / l_0_0;
        x(0) = x_0;
    }
};

template<typename T>
struct FusedCholeskySolver<T, 3> {
    template<typename M1, typename M2, typename M3>
    __forceinline__ __device__ 
    static void compute(
            M1 const& M,
            M2 const& b,
            M3 &x,
            MapV<char> const& mapping) {
        // element 0
        auto const real_0 = mapping(0);
        auto const l_0_0 = std::sqrt(M(real_0, real_0));
        auto const interm_0 = b(real_0) / l_0_0;

        // row 1
        auto const real_1 = mapping(1);
        auto const l_1_0 = M(real_1, real_0) / l_0_0;
        auto const l_1_1 = std::sqrt(M(real_1, real_1) - l_1_0*l_1_0);
        auto const interm_1 = (b(real_1) - interm_0 * l_1_0) / l_1_1;

        // row 2
        auto const real_2 = mapping(2);
        auto const l_2_0 = M(real_2, real_0) / l_0_0;
        auto const l_2_1 = (M(real_2, real_1) - l_2_0 * l_1_0) / l_1_1;
        auto const l_2_2 = std::sqrt(M(real_2, real_2) - l_2_0 * l_2_0 - l_2_1*l_2_1);
        auto const interm_2 = (b(real_2) - interm_0 * l_2_0 - interm_1 * l_2_1) / l_2_2;

        auto const x_2 = interm_2 / l_2_2;
        x(2) = x_2;
        auto const x_1 = (interm_1 - l_2_1 * x_2) / l_1_1;
        x(1) = x_1;
        auto const x_0 = (interm_0 - l_1_0 * x_1 - l_2_0 * x_2) / l_0_0;
        x(0) = x_0;
    }
};

template<typename T>
struct FusedCholeskyForwardSubst {
    template<typename M1, typename M2, typename M3, typename M4>
    __forceinline__ __device__ 
    static void compute(
            M1 const& M, 
            M2 const& b,
            M3 &L,
            M4 &intermediate,
            MapV<char> const& mapping,
            int view) {
        // compute element 0,0 for L
        auto const real_0 = mapping(0);
        auto const sqrtm_0_0 = std::sqrt(M(real_0, real_0));
        L(0, 0) = sqrtm_0_0;

        // compute solution for forward subst for element 0
        auto const interm_0 = b(real_0) / sqrtm_0_0;
        intermediate(0) = interm_0;

        for (int i=1; i<view; ++i) {
            // load the value to sub from
            auto const real_i = mapping(i);
            T total = b(real_i);

            // first compute elements to the left of the diagoanl
            T sumsq{static_cast<T>(0)};
            for (int j=0; j<i; ++j) {
                T sumsq2{static_cast<T>(0)};
                auto const real_j = mapping(j);
                auto const m_i_j = M(real_i, real_j);
                for (int k=0; k<j; ++k)
                    sumsq2 += L(i, k) * L(j, k);

                // comput the i,j : i>j, elements to the left of the diagonal
                auto const value_i_j = (m_i_j - sumsq2) / L(j, j);
                L(i, j) = value_i_j;

                // needed to compute diagonal element
                sumsq += value_i_j * value_i_j;

                total -= value_i_j * intermediate(j);
            }

            // second, compute the diagonal element
            auto const l_i_i = std::sqrt(M(real_i, real_i) - sumsq);
            L(i, i) = l_i_i;

            intermediate(i) = total / l_i_i;
        }
    }
};

template<typename T, int N>
struct FusedCholeskyForwardSubstUnrolled {
    template<typename M1, typename M2, typename M3, typename M4>
    __forceinline__ __device__ 
    static void compute(
            M1 const& M, 
            M2 const& b,
            M3 &L,
            M4 &intermediate,
            MapV<char> const& mapping) {
        // compute element 0,0 for L
        auto const real_0 = mapping(0);
        auto const sqrtm_0_0 = std::sqrt(M(real_0, real_0));
        L(0, 0) = sqrtm_0_0;

        // compute solution for forward subst for element 0
        auto const interm_0 = b(real_0) / sqrtm_0_0;
        intermediate(0) = interm_0;

        #pragma unroll
        for (int i=1; i<N; ++i) {
            // load the value to sub from
            auto const real_i = mapping(i);
            T total = b(real_i);

            // first compute elements to the left of the diagoanl
            T sumsq{static_cast<T>(0)};
            for (int j=0; j<i; ++j) {
                T sumsq2{static_cast<T>(0)};
                auto const real_j = mapping(j);
                auto const m_i_j = M(real_i, real_j);
                for (int k=0; k<j; ++k)
                    sumsq2 += L(i, k) * L(j, k);

                // comput the i,j : i>j, elements to the left of the diagonal
                auto const value_i_j = (m_i_j - sumsq2) / L(j, j);
                L(i, j) = value_i_j;

                // needed to compute diagonal element
                sumsq += value_i_j * value_i_j;

                total -= value_i_j * intermediate(j);
            }

            // second, compute the diagonal element
            auto const l_i_i = std::sqrt(M(real_i, real_i) - sumsq);
            L(i, i) = l_i_i;

            intermediate(i) = total / l_i_i;
        }
    }
};

}}

#endif // RecoLocalCalo_EcalRecAlgos_src_KernelHelpers_h
