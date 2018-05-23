#ifndef FWCore_Common_interface_CudaDecls_h
#define FWCore_Common_interface_CudaDecls_h

#ifdef __CUDACC__
    #define DEVICE __device__
    #define HOST __host__
#else
    #define DEVICE
    #define HOST
#endif 

#endif // FWCore_Common_interface_CudaDecls_h
