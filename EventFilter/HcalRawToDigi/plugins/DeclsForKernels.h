#ifndef EventFilter_HcalRawToDigi_interface_DeclsForKernels_h
#define EventFilter_HcalRawToDigi_interface_DeclsForKernels_h

namespace hcal { namespace raw {

constexpr int32_t empty_event_size = 32;

struct ConfigurationParameters {
    uint32_t maxChannels;
};

}}

#endif // EventFilter_HcalRawToDigi_interface_DeclsForKernels_h
