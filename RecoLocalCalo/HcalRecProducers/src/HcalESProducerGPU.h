#ifndef RecoLocalCalo_HcalRecProducers_src_HcalESProducerGPU_h
#define RecoLocalCalo_HcalRecProducers_src_HcalESProducerGPU_h

#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/typelookup.h"
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESProductHost.h"
#include "FWCore/Utilities/interface/ReusableObjectHolder.h"

#include <iostream>
#include <array>
#include <tuple>

template<typename Record, typename Target, typename Source>
class HcalESProducerGPU : public edm::ESProducer {
public:
    explicit HcalESProducerGPU(edm::ParameterSet const& ps) 
        : label_{ps.getParameter<std::string>("label")}
    {
        std::string name = ps.getParameter<std::string>("ComponentName");
        setWhatProduced(this, name);
    }
   
    std::unique_ptr<Target> produce(Record const& record) {
        // retrieve conditions in old format 
        edm::ESTransientHandle<Source> product;
        record.get(label_, product);

        return std::make_unique<Target>(*product);
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& confDesc) {
        edm::ParameterSetDescription desc;

        std::string label = Target::name() + "ESProducer";
        desc.add<std::string>("ComponentName", "");
        desc.add<std::string>("label", "")->setComment("Product Label");
        confDesc.add(label, desc);
    }

private:
    std::string label_;
};

template
<
    typename CombinedRecord,
    typename Target,
    typename... Dependencies
>
class HcalESProducerGPUWithDependencies;

template
<
    template<typename...> typename CombinedRecord, typename... DepsRecords,
    typename Target, 
    typename... Dependencies
>
class HcalESProducerGPUWithDependencies
<
    CombinedRecord<DepsRecords...>,
    Target,
    Dependencies...
> 
    : public edm::ESProducer {
public:
    static constexpr std::size_t nsources = sizeof...(Dependencies);
    /*using HostType = edm::ESProductHost<Target,
                                        DepsRecords...>;
                                        */

    explicit HcalESProducerGPUWithDependencies(edm::ParameterSet const& ps) {
        for (std::size_t i=0; i<labels_.size(); i++)
            labels_[i] = ps.getParameter<std::string>("label" + std::to_string(i));

        std::string name = ps.getParameter<std::string>("ComponentName");
        setWhatProduced(this, name);
    }

    std::unique_ptr<Target> produce(CombinedRecord<DepsRecords...> const& record) {
        auto handles = std::tuple<edm::ESTransientHandle<Dependencies>...>{};
        WalkAndCall<nsources-1, 
                    edm::ESTransientHandle<Dependencies>...>::iterate(
            record, handles, labels_);

        return std::apply(
            [](auto const&... handles) {
                return std::make_unique<Target>((*handles)...);
            },
            handles
        );
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& confDesc) {
        edm::ParameterSetDescription desc;

        std::string label = Target::name() + "ESProducerWithDependencies";
        desc.add<std::string>("ComponentName", "");
        for (std::size_t i=0; i<nsources; i++)
            desc.add<std::string>("label" + std::to_string(i), "")
                ->setComment("Product Label");
        confDesc.add(label, desc);
    }

private:
    template<std::size_t N, typename... Types>
    struct WalkAndCall {
        static void iterate(
                CombinedRecord<DepsRecords...> const& containingRecord, 
                std::tuple<Types...>& ts,
                std::array<std::string, nsources> const& labels) {
            using Record = 
                typename std::tuple_element<N, std::tuple<DepsRecords...>>::type;
            auto const& record = containingRecord.template getRecord<Record>();
            record.get(labels[N], std::get<N>(ts));
            WalkAndCall<N-1, Types...>::iterate(containingRecord, ts, labels);
        }
    };

    template<typename... Types>
    struct WalkAndCall<0, Types...> {
        static void iterate(
                CombinedRecord<DepsRecords...> const& containingRecord, 
                std::tuple<Types...>& ts,
                std::array<std::string, nsources> const& labels) {
            using Record = 
                typename std::tuple_element<0, std::tuple<DepsRecords...>>::type;
            auto const& record = containingRecord.template getRecord<Record>();
            record.get(labels[0], std::get<0>(ts));
        }
    };

private:
    std::array<std::string, nsources> labels_;
    //edm::ReusableObjectHolder<HostType> holder_;
};

#endif
