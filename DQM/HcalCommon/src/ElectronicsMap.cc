#include "DQM/HcalCommon/interface/ElectronicsMap.h"
#include <iomanip>

namespace hcaldqm
{
	namespace electronicsmap
	{
		void ElectronicsMap::initialize(HcalElectronicsMap const* emap,
			ElectronicsMapType etype/*=fHcalElectronicsMap*/)
		{
			_etype=etype;
			_emap = emap;
			//	if we actually use a HashMap then 
			if (_etype!=fHcalElectronicsMap)
			{
				if (_etype==fD2EHashMap)
				{
					//	HcalDetId -> HcalElectronicsId hash map
					std::vector<HcalGenericDetId> dids=emap->allPrecisionId();
					for (std::vector<HcalGenericDetId>::const_iterator it=
						dids.begin(); it!=dids.end(); ++it)
					{
						if (!it->isHcalDetId())
							continue;

						HcalElectronicsId eid = _emap->lookup(
							HcalDetId(it->rawId()));
						uint32_t hash = it->rawId();
						EMapType::iterator eit = _ids.find(hash);
						if (eit!=_ids.end())
							continue;

						_ids.insert(std::make_pair(hash, eid.rawId()));
					}
				}
				else if (_etype==fT2EHashMap)
				{
					//	HcalTrigTowerDetId -> HcalElectronicsId
					std::vector<HcalTrigTowerDetId> tids = 
						emap->allTriggerId();
					for (std::vector<HcalTrigTowerDetId>::const_iterator it=
						tids.begin(); it!=tids.end(); ++it)
					{
						HcalElectronicsId eid = _emap->lookupTrigger(*it);
						uint32_t hash = it->rawId();
						EMapType::iterator eit = _ids.find(hash);
						if (eit!=_ids.end())
							continue;
						
						_ids.insert(std::make_pair(hash, eid.rawId()));	
					}
				}
				else if (_etype==fE2DHashMap)
				{
					//	HcalElectronicId -> HcalDetId hash map
					std::vector<HcalElectronicsId> eids = 
						emap->allElectronicsIdPrecision();
					for (std::vector<HcalElectronicsId>::const_iterator
						it=eids.begin(); it!=eids.end(); ++it)
					{
						HcalDetId did = HcalDetId(_emap->lookup(*it));
						EMapType::iterator eit = _ids.find(it->rawId());
						if (eit!=_ids.end())
							continue;
						
						//	eid.rawId() -> did.rawId()
						_ids.insert(std::make_pair(it->rawId(), 
							did.rawId()));
					}
				}
				else if (_etype==fE2THashMap)
				{
					//	HcalElectronicId -> HcalDetId hash map
					std::vector<HcalElectronicsId> eids = 
						emap->allElectronicsIdTrigger();
					for (std::vector<HcalElectronicsId>::const_iterator
						it=eids.begin(); it!=eids.end(); ++it)
					{
						HcalTrigTowerDetId tid = 
							HcalTrigTowerDetId(_emap->lookupTrigger(*it));
						EMapType::iterator eit = _ids.find(it->rawId());
						if (eit!=_ids.end())
							continue;
						
						//	eid.rawId() -> tid.rawId()
						_ids.insert(std::make_pair(it->rawId(), 
							tid.rawId()));
					}
				}
			}
		}

		void ElectronicsMap::initialize(HcalElectronicsMap const* emap,
			ElectronicsMapType etype, filter::HashFilter const& filter)
		{
			_etype=etype;
			_emap = emap;
			
			//	note this initialization has iteration over electronics not 
			//	detector. 
			//	Filtering is done on Electronics id - possible to have 
			//	several electronics ids to 1 detid - not vice versa
			if (_etype!=fHcalElectronicsMap)
			{
				if (_etype==fD2EHashMap)
				{
					std::vector<HcalElectronicsId> eids = 
						emap->allElectronicsIdPrecision();
					for (std::vector<HcalElectronicsId>::const_iterator it=
						eids.begin(); it!=eids.end(); ++it)
					{
						HcalGenericDetId did = HcalGenericDetId(
							_emap->lookup(*it));
						if (filter.filter(*it))
							continue;
						if (!did.isHcalDetId())
							continue;

						_ids.insert(std::make_pair(did.rawId(), it->rawId()));
					}
				}
				else if (_etype==fT2EHashMap)
				{
					std::vector<HcalElectronicsId> eids=
						emap->allElectronicsIdTrigger();
					for (std::vector<HcalElectronicsId>::const_iterator it=
						eids.begin(); it!=eids.end(); ++it)
					{
						if (filter.filter(*it))
							continue;
						HcalTrigTowerDetId tid = emap->lookupTrigger(*it);
						_ids.insert(std::make_pair(tid.rawId(), it->rawId()));
					}
				}
				else if (_etype==fE2DHashMap)
				{
					std::vector<HcalElectronicsId> eids = 
						emap->allElectronicsIdPrecision();
					for (std::vector<HcalElectronicsId>::const_iterator it=
						eids.begin(); it!=eids.end(); ++it)
					{
						HcalGenericDetId did = HcalGenericDetId(
							_emap->lookup(*it));
						if (filter.filter(*it))
							continue;
						if (!did.isHcalDetId())
							continue;

						_ids.insert(std::make_pair(it->rawId(), did.rawId()));
					}
				}
				else if (_etype==fE2THashMap)
				{
					std::vector<HcalElectronicsId> eids=
						emap->allElectronicsIdTrigger();
					for (std::vector<HcalElectronicsId>::const_iterator it=
						eids.begin(); it!=eids.end(); ++it)
					{
						if (filter.filter(*it))
							continue;
						HcalTrigTowerDetId tid = emap->lookupTrigger(*it);
						_ids.insert(std::make_pair(it->rawId(), tid.rawId()));
					}
				}
			}
		}

		//	2 funcs below are only for 1->1 mappings
		uint32_t ElectronicsMap::lookup(DetId const &did)
		{
			uint32_t hash = did.rawId();
			return _etype==fHcalElectronicsMap? _emap->lookup(did).rawId(): 
				_ids[hash];
		}

		uint32_t ElectronicsMap::lookup(HcalElectronicsId const &did)
		{
			uint32_t hash = did.rawId();
			return _etype==fHcalElectronicsMap? _emap->lookup(did).rawId():
				_ids[hash];
		}

		void ElectronicsMap::print()
		{
			std::cout << "Electronics HashMap Type=" << _etype  << std::endl;
			BOOST_FOREACH(EMapType::value_type &v, _ids)
			{
				std::cout << std::hex << v.first 
					 << "  "<< v.second << std::dec << std::endl;
			}
		}
	}
}
