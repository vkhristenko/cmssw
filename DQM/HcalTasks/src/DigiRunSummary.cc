#include "DQM/HcalTasks/interface/DigiRunSummary.h"

namespace hcaldqm
{
	DigiRunSummary::DigiRunSummary(std::string const& name, 
		std::string const& taskname, edm::ParameterSet const& ps) :
		DQClient(name, taskname, ps)
	{}

	/* virtual */ void DigiRunSummary::beginRun(edm::Run const& r,
		edm::EventSetup const& es)
	{
		DQClient::beginRun(r,es);
	}

	/* virtual */ void DigiRunSummary::endLuminosityBlock(DQMStore::IBooker& ib,
		DQMStore::IGetter& ig, edm::LuminosityBlock const& lb,
		edm::EventSetup const& es)
	{
		DQClient::endLuminosityBlock(ib, ig, lb, es);
	}

	/* virtual */ std::vector<flag::Flag> DigiRunSummary::endJob(
		DQMStore::IBooker& ib, DQMStore::IGetter& ig)
	{
		ib.setCurrentFolder("Hcal/Folder1");
		MonitorElement *me = ib.book1D("MonitorElement", "MonitorElement",
			10, 0, 10);
		for (int i=0; i<10; i++)
			me->Fill(i);

		//	FILTERS, some useful vectors, hash maps
		std::vector<uint32_t> vhashFEDHF; std::vector<uint32_t> vhashVME;
		std::vector<uint32_t> vhashuTCA;
		vhashFEDHF.push_back(HcalElectronicsId(22, SLOT_uTCA_MIN,
			FIBER_uTCA_MIN1, FIBERCH_MIN, false).rawId());
		vhashFEDHF.push_back(HcalElectronicsId(29, SLOT_uTCA_MIN,
			FIBER_uTCA_MIN1, FIBERCH_MIN, false).rawId());
		vhashFEDHF.push_back(HcalElectronicsId(32, SLOT_uTCA_MIN,
			FIBER_uTCA_MIN1, FIBERCH_MIN, false).rawId());
		vhashVME.push_back(HcalElectronicsId(constants::FIBERCH_MIN,
			constants::FIBER_VME_MIN, SPIGOT_MIN, CRATE_VME_MIN).rawId());
		vhashuTCA.push_back(HcalElectronicsId(CRATE_uTCA_MIN, SLOT_uTCA_MIN,
			FIBER_uTCA_MIN1, FIBERCH_MIN, false).rawId());
		filter::HashFilter filter_FEDHF;
		filter_FEDHF.initialize(filter::fPreserver, hashfunctions::fFED,
			vhashFEDHF);	// preserve only HF FEDs
		filter::HashFilter filter_VME;
		filter::HashFilter filter_uTCA;
		filter_uTCA.initialize(filter::fFilter, hashfunctions::fElectronics,
			vhashuTCA);
		filter_VME.initialize(filter::fFilter, hashfunctions::fElectronics,
			vhashVME);
		electronicsmap::ElectronicsMap ehashmap;
		ehashmap.initialize(_emap, electronicsmap::fD2EHashMap);
		std::vector<flag::Flag> vflags; vflags.resize(nDigiFlag);
		vflags[fDead]=flag::Flag("Dead");
		vflags[fUniSlotHF]=flag::Flag("UniSlotHF");
		vflags[fDigiSize]=flag::Flag("DigiSize");

		// INITIALIZE 
		Container2D cOccupancy_depth;
		Container1D cDigiSize_FED;
		Container2D cOccupancyCut_depth;
		ContainerSingle2D cSummary;
		ContainerXXX<double> xDead, xDigiSize, xUniHF, xUni; 
		xDead.initialize(hashfunctions::fFED);
		xDigiSize.initialize(hashfunctions::fFED);
		xUni.initialize(hashfunctions::fFED);
		xUniHF.initialize(hashfunctions::fFEDSlot);
		cOccupancy_depth.initialize(_taskname, "Occupancy",
			hashfunctions::fdepth,
			new quantity::DetectorQuantity(quantity::fieta),
			new quantity::DetectorQuantity(quantity::fiphi),
			new quantity::ValueQuantity(quantity::fN));
		cOccupancyCut_depth.initialize(_taskname, "OccupancyCut",
			hashfunctions::fdepth,
			new quantity::DetectorQuantity(quantity::fieta),
			new quantity::DetectorQuantity(quantity::fiphi),
			new quantity::ValueQuantity(quantity::fN));
		cDigiSize_FED.initialize(_taskname, "DigiSize",
			hashfunctions::fFED,
			new quantity::ValueQuantity(quantity::fDigiSize),
			new quantity::ValueQuantity(quantity::fN));

		_cDead_depth.initialize(_name, "Dead",
			hashfunctions::fdepth,
			new quantity::DetectorQuantity(quantity::fieta),
			new quantity::DetectorQuantity(quantity::fiphi),
			new quantity::ValueQuantity(quantity::fN));
		_cDead_FEDVME.initialize(_name, "Dead",
			hashfunctions::fFED,
			new quantity::ElectronicsQuantity(quantity::fSpigot),
			new quantity::ElectronicsQuantity(quantity::fFiberVMEFiberCh),
			new quantity::ValueQuantity(quantity::fN));
		_cDead_FEDuTCA.initialize(_name, "Dead",
			hashfunctions::fFED,
			new quantity::ElectronicsQuantity(quantity::fSlotuTCA),
			new quantity::ElectronicsQuantity(quantity::fFiberuTCAFiberCh),
			new quantity::ValueQuantity(quantity::fN));
		cSummary.initialize(_name, "Summary",
			new quantity::FEDQuantity(_vFEDs),
			new quantity::FlagQuantity(vflags),
			new quantity::ValueQuantity(quantity::fState));

		//	BOOK
		xDead.book(_emap); xDigiSize.book(_emap); xUniHF.book(_emap);
		xUni.book(_emap, filter_FEDHF);

		//	LOAD
		cOccupancy_depth.load(ig, _emap, _subsystem);
		cOccupancyCut_depth.load(ig, _emap, _subsystem);
		cDigiSize_FED.load(ig, _emap, _subsystem);

		//	BOOK
		_cDead_depth.book(ib, _emap, _subsystem);
		_cDead_FEDVME.book(ib, _emap, filter_uTCA, _subsystem);
		_cDead_FEDuTCA.book(ib, _emap, filter_VME, _subsystem);
		cSummary.book(ib, _subsystem);

		//	iterate over all channels
		std::vector<HcalGenericDetId> gids = _emap->allPrecisionId();
		for (std::vector<HcalGenericDetId>::const_iterator it=gids.begin();
			it!=gids.end(); ++it)
		{
			if (!it->isHcalDetId())
				continue;

			HcalDetId did = HcalDetId(it->rawId());
			HcalElectronicsId eid = HcalElectronicsId(ehashmap.lookup(did));

			if (cOccupancy_depth.getBinContent(did)<1)
			{
				xDead.get(eid)++;
				_cDead_depth.fill(did);
				eid.isVMEid()?_cDead_FEDVME.fill(eid):_cDead_FEDuTCA.fill(eid);
			}
			else
				xDead.get(eid)+=0;
			if (did.subdet()==HcalForward)
				xUniHF.get(eid)+=cOccupancyCut_depth.getBinContent(did);
			cDigiSize_FED.getMean(eid)!=
				constants::DIGISIZE[did.subdet()-1]?
				xDigiSize.get(eid)++:xDigiSize.get(eid)+=0;
			cDigiSize_FED.getRMS(eid)!=0?
				xDigiSize.get(eid)++:xDigiSize.get(eid)+=0;
		}

		//	iterate over all slots in HF
		for (doubleCompactMap::const_iterator it=xUniHF.begin();
			it!=xUniHF.end(); ++it)
		{
			uint32_t hash1 = it->first;
			HcalElectronicsId eid1(hash1);
			double x1 = it->second;

			for (doubleCompactMap::const_iterator jt=xUniHF.begin();
				jt!=xUniHF.end(); ++jt)
			{
				if (jt==it)
					continue;

				double x2 = jt->second;
				if (x2==0)
					continue;
				if (x1/x2<0.2)
					xUni.get(eid1)++;
			}
		}

		//	iterate over all FEDs
		std::vector<flag::Flag> sumflags;
		for (std::vector<uint32_t>::const_iterator it=_vhashFEDs.begin();
			it!=_vhashFEDs.end(); ++it)
		{
			flag::Flag fSum("DIGI");
			HcalElectronicsId eid(*it);

			std::vector<uint32_t>::const_iterator cit=std::find(
				_vcdaqEids.begin(), _vcdaqEids.end(), *it);
			if (cit==_vcdaqEids.end())
			{
				//	not registered @cDAQ
				fSum._state = flag::fNCDAQ;
				sumflags.push_back(flag::Flag("DIGI", flag::fNCDAQ));
				continue;
			}

			//	registered @cDAQ
			if (utilities::isFEDHBHE(eid) || utilities::isFEDHF(eid) ||
				utilities::isFEDHO(eid))
			{
				if (xDead.get(eid)>0)
					vflags[fDead]._state = flag::fBAD;
				else
					vflags[fDead]._state = flag::fGOOD;
				if (utilities::isFEDHF(eid))
				{
					if (xUni.get(eid)>0)
						vflags[fUniSlotHF]._state = flag::fBAD;
					else
						vflags[fUniSlotHF]._state = flag::fGOOD;
				}
				if (xDigiSize.get(eid)>0)
					vflags[fDigiSize]._state = flag::fBAD;
				else
					vflags[fDigiSize]._state = flag::fGOOD;
			}

			//	compute the Summary Flag and push it
			//	also reset all the flags...
			int iflag=0;
			for (std::vector<flag::Flag>::iterator ft=vflags.begin();
				ft!=vflags.end(); ++ft)
			{
				cSummary.setBinContent(eid,iflag,ft->_state);
				fSum+=(*ft);

				//	reset
				ft->reset();
			}
			sumflags.push_back(fSum);
		}

		return sumflags;
	}
}
