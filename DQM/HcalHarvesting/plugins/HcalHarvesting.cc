#include "DQM/HcalHarvesting/interface/HcalHarvesting.h"

HcalHarvesting::HcalHarvesting(edm::ParameterSet const& ps) :
	DQHarvester(ps)
{
	_frawnames.push_back("EvnMsm");
	_frawnames.push_back("BcnMsm");
	_frawnames.push_back("BadQuality");
}

/* virtual */ void HcalHarvesting::_dqmEndLuminosityBlock(DQMStore::IBooker& ib,
	DQMStore::IGetter& ig, edm::LuminosityBlock const&, 
	edm::EventSetup const&)
{
	std::cout << "Harvesting LS" << _currentLS << std::endl;

	//	Save the Summary per LS
	ContainerSingle2D rawSummary;
	rawSummary.initialize("RawTask", "Summary",
		new quantity::FEDQuantity(_vFEDs),
		new quantity::FlagQuantity(_frawnames),
		new quantity::QualityQuantity());
	rawSummary.load(ig);
	ContainerSingle2D rawSummaryCopy;
	char name[20];
	sprintf(name, "LS%d", _currentLS);
	rawSummaryCopy.initialize("RawTask", "SummaryvsLS",
		new quantity::FEDQuantity(_vFEDs),
		new quantity::FlagQuantity(_frawnames),
		new quantity::QualityQuantity());
	rawSummaryCopy.book(ib, "Hcal", name);
	for (std::vector<uint32_t>::const_iterator it=_vhashFEDs.begin();
		it!=_vhashFEDs.end(); ++it)
	{
		HcalElectronicsId eid(*it);
		for (uint32_t f=0; f<_frawnames.size(); f++)
			rawSummaryCopy.setBinContent(eid, (int)f, 
				rawSummary.getBinContent(eid, (int)f));
	}
}

/* virtual */ void HcalHarvesting::_dqmEndJob(DQMStore::IBooker&,
	DQMStore::IGetter&)
{}

DEFINE_FWK_MODULE(HcalHarvesting);
