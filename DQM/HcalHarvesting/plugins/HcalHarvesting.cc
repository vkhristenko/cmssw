#include "DQM/HcalHarvesting/interface/HcalHarvesting.h"

HcalHarvesting::HcalHarvesting(edm::ParameterSet const& ps) :
	DQHarvester(ps)
{
	//	set labels for summaries
	_frawnames.push_back("EvnMsm");
	_frawnames.push_back("BcnMsm");
	_frawnames.push_back("BadQuality");

	_fdiginames.push_back("UniSlot");
	_fdiginames.push_back("Msn1LS");
	_fdiginames.push_back("CapIdRot");
	_fdiginames.push_back("DigiSize");

	_freconames.push_back("OcpUniSlot");
	_freconames.push_back("TimeUniSlot");
	_freconames.push_back("TCDS");
	_freconames.push_back("Msn1LS");

	_ftpnames.push_back("OcpUniSlotD");
	_ftpnames.push_back("OcpUniSlotE");
	_ftpnames.push_back("EtMsmUniSlot");
	_ftpnames.push_back("FGMsmUniSlot");
	_ftpnames.push_back("MsnUniSlotD");
	_ftpnames.push_back("MsnUniSlotE");
	_ftpnames.push_back("EtCorrRatio");
	_ftpnames.push_back("EtMsmRatio");
	_ftpnames.push_back("FGMsmNumber");
}

/* virtual */ void HcalHarvesting::_dqmEndLuminosityBlock(DQMStore::IBooker& ib,
	DQMStore::IGetter& ig, edm::LuminosityBlock const&, 
	edm::EventSetup const&)
{
	//	Initialize what you need
	ContainerSingle2D rawSummary;
	ContainerSingle2D digiSummary;
	ContainerSingle2D recoSummary;
	ContainerSingle2D tpSummary;
	rawSummary.initialize("RawTask", "Summary",
		new quantity::FEDQuantity(_vFEDs),
		new quantity::FlagQuantity(_frawnames),
		new quantity::QualityQuantity());
	digiSummary.initialize("DigiTask", "Summary",
		new quantity::FEDQuantity(_vFEDs),
		new quantity::FlagQuantity(_fdiginames),
		new quantity::QualityQuantity());
	recoSummary.initialize("RecHitTask", "Summary",
		new quantity::FEDQuantity(_vFEDs),
		new quantity::FlagQuantity(_freconames),
		new quantity::QualityQuantity());
	tpSummary.initialize("TPTask", "Summary",
		new quantity::FEDQuantity(_vFEDs),
		new quantity::FlagQuantity(_freconames),
		new quantity::QualityQuantity());

	ContainerSingle2D rawSummaryCopy;
	ContainerSingle2D digiSummaryCopy;
	ContainerSingle2D recoSummaryCopy;
	ContainerSingle2D tpSummaryCopy;
	char name[20];
	sprintf(name, "LS%d", _currentLS);
	rawSummaryCopy.initialize("RawTask", "SummaryvsLS",
		new quantity::FEDQuantity(_vFEDs),
		new quantity::FlagQuantity(_frawnames),
		new quantity::QualityQuantity());
	digiSummaryCopy.initialize("DigiTask", "SummaryvsLS",
		new quantity::FEDQuantity(_vFEDs),
		new quantity::FlagQuantity(_fdiginames),
		new quantity::QualityQuantity());
	recoSummaryCopy.initialize("RecHitTask", "SummaryvsLS",
		new quantity::FEDQuantity(_vFEDs),
		new quantity::FlagQuantity(_freconames),
		new quantity::QualityQuantity());
	tpSummaryCopy.initialize("TPTask", "SummaryvsLS",
		new quantity::FEDQuantity(_vFEDs),
		new quantity::FlagQuantity(_ftpnames),
		new quantity::QualityQuantity());

	//	book the new plots and load existing ones
	rawSummaryCopy.book(ib, "Hcal", name);
	digiSummaryCopy.book(ib, "Hcal", name);
	recoSummaryCopy.book(ib, "Hcal", name);
	tpSummaryCopy.book(ib, "Hcal", name);
	rawSummary.load(ig);
	digiSummary.load(ig);
	recoSummary.load(ig);
	tpSummary.load(ig);

	//	process: put the quality into the copy and set the reportSummaryMap
	//	contents
	MonitorElement *reportSummaryMap = ig.get(
		"Hcal/EventInfo/reportSummaryMap");
	int ifed = 0;
	for (std::vector<uint32_t>::const_iterator it=_vhashFEDs.begin();
		it!=_vhashFEDs.end(); ++it)
	{
		HcalElectronicsId eid(*it);
		//	RAW
		int counter = 0;
		for (uint32_t f=0; f<_frawnames.size(); f++)
		{
			quantity::Quality q = (quantity::Quality)
				((int)rawSummary.getBinContent(eid, (int)f));
			rawSummaryCopy.setBinContent(eid, (int)f, q);
			if (q>quantity::fGood)
				counter++;
		}
		counter>0?
			reportSummaryMap->setBinContent(ifed+1, 1, quantity::fLow):
			reportSummaryMap->setBinContent(ifed+1, 1, quantity::fGood);

		//	DIGI
		counter=0;
		for (uint32_t f=0; f<_fdiginames.size(); f++)
		{
			quantity::Quality q = (quantity::Quality)
				((int)digiSummary.getBinContent(eid, (int)f));
			digiSummaryCopy.setBinContent(eid, (int)f, q);
			if (q>quantity::fGood)
				counter++;
		}
		counter>0?
			reportSummaryMap->setBinContent(ifed+1, 2, quantity::fLow):
			reportSummaryMap->setBinContent(ifed+1, 2, quantity::fGood);

		//	RECO
		counter=0;
		for (uint32_t f=0; f<_freconames.size(); f++)
		{
			quantity::Quality q = (quantity::Quality)
				((int)recoSummary.getBinContent(eid, (int)f));
			recoSummaryCopy.setBinContent(eid, (int)f, q);
			if (q>quantity::fGood)
				counter++;
		}
		counter>0?
			reportSummaryMap->setBinContent(ifed+1, 3, quantity::fLow):
			reportSummaryMap->setBinContent(ifed+1, 3, quantity::fGood);

		//	TP
		counter=0;
		for (uint32_t f=0; f<_ftpnames.size(); f++)
		{
			quantity::Quality q = (quantity::Quality)
				((int)tpSummary.getBinContent(eid, (int)f));
			tpSummaryCopy.setBinContent(eid, (int)f, q);
			if (q>quantity::fGood)
				counter++;
		}
		counter>0?
			reportSummaryMap->setBinContent(ifed+1, 4, quantity::fLow):
			reportSummaryMap->setBinContent(ifed+1, 4, quantity::fGood);

		ifed++;
	}
}

/* virtual */ void HcalHarvesting::_dqmEndJob(DQMStore::IBooker&,
	DQMStore::IGetter&)
{}

DEFINE_FWK_MODULE(HcalHarvesting);
