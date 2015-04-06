//	cmssw includes
#include "DQM/HcalTasks/interface/HcalDigiTask.h"

//	system includes
#include <iostream>
#include <string>

HcalDigiTask::HcalDigiTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalDigiTask::~HcalDigiTask()
{}

/* virtual */ void HcalDigiTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
	edm::Handle<DcsStatusCollection>			cdcs;
	edm::Handle<edm::TriggerResults>			hltRes;
	edm::Handle<HBHEDigiCollection>				chbhe;
	edm::Handle<HODigiCollection>				cho;
	edm::Handle<HFDigiCollection>				chf;
	edm::Handle<HcalUnpackerReport>				report;
	edm::Handle<FEDRawDataCollection>			craw;

	INITCOLL(_labels["HBHEDigi"], chbhe);
//	INITCOLL(_labels["DCS"], cdcs);
	INITCOLL(_labels["HLTResults"], hltRes);
	INITCOLL(_labels["HODigi"], cho);
	INITCOLL(_labels["HFDigi"], chf);
//	INITCOLL(_labels["UnpackerReport"], report);
	INITCOLL(_labels["RAW"], craw);

	MonitorElement &mHF = _mes["HFDigiShape"];
	MonitorElement &mHBHE =  _mes["HEDigiShape"];
	MonitorElement &mHO = _mes["HODigiShape"];
	MonitorElement &mDigiSizeCheck = _mes["DigiSizeCheck"];

	//	Example
	//	-> Shapes
	//	-> DigiSizeCheck
	for (HFDigiCollection::const_iterator it=chf->begin(); 
			it!=chf->end(); ++it)
	{
		const HFDataFrame digi = (const HFDataFrame)(*it);
		for (int i=0; i<digi.size(); i++)
			mHF.Fill(i+1, digi.sample(i).adc());
		mDigiSizeCheck.Fill(0.5, digi.size());
	}
	for (HBHEDigiCollection::const_iterator it=chbhe->begin();
			it!=chbhe->end(); ++it)
	{
		const HBHEDataFrame digi = (const HBHEDataFrame)(*it);
		for (int i=0; i<digi.size(); i++)
			mHBHE.Fill(i+1, digi.sample(i).adc());
		mDigiSizeCheck.Fill(1.5, digi.size());
	}
	for (HODigiCollection::const_iterator it=cho->begin();
			it!=cho->end(); ++it)
	{
		const HODataFrame digi = (const HODataFrame)(*it);
		for (int i=0; i<digi.size(); i++)
			mHO.Fill(i+1, digi.sample(i).adc());
		mDigiSizeCheck.Fill(2.5, digi.size());
	}
}

DEFINE_FWK_MODULE(HcalDigiTask);



