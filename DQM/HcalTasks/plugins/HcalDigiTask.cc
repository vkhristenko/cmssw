//	cmssw includes
#include "DQM/HcalTasks/interface/HcalDigiTask.h"

//	system includes
#include <iostream>
#include <string>

HcalDigiTask::HcalDigiTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{
	_ornMsgTime = ps.getUntrackedParameter<int>("OrnMsgTime", 3559);
	for (int i=0; i<4; i++)
		_numDigis.push_back(0);
}

/* virtual */ HcalDigiTask::~HcalDigiTask()
{
}

/* virtual */ void HcalDigiTask::beginLuminosityBlock(
		edm::LuminosityBlock const& lb, edm::EventSetup const& es)
{
	HcalDQSource::beginLuminosityBlock(lb, es);
}

/* virtual */ void HcalDigiTask::endLuminosityBlock(
		edm::LuminosityBlock const& lb, edm::EventSetup const& es)
{
	HcalDQSource::endLuminosityBlock(lb, es);
}

/* virtual */ void HcalDigiTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
	edm::Handle<HBHEDigiCollection>			chbhe;
	edm::Handle<HODigiCollection>			cho;
	edm::Handle<HFDigiCollection>			chf;
	
	INITCOLL(_labels["HBHEDigi"], chbhe);
	INITCOLL(_labels["HODigi"], cho);
	INITCOLL(_labels["HFDigi"], chf);

	this->process(*chbhe, std::string("HB"));
	this->process(*chbhe, std::string("HE"));
	this->process(*cho, std::string("HO"));
	this->process(*chf, std::string("HF"));
	
	_mes["HBNumberDigis"].Fill(_numDigis[0]);
	_mes["HENumberDigis"].Fill(_numDigis[1]);
	_mes["HONumberDigis"].Fill(_numDigis[2]);
	_mes["HFNumberDigis"].Fill(_numDigis[3]);
}

/* virtual */ void HcalDigiTask::reset(int const periodflag)
{
	HcalDQSource::reset(periodflag);
	if (periodflag==0)
	{
		// each events reset
		for (unsigned int i=0; i<_numDigis.size(); i++)
			_numDigis[0]=0;
	}
	else if (periodflag==1)
	{
		// each LS reset
	}
}

//	specializer
template<typename Hit>
void HcalDigiTask::specialize(Hit const& hit, std::string const& nameRes)
{
	//	offset of BCN for this channel(digi) relative to the nominal 
	//	in the unpacker. range(-7, +7), -1000 indicates invalid data
	int offset = hit.fiberIdleOffset();
	if (offset!=INVALID_FIBER_IDLEOFFSET && _ornMsgTime>-1)
		_mes[nameRes+"bcnOffset"].Fill(offset);

	int iphi = hit.id().iphi();
	int ieta = hit.id().ieta();
	int depth = hit.id().depth();
	if (nameRes!="HO")
		_mes["HBHEHFOccupancyMapD" + 
			boost::lexical_cast<std::string>(depth)].Fill(ieta, iphi, 1);

	//	Get the digi size
	//	To be done better
	int subdet = 0;
	if (nameRes=="HF")
		subdet = SUBDET_HF;
	else if (nameRes=="HB" || nameRes=="HE")
		subdet = SUBDET_HBHE;
	else if (nameRes=="HO")
		subdet = SUBDET_HO;
	else 
		subdet = SUBDET_OTHER;
	_mes["DigiSizeCheck"].Fill(subdet, hit.size());
	
	//	Extract per digi Information
	for (int i=0; i<hit.size(); i++)
	{
		_mes[nameRes+"DigiShape"].Fill(i, hit.sample(i).nominal_fC());
		_mes[nameRes+"ADCCountPerTS"].Fill(hit.sample(i).adc());
		_mes[nameRes+"Presamples"].Fill(hit.presamples());
		_mes[nameRes+"CapId"].Fill(hit.sample(i).capid());
	}
}

DEFINE_FWK_MODULE(HcalDigiTask);



