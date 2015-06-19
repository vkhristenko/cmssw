//	cmssw includes
#include "DQM/HcalTasks/interface/HcalRawTask.h"

//	system includes
#include <iostream>
#include <string>

HcalRawTask::HcalRawTask(edm::ParameterSet const&ps):
	hcaldqm::HcalDQSource(ps)
{}

/* virtual */ HcalRawTask::~HcalRawTask()
{
}

/* virtual */ void HcalRawTask::beginLuminosityBlock(
		edm::LuminosityBlock const& lb, edm::EventSetup const& es)
{
	HcalDQSource::beginLuminosityBlock(lb, es);
}

/* virtual */ void HcalRawTask::endLuminosityBlock(
		edm::LuminosityBlock const& lb, edm::EventSetup const& es)
{
	HcalDQSource::endLuminosityBlock(lb, es);
}

/* virtual */ void HcalRawTask::doWork(edm::Event const& e,
		edm::EventSetup const& es)
{
	//	We need this module only for the Normal Calibration Mode
	if (_mi.currentCalibType>0)
		return;

	edm::Handle<FEDRawDataCollection> craw;
	INITCOLL(_labels["RAW"], craw);
	this->process(*craw);

	//	Per Event Fills
	_mes["NumFEDsUnpackedvsLS"].Fill(_mi.currentLS, _numFEDsUnpackedPerEvent);
}

/* virtyal */ void HcalRawTask::reset(int const periodflag)
{
	HcalDQSource::reset(periodflag);
	if (periodflag==0)
	{
		//	per event
		_numFEDsUnpackedPerEvent=0;
	}
	else if (periodflag==1)
	{
		//	per LS
	}
}

void HcalRawTask::specialize(FEDRawData const& raw, int ifed)
{
	if (this->isuTCA(ifed))
	{
		hcal::AMC13Header const* amc13h = (hcal::AMC13Header const*)(raw.data());
		if (!amc13h)
			return;
		amc13(amc13h, raw.size(), ifed);
		_numFEDsUnpackedPerEvent++;
		_mes["uTCA_FEDsUnpacked"].Fill(ifed);
	}
	else
	{
		HcalDCCHeader const* dcch = (HcalDCCHeader const*)(raw.data());
		if (!dcch)
			return;
		dcc(dcch, raw.size(), ifed);
		_numFEDsUnpackedPerEvent++;
		_mes["VME_FEDsUnpacked"].Fill(ifed);
	}
}

//	Some private functions
bool HcalRawTask::isuTCA(int const ifed) const
{
	if (ifed<FEDNumbering::MINHCALuTCAFEDID || ifed>FEDNumbering::MAXHCALuTCAFEDID)
		return false;
	else 
		return true;

	return false;
}

//	For AMC13/uTCA 
void HcalRawTask::amc13(hcal::AMC13Header const* amc13h, 
		unsigned int const size, int const ifed)
{
	//	Get The Info you need 
//	int				sourceId			= amc13h->sourceId();
	int				bx					= amc13h->bunchId();
	unsigned int	orn					= amc13h->orbitNumber();
	int				l1a					= amc13h->l1aNumber();
	int				namc				= amc13h->NAMC();
//	int				amc13version		= amc13h->AMC13FormatVersion();
	
	//	Iterate over all AMCs
	for (int iamc=0; iamc<namc; iamc++)
	{
		//	Get the info for that AMC13
		int slot		= amc13h->AMCSlot(iamc);
		int crate		= amc13h->AMCId(iamc)&0xFF;
//		int amcsize		= amc13h->AMCSize(iamc)/1000;

		_mes["uTCA_CratesVSslots"].Fill(slot, crate);
		HcalUHTRData uhtr(amc13h->AMCPayload(iamc), amc13h->AMCSize(iamc));
		for (HcalUHTRData::const_iterator iuhtr=uhtr.begin(); iuhtr!=uhtr.end();
				++iuhtr)
		{
			if (!iuhtr.isHeader())
				continue;

			//	Flavor determines what kind of data this uhtr contains
			//	tp, regular digi, upgrade qie digis, etc..
			if (iuhtr.flavor()==hcaldqm::constants::UTCA_DATAFLAVOR)
			{
				//	get the Info you need
				int fiber = (iuhtr.channelid()>>2)&0x1F;
				int fibchannel = iuhtr.channelid()&0x3;
				uint32_t	l1a_uhtr	= uhtr.l1ANumber();
				uint32_t	bx_uhtr		= uhtr.bunchNumber();
				uint32_t	orn_uhtr	= uhtr.orbitNumber();

				//	Fill 
				_mes["uTCA_C" + 
					boost::lexical_cast<std::string>(crate) + "S" + 
					boost::lexical_cast<std::string>(slot) + 
					"_Channels"].Fill(fiber, fibchannel);
//				_mes["uTCA_DataSize"].Fill(amcsize);
				_mes["uTCA_C" + 
					boost::lexical_cast<std::string>(crate) + "S" + 
					boost::lexical_cast<std::string>(slot) + 
					"_EvNComp"].Fill(l1a_uhtr-l1a);
				_mes["uTCA_C" + 
					boost::lexical_cast<std::string>(crate) + "S" + 
					boost::lexical_cast<std::string>(slot) + 
					"_ORNComp"].Fill(orn_uhtr-orn);
				_mes["uTCA_C" + 
					boost::lexical_cast<std::string>(crate) + "S" + 
					boost::lexical_cast<std::string>(slot) + 
					"_BcNComp"].Fill(bx_uhtr-bx);
			}
		}
	}
}

//	For DCC/VME
void HcalRawTask::dcc(HcalDCCHeader const* dcch, unsigned int const size, int const ifed)
{
	//	Get the Info you need
//	unsigned int		bytes			= dcch->getTotalLengthBytes();
	int					sourceId		= dcch->getSourceId();
//	int					bx				= dcch->getBunchId();
//	int					orn				= dcch->getOrbitNumber();
//	unsigned long		evn				= dcch->getDCCEventNumber();
	int					dccid			= sourceId-hcaldqm::constants::VME_DCC_OFFSET;
	HcalHTRData			htr;

	//	Iterate over all spigots
	for (int spigot=0; spigot<HcalDCCHeader::SPIGOT_COUNT; spigot++)
	{
		_mes["VME_DCCvsSpigots"].Fill(spigot, dccid);
//		dcch->getSpigotData(spigot, htr, size);
	}
}

DEFINE_FWK_MODULE(HcalRawTask);



