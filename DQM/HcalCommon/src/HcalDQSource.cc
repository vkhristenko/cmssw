

#include "DQM/HcalCommon/interface/HcalDQSource.h"

namespace hcaldqm
{
	HcalDQSource::HcalDQSource(edm::ParameterSet const& ps) :
		HcalDQMonitor(ps.getUntrackedParameterSet("moduleParameters")), 
		_mes(ps.getUntrackedParameterSet("MEs"), _mi.debug)
	{
		_si.currentCalibType = -1;
		_si.evsTotal = 0;
		_si.evsGood = 0;
		_si.evsPerLS = 0;
	}

	/* virtual */HcalDQSource::~HcalDQSource() 
	{
		this->debug_("Calling Destructor...");
	}

	//	Functions to be reimplemented from DQMEDAnalyzer
	/* virtual */ void HcalDQSource::analyze(edm::Event const &e,
			edm::EventSetup const& es)
	{
		this->reset(0);
		//	Update event counters;
		_si.evsTotal++; _mes["EventsProcessed"].Fill(_si.evsTotal);
		_si.evsGood++;
		_si.evsPerLS++; _mes["EventsProcessedPerLS"].Fill(_si.evsPerLS);

		this->extractCalibType(e);
		if (this->isAllowedCalibType()==false)
			return;
		this->doWork(e, es);
	}

	/* virtual */ void HcalDQSource::bookHistograms(DQMStore::IBooker &ib,
			edm::Run const& r, edm::EventSetup const& es)
	{
		_mes.book(ib);
	}

	/* virtual */ void HcalDQSource::dqmBeginRun(edm::Run const& r,
			edm::EventSetup const& es)
	{
		this->reset(0);
		this->reset(1);
	}

	/* virtual */ void HcalDQSource::beginLuminosityBlock(
			edm::LuminosityBlock const& lb, edm::EventSetup const& es)
	{
		this->reset(1);
		//	Reset things per LS.
		//	But at least 100 events in LS
//		if (_si.evsPerLS>100)
//			this->reset(1);
	}

	/* virtual */ void HcalDQSource::endLuminosityBlock(
			edm::LuminosityBlock const& lb, edm::EventSetup const& es)
	{}

	//	extract Event Calibration Type from FEDs
	void HcalDQSource::extractCalibType(edm::Event const&e)
	{
		edm::Handle<FEDRawDataCollection> craw;
		INITCOLL(_labels["RAW"], craw);

		//	for now
		int badFEDs = 0;
		std::vector<unsigned int> types(8,0);
		for (std::vector<int>::const_iterator it=_mi.feds.begin();
				it!=_mi.feds.end(); ++it)
		{
			FEDRawData const& fd = craw->FEDData(*it);
			if (fd.size()<RAWDATASIZE_CALIB)
			{badFEDs++; continue;}
			int cval = (int)((HcalDCCHeader const*)(fd.data()))->getCalibType();
			if (cval>MAXCALIBTYPE)
				warn_("Unexpected Calib Type in FED " + 
						boost::lexical_cast<std::string>(*it));
			types[cval]++;
		}

		unsigned int max = 0;
		for (unsigned int ic=0; ic<types.size(); ic++)
		{
			if (types[ic]>max)
			{_si.currentCalibType=ic; max=types[ic];}
			if (max==_mi.feds.size())
				break;
		}

		if (max!=(_mi.feds.size()-badFEDs))
			warn_("Conflictings Calibration Types found. Assigning " + 
					boost::lexical_cast<std::string>(_si.currentCalibType));
	}

	//	Check if calibration type set is among allowed
	bool HcalDQSource::isAllowedCalibType()
	{
		for (std::vector<int>::const_iterator it=_mi.calibTypesAllowed.begin();
				it!=_mi.calibTypesAllowed.end(); ++it)
			if (_si.currentCalibType==*it)
				return true;

		return false;
	}

	//	reset
	/* virtual */ void HcalDQSource::reset(int const periodflag)
	{
		//	Collection Class will determine itself who needs a reset and when
		_mes.reset(periodflag);

		if (periodflag==0)
		{
			//	each event reset
		}
		else if (periodflag==1)
		{
			//	each LS reset
			_si.evsPerLS=0;
		}
	}

}









