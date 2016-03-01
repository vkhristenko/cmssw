#ifndef LEDTask_h
#define LEDTask_h

/*
 *	file:			LEDTask.h
 *	Author:			Viktor Khristenko
 *	Date:			16.10.2015
 */

#include "DQM/HcalCommon/interface/DQTask.h"
#include "DQM/HcalCommon/interface/Utilities.h"
#include "DQM/HcalCommon/interface/ElectronicsMap.h"
#include "DQM/HcalCommon/interface/ContainerCompact.h"
#include "DQM/HcalCommon/interface/ContainerXXX.h"
#include "DQM/HcalCommon/interface/Container1D.h"
#include "DQM/HcalCommon/interface/Container2D.h"
#include "DQM/HcalCommon/interface/ContainerSingle2D.h"
#include "DQM/HcalCommon/interface/ContainerSingleProf2D.h"
#include "DQM/HcalCommon/interface/ContainerProf1D.h"
#include "DQM/HcalCommon/interface/ContainerProf2D.h"

using namespace hcaldqm;

class LEDTask : public DQTask
{
	public:
		LEDTask(edm::ParameterSet const&);
		virtual ~LEDTask()
		{}

		virtual void bookHistograms(DQMStore::IBooker&,
			edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&)
		{this->_dump();}

	protected:
		//	funcs
		virtual void _process(edm::Event const&, edm::EventSetup const&);
		virtual void _resetMonitors(UpdateFreq);
		virtual bool _isApplicable(edm::Event const&);
		virtual void _dump();

		//	tags and tokens
		edm::InputTag	_tagHBHE;
		edm::InputTag	_tagHO;
		edm::InputTag	_tagHF;
		edm::InputTag	_tagTrigger;
		edm::EDGetTokenT<HBHEDigiCollection> _tokHBHE;
		edm::EDGetTokenT<HODigiCollection> _tokHO;
		edm::EDGetTokenT<HFDigiCollection> _tokHF;
		edm::EDGetTokenT<HcalTBTriggerData> _tokTrigger;

		//	emap
		HcalElectronicsMap const* _emap;
		electronicsmap::ElectronicsMap _emaphash;

		//	Cuts
		double _lowHBHE;
		double _lowHO;
		double _lowHF;

		//	Compact
		ContainerXXX _cSignals_DChannel;
		ContainerXXX _cTiming_DChannel;

		//	1D
		Container1D		_cSignalMean_Subdet;
		Container1D		_cSignalRMS_Subdet;
		Container1D		_cTimingMean_Subdet;
		Container1D		_cTimingRMS_Subdet;

		//	Prof1D
		ContainerProf1D	_cShapeCut_FEDSlot;

		//	2D timing/signals
		ContainerProf2D		_cSignalMean_depth;
		ContainerProf2D		_cSignalRMS_depth;
		ContainerProf2D		_cTimingMean_depth;
		ContainerProf2D		_cTimingRMS_depth;

		ContainerSingleProf2D _cTimingVME;
		ContainerSingleProf2D _cSignalVME;
		ContainerSingle2D _cOccupancyVME;
		ContainerSingleProf2D _cTiminguTCA;
		ContainerSingleProf2D _cSignaluTCA;
		ContainerSingle2D _cOccupancyuTCA;

		//	Bad Quality and Missing Channels
		Container2D		_cMissing_depth;
};

#endif







