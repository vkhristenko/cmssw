#ifndef PedestalTask_h
#define PedestalTask_h

/*
 *	file:			PedestalTask.h
 *	Author:			Viktor Khristenko
 *	Date:			16.10.2015
 */

#include "DQM/HcalCommon/interface/DQTask.h"
#include "DQM/HcalCommon/interface/Utilities.h"
#include "DQM/HcalCommon/interface/Container1D.h"
#include "DQM/HcalCommon/interface/Container2D.h"
#include "DQM/HcalCommon/interface/ContainerProf1D.h"
#include "DQM/HcalCommon/interface/ContainerProf2D.h"
#include "DQM/HcalCommon/interface/ContainerCompact.h"
#include "DQM/HcalCommon/interface/ContainerXXX.h"

using namespace hcaldqm;
class PedestalTask : public DQTask
{
	public:
		PedestalTask(edm::ParameterSet const&);
		virtual ~PedestalTask()
		{}

		virtual void bookHistograms(DQMStore::IBooker&,
			edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);

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
		HcalElectronicsMap const*	_emap;
		HcalPedestal const*			_condPeds;

		ContainerXXX		_cPed;
		ContainerXXX		_cPedRef;

		//	1D Means/RMSs
		Container1D		_cMean_Subdet;
		Container1D		_cRMS_Subdet;

		//	1D Means/RMSs Conditions DB comparison
		Container1D		_cMeanDBRef_Subdet;
		Container1D		_cRMSDBRef_Subdet;

		//	2D
		ContainerProf2D		_cMean_depth;
		ContainerProf2D		_cRMS_depth;
		
		//	with DB Conditions comparison
		ContainerProf2D		_cMeanDBRef_depth;
		ContainerProf2D		_cRMSDBRef_depth;

		//	Missing + Bad Quality
		Container2D		_cMissing_depth;
		Container2D		_cMeanBad_depth;
		Container2D		_cRMSBad_depth;
};

#endif







