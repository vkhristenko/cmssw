#ifndef HCALHFPHASESCANTASK_H
#define HCALHFPHASESCANTASK_H

/*
 *	file:			HcalHFPhaseScanTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalHFPhaseScanTask : public hcaldqm::HcalDQSource
{
	public:
		HcalHFPhaseScanTask(edm::ParameterSet const&);
		virtual ~HcalHFPhaseScanTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
//		DEFPROCESSOR(HBHERecHitCollection, HBHERecHit);
//		DEFPROCESSOR(HORecHitCollection, HORecHit);
//		DEFPROCESSOR(HFRecHitCollection, HFRecHit);
//		DEFPROCESSOR(HBHEDigiCollection, HBHEDataFrame);
//		DEFPROCESSOR(HODigiCollection, HODataFrame);
//		DEFPROCESSOR(HFDigiCollection, HFDataFrame);

		template<typename Hit>
		void specialize(Hit const&, std::string const&);
};

#endif
