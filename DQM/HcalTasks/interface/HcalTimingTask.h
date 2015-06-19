#ifndef HCALTIMINGTASK_H
#define HCALTIMINGTASK_H

/*
 *	file:			HcalTimingTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalTimingTask : public hcaldqm::HcalDQSource
{
	public:
		HcalTimingTask(edm::ParameterSet const&);
		virtual ~HcalTimingTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

//	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
};

#endif
