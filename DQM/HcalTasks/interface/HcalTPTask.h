#ifndef HCALTPTASK_H
#define HCALTPTASK_H

/*
 *	file:			HcalTPTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalTPTask : public hcaldqm::HcalDQSource
{
	public:
		HcalTPTask(edm::ParameterSet const&);
		virtual ~HcalTPTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

//	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
};

#endif
