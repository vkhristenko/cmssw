#ifndef HCALLEDTASK_H
#define HCALLEDTASK_H

/*
 *	file:			HcalLEDTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalLEDTask : public hcaldqm::HcalDQSource
{
	public:
		HcalLEDTask(edm::ParameterSet const&);
		virtual ~HcalLEDTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

//	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
};

#endif
