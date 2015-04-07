#ifndef HCALRAWTASK_H
#define HCALRAWTASK_H

/*
 *	file:			HcalRawTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalRawTask : public hcaldqm::HcalDQSource
{
	public:
		HcalRawTask(edm::ParameterSet const&);
		virtual ~HcalRawTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

//	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
};

#endif
