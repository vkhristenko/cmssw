#ifndef HCALLASERTASK_H
#define HCALLASERTASK_H

/*
 *	file:			HcalLaserTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalLaserTask : public hcaldqm::HcalDQSource
{
	public:
		HcalLaserTask(edm::ParameterSet const&);
		virtual ~HcalLaserTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

//	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
};

#endif
