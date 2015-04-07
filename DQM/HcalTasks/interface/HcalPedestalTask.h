#ifndef HCALPEDESTALTASK_H
#define HCALPEDESTALTASK_H

/*
 *	file:			HcalPedestalTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalPedestalTask : public hcaldqm::HcalDQSource
{
	public:
		HcalPedestalTask(edm::ParameterSet const&);
		virtual ~HcalPedestalTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

//	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
};

#endif
