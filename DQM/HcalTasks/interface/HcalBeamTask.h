#ifndef HCALBEAMTASK_H
#define HCALBEAMTASK_H

/*
 *	file:			HcalBeamTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalBeamTask : public hcaldqm::HcalDQSource
{
	public:
		HcalBeamTask(edm::ParameterSet const&);
		virtual ~HcalBeamTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

//	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
};

#endif
