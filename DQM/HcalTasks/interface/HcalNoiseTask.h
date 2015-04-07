#ifndef HCALNOISETASK_H
#define HCALNOISETASK_H

/*
 *	file:			HcalNoiseTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalNoiseTask : public hcaldqm::HcalDQSource
{
	public:
		HcalNoiseTask(edm::ParameterSet const&);
		virtual ~HcalNoiseTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

//	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
};

#endif
