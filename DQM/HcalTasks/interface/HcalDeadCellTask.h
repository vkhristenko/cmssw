#ifndef HCALDEADCELLTASK_H
#define HCALDEADCELLTASK_H

/*
 *	file:			HcalDeadCellTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalDeadCellTask : public hcaldqm::HcalDQSource
{
	public:
		HcalDeadCellTask(edm::ParameterSet const&);
		virtual ~HcalDeadCellTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

//	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
};

#endif
