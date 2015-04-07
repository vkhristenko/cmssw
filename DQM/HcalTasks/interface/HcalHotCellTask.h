#ifndef HCALHOTCELLTASK_H
#define HCALHOTCELLTASK_H

/*
 *	file:			HcalHotCellTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalHotCellTask : public hcaldqm::HcalDQSource
{
	public:
		HcalHotCellTask(edm::ParameterSet const&);
		virtual ~HcalHotCellTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

//	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
};

#endif
