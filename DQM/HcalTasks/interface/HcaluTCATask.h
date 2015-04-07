#ifndef HCALUTCATASK_H
#define HCALUTCATASK_H

/*
 *	file:			HcaluTCATask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcaluTCATask : public hcaldqm::HcalDQSource
{
	public:
		HcaluTCATask(edm::ParameterSet const&);
		virtual ~HcaluTCATask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);

//	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
};

#endif
