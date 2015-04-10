#ifndef HCALDIGITASK_H
#define HCALDIGITASK_H

/*
 *	file:			HcalDigiTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalDigiTask : public hcaldqm::HcalDQSource
{
	public:
		HcalDigiTask(edm::ParameterSet const&);
		virtual ~HcalDigiTask();

		virtual void doWork(edm::Event const&e,
				edm::EventSetup const& es);
	
	private:
		//	declare the template for specializing
		template<typename Hit>
		void specialize(Hit const& hit, std::string const&);
		
		//	define and initialize the collection processors
		DEFPROCESSOR(HBHEDigiCollection, HBHEDataFrame);
		DEFPROCESSOR(HODigiCollection, HODataFrame);
		DEFPROCESSOR(HFDigiCollection, HFDataFrame);

	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
		int _ornMsgTime;
};

#endif
