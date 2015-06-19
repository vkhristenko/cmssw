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

		virtual void beginLuminosityBlock(edm::LuminosityBlock const&,
				edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&,
				edm::EventSetup const&);

		virtual void reset(int const);

	private:
		//	MEs Collection come from the base class
		//	Here, we only need module specific parameters
		template<typename Hit>
		void specialize(Hit const& hit1, Hit const& hit2, std::string const&);

		//	Define and Initialize Comparators
		DEFCOMPARATOR(HFDigiCollection, HFDataFrame);
};

#endif
