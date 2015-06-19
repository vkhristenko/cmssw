#ifndef HCALPHASESCANTASK_H
#define HCALPHASESCANTASK_H

/*
 *	file:			HcalPhaseScanTask.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQSource.h"

class HcalPhaseScanTask : public hcaldqm::HcalDQSource
{
	public:
		HcalPhaseScanTask(edm::ParameterSet const&);
		virtual ~HcalPhaseScanTask();

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
		void specialize(Hit const& hit, std::string const&,
				int const wtw=1);

		//	Define the processors
		DEFPROCESSOR(HFDigiCollection, HFDataFrame);

		//	Call HF Phase Scan
		void hf(HFDataFrame const&);
};

#endif
