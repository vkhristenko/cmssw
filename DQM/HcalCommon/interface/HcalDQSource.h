#ifndef HCALDQSOURCE_H
#define HCALDQSOURCE_H

/*
 *	file:			HcalDQSource.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 *
 *	TODO:
 *		1) Extracting the Calibration Type
 *		2) Other Source-specific functionality
 */

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"

#include "DQM/HcalCommon/interface/HcalMECollection.h"
#include "DQM/HcalCommon/interface/HcalDQMonitor.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <vector>

namespace hcaldqm
{
	/*
	 *	HcalDQSource Class - Base Class for DQSources
	 */
	class HcalDQSource : public DQMEDAnalyzer, public HcalDQMonitor
	{
		public:
			HcalDQSource(edm::ParameterSet const&);
			virtual ~HcalDQSource();

			//	Genetic doWork function for all DQSources
			//	Note: Different def from DQClients
			virtual void doWork(edm::Event const&e, 
					edm::EventSetup const& es) = 0;

			//	Functions which have to be reimplemented from DQMEDAnalyzer
			virtual void analyze(edm::Event const& e, edm::EventSetup const& es);
			virtual void bookHistograms(DQMStore::IBooker &ib, edm::Run const&,
					edm::EventSetup const&);
			virtual void dqmBeginRun(edm::Run const&, edm::EventSetup const&);

			virtual void beginLuminosityBlock(edm::LuminosityBlock const& ,
					edm::EventSetup const&);
			virtual void endLuminosityBlock(edm::LuminosityBlock const& ,
					edm::EventSetup const&);
			
		protected:
			//	Apply Reset/Update if neccessary
			//	Resets the contents of MEs
			//	periodflag: 0 for Event Reset and 1 for LS Reset
			virtual void reset(int const periodflag);

		protected:
			//	Functions specific for Sources, but generic to all of them
			void extractCalibType(edm::Event const&);
			bool isAllowedCalibType();

		protected:
			HcalMECollection		_mes;
	};
}

//	The use of this macro must be properly controlled!
//	Becaue the COLLECTIONTYPE and HITTYPE must go in accord with each other
#define DEFPROCESSOR(COLLECTIONTYPE, HITTYPE) \
	void process(COLLECTIONTYPE const& c, std::string const& nameRes) \
	{	\
		for (COLLECTIONTYPE::const_iterator it=c.begin(); it!=c.end(); ++it)	\
		{	\
			const HITTYPE hit = (const HITTYPE)(*it);	\
			if ((nameRes=="HB" && hit.id().subdet()!=HcalBarrel)	\
					|| (nameRes=="HE" && hit.id().subdet()!=HcalEndcap))	\
				continue;	\
			specialize<HITTYPE>(hit, nameRes);	\
		}	\
	}

//	Define a specializer
//#define DEFSPECIALIZER()

#endif

