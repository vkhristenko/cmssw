#ifndef DQM_HcalTasks_DigiRunSummary_h
#define DQM_HcalTasks_DigiRunSummary_h

#include "DQM/HcalCommon/interface/DQClient.h"
#include "DQM/HcalCommon/interface/ElectronicsMap.h"

namespace hcaldqm
{
	class DigiRunSummary : public DQClient
	{
		public:
			DigiRunSummary(std::string const&, std::string const&,
				edm::ParameterSet const&);
			virtual ~DigiRunSummary() {}

			virtual void beginRun(edm::Run const&, edm::EventSetup const&);
			virtual void endLuminosityBlock(DQMStore::IBooker&,
				DQMStore::IGetter&, edm::LuminosityBlock const&,
				edm::EventSetup const&);
			virtual std::vector<flag::Flag> endJob(
				DQMStore::IBooker&, DQMStore::IGetter&);

		protected:
			Container2D _cDead_depth;
			Container2D _cDead_FEDVME;
			Container2D _cDead_FEDuTCA;

			//	flag enum
			enum DigiFlag
			{
				fDead = 0,
				fUniSlotHF = 1,
				fDigiSize = 2,
				nDigiFlag = 3
			};
	};
}

#endif
