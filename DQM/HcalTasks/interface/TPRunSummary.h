#ifndef DQM_HcalTasks_TPRunSummary_h
#define DQM_HcalTasks_TPRunSummary_h

#include "DQM/HcalCommon/interface/DQClient.h"
#include "DQM/HcalCommon/interface/ElectronicsMap.h"

namespace hcaldqm
{
	class TPRunSummary : public DQClient
	{
		public:
			TPRunSummary(std::string const&, std::string const&,
				edm::ParameterSet const&);
			virtual ~TPRunSummary() {}

			virtual void beginRun(edm::Run const&, edm::EventSetup const&);
			virtual void endLuminosityBlock(DQMStore::IBooker&,
				DQMStore::IGetter&, edm::LuminosityBlock const&,
				edm::EventSetup const&);
			virtual std::vector<flag::Flag> endJob(
				DQMStore::IBooker&, DQMStore::IGetter&);

		protected:
	};
}

#endif
