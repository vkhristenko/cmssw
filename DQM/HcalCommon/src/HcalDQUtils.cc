//	cmssw includes
#include "DQM/HcalCommon/interface/HcalDQUtils.h"

//	system includes
#include <string>
#include <vector>

namespace hcaldqm
{
	//	Book all the Input Labels for Collections needed
	Labels::Labels(edm::ParameterSet const& ps)
	{
		std::vector<std::string> const& lNames = ps.getParameterNames();
		for (std::vector<std::string>::const_iterator it=lNames.begin();
				it!=lNames.end(); ++it)
		{
			edm::InputTag label = ps.getUntrackedParameter<edm::InputTag>(
					*it);
			std::pair<LabelMap::iterator, bool> r = 
				_labels.insert(std::make_pair(*it, label));
			if (r.second)
				continue;
		}
	}
}











