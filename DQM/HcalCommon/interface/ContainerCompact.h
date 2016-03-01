#ifndef ContainerCompact_h
#define ContainerCompact_h

/*
 *	file:			ContainerCompact.h
 *	Author:			Viktor Khristenko
 *
 *	Description:
 *		No Type Usage adjustment - only use in histo-like mode
 */

#include "DQM/HcalCommon/interface/Container1D.h"
#include "DQM/HcalCommon/interface/Compact.h"

namespace hcaldqm
{
	using namespace constants;
	using namespace compact;	
	class ContainerCompact
	{
		public:
			ContainerCompact()
			{}
			virtual ~ContainerCompact() {}

			//	fills
			virtual void fill(HcalDetId const&, double);
			virtual void dump(Container1D*, bool);

		protected:
			Compact			_data[SUBDET_NUM][IPHI_NUM][IETA_NUM][DEPTH_NUM];
	};
}

#endif


