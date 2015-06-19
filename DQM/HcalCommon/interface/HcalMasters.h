#ifndef HCALMASTERS_H
#define HCALMASTERS_H

/*
 *	file:			HcalMasters.h
 *	Author:			VK
 *	Description:
 *
 *	TODO:
 *
 */

#include "DQM/HcalCommon/interface/HcalCommonHeaders.h"

namespace hcaldqm
{
	/*
	 *	HcalDcsMaster Class: Helper Class to deal with anything 
	 *	related to the HCAL Detector:
	 *	1) DCS status
	 *	HcalMapMaster Class
	 *	2) Conditions if neccessary
	 *		- Logical Map
	 *		- Electronics Map
	 *	3) Channel Statuses
	 *	HcalTriggerMaster
	 *	4) Triggers
	 */
	class HcalDcsMaster
	{
		public:
			HcalDcsMaster() {}
			~HcalDcsMaster() {}
			
			void initDcs(DcsStatusCollection const&);
			inline bool getDcsStatus(int const i) const {return _dcs[i];}
			inline void setDcsStatus(int const i, bool status) 
			{_dcs[i]=status;}

		protected:
			bool _dcs[hcaldqm::constants::STD_NUMSUBS];
	};

	/*
	 *	HcalTriggerMaster Class: 
	 */
/*	class HcalTriggerMaster
	{
		public:
			HcalTriggerMaster();
			~HcalTriggerMaster();

			//	Initialize the Event Triggers
			void initL1GT(L1GlobalTriggerReadoutRecord const&);
			void initTBTrigger(HcalTBTriggerData const&);
			void initL1MuGMT(std::vector<L1MuGMTReadoutRecord> const&);

		protected:
			bool	_trPED;
			bool	_trLED;
			bool	_trLASER;
	};
*/

}

#endif








