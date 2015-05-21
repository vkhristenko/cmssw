#ifndef HCALDIGICLIENT_H
#define HCALDIGICLIENT_H

/*
 *	file:				HcalDigiClient.h
 *	Author:				Viktor Khristenko
 *	Start Date:			16/05/2015
 */

#include "DQM/HcalCommon/interface/HcalDQClient.h"
#include "DQM/HcalCommon/interface/HcalCommonHeaders.h"

class HcalDigiClient : public hcaldqm::HcalDQClient
{
	public:
		HcalDigiClient(edm::ParameterSet const&);
		virtual ~HcalDigiClient();

		virtual void doWork(DQMStore::IGetter&,
				edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void doWork(DQMStore::IBooker&, DQMStore::IGetter&);

	private:
		//	Digi Client specific things
};

#endif
