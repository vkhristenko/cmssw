#include "DQM/HcalClients/interface/HcalDigiClient.h"

HcalDigiClient::HcalDigiClient()
	:
{}

/* virtual */ HcalDigiClient::~HcalDigiClient()
{}

/* virtual */ void HcalDigiClient::doWork(DQMStore::IGetter& ig,
		edm::LuminosotyBlock const& ls, edm::EventSetup const& es)
{

}

/* virtual */ void HcalDigiClient::doWork(DQMStore::IBooker &ib,
		DQMStore::IGetter &ig)
{

}
