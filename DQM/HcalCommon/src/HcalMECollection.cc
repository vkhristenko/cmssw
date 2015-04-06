
#include "DQM/HcalCommon/interface/HcalMECollection.h"

namespace hcaldqm
{
	MEInfo::MEInfo(edm::ParameterSet const&ps) :
		meps(ps)
	{}

	HcalMECollection::HcalMECollection(edm::ParameterSet const& ps)
		: _ps(ps)
	{}

	HcalMECollection::~HcalMECollection()
	{
		//	Clean the ME map
	}

	//	for booking Monitor Elements based on PSet
	void HcalMECollection::book(DQMStore::IBooker &ib)
	{
		std::vector<std::string> const& meNames(_ps.getParameterNames());
		for (std::vector<std::string>::const_iterator it=meNames.begin();
				it!=meNames.end(); ++it)
		{
			std::string meName = *it;
			MEInfo meinfo(_ps.getUntrackedParameterSet(*it)); 
			meinfo.setName(meName);
			doBook(ib, meinfo);
		}
	}

	void HcalMECollection::doBook(DQMStore::IBooker &ib,
			MEInfo const& info)
	{
		std::string path = info.getPS().getUntrackedParameter<std::string>(
				"path");
		std::string kind = info.getPS().getUntrackedParameter<std::string>(
				"kind");
		MonitorElement *me;

		ib.setCurrentFolder(path);
		if  (kind=="INT")
//			case MonitorElement::DQM_KIND_INT:
			me = ib.bookInt(info.getName());
		else if (kind=="REAL")
//			case MonitorElement::DQM_KIND_REAL:
			me = ib.bookFloat(info.getName());
		else if (kind=="STRING")
//			case MonitorElement::DQM_KIND_STRING :
			me = ib.bookString(path, info.getName());
		else if (kind=="TH1D")
//			case MonitorElement::DQM_KIND_TH1D:
			me = this->create1D(ib, info);
		else if (kind=="TH2D")
//			case MonitorElement::DQM_KIND_TH2D:
			me = this->create2D(ib, info);
		else if (kind=="PROF")
//			case MonitorElement::DQM_KIND_TPROF:
			me = this->createProf(ib, info);
		else if (kind=="PROF2D")
//			case MonitorElement::DQM_KIND_TPROF2D:
			me = this->createProf2D(ib,info);
		else
			me = NULL;

		std::string key = info.getName();
		std::pair<MEMap::iterator, bool> r = _meMap.insert(key, me);
		if (r.second)
			return;
	}

	//	for retirieving Monitor Elements based on PSet
	void HcalMECollection::retrieve(DQMStore::IGetter &ig)
	{}

	MonitorElement& HcalMECollection::operator[](std::string name)
	{
		return _meMap[name];
	}
	
	MonitorElement* HcalMECollection::create1D(DQMStore::IBooker &ib, 
			MEInfo const& info)
	{
		edm::ParameterSet const& axisps = 
			info.getPS().getUntrackedParameterSet("xaxis");
		std::string desc = info.getPS().getUntrackedParameter<std::string>(
				"desc");
		MEAxis xaxis;
		xaxis.edges = GETPAR(axisps, bool, "edges");
		xaxis.nbins = GETPAR(axisps, int, "nbins");
		xaxis.title = GETPAR(axisps, std::string, "title");
		MonitorElement *me = NULL;
		if (!xaxis.edges)
		{
			xaxis.min = GETPAR(axisps, double, "min");
			xaxis.max = GETPAR(axisps, double, "max");
			me = ib.book1D(info.getName(), desc, xaxis.nbins, 
					xaxis.min, xaxis.max);
		}
		else
		{
			xaxis.bins = 
				axisps.getUntrackedParameter<std::vector<double> >("bins").data();
			me = ib.book1D(info.getName(), desc, xaxis.nbins, 
					(float const*)xaxis.bins);
		}
		
		me->setAxisTitle(xaxis.title);
		return me;
	}

	MonitorElement* HcalMECollection::createProf2D(DQMStore::IBooker &ib, 
			MEInfo const& info)
	{return NULL;}

	MonitorElement* HcalMECollection::createProf(DQMStore::IBooker &ib, 
			MEInfo const& info)
	{	
		edm::ParameterSet const &xaxisps = 
			info.getPS().getUntrackedParameterSet("xaxis");
		edm::ParameterSet const &yaxisps = 
			info.getPS().getUntrackedParameterSet("yaxis");
		std::string desc = info.getPS().getUntrackedParameter<std::string>(
				"desc");

		MEAxis xaxis, yaxis;
		xaxis.edges = GETPAR(xaxisps, bool, "edges");
		xaxis.nbins = GETPAR(xaxisps, int, "nbins");
		yaxis.wnbins = GETPAR(yaxisps, bool, "wnbins");
		yaxis.edges = false;
		xaxis.title = GETPAR(xaxisps, std::string, "title");
		yaxis.title = GETPAR(yaxisps, std::string, "title");
		MonitorElement *me = NULL;

		if (!xaxis.edges && yaxis.wnbins)
		{
			xaxis.min = GETPAR(xaxisps, double, "min");
			xaxis.max = GETPAR(xaxisps, double, "max");
			yaxis.min = GETPAR(yaxisps, double, "min");
			yaxis.max = GETPAR(yaxisps, double, "max");
			yaxis.nbins = GETPAR(yaxisps,int, "nbins");
			me = ib.bookProfile(info.getName(), desc, xaxis.nbins, 
					xaxis.min, xaxis.max, yaxis.nbins, yaxis.min, yaxis.max);
		}
		else if (!xaxis.edges && !yaxis.wnbins)
		{
			xaxis.min = GETPAR(xaxisps, double, "min");
			xaxis.max = GETPAR(xaxisps, double, "max");
			yaxis.min = GETPAR(yaxisps, double, "min");
			yaxis.max = GETPAR(yaxisps, double, "max");
			me = ib.bookProfile(info.getName(), desc, xaxis.nbins, xaxis.min,
					xaxis.max, yaxis.min, yaxis.max);
		}
		else if (xaxis.edges && yaxis.wnbins)
		{
			xaxis.bins = 
				xaxisps.getUntrackedParameter<std::vector<double> >(
						"bins").data();
			yaxis.nbins = GETPAR(yaxisps,int, "nbins");
			yaxis.min = GETPAR(yaxisps, double, "min");
			yaxis.max = GETPAR(yaxisps, double, "max");
			me = ib.bookProfile(info.getName(), desc, xaxis.nbins, 
					xaxis.bins,	yaxis.nbins, yaxis.min, yaxis.max);
		}
		else if (xaxis.edges && !yaxis.wnbins)
		{
			xaxis.bins = 
				xaxisps.getUntrackedParameter<std::vector<double> >(
						"bins").data();
			yaxis.max = GETPAR(yaxisps, double, "max");
			yaxis.min = GETPAR(yaxisps, double, "min");
			me = ib.bookProfile(info.getName(), desc, xaxis.nbins, 
					xaxis.bins,	yaxis.min, yaxis.max);
		}

		me->setAxisTitle(xaxis.title, 1);
		me->setAxisTitle(yaxis.title, 2);
		return me;
	}

	MonitorElement* HcalMECollection::create2D(DQMStore::IBooker &ib, 
			MEInfo const& info)
	{
		edm::ParameterSet const &xaxisps = 
			info.getPS().getUntrackedParameterSet("xaxis");
		edm::ParameterSet const &yaxisps = 
			info.getPS().getUntrackedParameterSet("yaxis");
		std::string desc = info.getPS().getUntrackedParameter<std::string>(
				"desc");

		MEAxis xaxis, yaxis;
		xaxis.edges = GETPAR(xaxisps, bool, "edges");
		xaxis.nbins = GETPAR(xaxisps, int, "nbins");
		yaxis.edges = GETPAR(yaxisps, bool, "edges");
		yaxis.nbins = GETPAR(yaxisps, int, "nbins");
		xaxis.wnbins = false;
		yaxis.wnbins = false;
		xaxis.title = GETPAR(xaxisps, std::string, "title");
		yaxis.title = GETPAR(yaxisps, std::string, "title");
		MonitorElement *me = NULL;

		if (!yaxis.edges && !xaxis.edges)
		{		
			xaxis.min = GETPAR(xaxisps, double, "min");
			xaxis.max = GETPAR(xaxisps, double, "max");
			yaxis.min = GETPAR(yaxisps, double, "min");
			yaxis.max = GETPAR(yaxisps, double, "max");
			me = ib.book2D(info.getName(), desc, xaxis.nbins, xaxis.min, 
					xaxis.max, yaxis.nbins, yaxis.min, yaxis.max);
		}
		else if (yaxis.edges && xaxis.edges)
		{
			xaxis.bins = 
				xaxisps.getUntrackedParameter<std::vector<double> >(
						"bins").data();
			yaxis.bins = 
				yaxisps.getUntrackedParameter<std::vector<double> >(
						"bins").data();
			me = ib.book2D(info.getName(), desc, xaxis.nbins, 
					(const float*)xaxis.bins, yaxis.nbins, 
					(const float*)yaxis.bins);
		}
		
		me->setAxisTitle(xaxis.title, 1);
		me->setAxisTitle(yaxis.title, 2);
		return me;
	}
}







