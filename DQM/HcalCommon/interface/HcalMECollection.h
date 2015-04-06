#ifndef HCALMECOLLECTION_H
#define HCALMECOLLECTION_H

/*
 *	file:			HcalMECollection.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalCommonHeaders.h"

#include <vector>
#include <string>
#include "boost/ptr_container/ptr_map.hpp"

namespace hcaldqm
{
	enum CollectionType
	{
		iBooker,
		iGetter,
		nCollectionType
	};

	class MEInfo
	{
		public:
			MEInfo(edm::ParameterSet const& ps);

			void setName(std::string n)				{name=n;}

			std::string const getName() const				{return name;}
			edm::ParameterSet const& getPS() const			{return meps;}

		private:
			std::string					name;
			edm::ParameterSet const&	meps;
	};

	//	Simple struct to detect the type of axis input
	struct MEAxis
	{
		bool edges, wnbins;
		std::string title;
		int nbins;
		double min,max;
		double *bins;
	};

	/*
	 *	HcalMECollection Class
	 *	Access to all Monitor Elements thru a dictionary
	 */
	class HcalMECollection
	{
		public:
			HcalMECollection(edm::ParameterSet const&);
			~HcalMECollection();

			//	Book MEs based on the PSet
			void book(DQMStore::IBooker&);

			//	Retrieve MEs based on PSet
			void retrieve(DQMStore::IGetter&);

			//	Simple getters
			MonitorElement& getME(std::string name) {return (*this)[name];}
			MonitorElement& operator[](std::string);

		private:
			//	do the actual Booking
			void doBook(DQMStore::IBooker&, MEInfo const&);

			//	for clarity - separate
			MonitorElement* create1D(DQMStore::IBooker&, MEInfo const&);
			MonitorElement* create2D(DQMStore::IBooker&, MEInfo const&);
			MonitorElement* createProf(DQMStore::IBooker&, MEInfo const&);
			MonitorElement* createProf2D(DQMStore::IBooker&, MEInfo const&);

		private:
			//	a Map: MEname -> ME*
			typedef boost::ptr_map<std::string, MonitorElement> MEMap;
			MEMap											_meMap;
			//	Parameter Set of MEs	
			edm::ParameterSet const&						_ps;
	};

}

#define GETPAR(PS, TYPE, NAME) \
	PS.getUntrackedParameter<TYPE>(NAME)

#endif




