#ifndef HCALDQUTILS_H
#define HCALDQUTILS_H

/*	
 *	file:				HcalDigiTask.h
 *	Author:				Viktor Khristenko
 *	Start Date:			03/04/2015
 */

#include "DQM/HcalCommon/interface/HcalCommonHeaders.h"

#include <vector>
#include <string>
#include "boost/ptr_container/ptr_map.hpp"
#include "boost/container/map.hpp"

namespace hcaldqm
{

	//	a struct of labels
	struct Labels
	{
		Labels(edm::ParameterSet const&);		
		edm::InputTag& operator[](std::string s) 
		{
			return _labels[s];
		}

		typedef boost::container::map<std::string, edm::InputTag> LabelMap;
		LabelMap _labels;
	};
}

#endif
