#ifndef HCALDQMONITOR_H
#define HCALDQMONITOR_H

/*
 *	file:			HcalDQMonitor.h
 *	Author:			Viktor Khristenko
 *	Start Date:		03/04/2015
 *
 *	TODO:
 *		1) Other module-generic parameters???
 *		2) Other module-genetic functionality???
 */

#include "DQM/HcalCommon/interface/HcalDQUtils.h"
#include "DQM/HcalCommon/interface/HcalCommonHeaders.h"

#include <string>


namespace hcaldqm
{
	/*
	 * Structs below are for the future
	 */
	enum EventType
	{
		iNormal,
		iCalibration,
		nEventType
	};

	enum RunType
	{
		iOnline,
		iOffline,
		iLocal,				//	a la DetDiag
		nRunType
	};

	enum ModuleType
	{
		iClient,
		iSource,
		nModuleType
	};

	struct ModuleInfo
	{
		std::string			type;
		std::string			runType;
		std::vector<int>	calibTypesAllowed;
		std::vector<int>	feds;
		std::string			name;
		int					debug;
	};

	/*
	 *	HcalDQMonitor Class: Commont Base Class for DQSources and DQClients
	 *	and all the modules to be used
	 */
	class HcalDQMonitor
	{
		public:
			HcalDQMonitor(edm::ParameterSet const&);
			virtual ~HcalDQMonitor();

			ModuleInfo info() const {return _mi;}

			inline void throw_(std::string const msg) const
			{
				throw cms::Exception("HcalDQM") << _mi.name << "::"	<< msg;
			}

			inline void warn_(std::string const msg) const
			{
				edm::LogWarning("HcalDQM") << _mi.name << "::" << msg;
			}

			//	For the case when msg and msg1 cannot be concatenated
			inline void warn_(std::string const msg, std::string const msg1) const
			{
				edm::LogWarning("HcalDQM") << _mi.name << "::" << msg << msg1;
			}

			inline void info_(std::string const msg) const
			{
				edm::LogInfo("HcalDQM") << _mi.name << "::" << msg;
			}

			inline void debug_(std::string const msg) const
			{
				std::cout << "%MSG" << std::endl;
				std::cout << "%MSG-d HcalDQM::" << _mi.name << "::" << msg;
				std::cout << std::endl;
			}

		protected:
			ModuleInfo		_mi;
			Labels			_labels;
	};

}

#define INITCOLL(LABEL, COLL) \
	if (!(e.getByLabel(LABEL, COLL))) \
		warn_("Collection" #COLL " not available", \
				"  " + LABEL.label()+ "  " + LABEL.instance())

#endif






