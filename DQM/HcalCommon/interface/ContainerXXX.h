#ifndef ContainerXXX_h
#define ContainerXXX_h

/*
 *	file:			ContainerXXX.h
 *	Author:			Viktor Khristenko
 *
 *	Description:
 *		1D Compact Container
 */

#include "DQM/HcalCommon/interface/Container1D.h"
#include "DQM/HcalCommon/interface/Compact.h"
#include "DQM/HcalCommon/interface/Logger.h"

namespace hcaldqm
{
	using namespace constants;
	using namespace compact;

	//	comparison and quality function typedefs
	//	plus some standard types
	typedef double (*comparison_function)(double, double);
	typedef bool (*quality_function)(double, double, double);
	double ratio(double, double);
	double diff(double, double);

	//	quality funcs - 2 values should be within 20%(LED,Laser) of each other
	//	configurable, 10%(PED)
	//	true is good, false is bad
	bool pedestal_quality(double, double, double);
	bool pedestal_quality2(double, double, double);
	bool led_quality(double, double, double);
	bool laser_quality(double, double, double);

	class QualityWrapper
	{
		public:
			QualityWrapper(quality_function f, double thr=0.2) : func(f), 
			_thr(thr){}
			virtual ~QualityWrapper() {}

			bool quality(double x, double y){return (*func)(x, y, _thr);}

		private:
			quality_function func;
			double	_thr;
	};

	class ContainerXXX
	{
		public:
			ContainerXXX() {}
			ContainerXXX(hashfunctions::HashType ht, CompactUsageType 
				ut=fHistogram)
				: _hashmap(ht), _usetype(ut)
			{}
			virtual ~ContainerXXX() {_cmap.clear();}

			//	initializer
			virtual void initialize(hashfunctions::HashType,
				CompactUsageType=fHistogram, int debug=0);

			//	book
			virtual void book(HcalElectronicsMap const*);

			//	fills - only in histo mode
			virtual void fill(HcalDetId const&, double);
			virtual void fill(HcalElectronicsId const&, double);
			virtual void fill(HcalTrigTowerDetId const&, double);

			//	sets - only in non-histo mode
			virtual void set(HcalDetId const&, double, double y=0);
			virtual void set(HcalElectronicsId const&, double, double y=0);
			virtual void set(HcalTrigTowerDetId const&, double, double y=0);

			//	get the number of entries
			virtual uint32_t getEntries(HcalDetId const&);
			virtual uint32_t getEntries(HcalElectronicsId const&);
			virtual uint32_t getEntries(HcalTrigTowerDetId const&);

			//	dump all the contents into the Container
			//	for mean bool=true
			virtual void dump(Container1D*, bool q=true);
			virtual void dump(std::vector<Container1D*> const&, bool q=true);
			//	dump the results + dump if absent
			virtual void dump(Container1D*, Container1D *, bool q=true);
			virtual void dump(std::vector<Container1D*> const&, 
				Container1D*, bool q=true);
			//virtual void dump(HcalElectronicsMap const*, Container1D*,
//				bool q=true);
			//virtual void dump(HcalElectronicsMap const*, ContainerSingle1D*,
//				bool q=true);
			//virtual void dump(HcalElectronicsMap const*, ContainerSingle2D*,
//				bool q=true);

			//	get all the contents from the Container
			virtual void load(Container1D*);
			virtual void load(Container1D*, Container1D*);

			//	loads from dbs
			virtual void loadPedestals(edm::ESHandle<HcalDbService> const&);

			//	compare 2 sets of data
			//	1-st Container1D is to where to dump the comparison
			//	2-nd Container1D is to where to dump the non-present channels
			//	3-rd COntainer1D is to where to dump the quality faulire channels.
			//	
			virtual void compare(ContainerXXX const&, Container1D*,
				Container1D*, Container1D*, comparison_function, 
				QualityWrapper, bool q=true);
			virtual void compare(ContainerXXX const&, 
				std::vector<Container1D*> const&, Container1D*,
				Container1D*, comparison_function, QualityWrapper, 
				bool q=true);

			//	reset
			virtual void reset();

			//	size
			virtual uint32_t size();

			//	print
			virtual void print();

		protected:
			typedef boost::unordered_map<uint32_t, Compact> CompactMap;
			CompactMap				_cmap;
			mapper::HashMapper		_hashmap;
			CompactUsageType		_usetype;
			Logger					_logger;
	};
}

#endif


