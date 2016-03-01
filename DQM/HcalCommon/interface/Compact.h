#ifndef Compact_h
#define Compact_h

/**
 *	file:		Compact.h
 *	Auhtor:		Viktor Khristenko
 */

#include "DQM/HcalCommon/interface/HcalCommonHeaders.h"
#include "DQM/HcalCommon/interface/Constants.h"

namespace hcaldqm
{
	namespace compact
	{
		enum CompactUsageType
		{
			fHistogram = 0, //	always fill and sum only (for mean/rms usage)
			fSingleValued = 1, // use only 1 value with set
			fDoubleValued = 2, // use both values with set
			nCompactUsageType = 3
		};

		struct Compact 
		{
			Compact(){_x1=0; _x2=0; _n=0;}
			void reset() {_x1=0; _x2=0; _n=0;}
			
			double mean() {return _n>0?_x1/_n:constants::GARBAGE_VALUE;}
			double rms() 
			{
				if (mean()==constants::GARBAGE_VALUE)
					return constants::GARBAGE_VALUE;

				double m = this->mean();
				return sqrt(_x2/_n-m*m);
			}
			std::pair<double, double> getValues(CompactUsageType utype)
			{
				std::pair<double, double> p;
				if (utype==fHistogram)
					p = std::make_pair<double, double>(mean(), rms());
				else
					p = std::make_pair<double, double>((double)_x1, (double)_x2);

				return p;
			}

			double		_x1;
			double		_x2;
			uint32_t	_n;
		};
		std::ostream& operator<<(std::ostream&, Compact const&);
	}
}

#endif
