#include "DQM/HcalCommon/interface/Compact.h"

namespace hcaldqm
{
	namespace compact
	{
		std::ostream& operator<<(std::ostream& o, Compact const& x)
		{
			return o << "sum="<< x._x1 << "  sum2=" 
				<< x._x2 << "  N=" << x._n;
		}
	}
}
