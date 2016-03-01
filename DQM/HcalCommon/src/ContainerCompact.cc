
#include "DQM/HcalCommon/interface/ContainerCompact.h"

namespace hcaldqm
{
	using namespace constants;

	/* virtual */ void ContainerCompact::fill(HcalDetId const& did, double x)
	{
		if (x==GARBAGE_VALUE)
			return;

		int isubdet = did.subdet()-1;
		int iiphi = did.iphi()-1;
		int iieta = did.ieta()<0 ? abs(did.ieta())-IETA_MIN : 
			did.ieta()-IETA_MIN+IETA_NUM/2;
		int idepth = did.depth()-1;

		_data[isubdet][iiphi][iieta][idepth]._x1 += x;
		_data[isubdet][iiphi][iieta][idepth]._x2 += x*x;
		_data[isubdet][iiphi][iieta][idepth]._n++;
	}

	/* virtual */ void ContainerCompact::dump(Container1D* c, bool q)
	{
		for (int idet=0; idet<SUBDET_NUM; idet++)
		{
			HcalSubdetector subd = HcalEmpty;
			if (idet+1==HB)
				subd = HcalBarrel;
			else if (idet+1==HE)
				subd = HcalEndcap;
			else if (idet+1==HO)
				subd = HcalOuter;
			else 
				subd = HcalForward;
			for (int iiphi=0; iiphi<IPHI_NUM; iiphi++)
				for (int iieta=0; iieta<IETA_NUM; iieta++)
					for (int id=0; id<DEPTH_NUM; id++)
					{
						Compact tmp(_data[idet][iiphi][iieta][id]);
						if (tmp._n<=0)
							continue;

						int iphi = iiphi+1;
						int ieta = iieta<IETA_NUM/2 ? -(iieta+IETA_MIN) :
							iieta-IETA_NUM/2+IETA_MIN;
						int d = id+1;
						HcalDetId did(subd, ieta, iphi, d);
						double mean= tmp._x1/tmp._n;
						double rms = sqrt(tmp._x2/tmp._n - 
							mean*mean);
						if (q)
							c->fill(did, mean);
						else
							c->fill(did, rms);
					}
		}
	}
}


