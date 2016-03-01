
#include "DQM/HcalCommon/interface/ContainerXXX.h"
#include <cmath>

namespace hcaldqm
{
	using namespace constants;

	double ratio(double x, double y)
	{return y==0?GARBAGE_VALUE:x/y;}
	double diff(double x, double y)
	{return x-y;}
	bool pedestal_quality(double x, double ref, double thr) 
	{return std::fabs(x-ref)/ref<=thr ? true : false;}
	bool pedestal_quality2(double x, double ref, double thr) 
	{return std::fabs(x-ref)<=thr ? true : false;}
	bool led_quality(double x, double ref, double thr) 
	{return std::fabs(x-ref)/ref<=thr?true:false;}
	bool laser_quality(double x, double ref, double thr) 
	{return std::fabs(x-ref)/ref<=thr?true:false;}

	/* virtual */ void ContainerXXX::initialize(hashfunctions::HashType ht,
		CompactUsageType ut, int debug)
	{
		_hashmap.initialize(ht);
		_usetype = ut;
		_logger.set("XXX", debug);
	}

	/* virtual */ void ContainerXXX::book(HcalElectronicsMap const* emap)
	{
		if (_hashmap.isDHash())
		{
			std::vector<HcalGenericDetId> dids = emap->allPrecisionId();
			for (std::vector<HcalGenericDetId>::const_iterator it=
				dids.begin(); it!=dids.end(); ++it)
			{
				if (!it->isHcalDetId())
					continue;

				uint32_t hash = it->rawId();
				HcalDetId did = HcalDetId(hash);
				_logger.debug(_hashmap.getName(did));
				CompactMap::iterator mit = _cmap.find(hash);
				if (mit!=_cmap.end())
					continue;

				_cmap.insert(
					std::make_pair(hash, Compact()));
			}
		}
		else if (_hashmap.isEHash())
		{
			std::vector<HcalElectronicsId> eids = 
				emap->allElectronicsIdPrecision();
			for (std::vector<HcalElectronicsId>::const_iterator it=
				eids.begin(); it!=eids.end(); ++it)
			{
				uint32_t hash = it->rawId();
				HcalElectronicsId eid = HcalElectronicsId(hash);
				_logger.debug(_hashmap.getName(eid));
				CompactMap::iterator mit = _cmap.find(hash);
				if (mit!=_cmap.end())
					continue;

				_cmap.insert(
					std::make_pair(hash, Compact()));
			}
		}
		else if (_hashmap.isTHash())
		{
			std::vector<HcalTrigTowerDetId> tids = emap->allTriggerId();
			for (std::vector<HcalTrigTowerDetId>::const_iterator it=
				tids.begin(); it!=tids.end(); ++it)
			{
				uint32_t hash = it->rawId();
				HcalTrigTowerDetId tid = HcalTrigTowerDetId(hash);
				_logger.debug(_hashmap.getName(tid));
				CompactMap::iterator mit = _cmap.find(hash);
				if (mit!=_cmap.end())
					continue;

				_cmap.insert(
					std::make_pair(hash, Compact()));
			}
		}
	}

	/* virtual */ void ContainerXXX::fill(HcalDetId const& did, double x)
	{
		if (_usetype!=fHistogram)
			return;

		Compact &c = _cmap[_hashmap.getHash(did)];
		c._x1 += x;
		c._x2 += x*x;
		c._n++;
	}

	/* virtual */ void ContainerXXX::fill(HcalElectronicsId const& did, 
		double x)
	{
		if (_usetype!=fHistogram)
			return;

		Compact &c = _cmap[_hashmap.getHash(did)];
		c._x1 += x;
		c._x2 += x*x;
		c._n++;
	}

	/* virtual */ void ContainerXXX::fill(HcalTrigTowerDetId const& did, 
		double x)
	{
		if (_usetype!=fHistogram)
			return;

		Compact &c = _cmap[_hashmap.getHash(did)];
		c._x1 += x;
		c._x2 += x*x;
		c._n++;
	}

	/* virtual */ void ContainerXXX::set(HcalDetId const& did, double x,
		double y)
	{
		if (_usetype==fHistogram)
			return;

		Compact &c = _cmap[_hashmap.getHash(did)];
		c._x1 = x;
		c._x2 = y;
		c._n=1;
	}

	/* virtual */ void ContainerXXX::set(HcalElectronicsId const& did, 
		double x, double y)
	{
		if (_usetype==fHistogram)
			return;

		Compact &c = _cmap[_hashmap.getHash(did)];
		c._x1 = x;
		c._x2 = y;
		c._n=1;
	}

	/* virtual */ void ContainerXXX::set(HcalTrigTowerDetId const& did, 
		double x, double y)
	{
		if (_usetype==fHistogram)
			return;

		Compact &c = _cmap[_hashmap.getHash(did)];
		c._x1 = x;
		c._x2 = y;
		c._n=1;
	}

	/* virtual */ uint32_t ContainerXXX::getEntries(HcalDetId const& did)
	{
		return _cmap[_hashmap.getHash(did)]._n;
	}

	/* virtual */ uint32_t ContainerXXX::getEntries(HcalElectronicsId const& 
		did)
	{
		return _cmap[_hashmap.getHash(did)]._n;
	}

	/* virtual */ uint32_t ContainerXXX::getEntries(HcalTrigTowerDetId const& 
		did)
	{
		return _cmap[_hashmap.getHash(did)]._n;
	}

	/* virtual */ uint32_t ContainerXXX::size()
	{
		return (uint32_t)(_cmap.size());
	}

	/* virtual */ void ContainerXXX::dump(Container1D* c, bool q)
	{
		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;
			if (x._n<=0)
				continue;

			std::pair<double, double> xpair = x.getValues(_usetype);
			q ? c->fill(hash, xpair.first) : 
				c->fill(hash, xpair.second);
		}
	}

	/* virtual */ 
	/*void ContainerXXX::dump(HcalElectronicsMap const* emap,
		Container1D* c, bool q)
	{
		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;
			if (x._n<=0)
				continue;

			uint32_t althash = 0;
			if (_hashmap.isDHash())	// then look up electronics
				hashToUse = emap->lookup(HcalDetId(hash)).rawId();
			else if (_hashmap.isTHash())
				hashToUse = emap->lookup(HcalTrigTowerDetId(hash)).rawId();
			else if (_hashmap.isEHash())
				hashToUse = emap->lookup(HcalElectronicsId(hash)).rawId();
			std::pair<double, double> xpair = x.getValues(_usetype);
			q ? c->fill(hashToUse, xpair.first) : 
				c->fill(hashToUse, xpair.second);
		}
	}*/

	/* virtual */ 
	/*void ContainerXXX::dump(HcalElectronicsMap const* emap,
		ContainerSingle1D* c, bool q)
	{
		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;
			if (x._n<=0)
				continue;

			uint32_t althash = 0;
			std::pair<double, double> xpair = x.getValues(_usetype);
			if (_hashmap.isDHash())	// then look up electronics
				c->fill(emap->lookup(HcalDetId(hash)), 
					q?xpair.first:xpair.second);
			else if (_hashmap.isTHash())
				c->fill(emap->lookupTrigger(HcalTrigTowerDetId(hash)), 
					q?xpair.first:xpair.second);
			else if (_hashmap.isEHash())
				c->fill(emap->lookupTrigger(HcalElectronicsId(hash)), 
					q?xpair.first:xpair.second);
		}
	}*/

	/* virtual */ 
	/*void ContainerXXX::dump(HcalElectronicsMap const* emap,
		ContainerSingle2D* c, bool q)
	{
		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;
			if (x._n<=0)
				continue;

			uint32_t althash = 0;
			std::pair<double, double> xpair = x.getValues(_usetype);
			if (_hashmap.isDHash())	// then look up electronics
			{
				HcalElectronicsId id=emap->lookup(HcalDetId(hash));
				if (id.isVMEid())
				c->fill(emap->lookup(HcalDetId(hash)), 
					q?xpair.first:xpair.second);
			}
			else if (_hashmap.isTHash())	// look Trigger Electronics
				c->fill(emap->lookupTrigger(HcalTrigTowerDetId(hash)), 
					q?xpair.first:xpair.second);
			else if (_hashmap.isEHash())	// look HcalDetId
				c->fill(emap->lookup(HcalElectronicsId(hash)), 
					q?xpair.first:xpair.second);
		}
	}*/

	/* virtual */ void ContainerXXX::dump(std::vector<Container1D*> const &vc, 
		bool q)
	{
		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;
			if (x._n<=0)
				continue;

			std::pair<double, double> xpair = x.getValues(_usetype);
			for (std::vector<Container1D*>::const_iterator it=vc.begin();
				it!=vc.end(); ++it)
				q ? (*it)->fill(hash, xpair.first) : 
					(*it)->fill(hash, xpair.second);
		}
	}

	/* virtual */ void ContainerXXX::dump(std::vector<Container1D*> const &vc, 
		Container1D* cmissing, bool q)
	{
		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;
			if (x._n<=0)
			{
				cmissing->fill(hash);
				continue;
			}

			std::pair<double, double> xpair = x.getValues(_usetype);
			for (std::vector<Container1D*>::const_iterator it=vc.begin();
				it!=vc.end(); ++it)
				q ? (*it)->fill(hash, xpair.first) : 
					(*it)->fill(hash, xpair.second);
		}
	}

	/* virtual */ void ContainerXXX::dump(Container1D* c, 
		Container1D *cmissing, bool q)
	{
		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;
			if (x._n<=0)
			{
				cmissing->fill(hash);
				continue;
			}

			std::pair<double, double> xpair = x.getValues(_usetype);
			q ? c->fill(hash, xpair.first) : 
				c->fill(hash, xpair.second);
		}
	}

	/* virtual */ void ContainerXXX::print()
	{
		std::cout << "Container by " << _hashmap.getHashTypeName() << std::endl;
		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			if (_hashmap.isDHash())
				std::cout << HcalDetId(p.first) << p.second << std::endl;
			else if (_hashmap.isEHash())
				std::cout << HcalElectronicsId(p.first) << p.second 
					<< std::endl;
			else if (_hashmap.isTHash())
				std::cout << HcalTrigTowerDetId(p.first) << p.second 
					<< std::endl;
		}
	}

	/* virtual */ void ContainerXXX::reset()
	{
		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			p.second.reset();
		}
	}

	/* virtual */ void ContainerXXX::load(Container1D* cont)
	{
		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;

			double value=0;
			if (_hashmap.isDHash())
				value = cont->getBinContent(HcalDetId(hash));
			else if (_hashmap.isEHash())
				value = cont->getBinContent(HcalElectronicsId(hash));
			else if (_hashmap.isTHash())
				value = cont->getBinContent(HcalTrigTowerDetId(hash));

			x._x1 += value;
			x._x2 += value*value;
			x._n++;
		}	
	}

	/* virtual */ void ContainerXXX::load(Container1D* c1, Container1D* c2)
	{
		//	only for value use cases
		if (_usetype==fHistogram)
			return;

		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;

			double m=0; double rms=0;
			if (_hashmap.isDHash())
			{
				m = c1->getBinContent(HcalDetId(hash));
				rms = c2->getBinContent(HcalDetId(hash));
			}
			else if (_hashmap.isEHash())
			{
				m = c1->getBinContent(HcalElectronicsId(hash));
				rms = c2->getBinContent(HcalElectronicsId(hash));
			}
			else if (_hashmap.isTHash())
			{
				m = c1->getBinContent(HcalTrigTowerDetId(hash));
				rms = c2->getBinContent(HcalTrigTowerDetId(hash));
			}

			x._x1 = m;
			x._x2 = rms;
			x._n=1;
		}	
	}

	/* virtual */ void ContainerXXX::loadPedestals(
		edm::ESHandle<HcalDbService> const& dbs)
	{
		if (!_hashmap.isDHash() || _usetype==fHistogram)
			return;

		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;
			HcalPedestal const* peds = dbs->getPedestal(
				HcalGenericDetId(hash));

			float const* means = peds->getValues();
			float const* widths = peds->getWidths();

			double msum=0; double rsum=0;
			for (uint32_t i=0; i<4; i++)
			{msum+=means[i];rsum+=widths[i];}
			msum/=4;rsum/=4; //	4 is the #caps - fixed even for QIE10...
			x._x1 = msum;
			x._x2 = rsum;
			x._n = 1;
		}
	}

	/* virtual */ void ContainerXXX::compare(ContainerXXX const& ctarget, 
		Container1D* cdump, Container1D* cnonpres, Container1D *cbad,
		comparison_function cfunc, QualityWrapper qwrap, bool q)
	{
		//	cannot compare sets with different hashtypes
		if (_hashmap.getHashType()!=ctarget._hashmap.getHashType())
			return;

		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;
			CompactMap::const_iterator it = ctarget._cmap.find(hash);

			//	skip if reference doesn't have that channel...
			if (x._n<=0)
				continue;

			//	if this channel is absent or didn't have entries - fill
			if (it==ctarget._cmap.end() || it->second._n==0)
			{
				cnonpres->fill(hash);
				continue;
			}

			//	both channels exist and do have entries
			Compact y = it->second;
			std::pair<double, double> p1 = x.getValues(_usetype);
			std::pair<double, double> p2 = y.getValues(ctarget._usetype);

			//	get the comparison value
			double xxx = q?p2.first:p2.second;
			double ref = q?p1.first:p1.second;
			double ccc = (*cfunc)(xxx, ref);

			//	check quality and fill the container for comparison
			if (!qwrap.quality(xxx,ref) || ccc==GARBAGE_VALUE)
				cbad->fill(hash);
			cdump->fill(hash, ccc);
		}
	}

	/* virtual */ void ContainerXXX::compare(ContainerXXX const& ctarget, 
		std::vector<Container1D*> const &vcdump, 
		Container1D* cnonpres, Container1D *cbad,
		comparison_function cfunc, QualityWrapper qwrap, bool q)
	{
		//	cannot compare sets with different hashtypes
		if (_hashmap.getHashType()!=ctarget._hashmap.getHashType())
			return;

		BOOST_FOREACH(CompactMap::value_type &p, _cmap)
		{
			Compact &x = p.second;
			uint32_t hash = p.first;
			CompactMap::const_iterator it = ctarget._cmap.find(hash);

			//	skip if reference doesn't have that channel...
			if (x._n<=0)
				continue;

			//	if this channel is absent or didn't have entries - fill
			if (it==ctarget._cmap.end() || it->second._n==0)
			{
				cnonpres->fill(hash);
				continue;
			}

			//	both channels exist and do have entries
			Compact y = it->second;
			std::pair<double, double> p1 = x.getValues(_usetype);
			std::pair<double, double> p2 = y.getValues(ctarget._usetype);

			//	get the comparison value
			double xxx = q?p2.first:p2.second;
			double ref = q?p1.first:p1.second;
			double ccc = (*cfunc)(xxx, ref);

			//	check quality and fill the container for comparison
			if (!qwrap.quality(xxx,ref) || ccc==GARBAGE_VALUE)
				cbad->fill(hash);
			for (std::vector<Container1D*>::const_iterator it=vcdump.begin();
				it!=vcdump.end(); ++it)
				(*it)->fill(hash, ccc);
		}
	}
}


