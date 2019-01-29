#ifndef ECALDETID_ECALTRIGTOWERDETID_H
#define ECALDETID_ECALTRIGTOWERDETID_H

#include <iosfwd>
#include <cassert>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "FWCore/Utilities/interface/Exception.h"

/** \class EcalTrigTowerDetId
    
   DetId for an Ecal Trigger tower

*/

// a workaround for std::abs not being a constexpr function
namespace ecal {

template<typename T>
constexpr T abs(T const& value) {
    return ::std::max(value, -value);
}

}

class EcalTrigTowerDetId : public DetId {
 public:
  /** Constructor of a null id */
  constexpr EcalTrigTowerDetId() {
  }
  /** Constructor from a raw value */
  constexpr EcalTrigTowerDetId(uint32_t rawid) : DetId(rawid) {
  }
  /** \brief Constructor from signed ieta, iphi
   */
  constexpr EcalTrigTowerDetId(int zside, EcalSubdetector subDet, int i, int j, int mode=SUBDETIJMODE)
    : DetId(Ecal,EcalTriggerTower) 
  {
    int tower_i=0;
    int tower_j=0;

    if (mode == SUBDETIJMODE)
    {
      tower_i=i;
      tower_j=j;
    }
    else if (mode == SUBDETDCCTTMODE)
    {
      throw cms::Exception("InvalidDetId") << "EcalTriggerTowerDetId:  Cannot create object. SUBDETDCCTTMODE not yet implemented.";   
    }
    else
      throw cms::Exception("InvalidDetId") << "EcalTriggerTowerDetId:  Cannot create object.  Unknown mode for (int, EcalSubdetector, int, int) constructor.";
  
    if (tower_i > MAX_I || tower_i < MIN_I  || tower_j > MAX_J || tower_j < MIN_J)
      throw cms::Exception("InvalidDetId") << "EcalTriggerTowerDetId:  Cannot create object.  Indexes out of bounds.";
  
    id_|= ((zside>0)?(0x8000):(0x0)) | ((subDet == EcalBarrel)?(0x4000):(0x0)) | (tower_i<<7) | (tower_j & 0x7F);

  }
  
  /** Constructor from a generic cell id */
  constexpr EcalTrigTowerDetId(const DetId& gen) 
  {
    if (!gen.null() && ( gen.det()!=Ecal || gen.subdetId()!=EcalTriggerTower )) {
      throw cms::Exception("InvalidDetId");  }
    id_=gen.rawId();
  }
  /** Assignment from a generic cell id */
  constexpr EcalTrigTowerDetId& operator=(const DetId& gen) {
    if (!gen.null() && ( gen.det()!=Ecal || gen.subdetId()!=EcalTriggerTower )) {
      throw cms::Exception("InvalidDetId");
    }
    id_=gen.rawId();
    return *this;
  }
  

  /// get the z-side of the tower (1/-1)
  constexpr int zside() const { return (id_&0x8000)?(1):(-1); }

  /// get the subDetector associated to the Trigger Tower
  constexpr EcalSubdetector subDet() const { return (id_&0x4000) ? EcalBarrel:EcalEndcap; }

  /// get the absolute value of the tower ieta 
  constexpr int ietaAbs() const 
   { 
     /*       if ( subDet() == EcalBarrel) */
     return (id_>>7)&0x7f; 
     /*       else */
     /* 	throw(std::runtime_error("EcalTrigTowerDetId: ietaAbs not applicable for this subDetector.")); */
   }  

  /// get the tower ieta 
  constexpr int ieta() const 
    { 
      /*       if ( subDet() == EcalBarrel) */
      return zside()*ietaAbs(); 
      /*       else */
      /* 	throw(std::runtime_error("EcalTrigTowerDetId: ieta not applicable for this subDetector.")); */
    } 
  
  /// get the tower iphi 
  constexpr int iphi() const 
    { 
      /*       if ( subDet() == EcalBarrel) */
      return id_&0x7F; 
      /*       else */
      /* 	throw(std::runtime_error("EcalTrigTowerDetId: iphi not applicable for this subDetector.")); */
      
    } 

  constexpr int iquadrant() const
  {
    if ( subDet() == EcalEndcap )
      return int((iphi()-1)/kEETowersInPhiPerQuadrant)+1;
    else
      throw cms::Exception("MethodNotApplicable") << "EcalTriggerTowerDetId: iquadrant not applicable";
  }  


   /// get the tower ix (Endcap case) */
   constexpr int ix() const  
     {  
       if ( subDet() == EcalEndcap) 
 	return (id_>>7)&0x7f;  
       else 
 	throw(std::runtime_error("EcalTrigTowerDetId: ix not applicable for this subDetector."));
     }  
  
   /// get the tower iy (Endcap case) */
   constexpr int iy() const
     { 
      if ( subDet() == EcalEndcap)
        return id_&0x7F; 
      else 
        throw(std::runtime_error("EcalTrigTowerDetId: ix not applicable for this subDetector.")); 
     }  
  

  /// get a compact index for arrays [TODO: NEEDS WORK]
  constexpr int 
  hashedIndex() const 
  {
    const unsigned int iea ( ietaAbs() ) ;
    const unsigned int iph ( iphi()    ) ;
    return ( subDet() == EcalBarrel  ? 
	    ( iDCC() - 1 )*kEBTowersPerSM + iTT() - 1 :
	    kEBTotalTowers + ( ( zside() + 1 )/2 )*kEETowersPerEndcap +
	    ( ( iea < 27 ? iea : 27 ) - kEEOuterEta )*kEETowersInPhiPerEndcap + 
	    ( iea < 27 ? iph : // for iphi=27,28 only half TT present, odd for EE-, even EE+
	      ( iea - 27 )*kEETowersInPhiPerEndcap/2 + ( iph + 1 )/2 ) - 1 ) ;
  }

  constexpr uint32_t denseIndex() const { return hashedIndex() ; }

  static constexpr bool validDenseIndex( uint32_t din ) { return ( din < kSizeForDenseIndexing ) ; }

  static constexpr EcalTrigTowerDetId 
  detIdFromDenseIndex( uint32_t di ) 
  {
    const EcalSubdetector sd ( di < kEBTotalTowers ? EcalBarrel : EcalEndcap ) ;
    const int iz ( di < kEBTotalTowers ? 
		  ( di < kEBHalfTowers ?  1 : -1 ) :
		  ( di - kEBTotalTowers < kEETowersPerEndcap ? -1 : 1 ) ) ;
    int i {0};
    int j {0};
    if( di < kEBTotalTowers ) // barrel
    {
      const unsigned int itt ( di%kEBTowersPerSM ) ;
      const unsigned int idc ( di/kEBTowersPerSM ) ;
      j = (idc%18)*kEBTowersInPhi + 
	 ( (1+iz)/2 )*kEBTowersInPhi - 
	 iz*(itt%kEBTowersInPhi)  + 1 - (1+iz)/2 - 2 ;
      if( j < 1 ) j += 72 ;
      i = 1 + itt/kEBTowersInPhi ;
    }
    else
    {
      const int eonly ( ( di - kEBTotalTowers )%kEETowersPerEndcap ) ;
      i = kEEOuterEta + eonly/kEETowersInPhiPerEndcap ;
      j = 1 + eonly%kEETowersInPhiPerEndcap ;
      if( 27 == i ) // last two rings have half of normal phi elementes
      {
	 if( j > kEETowersInPhiPerEndcap/2 )
	 {
	    ++i ; 
	    j -= kEETowersInPhiPerEndcap/2 ;
	 }
	 j = 2*j ;
	 if( 0 < iz ) --j ;
      }
    }
    assert( validDetId( iz, sd, i, j ) ) ;
    return EcalTrigTowerDetId( iz, sd, i, j ) ;
  }

  /// check if a valid index combination
  // TODO: make sure ecal::abs does not introduce any penalty
  static constexpr bool
  validDetId( int iz, EcalSubdetector sd , int i, int j )
  {
    return ( 1 == ecal::abs( iz )               &&
	    0 <  i                       &&
	    0 <  j                       &&
	    kEETowersInPhiPerEndcap >= j &&
	    ( ( EcalBarrel     == sd &&
		kEBTowersInEta >=  i     ) ||
	      ( EcalEndcap     == sd &&
		kEEOuterEta    <=  i &&
		kEEInnerEta    >=  i &&
		( 27 > i ||
		  (  ( 0 >  iz  &&
		       0 == j%2    ) ||
		     ( 0 <  iz  &&
		       1 == j%2         ) ) ) ) ) ) ;
	    
  }

  /// get the ECAL DCC id - in the  barrrel ism == iDCC
  //New SM numbering scheme. Avoids discontinuity in phi crossing \eta=0  
  constexpr int iDCC() const 
  {
    if ( subDet() == EcalBarrel )
    {
      //Correction since iphi is uniformized with HB convention 
      int iphi_simple = iphi() + 2 ;
      if (iphi_simple > 72 ) iphi_simple = iphi_simple % 72;
      int id = ( iphi_simple - 1 ) / kEBTowersInPhi + 1;
      if ( zside() < 0 ) id += 18;
      return id;
    }
    else
      throw cms::Exception("MethodNotImplemented") << "EcalTriggerTowerDetId: iDCC not yet implemented";
  }

  /// sequential index within one DCC
  constexpr int iTT() const 
  {
    if ( subDet() == EcalBarrel )
    {
      int ie = ietaAbs() -1;
      int ip{0};
      int iphi_simple = iphi() + 2 ;
      if (iphi_simple > 72 )  iphi_simple = iphi_simple % 72;
      if (zside() < 0) {
	ip = (( iphi_simple -1 ) % kEBTowersInPhi ) + 1;
      } else {
	ip = kEBTowersInPhi - ((iphi_simple -1 ) % kEBTowersInPhi );
      }
      
      return (ie * kEBTowersInPhi) + ip;
    }
    else
      throw cms::Exception("MethodNotImplemented") << "EcalTriggerTowerDetId: iTT not yet implemented";
  }



  static const int MIN_I = 1;
  static const int MIN_J = 1;
  static const int MAX_I = 127;
  static const int MAX_J = 127;

  static const int kEBTowersInPhi = 4; // per SM (in the Barrel)
  static const int kEBTowersPerSM = 68; // per SM (in the Barrel)
  static const int kEBTowersInEta = 17; // per SM (in the Barrel)
//  static const int kEETowersInEta = 11; // Endcap
  static const int kEETowersInPhiPerQuadrant = 18; // per Quadrant (in the Endcap)

  // function modes for (int, int) constructor
  static const int SUBDETIJMODE = 0;
  static const int SUBDETDCCTTMODE = 1;

      enum { kEETowersInPhiPerEndcap = 4* kEETowersInPhiPerQuadrant ,
	     kEEOuterEta             = 18 ,
	     kEEInnerEta             = 28 ,
	     kEETowersInEta          = ( kEEInnerEta - kEEOuterEta + 1 ) ,
	     kEBHalfTowers           = kEBTowersPerSM*18 ,
	     kEBTotalTowers          = kEBHalfTowers*2 ,
	     kEETowersPerEndcap      = kEETowersInEta*kEETowersInPhiPerEndcap - 72,
	     kEETotalTowers          = kEETowersPerEndcap*2,
	     kSizeForDenseIndexing   = kEBTotalTowers + kEETotalTowers } ;
};

std::ostream& operator<<(std::ostream&,const EcalTrigTowerDetId& id);

#endif
