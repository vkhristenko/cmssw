#ifndef CALIBCALORIMETRY_HCALALGOS_HCALTIMESLEW_H
#define CALIBCALORIMETRY_HCALALGOS_HCALTIMESLEW_H 1

#ifndef __CUDA_ARCH__
#include <vector>
#endif

/** \class HcalTimeSlew
  * 
  * Provides pulse delay as a function of amplitude for three choices
  * of QIE bias setting.  The "Medium" setting is used in HB and HE,
  * while the "Slow" (and lower noise) setting is used in HO.  All
  * data taken from bench measurements of the QIE and plotted in
  * Physics TDR Vol 1.
  * 
  * Not to be used for HF at this time (unlikely to have much effect)
  *
  * \original author        J. Mans   - Minnesota
  * \upgraded 2017          C. Madrid - Baylor
  */
class HcalTimeSlew {
 public:
  class HcalTimeSlewM2Parameters{
  public:
    //M2 Parameters
    double tzero; 
    double slope; 
    double  tmax; 

#ifndef __CUDA_ARCH__
    HcalTimeSlewM2Parameters(double t0, double m, double tmaximum):tzero(t0), slope(m), tmax(tmaximum){}
#endif
  };

  class HcalTimeSlewM3Parameters{
  public:
    //M3 Parameters
    double cap;            
    double tspar0;         
    double tspar1;         
    double tspar2;         
    double tspar0_siPM;    
    double tspar1_siPM;    
    double tspar2_siPM;    
  
#ifndef __CUDA_ARCH__
  HcalTimeSlewM3Parameters(double capCon, double tspar0Con, double tspar1Con, double tspar2Con, double tspar0_siPMCon, double tspar1_siPMCon, double tspar2_siPMCon):cap(capCon), tspar0(tspar0Con), tspar1(tspar1Con), tspar2(tspar2Con), tspar0_siPM(tspar0_siPMCon), tspar1_siPM(tspar1_siPMCon), tspar2_siPM(tspar2_siPMCon){} 
#endif
  };

#ifndef __CUDA_ARCH__
  HcalTimeSlew() {};
  ~HcalTimeSlew() {}

  void addM2ParameterSet(double tzero, double slope, double tmax);
  void addM3ParameterSet(double cap, double tspar0, double tspar1, double tspar2, double tspar0_siPM, double tspar1_siPM, double tspar2_siPM);
#endif

  enum ParaSource { TestStand=0, Data=1, MC=2, HBHE=3 };
  enum BiasSetting { Slow=0, Medium=1, Fast=2 };
  /** \brief Returns the amount (ns) by which a pulse of the given
   number of fC will be delayed by the timeslew effect, for the
   specified bias setting. */
#ifndef __CUDA_ARCH__
  double delay(double fC, BiasSetting bias=Medium) const; 
  double delay(double fC, ParaSource source=HBHE, BiasSetting bias=Medium, bool isHPD=true) const;

  std::vector<HcalTimeSlewM2Parameters> const& get_m2_params() const
  { return parametersM2_; }
  std::vector<HcalTimeSlewM3Parameters> const& get_m3_params() const
  { return parametersM3_; }

private:
  std::vector<HcalTimeSlewM2Parameters> parametersM2_;
  std::vector<HcalTimeSlewM3Parameters> parametersM3_;
#endif
};

#endif
