#ifndef RecoLocalCalo_HcalRecAlgos_MahiFit_gpu_test4cpu_HH
#define RecoLocalCalo_HcalRecAlgos_MahiFit_gpu_test4cpu_HH

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/EigenMatrixTypes.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/PulseShapeFunctor_gpu_test4cpu.h"

namespace hcal { namespace mahi { namespace test {

struct MahiNnlsWorkspace {

  unsigned int nPulseTot;
  unsigned int tsSize;
  unsigned int tsOffset;
  unsigned int fullTSOffset;
  int bxOffset;
  int maxoffset;
  double dt;

  //holds active bunch crossings
  BXVector bxs;  

  //holds data samples
  SampleVector amplitudes;

  //holds inverse covariance matrix
  SampleMatrix invCovMat;

  //holds diagonal noise terms
  SampleVector noiseTerms;

  //holds flat pedestal uncertainty
  SampleMatrix pedConstraint;
  
  //holds full covariance matrix for a pulse shape 
  //varied in time
  FullSampleMatrix pulseCovArray[MaxPVSize];

  //holds full pulse shape template
  FullSampleVector pulseShapeArray[MaxPVSize];

  //holds full pulse shape derivatives
  FullSampleVector pulseDerivArray[MaxPVSize];

  //holders for calculating pulse shape & covariance matrices
  double pulseN[MaxSVSize];
  double pulseM[MaxSVSize];
  double pulseP[MaxSVSize];

  //holds matrix of pulse shape templates for each BX
  SamplePulseMatrix pulseMat;

  //holds matrix of pulse shape derivatives for each BX
  SamplePulseMatrix pulseDerivMat;

  //holds residual vector
  PulseVector residuals;

  //for FNNLS algorithm
  unsigned int nP;
  PulseVector ampVec;

  PulseVector ampvecpermtest;

  SamplePulseMatrix invcovp;
  PulseMatrix aTaMat; // A-transpose A (matrix)
  PulseVector aTbVec; // A-transpose b (vector)
  PulseVector updateWork; // w (vector)

  SampleDecompLLT covDecomp;
  PulseDecompLDLT pulseDecomp;

};

struct MahiDebugInfo {

  int   nSamples;
  int   soi;

  bool  use3;

  float inTimeConst;
  float inDarkCurrent;
  float inPedAvg;
  float inGain;
  
  float inNoiseADC[MaxSVSize];
  float inNoiseDC[MaxSVSize];
  float inNoisePhoto[MaxSVSize];
  float inPedestal[MaxSVSize];

  float totalUCNoise[MaxSVSize];

  float mahiEnergy;
  float chiSq;
  float arrivalTime;

  float pEnergy;
  float nEnergy;
  float pedEnergy;

  float count[MaxSVSize];
  float inputTS[MaxSVSize];
  int inputTDC[MaxSVSize];
  float itPulse[MaxSVSize];
  float pPulse[MaxSVSize];
  float nPulse[MaxSVSize];
  
};

class MahiFit
{
 public:
   
  MahiFit(float const*);

  
  void phase1Apply(const HBHEChannelInfo& channelData, 
		   float& reconstructedEnergy, 
		   float& reconstructedTime, 
		   bool& useTriple,
		   float& chi2) const;

  
  void phase1Debug(const HBHEChannelInfo& channelData,
		   MahiDebugInfo& mdi) const;

  
  void doFit(float correctedOutput[3], const int nbx) const;

  
  //void setPulseShapeTemplate  (const HcalPulseShapes::Shape& ps,const HcalTimeSlew * hcalTimeSlewDelay);
  void setPulseShapeTemplate(float const* pshape);
/*  
  void resetPulseShapeTemplate(const HcalPulseShapes::Shape& ps);
  */

  typedef BXVector::Index Index;
  float const* pshape_;
  // TODO: resolve this
//  const HcalTimeSlew* hcalTimeSlewDelay_=nullptr;

 private:

  
  double minimize() const;
  
  void onePulseMinimize() const;
  
  void updateCov() const;
  
  void updatePulseShape(double itQ, FullSampleVector &pulseShape, 
			FullSampleVector &pulseDeriv,
			FullSampleMatrix &pulseCov) const;

  
  double calculateArrivalTime() const;
  
  double calculateChiSq() const;
  
  void nnls() const;
  
  void resetWorkspace() const;

  
  void nnlsUnconstrainParameter(Index idxp) const;
  
  void nnlsConstrainParameter(Index minratioidx) const;

  
  void solveSubmatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned nP) const;

  mutable MahiNnlsWorkspace nnlsWork_;

  //hard coded in initializer
  const unsigned int fullTSSize_;
  const unsigned int fullTSofInterest_;

  static constexpr int pedestalBX_ = 100;

  // used to restrict returned time value to a 25 ns window centered 
  // on the nominal arrival time
  static constexpr float timeLimit_ = 12.5;

  // Python-configurables
  bool dynamicPed_ {true};
  float ts4Thresh_ {0.0};
  float chiSqSwitch_{15.0}; 

  bool applyTimeSlew_{true}; 
  /* HcalTimeSlew::BiasSetting */ int slewFlavor_{1};
  double tsDelay1GeV_{10};

  float meanTime_{0.0};
  float timeSigmaHPD_{5.0}; 
  float timeSigmaSiPM_{2.5};

  int activeBXs_[3] = {-1, 0, 1};

  int nMaxItersMin_{500}; 
  int nMaxItersNNLS_{500}; 

  float deltaChiSqThresh_{1e-3}; 
  float nnlsThresh_{1e-11}; 

  unsigned int bxSizeConf_{3};
  int bxOffsetConf_{1}; // -(min_element({-1, 0, 1}))

  //for pulse shapes
  int cntsetPulseShape_;

  mutable FitterFuncs::PulseShapeFunctor functor_;
};

}}}

#endif
