#ifndef RecoLocalCalo_HcalRecAlgos_MahiAux_h
#define RecoLocalCalo_HcalRecAlgos_MahiAux_h

namespace HcalConst{

   constexpr int maxSamples = 10;
   constexpr int maxPSshapeBin = 256;
   constexpr int nsPerBX = 25;
   constexpr float iniTimeShift = 92.5f;
   constexpr double invertnsPerBx = 0.04;
   constexpr int shiftTS = 4;

}

namespace FitterFuncs{
  
   class MahiFunctor {
      public:
     // assume pulse has the above-assumed number of bins -> maxPSshapeBikn
     __device__ MahiFunctor();
     __device__ void assign(float const * pulse,bool iPedestalConstraint, bool iTimeConstraint,bool iAddPulseJitter,
		       double iPulseJitter,double iTimeMean,double iPedMean,
		       unsigned int nSamplesToFit);
     
     __device__ double EvalPulse(const double *pars, const unsigned nPar);
     
     __device__ double singlePulseShapeFunc( const double *x );
     __device__ double doublePulseShapeFunc( const double *x );
     __device__ double triplePulseShapeFunc( const double *x );

     __device__ void getPulseShape(double *fillPulseShape) { 
         for (int i=0; i<HcalConst::maxSamples; i++)
             fillPulseShape[i] = pulse_shape_[i];
     }
     
   private:
     float pulse_hist[HcalConst::maxPSshapeBin];
     
     int cntNANinfit;
     float acc25nsVec[HcalConst::maxPSshapeBin], 
           diff25nsItvlVec[HcalConst::maxPSshapeBin];
     float accVarLenIdxZEROVec[HcalConst::nsPerBX], 
           diffVarItvlIdxZEROVec[HcalConst::nsPerBX];
     float accVarLenIdxMinusOneVec[HcalConst::nsPerBX], 
           diffVarItvlIdxMinusOneVec[HcalConst::nsPerBX];

     //
     __device__
     void funcShape(double *ntmpbin, const double pulseTime, 
                    const double pulseHeight,const double slew);
     double psFit_x[HcalConst::maxSamples], psFit_y[HcalConst::maxSamples], psFit_erry[HcalConst::maxSamples], psFit_erry2[HcalConst::maxSamples], psFit_slew[HcalConst::maxSamples];
     
     unsigned nSamplesToFit_;
     bool pedestalConstraint_;
     bool timeConstraint_;
     bool addPulseJitter_;
     bool unConstrainedFit_;
     double pulseJitter_;
     double timeMean_;
     double timeSig_;
     double pedMean_;
     double timeShift_;

     double inverttimeSig2_;
     double invertpedSig2_;
     double pulse_shape_[HcalConst::maxSamples];
     double pulse_shape_sum_[HcalConst::maxSamples];

   };
   
}

#endif
