#ifndef RecoLocalCalo_HcalRecAlgos_PulseShapeFunctor_gpu_h
#define RecoLocalCalo_HcalRecAlgos_PulseShapeFunctor_gpu_h

namespace hcal { namespace mahi {

namespace HcalConst{

   constexpr int maxSamples = 10;
   constexpr int maxPSshapeBin = 256;
   constexpr int nsPerBX = 25;
   constexpr float iniTimeShift = 92.5f;
   constexpr double invertnsPerBx = 0.04;
   constexpr int shiftTS = 4;

}

namespace FitterFuncs {
  
   class PulseShapeFunctor {
   public:
     __device__
     PulseShapeFunctor(float const* pulse, bool iPedestalConstraint, 
                       bool iTimeConstraint,bool iAddPulseJitter,
		               double iPulseJitter,double iTimeMean,double iPedMean,
		               unsigned int nSamplesToFit);
    
     __device__
     double EvalPulse(const double *pars, const unsigned nPar);
     
     __device__
     void setDefaultcntNANinfit(){ cntNANinfit =0; }
     __device__
     int getcntNANinfit(){ return cntNANinfit; }
     
     __device__
     void setpsFitx(double *x )
     { for(int i=0; i<HcalConst::maxSamples; ++i) psFit_x[i] = x[i]; }
     __device__
     void setpsFity(double *y )
     { for(int i=0; i<HcalConst::maxSamples; ++i) psFit_y[i] = y[i]; }
     __device__
     void setpsFiterry (double *erry  )
     { for(int i=0; i<HcalConst::maxSamples; ++i) psFit_erry  [i] = erry [i]; }
     __device__
     void setpsFiterry2(double *erry2 )
     { for(int i=0; i<HcalConst::maxSamples; ++i) psFit_erry2 [i] = erry2[i]; }
     __device__
     void setpsFitslew (double *slew  )
     { for(int i=0; i<HcalConst::maxSamples; ++i) {psFit_slew [i] = slew [i]; } }
     __device__
     double getSiPMDarkCurrent(double darkCurrent, double fcByPE, double lambda);
     __device__
     void setinvertpedSig2(double x) { invertpedSig2_ = x; }
     __device__
     void setinverttimeSig2(double x) { inverttimeSig2_ = x; }

     __device__
     double singlePulseShapeFunc( const double *x );
     __device__
     double doublePulseShapeFunc( const double *x );
     __device__
     double triplePulseShapeFunc( const double *x );

     __device__
     void getPulseShape(double fillPulseShape[HcalConst::maxSamples]) { 
       for (unsigned int i=0; i<HcalConst::maxSamples; i++)
         fillPulseShape[i] = pulse_shape_[i];
     }
     
   private:
     float pulse_hist[HcalConst::maxPSshapeBin];
     
     int cntNANinfit;
     
     /*
     std::vector<float> acc25nsVec, diff25nsItvlVec;
     std::vector<float> accVarLenIdxZEROVec, diffVarItvlIdxZEROVec;
     std::vector<float> accVarLenIdxMinusOneVec, diffVarItvlIdxMinusOneVec;
     */
     float acc25nsVec[HcalConst::maxPSshapeBin], 
           diff25nsItvlVec[HcalConst::maxPSshapeBin];
     float accVarLenIdxZEROVec[HcalConst::nsPerBX],
           diffVarItvlIdxZEROVec[HcalConst::nsPerBX];
     float accVarLenIdxMinusOneVec[HcalConst::nsPerBX],
           diffVarItvlIdxMinusOneVec[HcalConst::nsPerBX];

     __device__
     void funcShape(double ntmpbin[HcalConst::maxSamples], 
                    const double pulseTime, const double pulseHeight,
                    const double slew);
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
     double pulse_shape_[HcalConst::maxSamples] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
     double pulse_shape_sum_[HcalConst::maxSamples]
         = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//     std::array<double,HcalConst::maxSamples> pulse_shape_;
//     std::array<double,HcalConst::maxSamples> pulse_shape_sum_;

   };
   
}

}}

#endif // PulseShapeFunctor_h
