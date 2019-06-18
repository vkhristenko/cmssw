/*
 *
 */
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalCorrectionFunctions.h"

const int num_bins_hbheho = 61;
const float wpksamp0_hbheho = 0.5;

#ifdef __CUDA_ARCH__
__constant__
#endif
const float actual_ns_hbheho[num_bins_hbheho] = {
-5.44000, // 0.500, 0.000-0.017
-4.84250, // 0.517, 0.017-0.033
-4.26500, // 0.533, 0.033-0.050
-3.71000, // 0.550, 0.050-0.067
-3.18000, // 0.567, 0.067-0.083
-2.66250, // 0.583, 0.083-0.100
-2.17250, // 0.600, 0.100-0.117
-1.69000, // 0.617, 0.117-0.133
-1.23000, // 0.633, 0.133-0.150
-0.78000, // 0.650, 0.150-0.167
-0.34250, // 0.667, 0.167-0.183
 0.08250, // 0.683, 0.183-0.200
 0.50250, // 0.700, 0.200-0.217
 0.90500, // 0.717, 0.217-0.233
 1.30500, // 0.733, 0.233-0.250
 1.69500, // 0.750, 0.250-0.267
 2.07750, // 0.767, 0.267-0.283
 2.45750, // 0.783, 0.283-0.300
 2.82500, // 0.800, 0.300-0.317
 3.19250, // 0.817, 0.317-0.333
 3.55750, // 0.833, 0.333-0.350
 3.91750, // 0.850, 0.350-0.367
 4.27500, // 0.867, 0.367-0.383
 4.63000, // 0.883, 0.383-0.400
 4.98500, // 0.900, 0.400-0.417
 5.33750, // 0.917, 0.417-0.433
 5.69500, // 0.933, 0.433-0.450
 6.05000, // 0.950, 0.450-0.467
 6.40500, // 0.967, 0.467-0.483
 6.77000, // 0.983, 0.483-0.500
 7.13500, // 1.000, 0.500-0.517
 7.50000, // 1.017, 0.517-0.533
 7.88250, // 1.033, 0.533-0.550
 8.26500, // 1.050, 0.550-0.567
 8.66000, // 1.067, 0.567-0.583
 9.07000, // 1.083, 0.583-0.600
 9.48250, // 1.100, 0.600-0.617
 9.92750, // 1.117, 0.617-0.633
10.37750, // 1.133, 0.633-0.650
10.87500, // 1.150, 0.650-0.667
11.38000, // 1.167, 0.667-0.683
11.95250, // 1.183, 0.683-0.700
12.55000, // 1.200, 0.700-0.717
13.22750, // 1.217, 0.717-0.733
13.98500, // 1.233, 0.733-0.750
14.81500, // 1.250, 0.750-0.767
15.71500, // 1.267, 0.767-0.783
16.63750, // 1.283, 0.783-0.800
17.53750, // 1.300, 0.800-0.817
18.38500, // 1.317, 0.817-0.833
19.16500, // 1.333, 0.833-0.850
19.89750, // 1.350, 0.850-0.867
20.59250, // 1.367, 0.867-0.883
21.24250, // 1.383, 0.883-0.900
21.85250, // 1.400, 0.900-0.917
22.44500, // 1.417, 0.917-0.933
22.99500, // 1.433, 0.933-0.950
23.53250, // 1.450, 0.950-0.967
24.03750, // 1.467, 0.967-0.983
24.53250, // 1.483, 0.983-1.000
25.00000  // 1.500, 1.000-1.017 - keep for interpolation
};

__host__ __device__ 
float timeshift_ns_hbheho(float wpksamp) {
    float flx = (num_bins_hbheho-1)*(wpksamp - wpksamp0_hbheho);
    int index = (int)flx;

    if      (index <    0)               return actual_ns_hbheho[0];
    else if (index >= num_bins_hbheho-1) return actual_ns_hbheho[num_bins_hbheho-1];

    // else interpolate:
    float y1 = actual_ns_hbheho[index];
    float y2 = actual_ns_hbheho[index+1];
    return y1 + (y2-y1)*(flx-(float)index);
}

/// Ugly hack to apply energy corrections to some HB- cells
__host__ __device__
float hbminus_special_ecorr(int ieta, int iphi, double energy, int runnum) {
    // return energy correction factor for HBM channels 
    // iphi=6 ieta=(-1,-15) and iphi=32 ieta=(-1,-7)
    // I.Vodopianov 28 Feb. 2011
    static const float low32[7]  = {0.741,0.721,0.730,0.698,0.708,0.751,0.861};
    static const float high32[7] = {0.973,0.925,0.900,0.897,0.950,0.935,1};
    static const float low6[15]  = {0.635,0.623,0.670,0.633,0.644,0.648,0.600,
                                0.570,0.595,0.554,0.505,0.513,0.515,0.561,0.579};
    static const float high6[15] = {0.875,0.937,0.942,0.900,0.922,0.925,0.901,
                                  0.850,0.852,0.818,0.731,0.717,0.782,0.853,0.778};

              
    double slope, mid, en; 
    double corr = 1.0;

    if (!(iphi==6 && ieta<0 && ieta>-16) && !(iphi==32 && ieta<0 && ieta>-8)) 
        return corr;

    int jeta = -ieta-1;
    double xeta = (double) ieta;
    if (energy > 0.) en=energy;
    else en = 0.; 

    if (iphi == 32) {
        slope = 0.2272;
        mid = 17.14 + 0.7147*xeta;
        if (en > 100.) corr = high32[jeta];
        else corr = low32[jeta]+(high32[jeta]-low32[jeta])/(1.0+exp(-(en-mid)*slope));
    } else if (iphi == 6 && runnum < 216091 ) { 
        slope = 0.1956;
        mid = 15.96 + 0.3075*xeta;
        if (en > 100.0) corr = high6[jeta];
        else corr = low6[jeta]+(high6[jeta]-low6[jeta])/(1.0+exp(-(en-mid)*slope));
    }

    //  std::cout << "HBHE cell:  ieta, iphi = " << ieta << "  " << iphi 
    //        << "  ->  energy = " << en << "   corr = " << corr << std::endl;

    return corr;
}
