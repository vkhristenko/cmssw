/* Sushil Dubey, Shashi Dugad, TIFR, July 2017
 *
 * File Name: RawToDigiGPU.cu
 * Description: It converts Raw data into Digi Format on GPU
 * then it converts adc -> electron and
 * applies the adc threshold to needed for clustering
 * Finaly the Output of RawToDigi data is given to pixelClusterizer
 *
**/

// System includes
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>

// CUDA runtime
#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>

#include "cudaCheck.h"
#include "EventInfoGPU.h"
#include "RawToDigiGPU.h"
#include "SiPixelFedCablingMapGPU.h"


context initDeviceMemory() {
  context c;

  // Number of words for all the feds
  constexpr uint32_t MAX_WORD_SIZE = MAX_FED * MAX_WORD  * sizeof(uint32_t);
  constexpr uint32_t FSIZE = 2*MAX_FED*sizeof(uint32_t)+sizeof(uint32_t);
  constexpr uint32_t MSIZE = NMODULE*sizeof(int)+sizeof(int);

  cudaCheck(cudaMalloc((void**) & c.word_d,        MAX_WORD_SIZE));
  cudaCheck(cudaMalloc((void**) & c.fedIndex_d,    FSIZE));
  cudaCheck(cudaMalloc((void**) & c.xx_d,          MAX_WORD_SIZE)); // to store the x and y coordinate
  cudaCheck(cudaMalloc((void**) & c.yy_d,          MAX_WORD_SIZE));
  cudaCheck(cudaMalloc((void**) & c.xx_adc,        MAX_WORD_SIZE)); // to store the x and y coordinate
  cudaCheck(cudaMalloc((void**) & c.yy_adc,        MAX_WORD_SIZE));
  cudaCheck(cudaMalloc((void**) & c.adc_d,         MAX_WORD_SIZE));
  cudaCheck(cudaMalloc((void**) & c.layer_d ,      MAX_WORD_SIZE));
  cudaCheck(cudaMalloc((void**) & c.rawIdArr_d,    MAX_WORD_SIZE));
  cudaCheck(cudaMalloc((void**) & c.errType_d,     MAX_WORD_SIZE));
  cudaCheck(cudaMalloc((void**) & c.errWord_d,     MAX_WORD_SIZE));
  cudaCheck(cudaMalloc((void**) & c.errFedID_d,    MAX_WORD_SIZE));
  cudaCheck(cudaMalloc((void**) & c.errRawID_d,    MAX_WORD_SIZE));
  cudaCheck(cudaMalloc((void**) & c.moduleId_d,    MAX_WORD_SIZE));
  cudaCheck(cudaMalloc((void**) & c.mIndexStart_d, MSIZE));
  cudaCheck(cudaMalloc((void**) & c.mIndexEnd_d,   MSIZE));

  // create a CUDA stream
  cudaCheck(cudaStreamCreate(&c.stream));

  return c;
}


void freeMemory(context & c) {
  // free the GPU memory
  cudaCheck(cudaFree(c.word_d));
  cudaCheck(cudaFree(c.fedIndex_d));
  cudaCheck(cudaFree(c.adc_d));
  cudaCheck(cudaFree(c.layer_d));
  cudaCheck(cudaFree(c.xx_d));
  cudaCheck(cudaFree(c.yy_d));
  cudaCheck(cudaFree(c.xx_adc));
  cudaCheck(cudaFree(c.yy_adc));
  cudaCheck(cudaFree(c.rawIdArr_d));
  cudaCheck(cudaFree(c.moduleId_d));
  cudaCheck(cudaFree(c.mIndexStart_d));
  cudaCheck(cudaFree(c.mIndexEnd_d));

  // destroy the CUDA stream
  cudaCheck(cudaStreamDestroy(c.stream));
}


__device__ uint32_t getLink(uint32_t ww)  {
  return ((ww >> LINK_shift) & LINK_mask);
}


__device__ uint32_t getRoc(uint32_t ww) {
  return ((ww >> ROC_shift ) & ROC_mask);
}


__device__ uint32_t getADC(uint32_t ww) {
  return ((ww >> ADC_shift) & ADC_mask);
}


__device__ bool isBarrel(uint32_t rawId) {
  return (1==((rawId>>25)&0x7));
}



__device__ DetIdGPU getRawId(const SiPixelFedCablingMapGPU * Map, uint32_t fed, uint32_t link, uint32_t roc) {
  uint32_t index = fed * MAX_LINK * MAX_ROC + (link-1) * MAX_ROC + roc;
  DetIdGPU detId = { Map->RawId[index], Map->rocInDet[index], Map->moduleId[index] };
  return detId;
}


//reference http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_2_0/doc/html/dd/d31/FrameConversion_8cc_source.html
//http://cmslxr.fnal.gov/source/CondFormats/SiPixelObjects/src/PixelROC.cc?v=CMSSW_9_2_0#0071
// Convert local pixel to global pixel
__device__ Pixel frameConversion(bool bpix, int side, uint32_t layer, uint32_t rocIdInDetUnit, Pixel local) {

  int slopeRow  = 0,  slopeCol = 0;
  int rowOffset = 0, colOffset = 0;

  if (bpix) {

    if (side == -1 && layer != 1) { // -Z side: 4 non-flipped modules oriented like 'dddd', except Layer 1
      if (rocIdInDetUnit < 8) {
        slopeRow = 1;
        slopeCol = -1;
        rowOffset = 0;
        colOffset = (8-rocIdInDetUnit)*numColsInRoc-1;
      }
      else {
        slopeRow  = -1;
        slopeCol  = 1;
        rowOffset = 2*numRowsInRoc-1;
        colOffset = (rocIdInDetUnit-8)*numColsInRoc;
      } // if roc
    }
    else { // +Z side: 4 non-flipped modules oriented like 'pppp', but all 8 in layer1
      if (rocIdInDetUnit < 8) {
        slopeRow  = -1;
        slopeCol  =  1;
        rowOffset = 2*numRowsInRoc-1;
        colOffset = rocIdInDetUnit * numColsInRoc;
      }
      else {
        slopeRow  = 1;
        slopeCol  = -1;
        rowOffset = 0;
        colOffset = (16-rocIdInDetUnit)*numColsInRoc-1;
      }
    }

  }
  else { // fpix
    if (side==-1) { // pannel 1
      if (rocIdInDetUnit < 8) {
        slopeRow = 1;
        slopeCol = -1;
        rowOffset = 0;
        colOffset = (8-rocIdInDetUnit)*numColsInRoc-1;
      }
      else {
        slopeRow = -1;
        slopeCol = 1;
        rowOffset = 2*numRowsInRoc-1;
        colOffset = (rocIdInDetUnit-8)*numColsInRoc;
      }
    }
    else { // pannel 2
      if (rocIdInDetUnit < 8) {
        slopeRow = 1;
        slopeCol = -1;
        rowOffset = 0;
        colOffset = (8-rocIdInDetUnit)*numColsInRoc-1;
      }
      else {
        slopeRow = -1;
        slopeCol = 1;
        rowOffset = 2*numRowsInRoc-1;
        colOffset = (rocIdInDetUnit-8)*numColsInRoc;
      }

    } // side

  }

  uint32_t gRow = rowOffset+slopeRow*local.row;
  uint32_t gCol = colOffset+slopeCol*local.col;
  //printf("Inside frameConversion row: %u, column: %u\n",gRow, gCol);
  Pixel global = {gRow, gCol};
  return global;
}


__device__ uint32_t conversionError(uint32_t fedId, uint32_t status, bool debug = false)
{

  uint32_t errorType = 0;

  switch (status) {
      case(1) : {
        if (debug) printf("Error in Fed: %i, invalid channel Id (errorType=35)", fedId );
        errorType = 35;
        break;
      }
      case(2) : {
        if (debug) printf("Error in Fed: %i, invalid ROC Id (errorType=36)", fedId);
        errorType = 36;
        break;
      }
      case(3) : {
        if (debug) printf("Error in Fed: %i, invalid dcol/pixel value (errorType=37)", fedId);
        errorType = 37;
        break;
      }
      case(4) : {
        if (debug) printf("Error in Fed: %i, dcol/pixel read out of order (errorType=38)", fedId);
        errorType = 38;
        break;
      }
      default: if (debug) printf("Cabling check returned unexpected result, status = %i", status);
  };

  return errorType;

}


__device__ bool rocRowColIsValid(uint32_t rocRow, uint32_t rocCol)
{
    uint32_t numRowsInRoc = 80;
    uint32_t numColsInRoc = 52;

    /// row and collumn in ROC representation
    return ( (rocRow < numRowsInRoc) & (rocCol < numColsInRoc) );
}


__device__ bool dcolIsValid(uint32_t dcol, uint32_t pxid)
{
    return ( (dcol < 26) &  (2 <= pxid) & (pxid < 162) );
}


__device__ uint32_t checkROC(uint32_t errorWord, uint32_t fedId, uint32_t link, const SiPixelFedCablingMapGPU *Map, bool debug = false)
{

 int errorType = (errorWord >> ROC_shift) & ERROR_mask;
 bool errorFound = false;

 switch (errorType) {
    case(25) : {
     uint32_t index = fedId * MAX_LINK * MAX_ROC + (link-1) * MAX_ROC + 1;
     if (index > 1 && index <= Map->size){
       if (!(link == Map->link[index] && 1 == Map->roc[index])) errorFound = false;
     }
     else{
       errorFound = true;
       if (debug) printf("Invalid ROC = 25 found (errorType=25)");
     }
     break;
   }
   case(26) : {
     if (debug) printf("Gap word found (errorType=26)");
     errorFound = true;
     break;
   }
   case(27) : {
     if (debug) printf("Dummy word found (errorType=27)");
     errorFound = true;
     break;
   }
   case(28) : {
     if (debug) printf("Error fifo nearly full (errorType=28)");
     errorFound = true;
     break;
   }
   case(29) : {
     if (debug) printf("Timeout on a channel (errorType=29)");
     if ((errorWord >> OMIT_ERR_shift) & OMIT_ERR_mask) {
       if (debug) printf("...first errorType=29 error, this gets masked out");
     }
     errorFound = true;
     break;
   }
   case(30) : {
     if (debug) printf("TBM error trailer (errorType=30)");
     int StateMatch_bits      = 4;
     int StateMatch_shift     = 8;
     uint32_t StateMatch_mask = ~(~uint32_t(0) << StateMatch_bits);
     int StateMatch = (errorWord >> StateMatch_shift) & StateMatch_mask;
     if ( StateMatch != 1 && StateMatch != 8 ) {
       if (debug) printf("FED error 30 with unexpected State Bits (errorType=30)");
     }
     if ( StateMatch==1 ) errorType = 40; // 1=Overflow -> 40, 8=number of ROCs -> 30
     errorFound = true;
     break;
   }
   case(31) : {
     if (debug) printf("Event number error (errorType=31)");
     errorFound = true;
     break;
   }
   default: errorFound = false;

 };

 return errorFound? errorType : 0;

}


__device__ uint32_t getErrRawID(uint32_t fedId, uint32_t errWord, uint32_t errorType, const SiPixelFedCablingMapGPU *Map, bool debug = false)
{

  uint32_t rID = 0xffffffff;

  switch (errorType) {
    case  25 : case  30 : case  31 : case  36 : case 40 : {
      // set dummy values for cabling just to get detId from link
      //cabling.dcol = 0;
      //cabling.pxid = 2;
      uint32_t roc  = 1;
      uint32_t link = (errWord >> LINK_shift) & LINK_mask;

      rID = getRawId(Map, fedId, link, roc).RawId;
      break;
    }
    case  29 : {
      int chanNmbr = 0;
      const int DB0_shift = 0;
      const int DB1_shift = DB0_shift + 1;
      const int DB2_shift = DB1_shift + 1;
      const int DB3_shift = DB2_shift + 1;
      const int DB4_shift = DB3_shift + 1;
      const uint32_t DataBit_mask = ~(~uint32_t(0) << 1);

      int CH1 = (errWord >> DB0_shift) & DataBit_mask;
      int CH2 = (errWord >> DB1_shift) & DataBit_mask;
      int CH3 = (errWord >> DB2_shift) & DataBit_mask;
      int CH4 = (errWord >> DB3_shift) & DataBit_mask;
      int CH5 = (errWord >> DB4_shift) & DataBit_mask;
      int BLOCK_bits      = 3;
      int BLOCK_shift     = 8;
      uint32_t BLOCK_mask = ~(~uint32_t(0) << BLOCK_bits);
      int BLOCK = (errWord >> BLOCK_shift) & BLOCK_mask;
      int localCH = 1*CH1+2*CH2+3*CH3+4*CH4+5*CH5;
      if (BLOCK%2==0) chanNmbr=(BLOCK/2)*9+localCH;
      else chanNmbr = ((BLOCK-1)/2)*9+4+localCH;
      if ((chanNmbr < 1)||(chanNmbr > 36)) break;  // signifies unexpected result

      // set dummy values for cabling just to get detId from link if in Barrel
      //cabling.dcol = 0;
      //cabling.pxid = 2;
      uint32_t roc  = 1;
      uint32_t link = chanNmbr;
      rID = getRawId(Map, fedId, link, roc).RawId;
      break;
    }
    case  37 : case  38: {
      //cabling.dcol = 0;
      //cabling.pxid = 2;
      uint32_t roc  = (errWord >> ROC_shift) & ROC_mask;
      uint32_t link = (errWord >> LINK_shift) & LINK_mask;
      rID = getRawId(Map, fedId, link, roc).RawId;
      break;
    }

  default : break;

  };

  return rID;

}


/*----------
* Name: applyADCthreshold_kernel()
* Desc: converts adc count to electrons and then applies the
* threshold on each channel.
* make pixel to 0 if it is below the threshold
* Input: xx_d[], yy_d[], layer_d[], wordCounter, adc[], ADCThreshold
*-----------
* Output: xx_adc[], yy_adc[] with pixel threshold applied
*/
// kernel to apply adc threshold on the channels


// Felice: gains and pedestals are not the same for each pixel. This code should be rewritten to take
// in account local gains/pedestals
// __global__ void applyADCthreshold_kernel(const uint32_t *xx_d, const uint32_t *yy_d, const uint32_t *layer_d, uint32_t *adc, const uint32_t wordCounter,
//  const ADCThreshold adcThreshold, uint32_t *xx_adc, uint32_t *yy_adc ) {
//   int tid = threadIdx.x;
//   int gIndex = blockDim.x*blockIdx.x+tid;
//   if (gIndex<wordCounter) {
//     uint32_t adcOld = adc[gIndex];
//     const float gain = adcThreshold.theElectronPerADCGain_; // default: 1 ADC = 135 electrons
//     const float pedestal = 0; //
//     int adcNew = int(adcOld*gain+pedestal);
//     // rare chance of entering into the if ()
//     if (layer_d[gIndex]>=adcThreshold.theFirstStack_) {
//       if (adcThreshold.theStackADC_==1 && adcOld==1) {
//         adcNew = int(255*135); // Arbitrarily use overflow value.
//       }
//       if (adcThreshold.theStackADC_ >1 && adcThreshold.theStackADC_!=255 && adcOld>=1){
//         adcNew = int((adcOld-1) * gain * 255/float(adcThreshold.theStackADC_-1));
//       }
//     }
//
//     if (adcNew >adcThreshold.thePixelThreshold ) {
//       xx_adc[gIndex]=xx_d[gIndex];
//       yy_adc[gIndex]=yy_d[gIndex];
//     }
//     else {
//       xx_adc[gIndex]=0; // 0: dead pixel
//       yy_adc[gIndex]=0;
//     }
//     adc[gIndex] = adcNew;
//   }
// }


// Kernel to perform Raw to Digi conversion
__global__ void RawToDigi_kernel(const SiPixelFedCablingMapGPU *Map, const uint32_t wordCounter, const uint32_t *Word, const uint32_t *fedIndex,
                                 uint32_t *XX, uint32_t *YY, uint32_t *moduleId, int *mIndexStart,
                                 int *mIndexEnd, uint32_t *ADC, uint32_t *layerArr, uint32_t *rawIdArr,
                                 uint32_t *errType, uint32_t *errWord, uint32_t *errFedID, uint32_t *errRawID,
                                 bool useQualityInfo, bool includeErrors, bool debug)
{
  uint32_t blockId  = blockIdx.x;

  //const uint32_t eventOffset  = eventIndex[eventno];
  uint32_t fedId     = fedIndex[blockId];
  uint32_t threadId  = threadIdx.x;

  uint32_t begin  = fedIndex[MAX_FED + blockId];
  uint32_t end    = fedIndex[MAX_FED + blockId + 1];

  if (blockIdx.x == gridDim.x - 1) {
    end = wordCounter; // for last fed to get the end index
  }

  uint32_t link = 0;
  uint32_t roc  = 0;

  bool skipROC = false;
  //if (threadId==0) printf("Event: %u blockId: %u start: %u end: %u\n", eventno, blockId, begin, end);
  int no_itr = (end - begin)/blockDim.x + 1; // to deal with number of hits greater than blockDim.x
  for (int i = 0; i < no_itr; i++) {
    auto gIndex = begin + threadId + i*blockDim.x;
    if (gIndex < end) {
      uint32_t ww = Word[gIndex]; // Array containing 32 bit raw data
      if (includeErrors) {
        errType[gIndex]  = 0;
        errWord[gIndex]  = ww;
        errFedID[gIndex] = fedId;
        errRawID[gIndex] = 0;
      }
      if (ww == 0) {
        //noise and dead channels are ignored
        XX[gIndex]    = 0;  // 0 is an indicator of a noise/dead channel
        YY[gIndex]    = 0; // skip these pixels during clusterization
        ADC[gIndex]   = 0;
        layerArr[gIndex] = 0;
        moduleId[gIndex] = 9999; //9999 is the indication of bad module, taken care later
        rawIdArr[gIndex] = 9999;
        continue ; // 0: bad word,
      }

      uint32_t nlink  = getLink(ww);            // Extract link
      uint32_t nroc   = getRoc(ww);             // Extract Roc in link
      if (!((nlink != link) | (nroc != roc))) continue;
      link = nlink;
      roc = nroc;
      DetIdGPU detId = getRawId(Map, fedId, link, roc);


      uint32_t errorType = checkROC(ww, fedId, link, Map, debug);
      skipROC = (roc < maxROCIndex) ? false : (errorType != 0);
      if (includeErrors and skipROC)
      {
        uint32_t rID = getErrRawID(fedId, ww, errorType, Map, debug); //write the function
        errType[gIndex]  = errorType;
        errWord[gIndex]  = ww;
        errFedID[gIndex] = fedId;
        errRawID[gIndex] = rID;
        continue;
      }


      uint32_t rawId  = detId.RawId;
      uint32_t rocIdInDetUnit = detId.rocInDet;

      bool barrel = isBarrel(rawId);

      uint32_t index = fedId * MAX_LINK * MAX_ROC + (link-1) * MAX_ROC + roc;
      if (useQualityInfo) {

          skipROC = Map->badRocs[index];
          if (skipROC) continue;

      }
      skipROC = Map->modToUnp[index];
      if (skipROC) continue;

      uint32_t layer = 0;//, ladder =0;
      int side = 0, panel = 0, module = 0;//disk = 0,blade = 0

      if (barrel)
      {
        layer  = (rawId >> layerStartBit_) & layerMask_;
        module = (rawId >> moduleStartBit_) & moduleMask_;
        side   = (module<5)? -1 : 1;
      }
      else {
        // endcap ids
        layer = 0;
        panel = (rawId >> panelStartBit_) & panelMask_;
        //disk  = (rawId >> diskStartBit_)  & diskMask_ ;
        side  = (panel==1)? -1 : 1;
        //blade = (rawId>>bladeStartBit_) & bladeMask_;
      }

      // ***special case of layer to 1 be handled here
      Pixel localPix;
      if (layer == 1) {
        uint32_t col = (ww >> COL_shift) & COL_mask;
        uint32_t row = (ww >> ROW_shift) & ROW_mask;
        localPix.row = row;
        localPix.col = col;
        if (includeErrors) {
          if (not rocRowColIsValid(row, col)) {
            uint32_t error = conversionError(fedId, 3, debug); //use the device function and fill the arrays
            errType[gIndex]  = error;
            errWord[gIndex]  = ww;
            errFedID[gIndex] = fedId;
            errRawID[gIndex] = rawId;
            printf("Error status: %i\n", error);
            continue;
          }
        }
      } else {
        // ***conversion rules for dcol and pxid
        uint32_t dcol = (ww >> DCOL_shift) & DCOL_mask;
        uint32_t pxid = (ww >> PXID_shift) & PXID_mask;
        uint32_t row  = numRowsInRoc - pxid/2;
        uint32_t col  = dcol*2 + pxid%2;
        localPix.row = row;
        localPix.col = col;
        if (includeErrors and not dcolIsValid(dcol, pxid)) {
          uint32_t error = conversionError(fedId, 3, debug);
          errType[gIndex] = error;
          errWord[gIndex] = ww;
          errFedID[gIndex] = fedId;
          errRawID[gIndex] = rawId;
          printf("Error status: %i\n", error);
          continue;
        }
      }

      Pixel globalPix = frameConversion(barrel, side, layer, rocIdInDetUnit, localPix);
      //printf("GPU side: %i, layer: %i, roc: %i, lrow: %i, lcol: %i, grow: %i, gcol: %i, word: %i\n", side, layer, rocIdInDetUnit, localPix.row, localPix.col, globalPix.row, globalPix.col, ww);
      XX[gIndex]    = globalPix.row  ; // origin shifting by 1 0-159
      YY[gIndex]    = globalPix.col ; // origin shifting by 1 0-415
      ADC[gIndex]   = getADC(ww);
      layerArr[gIndex] = layer;
      moduleId[gIndex] = detId.moduleId;
      rawIdArr[gIndex] = rawId;
      if (includeErrors) {
        //fill error type
        errType[gIndex]  = 0;
        errWord[gIndex]  = ww;
        errFedID[gIndex] = fedId;
        errRawID[gIndex] = rawId;
      }
    } // end of if (gIndex < end)
  } // end of for(int i =0;i<no_itr...)

  __syncthreads();

  // three cases possible
  // case 1: 21 21 21 22 21 22 22
  // pos   : 0  1  2  3  4  5  6
  // solution swap 21 with 22 : 21 21 21 21 22 22 22
  // atomicExch(address, value), set the variable at address to value.
  // do the swapping for above case and replace the 9999 with
  // valid moduleId

  for(int i = 0; i < no_itr; i++) {
    uint32_t gIndex = begin + threadId + i*blockDim.x;
    if (gIndex+2 < end) {
      //rare condition
      if (moduleId[gIndex] == moduleId[gIndex+2] and moduleId[gIndex] < moduleId[gIndex+1]) {
        atomicExch(&moduleId[gIndex+2], atomicExch(&moduleId[gIndex+1], moduleId[gIndex+2]));
        //*swap all the digi id
        atomicExch(&XX[gIndex+2], atomicExch(&XX[gIndex+1], XX[gIndex+2]));
        atomicExch(&YY[gIndex+2], atomicExch(&YY[gIndex+1], YY[gIndex+2]));
        atomicExch(&ADC[gIndex+2], atomicExch(&ADC[gIndex+1], ADC[gIndex+2]));
        atomicExch(&layerArr[gIndex+2], atomicExch(&layerArr[gIndex+1], layerArr[gIndex+2]));
        atomicExch(&rawIdArr[gIndex+2], atomicExch(&rawIdArr[gIndex+1], rawIdArr[gIndex+2]));
      }
      __syncthreads();

      //rarest condition
      // above condition fails at 361 361 361 363 362 363 363
      // here we need to swap 362 with previous 363
      if (moduleId[gIndex]==moduleId[gIndex+2] && moduleId[gIndex]>moduleId[gIndex+1]) {
        atomicExch(&moduleId[gIndex+1], atomicExch(&moduleId[gIndex], moduleId[gIndex+1]));
        //*swap all the digi id
        atomicExch(&XX[gIndex+1], atomicExch(&XX[gIndex], XX[gIndex+1]));
        atomicExch(&YY[gIndex+1], atomicExch(&YY[gIndex], YY[gIndex+1]));
        atomicExch(&ADC[gIndex+1], atomicExch(&ADC[gIndex], ADC[gIndex+1]));
        atomicExch(&layerArr[gIndex+1], atomicExch(&layerArr[gIndex], layerArr[gIndex+1]));
        atomicExch(&rawIdArr[gIndex+1], atomicExch(&rawIdArr[gIndex], rawIdArr[gIndex+1]));
      }

      // moduleId== 9999 then pixel is bad with x=y=layer=adc=0
      // this bad pixel will not affect the cluster, since for cluster
      // the origin is shifted at (1,1) so x=y=0 will be ignored
      // assign the previous valid moduleId to this pixel to remove 9999
      // so that we can get the start & end index of module easily.
      __syncthreads(); // let the swapping finish first

      if (moduleId[gIndex] == 9999) {
        int m=gIndex;
        while(moduleId[--m] == 9999) {} //skip till you get the valid module
        moduleId[gIndex] = moduleId[m];
      }
    } // end of if (gIndex<end)
  } //  end of for(int i=0;i<no_itr;...)
  __syncthreads();

  // mIndexStart stores starting index of module
  // mIndexEnd stores end index of module
  // both indexes are inclusive
  // check consecutive module numbers
  // for start of fed
  for(int i = 0; i < no_itr; i++) {
    uint32_t gIndex = begin + threadId + i*blockDim.x;
    uint32_t moduleOffset = NMODULE;
    //if (threadId==0) printf("moduleOffset: %u\n",moduleOffset );
    if (gIndex < end) {
      if (gIndex == begin) {
        mIndexStart[moduleOffset+moduleId[gIndex]] = gIndex;
      }
      // for end of the fed
      if (gIndex == (end-1)) {
        mIndexEnd[moduleOffset+moduleId[gIndex]] = gIndex;
      }
      // point to the gIndex where two consecutive moduleId varies
      if (gIndex!= begin && (gIndex<(end-1)) && moduleId[gIndex]!=9999) {
        if (moduleId[gIndex]<moduleId[gIndex+1] ) {
          mIndexEnd[moduleOffset + moduleId[gIndex]] = gIndex;
        }
        if (moduleId[gIndex] > moduleId[gIndex-1] ) {
          mIndexStart[moduleOffset+ moduleId[gIndex]] = gIndex;
        }
      } //end of if (gIndex!= begin && (gIndex<(end-1)) ...
    } //end of if (gIndex <end)
  }
} // end of Raw to Digi kernel


// kernel wrapper called from runRawToDigi_kernel
void RawToDigi_wrapper(
    context & c,
    const SiPixelFedCablingMapGPU* cablingMapDevice, const uint32_t wordCounter, uint32_t *word, const uint32_t fedCounter,  uint32_t *fedIndex,
    bool convertADCtoElectrons, uint32_t *xx_h, uint32_t *yy_h, uint32_t *adc_h, int *mIndexStart_h,
    int *mIndexEnd_h, uint32_t *rawIdArr_h, uint32_t *errType_h, uint32_t *errWord_h, uint32_t *errFedID_h, uint32_t *errRawID_h,
    bool useQualityInfo, bool includeErrors, bool debug)
{
  const int threadsPerBlock = 512;
  const int blocks = 108; // only 108 feds are present

  int MSIZE = NMODULE*sizeof(int)+sizeof(int);
  // initialize moduleStart & moduleEnd with some constant(-1)
  // just to check if it updated in kernel or not
  cudaCheck(cudaMemsetAsync(c.mIndexStart_d, -1, MSIZE, c.stream));
  cudaCheck(cudaMemsetAsync(c.mIndexEnd_d, -1, MSIZE, c.stream));

  int FSIZE     = (2*MAX_FED +1)*sizeof(uint32_t); // 0 to 150:fedId, 150:300: fedIndex
  // wordCounter is the total no of words in each event to be trasfered on device
  cudaCheck(cudaMemcpyAsync(&c.word_d[0],     &word[0],     wordCounter*sizeof(uint32_t), cudaMemcpyHostToDevice, c.stream));
  cudaCheck(cudaMemcpyAsync(&c.fedIndex_d[0], &fedIndex[0], FSIZE, cudaMemcpyHostToDevice, c.stream));

  // Launch rawToDigi kernel
  RawToDigi_kernel<<<blocks, threadsPerBlock, 0, c.stream>>>(
      cablingMapDevice,
      wordCounter,
      c.word_d,
      c.fedIndex_d,
      c.xx_d,
      c.yy_d,
      c.moduleId_d,
      c.mIndexStart_d,
      c.mIndexEnd_d,
      c.adc_d,c.layer_d,
      c.rawIdArr_d,
      c.errType_d,
      c.errWord_d,
      c.errFedID_d,
      c.errRawID_d,
      useQualityInfo,
      includeErrors,
      debug);
  cudaCheck(cudaGetLastError());

  // copy data to host variable
  // if you want to copy data after applying ADC threshold
  if (convertADCtoElectrons) {
    cudaCheck(cudaMemcpy(xx_h, c.xx_adc, wordCounter * sizeof(uint32_t), cudaMemcpyDeviceToHost));
    cudaCheck(cudaMemcpy(yy_h, c.yy_adc, wordCounter * sizeof(uint32_t), cudaMemcpyDeviceToHost));
  } else {
    cudaCheck(cudaMemcpy(xx_h, c.xx_d, wordCounter*sizeof(uint32_t), cudaMemcpyDeviceToHost));
    cudaCheck(cudaMemcpy(yy_h, c.yy_d, wordCounter*sizeof(uint32_t), cudaMemcpyDeviceToHost));
  }

  cudaCheck(cudaMemcpyAsync(adc_h, c.adc_d, wordCounter*sizeof(uint32_t), cudaMemcpyDeviceToHost, c.stream));
  cudaCheck(cudaMemcpyAsync(rawIdArr_h, c.rawIdArr_d, wordCounter*sizeof(uint32_t), cudaMemcpyDeviceToHost, c.stream));

  cudaCheck(cudaMemcpyAsync(mIndexStart_h, c.mIndexStart_d, NMODULE*sizeof(int), cudaMemcpyDeviceToHost, c.stream));
  cudaCheck(cudaMemcpyAsync(mIndexEnd_h, c.mIndexEnd_d, NMODULE*sizeof(int), cudaMemcpyDeviceToHost, c.stream));

  if (includeErrors) {
    cudaCheck(cudaMemcpyAsync(errType_h, c.errType_d, wordCounter*sizeof(uint32_t), cudaMemcpyDeviceToHost, c.stream));
    cudaCheck(cudaMemcpyAsync(errWord_h, c.errWord_d, wordCounter*sizeof(uint32_t), cudaMemcpyDeviceToHost, c.stream));
    cudaCheck(cudaMemcpyAsync(errFedID_h, c.errFedID_d, wordCounter*sizeof(uint32_t), cudaMemcpyDeviceToHost, c.stream));
    cudaCheck(cudaMemcpyAsync(errRawID_h, c.errRawID_d, wordCounter*sizeof(uint32_t), cudaMemcpyDeviceToHost, c.stream));
  }
  cudaStreamSynchronize(c.stream);
  // End  of Raw2Digi and passing data for cluserisation
  // PixelCluster_Wrapper(c.xx_adc , c.yy_adc, c.adc_d,wordCounter, c.mIndexStart_d, c.mIndexEnd_d);
}
