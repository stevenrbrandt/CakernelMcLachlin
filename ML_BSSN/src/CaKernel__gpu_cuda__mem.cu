#include"CaKernel__gpu_cuda__shared.h"
#include"CaKernel__gpu_cuda__mem.h"
/* Assume Piraha will generate this file and this file will be pushed here as well */

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#ifdef __CUDACC__ 
#include <cuda.h> 
#include <cuda_runtime.h> 
#endif

#include<assert.h>
#include<algorithm>
using namespace std;

 



Vars d_vars;
Pars d_pars;

void CaKernel_AllocDevVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  
}

void CaKernel_FreeDevVars(){
  
}

void CaKernel_AllocDevPars(CCTK_ARGUMENTS)
{
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int numvars =  13 ;
  size_t sizes[ 13  + 1] = {
    
      sizeof (CCTK_REAL), 
      sizeof (CCTK_REAL), 
      sizeof (CCTK_REAL), 
      sizeof (CCTK_REAL), 
      sizeof (CCTK_REAL), 
      sizeof (CCTK_INT), 
      sizeof (CCTK_REAL), 
      sizeof (CCTK_INT), 
      sizeof (CCTK_REAL), 
      sizeof (CCTK_REAL), 
      sizeof (CCTK_REAL), 
      sizeof (CCTK_INT), 
      sizeof (CCTK_REAL), 
      16 // to make sure that the "structure" is aligned to 16 bytes
    };   

  size_t elements[ 13  + 1] = {
      1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, 
      0
    };   

  int offsets[ 13  + 1];   
  
  offsets[0] = 0;
  for (int i = 1; i <= numvars; i++){
    int oldoffset = offsets[i - 1] + sizes[i - 1] * elements[i - 1];
    for(int b = 2; b <= 16 && b <= sizes[i]; b <<= 1)
      oldoffset += (b - (oldoffset % b)) % b;
    offsets[i] =  oldoffset;
  }  
  
   
    d_pars.MinimumLapse_offset = offsets[(int)0];  
    d_pars.LapseAdvectionCoeff_offset = offsets[(int)1.0];  
    d_pars.BetaDriver_offset = offsets[(int)2.0];  
    d_pars.EpsDiss_offset = offsets[(int)3.0];  
    d_pars.ShiftBCoeff_offset = offsets[(int)4.0];  
    d_pars.harmonicN_offset = offsets[(int)5.0];  
    d_pars.ShiftAdvectionCoeff_offset = offsets[(int)6.0];  
    d_pars.harmonicShift_offset = offsets[(int)7.0];  
    d_pars.harmonicF_offset = offsets[(int)8.0];  
    d_pars.ShiftGammaCoeff_offset = offsets[(int)9.0];  
    d_pars.LapseACoeff_offset = offsets[(int)10.0];  
    d_pars.conformalMethod_offset = offsets[(int)11.0];  
    d_pars.AlphaDriver_offset = offsets[(int)12.0]; 

//  d_pars = pars;
  char * tmpptr = (char *)malloc (offsets[numvars]);
  CUDA_SAFE_CALL (cudaMalloc ((void **) &(d_pars.ptr), offsets[numvars]));
  
   
  { 
    CCTK_REAL * ptr = (CCTK_REAL *) &tmpptr[d_pars.MinimumLapse_offset];
    *ptr = MinimumLapse;    
  }  
  { 
    CCTK_REAL * ptr = (CCTK_REAL *) &tmpptr[d_pars.LapseAdvectionCoeff_offset];
    *ptr = LapseAdvectionCoeff;    
  }  
  { 
    CCTK_REAL * ptr = (CCTK_REAL *) &tmpptr[d_pars.BetaDriver_offset];
    *ptr = BetaDriver;    
  }  
  { 
    CCTK_REAL * ptr = (CCTK_REAL *) &tmpptr[d_pars.EpsDiss_offset];
    *ptr = EpsDiss;    
  }  
  { 
    CCTK_REAL * ptr = (CCTK_REAL *) &tmpptr[d_pars.ShiftBCoeff_offset];
    *ptr = ShiftBCoeff;    
  }  
  { 
    CCTK_INT * ptr = (CCTK_INT *) &tmpptr[d_pars.harmonicN_offset];
    *ptr = harmonicN;    
  }  
  { 
    CCTK_REAL * ptr = (CCTK_REAL *) &tmpptr[d_pars.ShiftAdvectionCoeff_offset];
    *ptr = ShiftAdvectionCoeff;    
  }  
  { 
    CCTK_INT * ptr = (CCTK_INT *) &tmpptr[d_pars.harmonicShift_offset];
    *ptr = harmonicShift;    
  }  
  { 
    CCTK_REAL * ptr = (CCTK_REAL *) &tmpptr[d_pars.harmonicF_offset];
    *ptr = harmonicF;    
  }  
  { 
    CCTK_REAL * ptr = (CCTK_REAL *) &tmpptr[d_pars.ShiftGammaCoeff_offset];
    *ptr = ShiftGammaCoeff;    
  }  
  { 
    CCTK_REAL * ptr = (CCTK_REAL *) &tmpptr[d_pars.LapseACoeff_offset];
    *ptr = LapseACoeff;    
  }  
  { 
    CCTK_INT * ptr = (CCTK_INT *) &tmpptr[d_pars.conformalMethod_offset];
    *ptr = conformalMethod;    
  }  
  { 
    CCTK_REAL * ptr = (CCTK_REAL *) &tmpptr[d_pars.AlphaDriver_offset];
    *ptr = AlphaDriver;    
  } 
  
  CUDA_SAFE_CALL (cudaMemcpy(d_pars.ptr, tmpptr, offsets[numvars], cudaMemcpyHostToDevice));
  free(tmpptr);
  
}

void CaKernel_FreeDevPars()
{
  
//  free(pars.ptr);
  CUDA_SAFE_CALL(cudaFree ((void *) d_pars.ptr));
  
}

void CaKernel_AllocDevMem (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  int num_tl = 1;
  int vi = -1;
  
    assert((vi = CCTK_VarIndex("ML_BSSN::A"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::alpha"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::alpharhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::Arhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At11"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At11rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At12"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At12rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At13"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At13rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At22"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At22rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At23"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At23rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At33"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::At33rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::B1"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::B1rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::B2"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::B2rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::B3"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::B3rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::beta1"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::beta1rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::beta2"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::beta2rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::beta3"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::beta3rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt11"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt11rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt12"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt12rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt13"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt13rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt22"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt22rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt23"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt23rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt33"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::gt33rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::phi"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::phirhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::trK"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::trKrhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::Xt1"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::Xt1rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::Xt2"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::Xt2rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::Xt3"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  
    assert((vi = CCTK_VarIndex("ML_BSSN::Xt3rhs"))>=0);
    num_tl = CCTK_MaxTimeLevelsVI(vi);
    Device_RegisterMem(CCTK_PASS_CTOC, vi, num_tl); 
  

  
  CaKernel_AllocDevPars(CCTK_PASS_CTOC);
  CaKernel_AllocDevVars(CCTK_PASS_CTOC);
}

/**
 * Free the memory of grid variables on GPU devices.
 */

int CaKernel_FreeDevMem () {
  
  CaKernel_FreeDevPars();
  CaKernel_FreeDevVars();
}


void CaKernel_InitDevice(CCTK_ARGUMENTS){
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int num_dev; 
  CUDA_SAFE_CALL (cudaGetDeviceCount (&num_dev));
  if (num_dev == 0)
  {
    CCTK_WARN
      (0, "There are no CUDA devices available. (They may be busy, or the driver may not be installed.)");
    exit(-1);
  }

  int myproc, devidx;

  myproc = CCTK_MyProc(cctkGH);
  devidx = myproc % num_dev;

  /* we set device based on the number of devices available on each node */
  CUDA_SAFE_CALL (cudaSetDevice (devidx));
  CCTK_VInfo(CCTK_THORNSTRING, "number of device %d", num_dev);
  CCTK_VInfo(CCTK_THORNSTRING, "device %d is successfully assigned to process %d", devidx, myproc);
}



