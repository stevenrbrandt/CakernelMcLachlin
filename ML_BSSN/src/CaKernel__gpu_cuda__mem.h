#ifndef CAKERNEL__GPU_CUDA__MEM_H
#define CAKERNEL__GPU_CUDA__MEM_H

/* definition of CCTK_REAL */
#include "cctk.h"

#ifdef __CUDACC__

/* device pointers */
/* All the external variables will be declared here.
 * We will generate the this file from the parser directly but we
 * need to generate the memory allocation routines as well.
 */



struct Vars {
  char  * ptr;
  union
  {
    struct
    {
      
   
      
    
      
    
      unsigned int all;
      unsigned int sepinout;
      unsigned int sepinout_out;
    };
    unsigned int field[];
  };
} ;

extern Vars d_vars;

struct Pars {
  char  * ptr;

  
  unsigned int MinimumLapse_offset; 
  unsigned int LapseAdvectionCoeff_offset; 
  unsigned int BetaDriver_offset; 
  unsigned int EpsDiss_offset; 
  unsigned int ShiftBCoeff_offset; 
  unsigned int harmonicN_offset; 
  unsigned int ShiftAdvectionCoeff_offset; 
  unsigned int harmonicShift_offset; 
  unsigned int harmonicF_offset; 
  unsigned int ShiftGammaCoeff_offset; 
  unsigned int LapseACoeff_offset; 
  unsigned int conformalMethod_offset; 
  unsigned int AlphaDriver_offset; 
  unsigned int all;
} ;

extern Pars d_pars;

#endif

#endif /* CAKERNEL__GPU_CUDA__MEM_H */
