/*  File produced by user eschnett */
/*  Produced with Mathematica Version 7.0 for Mac OS X x86 (64-bit) (February 19, 2009) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_FOWaveToy_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_FOWaveToy::rho"),  CCTK_VarIndex("ML_FOWaveToy::rhorhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_FOWaveToy::u"),  CCTK_VarIndex("ML_FOWaveToy::urhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_FOWaveToy::v1"),  CCTK_VarIndex("ML_FOWaveToy::v1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_FOWaveToy::v2"),  CCTK_VarIndex("ML_FOWaveToy::v2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_FOWaveToy::v3"),  CCTK_VarIndex("ML_FOWaveToy::v3rhs"));
  return;
}
