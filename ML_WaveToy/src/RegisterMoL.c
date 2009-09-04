/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_WaveToy_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToy::rho"),  CCTK_VarIndex("ML_WaveToy::rhorhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToy::u"),  CCTK_VarIndex("ML_WaveToy::urhs"));
  return;
}
