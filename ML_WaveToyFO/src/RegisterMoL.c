/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_WaveToyFO_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToyFO::rho"),  CCTK_VarIndex("ML_WaveToyFO::rhorhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToyFO::u"),  CCTK_VarIndex("ML_WaveToyFO::urhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToyFO::v1"),  CCTK_VarIndex("ML_WaveToyFO::v1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToyFO::v2"),  CCTK_VarIndex("ML_WaveToyFO::v2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToyFO::v3"),  CCTK_VarIndex("ML_WaveToyFO::v3rhs"));
  return;
}
