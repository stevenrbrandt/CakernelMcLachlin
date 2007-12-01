/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (32-bit) (April 20, 2007) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void MLFO_WaveToy_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("MLFO_WaveToy::rho"),  CCTK_VarIndex("MLFO_WaveToy::rhorhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("MLFO_WaveToy::u"),  CCTK_VarIndex("MLFO_WaveToy::urhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("MLFO_WaveToy::v1"),  CCTK_VarIndex("MLFO_WaveToy::v1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("MLFO_WaveToy::v2"),  CCTK_VarIndex("MLFO_WaveToy::v2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("MLFO_WaveToy::v3"),  CCTK_VarIndex("MLFO_WaveToy::v3rhs"));
  return;
}
