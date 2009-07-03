/*  File produced by user eschnett */
/*  Produced with Mathematica Version 7.0 for Mac OS X x86 (64-bit) (February 19, 2009) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_hydro_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_hydro::ene"),  CCTK_VarIndex("ML_hydro::enerhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_hydro::mass"),  CCTK_VarIndex("ML_hydro::massrhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_hydro::mom1"),  CCTK_VarIndex("ML_hydro::mom1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_hydro::mom2"),  CCTK_VarIndex("ML_hydro::mom2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_hydro::mom3"),  CCTK_VarIndex("ML_hydro::mom3rhs"));
  return;
}
