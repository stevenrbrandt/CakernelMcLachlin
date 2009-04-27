/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (64-bit) (May 21, 2008) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_ADM_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K11"),  CCTK_VarIndex("ML_ADM::K11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K12"),  CCTK_VarIndex("ML_ADM::K12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K13"),  CCTK_VarIndex("ML_ADM::K13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K22"),  CCTK_VarIndex("ML_ADM::K22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K23"),  CCTK_VarIndex("ML_ADM::K23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::K33"),  CCTK_VarIndex("ML_ADM::K33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::alpha"),  CCTK_VarIndex("ML_ADM::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g11"),  CCTK_VarIndex("ML_ADM::g11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g12"),  CCTK_VarIndex("ML_ADM::g12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g13"),  CCTK_VarIndex("ML_ADM::g13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g22"),  CCTK_VarIndex("ML_ADM::g22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g23"),  CCTK_VarIndex("ML_ADM::g23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::g33"),  CCTK_VarIndex("ML_ADM::g33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::beta1"),  CCTK_VarIndex("ML_ADM::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::beta2"),  CCTK_VarIndex("ML_ADM::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_ADM::beta3"),  CCTK_VarIndex("ML_ADM::beta3rhs"));
  return;
}
