/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (64-bit) (May 21, 2008) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_BSSN6_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::At11"),  CCTK_VarIndex("ML_BSSN6::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::At12"),  CCTK_VarIndex("ML_BSSN6::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::At13"),  CCTK_VarIndex("ML_BSSN6::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::At22"),  CCTK_VarIndex("ML_BSSN6::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::At23"),  CCTK_VarIndex("ML_BSSN6::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::At33"),  CCTK_VarIndex("ML_BSSN6::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::A"),  CCTK_VarIndex("ML_BSSN6::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::B1"),  CCTK_VarIndex("ML_BSSN6::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::B2"),  CCTK_VarIndex("ML_BSSN6::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::B3"),  CCTK_VarIndex("ML_BSSN6::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::Xt1"),  CCTK_VarIndex("ML_BSSN6::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::Xt2"),  CCTK_VarIndex("ML_BSSN6::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::Xt3"),  CCTK_VarIndex("ML_BSSN6::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::alpha"),  CCTK_VarIndex("ML_BSSN6::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::phi"),  CCTK_VarIndex("ML_BSSN6::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::gt11"),  CCTK_VarIndex("ML_BSSN6::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::gt12"),  CCTK_VarIndex("ML_BSSN6::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::gt13"),  CCTK_VarIndex("ML_BSSN6::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::gt22"),  CCTK_VarIndex("ML_BSSN6::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::gt23"),  CCTK_VarIndex("ML_BSSN6::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::gt33"),  CCTK_VarIndex("ML_BSSN6::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::beta1"),  CCTK_VarIndex("ML_BSSN6::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::beta2"),  CCTK_VarIndex("ML_BSSN6::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::beta3"),  CCTK_VarIndex("ML_BSSN6::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN6::trK"),  CCTK_VarIndex("ML_BSSN6::trKrhs"));
  return;
}
