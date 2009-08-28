/*  File produced by user eschnett */
/*  Produced with Mathematica Version 7.0 for Mac OS X x86 (64-bit) (November 11, 2008) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_BSSN_M_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::At11"),  CCTK_VarIndex("ML_BSSN_M::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::At12"),  CCTK_VarIndex("ML_BSSN_M::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::At13"),  CCTK_VarIndex("ML_BSSN_M::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::At22"),  CCTK_VarIndex("ML_BSSN_M::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::At23"),  CCTK_VarIndex("ML_BSSN_M::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::At33"),  CCTK_VarIndex("ML_BSSN_M::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::A"),  CCTK_VarIndex("ML_BSSN_M::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::B1"),  CCTK_VarIndex("ML_BSSN_M::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::B2"),  CCTK_VarIndex("ML_BSSN_M::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::B3"),  CCTK_VarIndex("ML_BSSN_M::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::Xt1"),  CCTK_VarIndex("ML_BSSN_M::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::Xt2"),  CCTK_VarIndex("ML_BSSN_M::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::Xt3"),  CCTK_VarIndex("ML_BSSN_M::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::alpha"),  CCTK_VarIndex("ML_BSSN_M::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::phi"),  CCTK_VarIndex("ML_BSSN_M::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::gt11"),  CCTK_VarIndex("ML_BSSN_M::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::gt12"),  CCTK_VarIndex("ML_BSSN_M::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::gt13"),  CCTK_VarIndex("ML_BSSN_M::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::gt22"),  CCTK_VarIndex("ML_BSSN_M::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::gt23"),  CCTK_VarIndex("ML_BSSN_M::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::gt33"),  CCTK_VarIndex("ML_BSSN_M::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::beta1"),  CCTK_VarIndex("ML_BSSN_M::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::beta2"),  CCTK_VarIndex("ML_BSSN_M::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::beta3"),  CCTK_VarIndex("ML_BSSN_M::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_M::trK"),  CCTK_VarIndex("ML_BSSN_M::trKrhs"));
  return;
}
