/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_BSSN_O2_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::At11"),  CCTK_VarIndex("ML_BSSN_O2::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::At12"),  CCTK_VarIndex("ML_BSSN_O2::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::At13"),  CCTK_VarIndex("ML_BSSN_O2::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::At22"),  CCTK_VarIndex("ML_BSSN_O2::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::At23"),  CCTK_VarIndex("ML_BSSN_O2::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::At33"),  CCTK_VarIndex("ML_BSSN_O2::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::A"),  CCTK_VarIndex("ML_BSSN_O2::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::B1"),  CCTK_VarIndex("ML_BSSN_O2::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::B2"),  CCTK_VarIndex("ML_BSSN_O2::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::B3"),  CCTK_VarIndex("ML_BSSN_O2::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::Xt1"),  CCTK_VarIndex("ML_BSSN_O2::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::Xt2"),  CCTK_VarIndex("ML_BSSN_O2::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::Xt3"),  CCTK_VarIndex("ML_BSSN_O2::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::alpha"),  CCTK_VarIndex("ML_BSSN_O2::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::phi"),  CCTK_VarIndex("ML_BSSN_O2::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::gt11"),  CCTK_VarIndex("ML_BSSN_O2::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::gt12"),  CCTK_VarIndex("ML_BSSN_O2::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::gt13"),  CCTK_VarIndex("ML_BSSN_O2::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::gt22"),  CCTK_VarIndex("ML_BSSN_O2::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::gt23"),  CCTK_VarIndex("ML_BSSN_O2::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::gt33"),  CCTK_VarIndex("ML_BSSN_O2::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::beta1"),  CCTK_VarIndex("ML_BSSN_O2::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::beta2"),  CCTK_VarIndex("ML_BSSN_O2::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::beta3"),  CCTK_VarIndex("ML_BSSN_O2::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_O2::trK"),  CCTK_VarIndex("ML_BSSN_O2::trKrhs"));
  return;
}
