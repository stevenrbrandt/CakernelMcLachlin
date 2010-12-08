/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_BSSN_UPW_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::At11"),  CCTK_VarIndex("ML_BSSN_UPW::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::At12"),  CCTK_VarIndex("ML_BSSN_UPW::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::At13"),  CCTK_VarIndex("ML_BSSN_UPW::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::At22"),  CCTK_VarIndex("ML_BSSN_UPW::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::At23"),  CCTK_VarIndex("ML_BSSN_UPW::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::At33"),  CCTK_VarIndex("ML_BSSN_UPW::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::A"),  CCTK_VarIndex("ML_BSSN_UPW::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::B1"),  CCTK_VarIndex("ML_BSSN_UPW::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::B2"),  CCTK_VarIndex("ML_BSSN_UPW::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::B3"),  CCTK_VarIndex("ML_BSSN_UPW::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::Xt1"),  CCTK_VarIndex("ML_BSSN_UPW::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::Xt2"),  CCTK_VarIndex("ML_BSSN_UPW::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::Xt3"),  CCTK_VarIndex("ML_BSSN_UPW::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::alpha"),  CCTK_VarIndex("ML_BSSN_UPW::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::phi"),  CCTK_VarIndex("ML_BSSN_UPW::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::gt11"),  CCTK_VarIndex("ML_BSSN_UPW::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::gt12"),  CCTK_VarIndex("ML_BSSN_UPW::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::gt13"),  CCTK_VarIndex("ML_BSSN_UPW::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::gt22"),  CCTK_VarIndex("ML_BSSN_UPW::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::gt23"),  CCTK_VarIndex("ML_BSSN_UPW::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::gt33"),  CCTK_VarIndex("ML_BSSN_UPW::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::beta1"),  CCTK_VarIndex("ML_BSSN_UPW::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::beta2"),  CCTK_VarIndex("ML_BSSN_UPW::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::beta3"),  CCTK_VarIndex("ML_BSSN_UPW::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_UPW::trK"),  CCTK_VarIndex("ML_BSSN_UPW::trKrhs"));
  return;
}
