/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_MP_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::At11"),  CCTK_VarIndex("ML_BSSN_MP::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::At12"),  CCTK_VarIndex("ML_BSSN_MP::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::At13"),  CCTK_VarIndex("ML_BSSN_MP::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::At22"),  CCTK_VarIndex("ML_BSSN_MP::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::At23"),  CCTK_VarIndex("ML_BSSN_MP::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::At33"),  CCTK_VarIndex("ML_BSSN_MP::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::A"),  CCTK_VarIndex("ML_BSSN_MP::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::B1"),  CCTK_VarIndex("ML_BSSN_MP::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::B2"),  CCTK_VarIndex("ML_BSSN_MP::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::B3"),  CCTK_VarIndex("ML_BSSN_MP::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::Xt1"),  CCTK_VarIndex("ML_BSSN_MP::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::Xt2"),  CCTK_VarIndex("ML_BSSN_MP::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::Xt3"),  CCTK_VarIndex("ML_BSSN_MP::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::alpha"),  CCTK_VarIndex("ML_BSSN_MP::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::phi"),  CCTK_VarIndex("ML_BSSN_MP::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::gt11"),  CCTK_VarIndex("ML_BSSN_MP::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::gt12"),  CCTK_VarIndex("ML_BSSN_MP::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::gt13"),  CCTK_VarIndex("ML_BSSN_MP::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::gt22"),  CCTK_VarIndex("ML_BSSN_MP::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::gt23"),  CCTK_VarIndex("ML_BSSN_MP::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::gt33"),  CCTK_VarIndex("ML_BSSN_MP::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::beta1"),  CCTK_VarIndex("ML_BSSN_MP::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::beta2"),  CCTK_VarIndex("ML_BSSN_MP::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::beta3"),  CCTK_VarIndex("ML_BSSN_MP::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP::trK"),  CCTK_VarIndex("ML_BSSN_MP::trKrhs"));
  return;
}
