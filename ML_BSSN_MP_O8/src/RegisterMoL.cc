/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_MP_O8_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::At11"),  CCTK_VarIndex("ML_BSSN_MP_O8::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::At12"),  CCTK_VarIndex("ML_BSSN_MP_O8::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::At13"),  CCTK_VarIndex("ML_BSSN_MP_O8::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::At22"),  CCTK_VarIndex("ML_BSSN_MP_O8::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::At23"),  CCTK_VarIndex("ML_BSSN_MP_O8::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::At33"),  CCTK_VarIndex("ML_BSSN_MP_O8::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::A"),  CCTK_VarIndex("ML_BSSN_MP_O8::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::B1"),  CCTK_VarIndex("ML_BSSN_MP_O8::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::B2"),  CCTK_VarIndex("ML_BSSN_MP_O8::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::B3"),  CCTK_VarIndex("ML_BSSN_MP_O8::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::Xt1"),  CCTK_VarIndex("ML_BSSN_MP_O8::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::Xt2"),  CCTK_VarIndex("ML_BSSN_MP_O8::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::Xt3"),  CCTK_VarIndex("ML_BSSN_MP_O8::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::alpha"),  CCTK_VarIndex("ML_BSSN_MP_O8::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::phi"),  CCTK_VarIndex("ML_BSSN_MP_O8::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::gt11"),  CCTK_VarIndex("ML_BSSN_MP_O8::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::gt12"),  CCTK_VarIndex("ML_BSSN_MP_O8::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::gt13"),  CCTK_VarIndex("ML_BSSN_MP_O8::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::gt22"),  CCTK_VarIndex("ML_BSSN_MP_O8::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::gt23"),  CCTK_VarIndex("ML_BSSN_MP_O8::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::gt33"),  CCTK_VarIndex("ML_BSSN_MP_O8::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::beta1"),  CCTK_VarIndex("ML_BSSN_MP_O8::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::beta2"),  CCTK_VarIndex("ML_BSSN_MP_O8::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::beta3"),  CCTK_VarIndex("ML_BSSN_MP_O8::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_MP_O8::trK"),  CCTK_VarIndex("ML_BSSN_MP_O8::trKrhs"));
  
  /* Register all the evolved Array functions with MoL */
  return;
}
