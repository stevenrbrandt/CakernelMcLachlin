/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_Host_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::At11"),  CCTK_VarIndex("ML_BSSN_Host::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::At12"),  CCTK_VarIndex("ML_BSSN_Host::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::At13"),  CCTK_VarIndex("ML_BSSN_Host::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::At22"),  CCTK_VarIndex("ML_BSSN_Host::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::At23"),  CCTK_VarIndex("ML_BSSN_Host::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::At33"),  CCTK_VarIndex("ML_BSSN_Host::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::A"),  CCTK_VarIndex("ML_BSSN_Host::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::B1"),  CCTK_VarIndex("ML_BSSN_Host::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::B2"),  CCTK_VarIndex("ML_BSSN_Host::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::B3"),  CCTK_VarIndex("ML_BSSN_Host::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::Xt1"),  CCTK_VarIndex("ML_BSSN_Host::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::Xt2"),  CCTK_VarIndex("ML_BSSN_Host::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::Xt3"),  CCTK_VarIndex("ML_BSSN_Host::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::alpha"),  CCTK_VarIndex("ML_BSSN_Host::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::phi"),  CCTK_VarIndex("ML_BSSN_Host::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::gt11"),  CCTK_VarIndex("ML_BSSN_Host::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::gt12"),  CCTK_VarIndex("ML_BSSN_Host::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::gt13"),  CCTK_VarIndex("ML_BSSN_Host::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::gt22"),  CCTK_VarIndex("ML_BSSN_Host::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::gt23"),  CCTK_VarIndex("ML_BSSN_Host::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::gt33"),  CCTK_VarIndex("ML_BSSN_Host::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::beta1"),  CCTK_VarIndex("ML_BSSN_Host::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::beta2"),  CCTK_VarIndex("ML_BSSN_Host::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::beta3"),  CCTK_VarIndex("ML_BSSN_Host::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_Host::trK"),  CCTK_VarIndex("ML_BSSN_Host::trKrhs"));
  
  /* Register all the evolved Array functions with MoL */
  return;
}
