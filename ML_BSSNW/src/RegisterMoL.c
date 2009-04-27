/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (64-bit) (May 21, 2008) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_BSSNW_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::At11"),  CCTK_VarIndex("ML_BSSNW::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::At12"),  CCTK_VarIndex("ML_BSSNW::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::At13"),  CCTK_VarIndex("ML_BSSNW::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::At22"),  CCTK_VarIndex("ML_BSSNW::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::At23"),  CCTK_VarIndex("ML_BSSNW::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::At33"),  CCTK_VarIndex("ML_BSSNW::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::A"),  CCTK_VarIndex("ML_BSSNW::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::B1"),  CCTK_VarIndex("ML_BSSNW::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::B2"),  CCTK_VarIndex("ML_BSSNW::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::B3"),  CCTK_VarIndex("ML_BSSNW::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::Xt1"),  CCTK_VarIndex("ML_BSSNW::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::Xt2"),  CCTK_VarIndex("ML_BSSNW::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::Xt3"),  CCTK_VarIndex("ML_BSSNW::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::alpha"),  CCTK_VarIndex("ML_BSSNW::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::gt11"),  CCTK_VarIndex("ML_BSSNW::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::gt12"),  CCTK_VarIndex("ML_BSSNW::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::gt13"),  CCTK_VarIndex("ML_BSSNW::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::gt22"),  CCTK_VarIndex("ML_BSSNW::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::gt23"),  CCTK_VarIndex("ML_BSSNW::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::gt33"),  CCTK_VarIndex("ML_BSSNW::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::beta1"),  CCTK_VarIndex("ML_BSSNW::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::beta2"),  CCTK_VarIndex("ML_BSSNW::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::beta3"),  CCTK_VarIndex("ML_BSSNW::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::trK"),  CCTK_VarIndex("ML_BSSNW::trKrhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNW::W"),  CCTK_VarIndex("ML_BSSNW::Wrhs"));
  return;
}
