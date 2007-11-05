/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (32-bit) (April 20, 2007) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_BSSN_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At11"),  CCTK_VarIndex("ML_BSSN::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At21"),  CCTK_VarIndex("ML_BSSN::At21rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At31"),  CCTK_VarIndex("ML_BSSN::At31rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At22"),  CCTK_VarIndex("ML_BSSN::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At32"),  CCTK_VarIndex("ML_BSSN::At32rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At33"),  CCTK_VarIndex("ML_BSSN::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::dtalpha"),  CCTK_VarIndex("ML_BSSN::dtalpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::dtbeta1"),  CCTK_VarIndex("ML_BSSN::dtbeta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::dtbeta2"),  CCTK_VarIndex("ML_BSSN::dtbeta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::dtbeta3"),  CCTK_VarIndex("ML_BSSN::dtbeta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::Xt1"),  CCTK_VarIndex("ML_BSSN::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::Xt2"),  CCTK_VarIndex("ML_BSSN::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::Xt3"),  CCTK_VarIndex("ML_BSSN::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::alpha"),  CCTK_VarIndex("ML_BSSN::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::phi"),  CCTK_VarIndex("ML_BSSN::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt11"),  CCTK_VarIndex("ML_BSSN::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt21"),  CCTK_VarIndex("ML_BSSN::gt21rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt31"),  CCTK_VarIndex("ML_BSSN::gt31rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt22"),  CCTK_VarIndex("ML_BSSN::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt32"),  CCTK_VarIndex("ML_BSSN::gt32rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt33"),  CCTK_VarIndex("ML_BSSN::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::beta1"),  CCTK_VarIndex("ML_BSSN::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::beta2"),  CCTK_VarIndex("ML_BSSN::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::beta3"),  CCTK_VarIndex("ML_BSSN::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::trK"),  CCTK_VarIndex("ML_BSSN::trKrhs"));
  return;
}
