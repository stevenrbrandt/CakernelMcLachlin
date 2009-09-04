/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ML_BSSNUp_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::At11"),  CCTK_VarIndex("ML_BSSNUp::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::At12"),  CCTK_VarIndex("ML_BSSNUp::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::At13"),  CCTK_VarIndex("ML_BSSNUp::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::At22"),  CCTK_VarIndex("ML_BSSNUp::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::At23"),  CCTK_VarIndex("ML_BSSNUp::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::At33"),  CCTK_VarIndex("ML_BSSNUp::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::A"),  CCTK_VarIndex("ML_BSSNUp::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::B1"),  CCTK_VarIndex("ML_BSSNUp::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::B2"),  CCTK_VarIndex("ML_BSSNUp::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::B3"),  CCTK_VarIndex("ML_BSSNUp::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::Xt1"),  CCTK_VarIndex("ML_BSSNUp::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::Xt2"),  CCTK_VarIndex("ML_BSSNUp::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::Xt3"),  CCTK_VarIndex("ML_BSSNUp::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::alpha"),  CCTK_VarIndex("ML_BSSNUp::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::phi"),  CCTK_VarIndex("ML_BSSNUp::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::gt11"),  CCTK_VarIndex("ML_BSSNUp::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::gt12"),  CCTK_VarIndex("ML_BSSNUp::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::gt13"),  CCTK_VarIndex("ML_BSSNUp::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::gt22"),  CCTK_VarIndex("ML_BSSNUp::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::gt23"),  CCTK_VarIndex("ML_BSSNUp::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::gt33"),  CCTK_VarIndex("ML_BSSNUp::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::beta1"),  CCTK_VarIndex("ML_BSSNUp::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::beta2"),  CCTK_VarIndex("ML_BSSNUp::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::beta3"),  CCTK_VarIndex("ML_BSSNUp::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSNUp::trK"),  CCTK_VarIndex("ML_BSSNUp::trKrhs"));
  return;
}
