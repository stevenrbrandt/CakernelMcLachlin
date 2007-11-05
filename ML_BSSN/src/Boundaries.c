/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (32-bit) (April 20, 2007) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"


/* the boundary treatment is split into 3 steps:    */
/* 1. excision                                      */
/* 2. symmetries                                    */
/* 3. "other" boundary conditions, e.g. radiative */

/* to simplify scheduling and testing, the 3 steps  */
/* are currently applied in separate functions      */


void ML_BSSN_CheckBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  return;
}

void ML_BSSN_ApplyBoundConds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT ierr = 0;
  
  if (CCTK_EQUALS(curv_bound, "none"  ) ||
      CCTK_EQUALS(curv_bound, "static") ||
      CCTK_EQUALS(curv_bound, "flat"  ) ||
      CCTK_EQUALS(curv_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::curv", curv_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register curv_bound BC for ML_BSSN::curv!");
  }
  
  if (CCTK_EQUALS(dtlapse_bound, "none"  ) ||
      CCTK_EQUALS(dtlapse_bound, "static") ||
      CCTK_EQUALS(dtlapse_bound, "flat"  ) ||
      CCTK_EQUALS(dtlapse_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::dtlapse", dtlapse_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register dtlapse_bound BC for ML_BSSN::dtlapse!");
  }
  
  if (CCTK_EQUALS(dtshift_bound, "none"  ) ||
      CCTK_EQUALS(dtshift_bound, "static") ||
      CCTK_EQUALS(dtshift_bound, "flat"  ) ||
      CCTK_EQUALS(dtshift_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::dtshift", dtshift_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register dtshift_bound BC for ML_BSSN::dtshift!");
  }
  
  if (CCTK_EQUALS(Gamma_bound, "none"  ) ||
      CCTK_EQUALS(Gamma_bound, "static") ||
      CCTK_EQUALS(Gamma_bound, "flat"  ) ||
      CCTK_EQUALS(Gamma_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::Gamma", Gamma_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Gamma_bound BC for ML_BSSN::Gamma!");
  }
  
  if (CCTK_EQUALS(lapse_bound, "none"  ) ||
      CCTK_EQUALS(lapse_bound, "static") ||
      CCTK_EQUALS(lapse_bound, "flat"  ) ||
      CCTK_EQUALS(lapse_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::lapse", lapse_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register lapse_bound BC for ML_BSSN::lapse!");
  }
  
  if (CCTK_EQUALS(log_confac_bound, "none"  ) ||
      CCTK_EQUALS(log_confac_bound, "static") ||
      CCTK_EQUALS(log_confac_bound, "flat"  ) ||
      CCTK_EQUALS(log_confac_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::log_confac", log_confac_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register log_confac_bound BC for ML_BSSN::log_confac!");
  }
  
  if (CCTK_EQUALS(metric_bound, "none"  ) ||
      CCTK_EQUALS(metric_bound, "static") ||
      CCTK_EQUALS(metric_bound, "flat"  ) ||
      CCTK_EQUALS(metric_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::metric", metric_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register metric_bound BC for ML_BSSN::metric!");
  }
  
  if (CCTK_EQUALS(shift_bound, "none"  ) ||
      CCTK_EQUALS(shift_bound, "static") ||
      CCTK_EQUALS(shift_bound, "flat"  ) ||
      CCTK_EQUALS(shift_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::shift", shift_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register shift_bound BC for ML_BSSN::shift!");
  }
  
  if (CCTK_EQUALS(trace_curv_bound, "none"  ) ||
      CCTK_EQUALS(trace_curv_bound, "static") ||
      CCTK_EQUALS(trace_curv_bound, "flat"  ) ||
      CCTK_EQUALS(trace_curv_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::trace_curv", trace_curv_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register trace_curv_bound BC for ML_BSSN::trace_curv!");
  }
  
  if (CCTK_EQUALS(At11_bound, "none"  ) ||
      CCTK_EQUALS(At11_bound, "static") ||
      CCTK_EQUALS(At11_bound, "flat"  ) ||
      CCTK_EQUALS(At11_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::At11", At11_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register At11_bound BC for ML_BSSN::At11!");
  }
  
  if (CCTK_EQUALS(At21_bound, "none"  ) ||
      CCTK_EQUALS(At21_bound, "static") ||
      CCTK_EQUALS(At21_bound, "flat"  ) ||
      CCTK_EQUALS(At21_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::At21", At21_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register At21_bound BC for ML_BSSN::At21!");
  }
  
  if (CCTK_EQUALS(At31_bound, "none"  ) ||
      CCTK_EQUALS(At31_bound, "static") ||
      CCTK_EQUALS(At31_bound, "flat"  ) ||
      CCTK_EQUALS(At31_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::At31", At31_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register At31_bound BC for ML_BSSN::At31!");
  }
  
  if (CCTK_EQUALS(At22_bound, "none"  ) ||
      CCTK_EQUALS(At22_bound, "static") ||
      CCTK_EQUALS(At22_bound, "flat"  ) ||
      CCTK_EQUALS(At22_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::At22", At22_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register At22_bound BC for ML_BSSN::At22!");
  }
  
  if (CCTK_EQUALS(At32_bound, "none"  ) ||
      CCTK_EQUALS(At32_bound, "static") ||
      CCTK_EQUALS(At32_bound, "flat"  ) ||
      CCTK_EQUALS(At32_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::At32", At32_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register At32_bound BC for ML_BSSN::At32!");
  }
  
  if (CCTK_EQUALS(At33_bound, "none"  ) ||
      CCTK_EQUALS(At33_bound, "static") ||
      CCTK_EQUALS(At33_bound, "flat"  ) ||
      CCTK_EQUALS(At33_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::At33", At33_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register At33_bound BC for ML_BSSN::At33!");
  }
  
  if (CCTK_EQUALS(dtalpha_bound, "none"  ) ||
      CCTK_EQUALS(dtalpha_bound, "static") ||
      CCTK_EQUALS(dtalpha_bound, "flat"  ) ||
      CCTK_EQUALS(dtalpha_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::dtalpha", dtalpha_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register dtalpha_bound BC for ML_BSSN::dtalpha!");
  }
  
  if (CCTK_EQUALS(dtbeta1_bound, "none"  ) ||
      CCTK_EQUALS(dtbeta1_bound, "static") ||
      CCTK_EQUALS(dtbeta1_bound, "flat"  ) ||
      CCTK_EQUALS(dtbeta1_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::dtbeta1", dtbeta1_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register dtbeta1_bound BC for ML_BSSN::dtbeta1!");
  }
  
  if (CCTK_EQUALS(dtbeta2_bound, "none"  ) ||
      CCTK_EQUALS(dtbeta2_bound, "static") ||
      CCTK_EQUALS(dtbeta2_bound, "flat"  ) ||
      CCTK_EQUALS(dtbeta2_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::dtbeta2", dtbeta2_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register dtbeta2_bound BC for ML_BSSN::dtbeta2!");
  }
  
  if (CCTK_EQUALS(dtbeta3_bound, "none"  ) ||
      CCTK_EQUALS(dtbeta3_bound, "static") ||
      CCTK_EQUALS(dtbeta3_bound, "flat"  ) ||
      CCTK_EQUALS(dtbeta3_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::dtbeta3", dtbeta3_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register dtbeta3_bound BC for ML_BSSN::dtbeta3!");
  }
  
  if (CCTK_EQUALS(Xt1_bound, "none"  ) ||
      CCTK_EQUALS(Xt1_bound, "static") ||
      CCTK_EQUALS(Xt1_bound, "flat"  ) ||
      CCTK_EQUALS(Xt1_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::Xt1", Xt1_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Xt1_bound BC for ML_BSSN::Xt1!");
  }
  
  if (CCTK_EQUALS(Xt2_bound, "none"  ) ||
      CCTK_EQUALS(Xt2_bound, "static") ||
      CCTK_EQUALS(Xt2_bound, "flat"  ) ||
      CCTK_EQUALS(Xt2_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::Xt2", Xt2_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Xt2_bound BC for ML_BSSN::Xt2!");
  }
  
  if (CCTK_EQUALS(Xt3_bound, "none"  ) ||
      CCTK_EQUALS(Xt3_bound, "static") ||
      CCTK_EQUALS(Xt3_bound, "flat"  ) ||
      CCTK_EQUALS(Xt3_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::Xt3", Xt3_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Xt3_bound BC for ML_BSSN::Xt3!");
  }
  
  if (CCTK_EQUALS(alpha_bound, "none"  ) ||
      CCTK_EQUALS(alpha_bound, "static") ||
      CCTK_EQUALS(alpha_bound, "flat"  ) ||
      CCTK_EQUALS(alpha_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::alpha", alpha_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register alpha_bound BC for ML_BSSN::alpha!");
  }
  
  if (CCTK_EQUALS(phi_bound, "none"  ) ||
      CCTK_EQUALS(phi_bound, "static") ||
      CCTK_EQUALS(phi_bound, "flat"  ) ||
      CCTK_EQUALS(phi_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::phi", phi_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register phi_bound BC for ML_BSSN::phi!");
  }
  
  if (CCTK_EQUALS(gt11_bound, "none"  ) ||
      CCTK_EQUALS(gt11_bound, "static") ||
      CCTK_EQUALS(gt11_bound, "flat"  ) ||
      CCTK_EQUALS(gt11_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::gt11", gt11_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register gt11_bound BC for ML_BSSN::gt11!");
  }
  
  if (CCTK_EQUALS(gt21_bound, "none"  ) ||
      CCTK_EQUALS(gt21_bound, "static") ||
      CCTK_EQUALS(gt21_bound, "flat"  ) ||
      CCTK_EQUALS(gt21_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::gt21", gt21_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register gt21_bound BC for ML_BSSN::gt21!");
  }
  
  if (CCTK_EQUALS(gt31_bound, "none"  ) ||
      CCTK_EQUALS(gt31_bound, "static") ||
      CCTK_EQUALS(gt31_bound, "flat"  ) ||
      CCTK_EQUALS(gt31_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::gt31", gt31_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register gt31_bound BC for ML_BSSN::gt31!");
  }
  
  if (CCTK_EQUALS(gt22_bound, "none"  ) ||
      CCTK_EQUALS(gt22_bound, "static") ||
      CCTK_EQUALS(gt22_bound, "flat"  ) ||
      CCTK_EQUALS(gt22_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::gt22", gt22_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register gt22_bound BC for ML_BSSN::gt22!");
  }
  
  if (CCTK_EQUALS(gt32_bound, "none"  ) ||
      CCTK_EQUALS(gt32_bound, "static") ||
      CCTK_EQUALS(gt32_bound, "flat"  ) ||
      CCTK_EQUALS(gt32_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::gt32", gt32_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register gt32_bound BC for ML_BSSN::gt32!");
  }
  
  if (CCTK_EQUALS(gt33_bound, "none"  ) ||
      CCTK_EQUALS(gt33_bound, "static") ||
      CCTK_EQUALS(gt33_bound, "flat"  ) ||
      CCTK_EQUALS(gt33_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::gt33", gt33_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register gt33_bound BC for ML_BSSN::gt33!");
  }
  
  if (CCTK_EQUALS(beta1_bound, "none"  ) ||
      CCTK_EQUALS(beta1_bound, "static") ||
      CCTK_EQUALS(beta1_bound, "flat"  ) ||
      CCTK_EQUALS(beta1_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::beta1", beta1_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register beta1_bound BC for ML_BSSN::beta1!");
  }
  
  if (CCTK_EQUALS(beta2_bound, "none"  ) ||
      CCTK_EQUALS(beta2_bound, "static") ||
      CCTK_EQUALS(beta2_bound, "flat"  ) ||
      CCTK_EQUALS(beta2_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::beta2", beta2_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register beta2_bound BC for ML_BSSN::beta2!");
  }
  
  if (CCTK_EQUALS(beta3_bound, "none"  ) ||
      CCTK_EQUALS(beta3_bound, "static") ||
      CCTK_EQUALS(beta3_bound, "flat"  ) ||
      CCTK_EQUALS(beta3_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::beta3", beta3_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register beta3_bound BC for ML_BSSN::beta3!");
  }
  
  if (CCTK_EQUALS(trK_bound, "none"  ) ||
      CCTK_EQUALS(trK_bound, "static") ||
      CCTK_EQUALS(trK_bound, "flat"  ) ||
      CCTK_EQUALS(trK_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_BSSN::trK", trK_bound);
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register trK_bound BC for ML_BSSN::trK!");
  }
  
  if (CCTK_EQUALS(curv_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_curv_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_curv_bound , curv_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_curv_bound ,curv_bound_speed, "SPEED") < 0)
       CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_curv_bound, 
                      "ML_BSSN::curv", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::curv!");
  
  }
  
  if (CCTK_EQUALS(dtlapse_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_dtlapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtlapse_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtlapse_bound , dtlapse_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dtlapse_bound ,dtlapse_bound_speed, "SPEED") < 0)
       CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtlapse_bound, 
                      "ML_BSSN::dtlapse", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::dtlapse!");
  
  }
  
  if (CCTK_EQUALS(dtshift_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_dtshift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtshift_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtshift_bound , dtshift_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dtshift_bound ,dtshift_bound_speed, "SPEED") < 0)
       CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtshift_bound, 
                      "ML_BSSN::dtshift", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::dtshift!");
  
  }
  
  if (CCTK_EQUALS(Gamma_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_Gamma_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Gamma_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_Gamma_bound , Gamma_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Gamma_bound ,Gamma_bound_speed, "SPEED") < 0)
       CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Gamma_bound, 
                      "ML_BSSN::Gamma", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::Gamma!");
  
  }
  
  if (CCTK_EQUALS(lapse_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_lapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_lapse_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_lapse_bound , lapse_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_lapse_bound ,lapse_bound_speed, "SPEED") < 0)
       CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_lapse_bound, 
                      "ML_BSSN::lapse", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::lapse!");
  
  }
  
  if (CCTK_EQUALS(log_confac_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_log_confac_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_log_confac_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_log_confac_bound , log_confac_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_log_confac_bound ,log_confac_bound_speed, "SPEED") < 0)
       CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_log_confac_bound, 
                      "ML_BSSN::log_confac", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::log_confac!");
  
  }
  
  if (CCTK_EQUALS(metric_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_metric_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_metric_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_metric_bound , metric_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_metric_bound ,metric_bound_speed, "SPEED") < 0)
       CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_metric_bound, 
                      "ML_BSSN::metric", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::metric!");
  
  }
  
  if (CCTK_EQUALS(shift_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_shift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_shift_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_shift_bound , shift_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_shift_bound ,shift_bound_speed, "SPEED") < 0)
       CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_shift_bound, 
                      "ML_BSSN::shift", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::shift!");
  
  }
  
  if (CCTK_EQUALS(trace_curv_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_trace_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_trace_curv_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_trace_curv_bound , trace_curv_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_trace_curv_bound ,trace_curv_bound_speed, "SPEED") < 0)
       CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_trace_curv_bound, 
                      "ML_BSSN::trace_curv", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::trace_curv!");
  
  }
  
  if (CCTK_EQUALS(At11_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_At11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At11_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At11_bound , At11_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At11_bound ,At11_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At11_bound, 
                      "ML_BSSN::At11", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::At11!");
  
  }
  
  if (CCTK_EQUALS(At21_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_At21_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At21_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At21_bound , At21_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At21_bound ,At21_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At21_bound, 
                      "ML_BSSN::At21", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::At21!");
  
  }
  
  if (CCTK_EQUALS(At31_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_At31_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At31_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At31_bound , At31_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At31_bound ,At31_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At31_bound, 
                      "ML_BSSN::At31", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::At31!");
  
  }
  
  if (CCTK_EQUALS(At22_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_At22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At22_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At22_bound , At22_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At22_bound ,At22_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At22_bound, 
                      "ML_BSSN::At22", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::At22!");
  
  }
  
  if (CCTK_EQUALS(At32_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_At32_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At32_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At32_bound , At32_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At32_bound ,At32_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At32_bound, 
                      "ML_BSSN::At32", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::At32!");
  
  }
  
  if (CCTK_EQUALS(At33_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_At33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At33_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At33_bound , At33_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At33_bound ,At33_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At33_bound, 
                      "ML_BSSN::At33", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::At33!");
  
  }
  
  if (CCTK_EQUALS(dtalpha_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_dtalpha_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtalpha_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtalpha_bound , dtalpha_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dtalpha_bound ,dtalpha_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtalpha_bound, 
                      "ML_BSSN::dtalpha", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::dtalpha!");
  
  }
  
  if (CCTK_EQUALS(dtbeta1_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_dtbeta1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtbeta1_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtbeta1_bound , dtbeta1_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dtbeta1_bound ,dtbeta1_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtbeta1_bound, 
                      "ML_BSSN::dtbeta1", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::dtbeta1!");
  
  }
  
  if (CCTK_EQUALS(dtbeta2_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_dtbeta2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtbeta2_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtbeta2_bound , dtbeta2_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dtbeta2_bound ,dtbeta2_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtbeta2_bound, 
                      "ML_BSSN::dtbeta2", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::dtbeta2!");
  
  }
  
  if (CCTK_EQUALS(dtbeta3_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_dtbeta3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtbeta3_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtbeta3_bound , dtbeta3_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_dtbeta3_bound ,dtbeta3_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtbeta3_bound, 
                      "ML_BSSN::dtbeta3", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::dtbeta3!");
  
  }
  
  if (CCTK_EQUALS(Xt1_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_Xt1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt1_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_Xt1_bound , Xt1_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Xt1_bound ,Xt1_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt1_bound, 
                      "ML_BSSN::Xt1", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::Xt1!");
  
  }
  
  if (CCTK_EQUALS(Xt2_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_Xt2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt2_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_Xt2_bound , Xt2_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Xt2_bound ,Xt2_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt2_bound, 
                      "ML_BSSN::Xt2", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::Xt2!");
  
  }
  
  if (CCTK_EQUALS(Xt3_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_Xt3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt3_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_Xt3_bound , Xt3_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Xt3_bound ,Xt3_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt3_bound, 
                      "ML_BSSN::Xt3", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::Xt3!");
  
  }
  
  if (CCTK_EQUALS(alpha_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_alpha_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_alpha_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_alpha_bound , alpha_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_alpha_bound ,alpha_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_alpha_bound, 
                      "ML_BSSN::alpha", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::alpha!");
  
  }
  
  if (CCTK_EQUALS(phi_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_phi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phi_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_phi_bound , phi_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_phi_bound ,phi_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phi_bound, 
                      "ML_BSSN::phi", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::phi!");
  
  }
  
  if (CCTK_EQUALS(gt11_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_gt11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt11_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt11_bound , gt11_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt11_bound ,gt11_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt11_bound, 
                      "ML_BSSN::gt11", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::gt11!");
  
  }
  
  if (CCTK_EQUALS(gt21_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_gt21_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt21_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt21_bound , gt21_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt21_bound ,gt21_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt21_bound, 
                      "ML_BSSN::gt21", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::gt21!");
  
  }
  
  if (CCTK_EQUALS(gt31_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_gt31_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt31_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt31_bound , gt31_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt31_bound ,gt31_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt31_bound, 
                      "ML_BSSN::gt31", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::gt31!");
  
  }
  
  if (CCTK_EQUALS(gt22_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_gt22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt22_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt22_bound , gt22_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt22_bound ,gt22_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt22_bound, 
                      "ML_BSSN::gt22", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::gt22!");
  
  }
  
  if (CCTK_EQUALS(gt32_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_gt32_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt32_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt32_bound , gt32_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt32_bound ,gt32_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt32_bound, 
                      "ML_BSSN::gt32", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::gt32!");
  
  }
  
  if (CCTK_EQUALS(gt33_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_gt33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt33_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt33_bound , gt33_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt33_bound ,gt33_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt33_bound, 
                      "ML_BSSN::gt33", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::gt33!");
  
  }
  
  if (CCTK_EQUALS(beta1_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_beta1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta1_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_beta1_bound , beta1_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta1_bound ,beta1_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta1_bound, 
                      "ML_BSSN::beta1", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::beta1!");
  
  }
  
  if (CCTK_EQUALS(beta2_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_beta2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta2_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_beta2_bound , beta2_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta2_bound ,beta2_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta2_bound, 
                      "ML_BSSN::beta2", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::beta2!");
  
  }
  
  if (CCTK_EQUALS(beta3_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_beta3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta3_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_beta3_bound , beta3_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta3_bound ,beta3_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta3_bound, 
                      "ML_BSSN::beta3", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::beta3!");
  
  }
  
  if (CCTK_EQUALS(trK_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    CCTK_INT handle_trK_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_trK_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_trK_bound , trK_bound_limit, "LIMIT") < 0)
       CCTK_WARN(-1, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_trK_bound ,trK_bound_speed, "SPEED") < 0)
        CCTK_WARN(-1, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_trK_bound, 
                      "ML_BSSN::trK", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Radiation BC for ML_BSSN::trK!");
  
  }
  
  if (CCTK_EQUALS(curv_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_curv_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_curv_bound ,curv_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_curv_bound, 
                      "ML_BSSN::curv", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Scalar BC for ML_BSSN::curv!");
  
  }
  
  if (CCTK_EQUALS(dtlapse_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_dtlapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtlapse_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtlapse_bound ,dtlapse_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtlapse_bound, 
                      "ML_BSSN::dtlapse", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Scalar BC for ML_BSSN::dtlapse!");
  
  }
  
  if (CCTK_EQUALS(dtshift_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_dtshift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtshift_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtshift_bound ,dtshift_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtshift_bound, 
                      "ML_BSSN::dtshift", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Scalar BC for ML_BSSN::dtshift!");
  
  }
  
  if (CCTK_EQUALS(Gamma_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_Gamma_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Gamma_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_Gamma_bound ,Gamma_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Gamma_bound, 
                      "ML_BSSN::Gamma", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Scalar BC for ML_BSSN::Gamma!");
  
  }
  
  if (CCTK_EQUALS(lapse_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_lapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_lapse_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_lapse_bound ,lapse_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_lapse_bound, 
                      "ML_BSSN::lapse", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Scalar BC for ML_BSSN::lapse!");
  
  }
  
  if (CCTK_EQUALS(log_confac_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_log_confac_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_log_confac_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_log_confac_bound ,log_confac_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_log_confac_bound, 
                      "ML_BSSN::log_confac", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Scalar BC for ML_BSSN::log_confac!");
  
  }
  
  if (CCTK_EQUALS(metric_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_metric_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_metric_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_metric_bound ,metric_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_metric_bound, 
                      "ML_BSSN::metric", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Scalar BC for ML_BSSN::metric!");
  
  }
  
  if (CCTK_EQUALS(shift_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_shift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_shift_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_shift_bound ,shift_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_shift_bound, 
                      "ML_BSSN::shift", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Scalar BC for ML_BSSN::shift!");
  
  }
  
  if (CCTK_EQUALS(trace_curv_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_trace_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_trace_curv_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_trace_curv_bound ,trace_curv_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_trace_curv_bound, 
                      "ML_BSSN::trace_curv", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register Scalar BC for ML_BSSN::trace_curv!");
  
  }
  
  if (CCTK_EQUALS(At11_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_At11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At11_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At11_bound ,At11_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At11_bound, 
                      "ML_BSSN::At11", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::At11!");
  
  }
  
  if (CCTK_EQUALS(At21_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_At21_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At21_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At21_bound ,At21_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At21_bound, 
                      "ML_BSSN::At21", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::At21!");
  
  }
  
  if (CCTK_EQUALS(At31_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_At31_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At31_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At31_bound ,At31_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At31_bound, 
                      "ML_BSSN::At31", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::At31!");
  
  }
  
  if (CCTK_EQUALS(At22_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_At22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At22_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At22_bound ,At22_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At22_bound, 
                      "ML_BSSN::At22", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::At22!");
  
  }
  
  if (CCTK_EQUALS(At32_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_At32_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At32_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At32_bound ,At32_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At32_bound, 
                      "ML_BSSN::At32", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::At32!");
  
  }
  
  if (CCTK_EQUALS(At33_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_At33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At33_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_At33_bound ,At33_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At33_bound, 
                      "ML_BSSN::At33", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::At33!");
  
  }
  
  if (CCTK_EQUALS(dtalpha_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_dtalpha_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtalpha_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtalpha_bound ,dtalpha_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtalpha_bound, 
                      "ML_BSSN::dtalpha", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::dtalpha!");
  
  }
  
  if (CCTK_EQUALS(dtbeta1_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_dtbeta1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtbeta1_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtbeta1_bound ,dtbeta1_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtbeta1_bound, 
                      "ML_BSSN::dtbeta1", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::dtbeta1!");
  
  }
  
  if (CCTK_EQUALS(dtbeta2_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_dtbeta2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtbeta2_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtbeta2_bound ,dtbeta2_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtbeta2_bound, 
                      "ML_BSSN::dtbeta2", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::dtbeta2!");
  
  }
  
  if (CCTK_EQUALS(dtbeta3_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_dtbeta3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_dtbeta3_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_dtbeta3_bound ,dtbeta3_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_dtbeta3_bound, 
                      "ML_BSSN::dtbeta3", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::dtbeta3!");
  
  }
  
  if (CCTK_EQUALS(Xt1_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_Xt1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt1_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_Xt1_bound ,Xt1_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt1_bound, 
                      "ML_BSSN::Xt1", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::Xt1!");
  
  }
  
  if (CCTK_EQUALS(Xt2_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_Xt2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt2_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_Xt2_bound ,Xt2_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt2_bound, 
                      "ML_BSSN::Xt2", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::Xt2!");
  
  }
  
  if (CCTK_EQUALS(Xt3_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_Xt3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt3_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_Xt3_bound ,Xt3_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt3_bound, 
                      "ML_BSSN::Xt3", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::Xt3!");
  
  }
  
  if (CCTK_EQUALS(alpha_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_alpha_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_alpha_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_alpha_bound ,alpha_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_alpha_bound, 
                      "ML_BSSN::alpha", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::alpha!");
  
  }
  
  if (CCTK_EQUALS(phi_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_phi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phi_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_phi_bound ,phi_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phi_bound, 
                      "ML_BSSN::phi", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::phi!");
  
  }
  
  if (CCTK_EQUALS(gt11_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_gt11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt11_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt11_bound ,gt11_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt11_bound, 
                      "ML_BSSN::gt11", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::gt11!");
  
  }
  
  if (CCTK_EQUALS(gt21_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_gt21_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt21_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt21_bound ,gt21_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt21_bound, 
                      "ML_BSSN::gt21", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::gt21!");
  
  }
  
  if (CCTK_EQUALS(gt31_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_gt31_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt31_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt31_bound ,gt31_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt31_bound, 
                      "ML_BSSN::gt31", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::gt31!");
  
  }
  
  if (CCTK_EQUALS(gt22_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_gt22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt22_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt22_bound ,gt22_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt22_bound, 
                      "ML_BSSN::gt22", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::gt22!");
  
  }
  
  if (CCTK_EQUALS(gt32_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_gt32_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt32_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt32_bound ,gt32_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt32_bound, 
                      "ML_BSSN::gt32", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::gt32!");
  
  }
  
  if (CCTK_EQUALS(gt33_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_gt33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt33_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_gt33_bound ,gt33_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt33_bound, 
                      "ML_BSSN::gt33", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::gt33!");
  
  }
  
  if (CCTK_EQUALS(beta1_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_beta1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta1_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_beta1_bound ,beta1_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta1_bound, 
                      "ML_BSSN::beta1", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::beta1!");
  
  }
  
  if (CCTK_EQUALS(beta2_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_beta2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta2_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_beta2_bound ,beta2_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta2_bound, 
                      "ML_BSSN::beta2", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::beta2!");
  
  }
  
  if (CCTK_EQUALS(beta3_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_beta3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta3_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_beta3_bound ,beta3_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta3_bound, 
                      "ML_BSSN::beta3", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::beta3!");
  
  }
  
  if (CCTK_EQUALS(trK_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    CCTK_INT handle_trK_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_trK_bound < 0) CCTK_WARN(-1, "could not create table!");
    if (Util_TableSetReal(handle_trK_bound ,trK_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(-1, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_trK_bound, 
                      "ML_BSSN::trK", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(-1, "Error in registering Scalar BC for ML_BSSN::trK!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#ML_BSSN::curv_bound       = "skip"
#$bound$#ML_BSSN::curv_bound_speed = 1.0
#$bound$#ML_BSSN::curv_bound_limit = 0.0
#$bound$#ML_BSSN::curv_bound_scalar = 0.0

#$bound$#ML_BSSN::dtlapse_bound       = "skip"
#$bound$#ML_BSSN::dtlapse_bound_speed = 1.0
#$bound$#ML_BSSN::dtlapse_bound_limit = 0.0
#$bound$#ML_BSSN::dtlapse_bound_scalar = 0.0

#$bound$#ML_BSSN::dtshift_bound       = "skip"
#$bound$#ML_BSSN::dtshift_bound_speed = 1.0
#$bound$#ML_BSSN::dtshift_bound_limit = 0.0
#$bound$#ML_BSSN::dtshift_bound_scalar = 0.0

#$bound$#ML_BSSN::Gamma_bound       = "skip"
#$bound$#ML_BSSN::Gamma_bound_speed = 1.0
#$bound$#ML_BSSN::Gamma_bound_limit = 0.0
#$bound$#ML_BSSN::Gamma_bound_scalar = 0.0

#$bound$#ML_BSSN::lapse_bound       = "skip"
#$bound$#ML_BSSN::lapse_bound_speed = 1.0
#$bound$#ML_BSSN::lapse_bound_limit = 0.0
#$bound$#ML_BSSN::lapse_bound_scalar = 0.0

#$bound$#ML_BSSN::log_confac_bound       = "skip"
#$bound$#ML_BSSN::log_confac_bound_speed = 1.0
#$bound$#ML_BSSN::log_confac_bound_limit = 0.0
#$bound$#ML_BSSN::log_confac_bound_scalar = 0.0

#$bound$#ML_BSSN::metric_bound       = "skip"
#$bound$#ML_BSSN::metric_bound_speed = 1.0
#$bound$#ML_BSSN::metric_bound_limit = 0.0
#$bound$#ML_BSSN::metric_bound_scalar = 0.0

#$bound$#ML_BSSN::shift_bound       = "skip"
#$bound$#ML_BSSN::shift_bound_speed = 1.0
#$bound$#ML_BSSN::shift_bound_limit = 0.0
#$bound$#ML_BSSN::shift_bound_scalar = 0.0

#$bound$#ML_BSSN::trace_curv_bound       = "skip"
#$bound$#ML_BSSN::trace_curv_bound_speed = 1.0
#$bound$#ML_BSSN::trace_curv_bound_limit = 0.0
#$bound$#ML_BSSN::trace_curv_bound_scalar = 0.0

#$bound$#ML_BSSN::At11_bound       = "skip"
#$bound$#ML_BSSN::At11_bound_speed = 1.0
#$bound$#ML_BSSN::At11_bound_limit = 0.0
#$bound$#ML_BSSN::At11_bound_scalar = 0.0

#$bound$#ML_BSSN::At21_bound       = "skip"
#$bound$#ML_BSSN::At21_bound_speed = 1.0
#$bound$#ML_BSSN::At21_bound_limit = 0.0
#$bound$#ML_BSSN::At21_bound_scalar = 0.0

#$bound$#ML_BSSN::At31_bound       = "skip"
#$bound$#ML_BSSN::At31_bound_speed = 1.0
#$bound$#ML_BSSN::At31_bound_limit = 0.0
#$bound$#ML_BSSN::At31_bound_scalar = 0.0

#$bound$#ML_BSSN::At22_bound       = "skip"
#$bound$#ML_BSSN::At22_bound_speed = 1.0
#$bound$#ML_BSSN::At22_bound_limit = 0.0
#$bound$#ML_BSSN::At22_bound_scalar = 0.0

#$bound$#ML_BSSN::At32_bound       = "skip"
#$bound$#ML_BSSN::At32_bound_speed = 1.0
#$bound$#ML_BSSN::At32_bound_limit = 0.0
#$bound$#ML_BSSN::At32_bound_scalar = 0.0

#$bound$#ML_BSSN::At33_bound       = "skip"
#$bound$#ML_BSSN::At33_bound_speed = 1.0
#$bound$#ML_BSSN::At33_bound_limit = 0.0
#$bound$#ML_BSSN::At33_bound_scalar = 0.0

#$bound$#ML_BSSN::dtalpha_bound       = "skip"
#$bound$#ML_BSSN::dtalpha_bound_speed = 1.0
#$bound$#ML_BSSN::dtalpha_bound_limit = 0.0
#$bound$#ML_BSSN::dtalpha_bound_scalar = 0.0

#$bound$#ML_BSSN::dtbeta1_bound       = "skip"
#$bound$#ML_BSSN::dtbeta1_bound_speed = 1.0
#$bound$#ML_BSSN::dtbeta1_bound_limit = 0.0
#$bound$#ML_BSSN::dtbeta1_bound_scalar = 0.0

#$bound$#ML_BSSN::dtbeta2_bound       = "skip"
#$bound$#ML_BSSN::dtbeta2_bound_speed = 1.0
#$bound$#ML_BSSN::dtbeta2_bound_limit = 0.0
#$bound$#ML_BSSN::dtbeta2_bound_scalar = 0.0

#$bound$#ML_BSSN::dtbeta3_bound       = "skip"
#$bound$#ML_BSSN::dtbeta3_bound_speed = 1.0
#$bound$#ML_BSSN::dtbeta3_bound_limit = 0.0
#$bound$#ML_BSSN::dtbeta3_bound_scalar = 0.0

#$bound$#ML_BSSN::Xt1_bound       = "skip"
#$bound$#ML_BSSN::Xt1_bound_speed = 1.0
#$bound$#ML_BSSN::Xt1_bound_limit = 0.0
#$bound$#ML_BSSN::Xt1_bound_scalar = 0.0

#$bound$#ML_BSSN::Xt2_bound       = "skip"
#$bound$#ML_BSSN::Xt2_bound_speed = 1.0
#$bound$#ML_BSSN::Xt2_bound_limit = 0.0
#$bound$#ML_BSSN::Xt2_bound_scalar = 0.0

#$bound$#ML_BSSN::Xt3_bound       = "skip"
#$bound$#ML_BSSN::Xt3_bound_speed = 1.0
#$bound$#ML_BSSN::Xt3_bound_limit = 0.0
#$bound$#ML_BSSN::Xt3_bound_scalar = 0.0

#$bound$#ML_BSSN::alpha_bound       = "skip"
#$bound$#ML_BSSN::alpha_bound_speed = 1.0
#$bound$#ML_BSSN::alpha_bound_limit = 0.0
#$bound$#ML_BSSN::alpha_bound_scalar = 0.0

#$bound$#ML_BSSN::phi_bound       = "skip"
#$bound$#ML_BSSN::phi_bound_speed = 1.0
#$bound$#ML_BSSN::phi_bound_limit = 0.0
#$bound$#ML_BSSN::phi_bound_scalar = 0.0

#$bound$#ML_BSSN::gt11_bound       = "skip"
#$bound$#ML_BSSN::gt11_bound_speed = 1.0
#$bound$#ML_BSSN::gt11_bound_limit = 0.0
#$bound$#ML_BSSN::gt11_bound_scalar = 0.0

#$bound$#ML_BSSN::gt21_bound       = "skip"
#$bound$#ML_BSSN::gt21_bound_speed = 1.0
#$bound$#ML_BSSN::gt21_bound_limit = 0.0
#$bound$#ML_BSSN::gt21_bound_scalar = 0.0

#$bound$#ML_BSSN::gt31_bound       = "skip"
#$bound$#ML_BSSN::gt31_bound_speed = 1.0
#$bound$#ML_BSSN::gt31_bound_limit = 0.0
#$bound$#ML_BSSN::gt31_bound_scalar = 0.0

#$bound$#ML_BSSN::gt22_bound       = "skip"
#$bound$#ML_BSSN::gt22_bound_speed = 1.0
#$bound$#ML_BSSN::gt22_bound_limit = 0.0
#$bound$#ML_BSSN::gt22_bound_scalar = 0.0

#$bound$#ML_BSSN::gt32_bound       = "skip"
#$bound$#ML_BSSN::gt32_bound_speed = 1.0
#$bound$#ML_BSSN::gt32_bound_limit = 0.0
#$bound$#ML_BSSN::gt32_bound_scalar = 0.0

#$bound$#ML_BSSN::gt33_bound       = "skip"
#$bound$#ML_BSSN::gt33_bound_speed = 1.0
#$bound$#ML_BSSN::gt33_bound_limit = 0.0
#$bound$#ML_BSSN::gt33_bound_scalar = 0.0

#$bound$#ML_BSSN::beta1_bound       = "skip"
#$bound$#ML_BSSN::beta1_bound_speed = 1.0
#$bound$#ML_BSSN::beta1_bound_limit = 0.0
#$bound$#ML_BSSN::beta1_bound_scalar = 0.0

#$bound$#ML_BSSN::beta2_bound       = "skip"
#$bound$#ML_BSSN::beta2_bound_speed = 1.0
#$bound$#ML_BSSN::beta2_bound_limit = 0.0
#$bound$#ML_BSSN::beta2_bound_scalar = 0.0

#$bound$#ML_BSSN::beta3_bound       = "skip"
#$bound$#ML_BSSN::beta3_bound_speed = 1.0
#$bound$#ML_BSSN::beta3_bound_limit = 0.0
#$bound$#ML_BSSN::beta3_bound_scalar = 0.0

#$bound$#ML_BSSN::trK_bound       = "skip"
#$bound$#ML_BSSN::trK_bound_speed = 1.0
#$bound$#ML_BSSN::trK_bound_limit = 0.0
#$bound$#ML_BSSN::trK_bound_scalar = 0.0

*/

