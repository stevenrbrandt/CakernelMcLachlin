/*  File produced by Kranc */

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


void ML_hydro_CheckBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

void ML_hydro_ApplyBoundConds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  
  if (CCTK_EQUALS(ene_group_bound, "none"  ) ||
      CCTK_EQUALS(ene_group_bound, "static") ||
      CCTK_EQUALS(ene_group_bound, "flat"  ) ||
      CCTK_EQUALS(ene_group_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_hydro::ene_group", ene_group_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register ene_group_bound BC for ML_hydro::ene_group!");
  }
  
  if (CCTK_EQUALS(mass_group_bound, "none"  ) ||
      CCTK_EQUALS(mass_group_bound, "static") ||
      CCTK_EQUALS(mass_group_bound, "flat"  ) ||
      CCTK_EQUALS(mass_group_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_hydro::mass_group", mass_group_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register mass_group_bound BC for ML_hydro::mass_group!");
  }
  
  if (CCTK_EQUALS(mom_group_bound, "none"  ) ||
      CCTK_EQUALS(mom_group_bound, "static") ||
      CCTK_EQUALS(mom_group_bound, "flat"  ) ||
      CCTK_EQUALS(mom_group_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_hydro::mom_group", mom_group_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register mom_group_bound BC for ML_hydro::mom_group!");
  }
  
  if (CCTK_EQUALS(ene_bound, "none"  ) ||
      CCTK_EQUALS(ene_bound, "static") ||
      CCTK_EQUALS(ene_bound, "flat"  ) ||
      CCTK_EQUALS(ene_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_hydro::ene", ene_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register ene_bound BC for ML_hydro::ene!");
  }
  
  if (CCTK_EQUALS(mass_bound, "none"  ) ||
      CCTK_EQUALS(mass_bound, "static") ||
      CCTK_EQUALS(mass_bound, "flat"  ) ||
      CCTK_EQUALS(mass_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_hydro::mass", mass_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register mass_bound BC for ML_hydro::mass!");
  }
  
  if (CCTK_EQUALS(mom1_bound, "none"  ) ||
      CCTK_EQUALS(mom1_bound, "static") ||
      CCTK_EQUALS(mom1_bound, "flat"  ) ||
      CCTK_EQUALS(mom1_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_hydro::mom1", mom1_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register mom1_bound BC for ML_hydro::mom1!");
  }
  
  if (CCTK_EQUALS(mom2_bound, "none"  ) ||
      CCTK_EQUALS(mom2_bound, "static") ||
      CCTK_EQUALS(mom2_bound, "flat"  ) ||
      CCTK_EQUALS(mom2_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_hydro::mom2", mom2_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register mom2_bound BC for ML_hydro::mom2!");
  }
  
  if (CCTK_EQUALS(mom3_bound, "none"  ) ||
      CCTK_EQUALS(mom3_bound, "static") ||
      CCTK_EQUALS(mom3_bound, "flat"  ) ||
      CCTK_EQUALS(mom3_bound, "zero"  ) ) 
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                      "ML_hydro::mom3", mom3_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register mom3_bound BC for ML_hydro::mom3!");
  }
  
  if (CCTK_EQUALS(ene_group_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    static CCTK_INT handle_ene_group_bound = -1;
    if (handle_ene_group_bound < 0) handle_ene_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ene_group_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_ene_group_bound , ene_group_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ene_group_bound ,ene_group_bound_speed, "SPEED") < 0)
       CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ene_group_bound, 
                      "ML_hydro::ene_group", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for ML_hydro::ene_group!");
  
  }
  
  if (CCTK_EQUALS(mass_group_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    static CCTK_INT handle_mass_group_bound = -1;
    if (handle_mass_group_bound < 0) handle_mass_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mass_group_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mass_group_bound , mass_group_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_mass_group_bound ,mass_group_bound_speed, "SPEED") < 0)
       CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mass_group_bound, 
                      "ML_hydro::mass_group", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for ML_hydro::mass_group!");
  
  }
  
  if (CCTK_EQUALS(mom_group_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    static CCTK_INT handle_mom_group_bound = -1;
    if (handle_mom_group_bound < 0) handle_mom_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mom_group_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mom_group_bound , mom_group_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_mom_group_bound ,mom_group_bound_speed, "SPEED") < 0)
       CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mom_group_bound, 
                      "ML_hydro::mom_group", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for ML_hydro::mom_group!");
  
  }
  
  if (CCTK_EQUALS(ene_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    static CCTK_INT handle_ene_bound = -1;
    if (handle_ene_bound < 0) handle_ene_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ene_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_ene_bound , ene_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ene_bound ,ene_bound_speed, "SPEED") < 0)
        CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ene_bound, 
                      "ML_hydro::ene", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for ML_hydro::ene!");
  
  }
  
  if (CCTK_EQUALS(mass_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    static CCTK_INT handle_mass_bound = -1;
    if (handle_mass_bound < 0) handle_mass_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mass_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mass_bound , mass_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_mass_bound ,mass_bound_speed, "SPEED") < 0)
        CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mass_bound, 
                      "ML_hydro::mass", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for ML_hydro::mass!");
  
  }
  
  if (CCTK_EQUALS(mom1_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    static CCTK_INT handle_mom1_bound = -1;
    if (handle_mom1_bound < 0) handle_mom1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mom1_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mom1_bound , mom1_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_mom1_bound ,mom1_bound_speed, "SPEED") < 0)
        CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mom1_bound, 
                      "ML_hydro::mom1", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for ML_hydro::mom1!");
  
  }
  
  if (CCTK_EQUALS(mom2_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    static CCTK_INT handle_mom2_bound = -1;
    if (handle_mom2_bound < 0) handle_mom2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mom2_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mom2_bound , mom2_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_mom2_bound ,mom2_bound_speed, "SPEED") < 0)
        CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mom2_bound, 
                      "ML_hydro::mom2", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for ML_hydro::mom2!");
  
  }
  
  if (CCTK_EQUALS(mom3_bound, "radiative"))
  {
   /* apply radiation boundary condition */
    static CCTK_INT handle_mom3_bound = -1;
    if (handle_mom3_bound < 0) handle_mom3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mom3_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mom3_bound , mom3_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_mom3_bound ,mom3_bound_speed, "SPEED") < 0)
        CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mom3_bound, 
                      "ML_hydro::mom3", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for ML_hydro::mom3!");
  
  }
  
  if (CCTK_EQUALS(ene_group_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    static CCTK_INT handle_ene_group_bound = -1;
    if (handle_ene_group_bound < 0) handle_ene_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ene_group_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_ene_group_bound ,ene_group_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ene_group_bound, 
                      "ML_hydro::ene_group", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Scalar BC for ML_hydro::ene_group!");
  
  }
  
  if (CCTK_EQUALS(mass_group_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    static CCTK_INT handle_mass_group_bound = -1;
    if (handle_mass_group_bound < 0) handle_mass_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mass_group_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mass_group_bound ,mass_group_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mass_group_bound, 
                      "ML_hydro::mass_group", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Scalar BC for ML_hydro::mass_group!");
  
  }
  
  if (CCTK_EQUALS(mom_group_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    static CCTK_INT handle_mom_group_bound = -1;
    if (handle_mom_group_bound < 0) handle_mom_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mom_group_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mom_group_bound ,mom_group_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mom_group_bound, 
                      "ML_hydro::mom_group", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Scalar BC for ML_hydro::mom_group!");
  
  }
  
  if (CCTK_EQUALS(ene_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    static CCTK_INT handle_ene_bound = -1;
    if (handle_ene_bound < 0) handle_ene_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ene_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_ene_bound ,ene_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ene_bound, 
                      "ML_hydro::ene", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Error in registering Scalar BC for ML_hydro::ene!");
  
  }
  
  if (CCTK_EQUALS(mass_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    static CCTK_INT handle_mass_bound = -1;
    if (handle_mass_bound < 0) handle_mass_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mass_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mass_bound ,mass_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mass_bound, 
                      "ML_hydro::mass", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Error in registering Scalar BC for ML_hydro::mass!");
  
  }
  
  if (CCTK_EQUALS(mom1_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    static CCTK_INT handle_mom1_bound = -1;
    if (handle_mom1_bound < 0) handle_mom1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mom1_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mom1_bound ,mom1_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mom1_bound, 
                      "ML_hydro::mom1", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Error in registering Scalar BC for ML_hydro::mom1!");
  
  }
  
  if (CCTK_EQUALS(mom2_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    static CCTK_INT handle_mom2_bound = -1;
    if (handle_mom2_bound < 0) handle_mom2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mom2_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mom2_bound ,mom2_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mom2_bound, 
                      "ML_hydro::mom2", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Error in registering Scalar BC for ML_hydro::mom2!");
  
  }
  
  if (CCTK_EQUALS(mom3_bound, "scalar"))
  {
   /* apply scalar boundary condition */
    static CCTK_INT handle_mom3_bound = -1;
    if (handle_mom3_bound < 0) handle_mom3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_mom3_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_mom3_bound ,mom3_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_mom3_bound, 
                      "ML_hydro::mom3", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Error in registering Scalar BC for ML_hydro::mom3!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#ML_hydro::ene_group_bound       = "skip"
#$bound$#ML_hydro::ene_group_bound_speed = 1.0
#$bound$#ML_hydro::ene_group_bound_limit = 0.0
#$bound$#ML_hydro::ene_group_bound_scalar = 0.0

#$bound$#ML_hydro::mass_group_bound       = "skip"
#$bound$#ML_hydro::mass_group_bound_speed = 1.0
#$bound$#ML_hydro::mass_group_bound_limit = 0.0
#$bound$#ML_hydro::mass_group_bound_scalar = 0.0

#$bound$#ML_hydro::mom_group_bound       = "skip"
#$bound$#ML_hydro::mom_group_bound_speed = 1.0
#$bound$#ML_hydro::mom_group_bound_limit = 0.0
#$bound$#ML_hydro::mom_group_bound_scalar = 0.0

#$bound$#ML_hydro::ene_bound       = "skip"
#$bound$#ML_hydro::ene_bound_speed = 1.0
#$bound$#ML_hydro::ene_bound_limit = 0.0
#$bound$#ML_hydro::ene_bound_scalar = 0.0

#$bound$#ML_hydro::mass_bound       = "skip"
#$bound$#ML_hydro::mass_bound_speed = 1.0
#$bound$#ML_hydro::mass_bound_limit = 0.0
#$bound$#ML_hydro::mass_bound_scalar = 0.0

#$bound$#ML_hydro::mom1_bound       = "skip"
#$bound$#ML_hydro::mom1_bound_speed = 1.0
#$bound$#ML_hydro::mom1_bound_limit = 0.0
#$bound$#ML_hydro::mom1_bound_scalar = 0.0

#$bound$#ML_hydro::mom2_bound       = "skip"
#$bound$#ML_hydro::mom2_bound_speed = 1.0
#$bound$#ML_hydro::mom2_bound_limit = 0.0
#$bound$#ML_hydro::mom2_bound_scalar = 0.0

#$bound$#ML_hydro::mom3_bound       = "skip"
#$bound$#ML_hydro::mom3_bound_speed = 1.0
#$bound$#ML_hydro::mom3_bound_limit = 0.0
#$bound$#ML_hydro::mom3_bound_scalar = 0.0

*/

