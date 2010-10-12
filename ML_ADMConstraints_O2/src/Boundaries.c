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


void ML_ADMConstraints_O2_CheckBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

void ML_ADMConstraints_O2_SelectBoundConds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC (cctkGH, CCTK_ALL_FACES, boundary_width, -1, "ML_ADMConstraints_O2::ML_Ham", "scalar");
  if (ierr<0) CCTK_WARN (CCTK_WARN_ABORT, "Failed to select boundary condition for ML_ADMConstraints_O2::ML_Ham");

  ierr = Boundary_SelectGroupForBC (cctkGH, CCTK_ALL_FACES, boundary_width, -1, "ML_ADMConstraints_O2::ML_mom", "scalar");
  if (ierr<0) CCTK_WARN (CCTK_WARN_ABORT, "Failed to select boundary condition for ML_ADMConstraints_O2::ML_mom");
}



/* template for entries in parameter file:
*/

