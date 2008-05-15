#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

static void
select_bcs (cGH const * restrict cctkGH, char const * restrict gn);

void
ML_BSSN_MP_SelectBCsADMBase (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  select_bcs (cctkGH, "ADMBase::metric");
  select_bcs (cctkGH, "ADMBase::curv");
  select_bcs (cctkGH, "ADMBase::lapse");
  select_bcs (cctkGH, "ADMBase::shift");
  select_bcs (cctkGH, "ADMBase::dtlapse");
  select_bcs (cctkGH, "ADMBase::dtshift");
}

static void
select_bcs (cGH const * restrict const cctkGH, char const * restrict const gn)
{
  DECLARE_CCTK_PARAMETERS;
  
  Boundary_SelectGroupForBC
    (cctkGH, CCTK_ALL_FACES, boundary_width, -1, gn, "none");
}
