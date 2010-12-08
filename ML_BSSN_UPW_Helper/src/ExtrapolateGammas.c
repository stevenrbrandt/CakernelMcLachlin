#include <cctk.h>
#include <cctk_Arguments.h>

static void
extrap (cGH const * restrict cctkGH,
        CCTK_REAL * restrict var);

void
ML_BSSN_UPW_ExtrapolateGammas (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  extrap (cctkGH, Xt1);
  extrap (cctkGH, Xt2);
  extrap (cctkGH, Xt3);
  
  extrap (cctkGH, A);
  
  extrap (cctkGH, B1);
  extrap (cctkGH, B2);
  extrap (cctkGH, B3);
}

static void
extrap (cGH const * restrict const cctkGH,
        CCTK_REAL * restrict const var)
{
  ExtrapolateGammas (cctkGH, var);
}
