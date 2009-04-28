#include <stdlib.h>
#include <string.h>

#include <cctk.h>
#include <cctk_Arguments.h>

static void
copy (cGH const * restrict cctkGH,
      CCTK_REAL * restrict dst, CCTK_REAL const * restrict src);

void
ML_BSSNW_CopyADMBase (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  copy (cctkGH, gxx, gxx_p);
  copy (cctkGH, gxy, gxx_p);
  copy (cctkGH, gxz, gxx_p);
  copy (cctkGH, gyy, gxx_p);
  copy (cctkGH, gyz, gxx_p);
  copy (cctkGH, gzz, gxx_p);
  
  copy (cctkGH, kxx, gxx_p);
  copy (cctkGH, kxy, gxx_p);
  copy (cctkGH, kxz, gxx_p);
  copy (cctkGH, kyy, gxx_p);
  copy (cctkGH, kyz, gxx_p);
  copy (cctkGH, kzz, gxx_p);
  
  copy (cctkGH, alp, alp_p);
  
  copy (cctkGH, betax, betax_p);
  copy (cctkGH, betay, betay_p);
  copy (cctkGH, betaz, betaz_p);
  
  copy (cctkGH, dtalp, dtalp_p);
  
  copy (cctkGH, dtbetax, dtbetax_p);
  copy (cctkGH, dtbetay, dtbetay_p);
  copy (cctkGH, dtbetaz, dtbetaz_p);
}

static void
copy (cGH const * restrict const cctkGH,
      CCTK_REAL * restrict const dst, CCTK_REAL const * restrict const src)
{
  size_t const npoints =
    (size_t) cctkGH->cctk_lsh[0] * cctkGH->cctk_lsh[1] * cctkGH->cctk_lsh[2];
  memcpy (dst, src, npoints * sizeof *dst);
}
