/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_BSSN_O2_Startup(void)
{
  const char * banner = "ML_BSSN_O2";
  CCTK_RegisterBanner(banner);
  return 0;
}
