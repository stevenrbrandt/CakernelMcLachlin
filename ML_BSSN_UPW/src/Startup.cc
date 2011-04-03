/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_BSSN_UPW_Startup(void)
{
  const char * banner = "ML_BSSN_UPW";
  CCTK_RegisterBanner(banner);
  return 0;
}
