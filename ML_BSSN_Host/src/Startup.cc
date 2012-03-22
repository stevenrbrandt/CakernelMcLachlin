/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_BSSN_Host_Startup(void)
{
  const char * banner = "ML_BSSN_Host";
  CCTK_RegisterBanner(banner);
  return 0;
}
