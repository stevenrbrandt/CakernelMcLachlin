/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_ADM_Startup(void)
{
  const char * banner = "ML_ADM";
  CCTK_RegisterBanner(banner);
  return 0;
}
