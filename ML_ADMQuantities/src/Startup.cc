/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_ADMQuantities_Startup(void)
{
  const char * banner = "ML_ADMQuantities";
  CCTK_RegisterBanner(banner);
  return 0;
}
