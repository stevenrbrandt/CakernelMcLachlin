/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_WaveToyFO_Startup(void)
{
  const char * banner = "ML_WaveToyFO";
  CCTK_RegisterBanner(banner);
  return 0;
}
