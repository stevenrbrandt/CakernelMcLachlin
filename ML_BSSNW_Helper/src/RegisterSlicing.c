#include <cctk.h>

#include "CactusEinstein/CoordGauge/src/Slicing.h"

int
ML_BSSNW_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSNW");
  return 0;
}
