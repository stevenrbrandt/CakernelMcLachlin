#include <cctk.h>

#include "CactusEinstein/CoordGauge/src/Slicing.h"

int
ML_BSSN6_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSN6");
  return 0;
}
