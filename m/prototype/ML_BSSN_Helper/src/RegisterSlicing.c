#include <cctk.h>

#include "CactusEinstein/CoordGauge/src/Slicing.h"

int
ML_BSSN_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSN");
  return 0;
}
