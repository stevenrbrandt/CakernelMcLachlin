#include <cctk.h>

#include "CactusEinstein/CoordGauge/src/Slicing.h"

int
ML_BSSN_O2_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSN_O2");
  return 0;
}
