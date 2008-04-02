#include <cctk.h>

#include "CactusEinstein/CoordGauge/src/Slicing.h"

int
ML_BSSN_MP_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSN_MP");
  return 0;
}
