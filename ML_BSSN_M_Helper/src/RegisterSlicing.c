#include <cctk.h>

#include "CactusEinstein/CoordGauge/src/Slicing.h"

int
ML_BSSN_M_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSN_M");
  return 0;
}
