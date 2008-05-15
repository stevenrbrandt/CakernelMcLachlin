#include <cctk.h>

#include "CactusEinstein/CoordGauge/src/Slicing.h"

int
ML_BSSNUp_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSNUp");
  return 0;
}
