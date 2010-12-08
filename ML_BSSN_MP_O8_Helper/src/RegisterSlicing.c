#include <cctk.h>
#include <Slicing.h>

int
ML_BSSN_MP_O8_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSN_MP_O8");
  return 0;
}
