#include <cctk.h>
#include <Slicing.h>

int
ML_BSSN_O8_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSN_O8");
  return 0;
}
