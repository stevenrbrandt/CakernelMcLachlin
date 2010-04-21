#include <cctk.h>
#include <Slicing.h>

int
ML_BSSN_O2_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSN_O2");
  return 0;
}
