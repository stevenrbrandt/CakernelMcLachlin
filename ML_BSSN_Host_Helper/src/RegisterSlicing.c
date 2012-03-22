#include <cctk.h>
#include <Slicing.h>

int
ML_BSSN_Host_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSN_Host");
  return 0;
}
