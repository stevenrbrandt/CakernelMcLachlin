#include <cctk.h>
#include <Slicing.h>

int
ML_BSSN_UPW_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("ML_BSSN_UPW");
  return 0;
}
