#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>

#include <Symmetry.h>

static void
setcartsym (cGH const * restrict cctkGH, char const * restrict vn, int symdesc);

void
ML_BSSN_M_RegisterSymmetry (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  setcartsym (cctkGH, "ML_BSSN_M::phi"  ,  0);
  setcartsym (cctkGH, "ML_BSSN_M::gt11" , 11);
  setcartsym (cctkGH, "ML_BSSN_M::gt12" , 12);
  setcartsym (cctkGH, "ML_BSSN_M::gt13" , 13);
  setcartsym (cctkGH, "ML_BSSN_M::gt22" , 22);
  setcartsym (cctkGH, "ML_BSSN_M::gt23" , 23);
  setcartsym (cctkGH, "ML_BSSN_M::gt33" , 33);
  setcartsym (cctkGH, "ML_BSSN_M::trK"  ,  0);
  setcartsym (cctkGH, "ML_BSSN_M::At11" , 11);
  setcartsym (cctkGH, "ML_BSSN_M::At12" , 12);
  setcartsym (cctkGH, "ML_BSSN_M::At13" , 13);
  setcartsym (cctkGH, "ML_BSSN_M::At22" , 22);
  setcartsym (cctkGH, "ML_BSSN_M::At23" , 23);
  setcartsym (cctkGH, "ML_BSSN_M::At33" , 33);
  setcartsym (cctkGH, "ML_BSSN_M::Xt1"  ,  1);
  setcartsym (cctkGH, "ML_BSSN_M::Xt2"  ,  2);
  setcartsym (cctkGH, "ML_BSSN_M::Xt3"  ,  3);
  setcartsym (cctkGH, "ML_BSSN_M::alpha",  0);
  setcartsym (cctkGH, "ML_BSSN_M::beta1",  1);
  setcartsym (cctkGH, "ML_BSSN_M::beta2",  2);
  setcartsym (cctkGH, "ML_BSSN_M::beta3",  3);
  setcartsym (cctkGH, "ML_BSSN_M::A"    ,  0);
  setcartsym (cctkGH, "ML_BSSN_M::B1"   ,  1);
  setcartsym (cctkGH, "ML_BSSN_M::B2"   ,  2);
  setcartsym (cctkGH, "ML_BSSN_M::B3"   ,  3);
}

void
setcartsym (cGH const * restrict const cctkGH,
            char const * restrict const vn,
            int symdesc)
{
  int sym[3];
  
  for (int d=0; d<3; ++d) {
    sym[d] = +1;
  }
  
  assert (symdesc >= 0);
  while (symdesc > 0) {
    int const d = symdesc % 10;
    assert (d>=1 && d<=3);
    sym[d-1] *= -1;
    symdesc /= 10;
  }
  
  int const ierr = SetCartSymVN (cctkGH, sym, vn);
  if (ierr != 0) {
    CCTK_WARN (CCTK_WARN_ABORT, "internal error");
  }
}
