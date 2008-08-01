#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>

#include <Symmetry.h>

static void
setcartsym (cGH const * restrict cctkGH, char const * restrict vn, int symdesc);

void
ML_BSSN_RegisterSymmetry (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  setcartsym (cctkGH, "ML_BSSN::phi"  ,  0);
  setcartsym (cctkGH, "ML_BSSN::gt11" , 11);
  setcartsym (cctkGH, "ML_BSSN::gt12" , 12);
  setcartsym (cctkGH, "ML_BSSN::gt13" , 13);
  setcartsym (cctkGH, "ML_BSSN::gt22" , 22);
  setcartsym (cctkGH, "ML_BSSN::gt23" , 23);
  setcartsym (cctkGH, "ML_BSSN::gt33" , 33);
  setcartsym (cctkGH, "ML_BSSN::trK"  ,  0);
  setcartsym (cctkGH, "ML_BSSN::At11" , 11);
  setcartsym (cctkGH, "ML_BSSN::At12" , 12);
  setcartsym (cctkGH, "ML_BSSN::At13" , 13);
  setcartsym (cctkGH, "ML_BSSN::At22" , 22);
  setcartsym (cctkGH, "ML_BSSN::At23" , 23);
  setcartsym (cctkGH, "ML_BSSN::At33" , 33);
  setcartsym (cctkGH, "ML_BSSN::Xt1"  ,  1);
  setcartsym (cctkGH, "ML_BSSN::Xt2"  ,  2);
  setcartsym (cctkGH, "ML_BSSN::Xt3"  ,  3);
  setcartsym (cctkGH, "ML_BSSN::alpha",  0);
  setcartsym (cctkGH, "ML_BSSN::beta1",  1);
  setcartsym (cctkGH, "ML_BSSN::beta2",  2);
  setcartsym (cctkGH, "ML_BSSN::beta3",  3);
  setcartsym (cctkGH, "ML_BSSN::A"    ,  0);
  setcartsym (cctkGH, "ML_BSSN::B1"   ,  1);
  setcartsym (cctkGH, "ML_BSSN::B2"   ,  2);
  setcartsym (cctkGH, "ML_BSSN::B3"   ,  3);
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
