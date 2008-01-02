/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (32-bit) (April 20, 2007) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#define KRANC_C

#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"
#include "Differencing.h"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

void ML_BSSN_convertFromADMBase_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  
  /* Declare finite differencing variables */
  CCTK_REAL dx = INITVALUE, dy = INITVALUE, dz = INITVALUE;
  CCTK_REAL dxi = INITVALUE, dyi = INITVALUE, dzi = INITVALUE;
  CCTK_REAL khalf = INITVALUE, kthird = INITVALUE, ktwothird = INITVALUE, kfourthird = INITVALUE, keightthird = INITVALUE;
  CCTK_REAL hdxi = INITVALUE, hdyi = INITVALUE, hdzi = INITVALUE;
  
  
  /* Declare predefined quantities */
  CCTK_REAL p1o12dx = INITVALUE;
  CCTK_REAL p1o12dy = INITVALUE;
  CCTK_REAL p1o12dz = INITVALUE;
  CCTK_REAL p1o144dxdy = INITVALUE;
  CCTK_REAL p1o144dxdz = INITVALUE;
  CCTK_REAL p1o144dydz = INITVALUE;
  CCTK_REAL pm1o12dx2 = INITVALUE;
  CCTK_REAL pm1o12dy2 = INITVALUE;
  CCTK_REAL pm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_convertFromADMBase_Body");
  }
  
  if (cctk_iteration % ML_BSSN_convertFromADMBase_calc_every != ML_BSSN_convertFromADMBase_calc_offset)
  {
    return;
  }
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  dx = CCTK_DELTA_SPACE(0);
  dy = CCTK_DELTA_SPACE(1);
  dz = CCTK_DELTA_SPACE(2);
  dxi = 1.0 / dx;
  dyi = 1.0 / dy;
  dzi = 1.0 / dz;
  khalf = 0.5;
  kthird = 1/3.0;
  ktwothird = 2.0/3.0;
  kfourthird = 4.0/3.0;
  keightthird = 8.0/3.0;
  hdxi = 0.5 * dxi;
  hdyi = 0.5 * dyi;
  hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  p1o12dx = INV(dx)/12.;
  p1o12dy = INV(dy)/12.;
  p1o12dz = INV(dz)/12.;
  p1o144dxdy = (INV(dx)*INV(dy))/144.;
  p1o144dxdz = (INV(dx)*INV(dz))/144.;
  p1o144dydz = (INV(dy)*INV(dz))/144.;
  pm1o12dx2 = -pow(dx,-2)/12.;
  pm1o12dy2 = -pow(dy,-2)/12.;
  pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  _Pragma ("omp parallel")
  LC_LOOP3 (ML_BSSN_convertFromADMBase,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL detg = INITVALUE;
    CCTK_REAL em4phi = INITVALUE;
    CCTK_REAL g11 = INITVALUE, g12 = INITVALUE, g13 = INITVALUE, g22 = INITVALUE, g23 = INITVALUE, g33 = INITVALUE;
    CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    CCTK_REAL K11 = INITVALUE, K12 = INITVALUE, K13 = INITVALUE, K22 = INITVALUE, K23 = INITVALUE, K33 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL alpL = INITVALUE;
    CCTK_REAL alphaL = INITVALUE;
    CCTK_REAL At11L = INITVALUE, At12L = INITVALUE, At13L = INITVALUE, At22L = INITVALUE, At23L = INITVALUE, At33L = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    CCTK_REAL betaxL = INITVALUE;
    CCTK_REAL betayL = INITVALUE;
    CCTK_REAL betazL = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt12L = INITVALUE, gt13L = INITVALUE, gt22L = INITVALUE, gt23L = INITVALUE, gt33L = INITVALUE;
    CCTK_REAL gxxL = INITVALUE;
    CCTK_REAL gxyL = INITVALUE;
    CCTK_REAL gxzL = INITVALUE;
    CCTK_REAL gyyL = INITVALUE;
    CCTK_REAL gyzL = INITVALUE;
    CCTK_REAL gzzL = INITVALUE;
    CCTK_REAL kxxL = INITVALUE;
    CCTK_REAL kxyL = INITVALUE;
    CCTK_REAL kxzL = INITVALUE;
    CCTK_REAL kyyL = INITVALUE;
    CCTK_REAL kyzL = INITVALUE;
    CCTK_REAL kzzL = INITVALUE;
    CCTK_REAL phiL = INITVALUE;
    CCTK_REAL trKL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    alpL = alp[index];
    betaxL = betax[index];
    betayL = betay[index];
    betazL = betaz[index];
    gxxL = gxx[index];
    gxyL = gxy[index];
    gxzL = gxz[index];
    gyyL = gyy[index];
    gyzL = gyz[index];
    gzzL = gzz[index];
    kxxL = kxx[index];
    kxyL = kxy[index];
    kxzL = kxz[index];
    kyyL = kyy[index];
    kyzL = kyz[index];
    kzzL = kzz[index];
    phiL = phi[index];
    trKL = trK[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    g11  =  gxxL;
    
    g12  =  gxyL;
    
    g13  =  gxzL;
    
    g22  =  gyyL;
    
    g23  =  gyzL;
    
    g33  =  gzzL;
    
    detg  =  2*g12*g13*g23 + g33*(g11*g22 - SQR(g12)) - g22*SQR(g13) - g11*SQR(g23);
    
    gu11  =  INV(detg)*(g22*g33 - SQR(g23));
    
    gu21  =  (g13*g23 - g12*g33)*INV(detg);
    
    gu31  =  (-(g13*g22) + g12*g23)*INV(detg);
    
    gu22  =  INV(detg)*(g11*g33 - SQR(g13));
    
    gu32  =  (g12*g13 - g11*g23)*INV(detg);
    
    gu33  =  INV(detg)*(g11*g22 - SQR(g12));
    
    phiL  =  Log(detg)/12.;
    
    em4phi  =  exp(-4*phiL);
    
    gt11L  =  em4phi*g11;
    
    gt12L  =  em4phi*g12;
    
    gt13L  =  em4phi*g13;
    
    gt22L  =  em4phi*g22;
    
    gt23L  =  em4phi*g23;
    
    gt33L  =  em4phi*g33;
    
    K11  =  kxxL;
    
    K12  =  kxyL;
    
    K13  =  kxzL;
    
    K22  =  kyyL;
    
    K23  =  kyzL;
    
    K33  =  kzzL;
    
    trKL  =  gu11*K11 + gu22*K22 + 2*(gu21*K12 + gu31*K13 + gu32*K23) + gu33*K33;
    
    At11L  =  em4phi*(K11 - g11*kthird*trKL);
    
    At12L  =  em4phi*(K12 - g12*kthird*trKL);
    
    At13L  =  em4phi*(K13 - g13*kthird*trKL);
    
    At22L  =  em4phi*(K22 - g22*kthird*trKL);
    
    At23L  =  em4phi*(K23 - g23*kthird*trKL);
    
    At33L  =  em4phi*(K33 - g33*kthird*trKL);
    
    alphaL  =  alpL;
    
    beta1L  =  betaxL;
    
    beta2L  =  betayL;
    
    beta3L  =  betazL;
    
    
    /* Copy local copies back to grid functions */
    alpha[index] = alphaL;
    At11[index] = At11L;
    At12[index] = At12L;
    At13[index] = At13L;
    At22[index] = At22L;
    At23[index] = At23L;
    At33[index] = At33L;
    beta1[index] = beta1L;
    beta2[index] = beta2L;
    beta3[index] = beta3L;
    gt11[index] = gt11L;
    gt12[index] = gt12L;
    gt13[index] = gt13L;
    gt22[index] = gt22L;
    gt23[index] = gt23L;
    gt33[index] = gt33L;
    phi[index] = phiL;
    trK[index] = trKL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_BSSN_convertFromADMBase);
}

void ML_BSSN_convertFromADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverEverything(cctkGH, &ML_BSSN_convertFromADMBase_Body);
}
