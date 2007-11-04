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
  LC_LOOP3 (somename,
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
    CCTK_REAL g11 = INITVALUE, g21 = INITVALUE, g22 = INITVALUE, g31 = INITVALUE, g32 = INITVALUE, g33 = INITVALUE;
    CCTK_REAL gu11 = INITVALUE, gu12 = INITVALUE, gu13 = INITVALUE, gu22 = INITVALUE, gu23 = INITVALUE, gu33 = INITVALUE;
    CCTK_REAL K11 = INITVALUE, K21 = INITVALUE, K22 = INITVALUE, K31 = INITVALUE, K32 = INITVALUE, K33 = INITVALUE;
    CCTK_REAL Km11 = INITVALUE, Km22 = INITVALUE, Km33 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL alpL = INITVALUE;
    CCTK_REAL alphaL = INITVALUE;
    CCTK_REAL At11L = INITVALUE, At21L = INITVALUE, At22L = INITVALUE, At31L = INITVALUE, At32L = INITVALUE, At33L = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    CCTK_REAL betaxL = INITVALUE;
    CCTK_REAL betayL = INITVALUE;
    CCTK_REAL betazL = INITVALUE;
    CCTK_REAL dtalpL = INITVALUE;
    CCTK_REAL dtalphaL = INITVALUE;
    CCTK_REAL dtbeta1L = INITVALUE, dtbeta2L = INITVALUE, dtbeta3L = INITVALUE;
    CCTK_REAL dtbetaxL = INITVALUE;
    CCTK_REAL dtbetayL = INITVALUE;
    CCTK_REAL dtbetazL = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt21L = INITVALUE, gt22L = INITVALUE, gt31L = INITVALUE, gt32L = INITVALUE, gt33L = INITVALUE;
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
    dtalpL = dtalp[index];
    dtbetaxL = dtbetax[index];
    dtbetayL = dtbetay[index];
    dtbetazL = dtbetaz[index];
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
    trKL = trK[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    g11  =  gxxL;
    
    g21  =  gxyL;
    
    g31  =  gxzL;
    
    g22  =  gyyL;
    
    g32  =  gyzL;
    
    g33  =  gzzL;
    
    detg  =  2*g21*g31*g32 + g33*(g11*g22 - SQR(g21)) - g22*SQR(g31) - g11*SQR(g32);
    
    gu11  =  INV(detg)*(g22*g33 - SQR(g32));
    
    gu12  =  (g31*g32 - g21*g33)*INV(detg);
    
    gu13  =  (-(g22*g31) + g21*g32)*INV(detg);
    
    gu22  =  INV(detg)*(g11*g33 - SQR(g31));
    
    gu23  =  (g21*g31 - g11*g32)*INV(detg);
    
    gu33  =  INV(detg)*(g11*g22 - SQR(g21));
    
    em4phi  =  pow(detg,-3);
    
    phiL  =  Log(detg)/12.;
    
    gt11L  =  em4phi*g11;
    
    gt21L  =  em4phi*g21;
    
    gt31L  =  em4phi*g31;
    
    gt22L  =  em4phi*g22;
    
    gt32L  =  em4phi*g32;
    
    gt33L  =  em4phi*g33;
    
    K11  =  kxxL;
    
    K21  =  kxyL;
    
    K31  =  kxzL;
    
    K22  =  kyyL;
    
    K32  =  kyzL;
    
    K33  =  kzzL;
    
    Km11  =  gu11*K11 + gu12*K21 + gu13*K31;
    
    Km22  =  gu12*K21 + gu22*K22 + gu23*K32;
    
    Km33  =  gu13*K31 + gu23*K32 + gu33*K33;
    
    trKL  =  Km11 + Km22 + Km33;
    
    At11L  =  em4phi*(K11 - g11*kthird*trKL);
    
    At21L  =  em4phi*(K21 - g21*kthird*trKL);
    
    At31L  =  em4phi*(K31 - g31*kthird*trKL);
    
    At22L  =  em4phi*(K22 - g22*kthird*trKL);
    
    At32L  =  em4phi*(K32 - g32*kthird*trKL);
    
    At33L  =  em4phi*(K33 - g33*kthird*trKL);
    
    alphaL  =  alpL;
    
    dtalphaL  =  dtalpL;
    
    beta1L  =  betaxL;
    
    beta2L  =  betayL;
    
    beta3L  =  betazL;
    
    dtbeta1L  =  dtbetaxL;
    
    dtbeta2L  =  dtbetayL;
    
    dtbeta3L  =  dtbetazL;
    
    
    /* Copy local copies back to grid functions */
    alpha[index] = alphaL;
    At11[index] = At11L;
    At21[index] = At21L;
    At22[index] = At22L;
    At31[index] = At31L;
    At32[index] = At32L;
    At33[index] = At33L;
    beta1[index] = beta1L;
    beta2[index] = beta2L;
    beta3[index] = beta3L;
    dtalpha[index] = dtalphaL;
    dtbeta1[index] = dtbeta1L;
    dtbeta2[index] = dtbeta2L;
    dtbeta3[index] = dtbeta3L;
    gt11[index] = gt11L;
    gt21[index] = gt21L;
    gt22[index] = gt22L;
    gt31[index] = gt31L;
    gt32[index] = gt32L;
    gt33[index] = gt33L;
    phi[index] = phiL;
    trK[index] = trKL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (somename);
}

void ML_BSSN_convertFromADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverEverything(cctkGH, &ML_BSSN_convertFromADMBase_Body);
}
