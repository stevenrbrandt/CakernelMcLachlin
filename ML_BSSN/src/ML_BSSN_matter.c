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

void ML_BSSN_matter_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_matter_Body");
  }
  
  if (cctk_iteration % ML_BSSN_matter_calc_every != ML_BSSN_matter_calc_offset)
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
  LC_LOOP3 (ML_BSSN_matter,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL detgt = INITVALUE;
    CCTK_REAL e4phi = INITVALUE;
    CCTK_REAL em4phi = INITVALUE;
    CCTK_REAL g11 = INITVALUE, g12 = INITVALUE, g13 = INITVALUE, g22 = INITVALUE, g23 = INITVALUE, g33 = INITVALUE;
    CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    CCTK_REAL rho = INITVALUE;
    CCTK_REAL S1 = INITVALUE, S2 = INITVALUE, S3 = INITVALUE;
    CCTK_REAL T00 = INITVALUE, T01 = INITVALUE, T02 = INITVALUE, T03 = INITVALUE, T11 = INITVALUE, T12 = INITVALUE;
    CCTK_REAL T13 = INITVALUE, T22 = INITVALUE, T23 = INITVALUE, T33 = INITVALUE;
    CCTK_REAL trS = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL alphaL = INITVALUE;
    CCTK_REAL At11rhsL = INITVALUE, At12rhsL = INITVALUE, At13rhsL = INITVALUE, At22rhsL = INITVALUE, At23rhsL = INITVALUE, At33rhsL = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    CCTK_REAL eTttL = INITVALUE;
    CCTK_REAL eTtxL = INITVALUE;
    CCTK_REAL eTtyL = INITVALUE;
    CCTK_REAL eTtzL = INITVALUE;
    CCTK_REAL eTxxL = INITVALUE;
    CCTK_REAL eTxyL = INITVALUE;
    CCTK_REAL eTxzL = INITVALUE;
    CCTK_REAL eTyyL = INITVALUE;
    CCTK_REAL eTyzL = INITVALUE;
    CCTK_REAL eTzzL = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt12L = INITVALUE, gt13L = INITVALUE, gt22L = INITVALUE, gt23L = INITVALUE, gt33L = INITVALUE;
    CCTK_REAL phiL = INITVALUE;
    CCTK_REAL trKrhsL = INITVALUE;
    CCTK_REAL Xt1rhsL = INITVALUE, Xt2rhsL = INITVALUE, Xt3rhsL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    alphaL = alpha[index];
    At11rhsL = At11rhs[index];
    At12rhsL = At12rhs[index];
    At13rhsL = At13rhs[index];
    At22rhsL = At22rhs[index];
    At23rhsL = At23rhs[index];
    At33rhsL = At33rhs[index];
    beta1L = beta1[index];
    beta2L = beta2[index];
    beta3L = beta3[index];
    eTttL = eTtt[index];
    eTtxL = eTtx[index];
    eTtyL = eTty[index];
    eTtzL = eTtz[index];
    eTxxL = eTxx[index];
    eTxyL = eTxy[index];
    eTxzL = eTxz[index];
    eTyyL = eTyy[index];
    eTyzL = eTyz[index];
    eTzzL = eTzz[index];
    gt11L = gt11[index];
    gt12L = gt12[index];
    gt13L = gt13[index];
    gt22L = gt22[index];
    gt23L = gt23[index];
    gt33L = gt33[index];
    phiL = phi[index];
    trKrhsL = trKrhs[index];
    Xt1rhsL = Xt1rhs[index];
    Xt2rhsL = Xt2rhs[index];
    Xt3rhsL = Xt3rhs[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    T00  =  eTttL;
    
    T01  =  eTtxL;
    
    T02  =  eTtyL;
    
    T03  =  eTtzL;
    
    T11  =  eTxxL;
    
    T12  =  eTxyL;
    
    T13  =  eTxzL;
    
    T22  =  eTyyL;
    
    T23  =  eTyzL;
    
    T33  =  eTzzL;
    
    rho  =  pow(alphaL,-2)*(T00 - 2*(beta2L*T02 + beta3L*T03) + 
          2*(beta1L*(-T01 + beta2L*T12 + beta3L*T13) + beta2L*beta3L*T23) + T11*SQR(beta1L) + T22*SQR(beta2L) + 
          T33*SQR(beta3L));
    
    S1  =  (-T01 + beta1L*T11 + beta2L*T12 + beta3L*T13)*INV(alphaL);
    
    S2  =  (-T02 + beta1L*T12 + beta2L*T22 + beta3L*T23)*INV(alphaL);
    
    S3  =  (-T03 + beta1L*T13 + beta2L*T23 + beta3L*T33)*INV(alphaL);
    
    detgt  =  1;
    
    gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    gtu21  =  (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    gtu31  =  (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    gtu32  =  (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    e4phi  =  exp(4*phiL);
    
    em4phi  =  INV(e4phi);
    
    g11  =  e4phi*gt11L;
    
    g12  =  e4phi*gt12L;
    
    g13  =  e4phi*gt13L;
    
    g22  =  e4phi*gt22L;
    
    g23  =  e4phi*gt23L;
    
    g33  =  e4phi*gt33L;
    
    gu11  =  em4phi*gtu11;
    
    gu21  =  em4phi*gtu21;
    
    gu31  =  em4phi*gtu31;
    
    gu22  =  em4phi*gtu22;
    
    gu32  =  em4phi*gtu32;
    
    gu33  =  em4phi*gtu33;
    
    trS  =  gu11*T11 + gu22*T22 + 2*(gu21*T12 + gu31*T13 + gu32*T23) + gu33*T33;
    
    trKrhsL  =  trKrhsL + 12.56637061435917295385057353311801153679*alphaL*(rho + trS);
    
    At11rhsL  =  At11rhsL + alphaL*em4phi*(-25.13274122871834590770114706623602307358*T11 + 
           8.377580409572781969233715688745341024526*g11*trS);
    
    At12rhsL  =  At12rhsL + alphaL*em4phi*(-25.13274122871834590770114706623602307358*T12 + 
           8.377580409572781969233715688745341024526*g12*trS);
    
    At13rhsL  =  At13rhsL + alphaL*em4phi*(-25.13274122871834590770114706623602307358*T13 + 
           8.377580409572781969233715688745341024526*g13*trS);
    
    At22rhsL  =  At22rhsL + alphaL*em4phi*(-25.13274122871834590770114706623602307358*T22 + 
           8.377580409572781969233715688745341024526*g22*trS);
    
    At23rhsL  =  At23rhsL + alphaL*em4phi*(-25.13274122871834590770114706623602307358*T23 + 
           8.377580409572781969233715688745341024526*g23*trS);
    
    At33rhsL  =  At33rhsL + alphaL*em4phi*(-25.13274122871834590770114706623602307358*T33 + 
           8.377580409572781969233715688745341024526*g33*trS);
    
    Xt1rhsL  =  -50.26548245743669181540229413247204614715*alphaL*(gtu11*S1 + gtu21*S2 + gtu31*S3) + Xt1rhsL;
    
    Xt2rhsL  =  -50.26548245743669181540229413247204614715*alphaL*(gtu21*S1 + gtu22*S2 + gtu32*S3) + Xt2rhsL;
    
    Xt3rhsL  =  -50.26548245743669181540229413247204614715*alphaL*(gtu31*S1 + gtu32*S2 + gtu33*S3) + Xt3rhsL;
    
    
    /* Copy local copies back to grid functions */
    At11rhs[index] = At11rhsL;
    At12rhs[index] = At12rhsL;
    At13rhs[index] = At13rhsL;
    At22rhs[index] = At22rhsL;
    At23rhs[index] = At23rhsL;
    At33rhs[index] = At33rhsL;
    trKrhs[index] = trKrhsL;
    Xt1rhs[index] = Xt1rhsL;
    Xt2rhs[index] = Xt2rhsL;
    Xt3rhs[index] = Xt3rhsL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_BSSN_matter);
}

void ML_BSSN_matter(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_matter_Body);
}
