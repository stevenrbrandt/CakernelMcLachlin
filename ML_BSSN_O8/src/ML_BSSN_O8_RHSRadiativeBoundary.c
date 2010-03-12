/*  File produced by Kranc */

#define KRANC_C

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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

void ML_BSSN_O8_RHSRadiativeBoundary_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  /* Declare predefined quantities */
  // CCTK_REAL p1o5040dx2 = INITVALUE;
  // CCTK_REAL p1o5040dy2 = INITVALUE;
  // CCTK_REAL p1o5040dz2 = INITVALUE;
  // CCTK_REAL p1o705600dxdy = INITVALUE;
  // CCTK_REAL p1o705600dxdz = INITVALUE;
  // CCTK_REAL p1o705600dydz = INITVALUE;
  // CCTK_REAL p1o840dx = INITVALUE;
  // CCTK_REAL p1o840dy = INITVALUE;
  // CCTK_REAL p1o840dz = INITVALUE;
  // CCTK_REAL p1odx = INITVALUE;
  // CCTK_REAL p1ody = INITVALUE;
  // CCTK_REAL p1odz = INITVALUE;
  // CCTK_REAL pm1o840dx = INITVALUE;
  // CCTK_REAL pm1o840dy = INITVALUE;
  // CCTK_REAL pm1o840dz = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O8_RHSRadiativeBoundary_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O8_RHSRadiativeBoundary_calc_every != ML_BSSN_O8_RHSRadiativeBoundary_calc_offset)
  {
    return;
  }
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  CCTK_REAL const dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL const dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL const dz = CCTK_DELTA_SPACE(2);
  int const di = 1;
  int const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  int const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  CCTK_REAL const dxi = 1.0 / dx;
  CCTK_REAL const dyi = 1.0 / dy;
  CCTK_REAL const dzi = 1.0 / dz;
  CCTK_REAL const khalf = 0.5;
  CCTK_REAL const kthird = 1/3.0;
  CCTK_REAL const ktwothird = 2.0/3.0;
  CCTK_REAL const kfourthird = 4.0/3.0;
  CCTK_REAL const keightthird = 8.0/3.0;
  CCTK_REAL const hdxi = 0.5 * dxi;
  CCTK_REAL const hdyi = 0.5 * dyi;
  CCTK_REAL const hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  CCTK_REAL const p1o5040dx2 = pow(dx,-2)/5040.;
  CCTK_REAL const p1o5040dy2 = pow(dy,-2)/5040.;
  CCTK_REAL const p1o5040dz2 = pow(dz,-2)/5040.;
  CCTK_REAL const p1o705600dxdy = (INV(dx)*INV(dy))/705600.;
  CCTK_REAL const p1o705600dxdz = (INV(dx)*INV(dz))/705600.;
  CCTK_REAL const p1o705600dydz = (INV(dy)*INV(dz))/705600.;
  CCTK_REAL const p1o840dx = INV(dx)/840.;
  CCTK_REAL const p1o840dy = INV(dy)/840.;
  CCTK_REAL const p1o840dz = INV(dz)/840.;
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o840dx = -INV(dx)/840.;
  CCTK_REAL const pm1o840dy = -INV(dy)/840.;
  CCTK_REAL const pm1o840dz = -INV(dz)/840.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_O8_RHSRadiativeBoundary,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    // int subblock_index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int const subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    // CCTK_REAL detgt = INITVALUE;
    // CCTK_REAL dir1 = INITVALUE, dir2 = INITVALUE, dir3 = INITVALUE;
    // CCTK_REAL em4phi = INITVALUE;
    // CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    // CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    // CCTK_REAL nlen = INITVALUE, nlen2 = INITVALUE;
    // CCTK_REAL nn1 = INITVALUE, nn2 = INITVALUE, nn3 = INITVALUE;
    // CCTK_REAL nu1 = INITVALUE, nu2 = INITVALUE, nu3 = INITVALUE;
    // CCTK_REAL su1 = INITVALUE, su2 = INITVALUE, su3 = INITVALUE;
    // CCTK_REAL vg = INITVALUE;
    
    /* Declare local copies of grid functions */
    // CCTK_REAL AL = INITVALUE;
    // CCTK_REAL alphaL = INITVALUE, alpharhsL = INITVALUE;
    // CCTK_REAL ArhsL = INITVALUE;
    // CCTK_REAL At11L = INITVALUE, At11rhsL = INITVALUE, At12L = INITVALUE, At12rhsL = INITVALUE, At13L = INITVALUE, At13rhsL = INITVALUE;
    // CCTK_REAL At22L = INITVALUE, At22rhsL = INITVALUE, At23L = INITVALUE, At23rhsL = INITVALUE, At33L = INITVALUE, At33rhsL = INITVALUE;
    // CCTK_REAL B1L = INITVALUE, B1rhsL = INITVALUE, B2L = INITVALUE, B2rhsL = INITVALUE, B3L = INITVALUE, B3rhsL = INITVALUE;
    // CCTK_REAL beta1L = INITVALUE, beta1rhsL = INITVALUE, beta2L = INITVALUE, beta2rhsL = INITVALUE, beta3L = INITVALUE, beta3rhsL = INITVALUE;
    // CCTK_REAL gt11L = INITVALUE, gt11rhsL = INITVALUE, gt12L = INITVALUE, gt12rhsL = INITVALUE, gt13L = INITVALUE, gt13rhsL = INITVALUE;
    // CCTK_REAL gt22L = INITVALUE, gt22rhsL = INITVALUE, gt23L = INITVALUE, gt23rhsL = INITVALUE, gt33L = INITVALUE, gt33rhsL = INITVALUE;
    // CCTK_REAL phiL = INITVALUE, phirhsL = INITVALUE;
    // CCTK_REAL trKL = INITVALUE, trKrhsL = INITVALUE;
    // CCTK_REAL Xt1L = INITVALUE, Xt1rhsL = INITVALUE, Xt2L = INITVALUE, Xt2rhsL = INITVALUE, Xt3L = INITVALUE, Xt3rhsL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    CCTK_REAL const AL = A[index];
    CCTK_REAL const alphaL = alpha[index];
    CCTK_REAL const At11L = At11[index];
    CCTK_REAL const At12L = At12[index];
    CCTK_REAL const At13L = At13[index];
    CCTK_REAL const At22L = At22[index];
    CCTK_REAL const At23L = At23[index];
    CCTK_REAL const At33L = At33[index];
    CCTK_REAL const B1L = B1[index];
    CCTK_REAL const B2L = B2[index];
    CCTK_REAL const B3L = B3[index];
    CCTK_REAL const beta1L = beta1[index];
    CCTK_REAL const beta2L = beta2[index];
    CCTK_REAL const beta3L = beta3[index];
    CCTK_REAL const gt11L = gt11[index];
    CCTK_REAL const gt12L = gt12[index];
    CCTK_REAL const gt13L = gt13[index];
    CCTK_REAL const gt22L = gt22[index];
    CCTK_REAL const gt23L = gt23[index];
    CCTK_REAL const gt33L = gt33[index];
    CCTK_REAL const phiL = phi[index];
    CCTK_REAL const trKL = trK[index];
    CCTK_REAL const Xt1L = Xt1[index];
    CCTK_REAL const Xt2L = Xt2[index];
    CCTK_REAL const Xt3L = Xt3[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    int const dir1  =  Sign(normal[0]);
    
    int const dir2  =  Sign(normal[1]);
    
    int const dir3  =  Sign(normal[2]);
    
    CCTK_REAL const detgt  =  1;
    
    CCTK_REAL const gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL const gtu21  =  (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL const gtu31  =  (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL const gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL const gtu32  =  (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL const gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL const em4phi  =  IfThen(conformalMethod,SQR(phiL),exp(-4*phiL));
    
    CCTK_REAL const gu11  =  em4phi*gtu11;
    
    CCTK_REAL const gu21  =  em4phi*gtu21;
    
    CCTK_REAL const gu31  =  em4phi*gtu31;
    
    CCTK_REAL const gu22  =  em4phi*gtu22;
    
    CCTK_REAL const gu32  =  em4phi*gtu32;
    
    CCTK_REAL const gu33  =  em4phi*gtu33;
    
    CCTK_REAL const nn1  =  normal[0];
    
    CCTK_REAL const nn2  =  normal[1];
    
    CCTK_REAL const nn3  =  normal[2];
    
    CCTK_REAL const nu1  =  gu11*nn1 + gu21*nn2 + gu31*nn3;
    
    CCTK_REAL const nu2  =  gu21*nn1 + gu22*nn2 + gu32*nn3;
    
    CCTK_REAL const nu3  =  gu31*nn1 + gu32*nn2 + gu33*nn3;
    
    CCTK_REAL const nlen2  =  nn1*nu1 + nn2*nu2 + nn3*nu3;
    
    CCTK_REAL const nlen  =  pow(nlen2,0.5);
    
    CCTK_REAL const su1  =  nu1*INV(nlen);
    
    CCTK_REAL const su2  =  nu2*INV(nlen);
    
    CCTK_REAL const su3  =  nu3*INV(nlen);
    
    CCTK_REAL const vg  =  pow(harmonicF,0.5);
    
    CCTK_REAL const phirhsL  =  -((PDonesided1(phi, i, j, k)*su1 + PDonesided2(phi, i, j, k)*su2 + PDonesided3(phi, i, j, k)*su3)*vg);
    
    CCTK_REAL const gt11rhsL  =  -(PDonesided1(gt11, i, j, k)*su1) - PDonesided2(gt11, i, j, k)*su2 - PDonesided3(gt11, i, j, k)*su3;
    
    CCTK_REAL const gt12rhsL  =  -(PDonesided1(gt12, i, j, k)*su1) - PDonesided2(gt12, i, j, k)*su2 - PDonesided3(gt12, i, j, k)*su3;
    
    CCTK_REAL const gt13rhsL  =  -(PDonesided1(gt13, i, j, k)*su1) - PDonesided2(gt13, i, j, k)*su2 - PDonesided3(gt13, i, j, k)*su3;
    
    CCTK_REAL const gt22rhsL  =  -(PDonesided1(gt22, i, j, k)*su1) - PDonesided2(gt22, i, j, k)*su2 - PDonesided3(gt22, i, j, k)*su3;
    
    CCTK_REAL const gt23rhsL  =  -(PDonesided1(gt23, i, j, k)*su1) - PDonesided2(gt23, i, j, k)*su2 - PDonesided3(gt23, i, j, k)*su3;
    
    CCTK_REAL const gt33rhsL  =  -(PDonesided1(gt33, i, j, k)*su1) - PDonesided2(gt33, i, j, k)*su2 - PDonesided3(gt33, i, j, k)*su3;
    
    CCTK_REAL const trKrhsL  =  -((PDonesided1(trK, i, j, k)*su1 + PDonesided2(trK, i, j, k)*su2 + PDonesided3(trK, i, j, k)*su3)*vg);
    
    CCTK_REAL const At11rhsL  =  -(PDonesided1(At11, i, j, k)*su1) - PDonesided2(At11, i, j, k)*su2 - PDonesided3(At11, i, j, k)*su3;
    
    CCTK_REAL const At12rhsL  =  -(PDonesided1(At12, i, j, k)*su1) - PDonesided2(At12, i, j, k)*su2 - PDonesided3(At12, i, j, k)*su3;
    
    CCTK_REAL const At13rhsL  =  -(PDonesided1(At13, i, j, k)*su1) - PDonesided2(At13, i, j, k)*su2 - PDonesided3(At13, i, j, k)*su3;
    
    CCTK_REAL const At22rhsL  =  -(PDonesided1(At22, i, j, k)*su1) - PDonesided2(At22, i, j, k)*su2 - PDonesided3(At22, i, j, k)*su3;
    
    CCTK_REAL const At23rhsL  =  -(PDonesided1(At23, i, j, k)*su1) - PDonesided2(At23, i, j, k)*su2 - PDonesided3(At23, i, j, k)*su3;
    
    CCTK_REAL const At33rhsL  =  -(PDonesided1(At33, i, j, k)*su1) - PDonesided2(At33, i, j, k)*su2 - PDonesided3(At33, i, j, k)*su3;
    
    CCTK_REAL const Xt1rhsL  =  -(PDonesided1(Xt1, i, j, k)*su1) - PDonesided2(Xt1, i, j, k)*su2 - PDonesided3(Xt1, i, j, k)*su3;
    
    CCTK_REAL const Xt2rhsL  =  -(PDonesided1(Xt2, i, j, k)*su1) - PDonesided2(Xt2, i, j, k)*su2 - PDonesided3(Xt2, i, j, k)*su3;
    
    CCTK_REAL const Xt3rhsL  =  -(PDonesided1(Xt3, i, j, k)*su1) - PDonesided2(Xt3, i, j, k)*su2 - PDonesided3(Xt3, i, j, k)*su3;
    
    CCTK_REAL const alpharhsL  =  -((PDonesided1(alpha, i, j, k)*su1 + PDonesided2(alpha, i, j, k)*su2 + 
            PDonesided3(alpha, i, j, k)*su3)*vg);
    
    CCTK_REAL const ArhsL  =  -((PDonesided1(A, i, j, k)*su1 + PDonesided2(A, i, j, k)*su2 + PDonesided3(A, i, j, k)*su3)*vg);
    
    CCTK_REAL const beta1rhsL  =  -(PDonesided1(beta1, i, j, k)*su1) - PDonesided2(beta1, i, j, k)*su2 - 
        PDonesided3(beta1, i, j, k)*su3;
    
    CCTK_REAL const beta2rhsL  =  -(PDonesided1(beta2, i, j, k)*su1) - PDonesided2(beta2, i, j, k)*su2 - 
        PDonesided3(beta2, i, j, k)*su3;
    
    CCTK_REAL const beta3rhsL  =  -(PDonesided1(beta3, i, j, k)*su1) - PDonesided2(beta3, i, j, k)*su2 - 
        PDonesided3(beta3, i, j, k)*su3;
    
    CCTK_REAL const B1rhsL  =  -(PDonesided1(B1, i, j, k)*su1) - PDonesided2(B1, i, j, k)*su2 - PDonesided3(B1, i, j, k)*su3;
    
    CCTK_REAL const B2rhsL  =  -(PDonesided1(B2, i, j, k)*su1) - PDonesided2(B2, i, j, k)*su2 - PDonesided3(B2, i, j, k)*su3;
    
    CCTK_REAL const B3rhsL  =  -(PDonesided1(B3, i, j, k)*su1) - PDonesided2(B3, i, j, k)*su2 - PDonesided3(B3, i, j, k)*su3;
    
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    Arhs[index] = ArhsL;
    At11rhs[index] = At11rhsL;
    At12rhs[index] = At12rhsL;
    At13rhs[index] = At13rhsL;
    At22rhs[index] = At22rhsL;
    At23rhs[index] = At23rhsL;
    At33rhs[index] = At33rhsL;
    B1rhs[index] = B1rhsL;
    B2rhs[index] = B2rhsL;
    B3rhs[index] = B3rhsL;
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    gt11rhs[index] = gt11rhsL;
    gt12rhs[index] = gt12rhsL;
    gt13rhs[index] = gt13rhsL;
    gt22rhs[index] = gt22rhsL;
    gt23rhs[index] = gt23rhsL;
    gt33rhs[index] = gt33rhsL;
    phirhs[index] = phirhsL;
    trKrhs[index] = trKrhsL;
    Xt1rhs[index] = Xt1rhsL;
    Xt2rhs[index] = Xt2rhsL;
    Xt3rhs[index] = Xt3rhsL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_BSSN_O8_RHSRadiativeBoundary);
}

void ML_BSSN_O8_RHSRadiativeBoundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverBoundary(cctkGH, &ML_BSSN_O8_RHSRadiativeBoundary_Body);
}
