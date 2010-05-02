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
#include "Vectors.hh"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

static void ML_BSSN_O8_RHSRadiativeBoundary_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
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
  CCTK_REAL const p1o1024dx = INV(dx)/1024.;
  CCTK_REAL const p1o1024dy = INV(dy)/1024.;
  CCTK_REAL const p1o1024dz = INV(dz)/1024.;
  CCTK_REAL const p1o1680dx = INV(dx)/1680.;
  CCTK_REAL const p1o1680dy = INV(dy)/1680.;
  CCTK_REAL const p1o1680dz = INV(dz)/1680.;
  CCTK_REAL const p1o5040dx2 = pow(dx,-2)/5040.;
  CCTK_REAL const p1o5040dy2 = pow(dy,-2)/5040.;
  CCTK_REAL const p1o5040dz2 = pow(dz,-2)/5040.;
  CCTK_REAL const p1o560dx = INV(dx)/560.;
  CCTK_REAL const p1o560dy = INV(dy)/560.;
  CCTK_REAL const p1o560dz = INV(dz)/560.;
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
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    CCTK_REAL_VEC  AL = vec_load(A[index]);
    CCTK_REAL_VEC  alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC  At11L = vec_load(At11[index]);
    CCTK_REAL_VEC  At12L = vec_load(At12[index]);
    CCTK_REAL_VEC  At13L = vec_load(At13[index]);
    CCTK_REAL_VEC  At22L = vec_load(At22[index]);
    CCTK_REAL_VEC  At23L = vec_load(At23[index]);
    CCTK_REAL_VEC  At33L = vec_load(At33[index]);
    CCTK_REAL_VEC  B1L = vec_load(B1[index]);
    CCTK_REAL_VEC  B2L = vec_load(B2[index]);
    CCTK_REAL_VEC  B3L = vec_load(B3[index]);
    CCTK_REAL_VEC  beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC  beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC  beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC  gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC  gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC  gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC  gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC  gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC  gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC  phiL = vec_load(phi[index]);
    CCTK_REAL_VEC  trKL = vec_load(trK[index]);
    CCTK_REAL_VEC  Xt1L = vec_load(Xt1[index]);
    CCTK_REAL_VEC  Xt2L = vec_load(Xt2[index]);
    CCTK_REAL_VEC  Xt3L = vec_load(Xt3[index]);
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    int dir1 = Sign(normal[0]);
    
    int dir2 = Sign(normal[1]);
    
    int dir3 = Sign(normal[2]);
    
    CCTK_REAL_VEC detgt = 1;
    
    CCTK_REAL_VEC gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL_VEC gtu21 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL_VEC gtu31 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL_VEC gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL_VEC gtu32 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL_VEC gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL_VEC em4phi = IfThen(conformalMethod,SQR(phiL),exp(-4*phiL));
    
    CCTK_REAL_VEC gu11 = em4phi*gtu11;
    
    CCTK_REAL_VEC gu21 = em4phi*gtu21;
    
    CCTK_REAL_VEC gu31 = em4phi*gtu31;
    
    CCTK_REAL_VEC gu22 = em4phi*gtu22;
    
    CCTK_REAL_VEC gu32 = em4phi*gtu32;
    
    CCTK_REAL_VEC gu33 = em4phi*gtu33;
    
    CCTK_REAL_VEC nn1 = normal[0];
    
    CCTK_REAL_VEC nn2 = normal[1];
    
    CCTK_REAL_VEC nn3 = normal[2];
    
    CCTK_REAL_VEC nu1 = gu11*nn1 + gu21*nn2 + gu31*nn3;
    
    CCTK_REAL_VEC nu2 = gu21*nn1 + gu22*nn2 + gu32*nn3;
    
    CCTK_REAL_VEC nu3 = gu31*nn1 + gu32*nn2 + gu33*nn3;
    
    CCTK_REAL_VEC nlen2 = nn1*nu1 + nn2*nu2 + nn3*nu3;
    
    CCTK_REAL_VEC nlen = pow(nlen2,0.5);
    
    CCTK_REAL_VEC su1 = nu1*INV(nlen);
    
    CCTK_REAL_VEC su2 = nu2*INV(nlen);
    
    CCTK_REAL_VEC su3 = nu3*INV(nlen);
    
    CCTK_REAL_VEC vg = pow(harmonicF,0.5);
    
    CCTK_REAL_VEC phirhsL = -((PDonesided1(phi, i, j, k)*su1 + 
      PDonesided2(phi, i, j, k)*su2 + PDonesided3(phi, i, j, k)*su3)*vg);
    
    CCTK_REAL_VEC gt11rhsL = -(PDonesided1(gt11, i, j, k)*su1) - 
      PDonesided2(gt11, i, j, k)*su2 - PDonesided3(gt11, i, j, k)*su3;
    
    CCTK_REAL_VEC gt12rhsL = -(PDonesided1(gt12, i, j, k)*su1) - 
      PDonesided2(gt12, i, j, k)*su2 - PDonesided3(gt12, i, j, k)*su3;
    
    CCTK_REAL_VEC gt13rhsL = -(PDonesided1(gt13, i, j, k)*su1) - 
      PDonesided2(gt13, i, j, k)*su2 - PDonesided3(gt13, i, j, k)*su3;
    
    CCTK_REAL_VEC gt22rhsL = -(PDonesided1(gt22, i, j, k)*su1) - 
      PDonesided2(gt22, i, j, k)*su2 - PDonesided3(gt22, i, j, k)*su3;
    
    CCTK_REAL_VEC gt23rhsL = -(PDonesided1(gt23, i, j, k)*su1) - 
      PDonesided2(gt23, i, j, k)*su2 - PDonesided3(gt23, i, j, k)*su3;
    
    CCTK_REAL_VEC gt33rhsL = -(PDonesided1(gt33, i, j, k)*su1) - 
      PDonesided2(gt33, i, j, k)*su2 - PDonesided3(gt33, i, j, k)*su3;
    
    CCTK_REAL_VEC trKrhsL = -((PDonesided1(trK, i, j, k)*su1 + 
      PDonesided2(trK, i, j, k)*su2 + PDonesided3(trK, i, j, k)*su3)*vg);
    
    CCTK_REAL_VEC At11rhsL = -(PDonesided1(At11, i, j, k)*su1) - 
      PDonesided2(At11, i, j, k)*su2 - PDonesided3(At11, i, j, k)*su3;
    
    CCTK_REAL_VEC At12rhsL = -(PDonesided1(At12, i, j, k)*su1) - 
      PDonesided2(At12, i, j, k)*su2 - PDonesided3(At12, i, j, k)*su3;
    
    CCTK_REAL_VEC At13rhsL = -(PDonesided1(At13, i, j, k)*su1) - 
      PDonesided2(At13, i, j, k)*su2 - PDonesided3(At13, i, j, k)*su3;
    
    CCTK_REAL_VEC At22rhsL = -(PDonesided1(At22, i, j, k)*su1) - 
      PDonesided2(At22, i, j, k)*su2 - PDonesided3(At22, i, j, k)*su3;
    
    CCTK_REAL_VEC At23rhsL = -(PDonesided1(At23, i, j, k)*su1) - 
      PDonesided2(At23, i, j, k)*su2 - PDonesided3(At23, i, j, k)*su3;
    
    CCTK_REAL_VEC At33rhsL = -(PDonesided1(At33, i, j, k)*su1) - 
      PDonesided2(At33, i, j, k)*su2 - PDonesided3(At33, i, j, k)*su3;
    
    CCTK_REAL_VEC Xt1rhsL = -(PDonesided1(Xt1, i, j, k)*su1) - 
      PDonesided2(Xt1, i, j, k)*su2 - PDonesided3(Xt1, i, j, k)*su3;
    
    CCTK_REAL_VEC Xt2rhsL = -(PDonesided1(Xt2, i, j, k)*su1) - 
      PDonesided2(Xt2, i, j, k)*su2 - PDonesided3(Xt2, i, j, k)*su3;
    
    CCTK_REAL_VEC Xt3rhsL = -(PDonesided1(Xt3, i, j, k)*su1) - 
      PDonesided2(Xt3, i, j, k)*su2 - PDonesided3(Xt3, i, j, k)*su3;
    
    CCTK_REAL_VEC alpharhsL = -((PDonesided1(alpha, i, j, k)*su1 + 
      PDonesided2(alpha, i, j, k)*su2 + PDonesided3(alpha, i, j, 
      k)*su3)*vg);
    
    CCTK_REAL_VEC ArhsL = -((PDonesided1(A, i, j, k)*su1 + 
      PDonesided2(A, i, j, k)*su2 + PDonesided3(A, i, j, k)*su3)*vg);
    
    CCTK_REAL_VEC beta1rhsL = -(PDonesided1(beta1, i, j, k)*su1) - 
      PDonesided2(beta1, i, j, k)*su2 - PDonesided3(beta1, i, j, k)*su3;
    
    CCTK_REAL_VEC beta2rhsL = -(PDonesided1(beta2, i, j, k)*su1) - 
      PDonesided2(beta2, i, j, k)*su2 - PDonesided3(beta2, i, j, k)*su3;
    
    CCTK_REAL_VEC beta3rhsL = -(PDonesided1(beta3, i, j, k)*su1) - 
      PDonesided2(beta3, i, j, k)*su2 - PDonesided3(beta3, i, j, k)*su3;
    
    CCTK_REAL_VEC B1rhsL = -(PDonesided1(B1, i, j, k)*su1) - 
      PDonesided2(B1, i, j, k)*su2 - PDonesided3(B1, i, j, k)*su3;
    
    CCTK_REAL_VEC B2rhsL = -(PDonesided1(B2, i, j, k)*su1) - 
      PDonesided2(B2, i, j, k)*su2 - PDonesided3(B2, i, j, k)*su3;
    
    CCTK_REAL_VEC B3rhsL = -(PDonesided1(B3, i, j, k)*su1) - 
      PDonesided2(B3, i, j, k)*su2 - PDonesided3(B3, i, j, k)*su3;
    
    
    /* Copy local copies back to grid functions */
    vec_store_nta(alpharhs[index],alpharhsL);
    vec_store_nta(Arhs[index],ArhsL);
    vec_store_nta(At11rhs[index],At11rhsL);
    vec_store_nta(At12rhs[index],At12rhsL);
    vec_store_nta(At13rhs[index],At13rhsL);
    vec_store_nta(At22rhs[index],At22rhsL);
    vec_store_nta(At23rhs[index],At23rhsL);
    vec_store_nta(At33rhs[index],At33rhsL);
    vec_store_nta(B1rhs[index],B1rhsL);
    vec_store_nta(B2rhs[index],B2rhsL);
    vec_store_nta(B3rhs[index],B3rhsL);
    vec_store_nta(beta1rhs[index],beta1rhsL);
    vec_store_nta(beta2rhs[index],beta2rhsL);
    vec_store_nta(beta3rhs[index],beta3rhsL);
    vec_store_nta(gt11rhs[index],gt11rhsL);
    vec_store_nta(gt12rhs[index],gt12rhsL);
    vec_store_nta(gt13rhs[index],gt13rhsL);
    vec_store_nta(gt22rhs[index],gt22rhsL);
    vec_store_nta(gt23rhs[index],gt23rhsL);
    vec_store_nta(gt33rhs[index],gt33rhsL);
    vec_store_nta(phirhs[index],phirhsL);
    vec_store_nta(trKrhs[index],trKrhsL);
    vec_store_nta(Xt1rhs[index],Xt1rhsL);
    vec_store_nta(Xt2rhs[index],Xt2rhsL);
    vec_store_nta(Xt3rhs[index],Xt3rhsL);
  i += CCTK_REAL_VEC_SIZE-1;
  }
  LC_ENDLOOP3 (ML_BSSN_O8_RHSRadiativeBoundary);
}

extern "C" void ML_BSSN_O8_RHSRadiativeBoundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverBoundary(cctkGH, &ML_BSSN_O8_RHSRadiativeBoundary_Body);
}
