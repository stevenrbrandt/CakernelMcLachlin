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

void ML_BSSN_O2_RHSRadiativeBoundary_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O2::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O2::ML_curvrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O2::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O2::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O2::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O2::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O2::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O2::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O2::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O2::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O2::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O2::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O2::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O2::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O2::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O2::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O2::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O2::ML_trace_curvrhs.");
  return;
}

void ML_BSSN_O2_RHSRadiativeBoundary_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O2_RHSRadiativeBoundary_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O2_RHSRadiativeBoundary_calc_every != ML_BSSN_O2_RHSRadiativeBoundary_calc_offset)
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
  CCTK_REAL const p1o16dx = INV(dx)/16.;
  CCTK_REAL const p1o16dy = INV(dy)/16.;
  CCTK_REAL const p1o16dz = INV(dz)/16.;
  CCTK_REAL const p1o2dx = khalf*INV(dx);
  CCTK_REAL const p1o2dy = khalf*INV(dy);
  CCTK_REAL const p1o2dz = khalf*INV(dz);
  CCTK_REAL const p1o4dx = INV(dx)/4.;
  CCTK_REAL const p1o4dxdy = (INV(dx)*INV(dy))/4.;
  CCTK_REAL const p1o4dxdz = (INV(dx)*INV(dz))/4.;
  CCTK_REAL const p1o4dy = INV(dy)/4.;
  CCTK_REAL const p1o4dydz = (INV(dy)*INV(dz))/4.;
  CCTK_REAL const p1o4dz = INV(dz)/4.;
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1odx2 = pow(dx,-2);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1ody2 = pow(dy,-2);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const p1odz2 = pow(dz,-2);
  CCTK_REAL const pm1o2dx = -(khalf*INV(dx));
  CCTK_REAL const pm1o2dy = -(khalf*INV(dy));
  CCTK_REAL const pm1o2dz = -(khalf*INV(dz));
  CCTK_REAL const pm1o4dx = -INV(dx)/4.;
  CCTK_REAL const pm1o4dy = -INV(dy)/4.;
  CCTK_REAL const pm1o4dz = -INV(dz)/4.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_O2_RHSRadiativeBoundary,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    CCTK_REAL  AL = A[index];
    CCTK_REAL  alphaL = alpha[index];
    CCTK_REAL  At11L = At11[index];
    CCTK_REAL  At12L = At12[index];
    CCTK_REAL  At13L = At13[index];
    CCTK_REAL  At22L = At22[index];
    CCTK_REAL  At23L = At23[index];
    CCTK_REAL  At33L = At33[index];
    CCTK_REAL  B1L = B1[index];
    CCTK_REAL  B2L = B2[index];
    CCTK_REAL  B3L = B3[index];
    CCTK_REAL  beta1L = beta1[index];
    CCTK_REAL  beta2L = beta2[index];
    CCTK_REAL  beta3L = beta3[index];
    CCTK_REAL  gt11L = gt11[index];
    CCTK_REAL  gt12L = gt12[index];
    CCTK_REAL  gt13L = gt13[index];
    CCTK_REAL  gt22L = gt22[index];
    CCTK_REAL  gt23L = gt23[index];
    CCTK_REAL  gt33L = gt33[index];
    CCTK_REAL  phiL = phi[index];
    CCTK_REAL  trKL = trK[index];
    CCTK_REAL  Xt1L = Xt1[index];
    CCTK_REAL  Xt2L = Xt2[index];
    CCTK_REAL  Xt3L = Xt3[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    int dir1 = Sign(normal[0]);
    
    int dir2 = Sign(normal[1]);
    
    int dir3 = Sign(normal[2]);
    
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu12 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu13 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu23 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL em4phi = IfThen(conformalMethod,SQR(phiL),exp(-4*phiL));
    
    CCTK_REAL gu11 = em4phi*gtu11;
    
    CCTK_REAL gu12 = em4phi*gtu12;
    
    CCTK_REAL gu13 = em4phi*gtu13;
    
    CCTK_REAL gu22 = em4phi*gtu22;
    
    CCTK_REAL gu23 = em4phi*gtu23;
    
    CCTK_REAL gu33 = em4phi*gtu33;
    
    CCTK_REAL nn1 = normal[0];
    
    CCTK_REAL nn2 = normal[1];
    
    CCTK_REAL nn3 = normal[2];
    
    CCTK_REAL nu1 = gu11*nn1 + gu12*nn2 + gu13*nn3;
    
    CCTK_REAL nu2 = gu12*nn1 + gu22*nn2 + gu23*nn3;
    
    CCTK_REAL nu3 = gu13*nn1 + gu23*nn2 + gu33*nn3;
    
    CCTK_REAL nlen2 = nn1*nu1 + nn2*nu2 + nn3*nu3;
    
    CCTK_REAL nlen = pow(nlen2,0.5);
    
    CCTK_REAL su1 = nu1*INV(nlen);
    
    CCTK_REAL su2 = nu2*INV(nlen);
    
    CCTK_REAL su3 = nu3*INV(nlen);
    
    CCTK_REAL vg = pow(harmonicF,0.5);
    
    CCTK_REAL phirhsL = -((PDonesided1(phi, i, j, k)*su1 + 
      PDonesided2(phi, i, j, k)*su2 + PDonesided3(phi, i, j, k)*su3)*vg);
    
    CCTK_REAL gt11rhsL = -(PDonesided1(gt11, i, j, k)*su1) - 
      PDonesided2(gt11, i, j, k)*su2 - PDonesided3(gt11, i, j, k)*su3;
    
    CCTK_REAL gt12rhsL = -(PDonesided1(gt12, i, j, k)*su1) - 
      PDonesided2(gt12, i, j, k)*su2 - PDonesided3(gt12, i, j, k)*su3;
    
    CCTK_REAL gt13rhsL = -(PDonesided1(gt13, i, j, k)*su1) - 
      PDonesided2(gt13, i, j, k)*su2 - PDonesided3(gt13, i, j, k)*su3;
    
    CCTK_REAL gt22rhsL = -(PDonesided1(gt22, i, j, k)*su1) - 
      PDonesided2(gt22, i, j, k)*su2 - PDonesided3(gt22, i, j, k)*su3;
    
    CCTK_REAL gt23rhsL = -(PDonesided1(gt23, i, j, k)*su1) - 
      PDonesided2(gt23, i, j, k)*su2 - PDonesided3(gt23, i, j, k)*su3;
    
    CCTK_REAL gt33rhsL = -(PDonesided1(gt33, i, j, k)*su1) - 
      PDonesided2(gt33, i, j, k)*su2 - PDonesided3(gt33, i, j, k)*su3;
    
    CCTK_REAL trKrhsL = -((PDonesided1(trK, i, j, k)*su1 + 
      PDonesided2(trK, i, j, k)*su2 + PDonesided3(trK, i, j, k)*su3)*vg);
    
    CCTK_REAL At11rhsL = -(PDonesided1(At11, i, j, k)*su1) - 
      PDonesided2(At11, i, j, k)*su2 - PDonesided3(At11, i, j, k)*su3;
    
    CCTK_REAL At12rhsL = -(PDonesided1(At12, i, j, k)*su1) - 
      PDonesided2(At12, i, j, k)*su2 - PDonesided3(At12, i, j, k)*su3;
    
    CCTK_REAL At13rhsL = -(PDonesided1(At13, i, j, k)*su1) - 
      PDonesided2(At13, i, j, k)*su2 - PDonesided3(At13, i, j, k)*su3;
    
    CCTK_REAL At22rhsL = -(PDonesided1(At22, i, j, k)*su1) - 
      PDonesided2(At22, i, j, k)*su2 - PDonesided3(At22, i, j, k)*su3;
    
    CCTK_REAL At23rhsL = -(PDonesided1(At23, i, j, k)*su1) - 
      PDonesided2(At23, i, j, k)*su2 - PDonesided3(At23, i, j, k)*su3;
    
    CCTK_REAL At33rhsL = -(PDonesided1(At33, i, j, k)*su1) - 
      PDonesided2(At33, i, j, k)*su2 - PDonesided3(At33, i, j, k)*su3;
    
    CCTK_REAL Xt1rhsL = -(PDonesided1(Xt1, i, j, k)*su1) - 
      PDonesided2(Xt1, i, j, k)*su2 - PDonesided3(Xt1, i, j, k)*su3;
    
    CCTK_REAL Xt2rhsL = -(PDonesided1(Xt2, i, j, k)*su1) - 
      PDonesided2(Xt2, i, j, k)*su2 - PDonesided3(Xt2, i, j, k)*su3;
    
    CCTK_REAL Xt3rhsL = -(PDonesided1(Xt3, i, j, k)*su1) - 
      PDonesided2(Xt3, i, j, k)*su2 - PDonesided3(Xt3, i, j, k)*su3;
    
    CCTK_REAL alpharhsL = -((PDonesided1(alpha, i, j, k)*su1 + 
      PDonesided2(alpha, i, j, k)*su2 + PDonesided3(alpha, i, j, 
      k)*su3)*vg);
    
    CCTK_REAL ArhsL = -((PDonesided1(A, i, j, k)*su1 + PDonesided2(A, 
      i, j, k)*su2 + PDonesided3(A, i, j, k)*su3)*vg);
    
    CCTK_REAL beta1rhsL = -(PDonesided1(beta1, i, j, k)*su1) - 
      PDonesided2(beta1, i, j, k)*su2 - PDonesided3(beta1, i, j, k)*su3;
    
    CCTK_REAL beta2rhsL = -(PDonesided1(beta2, i, j, k)*su1) - 
      PDonesided2(beta2, i, j, k)*su2 - PDonesided3(beta2, i, j, k)*su3;
    
    CCTK_REAL beta3rhsL = -(PDonesided1(beta3, i, j, k)*su1) - 
      PDonesided2(beta3, i, j, k)*su2 - PDonesided3(beta3, i, j, k)*su3;
    
    CCTK_REAL B1rhsL = -(PDonesided1(B1, i, j, k)*su1) - 
      PDonesided2(B1, i, j, k)*su2 - PDonesided3(B1, i, j, k)*su3;
    
    CCTK_REAL B2rhsL = -(PDonesided1(B2, i, j, k)*su1) - 
      PDonesided2(B2, i, j, k)*su2 - PDonesided3(B2, i, j, k)*su3;
    
    CCTK_REAL B3rhsL = -(PDonesided1(B3, i, j, k)*su1) - 
      PDonesided2(B3, i, j, k)*su2 - PDonesided3(B3, i, j, k)*su3;
    
    
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
  }
  LC_ENDLOOP3 (ML_BSSN_O2_RHSRadiativeBoundary);
}

void ML_BSSN_O2_RHSRadiativeBoundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverBoundary(cctkGH, &ML_BSSN_O2_RHSRadiativeBoundary_Body);
}
