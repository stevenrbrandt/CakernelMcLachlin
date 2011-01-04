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

void ML_BSSN_constraints2_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_cons_detg","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_cons_detg.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_cons_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_cons_Gamma.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_cons_traceA","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_cons_traceA.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_mom","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_mom.");
  return;
}

void ML_BSSN_constraints2_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_constraints2_Body");
  }
  
  if (cctk_iteration % ML_BSSN_constraints2_calc_every != ML_BSSN_constraints2_calc_offset)
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
  CCTK_REAL const p1o12dx = INV(dx)/12.;
  CCTK_REAL const p1o12dy = INV(dy)/12.;
  CCTK_REAL const p1o12dz = INV(dz)/12.;
  CCTK_REAL const p1o144dxdy = (INV(dx)*INV(dy))/144.;
  CCTK_REAL const p1o144dxdz = (INV(dx)*INV(dz))/144.;
  CCTK_REAL const p1o144dydz = (INV(dy)*INV(dz))/144.;
  CCTK_REAL const p1o24dx = INV(dx)/24.;
  CCTK_REAL const p1o24dy = INV(dy)/24.;
  CCTK_REAL const p1o24dz = INV(dz)/24.;
  CCTK_REAL const p1o64dx = INV(dx)/64.;
  CCTK_REAL const p1o64dy = INV(dy)/64.;
  CCTK_REAL const p1o64dz = INV(dz)/64.;
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o12dx2 = -pow(dx,-2)/12.;
  CCTK_REAL const pm1o12dy2 = -pow(dy,-2)/12.;
  CCTK_REAL const pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_constraints2,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL PDstandardNth1At11 = INITVALUE;
    // CCTK_REAL PDstandardNth2At11 = INITVALUE;
    // CCTK_REAL PDstandardNth3At11 = INITVALUE;
    // CCTK_REAL PDstandardNth1At12 = INITVALUE;
    // CCTK_REAL PDstandardNth2At12 = INITVALUE;
    // CCTK_REAL PDstandardNth3At12 = INITVALUE;
    // CCTK_REAL PDstandardNth1At13 = INITVALUE;
    // CCTK_REAL PDstandardNth2At13 = INITVALUE;
    // CCTK_REAL PDstandardNth3At13 = INITVALUE;
    // CCTK_REAL PDstandardNth1At22 = INITVALUE;
    // CCTK_REAL PDstandardNth2At22 = INITVALUE;
    // CCTK_REAL PDstandardNth3At22 = INITVALUE;
    // CCTK_REAL PDstandardNth1At23 = INITVALUE;
    // CCTK_REAL PDstandardNth2At23 = INITVALUE;
    // CCTK_REAL PDstandardNth3At23 = INITVALUE;
    // CCTK_REAL PDstandardNth1At33 = INITVALUE;
    // CCTK_REAL PDstandardNth2At33 = INITVALUE;
    // CCTK_REAL PDstandardNth3At33 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth1phi = INITVALUE;
    // CCTK_REAL PDstandardNth2phi = INITVALUE;
    // CCTK_REAL PDstandardNth3phi = INITVALUE;
    // CCTK_REAL PDstandardNth1trK = INITVALUE;
    // CCTK_REAL PDstandardNth2trK = INITVALUE;
    // CCTK_REAL PDstandardNth3trK = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL  alphaL = alpha[index];
    CCTK_REAL  At11L = At11[index];
    CCTK_REAL  At12L = At12[index];
    CCTK_REAL  At13L = At13[index];
    CCTK_REAL  At22L = At22[index];
    CCTK_REAL  At23L = At23[index];
    CCTK_REAL  At33L = At33[index];
    CCTK_REAL  beta1L = beta1[index];
    CCTK_REAL  beta2L = beta2[index];
    CCTK_REAL  beta3L = beta3[index];
    CCTK_REAL  eTtxL = (*stress_energy_state) ? (eTtx[index]) : 0.0;
    CCTK_REAL  eTtyL = (*stress_energy_state) ? (eTty[index]) : 0.0;
    CCTK_REAL  eTtzL = (*stress_energy_state) ? (eTtz[index]) : 0.0;
    CCTK_REAL  eTxxL = (*stress_energy_state) ? (eTxx[index]) : 0.0;
    CCTK_REAL  eTxyL = (*stress_energy_state) ? (eTxy[index]) : 0.0;
    CCTK_REAL  eTxzL = (*stress_energy_state) ? (eTxz[index]) : 0.0;
    CCTK_REAL  eTyyL = (*stress_energy_state) ? (eTyy[index]) : 0.0;
    CCTK_REAL  eTyzL = (*stress_energy_state) ? (eTyz[index]) : 0.0;
    CCTK_REAL  eTzzL = (*stress_energy_state) ? (eTzz[index]) : 0.0;
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
    CCTK_REAL const PDstandardNth1At11 = PDstandardNth1(At11, i, j, k);
    CCTK_REAL const PDstandardNth2At11 = PDstandardNth2(At11, i, j, k);
    CCTK_REAL const PDstandardNth3At11 = PDstandardNth3(At11, i, j, k);
    CCTK_REAL const PDstandardNth1At12 = PDstandardNth1(At12, i, j, k);
    CCTK_REAL const PDstandardNth2At12 = PDstandardNth2(At12, i, j, k);
    CCTK_REAL const PDstandardNth3At12 = PDstandardNth3(At12, i, j, k);
    CCTK_REAL const PDstandardNth1At13 = PDstandardNth1(At13, i, j, k);
    CCTK_REAL const PDstandardNth2At13 = PDstandardNth2(At13, i, j, k);
    CCTK_REAL const PDstandardNth3At13 = PDstandardNth3(At13, i, j, k);
    CCTK_REAL const PDstandardNth1At22 = PDstandardNth1(At22, i, j, k);
    CCTK_REAL const PDstandardNth2At22 = PDstandardNth2(At22, i, j, k);
    CCTK_REAL const PDstandardNth3At22 = PDstandardNth3(At22, i, j, k);
    CCTK_REAL const PDstandardNth1At23 = PDstandardNth1(At23, i, j, k);
    CCTK_REAL const PDstandardNth2At23 = PDstandardNth2(At23, i, j, k);
    CCTK_REAL const PDstandardNth3At23 = PDstandardNth3(At23, i, j, k);
    CCTK_REAL const PDstandardNth1At33 = PDstandardNth1(At33, i, j, k);
    CCTK_REAL const PDstandardNth2At33 = PDstandardNth2(At33, i, j, k);
    CCTK_REAL const PDstandardNth3At33 = PDstandardNth3(At33, i, j, k);
    CCTK_REAL const PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    CCTK_REAL const PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    CCTK_REAL const PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    CCTK_REAL const PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    CCTK_REAL const PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    CCTK_REAL const PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    CCTK_REAL const PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    CCTK_REAL const PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    CCTK_REAL const PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    CCTK_REAL const PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    CCTK_REAL const PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    CCTK_REAL const PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    CCTK_REAL const PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    CCTK_REAL const PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    CCTK_REAL const PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    CCTK_REAL const PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    CCTK_REAL const PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    CCTK_REAL const PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    CCTK_REAL const PDstandardNth1phi = PDstandardNth1(phi, i, j, k);
    CCTK_REAL const PDstandardNth2phi = PDstandardNth2(phi, i, j, k);
    CCTK_REAL const PDstandardNth3phi = PDstandardNth3(phi, i, j, k);
    CCTK_REAL const PDstandardNth1trK = PDstandardNth1(trK, i, j, k);
    CCTK_REAL const PDstandardNth2trK = PDstandardNth2(trK, i, j, k);
    CCTK_REAL const PDstandardNth3trK = PDstandardNth3(trK, i, j, k);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu12 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu13 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu23 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL Gtl111 = khalf*PDstandardNth1gt11;
    
    CCTK_REAL Gtl112 = khalf*PDstandardNth2gt11;
    
    CCTK_REAL Gtl113 = khalf*PDstandardNth3gt11;
    
    CCTK_REAL Gtl122 = -(khalf*PDstandardNth1gt22) + PDstandardNth2gt12;
    
    CCTK_REAL Gtl123 = khalf*(-PDstandardNth1gt23 + PDstandardNth2gt13 + 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl133 = -(khalf*PDstandardNth1gt33) + PDstandardNth3gt13;
    
    CCTK_REAL Gtl211 = PDstandardNth1gt12 - khalf*PDstandardNth2gt11;
    
    CCTK_REAL Gtl212 = khalf*PDstandardNth1gt22;
    
    CCTK_REAL Gtl213 = khalf*(PDstandardNth1gt23 - PDstandardNth2gt13 + 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl222 = khalf*PDstandardNth2gt22;
    
    CCTK_REAL Gtl223 = khalf*PDstandardNth3gt22;
    
    CCTK_REAL Gtl233 = -(khalf*PDstandardNth2gt33) + PDstandardNth3gt23;
    
    CCTK_REAL Gtl311 = PDstandardNth1gt13 - khalf*PDstandardNth3gt11;
    
    CCTK_REAL Gtl312 = khalf*(PDstandardNth1gt23 + PDstandardNth2gt13 - 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl313 = khalf*PDstandardNth1gt33;
    
    CCTK_REAL Gtl322 = PDstandardNth2gt23 - khalf*PDstandardNth3gt22;
    
    CCTK_REAL Gtl323 = khalf*PDstandardNth2gt33;
    
    CCTK_REAL Gtl333 = khalf*PDstandardNth3gt33;
    
    CCTK_REAL Gt111 = Gtl111*gtu11 + Gtl211*gtu12 + Gtl311*gtu13;
    
    CCTK_REAL Gt211 = Gtl111*gtu12 + Gtl211*gtu22 + Gtl311*gtu23;
    
    CCTK_REAL Gt311 = Gtl111*gtu13 + Gtl211*gtu23 + Gtl311*gtu33;
    
    CCTK_REAL Gt112 = Gtl112*gtu11 + Gtl212*gtu12 + Gtl312*gtu13;
    
    CCTK_REAL Gt212 = Gtl112*gtu12 + Gtl212*gtu22 + Gtl312*gtu23;
    
    CCTK_REAL Gt312 = Gtl112*gtu13 + Gtl212*gtu23 + Gtl312*gtu33;
    
    CCTK_REAL Gt113 = Gtl113*gtu11 + Gtl213*gtu12 + Gtl313*gtu13;
    
    CCTK_REAL Gt213 = Gtl113*gtu12 + Gtl213*gtu22 + Gtl313*gtu23;
    
    CCTK_REAL Gt313 = Gtl113*gtu13 + Gtl213*gtu23 + Gtl313*gtu33;
    
    CCTK_REAL Gt122 = Gtl122*gtu11 + Gtl222*gtu12 + Gtl322*gtu13;
    
    CCTK_REAL Gt222 = Gtl122*gtu12 + Gtl222*gtu22 + Gtl322*gtu23;
    
    CCTK_REAL Gt322 = Gtl122*gtu13 + Gtl222*gtu23 + Gtl322*gtu33;
    
    CCTK_REAL Gt123 = Gtl123*gtu11 + Gtl223*gtu12 + Gtl323*gtu13;
    
    CCTK_REAL Gt223 = Gtl123*gtu12 + Gtl223*gtu22 + Gtl323*gtu23;
    
    CCTK_REAL Gt323 = Gtl123*gtu13 + Gtl223*gtu23 + Gtl323*gtu33;
    
    CCTK_REAL Gt133 = Gtl133*gtu11 + Gtl233*gtu12 + Gtl333*gtu13;
    
    CCTK_REAL Gt233 = Gtl133*gtu12 + Gtl233*gtu22 + Gtl333*gtu23;
    
    CCTK_REAL Gt333 = Gtl133*gtu13 + Gtl233*gtu23 + Gtl333*gtu33;
    
    CCTK_REAL fac1 = IfThen(conformalMethod,-(khalf*INV(phiL)),1);
    
    CCTK_REAL cdphi1 = fac1*PDstandardNth1phi;
    
    CCTK_REAL cdphi2 = fac1*PDstandardNth2phi;
    
    CCTK_REAL cdphi3 = fac1*PDstandardNth3phi;
    
    CCTK_REAL S1 = (-eTtxL + beta1L*eTxxL + beta2L*eTxyL + 
      beta3L*eTxzL)*INV(alphaL);
    
    CCTK_REAL S2 = (-eTtyL + beta1L*eTxyL + beta2L*eTyyL + 
      beta3L*eTyzL)*INV(alphaL);
    
    CCTK_REAL S3 = (-eTtzL + beta1L*eTxzL + beta2L*eTyzL + 
      beta3L*eTzzL)*INV(alphaL);
    
    CCTK_REAL M1L = -2.*((At12L*Gt211 + At13L*Gt311)*gtu11 + 
      At11L*Gt123*gtu23) + At12L*((6.*cdphi1 - 1.*Gt111)*gtu12 - 
      3.*Gt213*gtu13 + 6.*(cdphi2*gtu22 + cdphi3*gtu23) - 1.*Gt233*gtu33) - 
      1.*((At22L*Gt212 + At12L*(Gt112 + Gt222) + At23L*Gt312)*gtu22 + 
      (At13L*Gt112 + At12L*Gt113 + At23L*Gt212)*gtu23 + (At13L*Gt113 + 
      At23L*Gt213 + At33L*Gt313)*gtu33) + At11L*((6.*cdphi1 - 2.*Gt111)*gtu11 
      + 6.*(cdphi2*gtu12 + cdphi3*gtu13) - 1.*(Gt122*gtu22 + Gt133*gtu33)) + 
      At13L*((6.*cdphi1 - 1.*Gt111)*gtu13 + 6.*(cdphi2*gtu23 + cdphi3*gtu33) 
      - 1.*(Gt322*gtu22 + Gt333*gtu33)) + gtu11*PDstandardNth1At11 - 
      0.6666666666666666666666666666666666666667*PDstandardNth1trK + 
      gtu12*(-1.*(At22L*Gt211 + At23L*Gt311) - 3.*(At11L*Gt112 + At12L*Gt212 
      + At13L*Gt312) + PDstandardNth1At12 + PDstandardNth2At11) + 
      gtu22*PDstandardNth2At12 + gtu13*(-1.*(At23L*Gt211 + At33L*Gt311) - 
      3.*(At11L*Gt113 + At13L*Gt313) + PDstandardNth1At13 + 
      PDstandardNth3At11) + gtu23*(-1.*(At22L*Gt213 + At33L*Gt312 + 
      At23L*Gt313) - 2.*(At12L*Gt223 + At13L*Gt323) + PDstandardNth2At13 + 
      PDstandardNth3At12) + gtu33*PDstandardNth3At13 - 
      25.13274122871834590770114706623602307358*S1;
    
    CCTK_REAL M2L = At12L*((6.*cdphi1 - 1.*Gt111)*gtu11 + 6.*(cdphi2*gtu12 
      + cdphi3*gtu13) - 2.*Gt122*gtu22 - 3.*Gt123*gtu23 - 1.*Gt133*gtu33) + 
      At22L*((6.*cdphi2 - 2.*Gt222)*gtu22 + 6.*cdphi3*gtu23 - 1.*Gt233*gtu33) 
      + At23L*(-2.*Gt322*gtu22 - 1.*Gt333*gtu33 + 6.*(cdphi2*gtu23 + 
      cdphi3*gtu33)) - 1.*((At11L*Gt112 + At22L*Gt211 + At12L*Gt212 + 
      At23L*Gt311 + At13L*Gt312)*gtu11 + Gt122*(At11L*gtu12 + At13L*gtu23) + 
      (At23L*Gt223 + At33L*Gt323)*gtu33 + At13L*(Gt112*gtu13 + Gt123*gtu33)) 
      + gtu11*PDstandardNth1At12 + gtu12*(At22L*(6.*cdphi1 - 3.*Gt212) + 
      At12L*(-3.*Gt112 - 1.*Gt222) - 3.*At23L*Gt312 - 1.*At13L*Gt322 + 
      PDstandardNth1At22 + PDstandardNth2At12) + gtu22*PDstandardNth2At22 - 
      0.6666666666666666666666666666666666666667*PDstandardNth2trK + 
      gtu13*(-2.*(At12L*Gt113 + At22L*Gt213) + At23L*(6.*cdphi1 - 1.*Gt212 - 
      2.*Gt313) - 1.*(At11L*Gt123 + At12L*Gt223 + At33L*Gt312 + At13L*Gt323) 
      + PDstandardNth1At23 + PDstandardNth3At12) + gtu23*(-1.*(At23L*Gt222 + 
      At33L*Gt322) - 3.*(At22L*Gt223 + At23L*Gt323) + PDstandardNth2At23 + 
      PDstandardNth3At22) + gtu33*PDstandardNth3At23 - 
      25.13274122871834590770114706623602307358*S2;
    
    CCTK_REAL M3L = -1.*((At11L*Gt113 + At23L*Gt211 + At12L*Gt213 + 
      At33L*Gt311)*gtu11 + (At22L*Gt223 + At33L*Gt322 + At23L*Gt323)*gtu22 + 
      At12L*(Gt113*gtu12 + Gt123*gtu22) + Gt133*(At11L*gtu13 + At12L*gtu23)) 
      + At13L*((6.*cdphi1 - 1.*(Gt111 + Gt313))*gtu11 + 6.*(cdphi2*gtu12 + 
      cdphi3*gtu13) - 1.*Gt122*gtu22 - 3.*Gt123*gtu23 - 2.*Gt133*gtu33) + 
      At23L*((6.*cdphi2 - 1.*Gt222)*gtu22 + 6.*cdphi3*gtu23 - 2.*Gt233*gtu33) 
      + gtu11*PDstandardNth1At13 + gtu12*(-2.*(At13L*Gt112 + At33L*Gt312) + 
      At23L*(6.*cdphi1 - 2.*Gt212 - 1.*Gt313) - 1.*(At11L*Gt123 + At22L*Gt213 
      + At12L*Gt223 + At13L*Gt323) + PDstandardNth1At23 + PDstandardNth2At13) 
      + gtu22*PDstandardNth2At23 + gtu13*(-3.*(At13L*Gt113 + At23L*Gt213) + 
      At33L*(6.*cdphi1 - 3.*Gt313) - 1.*(At12L*Gt233 + At13L*Gt333) + 
      PDstandardNth1At33 + PDstandardNth3At13) + gtu23*(-1.*At22L*Gt233 + 
      At33L*(6.*cdphi2 - 3.*Gt323) + At23L*(-3.*Gt223 - 1.*Gt333) + 
      PDstandardNth2At33 + PDstandardNth3At23) + gtu33*(At33L*(6.*cdphi3 - 
      2.*Gt333) + PDstandardNth3At33) - 
      0.6666666666666666666666666666666666666667*PDstandardNth3trK - 
      25.13274122871834590770114706623602307358*S3;
    
    CCTK_REAL cSL = Log(detgt);
    
    CCTK_REAL cXt1L = Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu12 + 
      Gt113*gtu13 + Gt123*gtu23) + Gt133*gtu33 - Xt1L;
    
    CCTK_REAL cXt2L = Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu12 + 
      Gt213*gtu13 + Gt223*gtu23) + Gt233*gtu33 - Xt2L;
    
    CCTK_REAL cXt3L = Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu12 + 
      Gt313*gtu13 + Gt323*gtu23) + Gt333*gtu33 - Xt3L;
    
    CCTK_REAL cAL = At11L*gtu11 + At22L*gtu22 + 2*(At12L*gtu12 + 
      At13L*gtu13 + At23L*gtu23) + At33L*gtu33;
    
    
    /* Copy local copies back to grid functions */
    cA[index] = cAL;
    cS[index] = cSL;
    cXt1[index] = cXt1L;
    cXt2[index] = cXt2L;
    cXt3[index] = cXt3L;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
  }
  LC_ENDLOOP3 (ML_BSSN_constraints2);
}

void ML_BSSN_constraints2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_constraints2_Body);
}
