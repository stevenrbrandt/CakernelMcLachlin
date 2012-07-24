/*  File produced by Kranc */

#define KRANC_C

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"
#include "Differencing.h"
#include "cctk_Loop.h"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define INV(x) ((CCTK_REAL)1.0 / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * SQR(x))
#define QAD(x) (SQR(SQR(x)))

extern "C" void ML_ADMQuantities_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMQuantities::ML_Jadm","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMQuantities::ML_Jadm.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMQuantities::ML_Madm","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMQuantities::ML_Madm.");
  return;
}

static void ML_ADMQuantities_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  ptrdiff_t const di = 1;
  ptrdiff_t const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const cdi = sizeof(CCTK_REAL) * di;
  ptrdiff_t const cdj = sizeof(CCTK_REAL) * dj;
  ptrdiff_t const cdk = sizeof(CCTK_REAL) * dk;
  CCTK_REAL const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL const dz = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL const dt = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL const t = ToReal(cctk_time);
  CCTK_REAL const dxi = INV(dx);
  CCTK_REAL const dyi = INV(dy);
  CCTK_REAL const dzi = INV(dz);
  CCTK_REAL const khalf = 0.5;
  CCTK_REAL const kthird = 1/3.0;
  CCTK_REAL const ktwothird = 2.0/3.0;
  CCTK_REAL const kfourthird = 4.0/3.0;
  CCTK_REAL const keightthird = 8.0/3.0;
  CCTK_REAL const hdxi = 0.5 * dxi;
  CCTK_REAL const hdyi = 0.5 * dyi;
  CCTK_REAL const hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  CCTK_REAL const p1o12dx = 0.0833333333333333333333333333333*INV(dx);
  CCTK_REAL const p1o12dy = 0.0833333333333333333333333333333*INV(dy);
  CCTK_REAL const p1o12dz = 0.0833333333333333333333333333333*INV(dz);
  CCTK_REAL const p1o144dxdy = 0.00694444444444444444444444444444*INV(dx*dy);
  CCTK_REAL const p1o144dxdz = 0.00694444444444444444444444444444*INV(dx*dz);
  CCTK_REAL const p1o144dydz = 0.00694444444444444444444444444444*INV(dy*dz);
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o12dx2 = -0.0833333333333333333333333333333*INV(SQR(dx));
  CCTK_REAL const pm1o12dy2 = -0.0833333333333333333333333333333*INV(SQR(dy));
  CCTK_REAL const pm1o12dz2 = -0.0833333333333333333333333333333*INV(SQR(dz));
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  CCTK_LOOP3(ML_ADMQuantities,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL alphaL = alpha[index];
    CCTK_REAL At11L = At11[index];
    CCTK_REAL At12L = At12[index];
    CCTK_REAL At13L = At13[index];
    CCTK_REAL At22L = At22[index];
    CCTK_REAL At23L = At23[index];
    CCTK_REAL At33L = At33[index];
    CCTK_REAL beta1L = beta1[index];
    CCTK_REAL beta2L = beta2[index];
    CCTK_REAL beta3L = beta3[index];
    CCTK_REAL gt11L = gt11[index];
    CCTK_REAL gt12L = gt12[index];
    CCTK_REAL gt13L = gt13[index];
    CCTK_REAL gt22L = gt22[index];
    CCTK_REAL gt23L = gt23[index];
    CCTK_REAL gt33L = gt33[index];
    CCTK_REAL phiL = phi[index];
    CCTK_REAL trKL = trK[index];
    CCTK_REAL xL = x[index];
    CCTK_REAL Xt1L = Xt1[index];
    CCTK_REAL Xt2L = Xt2[index];
    CCTK_REAL Xt3L = Xt3[index];
    CCTK_REAL yL = y[index];
    CCTK_REAL zL = z[index];
    
    CCTK_REAL eTttL, eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL;
    
    if (*stress_energy_state)
    {
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
    }
    else
    {
      eTttL = ToReal(0.0);
      eTtxL = ToReal(0.0);
      eTtyL = ToReal(0.0);
      eTtzL = ToReal(0.0);
      eTxxL = ToReal(0.0);
      eTxyL = ToReal(0.0);
      eTxzL = ToReal(0.0);
      eTyyL = ToReal(0.0);
      eTyzL = ToReal(0.0);
      eTzzL = ToReal(0.0);
    }
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDstandardNth1gt11 = PDstandardNth1(&gt11[index]);
    CCTK_REAL const PDstandardNth2gt11 = PDstandardNth2(&gt11[index]);
    CCTK_REAL const PDstandardNth3gt11 = PDstandardNth3(&gt11[index]);
    CCTK_REAL const PDstandardNth11gt11 = PDstandardNth11(&gt11[index]);
    CCTK_REAL const PDstandardNth22gt11 = PDstandardNth22(&gt11[index]);
    CCTK_REAL const PDstandardNth33gt11 = PDstandardNth33(&gt11[index]);
    CCTK_REAL const PDstandardNth12gt11 = PDstandardNth12(&gt11[index]);
    CCTK_REAL const PDstandardNth13gt11 = PDstandardNth13(&gt11[index]);
    CCTK_REAL const PDstandardNth23gt11 = PDstandardNth23(&gt11[index]);
    CCTK_REAL const PDstandardNth1gt12 = PDstandardNth1(&gt12[index]);
    CCTK_REAL const PDstandardNth2gt12 = PDstandardNth2(&gt12[index]);
    CCTK_REAL const PDstandardNth3gt12 = PDstandardNth3(&gt12[index]);
    CCTK_REAL const PDstandardNth11gt12 = PDstandardNth11(&gt12[index]);
    CCTK_REAL const PDstandardNth22gt12 = PDstandardNth22(&gt12[index]);
    CCTK_REAL const PDstandardNth33gt12 = PDstandardNth33(&gt12[index]);
    CCTK_REAL const PDstandardNth12gt12 = PDstandardNth12(&gt12[index]);
    CCTK_REAL const PDstandardNth13gt12 = PDstandardNth13(&gt12[index]);
    CCTK_REAL const PDstandardNth23gt12 = PDstandardNth23(&gt12[index]);
    CCTK_REAL const PDstandardNth1gt13 = PDstandardNth1(&gt13[index]);
    CCTK_REAL const PDstandardNth2gt13 = PDstandardNth2(&gt13[index]);
    CCTK_REAL const PDstandardNth3gt13 = PDstandardNth3(&gt13[index]);
    CCTK_REAL const PDstandardNth11gt13 = PDstandardNth11(&gt13[index]);
    CCTK_REAL const PDstandardNth22gt13 = PDstandardNth22(&gt13[index]);
    CCTK_REAL const PDstandardNth33gt13 = PDstandardNth33(&gt13[index]);
    CCTK_REAL const PDstandardNth12gt13 = PDstandardNth12(&gt13[index]);
    CCTK_REAL const PDstandardNth13gt13 = PDstandardNth13(&gt13[index]);
    CCTK_REAL const PDstandardNth23gt13 = PDstandardNth23(&gt13[index]);
    CCTK_REAL const PDstandardNth1gt22 = PDstandardNth1(&gt22[index]);
    CCTK_REAL const PDstandardNth2gt22 = PDstandardNth2(&gt22[index]);
    CCTK_REAL const PDstandardNth3gt22 = PDstandardNth3(&gt22[index]);
    CCTK_REAL const PDstandardNth11gt22 = PDstandardNth11(&gt22[index]);
    CCTK_REAL const PDstandardNth22gt22 = PDstandardNth22(&gt22[index]);
    CCTK_REAL const PDstandardNth33gt22 = PDstandardNth33(&gt22[index]);
    CCTK_REAL const PDstandardNth12gt22 = PDstandardNth12(&gt22[index]);
    CCTK_REAL const PDstandardNth13gt22 = PDstandardNth13(&gt22[index]);
    CCTK_REAL const PDstandardNth23gt22 = PDstandardNth23(&gt22[index]);
    CCTK_REAL const PDstandardNth1gt23 = PDstandardNth1(&gt23[index]);
    CCTK_REAL const PDstandardNth2gt23 = PDstandardNth2(&gt23[index]);
    CCTK_REAL const PDstandardNth3gt23 = PDstandardNth3(&gt23[index]);
    CCTK_REAL const PDstandardNth11gt23 = PDstandardNth11(&gt23[index]);
    CCTK_REAL const PDstandardNth22gt23 = PDstandardNth22(&gt23[index]);
    CCTK_REAL const PDstandardNth33gt23 = PDstandardNth33(&gt23[index]);
    CCTK_REAL const PDstandardNth12gt23 = PDstandardNth12(&gt23[index]);
    CCTK_REAL const PDstandardNth13gt23 = PDstandardNth13(&gt23[index]);
    CCTK_REAL const PDstandardNth23gt23 = PDstandardNth23(&gt23[index]);
    CCTK_REAL const PDstandardNth1gt33 = PDstandardNth1(&gt33[index]);
    CCTK_REAL const PDstandardNth2gt33 = PDstandardNth2(&gt33[index]);
    CCTK_REAL const PDstandardNth3gt33 = PDstandardNth3(&gt33[index]);
    CCTK_REAL const PDstandardNth11gt33 = PDstandardNth11(&gt33[index]);
    CCTK_REAL const PDstandardNth22gt33 = PDstandardNth22(&gt33[index]);
    CCTK_REAL const PDstandardNth33gt33 = PDstandardNth33(&gt33[index]);
    CCTK_REAL const PDstandardNth12gt33 = PDstandardNth12(&gt33[index]);
    CCTK_REAL const PDstandardNth13gt33 = PDstandardNth13(&gt33[index]);
    CCTK_REAL const PDstandardNth23gt33 = PDstandardNth23(&gt33[index]);
    CCTK_REAL const PDstandardNth1trK = PDstandardNth1(&trK[index]);
    CCTK_REAL const PDstandardNth2trK = PDstandardNth2(&trK[index]);
    CCTK_REAL const PDstandardNth3trK = PDstandardNth3(&trK[index]);
    CCTK_REAL const PDstandardNth1Xt1 = PDstandardNth1(&Xt1[index]);
    CCTK_REAL const PDstandardNth2Xt1 = PDstandardNth2(&Xt1[index]);
    CCTK_REAL const PDstandardNth3Xt1 = PDstandardNth3(&Xt1[index]);
    CCTK_REAL const PDstandardNth1Xt2 = PDstandardNth1(&Xt2[index]);
    CCTK_REAL const PDstandardNth2Xt2 = PDstandardNth2(&Xt2[index]);
    CCTK_REAL const PDstandardNth3Xt2 = PDstandardNth3(&Xt2[index]);
    CCTK_REAL const PDstandardNth1Xt3 = PDstandardNth1(&Xt3[index]);
    CCTK_REAL const PDstandardNth2Xt3 = PDstandardNth2(&Xt3[index]);
    CCTK_REAL const PDstandardNth3Xt3 = PDstandardNth3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu21 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu31 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu32 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL dgtu111 = -2*(gtu11*gtu21*PDstandardNth1gt12 + 
      gtu11*gtu31*PDstandardNth1gt13 + gtu21*gtu31*PDstandardNth1gt23) - 
      PDstandardNth1gt11*SQR(gtu11) - PDstandardNth1gt22*SQR(gtu21) - 
      PDstandardNth1gt33*SQR(gtu31);
    
    CCTK_REAL dgtu211 = gtu21*(-(gtu11*PDstandardNth1gt11) - 
      gtu31*PDstandardNth1gt13 - gtu22*PDstandardNth1gt22) + 
      gtu32*(-(gtu11*PDstandardNth1gt13) - gtu21*PDstandardNth1gt23) + 
      gtu31*(-(gtu22*PDstandardNth1gt23) - gtu32*PDstandardNth1gt33) - 
      PDstandardNth1gt12*(gtu11*gtu22 + SQR(gtu21));
    
    CCTK_REAL dgtu311 = -((gtu21*gtu31 + gtu11*gtu32)*PDstandardNth1gt12) 
      + gtu11*(-(gtu31*PDstandardNth1gt11) - gtu33*PDstandardNth1gt13) + 
      gtu32*(-(gtu21*PDstandardNth1gt22) - gtu31*PDstandardNth1gt23) + 
      gtu33*(-(gtu21*PDstandardNth1gt23) - gtu31*PDstandardNth1gt33) - 
      PDstandardNth1gt13*SQR(gtu31);
    
    CCTK_REAL dgtu221 = -2*(gtu21*gtu22*PDstandardNth1gt12 + 
      gtu21*gtu32*PDstandardNth1gt13 + gtu22*gtu32*PDstandardNth1gt23) - 
      PDstandardNth1gt11*SQR(gtu21) - PDstandardNth1gt22*SQR(gtu22) - 
      PDstandardNth1gt33*SQR(gtu32);
    
    CCTK_REAL dgtu321 = -((gtu22*gtu31 + gtu21*gtu32)*PDstandardNth1gt12) 
      + gtu31*(-(gtu21*PDstandardNth1gt11) - gtu32*PDstandardNth1gt13) + 
      gtu33*(-(gtu21*PDstandardNth1gt13) - gtu22*PDstandardNth1gt23) + 
      gtu32*(-(gtu22*PDstandardNth1gt22) - gtu33*PDstandardNth1gt33) - 
      PDstandardNth1gt23*SQR(gtu32);
    
    CCTK_REAL dgtu331 = -2*(gtu31*gtu32*PDstandardNth1gt12 + 
      gtu31*gtu33*PDstandardNth1gt13 + gtu32*gtu33*PDstandardNth1gt23) - 
      PDstandardNth1gt11*SQR(gtu31) - PDstandardNth1gt22*SQR(gtu32) - 
      PDstandardNth1gt33*SQR(gtu33);
    
    CCTK_REAL dgtu112 = -2*(gtu11*gtu21*PDstandardNth2gt12 + 
      gtu11*gtu31*PDstandardNth2gt13 + gtu21*gtu31*PDstandardNth2gt23) - 
      PDstandardNth2gt11*SQR(gtu11) - PDstandardNth2gt22*SQR(gtu21) - 
      PDstandardNth2gt33*SQR(gtu31);
    
    CCTK_REAL dgtu212 = gtu21*(-(gtu11*PDstandardNth2gt11) - 
      gtu31*PDstandardNth2gt13 - gtu22*PDstandardNth2gt22) + 
      gtu32*(-(gtu11*PDstandardNth2gt13) - gtu21*PDstandardNth2gt23) + 
      gtu31*(-(gtu22*PDstandardNth2gt23) - gtu32*PDstandardNth2gt33) - 
      PDstandardNth2gt12*(gtu11*gtu22 + SQR(gtu21));
    
    CCTK_REAL dgtu312 = -((gtu21*gtu31 + gtu11*gtu32)*PDstandardNth2gt12) 
      + gtu11*(-(gtu31*PDstandardNth2gt11) - gtu33*PDstandardNth2gt13) + 
      gtu32*(-(gtu21*PDstandardNth2gt22) - gtu31*PDstandardNth2gt23) + 
      gtu33*(-(gtu21*PDstandardNth2gt23) - gtu31*PDstandardNth2gt33) - 
      PDstandardNth2gt13*SQR(gtu31);
    
    CCTK_REAL dgtu222 = -2*(gtu21*gtu22*PDstandardNth2gt12 + 
      gtu21*gtu32*PDstandardNth2gt13 + gtu22*gtu32*PDstandardNth2gt23) - 
      PDstandardNth2gt11*SQR(gtu21) - PDstandardNth2gt22*SQR(gtu22) - 
      PDstandardNth2gt33*SQR(gtu32);
    
    CCTK_REAL dgtu322 = -((gtu22*gtu31 + gtu21*gtu32)*PDstandardNth2gt12) 
      + gtu31*(-(gtu21*PDstandardNth2gt11) - gtu32*PDstandardNth2gt13) + 
      gtu33*(-(gtu21*PDstandardNth2gt13) - gtu22*PDstandardNth2gt23) + 
      gtu32*(-(gtu22*PDstandardNth2gt22) - gtu33*PDstandardNth2gt33) - 
      PDstandardNth2gt23*SQR(gtu32);
    
    CCTK_REAL dgtu332 = -2*(gtu31*gtu32*PDstandardNth2gt12 + 
      gtu31*gtu33*PDstandardNth2gt13 + gtu32*gtu33*PDstandardNth2gt23) - 
      PDstandardNth2gt11*SQR(gtu31) - PDstandardNth2gt22*SQR(gtu32) - 
      PDstandardNth2gt33*SQR(gtu33);
    
    CCTK_REAL dgtu113 = -2*(gtu11*gtu21*PDstandardNth3gt12 + 
      gtu11*gtu31*PDstandardNth3gt13 + gtu21*gtu31*PDstandardNth3gt23) - 
      PDstandardNth3gt11*SQR(gtu11) - PDstandardNth3gt22*SQR(gtu21) - 
      PDstandardNth3gt33*SQR(gtu31);
    
    CCTK_REAL dgtu213 = gtu21*(-(gtu11*PDstandardNth3gt11) - 
      gtu31*PDstandardNth3gt13 - gtu22*PDstandardNth3gt22) + 
      gtu32*(-(gtu11*PDstandardNth3gt13) - gtu21*PDstandardNth3gt23) + 
      gtu31*(-(gtu22*PDstandardNth3gt23) - gtu32*PDstandardNth3gt33) - 
      PDstandardNth3gt12*(gtu11*gtu22 + SQR(gtu21));
    
    CCTK_REAL dgtu313 = -((gtu21*gtu31 + gtu11*gtu32)*PDstandardNth3gt12) 
      + gtu11*(-(gtu31*PDstandardNth3gt11) - gtu33*PDstandardNth3gt13) + 
      gtu32*(-(gtu21*PDstandardNth3gt22) - gtu31*PDstandardNth3gt23) + 
      gtu33*(-(gtu21*PDstandardNth3gt23) - gtu31*PDstandardNth3gt33) - 
      PDstandardNth3gt13*SQR(gtu31);
    
    CCTK_REAL dgtu223 = -2*(gtu21*gtu22*PDstandardNth3gt12 + 
      gtu21*gtu32*PDstandardNth3gt13 + gtu22*gtu32*PDstandardNth3gt23) - 
      PDstandardNth3gt11*SQR(gtu21) - PDstandardNth3gt22*SQR(gtu22) - 
      PDstandardNth3gt33*SQR(gtu32);
    
    CCTK_REAL dgtu323 = -((gtu22*gtu31 + gtu21*gtu32)*PDstandardNth3gt12) 
      + gtu31*(-(gtu21*PDstandardNth3gt11) - gtu32*PDstandardNth3gt13) + 
      gtu33*(-(gtu21*PDstandardNth3gt13) - gtu22*PDstandardNth3gt23) + 
      gtu32*(-(gtu22*PDstandardNth3gt22) - gtu33*PDstandardNth3gt33) - 
      PDstandardNth3gt23*SQR(gtu32);
    
    CCTK_REAL dgtu333 = -2*(gtu31*gtu32*PDstandardNth3gt12 + 
      gtu31*gtu33*PDstandardNth3gt13 + gtu32*gtu33*PDstandardNth3gt23) - 
      PDstandardNth3gt11*SQR(gtu31) - PDstandardNth3gt22*SQR(gtu32) - 
      PDstandardNth3gt33*SQR(gtu33);
    
    CCTK_REAL Gtl111 = 0.5*PDstandardNth1gt11;
    
    CCTK_REAL Gtl112 = 0.5*PDstandardNth2gt11;
    
    CCTK_REAL Gtl113 = 0.5*PDstandardNth3gt11;
    
    CCTK_REAL Gtl122 = -0.5*PDstandardNth1gt22 + PDstandardNth2gt12;
    
    CCTK_REAL Gtl123 = 0.5*(-PDstandardNth1gt23 + PDstandardNth2gt13 + 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl133 = -0.5*PDstandardNth1gt33 + PDstandardNth3gt13;
    
    CCTK_REAL Gtl211 = PDstandardNth1gt12 - 0.5*PDstandardNth2gt11;
    
    CCTK_REAL Gtl212 = 0.5*PDstandardNth1gt22;
    
    CCTK_REAL Gtl213 = 0.5*(PDstandardNth1gt23 - PDstandardNth2gt13 + 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl222 = 0.5*PDstandardNth2gt22;
    
    CCTK_REAL Gtl223 = 0.5*PDstandardNth3gt22;
    
    CCTK_REAL Gtl233 = -0.5*PDstandardNth2gt33 + PDstandardNth3gt23;
    
    CCTK_REAL Gtl311 = PDstandardNth1gt13 - 0.5*PDstandardNth3gt11;
    
    CCTK_REAL Gtl312 = 0.5*(PDstandardNth1gt23 + PDstandardNth2gt13 - 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl313 = 0.5*PDstandardNth1gt33;
    
    CCTK_REAL Gtl322 = PDstandardNth2gt23 - 0.5*PDstandardNth3gt22;
    
    CCTK_REAL Gtl323 = 0.5*PDstandardNth2gt33;
    
    CCTK_REAL Gtl333 = 0.5*PDstandardNth3gt33;
    
    CCTK_REAL Gtlu111 = Gtl111*gtu11 + Gtl112*gtu21 + Gtl113*gtu31;
    
    CCTK_REAL Gtlu112 = Gtl111*gtu21 + Gtl112*gtu22 + Gtl113*gtu32;
    
    CCTK_REAL Gtlu113 = Gtl111*gtu31 + Gtl112*gtu32 + Gtl113*gtu33;
    
    CCTK_REAL Gtlu121 = Gtl112*gtu11 + Gtl122*gtu21 + Gtl123*gtu31;
    
    CCTK_REAL Gtlu122 = Gtl112*gtu21 + Gtl122*gtu22 + Gtl123*gtu32;
    
    CCTK_REAL Gtlu123 = Gtl112*gtu31 + Gtl122*gtu32 + Gtl123*gtu33;
    
    CCTK_REAL Gtlu131 = Gtl113*gtu11 + Gtl123*gtu21 + Gtl133*gtu31;
    
    CCTK_REAL Gtlu132 = Gtl113*gtu21 + Gtl123*gtu22 + Gtl133*gtu32;
    
    CCTK_REAL Gtlu133 = Gtl113*gtu31 + Gtl123*gtu32 + Gtl133*gtu33;
    
    CCTK_REAL Gtlu211 = Gtl211*gtu11 + Gtl212*gtu21 + Gtl213*gtu31;
    
    CCTK_REAL Gtlu212 = Gtl211*gtu21 + Gtl212*gtu22 + Gtl213*gtu32;
    
    CCTK_REAL Gtlu213 = Gtl211*gtu31 + Gtl212*gtu32 + Gtl213*gtu33;
    
    CCTK_REAL Gtlu221 = Gtl212*gtu11 + Gtl222*gtu21 + Gtl223*gtu31;
    
    CCTK_REAL Gtlu222 = Gtl212*gtu21 + Gtl222*gtu22 + Gtl223*gtu32;
    
    CCTK_REAL Gtlu223 = Gtl212*gtu31 + Gtl222*gtu32 + Gtl223*gtu33;
    
    CCTK_REAL Gtlu231 = Gtl213*gtu11 + Gtl223*gtu21 + Gtl233*gtu31;
    
    CCTK_REAL Gtlu232 = Gtl213*gtu21 + Gtl223*gtu22 + Gtl233*gtu32;
    
    CCTK_REAL Gtlu233 = Gtl213*gtu31 + Gtl223*gtu32 + Gtl233*gtu33;
    
    CCTK_REAL Gtlu311 = Gtl311*gtu11 + Gtl312*gtu21 + Gtl313*gtu31;
    
    CCTK_REAL Gtlu312 = Gtl311*gtu21 + Gtl312*gtu22 + Gtl313*gtu32;
    
    CCTK_REAL Gtlu313 = Gtl311*gtu31 + Gtl312*gtu32 + Gtl313*gtu33;
    
    CCTK_REAL Gtlu321 = Gtl312*gtu11 + Gtl322*gtu21 + Gtl323*gtu31;
    
    CCTK_REAL Gtlu322 = Gtl312*gtu21 + Gtl322*gtu22 + Gtl323*gtu32;
    
    CCTK_REAL Gtlu323 = Gtl312*gtu31 + Gtl322*gtu32 + Gtl323*gtu33;
    
    CCTK_REAL Gtlu331 = Gtl313*gtu11 + Gtl323*gtu21 + Gtl333*gtu31;
    
    CCTK_REAL Gtlu332 = Gtl313*gtu21 + Gtl323*gtu22 + Gtl333*gtu32;
    
    CCTK_REAL Gtlu333 = Gtl313*gtu31 + Gtl323*gtu32 + Gtl333*gtu33;
    
    CCTK_REAL Gt111 = Gtl111*gtu11 + Gtl211*gtu21 + Gtl311*gtu31;
    
    CCTK_REAL Gt211 = Gtl111*gtu21 + Gtl211*gtu22 + Gtl311*gtu32;
    
    CCTK_REAL Gt311 = Gtl111*gtu31 + Gtl211*gtu32 + Gtl311*gtu33;
    
    CCTK_REAL Gt112 = Gtl112*gtu11 + Gtl212*gtu21 + Gtl312*gtu31;
    
    CCTK_REAL Gt212 = Gtl112*gtu21 + Gtl212*gtu22 + Gtl312*gtu32;
    
    CCTK_REAL Gt312 = Gtl112*gtu31 + Gtl212*gtu32 + Gtl312*gtu33;
    
    CCTK_REAL Gt113 = Gtl113*gtu11 + Gtl213*gtu21 + Gtl313*gtu31;
    
    CCTK_REAL Gt213 = Gtl113*gtu21 + Gtl213*gtu22 + Gtl313*gtu32;
    
    CCTK_REAL Gt313 = Gtl113*gtu31 + Gtl213*gtu32 + Gtl313*gtu33;
    
    CCTK_REAL Gt122 = Gtl122*gtu11 + Gtl222*gtu21 + Gtl322*gtu31;
    
    CCTK_REAL Gt222 = Gtl122*gtu21 + Gtl222*gtu22 + Gtl322*gtu32;
    
    CCTK_REAL Gt322 = Gtl122*gtu31 + Gtl222*gtu32 + Gtl322*gtu33;
    
    CCTK_REAL Gt123 = Gtl123*gtu11 + Gtl223*gtu21 + Gtl323*gtu31;
    
    CCTK_REAL Gt223 = Gtl123*gtu21 + Gtl223*gtu22 + Gtl323*gtu32;
    
    CCTK_REAL Gt323 = Gtl123*gtu31 + Gtl223*gtu32 + Gtl323*gtu33;
    
    CCTK_REAL Gt133 = Gtl133*gtu11 + Gtl233*gtu21 + Gtl333*gtu31;
    
    CCTK_REAL Gt233 = Gtl133*gtu21 + Gtl233*gtu22 + Gtl333*gtu32;
    
    CCTK_REAL Gt333 = Gtl133*gtu31 + Gtl233*gtu32 + Gtl333*gtu33;
    
    CCTK_REAL Xtn1 = Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + 
      Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    CCTK_REAL Xtn2 = Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + 
      Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    CCTK_REAL Xtn3 = Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + 
      Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    CCTK_REAL Rt11 = 0.5*(6*(Gt111*Gtlu111 + Gt112*Gtlu112 + 
      Gt113*Gtlu113) + 4*(Gt211*Gtlu121 + Gt212*Gtlu122 + Gt213*Gtlu123 + 
      Gt311*Gtlu131 + Gt312*Gtlu132 + Gt313*Gtlu133) - 
      gtu11*PDstandardNth11gt11 - 2*gtu21*PDstandardNth12gt11 - 
      2*gtu31*PDstandardNth13gt11 + 2*(Gt211*Gtlu211 + Gt212*Gtlu212 + 
      Gt213*Gtlu213 + Gt311*Gtlu311 + Gt312*Gtlu312 + Gt313*Gtlu313 + 
      gt11L*PDstandardNth1Xt1) + 2*gt12L*PDstandardNth1Xt2 + 
      2*gt13L*PDstandardNth1Xt3 - gtu22*PDstandardNth22gt11 - 
      2*gtu32*PDstandardNth23gt11 - gtu33*PDstandardNth33gt11 + 2*Gtl111*Xtn1 
      + 2*Gtl112*Xtn2 + 2*Gtl113*Xtn3);
    
    CCTK_REAL Rt12 = 0.5*(4*(Gt211*Gtlu221 + Gt212*Gtlu222 + 
      Gt213*Gtlu223) + 2*(Gt112*Gtlu111 + Gt122*Gtlu112 + Gt123*Gtlu113 + 
      Gt111*Gtlu121 + Gt212*Gtlu121 + Gt112*Gtlu122 + Gt222*Gtlu122 + 
      Gt113*Gtlu123 + Gt223*Gtlu123 + Gt312*Gtlu131 + Gt322*Gtlu132 + 
      Gt323*Gtlu133 + Gt111*Gtlu211 + Gt112*Gtlu212 + Gt113*Gtlu213 + 
      Gt311*Gtlu231 + Gt312*Gtlu232 + Gt313*Gtlu233 + Gt311*Gtlu321 + 
      Gt312*Gtlu322 + Gt313*Gtlu323) - gtu11*PDstandardNth11gt12 - 
      2*gtu21*PDstandardNth12gt12 - 2*gtu31*PDstandardNth13gt12 + 
      gt12L*PDstandardNth1Xt1 + gt22L*PDstandardNth1Xt2 + 
      gt23L*PDstandardNth1Xt3 - gtu22*PDstandardNth22gt12 - 
      2*gtu32*PDstandardNth23gt12 + gt11L*PDstandardNth2Xt1 + 
      gt12L*PDstandardNth2Xt2 + gt13L*PDstandardNth2Xt3 - 
      gtu33*PDstandardNth33gt12 + Gtl112*Xtn1 + Gtl211*Xtn1 + Gtl122*Xtn2 + 
      Gtl212*Xtn2 + Gtl123*Xtn3 + Gtl213*Xtn3);
    
    CCTK_REAL Rt13 = 0.5*(2*(Gt113*Gtlu111 + Gt123*Gtlu112 + Gt133*Gtlu113 
      + Gt213*Gtlu121 + Gt223*Gtlu122 + Gt233*Gtlu123 + Gt111*Gtlu131 + 
      Gt313*Gtlu131 + Gt112*Gtlu132 + Gt323*Gtlu132 + Gt113*Gtlu133 + 
      Gt333*Gtlu133 + Gt211*Gtlu231 + Gt212*Gtlu232 + Gt213*Gtlu233 + 
      Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*Gtlu313 + Gt211*Gtlu321 + 
      Gt212*Gtlu322 + Gt213*Gtlu323) + 4*(Gt311*Gtlu331 + Gt312*Gtlu332 + 
      Gt313*Gtlu333) - gtu11*PDstandardNth11gt13 - 
      2*gtu21*PDstandardNth12gt13 - 2*gtu31*PDstandardNth13gt13 + 
      gt13L*PDstandardNth1Xt1 + gt23L*PDstandardNth1Xt2 + 
      gt33L*PDstandardNth1Xt3 - gtu22*PDstandardNth22gt13 - 
      2*gtu32*PDstandardNth23gt13 - gtu33*PDstandardNth33gt13 + 
      gt11L*PDstandardNth3Xt1 + gt12L*PDstandardNth3Xt2 + 
      gt13L*PDstandardNth3Xt3 + Gtl113*Xtn1 + Gtl311*Xtn1 + Gtl123*Xtn2 + 
      Gtl312*Xtn2 + Gtl133*Xtn3 + Gtl313*Xtn3);
    
    CCTK_REAL Rt22 = 0.5*(6*(Gt212*Gtlu221 + Gt222*Gtlu222 + 
      Gt223*Gtlu223) + 4*(Gt112*Gtlu211 + Gt122*Gtlu212 + Gt123*Gtlu213 + 
      Gt312*Gtlu231 + Gt322*Gtlu232 + Gt323*Gtlu233) - 
      gtu11*PDstandardNth11gt22 - 2*gtu21*PDstandardNth12gt22 - 
      2*gtu31*PDstandardNth13gt22 - gtu22*PDstandardNth22gt22 - 
      2*gtu32*PDstandardNth23gt22 + 2*(Gt112*Gtlu121 + Gt122*Gtlu122 + 
      Gt123*Gtlu123 + Gt312*Gtlu321 + Gt322*Gtlu322 + Gt323*Gtlu323 + 
      gt12L*PDstandardNth2Xt1) + 2*gt22L*PDstandardNth2Xt2 + 
      2*gt23L*PDstandardNth2Xt3 - gtu33*PDstandardNth33gt22 + 2*Gtl212*Xtn1 
      + 2*Gtl222*Xtn2 + 2*Gtl223*Xtn3);
    
    CCTK_REAL Rt23 = 0.5*(2*(Gt112*Gtlu131 + Gt122*Gtlu132 + Gt123*Gtlu133 
      + Gt113*Gtlu211 + Gt123*Gtlu212 + Gt133*Gtlu213 + Gt213*Gtlu221 + 
      Gt223*Gtlu222 + Gt233*Gtlu223 + Gt212*Gtlu231 + Gt313*Gtlu231 + 
      Gt222*Gtlu232 + Gt323*Gtlu232 + Gt223*Gtlu233 + Gt333*Gtlu233 + 
      Gt112*Gtlu311 + Gt122*Gtlu312 + Gt123*Gtlu313 + Gt212*Gtlu321 + 
      Gt222*Gtlu322 + Gt223*Gtlu323) + 4*(Gt312*Gtlu331 + Gt322*Gtlu332 + 
      Gt323*Gtlu333) - gtu11*PDstandardNth11gt23 - 
      2*gtu21*PDstandardNth12gt23 - 2*gtu31*PDstandardNth13gt23 - 
      gtu22*PDstandardNth22gt23 - 2*gtu32*PDstandardNth23gt23 + 
      gt13L*PDstandardNth2Xt1 + gt23L*PDstandardNth2Xt2 + 
      gt33L*PDstandardNth2Xt3 - gtu33*PDstandardNth33gt23 + 
      gt12L*PDstandardNth3Xt1 + gt22L*PDstandardNth3Xt2 + 
      gt23L*PDstandardNth3Xt3 + Gtl213*Xtn1 + Gtl312*Xtn1 + Gtl223*Xtn2 + 
      Gtl322*Xtn2 + Gtl233*Xtn3 + Gtl323*Xtn3);
    
    CCTK_REAL Rt33 = 0.5*(4*(Gt113*Gtlu311 + Gt123*Gtlu312 + Gt133*Gtlu313 
      + Gt213*Gtlu321 + Gt223*Gtlu322 + Gt233*Gtlu323) + 6*(Gt313*Gtlu331 + 
      Gt323*Gtlu332 + Gt333*Gtlu333) - gtu11*PDstandardNth11gt33 - 
      2*gtu21*PDstandardNth12gt33 - 2*gtu31*PDstandardNth13gt33 - 
      gtu22*PDstandardNth22gt33 - 2*gtu32*PDstandardNth23gt33 - 
      gtu33*PDstandardNth33gt33 + 2*(Gt113*Gtlu131 + Gt123*Gtlu132 + 
      Gt133*Gtlu133 + Gt213*Gtlu231 + Gt223*Gtlu232 + Gt233*Gtlu233 + 
      gt13L*PDstandardNth3Xt1) + 2*gt23L*PDstandardNth3Xt2 + 
      2*gt33L*PDstandardNth3Xt3 + 2*Gtl313*Xtn1 + 2*Gtl323*Xtn2 + 
      2*Gtl333*Xtn3);
    
    CCTK_REAL trRt = gtu11*Rt11 + gtu22*Rt22 + 2*(gtu21*Rt12 + gtu31*Rt13 
      + gtu32*Rt23) + gtu33*Rt33;
    
    CCTK_REAL ephi = 
      IfThen(conformalMethod,INV(sqrt(phiL)),exp(phiL));
    
    CCTK_REAL Atm11 = At11L*gtu11 + At12L*gtu21 + At13L*gtu31;
    
    CCTK_REAL Atm21 = At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    CCTK_REAL Atm31 = At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    CCTK_REAL Atm12 = At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    CCTK_REAL Atm22 = At12L*gtu21 + At22L*gtu22 + At23L*gtu32;
    
    CCTK_REAL Atm32 = At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    CCTK_REAL Atm13 = At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    CCTK_REAL Atm23 = At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    CCTK_REAL Atm33 = At13L*gtu31 + At23L*gtu32 + At33L*gtu33;
    
    CCTK_REAL rho = INV(SQR(alphaL))*(eTttL - 2*(beta2L*eTtyL + 
      beta3L*eTtzL) + 2*(beta1L*(-eTtxL + beta2L*eTxyL + 
      beta3L*eTxzL) + beta2L*beta3L*eTyzL) + eTxxL*SQR(beta1L) 
      + eTyyL*SQR(beta2L) + eTzzL*SQR(beta3L));
    
    CCTK_REAL S1 = (-eTtxL + beta1L*eTxxL + beta2L*eTxyL + 
      beta3L*eTxzL)*INV(alphaL);
    
    CCTK_REAL S2 = (-eTtyL + beta1L*eTxyL + beta2L*eTyyL + 
      beta3L*eTyzL)*INV(alphaL);
    
    CCTK_REAL S3 = (-eTtzL + beta1L*eTxzL + beta2L*eTyzL + 
      beta3L*eTzzL)*INV(alphaL);
    
    CCTK_REAL MadmL = -0.0625*INV(Pi)*((Gt111*Gtlu111 + Gt112*Gtlu112 + 
      Gt113*Gtlu113 + Gt211*Gtlu121 + Gt212*Gtlu122 + Gt213*Gtlu123 + 
      Gt311*Gtlu131 + Gt312*Gtlu132 + Gt313*Gtlu133)*gtu11 + (Gt122*Gtlu112 + 
      Gt123*Gtlu113 + Gt222*Gtlu122 + Gt223*Gtlu123 + Gt322*Gtlu132 + 
      Gt323*Gtlu133 + Gt111*Gtlu211 + Gt112*(Gtlu111 + Gtlu212) + 
      Gt113*Gtlu213 + Gt211*Gtlu221 + Gt212*(Gtlu121 + Gtlu222) + 
      Gt213*Gtlu223 + Gt311*Gtlu231 + Gt312*(Gtlu131 + Gtlu232) + 
      Gt313*Gtlu233)*gtu21 + (Gt112*Gtlu211 + Gt122*Gtlu212 + Gt123*Gtlu213 + 
      Gt212*Gtlu221 + Gt222*Gtlu222 + Gt223*Gtlu223 + Gt312*Gtlu231 + 
      Gt322*Gtlu232 + Gt323*Gtlu233)*gtu22 + (Gt123*Gtlu112 + Gt133*Gtlu113 + 
      Gt223*Gtlu122 + Gt233*Gtlu123 + Gt323*Gtlu132 + Gt333*Gtlu133 + 
      Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*(Gtlu111 + Gtlu313) + 
      Gt211*Gtlu321 + Gt212*Gtlu322 + Gt213*(Gtlu121 + Gtlu323) + 
      Gt311*Gtlu331 + Gt312*Gtlu332 + Gt313*(Gtlu131 + Gtlu333))*gtu31 + 
      (Gt113*Gtlu211 + Gt133*Gtlu213 + Gt213*Gtlu221 + Gt233*Gtlu223 + 
      Gt313*Gtlu231 + Gt333*Gtlu233 + Gt112*Gtlu311 + Gt122*Gtlu312 + 
      Gt123*(Gtlu212 + Gtlu313) + Gt212*Gtlu321 + Gt222*Gtlu322 + 
      Gt223*(Gtlu222 + Gtlu323) + Gt312*Gtlu331 + Gt322*Gtlu332 + 
      Gt323*(Gtlu232 + Gtlu333))*gtu32 + (Gt113*Gtlu311 + Gt123*Gtlu312 + 
      Gt133*Gtlu313 + Gt213*Gtlu321 + Gt223*Gtlu322 + Gt233*Gtlu323 + 
      Gt313*Gtlu331 + Gt323*Gtlu332 + Gt333*Gtlu333)*gtu33 + (-1 + ephi)*trRt 
      - pow(ephi,5)*(2*(Atm12*Atm21 + Atm13*Atm31 + Atm23*Atm32) + 16*Pi*rho 
      - 0.666666666666666666666666666667*SQR(trKL) + SQR(Atm11) + 
      SQR(Atm22) + SQR(Atm33)));
    
    CCTK_REAL Jadm1L = 
      0.0208333333333333333333333333333*(-4*zL*PDstandardNth2trK + 
      4*yL*PDstandardNth3trK + 3*(At11L*(zL*dgtu112 - yL*dgtu113) + 
      At22L*(zL*dgtu222 - yL*dgtu223) + 2*(Atm23 + 
      zL*(At12L*dgtu212 + At13L*dgtu312 + At23L*dgtu322)) - 2*(Atm32 
      + yL*(At12L*dgtu213 + At13L*dgtu313 + At23L*dgtu323)) + 
      At33L*(zL*dgtu332 - yL*dgtu333) + Pi*(-16*zL*S2 + 
      16*yL*S3)))*INV(Pi)*pow(ephi,6);
    
    CCTK_REAL Jadm2L = 
      0.0208333333333333333333333333333*(4*zL*PDstandardNth1trK - 
      4*xL*PDstandardNth3trK + 3*(At11L*(-(zL*dgtu111) + xL*dgtu113) 
      + At22L*(-(zL*dgtu221) + xL*dgtu223) - 2*(Atm13 + 
      zL*(At12L*dgtu211 + At13L*dgtu311 + At23L*dgtu321)) + 2*(Atm31 
      + xL*(At12L*dgtu213 + At13L*dgtu313 + At23L*dgtu323)) + 
      At33L*(-(zL*dgtu331) + xL*dgtu333) + Pi*(16*zL*S1 - 
      16*xL*S3)))*INV(Pi)*pow(ephi,6);
    
    CCTK_REAL Jadm3L = 
      0.0208333333333333333333333333333*(-4*yL*PDstandardNth1trK + 
      4*xL*PDstandardNth2trK + 3*(At11L*(yL*dgtu111 - xL*dgtu112) + 
      At22L*(yL*dgtu221 - xL*dgtu222) + 2*(Atm12 + 
      yL*(At12L*dgtu211 + At13L*dgtu311 + At23L*dgtu321)) - 2*(Atm21 
      + xL*(At12L*dgtu212 + At13L*dgtu312 + At23L*dgtu322)) + 
      At33L*(yL*dgtu331 - xL*dgtu332) + Pi*(-16*yL*S1 + 
      16*xL*S2)))*INV(Pi)*pow(ephi,6);
    
    /* Copy local copies back to grid functions */
    Jadm1[index] = Jadm1L;
    Jadm2[index] = Jadm2L;
    Jadm3[index] = Jadm3L;
    Madm[index] = MadmL;
  }
  CCTK_ENDLOOP3(ML_ADMQuantities);
}

extern "C" void ML_ADMQuantities(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMQuantities_Body");
  }
  
  if (cctk_iteration % ML_ADMQuantities_calc_every != ML_ADMQuantities_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "grid::coordinates",
    "Grid::coordinates",
    "ML_BSSN::ML_curv",
    "ML_BSSN::ML_Gamma",
    "ML_BSSN::ML_lapse",
    "ML_BSSN::ML_log_confac",
    "ML_BSSN::ML_metric",
    "ML_BSSN::ML_shift",
    "ML_BSSN::ML_trace_curv",
    "ML_ADMQuantities::ML_Jadm",
    "ML_ADMQuantities::ML_Madm"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADMQuantities", 11, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_ADMQuantities", 2, 2, 2);
  
  GenericFD_LoopOverInterior(cctkGH, ML_ADMQuantities_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_ADMQuantities_Body");
  }
}
