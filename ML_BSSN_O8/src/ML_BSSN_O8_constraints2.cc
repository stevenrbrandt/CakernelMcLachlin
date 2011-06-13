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
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define QAD(x) (SQR(SQR(x)))
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))

extern "C" void ML_BSSN_O8_constraints2_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_cons_detg","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_cons_detg.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_cons_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_cons_Gamma.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_cons_traceA","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_cons_traceA.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_mom","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_mom.");
  return;
}

static void ML_BSSN_O8_constraints2_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O8_constraints2_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O8_constraints2_calc_every != ML_BSSN_O8_constraints2_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_BSSN_O8::ML_cons_detg","ML_BSSN_O8::ML_cons_Gamma","ML_BSSN_O8::ML_cons_traceA","ML_BSSN_O8::ML_curv","ML_BSSN_O8::ML_Gamma","ML_BSSN_O8::ML_lapse","ML_BSSN_O8::ML_log_confac","ML_BSSN_O8::ML_metric","ML_BSSN_O8::ML_mom","ML_BSSN_O8::ML_shift","ML_BSSN_O8::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_O8_constraints2", 11, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_O8_constraints2", 4, 4, 4);
  
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
  CCTK_REAL const p1o1024dx = 0.0009765625*INV(dx);
  CCTK_REAL const p1o1024dy = 0.0009765625*INV(dy);
  CCTK_REAL const p1o1024dz = 0.0009765625*INV(dz);
  CCTK_REAL const p1o1680dx = 0.000595238095238095238095238095238*INV(dx);
  CCTK_REAL const p1o1680dy = 0.000595238095238095238095238095238*INV(dy);
  CCTK_REAL const p1o1680dz = 0.000595238095238095238095238095238*INV(dz);
  CCTK_REAL const p1o5040dx2 = 0.000198412698412698412698412698413*INV(SQR(dx));
  CCTK_REAL const p1o5040dy2 = 0.000198412698412698412698412698413*INV(SQR(dy));
  CCTK_REAL const p1o5040dz2 = 0.000198412698412698412698412698413*INV(SQR(dz));
  CCTK_REAL const p1o560dx = 0.00178571428571428571428571428571*INV(dx);
  CCTK_REAL const p1o560dy = 0.00178571428571428571428571428571*INV(dy);
  CCTK_REAL const p1o560dz = 0.00178571428571428571428571428571*INV(dz);
  CCTK_REAL const p1o705600dxdy = 1.41723356009070294784580498866e-6*INV(dx)*INV(dy);
  CCTK_REAL const p1o705600dxdz = 1.41723356009070294784580498866e-6*INV(dx)*INV(dz);
  CCTK_REAL const p1o705600dydz = 1.41723356009070294784580498866e-6*INV(dy)*INV(dz);
  CCTK_REAL const p1o840dx = 0.00119047619047619047619047619048*INV(dx);
  CCTK_REAL const p1o840dy = 0.00119047619047619047619047619048*INV(dy);
  CCTK_REAL const p1o840dz = 0.00119047619047619047619047619048*INV(dz);
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o840dx = -0.00119047619047619047619047619048*INV(dx);
  CCTK_REAL const pm1o840dy = -0.00119047619047619047619047619048*INV(dy);
  CCTK_REAL const pm1o840dz = -0.00119047619047619047619047619048*INV(dz);
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_O8_constraints2,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
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
    CCTK_REAL Xt1L = Xt1[index];
    CCTK_REAL Xt2L = Xt2[index];
    CCTK_REAL Xt3L = Xt3[index];
    
    CCTK_REAL eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL;
    
    if (*stress_energy_state)
    {
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
    CCTK_REAL const PDstandardNth1At11 = PDstandardNth1(&At11[index]);
    CCTK_REAL const PDstandardNth2At11 = PDstandardNth2(&At11[index]);
    CCTK_REAL const PDstandardNth3At11 = PDstandardNth3(&At11[index]);
    CCTK_REAL const PDstandardNth1At12 = PDstandardNth1(&At12[index]);
    CCTK_REAL const PDstandardNth2At12 = PDstandardNth2(&At12[index]);
    CCTK_REAL const PDstandardNth3At12 = PDstandardNth3(&At12[index]);
    CCTK_REAL const PDstandardNth1At13 = PDstandardNth1(&At13[index]);
    CCTK_REAL const PDstandardNth2At13 = PDstandardNth2(&At13[index]);
    CCTK_REAL const PDstandardNth3At13 = PDstandardNth3(&At13[index]);
    CCTK_REAL const PDstandardNth1At22 = PDstandardNth1(&At22[index]);
    CCTK_REAL const PDstandardNth2At22 = PDstandardNth2(&At22[index]);
    CCTK_REAL const PDstandardNth3At22 = PDstandardNth3(&At22[index]);
    CCTK_REAL const PDstandardNth1At23 = PDstandardNth1(&At23[index]);
    CCTK_REAL const PDstandardNth2At23 = PDstandardNth2(&At23[index]);
    CCTK_REAL const PDstandardNth3At23 = PDstandardNth3(&At23[index]);
    CCTK_REAL const PDstandardNth1At33 = PDstandardNth1(&At33[index]);
    CCTK_REAL const PDstandardNth2At33 = PDstandardNth2(&At33[index]);
    CCTK_REAL const PDstandardNth3At33 = PDstandardNth3(&At33[index]);
    CCTK_REAL const PDstandardNth1gt11 = PDstandardNth1(&gt11[index]);
    CCTK_REAL const PDstandardNth2gt11 = PDstandardNth2(&gt11[index]);
    CCTK_REAL const PDstandardNth3gt11 = PDstandardNth3(&gt11[index]);
    CCTK_REAL const PDstandardNth1gt12 = PDstandardNth1(&gt12[index]);
    CCTK_REAL const PDstandardNth2gt12 = PDstandardNth2(&gt12[index]);
    CCTK_REAL const PDstandardNth3gt12 = PDstandardNth3(&gt12[index]);
    CCTK_REAL const PDstandardNth1gt13 = PDstandardNth1(&gt13[index]);
    CCTK_REAL const PDstandardNth2gt13 = PDstandardNth2(&gt13[index]);
    CCTK_REAL const PDstandardNth3gt13 = PDstandardNth3(&gt13[index]);
    CCTK_REAL const PDstandardNth1gt22 = PDstandardNth1(&gt22[index]);
    CCTK_REAL const PDstandardNth2gt22 = PDstandardNth2(&gt22[index]);
    CCTK_REAL const PDstandardNth3gt22 = PDstandardNth3(&gt22[index]);
    CCTK_REAL const PDstandardNth1gt23 = PDstandardNth1(&gt23[index]);
    CCTK_REAL const PDstandardNth2gt23 = PDstandardNth2(&gt23[index]);
    CCTK_REAL const PDstandardNth3gt23 = PDstandardNth3(&gt23[index]);
    CCTK_REAL const PDstandardNth1gt33 = PDstandardNth1(&gt33[index]);
    CCTK_REAL const PDstandardNth2gt33 = PDstandardNth2(&gt33[index]);
    CCTK_REAL const PDstandardNth3gt33 = PDstandardNth3(&gt33[index]);
    CCTK_REAL const PDstandardNth1phi = PDstandardNth1(&phi[index]);
    CCTK_REAL const PDstandardNth2phi = PDstandardNth2(&phi[index]);
    CCTK_REAL const PDstandardNth3phi = PDstandardNth3(&phi[index]);
    CCTK_REAL const PDstandardNth1trK = PDstandardNth1(&trK[index]);
    CCTK_REAL const PDstandardNth2trK = PDstandardNth2(&trK[index]);
    CCTK_REAL const PDstandardNth3trK = PDstandardNth3(&trK[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu12 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu13 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu23 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
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
    
    CCTK_REAL fac1 = IfThen(conformalMethod,-0.5*INV(phiL),1);
    
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
  LC_ENDLOOP3 (ML_BSSN_O8_constraints2);
}

extern "C" void ML_BSSN_O8_constraints2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_O8_constraints2_Body);
}
