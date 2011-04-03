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
#define INITVALUE (42)
#define QAD(x) (SQR(SQR(x)))
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))

extern "C" void ML_BSSN_MP_O8_constraints2_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_cons_detg","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_cons_detg.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_cons_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_cons_Gamma.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_cons_traceA","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_cons_traceA.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_mom","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_mom.");
  return;
}

static void ML_BSSN_MP_O8_constraints2_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_O8_constraints2_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_O8_constraints2_calc_every != ML_BSSN_MP_O8_constraints2_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"Coordinates::jacobian","ML_BSSN_MP_O8::ML_cons_detg","ML_BSSN_MP_O8::ML_cons_Gamma","ML_BSSN_MP_O8::ML_cons_traceA","ML_BSSN_MP_O8::ML_curv","ML_BSSN_MP_O8::ML_Gamma","ML_BSSN_MP_O8::ML_lapse","ML_BSSN_MP_O8::ML_log_confac","ML_BSSN_MP_O8::ML_metric","ML_BSSN_MP_O8::ML_mom","ML_BSSN_MP_O8::ML_shift","ML_BSSN_MP_O8::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_O8_constraints2", 12, groups);
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  ptrdiff_t const di = 1;
  ptrdiff_t const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  CCTK_REAL const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL const dz = ToReal(CCTK_DELTA_SPACE(2));
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
  LC_LOOP3 (ML_BSSN_MP_O8_constraints2,
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
    CCTK_REAL eTtxL = (*stress_energy_state) ? eTtx[index] : ToReal(0.0);
    CCTK_REAL eTtyL = (*stress_energy_state) ? eTty[index] : ToReal(0.0);
    CCTK_REAL eTtzL = (*stress_energy_state) ? eTtz[index] : ToReal(0.0);
    CCTK_REAL eTxxL = (*stress_energy_state) ? eTxx[index] : ToReal(0.0);
    CCTK_REAL eTxyL = (*stress_energy_state) ? eTxy[index] : ToReal(0.0);
    CCTK_REAL eTxzL = (*stress_energy_state) ? eTxz[index] : ToReal(0.0);
    CCTK_REAL eTyyL = (*stress_energy_state) ? eTyy[index] : ToReal(0.0);
    CCTK_REAL eTyzL = (*stress_energy_state) ? eTyz[index] : ToReal(0.0);
    CCTK_REAL eTzzL = (*stress_energy_state) ? eTzz[index] : ToReal(0.0);
    CCTK_REAL gt11L = gt11[index];
    CCTK_REAL gt12L = gt12[index];
    CCTK_REAL gt13L = gt13[index];
    CCTK_REAL gt22L = gt22[index];
    CCTK_REAL gt23L = gt23[index];
    CCTK_REAL gt33L = gt33[index];
    CCTK_REAL J11L = J11[index];
    CCTK_REAL J12L = J12[index];
    CCTK_REAL J13L = J13[index];
    CCTK_REAL J21L = J21[index];
    CCTK_REAL J22L = J22[index];
    CCTK_REAL J23L = J23[index];
    CCTK_REAL J31L = J31[index];
    CCTK_REAL J32L = J32[index];
    CCTK_REAL J33L = J33[index];
    CCTK_REAL phiL = phi[index];
    CCTK_REAL trKL = trK[index];
    CCTK_REAL Xt1L = Xt1[index];
    CCTK_REAL Xt2L = Xt2[index];
    CCTK_REAL Xt3L = Xt3[index];
    
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
    
    CCTK_REAL Gtl111 = 0.5*(J11L*PDstandardNth1gt11 + 
      J21L*PDstandardNth2gt11 + J31L*PDstandardNth3gt11);
    
    CCTK_REAL Gtl112 = 0.5*(J12L*PDstandardNth1gt11 + 
      J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11);
    
    CCTK_REAL Gtl113 = 0.5*(J13L*PDstandardNth1gt11 + 
      J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11);
    
    CCTK_REAL Gtl122 = J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12 + 
      J32L*PDstandardNth3gt12 - 0.5*(J11L*PDstandardNth1gt22 + 
      J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22);
    
    CCTK_REAL Gtl123 = 0.5*(J13L*PDstandardNth1gt12 + 
      J12L*PDstandardNth1gt13 - J11L*PDstandardNth1gt23 + 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 - 
      J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 - J31L*PDstandardNth3gt23);
    
    CCTK_REAL Gtl133 = J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13 + 
      J33L*PDstandardNth3gt13 - 0.5*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33);
    
    CCTK_REAL Gtl211 = J11L*PDstandardNth1gt12 + J21L*PDstandardNth2gt12 - 
      0.5*(J12L*PDstandardNth1gt11 + J22L*PDstandardNth2gt11 + 
      J32L*PDstandardNth3gt11) + J31L*PDstandardNth3gt12;
    
    CCTK_REAL Gtl212 = 0.5*(J11L*PDstandardNth1gt22 + 
      J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22);
    
    CCTK_REAL Gtl213 = 0.5*(J13L*PDstandardNth1gt12 - 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 + 
      J23L*PDstandardNth2gt12 - J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 - 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23);
    
    CCTK_REAL Gtl222 = 0.5*(J12L*PDstandardNth1gt22 + 
      J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22);
    
    CCTK_REAL Gtl223 = 0.5*(J13L*PDstandardNth1gt22 + 
      J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22);
    
    CCTK_REAL Gtl233 = J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt23 - 0.5*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33);
    
    CCTK_REAL Gtl311 = J11L*PDstandardNth1gt13 + J21L*PDstandardNth2gt13 - 
      0.5*(J13L*PDstandardNth1gt11 + J23L*PDstandardNth2gt11 + 
      J33L*PDstandardNth3gt11) + J31L*PDstandardNth3gt13;
    
    CCTK_REAL Gtl312 = 0.5*(-(J13L*PDstandardNth1gt12) + 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23);
    
    CCTK_REAL Gtl313 = 0.5*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33);
    
    CCTK_REAL Gtl322 = J12L*PDstandardNth1gt23 + J22L*PDstandardNth2gt23 - 
      0.5*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
      J33L*PDstandardNth3gt22) + J32L*PDstandardNth3gt23;
    
    CCTK_REAL Gtl323 = 0.5*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33);
    
    CCTK_REAL Gtl333 = 0.5*(J13L*PDstandardNth1gt33 + 
      J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33);
    
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
    
    CCTK_REAL fac1 = IfThen(ToReal(conformalMethod),-0.5*INV(phiL),1);
    
    CCTK_REAL cdphi1 = fac1*(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + J31L*PDstandardNth3phi);
    
    CCTK_REAL cdphi2 = fac1*(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + J32L*PDstandardNth3phi);
    
    CCTK_REAL cdphi3 = fac1*(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi);
    
    CCTK_REAL S1 = (-eTtxL + beta1L*eTxxL + beta2L*eTxyL + 
      beta3L*eTxzL)*INV(alphaL);
    
    CCTK_REAL S2 = (-eTtyL + beta1L*eTxyL + beta2L*eTyyL + 
      beta3L*eTyzL)*INV(alphaL);
    
    CCTK_REAL S3 = (-eTtzL + beta1L*eTxzL + beta2L*eTyzL + 
      beta3L*eTzzL)*INV(alphaL);
    
    CCTK_REAL M1L = (-(At22L*Gt212) + At12L*(-Gt112 - Gt222) - 
      At23L*Gt312)*gtu22 + (6*At13L*cdphi2 - At23L*Gt212 - At33L*Gt312)*gtu23 
      + (6*At12L*cdphi3 - At22L*Gt213 - At23L*Gt313)*gtu23 + At12L*((6*cdphi1 
      - Gt111)*gtu12 - 3*Gt213*gtu13 + 6*cdphi2*gtu22 - Gt113*gtu23) - 
      2*(At12L*Gt211*gtu11 + At13L*Gt311*gtu11 + At11L*Gt123*gtu23) + 
      (6*At13L*cdphi3 - At23L*Gt213 - At33L*Gt313)*gtu33 + 
      At13L*(-2*Gt323*gtu23 - Gt113*gtu33) + At11L*(6*cdphi1*gtu11 - 
      2*Gt111*gtu11 + 6*cdphi2*gtu12 + 6*cdphi3*gtu13 - Gt122*gtu22 - 
      Gt133*gtu33) + At12L*(-2*Gt223*gtu23 - Gt233*gtu33) + At13L*((6*cdphi1 
      - Gt111)*gtu13 - Gt322*gtu22 - Gt112*gtu23 - Gt333*gtu33) + (gtu11*J11L 
      + gtu12*J12L + gtu13*J13L)*PDstandardNth1At11 + 
      gtu22*J12L*PDstandardNth1At12 + gtu23*J13L*PDstandardNth1At12 + 
      gtu12*(-3*At11L*Gt112 - At22L*Gt211 - 3*At12L*Gt212 - At23L*Gt311 - 
      3*At13L*Gt312 + J11L*PDstandardNth1At12) + 
      gtu23*J12L*PDstandardNth1At13 + gtu33*J13L*PDstandardNth1At13 + 
      gtu13*(-3*At11L*Gt113 - At23L*Gt211 - At33L*Gt311 - 3*At13L*Gt313 + 
      J11L*PDstandardNth1At13) - 
      0.666666666666666666666666666667*J11L*PDstandardNth1trK + (gtu11*J21L + 
      gtu12*J22L + gtu13*J23L)*PDstandardNth2At11 + 
      gtu12*J21L*PDstandardNth2At12 + gtu22*J22L*PDstandardNth2At12 + 
      gtu23*J23L*PDstandardNth2At12 + gtu13*J21L*PDstandardNth2At13 + 
      gtu23*J22L*PDstandardNth2At13 + gtu33*J23L*PDstandardNth2At13 - 
      0.666666666666666666666666666667*J21L*PDstandardNth2trK + 
      gtu11*J31L*PDstandardNth3At11 + gtu12*J32L*PDstandardNth3At11 + 
      gtu13*J33L*PDstandardNth3At11 + gtu12*J31L*PDstandardNth3At12 + 
      gtu22*J32L*PDstandardNth3At12 + gtu23*J33L*PDstandardNth3At12 + 
      gtu13*J31L*PDstandardNth3At13 + gtu23*J32L*PDstandardNth3At13 + 
      gtu33*J33L*PDstandardNth3At13 - 
      0.666666666666666666666666666667*J31L*PDstandardNth3trK - 
      25.13274122871834590770114706623602307358*S1;
    
    CCTK_REAL M2L = (-(At11L*Gt112) - At22L*Gt211 - At12L*Gt212 - 
      At23L*Gt311 - At13L*Gt312)*gtu11 + (6*At22L*cdphi1 - At11L*Gt122)*gtu12 
      + (6*At12L*cdphi2 - At13L*Gt322)*gtu12 + (6*At23L*cdphi1 - 
      At11L*Gt123)*gtu13 + (6*At12L*cdphi3 - At13L*Gt323)*gtu13 + 
      Gt112*(-3*At12L*gtu12 - At13L*gtu13) + Gt212*(-3*At22L*gtu12 - 
      At23L*gtu13) + Gt312*(-3*At23L*gtu12 - At33L*gtu13) + At12L*((6*cdphi1 
      - Gt111)*gtu11 - Gt222*gtu12 - Gt223*gtu13) + 6*At22L*cdphi2*gtu22 - 
      2*At22L*Gt222*gtu22 + 6*At22L*cdphi3*gtu23 - 3*At22L*Gt223*gtu23 + 
      Gt122*(-2*At12L*gtu22 - At13L*gtu23) + Gt322*(-2*At23L*gtu22 - 
      At33L*gtu23) + 6*At23L*cdphi3*gtu33 + Gt123*(-3*At12L*gtu23 - 
      At13L*gtu33) + Gt323*(-3*At23L*gtu23 - At33L*gtu33) + 
      At12L*(-2*Gt113*gtu13 - Gt133*gtu33) + At23L*(6*cdphi2*gtu23 - 
      Gt223*gtu33) + At22L*(-2*Gt213*gtu13 - Gt233*gtu33) + 
      At23L*(-2*Gt313*gtu13 - Gt222*gtu23 - Gt333*gtu33) + (gtu11*J11L + 
      gtu12*J12L + gtu13*J13L)*PDstandardNth1At12 + 
      gtu12*J11L*PDstandardNth1At22 + gtu22*J12L*PDstandardNth1At22 + 
      gtu23*J13L*PDstandardNth1At22 + gtu13*J11L*PDstandardNth1At23 + 
      gtu23*J12L*PDstandardNth1At23 + gtu33*J13L*PDstandardNth1At23 - 
      0.666666666666666666666666666667*J12L*PDstandardNth1trK + (gtu11*J21L + 
      gtu12*J22L + gtu13*J23L)*PDstandardNth2At12 + 
      gtu12*J21L*PDstandardNth2At22 + gtu22*J22L*PDstandardNth2At22 + 
      gtu23*J23L*PDstandardNth2At22 + gtu13*J21L*PDstandardNth2At23 + 
      gtu23*J22L*PDstandardNth2At23 + gtu33*J23L*PDstandardNth2At23 - 
      0.666666666666666666666666666667*J22L*PDstandardNth2trK + 
      gtu11*J31L*PDstandardNth3At12 + gtu12*J32L*PDstandardNth3At12 + 
      gtu13*J33L*PDstandardNth3At12 + gtu12*J31L*PDstandardNth3At22 + 
      gtu22*J32L*PDstandardNth3At22 + gtu23*J33L*PDstandardNth3At22 + 
      gtu13*J31L*PDstandardNth3At23 + gtu23*J32L*PDstandardNth3At23 + 
      gtu33*J33L*PDstandardNth3At23 - 
      0.666666666666666666666666666667*J32L*PDstandardNth3trK - 
      25.13274122871834590770114706623602307358*S2;
    
    CCTK_REAL M3L = (-(At11L*Gt113) - At23L*Gt211 - At12L*Gt213 - 
      At33L*Gt311 - At13L*Gt313)*gtu11 + (-2*At13L*Gt112 - At22L*Gt213)*gtu12 
      + (6*At23L*cdphi1 + At12L*(-Gt113 - Gt223))*gtu12 + (6*At13L*cdphi2 - 
      At11L*Gt123 - At23L*Gt313)*gtu12 + (6*At33L*cdphi1 - At11L*Gt133)*gtu13 
      + (6*At13L*cdphi3 - At12L*Gt233)*gtu13 - 3*At33L*Gt313*gtu13 + 
      At13L*((6*cdphi1 - Gt111)*gtu11 - Gt323*gtu12 - Gt333*gtu13) + 
      (6*At23L*cdphi2 - At12L*Gt123 - At22L*Gt223)*gtu22 + 
      At13L*(-3*Gt113*gtu13 - Gt122*gtu22) + At33L*(-2*Gt312*gtu12 - 
      Gt322*gtu22) + At23L*(-3*Gt213*gtu13 - Gt323*gtu22) - 
      3*At13L*Gt123*gtu23 + (6*At33L*cdphi2 - At12L*Gt133)*gtu23 - 
      3*At23L*Gt223*gtu23 + (6*At23L*cdphi3 - At22L*Gt233)*gtu23 - 
      3*At33L*Gt323*gtu23 + At23L*(-2*Gt212*gtu12 - Gt222*gtu22 - 
      Gt333*gtu23) + 6*At33L*cdphi3*gtu33 - 2*At13L*Gt133*gtu33 - 
      2*At23L*Gt233*gtu33 - 2*At33L*Gt333*gtu33 + (gtu11*J11L + gtu12*J12L + 
      gtu13*J13L)*PDstandardNth1At13 + gtu12*J11L*PDstandardNth1At23 + 
      gtu22*J12L*PDstandardNth1At23 + gtu23*J13L*PDstandardNth1At23 + 
      gtu13*J11L*PDstandardNth1At33 + gtu23*J12L*PDstandardNth1At33 + 
      gtu33*J13L*PDstandardNth1At33 - 
      0.666666666666666666666666666667*J13L*PDstandardNth1trK + (gtu11*J21L + 
      gtu12*J22L + gtu13*J23L)*PDstandardNth2At13 + 
      gtu12*J21L*PDstandardNth2At23 + gtu22*J22L*PDstandardNth2At23 + 
      gtu23*J23L*PDstandardNth2At23 + gtu13*J21L*PDstandardNth2At33 + 
      gtu23*J22L*PDstandardNth2At33 + gtu33*J23L*PDstandardNth2At33 - 
      0.666666666666666666666666666667*J23L*PDstandardNth2trK + 
      gtu11*J31L*PDstandardNth3At13 + gtu12*J32L*PDstandardNth3At13 + 
      gtu13*J33L*PDstandardNth3At13 + gtu12*J31L*PDstandardNth3At23 + 
      gtu22*J32L*PDstandardNth3At23 + gtu23*J33L*PDstandardNth3At23 + 
      gtu13*J31L*PDstandardNth3At33 + gtu23*J32L*PDstandardNth3At33 + 
      gtu33*J33L*PDstandardNth3At33 - 
      0.666666666666666666666666666667*J33L*PDstandardNth3trK - 
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
  LC_ENDLOOP3 (ML_BSSN_MP_O8_constraints2);
}

extern "C" void ML_BSSN_MP_O8_constraints2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_O8_constraints2_Body);
}
