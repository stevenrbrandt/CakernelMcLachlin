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
#define QAD(x) (SQR(SQR(x)))
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))

extern "C" void HOST_ML_BSSN_RHS_NonDerivatives1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_trace_curvrhs.");
  return;
}

static void HOST_ML_BSSN_RHS_NonDerivatives1_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_REAL const p1o705600dxdy = 1.41723356009070294784580498866e-6*INV(dx*dy);
  CCTK_REAL const p1o705600dxdz = 1.41723356009070294784580498866e-6*INV(dx*dz);
  CCTK_REAL const p1o705600dydz = 1.41723356009070294784580498866e-6*INV(dy*dz);
  CCTK_REAL const p1o840dx = 0.00119047619047619047619047619048*INV(dx);
  CCTK_REAL const p1o840dy = 0.00119047619047619047619047619048*INV(dy);
  CCTK_REAL const p1o840dz = 0.00119047619047619047619047619048*INV(dz);
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o840dx = -0.00119047619047619047619047619048*INV(dx);
  CCTK_REAL const pm1o840dy = -0.00119047619047619047619047619048*INV(dy);
  CCTK_REAL const pm1o840dz = -0.00119047619047619047619047619048*INV(dz);
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  CCTK_LOOP3(HOST_ML_BSSN_RHS_NonDerivatives1,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL AL = A[index];
    CCTK_REAL alphaL = alpha[index];
    CCTK_REAL At11L = At11[index];
    CCTK_REAL At12L = At12[index];
    CCTK_REAL At13L = At13[index];
    CCTK_REAL At22L = At22[index];
    CCTK_REAL At23L = At23[index];
    CCTK_REAL At33L = At33[index];
    CCTK_REAL B1L = B1[index];
    CCTK_REAL B2L = B2[index];
    CCTK_REAL B3L = B3[index];
    CCTK_REAL beta1L = beta1[index];
    CCTK_REAL beta2L = beta2[index];
    CCTK_REAL beta3L = beta3[index];
    CCTK_REAL DPDstandardNthalpha1L = DPDstandardNthalpha1[index];
    CCTK_REAL DPDstandardNthalpha11L = DPDstandardNthalpha11[index];
    CCTK_REAL DPDstandardNthalpha12L = DPDstandardNthalpha12[index];
    CCTK_REAL DPDstandardNthalpha13L = DPDstandardNthalpha13[index];
    CCTK_REAL DPDstandardNthalpha2L = DPDstandardNthalpha2[index];
    CCTK_REAL DPDstandardNthalpha22L = DPDstandardNthalpha22[index];
    CCTK_REAL DPDstandardNthalpha23L = DPDstandardNthalpha23[index];
    CCTK_REAL DPDstandardNthalpha3L = DPDstandardNthalpha3[index];
    CCTK_REAL DPDstandardNthalpha33L = DPDstandardNthalpha33[index];
    CCTK_REAL DPDstandardNthbeta11L = DPDstandardNthbeta11[index];
    CCTK_REAL DPDstandardNthbeta111L = DPDstandardNthbeta111[index];
    CCTK_REAL DPDstandardNthbeta112L = DPDstandardNthbeta112[index];
    CCTK_REAL DPDstandardNthbeta113L = DPDstandardNthbeta113[index];
    CCTK_REAL DPDstandardNthbeta12L = DPDstandardNthbeta12[index];
    CCTK_REAL DPDstandardNthbeta122L = DPDstandardNthbeta122[index];
    CCTK_REAL DPDstandardNthbeta123L = DPDstandardNthbeta123[index];
    CCTK_REAL DPDstandardNthbeta13L = DPDstandardNthbeta13[index];
    CCTK_REAL DPDstandardNthbeta133L = DPDstandardNthbeta133[index];
    CCTK_REAL DPDstandardNthbeta21L = DPDstandardNthbeta21[index];
    CCTK_REAL DPDstandardNthbeta211L = DPDstandardNthbeta211[index];
    CCTK_REAL DPDstandardNthbeta212L = DPDstandardNthbeta212[index];
    CCTK_REAL DPDstandardNthbeta213L = DPDstandardNthbeta213[index];
    CCTK_REAL DPDstandardNthbeta22L = DPDstandardNthbeta22[index];
    CCTK_REAL DPDstandardNthbeta222L = DPDstandardNthbeta222[index];
    CCTK_REAL DPDstandardNthbeta223L = DPDstandardNthbeta223[index];
    CCTK_REAL DPDstandardNthbeta23L = DPDstandardNthbeta23[index];
    CCTK_REAL DPDstandardNthbeta233L = DPDstandardNthbeta233[index];
    CCTK_REAL DPDstandardNthbeta31L = DPDstandardNthbeta31[index];
    CCTK_REAL DPDstandardNthbeta311L = DPDstandardNthbeta311[index];
    CCTK_REAL DPDstandardNthbeta312L = DPDstandardNthbeta312[index];
    CCTK_REAL DPDstandardNthbeta313L = DPDstandardNthbeta313[index];
    CCTK_REAL DPDstandardNthbeta32L = DPDstandardNthbeta32[index];
    CCTK_REAL DPDstandardNthbeta322L = DPDstandardNthbeta322[index];
    CCTK_REAL DPDstandardNthbeta323L = DPDstandardNthbeta323[index];
    CCTK_REAL DPDstandardNthbeta33L = DPDstandardNthbeta33[index];
    CCTK_REAL DPDstandardNthbeta333L = DPDstandardNthbeta333[index];
    CCTK_REAL DPDstandardNthgt111L = DPDstandardNthgt111[index];
    CCTK_REAL DPDstandardNthgt112L = DPDstandardNthgt112[index];
    CCTK_REAL DPDstandardNthgt113L = DPDstandardNthgt113[index];
    CCTK_REAL DPDstandardNthgt121L = DPDstandardNthgt121[index];
    CCTK_REAL DPDstandardNthgt122L = DPDstandardNthgt122[index];
    CCTK_REAL DPDstandardNthgt123L = DPDstandardNthgt123[index];
    CCTK_REAL DPDstandardNthgt131L = DPDstandardNthgt131[index];
    CCTK_REAL DPDstandardNthgt132L = DPDstandardNthgt132[index];
    CCTK_REAL DPDstandardNthgt133L = DPDstandardNthgt133[index];
    CCTK_REAL DPDstandardNthgt221L = DPDstandardNthgt221[index];
    CCTK_REAL DPDstandardNthgt222L = DPDstandardNthgt222[index];
    CCTK_REAL DPDstandardNthgt223L = DPDstandardNthgt223[index];
    CCTK_REAL DPDstandardNthgt231L = DPDstandardNthgt231[index];
    CCTK_REAL DPDstandardNthgt232L = DPDstandardNthgt232[index];
    CCTK_REAL DPDstandardNthgt233L = DPDstandardNthgt233[index];
    CCTK_REAL DPDstandardNthgt331L = DPDstandardNthgt331[index];
    CCTK_REAL DPDstandardNthgt332L = DPDstandardNthgt332[index];
    CCTK_REAL DPDstandardNthgt333L = DPDstandardNthgt333[index];
    CCTK_REAL DPDstandardNthphi1L = DPDstandardNthphi1[index];
    CCTK_REAL DPDstandardNthphi2L = DPDstandardNthphi2[index];
    CCTK_REAL DPDstandardNthphi3L = DPDstandardNthphi3[index];
    CCTK_REAL DPDstandardNthtrK1L = DPDstandardNthtrK1[index];
    CCTK_REAL DPDstandardNthtrK2L = DPDstandardNthtrK2[index];
    CCTK_REAL DPDstandardNthtrK3L = DPDstandardNthtrK3[index];
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
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu12 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu13 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu23 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL Gtl111 = 0.5*DPDstandardNthgt111L;
    
    CCTK_REAL Gtl112 = 0.5*DPDstandardNthgt112L;
    
    CCTK_REAL Gtl113 = 0.5*DPDstandardNthgt113L;
    
    CCTK_REAL Gtl122 = DPDstandardNthgt122L - 
      0.5*DPDstandardNthgt221L;
    
    CCTK_REAL Gtl123 = 0.5*(DPDstandardNthgt123L + 
      DPDstandardNthgt132L - DPDstandardNthgt231L);
    
    CCTK_REAL Gtl133 = DPDstandardNthgt133L - 
      0.5*DPDstandardNthgt331L;
    
    CCTK_REAL Gtl211 = -0.5*DPDstandardNthgt112L + 
      DPDstandardNthgt121L;
    
    CCTK_REAL Gtl212 = 0.5*DPDstandardNthgt221L;
    
    CCTK_REAL Gtl213 = 0.5*(DPDstandardNthgt123L - 
      DPDstandardNthgt132L + DPDstandardNthgt231L);
    
    CCTK_REAL Gtl222 = 0.5*DPDstandardNthgt222L;
    
    CCTK_REAL Gtl223 = 0.5*DPDstandardNthgt223L;
    
    CCTK_REAL Gtl233 = DPDstandardNthgt233L - 
      0.5*DPDstandardNthgt332L;
    
    CCTK_REAL Gtl311 = -0.5*DPDstandardNthgt113L + 
      DPDstandardNthgt131L;
    
    CCTK_REAL Gtl312 = 0.5*(-DPDstandardNthgt123L + 
      DPDstandardNthgt132L + DPDstandardNthgt231L);
    
    CCTK_REAL Gtl313 = 0.5*DPDstandardNthgt331L;
    
    CCTK_REAL Gtl322 = -0.5*DPDstandardNthgt223L + 
      DPDstandardNthgt232L;
    
    CCTK_REAL Gtl323 = 0.5*DPDstandardNthgt332L;
    
    CCTK_REAL Gtl333 = 0.5*DPDstandardNthgt333L;
    
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
    
    CCTK_REAL Xtn1 = Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu12 + 
      Gt113*gtu13 + Gt123*gtu23) + Gt133*gtu33;
    
    CCTK_REAL Xtn2 = Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu12 + 
      Gt213*gtu13 + Gt223*gtu23) + Gt233*gtu33;
    
    CCTK_REAL Xtn3 = Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu12 + 
      Gt313*gtu13 + Gt323*gtu23) + Gt333*gtu33;
    
    CCTK_REAL fac1 = IfThen(conformalMethod,-0.5*INV(phiL),1);
    
    CCTK_REAL cdphi1 = DPDstandardNthphi1L*fac1;
    
    CCTK_REAL cdphi2 = DPDstandardNthphi2L*fac1;
    
    CCTK_REAL cdphi3 = DPDstandardNthphi3L*fac1;
    
    CCTK_REAL Atm11 = At11L*gtu11 + At12L*gtu12 + At13L*gtu13;
    
    CCTK_REAL Atm21 = At11L*gtu12 + At12L*gtu22 + At13L*gtu23;
    
    CCTK_REAL Atm31 = At11L*gtu13 + At12L*gtu23 + At13L*gtu33;
    
    CCTK_REAL Atm12 = At12L*gtu11 + At22L*gtu12 + At23L*gtu13;
    
    CCTK_REAL Atm22 = At12L*gtu12 + At22L*gtu22 + At23L*gtu23;
    
    CCTK_REAL Atm32 = At12L*gtu13 + At22L*gtu23 + At23L*gtu33;
    
    CCTK_REAL Atm13 = At13L*gtu11 + At23L*gtu12 + At33L*gtu13;
    
    CCTK_REAL Atm23 = At13L*gtu12 + At23L*gtu22 + At33L*gtu23;
    
    CCTK_REAL Atm33 = At13L*gtu13 + At23L*gtu23 + At33L*gtu33;
    
    CCTK_REAL Atu11 = Atm11*gtu11 + Atm12*gtu12 + Atm13*gtu13;
    
    CCTK_REAL Atu12 = Atm11*gtu12 + Atm12*gtu22 + Atm13*gtu23;
    
    CCTK_REAL Atu13 = Atm11*gtu13 + Atm12*gtu23 + Atm13*gtu33;
    
    CCTK_REAL Atu22 = Atm21*gtu12 + Atm22*gtu22 + Atm23*gtu23;
    
    CCTK_REAL Atu23 = Atm21*gtu13 + Atm22*gtu23 + Atm23*gtu33;
    
    CCTK_REAL Atu33 = Atm31*gtu13 + Atm32*gtu23 + Atm33*gtu33;
    
    CCTK_REAL e4phi = 
      IfThen(conformalMethod,INV(SQR(phiL)),exp(4*phiL));
    
    CCTK_REAL em4phi = INV(e4phi);
    
    CCTK_REAL phirhsL = 
      IfThen(conformalMethod,phiL*(-0.333333333333333333333333333333*(DPDstandardNthbeta11L 
      + DPDstandardNthbeta22L + DPDstandardNthbeta33L) + 
      0.333333333333333333333333333333*alphaL*trKL),0.166666666666666666666666666667*(DPDstandardNthbeta11L 
      + DPDstandardNthbeta22L + DPDstandardNthbeta33L) - 
      0.166666666666666666666666666667*alphaL*trKL);
    
    CCTK_REAL gt11rhsL = 
      -0.666666666666666666666666666667*(3*alphaL*At11L + 
      (-2*DPDstandardNthbeta11L + DPDstandardNthbeta22L + 
      DPDstandardNthbeta33L)*gt11L - 3*(DPDstandardNthbeta21L*gt12L + 
      DPDstandardNthbeta31L*gt13L));
    
    CCTK_REAL gt12rhsL = 
      0.333333333333333333333333333333*(-6*alphaL*At12L + 
      (DPDstandardNthbeta11L + DPDstandardNthbeta22L - 
      2*DPDstandardNthbeta33L)*gt12L + 3*(DPDstandardNthbeta12L*gt11L 
      + DPDstandardNthbeta32L*gt13L + DPDstandardNthbeta21L*gt22L + 
      DPDstandardNthbeta31L*gt23L));
    
    CCTK_REAL gt13rhsL = 
      0.333333333333333333333333333333*(-6*alphaL*At13L + 
      (DPDstandardNthbeta11L - 2*DPDstandardNthbeta22L + 
      DPDstandardNthbeta33L)*gt13L + 3*(DPDstandardNthbeta13L*gt11L + 
      DPDstandardNthbeta23L*gt12L + DPDstandardNthbeta21L*gt23L + 
      DPDstandardNthbeta31L*gt33L));
    
    CCTK_REAL gt22rhsL = 
      -0.666666666666666666666666666667*(3*alphaL*At22L + 
      (DPDstandardNthbeta11L - 2*DPDstandardNthbeta22L + 
      DPDstandardNthbeta33L)*gt22L - 3*(DPDstandardNthbeta12L*gt12L + 
      DPDstandardNthbeta32L*gt23L));
    
    CCTK_REAL gt23rhsL = 
      0.333333333333333333333333333333*(-6*alphaL*At23L + 
      (-2*DPDstandardNthbeta11L + DPDstandardNthbeta22L + 
      DPDstandardNthbeta33L)*gt23L + 3*(DPDstandardNthbeta13L*gt12L + 
      DPDstandardNthbeta12L*gt13L + DPDstandardNthbeta23L*gt22L + 
      DPDstandardNthbeta32L*gt33L));
    
    CCTK_REAL gt33rhsL = 
      -0.666666666666666666666666666667*(3*alphaL*At33L - 
      3*(DPDstandardNthbeta13L*gt13L + DPDstandardNthbeta23L*gt23L) + 
      (DPDstandardNthbeta11L + DPDstandardNthbeta22L - 
      2*DPDstandardNthbeta33L)*gt33L);
    
    CCTK_REAL dotXt1 = 
      0.333333333333333333333333333333*(-6*(DPDstandardNthalpha1L*Atu11 + 
      DPDstandardNthalpha2L*Atu12 + DPDstandardNthalpha3L*Atu13) + 
      (4*DPDstandardNthbeta111L + DPDstandardNthbeta212L + 
      DPDstandardNthbeta313L)*gtu11 + (7*DPDstandardNthbeta112L + 
      DPDstandardNthbeta222L + DPDstandardNthbeta323L)*gtu12 + 
      (7*DPDstandardNthbeta113L + DPDstandardNthbeta223L + 
      DPDstandardNthbeta333L)*gtu13 + 6*DPDstandardNthbeta123L*gtu23 + 
      3*(DPDstandardNthbeta122L*gtu22 + DPDstandardNthbeta133L*gtu33) + 
      (-DPDstandardNthbeta11L + 2*DPDstandardNthbeta33L)*Xtn1 + 
      2*(alphaL*(Atu13*(18*cdphi3 + 6*Gt113) + 6*(Atu12*(3*cdphi2 + Gt112) 
      + Atu23*Gt123) + 3*(Atu11*(6*cdphi1 + Gt111) + Atu22*Gt122 + 
      Atu33*Gt133) - 2*(DPDstandardNthtrK1L*gtu11 + 
      DPDstandardNthtrK2L*gtu12 + DPDstandardNthtrK3L*gtu13)) + 
      DPDstandardNthbeta22L*Xtn1) - 3*(DPDstandardNthbeta12L*Xtn2 + 
      DPDstandardNthbeta13L*Xtn3));
    
    CCTK_REAL dotXt2 = 
      0.333333333333333333333333333333*(-6*(DPDstandardNthalpha1L*Atu12 + 
      DPDstandardNthalpha2L*Atu22 + DPDstandardNthalpha3L*Atu23) + 
      (DPDstandardNthbeta111L + 7*DPDstandardNthbeta212L + 
      DPDstandardNthbeta313L)*gtu12 + 6*DPDstandardNthbeta213L*gtu13 + 
      (DPDstandardNthbeta112L + 4*DPDstandardNthbeta222L + 
      DPDstandardNthbeta323L)*gtu22 + (DPDstandardNthbeta113L + 
      7*DPDstandardNthbeta223L + DPDstandardNthbeta333L)*gtu23 + 
      3*(DPDstandardNthbeta211L*gtu11 + DPDstandardNthbeta233L*gtu33) + 
      (-DPDstandardNthbeta22L + 2*DPDstandardNthbeta33L)*Xtn2 + 
      2*(alphaL*(6*(Atu12*(3*cdphi1 + Gt212) + Atu13*Gt213) + 
      Atu23*(18*cdphi3 + 6*Gt223) + 3*(Atu11*Gt211 + Atu22*(6*cdphi2 + Gt222) 
      + Atu33*Gt233) - 2*(DPDstandardNthtrK1L*gtu12 + 
      DPDstandardNthtrK2L*gtu22 + DPDstandardNthtrK3L*gtu23)) + 
      DPDstandardNthbeta11L*Xtn2) - 3*(DPDstandardNthbeta21L*Xtn1 + 
      DPDstandardNthbeta23L*Xtn3));
    
    CCTK_REAL dotXt3 = 
      0.333333333333333333333333333333*(-6*(DPDstandardNthalpha1L*Atu13 + 
      DPDstandardNthalpha2L*Atu23 + DPDstandardNthalpha3L*Atu33) + 
      6*DPDstandardNthbeta312L*gtu12 + (DPDstandardNthbeta111L + 
      DPDstandardNthbeta212L + 7*DPDstandardNthbeta313L)*gtu13 + 
      3*(DPDstandardNthbeta311L*gtu11 + DPDstandardNthbeta322L*gtu22) + 
      (DPDstandardNthbeta112L + DPDstandardNthbeta222L + 
      7*DPDstandardNthbeta323L)*gtu23 + (DPDstandardNthbeta113L + 
      DPDstandardNthbeta223L + 4*DPDstandardNthbeta333L)*gtu33 - 
      3*(DPDstandardNthbeta31L*Xtn1 + DPDstandardNthbeta32L*Xtn2) + 
      (2*DPDstandardNthbeta22L - DPDstandardNthbeta33L)*Xtn3 + 
      2*(alphaL*(3*(Atu11*Gt311 + Atu22*Gt322) + 6*(Atu12*Gt312 + 
      Atu13*(3*cdphi1 + Gt313) + Atu23*(3*cdphi2 + Gt323)) + Atu33*(18*cdphi3 
      + 3*Gt333) - 2*(DPDstandardNthtrK1L*gtu13 + 
      DPDstandardNthtrK2L*gtu23 + DPDstandardNthtrK3L*gtu33)) + 
      DPDstandardNthbeta11L*Xtn3));
    
    CCTK_REAL Xt1rhsL = dotXt1;
    
    CCTK_REAL Xt2rhsL = dotXt2;
    
    CCTK_REAL Xt3rhsL = dotXt3;
    
    CCTK_REAL dottrK = -(em4phi*(DPDstandardNthalpha11L*gtu11 + 
      DPDstandardNthalpha22L*gtu22 + DPDstandardNthalpha33L*gtu33 + 
      2*(DPDstandardNthalpha12L*gtu12 + 
      DPDstandardNthalpha1L*cdphi2*gtu12 + DPDstandardNthalpha13L*gtu13 + 
      DPDstandardNthalpha1L*cdphi3*gtu13 + 
      cdphi1*(DPDstandardNthalpha1L*gtu11 + DPDstandardNthalpha2L*gtu12 + 
      DPDstandardNthalpha3L*gtu13) + DPDstandardNthalpha2L*cdphi2*gtu22 + 
      DPDstandardNthalpha23L*gtu23 + DPDstandardNthalpha3L*cdphi2*gtu23 + 
      DPDstandardNthalpha2L*cdphi3*gtu23 + 
      DPDstandardNthalpha3L*cdphi3*gtu33) - DPDstandardNthalpha1L*Xtn1 - 
      DPDstandardNthalpha2L*Xtn2 - DPDstandardNthalpha3L*Xtn3)) + 
      alphaL*(2*(Atm12*Atm21 + Atm13*Atm31 + Atm23*Atm32) + 
      0.333333333333333333333333333333*SQR(trKL) + SQR(Atm11) + SQR(Atm22) 
      + SQR(Atm33));
    
    CCTK_REAL trKrhsL = dottrK;
    
    CCTK_REAL alpharhsL = 
      -(pow(alphaL,ToReal(harmonicN))*ToReal(harmonicF)*((trKL + (-1 + 
      alphaL)*ToReal(AlphaDriver))*(1 - ToReal(LapseACoeff)) + 
      AL*ToReal(LapseACoeff)));
    
    CCTK_REAL ArhsL = (dottrK - 
      AL*ToReal(AlphaDriver))*ToReal(LapseACoeff);
    
    CCTK_REAL eta = 1;
    
    CCTK_REAL theta = 1;
    
    CCTK_REAL beta1rhsL;
    CCTK_REAL beta2rhsL;
    CCTK_REAL beta3rhsL;
    
    if (harmonicShift)
    {
      beta1rhsL = 
        0.5*alphaL*phiL*(-2*phiL*(DPDstandardNthalpha1L*gtu11 + 
        DPDstandardNthalpha2L*gtu12 + DPDstandardNthalpha3L*gtu13) + 
        alphaL*(2*(DPDstandardNthphi2L*gtu12 + DPDstandardNthphi3L*gtu13 
        + gtu11*(DPDstandardNthphi1L + phiL*(DPDstandardNthgt122L*gtu22 + 
        DPDstandardNthgt123L*gtu23 + DPDstandardNthgt133L*gtu33)) + 
        DPDstandardNthgt221L*phiL*SQR(gtu12)) + 
        phiL*((DPDstandardNthgt222L*gtu12 + 
        2*DPDstandardNthgt232L*gtu13)*gtu22 + 
        (-2*DPDstandardNthgt231L*gtu11 + 
        2*DPDstandardNthgt223L*gtu12)*gtu23 + (2*DPDstandardNthgt233L - 
        DPDstandardNthgt332L)*gtu12*gtu33 + gtu11*((DPDstandardNthgt112L + 
        2*DPDstandardNthgt121L)*gtu12 + DPDstandardNthgt113L*gtu13 - 
        DPDstandardNthgt221L*gtu22 + 2*DPDstandardNthgt132L*gtu23 - 
        DPDstandardNthgt331L*gtu33) + gtu13*(4*DPDstandardNthgt231L*gtu12 - 
        DPDstandardNthgt223L*gtu22 + 2*(DPDstandardNthgt131L*gtu11 + 
        DPDstandardNthgt332L*gtu23) + DPDstandardNthgt333L*gtu33) + 
        DPDstandardNthgt111L*SQR(gtu11) + 
        2*DPDstandardNthgt331L*SQR(gtu13))));
      
      beta2rhsL = 
        0.5*alphaL*phiL*(-2*phiL*(DPDstandardNthalpha1L*gtu12 + 
        DPDstandardNthalpha2L*gtu22 + DPDstandardNthalpha3L*gtu23) + 
        alphaL*(gtu12*(phiL*(DPDstandardNthgt111L*gtu11 + 
        2*DPDstandardNthgt122L*gtu22 + 4*DPDstandardNthgt132L*gtu23) + 
        2*(DPDstandardNthphi1L + phiL*(DPDstandardNthgt113L*gtu13 + 
        DPDstandardNthgt133L*gtu33))) + 
        gtu22*(phiL*(-(DPDstandardNthgt112L*gtu11) + 
        DPDstandardNthgt221L*gtu12 + DPDstandardNthgt223L*gtu23) + 
        2*(DPDstandardNthphi2L + phiL*(DPDstandardNthgt121L*gtu11 + 
        DPDstandardNthgt123L*gtu13 + DPDstandardNthgt233L*gtu33))) + 
        phiL*(gtu13*(-2*DPDstandardNthgt132L*gtu22 + 
        2*DPDstandardNthgt231L*gtu22) + (-(DPDstandardNthgt331L*gtu12) - 
        DPDstandardNthgt332L*gtu22)*gtu33 + 
        gtu23*(-(DPDstandardNthgt113L*gtu11) + DPDstandardNthgt333L*gtu33) 
        + DPDstandardNthgt222L*SQR(gtu22)) + 2*((DPDstandardNthphi3L + 
        phiL*(DPDstandardNthgt131L*gtu11 + DPDstandardNthgt331L*gtu13 + 
        DPDstandardNthgt232L*gtu22))*gtu23 + 
        phiL*(DPDstandardNthgt112L*SQR(gtu12) + 
        DPDstandardNthgt332L*SQR(gtu23)))));
      
      beta3rhsL = 
        0.5*alphaL*phiL*(-2*phiL*(DPDstandardNthalpha1L*gtu13 + 
        DPDstandardNthalpha2L*gtu23 + DPDstandardNthalpha3L*gtu33) + 
        alphaL*((2*(DPDstandardNthphi3L + 
        phiL*(DPDstandardNthgt131L*gtu11 + DPDstandardNthgt232L*gtu22)) + 
        phiL*(-(DPDstandardNthgt113L*gtu11) - 
        2*DPDstandardNthgt123L*gtu12 + DPDstandardNthgt331L*gtu13 - 
        DPDstandardNthgt223L*gtu22 + DPDstandardNthgt332L*gtu23))*gtu33 + 
        gtu13*(phiL*(DPDstandardNthgt111L*gtu11 + (2*DPDstandardNthgt122L 
        - DPDstandardNthgt221L)*gtu22 + 4*DPDstandardNthgt123L*gtu23) + 
        2*(DPDstandardNthphi1L + phiL*(DPDstandardNthgt112L*gtu12 + 
        DPDstandardNthgt133L*gtu33))) + 
        gtu23*(phiL*(-(DPDstandardNthgt112L*gtu11) + 
        2*DPDstandardNthgt221L*gtu12 + DPDstandardNthgt222L*gtu22) + 
        2*(DPDstandardNthphi2L + phiL*(DPDstandardNthgt121L*gtu11 + 
        DPDstandardNthgt233L*gtu33))) + phiL*(2*((DPDstandardNthgt132L + 
        DPDstandardNthgt231L)*gtu12*gtu33 + DPDstandardNthgt113L*SQR(gtu13) 
        + DPDstandardNthgt223L*SQR(gtu23)) + 
        DPDstandardNthgt333L*SQR(gtu33))));
    }
    else
    {
      beta1rhsL = theta*(Xt1L + beta1L*eta*ToReal(BetaDriver)*(-1 + 
        ToReal(ShiftBCoeff)) + (B1L - 
        Xt1L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
      
      beta2rhsL = theta*(Xt2L + beta2L*eta*ToReal(BetaDriver)*(-1 + 
        ToReal(ShiftBCoeff)) + (B2L - 
        Xt2L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
      
      beta3rhsL = theta*(Xt3L + beta3L*eta*ToReal(BetaDriver)*(-1 + 
        ToReal(ShiftBCoeff)) + (B3L - 
        Xt3L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    }
    
    CCTK_REAL B1rhsL = (dotXt1 - 
      B1L*eta*ToReal(BetaDriver))*ToReal(ShiftBCoeff);
    
    CCTK_REAL B2rhsL = (dotXt2 - 
      B2L*eta*ToReal(BetaDriver))*ToReal(ShiftBCoeff);
    
    CCTK_REAL B3rhsL = (dotXt3 - 
      B3L*eta*ToReal(BetaDriver))*ToReal(ShiftBCoeff);
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    Arhs[index] = ArhsL;
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
  CCTK_ENDLOOP3(HOST_ML_BSSN_RHS_NonDerivatives1);
}

extern "C" void HOST_ML_BSSN_RHS_NonDerivatives1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering HOST_ML_BSSN_RHS_NonDerivatives1_Body");
  }
  
  if (cctk_iteration % HOST_ML_BSSN_RHS_NonDerivatives1_calc_every != HOST_ML_BSSN_RHS_NonDerivatives1_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "ML_BSSN::DPDstandardNthalpha11_group",
    "ML_BSSN::DPDstandardNthalpha12_group",
    "ML_BSSN::DPDstandardNthalpha13_group",
    "ML_BSSN::DPDstandardNthalpha1_group",
    "ML_BSSN::DPDstandardNthalpha22_group",
    "ML_BSSN::DPDstandardNthalpha23_group",
    "ML_BSSN::DPDstandardNthalpha2_group",
    "ML_BSSN::DPDstandardNthalpha33_group",
    "ML_BSSN::DPDstandardNthalpha3_group",
    "ML_BSSN::DPDstandardNthbeta111_group",
    "ML_BSSN::DPDstandardNthbeta112_group",
    "ML_BSSN::DPDstandardNthbeta113_group",
    "ML_BSSN::DPDstandardNthbeta11_group",
    "ML_BSSN::DPDstandardNthbeta122_group",
    "ML_BSSN::DPDstandardNthbeta123_group",
    "ML_BSSN::DPDstandardNthbeta12_group",
    "ML_BSSN::DPDstandardNthbeta133_group",
    "ML_BSSN::DPDstandardNthbeta13_group",
    "ML_BSSN::DPDstandardNthbeta211_group",
    "ML_BSSN::DPDstandardNthbeta212_group",
    "ML_BSSN::DPDstandardNthbeta213_group",
    "ML_BSSN::DPDstandardNthbeta21_group",
    "ML_BSSN::DPDstandardNthbeta222_group",
    "ML_BSSN::DPDstandardNthbeta223_group",
    "ML_BSSN::DPDstandardNthbeta22_group",
    "ML_BSSN::DPDstandardNthbeta233_group",
    "ML_BSSN::DPDstandardNthbeta23_group",
    "ML_BSSN::DPDstandardNthbeta311_group",
    "ML_BSSN::DPDstandardNthbeta312_group",
    "ML_BSSN::DPDstandardNthbeta313_group",
    "ML_BSSN::DPDstandardNthbeta31_group",
    "ML_BSSN::DPDstandardNthbeta322_group",
    "ML_BSSN::DPDstandardNthbeta323_group",
    "ML_BSSN::DPDstandardNthbeta32_group",
    "ML_BSSN::DPDstandardNthbeta333_group",
    "ML_BSSN::DPDstandardNthbeta33_group",
    "ML_BSSN::DPDstandardNthgt111_group",
    "ML_BSSN::DPDstandardNthgt112_group",
    "ML_BSSN::DPDstandardNthgt113_group",
    "ML_BSSN::DPDstandardNthgt121_group",
    "ML_BSSN::DPDstandardNthgt122_group",
    "ML_BSSN::DPDstandardNthgt123_group",
    "ML_BSSN::DPDstandardNthgt131_group",
    "ML_BSSN::DPDstandardNthgt132_group",
    "ML_BSSN::DPDstandardNthgt133_group",
    "ML_BSSN::DPDstandardNthgt221_group",
    "ML_BSSN::DPDstandardNthgt222_group",
    "ML_BSSN::DPDstandardNthgt223_group",
    "ML_BSSN::DPDstandardNthgt231_group",
    "ML_BSSN::DPDstandardNthgt232_group",
    "ML_BSSN::DPDstandardNthgt233_group",
    "ML_BSSN::DPDstandardNthgt331_group",
    "ML_BSSN::DPDstandardNthgt332_group",
    "ML_BSSN::DPDstandardNthgt333_group",
    "ML_BSSN::DPDstandardNthphi1_group",
    "ML_BSSN::DPDstandardNthphi2_group",
    "ML_BSSN::DPDstandardNthphi3_group",
    "ML_BSSN::DPDstandardNthtrK1_group",
    "ML_BSSN::DPDstandardNthtrK2_group",
    "ML_BSSN::DPDstandardNthtrK3_group",
    "ML_BSSN::ML_curv",
    "ML_BSSN::ML_dtlapse",
    "ML_BSSN::ML_dtlapserhs",
    "ML_BSSN::ML_dtshift",
    "ML_BSSN::ML_dtshiftrhs",
    "ML_BSSN::ML_Gamma",
    "ML_BSSN::ML_Gammarhs",
    "ML_BSSN::ML_lapse",
    "ML_BSSN::ML_lapserhs",
    "ML_BSSN::ML_log_confac",
    "ML_BSSN::ML_log_confacrhs",
    "ML_BSSN::ML_metric",
    "ML_BSSN::ML_metricrhs",
    "ML_BSSN::ML_shift",
    "ML_BSSN::ML_shiftrhs",
    "ML_BSSN::ML_trace_curv",
    "ML_BSSN::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "HOST_ML_BSSN_RHS_NonDerivatives1", 77, groups);
  
  
  GenericFD_LoopOverInterior(cctkGH, HOST_ML_BSSN_RHS_NonDerivatives1_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving HOST_ML_BSSN_RHS_NonDerivatives1_Body");
  }
}
