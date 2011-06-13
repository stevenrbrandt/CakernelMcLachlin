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

extern "C" void ML_BSSN_O2_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
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

static void ML_BSSN_O2_RHS1_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O2_RHS1_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O2_RHS1_calc_every != ML_BSSN_O2_RHS1_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"grid::coordinates","Grid::coordinates","ML_BSSN_O2::ML_curv","ML_BSSN_O2::ML_dtlapse","ML_BSSN_O2::ML_dtlapserhs","ML_BSSN_O2::ML_dtshift","ML_BSSN_O2::ML_dtshiftrhs","ML_BSSN_O2::ML_Gamma","ML_BSSN_O2::ML_Gammarhs","ML_BSSN_O2::ML_lapse","ML_BSSN_O2::ML_lapserhs","ML_BSSN_O2::ML_log_confac","ML_BSSN_O2::ML_log_confacrhs","ML_BSSN_O2::ML_metric","ML_BSSN_O2::ML_metricrhs","ML_BSSN_O2::ML_shift","ML_BSSN_O2::ML_shiftrhs","ML_BSSN_O2::ML_trace_curv","ML_BSSN_O2::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_O2_RHS1", 19, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_O2_RHS1", 1, 1, 1);
  
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
  CCTK_REAL const p1o16dx = 0.0625*INV(dx);
  CCTK_REAL const p1o16dy = 0.0625*INV(dy);
  CCTK_REAL const p1o16dz = 0.0625*INV(dz);
  CCTK_REAL const p1o2dx = 0.5*INV(dx);
  CCTK_REAL const p1o2dy = 0.5*INV(dy);
  CCTK_REAL const p1o2dz = 0.5*INV(dz);
  CCTK_REAL const p1o4dx = 0.25*INV(dx);
  CCTK_REAL const p1o4dxdy = 0.25*INV(dx)*INV(dy);
  CCTK_REAL const p1o4dxdz = 0.25*INV(dx)*INV(dz);
  CCTK_REAL const p1o4dy = 0.25*INV(dy);
  CCTK_REAL const p1o4dydz = 0.25*INV(dy)*INV(dz);
  CCTK_REAL const p1o4dz = 0.25*INV(dz);
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1odx2 = INV(SQR(dx));
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1ody2 = INV(SQR(dy));
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const p1odz2 = INV(SQR(dz));
  CCTK_REAL const pm1o2dx = -0.5*INV(dx);
  CCTK_REAL const pm1o2dy = -0.5*INV(dy);
  CCTK_REAL const pm1o2dz = -0.5*INV(dz);
  CCTK_REAL const pm1o4dx = -0.25*INV(dx);
  CCTK_REAL const pm1o4dy = -0.25*INV(dy);
  CCTK_REAL const pm1o4dz = -0.25*INV(dz);
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_O2_RHS1,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
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
    CCTK_REAL gt11L = gt11[index];
    CCTK_REAL gt12L = gt12[index];
    CCTK_REAL gt13L = gt13[index];
    CCTK_REAL gt22L = gt22[index];
    CCTK_REAL gt23L = gt23[index];
    CCTK_REAL gt33L = gt33[index];
    CCTK_REAL phiL = phi[index];
    CCTK_REAL rL = r[index];
    CCTK_REAL trKL = trK[index];
    CCTK_REAL Xt1L = Xt1[index];
    CCTK_REAL Xt2L = Xt2[index];
    CCTK_REAL Xt3L = Xt3[index];
    
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
    CCTK_REAL const PDstandardNth1alpha = PDstandardNth1(&alpha[index]);
    CCTK_REAL const PDstandardNth2alpha = PDstandardNth2(&alpha[index]);
    CCTK_REAL const PDstandardNth3alpha = PDstandardNth3(&alpha[index]);
    CCTK_REAL const PDstandardNth11alpha = PDstandardNth11(&alpha[index]);
    CCTK_REAL const PDstandardNth22alpha = PDstandardNth22(&alpha[index]);
    CCTK_REAL const PDstandardNth33alpha = PDstandardNth33(&alpha[index]);
    CCTK_REAL const PDstandardNth12alpha = PDstandardNth12(&alpha[index]);
    CCTK_REAL const PDstandardNth13alpha = PDstandardNth13(&alpha[index]);
    CCTK_REAL const PDstandardNth23alpha = PDstandardNth23(&alpha[index]);
    CCTK_REAL const PDstandardNth1beta1 = PDstandardNth1(&beta1[index]);
    CCTK_REAL const PDstandardNth2beta1 = PDstandardNth2(&beta1[index]);
    CCTK_REAL const PDstandardNth3beta1 = PDstandardNth3(&beta1[index]);
    CCTK_REAL const PDstandardNth11beta1 = PDstandardNth11(&beta1[index]);
    CCTK_REAL const PDstandardNth22beta1 = PDstandardNth22(&beta1[index]);
    CCTK_REAL const PDstandardNth33beta1 = PDstandardNth33(&beta1[index]);
    CCTK_REAL const PDstandardNth12beta1 = PDstandardNth12(&beta1[index]);
    CCTK_REAL const PDstandardNth13beta1 = PDstandardNth13(&beta1[index]);
    CCTK_REAL const PDstandardNth23beta1 = PDstandardNth23(&beta1[index]);
    CCTK_REAL const PDstandardNth1beta2 = PDstandardNth1(&beta2[index]);
    CCTK_REAL const PDstandardNth2beta2 = PDstandardNth2(&beta2[index]);
    CCTK_REAL const PDstandardNth3beta2 = PDstandardNth3(&beta2[index]);
    CCTK_REAL const PDstandardNth11beta2 = PDstandardNth11(&beta2[index]);
    CCTK_REAL const PDstandardNth22beta2 = PDstandardNth22(&beta2[index]);
    CCTK_REAL const PDstandardNth33beta2 = PDstandardNth33(&beta2[index]);
    CCTK_REAL const PDstandardNth12beta2 = PDstandardNth12(&beta2[index]);
    CCTK_REAL const PDstandardNth13beta2 = PDstandardNth13(&beta2[index]);
    CCTK_REAL const PDstandardNth23beta2 = PDstandardNth23(&beta2[index]);
    CCTK_REAL const PDstandardNth1beta3 = PDstandardNth1(&beta3[index]);
    CCTK_REAL const PDstandardNth2beta3 = PDstandardNth2(&beta3[index]);
    CCTK_REAL const PDstandardNth3beta3 = PDstandardNth3(&beta3[index]);
    CCTK_REAL const PDstandardNth11beta3 = PDstandardNth11(&beta3[index]);
    CCTK_REAL const PDstandardNth22beta3 = PDstandardNth22(&beta3[index]);
    CCTK_REAL const PDstandardNth33beta3 = PDstandardNth33(&beta3[index]);
    CCTK_REAL const PDstandardNth12beta3 = PDstandardNth12(&beta3[index]);
    CCTK_REAL const PDstandardNth13beta3 = PDstandardNth13(&beta3[index]);
    CCTK_REAL const PDstandardNth23beta3 = PDstandardNth23(&beta3[index]);
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
    
    CCTK_REAL Xtn1 = Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu12 + 
      Gt113*gtu13 + Gt123*gtu23) + Gt133*gtu33;
    
    CCTK_REAL Xtn2 = Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu12 + 
      Gt213*gtu13 + Gt223*gtu23) + Gt233*gtu33;
    
    CCTK_REAL Xtn3 = Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu12 + 
      Gt313*gtu13 + Gt323*gtu23) + Gt333*gtu33;
    
    CCTK_REAL fac1 = IfThen(conformalMethod,-0.5*INV(phiL),1);
    
    CCTK_REAL cdphi1 = fac1*PDstandardNth1phi;
    
    CCTK_REAL cdphi2 = fac1*PDstandardNth2phi;
    
    CCTK_REAL cdphi3 = fac1*PDstandardNth3phi;
    
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
    
    CCTK_REAL e4phi = IfThen(conformalMethod,INV(SQR(phiL)),exp(4*phiL));
    
    CCTK_REAL em4phi = INV(e4phi);
    
    CCTK_REAL rho = INV(SQR(alphaL))*(eTttL - 2*(beta2L*eTtyL + 
      beta3L*eTtzL) + 2*(beta1L*(-eTtxL + beta2L*eTxyL + beta3L*eTxzL) + 
      beta2L*beta3L*eTyzL) + eTxxL*SQR(beta1L) + eTyyL*SQR(beta2L) + 
      eTzzL*SQR(beta3L));
    
    CCTK_REAL S1 = (-eTtxL + beta1L*eTxxL + beta2L*eTxyL + 
      beta3L*eTxzL)*INV(alphaL);
    
    CCTK_REAL S2 = (-eTtyL + beta1L*eTxyL + beta2L*eTyyL + 
      beta3L*eTyzL)*INV(alphaL);
    
    CCTK_REAL S3 = (-eTtzL + beta1L*eTxzL + beta2L*eTyzL + 
      beta3L*eTzzL)*INV(alphaL);
    
    CCTK_REAL trS = em4phi*(eTxxL*gtu11 + eTyyL*gtu22 + 2*(eTxyL*gtu12 + 
      eTxzL*gtu13 + eTyzL*gtu23) + eTzzL*gtu33);
    
    CCTK_REAL phirhsL = 
      IfThen(conformalMethod,phiL*(-0.333333333333333333333333333333*(PDstandardNth1beta1 
      + PDstandardNth2beta2 + PDstandardNth3beta3) + 
      0.333333333333333333333333333333*alphaL*trKL),0.166666666666666666666666666667*(PDstandardNth1beta1 
      + PDstandardNth2beta2 + PDstandardNth3beta3) - 
      0.166666666666666666666666666667*alphaL*trKL);
    
    CCTK_REAL gt11rhsL = -0.666666666666666666666666666667*(3*alphaL*At11L 
      - 3*(gt12L*PDstandardNth1beta2 + gt13L*PDstandardNth1beta3) + 
      gt11L*(-2*PDstandardNth1beta1 + PDstandardNth2beta2 + 
      PDstandardNth3beta3));
    
    CCTK_REAL gt12rhsL = 0.333333333333333333333333333333*(-6*alphaL*At12L 
      + 3*(gt22L*PDstandardNth1beta2 + gt23L*PDstandardNth1beta3 + 
      gt11L*PDstandardNth2beta1 + gt13L*PDstandardNth2beta3) + 
      gt12L*(PDstandardNth1beta1 + PDstandardNth2beta2 - 
      2*PDstandardNth3beta3));
    
    CCTK_REAL gt13rhsL = 0.333333333333333333333333333333*(-6*alphaL*At13L 
      + 3*(gt23L*PDstandardNth1beta2 + gt33L*PDstandardNth1beta3 + 
      gt11L*PDstandardNth3beta1 + gt12L*PDstandardNth3beta2) + 
      gt13L*(PDstandardNth1beta1 - 2*PDstandardNth2beta2 + 
      PDstandardNth3beta3));
    
    CCTK_REAL gt22rhsL = -0.666666666666666666666666666667*(3*alphaL*At22L 
      - 3*(gt12L*PDstandardNth2beta1 + gt23L*PDstandardNth2beta3) + 
      gt22L*(PDstandardNth1beta1 - 2*PDstandardNth2beta2 + 
      PDstandardNth3beta3));
    
    CCTK_REAL gt23rhsL = 0.333333333333333333333333333333*(-6*alphaL*At23L 
      + 3*(gt13L*PDstandardNth2beta1 + gt33L*PDstandardNth2beta3 + 
      gt12L*PDstandardNth3beta1 + gt22L*PDstandardNth3beta2) + 
      gt23L*(-2*PDstandardNth1beta1 + PDstandardNth2beta2 + 
      PDstandardNth3beta3));
    
    CCTK_REAL gt33rhsL = -0.666666666666666666666666666667*(3*alphaL*At33L 
      - 3*(gt13L*PDstandardNth3beta1 + gt23L*PDstandardNth3beta2) + 
      gt33L*(PDstandardNth1beta1 + PDstandardNth2beta2 - 
      2*PDstandardNth3beta3));
    
    CCTK_REAL dotXt1 = 
      0.333333333333333333333333333333*(7*(gtu12*PDstandardNth12beta1 + 
      gtu13*PDstandardNth13beta1) + 6*gtu23*PDstandardNth23beta1 + 
      3*(gtu22*PDstandardNth22beta1 + gtu33*PDstandardNth33beta1) - 
      6*(Atu11*PDstandardNth1alpha + Atu12*PDstandardNth2alpha + 
      Atu13*PDstandardNth3alpha) + gtu11*(4*PDstandardNth11beta1 + 
      PDstandardNth12beta2 + PDstandardNth13beta3 - 
      150.7964473723100754462068823974161384415*alphaL*S1) + 
      gtu12*(PDstandardNth22beta2 + PDstandardNth23beta3 - 
      150.7964473723100754462068823974161384415*alphaL*S2) + 
      gtu13*(PDstandardNth23beta2 + PDstandardNth33beta3 - 
      150.7964473723100754462068823974161384415*alphaL*S3) + 
      (-PDstandardNth1beta1 + 2*PDstandardNth3beta3)*Xtn1 + 
      2*(alphaL*(18*(Atu11*cdphi1 + Atu12*cdphi2 + Atu13*cdphi3) + 
      6*(Atu12*Gt112 + Atu13*Gt113 + Atu23*Gt123) + 3*(Atu11*Gt111 + 
      Atu22*Gt122 + Atu33*Gt133) - 2*(gtu11*PDstandardNth1trK + 
      gtu12*PDstandardNth2trK + gtu13*PDstandardNth3trK)) + 
      PDstandardNth2beta2*Xtn1) - 3*(PDstandardNth2beta1*Xtn2 + 
      PDstandardNth3beta1*Xtn3));
    
    CCTK_REAL dotXt2 = 
      0.333333333333333333333333333333*(6*gtu13*PDstandardNth13beta2 + 
      3*(gtu11*PDstandardNth11beta2 + gtu33*PDstandardNth33beta2) - 
      6*(Atu12*PDstandardNth1alpha + Atu22*PDstandardNth2alpha + 
      Atu23*PDstandardNth3alpha) + gtu12*(PDstandardNth11beta1 + 
      7*PDstandardNth12beta2 + PDstandardNth13beta3 - 
      150.7964473723100754462068823974161384415*alphaL*S1) + 
      gtu22*(PDstandardNth12beta1 + 4*PDstandardNth22beta2 + 
      PDstandardNth23beta3 - 
      150.7964473723100754462068823974161384415*alphaL*S2) + 
      gtu23*(PDstandardNth13beta1 + 7*PDstandardNth23beta2 + 
      PDstandardNth33beta3 - 
      150.7964473723100754462068823974161384415*alphaL*S3) + 
      (-PDstandardNth2beta2 + 2*PDstandardNth3beta3)*Xtn2 + 
      2*(alphaL*(18*(Atu12*cdphi1 + Atu22*cdphi2 + Atu23*cdphi3) + 
      6*(Atu12*Gt212 + Atu13*Gt213 + Atu23*Gt223) + 3*(Atu11*Gt211 + 
      Atu22*Gt222 + Atu33*Gt233) - 2*(gtu12*PDstandardNth1trK + 
      gtu22*PDstandardNth2trK + gtu23*PDstandardNth3trK)) + 
      PDstandardNth1beta1*Xtn2) - 3*(PDstandardNth1beta2*Xtn1 + 
      PDstandardNth3beta2*Xtn3));
    
    CCTK_REAL dotXt3 = 
      0.333333333333333333333333333333*(6*gtu12*PDstandardNth12beta3 + 
      3*(gtu11*PDstandardNth11beta3 + gtu22*PDstandardNth22beta3) - 
      6*(Atu13*PDstandardNth1alpha + Atu23*PDstandardNth2alpha + 
      Atu33*PDstandardNth3alpha) + gtu13*(PDstandardNth11beta1 + 
      PDstandardNth12beta2 + 7*PDstandardNth13beta3 - 
      150.7964473723100754462068823974161384415*alphaL*S1) + 
      gtu23*(PDstandardNth12beta1 + PDstandardNth22beta2 + 
      7*PDstandardNth23beta3 - 
      150.7964473723100754462068823974161384415*alphaL*S2) + 
      gtu33*(PDstandardNth13beta1 + PDstandardNth23beta2 + 
      4*PDstandardNth33beta3 - 
      150.7964473723100754462068823974161384415*alphaL*S3) - 
      3*(PDstandardNth1beta3*Xtn1 + PDstandardNth2beta3*Xtn2) + 
      (2*PDstandardNth2beta2 - PDstandardNth3beta3)*Xtn3 + 
      2*(alphaL*(18*(Atu13*cdphi1 + Atu23*cdphi2 + Atu33*cdphi3) + 
      6*(Atu12*Gt312 + Atu13*Gt313 + Atu23*Gt323) + 3*(Atu11*Gt311 + 
      Atu22*Gt322 + Atu33*Gt333) - 2*(gtu13*PDstandardNth1trK + 
      gtu23*PDstandardNth2trK + gtu33*PDstandardNth3trK)) + 
      PDstandardNth1beta1*Xtn3));
    
    CCTK_REAL Xt1rhsL = dotXt1;
    
    CCTK_REAL Xt2rhsL = dotXt2;
    
    CCTK_REAL Xt3rhsL = dotXt3;
    
    CCTK_REAL dottrK = -(em4phi*(gtu11*PDstandardNth11alpha + 
      gtu22*PDstandardNth22alpha + gtu33*(PDstandardNth33alpha + 
      2*cdphi3*PDstandardNth3alpha) + 2*(gtu12*PDstandardNth12alpha + 
      gtu13*(PDstandardNth13alpha + cdphi1*PDstandardNth3alpha) + 
      gtu23*(PDstandardNth23alpha + cdphi2*PDstandardNth3alpha)) + 
      PDstandardNth1alpha*(2*(cdphi1*gtu11 + cdphi2*gtu12 + cdphi3*gtu13) - 
      Xtn1) + PDstandardNth2alpha*(2*(cdphi1*gtu12 + cdphi2*gtu22 + 
      cdphi3*gtu23) - Xtn2) - PDstandardNth3alpha*Xtn3)) + 
      alphaL*(2*(Atm12*Atm21 + Atm13*Atm31 + Atm23*Atm32) + 
      12.56637061435917295385057353311801153679*(rho + trS) + SQR(Atm11) + 
      SQR(Atm22) + SQR(Atm33) + 0.333333333333333333333333333333*SQR(trKL));
    
    CCTK_REAL trKrhsL = dottrK;
    
    CCTK_REAL alpharhsL = 
      -(pow(alphaL,ToReal(harmonicN))*ToReal(harmonicF)*(trKL + (AL - 
      trKL)*ToReal(LapseACoeff)));
    
    CCTK_REAL ArhsL = (dottrK - 
      AL*ToReal(AlphaDriver))*ToReal(LapseACoeff);
    
    CCTK_REAL eta = fmin(1,INV(rL)*ToReal(SpatialBetaDriverRadius));
    
    CCTK_REAL theta = fmin(1,exp(1 - 
      rL*INV(ToReal(SpatialShiftGammaCoeffRadius))));
    
    CCTK_REAL beta1rhsL = theta*(Xt1L + beta1L*eta*ToReal(BetaDriver)*(-1 
      + ToReal(ShiftBCoeff)) + (B1L - 
      Xt1L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    CCTK_REAL beta2rhsL = theta*(Xt2L + beta2L*eta*ToReal(BetaDriver)*(-1 
      + ToReal(ShiftBCoeff)) + (B2L - 
      Xt2L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    CCTK_REAL beta3rhsL = theta*(Xt3L + beta3L*eta*ToReal(BetaDriver)*(-1 
      + ToReal(ShiftBCoeff)) + (B3L - 
      Xt3L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
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
  LC_ENDLOOP3 (ML_BSSN_O2_RHS1);
}

extern "C" void ML_BSSN_O2_RHS1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_O2_RHS1_Body);
}
