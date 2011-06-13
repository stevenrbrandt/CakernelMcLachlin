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

extern "C" void ML_BSSN_MP_RHS2_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_curvrhs.");
  return;
}

static void ML_BSSN_MP_RHS2_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_RHS2_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_RHS2_calc_every != ML_BSSN_MP_RHS2_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"Coordinates::jacobian","Coordinates::jacobian2","ML_BSSN_MP::ML_curv","ML_BSSN_MP::ML_curvrhs","ML_BSSN_MP::ML_Gamma","ML_BSSN_MP::ML_lapse","ML_BSSN_MP::ML_log_confac","ML_BSSN_MP::ML_metric","ML_BSSN_MP::ML_shift","ML_BSSN_MP::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_RHS2", 10, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_MP_RHS2", 3, 3, 3);
  
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
  CCTK_REAL const p1o12dx = 0.0833333333333333333333333333333*INV(dx);
  CCTK_REAL const p1o12dy = 0.0833333333333333333333333333333*INV(dy);
  CCTK_REAL const p1o12dz = 0.0833333333333333333333333333333*INV(dz);
  CCTK_REAL const p1o144dxdy = 0.00694444444444444444444444444444*INV(dx)*INV(dy);
  CCTK_REAL const p1o144dxdz = 0.00694444444444444444444444444444*INV(dx)*INV(dz);
  CCTK_REAL const p1o144dydz = 0.00694444444444444444444444444444*INV(dy)*INV(dz);
  CCTK_REAL const p1o24dx = 0.0416666666666666666666666666667*INV(dx);
  CCTK_REAL const p1o24dy = 0.0416666666666666666666666666667*INV(dy);
  CCTK_REAL const p1o24dz = 0.0416666666666666666666666666667*INV(dz);
  CCTK_REAL const p1o64dx = 0.015625*INV(dx);
  CCTK_REAL const p1o64dy = 0.015625*INV(dy);
  CCTK_REAL const p1o64dz = 0.015625*INV(dz);
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o12dx2 = -0.0833333333333333333333333333333*INV(SQR(dx));
  CCTK_REAL const pm1o12dy2 = -0.0833333333333333333333333333333*INV(SQR(dy));
  CCTK_REAL const pm1o12dz2 = -0.0833333333333333333333333333333*INV(SQR(dz));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_MP_RHS2,
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
    CCTK_REAL dJ111L = dJ111[index];
    CCTK_REAL dJ112L = dJ112[index];
    CCTK_REAL dJ113L = dJ113[index];
    CCTK_REAL dJ122L = dJ122[index];
    CCTK_REAL dJ123L = dJ123[index];
    CCTK_REAL dJ133L = dJ133[index];
    CCTK_REAL dJ211L = dJ211[index];
    CCTK_REAL dJ212L = dJ212[index];
    CCTK_REAL dJ213L = dJ213[index];
    CCTK_REAL dJ222L = dJ222[index];
    CCTK_REAL dJ223L = dJ223[index];
    CCTK_REAL dJ233L = dJ233[index];
    CCTK_REAL dJ311L = dJ311[index];
    CCTK_REAL dJ312L = dJ312[index];
    CCTK_REAL dJ313L = dJ313[index];
    CCTK_REAL dJ322L = dJ322[index];
    CCTK_REAL dJ323L = dJ323[index];
    CCTK_REAL dJ333L = dJ333[index];
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
    
    CCTK_REAL eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL;
    
    if (*stress_energy_state)
    {
      eTxxL = eTxx[index];
      eTxyL = eTxy[index];
      eTxzL = eTxz[index];
      eTyyL = eTyy[index];
      eTyzL = eTyz[index];
      eTzzL = eTzz[index];
    }
    else
    {
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
    CCTK_REAL const PDdissipationNth1At11 = PDdissipationNth1(&At11[index]);
    CCTK_REAL const PDdissipationNth2At11 = PDdissipationNth2(&At11[index]);
    CCTK_REAL const PDdissipationNth3At11 = PDdissipationNth3(&At11[index]);
    CCTK_REAL const PDupwindNthAnti1At11 = PDupwindNthAnti1(&At11[index]);
    CCTK_REAL const PDupwindNthSymm1At11 = PDupwindNthSymm1(&At11[index]);
    CCTK_REAL const PDupwindNthAnti2At11 = PDupwindNthAnti2(&At11[index]);
    CCTK_REAL const PDupwindNthSymm2At11 = PDupwindNthSymm2(&At11[index]);
    CCTK_REAL const PDupwindNthAnti3At11 = PDupwindNthAnti3(&At11[index]);
    CCTK_REAL const PDupwindNthSymm3At11 = PDupwindNthSymm3(&At11[index]);
    CCTK_REAL const PDdissipationNth1At12 = PDdissipationNth1(&At12[index]);
    CCTK_REAL const PDdissipationNth2At12 = PDdissipationNth2(&At12[index]);
    CCTK_REAL const PDdissipationNth3At12 = PDdissipationNth3(&At12[index]);
    CCTK_REAL const PDupwindNthAnti1At12 = PDupwindNthAnti1(&At12[index]);
    CCTK_REAL const PDupwindNthSymm1At12 = PDupwindNthSymm1(&At12[index]);
    CCTK_REAL const PDupwindNthAnti2At12 = PDupwindNthAnti2(&At12[index]);
    CCTK_REAL const PDupwindNthSymm2At12 = PDupwindNthSymm2(&At12[index]);
    CCTK_REAL const PDupwindNthAnti3At12 = PDupwindNthAnti3(&At12[index]);
    CCTK_REAL const PDupwindNthSymm3At12 = PDupwindNthSymm3(&At12[index]);
    CCTK_REAL const PDdissipationNth1At13 = PDdissipationNth1(&At13[index]);
    CCTK_REAL const PDdissipationNth2At13 = PDdissipationNth2(&At13[index]);
    CCTK_REAL const PDdissipationNth3At13 = PDdissipationNth3(&At13[index]);
    CCTK_REAL const PDupwindNthAnti1At13 = PDupwindNthAnti1(&At13[index]);
    CCTK_REAL const PDupwindNthSymm1At13 = PDupwindNthSymm1(&At13[index]);
    CCTK_REAL const PDupwindNthAnti2At13 = PDupwindNthAnti2(&At13[index]);
    CCTK_REAL const PDupwindNthSymm2At13 = PDupwindNthSymm2(&At13[index]);
    CCTK_REAL const PDupwindNthAnti3At13 = PDupwindNthAnti3(&At13[index]);
    CCTK_REAL const PDupwindNthSymm3At13 = PDupwindNthSymm3(&At13[index]);
    CCTK_REAL const PDdissipationNth1At22 = PDdissipationNth1(&At22[index]);
    CCTK_REAL const PDdissipationNth2At22 = PDdissipationNth2(&At22[index]);
    CCTK_REAL const PDdissipationNth3At22 = PDdissipationNth3(&At22[index]);
    CCTK_REAL const PDupwindNthAnti1At22 = PDupwindNthAnti1(&At22[index]);
    CCTK_REAL const PDupwindNthSymm1At22 = PDupwindNthSymm1(&At22[index]);
    CCTK_REAL const PDupwindNthAnti2At22 = PDupwindNthAnti2(&At22[index]);
    CCTK_REAL const PDupwindNthSymm2At22 = PDupwindNthSymm2(&At22[index]);
    CCTK_REAL const PDupwindNthAnti3At22 = PDupwindNthAnti3(&At22[index]);
    CCTK_REAL const PDupwindNthSymm3At22 = PDupwindNthSymm3(&At22[index]);
    CCTK_REAL const PDdissipationNth1At23 = PDdissipationNth1(&At23[index]);
    CCTK_REAL const PDdissipationNth2At23 = PDdissipationNth2(&At23[index]);
    CCTK_REAL const PDdissipationNth3At23 = PDdissipationNth3(&At23[index]);
    CCTK_REAL const PDupwindNthAnti1At23 = PDupwindNthAnti1(&At23[index]);
    CCTK_REAL const PDupwindNthSymm1At23 = PDupwindNthSymm1(&At23[index]);
    CCTK_REAL const PDupwindNthAnti2At23 = PDupwindNthAnti2(&At23[index]);
    CCTK_REAL const PDupwindNthSymm2At23 = PDupwindNthSymm2(&At23[index]);
    CCTK_REAL const PDupwindNthAnti3At23 = PDupwindNthAnti3(&At23[index]);
    CCTK_REAL const PDupwindNthSymm3At23 = PDupwindNthSymm3(&At23[index]);
    CCTK_REAL const PDdissipationNth1At33 = PDdissipationNth1(&At33[index]);
    CCTK_REAL const PDdissipationNth2At33 = PDdissipationNth2(&At33[index]);
    CCTK_REAL const PDdissipationNth3At33 = PDdissipationNth3(&At33[index]);
    CCTK_REAL const PDupwindNthAnti1At33 = PDupwindNthAnti1(&At33[index]);
    CCTK_REAL const PDupwindNthSymm1At33 = PDupwindNthSymm1(&At33[index]);
    CCTK_REAL const PDupwindNthAnti2At33 = PDupwindNthAnti2(&At33[index]);
    CCTK_REAL const PDupwindNthSymm2At33 = PDupwindNthSymm2(&At33[index]);
    CCTK_REAL const PDupwindNthAnti3At33 = PDupwindNthAnti3(&At33[index]);
    CCTK_REAL const PDupwindNthSymm3At33 = PDupwindNthSymm3(&At33[index]);
    CCTK_REAL const PDstandardNth1beta1 = PDstandardNth1(&beta1[index]);
    CCTK_REAL const PDstandardNth2beta1 = PDstandardNth2(&beta1[index]);
    CCTK_REAL const PDstandardNth3beta1 = PDstandardNth3(&beta1[index]);
    CCTK_REAL const PDstandardNth1beta2 = PDstandardNth1(&beta2[index]);
    CCTK_REAL const PDstandardNth2beta2 = PDstandardNth2(&beta2[index]);
    CCTK_REAL const PDstandardNth3beta2 = PDstandardNth3(&beta2[index]);
    CCTK_REAL const PDstandardNth1beta3 = PDstandardNth1(&beta3[index]);
    CCTK_REAL const PDstandardNth2beta3 = PDstandardNth2(&beta3[index]);
    CCTK_REAL const PDstandardNth3beta3 = PDstandardNth3(&beta3[index]);
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
    CCTK_REAL const PDstandardNth1phi = PDstandardNth1(&phi[index]);
    CCTK_REAL const PDstandardNth2phi = PDstandardNth2(&phi[index]);
    CCTK_REAL const PDstandardNth3phi = PDstandardNth3(&phi[index]);
    CCTK_REAL const PDstandardNth11phi = PDstandardNth11(&phi[index]);
    CCTK_REAL const PDstandardNth22phi = PDstandardNth22(&phi[index]);
    CCTK_REAL const PDstandardNth33phi = PDstandardNth33(&phi[index]);
    CCTK_REAL const PDstandardNth12phi = PDstandardNth12(&phi[index]);
    CCTK_REAL const PDstandardNth13phi = PDstandardNth13(&phi[index]);
    CCTK_REAL const PDstandardNth23phi = PDstandardNth23(&phi[index]);
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
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    CCTK_REAL epsdiss1 = ToReal(EpsDiss);
    
    CCTK_REAL epsdiss2 = ToReal(EpsDiss);
    
    CCTK_REAL epsdiss3 = ToReal(EpsDiss);
    
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
    
    CCTK_REAL Gtlu111 = Gtl111*gtu11 + Gtl112*gtu12 + Gtl113*gtu13;
    
    CCTK_REAL Gtlu112 = Gtl111*gtu12 + Gtl112*gtu22 + Gtl113*gtu23;
    
    CCTK_REAL Gtlu113 = Gtl111*gtu13 + Gtl112*gtu23 + Gtl113*gtu33;
    
    CCTK_REAL Gtlu121 = Gtl112*gtu11 + Gtl122*gtu12 + Gtl123*gtu13;
    
    CCTK_REAL Gtlu122 = Gtl112*gtu12 + Gtl122*gtu22 + Gtl123*gtu23;
    
    CCTK_REAL Gtlu123 = Gtl112*gtu13 + Gtl122*gtu23 + Gtl123*gtu33;
    
    CCTK_REAL Gtlu131 = Gtl113*gtu11 + Gtl123*gtu12 + Gtl133*gtu13;
    
    CCTK_REAL Gtlu132 = Gtl113*gtu12 + Gtl123*gtu22 + Gtl133*gtu23;
    
    CCTK_REAL Gtlu133 = Gtl113*gtu13 + Gtl123*gtu23 + Gtl133*gtu33;
    
    CCTK_REAL Gtlu211 = Gtl211*gtu11 + Gtl212*gtu12 + Gtl213*gtu13;
    
    CCTK_REAL Gtlu212 = Gtl211*gtu12 + Gtl212*gtu22 + Gtl213*gtu23;
    
    CCTK_REAL Gtlu213 = Gtl211*gtu13 + Gtl212*gtu23 + Gtl213*gtu33;
    
    CCTK_REAL Gtlu221 = Gtl212*gtu11 + Gtl222*gtu12 + Gtl223*gtu13;
    
    CCTK_REAL Gtlu222 = Gtl212*gtu12 + Gtl222*gtu22 + Gtl223*gtu23;
    
    CCTK_REAL Gtlu223 = Gtl212*gtu13 + Gtl222*gtu23 + Gtl223*gtu33;
    
    CCTK_REAL Gtlu231 = Gtl213*gtu11 + Gtl223*gtu12 + Gtl233*gtu13;
    
    CCTK_REAL Gtlu232 = Gtl213*gtu12 + Gtl223*gtu22 + Gtl233*gtu23;
    
    CCTK_REAL Gtlu233 = Gtl213*gtu13 + Gtl223*gtu23 + Gtl233*gtu33;
    
    CCTK_REAL Gtlu311 = Gtl311*gtu11 + Gtl312*gtu12 + Gtl313*gtu13;
    
    CCTK_REAL Gtlu312 = Gtl311*gtu12 + Gtl312*gtu22 + Gtl313*gtu23;
    
    CCTK_REAL Gtlu313 = Gtl311*gtu13 + Gtl312*gtu23 + Gtl313*gtu33;
    
    CCTK_REAL Gtlu321 = Gtl312*gtu11 + Gtl322*gtu12 + Gtl323*gtu13;
    
    CCTK_REAL Gtlu322 = Gtl312*gtu12 + Gtl322*gtu22 + Gtl323*gtu23;
    
    CCTK_REAL Gtlu323 = Gtl312*gtu13 + Gtl322*gtu23 + Gtl323*gtu33;
    
    CCTK_REAL Gtlu331 = Gtl313*gtu11 + Gtl323*gtu12 + Gtl333*gtu13;
    
    CCTK_REAL Gtlu332 = Gtl313*gtu12 + Gtl323*gtu22 + Gtl333*gtu23;
    
    CCTK_REAL Gtlu333 = Gtl313*gtu13 + Gtl323*gtu23 + Gtl333*gtu33;
    
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
    
    CCTK_REAL Rt11 = 3*(Gt111*Gtlu111 + Gt112*Gtlu112 + Gt113*Gtlu113) + 
      2*(Gt211*Gtlu121 + Gt212*Gtlu122 + Gt213*Gtlu123 + Gt311*Gtlu131 + 
      Gt312*Gtlu132 + Gt313*Gtlu133) + Gt211*Gtlu211 + Gt212*Gtlu212 + 
      Gt213*Gtlu213 + Gt311*Gtlu311 + Gt312*Gtlu312 + Gt313*Gtlu313 + 
      J11L*(gt11L*PDstandardNth1Xt1 + gt12L*PDstandardNth1Xt2 + 
      gt13L*PDstandardNth1Xt3) + J21L*(gt11L*PDstandardNth2Xt1 + 
      gt12L*PDstandardNth2Xt2 + gt13L*PDstandardNth2Xt3) + 
      J31L*(gt11L*PDstandardNth3Xt1 + gt12L*PDstandardNth3Xt2 + 
      gt13L*PDstandardNth3Xt3) + Gtl111*Xtn1 + Gtl112*Xtn2 + Gtl113*Xtn3 + 
      0.5*(-2*(gtu12*(J11L*J12L*PDstandardNth11gt11 + 
      J12L*J21L*PDstandardNth12gt11 + J11L*J22L*PDstandardNth12gt11 + 
      J12L*J31L*PDstandardNth13gt11 + J11L*J32L*PDstandardNth13gt11 + 
      dJ112L*PDstandardNth1gt11 + J21L*J22L*PDstandardNth22gt11 + 
      J22L*J31L*PDstandardNth23gt11 + J21L*J32L*PDstandardNth23gt11 + 
      dJ212L*PDstandardNth2gt11 + J31L*J32L*PDstandardNth33gt11 + 
      dJ312L*PDstandardNth3gt11) + gtu13*(J11L*J13L*PDstandardNth11gt11 + 
      J13L*J21L*PDstandardNth12gt11 + J11L*J23L*PDstandardNth12gt11 + 
      J13L*J31L*PDstandardNth13gt11 + J11L*J33L*PDstandardNth13gt11 + 
      dJ113L*PDstandardNth1gt11 + J21L*J23L*PDstandardNth22gt11 + 
      J23L*J31L*PDstandardNth23gt11 + J21L*J33L*PDstandardNth23gt11 + 
      dJ213L*PDstandardNth2gt11 + J31L*J33L*PDstandardNth33gt11 + 
      dJ313L*PDstandardNth3gt11) + gtu23*(J12L*J13L*PDstandardNth11gt11 + 
      J13L*J22L*PDstandardNth12gt11 + J12L*J23L*PDstandardNth12gt11 + 
      J13L*J32L*PDstandardNth13gt11 + J12L*J33L*PDstandardNth13gt11 + 
      dJ123L*PDstandardNth1gt11 + J22L*J23L*PDstandardNth22gt11 + 
      J23L*J32L*PDstandardNth23gt11 + J22L*J33L*PDstandardNth23gt11 + 
      dJ223L*PDstandardNth2gt11 + J32L*J33L*PDstandardNth33gt11 + 
      dJ323L*PDstandardNth3gt11)) - gtu11*(2*J11L*J21L*PDstandardNth12gt11 + 
      2*J11L*J31L*PDstandardNth13gt11 + dJ111L*PDstandardNth1gt11 + 
      2*J21L*J31L*PDstandardNth23gt11 + dJ211L*PDstandardNth2gt11 + 
      dJ311L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J11L) + 
      PDstandardNth22gt11*SQR(J21L) + PDstandardNth33gt11*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt11 + 
      2*J12L*J32L*PDstandardNth13gt11 + dJ122L*PDstandardNth1gt11 + 
      2*J22L*J32L*PDstandardNth23gt11 + dJ222L*PDstandardNth2gt11 + 
      dJ322L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J12L) + 
      PDstandardNth22gt11*SQR(J22L) + PDstandardNth33gt11*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt11 + 
      2*J13L*J33L*PDstandardNth13gt11 + dJ133L*PDstandardNth1gt11 + 
      2*J23L*J33L*PDstandardNth23gt11 + dJ233L*PDstandardNth2gt11 + 
      dJ333L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J13L) + 
      PDstandardNth22gt11*SQR(J23L) + PDstandardNth33gt11*SQR(J33L)));
    
    CCTK_REAL Rt12 = Gt122*Gtlu112 + Gt123*Gtlu113 + (Gt111 + 
      Gt212)*Gtlu121 + Gt222*Gtlu122 + (Gt113 + Gt223)*Gtlu123 + 
      Gt322*Gtlu132 + Gt323*Gtlu133 + Gt111*Gtlu211 + Gt112*(Gtlu111 + 
      Gtlu122 + Gtlu212) + Gt113*Gtlu213 + 2*(Gt211*Gtlu221 + Gt212*Gtlu222 + 
      Gt213*Gtlu223) + Gt311*(Gtlu231 + Gtlu321) + Gt312*(Gtlu131 + Gtlu232 + 
      Gtlu322) + Gt313*(Gtlu233 + Gtlu323) + 0.5*((gt12L*J11L + 
      gt11L*J12L)*PDstandardNth1Xt1 + (gt22L*J11L + 
      gt12L*J12L)*PDstandardNth1Xt2 + (gt23L*J11L + 
      gt13L*J12L)*PDstandardNth1Xt3 + (gt12L*J21L + 
      gt11L*J22L)*PDstandardNth2Xt1 + (gt22L*J21L + 
      gt12L*J22L)*PDstandardNth2Xt2 + (gt23L*J21L + 
      gt13L*J22L)*PDstandardNth2Xt3 - 2*(gtu12*(J11L*J12L*PDstandardNth11gt12 
      + J12L*J21L*PDstandardNth12gt12 + J11L*J22L*PDstandardNth12gt12 + 
      J12L*J31L*PDstandardNth13gt12 + J11L*J32L*PDstandardNth13gt12 + 
      dJ112L*PDstandardNth1gt12 + J21L*J22L*PDstandardNth22gt12 + 
      J22L*J31L*PDstandardNth23gt12 + J21L*J32L*PDstandardNth23gt12 + 
      dJ212L*PDstandardNth2gt12 + J31L*J32L*PDstandardNth33gt12 + 
      dJ312L*PDstandardNth3gt12) + gtu13*(J11L*J13L*PDstandardNth11gt12 + 
      J13L*J21L*PDstandardNth12gt12 + J11L*J23L*PDstandardNth12gt12 + 
      J13L*J31L*PDstandardNth13gt12 + J11L*J33L*PDstandardNth13gt12 + 
      dJ113L*PDstandardNth1gt12 + J21L*J23L*PDstandardNth22gt12 + 
      J23L*J31L*PDstandardNth23gt12 + J21L*J33L*PDstandardNth23gt12 + 
      dJ213L*PDstandardNth2gt12 + J31L*J33L*PDstandardNth33gt12 + 
      dJ313L*PDstandardNth3gt12) + gtu23*(J12L*J13L*PDstandardNth11gt12 + 
      J13L*J22L*PDstandardNth12gt12 + J12L*J23L*PDstandardNth12gt12 + 
      J13L*J32L*PDstandardNth13gt12 + J12L*J33L*PDstandardNth13gt12 + 
      dJ123L*PDstandardNth1gt12 + J22L*J23L*PDstandardNth22gt12 + 
      J23L*J32L*PDstandardNth23gt12 + J22L*J33L*PDstandardNth23gt12 + 
      dJ223L*PDstandardNth2gt12 + J32L*J33L*PDstandardNth33gt12 + 
      dJ323L*PDstandardNth3gt12)) + (gt12L*J31L + 
      gt11L*J32L)*PDstandardNth3Xt1 + (gt22L*J31L + 
      gt12L*J32L)*PDstandardNth3Xt2 + (gt23L*J31L + 
      gt13L*J32L)*PDstandardNth3Xt3 + (Gtl112 + Gtl211)*Xtn1 + (Gtl122 + 
      Gtl212)*Xtn2 + (Gtl123 + Gtl213)*Xtn3 - 
      gtu11*(2*J11L*J21L*PDstandardNth12gt12 + 
      2*J11L*J31L*PDstandardNth13gt12 + dJ111L*PDstandardNth1gt12 + 
      2*J21L*J31L*PDstandardNth23gt12 + dJ211L*PDstandardNth2gt12 + 
      dJ311L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J11L) + 
      PDstandardNth22gt12*SQR(J21L) + PDstandardNth33gt12*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt12 + 
      2*J12L*J32L*PDstandardNth13gt12 + dJ122L*PDstandardNth1gt12 + 
      2*J22L*J32L*PDstandardNth23gt12 + dJ222L*PDstandardNth2gt12 + 
      dJ322L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J12L) + 
      PDstandardNth22gt12*SQR(J22L) + PDstandardNth33gt12*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt12 + 
      2*J13L*J33L*PDstandardNth13gt12 + dJ133L*PDstandardNth1gt12 + 
      2*J23L*J33L*PDstandardNth23gt12 + dJ233L*PDstandardNth2gt12 + 
      dJ333L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J13L) + 
      PDstandardNth22gt12*SQR(J23L) + PDstandardNth33gt12*SQR(J33L)));
    
    CCTK_REAL Rt13 = Gt123*Gtlu112 + Gt133*Gtlu113 + Gt223*Gtlu122 + 
      Gt233*Gtlu123 + (Gt111 + Gt313)*Gtlu131 + (Gt112 + Gt323)*Gtlu132 + 
      Gt333*Gtlu133 + Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*(Gtlu111 + 
      Gtlu133 + Gtlu313) + Gt211*(Gtlu231 + Gtlu321) + Gt212*(Gtlu232 + 
      Gtlu322) + Gt213*(Gtlu121 + Gtlu233 + Gtlu323) + 2*(Gt311*Gtlu331 + 
      Gt312*Gtlu332 + Gt313*Gtlu333) + 0.5*((gt13L*J11L + 
      gt11L*J13L)*PDstandardNth1Xt1 + (gt23L*J11L + 
      gt12L*J13L)*PDstandardNth1Xt2 + (gt33L*J11L + 
      gt13L*J13L)*PDstandardNth1Xt3 + (gt13L*J21L + 
      gt11L*J23L)*PDstandardNth2Xt1 + (gt23L*J21L + 
      gt12L*J23L)*PDstandardNth2Xt2 + (gt33L*J21L + 
      gt13L*J23L)*PDstandardNth2Xt3 - 2*(gtu12*(J11L*J12L*PDstandardNth11gt13 
      + J12L*J21L*PDstandardNth12gt13 + J11L*J22L*PDstandardNth12gt13 + 
      J12L*J31L*PDstandardNth13gt13 + J11L*J32L*PDstandardNth13gt13 + 
      dJ112L*PDstandardNth1gt13 + J21L*J22L*PDstandardNth22gt13 + 
      J22L*J31L*PDstandardNth23gt13 + J21L*J32L*PDstandardNth23gt13 + 
      dJ212L*PDstandardNth2gt13 + J31L*J32L*PDstandardNth33gt13 + 
      dJ312L*PDstandardNth3gt13) + gtu13*(J11L*J13L*PDstandardNth11gt13 + 
      J13L*J21L*PDstandardNth12gt13 + J11L*J23L*PDstandardNth12gt13 + 
      J13L*J31L*PDstandardNth13gt13 + J11L*J33L*PDstandardNth13gt13 + 
      dJ113L*PDstandardNth1gt13 + J21L*J23L*PDstandardNth22gt13 + 
      J23L*J31L*PDstandardNth23gt13 + J21L*J33L*PDstandardNth23gt13 + 
      dJ213L*PDstandardNth2gt13 + J31L*J33L*PDstandardNth33gt13 + 
      dJ313L*PDstandardNth3gt13) + gtu23*(J12L*J13L*PDstandardNth11gt13 + 
      J13L*J22L*PDstandardNth12gt13 + J12L*J23L*PDstandardNth12gt13 + 
      J13L*J32L*PDstandardNth13gt13 + J12L*J33L*PDstandardNth13gt13 + 
      dJ123L*PDstandardNth1gt13 + J22L*J23L*PDstandardNth22gt13 + 
      J23L*J32L*PDstandardNth23gt13 + J22L*J33L*PDstandardNth23gt13 + 
      dJ223L*PDstandardNth2gt13 + J32L*J33L*PDstandardNth33gt13 + 
      dJ323L*PDstandardNth3gt13)) + (gt13L*J31L + 
      gt11L*J33L)*PDstandardNth3Xt1 + (gt23L*J31L + 
      gt12L*J33L)*PDstandardNth3Xt2 + (gt33L*J31L + 
      gt13L*J33L)*PDstandardNth3Xt3 + (Gtl113 + Gtl311)*Xtn1 + (Gtl123 + 
      Gtl312)*Xtn2 + (Gtl133 + Gtl313)*Xtn3 - 
      gtu11*(2*J11L*J21L*PDstandardNth12gt13 + 
      2*J11L*J31L*PDstandardNth13gt13 + dJ111L*PDstandardNth1gt13 + 
      2*J21L*J31L*PDstandardNth23gt13 + dJ211L*PDstandardNth2gt13 + 
      dJ311L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J11L) + 
      PDstandardNth22gt13*SQR(J21L) + PDstandardNth33gt13*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt13 + 
      2*J12L*J32L*PDstandardNth13gt13 + dJ122L*PDstandardNth1gt13 + 
      2*J22L*J32L*PDstandardNth23gt13 + dJ222L*PDstandardNth2gt13 + 
      dJ322L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J12L) + 
      PDstandardNth22gt13*SQR(J22L) + PDstandardNth33gt13*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt13 + 
      2*J13L*J33L*PDstandardNth13gt13 + dJ133L*PDstandardNth1gt13 + 
      2*J23L*J33L*PDstandardNth23gt13 + dJ233L*PDstandardNth2gt13 + 
      dJ333L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J13L) + 
      PDstandardNth22gt13*SQR(J23L) + PDstandardNth33gt13*SQR(J33L)));
    
    CCTK_REAL Rt22 = Gt112*(Gtlu121 + 2*Gtlu211) + Gt122*(Gtlu122 + 
      2*Gtlu212) + Gt123*(Gtlu123 + 2*Gtlu213) + 3*(Gt212*Gtlu221 + 
      Gt222*Gtlu222 + Gt223*Gtlu223) + 2*(Gt312*Gtlu231 + Gt322*Gtlu232 + 
      Gt323*Gtlu233) + Gt312*Gtlu321 + Gt322*Gtlu322 + Gt323*Gtlu323 + 
      J12L*(gt12L*PDstandardNth1Xt1 + gt22L*PDstandardNth1Xt2 + 
      gt23L*PDstandardNth1Xt3) + J22L*(gt12L*PDstandardNth2Xt1 + 
      gt22L*PDstandardNth2Xt2 + gt23L*PDstandardNth2Xt3) + 
      J32L*(gt12L*PDstandardNth3Xt1 + gt22L*PDstandardNth3Xt2 + 
      gt23L*PDstandardNth3Xt3) + Gtl212*Xtn1 + Gtl222*Xtn2 + Gtl223*Xtn3 + 
      0.5*(-2*(gtu12*(J11L*J12L*PDstandardNth11gt22 + 
      J12L*J21L*PDstandardNth12gt22 + J11L*J22L*PDstandardNth12gt22 + 
      J12L*J31L*PDstandardNth13gt22 + J11L*J32L*PDstandardNth13gt22 + 
      dJ112L*PDstandardNth1gt22 + J21L*J22L*PDstandardNth22gt22 + 
      J22L*J31L*PDstandardNth23gt22 + J21L*J32L*PDstandardNth23gt22 + 
      dJ212L*PDstandardNth2gt22 + J31L*J32L*PDstandardNth33gt22 + 
      dJ312L*PDstandardNth3gt22) + gtu13*(J11L*J13L*PDstandardNth11gt22 + 
      J13L*J21L*PDstandardNth12gt22 + J11L*J23L*PDstandardNth12gt22 + 
      J13L*J31L*PDstandardNth13gt22 + J11L*J33L*PDstandardNth13gt22 + 
      dJ113L*PDstandardNth1gt22 + J21L*J23L*PDstandardNth22gt22 + 
      J23L*J31L*PDstandardNth23gt22 + J21L*J33L*PDstandardNth23gt22 + 
      dJ213L*PDstandardNth2gt22 + J31L*J33L*PDstandardNth33gt22 + 
      dJ313L*PDstandardNth3gt22) + gtu23*(J12L*J13L*PDstandardNth11gt22 + 
      J13L*J22L*PDstandardNth12gt22 + J12L*J23L*PDstandardNth12gt22 + 
      J13L*J32L*PDstandardNth13gt22 + J12L*J33L*PDstandardNth13gt22 + 
      dJ123L*PDstandardNth1gt22 + J22L*J23L*PDstandardNth22gt22 + 
      J23L*J32L*PDstandardNth23gt22 + J22L*J33L*PDstandardNth23gt22 + 
      dJ223L*PDstandardNth2gt22 + J32L*J33L*PDstandardNth33gt22 + 
      dJ323L*PDstandardNth3gt22)) - gtu11*(2*J11L*J21L*PDstandardNth12gt22 + 
      2*J11L*J31L*PDstandardNth13gt22 + dJ111L*PDstandardNth1gt22 + 
      2*J21L*J31L*PDstandardNth23gt22 + dJ211L*PDstandardNth2gt22 + 
      dJ311L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J11L) + 
      PDstandardNth22gt22*SQR(J21L) + PDstandardNth33gt22*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt22 + 
      2*J12L*J32L*PDstandardNth13gt22 + dJ122L*PDstandardNth1gt22 + 
      2*J22L*J32L*PDstandardNth23gt22 + dJ222L*PDstandardNth2gt22 + 
      dJ322L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J12L) + 
      PDstandardNth22gt22*SQR(J22L) + PDstandardNth33gt22*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt22 + 
      2*J13L*J33L*PDstandardNth13gt22 + dJ133L*PDstandardNth1gt22 + 
      2*J23L*J33L*PDstandardNth23gt22 + dJ233L*PDstandardNth2gt22 + 
      dJ333L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J13L) + 
      PDstandardNth22gt22*SQR(J23L) + PDstandardNth33gt22*SQR(J33L)));
    
    CCTK_REAL Rt23 = Gt113*Gtlu211 + Gt133*Gtlu213 + Gt213*Gtlu221 + 
      Gt233*Gtlu223 + (Gt212 + Gt313)*Gtlu231 + (Gt222 + Gt323)*Gtlu232 + 
      Gt333*Gtlu233 + Gt112*(Gtlu131 + Gtlu311) + Gt122*(Gtlu132 + Gtlu312) + 
      Gt123*(Gtlu133 + Gtlu212 + Gtlu313) + Gt212*Gtlu321 + Gt222*Gtlu322 + 
      Gt223*(Gtlu222 + Gtlu233 + Gtlu323) + 2*(Gt312*Gtlu331 + Gt322*Gtlu332 
      + Gt323*Gtlu333) + 0.5*((gt13L*J12L + gt12L*J13L)*PDstandardNth1Xt1 + 
      (gt23L*J12L + gt22L*J13L)*PDstandardNth1Xt2 + (gt33L*J12L + 
      gt23L*J13L)*PDstandardNth1Xt3 + (gt13L*J22L + 
      gt12L*J23L)*PDstandardNth2Xt1 + (gt23L*J22L + 
      gt22L*J23L)*PDstandardNth2Xt2 + (gt33L*J22L + 
      gt23L*J23L)*PDstandardNth2Xt3 - 2*(gtu12*(J11L*J12L*PDstandardNth11gt23 
      + J12L*J21L*PDstandardNth12gt23 + J11L*J22L*PDstandardNth12gt23 + 
      J12L*J31L*PDstandardNth13gt23 + J11L*J32L*PDstandardNth13gt23 + 
      dJ112L*PDstandardNth1gt23 + J21L*J22L*PDstandardNth22gt23 + 
      J22L*J31L*PDstandardNth23gt23 + J21L*J32L*PDstandardNth23gt23 + 
      dJ212L*PDstandardNth2gt23 + J31L*J32L*PDstandardNth33gt23 + 
      dJ312L*PDstandardNth3gt23) + gtu13*(J11L*J13L*PDstandardNth11gt23 + 
      J13L*J21L*PDstandardNth12gt23 + J11L*J23L*PDstandardNth12gt23 + 
      J13L*J31L*PDstandardNth13gt23 + J11L*J33L*PDstandardNth13gt23 + 
      dJ113L*PDstandardNth1gt23 + J21L*J23L*PDstandardNth22gt23 + 
      J23L*J31L*PDstandardNth23gt23 + J21L*J33L*PDstandardNth23gt23 + 
      dJ213L*PDstandardNth2gt23 + J31L*J33L*PDstandardNth33gt23 + 
      dJ313L*PDstandardNth3gt23) + gtu23*(J12L*J13L*PDstandardNth11gt23 + 
      J13L*J22L*PDstandardNth12gt23 + J12L*J23L*PDstandardNth12gt23 + 
      J13L*J32L*PDstandardNth13gt23 + J12L*J33L*PDstandardNth13gt23 + 
      dJ123L*PDstandardNth1gt23 + J22L*J23L*PDstandardNth22gt23 + 
      J23L*J32L*PDstandardNth23gt23 + J22L*J33L*PDstandardNth23gt23 + 
      dJ223L*PDstandardNth2gt23 + J32L*J33L*PDstandardNth33gt23 + 
      dJ323L*PDstandardNth3gt23)) + (gt13L*J32L + 
      gt12L*J33L)*PDstandardNth3Xt1 + (gt23L*J32L + 
      gt22L*J33L)*PDstandardNth3Xt2 + (gt33L*J32L + 
      gt23L*J33L)*PDstandardNth3Xt3 + (Gtl213 + Gtl312)*Xtn1 + (Gtl223 + 
      Gtl322)*Xtn2 + (Gtl233 + Gtl323)*Xtn3 - 
      gtu11*(2*J11L*J21L*PDstandardNth12gt23 + 
      2*J11L*J31L*PDstandardNth13gt23 + dJ111L*PDstandardNth1gt23 + 
      2*J21L*J31L*PDstandardNth23gt23 + dJ211L*PDstandardNth2gt23 + 
      dJ311L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J11L) + 
      PDstandardNth22gt23*SQR(J21L) + PDstandardNth33gt23*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt23 + 
      2*J12L*J32L*PDstandardNth13gt23 + dJ122L*PDstandardNth1gt23 + 
      2*J22L*J32L*PDstandardNth23gt23 + dJ222L*PDstandardNth2gt23 + 
      dJ322L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J12L) + 
      PDstandardNth22gt23*SQR(J22L) + PDstandardNth33gt23*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt23 + 
      2*J13L*J33L*PDstandardNth13gt23 + dJ133L*PDstandardNth1gt23 + 
      2*J23L*J33L*PDstandardNth23gt23 + dJ233L*PDstandardNth2gt23 + 
      dJ333L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J13L) + 
      PDstandardNth22gt23*SQR(J23L) + PDstandardNth33gt23*SQR(J33L)));
    
    CCTK_REAL Rt33 = Gt113*(Gtlu131 + 2*Gtlu311) + Gt123*(Gtlu132 + 
      2*Gtlu312) + Gt133*(Gtlu133 + 2*Gtlu313) + Gt213*(Gtlu231 + 2*Gtlu321) 
      + Gt223*(Gtlu232 + 2*Gtlu322) + Gt233*(Gtlu233 + 2*Gtlu323) + 
      3*(Gt313*Gtlu331 + Gt323*Gtlu332 + Gt333*Gtlu333) + 
      J13L*(gt13L*PDstandardNth1Xt1 + gt23L*PDstandardNth1Xt2 + 
      gt33L*PDstandardNth1Xt3) + J23L*(gt13L*PDstandardNth2Xt1 + 
      gt23L*PDstandardNth2Xt2 + gt33L*PDstandardNth2Xt3) + 
      J33L*(gt13L*PDstandardNth3Xt1 + gt23L*PDstandardNth3Xt2 + 
      gt33L*PDstandardNth3Xt3) + Gtl313*Xtn1 + Gtl323*Xtn2 + Gtl333*Xtn3 + 
      0.5*(-2*(gtu12*(J11L*J12L*PDstandardNth11gt33 + 
      J12L*J21L*PDstandardNth12gt33 + J11L*J22L*PDstandardNth12gt33 + 
      J12L*J31L*PDstandardNth13gt33 + J11L*J32L*PDstandardNth13gt33 + 
      dJ112L*PDstandardNth1gt33 + J21L*J22L*PDstandardNth22gt33 + 
      J22L*J31L*PDstandardNth23gt33 + J21L*J32L*PDstandardNth23gt33 + 
      dJ212L*PDstandardNth2gt33 + J31L*J32L*PDstandardNth33gt33 + 
      dJ312L*PDstandardNth3gt33) + gtu13*(J11L*J13L*PDstandardNth11gt33 + 
      J13L*J21L*PDstandardNth12gt33 + J11L*J23L*PDstandardNth12gt33 + 
      J13L*J31L*PDstandardNth13gt33 + J11L*J33L*PDstandardNth13gt33 + 
      dJ113L*PDstandardNth1gt33 + J21L*J23L*PDstandardNth22gt33 + 
      J23L*J31L*PDstandardNth23gt33 + J21L*J33L*PDstandardNth23gt33 + 
      dJ213L*PDstandardNth2gt33 + J31L*J33L*PDstandardNth33gt33 + 
      dJ313L*PDstandardNth3gt33) + gtu23*(J12L*J13L*PDstandardNth11gt33 + 
      J13L*J22L*PDstandardNth12gt33 + J12L*J23L*PDstandardNth12gt33 + 
      J13L*J32L*PDstandardNth13gt33 + J12L*J33L*PDstandardNth13gt33 + 
      dJ123L*PDstandardNth1gt33 + J22L*J23L*PDstandardNth22gt33 + 
      J23L*J32L*PDstandardNth23gt33 + J22L*J33L*PDstandardNth23gt33 + 
      dJ223L*PDstandardNth2gt33 + J32L*J33L*PDstandardNth33gt33 + 
      dJ323L*PDstandardNth3gt33)) - gtu11*(2*J11L*J21L*PDstandardNth12gt33 + 
      2*J11L*J31L*PDstandardNth13gt33 + dJ111L*PDstandardNth1gt33 + 
      2*J21L*J31L*PDstandardNth23gt33 + dJ211L*PDstandardNth2gt33 + 
      dJ311L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J11L) + 
      PDstandardNth22gt33*SQR(J21L) + PDstandardNth33gt33*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt33 + 
      2*J12L*J32L*PDstandardNth13gt33 + dJ122L*PDstandardNth1gt33 + 
      2*J22L*J32L*PDstandardNth23gt33 + dJ222L*PDstandardNth2gt33 + 
      dJ322L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J12L) + 
      PDstandardNth22gt33*SQR(J22L) + PDstandardNth33gt33*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt33 + 
      2*J13L*J33L*PDstandardNth13gt33 + dJ133L*PDstandardNth1gt33 + 
      2*J23L*J33L*PDstandardNth23gt33 + dJ233L*PDstandardNth2gt33 + 
      dJ333L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J13L) + 
      PDstandardNth22gt33*SQR(J23L) + PDstandardNth33gt33*SQR(J33L)));
    
    CCTK_REAL fac1 = IfThen(ToReal(conformalMethod),-0.5*INV(phiL),1);
    
    CCTK_REAL cdphi1 = fac1*(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + J31L*PDstandardNth3phi);
    
    CCTK_REAL cdphi2 = fac1*(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + J32L*PDstandardNth3phi);
    
    CCTK_REAL cdphi3 = fac1*(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi);
    
    CCTK_REAL fac2 = IfThen(ToReal(conformalMethod),0.5*INV(SQR(phiL)),0);
    
    CCTK_REAL cdphi211 = fac1*((dJ111L - Gt111*J11L - Gt211*J12L - 
      Gt311*J13L)*PDstandardNth1phi + 2*(J11L*(J21L*PDstandardNth12phi + 
      J31L*PDstandardNth13phi) + J21L*J31L*PDstandardNth23phi) + (dJ211L - 
      Gt111*J21L - Gt211*J22L - Gt311*J23L)*PDstandardNth2phi + (dJ311L - 
      Gt111*J31L - Gt211*J32L - Gt311*J33L)*PDstandardNth3phi + 
      PDstandardNth11phi*SQR(J11L) + PDstandardNth22phi*SQR(J21L) + 
      PDstandardNth33phi*SQR(J31L)) + fac2*SQR(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + J31L*PDstandardNth3phi);
    
    CCTK_REAL cdphi212 = fac2*(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + 
      J31L*PDstandardNth3phi)*(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + J32L*PDstandardNth3phi) + 
      fac1*(J12L*(J11L*PDstandardNth11phi + J21L*PDstandardNth12phi + 
      J31L*PDstandardNth13phi) + J11L*(J22L*PDstandardNth12phi + 
      J32L*PDstandardNth13phi) + (dJ112L - Gt112*J11L - Gt212*J12L - 
      Gt312*J13L)*PDstandardNth1phi + J22L*(J21L*PDstandardNth22phi + 
      J31L*PDstandardNth23phi) + (dJ212L - Gt112*J21L - Gt212*J22L - 
      Gt312*J23L)*PDstandardNth2phi + J32L*(J21L*PDstandardNth23phi + 
      J31L*PDstandardNth33phi) + (dJ312L - Gt112*J31L - Gt212*J32L - 
      Gt312*J33L)*PDstandardNth3phi);
    
    CCTK_REAL cdphi213 = fac2*(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + 
      J31L*PDstandardNth3phi)*(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi) + 
      fac1*(J13L*(J11L*PDstandardNth11phi + J21L*PDstandardNth12phi + 
      J31L*PDstandardNth13phi) + J11L*(J23L*PDstandardNth12phi + 
      J33L*PDstandardNth13phi) + (dJ113L - Gt113*J11L - Gt213*J12L - 
      Gt313*J13L)*PDstandardNth1phi + J23L*(J21L*PDstandardNth22phi + 
      J31L*PDstandardNth23phi) + (dJ213L - Gt113*J21L - Gt213*J22L - 
      Gt313*J23L)*PDstandardNth2phi + J33L*(J21L*PDstandardNth23phi + 
      J31L*PDstandardNth33phi) + (dJ313L - Gt113*J31L - Gt213*J32L - 
      Gt313*J33L)*PDstandardNth3phi);
    
    CCTK_REAL cdphi222 = fac1*((dJ122L - Gt122*J11L - Gt222*J12L - 
      Gt322*J13L)*PDstandardNth1phi + 2*(J12L*(J22L*PDstandardNth12phi + 
      J32L*PDstandardNth13phi) + J22L*J32L*PDstandardNth23phi) + (dJ222L - 
      Gt122*J21L - Gt222*J22L - Gt322*J23L)*PDstandardNth2phi + (dJ322L - 
      Gt122*J31L - Gt222*J32L - Gt322*J33L)*PDstandardNth3phi + 
      PDstandardNth11phi*SQR(J12L) + PDstandardNth22phi*SQR(J22L) + 
      PDstandardNth33phi*SQR(J32L)) + fac2*SQR(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + J32L*PDstandardNth3phi);
    
    CCTK_REAL cdphi223 = fac2*(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + 
      J32L*PDstandardNth3phi)*(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi) + 
      fac1*(J13L*(J12L*PDstandardNth11phi + J22L*PDstandardNth12phi + 
      J32L*PDstandardNth13phi) + J12L*(J23L*PDstandardNth12phi + 
      J33L*PDstandardNth13phi) + (dJ123L - Gt123*J11L - Gt223*J12L - 
      Gt323*J13L)*PDstandardNth1phi + J23L*(J22L*PDstandardNth22phi + 
      J32L*PDstandardNth23phi) + (dJ223L - Gt123*J21L - Gt223*J22L - 
      Gt323*J23L)*PDstandardNth2phi + J33L*(J22L*PDstandardNth23phi + 
      J32L*PDstandardNth33phi) + (dJ323L - Gt123*J31L - Gt223*J32L - 
      Gt323*J33L)*PDstandardNth3phi);
    
    CCTK_REAL cdphi233 = fac1*((dJ133L - Gt133*J11L - Gt233*J12L - 
      Gt333*J13L)*PDstandardNth1phi + 2*(J13L*(J23L*PDstandardNth12phi + 
      J33L*PDstandardNth13phi) + J23L*J33L*PDstandardNth23phi) + (dJ233L - 
      Gt133*J21L - Gt233*J22L - Gt333*J23L)*PDstandardNth2phi + (dJ333L - 
      Gt133*J31L - Gt233*J32L - Gt333*J33L)*PDstandardNth3phi + 
      PDstandardNth11phi*SQR(J13L) + PDstandardNth22phi*SQR(J23L) + 
      PDstandardNth33phi*SQR(J33L)) + fac2*SQR(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi);
    
    CCTK_REAL Rphi11 = -2*(cdphi211 + 2*(-1 + gt11L*gtu11)*SQR(cdphi1) + 
      gt11L*(cdphi211*gtu11 + 4*(cdphi1*(cdphi2*gtu12 + cdphi3*gtu13) + 
      cdphi2*cdphi3*gtu23) + cdphi233*gtu33 + gtu22*(cdphi222 + 
      2*SQR(cdphi2)) + 2*(cdphi212*gtu12 + cdphi213*gtu13 + cdphi223*gtu23 + 
      gtu33*SQR(cdphi3))));
    
    CCTK_REAL Rphi12 = -2*(cdphi212 + cdphi1*(cdphi2*(-2 + 4*gt12L*gtu12) 
      + 4*cdphi3*gt12L*gtu13) + gt12L*(cdphi211*gtu11 + 4*cdphi2*cdphi3*gtu23 
      + 2*(cdphi212*gtu12 + cdphi213*gtu13 + cdphi223*gtu23 + 
      gtu11*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2)) + gtu33*(cdphi233 
      + 2*SQR(cdphi3))));
    
    CCTK_REAL Rphi13 = -2*(cdphi213 + cdphi1*(4*cdphi2*gt13L*gtu12 + 
      cdphi3*(-2 + 4*gt13L*gtu13)) + gt13L*(cdphi211*gtu11 + 
      4*cdphi2*cdphi3*gtu23 + 2*(cdphi212*gtu12 + cdphi213*gtu13 + 
      cdphi223*gtu23 + gtu11*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2)) 
      + gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL Rphi22 = -2*(cdphi222 + 2*(-1 + gt22L*gtu22)*SQR(cdphi2) + 
      gt22L*(cdphi222*gtu22 + 4*(cdphi1*cdphi3*gtu13 + cdphi2*(cdphi1*gtu12 + 
      cdphi3*gtu23)) + cdphi233*gtu33 + gtu11*(cdphi211 + 2*SQR(cdphi1)) + 
      2*(cdphi212*gtu12 + cdphi213*gtu13 + cdphi223*gtu23 + 
      gtu33*SQR(cdphi3))));
    
    CCTK_REAL Rphi23 = -2*(cdphi223 + cdphi2*(4*cdphi1*gt23L*gtu12 + 
      cdphi3*(-2 + 4*gt23L*gtu23)) + gt23L*(4*cdphi1*cdphi3*gtu13 + 
      cdphi222*gtu22 + gtu11*(cdphi211 + 2*SQR(cdphi1)) + 2*(cdphi212*gtu12 + 
      cdphi213*gtu13 + cdphi223*gtu23 + gtu22*SQR(cdphi2)) + gtu33*(cdphi233 
      + 2*SQR(cdphi3))));
    
    CCTK_REAL Rphi33 = -2*(cdphi233 + gt33L*((4*cdphi1*cdphi2 + 
      2*cdphi212)*gtu12 + 4*cdphi3*(cdphi1*gtu13 + cdphi2*gtu23) + 
      2*(cdphi213*gtu13 + cdphi223*gtu23) + cdphi233*gtu33 + gtu11*(cdphi211 
      + 2*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2))) + 2*(-1 + 
      gt33L*gtu33)*SQR(cdphi3));
    
    CCTK_REAL Atm11 = At11L*gtu11 + At12L*gtu12 + At13L*gtu13;
    
    CCTK_REAL Atm21 = At11L*gtu12 + At12L*gtu22 + At13L*gtu23;
    
    CCTK_REAL Atm31 = At11L*gtu13 + At12L*gtu23 + At13L*gtu33;
    
    CCTK_REAL Atm12 = At12L*gtu11 + At22L*gtu12 + At23L*gtu13;
    
    CCTK_REAL Atm22 = At12L*gtu12 + At22L*gtu22 + At23L*gtu23;
    
    CCTK_REAL Atm32 = At12L*gtu13 + At22L*gtu23 + At23L*gtu33;
    
    CCTK_REAL Atm13 = At13L*gtu11 + At23L*gtu12 + At33L*gtu13;
    
    CCTK_REAL Atm23 = At13L*gtu12 + At23L*gtu22 + At33L*gtu23;
    
    CCTK_REAL Atm33 = At13L*gtu13 + At23L*gtu23 + At33L*gtu33;
    
    CCTK_REAL e4phi = 
      IfThen(ToReal(conformalMethod),INV(SQR(phiL)),exp(4*phiL));
    
    CCTK_REAL em4phi = INV(e4phi);
    
    CCTK_REAL g11 = e4phi*gt11L;
    
    CCTK_REAL g12 = e4phi*gt12L;
    
    CCTK_REAL g13 = e4phi*gt13L;
    
    CCTK_REAL g22 = e4phi*gt22L;
    
    CCTK_REAL g23 = e4phi*gt23L;
    
    CCTK_REAL g33 = e4phi*gt33L;
    
    CCTK_REAL gu11 = em4phi*gtu11;
    
    CCTK_REAL gu12 = em4phi*gtu12;
    
    CCTK_REAL gu13 = em4phi*gtu13;
    
    CCTK_REAL gu22 = em4phi*gtu22;
    
    CCTK_REAL gu23 = em4phi*gtu23;
    
    CCTK_REAL gu33 = em4phi*gtu33;
    
    CCTK_REAL R11 = Rphi11 + Rt11;
    
    CCTK_REAL R12 = Rphi12 + Rt12;
    
    CCTK_REAL R13 = Rphi13 + Rt13;
    
    CCTK_REAL R22 = Rphi22 + Rt22;
    
    CCTK_REAL R23 = Rphi23 + Rt23;
    
    CCTK_REAL R33 = Rphi33 + Rt33;
    
    CCTK_REAL trS = em4phi*(eTxxL*gtu11 + eTyyL*gtu22 + 2*(eTxyL*gtu12 + 
      eTxzL*gtu13 + eTyzL*gtu23) + eTzzL*gtu33);
    
    CCTK_REAL Ats11 = (-dJ111L + 4*cdphi1*J11L + Gt111*J11L + Gt211*J12L + 
      Gt311*J13L)*PDstandardNth1alpha - 2*(J11L*J21L*PDstandardNth12alpha + 
      J11L*J31L*PDstandardNth13alpha + J21L*J31L*PDstandardNth23alpha) + 
      (-dJ211L + 4*cdphi1*J21L + Gt111*J21L + Gt211*J22L + 
      Gt311*J23L)*PDstandardNth2alpha + (-dJ311L + 4*cdphi1*J31L + Gt111*J31L 
      + Gt211*J32L + Gt311*J33L)*PDstandardNth3alpha + alphaL*R11 - 
      PDstandardNth11alpha*SQR(J11L) - PDstandardNth22alpha*SQR(J21L) - 
      PDstandardNth33alpha*SQR(J31L);
    
    CCTK_REAL Ats12 = J12L*(-(J11L*PDstandardNth11alpha) - 
      J21L*PDstandardNth12alpha - J31L*PDstandardNth13alpha) + 
      J11L*(-(J22L*PDstandardNth12alpha) - J32L*PDstandardNth13alpha) + 
      (-dJ112L + 2*cdphi2*J11L + Gt112*J11L + 2*cdphi1*J12L + Gt212*J12L + 
      Gt312*J13L)*PDstandardNth1alpha + J22L*(-(J21L*PDstandardNth22alpha) - 
      J31L*PDstandardNth23alpha) + (-dJ212L + 2*cdphi2*J21L + Gt112*J21L + 
      2*cdphi1*J22L + Gt212*J22L + Gt312*J23L)*PDstandardNth2alpha + 
      J32L*(-(J21L*PDstandardNth23alpha) - J31L*PDstandardNth33alpha) + 
      Gt112*J31L*PDstandardNth3alpha + (-dJ312L + 
      2*cdphi2*J31L)*PDstandardNth3alpha + 2*cdphi1*J32L*PDstandardNth3alpha 
      + Gt212*J32L*PDstandardNth3alpha + Gt312*J33L*PDstandardNth3alpha + 
      alphaL*R12;
    
    CCTK_REAL Ats13 = J13L*(-(J11L*PDstandardNth11alpha) - 
      J21L*PDstandardNth12alpha - J31L*PDstandardNth13alpha) + 
      J11L*(-(J23L*PDstandardNth12alpha) - J33L*PDstandardNth13alpha) + 
      (-dJ113L + 2*cdphi3*J11L + Gt113*J11L + Gt213*J12L + 2*cdphi1*J13L + 
      Gt313*J13L)*PDstandardNth1alpha + J23L*(-(J21L*PDstandardNth22alpha) - 
      J31L*PDstandardNth23alpha) + (-dJ213L + 2*cdphi3*J21L + Gt113*J21L + 
      Gt213*J22L + 2*cdphi1*J23L + Gt313*J23L)*PDstandardNth2alpha + 
      J33L*(-(J21L*PDstandardNth23alpha) - J31L*PDstandardNth33alpha) + 
      Gt113*J31L*PDstandardNth3alpha + (-dJ313L + 
      2*cdphi3*J31L)*PDstandardNth3alpha + Gt213*J32L*PDstandardNth3alpha + 
      2*cdphi1*J33L*PDstandardNth3alpha + Gt313*J33L*PDstandardNth3alpha + 
      alphaL*R13;
    
    CCTK_REAL Ats22 = (-dJ122L + Gt122*J11L + 4*cdphi2*J12L + Gt222*J12L + 
      Gt322*J13L)*PDstandardNth1alpha - 2*(J12L*J22L*PDstandardNth12alpha + 
      J12L*J32L*PDstandardNth13alpha + J22L*J32L*PDstandardNth23alpha) + 
      (-dJ222L + Gt122*J21L + 4*cdphi2*J22L + Gt222*J22L + 
      Gt322*J23L)*PDstandardNth2alpha + (-dJ322L + Gt122*J31L + 4*cdphi2*J32L 
      + Gt222*J32L + Gt322*J33L)*PDstandardNth3alpha + alphaL*R22 - 
      PDstandardNth11alpha*SQR(J12L) - PDstandardNth22alpha*SQR(J22L) - 
      PDstandardNth33alpha*SQR(J32L);
    
    CCTK_REAL Ats23 = J13L*(-(J12L*PDstandardNth11alpha) - 
      J22L*PDstandardNth12alpha - J32L*PDstandardNth13alpha) + 
      J12L*(-(J23L*PDstandardNth12alpha) - J33L*PDstandardNth13alpha) + 
      (-dJ123L + Gt123*J11L + 2*cdphi3*J12L + Gt223*J12L + 2*cdphi2*J13L + 
      Gt323*J13L)*PDstandardNth1alpha + J23L*(-(J22L*PDstandardNth22alpha) - 
      J32L*PDstandardNth23alpha) + (-dJ223L + Gt123*J21L + 2*cdphi3*J22L + 
      Gt223*J22L + 2*cdphi2*J23L + Gt323*J23L)*PDstandardNth2alpha + 
      J33L*(-(J22L*PDstandardNth23alpha) - J32L*PDstandardNth33alpha) + 
      (-dJ323L + Gt123*J31L)*PDstandardNth3alpha + 
      2*cdphi3*J32L*PDstandardNth3alpha + Gt223*J32L*PDstandardNth3alpha + 
      2*cdphi2*J33L*PDstandardNth3alpha + Gt323*J33L*PDstandardNth3alpha + 
      alphaL*R23;
    
    CCTK_REAL Ats33 = (-dJ133L + Gt133*J11L + Gt233*J12L + 4*cdphi3*J13L + 
      Gt333*J13L)*PDstandardNth1alpha - 2*(J13L*J23L*PDstandardNth12alpha + 
      J13L*J33L*PDstandardNth13alpha + J23L*J33L*PDstandardNth23alpha) + 
      (-dJ233L + Gt133*J21L + Gt233*J22L + 4*cdphi3*J23L + 
      Gt333*J23L)*PDstandardNth2alpha + (-dJ333L + Gt133*J31L + Gt233*J32L + 
      4*cdphi3*J33L + Gt333*J33L)*PDstandardNth3alpha + alphaL*R33 - 
      PDstandardNth11alpha*SQR(J13L) - PDstandardNth22alpha*SQR(J23L) - 
      PDstandardNth33alpha*SQR(J33L);
    
    CCTK_REAL trAts = Ats11*gu11 + Ats22*gu22 + 2*(Ats12*gu12 + Ats13*gu13 
      + Ats23*gu23) + Ats33*gu33;
    
    CCTK_REAL At11rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At11 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At11 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At11 - 
      0.666666666666666666666666666667*At11L*(J11L*PDstandardNth1beta1 + 
      J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
      J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + 
      J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
      J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + 
      2*(J11L*(At11L*PDstandardNth1beta1 + At12L*PDstandardNth1beta2 + 
      At13L*PDstandardNth1beta3) + J21L*(At11L*PDstandardNth2beta1 + 
      At12L*PDstandardNth2beta2 + At13L*PDstandardNth2beta3) + 
      J31L*(At11L*PDstandardNth3beta1 + At12L*PDstandardNth3beta2 + 
      At13L*PDstandardNth3beta3)) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1At11 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2At11 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3At11 + alphaL*(-2*(At11L*Atm11 + 
      At12L*Atm21 + At13L*Atm31) + At11L*trKL) + em4phi*(Ats11 - 
      0.333333333333333333333333333333*g11*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTxxL - 
      0.333333333333333333333333333333*g11*trS)) + (J11L*PDupwindNthSymm1At11 
      + J21L*PDupwindNthSymm2At11 + J31L*PDupwindNthSymm3At11)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At11 + J22L*PDupwindNthSymm2At11 + 
      J32L*PDupwindNthSymm3At11)*Abs(beta2L) + (J13L*PDupwindNthSymm1At11 + 
      J23L*PDupwindNthSymm2At11 + J33L*PDupwindNthSymm3At11)*Abs(beta3L);
    
    CCTK_REAL At12rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At12 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At12 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At12 + (At12L*J11L + 
      At11L*J12L)*PDstandardNth1beta1 + (At22L*J11L + 
      At12L*J12L)*PDstandardNth1beta2 + (At23L*J11L + 
      At13L*J12L)*PDstandardNth1beta3 + (At12L*J21L + 
      At11L*J22L)*PDstandardNth2beta1 + (At22L*J21L + 
      At12L*J22L)*PDstandardNth2beta2 + (At23L*J21L + 
      At13L*J22L)*PDstandardNth2beta3 + (At12L*J31L + 
      At11L*J32L)*PDstandardNth3beta1 + (At22L*J31L + 
      At12L*J32L)*PDstandardNth3beta2 + (At23L*J31L + 
      At13L*J32L)*PDstandardNth3beta3 - 
      0.666666666666666666666666666667*At12L*(J11L*PDstandardNth1beta1 + 
      J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
      J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + 
      J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
      J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + (beta1L*J11L + 
      beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1At12 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2At12 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3At12 + 
      alphaL*(-2*(At11L*Atm12 + At12L*Atm22 + At13L*Atm32) + At12L*trKL) + 
      em4phi*(Ats12 - 0.333333333333333333333333333333*g12*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTxyL - 
      0.333333333333333333333333333333*g12*trS)) + (J11L*PDupwindNthSymm1At12 
      + J21L*PDupwindNthSymm2At12 + J31L*PDupwindNthSymm3At12)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At12 + J22L*PDupwindNthSymm2At12 + 
      J32L*PDupwindNthSymm3At12)*Abs(beta2L) + (J13L*PDupwindNthSymm1At12 + 
      J23L*PDupwindNthSymm2At12 + J33L*PDupwindNthSymm3At12)*Abs(beta3L);
    
    CCTK_REAL At13rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At13 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At13 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At13 + (At13L*J11L + 
      At11L*J13L)*PDstandardNth1beta1 + (At23L*J11L + 
      At12L*J13L)*PDstandardNth1beta2 + (At33L*J11L + 
      At13L*J13L)*PDstandardNth1beta3 + (At13L*J21L + 
      At11L*J23L)*PDstandardNth2beta1 + (At23L*J21L + 
      At12L*J23L)*PDstandardNth2beta2 + (At33L*J21L + 
      At13L*J23L)*PDstandardNth2beta3 + (At13L*J31L + 
      At11L*J33L)*PDstandardNth3beta1 + (At23L*J31L + 
      At12L*J33L)*PDstandardNth3beta2 + (At33L*J31L + 
      At13L*J33L)*PDstandardNth3beta3 - 
      0.666666666666666666666666666667*At13L*(J11L*PDstandardNth1beta1 + 
      J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
      J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + 
      J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
      J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + (beta1L*J11L + 
      beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1At13 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2At13 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3At13 + 
      alphaL*(-2*(At11L*Atm13 + At12L*Atm23 + At13L*Atm33) + At13L*trKL) + 
      em4phi*(Ats13 - 0.333333333333333333333333333333*g13*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTxzL - 
      0.333333333333333333333333333333*g13*trS)) + (J11L*PDupwindNthSymm1At13 
      + J21L*PDupwindNthSymm2At13 + J31L*PDupwindNthSymm3At13)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At13 + J22L*PDupwindNthSymm2At13 + 
      J32L*PDupwindNthSymm3At13)*Abs(beta2L) + (J13L*PDupwindNthSymm1At13 + 
      J23L*PDupwindNthSymm2At13 + J33L*PDupwindNthSymm3At13)*Abs(beta3L);
    
    CCTK_REAL At22rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At22 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At22 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At22 - 
      0.666666666666666666666666666667*At22L*(J11L*PDstandardNth1beta1 + 
      J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
      J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + 
      J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
      J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + 
      2*(J12L*(At12L*PDstandardNth1beta1 + At22L*PDstandardNth1beta2 + 
      At23L*PDstandardNth1beta3) + J22L*(At12L*PDstandardNth2beta1 + 
      At22L*PDstandardNth2beta2 + At23L*PDstandardNth2beta3) + 
      J32L*(At12L*PDstandardNth3beta1 + At22L*PDstandardNth3beta2 + 
      At23L*PDstandardNth3beta3)) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1At22 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2At22 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3At22 + alphaL*(-2*(At12L*Atm12 + 
      At22L*Atm22 + At23L*Atm32) + At22L*trKL) + em4phi*(Ats22 - 
      0.333333333333333333333333333333*g22*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTyyL - 
      0.333333333333333333333333333333*g22*trS)) + (J11L*PDupwindNthSymm1At22 
      + J21L*PDupwindNthSymm2At22 + J31L*PDupwindNthSymm3At22)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At22 + J22L*PDupwindNthSymm2At22 + 
      J32L*PDupwindNthSymm3At22)*Abs(beta2L) + (J13L*PDupwindNthSymm1At22 + 
      J23L*PDupwindNthSymm2At22 + J33L*PDupwindNthSymm3At22)*Abs(beta3L);
    
    CCTK_REAL At23rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At23 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At23 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At23 + (At13L*J12L + 
      At12L*J13L)*PDstandardNth1beta1 + (At23L*J12L + 
      At22L*J13L)*PDstandardNth1beta2 + (At33L*J12L + 
      At23L*J13L)*PDstandardNth1beta3 + (At13L*J22L + 
      At12L*J23L)*PDstandardNth2beta1 + (At23L*J22L + 
      At22L*J23L)*PDstandardNth2beta2 + (At33L*J22L + 
      At23L*J23L)*PDstandardNth2beta3 + (At13L*J32L + 
      At12L*J33L)*PDstandardNth3beta1 + (At23L*J32L + 
      At22L*J33L)*PDstandardNth3beta2 + (At33L*J32L + 
      At23L*J33L)*PDstandardNth3beta3 - 
      0.666666666666666666666666666667*At23L*(J11L*PDstandardNth1beta1 + 
      J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
      J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + 
      J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
      J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + (beta1L*J11L + 
      beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1At23 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2At23 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3At23 + 
      alphaL*(-2*(At12L*Atm13 + At22L*Atm23 + At23L*Atm33) + At23L*trKL) + 
      em4phi*(Ats23 - 0.333333333333333333333333333333*g23*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTyzL - 
      0.333333333333333333333333333333*g23*trS)) + (J11L*PDupwindNthSymm1At23 
      + J21L*PDupwindNthSymm2At23 + J31L*PDupwindNthSymm3At23)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At23 + J22L*PDupwindNthSymm2At23 + 
      J32L*PDupwindNthSymm3At23)*Abs(beta2L) + (J13L*PDupwindNthSymm1At23 + 
      J23L*PDupwindNthSymm2At23 + J33L*PDupwindNthSymm3At23)*Abs(beta3L);
    
    CCTK_REAL At33rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At33 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At33 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At33 - 
      0.666666666666666666666666666667*At33L*(J11L*PDstandardNth1beta1 + 
      J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
      J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + 
      J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
      J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3) + 
      2*(J13L*(At13L*PDstandardNth1beta1 + At23L*PDstandardNth1beta2 + 
      At33L*PDstandardNth1beta3) + J23L*(At13L*PDstandardNth2beta1 + 
      At23L*PDstandardNth2beta2 + At33L*PDstandardNth2beta3) + 
      J33L*(At13L*PDstandardNth3beta1 + At23L*PDstandardNth3beta2 + 
      At33L*PDstandardNth3beta3)) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1At33 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2At33 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3At33 + alphaL*(-2*(At13L*Atm13 + 
      At23L*Atm23 + At33L*Atm33) + At33L*trKL) + em4phi*(Ats33 - 
      0.333333333333333333333333333333*g33*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTzzL - 
      0.333333333333333333333333333333*g33*trS)) + (J11L*PDupwindNthSymm1At33 
      + J21L*PDupwindNthSymm2At33 + J31L*PDupwindNthSymm3At33)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At33 + J22L*PDupwindNthSymm2At33 + 
      J32L*PDupwindNthSymm3At33)*Abs(beta2L) + (J13L*PDupwindNthSymm1At33 + 
      J23L*PDupwindNthSymm2At33 + J33L*PDupwindNthSymm3At33)*Abs(beta3L);
    
    /* Copy local copies back to grid functions */
    At11rhs[index] = At11rhsL;
    At12rhs[index] = At12rhsL;
    At13rhs[index] = At13rhsL;
    At22rhs[index] = At22rhsL;
    At23rhs[index] = At23rhsL;
    At33rhs[index] = At33rhsL;
  }
  LC_ENDLOOP3 (ML_BSSN_MP_RHS2);
}

extern "C" void ML_BSSN_MP_RHS2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_RHS2_Body);
}
