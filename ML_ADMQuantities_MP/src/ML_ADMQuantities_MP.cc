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

extern "C" void ML_ADMQuantities_MP_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMQuantities_MP::ML_Jadm","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMQuantities_MP::ML_Jadm.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMQuantities_MP::ML_Madm","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMQuantities_MP::ML_Madm.");
  return;
}

static void ML_ADMQuantities_MP_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  
  /* Jacobian variable pointers */
  bool const use_jacobian = (!CCTK_IsFunctionAliased("MultiPatch_GetMap") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)
                       && strlen(jacobian_group) > 0;
  if (use_jacobian && strlen(jacobian_derivative_group) == 0)
  {
    CCTK_WARN (1, "GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names");
  }
  
  CCTK_REAL const *restrict jacobian_ptrs[9];
  if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_group,
                                                9, jacobian_ptrs);
  
  CCTK_REAL const *restrict const J11 = use_jacobian ? jacobian_ptrs[0] : 0;
  CCTK_REAL const *restrict const J12 = use_jacobian ? jacobian_ptrs[1] : 0;
  CCTK_REAL const *restrict const J13 = use_jacobian ? jacobian_ptrs[2] : 0;
  CCTK_REAL const *restrict const J21 = use_jacobian ? jacobian_ptrs[3] : 0;
  CCTK_REAL const *restrict const J22 = use_jacobian ? jacobian_ptrs[4] : 0;
  CCTK_REAL const *restrict const J23 = use_jacobian ? jacobian_ptrs[5] : 0;
  CCTK_REAL const *restrict const J31 = use_jacobian ? jacobian_ptrs[6] : 0;
  CCTK_REAL const *restrict const J32 = use_jacobian ? jacobian_ptrs[7] : 0;
  CCTK_REAL const *restrict const J33 = use_jacobian ? jacobian_ptrs[8] : 0;
  
  CCTK_REAL const *restrict jacobian_derivative_ptrs[18];
  if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_derivative_group,
                                                18, jacobian_derivative_ptrs);
  
  CCTK_REAL const *restrict const dJ111 = use_jacobian ? jacobian_derivative_ptrs[0] : 0;
  CCTK_REAL const *restrict const dJ112 = use_jacobian ? jacobian_derivative_ptrs[1] : 0;
  CCTK_REAL const *restrict const dJ113 = use_jacobian ? jacobian_derivative_ptrs[2] : 0;
  CCTK_REAL const *restrict const dJ122 = use_jacobian ? jacobian_derivative_ptrs[3] : 0;
  CCTK_REAL const *restrict const dJ123 = use_jacobian ? jacobian_derivative_ptrs[4] : 0;
  CCTK_REAL const *restrict const dJ133 = use_jacobian ? jacobian_derivative_ptrs[5] : 0;
  CCTK_REAL const *restrict const dJ211 = use_jacobian ? jacobian_derivative_ptrs[6] : 0;
  CCTK_REAL const *restrict const dJ212 = use_jacobian ? jacobian_derivative_ptrs[7] : 0;
  CCTK_REAL const *restrict const dJ213 = use_jacobian ? jacobian_derivative_ptrs[8] : 0;
  CCTK_REAL const *restrict const dJ222 = use_jacobian ? jacobian_derivative_ptrs[9] : 0;
  CCTK_REAL const *restrict const dJ223 = use_jacobian ? jacobian_derivative_ptrs[10] : 0;
  CCTK_REAL const *restrict const dJ233 = use_jacobian ? jacobian_derivative_ptrs[11] : 0;
  CCTK_REAL const *restrict const dJ311 = use_jacobian ? jacobian_derivative_ptrs[12] : 0;
  CCTK_REAL const *restrict const dJ312 = use_jacobian ? jacobian_derivative_ptrs[13] : 0;
  CCTK_REAL const *restrict const dJ313 = use_jacobian ? jacobian_derivative_ptrs[14] : 0;
  CCTK_REAL const *restrict const dJ322 = use_jacobian ? jacobian_derivative_ptrs[15] : 0;
  CCTK_REAL const *restrict const dJ323 = use_jacobian ? jacobian_derivative_ptrs[16] : 0;
  CCTK_REAL const *restrict const dJ333 = use_jacobian ? jacobian_derivative_ptrs[17] : 0;
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  CCTK_LOOP3(ML_ADMQuantities_MP,
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
    
    CCTK_REAL dJ111L, dJ112L, dJ113L, dJ122L, dJ123L, dJ133L, dJ211L, dJ212L, dJ213L, dJ222L, dJ223L, dJ233L, dJ311L, dJ312L, dJ313L, dJ322L, dJ323L, dJ333L, J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L;
    
    if (use_jacobian)
    {
      dJ111L = dJ111[index];
      dJ112L = dJ112[index];
      dJ113L = dJ113[index];
      dJ122L = dJ122[index];
      dJ123L = dJ123[index];
      dJ133L = dJ133[index];
      dJ211L = dJ211[index];
      dJ212L = dJ212[index];
      dJ213L = dJ213[index];
      dJ222L = dJ222[index];
      dJ223L = dJ223[index];
      dJ233L = dJ233[index];
      dJ311L = dJ311[index];
      dJ312L = dJ312[index];
      dJ313L = dJ313[index];
      dJ322L = dJ322[index];
      dJ323L = dJ323[index];
      dJ333L = dJ333[index];
      J11L = J11[index];
      J12L = J12[index];
      J13L = J13[index];
      J21L = J21[index];
      J22L = J22[index];
      J23L = J23[index];
      J31L = J31[index];
      J32L = J32[index];
      J33L = J33[index];
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
    CCTK_REAL JacPDstandardNth11gt11;
    CCTK_REAL JacPDstandardNth11gt12;
    CCTK_REAL JacPDstandardNth11gt13;
    CCTK_REAL JacPDstandardNth11gt22;
    CCTK_REAL JacPDstandardNth11gt23;
    CCTK_REAL JacPDstandardNth11gt33;
    CCTK_REAL JacPDstandardNth12gt11;
    CCTK_REAL JacPDstandardNth12gt12;
    CCTK_REAL JacPDstandardNth12gt13;
    CCTK_REAL JacPDstandardNth12gt22;
    CCTK_REAL JacPDstandardNth12gt23;
    CCTK_REAL JacPDstandardNth12gt33;
    CCTK_REAL JacPDstandardNth13gt11;
    CCTK_REAL JacPDstandardNth13gt12;
    CCTK_REAL JacPDstandardNth13gt13;
    CCTK_REAL JacPDstandardNth13gt22;
    CCTK_REAL JacPDstandardNth13gt23;
    CCTK_REAL JacPDstandardNth13gt33;
    CCTK_REAL JacPDstandardNth1gt11;
    CCTK_REAL JacPDstandardNth1gt12;
    CCTK_REAL JacPDstandardNth1gt13;
    CCTK_REAL JacPDstandardNth1gt22;
    CCTK_REAL JacPDstandardNth1gt23;
    CCTK_REAL JacPDstandardNth1gt33;
    CCTK_REAL JacPDstandardNth1trK;
    CCTK_REAL JacPDstandardNth1Xt1;
    CCTK_REAL JacPDstandardNth1Xt2;
    CCTK_REAL JacPDstandardNth1Xt3;
    CCTK_REAL JacPDstandardNth21gt11;
    CCTK_REAL JacPDstandardNth21gt12;
    CCTK_REAL JacPDstandardNth21gt13;
    CCTK_REAL JacPDstandardNth21gt22;
    CCTK_REAL JacPDstandardNth21gt23;
    CCTK_REAL JacPDstandardNth21gt33;
    CCTK_REAL JacPDstandardNth22gt11;
    CCTK_REAL JacPDstandardNth22gt12;
    CCTK_REAL JacPDstandardNth22gt13;
    CCTK_REAL JacPDstandardNth22gt22;
    CCTK_REAL JacPDstandardNth22gt23;
    CCTK_REAL JacPDstandardNth22gt33;
    CCTK_REAL JacPDstandardNth23gt11;
    CCTK_REAL JacPDstandardNth23gt12;
    CCTK_REAL JacPDstandardNth23gt13;
    CCTK_REAL JacPDstandardNth23gt22;
    CCTK_REAL JacPDstandardNth23gt23;
    CCTK_REAL JacPDstandardNth23gt33;
    CCTK_REAL JacPDstandardNth2gt11;
    CCTK_REAL JacPDstandardNth2gt12;
    CCTK_REAL JacPDstandardNth2gt13;
    CCTK_REAL JacPDstandardNth2gt22;
    CCTK_REAL JacPDstandardNth2gt23;
    CCTK_REAL JacPDstandardNth2gt33;
    CCTK_REAL JacPDstandardNth2trK;
    CCTK_REAL JacPDstandardNth2Xt1;
    CCTK_REAL JacPDstandardNth2Xt2;
    CCTK_REAL JacPDstandardNth2Xt3;
    CCTK_REAL JacPDstandardNth31gt11;
    CCTK_REAL JacPDstandardNth31gt12;
    CCTK_REAL JacPDstandardNth31gt13;
    CCTK_REAL JacPDstandardNth31gt22;
    CCTK_REAL JacPDstandardNth31gt23;
    CCTK_REAL JacPDstandardNth31gt33;
    CCTK_REAL JacPDstandardNth32gt11;
    CCTK_REAL JacPDstandardNth32gt12;
    CCTK_REAL JacPDstandardNth32gt13;
    CCTK_REAL JacPDstandardNth32gt22;
    CCTK_REAL JacPDstandardNth32gt23;
    CCTK_REAL JacPDstandardNth32gt33;
    CCTK_REAL JacPDstandardNth33gt11;
    CCTK_REAL JacPDstandardNth33gt12;
    CCTK_REAL JacPDstandardNth33gt13;
    CCTK_REAL JacPDstandardNth33gt22;
    CCTK_REAL JacPDstandardNth33gt23;
    CCTK_REAL JacPDstandardNth33gt33;
    CCTK_REAL JacPDstandardNth3gt11;
    CCTK_REAL JacPDstandardNth3gt12;
    CCTK_REAL JacPDstandardNth3gt13;
    CCTK_REAL JacPDstandardNth3gt22;
    CCTK_REAL JacPDstandardNth3gt23;
    CCTK_REAL JacPDstandardNth3gt33;
    CCTK_REAL JacPDstandardNth3trK;
    CCTK_REAL JacPDstandardNth3Xt1;
    CCTK_REAL JacPDstandardNth3Xt2;
    CCTK_REAL JacPDstandardNth3Xt3;
    
    if (use_jacobian)
    {
      JacPDstandardNth1gt11 = J11L*PDstandardNth1gt11 + 
        J21L*PDstandardNth2gt11 + J31L*PDstandardNth3gt11;
      
      JacPDstandardNth1gt12 = J11L*PDstandardNth1gt12 + 
        J21L*PDstandardNth2gt12 + J31L*PDstandardNth3gt12;
      
      JacPDstandardNth1gt13 = J11L*PDstandardNth1gt13 + 
        J21L*PDstandardNth2gt13 + J31L*PDstandardNth3gt13;
      
      JacPDstandardNth1gt22 = J11L*PDstandardNth1gt22 + 
        J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22;
      
      JacPDstandardNth1gt23 = J11L*PDstandardNth1gt23 + 
        J21L*PDstandardNth2gt23 + J31L*PDstandardNth3gt23;
      
      JacPDstandardNth1gt33 = J11L*PDstandardNth1gt33 + 
        J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33;
      
      JacPDstandardNth1trK = J11L*PDstandardNth1trK + 
        J21L*PDstandardNth2trK + J31L*PDstandardNth3trK;
      
      JacPDstandardNth1Xt1 = J11L*PDstandardNth1Xt1 + 
        J21L*PDstandardNth2Xt1 + J31L*PDstandardNth3Xt1;
      
      JacPDstandardNth1Xt2 = J11L*PDstandardNth1Xt2 + 
        J21L*PDstandardNth2Xt2 + J31L*PDstandardNth3Xt2;
      
      JacPDstandardNth1Xt3 = J11L*PDstandardNth1Xt3 + 
        J21L*PDstandardNth2Xt3 + J31L*PDstandardNth3Xt3;
      
      JacPDstandardNth2gt11 = J12L*PDstandardNth1gt11 + 
        J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11;
      
      JacPDstandardNth2gt12 = J12L*PDstandardNth1gt12 + 
        J22L*PDstandardNth2gt12 + J32L*PDstandardNth3gt12;
      
      JacPDstandardNth2gt13 = J12L*PDstandardNth1gt13 + 
        J22L*PDstandardNth2gt13 + J32L*PDstandardNth3gt13;
      
      JacPDstandardNth2gt22 = J12L*PDstandardNth1gt22 + 
        J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22;
      
      JacPDstandardNth2gt23 = J12L*PDstandardNth1gt23 + 
        J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23;
      
      JacPDstandardNth2gt33 = J12L*PDstandardNth1gt33 + 
        J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33;
      
      JacPDstandardNth2trK = J12L*PDstandardNth1trK + 
        J22L*PDstandardNth2trK + J32L*PDstandardNth3trK;
      
      JacPDstandardNth2Xt1 = J12L*PDstandardNth1Xt1 + 
        J22L*PDstandardNth2Xt1 + J32L*PDstandardNth3Xt1;
      
      JacPDstandardNth2Xt2 = J12L*PDstandardNth1Xt2 + 
        J22L*PDstandardNth2Xt2 + J32L*PDstandardNth3Xt2;
      
      JacPDstandardNth2Xt3 = J12L*PDstandardNth1Xt3 + 
        J22L*PDstandardNth2Xt3 + J32L*PDstandardNth3Xt3;
      
      JacPDstandardNth3gt11 = J13L*PDstandardNth1gt11 + 
        J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = J13L*PDstandardNth1gt12 + 
        J23L*PDstandardNth2gt12 + J33L*PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = J13L*PDstandardNth1gt13 + 
        J23L*PDstandardNth2gt13 + J33L*PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = J13L*PDstandardNth1gt22 + 
        J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = J13L*PDstandardNth1gt23 + 
        J23L*PDstandardNth2gt23 + J33L*PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = J13L*PDstandardNth1gt33 + 
        J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33;
      
      JacPDstandardNth3trK = J13L*PDstandardNth1trK + 
        J23L*PDstandardNth2trK + J33L*PDstandardNth3trK;
      
      JacPDstandardNth3Xt1 = J13L*PDstandardNth1Xt1 + 
        J23L*PDstandardNth2Xt1 + J33L*PDstandardNth3Xt1;
      
      JacPDstandardNth3Xt2 = J13L*PDstandardNth1Xt2 + 
        J23L*PDstandardNth2Xt2 + J33L*PDstandardNth3Xt2;
      
      JacPDstandardNth3Xt3 = J13L*PDstandardNth1Xt3 + 
        J23L*PDstandardNth2Xt3 + J33L*PDstandardNth3Xt3;
      
      JacPDstandardNth11gt11 = dJ111L*PDstandardNth1gt11 + 
        2*(J11L*(J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J21L*J31L*PDstandardNth23gt11) + dJ211L*PDstandardNth2gt11 + 
        dJ311L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J11L) + 
        PDstandardNth22gt11*SQR(J21L) + PDstandardNth33gt11*SQR(J31L);
      
      JacPDstandardNth11gt12 = dJ111L*PDstandardNth1gt12 + 
        2*(J11L*(J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J21L*J31L*PDstandardNth23gt12) + dJ211L*PDstandardNth2gt12 + 
        dJ311L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J11L) + 
        PDstandardNth22gt12*SQR(J21L) + PDstandardNth33gt12*SQR(J31L);
      
      JacPDstandardNth11gt13 = dJ111L*PDstandardNth1gt13 + 
        2*(J11L*(J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J21L*J31L*PDstandardNth23gt13) + dJ211L*PDstandardNth2gt13 + 
        dJ311L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J11L) + 
        PDstandardNth22gt13*SQR(J21L) + PDstandardNth33gt13*SQR(J31L);
      
      JacPDstandardNth11gt22 = dJ111L*PDstandardNth1gt22 + 
        2*(J11L*(J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J21L*J31L*PDstandardNth23gt22) + dJ211L*PDstandardNth2gt22 + 
        dJ311L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J11L) + 
        PDstandardNth22gt22*SQR(J21L) + PDstandardNth33gt22*SQR(J31L);
      
      JacPDstandardNth11gt23 = dJ111L*PDstandardNth1gt23 + 
        2*(J11L*(J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J21L*J31L*PDstandardNth23gt23) + dJ211L*PDstandardNth2gt23 + 
        dJ311L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J11L) + 
        PDstandardNth22gt23*SQR(J21L) + PDstandardNth33gt23*SQR(J31L);
      
      JacPDstandardNth11gt33 = dJ111L*PDstandardNth1gt33 + 
        2*(J11L*(J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J21L*J31L*PDstandardNth23gt33) + dJ211L*PDstandardNth2gt33 + 
        dJ311L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J11L) + 
        PDstandardNth22gt33*SQR(J21L) + PDstandardNth33gt33*SQR(J31L);
      
      JacPDstandardNth22gt11 = dJ122L*PDstandardNth1gt11 + 
        2*(J12L*(J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        J22L*J32L*PDstandardNth23gt11) + dJ222L*PDstandardNth2gt11 + 
        dJ322L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J12L) + 
        PDstandardNth22gt11*SQR(J22L) + PDstandardNth33gt11*SQR(J32L);
      
      JacPDstandardNth22gt12 = dJ122L*PDstandardNth1gt12 + 
        2*(J12L*(J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        J22L*J32L*PDstandardNth23gt12) + dJ222L*PDstandardNth2gt12 + 
        dJ322L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J12L) + 
        PDstandardNth22gt12*SQR(J22L) + PDstandardNth33gt12*SQR(J32L);
      
      JacPDstandardNth22gt13 = dJ122L*PDstandardNth1gt13 + 
        2*(J12L*(J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        J22L*J32L*PDstandardNth23gt13) + dJ222L*PDstandardNth2gt13 + 
        dJ322L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J12L) + 
        PDstandardNth22gt13*SQR(J22L) + PDstandardNth33gt13*SQR(J32L);
      
      JacPDstandardNth22gt22 = dJ122L*PDstandardNth1gt22 + 
        2*(J12L*(J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        J22L*J32L*PDstandardNth23gt22) + dJ222L*PDstandardNth2gt22 + 
        dJ322L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J12L) + 
        PDstandardNth22gt22*SQR(J22L) + PDstandardNth33gt22*SQR(J32L);
      
      JacPDstandardNth22gt23 = dJ122L*PDstandardNth1gt23 + 
        2*(J12L*(J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        J22L*J32L*PDstandardNth23gt23) + dJ222L*PDstandardNth2gt23 + 
        dJ322L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J12L) + 
        PDstandardNth22gt23*SQR(J22L) + PDstandardNth33gt23*SQR(J32L);
      
      JacPDstandardNth22gt33 = dJ122L*PDstandardNth1gt33 + 
        2*(J12L*(J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        J22L*J32L*PDstandardNth23gt33) + dJ222L*PDstandardNth2gt33 + 
        dJ322L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J12L) + 
        PDstandardNth22gt33*SQR(J22L) + PDstandardNth33gt33*SQR(J32L);
      
      JacPDstandardNth33gt11 = dJ133L*PDstandardNth1gt11 + 
        2*(J13L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        J23L*J33L*PDstandardNth23gt11) + dJ233L*PDstandardNth2gt11 + 
        dJ333L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J13L) + 
        PDstandardNth22gt11*SQR(J23L) + PDstandardNth33gt11*SQR(J33L);
      
      JacPDstandardNth33gt12 = dJ133L*PDstandardNth1gt12 + 
        2*(J13L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        J23L*J33L*PDstandardNth23gt12) + dJ233L*PDstandardNth2gt12 + 
        dJ333L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J13L) + 
        PDstandardNth22gt12*SQR(J23L) + PDstandardNth33gt12*SQR(J33L);
      
      JacPDstandardNth33gt13 = dJ133L*PDstandardNth1gt13 + 
        2*(J13L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        J23L*J33L*PDstandardNth23gt13) + dJ233L*PDstandardNth2gt13 + 
        dJ333L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J13L) + 
        PDstandardNth22gt13*SQR(J23L) + PDstandardNth33gt13*SQR(J33L);
      
      JacPDstandardNth33gt22 = dJ133L*PDstandardNth1gt22 + 
        2*(J13L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        J23L*J33L*PDstandardNth23gt22) + dJ233L*PDstandardNth2gt22 + 
        dJ333L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J13L) + 
        PDstandardNth22gt22*SQR(J23L) + PDstandardNth33gt22*SQR(J33L);
      
      JacPDstandardNth33gt23 = dJ133L*PDstandardNth1gt23 + 
        2*(J13L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        J23L*J33L*PDstandardNth23gt23) + dJ233L*PDstandardNth2gt23 + 
        dJ333L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J13L) + 
        PDstandardNth22gt23*SQR(J23L) + PDstandardNth33gt23*SQR(J33L);
      
      JacPDstandardNth33gt33 = dJ133L*PDstandardNth1gt33 + 
        2*(J13L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        J23L*J33L*PDstandardNth23gt33) + dJ233L*PDstandardNth2gt33 + 
        dJ333L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J13L) + 
        PDstandardNth22gt33*SQR(J23L) + PDstandardNth33gt33*SQR(J33L);
      
      JacPDstandardNth12gt11 = J12L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        dJ112L*PDstandardNth1gt11 + J22L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ212L*PDstandardNth2gt11 + 
        J32L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ312L*PDstandardNth3gt11;
      
      JacPDstandardNth12gt12 = J12L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        dJ112L*PDstandardNth1gt12 + J22L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ212L*PDstandardNth2gt12 + 
        J32L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ312L*PDstandardNth3gt12;
      
      JacPDstandardNth12gt13 = J12L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        dJ112L*PDstandardNth1gt13 + J22L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ212L*PDstandardNth2gt13 + 
        J32L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ312L*PDstandardNth3gt13;
      
      JacPDstandardNth12gt22 = J12L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        dJ112L*PDstandardNth1gt22 + J22L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ212L*PDstandardNth2gt22 + 
        J32L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ312L*PDstandardNth3gt22;
      
      JacPDstandardNth12gt23 = J12L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        dJ112L*PDstandardNth1gt23 + J22L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ212L*PDstandardNth2gt23 + 
        J32L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ312L*PDstandardNth3gt23;
      
      JacPDstandardNth12gt33 = J12L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        dJ112L*PDstandardNth1gt33 + J22L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ212L*PDstandardNth2gt33 + 
        J32L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ312L*PDstandardNth3gt33;
      
      JacPDstandardNth13gt11 = J13L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ113L*PDstandardNth1gt11 + J23L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ213L*PDstandardNth2gt11 + 
        J33L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ313L*PDstandardNth3gt11;
      
      JacPDstandardNth13gt12 = J13L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ113L*PDstandardNth1gt12 + J23L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ213L*PDstandardNth2gt12 + 
        J33L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ313L*PDstandardNth3gt12;
      
      JacPDstandardNth13gt13 = J13L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ113L*PDstandardNth1gt13 + J23L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ213L*PDstandardNth2gt13 + 
        J33L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ313L*PDstandardNth3gt13;
      
      JacPDstandardNth13gt22 = J13L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ113L*PDstandardNth1gt22 + J23L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ213L*PDstandardNth2gt22 + 
        J33L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ313L*PDstandardNth3gt22;
      
      JacPDstandardNth13gt23 = J13L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ113L*PDstandardNth1gt23 + J23L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ213L*PDstandardNth2gt23 + 
        J33L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ313L*PDstandardNth3gt23;
      
      JacPDstandardNth13gt33 = J13L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ113L*PDstandardNth1gt33 + J23L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ213L*PDstandardNth2gt33 + 
        J33L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ313L*PDstandardNth3gt33;
      
      JacPDstandardNth21gt11 = J12L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        dJ112L*PDstandardNth1gt11 + J22L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ212L*PDstandardNth2gt11 + 
        J32L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ312L*PDstandardNth3gt11;
      
      JacPDstandardNth21gt12 = J12L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        dJ112L*PDstandardNth1gt12 + J22L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ212L*PDstandardNth2gt12 + 
        J32L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ312L*PDstandardNth3gt12;
      
      JacPDstandardNth21gt13 = J12L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        dJ112L*PDstandardNth1gt13 + J22L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ212L*PDstandardNth2gt13 + 
        J32L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ312L*PDstandardNth3gt13;
      
      JacPDstandardNth21gt22 = J12L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        dJ112L*PDstandardNth1gt22 + J22L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ212L*PDstandardNth2gt22 + 
        J32L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ312L*PDstandardNth3gt22;
      
      JacPDstandardNth21gt23 = J12L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        dJ112L*PDstandardNth1gt23 + J22L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ212L*PDstandardNth2gt23 + 
        J32L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ312L*PDstandardNth3gt23;
      
      JacPDstandardNth21gt33 = J12L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        dJ112L*PDstandardNth1gt33 + J22L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ212L*PDstandardNth2gt33 + 
        J32L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ312L*PDstandardNth3gt33;
      
      JacPDstandardNth23gt11 = J13L*(J12L*PDstandardNth11gt11 + 
        J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        J12L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ123L*PDstandardNth1gt11 + J23L*(J22L*PDstandardNth22gt11 + 
        J32L*PDstandardNth23gt11) + dJ223L*PDstandardNth2gt11 + 
        J33L*(J22L*PDstandardNth23gt11 + J32L*PDstandardNth33gt11) + 
        dJ323L*PDstandardNth3gt11;
      
      JacPDstandardNth23gt12 = J13L*(J12L*PDstandardNth11gt12 + 
        J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        J12L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ123L*PDstandardNth1gt12 + J23L*(J22L*PDstandardNth22gt12 + 
        J32L*PDstandardNth23gt12) + dJ223L*PDstandardNth2gt12 + 
        J33L*(J22L*PDstandardNth23gt12 + J32L*PDstandardNth33gt12) + 
        dJ323L*PDstandardNth3gt12;
      
      JacPDstandardNth23gt13 = J13L*(J12L*PDstandardNth11gt13 + 
        J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        J12L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ123L*PDstandardNth1gt13 + J23L*(J22L*PDstandardNth22gt13 + 
        J32L*PDstandardNth23gt13) + dJ223L*PDstandardNth2gt13 + 
        J33L*(J22L*PDstandardNth23gt13 + J32L*PDstandardNth33gt13) + 
        dJ323L*PDstandardNth3gt13;
      
      JacPDstandardNth23gt22 = J13L*(J12L*PDstandardNth11gt22 + 
        J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        J12L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ123L*PDstandardNth1gt22 + J23L*(J22L*PDstandardNth22gt22 + 
        J32L*PDstandardNth23gt22) + dJ223L*PDstandardNth2gt22 + 
        J33L*(J22L*PDstandardNth23gt22 + J32L*PDstandardNth33gt22) + 
        dJ323L*PDstandardNth3gt22;
      
      JacPDstandardNth23gt23 = J13L*(J12L*PDstandardNth11gt23 + 
        J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        J12L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ123L*PDstandardNth1gt23 + J23L*(J22L*PDstandardNth22gt23 + 
        J32L*PDstandardNth23gt23) + dJ223L*PDstandardNth2gt23 + 
        J33L*(J22L*PDstandardNth23gt23 + J32L*PDstandardNth33gt23) + 
        dJ323L*PDstandardNth3gt23;
      
      JacPDstandardNth23gt33 = J13L*(J12L*PDstandardNth11gt33 + 
        J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        J12L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ123L*PDstandardNth1gt33 + J23L*(J22L*PDstandardNth22gt33 + 
        J32L*PDstandardNth23gt33) + dJ223L*PDstandardNth2gt33 + 
        J33L*(J22L*PDstandardNth23gt33 + J32L*PDstandardNth33gt33) + 
        dJ323L*PDstandardNth3gt33;
      
      JacPDstandardNth31gt11 = J13L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ113L*PDstandardNth1gt11 + J23L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ213L*PDstandardNth2gt11 + 
        J33L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ313L*PDstandardNth3gt11;
      
      JacPDstandardNth31gt12 = J13L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ113L*PDstandardNth1gt12 + J23L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ213L*PDstandardNth2gt12 + 
        J33L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ313L*PDstandardNth3gt12;
      
      JacPDstandardNth31gt13 = J13L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ113L*PDstandardNth1gt13 + J23L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ213L*PDstandardNth2gt13 + 
        J33L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ313L*PDstandardNth3gt13;
      
      JacPDstandardNth31gt22 = J13L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ113L*PDstandardNth1gt22 + J23L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ213L*PDstandardNth2gt22 + 
        J33L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ313L*PDstandardNth3gt22;
      
      JacPDstandardNth31gt23 = J13L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ113L*PDstandardNth1gt23 + J23L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ213L*PDstandardNth2gt23 + 
        J33L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ313L*PDstandardNth3gt23;
      
      JacPDstandardNth31gt33 = J13L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ113L*PDstandardNth1gt33 + J23L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ213L*PDstandardNth2gt33 + 
        J33L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ313L*PDstandardNth3gt33;
      
      JacPDstandardNth32gt11 = J13L*(J12L*PDstandardNth11gt11 + 
        J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        J12L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ123L*PDstandardNth1gt11 + J23L*(J22L*PDstandardNth22gt11 + 
        J32L*PDstandardNth23gt11) + dJ223L*PDstandardNth2gt11 + 
        J33L*(J22L*PDstandardNth23gt11 + J32L*PDstandardNth33gt11) + 
        dJ323L*PDstandardNth3gt11;
      
      JacPDstandardNth32gt12 = J13L*(J12L*PDstandardNth11gt12 + 
        J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        J12L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ123L*PDstandardNth1gt12 + J23L*(J22L*PDstandardNth22gt12 + 
        J32L*PDstandardNth23gt12) + dJ223L*PDstandardNth2gt12 + 
        J33L*(J22L*PDstandardNth23gt12 + J32L*PDstandardNth33gt12) + 
        dJ323L*PDstandardNth3gt12;
      
      JacPDstandardNth32gt13 = J13L*(J12L*PDstandardNth11gt13 + 
        J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        J12L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ123L*PDstandardNth1gt13 + J23L*(J22L*PDstandardNth22gt13 + 
        J32L*PDstandardNth23gt13) + dJ223L*PDstandardNth2gt13 + 
        J33L*(J22L*PDstandardNth23gt13 + J32L*PDstandardNth33gt13) + 
        dJ323L*PDstandardNth3gt13;
      
      JacPDstandardNth32gt22 = J13L*(J12L*PDstandardNth11gt22 + 
        J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        J12L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ123L*PDstandardNth1gt22 + J23L*(J22L*PDstandardNth22gt22 + 
        J32L*PDstandardNth23gt22) + dJ223L*PDstandardNth2gt22 + 
        J33L*(J22L*PDstandardNth23gt22 + J32L*PDstandardNth33gt22) + 
        dJ323L*PDstandardNth3gt22;
      
      JacPDstandardNth32gt23 = J13L*(J12L*PDstandardNth11gt23 + 
        J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        J12L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ123L*PDstandardNth1gt23 + J23L*(J22L*PDstandardNth22gt23 + 
        J32L*PDstandardNth23gt23) + dJ223L*PDstandardNth2gt23 + 
        J33L*(J22L*PDstandardNth23gt23 + J32L*PDstandardNth33gt23) + 
        dJ323L*PDstandardNth3gt23;
      
      JacPDstandardNth32gt33 = J13L*(J12L*PDstandardNth11gt33 + 
        J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        J12L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ123L*PDstandardNth1gt33 + J23L*(J22L*PDstandardNth22gt33 + 
        J32L*PDstandardNth23gt33) + dJ223L*PDstandardNth2gt33 + 
        J33L*(J22L*PDstandardNth23gt33 + J32L*PDstandardNth33gt33) + 
        dJ323L*PDstandardNth3gt33;
    }
    else
    {
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth1trK = PDstandardNth1trK;
      
      JacPDstandardNth1Xt1 = PDstandardNth1Xt1;
      
      JacPDstandardNth1Xt2 = PDstandardNth1Xt2;
      
      JacPDstandardNth1Xt3 = PDstandardNth1Xt3;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth2trK = PDstandardNth2trK;
      
      JacPDstandardNth2Xt1 = PDstandardNth2Xt1;
      
      JacPDstandardNth2Xt2 = PDstandardNth2Xt2;
      
      JacPDstandardNth2Xt3 = PDstandardNth2Xt3;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDstandardNth3trK = PDstandardNth3trK;
      
      JacPDstandardNth3Xt1 = PDstandardNth3Xt1;
      
      JacPDstandardNth3Xt2 = PDstandardNth3Xt2;
      
      JacPDstandardNth3Xt3 = PDstandardNth3Xt3;
      
      JacPDstandardNth11gt11 = PDstandardNth11gt11;
      
      JacPDstandardNth11gt12 = PDstandardNth11gt12;
      
      JacPDstandardNth11gt13 = PDstandardNth11gt13;
      
      JacPDstandardNth11gt22 = PDstandardNth11gt22;
      
      JacPDstandardNth11gt23 = PDstandardNth11gt23;
      
      JacPDstandardNth11gt33 = PDstandardNth11gt33;
      
      JacPDstandardNth22gt11 = PDstandardNth22gt11;
      
      JacPDstandardNth22gt12 = PDstandardNth22gt12;
      
      JacPDstandardNth22gt13 = PDstandardNth22gt13;
      
      JacPDstandardNth22gt22 = PDstandardNth22gt22;
      
      JacPDstandardNth22gt23 = PDstandardNth22gt23;
      
      JacPDstandardNth22gt33 = PDstandardNth22gt33;
      
      JacPDstandardNth33gt11 = PDstandardNth33gt11;
      
      JacPDstandardNth33gt12 = PDstandardNth33gt12;
      
      JacPDstandardNth33gt13 = PDstandardNth33gt13;
      
      JacPDstandardNth33gt22 = PDstandardNth33gt22;
      
      JacPDstandardNth33gt23 = PDstandardNth33gt23;
      
      JacPDstandardNth33gt33 = PDstandardNth33gt33;
      
      JacPDstandardNth12gt11 = PDstandardNth12gt11;
      
      JacPDstandardNth12gt12 = PDstandardNth12gt12;
      
      JacPDstandardNth12gt13 = PDstandardNth12gt13;
      
      JacPDstandardNth12gt22 = PDstandardNth12gt22;
      
      JacPDstandardNth12gt23 = PDstandardNth12gt23;
      
      JacPDstandardNth12gt33 = PDstandardNth12gt33;
      
      JacPDstandardNth13gt11 = PDstandardNth13gt11;
      
      JacPDstandardNth13gt12 = PDstandardNth13gt12;
      
      JacPDstandardNth13gt13 = PDstandardNth13gt13;
      
      JacPDstandardNth13gt22 = PDstandardNth13gt22;
      
      JacPDstandardNth13gt23 = PDstandardNth13gt23;
      
      JacPDstandardNth13gt33 = PDstandardNth13gt33;
      
      JacPDstandardNth21gt11 = PDstandardNth12gt11;
      
      JacPDstandardNth21gt12 = PDstandardNth12gt12;
      
      JacPDstandardNth21gt13 = PDstandardNth12gt13;
      
      JacPDstandardNth21gt22 = PDstandardNth12gt22;
      
      JacPDstandardNth21gt23 = PDstandardNth12gt23;
      
      JacPDstandardNth21gt33 = PDstandardNth12gt33;
      
      JacPDstandardNth23gt11 = PDstandardNth23gt11;
      
      JacPDstandardNth23gt12 = PDstandardNth23gt12;
      
      JacPDstandardNth23gt13 = PDstandardNth23gt13;
      
      JacPDstandardNth23gt22 = PDstandardNth23gt22;
      
      JacPDstandardNth23gt23 = PDstandardNth23gt23;
      
      JacPDstandardNth23gt33 = PDstandardNth23gt33;
      
      JacPDstandardNth31gt11 = PDstandardNth13gt11;
      
      JacPDstandardNth31gt12 = PDstandardNth13gt12;
      
      JacPDstandardNth31gt13 = PDstandardNth13gt13;
      
      JacPDstandardNth31gt22 = PDstandardNth13gt22;
      
      JacPDstandardNth31gt23 = PDstandardNth13gt23;
      
      JacPDstandardNth31gt33 = PDstandardNth13gt33;
      
      JacPDstandardNth32gt11 = PDstandardNth23gt11;
      
      JacPDstandardNth32gt12 = PDstandardNth23gt12;
      
      JacPDstandardNth32gt13 = PDstandardNth23gt13;
      
      JacPDstandardNth32gt22 = PDstandardNth23gt22;
      
      JacPDstandardNth32gt23 = PDstandardNth23gt23;
      
      JacPDstandardNth32gt33 = PDstandardNth23gt33;
    }
    
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu21 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu31 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu32 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL dgtu111 = -2*(gtu11*(gtu21*JacPDstandardNth1gt12 + 
      gtu31*JacPDstandardNth1gt13) + gtu21*gtu31*JacPDstandardNth1gt23) - 
      JacPDstandardNth1gt11*SQR(gtu11) - JacPDstandardNth1gt22*SQR(gtu21) - 
      JacPDstandardNth1gt33*SQR(gtu31);
    
    CCTK_REAL dgtu211 = -(gtu11*(gtu21*JacPDstandardNth1gt11 + 
      gtu22*JacPDstandardNth1gt12 + gtu32*JacPDstandardNth1gt13)) - 
      gtu21*(gtu31*JacPDstandardNth1gt13 + gtu22*JacPDstandardNth1gt22 + 
      gtu32*JacPDstandardNth1gt23) - gtu31*(gtu22*JacPDstandardNth1gt23 + 
      gtu32*JacPDstandardNth1gt33) - JacPDstandardNth1gt12*SQR(gtu21);
    
    CCTK_REAL dgtu311 = -(gtu11*(gtu31*JacPDstandardNth1gt11 + 
      gtu32*JacPDstandardNth1gt12 + gtu33*JacPDstandardNth1gt13)) - 
      gtu21*(gtu31*JacPDstandardNth1gt12 + gtu32*JacPDstandardNth1gt22 + 
      gtu33*JacPDstandardNth1gt23) - gtu31*(gtu31*JacPDstandardNth1gt13 + 
      gtu32*JacPDstandardNth1gt23 + gtu33*JacPDstandardNth1gt33);
    
    CCTK_REAL dgtu221 = -2*(gtu21*(gtu22*JacPDstandardNth1gt12 + 
      gtu32*JacPDstandardNth1gt13) + gtu22*gtu32*JacPDstandardNth1gt23) - 
      JacPDstandardNth1gt11*SQR(gtu21) - JacPDstandardNth1gt22*SQR(gtu22) - 
      JacPDstandardNth1gt33*SQR(gtu32);
    
    CCTK_REAL dgtu321 = -(gtu21*(gtu31*JacPDstandardNth1gt11 + 
      gtu32*JacPDstandardNth1gt12 + gtu33*JacPDstandardNth1gt13)) - 
      gtu22*(gtu31*JacPDstandardNth1gt12 + gtu32*JacPDstandardNth1gt22 + 
      gtu33*JacPDstandardNth1gt23) - gtu32*(gtu31*JacPDstandardNth1gt13 + 
      gtu32*JacPDstandardNth1gt23 + gtu33*JacPDstandardNth1gt33);
    
    CCTK_REAL dgtu331 = -2*(gtu31*(gtu32*JacPDstandardNth1gt12 + 
      gtu33*JacPDstandardNth1gt13) + gtu32*gtu33*JacPDstandardNth1gt23) - 
      JacPDstandardNth1gt11*SQR(gtu31) - JacPDstandardNth1gt22*SQR(gtu32) - 
      JacPDstandardNth1gt33*SQR(gtu33);
    
    CCTK_REAL dgtu112 = -2*(gtu11*(gtu21*JacPDstandardNth2gt12 + 
      gtu31*JacPDstandardNth2gt13) + gtu21*gtu31*JacPDstandardNth2gt23) - 
      JacPDstandardNth2gt11*SQR(gtu11) - JacPDstandardNth2gt22*SQR(gtu21) - 
      JacPDstandardNth2gt33*SQR(gtu31);
    
    CCTK_REAL dgtu212 = -(gtu11*(gtu21*JacPDstandardNth2gt11 + 
      gtu22*JacPDstandardNth2gt12 + gtu32*JacPDstandardNth2gt13)) - 
      gtu21*(gtu31*JacPDstandardNth2gt13 + gtu22*JacPDstandardNth2gt22 + 
      gtu32*JacPDstandardNth2gt23) - gtu31*(gtu22*JacPDstandardNth2gt23 + 
      gtu32*JacPDstandardNth2gt33) - JacPDstandardNth2gt12*SQR(gtu21);
    
    CCTK_REAL dgtu312 = -(gtu11*(gtu31*JacPDstandardNth2gt11 + 
      gtu32*JacPDstandardNth2gt12 + gtu33*JacPDstandardNth2gt13)) - 
      gtu21*(gtu31*JacPDstandardNth2gt12 + gtu32*JacPDstandardNth2gt22 + 
      gtu33*JacPDstandardNth2gt23) - gtu31*(gtu31*JacPDstandardNth2gt13 + 
      gtu32*JacPDstandardNth2gt23 + gtu33*JacPDstandardNth2gt33);
    
    CCTK_REAL dgtu222 = -2*(gtu21*(gtu22*JacPDstandardNth2gt12 + 
      gtu32*JacPDstandardNth2gt13) + gtu22*gtu32*JacPDstandardNth2gt23) - 
      JacPDstandardNth2gt11*SQR(gtu21) - JacPDstandardNth2gt22*SQR(gtu22) - 
      JacPDstandardNth2gt33*SQR(gtu32);
    
    CCTK_REAL dgtu322 = -(gtu21*(gtu31*JacPDstandardNth2gt11 + 
      gtu32*JacPDstandardNth2gt12 + gtu33*JacPDstandardNth2gt13)) - 
      gtu22*(gtu31*JacPDstandardNth2gt12 + gtu32*JacPDstandardNth2gt22 + 
      gtu33*JacPDstandardNth2gt23) - gtu32*(gtu31*JacPDstandardNth2gt13 + 
      gtu32*JacPDstandardNth2gt23 + gtu33*JacPDstandardNth2gt33);
    
    CCTK_REAL dgtu332 = -2*(gtu31*(gtu32*JacPDstandardNth2gt12 + 
      gtu33*JacPDstandardNth2gt13) + gtu32*gtu33*JacPDstandardNth2gt23) - 
      JacPDstandardNth2gt11*SQR(gtu31) - JacPDstandardNth2gt22*SQR(gtu32) - 
      JacPDstandardNth2gt33*SQR(gtu33);
    
    CCTK_REAL dgtu113 = -2*(gtu11*(gtu21*JacPDstandardNth3gt12 + 
      gtu31*JacPDstandardNth3gt13) + gtu21*gtu31*JacPDstandardNth3gt23) - 
      JacPDstandardNth3gt11*SQR(gtu11) - JacPDstandardNth3gt22*SQR(gtu21) - 
      JacPDstandardNth3gt33*SQR(gtu31);
    
    CCTK_REAL dgtu213 = -(gtu11*(gtu21*JacPDstandardNth3gt11 + 
      gtu22*JacPDstandardNth3gt12 + gtu32*JacPDstandardNth3gt13)) - 
      gtu21*(gtu31*JacPDstandardNth3gt13 + gtu22*JacPDstandardNth3gt22 + 
      gtu32*JacPDstandardNth3gt23) - gtu31*(gtu22*JacPDstandardNth3gt23 + 
      gtu32*JacPDstandardNth3gt33) - JacPDstandardNth3gt12*SQR(gtu21);
    
    CCTK_REAL dgtu313 = -(gtu11*(gtu31*JacPDstandardNth3gt11 + 
      gtu32*JacPDstandardNth3gt12 + gtu33*JacPDstandardNth3gt13)) - 
      gtu21*(gtu31*JacPDstandardNth3gt12 + gtu32*JacPDstandardNth3gt22 + 
      gtu33*JacPDstandardNth3gt23) - gtu31*(gtu31*JacPDstandardNth3gt13 + 
      gtu32*JacPDstandardNth3gt23 + gtu33*JacPDstandardNth3gt33);
    
    CCTK_REAL dgtu223 = -2*(gtu21*(gtu22*JacPDstandardNth3gt12 + 
      gtu32*JacPDstandardNth3gt13) + gtu22*gtu32*JacPDstandardNth3gt23) - 
      JacPDstandardNth3gt11*SQR(gtu21) - JacPDstandardNth3gt22*SQR(gtu22) - 
      JacPDstandardNth3gt33*SQR(gtu32);
    
    CCTK_REAL dgtu323 = -(gtu21*(gtu31*JacPDstandardNth3gt11 + 
      gtu32*JacPDstandardNth3gt12 + gtu33*JacPDstandardNth3gt13)) - 
      gtu22*(gtu31*JacPDstandardNth3gt12 + gtu32*JacPDstandardNth3gt22 + 
      gtu33*JacPDstandardNth3gt23) - gtu32*(gtu31*JacPDstandardNth3gt13 + 
      gtu32*JacPDstandardNth3gt23 + gtu33*JacPDstandardNth3gt33);
    
    CCTK_REAL dgtu333 = -2*(gtu31*(gtu32*JacPDstandardNth3gt12 + 
      gtu33*JacPDstandardNth3gt13) + gtu32*gtu33*JacPDstandardNth3gt23) - 
      JacPDstandardNth3gt11*SQR(gtu31) - JacPDstandardNth3gt22*SQR(gtu32) - 
      JacPDstandardNth3gt33*SQR(gtu33);
    
    CCTK_REAL Gtl111 = 0.5*JacPDstandardNth1gt11;
    
    CCTK_REAL Gtl112 = 0.5*JacPDstandardNth2gt11;
    
    CCTK_REAL Gtl113 = 0.5*JacPDstandardNth3gt11;
    
    CCTK_REAL Gtl122 = -0.5*JacPDstandardNth1gt22 + JacPDstandardNth2gt12;
    
    CCTK_REAL Gtl123 = 0.5*(-JacPDstandardNth1gt23 + JacPDstandardNth2gt13 
      + JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl133 = -0.5*JacPDstandardNth1gt33 + JacPDstandardNth3gt13;
    
    CCTK_REAL Gtl211 = JacPDstandardNth1gt12 - 0.5*JacPDstandardNth2gt11;
    
    CCTK_REAL Gtl212 = 0.5*JacPDstandardNth1gt22;
    
    CCTK_REAL Gtl213 = 0.5*(JacPDstandardNth1gt23 - JacPDstandardNth2gt13 
      + JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl222 = 0.5*JacPDstandardNth2gt22;
    
    CCTK_REAL Gtl223 = 0.5*JacPDstandardNth3gt22;
    
    CCTK_REAL Gtl233 = -0.5*JacPDstandardNth2gt33 + JacPDstandardNth3gt23;
    
    CCTK_REAL Gtl311 = JacPDstandardNth1gt13 - 0.5*JacPDstandardNth3gt11;
    
    CCTK_REAL Gtl312 = 0.5*(JacPDstandardNth1gt23 + JacPDstandardNth2gt13 
      - JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl313 = 0.5*JacPDstandardNth1gt33;
    
    CCTK_REAL Gtl322 = JacPDstandardNth2gt23 - 0.5*JacPDstandardNth3gt22;
    
    CCTK_REAL Gtl323 = 0.5*JacPDstandardNth2gt33;
    
    CCTK_REAL Gtl333 = 0.5*JacPDstandardNth3gt33;
    
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
      gtu11*JacPDstandardNth11gt11 + 2*(Gt211*Gtlu211 + Gt212*Gtlu212 + 
      Gt213*Gtlu213 + Gt311*Gtlu311 + Gt312*Gtlu312 + Gt313*Gtlu313 + 
      gt11L*JacPDstandardNth1Xt1) + 2*gt12L*JacPDstandardNth1Xt2 + 
      2*gt13L*JacPDstandardNth1Xt3 + gtu21*(-JacPDstandardNth12gt11 - 
      JacPDstandardNth21gt11) - gtu22*JacPDstandardNth22gt11 + 
      gtu31*(-JacPDstandardNth13gt11 - JacPDstandardNth31gt11) + 
      gtu32*(-JacPDstandardNth23gt11 - JacPDstandardNth32gt11) - 
      gtu33*JacPDstandardNth33gt11 + 2*Gtl111*Xtn1 + 2*Gtl112*Xtn2 + 
      2*Gtl113*Xtn3);
    
    CCTK_REAL Rt12 = 0.5*(4*(Gt211*Gtlu221 + Gt212*Gtlu222 + 
      Gt213*Gtlu223) + 2*(Gt122*Gtlu112 + Gt123*Gtlu113 + Gt111*Gtlu121 + 
      Gt212*Gtlu121 + Gt222*Gtlu122 + Gt113*Gtlu123 + Gt223*Gtlu123 + 
      Gt312*Gtlu131 + Gt322*Gtlu132 + Gt323*Gtlu133 + Gt111*Gtlu211 + 
      Gt112*(Gtlu111 + Gtlu122 + Gtlu212) + Gt113*Gtlu213 + Gt311*Gtlu231 + 
      Gt312*Gtlu232 + Gt313*Gtlu233 + Gt311*Gtlu321 + Gt312*Gtlu322 + 
      Gt313*Gtlu323) - gtu11*JacPDstandardNth11gt12 + 
      gt12L*JacPDstandardNth1Xt1 + gt22L*JacPDstandardNth1Xt2 + 
      gt23L*JacPDstandardNth1Xt3 + gtu21*(-JacPDstandardNth12gt12 - 
      JacPDstandardNth21gt12) - gtu22*JacPDstandardNth22gt12 + 
      gt11L*JacPDstandardNth2Xt1 + gt12L*JacPDstandardNth2Xt2 + 
      gt13L*JacPDstandardNth2Xt3 + gtu31*(-JacPDstandardNth13gt12 - 
      JacPDstandardNth31gt12) + gtu32*(-JacPDstandardNth23gt12 - 
      JacPDstandardNth32gt12) - gtu33*JacPDstandardNth33gt12 + Gtl112*Xtn1 + 
      Gtl211*Xtn1 + Gtl122*Xtn2 + Gtl212*Xtn2 + Gtl123*Xtn3 + Gtl213*Xtn3);
    
    CCTK_REAL Rt13 = 0.5*(2*(Gt123*Gtlu112 + Gt133*Gtlu113 + Gt213*Gtlu121 
      + Gt223*Gtlu122 + Gt233*Gtlu123 + Gt111*Gtlu131 + Gt313*Gtlu131 + 
      Gt112*Gtlu132 + Gt323*Gtlu132 + Gt333*Gtlu133 + Gt211*Gtlu231 + 
      Gt212*Gtlu232 + Gt213*Gtlu233 + Gt111*Gtlu311 + Gt112*Gtlu312 + 
      Gt113*(Gtlu111 + Gtlu133 + Gtlu313) + Gt211*Gtlu321 + Gt212*Gtlu322 + 
      Gt213*Gtlu323) + 4*(Gt311*Gtlu331 + Gt312*Gtlu332 + Gt313*Gtlu333) - 
      gtu11*JacPDstandardNth11gt13 + gt13L*JacPDstandardNth1Xt1 + 
      gt23L*JacPDstandardNth1Xt2 + gt33L*JacPDstandardNth1Xt3 + 
      gtu21*(-JacPDstandardNth12gt13 - JacPDstandardNth21gt13) - 
      gtu22*JacPDstandardNth22gt13 + gtu31*(-JacPDstandardNth13gt13 - 
      JacPDstandardNth31gt13) + gtu32*(-JacPDstandardNth23gt13 - 
      JacPDstandardNth32gt13) - gtu33*JacPDstandardNth33gt13 + 
      gt11L*JacPDstandardNth3Xt1 + gt12L*JacPDstandardNth3Xt2 + 
      gt13L*JacPDstandardNth3Xt3 + Gtl113*Xtn1 + Gtl311*Xtn1 + Gtl123*Xtn2 
      + Gtl312*Xtn2 + Gtl133*Xtn3 + Gtl313*Xtn3);
    
    CCTK_REAL Rt22 = 0.5*(6*(Gt212*Gtlu221 + Gt222*Gtlu222 + 
      Gt223*Gtlu223) + 4*(Gt123*Gtlu213 + Gt312*Gtlu231 + Gt322*Gtlu232 + 
      Gt323*Gtlu233) - gtu11*JacPDstandardNth11gt22 + 
      gtu21*(-JacPDstandardNth12gt22 - JacPDstandardNth21gt22) - 
      gtu22*JacPDstandardNth22gt22 + 2*(Gt123*Gtlu123 + Gt112*(Gtlu121 + 
      2*Gtlu211) + Gt122*(Gtlu122 + 2*Gtlu212) + Gt312*Gtlu321 + 
      Gt322*Gtlu322 + Gt323*Gtlu323 + gt12L*JacPDstandardNth2Xt1) + 
      2*gt22L*JacPDstandardNth2Xt2 + 2*gt23L*JacPDstandardNth2Xt3 + 
      gtu31*(-JacPDstandardNth13gt22 - JacPDstandardNth31gt22) + 
      gtu32*(-JacPDstandardNth23gt22 - JacPDstandardNth32gt22) - 
      gtu33*JacPDstandardNth33gt22 + 2*Gtl212*Xtn1 + 2*Gtl222*Xtn2 + 
      2*Gtl223*Xtn3);
    
    CCTK_REAL Rt23 = 0.5*(2*(Gt123*Gtlu133 + Gt113*Gtlu211 + Gt123*Gtlu212 
      + Gt133*Gtlu213 + Gt213*Gtlu221 + Gt223*Gtlu222 + Gt233*Gtlu223 + 
      Gt212*Gtlu231 + Gt313*Gtlu231 + Gt222*Gtlu232 + Gt323*Gtlu232 + 
      Gt223*Gtlu233 + Gt333*Gtlu233 + Gt112*(Gtlu131 + Gtlu311) + 
      Gt122*(Gtlu132 + Gtlu312) + Gt123*Gtlu313 + Gt212*Gtlu321 + 
      Gt222*Gtlu322 + Gt223*Gtlu323) + 4*(Gt312*Gtlu331 + Gt322*Gtlu332 + 
      Gt323*Gtlu333) - gtu11*JacPDstandardNth11gt23 + 
      gtu21*(-JacPDstandardNth12gt23 - JacPDstandardNth21gt23) - 
      gtu22*JacPDstandardNth22gt23 + gt13L*JacPDstandardNth2Xt1 + 
      gt23L*JacPDstandardNth2Xt2 + gt33L*JacPDstandardNth2Xt3 + 
      gtu31*(-JacPDstandardNth13gt23 - JacPDstandardNth31gt23) + 
      gtu32*(-JacPDstandardNth23gt23 - JacPDstandardNth32gt23) - 
      gtu33*JacPDstandardNth33gt23 + gt12L*JacPDstandardNth3Xt1 + 
      gt22L*JacPDstandardNth3Xt2 + gt23L*JacPDstandardNth3Xt3 + 
      Gtl213*Xtn1 + Gtl312*Xtn1 + Gtl223*Xtn2 + Gtl322*Xtn2 + Gtl233*Xtn3 + 
      Gtl323*Xtn3);
    
    CCTK_REAL Rt33 = 0.5*(4*(Gt133*Gtlu313 + Gt213*Gtlu321 + Gt223*Gtlu322 
      + Gt233*Gtlu323) + 6*(Gt313*Gtlu331 + Gt323*Gtlu332 + Gt333*Gtlu333) - 
      gtu11*JacPDstandardNth11gt33 + gtu21*(-JacPDstandardNth12gt33 - 
      JacPDstandardNth21gt33) - gtu22*JacPDstandardNth22gt33 + 
      gtu31*(-JacPDstandardNth13gt33 - JacPDstandardNth31gt33) + 
      gtu32*(-JacPDstandardNth23gt33 - JacPDstandardNth32gt33) - 
      gtu33*JacPDstandardNth33gt33 + 2*(Gt133*Gtlu133 + Gt213*Gtlu231 + 
      Gt223*Gtlu232 + Gt233*Gtlu233 + Gt113*(Gtlu131 + 2*Gtlu311) + 
      Gt123*(Gtlu132 + 2*Gtlu312) + gt13L*JacPDstandardNth3Xt1) + 
      2*gt23L*JacPDstandardNth3Xt2 + 2*gt33L*JacPDstandardNth3Xt3 + 
      2*Gtl313*Xtn1 + 2*Gtl323*Xtn2 + 2*Gtl333*Xtn3);
    
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
      0.0208333333333333333333333333333*(At11L*(3*zL*dgtu112 - 
      3*yL*dgtu113) + At22L*(3*zL*dgtu222 - 3*yL*dgtu223) + 6*(Atm23 
      + zL*(At12L*dgtu212 + At13L*dgtu312 + At23L*dgtu322)) - 
      6*(Atm32 + yL*(At12L*dgtu213 + At13L*dgtu313 + At23L*dgtu323)) 
      + At33L*(3*zL*dgtu332 - 3*yL*dgtu333) + 
      zL*(-4*JacPDstandardNth2trK - 48*Pi*S2) + 
      yL*(4*JacPDstandardNth3trK + 48*Pi*S3))*INV(Pi)*pow(ephi,6);
    
    CCTK_REAL Jadm2L = 
      0.0208333333333333333333333333333*(At11L*(-3*zL*dgtu111 + 
      3*xL*dgtu113) + At22L*(-3*zL*dgtu221 + 3*xL*dgtu223) - 6*(Atm13 
      + zL*(At12L*dgtu211 + At13L*dgtu311 + At23L*dgtu321)) + 
      6*(Atm31 + xL*(At12L*dgtu213 + At13L*dgtu313 + At23L*dgtu323)) 
      + At33L*(-3*zL*dgtu331 + 3*xL*dgtu333) + 
      zL*(4*JacPDstandardNth1trK + 48*Pi*S1) + 
      xL*(-4*JacPDstandardNth3trK - 48*Pi*S3))*INV(Pi)*pow(ephi,6);
    
    CCTK_REAL Jadm3L = 
      0.0208333333333333333333333333333*(At11L*(3*yL*dgtu111 - 
      3*xL*dgtu112) + At22L*(3*yL*dgtu221 - 3*xL*dgtu222) + 6*(Atm12 
      + yL*(At12L*dgtu211 + At13L*dgtu311 + At23L*dgtu321)) - 
      6*(Atm21 + xL*(At12L*dgtu212 + At13L*dgtu312 + At23L*dgtu322)) 
      + At33L*(3*yL*dgtu331 - 3*xL*dgtu332) + 
      yL*(-4*JacPDstandardNth1trK - 48*Pi*S1) + 
      xL*(4*JacPDstandardNth2trK + 48*Pi*S2))*INV(Pi)*pow(ephi,6);
    
    /* Copy local copies back to grid functions */
    Jadm1[index] = Jadm1L;
    Jadm2[index] = Jadm2L;
    Jadm3[index] = Jadm3L;
    Madm[index] = MadmL;
  }
  CCTK_ENDLOOP3(ML_ADMQuantities_MP);
}

extern "C" void ML_ADMQuantities_MP(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMQuantities_MP_Body");
  }
  
  if (cctk_iteration % ML_ADMQuantities_MP_calc_every != ML_ADMQuantities_MP_calc_offset)
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
    "ML_ADMQuantities_MP::ML_Jadm",
    "ML_ADMQuantities_MP::ML_Madm"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADMQuantities_MP", 11, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_ADMQuantities_MP", 2, 2, 2);
  
  GenericFD_LoopOverInterior(cctkGH, ML_ADMQuantities_MP_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_ADMQuantities_MP_Body");
  }
}
