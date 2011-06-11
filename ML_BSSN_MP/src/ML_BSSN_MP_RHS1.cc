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

extern "C" void ML_BSSN_MP_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_trace_curvrhs.");
  return;
}

static void ML_BSSN_MP_RHS1_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_RHS1_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_RHS1_calc_every != ML_BSSN_MP_RHS1_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"grid::coordinates","Grid::coordinates","ML_BSSN_MP::ML_curv","ML_BSSN_MP::ML_dtlapse","ML_BSSN_MP::ML_dtlapserhs","ML_BSSN_MP::ML_dtshift","ML_BSSN_MP::ML_dtshiftrhs","ML_BSSN_MP::ML_Gamma","ML_BSSN_MP::ML_Gammarhs","ML_BSSN_MP::ML_lapse","ML_BSSN_MP::ML_lapserhs","ML_BSSN_MP::ML_log_confac","ML_BSSN_MP::ML_log_confacrhs","ML_BSSN_MP::ML_metric","ML_BSSN_MP::ML_metricrhs","ML_BSSN_MP::ML_shift","ML_BSSN_MP::ML_shiftrhs","ML_BSSN_MP::ML_trace_curv","ML_BSSN_MP::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_RHS1", 19, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_MP_RHS1", 2, 2, 2);
  
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
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_MP_RHS1,
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
    CCTK_REAL JacPDstandardNth11alpha;
    CCTK_REAL JacPDstandardNth11beta1;
    CCTK_REAL JacPDstandardNth11beta2;
    CCTK_REAL JacPDstandardNth11beta3;
    CCTK_REAL JacPDstandardNth12alpha;
    CCTK_REAL JacPDstandardNth12beta1;
    CCTK_REAL JacPDstandardNth12beta2;
    CCTK_REAL JacPDstandardNth12beta3;
    CCTK_REAL JacPDstandardNth13alpha;
    CCTK_REAL JacPDstandardNth13beta1;
    CCTK_REAL JacPDstandardNth13beta2;
    CCTK_REAL JacPDstandardNth13beta3;
    CCTK_REAL JacPDstandardNth1alpha;
    CCTK_REAL JacPDstandardNth1beta1;
    CCTK_REAL JacPDstandardNth1beta2;
    CCTK_REAL JacPDstandardNth1beta3;
    CCTK_REAL JacPDstandardNth1gt11;
    CCTK_REAL JacPDstandardNth1gt12;
    CCTK_REAL JacPDstandardNth1gt13;
    CCTK_REAL JacPDstandardNth1gt22;
    CCTK_REAL JacPDstandardNth1gt23;
    CCTK_REAL JacPDstandardNth1gt33;
    CCTK_REAL JacPDstandardNth1phi;
    CCTK_REAL JacPDstandardNth1trK;
    CCTK_REAL JacPDstandardNth21alpha;
    CCTK_REAL JacPDstandardNth21beta1;
    CCTK_REAL JacPDstandardNth21beta2;
    CCTK_REAL JacPDstandardNth21beta3;
    CCTK_REAL JacPDstandardNth22alpha;
    CCTK_REAL JacPDstandardNth22beta1;
    CCTK_REAL JacPDstandardNth22beta2;
    CCTK_REAL JacPDstandardNth22beta3;
    CCTK_REAL JacPDstandardNth23alpha;
    CCTK_REAL JacPDstandardNth23beta1;
    CCTK_REAL JacPDstandardNth23beta2;
    CCTK_REAL JacPDstandardNth23beta3;
    CCTK_REAL JacPDstandardNth2alpha;
    CCTK_REAL JacPDstandardNth2beta1;
    CCTK_REAL JacPDstandardNth2beta2;
    CCTK_REAL JacPDstandardNth2beta3;
    CCTK_REAL JacPDstandardNth2gt11;
    CCTK_REAL JacPDstandardNth2gt12;
    CCTK_REAL JacPDstandardNth2gt13;
    CCTK_REAL JacPDstandardNth2gt22;
    CCTK_REAL JacPDstandardNth2gt23;
    CCTK_REAL JacPDstandardNth2gt33;
    CCTK_REAL JacPDstandardNth2phi;
    CCTK_REAL JacPDstandardNth2trK;
    CCTK_REAL JacPDstandardNth31alpha;
    CCTK_REAL JacPDstandardNth31beta1;
    CCTK_REAL JacPDstandardNth31beta2;
    CCTK_REAL JacPDstandardNth31beta3;
    CCTK_REAL JacPDstandardNth32alpha;
    CCTK_REAL JacPDstandardNth32beta1;
    CCTK_REAL JacPDstandardNth32beta2;
    CCTK_REAL JacPDstandardNth32beta3;
    CCTK_REAL JacPDstandardNth33alpha;
    CCTK_REAL JacPDstandardNth33beta1;
    CCTK_REAL JacPDstandardNth33beta2;
    CCTK_REAL JacPDstandardNth33beta3;
    CCTK_REAL JacPDstandardNth3alpha;
    CCTK_REAL JacPDstandardNth3beta1;
    CCTK_REAL JacPDstandardNth3beta2;
    CCTK_REAL JacPDstandardNth3beta3;
    CCTK_REAL JacPDstandardNth3gt11;
    CCTK_REAL JacPDstandardNth3gt12;
    CCTK_REAL JacPDstandardNth3gt13;
    CCTK_REAL JacPDstandardNth3gt22;
    CCTK_REAL JacPDstandardNth3gt23;
    CCTK_REAL JacPDstandardNth3gt33;
    CCTK_REAL JacPDstandardNth3phi;
    CCTK_REAL JacPDstandardNth3trK;
    
    if (use_jacobian)
    {
      JacPDstandardNth1alpha = J11L*PDstandardNth1alpha + 
        J21L*PDstandardNth2alpha + J31L*PDstandardNth3alpha;
      
      JacPDstandardNth1beta1 = J11L*PDstandardNth1beta1 + 
        J21L*PDstandardNth2beta1 + J31L*PDstandardNth3beta1;
      
      JacPDstandardNth1beta2 = J11L*PDstandardNth1beta2 + 
        J21L*PDstandardNth2beta2 + J31L*PDstandardNth3beta2;
      
      JacPDstandardNth1beta3 = J11L*PDstandardNth1beta3 + 
        J21L*PDstandardNth2beta3 + J31L*PDstandardNth3beta3;
      
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
      
      JacPDstandardNth1phi = J11L*PDstandardNth1phi + J21L*PDstandardNth2phi 
        + J31L*PDstandardNth3phi;
      
      JacPDstandardNth1trK = J11L*PDstandardNth1trK + J21L*PDstandardNth2trK 
        + J31L*PDstandardNth3trK;
      
      JacPDstandardNth2alpha = J12L*PDstandardNth1alpha + 
        J22L*PDstandardNth2alpha + J32L*PDstandardNth3alpha;
      
      JacPDstandardNth2beta1 = J12L*PDstandardNth1beta1 + 
        J22L*PDstandardNth2beta1 + J32L*PDstandardNth3beta1;
      
      JacPDstandardNth2beta2 = J12L*PDstandardNth1beta2 + 
        J22L*PDstandardNth2beta2 + J32L*PDstandardNth3beta2;
      
      JacPDstandardNth2beta3 = J12L*PDstandardNth1beta3 + 
        J22L*PDstandardNth2beta3 + J32L*PDstandardNth3beta3;
      
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
      
      JacPDstandardNth2phi = J12L*PDstandardNth1phi + J22L*PDstandardNth2phi 
        + J32L*PDstandardNth3phi;
      
      JacPDstandardNth2trK = J12L*PDstandardNth1trK + J22L*PDstandardNth2trK 
        + J32L*PDstandardNth3trK;
      
      JacPDstandardNth3alpha = J13L*PDstandardNth1alpha + 
        J23L*PDstandardNth2alpha + J33L*PDstandardNth3alpha;
      
      JacPDstandardNth3beta1 = J13L*PDstandardNth1beta1 + 
        J23L*PDstandardNth2beta1 + J33L*PDstandardNth3beta1;
      
      JacPDstandardNth3beta2 = J13L*PDstandardNth1beta2 + 
        J23L*PDstandardNth2beta2 + J33L*PDstandardNth3beta2;
      
      JacPDstandardNth3beta3 = J13L*PDstandardNth1beta3 + 
        J23L*PDstandardNth2beta3 + J33L*PDstandardNth3beta3;
      
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
      
      JacPDstandardNth3phi = J13L*PDstandardNth1phi + J23L*PDstandardNth2phi 
        + J33L*PDstandardNth3phi;
      
      JacPDstandardNth3trK = J13L*PDstandardNth1trK + J23L*PDstandardNth2trK 
        + J33L*PDstandardNth3trK;
      
      JacPDstandardNth11alpha = dJ111L*PDstandardNth1alpha + 
        2*(J11L*(J21L*PDstandardNth12alpha + J31L*PDstandardNth13alpha) + 
        J21L*J31L*PDstandardNth23alpha) + dJ211L*PDstandardNth2alpha + 
        dJ311L*PDstandardNth3alpha + PDstandardNth11alpha*SQR(J11L) + 
        PDstandardNth22alpha*SQR(J21L) + PDstandardNth33alpha*SQR(J31L);
      
      JacPDstandardNth11beta1 = dJ111L*PDstandardNth1beta1 + 
        2*(J11L*(J21L*PDstandardNth12beta1 + J31L*PDstandardNth13beta1) + 
        J21L*J31L*PDstandardNth23beta1) + dJ211L*PDstandardNth2beta1 + 
        dJ311L*PDstandardNth3beta1 + PDstandardNth11beta1*SQR(J11L) + 
        PDstandardNth22beta1*SQR(J21L) + PDstandardNth33beta1*SQR(J31L);
      
      JacPDstandardNth11beta2 = dJ111L*PDstandardNth1beta2 + 
        2*(J11L*(J21L*PDstandardNth12beta2 + J31L*PDstandardNth13beta2) + 
        J21L*J31L*PDstandardNth23beta2) + dJ211L*PDstandardNth2beta2 + 
        dJ311L*PDstandardNth3beta2 + PDstandardNth11beta2*SQR(J11L) + 
        PDstandardNth22beta2*SQR(J21L) + PDstandardNth33beta2*SQR(J31L);
      
      JacPDstandardNth11beta3 = dJ111L*PDstandardNth1beta3 + 
        2*(J11L*(J21L*PDstandardNth12beta3 + J31L*PDstandardNth13beta3) + 
        J21L*J31L*PDstandardNth23beta3) + dJ211L*PDstandardNth2beta3 + 
        dJ311L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J11L) + 
        PDstandardNth22beta3*SQR(J21L) + PDstandardNth33beta3*SQR(J31L);
      
      JacPDstandardNth22alpha = dJ122L*PDstandardNth1alpha + 
        2*(J12L*(J22L*PDstandardNth12alpha + J32L*PDstandardNth13alpha) + 
        J22L*J32L*PDstandardNth23alpha) + dJ222L*PDstandardNth2alpha + 
        dJ322L*PDstandardNth3alpha + PDstandardNth11alpha*SQR(J12L) + 
        PDstandardNth22alpha*SQR(J22L) + PDstandardNth33alpha*SQR(J32L);
      
      JacPDstandardNth22beta1 = dJ122L*PDstandardNth1beta1 + 
        2*(J12L*(J22L*PDstandardNth12beta1 + J32L*PDstandardNth13beta1) + 
        J22L*J32L*PDstandardNth23beta1) + dJ222L*PDstandardNth2beta1 + 
        dJ322L*PDstandardNth3beta1 + PDstandardNth11beta1*SQR(J12L) + 
        PDstandardNth22beta1*SQR(J22L) + PDstandardNth33beta1*SQR(J32L);
      
      JacPDstandardNth22beta2 = dJ122L*PDstandardNth1beta2 + 
        2*(J12L*(J22L*PDstandardNth12beta2 + J32L*PDstandardNth13beta2) + 
        J22L*J32L*PDstandardNth23beta2) + dJ222L*PDstandardNth2beta2 + 
        dJ322L*PDstandardNth3beta2 + PDstandardNth11beta2*SQR(J12L) + 
        PDstandardNth22beta2*SQR(J22L) + PDstandardNth33beta2*SQR(J32L);
      
      JacPDstandardNth22beta3 = dJ122L*PDstandardNth1beta3 + 
        2*(J12L*(J22L*PDstandardNth12beta3 + J32L*PDstandardNth13beta3) + 
        J22L*J32L*PDstandardNth23beta3) + dJ222L*PDstandardNth2beta3 + 
        dJ322L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J12L) + 
        PDstandardNth22beta3*SQR(J22L) + PDstandardNth33beta3*SQR(J32L);
      
      JacPDstandardNth33alpha = dJ133L*PDstandardNth1alpha + 
        2*(J13L*(J23L*PDstandardNth12alpha + J33L*PDstandardNth13alpha) + 
        J23L*J33L*PDstandardNth23alpha) + dJ233L*PDstandardNth2alpha + 
        dJ333L*PDstandardNth3alpha + PDstandardNth11alpha*SQR(J13L) + 
        PDstandardNth22alpha*SQR(J23L) + PDstandardNth33alpha*SQR(J33L);
      
      JacPDstandardNth33beta1 = dJ133L*PDstandardNth1beta1 + 
        2*(J13L*(J23L*PDstandardNth12beta1 + J33L*PDstandardNth13beta1) + 
        J23L*J33L*PDstandardNth23beta1) + dJ233L*PDstandardNth2beta1 + 
        dJ333L*PDstandardNth3beta1 + PDstandardNth11beta1*SQR(J13L) + 
        PDstandardNth22beta1*SQR(J23L) + PDstandardNth33beta1*SQR(J33L);
      
      JacPDstandardNth33beta2 = dJ133L*PDstandardNth1beta2 + 
        2*(J13L*(J23L*PDstandardNth12beta2 + J33L*PDstandardNth13beta2) + 
        J23L*J33L*PDstandardNth23beta2) + dJ233L*PDstandardNth2beta2 + 
        dJ333L*PDstandardNth3beta2 + PDstandardNth11beta2*SQR(J13L) + 
        PDstandardNth22beta2*SQR(J23L) + PDstandardNth33beta2*SQR(J33L);
      
      JacPDstandardNth33beta3 = dJ133L*PDstandardNth1beta3 + 
        2*(J13L*(J23L*PDstandardNth12beta3 + J33L*PDstandardNth13beta3) + 
        J23L*J33L*PDstandardNth23beta3) + dJ233L*PDstandardNth2beta3 + 
        dJ333L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J13L) + 
        PDstandardNth22beta3*SQR(J23L) + PDstandardNth33beta3*SQR(J33L);
      
      JacPDstandardNth12alpha = J12L*(J11L*PDstandardNth11alpha + 
        J21L*PDstandardNth12alpha + J31L*PDstandardNth13alpha) + 
        J11L*(J22L*PDstandardNth12alpha + J32L*PDstandardNth13alpha) + 
        dJ112L*PDstandardNth1alpha + J22L*(J21L*PDstandardNth22alpha + 
        J31L*PDstandardNth23alpha) + dJ212L*PDstandardNth2alpha + 
        J32L*(J21L*PDstandardNth23alpha + J31L*PDstandardNth33alpha) + 
        dJ312L*PDstandardNth3alpha;
      
      JacPDstandardNth12beta1 = J12L*(J11L*PDstandardNth11beta1 + 
        J21L*PDstandardNth12beta1 + J31L*PDstandardNth13beta1) + 
        J11L*(J22L*PDstandardNth12beta1 + J32L*PDstandardNth13beta1) + 
        dJ112L*PDstandardNth1beta1 + J22L*(J21L*PDstandardNth22beta1 + 
        J31L*PDstandardNth23beta1) + dJ212L*PDstandardNth2beta1 + 
        J32L*(J21L*PDstandardNth23beta1 + J31L*PDstandardNth33beta1) + 
        dJ312L*PDstandardNth3beta1;
      
      JacPDstandardNth12beta2 = J12L*(J11L*PDstandardNth11beta2 + 
        J21L*PDstandardNth12beta2 + J31L*PDstandardNth13beta2) + 
        J11L*(J22L*PDstandardNth12beta2 + J32L*PDstandardNth13beta2) + 
        dJ112L*PDstandardNth1beta2 + J22L*(J21L*PDstandardNth22beta2 + 
        J31L*PDstandardNth23beta2) + dJ212L*PDstandardNth2beta2 + 
        J32L*(J21L*PDstandardNth23beta2 + J31L*PDstandardNth33beta2) + 
        dJ312L*PDstandardNth3beta2;
      
      JacPDstandardNth12beta3 = J12L*(J11L*PDstandardNth11beta3 + 
        J21L*PDstandardNth12beta3 + J31L*PDstandardNth13beta3) + 
        J11L*(J22L*PDstandardNth12beta3 + J32L*PDstandardNth13beta3) + 
        dJ112L*PDstandardNth1beta3 + J22L*(J21L*PDstandardNth22beta3 + 
        J31L*PDstandardNth23beta3) + dJ212L*PDstandardNth2beta3 + 
        J32L*(J21L*PDstandardNth23beta3 + J31L*PDstandardNth33beta3) + 
        dJ312L*PDstandardNth3beta3;
      
      JacPDstandardNth13alpha = J13L*(J11L*PDstandardNth11alpha + 
        J21L*PDstandardNth12alpha + J31L*PDstandardNth13alpha) + 
        J11L*(J23L*PDstandardNth12alpha + J33L*PDstandardNth13alpha) + 
        dJ113L*PDstandardNth1alpha + J23L*(J21L*PDstandardNth22alpha + 
        J31L*PDstandardNth23alpha) + dJ213L*PDstandardNth2alpha + 
        J33L*(J21L*PDstandardNth23alpha + J31L*PDstandardNth33alpha) + 
        dJ313L*PDstandardNth3alpha;
      
      JacPDstandardNth13beta1 = J13L*(J11L*PDstandardNth11beta1 + 
        J21L*PDstandardNth12beta1 + J31L*PDstandardNth13beta1) + 
        J11L*(J23L*PDstandardNth12beta1 + J33L*PDstandardNth13beta1) + 
        dJ113L*PDstandardNth1beta1 + J23L*(J21L*PDstandardNth22beta1 + 
        J31L*PDstandardNth23beta1) + dJ213L*PDstandardNth2beta1 + 
        J33L*(J21L*PDstandardNth23beta1 + J31L*PDstandardNth33beta1) + 
        dJ313L*PDstandardNth3beta1;
      
      JacPDstandardNth13beta2 = J13L*(J11L*PDstandardNth11beta2 + 
        J21L*PDstandardNth12beta2 + J31L*PDstandardNth13beta2) + 
        J11L*(J23L*PDstandardNth12beta2 + J33L*PDstandardNth13beta2) + 
        dJ113L*PDstandardNth1beta2 + J23L*(J21L*PDstandardNth22beta2 + 
        J31L*PDstandardNth23beta2) + dJ213L*PDstandardNth2beta2 + 
        J33L*(J21L*PDstandardNth23beta2 + J31L*PDstandardNth33beta2) + 
        dJ313L*PDstandardNth3beta2;
      
      JacPDstandardNth13beta3 = J13L*(J11L*PDstandardNth11beta3 + 
        J21L*PDstandardNth12beta3 + J31L*PDstandardNth13beta3) + 
        J11L*(J23L*PDstandardNth12beta3 + J33L*PDstandardNth13beta3) + 
        dJ113L*PDstandardNth1beta3 + J23L*(J21L*PDstandardNth22beta3 + 
        J31L*PDstandardNth23beta3) + dJ213L*PDstandardNth2beta3 + 
        J33L*(J21L*PDstandardNth23beta3 + J31L*PDstandardNth33beta3) + 
        dJ313L*PDstandardNth3beta3;
      
      JacPDstandardNth21alpha = J12L*(J11L*PDstandardNth11alpha + 
        J21L*PDstandardNth12alpha + J31L*PDstandardNth13alpha) + 
        J11L*(J22L*PDstandardNth12alpha + J32L*PDstandardNth13alpha) + 
        dJ112L*PDstandardNth1alpha + J22L*(J21L*PDstandardNth22alpha + 
        J31L*PDstandardNth23alpha) + dJ212L*PDstandardNth2alpha + 
        J32L*(J21L*PDstandardNth23alpha + J31L*PDstandardNth33alpha) + 
        dJ312L*PDstandardNth3alpha;
      
      JacPDstandardNth21beta1 = J12L*(J11L*PDstandardNth11beta1 + 
        J21L*PDstandardNth12beta1 + J31L*PDstandardNth13beta1) + 
        J11L*(J22L*PDstandardNth12beta1 + J32L*PDstandardNth13beta1) + 
        dJ112L*PDstandardNth1beta1 + J22L*(J21L*PDstandardNth22beta1 + 
        J31L*PDstandardNth23beta1) + dJ212L*PDstandardNth2beta1 + 
        J32L*(J21L*PDstandardNth23beta1 + J31L*PDstandardNth33beta1) + 
        dJ312L*PDstandardNth3beta1;
      
      JacPDstandardNth21beta2 = J12L*(J11L*PDstandardNth11beta2 + 
        J21L*PDstandardNth12beta2 + J31L*PDstandardNth13beta2) + 
        J11L*(J22L*PDstandardNth12beta2 + J32L*PDstandardNth13beta2) + 
        dJ112L*PDstandardNth1beta2 + J22L*(J21L*PDstandardNth22beta2 + 
        J31L*PDstandardNth23beta2) + dJ212L*PDstandardNth2beta2 + 
        J32L*(J21L*PDstandardNth23beta2 + J31L*PDstandardNth33beta2) + 
        dJ312L*PDstandardNth3beta2;
      
      JacPDstandardNth21beta3 = J12L*(J11L*PDstandardNth11beta3 + 
        J21L*PDstandardNth12beta3 + J31L*PDstandardNth13beta3) + 
        J11L*(J22L*PDstandardNth12beta3 + J32L*PDstandardNth13beta3) + 
        dJ112L*PDstandardNth1beta3 + J22L*(J21L*PDstandardNth22beta3 + 
        J31L*PDstandardNth23beta3) + dJ212L*PDstandardNth2beta3 + 
        J32L*(J21L*PDstandardNth23beta3 + J31L*PDstandardNth33beta3) + 
        dJ312L*PDstandardNth3beta3;
      
      JacPDstandardNth23alpha = J13L*(J12L*PDstandardNth11alpha + 
        J22L*PDstandardNth12alpha + J32L*PDstandardNth13alpha) + 
        J12L*(J23L*PDstandardNth12alpha + J33L*PDstandardNth13alpha) + 
        dJ123L*PDstandardNth1alpha + J23L*(J22L*PDstandardNth22alpha + 
        J32L*PDstandardNth23alpha) + dJ223L*PDstandardNth2alpha + 
        J33L*(J22L*PDstandardNth23alpha + J32L*PDstandardNth33alpha) + 
        dJ323L*PDstandardNth3alpha;
      
      JacPDstandardNth23beta1 = J13L*(J12L*PDstandardNth11beta1 + 
        J22L*PDstandardNth12beta1 + J32L*PDstandardNth13beta1) + 
        J12L*(J23L*PDstandardNth12beta1 + J33L*PDstandardNth13beta1) + 
        dJ123L*PDstandardNth1beta1 + J23L*(J22L*PDstandardNth22beta1 + 
        J32L*PDstandardNth23beta1) + dJ223L*PDstandardNth2beta1 + 
        J33L*(J22L*PDstandardNth23beta1 + J32L*PDstandardNth33beta1) + 
        dJ323L*PDstandardNth3beta1;
      
      JacPDstandardNth23beta2 = J13L*(J12L*PDstandardNth11beta2 + 
        J22L*PDstandardNth12beta2 + J32L*PDstandardNth13beta2) + 
        J12L*(J23L*PDstandardNth12beta2 + J33L*PDstandardNth13beta2) + 
        dJ123L*PDstandardNth1beta2 + J23L*(J22L*PDstandardNth22beta2 + 
        J32L*PDstandardNth23beta2) + dJ223L*PDstandardNth2beta2 + 
        J33L*(J22L*PDstandardNth23beta2 + J32L*PDstandardNth33beta2) + 
        dJ323L*PDstandardNth3beta2;
      
      JacPDstandardNth23beta3 = J13L*(J12L*PDstandardNth11beta3 + 
        J22L*PDstandardNth12beta3 + J32L*PDstandardNth13beta3) + 
        J12L*(J23L*PDstandardNth12beta3 + J33L*PDstandardNth13beta3) + 
        dJ123L*PDstandardNth1beta3 + J23L*(J22L*PDstandardNth22beta3 + 
        J32L*PDstandardNth23beta3) + dJ223L*PDstandardNth2beta3 + 
        J33L*(J22L*PDstandardNth23beta3 + J32L*PDstandardNth33beta3) + 
        dJ323L*PDstandardNth3beta3;
      
      JacPDstandardNth31alpha = J13L*(J11L*PDstandardNth11alpha + 
        J21L*PDstandardNth12alpha + J31L*PDstandardNth13alpha) + 
        J11L*(J23L*PDstandardNth12alpha + J33L*PDstandardNth13alpha) + 
        dJ113L*PDstandardNth1alpha + J23L*(J21L*PDstandardNth22alpha + 
        J31L*PDstandardNth23alpha) + dJ213L*PDstandardNth2alpha + 
        J33L*(J21L*PDstandardNth23alpha + J31L*PDstandardNth33alpha) + 
        dJ313L*PDstandardNth3alpha;
      
      JacPDstandardNth31beta1 = J13L*(J11L*PDstandardNth11beta1 + 
        J21L*PDstandardNth12beta1 + J31L*PDstandardNth13beta1) + 
        J11L*(J23L*PDstandardNth12beta1 + J33L*PDstandardNth13beta1) + 
        dJ113L*PDstandardNth1beta1 + J23L*(J21L*PDstandardNth22beta1 + 
        J31L*PDstandardNth23beta1) + dJ213L*PDstandardNth2beta1 + 
        J33L*(J21L*PDstandardNth23beta1 + J31L*PDstandardNth33beta1) + 
        dJ313L*PDstandardNth3beta1;
      
      JacPDstandardNth31beta2 = J13L*(J11L*PDstandardNth11beta2 + 
        J21L*PDstandardNth12beta2 + J31L*PDstandardNth13beta2) + 
        J11L*(J23L*PDstandardNth12beta2 + J33L*PDstandardNth13beta2) + 
        dJ113L*PDstandardNth1beta2 + J23L*(J21L*PDstandardNth22beta2 + 
        J31L*PDstandardNth23beta2) + dJ213L*PDstandardNth2beta2 + 
        J33L*(J21L*PDstandardNth23beta2 + J31L*PDstandardNth33beta2) + 
        dJ313L*PDstandardNth3beta2;
      
      JacPDstandardNth31beta3 = J13L*(J11L*PDstandardNth11beta3 + 
        J21L*PDstandardNth12beta3 + J31L*PDstandardNth13beta3) + 
        J11L*(J23L*PDstandardNth12beta3 + J33L*PDstandardNth13beta3) + 
        dJ113L*PDstandardNth1beta3 + J23L*(J21L*PDstandardNth22beta3 + 
        J31L*PDstandardNth23beta3) + dJ213L*PDstandardNth2beta3 + 
        J33L*(J21L*PDstandardNth23beta3 + J31L*PDstandardNth33beta3) + 
        dJ313L*PDstandardNth3beta3;
      
      JacPDstandardNth32alpha = J13L*(J12L*PDstandardNth11alpha + 
        J22L*PDstandardNth12alpha + J32L*PDstandardNth13alpha) + 
        J12L*(J23L*PDstandardNth12alpha + J33L*PDstandardNth13alpha) + 
        dJ123L*PDstandardNth1alpha + J23L*(J22L*PDstandardNth22alpha + 
        J32L*PDstandardNth23alpha) + dJ223L*PDstandardNth2alpha + 
        J33L*(J22L*PDstandardNth23alpha + J32L*PDstandardNth33alpha) + 
        dJ323L*PDstandardNth3alpha;
      
      JacPDstandardNth32beta1 = J13L*(J12L*PDstandardNth11beta1 + 
        J22L*PDstandardNth12beta1 + J32L*PDstandardNth13beta1) + 
        J12L*(J23L*PDstandardNth12beta1 + J33L*PDstandardNth13beta1) + 
        dJ123L*PDstandardNth1beta1 + J23L*(J22L*PDstandardNth22beta1 + 
        J32L*PDstandardNth23beta1) + dJ223L*PDstandardNth2beta1 + 
        J33L*(J22L*PDstandardNth23beta1 + J32L*PDstandardNth33beta1) + 
        dJ323L*PDstandardNth3beta1;
      
      JacPDstandardNth32beta2 = J13L*(J12L*PDstandardNth11beta2 + 
        J22L*PDstandardNth12beta2 + J32L*PDstandardNth13beta2) + 
        J12L*(J23L*PDstandardNth12beta2 + J33L*PDstandardNth13beta2) + 
        dJ123L*PDstandardNth1beta2 + J23L*(J22L*PDstandardNth22beta2 + 
        J32L*PDstandardNth23beta2) + dJ223L*PDstandardNth2beta2 + 
        J33L*(J22L*PDstandardNth23beta2 + J32L*PDstandardNth33beta2) + 
        dJ323L*PDstandardNth3beta2;
      
      JacPDstandardNth32beta3 = J13L*(J12L*PDstandardNth11beta3 + 
        J22L*PDstandardNth12beta3 + J32L*PDstandardNth13beta3) + 
        J12L*(J23L*PDstandardNth12beta3 + J33L*PDstandardNth13beta3) + 
        dJ123L*PDstandardNth1beta3 + J23L*(J22L*PDstandardNth22beta3 + 
        J32L*PDstandardNth23beta3) + dJ223L*PDstandardNth2beta3 + 
        J33L*(J22L*PDstandardNth23beta3 + J32L*PDstandardNth33beta3) + 
        dJ323L*PDstandardNth3beta3;
    }
    else
    {
      JacPDstandardNth1alpha = PDstandardNth1alpha;
      
      JacPDstandardNth1beta1 = PDstandardNth1beta1;
      
      JacPDstandardNth1beta2 = PDstandardNth1beta2;
      
      JacPDstandardNth1beta3 = PDstandardNth1beta3;
      
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth1phi = PDstandardNth1phi;
      
      JacPDstandardNth1trK = PDstandardNth1trK;
      
      JacPDstandardNth2alpha = PDstandardNth2alpha;
      
      JacPDstandardNth2beta1 = PDstandardNth2beta1;
      
      JacPDstandardNth2beta2 = PDstandardNth2beta2;
      
      JacPDstandardNth2beta3 = PDstandardNth2beta3;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth2phi = PDstandardNth2phi;
      
      JacPDstandardNth2trK = PDstandardNth2trK;
      
      JacPDstandardNth3alpha = PDstandardNth3alpha;
      
      JacPDstandardNth3beta1 = PDstandardNth3beta1;
      
      JacPDstandardNth3beta2 = PDstandardNth3beta2;
      
      JacPDstandardNth3beta3 = PDstandardNth3beta3;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDstandardNth3phi = PDstandardNth3phi;
      
      JacPDstandardNth3trK = PDstandardNth3trK;
      
      JacPDstandardNth11alpha = PDstandardNth11alpha;
      
      JacPDstandardNth11beta1 = PDstandardNth11beta1;
      
      JacPDstandardNth11beta2 = PDstandardNth11beta2;
      
      JacPDstandardNth11beta3 = PDstandardNth11beta3;
      
      JacPDstandardNth22alpha = PDstandardNth22alpha;
      
      JacPDstandardNth22beta1 = PDstandardNth22beta1;
      
      JacPDstandardNth22beta2 = PDstandardNth22beta2;
      
      JacPDstandardNth22beta3 = PDstandardNth22beta3;
      
      JacPDstandardNth33alpha = PDstandardNth33alpha;
      
      JacPDstandardNth33beta1 = PDstandardNth33beta1;
      
      JacPDstandardNth33beta2 = PDstandardNth33beta2;
      
      JacPDstandardNth33beta3 = PDstandardNth33beta3;
      
      JacPDstandardNth12alpha = PDstandardNth12alpha;
      
      JacPDstandardNth12beta1 = PDstandardNth12beta1;
      
      JacPDstandardNth12beta2 = PDstandardNth12beta2;
      
      JacPDstandardNth12beta3 = PDstandardNth12beta3;
      
      JacPDstandardNth13alpha = PDstandardNth13alpha;
      
      JacPDstandardNth13beta1 = PDstandardNth13beta1;
      
      JacPDstandardNth13beta2 = PDstandardNth13beta2;
      
      JacPDstandardNth13beta3 = PDstandardNth13beta3;
      
      JacPDstandardNth21alpha = PDstandardNth12alpha;
      
      JacPDstandardNth21beta1 = PDstandardNth12beta1;
      
      JacPDstandardNth21beta2 = PDstandardNth12beta2;
      
      JacPDstandardNth21beta3 = PDstandardNth12beta3;
      
      JacPDstandardNth23alpha = PDstandardNth23alpha;
      
      JacPDstandardNth23beta1 = PDstandardNth23beta1;
      
      JacPDstandardNth23beta2 = PDstandardNth23beta2;
      
      JacPDstandardNth23beta3 = PDstandardNth23beta3;
      
      JacPDstandardNth31alpha = PDstandardNth13alpha;
      
      JacPDstandardNth31beta1 = PDstandardNth13beta1;
      
      JacPDstandardNth31beta2 = PDstandardNth13beta2;
      
      JacPDstandardNth31beta3 = PDstandardNth13beta3;
      
      JacPDstandardNth32alpha = PDstandardNth23alpha;
      
      JacPDstandardNth32beta1 = PDstandardNth23beta1;
      
      JacPDstandardNth32beta2 = PDstandardNth23beta2;
      
      JacPDstandardNth32beta3 = PDstandardNth23beta3;
    }
    
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
    
    CCTK_REAL cdphi1 = fac1*JacPDstandardNth1phi;
    
    CCTK_REAL cdphi2 = fac1*JacPDstandardNth2phi;
    
    CCTK_REAL cdphi3 = fac1*JacPDstandardNth3phi;
    
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
      IfThen(conformalMethod,phiL*(-0.333333333333333333333333333333*(JacPDstandardNth1beta1 
      + JacPDstandardNth2beta2 + JacPDstandardNth3beta3) + 
      0.333333333333333333333333333333*alphaL*trKL),0.166666666666666666666666666667*(JacPDstandardNth1beta1 
      + JacPDstandardNth2beta2 + JacPDstandardNth3beta3) - 
      0.166666666666666666666666666667*alphaL*trKL);
    
    CCTK_REAL gt11rhsL = -0.666666666666666666666666666667*(3*alphaL*At11L 
      - 3*(gt12L*JacPDstandardNth1beta2 + gt13L*JacPDstandardNth1beta3) + 
      gt11L*(-2*JacPDstandardNth1beta1 + JacPDstandardNth2beta2 + 
      JacPDstandardNth3beta3));
    
    CCTK_REAL gt12rhsL = 0.333333333333333333333333333333*(-6*alphaL*At12L 
      + 3*(gt22L*JacPDstandardNth1beta2 + gt23L*JacPDstandardNth1beta3 + 
      gt11L*JacPDstandardNth2beta1 + gt13L*JacPDstandardNth2beta3) + 
      gt12L*(JacPDstandardNth1beta1 + JacPDstandardNth2beta2 - 
      2*JacPDstandardNth3beta3));
    
    CCTK_REAL gt13rhsL = 0.333333333333333333333333333333*(-6*alphaL*At13L 
      + 3*(gt23L*JacPDstandardNth1beta2 + gt33L*JacPDstandardNth1beta3 + 
      gt11L*JacPDstandardNth3beta1 + gt12L*JacPDstandardNth3beta2) + 
      gt13L*(JacPDstandardNth1beta1 - 2*JacPDstandardNth2beta2 + 
      JacPDstandardNth3beta3));
    
    CCTK_REAL gt22rhsL = -0.666666666666666666666666666667*(3*alphaL*At22L 
      - 3*(gt12L*JacPDstandardNth2beta1 + gt23L*JacPDstandardNth2beta3) + 
      gt22L*(JacPDstandardNth1beta1 - 2*JacPDstandardNth2beta2 + 
      JacPDstandardNth3beta3));
    
    CCTK_REAL gt23rhsL = 0.333333333333333333333333333333*(-6*alphaL*At23L 
      + 3*(gt13L*JacPDstandardNth2beta1 + gt33L*JacPDstandardNth2beta3 + 
      gt12L*JacPDstandardNth3beta1 + gt22L*JacPDstandardNth3beta2) + 
      gt23L*(-2*JacPDstandardNth1beta1 + JacPDstandardNth2beta2 + 
      JacPDstandardNth3beta3));
    
    CCTK_REAL gt33rhsL = -0.666666666666666666666666666667*(3*alphaL*At33L 
      - 3*(gt13L*JacPDstandardNth3beta1 + gt23L*JacPDstandardNth3beta2) + 
      gt33L*(JacPDstandardNth1beta1 + JacPDstandardNth2beta2 - 
      2*JacPDstandardNth3beta3));
    
    CCTK_REAL dotXt1 = 
      0.333333333333333333333333333333*(gtu11*(4*JacPDstandardNth11beta1 + 
      JacPDstandardNth12beta2 + JacPDstandardNth13beta3) + 
      gtu12*(4*JacPDstandardNth21beta1 + JacPDstandardNth22beta2 + 
      JacPDstandardNth23beta3) + 3*(gtu12*JacPDstandardNth12beta1 + 
      gtu13*JacPDstandardNth13beta1 + gtu22*JacPDstandardNth22beta1 + 
      gtu23*(JacPDstandardNth23beta1 + JacPDstandardNth32beta1) + 
      gtu33*JacPDstandardNth33beta1) + gtu13*(4*JacPDstandardNth31beta1 + 
      JacPDstandardNth32beta2 + JacPDstandardNth33beta3) - 
      6*(Atu11*JacPDstandardNth1alpha + Atu12*JacPDstandardNth2alpha + 
      Atu13*JacPDstandardNth3alpha) + alphaL*(Atu13*(36*cdphi3 + 12*Gt113) + 
      12*(Atu12*(3*cdphi2 + Gt112) + Atu23*Gt123) + 6*(Atu11*(6*cdphi1 + 
      Gt111) + Atu22*Gt122 + Atu33*Gt133) - 4*(gtu11*JacPDstandardNth1trK + 
      gtu12*JacPDstandardNth2trK + gtu13*JacPDstandardNth3trK) - 
      150.7964473723100754462068823974161384415*(gtu11*S1 + gtu12*S2 + 
      gtu13*S3)) + (-JacPDstandardNth1beta1 + 2*(JacPDstandardNth2beta2 + 
      JacPDstandardNth3beta3))*Xtn1 - 3*(JacPDstandardNth2beta1*Xtn2 + 
      JacPDstandardNth3beta1*Xtn3));
    
    CCTK_REAL dotXt2 = 
      0.333333333333333333333333333333*(gtu12*(JacPDstandardNth11beta1 + 
      4*JacPDstandardNth12beta2 + JacPDstandardNth13beta3 + 
      3*JacPDstandardNth21beta2) + gtu22*(JacPDstandardNth21beta1 + 
      4*JacPDstandardNth22beta2 + JacPDstandardNth23beta3) + 
      3*(gtu11*JacPDstandardNth11beta2 + gtu23*JacPDstandardNth23beta2 + 
      gtu13*(JacPDstandardNth13beta2 + JacPDstandardNth31beta2) + 
      gtu33*JacPDstandardNth33beta2) + gtu23*(JacPDstandardNth31beta1 + 
      4*JacPDstandardNth32beta2 + JacPDstandardNth33beta3) - 
      6*(Atu12*JacPDstandardNth1alpha + Atu22*JacPDstandardNth2alpha + 
      Atu23*JacPDstandardNth3alpha) + alphaL*(12*(Atu12*(3*cdphi1 + Gt212) + 
      Atu13*Gt213) + Atu23*(36*cdphi3 + 12*Gt223) + 6*(Atu11*Gt211 + 
      Atu22*(6*cdphi2 + Gt222) + Atu33*Gt233) - 4*(gtu12*JacPDstandardNth1trK 
      + gtu22*JacPDstandardNth2trK + gtu23*JacPDstandardNth3trK) - 
      150.7964473723100754462068823974161384415*(gtu12*S1 + gtu22*S2 + 
      gtu23*S3)) + (-JacPDstandardNth2beta2 + 2*(JacPDstandardNth1beta1 + 
      JacPDstandardNth3beta3))*Xtn2 - 3*(JacPDstandardNth1beta2*Xtn1 + 
      JacPDstandardNth3beta2*Xtn3));
    
    CCTK_REAL dotXt3 = 
      0.333333333333333333333333333333*(gtu23*(JacPDstandardNth21beta1 + 
      JacPDstandardNth22beta2 + 4*JacPDstandardNth23beta3) + 
      gtu13*(JacPDstandardNth11beta1 + JacPDstandardNth12beta2 + 
      4*JacPDstandardNth13beta3 + 3*JacPDstandardNth31beta3) + 
      3*(gtu11*JacPDstandardNth11beta3 + gtu12*(JacPDstandardNth12beta3 + 
      JacPDstandardNth21beta3) + gtu22*JacPDstandardNth22beta3 + 
      gtu23*JacPDstandardNth32beta3) + gtu33*(JacPDstandardNth31beta1 + 
      JacPDstandardNth32beta2 + 4*JacPDstandardNth33beta3) - 
      6*(Atu13*JacPDstandardNth1alpha + Atu23*JacPDstandardNth2alpha + 
      Atu33*JacPDstandardNth3alpha) + alphaL*(6*(Atu11*Gt311 + Atu22*Gt322) + 
      12*(Atu12*Gt312 + Atu13*(3*cdphi1 + Gt313) + Atu23*(3*cdphi2 + Gt323)) 
      + Atu33*(36*cdphi3 + 6*Gt333) - 4*(gtu13*JacPDstandardNth1trK + 
      gtu23*JacPDstandardNth2trK + gtu33*JacPDstandardNth3trK) - 
      150.7964473723100754462068823974161384415*(gtu13*S1 + gtu23*S2 + 
      gtu33*S3)) - 3*(JacPDstandardNth1beta3*Xtn1 + 
      JacPDstandardNth2beta3*Xtn2) + (2*(JacPDstandardNth1beta1 + 
      JacPDstandardNth2beta2) - JacPDstandardNth3beta3)*Xtn3);
    
    CCTK_REAL Xt1rhsL = dotXt1;
    
    CCTK_REAL Xt2rhsL = dotXt2;
    
    CCTK_REAL Xt3rhsL = dotXt3;
    
    CCTK_REAL dottrK = -(em4phi*(gtu11*(JacPDstandardNth11alpha + 
      2*cdphi1*JacPDstandardNth1alpha) + gtu12*(JacPDstandardNth12alpha + 
      2*cdphi2*JacPDstandardNth1alpha + JacPDstandardNth21alpha + 
      2*cdphi1*JacPDstandardNth2alpha) + gtu22*(JacPDstandardNth22alpha + 
      2*cdphi2*JacPDstandardNth2alpha) + gtu13*(JacPDstandardNth13alpha + 
      2*cdphi3*JacPDstandardNth1alpha + JacPDstandardNth31alpha + 
      2*cdphi1*JacPDstandardNth3alpha) + gtu23*(JacPDstandardNth23alpha + 
      2*cdphi3*JacPDstandardNth2alpha + JacPDstandardNth32alpha + 
      2*cdphi2*JacPDstandardNth3alpha) + gtu33*(JacPDstandardNth33alpha + 
      2*cdphi3*JacPDstandardNth3alpha) - JacPDstandardNth1alpha*Xtn1 - 
      JacPDstandardNth2alpha*Xtn2 - JacPDstandardNth3alpha*Xtn3)) + 
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
  LC_ENDLOOP3 (ML_BSSN_MP_RHS1);
}

extern "C" void ML_BSSN_MP_RHS1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_RHS1_Body);
}
