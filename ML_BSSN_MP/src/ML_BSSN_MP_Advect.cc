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

extern "C" void ML_BSSN_MP_Advect_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_curvrhs.");
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

static void ML_BSSN_MP_Advect_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_Advect_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_Advect_calc_every != ML_BSSN_MP_Advect_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_BSSN_MP::ML_curv","ML_BSSN_MP::ML_curvrhs","ML_BSSN_MP::ML_dtlapse","ML_BSSN_MP::ML_dtlapserhs","ML_BSSN_MP::ML_dtshift","ML_BSSN_MP::ML_dtshiftrhs","ML_BSSN_MP::ML_Gamma","ML_BSSN_MP::ML_Gammarhs","ML_BSSN_MP::ML_lapse","ML_BSSN_MP::ML_lapserhs","ML_BSSN_MP::ML_log_confac","ML_BSSN_MP::ML_log_confacrhs","ML_BSSN_MP::ML_metric","ML_BSSN_MP::ML_metricrhs","ML_BSSN_MP::ML_shift","ML_BSSN_MP::ML_shiftrhs","ML_BSSN_MP::ML_trace_curv","ML_BSSN_MP::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_Advect", 18, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_MP_Advect", 3, 3, 3);
  
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
  LC_LOOP3 (ML_BSSN_MP_Advect,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL AL = A[index];
    CCTK_REAL alphaL = alpha[index];
    CCTK_REAL alpharhsL = alpharhs[index];
    CCTK_REAL ArhsL = Arhs[index];
    CCTK_REAL At11L = At11[index];
    CCTK_REAL At11rhsL = At11rhs[index];
    CCTK_REAL At12L = At12[index];
    CCTK_REAL At12rhsL = At12rhs[index];
    CCTK_REAL At13L = At13[index];
    CCTK_REAL At13rhsL = At13rhs[index];
    CCTK_REAL At22L = At22[index];
    CCTK_REAL At22rhsL = At22rhs[index];
    CCTK_REAL At23L = At23[index];
    CCTK_REAL At23rhsL = At23rhs[index];
    CCTK_REAL At33L = At33[index];
    CCTK_REAL At33rhsL = At33rhs[index];
    CCTK_REAL B1L = B1[index];
    CCTK_REAL B1rhsL = B1rhs[index];
    CCTK_REAL B2L = B2[index];
    CCTK_REAL B2rhsL = B2rhs[index];
    CCTK_REAL B3L = B3[index];
    CCTK_REAL B3rhsL = B3rhs[index];
    CCTK_REAL beta1L = beta1[index];
    CCTK_REAL beta1rhsL = beta1rhs[index];
    CCTK_REAL beta2L = beta2[index];
    CCTK_REAL beta2rhsL = beta2rhs[index];
    CCTK_REAL beta3L = beta3[index];
    CCTK_REAL beta3rhsL = beta3rhs[index];
    CCTK_REAL gt11L = gt11[index];
    CCTK_REAL gt11rhsL = gt11rhs[index];
    CCTK_REAL gt12L = gt12[index];
    CCTK_REAL gt12rhsL = gt12rhs[index];
    CCTK_REAL gt13L = gt13[index];
    CCTK_REAL gt13rhsL = gt13rhs[index];
    CCTK_REAL gt22L = gt22[index];
    CCTK_REAL gt22rhsL = gt22rhs[index];
    CCTK_REAL gt23L = gt23[index];
    CCTK_REAL gt23rhsL = gt23rhs[index];
    CCTK_REAL gt33L = gt33[index];
    CCTK_REAL gt33rhsL = gt33rhs[index];
    CCTK_REAL phiL = phi[index];
    CCTK_REAL phirhsL = phirhs[index];
    CCTK_REAL trKL = trK[index];
    CCTK_REAL trKrhsL = trKrhs[index];
    CCTK_REAL Xt1L = Xt1[index];
    CCTK_REAL Xt1rhsL = Xt1rhs[index];
    CCTK_REAL Xt2L = Xt2[index];
    CCTK_REAL Xt2rhsL = Xt2rhs[index];
    CCTK_REAL Xt3L = Xt3[index];
    CCTK_REAL Xt3rhsL = Xt3rhs[index];
    
    
    CCTK_REAL J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L;
    
    if (use_jacobian)
    {
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
    CCTK_REAL const PDupwindNthAnti1A = PDupwindNthAnti1(&A[index]);
    CCTK_REAL const PDupwindNthSymm1A = PDupwindNthSymm1(&A[index]);
    CCTK_REAL const PDupwindNthAnti2A = PDupwindNthAnti2(&A[index]);
    CCTK_REAL const PDupwindNthSymm2A = PDupwindNthSymm2(&A[index]);
    CCTK_REAL const PDupwindNthAnti3A = PDupwindNthAnti3(&A[index]);
    CCTK_REAL const PDupwindNthSymm3A = PDupwindNthSymm3(&A[index]);
    CCTK_REAL const PDupwindNthAnti1alpha = PDupwindNthAnti1(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm1alpha = PDupwindNthSymm1(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti2alpha = PDupwindNthAnti2(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm2alpha = PDupwindNthSymm2(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti3alpha = PDupwindNthAnti3(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm3alpha = PDupwindNthSymm3(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti1At11 = PDupwindNthAnti1(&At11[index]);
    CCTK_REAL const PDupwindNthSymm1At11 = PDupwindNthSymm1(&At11[index]);
    CCTK_REAL const PDupwindNthAnti2At11 = PDupwindNthAnti2(&At11[index]);
    CCTK_REAL const PDupwindNthSymm2At11 = PDupwindNthSymm2(&At11[index]);
    CCTK_REAL const PDupwindNthAnti3At11 = PDupwindNthAnti3(&At11[index]);
    CCTK_REAL const PDupwindNthSymm3At11 = PDupwindNthSymm3(&At11[index]);
    CCTK_REAL const PDupwindNthAnti1At12 = PDupwindNthAnti1(&At12[index]);
    CCTK_REAL const PDupwindNthSymm1At12 = PDupwindNthSymm1(&At12[index]);
    CCTK_REAL const PDupwindNthAnti2At12 = PDupwindNthAnti2(&At12[index]);
    CCTK_REAL const PDupwindNthSymm2At12 = PDupwindNthSymm2(&At12[index]);
    CCTK_REAL const PDupwindNthAnti3At12 = PDupwindNthAnti3(&At12[index]);
    CCTK_REAL const PDupwindNthSymm3At12 = PDupwindNthSymm3(&At12[index]);
    CCTK_REAL const PDupwindNthAnti1At13 = PDupwindNthAnti1(&At13[index]);
    CCTK_REAL const PDupwindNthSymm1At13 = PDupwindNthSymm1(&At13[index]);
    CCTK_REAL const PDupwindNthAnti2At13 = PDupwindNthAnti2(&At13[index]);
    CCTK_REAL const PDupwindNthSymm2At13 = PDupwindNthSymm2(&At13[index]);
    CCTK_REAL const PDupwindNthAnti3At13 = PDupwindNthAnti3(&At13[index]);
    CCTK_REAL const PDupwindNthSymm3At13 = PDupwindNthSymm3(&At13[index]);
    CCTK_REAL const PDupwindNthAnti1At22 = PDupwindNthAnti1(&At22[index]);
    CCTK_REAL const PDupwindNthSymm1At22 = PDupwindNthSymm1(&At22[index]);
    CCTK_REAL const PDupwindNthAnti2At22 = PDupwindNthAnti2(&At22[index]);
    CCTK_REAL const PDupwindNthSymm2At22 = PDupwindNthSymm2(&At22[index]);
    CCTK_REAL const PDupwindNthAnti3At22 = PDupwindNthAnti3(&At22[index]);
    CCTK_REAL const PDupwindNthSymm3At22 = PDupwindNthSymm3(&At22[index]);
    CCTK_REAL const PDupwindNthAnti1At23 = PDupwindNthAnti1(&At23[index]);
    CCTK_REAL const PDupwindNthSymm1At23 = PDupwindNthSymm1(&At23[index]);
    CCTK_REAL const PDupwindNthAnti2At23 = PDupwindNthAnti2(&At23[index]);
    CCTK_REAL const PDupwindNthSymm2At23 = PDupwindNthSymm2(&At23[index]);
    CCTK_REAL const PDupwindNthAnti3At23 = PDupwindNthAnti3(&At23[index]);
    CCTK_REAL const PDupwindNthSymm3At23 = PDupwindNthSymm3(&At23[index]);
    CCTK_REAL const PDupwindNthAnti1At33 = PDupwindNthAnti1(&At33[index]);
    CCTK_REAL const PDupwindNthSymm1At33 = PDupwindNthSymm1(&At33[index]);
    CCTK_REAL const PDupwindNthAnti2At33 = PDupwindNthAnti2(&At33[index]);
    CCTK_REAL const PDupwindNthSymm2At33 = PDupwindNthSymm2(&At33[index]);
    CCTK_REAL const PDupwindNthAnti3At33 = PDupwindNthAnti3(&At33[index]);
    CCTK_REAL const PDupwindNthSymm3At33 = PDupwindNthSymm3(&At33[index]);
    CCTK_REAL const PDupwindNthAnti1B1 = PDupwindNthAnti1(&B1[index]);
    CCTK_REAL const PDupwindNthSymm1B1 = PDupwindNthSymm1(&B1[index]);
    CCTK_REAL const PDupwindNthAnti2B1 = PDupwindNthAnti2(&B1[index]);
    CCTK_REAL const PDupwindNthSymm2B1 = PDupwindNthSymm2(&B1[index]);
    CCTK_REAL const PDupwindNthAnti3B1 = PDupwindNthAnti3(&B1[index]);
    CCTK_REAL const PDupwindNthSymm3B1 = PDupwindNthSymm3(&B1[index]);
    CCTK_REAL const PDupwindNthAnti1B2 = PDupwindNthAnti1(&B2[index]);
    CCTK_REAL const PDupwindNthSymm1B2 = PDupwindNthSymm1(&B2[index]);
    CCTK_REAL const PDupwindNthAnti2B2 = PDupwindNthAnti2(&B2[index]);
    CCTK_REAL const PDupwindNthSymm2B2 = PDupwindNthSymm2(&B2[index]);
    CCTK_REAL const PDupwindNthAnti3B2 = PDupwindNthAnti3(&B2[index]);
    CCTK_REAL const PDupwindNthSymm3B2 = PDupwindNthSymm3(&B2[index]);
    CCTK_REAL const PDupwindNthAnti1B3 = PDupwindNthAnti1(&B3[index]);
    CCTK_REAL const PDupwindNthSymm1B3 = PDupwindNthSymm1(&B3[index]);
    CCTK_REAL const PDupwindNthAnti2B3 = PDupwindNthAnti2(&B3[index]);
    CCTK_REAL const PDupwindNthSymm2B3 = PDupwindNthSymm2(&B3[index]);
    CCTK_REAL const PDupwindNthAnti3B3 = PDupwindNthAnti3(&B3[index]);
    CCTK_REAL const PDupwindNthSymm3B3 = PDupwindNthSymm3(&B3[index]);
    CCTK_REAL const PDupwindNthAnti1beta1 = PDupwindNthAnti1(&beta1[index]);
    CCTK_REAL const PDupwindNthSymm1beta1 = PDupwindNthSymm1(&beta1[index]);
    CCTK_REAL const PDupwindNthAnti2beta1 = PDupwindNthAnti2(&beta1[index]);
    CCTK_REAL const PDupwindNthSymm2beta1 = PDupwindNthSymm2(&beta1[index]);
    CCTK_REAL const PDupwindNthAnti3beta1 = PDupwindNthAnti3(&beta1[index]);
    CCTK_REAL const PDupwindNthSymm3beta1 = PDupwindNthSymm3(&beta1[index]);
    CCTK_REAL const PDupwindNthAnti1beta2 = PDupwindNthAnti1(&beta2[index]);
    CCTK_REAL const PDupwindNthSymm1beta2 = PDupwindNthSymm1(&beta2[index]);
    CCTK_REAL const PDupwindNthAnti2beta2 = PDupwindNthAnti2(&beta2[index]);
    CCTK_REAL const PDupwindNthSymm2beta2 = PDupwindNthSymm2(&beta2[index]);
    CCTK_REAL const PDupwindNthAnti3beta2 = PDupwindNthAnti3(&beta2[index]);
    CCTK_REAL const PDupwindNthSymm3beta2 = PDupwindNthSymm3(&beta2[index]);
    CCTK_REAL const PDupwindNthAnti1beta3 = PDupwindNthAnti1(&beta3[index]);
    CCTK_REAL const PDupwindNthSymm1beta3 = PDupwindNthSymm1(&beta3[index]);
    CCTK_REAL const PDupwindNthAnti2beta3 = PDupwindNthAnti2(&beta3[index]);
    CCTK_REAL const PDupwindNthSymm2beta3 = PDupwindNthSymm2(&beta3[index]);
    CCTK_REAL const PDupwindNthAnti3beta3 = PDupwindNthAnti3(&beta3[index]);
    CCTK_REAL const PDupwindNthSymm3beta3 = PDupwindNthSymm3(&beta3[index]);
    CCTK_REAL const PDupwindNthAnti1gt11 = PDupwindNthAnti1(&gt11[index]);
    CCTK_REAL const PDupwindNthSymm1gt11 = PDupwindNthSymm1(&gt11[index]);
    CCTK_REAL const PDupwindNthAnti2gt11 = PDupwindNthAnti2(&gt11[index]);
    CCTK_REAL const PDupwindNthSymm2gt11 = PDupwindNthSymm2(&gt11[index]);
    CCTK_REAL const PDupwindNthAnti3gt11 = PDupwindNthAnti3(&gt11[index]);
    CCTK_REAL const PDupwindNthSymm3gt11 = PDupwindNthSymm3(&gt11[index]);
    CCTK_REAL const PDupwindNthAnti1gt12 = PDupwindNthAnti1(&gt12[index]);
    CCTK_REAL const PDupwindNthSymm1gt12 = PDupwindNthSymm1(&gt12[index]);
    CCTK_REAL const PDupwindNthAnti2gt12 = PDupwindNthAnti2(&gt12[index]);
    CCTK_REAL const PDupwindNthSymm2gt12 = PDupwindNthSymm2(&gt12[index]);
    CCTK_REAL const PDupwindNthAnti3gt12 = PDupwindNthAnti3(&gt12[index]);
    CCTK_REAL const PDupwindNthSymm3gt12 = PDupwindNthSymm3(&gt12[index]);
    CCTK_REAL const PDupwindNthAnti1gt13 = PDupwindNthAnti1(&gt13[index]);
    CCTK_REAL const PDupwindNthSymm1gt13 = PDupwindNthSymm1(&gt13[index]);
    CCTK_REAL const PDupwindNthAnti2gt13 = PDupwindNthAnti2(&gt13[index]);
    CCTK_REAL const PDupwindNthSymm2gt13 = PDupwindNthSymm2(&gt13[index]);
    CCTK_REAL const PDupwindNthAnti3gt13 = PDupwindNthAnti3(&gt13[index]);
    CCTK_REAL const PDupwindNthSymm3gt13 = PDupwindNthSymm3(&gt13[index]);
    CCTK_REAL const PDupwindNthAnti1gt22 = PDupwindNthAnti1(&gt22[index]);
    CCTK_REAL const PDupwindNthSymm1gt22 = PDupwindNthSymm1(&gt22[index]);
    CCTK_REAL const PDupwindNthAnti2gt22 = PDupwindNthAnti2(&gt22[index]);
    CCTK_REAL const PDupwindNthSymm2gt22 = PDupwindNthSymm2(&gt22[index]);
    CCTK_REAL const PDupwindNthAnti3gt22 = PDupwindNthAnti3(&gt22[index]);
    CCTK_REAL const PDupwindNthSymm3gt22 = PDupwindNthSymm3(&gt22[index]);
    CCTK_REAL const PDupwindNthAnti1gt23 = PDupwindNthAnti1(&gt23[index]);
    CCTK_REAL const PDupwindNthSymm1gt23 = PDupwindNthSymm1(&gt23[index]);
    CCTK_REAL const PDupwindNthAnti2gt23 = PDupwindNthAnti2(&gt23[index]);
    CCTK_REAL const PDupwindNthSymm2gt23 = PDupwindNthSymm2(&gt23[index]);
    CCTK_REAL const PDupwindNthAnti3gt23 = PDupwindNthAnti3(&gt23[index]);
    CCTK_REAL const PDupwindNthSymm3gt23 = PDupwindNthSymm3(&gt23[index]);
    CCTK_REAL const PDupwindNthAnti1gt33 = PDupwindNthAnti1(&gt33[index]);
    CCTK_REAL const PDupwindNthSymm1gt33 = PDupwindNthSymm1(&gt33[index]);
    CCTK_REAL const PDupwindNthAnti2gt33 = PDupwindNthAnti2(&gt33[index]);
    CCTK_REAL const PDupwindNthSymm2gt33 = PDupwindNthSymm2(&gt33[index]);
    CCTK_REAL const PDupwindNthAnti3gt33 = PDupwindNthAnti3(&gt33[index]);
    CCTK_REAL const PDupwindNthSymm3gt33 = PDupwindNthSymm3(&gt33[index]);
    CCTK_REAL const PDupwindNthAnti1phi = PDupwindNthAnti1(&phi[index]);
    CCTK_REAL const PDupwindNthSymm1phi = PDupwindNthSymm1(&phi[index]);
    CCTK_REAL const PDupwindNthAnti2phi = PDupwindNthAnti2(&phi[index]);
    CCTK_REAL const PDupwindNthSymm2phi = PDupwindNthSymm2(&phi[index]);
    CCTK_REAL const PDupwindNthAnti3phi = PDupwindNthAnti3(&phi[index]);
    CCTK_REAL const PDupwindNthSymm3phi = PDupwindNthSymm3(&phi[index]);
    CCTK_REAL const PDupwindNthAnti1trK = PDupwindNthAnti1(&trK[index]);
    CCTK_REAL const PDupwindNthSymm1trK = PDupwindNthSymm1(&trK[index]);
    CCTK_REAL const PDupwindNthAnti2trK = PDupwindNthAnti2(&trK[index]);
    CCTK_REAL const PDupwindNthSymm2trK = PDupwindNthSymm2(&trK[index]);
    CCTK_REAL const PDupwindNthAnti3trK = PDupwindNthAnti3(&trK[index]);
    CCTK_REAL const PDupwindNthSymm3trK = PDupwindNthSymm3(&trK[index]);
    CCTK_REAL const PDupwindNthAnti1Xt1 = PDupwindNthAnti1(&Xt1[index]);
    CCTK_REAL const PDupwindNthSymm1Xt1 = PDupwindNthSymm1(&Xt1[index]);
    CCTK_REAL const PDupwindNthAnti2Xt1 = PDupwindNthAnti2(&Xt1[index]);
    CCTK_REAL const PDupwindNthSymm2Xt1 = PDupwindNthSymm2(&Xt1[index]);
    CCTK_REAL const PDupwindNthAnti3Xt1 = PDupwindNthAnti3(&Xt1[index]);
    CCTK_REAL const PDupwindNthSymm3Xt1 = PDupwindNthSymm3(&Xt1[index]);
    CCTK_REAL const PDupwindNthAnti1Xt2 = PDupwindNthAnti1(&Xt2[index]);
    CCTK_REAL const PDupwindNthSymm1Xt2 = PDupwindNthSymm1(&Xt2[index]);
    CCTK_REAL const PDupwindNthAnti2Xt2 = PDupwindNthAnti2(&Xt2[index]);
    CCTK_REAL const PDupwindNthSymm2Xt2 = PDupwindNthSymm2(&Xt2[index]);
    CCTK_REAL const PDupwindNthAnti3Xt2 = PDupwindNthAnti3(&Xt2[index]);
    CCTK_REAL const PDupwindNthSymm3Xt2 = PDupwindNthSymm3(&Xt2[index]);
    CCTK_REAL const PDupwindNthAnti1Xt3 = PDupwindNthAnti1(&Xt3[index]);
    CCTK_REAL const PDupwindNthSymm1Xt3 = PDupwindNthSymm1(&Xt3[index]);
    CCTK_REAL const PDupwindNthAnti2Xt3 = PDupwindNthAnti2(&Xt3[index]);
    CCTK_REAL const PDupwindNthSymm2Xt3 = PDupwindNthSymm2(&Xt3[index]);
    CCTK_REAL const PDupwindNthAnti3Xt3 = PDupwindNthAnti3(&Xt3[index]);
    CCTK_REAL const PDupwindNthSymm3Xt3 = PDupwindNthSymm3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL JacPDupwindNthAnti1A;
    CCTK_REAL JacPDupwindNthAnti1alpha;
    CCTK_REAL JacPDupwindNthAnti1At11;
    CCTK_REAL JacPDupwindNthAnti1At12;
    CCTK_REAL JacPDupwindNthAnti1At13;
    CCTK_REAL JacPDupwindNthAnti1At22;
    CCTK_REAL JacPDupwindNthAnti1At23;
    CCTK_REAL JacPDupwindNthAnti1At33;
    CCTK_REAL JacPDupwindNthAnti1B1;
    CCTK_REAL JacPDupwindNthAnti1B2;
    CCTK_REAL JacPDupwindNthAnti1B3;
    CCTK_REAL JacPDupwindNthAnti1beta1;
    CCTK_REAL JacPDupwindNthAnti1beta2;
    CCTK_REAL JacPDupwindNthAnti1beta3;
    CCTK_REAL JacPDupwindNthAnti1gt11;
    CCTK_REAL JacPDupwindNthAnti1gt12;
    CCTK_REAL JacPDupwindNthAnti1gt13;
    CCTK_REAL JacPDupwindNthAnti1gt22;
    CCTK_REAL JacPDupwindNthAnti1gt23;
    CCTK_REAL JacPDupwindNthAnti1gt33;
    CCTK_REAL JacPDupwindNthAnti1phi;
    CCTK_REAL JacPDupwindNthAnti1trK;
    CCTK_REAL JacPDupwindNthAnti1Xt1;
    CCTK_REAL JacPDupwindNthAnti1Xt2;
    CCTK_REAL JacPDupwindNthAnti1Xt3;
    CCTK_REAL JacPDupwindNthAnti2A;
    CCTK_REAL JacPDupwindNthAnti2alpha;
    CCTK_REAL JacPDupwindNthAnti2At11;
    CCTK_REAL JacPDupwindNthAnti2At12;
    CCTK_REAL JacPDupwindNthAnti2At13;
    CCTK_REAL JacPDupwindNthAnti2At22;
    CCTK_REAL JacPDupwindNthAnti2At23;
    CCTK_REAL JacPDupwindNthAnti2At33;
    CCTK_REAL JacPDupwindNthAnti2B1;
    CCTK_REAL JacPDupwindNthAnti2B2;
    CCTK_REAL JacPDupwindNthAnti2B3;
    CCTK_REAL JacPDupwindNthAnti2beta1;
    CCTK_REAL JacPDupwindNthAnti2beta2;
    CCTK_REAL JacPDupwindNthAnti2beta3;
    CCTK_REAL JacPDupwindNthAnti2gt11;
    CCTK_REAL JacPDupwindNthAnti2gt12;
    CCTK_REAL JacPDupwindNthAnti2gt13;
    CCTK_REAL JacPDupwindNthAnti2gt22;
    CCTK_REAL JacPDupwindNthAnti2gt23;
    CCTK_REAL JacPDupwindNthAnti2gt33;
    CCTK_REAL JacPDupwindNthAnti2phi;
    CCTK_REAL JacPDupwindNthAnti2trK;
    CCTK_REAL JacPDupwindNthAnti2Xt1;
    CCTK_REAL JacPDupwindNthAnti2Xt2;
    CCTK_REAL JacPDupwindNthAnti2Xt3;
    CCTK_REAL JacPDupwindNthAnti3A;
    CCTK_REAL JacPDupwindNthAnti3alpha;
    CCTK_REAL JacPDupwindNthAnti3At11;
    CCTK_REAL JacPDupwindNthAnti3At12;
    CCTK_REAL JacPDupwindNthAnti3At13;
    CCTK_REAL JacPDupwindNthAnti3At22;
    CCTK_REAL JacPDupwindNthAnti3At23;
    CCTK_REAL JacPDupwindNthAnti3At33;
    CCTK_REAL JacPDupwindNthAnti3B1;
    CCTK_REAL JacPDupwindNthAnti3B2;
    CCTK_REAL JacPDupwindNthAnti3B3;
    CCTK_REAL JacPDupwindNthAnti3beta1;
    CCTK_REAL JacPDupwindNthAnti3beta2;
    CCTK_REAL JacPDupwindNthAnti3beta3;
    CCTK_REAL JacPDupwindNthAnti3gt11;
    CCTK_REAL JacPDupwindNthAnti3gt12;
    CCTK_REAL JacPDupwindNthAnti3gt13;
    CCTK_REAL JacPDupwindNthAnti3gt22;
    CCTK_REAL JacPDupwindNthAnti3gt23;
    CCTK_REAL JacPDupwindNthAnti3gt33;
    CCTK_REAL JacPDupwindNthAnti3phi;
    CCTK_REAL JacPDupwindNthAnti3trK;
    CCTK_REAL JacPDupwindNthAnti3Xt1;
    CCTK_REAL JacPDupwindNthAnti3Xt2;
    CCTK_REAL JacPDupwindNthAnti3Xt3;
    CCTK_REAL JacPDupwindNthSymm1A;
    CCTK_REAL JacPDupwindNthSymm1alpha;
    CCTK_REAL JacPDupwindNthSymm1At11;
    CCTK_REAL JacPDupwindNthSymm1At12;
    CCTK_REAL JacPDupwindNthSymm1At13;
    CCTK_REAL JacPDupwindNthSymm1At22;
    CCTK_REAL JacPDupwindNthSymm1At23;
    CCTK_REAL JacPDupwindNthSymm1At33;
    CCTK_REAL JacPDupwindNthSymm1B1;
    CCTK_REAL JacPDupwindNthSymm1B2;
    CCTK_REAL JacPDupwindNthSymm1B3;
    CCTK_REAL JacPDupwindNthSymm1beta1;
    CCTK_REAL JacPDupwindNthSymm1beta2;
    CCTK_REAL JacPDupwindNthSymm1beta3;
    CCTK_REAL JacPDupwindNthSymm1gt11;
    CCTK_REAL JacPDupwindNthSymm1gt12;
    CCTK_REAL JacPDupwindNthSymm1gt13;
    CCTK_REAL JacPDupwindNthSymm1gt22;
    CCTK_REAL JacPDupwindNthSymm1gt23;
    CCTK_REAL JacPDupwindNthSymm1gt33;
    CCTK_REAL JacPDupwindNthSymm1phi;
    CCTK_REAL JacPDupwindNthSymm1trK;
    CCTK_REAL JacPDupwindNthSymm1Xt1;
    CCTK_REAL JacPDupwindNthSymm1Xt2;
    CCTK_REAL JacPDupwindNthSymm1Xt3;
    CCTK_REAL JacPDupwindNthSymm2A;
    CCTK_REAL JacPDupwindNthSymm2alpha;
    CCTK_REAL JacPDupwindNthSymm2At11;
    CCTK_REAL JacPDupwindNthSymm2At12;
    CCTK_REAL JacPDupwindNthSymm2At13;
    CCTK_REAL JacPDupwindNthSymm2At22;
    CCTK_REAL JacPDupwindNthSymm2At23;
    CCTK_REAL JacPDupwindNthSymm2At33;
    CCTK_REAL JacPDupwindNthSymm2B1;
    CCTK_REAL JacPDupwindNthSymm2B2;
    CCTK_REAL JacPDupwindNthSymm2B3;
    CCTK_REAL JacPDupwindNthSymm2beta1;
    CCTK_REAL JacPDupwindNthSymm2beta2;
    CCTK_REAL JacPDupwindNthSymm2beta3;
    CCTK_REAL JacPDupwindNthSymm2gt11;
    CCTK_REAL JacPDupwindNthSymm2gt12;
    CCTK_REAL JacPDupwindNthSymm2gt13;
    CCTK_REAL JacPDupwindNthSymm2gt22;
    CCTK_REAL JacPDupwindNthSymm2gt23;
    CCTK_REAL JacPDupwindNthSymm2gt33;
    CCTK_REAL JacPDupwindNthSymm2phi;
    CCTK_REAL JacPDupwindNthSymm2trK;
    CCTK_REAL JacPDupwindNthSymm2Xt1;
    CCTK_REAL JacPDupwindNthSymm2Xt2;
    CCTK_REAL JacPDupwindNthSymm2Xt3;
    CCTK_REAL JacPDupwindNthSymm3A;
    CCTK_REAL JacPDupwindNthSymm3alpha;
    CCTK_REAL JacPDupwindNthSymm3At11;
    CCTK_REAL JacPDupwindNthSymm3At12;
    CCTK_REAL JacPDupwindNthSymm3At13;
    CCTK_REAL JacPDupwindNthSymm3At22;
    CCTK_REAL JacPDupwindNthSymm3At23;
    CCTK_REAL JacPDupwindNthSymm3At33;
    CCTK_REAL JacPDupwindNthSymm3B1;
    CCTK_REAL JacPDupwindNthSymm3B2;
    CCTK_REAL JacPDupwindNthSymm3B3;
    CCTK_REAL JacPDupwindNthSymm3beta1;
    CCTK_REAL JacPDupwindNthSymm3beta2;
    CCTK_REAL JacPDupwindNthSymm3beta3;
    CCTK_REAL JacPDupwindNthSymm3gt11;
    CCTK_REAL JacPDupwindNthSymm3gt12;
    CCTK_REAL JacPDupwindNthSymm3gt13;
    CCTK_REAL JacPDupwindNthSymm3gt22;
    CCTK_REAL JacPDupwindNthSymm3gt23;
    CCTK_REAL JacPDupwindNthSymm3gt33;
    CCTK_REAL JacPDupwindNthSymm3phi;
    CCTK_REAL JacPDupwindNthSymm3trK;
    CCTK_REAL JacPDupwindNthSymm3Xt1;
    CCTK_REAL JacPDupwindNthSymm3Xt2;
    CCTK_REAL JacPDupwindNthSymm3Xt3;
    
    if (use_jacobian)
    {
      JacPDupwindNthAnti1A = J11L*PDupwindNthAnti1A + J21L*PDupwindNthAnti2A 
        + J31L*PDupwindNthAnti3A;
      
      JacPDupwindNthAnti1alpha = J11L*PDupwindNthAnti1alpha + 
        J21L*PDupwindNthAnti2alpha + J31L*PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti1At11 = J11L*PDupwindNthAnti1At11 + 
        J21L*PDupwindNthAnti2At11 + J31L*PDupwindNthAnti3At11;
      
      JacPDupwindNthAnti1At12 = J11L*PDupwindNthAnti1At12 + 
        J21L*PDupwindNthAnti2At12 + J31L*PDupwindNthAnti3At12;
      
      JacPDupwindNthAnti1At13 = J11L*PDupwindNthAnti1At13 + 
        J21L*PDupwindNthAnti2At13 + J31L*PDupwindNthAnti3At13;
      
      JacPDupwindNthAnti1At22 = J11L*PDupwindNthAnti1At22 + 
        J21L*PDupwindNthAnti2At22 + J31L*PDupwindNthAnti3At22;
      
      JacPDupwindNthAnti1At23 = J11L*PDupwindNthAnti1At23 + 
        J21L*PDupwindNthAnti2At23 + J31L*PDupwindNthAnti3At23;
      
      JacPDupwindNthAnti1At33 = J11L*PDupwindNthAnti1At33 + 
        J21L*PDupwindNthAnti2At33 + J31L*PDupwindNthAnti3At33;
      
      JacPDupwindNthAnti1B1 = J11L*PDupwindNthAnti1B1 + 
        J21L*PDupwindNthAnti2B1 + J31L*PDupwindNthAnti3B1;
      
      JacPDupwindNthAnti1B2 = J11L*PDupwindNthAnti1B2 + 
        J21L*PDupwindNthAnti2B2 + J31L*PDupwindNthAnti3B2;
      
      JacPDupwindNthAnti1B3 = J11L*PDupwindNthAnti1B3 + 
        J21L*PDupwindNthAnti2B3 + J31L*PDupwindNthAnti3B3;
      
      JacPDupwindNthAnti1beta1 = J11L*PDupwindNthAnti1beta1 + 
        J21L*PDupwindNthAnti2beta1 + J31L*PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti1beta2 = J11L*PDupwindNthAnti1beta2 + 
        J21L*PDupwindNthAnti2beta2 + J31L*PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti1beta3 = J11L*PDupwindNthAnti1beta3 + 
        J21L*PDupwindNthAnti2beta3 + J31L*PDupwindNthAnti3beta3;
      
      JacPDupwindNthAnti1gt11 = J11L*PDupwindNthAnti1gt11 + 
        J21L*PDupwindNthAnti2gt11 + J31L*PDupwindNthAnti3gt11;
      
      JacPDupwindNthAnti1gt12 = J11L*PDupwindNthAnti1gt12 + 
        J21L*PDupwindNthAnti2gt12 + J31L*PDupwindNthAnti3gt12;
      
      JacPDupwindNthAnti1gt13 = J11L*PDupwindNthAnti1gt13 + 
        J21L*PDupwindNthAnti2gt13 + J31L*PDupwindNthAnti3gt13;
      
      JacPDupwindNthAnti1gt22 = J11L*PDupwindNthAnti1gt22 + 
        J21L*PDupwindNthAnti2gt22 + J31L*PDupwindNthAnti3gt22;
      
      JacPDupwindNthAnti1gt23 = J11L*PDupwindNthAnti1gt23 + 
        J21L*PDupwindNthAnti2gt23 + J31L*PDupwindNthAnti3gt23;
      
      JacPDupwindNthAnti1gt33 = J11L*PDupwindNthAnti1gt33 + 
        J21L*PDupwindNthAnti2gt33 + J31L*PDupwindNthAnti3gt33;
      
      JacPDupwindNthAnti1phi = J11L*PDupwindNthAnti1phi + 
        J21L*PDupwindNthAnti2phi + J31L*PDupwindNthAnti3phi;
      
      JacPDupwindNthAnti1trK = J11L*PDupwindNthAnti1trK + 
        J21L*PDupwindNthAnti2trK + J31L*PDupwindNthAnti3trK;
      
      JacPDupwindNthAnti1Xt1 = J11L*PDupwindNthAnti1Xt1 + 
        J21L*PDupwindNthAnti2Xt1 + J31L*PDupwindNthAnti3Xt1;
      
      JacPDupwindNthAnti1Xt2 = J11L*PDupwindNthAnti1Xt2 + 
        J21L*PDupwindNthAnti2Xt2 + J31L*PDupwindNthAnti3Xt2;
      
      JacPDupwindNthAnti1Xt3 = J11L*PDupwindNthAnti1Xt3 + 
        J21L*PDupwindNthAnti2Xt3 + J31L*PDupwindNthAnti3Xt3;
      
      JacPDupwindNthSymm1A = J11L*PDupwindNthSymm1A + J21L*PDupwindNthSymm2A 
        + J31L*PDupwindNthSymm3A;
      
      JacPDupwindNthSymm1alpha = J11L*PDupwindNthSymm1alpha + 
        J21L*PDupwindNthSymm2alpha + J31L*PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm1At11 = J11L*PDupwindNthSymm1At11 + 
        J21L*PDupwindNthSymm2At11 + J31L*PDupwindNthSymm3At11;
      
      JacPDupwindNthSymm1At12 = J11L*PDupwindNthSymm1At12 + 
        J21L*PDupwindNthSymm2At12 + J31L*PDupwindNthSymm3At12;
      
      JacPDupwindNthSymm1At13 = J11L*PDupwindNthSymm1At13 + 
        J21L*PDupwindNthSymm2At13 + J31L*PDupwindNthSymm3At13;
      
      JacPDupwindNthSymm1At22 = J11L*PDupwindNthSymm1At22 + 
        J21L*PDupwindNthSymm2At22 + J31L*PDupwindNthSymm3At22;
      
      JacPDupwindNthSymm1At23 = J11L*PDupwindNthSymm1At23 + 
        J21L*PDupwindNthSymm2At23 + J31L*PDupwindNthSymm3At23;
      
      JacPDupwindNthSymm1At33 = J11L*PDupwindNthSymm1At33 + 
        J21L*PDupwindNthSymm2At33 + J31L*PDupwindNthSymm3At33;
      
      JacPDupwindNthSymm1B1 = J11L*PDupwindNthSymm1B1 + 
        J21L*PDupwindNthSymm2B1 + J31L*PDupwindNthSymm3B1;
      
      JacPDupwindNthSymm1B2 = J11L*PDupwindNthSymm1B2 + 
        J21L*PDupwindNthSymm2B2 + J31L*PDupwindNthSymm3B2;
      
      JacPDupwindNthSymm1B3 = J11L*PDupwindNthSymm1B3 + 
        J21L*PDupwindNthSymm2B3 + J31L*PDupwindNthSymm3B3;
      
      JacPDupwindNthSymm1beta1 = J11L*PDupwindNthSymm1beta1 + 
        J21L*PDupwindNthSymm2beta1 + J31L*PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm1beta2 = J11L*PDupwindNthSymm1beta2 + 
        J21L*PDupwindNthSymm2beta2 + J31L*PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm1beta3 = J11L*PDupwindNthSymm1beta3 + 
        J21L*PDupwindNthSymm2beta3 + J31L*PDupwindNthSymm3beta3;
      
      JacPDupwindNthSymm1gt11 = J11L*PDupwindNthSymm1gt11 + 
        J21L*PDupwindNthSymm2gt11 + J31L*PDupwindNthSymm3gt11;
      
      JacPDupwindNthSymm1gt12 = J11L*PDupwindNthSymm1gt12 + 
        J21L*PDupwindNthSymm2gt12 + J31L*PDupwindNthSymm3gt12;
      
      JacPDupwindNthSymm1gt13 = J11L*PDupwindNthSymm1gt13 + 
        J21L*PDupwindNthSymm2gt13 + J31L*PDupwindNthSymm3gt13;
      
      JacPDupwindNthSymm1gt22 = J11L*PDupwindNthSymm1gt22 + 
        J21L*PDupwindNthSymm2gt22 + J31L*PDupwindNthSymm3gt22;
      
      JacPDupwindNthSymm1gt23 = J11L*PDupwindNthSymm1gt23 + 
        J21L*PDupwindNthSymm2gt23 + J31L*PDupwindNthSymm3gt23;
      
      JacPDupwindNthSymm1gt33 = J11L*PDupwindNthSymm1gt33 + 
        J21L*PDupwindNthSymm2gt33 + J31L*PDupwindNthSymm3gt33;
      
      JacPDupwindNthSymm1phi = J11L*PDupwindNthSymm1phi + 
        J21L*PDupwindNthSymm2phi + J31L*PDupwindNthSymm3phi;
      
      JacPDupwindNthSymm1trK = J11L*PDupwindNthSymm1trK + 
        J21L*PDupwindNthSymm2trK + J31L*PDupwindNthSymm3trK;
      
      JacPDupwindNthSymm1Xt1 = J11L*PDupwindNthSymm1Xt1 + 
        J21L*PDupwindNthSymm2Xt1 + J31L*PDupwindNthSymm3Xt1;
      
      JacPDupwindNthSymm1Xt2 = J11L*PDupwindNthSymm1Xt2 + 
        J21L*PDupwindNthSymm2Xt2 + J31L*PDupwindNthSymm3Xt2;
      
      JacPDupwindNthSymm1Xt3 = J11L*PDupwindNthSymm1Xt3 + 
        J21L*PDupwindNthSymm2Xt3 + J31L*PDupwindNthSymm3Xt3;
      
      JacPDupwindNthAnti2A = J12L*PDupwindNthAnti1A + J22L*PDupwindNthAnti2A 
        + J32L*PDupwindNthAnti3A;
      
      JacPDupwindNthAnti2alpha = J12L*PDupwindNthAnti1alpha + 
        J22L*PDupwindNthAnti2alpha + J32L*PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti2At11 = J12L*PDupwindNthAnti1At11 + 
        J22L*PDupwindNthAnti2At11 + J32L*PDupwindNthAnti3At11;
      
      JacPDupwindNthAnti2At12 = J12L*PDupwindNthAnti1At12 + 
        J22L*PDupwindNthAnti2At12 + J32L*PDupwindNthAnti3At12;
      
      JacPDupwindNthAnti2At13 = J12L*PDupwindNthAnti1At13 + 
        J22L*PDupwindNthAnti2At13 + J32L*PDupwindNthAnti3At13;
      
      JacPDupwindNthAnti2At22 = J12L*PDupwindNthAnti1At22 + 
        J22L*PDupwindNthAnti2At22 + J32L*PDupwindNthAnti3At22;
      
      JacPDupwindNthAnti2At23 = J12L*PDupwindNthAnti1At23 + 
        J22L*PDupwindNthAnti2At23 + J32L*PDupwindNthAnti3At23;
      
      JacPDupwindNthAnti2At33 = J12L*PDupwindNthAnti1At33 + 
        J22L*PDupwindNthAnti2At33 + J32L*PDupwindNthAnti3At33;
      
      JacPDupwindNthAnti2B1 = J12L*PDupwindNthAnti1B1 + 
        J22L*PDupwindNthAnti2B1 + J32L*PDupwindNthAnti3B1;
      
      JacPDupwindNthAnti2B2 = J12L*PDupwindNthAnti1B2 + 
        J22L*PDupwindNthAnti2B2 + J32L*PDupwindNthAnti3B2;
      
      JacPDupwindNthAnti2B3 = J12L*PDupwindNthAnti1B3 + 
        J22L*PDupwindNthAnti2B3 + J32L*PDupwindNthAnti3B3;
      
      JacPDupwindNthAnti2beta1 = J12L*PDupwindNthAnti1beta1 + 
        J22L*PDupwindNthAnti2beta1 + J32L*PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti2beta2 = J12L*PDupwindNthAnti1beta2 + 
        J22L*PDupwindNthAnti2beta2 + J32L*PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti2beta3 = J12L*PDupwindNthAnti1beta3 + 
        J22L*PDupwindNthAnti2beta3 + J32L*PDupwindNthAnti3beta3;
      
      JacPDupwindNthAnti2gt11 = J12L*PDupwindNthAnti1gt11 + 
        J22L*PDupwindNthAnti2gt11 + J32L*PDupwindNthAnti3gt11;
      
      JacPDupwindNthAnti2gt12 = J12L*PDupwindNthAnti1gt12 + 
        J22L*PDupwindNthAnti2gt12 + J32L*PDupwindNthAnti3gt12;
      
      JacPDupwindNthAnti2gt13 = J12L*PDupwindNthAnti1gt13 + 
        J22L*PDupwindNthAnti2gt13 + J32L*PDupwindNthAnti3gt13;
      
      JacPDupwindNthAnti2gt22 = J12L*PDupwindNthAnti1gt22 + 
        J22L*PDupwindNthAnti2gt22 + J32L*PDupwindNthAnti3gt22;
      
      JacPDupwindNthAnti2gt23 = J12L*PDupwindNthAnti1gt23 + 
        J22L*PDupwindNthAnti2gt23 + J32L*PDupwindNthAnti3gt23;
      
      JacPDupwindNthAnti2gt33 = J12L*PDupwindNthAnti1gt33 + 
        J22L*PDupwindNthAnti2gt33 + J32L*PDupwindNthAnti3gt33;
      
      JacPDupwindNthAnti2phi = J12L*PDupwindNthAnti1phi + 
        J22L*PDupwindNthAnti2phi + J32L*PDupwindNthAnti3phi;
      
      JacPDupwindNthAnti2trK = J12L*PDupwindNthAnti1trK + 
        J22L*PDupwindNthAnti2trK + J32L*PDupwindNthAnti3trK;
      
      JacPDupwindNthAnti2Xt1 = J12L*PDupwindNthAnti1Xt1 + 
        J22L*PDupwindNthAnti2Xt1 + J32L*PDupwindNthAnti3Xt1;
      
      JacPDupwindNthAnti2Xt2 = J12L*PDupwindNthAnti1Xt2 + 
        J22L*PDupwindNthAnti2Xt2 + J32L*PDupwindNthAnti3Xt2;
      
      JacPDupwindNthAnti2Xt3 = J12L*PDupwindNthAnti1Xt3 + 
        J22L*PDupwindNthAnti2Xt3 + J32L*PDupwindNthAnti3Xt3;
      
      JacPDupwindNthSymm2A = J12L*PDupwindNthSymm1A + J22L*PDupwindNthSymm2A 
        + J32L*PDupwindNthSymm3A;
      
      JacPDupwindNthSymm2alpha = J12L*PDupwindNthSymm1alpha + 
        J22L*PDupwindNthSymm2alpha + J32L*PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm2At11 = J12L*PDupwindNthSymm1At11 + 
        J22L*PDupwindNthSymm2At11 + J32L*PDupwindNthSymm3At11;
      
      JacPDupwindNthSymm2At12 = J12L*PDupwindNthSymm1At12 + 
        J22L*PDupwindNthSymm2At12 + J32L*PDupwindNthSymm3At12;
      
      JacPDupwindNthSymm2At13 = J12L*PDupwindNthSymm1At13 + 
        J22L*PDupwindNthSymm2At13 + J32L*PDupwindNthSymm3At13;
      
      JacPDupwindNthSymm2At22 = J12L*PDupwindNthSymm1At22 + 
        J22L*PDupwindNthSymm2At22 + J32L*PDupwindNthSymm3At22;
      
      JacPDupwindNthSymm2At23 = J12L*PDupwindNthSymm1At23 + 
        J22L*PDupwindNthSymm2At23 + J32L*PDupwindNthSymm3At23;
      
      JacPDupwindNthSymm2At33 = J12L*PDupwindNthSymm1At33 + 
        J22L*PDupwindNthSymm2At33 + J32L*PDupwindNthSymm3At33;
      
      JacPDupwindNthSymm2B1 = J12L*PDupwindNthSymm1B1 + 
        J22L*PDupwindNthSymm2B1 + J32L*PDupwindNthSymm3B1;
      
      JacPDupwindNthSymm2B2 = J12L*PDupwindNthSymm1B2 + 
        J22L*PDupwindNthSymm2B2 + J32L*PDupwindNthSymm3B2;
      
      JacPDupwindNthSymm2B3 = J12L*PDupwindNthSymm1B3 + 
        J22L*PDupwindNthSymm2B3 + J32L*PDupwindNthSymm3B3;
      
      JacPDupwindNthSymm2beta1 = J12L*PDupwindNthSymm1beta1 + 
        J22L*PDupwindNthSymm2beta1 + J32L*PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm2beta2 = J12L*PDupwindNthSymm1beta2 + 
        J22L*PDupwindNthSymm2beta2 + J32L*PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm2beta3 = J12L*PDupwindNthSymm1beta3 + 
        J22L*PDupwindNthSymm2beta3 + J32L*PDupwindNthSymm3beta3;
      
      JacPDupwindNthSymm2gt11 = J12L*PDupwindNthSymm1gt11 + 
        J22L*PDupwindNthSymm2gt11 + J32L*PDupwindNthSymm3gt11;
      
      JacPDupwindNthSymm2gt12 = J12L*PDupwindNthSymm1gt12 + 
        J22L*PDupwindNthSymm2gt12 + J32L*PDupwindNthSymm3gt12;
      
      JacPDupwindNthSymm2gt13 = J12L*PDupwindNthSymm1gt13 + 
        J22L*PDupwindNthSymm2gt13 + J32L*PDupwindNthSymm3gt13;
      
      JacPDupwindNthSymm2gt22 = J12L*PDupwindNthSymm1gt22 + 
        J22L*PDupwindNthSymm2gt22 + J32L*PDupwindNthSymm3gt22;
      
      JacPDupwindNthSymm2gt23 = J12L*PDupwindNthSymm1gt23 + 
        J22L*PDupwindNthSymm2gt23 + J32L*PDupwindNthSymm3gt23;
      
      JacPDupwindNthSymm2gt33 = J12L*PDupwindNthSymm1gt33 + 
        J22L*PDupwindNthSymm2gt33 + J32L*PDupwindNthSymm3gt33;
      
      JacPDupwindNthSymm2phi = J12L*PDupwindNthSymm1phi + 
        J22L*PDupwindNthSymm2phi + J32L*PDupwindNthSymm3phi;
      
      JacPDupwindNthSymm2trK = J12L*PDupwindNthSymm1trK + 
        J22L*PDupwindNthSymm2trK + J32L*PDupwindNthSymm3trK;
      
      JacPDupwindNthSymm2Xt1 = J12L*PDupwindNthSymm1Xt1 + 
        J22L*PDupwindNthSymm2Xt1 + J32L*PDupwindNthSymm3Xt1;
      
      JacPDupwindNthSymm2Xt2 = J12L*PDupwindNthSymm1Xt2 + 
        J22L*PDupwindNthSymm2Xt2 + J32L*PDupwindNthSymm3Xt2;
      
      JacPDupwindNthSymm2Xt3 = J12L*PDupwindNthSymm1Xt3 + 
        J22L*PDupwindNthSymm2Xt3 + J32L*PDupwindNthSymm3Xt3;
      
      JacPDupwindNthAnti3A = J13L*PDupwindNthAnti1A + J23L*PDupwindNthAnti2A 
        + J33L*PDupwindNthAnti3A;
      
      JacPDupwindNthAnti3alpha = J13L*PDupwindNthAnti1alpha + 
        J23L*PDupwindNthAnti2alpha + J33L*PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti3At11 = J13L*PDupwindNthAnti1At11 + 
        J23L*PDupwindNthAnti2At11 + J33L*PDupwindNthAnti3At11;
      
      JacPDupwindNthAnti3At12 = J13L*PDupwindNthAnti1At12 + 
        J23L*PDupwindNthAnti2At12 + J33L*PDupwindNthAnti3At12;
      
      JacPDupwindNthAnti3At13 = J13L*PDupwindNthAnti1At13 + 
        J23L*PDupwindNthAnti2At13 + J33L*PDupwindNthAnti3At13;
      
      JacPDupwindNthAnti3At22 = J13L*PDupwindNthAnti1At22 + 
        J23L*PDupwindNthAnti2At22 + J33L*PDupwindNthAnti3At22;
      
      JacPDupwindNthAnti3At23 = J13L*PDupwindNthAnti1At23 + 
        J23L*PDupwindNthAnti2At23 + J33L*PDupwindNthAnti3At23;
      
      JacPDupwindNthAnti3At33 = J13L*PDupwindNthAnti1At33 + 
        J23L*PDupwindNthAnti2At33 + J33L*PDupwindNthAnti3At33;
      
      JacPDupwindNthAnti3B1 = J13L*PDupwindNthAnti1B1 + 
        J23L*PDupwindNthAnti2B1 + J33L*PDupwindNthAnti3B1;
      
      JacPDupwindNthAnti3B2 = J13L*PDupwindNthAnti1B2 + 
        J23L*PDupwindNthAnti2B2 + J33L*PDupwindNthAnti3B2;
      
      JacPDupwindNthAnti3B3 = J13L*PDupwindNthAnti1B3 + 
        J23L*PDupwindNthAnti2B3 + J33L*PDupwindNthAnti3B3;
      
      JacPDupwindNthAnti3beta1 = J13L*PDupwindNthAnti1beta1 + 
        J23L*PDupwindNthAnti2beta1 + J33L*PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti3beta2 = J13L*PDupwindNthAnti1beta2 + 
        J23L*PDupwindNthAnti2beta2 + J33L*PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti3beta3 = J13L*PDupwindNthAnti1beta3 + 
        J23L*PDupwindNthAnti2beta3 + J33L*PDupwindNthAnti3beta3;
      
      JacPDupwindNthAnti3gt11 = J13L*PDupwindNthAnti1gt11 + 
        J23L*PDupwindNthAnti2gt11 + J33L*PDupwindNthAnti3gt11;
      
      JacPDupwindNthAnti3gt12 = J13L*PDupwindNthAnti1gt12 + 
        J23L*PDupwindNthAnti2gt12 + J33L*PDupwindNthAnti3gt12;
      
      JacPDupwindNthAnti3gt13 = J13L*PDupwindNthAnti1gt13 + 
        J23L*PDupwindNthAnti2gt13 + J33L*PDupwindNthAnti3gt13;
      
      JacPDupwindNthAnti3gt22 = J13L*PDupwindNthAnti1gt22 + 
        J23L*PDupwindNthAnti2gt22 + J33L*PDupwindNthAnti3gt22;
      
      JacPDupwindNthAnti3gt23 = J13L*PDupwindNthAnti1gt23 + 
        J23L*PDupwindNthAnti2gt23 + J33L*PDupwindNthAnti3gt23;
      
      JacPDupwindNthAnti3gt33 = J13L*PDupwindNthAnti1gt33 + 
        J23L*PDupwindNthAnti2gt33 + J33L*PDupwindNthAnti3gt33;
      
      JacPDupwindNthAnti3phi = J13L*PDupwindNthAnti1phi + 
        J23L*PDupwindNthAnti2phi + J33L*PDupwindNthAnti3phi;
      
      JacPDupwindNthAnti3trK = J13L*PDupwindNthAnti1trK + 
        J23L*PDupwindNthAnti2trK + J33L*PDupwindNthAnti3trK;
      
      JacPDupwindNthAnti3Xt1 = J13L*PDupwindNthAnti1Xt1 + 
        J23L*PDupwindNthAnti2Xt1 + J33L*PDupwindNthAnti3Xt1;
      
      JacPDupwindNthAnti3Xt2 = J13L*PDupwindNthAnti1Xt2 + 
        J23L*PDupwindNthAnti2Xt2 + J33L*PDupwindNthAnti3Xt2;
      
      JacPDupwindNthAnti3Xt3 = J13L*PDupwindNthAnti1Xt3 + 
        J23L*PDupwindNthAnti2Xt3 + J33L*PDupwindNthAnti3Xt3;
      
      JacPDupwindNthSymm3A = J13L*PDupwindNthSymm1A + J23L*PDupwindNthSymm2A 
        + J33L*PDupwindNthSymm3A;
      
      JacPDupwindNthSymm3alpha = J13L*PDupwindNthSymm1alpha + 
        J23L*PDupwindNthSymm2alpha + J33L*PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm3At11 = J13L*PDupwindNthSymm1At11 + 
        J23L*PDupwindNthSymm2At11 + J33L*PDupwindNthSymm3At11;
      
      JacPDupwindNthSymm3At12 = J13L*PDupwindNthSymm1At12 + 
        J23L*PDupwindNthSymm2At12 + J33L*PDupwindNthSymm3At12;
      
      JacPDupwindNthSymm3At13 = J13L*PDupwindNthSymm1At13 + 
        J23L*PDupwindNthSymm2At13 + J33L*PDupwindNthSymm3At13;
      
      JacPDupwindNthSymm3At22 = J13L*PDupwindNthSymm1At22 + 
        J23L*PDupwindNthSymm2At22 + J33L*PDupwindNthSymm3At22;
      
      JacPDupwindNthSymm3At23 = J13L*PDupwindNthSymm1At23 + 
        J23L*PDupwindNthSymm2At23 + J33L*PDupwindNthSymm3At23;
      
      JacPDupwindNthSymm3At33 = J13L*PDupwindNthSymm1At33 + 
        J23L*PDupwindNthSymm2At33 + J33L*PDupwindNthSymm3At33;
      
      JacPDupwindNthSymm3B1 = J13L*PDupwindNthSymm1B1 + 
        J23L*PDupwindNthSymm2B1 + J33L*PDupwindNthSymm3B1;
      
      JacPDupwindNthSymm3B2 = J13L*PDupwindNthSymm1B2 + 
        J23L*PDupwindNthSymm2B2 + J33L*PDupwindNthSymm3B2;
      
      JacPDupwindNthSymm3B3 = J13L*PDupwindNthSymm1B3 + 
        J23L*PDupwindNthSymm2B3 + J33L*PDupwindNthSymm3B3;
      
      JacPDupwindNthSymm3beta1 = J13L*PDupwindNthSymm1beta1 + 
        J23L*PDupwindNthSymm2beta1 + J33L*PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm3beta2 = J13L*PDupwindNthSymm1beta2 + 
        J23L*PDupwindNthSymm2beta2 + J33L*PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm3beta3 = J13L*PDupwindNthSymm1beta3 + 
        J23L*PDupwindNthSymm2beta3 + J33L*PDupwindNthSymm3beta3;
      
      JacPDupwindNthSymm3gt11 = J13L*PDupwindNthSymm1gt11 + 
        J23L*PDupwindNthSymm2gt11 + J33L*PDupwindNthSymm3gt11;
      
      JacPDupwindNthSymm3gt12 = J13L*PDupwindNthSymm1gt12 + 
        J23L*PDupwindNthSymm2gt12 + J33L*PDupwindNthSymm3gt12;
      
      JacPDupwindNthSymm3gt13 = J13L*PDupwindNthSymm1gt13 + 
        J23L*PDupwindNthSymm2gt13 + J33L*PDupwindNthSymm3gt13;
      
      JacPDupwindNthSymm3gt22 = J13L*PDupwindNthSymm1gt22 + 
        J23L*PDupwindNthSymm2gt22 + J33L*PDupwindNthSymm3gt22;
      
      JacPDupwindNthSymm3gt23 = J13L*PDupwindNthSymm1gt23 + 
        J23L*PDupwindNthSymm2gt23 + J33L*PDupwindNthSymm3gt23;
      
      JacPDupwindNthSymm3gt33 = J13L*PDupwindNthSymm1gt33 + 
        J23L*PDupwindNthSymm2gt33 + J33L*PDupwindNthSymm3gt33;
      
      JacPDupwindNthSymm3phi = J13L*PDupwindNthSymm1phi + 
        J23L*PDupwindNthSymm2phi + J33L*PDupwindNthSymm3phi;
      
      JacPDupwindNthSymm3trK = J13L*PDupwindNthSymm1trK + 
        J23L*PDupwindNthSymm2trK + J33L*PDupwindNthSymm3trK;
      
      JacPDupwindNthSymm3Xt1 = J13L*PDupwindNthSymm1Xt1 + 
        J23L*PDupwindNthSymm2Xt1 + J33L*PDupwindNthSymm3Xt1;
      
      JacPDupwindNthSymm3Xt2 = J13L*PDupwindNthSymm1Xt2 + 
        J23L*PDupwindNthSymm2Xt2 + J33L*PDupwindNthSymm3Xt2;
      
      JacPDupwindNthSymm3Xt3 = J13L*PDupwindNthSymm1Xt3 + 
        J23L*PDupwindNthSymm2Xt3 + J33L*PDupwindNthSymm3Xt3;
    }
    else
    {
      JacPDupwindNthAnti1A = PDupwindNthAnti1A;
      
      JacPDupwindNthAnti1alpha = PDupwindNthAnti1alpha;
      
      JacPDupwindNthAnti1At11 = PDupwindNthAnti1At11;
      
      JacPDupwindNthAnti1At12 = PDupwindNthAnti1At12;
      
      JacPDupwindNthAnti1At13 = PDupwindNthAnti1At13;
      
      JacPDupwindNthAnti1At22 = PDupwindNthAnti1At22;
      
      JacPDupwindNthAnti1At23 = PDupwindNthAnti1At23;
      
      JacPDupwindNthAnti1At33 = PDupwindNthAnti1At33;
      
      JacPDupwindNthAnti1B1 = PDupwindNthAnti1B1;
      
      JacPDupwindNthAnti1B2 = PDupwindNthAnti1B2;
      
      JacPDupwindNthAnti1B3 = PDupwindNthAnti1B3;
      
      JacPDupwindNthAnti1beta1 = PDupwindNthAnti1beta1;
      
      JacPDupwindNthAnti1beta2 = PDupwindNthAnti1beta2;
      
      JacPDupwindNthAnti1beta3 = PDupwindNthAnti1beta3;
      
      JacPDupwindNthAnti1gt11 = PDupwindNthAnti1gt11;
      
      JacPDupwindNthAnti1gt12 = PDupwindNthAnti1gt12;
      
      JacPDupwindNthAnti1gt13 = PDupwindNthAnti1gt13;
      
      JacPDupwindNthAnti1gt22 = PDupwindNthAnti1gt22;
      
      JacPDupwindNthAnti1gt23 = PDupwindNthAnti1gt23;
      
      JacPDupwindNthAnti1gt33 = PDupwindNthAnti1gt33;
      
      JacPDupwindNthAnti1phi = PDupwindNthAnti1phi;
      
      JacPDupwindNthAnti1trK = PDupwindNthAnti1trK;
      
      JacPDupwindNthAnti1Xt1 = PDupwindNthAnti1Xt1;
      
      JacPDupwindNthAnti1Xt2 = PDupwindNthAnti1Xt2;
      
      JacPDupwindNthAnti1Xt3 = PDupwindNthAnti1Xt3;
      
      JacPDupwindNthSymm1A = PDupwindNthSymm1A;
      
      JacPDupwindNthSymm1alpha = PDupwindNthSymm1alpha;
      
      JacPDupwindNthSymm1At11 = PDupwindNthSymm1At11;
      
      JacPDupwindNthSymm1At12 = PDupwindNthSymm1At12;
      
      JacPDupwindNthSymm1At13 = PDupwindNthSymm1At13;
      
      JacPDupwindNthSymm1At22 = PDupwindNthSymm1At22;
      
      JacPDupwindNthSymm1At23 = PDupwindNthSymm1At23;
      
      JacPDupwindNthSymm1At33 = PDupwindNthSymm1At33;
      
      JacPDupwindNthSymm1B1 = PDupwindNthSymm1B1;
      
      JacPDupwindNthSymm1B2 = PDupwindNthSymm1B2;
      
      JacPDupwindNthSymm1B3 = PDupwindNthSymm1B3;
      
      JacPDupwindNthSymm1beta1 = PDupwindNthSymm1beta1;
      
      JacPDupwindNthSymm1beta2 = PDupwindNthSymm1beta2;
      
      JacPDupwindNthSymm1beta3 = PDupwindNthSymm1beta3;
      
      JacPDupwindNthSymm1gt11 = PDupwindNthSymm1gt11;
      
      JacPDupwindNthSymm1gt12 = PDupwindNthSymm1gt12;
      
      JacPDupwindNthSymm1gt13 = PDupwindNthSymm1gt13;
      
      JacPDupwindNthSymm1gt22 = PDupwindNthSymm1gt22;
      
      JacPDupwindNthSymm1gt23 = PDupwindNthSymm1gt23;
      
      JacPDupwindNthSymm1gt33 = PDupwindNthSymm1gt33;
      
      JacPDupwindNthSymm1phi = PDupwindNthSymm1phi;
      
      JacPDupwindNthSymm1trK = PDupwindNthSymm1trK;
      
      JacPDupwindNthSymm1Xt1 = PDupwindNthSymm1Xt1;
      
      JacPDupwindNthSymm1Xt2 = PDupwindNthSymm1Xt2;
      
      JacPDupwindNthSymm1Xt3 = PDupwindNthSymm1Xt3;
      
      JacPDupwindNthAnti2A = PDupwindNthAnti2A;
      
      JacPDupwindNthAnti2alpha = PDupwindNthAnti2alpha;
      
      JacPDupwindNthAnti2At11 = PDupwindNthAnti2At11;
      
      JacPDupwindNthAnti2At12 = PDupwindNthAnti2At12;
      
      JacPDupwindNthAnti2At13 = PDupwindNthAnti2At13;
      
      JacPDupwindNthAnti2At22 = PDupwindNthAnti2At22;
      
      JacPDupwindNthAnti2At23 = PDupwindNthAnti2At23;
      
      JacPDupwindNthAnti2At33 = PDupwindNthAnti2At33;
      
      JacPDupwindNthAnti2B1 = PDupwindNthAnti2B1;
      
      JacPDupwindNthAnti2B2 = PDupwindNthAnti2B2;
      
      JacPDupwindNthAnti2B3 = PDupwindNthAnti2B3;
      
      JacPDupwindNthAnti2beta1 = PDupwindNthAnti2beta1;
      
      JacPDupwindNthAnti2beta2 = PDupwindNthAnti2beta2;
      
      JacPDupwindNthAnti2beta3 = PDupwindNthAnti2beta3;
      
      JacPDupwindNthAnti2gt11 = PDupwindNthAnti2gt11;
      
      JacPDupwindNthAnti2gt12 = PDupwindNthAnti2gt12;
      
      JacPDupwindNthAnti2gt13 = PDupwindNthAnti2gt13;
      
      JacPDupwindNthAnti2gt22 = PDupwindNthAnti2gt22;
      
      JacPDupwindNthAnti2gt23 = PDupwindNthAnti2gt23;
      
      JacPDupwindNthAnti2gt33 = PDupwindNthAnti2gt33;
      
      JacPDupwindNthAnti2phi = PDupwindNthAnti2phi;
      
      JacPDupwindNthAnti2trK = PDupwindNthAnti2trK;
      
      JacPDupwindNthAnti2Xt1 = PDupwindNthAnti2Xt1;
      
      JacPDupwindNthAnti2Xt2 = PDupwindNthAnti2Xt2;
      
      JacPDupwindNthAnti2Xt3 = PDupwindNthAnti2Xt3;
      
      JacPDupwindNthSymm2A = PDupwindNthSymm2A;
      
      JacPDupwindNthSymm2alpha = PDupwindNthSymm2alpha;
      
      JacPDupwindNthSymm2At11 = PDupwindNthSymm2At11;
      
      JacPDupwindNthSymm2At12 = PDupwindNthSymm2At12;
      
      JacPDupwindNthSymm2At13 = PDupwindNthSymm2At13;
      
      JacPDupwindNthSymm2At22 = PDupwindNthSymm2At22;
      
      JacPDupwindNthSymm2At23 = PDupwindNthSymm2At23;
      
      JacPDupwindNthSymm2At33 = PDupwindNthSymm2At33;
      
      JacPDupwindNthSymm2B1 = PDupwindNthSymm2B1;
      
      JacPDupwindNthSymm2B2 = PDupwindNthSymm2B2;
      
      JacPDupwindNthSymm2B3 = PDupwindNthSymm2B3;
      
      JacPDupwindNthSymm2beta1 = PDupwindNthSymm2beta1;
      
      JacPDupwindNthSymm2beta2 = PDupwindNthSymm2beta2;
      
      JacPDupwindNthSymm2beta3 = PDupwindNthSymm2beta3;
      
      JacPDupwindNthSymm2gt11 = PDupwindNthSymm2gt11;
      
      JacPDupwindNthSymm2gt12 = PDupwindNthSymm2gt12;
      
      JacPDupwindNthSymm2gt13 = PDupwindNthSymm2gt13;
      
      JacPDupwindNthSymm2gt22 = PDupwindNthSymm2gt22;
      
      JacPDupwindNthSymm2gt23 = PDupwindNthSymm2gt23;
      
      JacPDupwindNthSymm2gt33 = PDupwindNthSymm2gt33;
      
      JacPDupwindNthSymm2phi = PDupwindNthSymm2phi;
      
      JacPDupwindNthSymm2trK = PDupwindNthSymm2trK;
      
      JacPDupwindNthSymm2Xt1 = PDupwindNthSymm2Xt1;
      
      JacPDupwindNthSymm2Xt2 = PDupwindNthSymm2Xt2;
      
      JacPDupwindNthSymm2Xt3 = PDupwindNthSymm2Xt3;
      
      JacPDupwindNthAnti3A = PDupwindNthAnti3A;
      
      JacPDupwindNthAnti3alpha = PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti3At11 = PDupwindNthAnti3At11;
      
      JacPDupwindNthAnti3At12 = PDupwindNthAnti3At12;
      
      JacPDupwindNthAnti3At13 = PDupwindNthAnti3At13;
      
      JacPDupwindNthAnti3At22 = PDupwindNthAnti3At22;
      
      JacPDupwindNthAnti3At23 = PDupwindNthAnti3At23;
      
      JacPDupwindNthAnti3At33 = PDupwindNthAnti3At33;
      
      JacPDupwindNthAnti3B1 = PDupwindNthAnti3B1;
      
      JacPDupwindNthAnti3B2 = PDupwindNthAnti3B2;
      
      JacPDupwindNthAnti3B3 = PDupwindNthAnti3B3;
      
      JacPDupwindNthAnti3beta1 = PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti3beta2 = PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti3beta3 = PDupwindNthAnti3beta3;
      
      JacPDupwindNthAnti3gt11 = PDupwindNthAnti3gt11;
      
      JacPDupwindNthAnti3gt12 = PDupwindNthAnti3gt12;
      
      JacPDupwindNthAnti3gt13 = PDupwindNthAnti3gt13;
      
      JacPDupwindNthAnti3gt22 = PDupwindNthAnti3gt22;
      
      JacPDupwindNthAnti3gt23 = PDupwindNthAnti3gt23;
      
      JacPDupwindNthAnti3gt33 = PDupwindNthAnti3gt33;
      
      JacPDupwindNthAnti3phi = PDupwindNthAnti3phi;
      
      JacPDupwindNthAnti3trK = PDupwindNthAnti3trK;
      
      JacPDupwindNthAnti3Xt1 = PDupwindNthAnti3Xt1;
      
      JacPDupwindNthAnti3Xt2 = PDupwindNthAnti3Xt2;
      
      JacPDupwindNthAnti3Xt3 = PDupwindNthAnti3Xt3;
      
      JacPDupwindNthSymm3A = PDupwindNthSymm3A;
      
      JacPDupwindNthSymm3alpha = PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm3At11 = PDupwindNthSymm3At11;
      
      JacPDupwindNthSymm3At12 = PDupwindNthSymm3At12;
      
      JacPDupwindNthSymm3At13 = PDupwindNthSymm3At13;
      
      JacPDupwindNthSymm3At22 = PDupwindNthSymm3At22;
      
      JacPDupwindNthSymm3At23 = PDupwindNthSymm3At23;
      
      JacPDupwindNthSymm3At33 = PDupwindNthSymm3At33;
      
      JacPDupwindNthSymm3B1 = PDupwindNthSymm3B1;
      
      JacPDupwindNthSymm3B2 = PDupwindNthSymm3B2;
      
      JacPDupwindNthSymm3B3 = PDupwindNthSymm3B3;
      
      JacPDupwindNthSymm3beta1 = PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm3beta2 = PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm3beta3 = PDupwindNthSymm3beta3;
      
      JacPDupwindNthSymm3gt11 = PDupwindNthSymm3gt11;
      
      JacPDupwindNthSymm3gt12 = PDupwindNthSymm3gt12;
      
      JacPDupwindNthSymm3gt13 = PDupwindNthSymm3gt13;
      
      JacPDupwindNthSymm3gt22 = PDupwindNthSymm3gt22;
      
      JacPDupwindNthSymm3gt23 = PDupwindNthSymm3gt23;
      
      JacPDupwindNthSymm3gt33 = PDupwindNthSymm3gt33;
      
      JacPDupwindNthSymm3phi = PDupwindNthSymm3phi;
      
      JacPDupwindNthSymm3trK = PDupwindNthSymm3trK;
      
      JacPDupwindNthSymm3Xt1 = PDupwindNthSymm3Xt1;
      
      JacPDupwindNthSymm3Xt2 = PDupwindNthSymm3Xt2;
      
      JacPDupwindNthSymm3Xt3 = PDupwindNthSymm3Xt3;
    }
    
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    phirhsL = beta1L*JacPDupwindNthAnti1phi + 
      beta2L*JacPDupwindNthAnti2phi + beta3L*JacPDupwindNthAnti3phi + phirhsL 
      + JacPDupwindNthSymm1phi*Abs(beta1L) + 
      JacPDupwindNthSymm2phi*Abs(beta2L) + 
      JacPDupwindNthSymm3phi*Abs(beta3L);
    
    gt11rhsL = gt11rhsL + beta1L*JacPDupwindNthAnti1gt11 + 
      beta2L*JacPDupwindNthAnti2gt11 + beta3L*JacPDupwindNthAnti3gt11 + 
      JacPDupwindNthSymm1gt11*Abs(beta1L) + 
      JacPDupwindNthSymm2gt11*Abs(beta2L) + 
      JacPDupwindNthSymm3gt11*Abs(beta3L);
    
    gt12rhsL = gt12rhsL + beta1L*JacPDupwindNthAnti1gt12 + 
      beta2L*JacPDupwindNthAnti2gt12 + beta3L*JacPDupwindNthAnti3gt12 + 
      JacPDupwindNthSymm1gt12*Abs(beta1L) + 
      JacPDupwindNthSymm2gt12*Abs(beta2L) + 
      JacPDupwindNthSymm3gt12*Abs(beta3L);
    
    gt13rhsL = gt13rhsL + beta1L*JacPDupwindNthAnti1gt13 + 
      beta2L*JacPDupwindNthAnti2gt13 + beta3L*JacPDupwindNthAnti3gt13 + 
      JacPDupwindNthSymm1gt13*Abs(beta1L) + 
      JacPDupwindNthSymm2gt13*Abs(beta2L) + 
      JacPDupwindNthSymm3gt13*Abs(beta3L);
    
    gt22rhsL = gt22rhsL + beta1L*JacPDupwindNthAnti1gt22 + 
      beta2L*JacPDupwindNthAnti2gt22 + beta3L*JacPDupwindNthAnti3gt22 + 
      JacPDupwindNthSymm1gt22*Abs(beta1L) + 
      JacPDupwindNthSymm2gt22*Abs(beta2L) + 
      JacPDupwindNthSymm3gt22*Abs(beta3L);
    
    gt23rhsL = gt23rhsL + beta1L*JacPDupwindNthAnti1gt23 + 
      beta2L*JacPDupwindNthAnti2gt23 + beta3L*JacPDupwindNthAnti3gt23 + 
      JacPDupwindNthSymm1gt23*Abs(beta1L) + 
      JacPDupwindNthSymm2gt23*Abs(beta2L) + 
      JacPDupwindNthSymm3gt23*Abs(beta3L);
    
    gt33rhsL = gt33rhsL + beta1L*JacPDupwindNthAnti1gt33 + 
      beta2L*JacPDupwindNthAnti2gt33 + beta3L*JacPDupwindNthAnti3gt33 + 
      JacPDupwindNthSymm1gt33*Abs(beta1L) + 
      JacPDupwindNthSymm2gt33*Abs(beta2L) + 
      JacPDupwindNthSymm3gt33*Abs(beta3L);
    
    Xt1rhsL = beta1L*JacPDupwindNthAnti1Xt1 + 
      beta2L*JacPDupwindNthAnti2Xt1 + beta3L*JacPDupwindNthAnti3Xt1 + Xt1rhsL 
      + JacPDupwindNthSymm1Xt1*Abs(beta1L) + 
      JacPDupwindNthSymm2Xt1*Abs(beta2L) + 
      JacPDupwindNthSymm3Xt1*Abs(beta3L);
    
    Xt2rhsL = beta1L*JacPDupwindNthAnti1Xt2 + 
      beta2L*JacPDupwindNthAnti2Xt2 + beta3L*JacPDupwindNthAnti3Xt2 + Xt2rhsL 
      + JacPDupwindNthSymm1Xt2*Abs(beta1L) + 
      JacPDupwindNthSymm2Xt2*Abs(beta2L) + 
      JacPDupwindNthSymm3Xt2*Abs(beta3L);
    
    Xt3rhsL = beta1L*JacPDupwindNthAnti1Xt3 + 
      beta2L*JacPDupwindNthAnti2Xt3 + beta3L*JacPDupwindNthAnti3Xt3 + Xt3rhsL 
      + JacPDupwindNthSymm1Xt3*Abs(beta1L) + 
      JacPDupwindNthSymm2Xt3*Abs(beta2L) + 
      JacPDupwindNthSymm3Xt3*Abs(beta3L);
    
    trKrhsL = beta1L*JacPDupwindNthAnti1trK + 
      beta2L*JacPDupwindNthAnti2trK + beta3L*JacPDupwindNthAnti3trK + trKrhsL 
      + JacPDupwindNthSymm1trK*Abs(beta1L) + 
      JacPDupwindNthSymm2trK*Abs(beta2L) + 
      JacPDupwindNthSymm3trK*Abs(beta3L);
    
    At11rhsL = At11rhsL + beta1L*JacPDupwindNthAnti1At11 + 
      beta2L*JacPDupwindNthAnti2At11 + beta3L*JacPDupwindNthAnti3At11 + 
      JacPDupwindNthSymm1At11*Abs(beta1L) + 
      JacPDupwindNthSymm2At11*Abs(beta2L) + 
      JacPDupwindNthSymm3At11*Abs(beta3L);
    
    At12rhsL = At12rhsL + beta1L*JacPDupwindNthAnti1At12 + 
      beta2L*JacPDupwindNthAnti2At12 + beta3L*JacPDupwindNthAnti3At12 + 
      JacPDupwindNthSymm1At12*Abs(beta1L) + 
      JacPDupwindNthSymm2At12*Abs(beta2L) + 
      JacPDupwindNthSymm3At12*Abs(beta3L);
    
    At13rhsL = At13rhsL + beta1L*JacPDupwindNthAnti1At13 + 
      beta2L*JacPDupwindNthAnti2At13 + beta3L*JacPDupwindNthAnti3At13 + 
      JacPDupwindNthSymm1At13*Abs(beta1L) + 
      JacPDupwindNthSymm2At13*Abs(beta2L) + 
      JacPDupwindNthSymm3At13*Abs(beta3L);
    
    At22rhsL = At22rhsL + beta1L*JacPDupwindNthAnti1At22 + 
      beta2L*JacPDupwindNthAnti2At22 + beta3L*JacPDupwindNthAnti3At22 + 
      JacPDupwindNthSymm1At22*Abs(beta1L) + 
      JacPDupwindNthSymm2At22*Abs(beta2L) + 
      JacPDupwindNthSymm3At22*Abs(beta3L);
    
    At23rhsL = At23rhsL + beta1L*JacPDupwindNthAnti1At23 + 
      beta2L*JacPDupwindNthAnti2At23 + beta3L*JacPDupwindNthAnti3At23 + 
      JacPDupwindNthSymm1At23*Abs(beta1L) + 
      JacPDupwindNthSymm2At23*Abs(beta2L) + 
      JacPDupwindNthSymm3At23*Abs(beta3L);
    
    At33rhsL = At33rhsL + beta1L*JacPDupwindNthAnti1At33 + 
      beta2L*JacPDupwindNthAnti2At33 + beta3L*JacPDupwindNthAnti3At33 + 
      JacPDupwindNthSymm1At33*Abs(beta1L) + 
      JacPDupwindNthSymm2At33*Abs(beta2L) + 
      JacPDupwindNthSymm3At33*Abs(beta3L);
    
    alpharhsL = alpharhsL + (beta1L*JacPDupwindNthAnti1alpha + 
      beta2L*JacPDupwindNthAnti2alpha + beta3L*JacPDupwindNthAnti3alpha + 
      JacPDupwindNthSymm1alpha*Abs(beta1L) + 
      JacPDupwindNthSymm2alpha*Abs(beta2L) + 
      JacPDupwindNthSymm3alpha*Abs(beta3L))*ToReal(LapseAdvectionCoeff);
    
    ArhsL = ArhsL + (beta1L*JacPDupwindNthAnti1A + 
      beta2L*JacPDupwindNthAnti2A + beta3L*JacPDupwindNthAnti3A + 
      JacPDupwindNthSymm1A*Abs(beta1L) + JacPDupwindNthSymm2A*Abs(beta2L) + 
      JacPDupwindNthSymm3A*Abs(beta3L))*ToReal(LapseAdvectionCoeff);
    
    beta1rhsL = beta1rhsL + (beta1L*JacPDupwindNthAnti1beta1 + 
      beta2L*JacPDupwindNthAnti2beta1 + beta3L*JacPDupwindNthAnti3beta1 + 
      JacPDupwindNthSymm1beta1*Abs(beta1L) + 
      JacPDupwindNthSymm2beta1*Abs(beta2L) + 
      JacPDupwindNthSymm3beta1*Abs(beta3L))*ToReal(ShiftAdvectionCoeff);
    
    beta2rhsL = beta2rhsL + (beta1L*JacPDupwindNthAnti1beta2 + 
      beta2L*JacPDupwindNthAnti2beta2 + beta3L*JacPDupwindNthAnti3beta2 + 
      JacPDupwindNthSymm1beta2*Abs(beta1L) + 
      JacPDupwindNthSymm2beta2*Abs(beta2L) + 
      JacPDupwindNthSymm3beta2*Abs(beta3L))*ToReal(ShiftAdvectionCoeff);
    
    beta3rhsL = beta3rhsL + (beta1L*JacPDupwindNthAnti1beta3 + 
      beta2L*JacPDupwindNthAnti2beta3 + beta3L*JacPDupwindNthAnti3beta3 + 
      JacPDupwindNthSymm1beta3*Abs(beta1L) + 
      JacPDupwindNthSymm2beta3*Abs(beta2L) + 
      JacPDupwindNthSymm3beta3*Abs(beta3L))*ToReal(ShiftAdvectionCoeff);
    
    B1rhsL = B1rhsL + (beta1L*(JacPDupwindNthAnti1B1 - 
      JacPDupwindNthAnti1Xt1) + beta2L*(JacPDupwindNthAnti2B1 - 
      JacPDupwindNthAnti2Xt1) + beta3L*(JacPDupwindNthAnti3B1 - 
      JacPDupwindNthAnti3Xt1) + (JacPDupwindNthSymm1B1 - 
      JacPDupwindNthSymm1Xt1)*Abs(beta1L) + (JacPDupwindNthSymm2B1 - 
      JacPDupwindNthSymm2Xt1)*Abs(beta2L) + (JacPDupwindNthSymm3B1 - 
      JacPDupwindNthSymm3Xt1)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      (beta1L*JacPDupwindNthAnti1Xt1 + beta2L*JacPDupwindNthAnti2Xt1 + 
      beta3L*JacPDupwindNthAnti3Xt1 + JacPDupwindNthSymm1Xt1*Abs(beta1L) + 
      JacPDupwindNthSymm2Xt1*Abs(beta2L) + 
      JacPDupwindNthSymm3Xt1*Abs(beta3L))*ToReal(ShiftBCoeff);
    
    B2rhsL = B2rhsL + (beta1L*(JacPDupwindNthAnti1B2 - 
      JacPDupwindNthAnti1Xt2) + beta2L*(JacPDupwindNthAnti2B2 - 
      JacPDupwindNthAnti2Xt2) + beta3L*(JacPDupwindNthAnti3B2 - 
      JacPDupwindNthAnti3Xt2) + (JacPDupwindNthSymm1B2 - 
      JacPDupwindNthSymm1Xt2)*Abs(beta1L) + (JacPDupwindNthSymm2B2 - 
      JacPDupwindNthSymm2Xt2)*Abs(beta2L) + (JacPDupwindNthSymm3B2 - 
      JacPDupwindNthSymm3Xt2)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      (beta1L*JacPDupwindNthAnti1Xt2 + beta2L*JacPDupwindNthAnti2Xt2 + 
      beta3L*JacPDupwindNthAnti3Xt2 + JacPDupwindNthSymm1Xt2*Abs(beta1L) + 
      JacPDupwindNthSymm2Xt2*Abs(beta2L) + 
      JacPDupwindNthSymm3Xt2*Abs(beta3L))*ToReal(ShiftBCoeff);
    
    B3rhsL = B3rhsL + (beta1L*(JacPDupwindNthAnti1B3 - 
      JacPDupwindNthAnti1Xt3) + beta2L*(JacPDupwindNthAnti2B3 - 
      JacPDupwindNthAnti2Xt3) + beta3L*(JacPDupwindNthAnti3B3 - 
      JacPDupwindNthAnti3Xt3) + (JacPDupwindNthSymm1B3 - 
      JacPDupwindNthSymm1Xt3)*Abs(beta1L) + (JacPDupwindNthSymm2B3 - 
      JacPDupwindNthSymm2Xt3)*Abs(beta2L) + (JacPDupwindNthSymm3B3 - 
      JacPDupwindNthSymm3Xt3)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      (beta1L*JacPDupwindNthAnti1Xt3 + beta2L*JacPDupwindNthAnti2Xt3 + 
      beta3L*JacPDupwindNthAnti3Xt3 + JacPDupwindNthSymm1Xt3*Abs(beta1L) + 
      JacPDupwindNthSymm2Xt3*Abs(beta2L) + 
      JacPDupwindNthSymm3Xt3*Abs(beta3L))*ToReal(ShiftBCoeff);
    
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
  LC_ENDLOOP3 (ML_BSSN_MP_Advect);
}

extern "C" void ML_BSSN_MP_Advect(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_Advect_Body);
}
