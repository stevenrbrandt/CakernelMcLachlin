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

extern "C" void ML_BSSN_MP_O8_Dissipation_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_curvrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_trace_curvrhs.");
  return;
}

static void ML_BSSN_MP_O8_Dissipation_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_O8_Dissipation_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_O8_Dissipation_calc_every != ML_BSSN_MP_O8_Dissipation_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_BSSN_MP_O8::ML_curv","ML_BSSN_MP_O8::ML_curvrhs","ML_BSSN_MP_O8::ML_dtlapse","ML_BSSN_MP_O8::ML_dtlapserhs","ML_BSSN_MP_O8::ML_dtshift","ML_BSSN_MP_O8::ML_dtshiftrhs","ML_BSSN_MP_O8::ML_Gamma","ML_BSSN_MP_O8::ML_Gammarhs","ML_BSSN_MP_O8::ML_lapse","ML_BSSN_MP_O8::ML_lapserhs","ML_BSSN_MP_O8::ML_log_confac","ML_BSSN_MP_O8::ML_log_confacrhs","ML_BSSN_MP_O8::ML_metric","ML_BSSN_MP_O8::ML_metricrhs","ML_BSSN_MP_O8::ML_shift","ML_BSSN_MP_O8::ML_shiftrhs","ML_BSSN_MP_O8::ML_trace_curv","ML_BSSN_MP_O8::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_O8_Dissipation", 18, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_MP_O8_Dissipation", 5, 5, 5);
  
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
  LC_LOOP3 (ML_BSSN_MP_O8_Dissipation,
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
    CCTK_REAL const PDdissipationNth1A = PDdissipationNth1(&A[index]);
    CCTK_REAL const PDdissipationNth2A = PDdissipationNth2(&A[index]);
    CCTK_REAL const PDdissipationNth3A = PDdissipationNth3(&A[index]);
    CCTK_REAL const PDdissipationNth1alpha = PDdissipationNth1(&alpha[index]);
    CCTK_REAL const PDdissipationNth2alpha = PDdissipationNth2(&alpha[index]);
    CCTK_REAL const PDdissipationNth3alpha = PDdissipationNth3(&alpha[index]);
    CCTK_REAL const PDdissipationNth1At11 = PDdissipationNth1(&At11[index]);
    CCTK_REAL const PDdissipationNth2At11 = PDdissipationNth2(&At11[index]);
    CCTK_REAL const PDdissipationNth3At11 = PDdissipationNth3(&At11[index]);
    CCTK_REAL const PDdissipationNth1At12 = PDdissipationNth1(&At12[index]);
    CCTK_REAL const PDdissipationNth2At12 = PDdissipationNth2(&At12[index]);
    CCTK_REAL const PDdissipationNth3At12 = PDdissipationNth3(&At12[index]);
    CCTK_REAL const PDdissipationNth1At13 = PDdissipationNth1(&At13[index]);
    CCTK_REAL const PDdissipationNth2At13 = PDdissipationNth2(&At13[index]);
    CCTK_REAL const PDdissipationNth3At13 = PDdissipationNth3(&At13[index]);
    CCTK_REAL const PDdissipationNth1At22 = PDdissipationNth1(&At22[index]);
    CCTK_REAL const PDdissipationNth2At22 = PDdissipationNth2(&At22[index]);
    CCTK_REAL const PDdissipationNth3At22 = PDdissipationNth3(&At22[index]);
    CCTK_REAL const PDdissipationNth1At23 = PDdissipationNth1(&At23[index]);
    CCTK_REAL const PDdissipationNth2At23 = PDdissipationNth2(&At23[index]);
    CCTK_REAL const PDdissipationNth3At23 = PDdissipationNth3(&At23[index]);
    CCTK_REAL const PDdissipationNth1At33 = PDdissipationNth1(&At33[index]);
    CCTK_REAL const PDdissipationNth2At33 = PDdissipationNth2(&At33[index]);
    CCTK_REAL const PDdissipationNth3At33 = PDdissipationNth3(&At33[index]);
    CCTK_REAL const PDdissipationNth1B1 = PDdissipationNth1(&B1[index]);
    CCTK_REAL const PDdissipationNth2B1 = PDdissipationNth2(&B1[index]);
    CCTK_REAL const PDdissipationNth3B1 = PDdissipationNth3(&B1[index]);
    CCTK_REAL const PDdissipationNth1B2 = PDdissipationNth1(&B2[index]);
    CCTK_REAL const PDdissipationNth2B2 = PDdissipationNth2(&B2[index]);
    CCTK_REAL const PDdissipationNth3B2 = PDdissipationNth3(&B2[index]);
    CCTK_REAL const PDdissipationNth1B3 = PDdissipationNth1(&B3[index]);
    CCTK_REAL const PDdissipationNth2B3 = PDdissipationNth2(&B3[index]);
    CCTK_REAL const PDdissipationNth3B3 = PDdissipationNth3(&B3[index]);
    CCTK_REAL const PDdissipationNth1beta1 = PDdissipationNth1(&beta1[index]);
    CCTK_REAL const PDdissipationNth2beta1 = PDdissipationNth2(&beta1[index]);
    CCTK_REAL const PDdissipationNth3beta1 = PDdissipationNth3(&beta1[index]);
    CCTK_REAL const PDdissipationNth1beta2 = PDdissipationNth1(&beta2[index]);
    CCTK_REAL const PDdissipationNth2beta2 = PDdissipationNth2(&beta2[index]);
    CCTK_REAL const PDdissipationNth3beta2 = PDdissipationNth3(&beta2[index]);
    CCTK_REAL const PDdissipationNth1beta3 = PDdissipationNth1(&beta3[index]);
    CCTK_REAL const PDdissipationNth2beta3 = PDdissipationNth2(&beta3[index]);
    CCTK_REAL const PDdissipationNth3beta3 = PDdissipationNth3(&beta3[index]);
    CCTK_REAL const PDdissipationNth1gt11 = PDdissipationNth1(&gt11[index]);
    CCTK_REAL const PDdissipationNth2gt11 = PDdissipationNth2(&gt11[index]);
    CCTK_REAL const PDdissipationNth3gt11 = PDdissipationNth3(&gt11[index]);
    CCTK_REAL const PDdissipationNth1gt12 = PDdissipationNth1(&gt12[index]);
    CCTK_REAL const PDdissipationNth2gt12 = PDdissipationNth2(&gt12[index]);
    CCTK_REAL const PDdissipationNth3gt12 = PDdissipationNth3(&gt12[index]);
    CCTK_REAL const PDdissipationNth1gt13 = PDdissipationNth1(&gt13[index]);
    CCTK_REAL const PDdissipationNth2gt13 = PDdissipationNth2(&gt13[index]);
    CCTK_REAL const PDdissipationNth3gt13 = PDdissipationNth3(&gt13[index]);
    CCTK_REAL const PDdissipationNth1gt22 = PDdissipationNth1(&gt22[index]);
    CCTK_REAL const PDdissipationNth2gt22 = PDdissipationNth2(&gt22[index]);
    CCTK_REAL const PDdissipationNth3gt22 = PDdissipationNth3(&gt22[index]);
    CCTK_REAL const PDdissipationNth1gt23 = PDdissipationNth1(&gt23[index]);
    CCTK_REAL const PDdissipationNth2gt23 = PDdissipationNth2(&gt23[index]);
    CCTK_REAL const PDdissipationNth3gt23 = PDdissipationNth3(&gt23[index]);
    CCTK_REAL const PDdissipationNth1gt33 = PDdissipationNth1(&gt33[index]);
    CCTK_REAL const PDdissipationNth2gt33 = PDdissipationNth2(&gt33[index]);
    CCTK_REAL const PDdissipationNth3gt33 = PDdissipationNth3(&gt33[index]);
    CCTK_REAL const PDdissipationNth1phi = PDdissipationNth1(&phi[index]);
    CCTK_REAL const PDdissipationNth2phi = PDdissipationNth2(&phi[index]);
    CCTK_REAL const PDdissipationNth3phi = PDdissipationNth3(&phi[index]);
    CCTK_REAL const PDdissipationNth1trK = PDdissipationNth1(&trK[index]);
    CCTK_REAL const PDdissipationNth2trK = PDdissipationNth2(&trK[index]);
    CCTK_REAL const PDdissipationNth3trK = PDdissipationNth3(&trK[index]);
    CCTK_REAL const PDdissipationNth1Xt1 = PDdissipationNth1(&Xt1[index]);
    CCTK_REAL const PDdissipationNth2Xt1 = PDdissipationNth2(&Xt1[index]);
    CCTK_REAL const PDdissipationNth3Xt1 = PDdissipationNth3(&Xt1[index]);
    CCTK_REAL const PDdissipationNth1Xt2 = PDdissipationNth1(&Xt2[index]);
    CCTK_REAL const PDdissipationNth2Xt2 = PDdissipationNth2(&Xt2[index]);
    CCTK_REAL const PDdissipationNth3Xt2 = PDdissipationNth3(&Xt2[index]);
    CCTK_REAL const PDdissipationNth1Xt3 = PDdissipationNth1(&Xt3[index]);
    CCTK_REAL const PDdissipationNth2Xt3 = PDdissipationNth2(&Xt3[index]);
    CCTK_REAL const PDdissipationNth3Xt3 = PDdissipationNth3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL JacPDdissipationNth1A;
    CCTK_REAL JacPDdissipationNth1alpha;
    CCTK_REAL JacPDdissipationNth1At11;
    CCTK_REAL JacPDdissipationNth1At12;
    CCTK_REAL JacPDdissipationNth1At13;
    CCTK_REAL JacPDdissipationNth1At22;
    CCTK_REAL JacPDdissipationNth1At23;
    CCTK_REAL JacPDdissipationNth1At33;
    CCTK_REAL JacPDdissipationNth1B1;
    CCTK_REAL JacPDdissipationNth1B2;
    CCTK_REAL JacPDdissipationNth1B3;
    CCTK_REAL JacPDdissipationNth1beta1;
    CCTK_REAL JacPDdissipationNth1beta2;
    CCTK_REAL JacPDdissipationNth1beta3;
    CCTK_REAL JacPDdissipationNth1gt11;
    CCTK_REAL JacPDdissipationNth1gt12;
    CCTK_REAL JacPDdissipationNth1gt13;
    CCTK_REAL JacPDdissipationNth1gt22;
    CCTK_REAL JacPDdissipationNth1gt23;
    CCTK_REAL JacPDdissipationNth1gt33;
    CCTK_REAL JacPDdissipationNth1phi;
    CCTK_REAL JacPDdissipationNth1trK;
    CCTK_REAL JacPDdissipationNth1Xt1;
    CCTK_REAL JacPDdissipationNth1Xt2;
    CCTK_REAL JacPDdissipationNth1Xt3;
    CCTK_REAL JacPDdissipationNth2A;
    CCTK_REAL JacPDdissipationNth2alpha;
    CCTK_REAL JacPDdissipationNth2At11;
    CCTK_REAL JacPDdissipationNth2At12;
    CCTK_REAL JacPDdissipationNth2At13;
    CCTK_REAL JacPDdissipationNth2At22;
    CCTK_REAL JacPDdissipationNth2At23;
    CCTK_REAL JacPDdissipationNth2At33;
    CCTK_REAL JacPDdissipationNth2B1;
    CCTK_REAL JacPDdissipationNth2B2;
    CCTK_REAL JacPDdissipationNth2B3;
    CCTK_REAL JacPDdissipationNth2beta1;
    CCTK_REAL JacPDdissipationNth2beta2;
    CCTK_REAL JacPDdissipationNth2beta3;
    CCTK_REAL JacPDdissipationNth2gt11;
    CCTK_REAL JacPDdissipationNth2gt12;
    CCTK_REAL JacPDdissipationNth2gt13;
    CCTK_REAL JacPDdissipationNth2gt22;
    CCTK_REAL JacPDdissipationNth2gt23;
    CCTK_REAL JacPDdissipationNth2gt33;
    CCTK_REAL JacPDdissipationNth2phi;
    CCTK_REAL JacPDdissipationNth2trK;
    CCTK_REAL JacPDdissipationNth2Xt1;
    CCTK_REAL JacPDdissipationNth2Xt2;
    CCTK_REAL JacPDdissipationNth2Xt3;
    CCTK_REAL JacPDdissipationNth3A;
    CCTK_REAL JacPDdissipationNth3alpha;
    CCTK_REAL JacPDdissipationNth3At11;
    CCTK_REAL JacPDdissipationNth3At12;
    CCTK_REAL JacPDdissipationNth3At13;
    CCTK_REAL JacPDdissipationNth3At22;
    CCTK_REAL JacPDdissipationNth3At23;
    CCTK_REAL JacPDdissipationNth3At33;
    CCTK_REAL JacPDdissipationNth3B1;
    CCTK_REAL JacPDdissipationNth3B2;
    CCTK_REAL JacPDdissipationNth3B3;
    CCTK_REAL JacPDdissipationNth3beta1;
    CCTK_REAL JacPDdissipationNth3beta2;
    CCTK_REAL JacPDdissipationNth3beta3;
    CCTK_REAL JacPDdissipationNth3gt11;
    CCTK_REAL JacPDdissipationNth3gt12;
    CCTK_REAL JacPDdissipationNth3gt13;
    CCTK_REAL JacPDdissipationNth3gt22;
    CCTK_REAL JacPDdissipationNth3gt23;
    CCTK_REAL JacPDdissipationNth3gt33;
    CCTK_REAL JacPDdissipationNth3phi;
    CCTK_REAL JacPDdissipationNth3trK;
    CCTK_REAL JacPDdissipationNth3Xt1;
    CCTK_REAL JacPDdissipationNth3Xt2;
    CCTK_REAL JacPDdissipationNth3Xt3;
    
    if (use_jacobian)
    {
      JacPDdissipationNth1A = J11L*PDdissipationNth1A + 
        J21L*PDdissipationNth2A + J31L*PDdissipationNth3A;
      
      JacPDdissipationNth1alpha = J11L*PDdissipationNth1alpha + 
        J21L*PDdissipationNth2alpha + J31L*PDdissipationNth3alpha;
      
      JacPDdissipationNth1At11 = J11L*PDdissipationNth1At11 + 
        J21L*PDdissipationNth2At11 + J31L*PDdissipationNth3At11;
      
      JacPDdissipationNth1At12 = J11L*PDdissipationNth1At12 + 
        J21L*PDdissipationNth2At12 + J31L*PDdissipationNth3At12;
      
      JacPDdissipationNth1At13 = J11L*PDdissipationNth1At13 + 
        J21L*PDdissipationNth2At13 + J31L*PDdissipationNth3At13;
      
      JacPDdissipationNth1At22 = J11L*PDdissipationNth1At22 + 
        J21L*PDdissipationNth2At22 + J31L*PDdissipationNth3At22;
      
      JacPDdissipationNth1At23 = J11L*PDdissipationNth1At23 + 
        J21L*PDdissipationNth2At23 + J31L*PDdissipationNth3At23;
      
      JacPDdissipationNth1At33 = J11L*PDdissipationNth1At33 + 
        J21L*PDdissipationNth2At33 + J31L*PDdissipationNth3At33;
      
      JacPDdissipationNth1B1 = J11L*PDdissipationNth1B1 + 
        J21L*PDdissipationNth2B1 + J31L*PDdissipationNth3B1;
      
      JacPDdissipationNth1B2 = J11L*PDdissipationNth1B2 + 
        J21L*PDdissipationNth2B2 + J31L*PDdissipationNth3B2;
      
      JacPDdissipationNth1B3 = J11L*PDdissipationNth1B3 + 
        J21L*PDdissipationNth2B3 + J31L*PDdissipationNth3B3;
      
      JacPDdissipationNth1beta1 = J11L*PDdissipationNth1beta1 + 
        J21L*PDdissipationNth2beta1 + J31L*PDdissipationNth3beta1;
      
      JacPDdissipationNth1beta2 = J11L*PDdissipationNth1beta2 + 
        J21L*PDdissipationNth2beta2 + J31L*PDdissipationNth3beta2;
      
      JacPDdissipationNth1beta3 = J11L*PDdissipationNth1beta3 + 
        J21L*PDdissipationNth2beta3 + J31L*PDdissipationNth3beta3;
      
      JacPDdissipationNth1gt11 = J11L*PDdissipationNth1gt11 + 
        J21L*PDdissipationNth2gt11 + J31L*PDdissipationNth3gt11;
      
      JacPDdissipationNth1gt12 = J11L*PDdissipationNth1gt12 + 
        J21L*PDdissipationNth2gt12 + J31L*PDdissipationNth3gt12;
      
      JacPDdissipationNth1gt13 = J11L*PDdissipationNth1gt13 + 
        J21L*PDdissipationNth2gt13 + J31L*PDdissipationNth3gt13;
      
      JacPDdissipationNth1gt22 = J11L*PDdissipationNth1gt22 + 
        J21L*PDdissipationNth2gt22 + J31L*PDdissipationNth3gt22;
      
      JacPDdissipationNth1gt23 = J11L*PDdissipationNth1gt23 + 
        J21L*PDdissipationNth2gt23 + J31L*PDdissipationNth3gt23;
      
      JacPDdissipationNth1gt33 = J11L*PDdissipationNth1gt33 + 
        J21L*PDdissipationNth2gt33 + J31L*PDdissipationNth3gt33;
      
      JacPDdissipationNth1phi = J11L*PDdissipationNth1phi + 
        J21L*PDdissipationNth2phi + J31L*PDdissipationNth3phi;
      
      JacPDdissipationNth1trK = J11L*PDdissipationNth1trK + 
        J21L*PDdissipationNth2trK + J31L*PDdissipationNth3trK;
      
      JacPDdissipationNth1Xt1 = J11L*PDdissipationNth1Xt1 + 
        J21L*PDdissipationNth2Xt1 + J31L*PDdissipationNth3Xt1;
      
      JacPDdissipationNth1Xt2 = J11L*PDdissipationNth1Xt2 + 
        J21L*PDdissipationNth2Xt2 + J31L*PDdissipationNth3Xt2;
      
      JacPDdissipationNth1Xt3 = J11L*PDdissipationNth1Xt3 + 
        J21L*PDdissipationNth2Xt3 + J31L*PDdissipationNth3Xt3;
      
      JacPDdissipationNth2A = J12L*PDdissipationNth1A + 
        J22L*PDdissipationNth2A + J32L*PDdissipationNth3A;
      
      JacPDdissipationNth2alpha = J12L*PDdissipationNth1alpha + 
        J22L*PDdissipationNth2alpha + J32L*PDdissipationNth3alpha;
      
      JacPDdissipationNth2At11 = J12L*PDdissipationNth1At11 + 
        J22L*PDdissipationNth2At11 + J32L*PDdissipationNth3At11;
      
      JacPDdissipationNth2At12 = J12L*PDdissipationNth1At12 + 
        J22L*PDdissipationNth2At12 + J32L*PDdissipationNth3At12;
      
      JacPDdissipationNth2At13 = J12L*PDdissipationNth1At13 + 
        J22L*PDdissipationNth2At13 + J32L*PDdissipationNth3At13;
      
      JacPDdissipationNth2At22 = J12L*PDdissipationNth1At22 + 
        J22L*PDdissipationNth2At22 + J32L*PDdissipationNth3At22;
      
      JacPDdissipationNth2At23 = J12L*PDdissipationNth1At23 + 
        J22L*PDdissipationNth2At23 + J32L*PDdissipationNth3At23;
      
      JacPDdissipationNth2At33 = J12L*PDdissipationNth1At33 + 
        J22L*PDdissipationNth2At33 + J32L*PDdissipationNth3At33;
      
      JacPDdissipationNth2B1 = J12L*PDdissipationNth1B1 + 
        J22L*PDdissipationNth2B1 + J32L*PDdissipationNth3B1;
      
      JacPDdissipationNth2B2 = J12L*PDdissipationNth1B2 + 
        J22L*PDdissipationNth2B2 + J32L*PDdissipationNth3B2;
      
      JacPDdissipationNth2B3 = J12L*PDdissipationNth1B3 + 
        J22L*PDdissipationNth2B3 + J32L*PDdissipationNth3B3;
      
      JacPDdissipationNth2beta1 = J12L*PDdissipationNth1beta1 + 
        J22L*PDdissipationNth2beta1 + J32L*PDdissipationNth3beta1;
      
      JacPDdissipationNth2beta2 = J12L*PDdissipationNth1beta2 + 
        J22L*PDdissipationNth2beta2 + J32L*PDdissipationNth3beta2;
      
      JacPDdissipationNth2beta3 = J12L*PDdissipationNth1beta3 + 
        J22L*PDdissipationNth2beta3 + J32L*PDdissipationNth3beta3;
      
      JacPDdissipationNth2gt11 = J12L*PDdissipationNth1gt11 + 
        J22L*PDdissipationNth2gt11 + J32L*PDdissipationNth3gt11;
      
      JacPDdissipationNth2gt12 = J12L*PDdissipationNth1gt12 + 
        J22L*PDdissipationNth2gt12 + J32L*PDdissipationNth3gt12;
      
      JacPDdissipationNth2gt13 = J12L*PDdissipationNth1gt13 + 
        J22L*PDdissipationNth2gt13 + J32L*PDdissipationNth3gt13;
      
      JacPDdissipationNth2gt22 = J12L*PDdissipationNth1gt22 + 
        J22L*PDdissipationNth2gt22 + J32L*PDdissipationNth3gt22;
      
      JacPDdissipationNth2gt23 = J12L*PDdissipationNth1gt23 + 
        J22L*PDdissipationNth2gt23 + J32L*PDdissipationNth3gt23;
      
      JacPDdissipationNth2gt33 = J12L*PDdissipationNth1gt33 + 
        J22L*PDdissipationNth2gt33 + J32L*PDdissipationNth3gt33;
      
      JacPDdissipationNth2phi = J12L*PDdissipationNth1phi + 
        J22L*PDdissipationNth2phi + J32L*PDdissipationNth3phi;
      
      JacPDdissipationNth2trK = J12L*PDdissipationNth1trK + 
        J22L*PDdissipationNth2trK + J32L*PDdissipationNth3trK;
      
      JacPDdissipationNth2Xt1 = J12L*PDdissipationNth1Xt1 + 
        J22L*PDdissipationNth2Xt1 + J32L*PDdissipationNth3Xt1;
      
      JacPDdissipationNth2Xt2 = J12L*PDdissipationNth1Xt2 + 
        J22L*PDdissipationNth2Xt2 + J32L*PDdissipationNth3Xt2;
      
      JacPDdissipationNth2Xt3 = J12L*PDdissipationNth1Xt3 + 
        J22L*PDdissipationNth2Xt3 + J32L*PDdissipationNth3Xt3;
      
      JacPDdissipationNth3A = J13L*PDdissipationNth1A + 
        J23L*PDdissipationNth2A + J33L*PDdissipationNth3A;
      
      JacPDdissipationNth3alpha = J13L*PDdissipationNth1alpha + 
        J23L*PDdissipationNth2alpha + J33L*PDdissipationNth3alpha;
      
      JacPDdissipationNth3At11 = J13L*PDdissipationNth1At11 + 
        J23L*PDdissipationNth2At11 + J33L*PDdissipationNth3At11;
      
      JacPDdissipationNth3At12 = J13L*PDdissipationNth1At12 + 
        J23L*PDdissipationNth2At12 + J33L*PDdissipationNth3At12;
      
      JacPDdissipationNth3At13 = J13L*PDdissipationNth1At13 + 
        J23L*PDdissipationNth2At13 + J33L*PDdissipationNth3At13;
      
      JacPDdissipationNth3At22 = J13L*PDdissipationNth1At22 + 
        J23L*PDdissipationNth2At22 + J33L*PDdissipationNth3At22;
      
      JacPDdissipationNth3At23 = J13L*PDdissipationNth1At23 + 
        J23L*PDdissipationNth2At23 + J33L*PDdissipationNth3At23;
      
      JacPDdissipationNth3At33 = J13L*PDdissipationNth1At33 + 
        J23L*PDdissipationNth2At33 + J33L*PDdissipationNth3At33;
      
      JacPDdissipationNth3B1 = J13L*PDdissipationNth1B1 + 
        J23L*PDdissipationNth2B1 + J33L*PDdissipationNth3B1;
      
      JacPDdissipationNth3B2 = J13L*PDdissipationNth1B2 + 
        J23L*PDdissipationNth2B2 + J33L*PDdissipationNth3B2;
      
      JacPDdissipationNth3B3 = J13L*PDdissipationNth1B3 + 
        J23L*PDdissipationNth2B3 + J33L*PDdissipationNth3B3;
      
      JacPDdissipationNth3beta1 = J13L*PDdissipationNth1beta1 + 
        J23L*PDdissipationNth2beta1 + J33L*PDdissipationNth3beta1;
      
      JacPDdissipationNth3beta2 = J13L*PDdissipationNth1beta2 + 
        J23L*PDdissipationNth2beta2 + J33L*PDdissipationNth3beta2;
      
      JacPDdissipationNth3beta3 = J13L*PDdissipationNth1beta3 + 
        J23L*PDdissipationNth2beta3 + J33L*PDdissipationNth3beta3;
      
      JacPDdissipationNth3gt11 = J13L*PDdissipationNth1gt11 + 
        J23L*PDdissipationNth2gt11 + J33L*PDdissipationNth3gt11;
      
      JacPDdissipationNth3gt12 = J13L*PDdissipationNth1gt12 + 
        J23L*PDdissipationNth2gt12 + J33L*PDdissipationNth3gt12;
      
      JacPDdissipationNth3gt13 = J13L*PDdissipationNth1gt13 + 
        J23L*PDdissipationNth2gt13 + J33L*PDdissipationNth3gt13;
      
      JacPDdissipationNth3gt22 = J13L*PDdissipationNth1gt22 + 
        J23L*PDdissipationNth2gt22 + J33L*PDdissipationNth3gt22;
      
      JacPDdissipationNth3gt23 = J13L*PDdissipationNth1gt23 + 
        J23L*PDdissipationNth2gt23 + J33L*PDdissipationNth3gt23;
      
      JacPDdissipationNth3gt33 = J13L*PDdissipationNth1gt33 + 
        J23L*PDdissipationNth2gt33 + J33L*PDdissipationNth3gt33;
      
      JacPDdissipationNth3phi = J13L*PDdissipationNth1phi + 
        J23L*PDdissipationNth2phi + J33L*PDdissipationNth3phi;
      
      JacPDdissipationNth3trK = J13L*PDdissipationNth1trK + 
        J23L*PDdissipationNth2trK + J33L*PDdissipationNth3trK;
      
      JacPDdissipationNth3Xt1 = J13L*PDdissipationNth1Xt1 + 
        J23L*PDdissipationNth2Xt1 + J33L*PDdissipationNth3Xt1;
      
      JacPDdissipationNth3Xt2 = J13L*PDdissipationNth1Xt2 + 
        J23L*PDdissipationNth2Xt2 + J33L*PDdissipationNth3Xt2;
      
      JacPDdissipationNth3Xt3 = J13L*PDdissipationNth1Xt3 + 
        J23L*PDdissipationNth2Xt3 + J33L*PDdissipationNth3Xt3;
    }
    else
    {
      JacPDdissipationNth1A = PDdissipationNth1A;
      
      JacPDdissipationNth1alpha = PDdissipationNth1alpha;
      
      JacPDdissipationNth1At11 = PDdissipationNth1At11;
      
      JacPDdissipationNth1At12 = PDdissipationNth1At12;
      
      JacPDdissipationNth1At13 = PDdissipationNth1At13;
      
      JacPDdissipationNth1At22 = PDdissipationNth1At22;
      
      JacPDdissipationNth1At23 = PDdissipationNth1At23;
      
      JacPDdissipationNth1At33 = PDdissipationNth1At33;
      
      JacPDdissipationNth1B1 = PDdissipationNth1B1;
      
      JacPDdissipationNth1B2 = PDdissipationNth1B2;
      
      JacPDdissipationNth1B3 = PDdissipationNth1B3;
      
      JacPDdissipationNth1beta1 = PDdissipationNth1beta1;
      
      JacPDdissipationNth1beta2 = PDdissipationNth1beta2;
      
      JacPDdissipationNth1beta3 = PDdissipationNth1beta3;
      
      JacPDdissipationNth1gt11 = PDdissipationNth1gt11;
      
      JacPDdissipationNth1gt12 = PDdissipationNth1gt12;
      
      JacPDdissipationNth1gt13 = PDdissipationNth1gt13;
      
      JacPDdissipationNth1gt22 = PDdissipationNth1gt22;
      
      JacPDdissipationNth1gt23 = PDdissipationNth1gt23;
      
      JacPDdissipationNth1gt33 = PDdissipationNth1gt33;
      
      JacPDdissipationNth1phi = PDdissipationNth1phi;
      
      JacPDdissipationNth1trK = PDdissipationNth1trK;
      
      JacPDdissipationNth1Xt1 = PDdissipationNth1Xt1;
      
      JacPDdissipationNth1Xt2 = PDdissipationNth1Xt2;
      
      JacPDdissipationNth1Xt3 = PDdissipationNth1Xt3;
      
      JacPDdissipationNth2A = PDdissipationNth2A;
      
      JacPDdissipationNth2alpha = PDdissipationNth2alpha;
      
      JacPDdissipationNth2At11 = PDdissipationNth2At11;
      
      JacPDdissipationNth2At12 = PDdissipationNth2At12;
      
      JacPDdissipationNth2At13 = PDdissipationNth2At13;
      
      JacPDdissipationNth2At22 = PDdissipationNth2At22;
      
      JacPDdissipationNth2At23 = PDdissipationNth2At23;
      
      JacPDdissipationNth2At33 = PDdissipationNth2At33;
      
      JacPDdissipationNth2B1 = PDdissipationNth2B1;
      
      JacPDdissipationNth2B2 = PDdissipationNth2B2;
      
      JacPDdissipationNth2B3 = PDdissipationNth2B3;
      
      JacPDdissipationNth2beta1 = PDdissipationNth2beta1;
      
      JacPDdissipationNth2beta2 = PDdissipationNth2beta2;
      
      JacPDdissipationNth2beta3 = PDdissipationNth2beta3;
      
      JacPDdissipationNth2gt11 = PDdissipationNth2gt11;
      
      JacPDdissipationNth2gt12 = PDdissipationNth2gt12;
      
      JacPDdissipationNth2gt13 = PDdissipationNth2gt13;
      
      JacPDdissipationNth2gt22 = PDdissipationNth2gt22;
      
      JacPDdissipationNth2gt23 = PDdissipationNth2gt23;
      
      JacPDdissipationNth2gt33 = PDdissipationNth2gt33;
      
      JacPDdissipationNth2phi = PDdissipationNth2phi;
      
      JacPDdissipationNth2trK = PDdissipationNth2trK;
      
      JacPDdissipationNth2Xt1 = PDdissipationNth2Xt1;
      
      JacPDdissipationNth2Xt2 = PDdissipationNth2Xt2;
      
      JacPDdissipationNth2Xt3 = PDdissipationNth2Xt3;
      
      JacPDdissipationNth3A = PDdissipationNth3A;
      
      JacPDdissipationNth3alpha = PDdissipationNth3alpha;
      
      JacPDdissipationNth3At11 = PDdissipationNth3At11;
      
      JacPDdissipationNth3At12 = PDdissipationNth3At12;
      
      JacPDdissipationNth3At13 = PDdissipationNth3At13;
      
      JacPDdissipationNth3At22 = PDdissipationNth3At22;
      
      JacPDdissipationNth3At23 = PDdissipationNth3At23;
      
      JacPDdissipationNth3At33 = PDdissipationNth3At33;
      
      JacPDdissipationNth3B1 = PDdissipationNth3B1;
      
      JacPDdissipationNth3B2 = PDdissipationNth3B2;
      
      JacPDdissipationNth3B3 = PDdissipationNth3B3;
      
      JacPDdissipationNth3beta1 = PDdissipationNth3beta1;
      
      JacPDdissipationNth3beta2 = PDdissipationNth3beta2;
      
      JacPDdissipationNth3beta3 = PDdissipationNth3beta3;
      
      JacPDdissipationNth3gt11 = PDdissipationNth3gt11;
      
      JacPDdissipationNth3gt12 = PDdissipationNth3gt12;
      
      JacPDdissipationNth3gt13 = PDdissipationNth3gt13;
      
      JacPDdissipationNth3gt22 = PDdissipationNth3gt22;
      
      JacPDdissipationNth3gt23 = PDdissipationNth3gt23;
      
      JacPDdissipationNth3gt33 = PDdissipationNth3gt33;
      
      JacPDdissipationNth3phi = PDdissipationNth3phi;
      
      JacPDdissipationNth3trK = PDdissipationNth3trK;
      
      JacPDdissipationNth3Xt1 = PDdissipationNth3Xt1;
      
      JacPDdissipationNth3Xt2 = PDdissipationNth3Xt2;
      
      JacPDdissipationNth3Xt3 = PDdissipationNth3Xt3;
    }
    
    CCTK_REAL epsdiss1 = ToReal(EpsDiss);
    
    CCTK_REAL epsdiss2 = ToReal(EpsDiss);
    
    CCTK_REAL epsdiss3 = ToReal(EpsDiss);
    
    phirhsL = epsdiss1*JacPDdissipationNth1phi + 
      epsdiss2*JacPDdissipationNth2phi + epsdiss3*JacPDdissipationNth3phi + 
      phirhsL;
    
    gt11rhsL = gt11rhsL + epsdiss1*JacPDdissipationNth1gt11 + 
      epsdiss2*JacPDdissipationNth2gt11 + epsdiss3*JacPDdissipationNth3gt11;
    
    gt12rhsL = gt12rhsL + epsdiss1*JacPDdissipationNth1gt12 + 
      epsdiss2*JacPDdissipationNth2gt12 + epsdiss3*JacPDdissipationNth3gt12;
    
    gt13rhsL = gt13rhsL + epsdiss1*JacPDdissipationNth1gt13 + 
      epsdiss2*JacPDdissipationNth2gt13 + epsdiss3*JacPDdissipationNth3gt13;
    
    gt22rhsL = gt22rhsL + epsdiss1*JacPDdissipationNth1gt22 + 
      epsdiss2*JacPDdissipationNth2gt22 + epsdiss3*JacPDdissipationNth3gt22;
    
    gt23rhsL = gt23rhsL + epsdiss1*JacPDdissipationNth1gt23 + 
      epsdiss2*JacPDdissipationNth2gt23 + epsdiss3*JacPDdissipationNth3gt23;
    
    gt33rhsL = gt33rhsL + epsdiss1*JacPDdissipationNth1gt33 + 
      epsdiss2*JacPDdissipationNth2gt33 + epsdiss3*JacPDdissipationNth3gt33;
    
    Xt1rhsL = epsdiss1*JacPDdissipationNth1Xt1 + 
      epsdiss2*JacPDdissipationNth2Xt1 + epsdiss3*JacPDdissipationNth3Xt1 + 
      Xt1rhsL;
    
    Xt2rhsL = epsdiss1*JacPDdissipationNth1Xt2 + 
      epsdiss2*JacPDdissipationNth2Xt2 + epsdiss3*JacPDdissipationNth3Xt2 + 
      Xt2rhsL;
    
    Xt3rhsL = epsdiss1*JacPDdissipationNth1Xt3 + 
      epsdiss2*JacPDdissipationNth2Xt3 + epsdiss3*JacPDdissipationNth3Xt3 + 
      Xt3rhsL;
    
    trKrhsL = epsdiss1*JacPDdissipationNth1trK + 
      epsdiss2*JacPDdissipationNth2trK + epsdiss3*JacPDdissipationNth3trK + 
      trKrhsL;
    
    At11rhsL = At11rhsL + epsdiss1*JacPDdissipationNth1At11 + 
      epsdiss2*JacPDdissipationNth2At11 + epsdiss3*JacPDdissipationNth3At11;
    
    At12rhsL = At12rhsL + epsdiss1*JacPDdissipationNth1At12 + 
      epsdiss2*JacPDdissipationNth2At12 + epsdiss3*JacPDdissipationNth3At12;
    
    At13rhsL = At13rhsL + epsdiss1*JacPDdissipationNth1At13 + 
      epsdiss2*JacPDdissipationNth2At13 + epsdiss3*JacPDdissipationNth3At13;
    
    At22rhsL = At22rhsL + epsdiss1*JacPDdissipationNth1At22 + 
      epsdiss2*JacPDdissipationNth2At22 + epsdiss3*JacPDdissipationNth3At22;
    
    At23rhsL = At23rhsL + epsdiss1*JacPDdissipationNth1At23 + 
      epsdiss2*JacPDdissipationNth2At23 + epsdiss3*JacPDdissipationNth3At23;
    
    At33rhsL = At33rhsL + epsdiss1*JacPDdissipationNth1At33 + 
      epsdiss2*JacPDdissipationNth2At33 + epsdiss3*JacPDdissipationNth3At33;
    
    alpharhsL = alpharhsL + epsdiss1*JacPDdissipationNth1alpha + 
      epsdiss2*JacPDdissipationNth2alpha + 
      epsdiss3*JacPDdissipationNth3alpha;
    
    ArhsL = ArhsL + epsdiss1*JacPDdissipationNth1A + 
      epsdiss2*JacPDdissipationNth2A + epsdiss3*JacPDdissipationNth3A;
    
    beta1rhsL = beta1rhsL + epsdiss1*JacPDdissipationNth1beta1 + 
      epsdiss2*JacPDdissipationNth2beta1 + 
      epsdiss3*JacPDdissipationNth3beta1;
    
    beta2rhsL = beta2rhsL + epsdiss1*JacPDdissipationNth1beta2 + 
      epsdiss2*JacPDdissipationNth2beta2 + 
      epsdiss3*JacPDdissipationNth3beta2;
    
    beta3rhsL = beta3rhsL + epsdiss1*JacPDdissipationNth1beta3 + 
      epsdiss2*JacPDdissipationNth2beta3 + 
      epsdiss3*JacPDdissipationNth3beta3;
    
    B1rhsL = B1rhsL + epsdiss1*JacPDdissipationNth1B1 + 
      epsdiss2*JacPDdissipationNth2B1 + epsdiss3*JacPDdissipationNth3B1;
    
    B2rhsL = B2rhsL + epsdiss1*JacPDdissipationNth1B2 + 
      epsdiss2*JacPDdissipationNth2B2 + epsdiss3*JacPDdissipationNth3B2;
    
    B3rhsL = B3rhsL + epsdiss1*JacPDdissipationNth1B3 + 
      epsdiss2*JacPDdissipationNth2B3 + epsdiss3*JacPDdissipationNth3B3;
    
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
  LC_ENDLOOP3 (ML_BSSN_MP_O8_Dissipation);
}

extern "C" void ML_BSSN_MP_O8_Dissipation(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_O8_Dissipation_Body);
}
