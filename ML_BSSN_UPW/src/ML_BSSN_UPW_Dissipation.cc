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
#include "vectors.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define QAD(x) (SQR(SQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))

extern "C" void ML_BSSN_UPW_Dissipation_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_curvrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_trace_curvrhs.");
  return;
}

static void ML_BSSN_UPW_Dissipation_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  ptrdiff_t const di = 1;
  ptrdiff_t const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const cdi = sizeof(CCTK_REAL) * di;
  ptrdiff_t const cdj = sizeof(CCTK_REAL) * dj;
  ptrdiff_t const cdk = sizeof(CCTK_REAL) * dk;
  CCTK_REAL_VEC const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL_VEC const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL_VEC const dz = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL_VEC const dt = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL_VEC const t = ToReal(cctk_time);
  CCTK_REAL_VEC const dxi = INV(dx);
  CCTK_REAL_VEC const dyi = INV(dy);
  CCTK_REAL_VEC const dzi = INV(dz);
  CCTK_REAL_VEC const khalf = ToReal(0.5);
  CCTK_REAL_VEC const kthird = ToReal(1.0/3.0);
  CCTK_REAL_VEC const ktwothird = ToReal(2.0/3.0);
  CCTK_REAL_VEC const kfourthird = ToReal(4.0/3.0);
  CCTK_REAL_VEC const keightthird = ToReal(8.0/3.0);
  CCTK_REAL_VEC const hdxi = kmul(ToReal(0.5), dxi);
  CCTK_REAL_VEC const hdyi = kmul(ToReal(0.5), dyi);
  CCTK_REAL_VEC const hdzi = kmul(ToReal(0.5), dzi);
  
  /* Initialize predefined quantities */
  CCTK_REAL_VEC const p1o1024dx = kmul(INV(dx),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1024dy = kmul(INV(dy),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1024dz = kmul(INV(dz),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o120dx = kmul(INV(dx),ToReal(0.00833333333333333333333333333333));
  CCTK_REAL_VEC const p1o120dy = kmul(INV(dy),ToReal(0.00833333333333333333333333333333));
  CCTK_REAL_VEC const p1o120dz = kmul(INV(dz),ToReal(0.00833333333333333333333333333333));
  CCTK_REAL_VEC const p1o12dx = kmul(INV(dx),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o12dy = kmul(INV(dy),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o12dz = kmul(INV(dz),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o144dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o144dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o144dydz = kmul(INV(dy),kmul(INV(dz),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o1680dx = kmul(INV(dx),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o1680dy = kmul(INV(dy),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o1680dz = kmul(INV(dz),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o16dx = kmul(INV(dx),ToReal(0.0625));
  CCTK_REAL_VEC const p1o16dy = kmul(INV(dy),ToReal(0.0625));
  CCTK_REAL_VEC const p1o16dz = kmul(INV(dz),ToReal(0.0625));
  CCTK_REAL_VEC const p1o180dx2 = kmul(INV(SQR(dx)),ToReal(0.00555555555555555555555555555556));
  CCTK_REAL_VEC const p1o180dy2 = kmul(INV(SQR(dy)),ToReal(0.00555555555555555555555555555556));
  CCTK_REAL_VEC const p1o180dz2 = kmul(INV(SQR(dz)),ToReal(0.00555555555555555555555555555556));
  CCTK_REAL_VEC const p1o24dx = kmul(INV(dx),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o24dy = kmul(INV(dy),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o24dz = kmul(INV(dz),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o256dx = kmul(INV(dx),ToReal(0.00390625));
  CCTK_REAL_VEC const p1o256dy = kmul(INV(dy),ToReal(0.00390625));
  CCTK_REAL_VEC const p1o256dz = kmul(INV(dz),ToReal(0.00390625));
  CCTK_REAL_VEC const p1o2dx = kmul(INV(dx),ToReal(0.5));
  CCTK_REAL_VEC const p1o2dy = kmul(INV(dy),ToReal(0.5));
  CCTK_REAL_VEC const p1o2dz = kmul(INV(dz),ToReal(0.5));
  CCTK_REAL_VEC const p1o3600dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(0.000277777777777777777777777777778)));
  CCTK_REAL_VEC const p1o3600dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(0.000277777777777777777777777777778)));
  CCTK_REAL_VEC const p1o3600dydz = kmul(INV(dy),kmul(INV(dz),ToReal(0.000277777777777777777777777777778)));
  CCTK_REAL_VEC const p1o4dx = kmul(INV(dx),ToReal(0.25));
  CCTK_REAL_VEC const p1o4dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(0.25)));
  CCTK_REAL_VEC const p1o4dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(0.25)));
  CCTK_REAL_VEC const p1o4dy = kmul(INV(dy),ToReal(0.25));
  CCTK_REAL_VEC const p1o4dydz = kmul(INV(dy),kmul(INV(dz),ToReal(0.25)));
  CCTK_REAL_VEC const p1o4dz = kmul(INV(dz),ToReal(0.25));
  CCTK_REAL_VEC const p1o5040dx2 = kmul(INV(SQR(dx)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dy2 = kmul(INV(SQR(dy)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dz2 = kmul(INV(SQR(dz)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o560dx = kmul(INV(dx),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o560dy = kmul(INV(dy),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o560dz = kmul(INV(dz),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o60dx = kmul(INV(dx),ToReal(0.0166666666666666666666666666667));
  CCTK_REAL_VEC const p1o60dy = kmul(INV(dy),ToReal(0.0166666666666666666666666666667));
  CCTK_REAL_VEC const p1o60dz = kmul(INV(dz),ToReal(0.0166666666666666666666666666667));
  CCTK_REAL_VEC const p1o64dx = kmul(INV(dx),ToReal(0.015625));
  CCTK_REAL_VEC const p1o64dy = kmul(INV(dy),ToReal(0.015625));
  CCTK_REAL_VEC const p1o64dz = kmul(INV(dz),ToReal(0.015625));
  CCTK_REAL_VEC const p1o705600dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o705600dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o705600dydz = kmul(INV(dy),kmul(INV(dz),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o840dx = kmul(INV(dx),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dy = kmul(INV(dy),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dz = kmul(INV(dz),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1odx = INV(dx);
  CCTK_REAL_VEC const p1odx2 = INV(SQR(dx));
  CCTK_REAL_VEC const p1ody = INV(dy);
  CCTK_REAL_VEC const p1ody2 = INV(SQR(dy));
  CCTK_REAL_VEC const p1odz = INV(dz);
  CCTK_REAL_VEC const p1odz2 = INV(SQR(dz));
  CCTK_REAL_VEC const pm1o120dx = kmul(INV(dx),ToReal(-0.00833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o120dy = kmul(INV(dy),ToReal(-0.00833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o120dz = kmul(INV(dz),ToReal(-0.00833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dx2 = kmul(INV(SQR(dx)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dy2 = kmul(INV(SQR(dy)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dz2 = kmul(INV(SQR(dz)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o2dx = kmul(INV(dx),ToReal(-0.5));
  CCTK_REAL_VEC const pm1o2dy = kmul(INV(dy),ToReal(-0.5));
  CCTK_REAL_VEC const pm1o2dz = kmul(INV(dz),ToReal(-0.5));
  CCTK_REAL_VEC const pm1o4dx = kmul(INV(dx),ToReal(-0.25));
  CCTK_REAL_VEC const pm1o4dy = kmul(INV(dy),ToReal(-0.25));
  CCTK_REAL_VEC const pm1o4dz = kmul(INV(dz),ToReal(-0.25));
  CCTK_REAL_VEC const pm1o60dx = kmul(INV(dx),ToReal(-0.0166666666666666666666666666667));
  CCTK_REAL_VEC const pm1o60dy = kmul(INV(dy),ToReal(-0.0166666666666666666666666666667));
  CCTK_REAL_VEC const pm1o60dz = kmul(INV(dz),ToReal(-0.0166666666666666666666666666667));
  
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
  LC_LOOP3VEC (ML_BSSN_UPW_Dissipation,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC AL = vec_load(A[index]);
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC alpharhsL = vec_load(alpharhs[index]);
    CCTK_REAL_VEC ArhsL = vec_load(Arhs[index]);
    CCTK_REAL_VEC At11L = vec_load(At11[index]);
    CCTK_REAL_VEC At11rhsL = vec_load(At11rhs[index]);
    CCTK_REAL_VEC At12L = vec_load(At12[index]);
    CCTK_REAL_VEC At12rhsL = vec_load(At12rhs[index]);
    CCTK_REAL_VEC At13L = vec_load(At13[index]);
    CCTK_REAL_VEC At13rhsL = vec_load(At13rhs[index]);
    CCTK_REAL_VEC At22L = vec_load(At22[index]);
    CCTK_REAL_VEC At22rhsL = vec_load(At22rhs[index]);
    CCTK_REAL_VEC At23L = vec_load(At23[index]);
    CCTK_REAL_VEC At23rhsL = vec_load(At23rhs[index]);
    CCTK_REAL_VEC At33L = vec_load(At33[index]);
    CCTK_REAL_VEC At33rhsL = vec_load(At33rhs[index]);
    CCTK_REAL_VEC B1L = vec_load(B1[index]);
    CCTK_REAL_VEC B1rhsL = vec_load(B1rhs[index]);
    CCTK_REAL_VEC B2L = vec_load(B2[index]);
    CCTK_REAL_VEC B2rhsL = vec_load(B2rhs[index]);
    CCTK_REAL_VEC B3L = vec_load(B3[index]);
    CCTK_REAL_VEC B3rhsL = vec_load(B3rhs[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta1rhsL = vec_load(beta1rhs[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta2rhsL = vec_load(beta2rhs[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC beta3rhsL = vec_load(beta3rhs[index]);
    CCTK_REAL_VEC gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC gt11rhsL = vec_load(gt11rhs[index]);
    CCTK_REAL_VEC gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC gt12rhsL = vec_load(gt12rhs[index]);
    CCTK_REAL_VEC gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC gt13rhsL = vec_load(gt13rhs[index]);
    CCTK_REAL_VEC gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC gt22rhsL = vec_load(gt22rhs[index]);
    CCTK_REAL_VEC gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC gt23rhsL = vec_load(gt23rhs[index]);
    CCTK_REAL_VEC gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC gt33rhsL = vec_load(gt33rhs[index]);
    CCTK_REAL_VEC phiL = vec_load(phi[index]);
    CCTK_REAL_VEC phirhsL = vec_load(phirhs[index]);
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    CCTK_REAL_VEC trKrhsL = vec_load(trKrhs[index]);
    CCTK_REAL_VEC Xt1L = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt1rhsL = vec_load(Xt1rhs[index]);
    CCTK_REAL_VEC Xt2L = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt2rhsL = vec_load(Xt2rhs[index]);
    CCTK_REAL_VEC Xt3L = vec_load(Xt3[index]);
    CCTK_REAL_VEC Xt3rhsL = vec_load(Xt3rhs[index]);
    
    
    CCTK_REAL_VEC J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L;
    
    if (use_jacobian)
    {
      J11L = vec_load(J11[index]);
      J12L = vec_load(J12[index]);
      J13L = vec_load(J13[index]);
      J21L = vec_load(J21[index]);
      J22L = vec_load(J22[index]);
      J23L = vec_load(J23[index]);
      J31L = vec_load(J31[index]);
      J32L = vec_load(J32[index]);
      J33L = vec_load(J33[index]);
    }
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC PDdissipationNth1A;
    CCTK_REAL_VEC PDdissipationNth2A;
    CCTK_REAL_VEC PDdissipationNth3A;
    CCTK_REAL_VEC PDdissipationNth1alpha;
    CCTK_REAL_VEC PDdissipationNth2alpha;
    CCTK_REAL_VEC PDdissipationNth3alpha;
    CCTK_REAL_VEC PDdissipationNth1At11;
    CCTK_REAL_VEC PDdissipationNth2At11;
    CCTK_REAL_VEC PDdissipationNth3At11;
    CCTK_REAL_VEC PDdissipationNth1At12;
    CCTK_REAL_VEC PDdissipationNth2At12;
    CCTK_REAL_VEC PDdissipationNth3At12;
    CCTK_REAL_VEC PDdissipationNth1At13;
    CCTK_REAL_VEC PDdissipationNth2At13;
    CCTK_REAL_VEC PDdissipationNth3At13;
    CCTK_REAL_VEC PDdissipationNth1At22;
    CCTK_REAL_VEC PDdissipationNth2At22;
    CCTK_REAL_VEC PDdissipationNth3At22;
    CCTK_REAL_VEC PDdissipationNth1At23;
    CCTK_REAL_VEC PDdissipationNth2At23;
    CCTK_REAL_VEC PDdissipationNth3At23;
    CCTK_REAL_VEC PDdissipationNth1At33;
    CCTK_REAL_VEC PDdissipationNth2At33;
    CCTK_REAL_VEC PDdissipationNth3At33;
    CCTK_REAL_VEC PDdissipationNth1B1;
    CCTK_REAL_VEC PDdissipationNth2B1;
    CCTK_REAL_VEC PDdissipationNth3B1;
    CCTK_REAL_VEC PDdissipationNth1B2;
    CCTK_REAL_VEC PDdissipationNth2B2;
    CCTK_REAL_VEC PDdissipationNth3B2;
    CCTK_REAL_VEC PDdissipationNth1B3;
    CCTK_REAL_VEC PDdissipationNth2B3;
    CCTK_REAL_VEC PDdissipationNth3B3;
    CCTK_REAL_VEC PDdissipationNth1beta1;
    CCTK_REAL_VEC PDdissipationNth2beta1;
    CCTK_REAL_VEC PDdissipationNth3beta1;
    CCTK_REAL_VEC PDdissipationNth1beta2;
    CCTK_REAL_VEC PDdissipationNth2beta2;
    CCTK_REAL_VEC PDdissipationNth3beta2;
    CCTK_REAL_VEC PDdissipationNth1beta3;
    CCTK_REAL_VEC PDdissipationNth2beta3;
    CCTK_REAL_VEC PDdissipationNth3beta3;
    CCTK_REAL_VEC PDdissipationNth1gt11;
    CCTK_REAL_VEC PDdissipationNth2gt11;
    CCTK_REAL_VEC PDdissipationNth3gt11;
    CCTK_REAL_VEC PDdissipationNth1gt12;
    CCTK_REAL_VEC PDdissipationNth2gt12;
    CCTK_REAL_VEC PDdissipationNth3gt12;
    CCTK_REAL_VEC PDdissipationNth1gt13;
    CCTK_REAL_VEC PDdissipationNth2gt13;
    CCTK_REAL_VEC PDdissipationNth3gt13;
    CCTK_REAL_VEC PDdissipationNth1gt22;
    CCTK_REAL_VEC PDdissipationNth2gt22;
    CCTK_REAL_VEC PDdissipationNth3gt22;
    CCTK_REAL_VEC PDdissipationNth1gt23;
    CCTK_REAL_VEC PDdissipationNth2gt23;
    CCTK_REAL_VEC PDdissipationNth3gt23;
    CCTK_REAL_VEC PDdissipationNth1gt33;
    CCTK_REAL_VEC PDdissipationNth2gt33;
    CCTK_REAL_VEC PDdissipationNth3gt33;
    CCTK_REAL_VEC PDdissipationNth1phi;
    CCTK_REAL_VEC PDdissipationNth2phi;
    CCTK_REAL_VEC PDdissipationNth3phi;
    CCTK_REAL_VEC PDdissipationNth1trK;
    CCTK_REAL_VEC PDdissipationNth2trK;
    CCTK_REAL_VEC PDdissipationNth3trK;
    CCTK_REAL_VEC PDdissipationNth1Xt1;
    CCTK_REAL_VEC PDdissipationNth2Xt1;
    CCTK_REAL_VEC PDdissipationNth3Xt1;
    CCTK_REAL_VEC PDdissipationNth1Xt2;
    CCTK_REAL_VEC PDdissipationNth2Xt2;
    CCTK_REAL_VEC PDdissipationNth3Xt2;
    CCTK_REAL_VEC PDdissipationNth1Xt3;
    CCTK_REAL_VEC PDdissipationNth2Xt3;
    CCTK_REAL_VEC PDdissipationNth3Xt3;
    
    switch(fdOrder)
    {
      case 2:
        PDdissipationNth1A = PDdissipationNthfdOrder21(&A[index]);
        PDdissipationNth2A = PDdissipationNthfdOrder22(&A[index]);
        PDdissipationNth3A = PDdissipationNthfdOrder23(&A[index]);
        PDdissipationNth1alpha = PDdissipationNthfdOrder21(&alpha[index]);
        PDdissipationNth2alpha = PDdissipationNthfdOrder22(&alpha[index]);
        PDdissipationNth3alpha = PDdissipationNthfdOrder23(&alpha[index]);
        PDdissipationNth1At11 = PDdissipationNthfdOrder21(&At11[index]);
        PDdissipationNth2At11 = PDdissipationNthfdOrder22(&At11[index]);
        PDdissipationNth3At11 = PDdissipationNthfdOrder23(&At11[index]);
        PDdissipationNth1At12 = PDdissipationNthfdOrder21(&At12[index]);
        PDdissipationNth2At12 = PDdissipationNthfdOrder22(&At12[index]);
        PDdissipationNth3At12 = PDdissipationNthfdOrder23(&At12[index]);
        PDdissipationNth1At13 = PDdissipationNthfdOrder21(&At13[index]);
        PDdissipationNth2At13 = PDdissipationNthfdOrder22(&At13[index]);
        PDdissipationNth3At13 = PDdissipationNthfdOrder23(&At13[index]);
        PDdissipationNth1At22 = PDdissipationNthfdOrder21(&At22[index]);
        PDdissipationNth2At22 = PDdissipationNthfdOrder22(&At22[index]);
        PDdissipationNth3At22 = PDdissipationNthfdOrder23(&At22[index]);
        PDdissipationNth1At23 = PDdissipationNthfdOrder21(&At23[index]);
        PDdissipationNth2At23 = PDdissipationNthfdOrder22(&At23[index]);
        PDdissipationNth3At23 = PDdissipationNthfdOrder23(&At23[index]);
        PDdissipationNth1At33 = PDdissipationNthfdOrder21(&At33[index]);
        PDdissipationNth2At33 = PDdissipationNthfdOrder22(&At33[index]);
        PDdissipationNth3At33 = PDdissipationNthfdOrder23(&At33[index]);
        PDdissipationNth1B1 = PDdissipationNthfdOrder21(&B1[index]);
        PDdissipationNth2B1 = PDdissipationNthfdOrder22(&B1[index]);
        PDdissipationNth3B1 = PDdissipationNthfdOrder23(&B1[index]);
        PDdissipationNth1B2 = PDdissipationNthfdOrder21(&B2[index]);
        PDdissipationNth2B2 = PDdissipationNthfdOrder22(&B2[index]);
        PDdissipationNth3B2 = PDdissipationNthfdOrder23(&B2[index]);
        PDdissipationNth1B3 = PDdissipationNthfdOrder21(&B3[index]);
        PDdissipationNth2B3 = PDdissipationNthfdOrder22(&B3[index]);
        PDdissipationNth3B3 = PDdissipationNthfdOrder23(&B3[index]);
        PDdissipationNth1beta1 = PDdissipationNthfdOrder21(&beta1[index]);
        PDdissipationNth2beta1 = PDdissipationNthfdOrder22(&beta1[index]);
        PDdissipationNth3beta1 = PDdissipationNthfdOrder23(&beta1[index]);
        PDdissipationNth1beta2 = PDdissipationNthfdOrder21(&beta2[index]);
        PDdissipationNth2beta2 = PDdissipationNthfdOrder22(&beta2[index]);
        PDdissipationNth3beta2 = PDdissipationNthfdOrder23(&beta2[index]);
        PDdissipationNth1beta3 = PDdissipationNthfdOrder21(&beta3[index]);
        PDdissipationNth2beta3 = PDdissipationNthfdOrder22(&beta3[index]);
        PDdissipationNth3beta3 = PDdissipationNthfdOrder23(&beta3[index]);
        PDdissipationNth1gt11 = PDdissipationNthfdOrder21(&gt11[index]);
        PDdissipationNth2gt11 = PDdissipationNthfdOrder22(&gt11[index]);
        PDdissipationNth3gt11 = PDdissipationNthfdOrder23(&gt11[index]);
        PDdissipationNth1gt12 = PDdissipationNthfdOrder21(&gt12[index]);
        PDdissipationNth2gt12 = PDdissipationNthfdOrder22(&gt12[index]);
        PDdissipationNth3gt12 = PDdissipationNthfdOrder23(&gt12[index]);
        PDdissipationNth1gt13 = PDdissipationNthfdOrder21(&gt13[index]);
        PDdissipationNth2gt13 = PDdissipationNthfdOrder22(&gt13[index]);
        PDdissipationNth3gt13 = PDdissipationNthfdOrder23(&gt13[index]);
        PDdissipationNth1gt22 = PDdissipationNthfdOrder21(&gt22[index]);
        PDdissipationNth2gt22 = PDdissipationNthfdOrder22(&gt22[index]);
        PDdissipationNth3gt22 = PDdissipationNthfdOrder23(&gt22[index]);
        PDdissipationNth1gt23 = PDdissipationNthfdOrder21(&gt23[index]);
        PDdissipationNth2gt23 = PDdissipationNthfdOrder22(&gt23[index]);
        PDdissipationNth3gt23 = PDdissipationNthfdOrder23(&gt23[index]);
        PDdissipationNth1gt33 = PDdissipationNthfdOrder21(&gt33[index]);
        PDdissipationNth2gt33 = PDdissipationNthfdOrder22(&gt33[index]);
        PDdissipationNth3gt33 = PDdissipationNthfdOrder23(&gt33[index]);
        PDdissipationNth1phi = PDdissipationNthfdOrder21(&phi[index]);
        PDdissipationNth2phi = PDdissipationNthfdOrder22(&phi[index]);
        PDdissipationNth3phi = PDdissipationNthfdOrder23(&phi[index]);
        PDdissipationNth1trK = PDdissipationNthfdOrder21(&trK[index]);
        PDdissipationNth2trK = PDdissipationNthfdOrder22(&trK[index]);
        PDdissipationNth3trK = PDdissipationNthfdOrder23(&trK[index]);
        PDdissipationNth1Xt1 = PDdissipationNthfdOrder21(&Xt1[index]);
        PDdissipationNth2Xt1 = PDdissipationNthfdOrder22(&Xt1[index]);
        PDdissipationNth3Xt1 = PDdissipationNthfdOrder23(&Xt1[index]);
        PDdissipationNth1Xt2 = PDdissipationNthfdOrder21(&Xt2[index]);
        PDdissipationNth2Xt2 = PDdissipationNthfdOrder22(&Xt2[index]);
        PDdissipationNth3Xt2 = PDdissipationNthfdOrder23(&Xt2[index]);
        PDdissipationNth1Xt3 = PDdissipationNthfdOrder21(&Xt3[index]);
        PDdissipationNth2Xt3 = PDdissipationNthfdOrder22(&Xt3[index]);
        PDdissipationNth3Xt3 = PDdissipationNthfdOrder23(&Xt3[index]);
        break;
      
      case 4:
        PDdissipationNth1A = PDdissipationNthfdOrder41(&A[index]);
        PDdissipationNth2A = PDdissipationNthfdOrder42(&A[index]);
        PDdissipationNth3A = PDdissipationNthfdOrder43(&A[index]);
        PDdissipationNth1alpha = PDdissipationNthfdOrder41(&alpha[index]);
        PDdissipationNth2alpha = PDdissipationNthfdOrder42(&alpha[index]);
        PDdissipationNth3alpha = PDdissipationNthfdOrder43(&alpha[index]);
        PDdissipationNth1At11 = PDdissipationNthfdOrder41(&At11[index]);
        PDdissipationNth2At11 = PDdissipationNthfdOrder42(&At11[index]);
        PDdissipationNth3At11 = PDdissipationNthfdOrder43(&At11[index]);
        PDdissipationNth1At12 = PDdissipationNthfdOrder41(&At12[index]);
        PDdissipationNth2At12 = PDdissipationNthfdOrder42(&At12[index]);
        PDdissipationNth3At12 = PDdissipationNthfdOrder43(&At12[index]);
        PDdissipationNth1At13 = PDdissipationNthfdOrder41(&At13[index]);
        PDdissipationNth2At13 = PDdissipationNthfdOrder42(&At13[index]);
        PDdissipationNth3At13 = PDdissipationNthfdOrder43(&At13[index]);
        PDdissipationNth1At22 = PDdissipationNthfdOrder41(&At22[index]);
        PDdissipationNth2At22 = PDdissipationNthfdOrder42(&At22[index]);
        PDdissipationNth3At22 = PDdissipationNthfdOrder43(&At22[index]);
        PDdissipationNth1At23 = PDdissipationNthfdOrder41(&At23[index]);
        PDdissipationNth2At23 = PDdissipationNthfdOrder42(&At23[index]);
        PDdissipationNth3At23 = PDdissipationNthfdOrder43(&At23[index]);
        PDdissipationNth1At33 = PDdissipationNthfdOrder41(&At33[index]);
        PDdissipationNth2At33 = PDdissipationNthfdOrder42(&At33[index]);
        PDdissipationNth3At33 = PDdissipationNthfdOrder43(&At33[index]);
        PDdissipationNth1B1 = PDdissipationNthfdOrder41(&B1[index]);
        PDdissipationNth2B1 = PDdissipationNthfdOrder42(&B1[index]);
        PDdissipationNth3B1 = PDdissipationNthfdOrder43(&B1[index]);
        PDdissipationNth1B2 = PDdissipationNthfdOrder41(&B2[index]);
        PDdissipationNth2B2 = PDdissipationNthfdOrder42(&B2[index]);
        PDdissipationNth3B2 = PDdissipationNthfdOrder43(&B2[index]);
        PDdissipationNth1B3 = PDdissipationNthfdOrder41(&B3[index]);
        PDdissipationNth2B3 = PDdissipationNthfdOrder42(&B3[index]);
        PDdissipationNth3B3 = PDdissipationNthfdOrder43(&B3[index]);
        PDdissipationNth1beta1 = PDdissipationNthfdOrder41(&beta1[index]);
        PDdissipationNth2beta1 = PDdissipationNthfdOrder42(&beta1[index]);
        PDdissipationNth3beta1 = PDdissipationNthfdOrder43(&beta1[index]);
        PDdissipationNth1beta2 = PDdissipationNthfdOrder41(&beta2[index]);
        PDdissipationNth2beta2 = PDdissipationNthfdOrder42(&beta2[index]);
        PDdissipationNth3beta2 = PDdissipationNthfdOrder43(&beta2[index]);
        PDdissipationNth1beta3 = PDdissipationNthfdOrder41(&beta3[index]);
        PDdissipationNth2beta3 = PDdissipationNthfdOrder42(&beta3[index]);
        PDdissipationNth3beta3 = PDdissipationNthfdOrder43(&beta3[index]);
        PDdissipationNth1gt11 = PDdissipationNthfdOrder41(&gt11[index]);
        PDdissipationNth2gt11 = PDdissipationNthfdOrder42(&gt11[index]);
        PDdissipationNth3gt11 = PDdissipationNthfdOrder43(&gt11[index]);
        PDdissipationNth1gt12 = PDdissipationNthfdOrder41(&gt12[index]);
        PDdissipationNth2gt12 = PDdissipationNthfdOrder42(&gt12[index]);
        PDdissipationNth3gt12 = PDdissipationNthfdOrder43(&gt12[index]);
        PDdissipationNth1gt13 = PDdissipationNthfdOrder41(&gt13[index]);
        PDdissipationNth2gt13 = PDdissipationNthfdOrder42(&gt13[index]);
        PDdissipationNth3gt13 = PDdissipationNthfdOrder43(&gt13[index]);
        PDdissipationNth1gt22 = PDdissipationNthfdOrder41(&gt22[index]);
        PDdissipationNth2gt22 = PDdissipationNthfdOrder42(&gt22[index]);
        PDdissipationNth3gt22 = PDdissipationNthfdOrder43(&gt22[index]);
        PDdissipationNth1gt23 = PDdissipationNthfdOrder41(&gt23[index]);
        PDdissipationNth2gt23 = PDdissipationNthfdOrder42(&gt23[index]);
        PDdissipationNth3gt23 = PDdissipationNthfdOrder43(&gt23[index]);
        PDdissipationNth1gt33 = PDdissipationNthfdOrder41(&gt33[index]);
        PDdissipationNth2gt33 = PDdissipationNthfdOrder42(&gt33[index]);
        PDdissipationNth3gt33 = PDdissipationNthfdOrder43(&gt33[index]);
        PDdissipationNth1phi = PDdissipationNthfdOrder41(&phi[index]);
        PDdissipationNth2phi = PDdissipationNthfdOrder42(&phi[index]);
        PDdissipationNth3phi = PDdissipationNthfdOrder43(&phi[index]);
        PDdissipationNth1trK = PDdissipationNthfdOrder41(&trK[index]);
        PDdissipationNth2trK = PDdissipationNthfdOrder42(&trK[index]);
        PDdissipationNth3trK = PDdissipationNthfdOrder43(&trK[index]);
        PDdissipationNth1Xt1 = PDdissipationNthfdOrder41(&Xt1[index]);
        PDdissipationNth2Xt1 = PDdissipationNthfdOrder42(&Xt1[index]);
        PDdissipationNth3Xt1 = PDdissipationNthfdOrder43(&Xt1[index]);
        PDdissipationNth1Xt2 = PDdissipationNthfdOrder41(&Xt2[index]);
        PDdissipationNth2Xt2 = PDdissipationNthfdOrder42(&Xt2[index]);
        PDdissipationNth3Xt2 = PDdissipationNthfdOrder43(&Xt2[index]);
        PDdissipationNth1Xt3 = PDdissipationNthfdOrder41(&Xt3[index]);
        PDdissipationNth2Xt3 = PDdissipationNthfdOrder42(&Xt3[index]);
        PDdissipationNth3Xt3 = PDdissipationNthfdOrder43(&Xt3[index]);
        break;
      
      case 6:
        PDdissipationNth1A = PDdissipationNthfdOrder61(&A[index]);
        PDdissipationNth2A = PDdissipationNthfdOrder62(&A[index]);
        PDdissipationNth3A = PDdissipationNthfdOrder63(&A[index]);
        PDdissipationNth1alpha = PDdissipationNthfdOrder61(&alpha[index]);
        PDdissipationNth2alpha = PDdissipationNthfdOrder62(&alpha[index]);
        PDdissipationNth3alpha = PDdissipationNthfdOrder63(&alpha[index]);
        PDdissipationNth1At11 = PDdissipationNthfdOrder61(&At11[index]);
        PDdissipationNth2At11 = PDdissipationNthfdOrder62(&At11[index]);
        PDdissipationNth3At11 = PDdissipationNthfdOrder63(&At11[index]);
        PDdissipationNth1At12 = PDdissipationNthfdOrder61(&At12[index]);
        PDdissipationNth2At12 = PDdissipationNthfdOrder62(&At12[index]);
        PDdissipationNth3At12 = PDdissipationNthfdOrder63(&At12[index]);
        PDdissipationNth1At13 = PDdissipationNthfdOrder61(&At13[index]);
        PDdissipationNth2At13 = PDdissipationNthfdOrder62(&At13[index]);
        PDdissipationNth3At13 = PDdissipationNthfdOrder63(&At13[index]);
        PDdissipationNth1At22 = PDdissipationNthfdOrder61(&At22[index]);
        PDdissipationNth2At22 = PDdissipationNthfdOrder62(&At22[index]);
        PDdissipationNth3At22 = PDdissipationNthfdOrder63(&At22[index]);
        PDdissipationNth1At23 = PDdissipationNthfdOrder61(&At23[index]);
        PDdissipationNth2At23 = PDdissipationNthfdOrder62(&At23[index]);
        PDdissipationNth3At23 = PDdissipationNthfdOrder63(&At23[index]);
        PDdissipationNth1At33 = PDdissipationNthfdOrder61(&At33[index]);
        PDdissipationNth2At33 = PDdissipationNthfdOrder62(&At33[index]);
        PDdissipationNth3At33 = PDdissipationNthfdOrder63(&At33[index]);
        PDdissipationNth1B1 = PDdissipationNthfdOrder61(&B1[index]);
        PDdissipationNth2B1 = PDdissipationNthfdOrder62(&B1[index]);
        PDdissipationNth3B1 = PDdissipationNthfdOrder63(&B1[index]);
        PDdissipationNth1B2 = PDdissipationNthfdOrder61(&B2[index]);
        PDdissipationNth2B2 = PDdissipationNthfdOrder62(&B2[index]);
        PDdissipationNth3B2 = PDdissipationNthfdOrder63(&B2[index]);
        PDdissipationNth1B3 = PDdissipationNthfdOrder61(&B3[index]);
        PDdissipationNth2B3 = PDdissipationNthfdOrder62(&B3[index]);
        PDdissipationNth3B3 = PDdissipationNthfdOrder63(&B3[index]);
        PDdissipationNth1beta1 = PDdissipationNthfdOrder61(&beta1[index]);
        PDdissipationNth2beta1 = PDdissipationNthfdOrder62(&beta1[index]);
        PDdissipationNth3beta1 = PDdissipationNthfdOrder63(&beta1[index]);
        PDdissipationNth1beta2 = PDdissipationNthfdOrder61(&beta2[index]);
        PDdissipationNth2beta2 = PDdissipationNthfdOrder62(&beta2[index]);
        PDdissipationNth3beta2 = PDdissipationNthfdOrder63(&beta2[index]);
        PDdissipationNth1beta3 = PDdissipationNthfdOrder61(&beta3[index]);
        PDdissipationNth2beta3 = PDdissipationNthfdOrder62(&beta3[index]);
        PDdissipationNth3beta3 = PDdissipationNthfdOrder63(&beta3[index]);
        PDdissipationNth1gt11 = PDdissipationNthfdOrder61(&gt11[index]);
        PDdissipationNth2gt11 = PDdissipationNthfdOrder62(&gt11[index]);
        PDdissipationNth3gt11 = PDdissipationNthfdOrder63(&gt11[index]);
        PDdissipationNth1gt12 = PDdissipationNthfdOrder61(&gt12[index]);
        PDdissipationNth2gt12 = PDdissipationNthfdOrder62(&gt12[index]);
        PDdissipationNth3gt12 = PDdissipationNthfdOrder63(&gt12[index]);
        PDdissipationNth1gt13 = PDdissipationNthfdOrder61(&gt13[index]);
        PDdissipationNth2gt13 = PDdissipationNthfdOrder62(&gt13[index]);
        PDdissipationNth3gt13 = PDdissipationNthfdOrder63(&gt13[index]);
        PDdissipationNth1gt22 = PDdissipationNthfdOrder61(&gt22[index]);
        PDdissipationNth2gt22 = PDdissipationNthfdOrder62(&gt22[index]);
        PDdissipationNth3gt22 = PDdissipationNthfdOrder63(&gt22[index]);
        PDdissipationNth1gt23 = PDdissipationNthfdOrder61(&gt23[index]);
        PDdissipationNth2gt23 = PDdissipationNthfdOrder62(&gt23[index]);
        PDdissipationNth3gt23 = PDdissipationNthfdOrder63(&gt23[index]);
        PDdissipationNth1gt33 = PDdissipationNthfdOrder61(&gt33[index]);
        PDdissipationNth2gt33 = PDdissipationNthfdOrder62(&gt33[index]);
        PDdissipationNth3gt33 = PDdissipationNthfdOrder63(&gt33[index]);
        PDdissipationNth1phi = PDdissipationNthfdOrder61(&phi[index]);
        PDdissipationNth2phi = PDdissipationNthfdOrder62(&phi[index]);
        PDdissipationNth3phi = PDdissipationNthfdOrder63(&phi[index]);
        PDdissipationNth1trK = PDdissipationNthfdOrder61(&trK[index]);
        PDdissipationNth2trK = PDdissipationNthfdOrder62(&trK[index]);
        PDdissipationNth3trK = PDdissipationNthfdOrder63(&trK[index]);
        PDdissipationNth1Xt1 = PDdissipationNthfdOrder61(&Xt1[index]);
        PDdissipationNth2Xt1 = PDdissipationNthfdOrder62(&Xt1[index]);
        PDdissipationNth3Xt1 = PDdissipationNthfdOrder63(&Xt1[index]);
        PDdissipationNth1Xt2 = PDdissipationNthfdOrder61(&Xt2[index]);
        PDdissipationNth2Xt2 = PDdissipationNthfdOrder62(&Xt2[index]);
        PDdissipationNth3Xt2 = PDdissipationNthfdOrder63(&Xt2[index]);
        PDdissipationNth1Xt3 = PDdissipationNthfdOrder61(&Xt3[index]);
        PDdissipationNth2Xt3 = PDdissipationNthfdOrder62(&Xt3[index]);
        PDdissipationNth3Xt3 = PDdissipationNthfdOrder63(&Xt3[index]);
        break;
      
      case 8:
        PDdissipationNth1A = PDdissipationNthfdOrder81(&A[index]);
        PDdissipationNth2A = PDdissipationNthfdOrder82(&A[index]);
        PDdissipationNth3A = PDdissipationNthfdOrder83(&A[index]);
        PDdissipationNth1alpha = PDdissipationNthfdOrder81(&alpha[index]);
        PDdissipationNth2alpha = PDdissipationNthfdOrder82(&alpha[index]);
        PDdissipationNth3alpha = PDdissipationNthfdOrder83(&alpha[index]);
        PDdissipationNth1At11 = PDdissipationNthfdOrder81(&At11[index]);
        PDdissipationNth2At11 = PDdissipationNthfdOrder82(&At11[index]);
        PDdissipationNth3At11 = PDdissipationNthfdOrder83(&At11[index]);
        PDdissipationNth1At12 = PDdissipationNthfdOrder81(&At12[index]);
        PDdissipationNth2At12 = PDdissipationNthfdOrder82(&At12[index]);
        PDdissipationNth3At12 = PDdissipationNthfdOrder83(&At12[index]);
        PDdissipationNth1At13 = PDdissipationNthfdOrder81(&At13[index]);
        PDdissipationNth2At13 = PDdissipationNthfdOrder82(&At13[index]);
        PDdissipationNth3At13 = PDdissipationNthfdOrder83(&At13[index]);
        PDdissipationNth1At22 = PDdissipationNthfdOrder81(&At22[index]);
        PDdissipationNth2At22 = PDdissipationNthfdOrder82(&At22[index]);
        PDdissipationNth3At22 = PDdissipationNthfdOrder83(&At22[index]);
        PDdissipationNth1At23 = PDdissipationNthfdOrder81(&At23[index]);
        PDdissipationNth2At23 = PDdissipationNthfdOrder82(&At23[index]);
        PDdissipationNth3At23 = PDdissipationNthfdOrder83(&At23[index]);
        PDdissipationNth1At33 = PDdissipationNthfdOrder81(&At33[index]);
        PDdissipationNth2At33 = PDdissipationNthfdOrder82(&At33[index]);
        PDdissipationNth3At33 = PDdissipationNthfdOrder83(&At33[index]);
        PDdissipationNth1B1 = PDdissipationNthfdOrder81(&B1[index]);
        PDdissipationNth2B1 = PDdissipationNthfdOrder82(&B1[index]);
        PDdissipationNth3B1 = PDdissipationNthfdOrder83(&B1[index]);
        PDdissipationNth1B2 = PDdissipationNthfdOrder81(&B2[index]);
        PDdissipationNth2B2 = PDdissipationNthfdOrder82(&B2[index]);
        PDdissipationNth3B2 = PDdissipationNthfdOrder83(&B2[index]);
        PDdissipationNth1B3 = PDdissipationNthfdOrder81(&B3[index]);
        PDdissipationNth2B3 = PDdissipationNthfdOrder82(&B3[index]);
        PDdissipationNth3B3 = PDdissipationNthfdOrder83(&B3[index]);
        PDdissipationNth1beta1 = PDdissipationNthfdOrder81(&beta1[index]);
        PDdissipationNth2beta1 = PDdissipationNthfdOrder82(&beta1[index]);
        PDdissipationNth3beta1 = PDdissipationNthfdOrder83(&beta1[index]);
        PDdissipationNth1beta2 = PDdissipationNthfdOrder81(&beta2[index]);
        PDdissipationNth2beta2 = PDdissipationNthfdOrder82(&beta2[index]);
        PDdissipationNth3beta2 = PDdissipationNthfdOrder83(&beta2[index]);
        PDdissipationNth1beta3 = PDdissipationNthfdOrder81(&beta3[index]);
        PDdissipationNth2beta3 = PDdissipationNthfdOrder82(&beta3[index]);
        PDdissipationNth3beta3 = PDdissipationNthfdOrder83(&beta3[index]);
        PDdissipationNth1gt11 = PDdissipationNthfdOrder81(&gt11[index]);
        PDdissipationNth2gt11 = PDdissipationNthfdOrder82(&gt11[index]);
        PDdissipationNth3gt11 = PDdissipationNthfdOrder83(&gt11[index]);
        PDdissipationNth1gt12 = PDdissipationNthfdOrder81(&gt12[index]);
        PDdissipationNth2gt12 = PDdissipationNthfdOrder82(&gt12[index]);
        PDdissipationNth3gt12 = PDdissipationNthfdOrder83(&gt12[index]);
        PDdissipationNth1gt13 = PDdissipationNthfdOrder81(&gt13[index]);
        PDdissipationNth2gt13 = PDdissipationNthfdOrder82(&gt13[index]);
        PDdissipationNth3gt13 = PDdissipationNthfdOrder83(&gt13[index]);
        PDdissipationNth1gt22 = PDdissipationNthfdOrder81(&gt22[index]);
        PDdissipationNth2gt22 = PDdissipationNthfdOrder82(&gt22[index]);
        PDdissipationNth3gt22 = PDdissipationNthfdOrder83(&gt22[index]);
        PDdissipationNth1gt23 = PDdissipationNthfdOrder81(&gt23[index]);
        PDdissipationNth2gt23 = PDdissipationNthfdOrder82(&gt23[index]);
        PDdissipationNth3gt23 = PDdissipationNthfdOrder83(&gt23[index]);
        PDdissipationNth1gt33 = PDdissipationNthfdOrder81(&gt33[index]);
        PDdissipationNth2gt33 = PDdissipationNthfdOrder82(&gt33[index]);
        PDdissipationNth3gt33 = PDdissipationNthfdOrder83(&gt33[index]);
        PDdissipationNth1phi = PDdissipationNthfdOrder81(&phi[index]);
        PDdissipationNth2phi = PDdissipationNthfdOrder82(&phi[index]);
        PDdissipationNth3phi = PDdissipationNthfdOrder83(&phi[index]);
        PDdissipationNth1trK = PDdissipationNthfdOrder81(&trK[index]);
        PDdissipationNth2trK = PDdissipationNthfdOrder82(&trK[index]);
        PDdissipationNth3trK = PDdissipationNthfdOrder83(&trK[index]);
        PDdissipationNth1Xt1 = PDdissipationNthfdOrder81(&Xt1[index]);
        PDdissipationNth2Xt1 = PDdissipationNthfdOrder82(&Xt1[index]);
        PDdissipationNth3Xt1 = PDdissipationNthfdOrder83(&Xt1[index]);
        PDdissipationNth1Xt2 = PDdissipationNthfdOrder81(&Xt2[index]);
        PDdissipationNth2Xt2 = PDdissipationNthfdOrder82(&Xt2[index]);
        PDdissipationNth3Xt2 = PDdissipationNthfdOrder83(&Xt2[index]);
        PDdissipationNth1Xt3 = PDdissipationNthfdOrder81(&Xt3[index]);
        PDdissipationNth2Xt3 = PDdissipationNthfdOrder82(&Xt3[index]);
        PDdissipationNth3Xt3 = PDdissipationNthfdOrder83(&Xt3[index]);
        break;
    }
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC JacPDdissipationNth1A;
    CCTK_REAL_VEC JacPDdissipationNth1alpha;
    CCTK_REAL_VEC JacPDdissipationNth1At11;
    CCTK_REAL_VEC JacPDdissipationNth1At12;
    CCTK_REAL_VEC JacPDdissipationNth1At13;
    CCTK_REAL_VEC JacPDdissipationNth1At22;
    CCTK_REAL_VEC JacPDdissipationNth1At23;
    CCTK_REAL_VEC JacPDdissipationNth1At33;
    CCTK_REAL_VEC JacPDdissipationNth1B1;
    CCTK_REAL_VEC JacPDdissipationNth1B2;
    CCTK_REAL_VEC JacPDdissipationNth1B3;
    CCTK_REAL_VEC JacPDdissipationNth1beta1;
    CCTK_REAL_VEC JacPDdissipationNth1beta2;
    CCTK_REAL_VEC JacPDdissipationNth1beta3;
    CCTK_REAL_VEC JacPDdissipationNth1gt11;
    CCTK_REAL_VEC JacPDdissipationNth1gt12;
    CCTK_REAL_VEC JacPDdissipationNth1gt13;
    CCTK_REAL_VEC JacPDdissipationNth1gt22;
    CCTK_REAL_VEC JacPDdissipationNth1gt23;
    CCTK_REAL_VEC JacPDdissipationNth1gt33;
    CCTK_REAL_VEC JacPDdissipationNth1phi;
    CCTK_REAL_VEC JacPDdissipationNth1trK;
    CCTK_REAL_VEC JacPDdissipationNth1Xt1;
    CCTK_REAL_VEC JacPDdissipationNth1Xt2;
    CCTK_REAL_VEC JacPDdissipationNth1Xt3;
    CCTK_REAL_VEC JacPDdissipationNth2A;
    CCTK_REAL_VEC JacPDdissipationNth2alpha;
    CCTK_REAL_VEC JacPDdissipationNth2At11;
    CCTK_REAL_VEC JacPDdissipationNth2At12;
    CCTK_REAL_VEC JacPDdissipationNth2At13;
    CCTK_REAL_VEC JacPDdissipationNth2At22;
    CCTK_REAL_VEC JacPDdissipationNth2At23;
    CCTK_REAL_VEC JacPDdissipationNth2At33;
    CCTK_REAL_VEC JacPDdissipationNth2B1;
    CCTK_REAL_VEC JacPDdissipationNth2B2;
    CCTK_REAL_VEC JacPDdissipationNth2B3;
    CCTK_REAL_VEC JacPDdissipationNth2beta1;
    CCTK_REAL_VEC JacPDdissipationNth2beta2;
    CCTK_REAL_VEC JacPDdissipationNth2beta3;
    CCTK_REAL_VEC JacPDdissipationNth2gt11;
    CCTK_REAL_VEC JacPDdissipationNth2gt12;
    CCTK_REAL_VEC JacPDdissipationNth2gt13;
    CCTK_REAL_VEC JacPDdissipationNth2gt22;
    CCTK_REAL_VEC JacPDdissipationNth2gt23;
    CCTK_REAL_VEC JacPDdissipationNth2gt33;
    CCTK_REAL_VEC JacPDdissipationNth2phi;
    CCTK_REAL_VEC JacPDdissipationNth2trK;
    CCTK_REAL_VEC JacPDdissipationNth2Xt1;
    CCTK_REAL_VEC JacPDdissipationNth2Xt2;
    CCTK_REAL_VEC JacPDdissipationNth2Xt3;
    CCTK_REAL_VEC JacPDdissipationNth3A;
    CCTK_REAL_VEC JacPDdissipationNth3alpha;
    CCTK_REAL_VEC JacPDdissipationNth3At11;
    CCTK_REAL_VEC JacPDdissipationNth3At12;
    CCTK_REAL_VEC JacPDdissipationNth3At13;
    CCTK_REAL_VEC JacPDdissipationNth3At22;
    CCTK_REAL_VEC JacPDdissipationNth3At23;
    CCTK_REAL_VEC JacPDdissipationNth3At33;
    CCTK_REAL_VEC JacPDdissipationNth3B1;
    CCTK_REAL_VEC JacPDdissipationNth3B2;
    CCTK_REAL_VEC JacPDdissipationNth3B3;
    CCTK_REAL_VEC JacPDdissipationNth3beta1;
    CCTK_REAL_VEC JacPDdissipationNth3beta2;
    CCTK_REAL_VEC JacPDdissipationNth3beta3;
    CCTK_REAL_VEC JacPDdissipationNth3gt11;
    CCTK_REAL_VEC JacPDdissipationNth3gt12;
    CCTK_REAL_VEC JacPDdissipationNth3gt13;
    CCTK_REAL_VEC JacPDdissipationNth3gt22;
    CCTK_REAL_VEC JacPDdissipationNth3gt23;
    CCTK_REAL_VEC JacPDdissipationNth3gt33;
    CCTK_REAL_VEC JacPDdissipationNth3phi;
    CCTK_REAL_VEC JacPDdissipationNth3trK;
    CCTK_REAL_VEC JacPDdissipationNth3Xt1;
    CCTK_REAL_VEC JacPDdissipationNth3Xt2;
    CCTK_REAL_VEC JacPDdissipationNth3Xt3;
    
    if (use_jacobian)
    {
      JacPDdissipationNth1A = 
        kmadd(J11L,PDdissipationNth1A,kmadd(J21L,PDdissipationNth2A,kmul(J31L,PDdissipationNth3A)));
      
      JacPDdissipationNth1alpha = 
        kmadd(J11L,PDdissipationNth1alpha,kmadd(J21L,PDdissipationNth2alpha,kmul(J31L,PDdissipationNth3alpha)));
      
      JacPDdissipationNth1At11 = 
        kmadd(J11L,PDdissipationNth1At11,kmadd(J21L,PDdissipationNth2At11,kmul(J31L,PDdissipationNth3At11)));
      
      JacPDdissipationNth1At12 = 
        kmadd(J11L,PDdissipationNth1At12,kmadd(J21L,PDdissipationNth2At12,kmul(J31L,PDdissipationNth3At12)));
      
      JacPDdissipationNth1At13 = 
        kmadd(J11L,PDdissipationNth1At13,kmadd(J21L,PDdissipationNth2At13,kmul(J31L,PDdissipationNth3At13)));
      
      JacPDdissipationNth1At22 = 
        kmadd(J11L,PDdissipationNth1At22,kmadd(J21L,PDdissipationNth2At22,kmul(J31L,PDdissipationNth3At22)));
      
      JacPDdissipationNth1At23 = 
        kmadd(J11L,PDdissipationNth1At23,kmadd(J21L,PDdissipationNth2At23,kmul(J31L,PDdissipationNth3At23)));
      
      JacPDdissipationNth1At33 = 
        kmadd(J11L,PDdissipationNth1At33,kmadd(J21L,PDdissipationNth2At33,kmul(J31L,PDdissipationNth3At33)));
      
      JacPDdissipationNth1B1 = 
        kmadd(J11L,PDdissipationNth1B1,kmadd(J21L,PDdissipationNth2B1,kmul(J31L,PDdissipationNth3B1)));
      
      JacPDdissipationNth1B2 = 
        kmadd(J11L,PDdissipationNth1B2,kmadd(J21L,PDdissipationNth2B2,kmul(J31L,PDdissipationNth3B2)));
      
      JacPDdissipationNth1B3 = 
        kmadd(J11L,PDdissipationNth1B3,kmadd(J21L,PDdissipationNth2B3,kmul(J31L,PDdissipationNth3B3)));
      
      JacPDdissipationNth1beta1 = 
        kmadd(J11L,PDdissipationNth1beta1,kmadd(J21L,PDdissipationNth2beta1,kmul(J31L,PDdissipationNth3beta1)));
      
      JacPDdissipationNth1beta2 = 
        kmadd(J11L,PDdissipationNth1beta2,kmadd(J21L,PDdissipationNth2beta2,kmul(J31L,PDdissipationNth3beta2)));
      
      JacPDdissipationNth1beta3 = 
        kmadd(J11L,PDdissipationNth1beta3,kmadd(J21L,PDdissipationNth2beta3,kmul(J31L,PDdissipationNth3beta3)));
      
      JacPDdissipationNth1gt11 = 
        kmadd(J11L,PDdissipationNth1gt11,kmadd(J21L,PDdissipationNth2gt11,kmul(J31L,PDdissipationNth3gt11)));
      
      JacPDdissipationNth1gt12 = 
        kmadd(J11L,PDdissipationNth1gt12,kmadd(J21L,PDdissipationNth2gt12,kmul(J31L,PDdissipationNth3gt12)));
      
      JacPDdissipationNth1gt13 = 
        kmadd(J11L,PDdissipationNth1gt13,kmadd(J21L,PDdissipationNth2gt13,kmul(J31L,PDdissipationNth3gt13)));
      
      JacPDdissipationNth1gt22 = 
        kmadd(J11L,PDdissipationNth1gt22,kmadd(J21L,PDdissipationNth2gt22,kmul(J31L,PDdissipationNth3gt22)));
      
      JacPDdissipationNth1gt23 = 
        kmadd(J11L,PDdissipationNth1gt23,kmadd(J21L,PDdissipationNth2gt23,kmul(J31L,PDdissipationNth3gt23)));
      
      JacPDdissipationNth1gt33 = 
        kmadd(J11L,PDdissipationNth1gt33,kmadd(J21L,PDdissipationNth2gt33,kmul(J31L,PDdissipationNth3gt33)));
      
      JacPDdissipationNth1phi = 
        kmadd(J11L,PDdissipationNth1phi,kmadd(J21L,PDdissipationNth2phi,kmul(J31L,PDdissipationNth3phi)));
      
      JacPDdissipationNth1trK = 
        kmadd(J11L,PDdissipationNth1trK,kmadd(J21L,PDdissipationNth2trK,kmul(J31L,PDdissipationNth3trK)));
      
      JacPDdissipationNth1Xt1 = 
        kmadd(J11L,PDdissipationNth1Xt1,kmadd(J21L,PDdissipationNth2Xt1,kmul(J31L,PDdissipationNth3Xt1)));
      
      JacPDdissipationNth1Xt2 = 
        kmadd(J11L,PDdissipationNth1Xt2,kmadd(J21L,PDdissipationNth2Xt2,kmul(J31L,PDdissipationNth3Xt2)));
      
      JacPDdissipationNth1Xt3 = 
        kmadd(J11L,PDdissipationNth1Xt3,kmadd(J21L,PDdissipationNth2Xt3,kmul(J31L,PDdissipationNth3Xt3)));
      
      JacPDdissipationNth2A = 
        kmadd(J12L,PDdissipationNth1A,kmadd(J22L,PDdissipationNth2A,kmul(J32L,PDdissipationNth3A)));
      
      JacPDdissipationNth2alpha = 
        kmadd(J12L,PDdissipationNth1alpha,kmadd(J22L,PDdissipationNth2alpha,kmul(J32L,PDdissipationNth3alpha)));
      
      JacPDdissipationNth2At11 = 
        kmadd(J12L,PDdissipationNth1At11,kmadd(J22L,PDdissipationNth2At11,kmul(J32L,PDdissipationNth3At11)));
      
      JacPDdissipationNth2At12 = 
        kmadd(J12L,PDdissipationNth1At12,kmadd(J22L,PDdissipationNth2At12,kmul(J32L,PDdissipationNth3At12)));
      
      JacPDdissipationNth2At13 = 
        kmadd(J12L,PDdissipationNth1At13,kmadd(J22L,PDdissipationNth2At13,kmul(J32L,PDdissipationNth3At13)));
      
      JacPDdissipationNth2At22 = 
        kmadd(J12L,PDdissipationNth1At22,kmadd(J22L,PDdissipationNth2At22,kmul(J32L,PDdissipationNth3At22)));
      
      JacPDdissipationNth2At23 = 
        kmadd(J12L,PDdissipationNth1At23,kmadd(J22L,PDdissipationNth2At23,kmul(J32L,PDdissipationNth3At23)));
      
      JacPDdissipationNth2At33 = 
        kmadd(J12L,PDdissipationNth1At33,kmadd(J22L,PDdissipationNth2At33,kmul(J32L,PDdissipationNth3At33)));
      
      JacPDdissipationNth2B1 = 
        kmadd(J12L,PDdissipationNth1B1,kmadd(J22L,PDdissipationNth2B1,kmul(J32L,PDdissipationNth3B1)));
      
      JacPDdissipationNth2B2 = 
        kmadd(J12L,PDdissipationNth1B2,kmadd(J22L,PDdissipationNth2B2,kmul(J32L,PDdissipationNth3B2)));
      
      JacPDdissipationNth2B3 = 
        kmadd(J12L,PDdissipationNth1B3,kmadd(J22L,PDdissipationNth2B3,kmul(J32L,PDdissipationNth3B3)));
      
      JacPDdissipationNth2beta1 = 
        kmadd(J12L,PDdissipationNth1beta1,kmadd(J22L,PDdissipationNth2beta1,kmul(J32L,PDdissipationNth3beta1)));
      
      JacPDdissipationNth2beta2 = 
        kmadd(J12L,PDdissipationNth1beta2,kmadd(J22L,PDdissipationNth2beta2,kmul(J32L,PDdissipationNth3beta2)));
      
      JacPDdissipationNth2beta3 = 
        kmadd(J12L,PDdissipationNth1beta3,kmadd(J22L,PDdissipationNth2beta3,kmul(J32L,PDdissipationNth3beta3)));
      
      JacPDdissipationNth2gt11 = 
        kmadd(J12L,PDdissipationNth1gt11,kmadd(J22L,PDdissipationNth2gt11,kmul(J32L,PDdissipationNth3gt11)));
      
      JacPDdissipationNth2gt12 = 
        kmadd(J12L,PDdissipationNth1gt12,kmadd(J22L,PDdissipationNth2gt12,kmul(J32L,PDdissipationNth3gt12)));
      
      JacPDdissipationNth2gt13 = 
        kmadd(J12L,PDdissipationNth1gt13,kmadd(J22L,PDdissipationNth2gt13,kmul(J32L,PDdissipationNth3gt13)));
      
      JacPDdissipationNth2gt22 = 
        kmadd(J12L,PDdissipationNth1gt22,kmadd(J22L,PDdissipationNth2gt22,kmul(J32L,PDdissipationNth3gt22)));
      
      JacPDdissipationNth2gt23 = 
        kmadd(J12L,PDdissipationNth1gt23,kmadd(J22L,PDdissipationNth2gt23,kmul(J32L,PDdissipationNth3gt23)));
      
      JacPDdissipationNth2gt33 = 
        kmadd(J12L,PDdissipationNth1gt33,kmadd(J22L,PDdissipationNth2gt33,kmul(J32L,PDdissipationNth3gt33)));
      
      JacPDdissipationNth2phi = 
        kmadd(J12L,PDdissipationNth1phi,kmadd(J22L,PDdissipationNth2phi,kmul(J32L,PDdissipationNth3phi)));
      
      JacPDdissipationNth2trK = 
        kmadd(J12L,PDdissipationNth1trK,kmadd(J22L,PDdissipationNth2trK,kmul(J32L,PDdissipationNth3trK)));
      
      JacPDdissipationNth2Xt1 = 
        kmadd(J12L,PDdissipationNth1Xt1,kmadd(J22L,PDdissipationNth2Xt1,kmul(J32L,PDdissipationNth3Xt1)));
      
      JacPDdissipationNth2Xt2 = 
        kmadd(J12L,PDdissipationNth1Xt2,kmadd(J22L,PDdissipationNth2Xt2,kmul(J32L,PDdissipationNth3Xt2)));
      
      JacPDdissipationNth2Xt3 = 
        kmadd(J12L,PDdissipationNth1Xt3,kmadd(J22L,PDdissipationNth2Xt3,kmul(J32L,PDdissipationNth3Xt3)));
      
      JacPDdissipationNth3A = 
        kmadd(J13L,PDdissipationNth1A,kmadd(J23L,PDdissipationNth2A,kmul(J33L,PDdissipationNth3A)));
      
      JacPDdissipationNth3alpha = 
        kmadd(J13L,PDdissipationNth1alpha,kmadd(J23L,PDdissipationNth2alpha,kmul(J33L,PDdissipationNth3alpha)));
      
      JacPDdissipationNth3At11 = 
        kmadd(J13L,PDdissipationNth1At11,kmadd(J23L,PDdissipationNth2At11,kmul(J33L,PDdissipationNth3At11)));
      
      JacPDdissipationNth3At12 = 
        kmadd(J13L,PDdissipationNth1At12,kmadd(J23L,PDdissipationNth2At12,kmul(J33L,PDdissipationNth3At12)));
      
      JacPDdissipationNth3At13 = 
        kmadd(J13L,PDdissipationNth1At13,kmadd(J23L,PDdissipationNth2At13,kmul(J33L,PDdissipationNth3At13)));
      
      JacPDdissipationNth3At22 = 
        kmadd(J13L,PDdissipationNth1At22,kmadd(J23L,PDdissipationNth2At22,kmul(J33L,PDdissipationNth3At22)));
      
      JacPDdissipationNth3At23 = 
        kmadd(J13L,PDdissipationNth1At23,kmadd(J23L,PDdissipationNth2At23,kmul(J33L,PDdissipationNth3At23)));
      
      JacPDdissipationNth3At33 = 
        kmadd(J13L,PDdissipationNth1At33,kmadd(J23L,PDdissipationNth2At33,kmul(J33L,PDdissipationNth3At33)));
      
      JacPDdissipationNth3B1 = 
        kmadd(J13L,PDdissipationNth1B1,kmadd(J23L,PDdissipationNth2B1,kmul(J33L,PDdissipationNth3B1)));
      
      JacPDdissipationNth3B2 = 
        kmadd(J13L,PDdissipationNth1B2,kmadd(J23L,PDdissipationNth2B2,kmul(J33L,PDdissipationNth3B2)));
      
      JacPDdissipationNth3B3 = 
        kmadd(J13L,PDdissipationNth1B3,kmadd(J23L,PDdissipationNth2B3,kmul(J33L,PDdissipationNth3B3)));
      
      JacPDdissipationNth3beta1 = 
        kmadd(J13L,PDdissipationNth1beta1,kmadd(J23L,PDdissipationNth2beta1,kmul(J33L,PDdissipationNth3beta1)));
      
      JacPDdissipationNth3beta2 = 
        kmadd(J13L,PDdissipationNth1beta2,kmadd(J23L,PDdissipationNth2beta2,kmul(J33L,PDdissipationNth3beta2)));
      
      JacPDdissipationNth3beta3 = 
        kmadd(J13L,PDdissipationNth1beta3,kmadd(J23L,PDdissipationNth2beta3,kmul(J33L,PDdissipationNth3beta3)));
      
      JacPDdissipationNth3gt11 = 
        kmadd(J13L,PDdissipationNth1gt11,kmadd(J23L,PDdissipationNth2gt11,kmul(J33L,PDdissipationNth3gt11)));
      
      JacPDdissipationNth3gt12 = 
        kmadd(J13L,PDdissipationNth1gt12,kmadd(J23L,PDdissipationNth2gt12,kmul(J33L,PDdissipationNth3gt12)));
      
      JacPDdissipationNth3gt13 = 
        kmadd(J13L,PDdissipationNth1gt13,kmadd(J23L,PDdissipationNth2gt13,kmul(J33L,PDdissipationNth3gt13)));
      
      JacPDdissipationNth3gt22 = 
        kmadd(J13L,PDdissipationNth1gt22,kmadd(J23L,PDdissipationNth2gt22,kmul(J33L,PDdissipationNth3gt22)));
      
      JacPDdissipationNth3gt23 = 
        kmadd(J13L,PDdissipationNth1gt23,kmadd(J23L,PDdissipationNth2gt23,kmul(J33L,PDdissipationNth3gt23)));
      
      JacPDdissipationNth3gt33 = 
        kmadd(J13L,PDdissipationNth1gt33,kmadd(J23L,PDdissipationNth2gt33,kmul(J33L,PDdissipationNth3gt33)));
      
      JacPDdissipationNth3phi = 
        kmadd(J13L,PDdissipationNth1phi,kmadd(J23L,PDdissipationNth2phi,kmul(J33L,PDdissipationNth3phi)));
      
      JacPDdissipationNth3trK = 
        kmadd(J13L,PDdissipationNth1trK,kmadd(J23L,PDdissipationNth2trK,kmul(J33L,PDdissipationNth3trK)));
      
      JacPDdissipationNth3Xt1 = 
        kmadd(J13L,PDdissipationNth1Xt1,kmadd(J23L,PDdissipationNth2Xt1,kmul(J33L,PDdissipationNth3Xt1)));
      
      JacPDdissipationNth3Xt2 = 
        kmadd(J13L,PDdissipationNth1Xt2,kmadd(J23L,PDdissipationNth2Xt2,kmul(J33L,PDdissipationNth3Xt2)));
      
      JacPDdissipationNth3Xt3 = 
        kmadd(J13L,PDdissipationNth1Xt3,kmadd(J23L,PDdissipationNth2Xt3,kmul(J33L,PDdissipationNth3Xt3)));
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
    
    CCTK_REAL_VEC epsdiss1 = ToReal(EpsDiss);
    
    CCTK_REAL_VEC epsdiss2 = ToReal(EpsDiss);
    
    CCTK_REAL_VEC epsdiss3 = ToReal(EpsDiss);
    
    phirhsL = 
      kadd(phirhsL,kmadd(epsdiss1,JacPDdissipationNth1phi,kmadd(epsdiss2,JacPDdissipationNth2phi,kmul(epsdiss3,JacPDdissipationNth3phi))));
    
    gt11rhsL = 
      kadd(gt11rhsL,kmadd(epsdiss1,JacPDdissipationNth1gt11,kmadd(epsdiss2,JacPDdissipationNth2gt11,kmul(epsdiss3,JacPDdissipationNth3gt11))));
    
    gt12rhsL = 
      kadd(gt12rhsL,kmadd(epsdiss1,JacPDdissipationNth1gt12,kmadd(epsdiss2,JacPDdissipationNth2gt12,kmul(epsdiss3,JacPDdissipationNth3gt12))));
    
    gt13rhsL = 
      kadd(gt13rhsL,kmadd(epsdiss1,JacPDdissipationNth1gt13,kmadd(epsdiss2,JacPDdissipationNth2gt13,kmul(epsdiss3,JacPDdissipationNth3gt13))));
    
    gt22rhsL = 
      kadd(gt22rhsL,kmadd(epsdiss1,JacPDdissipationNth1gt22,kmadd(epsdiss2,JacPDdissipationNth2gt22,kmul(epsdiss3,JacPDdissipationNth3gt22))));
    
    gt23rhsL = 
      kadd(gt23rhsL,kmadd(epsdiss1,JacPDdissipationNth1gt23,kmadd(epsdiss2,JacPDdissipationNth2gt23,kmul(epsdiss3,JacPDdissipationNth3gt23))));
    
    gt33rhsL = 
      kadd(gt33rhsL,kmadd(epsdiss1,JacPDdissipationNth1gt33,kmadd(epsdiss2,JacPDdissipationNth2gt33,kmul(epsdiss3,JacPDdissipationNth3gt33))));
    
    Xt1rhsL = 
      kadd(Xt1rhsL,kmadd(epsdiss1,JacPDdissipationNth1Xt1,kmadd(epsdiss2,JacPDdissipationNth2Xt1,kmul(epsdiss3,JacPDdissipationNth3Xt1))));
    
    Xt2rhsL = 
      kadd(Xt2rhsL,kmadd(epsdiss1,JacPDdissipationNth1Xt2,kmadd(epsdiss2,JacPDdissipationNth2Xt2,kmul(epsdiss3,JacPDdissipationNth3Xt2))));
    
    Xt3rhsL = 
      kadd(Xt3rhsL,kmadd(epsdiss1,JacPDdissipationNth1Xt3,kmadd(epsdiss2,JacPDdissipationNth2Xt3,kmul(epsdiss3,JacPDdissipationNth3Xt3))));
    
    trKrhsL = 
      kadd(trKrhsL,kmadd(epsdiss1,JacPDdissipationNth1trK,kmadd(epsdiss2,JacPDdissipationNth2trK,kmul(epsdiss3,JacPDdissipationNth3trK))));
    
    At11rhsL = 
      kadd(At11rhsL,kmadd(epsdiss1,JacPDdissipationNth1At11,kmadd(epsdiss2,JacPDdissipationNth2At11,kmul(epsdiss3,JacPDdissipationNth3At11))));
    
    At12rhsL = 
      kadd(At12rhsL,kmadd(epsdiss1,JacPDdissipationNth1At12,kmadd(epsdiss2,JacPDdissipationNth2At12,kmul(epsdiss3,JacPDdissipationNth3At12))));
    
    At13rhsL = 
      kadd(At13rhsL,kmadd(epsdiss1,JacPDdissipationNth1At13,kmadd(epsdiss2,JacPDdissipationNth2At13,kmul(epsdiss3,JacPDdissipationNth3At13))));
    
    At22rhsL = 
      kadd(At22rhsL,kmadd(epsdiss1,JacPDdissipationNth1At22,kmadd(epsdiss2,JacPDdissipationNth2At22,kmul(epsdiss3,JacPDdissipationNth3At22))));
    
    At23rhsL = 
      kadd(At23rhsL,kmadd(epsdiss1,JacPDdissipationNth1At23,kmadd(epsdiss2,JacPDdissipationNth2At23,kmul(epsdiss3,JacPDdissipationNth3At23))));
    
    At33rhsL = 
      kadd(At33rhsL,kmadd(epsdiss1,JacPDdissipationNth1At33,kmadd(epsdiss2,JacPDdissipationNth2At33,kmul(epsdiss3,JacPDdissipationNth3At33))));
    
    alpharhsL = 
      kadd(alpharhsL,kmadd(epsdiss1,JacPDdissipationNth1alpha,kmadd(epsdiss2,JacPDdissipationNth2alpha,kmul(epsdiss3,JacPDdissipationNth3alpha))));
    
    ArhsL = 
      kadd(ArhsL,kmadd(epsdiss1,JacPDdissipationNth1A,kmadd(epsdiss2,JacPDdissipationNth2A,kmul(epsdiss3,JacPDdissipationNth3A))));
    
    beta1rhsL = 
      kadd(beta1rhsL,kmadd(epsdiss1,JacPDdissipationNth1beta1,kmadd(epsdiss2,JacPDdissipationNth2beta1,kmul(epsdiss3,JacPDdissipationNth3beta1))));
    
    beta2rhsL = 
      kadd(beta2rhsL,kmadd(epsdiss1,JacPDdissipationNth1beta2,kmadd(epsdiss2,JacPDdissipationNth2beta2,kmul(epsdiss3,JacPDdissipationNth3beta2))));
    
    beta3rhsL = 
      kadd(beta3rhsL,kmadd(epsdiss1,JacPDdissipationNth1beta3,kmadd(epsdiss2,JacPDdissipationNth2beta3,kmul(epsdiss3,JacPDdissipationNth3beta3))));
    
    B1rhsL = 
      kadd(B1rhsL,kmadd(epsdiss1,JacPDdissipationNth1B1,kmadd(epsdiss2,JacPDdissipationNth2B1,kmul(epsdiss3,JacPDdissipationNth3B1))));
    
    B2rhsL = 
      kadd(B2rhsL,kmadd(epsdiss1,JacPDdissipationNth1B2,kmadd(epsdiss2,JacPDdissipationNth2B2,kmul(epsdiss3,JacPDdissipationNth3B2))));
    
    B3rhsL = 
      kadd(B3rhsL,kmadd(epsdiss1,JacPDdissipationNth1B3,kmadd(epsdiss2,JacPDdissipationNth2B3,kmul(epsdiss3,JacPDdissipationNth3B3))));
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(alpharhs[index],alpharhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Arhs[index],ArhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At11rhs[index],At11rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At12rhs[index],At12rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At13rhs[index],At13rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At22rhs[index],At22rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At23rhs[index],At23rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At33rhs[index],At33rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B1rhs[index],B1rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B2rhs[index],B2rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B3rhs[index],B3rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta1rhs[index],beta1rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta2rhs[index],beta2rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta3rhs[index],beta3rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt11rhs[index],gt11rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt12rhs[index],gt12rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt13rhs[index],gt13rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt22rhs[index],gt22rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt23rhs[index],gt23rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt33rhs[index],gt33rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(phirhs[index],phirhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(trKrhs[index],trKrhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt1rhs[index],Xt1rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt2rhs[index],Xt2rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt3rhs[index],Xt3rhsL,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(alpharhs[index],alpharhsL,elt_count);
      vec_store_nta_partial_hi(Arhs[index],ArhsL,elt_count);
      vec_store_nta_partial_hi(At11rhs[index],At11rhsL,elt_count);
      vec_store_nta_partial_hi(At12rhs[index],At12rhsL,elt_count);
      vec_store_nta_partial_hi(At13rhs[index],At13rhsL,elt_count);
      vec_store_nta_partial_hi(At22rhs[index],At22rhsL,elt_count);
      vec_store_nta_partial_hi(At23rhs[index],At23rhsL,elt_count);
      vec_store_nta_partial_hi(At33rhs[index],At33rhsL,elt_count);
      vec_store_nta_partial_hi(B1rhs[index],B1rhsL,elt_count);
      vec_store_nta_partial_hi(B2rhs[index],B2rhsL,elt_count);
      vec_store_nta_partial_hi(B3rhs[index],B3rhsL,elt_count);
      vec_store_nta_partial_hi(beta1rhs[index],beta1rhsL,elt_count);
      vec_store_nta_partial_hi(beta2rhs[index],beta2rhsL,elt_count);
      vec_store_nta_partial_hi(beta3rhs[index],beta3rhsL,elt_count);
      vec_store_nta_partial_hi(gt11rhs[index],gt11rhsL,elt_count);
      vec_store_nta_partial_hi(gt12rhs[index],gt12rhsL,elt_count);
      vec_store_nta_partial_hi(gt13rhs[index],gt13rhsL,elt_count);
      vec_store_nta_partial_hi(gt22rhs[index],gt22rhsL,elt_count);
      vec_store_nta_partial_hi(gt23rhs[index],gt23rhsL,elt_count);
      vec_store_nta_partial_hi(gt33rhs[index],gt33rhsL,elt_count);
      vec_store_nta_partial_hi(phirhs[index],phirhsL,elt_count);
      vec_store_nta_partial_hi(trKrhs[index],trKrhsL,elt_count);
      vec_store_nta_partial_hi(Xt1rhs[index],Xt1rhsL,elt_count);
      vec_store_nta_partial_hi(Xt2rhs[index],Xt2rhsL,elt_count);
      vec_store_nta_partial_hi(Xt3rhs[index],Xt3rhsL,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(alpharhs[index],alpharhsL,elt_count);
      vec_store_nta_partial_lo(Arhs[index],ArhsL,elt_count);
      vec_store_nta_partial_lo(At11rhs[index],At11rhsL,elt_count);
      vec_store_nta_partial_lo(At12rhs[index],At12rhsL,elt_count);
      vec_store_nta_partial_lo(At13rhs[index],At13rhsL,elt_count);
      vec_store_nta_partial_lo(At22rhs[index],At22rhsL,elt_count);
      vec_store_nta_partial_lo(At23rhs[index],At23rhsL,elt_count);
      vec_store_nta_partial_lo(At33rhs[index],At33rhsL,elt_count);
      vec_store_nta_partial_lo(B1rhs[index],B1rhsL,elt_count);
      vec_store_nta_partial_lo(B2rhs[index],B2rhsL,elt_count);
      vec_store_nta_partial_lo(B3rhs[index],B3rhsL,elt_count);
      vec_store_nta_partial_lo(beta1rhs[index],beta1rhsL,elt_count);
      vec_store_nta_partial_lo(beta2rhs[index],beta2rhsL,elt_count);
      vec_store_nta_partial_lo(beta3rhs[index],beta3rhsL,elt_count);
      vec_store_nta_partial_lo(gt11rhs[index],gt11rhsL,elt_count);
      vec_store_nta_partial_lo(gt12rhs[index],gt12rhsL,elt_count);
      vec_store_nta_partial_lo(gt13rhs[index],gt13rhsL,elt_count);
      vec_store_nta_partial_lo(gt22rhs[index],gt22rhsL,elt_count);
      vec_store_nta_partial_lo(gt23rhs[index],gt23rhsL,elt_count);
      vec_store_nta_partial_lo(gt33rhs[index],gt33rhsL,elt_count);
      vec_store_nta_partial_lo(phirhs[index],phirhsL,elt_count);
      vec_store_nta_partial_lo(trKrhs[index],trKrhsL,elt_count);
      vec_store_nta_partial_lo(Xt1rhs[index],Xt1rhsL,elt_count);
      vec_store_nta_partial_lo(Xt2rhs[index],Xt2rhsL,elt_count);
      vec_store_nta_partial_lo(Xt3rhs[index],Xt3rhsL,elt_count);
      break;
    }
    vec_store_nta(alpharhs[index],alpharhsL);
    vec_store_nta(Arhs[index],ArhsL);
    vec_store_nta(At11rhs[index],At11rhsL);
    vec_store_nta(At12rhs[index],At12rhsL);
    vec_store_nta(At13rhs[index],At13rhsL);
    vec_store_nta(At22rhs[index],At22rhsL);
    vec_store_nta(At23rhs[index],At23rhsL);
    vec_store_nta(At33rhs[index],At33rhsL);
    vec_store_nta(B1rhs[index],B1rhsL);
    vec_store_nta(B2rhs[index],B2rhsL);
    vec_store_nta(B3rhs[index],B3rhsL);
    vec_store_nta(beta1rhs[index],beta1rhsL);
    vec_store_nta(beta2rhs[index],beta2rhsL);
    vec_store_nta(beta3rhs[index],beta3rhsL);
    vec_store_nta(gt11rhs[index],gt11rhsL);
    vec_store_nta(gt12rhs[index],gt12rhsL);
    vec_store_nta(gt13rhs[index],gt13rhsL);
    vec_store_nta(gt22rhs[index],gt22rhsL);
    vec_store_nta(gt23rhs[index],gt23rhsL);
    vec_store_nta(gt33rhs[index],gt33rhsL);
    vec_store_nta(phirhs[index],phirhsL);
    vec_store_nta(trKrhs[index],trKrhsL);
    vec_store_nta(Xt1rhs[index],Xt1rhsL);
    vec_store_nta(Xt2rhs[index],Xt2rhsL);
    vec_store_nta(Xt3rhs[index],Xt3rhsL);
  }
  LC_ENDLOOP3VEC (ML_BSSN_UPW_Dissipation);
}

extern "C" void ML_BSSN_UPW_Dissipation(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_UPW_Dissipation_Body");
  }
  
  if (cctk_iteration % ML_BSSN_UPW_Dissipation_calc_every != ML_BSSN_UPW_Dissipation_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_BSSN_UPW::ML_curv","ML_BSSN_UPW::ML_curvrhs","ML_BSSN_UPW::ML_dtlapse","ML_BSSN_UPW::ML_dtlapserhs","ML_BSSN_UPW::ML_dtshift","ML_BSSN_UPW::ML_dtshiftrhs","ML_BSSN_UPW::ML_Gamma","ML_BSSN_UPW::ML_Gammarhs","ML_BSSN_UPW::ML_lapse","ML_BSSN_UPW::ML_lapserhs","ML_BSSN_UPW::ML_log_confac","ML_BSSN_UPW::ML_log_confacrhs","ML_BSSN_UPW::ML_metric","ML_BSSN_UPW::ML_metricrhs","ML_BSSN_UPW::ML_shift","ML_BSSN_UPW::ML_shiftrhs","ML_BSSN_UPW::ML_trace_curv","ML_BSSN_UPW::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_UPW_Dissipation", 18, groups);
  
  switch(fdOrder)
  {
    case 2:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_Dissipation", 2, 2, 2);
      break;
    
    case 4:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_Dissipation", 3, 3, 3);
      break;
    
    case 6:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_Dissipation", 4, 4, 4);
      break;
    
    case 8:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_Dissipation", 5, 5, 5);
      break;
  }
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_UPW_Dissipation_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_UPW_Dissipation_Body");
  }
}
