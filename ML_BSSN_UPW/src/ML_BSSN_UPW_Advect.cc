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
#include "vectors.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define QAD(x) (SQR(SQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))

extern "C" void ML_BSSN_UPW_Advect_SelectBCs(CCTK_ARGUMENTS)
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

static void ML_BSSN_UPW_Advect_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  LC_LOOP3VEC (ML_BSSN_UPW_Advect,
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
    CCTK_REAL_VEC PDupwindNth1A;
    CCTK_REAL_VEC PDupwindNth2A;
    CCTK_REAL_VEC PDupwindNth3A;
    CCTK_REAL_VEC PDupwindNth1alpha;
    CCTK_REAL_VEC PDupwindNth2alpha;
    CCTK_REAL_VEC PDupwindNth3alpha;
    CCTK_REAL_VEC PDupwindNth1At11;
    CCTK_REAL_VEC PDupwindNth2At11;
    CCTK_REAL_VEC PDupwindNth3At11;
    CCTK_REAL_VEC PDupwindNth1At12;
    CCTK_REAL_VEC PDupwindNth2At12;
    CCTK_REAL_VEC PDupwindNth3At12;
    CCTK_REAL_VEC PDupwindNth1At13;
    CCTK_REAL_VEC PDupwindNth2At13;
    CCTK_REAL_VEC PDupwindNth3At13;
    CCTK_REAL_VEC PDupwindNth1At22;
    CCTK_REAL_VEC PDupwindNth2At22;
    CCTK_REAL_VEC PDupwindNth3At22;
    CCTK_REAL_VEC PDupwindNth1At23;
    CCTK_REAL_VEC PDupwindNth2At23;
    CCTK_REAL_VEC PDupwindNth3At23;
    CCTK_REAL_VEC PDupwindNth1At33;
    CCTK_REAL_VEC PDupwindNth2At33;
    CCTK_REAL_VEC PDupwindNth3At33;
    CCTK_REAL_VEC PDupwindNth1B1;
    CCTK_REAL_VEC PDupwindNth2B1;
    CCTK_REAL_VEC PDupwindNth3B1;
    CCTK_REAL_VEC PDupwindNth1B2;
    CCTK_REAL_VEC PDupwindNth2B2;
    CCTK_REAL_VEC PDupwindNth3B2;
    CCTK_REAL_VEC PDupwindNth1B3;
    CCTK_REAL_VEC PDupwindNth2B3;
    CCTK_REAL_VEC PDupwindNth3B3;
    CCTK_REAL_VEC PDupwindNth1beta1;
    CCTK_REAL_VEC PDupwindNth2beta1;
    CCTK_REAL_VEC PDupwindNth3beta1;
    CCTK_REAL_VEC PDupwindNth1beta2;
    CCTK_REAL_VEC PDupwindNth2beta2;
    CCTK_REAL_VEC PDupwindNth3beta2;
    CCTK_REAL_VEC PDupwindNth1beta3;
    CCTK_REAL_VEC PDupwindNth2beta3;
    CCTK_REAL_VEC PDupwindNth3beta3;
    CCTK_REAL_VEC PDupwindNth1gt11;
    CCTK_REAL_VEC PDupwindNth2gt11;
    CCTK_REAL_VEC PDupwindNth3gt11;
    CCTK_REAL_VEC PDupwindNth1gt12;
    CCTK_REAL_VEC PDupwindNth2gt12;
    CCTK_REAL_VEC PDupwindNth3gt12;
    CCTK_REAL_VEC PDupwindNth1gt13;
    CCTK_REAL_VEC PDupwindNth2gt13;
    CCTK_REAL_VEC PDupwindNth3gt13;
    CCTK_REAL_VEC PDupwindNth1gt22;
    CCTK_REAL_VEC PDupwindNth2gt22;
    CCTK_REAL_VEC PDupwindNth3gt22;
    CCTK_REAL_VEC PDupwindNth1gt23;
    CCTK_REAL_VEC PDupwindNth2gt23;
    CCTK_REAL_VEC PDupwindNth3gt23;
    CCTK_REAL_VEC PDupwindNth1gt33;
    CCTK_REAL_VEC PDupwindNth2gt33;
    CCTK_REAL_VEC PDupwindNth3gt33;
    CCTK_REAL_VEC PDupwindNth1phi;
    CCTK_REAL_VEC PDupwindNth2phi;
    CCTK_REAL_VEC PDupwindNth3phi;
    CCTK_REAL_VEC PDupwindNth1trK;
    CCTK_REAL_VEC PDupwindNth2trK;
    CCTK_REAL_VEC PDupwindNth3trK;
    CCTK_REAL_VEC PDupwindNth1Xt1;
    CCTK_REAL_VEC PDupwindNth2Xt1;
    CCTK_REAL_VEC PDupwindNth3Xt1;
    CCTK_REAL_VEC PDupwindNth1Xt2;
    CCTK_REAL_VEC PDupwindNth2Xt2;
    CCTK_REAL_VEC PDupwindNth3Xt2;
    CCTK_REAL_VEC PDupwindNth1Xt3;
    CCTK_REAL_VEC PDupwindNth2Xt3;
    CCTK_REAL_VEC PDupwindNth3Xt3;
    
    switch(fdOrder)
    {
      case 2:
        PDupwindNth1A = PDupwindNthfdOrder21(&A[index]);
        PDupwindNth2A = PDupwindNthfdOrder22(&A[index]);
        PDupwindNth3A = PDupwindNthfdOrder23(&A[index]);
        PDupwindNth1alpha = PDupwindNthfdOrder21(&alpha[index]);
        PDupwindNth2alpha = PDupwindNthfdOrder22(&alpha[index]);
        PDupwindNth3alpha = PDupwindNthfdOrder23(&alpha[index]);
        PDupwindNth1At11 = PDupwindNthfdOrder21(&At11[index]);
        PDupwindNth2At11 = PDupwindNthfdOrder22(&At11[index]);
        PDupwindNth3At11 = PDupwindNthfdOrder23(&At11[index]);
        PDupwindNth1At12 = PDupwindNthfdOrder21(&At12[index]);
        PDupwindNth2At12 = PDupwindNthfdOrder22(&At12[index]);
        PDupwindNth3At12 = PDupwindNthfdOrder23(&At12[index]);
        PDupwindNth1At13 = PDupwindNthfdOrder21(&At13[index]);
        PDupwindNth2At13 = PDupwindNthfdOrder22(&At13[index]);
        PDupwindNth3At13 = PDupwindNthfdOrder23(&At13[index]);
        PDupwindNth1At22 = PDupwindNthfdOrder21(&At22[index]);
        PDupwindNth2At22 = PDupwindNthfdOrder22(&At22[index]);
        PDupwindNth3At22 = PDupwindNthfdOrder23(&At22[index]);
        PDupwindNth1At23 = PDupwindNthfdOrder21(&At23[index]);
        PDupwindNth2At23 = PDupwindNthfdOrder22(&At23[index]);
        PDupwindNth3At23 = PDupwindNthfdOrder23(&At23[index]);
        PDupwindNth1At33 = PDupwindNthfdOrder21(&At33[index]);
        PDupwindNth2At33 = PDupwindNthfdOrder22(&At33[index]);
        PDupwindNth3At33 = PDupwindNthfdOrder23(&At33[index]);
        PDupwindNth1B1 = PDupwindNthfdOrder21(&B1[index]);
        PDupwindNth2B1 = PDupwindNthfdOrder22(&B1[index]);
        PDupwindNth3B1 = PDupwindNthfdOrder23(&B1[index]);
        PDupwindNth1B2 = PDupwindNthfdOrder21(&B2[index]);
        PDupwindNth2B2 = PDupwindNthfdOrder22(&B2[index]);
        PDupwindNth3B2 = PDupwindNthfdOrder23(&B2[index]);
        PDupwindNth1B3 = PDupwindNthfdOrder21(&B3[index]);
        PDupwindNth2B3 = PDupwindNthfdOrder22(&B3[index]);
        PDupwindNth3B3 = PDupwindNthfdOrder23(&B3[index]);
        PDupwindNth1beta1 = PDupwindNthfdOrder21(&beta1[index]);
        PDupwindNth2beta1 = PDupwindNthfdOrder22(&beta1[index]);
        PDupwindNth3beta1 = PDupwindNthfdOrder23(&beta1[index]);
        PDupwindNth1beta2 = PDupwindNthfdOrder21(&beta2[index]);
        PDupwindNth2beta2 = PDupwindNthfdOrder22(&beta2[index]);
        PDupwindNth3beta2 = PDupwindNthfdOrder23(&beta2[index]);
        PDupwindNth1beta3 = PDupwindNthfdOrder21(&beta3[index]);
        PDupwindNth2beta3 = PDupwindNthfdOrder22(&beta3[index]);
        PDupwindNth3beta3 = PDupwindNthfdOrder23(&beta3[index]);
        PDupwindNth1gt11 = PDupwindNthfdOrder21(&gt11[index]);
        PDupwindNth2gt11 = PDupwindNthfdOrder22(&gt11[index]);
        PDupwindNth3gt11 = PDupwindNthfdOrder23(&gt11[index]);
        PDupwindNth1gt12 = PDupwindNthfdOrder21(&gt12[index]);
        PDupwindNth2gt12 = PDupwindNthfdOrder22(&gt12[index]);
        PDupwindNth3gt12 = PDupwindNthfdOrder23(&gt12[index]);
        PDupwindNth1gt13 = PDupwindNthfdOrder21(&gt13[index]);
        PDupwindNth2gt13 = PDupwindNthfdOrder22(&gt13[index]);
        PDupwindNth3gt13 = PDupwindNthfdOrder23(&gt13[index]);
        PDupwindNth1gt22 = PDupwindNthfdOrder21(&gt22[index]);
        PDupwindNth2gt22 = PDupwindNthfdOrder22(&gt22[index]);
        PDupwindNth3gt22 = PDupwindNthfdOrder23(&gt22[index]);
        PDupwindNth1gt23 = PDupwindNthfdOrder21(&gt23[index]);
        PDupwindNth2gt23 = PDupwindNthfdOrder22(&gt23[index]);
        PDupwindNth3gt23 = PDupwindNthfdOrder23(&gt23[index]);
        PDupwindNth1gt33 = PDupwindNthfdOrder21(&gt33[index]);
        PDupwindNth2gt33 = PDupwindNthfdOrder22(&gt33[index]);
        PDupwindNth3gt33 = PDupwindNthfdOrder23(&gt33[index]);
        PDupwindNth1phi = PDupwindNthfdOrder21(&phi[index]);
        PDupwindNth2phi = PDupwindNthfdOrder22(&phi[index]);
        PDupwindNth3phi = PDupwindNthfdOrder23(&phi[index]);
        PDupwindNth1trK = PDupwindNthfdOrder21(&trK[index]);
        PDupwindNth2trK = PDupwindNthfdOrder22(&trK[index]);
        PDupwindNth3trK = PDupwindNthfdOrder23(&trK[index]);
        PDupwindNth1Xt1 = PDupwindNthfdOrder21(&Xt1[index]);
        PDupwindNth2Xt1 = PDupwindNthfdOrder22(&Xt1[index]);
        PDupwindNth3Xt1 = PDupwindNthfdOrder23(&Xt1[index]);
        PDupwindNth1Xt2 = PDupwindNthfdOrder21(&Xt2[index]);
        PDupwindNth2Xt2 = PDupwindNthfdOrder22(&Xt2[index]);
        PDupwindNth3Xt2 = PDupwindNthfdOrder23(&Xt2[index]);
        PDupwindNth1Xt3 = PDupwindNthfdOrder21(&Xt3[index]);
        PDupwindNth2Xt3 = PDupwindNthfdOrder22(&Xt3[index]);
        PDupwindNth3Xt3 = PDupwindNthfdOrder23(&Xt3[index]);
        break;
      
      case 4:
        PDupwindNth1A = PDupwindNthfdOrder41(&A[index]);
        PDupwindNth2A = PDupwindNthfdOrder42(&A[index]);
        PDupwindNth3A = PDupwindNthfdOrder43(&A[index]);
        PDupwindNth1alpha = PDupwindNthfdOrder41(&alpha[index]);
        PDupwindNth2alpha = PDupwindNthfdOrder42(&alpha[index]);
        PDupwindNth3alpha = PDupwindNthfdOrder43(&alpha[index]);
        PDupwindNth1At11 = PDupwindNthfdOrder41(&At11[index]);
        PDupwindNth2At11 = PDupwindNthfdOrder42(&At11[index]);
        PDupwindNth3At11 = PDupwindNthfdOrder43(&At11[index]);
        PDupwindNth1At12 = PDupwindNthfdOrder41(&At12[index]);
        PDupwindNth2At12 = PDupwindNthfdOrder42(&At12[index]);
        PDupwindNth3At12 = PDupwindNthfdOrder43(&At12[index]);
        PDupwindNth1At13 = PDupwindNthfdOrder41(&At13[index]);
        PDupwindNth2At13 = PDupwindNthfdOrder42(&At13[index]);
        PDupwindNth3At13 = PDupwindNthfdOrder43(&At13[index]);
        PDupwindNth1At22 = PDupwindNthfdOrder41(&At22[index]);
        PDupwindNth2At22 = PDupwindNthfdOrder42(&At22[index]);
        PDupwindNth3At22 = PDupwindNthfdOrder43(&At22[index]);
        PDupwindNth1At23 = PDupwindNthfdOrder41(&At23[index]);
        PDupwindNth2At23 = PDupwindNthfdOrder42(&At23[index]);
        PDupwindNth3At23 = PDupwindNthfdOrder43(&At23[index]);
        PDupwindNth1At33 = PDupwindNthfdOrder41(&At33[index]);
        PDupwindNth2At33 = PDupwindNthfdOrder42(&At33[index]);
        PDupwindNth3At33 = PDupwindNthfdOrder43(&At33[index]);
        PDupwindNth1B1 = PDupwindNthfdOrder41(&B1[index]);
        PDupwindNth2B1 = PDupwindNthfdOrder42(&B1[index]);
        PDupwindNth3B1 = PDupwindNthfdOrder43(&B1[index]);
        PDupwindNth1B2 = PDupwindNthfdOrder41(&B2[index]);
        PDupwindNth2B2 = PDupwindNthfdOrder42(&B2[index]);
        PDupwindNth3B2 = PDupwindNthfdOrder43(&B2[index]);
        PDupwindNth1B3 = PDupwindNthfdOrder41(&B3[index]);
        PDupwindNth2B3 = PDupwindNthfdOrder42(&B3[index]);
        PDupwindNth3B3 = PDupwindNthfdOrder43(&B3[index]);
        PDupwindNth1beta1 = PDupwindNthfdOrder41(&beta1[index]);
        PDupwindNth2beta1 = PDupwindNthfdOrder42(&beta1[index]);
        PDupwindNth3beta1 = PDupwindNthfdOrder43(&beta1[index]);
        PDupwindNth1beta2 = PDupwindNthfdOrder41(&beta2[index]);
        PDupwindNth2beta2 = PDupwindNthfdOrder42(&beta2[index]);
        PDupwindNth3beta2 = PDupwindNthfdOrder43(&beta2[index]);
        PDupwindNth1beta3 = PDupwindNthfdOrder41(&beta3[index]);
        PDupwindNth2beta3 = PDupwindNthfdOrder42(&beta3[index]);
        PDupwindNth3beta3 = PDupwindNthfdOrder43(&beta3[index]);
        PDupwindNth1gt11 = PDupwindNthfdOrder41(&gt11[index]);
        PDupwindNth2gt11 = PDupwindNthfdOrder42(&gt11[index]);
        PDupwindNth3gt11 = PDupwindNthfdOrder43(&gt11[index]);
        PDupwindNth1gt12 = PDupwindNthfdOrder41(&gt12[index]);
        PDupwindNth2gt12 = PDupwindNthfdOrder42(&gt12[index]);
        PDupwindNth3gt12 = PDupwindNthfdOrder43(&gt12[index]);
        PDupwindNth1gt13 = PDupwindNthfdOrder41(&gt13[index]);
        PDupwindNth2gt13 = PDupwindNthfdOrder42(&gt13[index]);
        PDupwindNth3gt13 = PDupwindNthfdOrder43(&gt13[index]);
        PDupwindNth1gt22 = PDupwindNthfdOrder41(&gt22[index]);
        PDupwindNth2gt22 = PDupwindNthfdOrder42(&gt22[index]);
        PDupwindNth3gt22 = PDupwindNthfdOrder43(&gt22[index]);
        PDupwindNth1gt23 = PDupwindNthfdOrder41(&gt23[index]);
        PDupwindNth2gt23 = PDupwindNthfdOrder42(&gt23[index]);
        PDupwindNth3gt23 = PDupwindNthfdOrder43(&gt23[index]);
        PDupwindNth1gt33 = PDupwindNthfdOrder41(&gt33[index]);
        PDupwindNth2gt33 = PDupwindNthfdOrder42(&gt33[index]);
        PDupwindNth3gt33 = PDupwindNthfdOrder43(&gt33[index]);
        PDupwindNth1phi = PDupwindNthfdOrder41(&phi[index]);
        PDupwindNth2phi = PDupwindNthfdOrder42(&phi[index]);
        PDupwindNth3phi = PDupwindNthfdOrder43(&phi[index]);
        PDupwindNth1trK = PDupwindNthfdOrder41(&trK[index]);
        PDupwindNth2trK = PDupwindNthfdOrder42(&trK[index]);
        PDupwindNth3trK = PDupwindNthfdOrder43(&trK[index]);
        PDupwindNth1Xt1 = PDupwindNthfdOrder41(&Xt1[index]);
        PDupwindNth2Xt1 = PDupwindNthfdOrder42(&Xt1[index]);
        PDupwindNth3Xt1 = PDupwindNthfdOrder43(&Xt1[index]);
        PDupwindNth1Xt2 = PDupwindNthfdOrder41(&Xt2[index]);
        PDupwindNth2Xt2 = PDupwindNthfdOrder42(&Xt2[index]);
        PDupwindNth3Xt2 = PDupwindNthfdOrder43(&Xt2[index]);
        PDupwindNth1Xt3 = PDupwindNthfdOrder41(&Xt3[index]);
        PDupwindNth2Xt3 = PDupwindNthfdOrder42(&Xt3[index]);
        PDupwindNth3Xt3 = PDupwindNthfdOrder43(&Xt3[index]);
        break;
      
      case 6:
        PDupwindNth1A = PDupwindNthfdOrder61(&A[index]);
        PDupwindNth2A = PDupwindNthfdOrder62(&A[index]);
        PDupwindNth3A = PDupwindNthfdOrder63(&A[index]);
        PDupwindNth1alpha = PDupwindNthfdOrder61(&alpha[index]);
        PDupwindNth2alpha = PDupwindNthfdOrder62(&alpha[index]);
        PDupwindNth3alpha = PDupwindNthfdOrder63(&alpha[index]);
        PDupwindNth1At11 = PDupwindNthfdOrder61(&At11[index]);
        PDupwindNth2At11 = PDupwindNthfdOrder62(&At11[index]);
        PDupwindNth3At11 = PDupwindNthfdOrder63(&At11[index]);
        PDupwindNth1At12 = PDupwindNthfdOrder61(&At12[index]);
        PDupwindNth2At12 = PDupwindNthfdOrder62(&At12[index]);
        PDupwindNth3At12 = PDupwindNthfdOrder63(&At12[index]);
        PDupwindNth1At13 = PDupwindNthfdOrder61(&At13[index]);
        PDupwindNth2At13 = PDupwindNthfdOrder62(&At13[index]);
        PDupwindNth3At13 = PDupwindNthfdOrder63(&At13[index]);
        PDupwindNth1At22 = PDupwindNthfdOrder61(&At22[index]);
        PDupwindNth2At22 = PDupwindNthfdOrder62(&At22[index]);
        PDupwindNth3At22 = PDupwindNthfdOrder63(&At22[index]);
        PDupwindNth1At23 = PDupwindNthfdOrder61(&At23[index]);
        PDupwindNth2At23 = PDupwindNthfdOrder62(&At23[index]);
        PDupwindNth3At23 = PDupwindNthfdOrder63(&At23[index]);
        PDupwindNth1At33 = PDupwindNthfdOrder61(&At33[index]);
        PDupwindNth2At33 = PDupwindNthfdOrder62(&At33[index]);
        PDupwindNth3At33 = PDupwindNthfdOrder63(&At33[index]);
        PDupwindNth1B1 = PDupwindNthfdOrder61(&B1[index]);
        PDupwindNth2B1 = PDupwindNthfdOrder62(&B1[index]);
        PDupwindNth3B1 = PDupwindNthfdOrder63(&B1[index]);
        PDupwindNth1B2 = PDupwindNthfdOrder61(&B2[index]);
        PDupwindNth2B2 = PDupwindNthfdOrder62(&B2[index]);
        PDupwindNth3B2 = PDupwindNthfdOrder63(&B2[index]);
        PDupwindNth1B3 = PDupwindNthfdOrder61(&B3[index]);
        PDupwindNth2B3 = PDupwindNthfdOrder62(&B3[index]);
        PDupwindNth3B3 = PDupwindNthfdOrder63(&B3[index]);
        PDupwindNth1beta1 = PDupwindNthfdOrder61(&beta1[index]);
        PDupwindNth2beta1 = PDupwindNthfdOrder62(&beta1[index]);
        PDupwindNth3beta1 = PDupwindNthfdOrder63(&beta1[index]);
        PDupwindNth1beta2 = PDupwindNthfdOrder61(&beta2[index]);
        PDupwindNth2beta2 = PDupwindNthfdOrder62(&beta2[index]);
        PDupwindNth3beta2 = PDupwindNthfdOrder63(&beta2[index]);
        PDupwindNth1beta3 = PDupwindNthfdOrder61(&beta3[index]);
        PDupwindNth2beta3 = PDupwindNthfdOrder62(&beta3[index]);
        PDupwindNth3beta3 = PDupwindNthfdOrder63(&beta3[index]);
        PDupwindNth1gt11 = PDupwindNthfdOrder61(&gt11[index]);
        PDupwindNth2gt11 = PDupwindNthfdOrder62(&gt11[index]);
        PDupwindNth3gt11 = PDupwindNthfdOrder63(&gt11[index]);
        PDupwindNth1gt12 = PDupwindNthfdOrder61(&gt12[index]);
        PDupwindNth2gt12 = PDupwindNthfdOrder62(&gt12[index]);
        PDupwindNth3gt12 = PDupwindNthfdOrder63(&gt12[index]);
        PDupwindNth1gt13 = PDupwindNthfdOrder61(&gt13[index]);
        PDupwindNth2gt13 = PDupwindNthfdOrder62(&gt13[index]);
        PDupwindNth3gt13 = PDupwindNthfdOrder63(&gt13[index]);
        PDupwindNth1gt22 = PDupwindNthfdOrder61(&gt22[index]);
        PDupwindNth2gt22 = PDupwindNthfdOrder62(&gt22[index]);
        PDupwindNth3gt22 = PDupwindNthfdOrder63(&gt22[index]);
        PDupwindNth1gt23 = PDupwindNthfdOrder61(&gt23[index]);
        PDupwindNth2gt23 = PDupwindNthfdOrder62(&gt23[index]);
        PDupwindNth3gt23 = PDupwindNthfdOrder63(&gt23[index]);
        PDupwindNth1gt33 = PDupwindNthfdOrder61(&gt33[index]);
        PDupwindNth2gt33 = PDupwindNthfdOrder62(&gt33[index]);
        PDupwindNth3gt33 = PDupwindNthfdOrder63(&gt33[index]);
        PDupwindNth1phi = PDupwindNthfdOrder61(&phi[index]);
        PDupwindNth2phi = PDupwindNthfdOrder62(&phi[index]);
        PDupwindNth3phi = PDupwindNthfdOrder63(&phi[index]);
        PDupwindNth1trK = PDupwindNthfdOrder61(&trK[index]);
        PDupwindNth2trK = PDupwindNthfdOrder62(&trK[index]);
        PDupwindNth3trK = PDupwindNthfdOrder63(&trK[index]);
        PDupwindNth1Xt1 = PDupwindNthfdOrder61(&Xt1[index]);
        PDupwindNth2Xt1 = PDupwindNthfdOrder62(&Xt1[index]);
        PDupwindNth3Xt1 = PDupwindNthfdOrder63(&Xt1[index]);
        PDupwindNth1Xt2 = PDupwindNthfdOrder61(&Xt2[index]);
        PDupwindNth2Xt2 = PDupwindNthfdOrder62(&Xt2[index]);
        PDupwindNth3Xt2 = PDupwindNthfdOrder63(&Xt2[index]);
        PDupwindNth1Xt3 = PDupwindNthfdOrder61(&Xt3[index]);
        PDupwindNth2Xt3 = PDupwindNthfdOrder62(&Xt3[index]);
        PDupwindNth3Xt3 = PDupwindNthfdOrder63(&Xt3[index]);
        break;
      
      case 8:
        PDupwindNth1A = PDupwindNthfdOrder81(&A[index]);
        PDupwindNth2A = PDupwindNthfdOrder82(&A[index]);
        PDupwindNth3A = PDupwindNthfdOrder83(&A[index]);
        PDupwindNth1alpha = PDupwindNthfdOrder81(&alpha[index]);
        PDupwindNth2alpha = PDupwindNthfdOrder82(&alpha[index]);
        PDupwindNth3alpha = PDupwindNthfdOrder83(&alpha[index]);
        PDupwindNth1At11 = PDupwindNthfdOrder81(&At11[index]);
        PDupwindNth2At11 = PDupwindNthfdOrder82(&At11[index]);
        PDupwindNth3At11 = PDupwindNthfdOrder83(&At11[index]);
        PDupwindNth1At12 = PDupwindNthfdOrder81(&At12[index]);
        PDupwindNth2At12 = PDupwindNthfdOrder82(&At12[index]);
        PDupwindNth3At12 = PDupwindNthfdOrder83(&At12[index]);
        PDupwindNth1At13 = PDupwindNthfdOrder81(&At13[index]);
        PDupwindNth2At13 = PDupwindNthfdOrder82(&At13[index]);
        PDupwindNth3At13 = PDupwindNthfdOrder83(&At13[index]);
        PDupwindNth1At22 = PDupwindNthfdOrder81(&At22[index]);
        PDupwindNth2At22 = PDupwindNthfdOrder82(&At22[index]);
        PDupwindNth3At22 = PDupwindNthfdOrder83(&At22[index]);
        PDupwindNth1At23 = PDupwindNthfdOrder81(&At23[index]);
        PDupwindNth2At23 = PDupwindNthfdOrder82(&At23[index]);
        PDupwindNth3At23 = PDupwindNthfdOrder83(&At23[index]);
        PDupwindNth1At33 = PDupwindNthfdOrder81(&At33[index]);
        PDupwindNth2At33 = PDupwindNthfdOrder82(&At33[index]);
        PDupwindNth3At33 = PDupwindNthfdOrder83(&At33[index]);
        PDupwindNth1B1 = PDupwindNthfdOrder81(&B1[index]);
        PDupwindNth2B1 = PDupwindNthfdOrder82(&B1[index]);
        PDupwindNth3B1 = PDupwindNthfdOrder83(&B1[index]);
        PDupwindNth1B2 = PDupwindNthfdOrder81(&B2[index]);
        PDupwindNth2B2 = PDupwindNthfdOrder82(&B2[index]);
        PDupwindNth3B2 = PDupwindNthfdOrder83(&B2[index]);
        PDupwindNth1B3 = PDupwindNthfdOrder81(&B3[index]);
        PDupwindNth2B3 = PDupwindNthfdOrder82(&B3[index]);
        PDupwindNth3B3 = PDupwindNthfdOrder83(&B3[index]);
        PDupwindNth1beta1 = PDupwindNthfdOrder81(&beta1[index]);
        PDupwindNth2beta1 = PDupwindNthfdOrder82(&beta1[index]);
        PDupwindNth3beta1 = PDupwindNthfdOrder83(&beta1[index]);
        PDupwindNth1beta2 = PDupwindNthfdOrder81(&beta2[index]);
        PDupwindNth2beta2 = PDupwindNthfdOrder82(&beta2[index]);
        PDupwindNth3beta2 = PDupwindNthfdOrder83(&beta2[index]);
        PDupwindNth1beta3 = PDupwindNthfdOrder81(&beta3[index]);
        PDupwindNth2beta3 = PDupwindNthfdOrder82(&beta3[index]);
        PDupwindNth3beta3 = PDupwindNthfdOrder83(&beta3[index]);
        PDupwindNth1gt11 = PDupwindNthfdOrder81(&gt11[index]);
        PDupwindNth2gt11 = PDupwindNthfdOrder82(&gt11[index]);
        PDupwindNth3gt11 = PDupwindNthfdOrder83(&gt11[index]);
        PDupwindNth1gt12 = PDupwindNthfdOrder81(&gt12[index]);
        PDupwindNth2gt12 = PDupwindNthfdOrder82(&gt12[index]);
        PDupwindNth3gt12 = PDupwindNthfdOrder83(&gt12[index]);
        PDupwindNth1gt13 = PDupwindNthfdOrder81(&gt13[index]);
        PDupwindNth2gt13 = PDupwindNthfdOrder82(&gt13[index]);
        PDupwindNth3gt13 = PDupwindNthfdOrder83(&gt13[index]);
        PDupwindNth1gt22 = PDupwindNthfdOrder81(&gt22[index]);
        PDupwindNth2gt22 = PDupwindNthfdOrder82(&gt22[index]);
        PDupwindNth3gt22 = PDupwindNthfdOrder83(&gt22[index]);
        PDupwindNth1gt23 = PDupwindNthfdOrder81(&gt23[index]);
        PDupwindNth2gt23 = PDupwindNthfdOrder82(&gt23[index]);
        PDupwindNth3gt23 = PDupwindNthfdOrder83(&gt23[index]);
        PDupwindNth1gt33 = PDupwindNthfdOrder81(&gt33[index]);
        PDupwindNth2gt33 = PDupwindNthfdOrder82(&gt33[index]);
        PDupwindNth3gt33 = PDupwindNthfdOrder83(&gt33[index]);
        PDupwindNth1phi = PDupwindNthfdOrder81(&phi[index]);
        PDupwindNth2phi = PDupwindNthfdOrder82(&phi[index]);
        PDupwindNth3phi = PDupwindNthfdOrder83(&phi[index]);
        PDupwindNth1trK = PDupwindNthfdOrder81(&trK[index]);
        PDupwindNth2trK = PDupwindNthfdOrder82(&trK[index]);
        PDupwindNth3trK = PDupwindNthfdOrder83(&trK[index]);
        PDupwindNth1Xt1 = PDupwindNthfdOrder81(&Xt1[index]);
        PDupwindNth2Xt1 = PDupwindNthfdOrder82(&Xt1[index]);
        PDupwindNth3Xt1 = PDupwindNthfdOrder83(&Xt1[index]);
        PDupwindNth1Xt2 = PDupwindNthfdOrder81(&Xt2[index]);
        PDupwindNth2Xt2 = PDupwindNthfdOrder82(&Xt2[index]);
        PDupwindNth3Xt2 = PDupwindNthfdOrder83(&Xt2[index]);
        PDupwindNth1Xt3 = PDupwindNthfdOrder81(&Xt3[index]);
        PDupwindNth2Xt3 = PDupwindNthfdOrder82(&Xt3[index]);
        PDupwindNth3Xt3 = PDupwindNthfdOrder83(&Xt3[index]);
        break;
    }
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    CCTK_REAL_VEC JacPDupwindNth1A;
    CCTK_REAL_VEC JacPDupwindNth1alpha;
    CCTK_REAL_VEC JacPDupwindNth1At11;
    CCTK_REAL_VEC JacPDupwindNth1At12;
    CCTK_REAL_VEC JacPDupwindNth1At13;
    CCTK_REAL_VEC JacPDupwindNth1At22;
    CCTK_REAL_VEC JacPDupwindNth1At23;
    CCTK_REAL_VEC JacPDupwindNth1At33;
    CCTK_REAL_VEC JacPDupwindNth1B1;
    CCTK_REAL_VEC JacPDupwindNth1B2;
    CCTK_REAL_VEC JacPDupwindNth1B3;
    CCTK_REAL_VEC JacPDupwindNth1beta1;
    CCTK_REAL_VEC JacPDupwindNth1beta2;
    CCTK_REAL_VEC JacPDupwindNth1beta3;
    CCTK_REAL_VEC JacPDupwindNth1gt11;
    CCTK_REAL_VEC JacPDupwindNth1gt12;
    CCTK_REAL_VEC JacPDupwindNth1gt13;
    CCTK_REAL_VEC JacPDupwindNth1gt22;
    CCTK_REAL_VEC JacPDupwindNth1gt23;
    CCTK_REAL_VEC JacPDupwindNth1gt33;
    CCTK_REAL_VEC JacPDupwindNth1phi;
    CCTK_REAL_VEC JacPDupwindNth1trK;
    CCTK_REAL_VEC JacPDupwindNth1Xt1;
    CCTK_REAL_VEC JacPDupwindNth1Xt2;
    CCTK_REAL_VEC JacPDupwindNth1Xt3;
    CCTK_REAL_VEC JacPDupwindNth2A;
    CCTK_REAL_VEC JacPDupwindNth2alpha;
    CCTK_REAL_VEC JacPDupwindNth2At11;
    CCTK_REAL_VEC JacPDupwindNth2At12;
    CCTK_REAL_VEC JacPDupwindNth2At13;
    CCTK_REAL_VEC JacPDupwindNth2At22;
    CCTK_REAL_VEC JacPDupwindNth2At23;
    CCTK_REAL_VEC JacPDupwindNth2At33;
    CCTK_REAL_VEC JacPDupwindNth2B1;
    CCTK_REAL_VEC JacPDupwindNth2B2;
    CCTK_REAL_VEC JacPDupwindNth2B3;
    CCTK_REAL_VEC JacPDupwindNth2beta1;
    CCTK_REAL_VEC JacPDupwindNth2beta2;
    CCTK_REAL_VEC JacPDupwindNth2beta3;
    CCTK_REAL_VEC JacPDupwindNth2gt11;
    CCTK_REAL_VEC JacPDupwindNth2gt12;
    CCTK_REAL_VEC JacPDupwindNth2gt13;
    CCTK_REAL_VEC JacPDupwindNth2gt22;
    CCTK_REAL_VEC JacPDupwindNth2gt23;
    CCTK_REAL_VEC JacPDupwindNth2gt33;
    CCTK_REAL_VEC JacPDupwindNth2phi;
    CCTK_REAL_VEC JacPDupwindNth2trK;
    CCTK_REAL_VEC JacPDupwindNth2Xt1;
    CCTK_REAL_VEC JacPDupwindNth2Xt2;
    CCTK_REAL_VEC JacPDupwindNth2Xt3;
    CCTK_REAL_VEC JacPDupwindNth3A;
    CCTK_REAL_VEC JacPDupwindNth3alpha;
    CCTK_REAL_VEC JacPDupwindNth3At11;
    CCTK_REAL_VEC JacPDupwindNth3At12;
    CCTK_REAL_VEC JacPDupwindNth3At13;
    CCTK_REAL_VEC JacPDupwindNth3At22;
    CCTK_REAL_VEC JacPDupwindNth3At23;
    CCTK_REAL_VEC JacPDupwindNth3At33;
    CCTK_REAL_VEC JacPDupwindNth3B1;
    CCTK_REAL_VEC JacPDupwindNth3B2;
    CCTK_REAL_VEC JacPDupwindNth3B3;
    CCTK_REAL_VEC JacPDupwindNth3beta1;
    CCTK_REAL_VEC JacPDupwindNth3beta2;
    CCTK_REAL_VEC JacPDupwindNth3beta3;
    CCTK_REAL_VEC JacPDupwindNth3gt11;
    CCTK_REAL_VEC JacPDupwindNth3gt12;
    CCTK_REAL_VEC JacPDupwindNth3gt13;
    CCTK_REAL_VEC JacPDupwindNth3gt22;
    CCTK_REAL_VEC JacPDupwindNth3gt23;
    CCTK_REAL_VEC JacPDupwindNth3gt33;
    CCTK_REAL_VEC JacPDupwindNth3phi;
    CCTK_REAL_VEC JacPDupwindNth3trK;
    CCTK_REAL_VEC JacPDupwindNth3Xt1;
    CCTK_REAL_VEC JacPDupwindNth3Xt2;
    CCTK_REAL_VEC JacPDupwindNth3Xt3;
    
    if (use_jacobian)
    {
      JacPDupwindNth1A = 
        kmadd(J11L,PDupwindNth1A,kmadd(J21L,PDupwindNth2A,kmul(J31L,PDupwindNth3A)));
      
      JacPDupwindNth1alpha = 
        kmadd(J11L,PDupwindNth1alpha,kmadd(J21L,PDupwindNth2alpha,kmul(J31L,PDupwindNth3alpha)));
      
      JacPDupwindNth1At11 = 
        kmadd(J11L,PDupwindNth1At11,kmadd(J21L,PDupwindNth2At11,kmul(J31L,PDupwindNth3At11)));
      
      JacPDupwindNth1At12 = 
        kmadd(J11L,PDupwindNth1At12,kmadd(J21L,PDupwindNth2At12,kmul(J31L,PDupwindNth3At12)));
      
      JacPDupwindNth1At13 = 
        kmadd(J11L,PDupwindNth1At13,kmadd(J21L,PDupwindNth2At13,kmul(J31L,PDupwindNth3At13)));
      
      JacPDupwindNth1At22 = 
        kmadd(J11L,PDupwindNth1At22,kmadd(J21L,PDupwindNth2At22,kmul(J31L,PDupwindNth3At22)));
      
      JacPDupwindNth1At23 = 
        kmadd(J11L,PDupwindNth1At23,kmadd(J21L,PDupwindNth2At23,kmul(J31L,PDupwindNth3At23)));
      
      JacPDupwindNth1At33 = 
        kmadd(J11L,PDupwindNth1At33,kmadd(J21L,PDupwindNth2At33,kmul(J31L,PDupwindNth3At33)));
      
      JacPDupwindNth1B1 = 
        kmadd(J11L,PDupwindNth1B1,kmadd(J21L,PDupwindNth2B1,kmul(J31L,PDupwindNth3B1)));
      
      JacPDupwindNth1B2 = 
        kmadd(J11L,PDupwindNth1B2,kmadd(J21L,PDupwindNth2B2,kmul(J31L,PDupwindNth3B2)));
      
      JacPDupwindNth1B3 = 
        kmadd(J11L,PDupwindNth1B3,kmadd(J21L,PDupwindNth2B3,kmul(J31L,PDupwindNth3B3)));
      
      JacPDupwindNth1beta1 = 
        kmadd(J11L,PDupwindNth1beta1,kmadd(J21L,PDupwindNth2beta1,kmul(J31L,PDupwindNth3beta1)));
      
      JacPDupwindNth1beta2 = 
        kmadd(J11L,PDupwindNth1beta2,kmadd(J21L,PDupwindNth2beta2,kmul(J31L,PDupwindNth3beta2)));
      
      JacPDupwindNth1beta3 = 
        kmadd(J11L,PDupwindNth1beta3,kmadd(J21L,PDupwindNth2beta3,kmul(J31L,PDupwindNth3beta3)));
      
      JacPDupwindNth1gt11 = 
        kmadd(J11L,PDupwindNth1gt11,kmadd(J21L,PDupwindNth2gt11,kmul(J31L,PDupwindNth3gt11)));
      
      JacPDupwindNth1gt12 = 
        kmadd(J11L,PDupwindNth1gt12,kmadd(J21L,PDupwindNth2gt12,kmul(J31L,PDupwindNth3gt12)));
      
      JacPDupwindNth1gt13 = 
        kmadd(J11L,PDupwindNth1gt13,kmadd(J21L,PDupwindNth2gt13,kmul(J31L,PDupwindNth3gt13)));
      
      JacPDupwindNth1gt22 = 
        kmadd(J11L,PDupwindNth1gt22,kmadd(J21L,PDupwindNth2gt22,kmul(J31L,PDupwindNth3gt22)));
      
      JacPDupwindNth1gt23 = 
        kmadd(J11L,PDupwindNth1gt23,kmadd(J21L,PDupwindNth2gt23,kmul(J31L,PDupwindNth3gt23)));
      
      JacPDupwindNth1gt33 = 
        kmadd(J11L,PDupwindNth1gt33,kmadd(J21L,PDupwindNth2gt33,kmul(J31L,PDupwindNth3gt33)));
      
      JacPDupwindNth1phi = 
        kmadd(J11L,PDupwindNth1phi,kmadd(J21L,PDupwindNth2phi,kmul(J31L,PDupwindNth3phi)));
      
      JacPDupwindNth1trK = 
        kmadd(J11L,PDupwindNth1trK,kmadd(J21L,PDupwindNth2trK,kmul(J31L,PDupwindNth3trK)));
      
      JacPDupwindNth1Xt1 = 
        kmadd(J11L,PDupwindNth1Xt1,kmadd(J21L,PDupwindNth2Xt1,kmul(J31L,PDupwindNth3Xt1)));
      
      JacPDupwindNth1Xt2 = 
        kmadd(J11L,PDupwindNth1Xt2,kmadd(J21L,PDupwindNth2Xt2,kmul(J31L,PDupwindNth3Xt2)));
      
      JacPDupwindNth1Xt3 = 
        kmadd(J11L,PDupwindNth1Xt3,kmadd(J21L,PDupwindNth2Xt3,kmul(J31L,PDupwindNth3Xt3)));
      
      JacPDupwindNth2A = 
        kmadd(J12L,PDupwindNth1A,kmadd(J22L,PDupwindNth2A,kmul(J32L,PDupwindNth3A)));
      
      JacPDupwindNth2alpha = 
        kmadd(J12L,PDupwindNth1alpha,kmadd(J22L,PDupwindNth2alpha,kmul(J32L,PDupwindNth3alpha)));
      
      JacPDupwindNth2At11 = 
        kmadd(J12L,PDupwindNth1At11,kmadd(J22L,PDupwindNth2At11,kmul(J32L,PDupwindNth3At11)));
      
      JacPDupwindNth2At12 = 
        kmadd(J12L,PDupwindNth1At12,kmadd(J22L,PDupwindNth2At12,kmul(J32L,PDupwindNth3At12)));
      
      JacPDupwindNth2At13 = 
        kmadd(J12L,PDupwindNth1At13,kmadd(J22L,PDupwindNth2At13,kmul(J32L,PDupwindNth3At13)));
      
      JacPDupwindNth2At22 = 
        kmadd(J12L,PDupwindNth1At22,kmadd(J22L,PDupwindNth2At22,kmul(J32L,PDupwindNth3At22)));
      
      JacPDupwindNth2At23 = 
        kmadd(J12L,PDupwindNth1At23,kmadd(J22L,PDupwindNth2At23,kmul(J32L,PDupwindNth3At23)));
      
      JacPDupwindNth2At33 = 
        kmadd(J12L,PDupwindNth1At33,kmadd(J22L,PDupwindNth2At33,kmul(J32L,PDupwindNth3At33)));
      
      JacPDupwindNth2B1 = 
        kmadd(J12L,PDupwindNth1B1,kmadd(J22L,PDupwindNth2B1,kmul(J32L,PDupwindNth3B1)));
      
      JacPDupwindNth2B2 = 
        kmadd(J12L,PDupwindNth1B2,kmadd(J22L,PDupwindNth2B2,kmul(J32L,PDupwindNth3B2)));
      
      JacPDupwindNth2B3 = 
        kmadd(J12L,PDupwindNth1B3,kmadd(J22L,PDupwindNth2B3,kmul(J32L,PDupwindNth3B3)));
      
      JacPDupwindNth2beta1 = 
        kmadd(J12L,PDupwindNth1beta1,kmadd(J22L,PDupwindNth2beta1,kmul(J32L,PDupwindNth3beta1)));
      
      JacPDupwindNth2beta2 = 
        kmadd(J12L,PDupwindNth1beta2,kmadd(J22L,PDupwindNth2beta2,kmul(J32L,PDupwindNth3beta2)));
      
      JacPDupwindNth2beta3 = 
        kmadd(J12L,PDupwindNth1beta3,kmadd(J22L,PDupwindNth2beta3,kmul(J32L,PDupwindNth3beta3)));
      
      JacPDupwindNth2gt11 = 
        kmadd(J12L,PDupwindNth1gt11,kmadd(J22L,PDupwindNth2gt11,kmul(J32L,PDupwindNth3gt11)));
      
      JacPDupwindNth2gt12 = 
        kmadd(J12L,PDupwindNth1gt12,kmadd(J22L,PDupwindNth2gt12,kmul(J32L,PDupwindNth3gt12)));
      
      JacPDupwindNth2gt13 = 
        kmadd(J12L,PDupwindNth1gt13,kmadd(J22L,PDupwindNth2gt13,kmul(J32L,PDupwindNth3gt13)));
      
      JacPDupwindNth2gt22 = 
        kmadd(J12L,PDupwindNth1gt22,kmadd(J22L,PDupwindNth2gt22,kmul(J32L,PDupwindNth3gt22)));
      
      JacPDupwindNth2gt23 = 
        kmadd(J12L,PDupwindNth1gt23,kmadd(J22L,PDupwindNth2gt23,kmul(J32L,PDupwindNth3gt23)));
      
      JacPDupwindNth2gt33 = 
        kmadd(J12L,PDupwindNth1gt33,kmadd(J22L,PDupwindNth2gt33,kmul(J32L,PDupwindNth3gt33)));
      
      JacPDupwindNth2phi = 
        kmadd(J12L,PDupwindNth1phi,kmadd(J22L,PDupwindNth2phi,kmul(J32L,PDupwindNth3phi)));
      
      JacPDupwindNth2trK = 
        kmadd(J12L,PDupwindNth1trK,kmadd(J22L,PDupwindNth2trK,kmul(J32L,PDupwindNth3trK)));
      
      JacPDupwindNth2Xt1 = 
        kmadd(J12L,PDupwindNth1Xt1,kmadd(J22L,PDupwindNth2Xt1,kmul(J32L,PDupwindNth3Xt1)));
      
      JacPDupwindNth2Xt2 = 
        kmadd(J12L,PDupwindNth1Xt2,kmadd(J22L,PDupwindNth2Xt2,kmul(J32L,PDupwindNth3Xt2)));
      
      JacPDupwindNth2Xt3 = 
        kmadd(J12L,PDupwindNth1Xt3,kmadd(J22L,PDupwindNth2Xt3,kmul(J32L,PDupwindNth3Xt3)));
      
      JacPDupwindNth3A = 
        kmadd(J13L,PDupwindNth1A,kmadd(J23L,PDupwindNth2A,kmul(J33L,PDupwindNth3A)));
      
      JacPDupwindNth3alpha = 
        kmadd(J13L,PDupwindNth1alpha,kmadd(J23L,PDupwindNth2alpha,kmul(J33L,PDupwindNth3alpha)));
      
      JacPDupwindNth3At11 = 
        kmadd(J13L,PDupwindNth1At11,kmadd(J23L,PDupwindNth2At11,kmul(J33L,PDupwindNth3At11)));
      
      JacPDupwindNth3At12 = 
        kmadd(J13L,PDupwindNth1At12,kmadd(J23L,PDupwindNth2At12,kmul(J33L,PDupwindNth3At12)));
      
      JacPDupwindNth3At13 = 
        kmadd(J13L,PDupwindNth1At13,kmadd(J23L,PDupwindNth2At13,kmul(J33L,PDupwindNth3At13)));
      
      JacPDupwindNth3At22 = 
        kmadd(J13L,PDupwindNth1At22,kmadd(J23L,PDupwindNth2At22,kmul(J33L,PDupwindNth3At22)));
      
      JacPDupwindNth3At23 = 
        kmadd(J13L,PDupwindNth1At23,kmadd(J23L,PDupwindNth2At23,kmul(J33L,PDupwindNth3At23)));
      
      JacPDupwindNth3At33 = 
        kmadd(J13L,PDupwindNth1At33,kmadd(J23L,PDupwindNth2At33,kmul(J33L,PDupwindNth3At33)));
      
      JacPDupwindNth3B1 = 
        kmadd(J13L,PDupwindNth1B1,kmadd(J23L,PDupwindNth2B1,kmul(J33L,PDupwindNth3B1)));
      
      JacPDupwindNth3B2 = 
        kmadd(J13L,PDupwindNth1B2,kmadd(J23L,PDupwindNth2B2,kmul(J33L,PDupwindNth3B2)));
      
      JacPDupwindNth3B3 = 
        kmadd(J13L,PDupwindNth1B3,kmadd(J23L,PDupwindNth2B3,kmul(J33L,PDupwindNth3B3)));
      
      JacPDupwindNth3beta1 = 
        kmadd(J13L,PDupwindNth1beta1,kmadd(J23L,PDupwindNth2beta1,kmul(J33L,PDupwindNth3beta1)));
      
      JacPDupwindNth3beta2 = 
        kmadd(J13L,PDupwindNth1beta2,kmadd(J23L,PDupwindNth2beta2,kmul(J33L,PDupwindNth3beta2)));
      
      JacPDupwindNth3beta3 = 
        kmadd(J13L,PDupwindNth1beta3,kmadd(J23L,PDupwindNth2beta3,kmul(J33L,PDupwindNth3beta3)));
      
      JacPDupwindNth3gt11 = 
        kmadd(J13L,PDupwindNth1gt11,kmadd(J23L,PDupwindNth2gt11,kmul(J33L,PDupwindNth3gt11)));
      
      JacPDupwindNth3gt12 = 
        kmadd(J13L,PDupwindNth1gt12,kmadd(J23L,PDupwindNth2gt12,kmul(J33L,PDupwindNth3gt12)));
      
      JacPDupwindNth3gt13 = 
        kmadd(J13L,PDupwindNth1gt13,kmadd(J23L,PDupwindNth2gt13,kmul(J33L,PDupwindNth3gt13)));
      
      JacPDupwindNth3gt22 = 
        kmadd(J13L,PDupwindNth1gt22,kmadd(J23L,PDupwindNth2gt22,kmul(J33L,PDupwindNth3gt22)));
      
      JacPDupwindNth3gt23 = 
        kmadd(J13L,PDupwindNth1gt23,kmadd(J23L,PDupwindNth2gt23,kmul(J33L,PDupwindNth3gt23)));
      
      JacPDupwindNth3gt33 = 
        kmadd(J13L,PDupwindNth1gt33,kmadd(J23L,PDupwindNth2gt33,kmul(J33L,PDupwindNth3gt33)));
      
      JacPDupwindNth3phi = 
        kmadd(J13L,PDupwindNth1phi,kmadd(J23L,PDupwindNth2phi,kmul(J33L,PDupwindNth3phi)));
      
      JacPDupwindNth3trK = 
        kmadd(J13L,PDupwindNth1trK,kmadd(J23L,PDupwindNth2trK,kmul(J33L,PDupwindNth3trK)));
      
      JacPDupwindNth3Xt1 = 
        kmadd(J13L,PDupwindNth1Xt1,kmadd(J23L,PDupwindNth2Xt1,kmul(J33L,PDupwindNth3Xt1)));
      
      JacPDupwindNth3Xt2 = 
        kmadd(J13L,PDupwindNth1Xt2,kmadd(J23L,PDupwindNth2Xt2,kmul(J33L,PDupwindNth3Xt2)));
      
      JacPDupwindNth3Xt3 = 
        kmadd(J13L,PDupwindNth1Xt3,kmadd(J23L,PDupwindNth2Xt3,kmul(J33L,PDupwindNth3Xt3)));
    }
    else
    {
      JacPDupwindNth1A = PDupwindNth1A;
      
      JacPDupwindNth1alpha = PDupwindNth1alpha;
      
      JacPDupwindNth1At11 = PDupwindNth1At11;
      
      JacPDupwindNth1At12 = PDupwindNth1At12;
      
      JacPDupwindNth1At13 = PDupwindNth1At13;
      
      JacPDupwindNth1At22 = PDupwindNth1At22;
      
      JacPDupwindNth1At23 = PDupwindNth1At23;
      
      JacPDupwindNth1At33 = PDupwindNth1At33;
      
      JacPDupwindNth1B1 = PDupwindNth1B1;
      
      JacPDupwindNth1B2 = PDupwindNth1B2;
      
      JacPDupwindNth1B3 = PDupwindNth1B3;
      
      JacPDupwindNth1beta1 = PDupwindNth1beta1;
      
      JacPDupwindNth1beta2 = PDupwindNth1beta2;
      
      JacPDupwindNth1beta3 = PDupwindNth1beta3;
      
      JacPDupwindNth1gt11 = PDupwindNth1gt11;
      
      JacPDupwindNth1gt12 = PDupwindNth1gt12;
      
      JacPDupwindNth1gt13 = PDupwindNth1gt13;
      
      JacPDupwindNth1gt22 = PDupwindNth1gt22;
      
      JacPDupwindNth1gt23 = PDupwindNth1gt23;
      
      JacPDupwindNth1gt33 = PDupwindNth1gt33;
      
      JacPDupwindNth1phi = PDupwindNth1phi;
      
      JacPDupwindNth1trK = PDupwindNth1trK;
      
      JacPDupwindNth1Xt1 = PDupwindNth1Xt1;
      
      JacPDupwindNth1Xt2 = PDupwindNth1Xt2;
      
      JacPDupwindNth1Xt3 = PDupwindNth1Xt3;
      
      JacPDupwindNth2A = PDupwindNth2A;
      
      JacPDupwindNth2alpha = PDupwindNth2alpha;
      
      JacPDupwindNth2At11 = PDupwindNth2At11;
      
      JacPDupwindNth2At12 = PDupwindNth2At12;
      
      JacPDupwindNth2At13 = PDupwindNth2At13;
      
      JacPDupwindNth2At22 = PDupwindNth2At22;
      
      JacPDupwindNth2At23 = PDupwindNth2At23;
      
      JacPDupwindNth2At33 = PDupwindNth2At33;
      
      JacPDupwindNth2B1 = PDupwindNth2B1;
      
      JacPDupwindNth2B2 = PDupwindNth2B2;
      
      JacPDupwindNth2B3 = PDupwindNth2B3;
      
      JacPDupwindNth2beta1 = PDupwindNth2beta1;
      
      JacPDupwindNth2beta2 = PDupwindNth2beta2;
      
      JacPDupwindNth2beta3 = PDupwindNth2beta3;
      
      JacPDupwindNth2gt11 = PDupwindNth2gt11;
      
      JacPDupwindNth2gt12 = PDupwindNth2gt12;
      
      JacPDupwindNth2gt13 = PDupwindNth2gt13;
      
      JacPDupwindNth2gt22 = PDupwindNth2gt22;
      
      JacPDupwindNth2gt23 = PDupwindNth2gt23;
      
      JacPDupwindNth2gt33 = PDupwindNth2gt33;
      
      JacPDupwindNth2phi = PDupwindNth2phi;
      
      JacPDupwindNth2trK = PDupwindNth2trK;
      
      JacPDupwindNth2Xt1 = PDupwindNth2Xt1;
      
      JacPDupwindNth2Xt2 = PDupwindNth2Xt2;
      
      JacPDupwindNth2Xt3 = PDupwindNth2Xt3;
      
      JacPDupwindNth3A = PDupwindNth3A;
      
      JacPDupwindNth3alpha = PDupwindNth3alpha;
      
      JacPDupwindNth3At11 = PDupwindNth3At11;
      
      JacPDupwindNth3At12 = PDupwindNth3At12;
      
      JacPDupwindNth3At13 = PDupwindNth3At13;
      
      JacPDupwindNth3At22 = PDupwindNth3At22;
      
      JacPDupwindNth3At23 = PDupwindNth3At23;
      
      JacPDupwindNth3At33 = PDupwindNth3At33;
      
      JacPDupwindNth3B1 = PDupwindNth3B1;
      
      JacPDupwindNth3B2 = PDupwindNth3B2;
      
      JacPDupwindNth3B3 = PDupwindNth3B3;
      
      JacPDupwindNth3beta1 = PDupwindNth3beta1;
      
      JacPDupwindNth3beta2 = PDupwindNth3beta2;
      
      JacPDupwindNth3beta3 = PDupwindNth3beta3;
      
      JacPDupwindNth3gt11 = PDupwindNth3gt11;
      
      JacPDupwindNth3gt12 = PDupwindNth3gt12;
      
      JacPDupwindNth3gt13 = PDupwindNth3gt13;
      
      JacPDupwindNth3gt22 = PDupwindNth3gt22;
      
      JacPDupwindNth3gt23 = PDupwindNth3gt23;
      
      JacPDupwindNth3gt33 = PDupwindNth3gt33;
      
      JacPDupwindNth3phi = PDupwindNth3phi;
      
      JacPDupwindNth3trK = PDupwindNth3trK;
      
      JacPDupwindNth3Xt1 = PDupwindNth3Xt1;
      
      JacPDupwindNth3Xt2 = PDupwindNth3Xt2;
      
      JacPDupwindNth3Xt3 = PDupwindNth3Xt3;
    }
    
    phirhsL = 
      kadd(phirhsL,kmadd(beta1L,JacPDupwindNth1phi,kmadd(beta2L,JacPDupwindNth2phi,kmul(beta3L,JacPDupwindNth3phi))));
    
    gt11rhsL = 
      kadd(gt11rhsL,kmadd(beta1L,JacPDupwindNth1gt11,kmadd(beta2L,JacPDupwindNth2gt11,kmul(beta3L,JacPDupwindNth3gt11))));
    
    gt12rhsL = 
      kadd(gt12rhsL,kmadd(beta1L,JacPDupwindNth1gt12,kmadd(beta2L,JacPDupwindNth2gt12,kmul(beta3L,JacPDupwindNth3gt12))));
    
    gt13rhsL = 
      kadd(gt13rhsL,kmadd(beta1L,JacPDupwindNth1gt13,kmadd(beta2L,JacPDupwindNth2gt13,kmul(beta3L,JacPDupwindNth3gt13))));
    
    gt22rhsL = 
      kadd(gt22rhsL,kmadd(beta1L,JacPDupwindNth1gt22,kmadd(beta2L,JacPDupwindNth2gt22,kmul(beta3L,JacPDupwindNth3gt22))));
    
    gt23rhsL = 
      kadd(gt23rhsL,kmadd(beta1L,JacPDupwindNth1gt23,kmadd(beta2L,JacPDupwindNth2gt23,kmul(beta3L,JacPDupwindNth3gt23))));
    
    gt33rhsL = 
      kadd(gt33rhsL,kmadd(beta1L,JacPDupwindNth1gt33,kmadd(beta2L,JacPDupwindNth2gt33,kmul(beta3L,JacPDupwindNth3gt33))));
    
    Xt1rhsL = 
      kadd(Xt1rhsL,kmadd(beta1L,JacPDupwindNth1Xt1,kmadd(beta2L,JacPDupwindNth2Xt1,kmul(beta3L,JacPDupwindNth3Xt1))));
    
    Xt2rhsL = 
      kadd(Xt2rhsL,kmadd(beta1L,JacPDupwindNth1Xt2,kmadd(beta2L,JacPDupwindNth2Xt2,kmul(beta3L,JacPDupwindNth3Xt2))));
    
    Xt3rhsL = 
      kadd(Xt3rhsL,kmadd(beta1L,JacPDupwindNth1Xt3,kmadd(beta2L,JacPDupwindNth2Xt3,kmul(beta3L,JacPDupwindNth3Xt3))));
    
    trKrhsL = 
      kadd(trKrhsL,kmadd(beta1L,JacPDupwindNth1trK,kmadd(beta2L,JacPDupwindNth2trK,kmul(beta3L,JacPDupwindNth3trK))));
    
    At11rhsL = 
      kadd(At11rhsL,kmadd(beta1L,JacPDupwindNth1At11,kmadd(beta2L,JacPDupwindNth2At11,kmul(beta3L,JacPDupwindNth3At11))));
    
    At12rhsL = 
      kadd(At12rhsL,kmadd(beta1L,JacPDupwindNth1At12,kmadd(beta2L,JacPDupwindNth2At12,kmul(beta3L,JacPDupwindNth3At12))));
    
    At13rhsL = 
      kadd(At13rhsL,kmadd(beta1L,JacPDupwindNth1At13,kmadd(beta2L,JacPDupwindNth2At13,kmul(beta3L,JacPDupwindNth3At13))));
    
    At22rhsL = 
      kadd(At22rhsL,kmadd(beta1L,JacPDupwindNth1At22,kmadd(beta2L,JacPDupwindNth2At22,kmul(beta3L,JacPDupwindNth3At22))));
    
    At23rhsL = 
      kadd(At23rhsL,kmadd(beta1L,JacPDupwindNth1At23,kmadd(beta2L,JacPDupwindNth2At23,kmul(beta3L,JacPDupwindNth3At23))));
    
    At33rhsL = 
      kadd(At33rhsL,kmadd(beta1L,JacPDupwindNth1At33,kmadd(beta2L,JacPDupwindNth2At33,kmul(beta3L,JacPDupwindNth3At33))));
    
    alpharhsL = 
      kmadd(kmadd(beta1L,JacPDupwindNth1alpha,kmadd(beta2L,JacPDupwindNth2alpha,kmul(beta3L,JacPDupwindNth3alpha))),ToReal(LapseAdvectionCoeff),alpharhsL);
    
    ArhsL = 
      kmadd(ToReal(LapseACoeff),kmadd(beta1L,kmadd(ksub(JacPDupwindNth1A,JacPDupwindNth1trK),ToReal(LapseAdvectionCoeff),JacPDupwindNth1trK),kmadd(beta2L,kmadd(ksub(JacPDupwindNth2A,JacPDupwindNth2trK),ToReal(LapseAdvectionCoeff),JacPDupwindNth2trK),kmul(beta3L,kmadd(ksub(JacPDupwindNth3A,JacPDupwindNth3trK),ToReal(LapseAdvectionCoeff),JacPDupwindNth3trK)))),ArhsL);
    
    beta1rhsL = 
      kmadd(kmadd(beta1L,JacPDupwindNth1beta1,kmadd(beta2L,JacPDupwindNth2beta1,kmul(beta3L,JacPDupwindNth3beta1))),ToReal(ShiftAdvectionCoeff),beta1rhsL);
    
    beta2rhsL = 
      kmadd(kmadd(beta1L,JacPDupwindNth1beta2,kmadd(beta2L,JacPDupwindNth2beta2,kmul(beta3L,JacPDupwindNth3beta2))),ToReal(ShiftAdvectionCoeff),beta2rhsL);
    
    beta3rhsL = 
      kmadd(kmadd(beta1L,JacPDupwindNth1beta3,kmadd(beta2L,JacPDupwindNth2beta3,kmul(beta3L,JacPDupwindNth3beta3))),ToReal(ShiftAdvectionCoeff),beta3rhsL);
    
    B1rhsL = 
      kmadd(kmadd(beta1L,kmadd(ksub(JacPDupwindNth1B1,JacPDupwindNth1Xt1),ToReal(ShiftAdvectionCoeff),JacPDupwindNth1Xt1),kmadd(beta2L,kmadd(ksub(JacPDupwindNth2B1,JacPDupwindNth2Xt1),ToReal(ShiftAdvectionCoeff),JacPDupwindNth2Xt1),kmul(beta3L,kmadd(ksub(JacPDupwindNth3B1,JacPDupwindNth3Xt1),ToReal(ShiftAdvectionCoeff),JacPDupwindNth3Xt1)))),ToReal(ShiftBCoeff),B1rhsL);
    
    B2rhsL = 
      kmadd(kmadd(beta1L,kmadd(ksub(JacPDupwindNth1B2,JacPDupwindNth1Xt2),ToReal(ShiftAdvectionCoeff),JacPDupwindNth1Xt2),kmadd(beta2L,kmadd(ksub(JacPDupwindNth2B2,JacPDupwindNth2Xt2),ToReal(ShiftAdvectionCoeff),JacPDupwindNth2Xt2),kmul(beta3L,kmadd(ksub(JacPDupwindNth3B2,JacPDupwindNth3Xt2),ToReal(ShiftAdvectionCoeff),JacPDupwindNth3Xt2)))),ToReal(ShiftBCoeff),B2rhsL);
    
    B3rhsL = 
      kmadd(kmadd(beta1L,kmadd(ksub(JacPDupwindNth1B3,JacPDupwindNth1Xt3),ToReal(ShiftAdvectionCoeff),JacPDupwindNth1Xt3),kmadd(beta2L,kmadd(ksub(JacPDupwindNth2B3,JacPDupwindNth2Xt3),ToReal(ShiftAdvectionCoeff),JacPDupwindNth2Xt3),kmul(beta3L,kmadd(ksub(JacPDupwindNth3B3,JacPDupwindNth3Xt3),ToReal(ShiftAdvectionCoeff),JacPDupwindNth3Xt3)))),ToReal(ShiftBCoeff),B3rhsL);
    
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
  LC_ENDLOOP3VEC (ML_BSSN_UPW_Advect);
}

extern "C" void ML_BSSN_UPW_Advect(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_UPW_Advect_Body");
  }
  
  if (cctk_iteration % ML_BSSN_UPW_Advect_calc_every != ML_BSSN_UPW_Advect_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_BSSN_UPW::ML_curv","ML_BSSN_UPW::ML_curvrhs","ML_BSSN_UPW::ML_dtlapse","ML_BSSN_UPW::ML_dtlapserhs","ML_BSSN_UPW::ML_dtshift","ML_BSSN_UPW::ML_dtshiftrhs","ML_BSSN_UPW::ML_Gamma","ML_BSSN_UPW::ML_Gammarhs","ML_BSSN_UPW::ML_lapse","ML_BSSN_UPW::ML_lapserhs","ML_BSSN_UPW::ML_log_confac","ML_BSSN_UPW::ML_log_confacrhs","ML_BSSN_UPW::ML_metric","ML_BSSN_UPW::ML_metricrhs","ML_BSSN_UPW::ML_shift","ML_BSSN_UPW::ML_shiftrhs","ML_BSSN_UPW::ML_trace_curv","ML_BSSN_UPW::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_UPW_Advect", 18, groups);
  
  switch(fdOrder)
  {
    case 2:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_Advect", 2, 2, 2);
      break;
    
    case 4:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_Advect", 3, 3, 3);
      break;
    
    case 6:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_Advect", 4, 4, 4);
      break;
    
    case 8:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_Advect", 5, 5, 5);
      break;
  }
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_UPW_Advect_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_UPW_Advect_Body");
  }
}
