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

extern "C" void ML_BSSN_Advect_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_curvrhs.");
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

static void ML_BSSN_Advect_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  LC_LOOP3VEC (ML_BSSN_Advect,
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
    CCTK_REAL_VEC PDupwindNthAnti1A;
    CCTK_REAL_VEC PDupwindNthSymm1A;
    CCTK_REAL_VEC PDupwindNthAnti2A;
    CCTK_REAL_VEC PDupwindNthSymm2A;
    CCTK_REAL_VEC PDupwindNthAnti3A;
    CCTK_REAL_VEC PDupwindNthSymm3A;
    CCTK_REAL_VEC PDupwindNthAnti1alpha;
    CCTK_REAL_VEC PDupwindNthSymm1alpha;
    CCTK_REAL_VEC PDupwindNthAnti2alpha;
    CCTK_REAL_VEC PDupwindNthSymm2alpha;
    CCTK_REAL_VEC PDupwindNthAnti3alpha;
    CCTK_REAL_VEC PDupwindNthSymm3alpha;
    CCTK_REAL_VEC PDupwindNthAnti1At11;
    CCTK_REAL_VEC PDupwindNthSymm1At11;
    CCTK_REAL_VEC PDupwindNthAnti2At11;
    CCTK_REAL_VEC PDupwindNthSymm2At11;
    CCTK_REAL_VEC PDupwindNthAnti3At11;
    CCTK_REAL_VEC PDupwindNthSymm3At11;
    CCTK_REAL_VEC PDupwindNthAnti1At12;
    CCTK_REAL_VEC PDupwindNthSymm1At12;
    CCTK_REAL_VEC PDupwindNthAnti2At12;
    CCTK_REAL_VEC PDupwindNthSymm2At12;
    CCTK_REAL_VEC PDupwindNthAnti3At12;
    CCTK_REAL_VEC PDupwindNthSymm3At12;
    CCTK_REAL_VEC PDupwindNthAnti1At13;
    CCTK_REAL_VEC PDupwindNthSymm1At13;
    CCTK_REAL_VEC PDupwindNthAnti2At13;
    CCTK_REAL_VEC PDupwindNthSymm2At13;
    CCTK_REAL_VEC PDupwindNthAnti3At13;
    CCTK_REAL_VEC PDupwindNthSymm3At13;
    CCTK_REAL_VEC PDupwindNthAnti1At22;
    CCTK_REAL_VEC PDupwindNthSymm1At22;
    CCTK_REAL_VEC PDupwindNthAnti2At22;
    CCTK_REAL_VEC PDupwindNthSymm2At22;
    CCTK_REAL_VEC PDupwindNthAnti3At22;
    CCTK_REAL_VEC PDupwindNthSymm3At22;
    CCTK_REAL_VEC PDupwindNthAnti1At23;
    CCTK_REAL_VEC PDupwindNthSymm1At23;
    CCTK_REAL_VEC PDupwindNthAnti2At23;
    CCTK_REAL_VEC PDupwindNthSymm2At23;
    CCTK_REAL_VEC PDupwindNthAnti3At23;
    CCTK_REAL_VEC PDupwindNthSymm3At23;
    CCTK_REAL_VEC PDupwindNthAnti1At33;
    CCTK_REAL_VEC PDupwindNthSymm1At33;
    CCTK_REAL_VEC PDupwindNthAnti2At33;
    CCTK_REAL_VEC PDupwindNthSymm2At33;
    CCTK_REAL_VEC PDupwindNthAnti3At33;
    CCTK_REAL_VEC PDupwindNthSymm3At33;
    CCTK_REAL_VEC PDupwindNthAnti1B1;
    CCTK_REAL_VEC PDupwindNthSymm1B1;
    CCTK_REAL_VEC PDupwindNthAnti2B1;
    CCTK_REAL_VEC PDupwindNthSymm2B1;
    CCTK_REAL_VEC PDupwindNthAnti3B1;
    CCTK_REAL_VEC PDupwindNthSymm3B1;
    CCTK_REAL_VEC PDupwindNthAnti1B2;
    CCTK_REAL_VEC PDupwindNthSymm1B2;
    CCTK_REAL_VEC PDupwindNthAnti2B2;
    CCTK_REAL_VEC PDupwindNthSymm2B2;
    CCTK_REAL_VEC PDupwindNthAnti3B2;
    CCTK_REAL_VEC PDupwindNthSymm3B2;
    CCTK_REAL_VEC PDupwindNthAnti1B3;
    CCTK_REAL_VEC PDupwindNthSymm1B3;
    CCTK_REAL_VEC PDupwindNthAnti2B3;
    CCTK_REAL_VEC PDupwindNthSymm2B3;
    CCTK_REAL_VEC PDupwindNthAnti3B3;
    CCTK_REAL_VEC PDupwindNthSymm3B3;
    CCTK_REAL_VEC PDupwindNthAnti1beta1;
    CCTK_REAL_VEC PDupwindNthSymm1beta1;
    CCTK_REAL_VEC PDupwindNthAnti2beta1;
    CCTK_REAL_VEC PDupwindNthSymm2beta1;
    CCTK_REAL_VEC PDupwindNthAnti3beta1;
    CCTK_REAL_VEC PDupwindNthSymm3beta1;
    CCTK_REAL_VEC PDupwindNthAnti1beta2;
    CCTK_REAL_VEC PDupwindNthSymm1beta2;
    CCTK_REAL_VEC PDupwindNthAnti2beta2;
    CCTK_REAL_VEC PDupwindNthSymm2beta2;
    CCTK_REAL_VEC PDupwindNthAnti3beta2;
    CCTK_REAL_VEC PDupwindNthSymm3beta2;
    CCTK_REAL_VEC PDupwindNthAnti1beta3;
    CCTK_REAL_VEC PDupwindNthSymm1beta3;
    CCTK_REAL_VEC PDupwindNthAnti2beta3;
    CCTK_REAL_VEC PDupwindNthSymm2beta3;
    CCTK_REAL_VEC PDupwindNthAnti3beta3;
    CCTK_REAL_VEC PDupwindNthSymm3beta3;
    CCTK_REAL_VEC PDupwindNthAnti1gt11;
    CCTK_REAL_VEC PDupwindNthSymm1gt11;
    CCTK_REAL_VEC PDupwindNthAnti2gt11;
    CCTK_REAL_VEC PDupwindNthSymm2gt11;
    CCTK_REAL_VEC PDupwindNthAnti3gt11;
    CCTK_REAL_VEC PDupwindNthSymm3gt11;
    CCTK_REAL_VEC PDupwindNthAnti1gt12;
    CCTK_REAL_VEC PDupwindNthSymm1gt12;
    CCTK_REAL_VEC PDupwindNthAnti2gt12;
    CCTK_REAL_VEC PDupwindNthSymm2gt12;
    CCTK_REAL_VEC PDupwindNthAnti3gt12;
    CCTK_REAL_VEC PDupwindNthSymm3gt12;
    CCTK_REAL_VEC PDupwindNthAnti1gt13;
    CCTK_REAL_VEC PDupwindNthSymm1gt13;
    CCTK_REAL_VEC PDupwindNthAnti2gt13;
    CCTK_REAL_VEC PDupwindNthSymm2gt13;
    CCTK_REAL_VEC PDupwindNthAnti3gt13;
    CCTK_REAL_VEC PDupwindNthSymm3gt13;
    CCTK_REAL_VEC PDupwindNthAnti1gt22;
    CCTK_REAL_VEC PDupwindNthSymm1gt22;
    CCTK_REAL_VEC PDupwindNthAnti2gt22;
    CCTK_REAL_VEC PDupwindNthSymm2gt22;
    CCTK_REAL_VEC PDupwindNthAnti3gt22;
    CCTK_REAL_VEC PDupwindNthSymm3gt22;
    CCTK_REAL_VEC PDupwindNthAnti1gt23;
    CCTK_REAL_VEC PDupwindNthSymm1gt23;
    CCTK_REAL_VEC PDupwindNthAnti2gt23;
    CCTK_REAL_VEC PDupwindNthSymm2gt23;
    CCTK_REAL_VEC PDupwindNthAnti3gt23;
    CCTK_REAL_VEC PDupwindNthSymm3gt23;
    CCTK_REAL_VEC PDupwindNthAnti1gt33;
    CCTK_REAL_VEC PDupwindNthSymm1gt33;
    CCTK_REAL_VEC PDupwindNthAnti2gt33;
    CCTK_REAL_VEC PDupwindNthSymm2gt33;
    CCTK_REAL_VEC PDupwindNthAnti3gt33;
    CCTK_REAL_VEC PDupwindNthSymm3gt33;
    CCTK_REAL_VEC PDupwindNthAnti1phi;
    CCTK_REAL_VEC PDupwindNthSymm1phi;
    CCTK_REAL_VEC PDupwindNthAnti2phi;
    CCTK_REAL_VEC PDupwindNthSymm2phi;
    CCTK_REAL_VEC PDupwindNthAnti3phi;
    CCTK_REAL_VEC PDupwindNthSymm3phi;
    CCTK_REAL_VEC PDupwindNthAnti1trK;
    CCTK_REAL_VEC PDupwindNthSymm1trK;
    CCTK_REAL_VEC PDupwindNthAnti2trK;
    CCTK_REAL_VEC PDupwindNthSymm2trK;
    CCTK_REAL_VEC PDupwindNthAnti3trK;
    CCTK_REAL_VEC PDupwindNthSymm3trK;
    CCTK_REAL_VEC PDupwindNthAnti1Xt1;
    CCTK_REAL_VEC PDupwindNthSymm1Xt1;
    CCTK_REAL_VEC PDupwindNthAnti2Xt1;
    CCTK_REAL_VEC PDupwindNthSymm2Xt1;
    CCTK_REAL_VEC PDupwindNthAnti3Xt1;
    CCTK_REAL_VEC PDupwindNthSymm3Xt1;
    CCTK_REAL_VEC PDupwindNthAnti1Xt2;
    CCTK_REAL_VEC PDupwindNthSymm1Xt2;
    CCTK_REAL_VEC PDupwindNthAnti2Xt2;
    CCTK_REAL_VEC PDupwindNthSymm2Xt2;
    CCTK_REAL_VEC PDupwindNthAnti3Xt2;
    CCTK_REAL_VEC PDupwindNthSymm3Xt2;
    CCTK_REAL_VEC PDupwindNthAnti1Xt3;
    CCTK_REAL_VEC PDupwindNthSymm1Xt3;
    CCTK_REAL_VEC PDupwindNthAnti2Xt3;
    CCTK_REAL_VEC PDupwindNthSymm2Xt3;
    CCTK_REAL_VEC PDupwindNthAnti3Xt3;
    CCTK_REAL_VEC PDupwindNthSymm3Xt3;
    
    switch(fdOrder)
    {
      case 2:
        PDupwindNthAnti1A = PDupwindNthAntifdOrder21(&A[index]);
        PDupwindNthSymm1A = PDupwindNthSymmfdOrder21(&A[index]);
        PDupwindNthAnti2A = PDupwindNthAntifdOrder22(&A[index]);
        PDupwindNthSymm2A = PDupwindNthSymmfdOrder22(&A[index]);
        PDupwindNthAnti3A = PDupwindNthAntifdOrder23(&A[index]);
        PDupwindNthSymm3A = PDupwindNthSymmfdOrder23(&A[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder21(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder21(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder22(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder22(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder23(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder23(&alpha[index]);
        PDupwindNthAnti1At11 = PDupwindNthAntifdOrder21(&At11[index]);
        PDupwindNthSymm1At11 = PDupwindNthSymmfdOrder21(&At11[index]);
        PDupwindNthAnti2At11 = PDupwindNthAntifdOrder22(&At11[index]);
        PDupwindNthSymm2At11 = PDupwindNthSymmfdOrder22(&At11[index]);
        PDupwindNthAnti3At11 = PDupwindNthAntifdOrder23(&At11[index]);
        PDupwindNthSymm3At11 = PDupwindNthSymmfdOrder23(&At11[index]);
        PDupwindNthAnti1At12 = PDupwindNthAntifdOrder21(&At12[index]);
        PDupwindNthSymm1At12 = PDupwindNthSymmfdOrder21(&At12[index]);
        PDupwindNthAnti2At12 = PDupwindNthAntifdOrder22(&At12[index]);
        PDupwindNthSymm2At12 = PDupwindNthSymmfdOrder22(&At12[index]);
        PDupwindNthAnti3At12 = PDupwindNthAntifdOrder23(&At12[index]);
        PDupwindNthSymm3At12 = PDupwindNthSymmfdOrder23(&At12[index]);
        PDupwindNthAnti1At13 = PDupwindNthAntifdOrder21(&At13[index]);
        PDupwindNthSymm1At13 = PDupwindNthSymmfdOrder21(&At13[index]);
        PDupwindNthAnti2At13 = PDupwindNthAntifdOrder22(&At13[index]);
        PDupwindNthSymm2At13 = PDupwindNthSymmfdOrder22(&At13[index]);
        PDupwindNthAnti3At13 = PDupwindNthAntifdOrder23(&At13[index]);
        PDupwindNthSymm3At13 = PDupwindNthSymmfdOrder23(&At13[index]);
        PDupwindNthAnti1At22 = PDupwindNthAntifdOrder21(&At22[index]);
        PDupwindNthSymm1At22 = PDupwindNthSymmfdOrder21(&At22[index]);
        PDupwindNthAnti2At22 = PDupwindNthAntifdOrder22(&At22[index]);
        PDupwindNthSymm2At22 = PDupwindNthSymmfdOrder22(&At22[index]);
        PDupwindNthAnti3At22 = PDupwindNthAntifdOrder23(&At22[index]);
        PDupwindNthSymm3At22 = PDupwindNthSymmfdOrder23(&At22[index]);
        PDupwindNthAnti1At23 = PDupwindNthAntifdOrder21(&At23[index]);
        PDupwindNthSymm1At23 = PDupwindNthSymmfdOrder21(&At23[index]);
        PDupwindNthAnti2At23 = PDupwindNthAntifdOrder22(&At23[index]);
        PDupwindNthSymm2At23 = PDupwindNthSymmfdOrder22(&At23[index]);
        PDupwindNthAnti3At23 = PDupwindNthAntifdOrder23(&At23[index]);
        PDupwindNthSymm3At23 = PDupwindNthSymmfdOrder23(&At23[index]);
        PDupwindNthAnti1At33 = PDupwindNthAntifdOrder21(&At33[index]);
        PDupwindNthSymm1At33 = PDupwindNthSymmfdOrder21(&At33[index]);
        PDupwindNthAnti2At33 = PDupwindNthAntifdOrder22(&At33[index]);
        PDupwindNthSymm2At33 = PDupwindNthSymmfdOrder22(&At33[index]);
        PDupwindNthAnti3At33 = PDupwindNthAntifdOrder23(&At33[index]);
        PDupwindNthSymm3At33 = PDupwindNthSymmfdOrder23(&At33[index]);
        PDupwindNthAnti1B1 = PDupwindNthAntifdOrder21(&B1[index]);
        PDupwindNthSymm1B1 = PDupwindNthSymmfdOrder21(&B1[index]);
        PDupwindNthAnti2B1 = PDupwindNthAntifdOrder22(&B1[index]);
        PDupwindNthSymm2B1 = PDupwindNthSymmfdOrder22(&B1[index]);
        PDupwindNthAnti3B1 = PDupwindNthAntifdOrder23(&B1[index]);
        PDupwindNthSymm3B1 = PDupwindNthSymmfdOrder23(&B1[index]);
        PDupwindNthAnti1B2 = PDupwindNthAntifdOrder21(&B2[index]);
        PDupwindNthSymm1B2 = PDupwindNthSymmfdOrder21(&B2[index]);
        PDupwindNthAnti2B2 = PDupwindNthAntifdOrder22(&B2[index]);
        PDupwindNthSymm2B2 = PDupwindNthSymmfdOrder22(&B2[index]);
        PDupwindNthAnti3B2 = PDupwindNthAntifdOrder23(&B2[index]);
        PDupwindNthSymm3B2 = PDupwindNthSymmfdOrder23(&B2[index]);
        PDupwindNthAnti1B3 = PDupwindNthAntifdOrder21(&B3[index]);
        PDupwindNthSymm1B3 = PDupwindNthSymmfdOrder21(&B3[index]);
        PDupwindNthAnti2B3 = PDupwindNthAntifdOrder22(&B3[index]);
        PDupwindNthSymm2B3 = PDupwindNthSymmfdOrder22(&B3[index]);
        PDupwindNthAnti3B3 = PDupwindNthAntifdOrder23(&B3[index]);
        PDupwindNthSymm3B3 = PDupwindNthSymmfdOrder23(&B3[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder21(&beta1[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder21(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder22(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder22(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder23(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder23(&beta1[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder21(&beta2[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder21(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder22(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder22(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder23(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder23(&beta2[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder21(&beta3[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder21(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder22(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder22(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder23(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder23(&beta3[index]);
        PDupwindNthAnti1gt11 = PDupwindNthAntifdOrder21(&gt11[index]);
        PDupwindNthSymm1gt11 = PDupwindNthSymmfdOrder21(&gt11[index]);
        PDupwindNthAnti2gt11 = PDupwindNthAntifdOrder22(&gt11[index]);
        PDupwindNthSymm2gt11 = PDupwindNthSymmfdOrder22(&gt11[index]);
        PDupwindNthAnti3gt11 = PDupwindNthAntifdOrder23(&gt11[index]);
        PDupwindNthSymm3gt11 = PDupwindNthSymmfdOrder23(&gt11[index]);
        PDupwindNthAnti1gt12 = PDupwindNthAntifdOrder21(&gt12[index]);
        PDupwindNthSymm1gt12 = PDupwindNthSymmfdOrder21(&gt12[index]);
        PDupwindNthAnti2gt12 = PDupwindNthAntifdOrder22(&gt12[index]);
        PDupwindNthSymm2gt12 = PDupwindNthSymmfdOrder22(&gt12[index]);
        PDupwindNthAnti3gt12 = PDupwindNthAntifdOrder23(&gt12[index]);
        PDupwindNthSymm3gt12 = PDupwindNthSymmfdOrder23(&gt12[index]);
        PDupwindNthAnti1gt13 = PDupwindNthAntifdOrder21(&gt13[index]);
        PDupwindNthSymm1gt13 = PDupwindNthSymmfdOrder21(&gt13[index]);
        PDupwindNthAnti2gt13 = PDupwindNthAntifdOrder22(&gt13[index]);
        PDupwindNthSymm2gt13 = PDupwindNthSymmfdOrder22(&gt13[index]);
        PDupwindNthAnti3gt13 = PDupwindNthAntifdOrder23(&gt13[index]);
        PDupwindNthSymm3gt13 = PDupwindNthSymmfdOrder23(&gt13[index]);
        PDupwindNthAnti1gt22 = PDupwindNthAntifdOrder21(&gt22[index]);
        PDupwindNthSymm1gt22 = PDupwindNthSymmfdOrder21(&gt22[index]);
        PDupwindNthAnti2gt22 = PDupwindNthAntifdOrder22(&gt22[index]);
        PDupwindNthSymm2gt22 = PDupwindNthSymmfdOrder22(&gt22[index]);
        PDupwindNthAnti3gt22 = PDupwindNthAntifdOrder23(&gt22[index]);
        PDupwindNthSymm3gt22 = PDupwindNthSymmfdOrder23(&gt22[index]);
        PDupwindNthAnti1gt23 = PDupwindNthAntifdOrder21(&gt23[index]);
        PDupwindNthSymm1gt23 = PDupwindNthSymmfdOrder21(&gt23[index]);
        PDupwindNthAnti2gt23 = PDupwindNthAntifdOrder22(&gt23[index]);
        PDupwindNthSymm2gt23 = PDupwindNthSymmfdOrder22(&gt23[index]);
        PDupwindNthAnti3gt23 = PDupwindNthAntifdOrder23(&gt23[index]);
        PDupwindNthSymm3gt23 = PDupwindNthSymmfdOrder23(&gt23[index]);
        PDupwindNthAnti1gt33 = PDupwindNthAntifdOrder21(&gt33[index]);
        PDupwindNthSymm1gt33 = PDupwindNthSymmfdOrder21(&gt33[index]);
        PDupwindNthAnti2gt33 = PDupwindNthAntifdOrder22(&gt33[index]);
        PDupwindNthSymm2gt33 = PDupwindNthSymmfdOrder22(&gt33[index]);
        PDupwindNthAnti3gt33 = PDupwindNthAntifdOrder23(&gt33[index]);
        PDupwindNthSymm3gt33 = PDupwindNthSymmfdOrder23(&gt33[index]);
        PDupwindNthAnti1phi = PDupwindNthAntifdOrder21(&phi[index]);
        PDupwindNthSymm1phi = PDupwindNthSymmfdOrder21(&phi[index]);
        PDupwindNthAnti2phi = PDupwindNthAntifdOrder22(&phi[index]);
        PDupwindNthSymm2phi = PDupwindNthSymmfdOrder22(&phi[index]);
        PDupwindNthAnti3phi = PDupwindNthAntifdOrder23(&phi[index]);
        PDupwindNthSymm3phi = PDupwindNthSymmfdOrder23(&phi[index]);
        PDupwindNthAnti1trK = PDupwindNthAntifdOrder21(&trK[index]);
        PDupwindNthSymm1trK = PDupwindNthSymmfdOrder21(&trK[index]);
        PDupwindNthAnti2trK = PDupwindNthAntifdOrder22(&trK[index]);
        PDupwindNthSymm2trK = PDupwindNthSymmfdOrder22(&trK[index]);
        PDupwindNthAnti3trK = PDupwindNthAntifdOrder23(&trK[index]);
        PDupwindNthSymm3trK = PDupwindNthSymmfdOrder23(&trK[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder21(&Xt1[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder21(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder22(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder22(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder23(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder23(&Xt1[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder21(&Xt2[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder21(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder22(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder22(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder23(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder23(&Xt2[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder21(&Xt3[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder21(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder22(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder22(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder23(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder23(&Xt3[index]);
        break;
      
      case 4:
        PDupwindNthAnti1A = PDupwindNthAntifdOrder41(&A[index]);
        PDupwindNthSymm1A = PDupwindNthSymmfdOrder41(&A[index]);
        PDupwindNthAnti2A = PDupwindNthAntifdOrder42(&A[index]);
        PDupwindNthSymm2A = PDupwindNthSymmfdOrder42(&A[index]);
        PDupwindNthAnti3A = PDupwindNthAntifdOrder43(&A[index]);
        PDupwindNthSymm3A = PDupwindNthSymmfdOrder43(&A[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder41(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder41(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder42(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder42(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder43(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder43(&alpha[index]);
        PDupwindNthAnti1At11 = PDupwindNthAntifdOrder41(&At11[index]);
        PDupwindNthSymm1At11 = PDupwindNthSymmfdOrder41(&At11[index]);
        PDupwindNthAnti2At11 = PDupwindNthAntifdOrder42(&At11[index]);
        PDupwindNthSymm2At11 = PDupwindNthSymmfdOrder42(&At11[index]);
        PDupwindNthAnti3At11 = PDupwindNthAntifdOrder43(&At11[index]);
        PDupwindNthSymm3At11 = PDupwindNthSymmfdOrder43(&At11[index]);
        PDupwindNthAnti1At12 = PDupwindNthAntifdOrder41(&At12[index]);
        PDupwindNthSymm1At12 = PDupwindNthSymmfdOrder41(&At12[index]);
        PDupwindNthAnti2At12 = PDupwindNthAntifdOrder42(&At12[index]);
        PDupwindNthSymm2At12 = PDupwindNthSymmfdOrder42(&At12[index]);
        PDupwindNthAnti3At12 = PDupwindNthAntifdOrder43(&At12[index]);
        PDupwindNthSymm3At12 = PDupwindNthSymmfdOrder43(&At12[index]);
        PDupwindNthAnti1At13 = PDupwindNthAntifdOrder41(&At13[index]);
        PDupwindNthSymm1At13 = PDupwindNthSymmfdOrder41(&At13[index]);
        PDupwindNthAnti2At13 = PDupwindNthAntifdOrder42(&At13[index]);
        PDupwindNthSymm2At13 = PDupwindNthSymmfdOrder42(&At13[index]);
        PDupwindNthAnti3At13 = PDupwindNthAntifdOrder43(&At13[index]);
        PDupwindNthSymm3At13 = PDupwindNthSymmfdOrder43(&At13[index]);
        PDupwindNthAnti1At22 = PDupwindNthAntifdOrder41(&At22[index]);
        PDupwindNthSymm1At22 = PDupwindNthSymmfdOrder41(&At22[index]);
        PDupwindNthAnti2At22 = PDupwindNthAntifdOrder42(&At22[index]);
        PDupwindNthSymm2At22 = PDupwindNthSymmfdOrder42(&At22[index]);
        PDupwindNthAnti3At22 = PDupwindNthAntifdOrder43(&At22[index]);
        PDupwindNthSymm3At22 = PDupwindNthSymmfdOrder43(&At22[index]);
        PDupwindNthAnti1At23 = PDupwindNthAntifdOrder41(&At23[index]);
        PDupwindNthSymm1At23 = PDupwindNthSymmfdOrder41(&At23[index]);
        PDupwindNthAnti2At23 = PDupwindNthAntifdOrder42(&At23[index]);
        PDupwindNthSymm2At23 = PDupwindNthSymmfdOrder42(&At23[index]);
        PDupwindNthAnti3At23 = PDupwindNthAntifdOrder43(&At23[index]);
        PDupwindNthSymm3At23 = PDupwindNthSymmfdOrder43(&At23[index]);
        PDupwindNthAnti1At33 = PDupwindNthAntifdOrder41(&At33[index]);
        PDupwindNthSymm1At33 = PDupwindNthSymmfdOrder41(&At33[index]);
        PDupwindNthAnti2At33 = PDupwindNthAntifdOrder42(&At33[index]);
        PDupwindNthSymm2At33 = PDupwindNthSymmfdOrder42(&At33[index]);
        PDupwindNthAnti3At33 = PDupwindNthAntifdOrder43(&At33[index]);
        PDupwindNthSymm3At33 = PDupwindNthSymmfdOrder43(&At33[index]);
        PDupwindNthAnti1B1 = PDupwindNthAntifdOrder41(&B1[index]);
        PDupwindNthSymm1B1 = PDupwindNthSymmfdOrder41(&B1[index]);
        PDupwindNthAnti2B1 = PDupwindNthAntifdOrder42(&B1[index]);
        PDupwindNthSymm2B1 = PDupwindNthSymmfdOrder42(&B1[index]);
        PDupwindNthAnti3B1 = PDupwindNthAntifdOrder43(&B1[index]);
        PDupwindNthSymm3B1 = PDupwindNthSymmfdOrder43(&B1[index]);
        PDupwindNthAnti1B2 = PDupwindNthAntifdOrder41(&B2[index]);
        PDupwindNthSymm1B2 = PDupwindNthSymmfdOrder41(&B2[index]);
        PDupwindNthAnti2B2 = PDupwindNthAntifdOrder42(&B2[index]);
        PDupwindNthSymm2B2 = PDupwindNthSymmfdOrder42(&B2[index]);
        PDupwindNthAnti3B2 = PDupwindNthAntifdOrder43(&B2[index]);
        PDupwindNthSymm3B2 = PDupwindNthSymmfdOrder43(&B2[index]);
        PDupwindNthAnti1B3 = PDupwindNthAntifdOrder41(&B3[index]);
        PDupwindNthSymm1B3 = PDupwindNthSymmfdOrder41(&B3[index]);
        PDupwindNthAnti2B3 = PDupwindNthAntifdOrder42(&B3[index]);
        PDupwindNthSymm2B3 = PDupwindNthSymmfdOrder42(&B3[index]);
        PDupwindNthAnti3B3 = PDupwindNthAntifdOrder43(&B3[index]);
        PDupwindNthSymm3B3 = PDupwindNthSymmfdOrder43(&B3[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder41(&beta1[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder41(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder42(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder42(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder43(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder43(&beta1[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder41(&beta2[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder41(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder42(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder42(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder43(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder43(&beta2[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder41(&beta3[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder41(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder42(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder42(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder43(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder43(&beta3[index]);
        PDupwindNthAnti1gt11 = PDupwindNthAntifdOrder41(&gt11[index]);
        PDupwindNthSymm1gt11 = PDupwindNthSymmfdOrder41(&gt11[index]);
        PDupwindNthAnti2gt11 = PDupwindNthAntifdOrder42(&gt11[index]);
        PDupwindNthSymm2gt11 = PDupwindNthSymmfdOrder42(&gt11[index]);
        PDupwindNthAnti3gt11 = PDupwindNthAntifdOrder43(&gt11[index]);
        PDupwindNthSymm3gt11 = PDupwindNthSymmfdOrder43(&gt11[index]);
        PDupwindNthAnti1gt12 = PDupwindNthAntifdOrder41(&gt12[index]);
        PDupwindNthSymm1gt12 = PDupwindNthSymmfdOrder41(&gt12[index]);
        PDupwindNthAnti2gt12 = PDupwindNthAntifdOrder42(&gt12[index]);
        PDupwindNthSymm2gt12 = PDupwindNthSymmfdOrder42(&gt12[index]);
        PDupwindNthAnti3gt12 = PDupwindNthAntifdOrder43(&gt12[index]);
        PDupwindNthSymm3gt12 = PDupwindNthSymmfdOrder43(&gt12[index]);
        PDupwindNthAnti1gt13 = PDupwindNthAntifdOrder41(&gt13[index]);
        PDupwindNthSymm1gt13 = PDupwindNthSymmfdOrder41(&gt13[index]);
        PDupwindNthAnti2gt13 = PDupwindNthAntifdOrder42(&gt13[index]);
        PDupwindNthSymm2gt13 = PDupwindNthSymmfdOrder42(&gt13[index]);
        PDupwindNthAnti3gt13 = PDupwindNthAntifdOrder43(&gt13[index]);
        PDupwindNthSymm3gt13 = PDupwindNthSymmfdOrder43(&gt13[index]);
        PDupwindNthAnti1gt22 = PDupwindNthAntifdOrder41(&gt22[index]);
        PDupwindNthSymm1gt22 = PDupwindNthSymmfdOrder41(&gt22[index]);
        PDupwindNthAnti2gt22 = PDupwindNthAntifdOrder42(&gt22[index]);
        PDupwindNthSymm2gt22 = PDupwindNthSymmfdOrder42(&gt22[index]);
        PDupwindNthAnti3gt22 = PDupwindNthAntifdOrder43(&gt22[index]);
        PDupwindNthSymm3gt22 = PDupwindNthSymmfdOrder43(&gt22[index]);
        PDupwindNthAnti1gt23 = PDupwindNthAntifdOrder41(&gt23[index]);
        PDupwindNthSymm1gt23 = PDupwindNthSymmfdOrder41(&gt23[index]);
        PDupwindNthAnti2gt23 = PDupwindNthAntifdOrder42(&gt23[index]);
        PDupwindNthSymm2gt23 = PDupwindNthSymmfdOrder42(&gt23[index]);
        PDupwindNthAnti3gt23 = PDupwindNthAntifdOrder43(&gt23[index]);
        PDupwindNthSymm3gt23 = PDupwindNthSymmfdOrder43(&gt23[index]);
        PDupwindNthAnti1gt33 = PDupwindNthAntifdOrder41(&gt33[index]);
        PDupwindNthSymm1gt33 = PDupwindNthSymmfdOrder41(&gt33[index]);
        PDupwindNthAnti2gt33 = PDupwindNthAntifdOrder42(&gt33[index]);
        PDupwindNthSymm2gt33 = PDupwindNthSymmfdOrder42(&gt33[index]);
        PDupwindNthAnti3gt33 = PDupwindNthAntifdOrder43(&gt33[index]);
        PDupwindNthSymm3gt33 = PDupwindNthSymmfdOrder43(&gt33[index]);
        PDupwindNthAnti1phi = PDupwindNthAntifdOrder41(&phi[index]);
        PDupwindNthSymm1phi = PDupwindNthSymmfdOrder41(&phi[index]);
        PDupwindNthAnti2phi = PDupwindNthAntifdOrder42(&phi[index]);
        PDupwindNthSymm2phi = PDupwindNthSymmfdOrder42(&phi[index]);
        PDupwindNthAnti3phi = PDupwindNthAntifdOrder43(&phi[index]);
        PDupwindNthSymm3phi = PDupwindNthSymmfdOrder43(&phi[index]);
        PDupwindNthAnti1trK = PDupwindNthAntifdOrder41(&trK[index]);
        PDupwindNthSymm1trK = PDupwindNthSymmfdOrder41(&trK[index]);
        PDupwindNthAnti2trK = PDupwindNthAntifdOrder42(&trK[index]);
        PDupwindNthSymm2trK = PDupwindNthSymmfdOrder42(&trK[index]);
        PDupwindNthAnti3trK = PDupwindNthAntifdOrder43(&trK[index]);
        PDupwindNthSymm3trK = PDupwindNthSymmfdOrder43(&trK[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder41(&Xt1[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder41(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder42(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder42(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder43(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder43(&Xt1[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder41(&Xt2[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder41(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder42(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder42(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder43(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder43(&Xt2[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder41(&Xt3[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder41(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder42(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder42(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder43(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder43(&Xt3[index]);
        break;
      
      case 6:
        PDupwindNthAnti1A = PDupwindNthAntifdOrder61(&A[index]);
        PDupwindNthSymm1A = PDupwindNthSymmfdOrder61(&A[index]);
        PDupwindNthAnti2A = PDupwindNthAntifdOrder62(&A[index]);
        PDupwindNthSymm2A = PDupwindNthSymmfdOrder62(&A[index]);
        PDupwindNthAnti3A = PDupwindNthAntifdOrder63(&A[index]);
        PDupwindNthSymm3A = PDupwindNthSymmfdOrder63(&A[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder61(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder61(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder62(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder62(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder63(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder63(&alpha[index]);
        PDupwindNthAnti1At11 = PDupwindNthAntifdOrder61(&At11[index]);
        PDupwindNthSymm1At11 = PDupwindNthSymmfdOrder61(&At11[index]);
        PDupwindNthAnti2At11 = PDupwindNthAntifdOrder62(&At11[index]);
        PDupwindNthSymm2At11 = PDupwindNthSymmfdOrder62(&At11[index]);
        PDupwindNthAnti3At11 = PDupwindNthAntifdOrder63(&At11[index]);
        PDupwindNthSymm3At11 = PDupwindNthSymmfdOrder63(&At11[index]);
        PDupwindNthAnti1At12 = PDupwindNthAntifdOrder61(&At12[index]);
        PDupwindNthSymm1At12 = PDupwindNthSymmfdOrder61(&At12[index]);
        PDupwindNthAnti2At12 = PDupwindNthAntifdOrder62(&At12[index]);
        PDupwindNthSymm2At12 = PDupwindNthSymmfdOrder62(&At12[index]);
        PDupwindNthAnti3At12 = PDupwindNthAntifdOrder63(&At12[index]);
        PDupwindNthSymm3At12 = PDupwindNthSymmfdOrder63(&At12[index]);
        PDupwindNthAnti1At13 = PDupwindNthAntifdOrder61(&At13[index]);
        PDupwindNthSymm1At13 = PDupwindNthSymmfdOrder61(&At13[index]);
        PDupwindNthAnti2At13 = PDupwindNthAntifdOrder62(&At13[index]);
        PDupwindNthSymm2At13 = PDupwindNthSymmfdOrder62(&At13[index]);
        PDupwindNthAnti3At13 = PDupwindNthAntifdOrder63(&At13[index]);
        PDupwindNthSymm3At13 = PDupwindNthSymmfdOrder63(&At13[index]);
        PDupwindNthAnti1At22 = PDupwindNthAntifdOrder61(&At22[index]);
        PDupwindNthSymm1At22 = PDupwindNthSymmfdOrder61(&At22[index]);
        PDupwindNthAnti2At22 = PDupwindNthAntifdOrder62(&At22[index]);
        PDupwindNthSymm2At22 = PDupwindNthSymmfdOrder62(&At22[index]);
        PDupwindNthAnti3At22 = PDupwindNthAntifdOrder63(&At22[index]);
        PDupwindNthSymm3At22 = PDupwindNthSymmfdOrder63(&At22[index]);
        PDupwindNthAnti1At23 = PDupwindNthAntifdOrder61(&At23[index]);
        PDupwindNthSymm1At23 = PDupwindNthSymmfdOrder61(&At23[index]);
        PDupwindNthAnti2At23 = PDupwindNthAntifdOrder62(&At23[index]);
        PDupwindNthSymm2At23 = PDupwindNthSymmfdOrder62(&At23[index]);
        PDupwindNthAnti3At23 = PDupwindNthAntifdOrder63(&At23[index]);
        PDupwindNthSymm3At23 = PDupwindNthSymmfdOrder63(&At23[index]);
        PDupwindNthAnti1At33 = PDupwindNthAntifdOrder61(&At33[index]);
        PDupwindNthSymm1At33 = PDupwindNthSymmfdOrder61(&At33[index]);
        PDupwindNthAnti2At33 = PDupwindNthAntifdOrder62(&At33[index]);
        PDupwindNthSymm2At33 = PDupwindNthSymmfdOrder62(&At33[index]);
        PDupwindNthAnti3At33 = PDupwindNthAntifdOrder63(&At33[index]);
        PDupwindNthSymm3At33 = PDupwindNthSymmfdOrder63(&At33[index]);
        PDupwindNthAnti1B1 = PDupwindNthAntifdOrder61(&B1[index]);
        PDupwindNthSymm1B1 = PDupwindNthSymmfdOrder61(&B1[index]);
        PDupwindNthAnti2B1 = PDupwindNthAntifdOrder62(&B1[index]);
        PDupwindNthSymm2B1 = PDupwindNthSymmfdOrder62(&B1[index]);
        PDupwindNthAnti3B1 = PDupwindNthAntifdOrder63(&B1[index]);
        PDupwindNthSymm3B1 = PDupwindNthSymmfdOrder63(&B1[index]);
        PDupwindNthAnti1B2 = PDupwindNthAntifdOrder61(&B2[index]);
        PDupwindNthSymm1B2 = PDupwindNthSymmfdOrder61(&B2[index]);
        PDupwindNthAnti2B2 = PDupwindNthAntifdOrder62(&B2[index]);
        PDupwindNthSymm2B2 = PDupwindNthSymmfdOrder62(&B2[index]);
        PDupwindNthAnti3B2 = PDupwindNthAntifdOrder63(&B2[index]);
        PDupwindNthSymm3B2 = PDupwindNthSymmfdOrder63(&B2[index]);
        PDupwindNthAnti1B3 = PDupwindNthAntifdOrder61(&B3[index]);
        PDupwindNthSymm1B3 = PDupwindNthSymmfdOrder61(&B3[index]);
        PDupwindNthAnti2B3 = PDupwindNthAntifdOrder62(&B3[index]);
        PDupwindNthSymm2B3 = PDupwindNthSymmfdOrder62(&B3[index]);
        PDupwindNthAnti3B3 = PDupwindNthAntifdOrder63(&B3[index]);
        PDupwindNthSymm3B3 = PDupwindNthSymmfdOrder63(&B3[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder61(&beta1[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder61(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder62(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder62(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder63(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder63(&beta1[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder61(&beta2[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder61(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder62(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder62(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder63(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder63(&beta2[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder61(&beta3[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder61(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder62(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder62(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder63(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder63(&beta3[index]);
        PDupwindNthAnti1gt11 = PDupwindNthAntifdOrder61(&gt11[index]);
        PDupwindNthSymm1gt11 = PDupwindNthSymmfdOrder61(&gt11[index]);
        PDupwindNthAnti2gt11 = PDupwindNthAntifdOrder62(&gt11[index]);
        PDupwindNthSymm2gt11 = PDupwindNthSymmfdOrder62(&gt11[index]);
        PDupwindNthAnti3gt11 = PDupwindNthAntifdOrder63(&gt11[index]);
        PDupwindNthSymm3gt11 = PDupwindNthSymmfdOrder63(&gt11[index]);
        PDupwindNthAnti1gt12 = PDupwindNthAntifdOrder61(&gt12[index]);
        PDupwindNthSymm1gt12 = PDupwindNthSymmfdOrder61(&gt12[index]);
        PDupwindNthAnti2gt12 = PDupwindNthAntifdOrder62(&gt12[index]);
        PDupwindNthSymm2gt12 = PDupwindNthSymmfdOrder62(&gt12[index]);
        PDupwindNthAnti3gt12 = PDupwindNthAntifdOrder63(&gt12[index]);
        PDupwindNthSymm3gt12 = PDupwindNthSymmfdOrder63(&gt12[index]);
        PDupwindNthAnti1gt13 = PDupwindNthAntifdOrder61(&gt13[index]);
        PDupwindNthSymm1gt13 = PDupwindNthSymmfdOrder61(&gt13[index]);
        PDupwindNthAnti2gt13 = PDupwindNthAntifdOrder62(&gt13[index]);
        PDupwindNthSymm2gt13 = PDupwindNthSymmfdOrder62(&gt13[index]);
        PDupwindNthAnti3gt13 = PDupwindNthAntifdOrder63(&gt13[index]);
        PDupwindNthSymm3gt13 = PDupwindNthSymmfdOrder63(&gt13[index]);
        PDupwindNthAnti1gt22 = PDupwindNthAntifdOrder61(&gt22[index]);
        PDupwindNthSymm1gt22 = PDupwindNthSymmfdOrder61(&gt22[index]);
        PDupwindNthAnti2gt22 = PDupwindNthAntifdOrder62(&gt22[index]);
        PDupwindNthSymm2gt22 = PDupwindNthSymmfdOrder62(&gt22[index]);
        PDupwindNthAnti3gt22 = PDupwindNthAntifdOrder63(&gt22[index]);
        PDupwindNthSymm3gt22 = PDupwindNthSymmfdOrder63(&gt22[index]);
        PDupwindNthAnti1gt23 = PDupwindNthAntifdOrder61(&gt23[index]);
        PDupwindNthSymm1gt23 = PDupwindNthSymmfdOrder61(&gt23[index]);
        PDupwindNthAnti2gt23 = PDupwindNthAntifdOrder62(&gt23[index]);
        PDupwindNthSymm2gt23 = PDupwindNthSymmfdOrder62(&gt23[index]);
        PDupwindNthAnti3gt23 = PDupwindNthAntifdOrder63(&gt23[index]);
        PDupwindNthSymm3gt23 = PDupwindNthSymmfdOrder63(&gt23[index]);
        PDupwindNthAnti1gt33 = PDupwindNthAntifdOrder61(&gt33[index]);
        PDupwindNthSymm1gt33 = PDupwindNthSymmfdOrder61(&gt33[index]);
        PDupwindNthAnti2gt33 = PDupwindNthAntifdOrder62(&gt33[index]);
        PDupwindNthSymm2gt33 = PDupwindNthSymmfdOrder62(&gt33[index]);
        PDupwindNthAnti3gt33 = PDupwindNthAntifdOrder63(&gt33[index]);
        PDupwindNthSymm3gt33 = PDupwindNthSymmfdOrder63(&gt33[index]);
        PDupwindNthAnti1phi = PDupwindNthAntifdOrder61(&phi[index]);
        PDupwindNthSymm1phi = PDupwindNthSymmfdOrder61(&phi[index]);
        PDupwindNthAnti2phi = PDupwindNthAntifdOrder62(&phi[index]);
        PDupwindNthSymm2phi = PDupwindNthSymmfdOrder62(&phi[index]);
        PDupwindNthAnti3phi = PDupwindNthAntifdOrder63(&phi[index]);
        PDupwindNthSymm3phi = PDupwindNthSymmfdOrder63(&phi[index]);
        PDupwindNthAnti1trK = PDupwindNthAntifdOrder61(&trK[index]);
        PDupwindNthSymm1trK = PDupwindNthSymmfdOrder61(&trK[index]);
        PDupwindNthAnti2trK = PDupwindNthAntifdOrder62(&trK[index]);
        PDupwindNthSymm2trK = PDupwindNthSymmfdOrder62(&trK[index]);
        PDupwindNthAnti3trK = PDupwindNthAntifdOrder63(&trK[index]);
        PDupwindNthSymm3trK = PDupwindNthSymmfdOrder63(&trK[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder61(&Xt1[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder61(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder62(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder62(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder63(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder63(&Xt1[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder61(&Xt2[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder61(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder62(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder62(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder63(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder63(&Xt2[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder61(&Xt3[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder61(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder62(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder62(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder63(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder63(&Xt3[index]);
        break;
      
      case 8:
        PDupwindNthAnti1A = PDupwindNthAntifdOrder81(&A[index]);
        PDupwindNthSymm1A = PDupwindNthSymmfdOrder81(&A[index]);
        PDupwindNthAnti2A = PDupwindNthAntifdOrder82(&A[index]);
        PDupwindNthSymm2A = PDupwindNthSymmfdOrder82(&A[index]);
        PDupwindNthAnti3A = PDupwindNthAntifdOrder83(&A[index]);
        PDupwindNthSymm3A = PDupwindNthSymmfdOrder83(&A[index]);
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder81(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder81(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder82(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder82(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder83(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder83(&alpha[index]);
        PDupwindNthAnti1At11 = PDupwindNthAntifdOrder81(&At11[index]);
        PDupwindNthSymm1At11 = PDupwindNthSymmfdOrder81(&At11[index]);
        PDupwindNthAnti2At11 = PDupwindNthAntifdOrder82(&At11[index]);
        PDupwindNthSymm2At11 = PDupwindNthSymmfdOrder82(&At11[index]);
        PDupwindNthAnti3At11 = PDupwindNthAntifdOrder83(&At11[index]);
        PDupwindNthSymm3At11 = PDupwindNthSymmfdOrder83(&At11[index]);
        PDupwindNthAnti1At12 = PDupwindNthAntifdOrder81(&At12[index]);
        PDupwindNthSymm1At12 = PDupwindNthSymmfdOrder81(&At12[index]);
        PDupwindNthAnti2At12 = PDupwindNthAntifdOrder82(&At12[index]);
        PDupwindNthSymm2At12 = PDupwindNthSymmfdOrder82(&At12[index]);
        PDupwindNthAnti3At12 = PDupwindNthAntifdOrder83(&At12[index]);
        PDupwindNthSymm3At12 = PDupwindNthSymmfdOrder83(&At12[index]);
        PDupwindNthAnti1At13 = PDupwindNthAntifdOrder81(&At13[index]);
        PDupwindNthSymm1At13 = PDupwindNthSymmfdOrder81(&At13[index]);
        PDupwindNthAnti2At13 = PDupwindNthAntifdOrder82(&At13[index]);
        PDupwindNthSymm2At13 = PDupwindNthSymmfdOrder82(&At13[index]);
        PDupwindNthAnti3At13 = PDupwindNthAntifdOrder83(&At13[index]);
        PDupwindNthSymm3At13 = PDupwindNthSymmfdOrder83(&At13[index]);
        PDupwindNthAnti1At22 = PDupwindNthAntifdOrder81(&At22[index]);
        PDupwindNthSymm1At22 = PDupwindNthSymmfdOrder81(&At22[index]);
        PDupwindNthAnti2At22 = PDupwindNthAntifdOrder82(&At22[index]);
        PDupwindNthSymm2At22 = PDupwindNthSymmfdOrder82(&At22[index]);
        PDupwindNthAnti3At22 = PDupwindNthAntifdOrder83(&At22[index]);
        PDupwindNthSymm3At22 = PDupwindNthSymmfdOrder83(&At22[index]);
        PDupwindNthAnti1At23 = PDupwindNthAntifdOrder81(&At23[index]);
        PDupwindNthSymm1At23 = PDupwindNthSymmfdOrder81(&At23[index]);
        PDupwindNthAnti2At23 = PDupwindNthAntifdOrder82(&At23[index]);
        PDupwindNthSymm2At23 = PDupwindNthSymmfdOrder82(&At23[index]);
        PDupwindNthAnti3At23 = PDupwindNthAntifdOrder83(&At23[index]);
        PDupwindNthSymm3At23 = PDupwindNthSymmfdOrder83(&At23[index]);
        PDupwindNthAnti1At33 = PDupwindNthAntifdOrder81(&At33[index]);
        PDupwindNthSymm1At33 = PDupwindNthSymmfdOrder81(&At33[index]);
        PDupwindNthAnti2At33 = PDupwindNthAntifdOrder82(&At33[index]);
        PDupwindNthSymm2At33 = PDupwindNthSymmfdOrder82(&At33[index]);
        PDupwindNthAnti3At33 = PDupwindNthAntifdOrder83(&At33[index]);
        PDupwindNthSymm3At33 = PDupwindNthSymmfdOrder83(&At33[index]);
        PDupwindNthAnti1B1 = PDupwindNthAntifdOrder81(&B1[index]);
        PDupwindNthSymm1B1 = PDupwindNthSymmfdOrder81(&B1[index]);
        PDupwindNthAnti2B1 = PDupwindNthAntifdOrder82(&B1[index]);
        PDupwindNthSymm2B1 = PDupwindNthSymmfdOrder82(&B1[index]);
        PDupwindNthAnti3B1 = PDupwindNthAntifdOrder83(&B1[index]);
        PDupwindNthSymm3B1 = PDupwindNthSymmfdOrder83(&B1[index]);
        PDupwindNthAnti1B2 = PDupwindNthAntifdOrder81(&B2[index]);
        PDupwindNthSymm1B2 = PDupwindNthSymmfdOrder81(&B2[index]);
        PDupwindNthAnti2B2 = PDupwindNthAntifdOrder82(&B2[index]);
        PDupwindNthSymm2B2 = PDupwindNthSymmfdOrder82(&B2[index]);
        PDupwindNthAnti3B2 = PDupwindNthAntifdOrder83(&B2[index]);
        PDupwindNthSymm3B2 = PDupwindNthSymmfdOrder83(&B2[index]);
        PDupwindNthAnti1B3 = PDupwindNthAntifdOrder81(&B3[index]);
        PDupwindNthSymm1B3 = PDupwindNthSymmfdOrder81(&B3[index]);
        PDupwindNthAnti2B3 = PDupwindNthAntifdOrder82(&B3[index]);
        PDupwindNthSymm2B3 = PDupwindNthSymmfdOrder82(&B3[index]);
        PDupwindNthAnti3B3 = PDupwindNthAntifdOrder83(&B3[index]);
        PDupwindNthSymm3B3 = PDupwindNthSymmfdOrder83(&B3[index]);
        PDupwindNthAnti1beta1 = PDupwindNthAntifdOrder81(&beta1[index]);
        PDupwindNthSymm1beta1 = PDupwindNthSymmfdOrder81(&beta1[index]);
        PDupwindNthAnti2beta1 = PDupwindNthAntifdOrder82(&beta1[index]);
        PDupwindNthSymm2beta1 = PDupwindNthSymmfdOrder82(&beta1[index]);
        PDupwindNthAnti3beta1 = PDupwindNthAntifdOrder83(&beta1[index]);
        PDupwindNthSymm3beta1 = PDupwindNthSymmfdOrder83(&beta1[index]);
        PDupwindNthAnti1beta2 = PDupwindNthAntifdOrder81(&beta2[index]);
        PDupwindNthSymm1beta2 = PDupwindNthSymmfdOrder81(&beta2[index]);
        PDupwindNthAnti2beta2 = PDupwindNthAntifdOrder82(&beta2[index]);
        PDupwindNthSymm2beta2 = PDupwindNthSymmfdOrder82(&beta2[index]);
        PDupwindNthAnti3beta2 = PDupwindNthAntifdOrder83(&beta2[index]);
        PDupwindNthSymm3beta2 = PDupwindNthSymmfdOrder83(&beta2[index]);
        PDupwindNthAnti1beta3 = PDupwindNthAntifdOrder81(&beta3[index]);
        PDupwindNthSymm1beta3 = PDupwindNthSymmfdOrder81(&beta3[index]);
        PDupwindNthAnti2beta3 = PDupwindNthAntifdOrder82(&beta3[index]);
        PDupwindNthSymm2beta3 = PDupwindNthSymmfdOrder82(&beta3[index]);
        PDupwindNthAnti3beta3 = PDupwindNthAntifdOrder83(&beta3[index]);
        PDupwindNthSymm3beta3 = PDupwindNthSymmfdOrder83(&beta3[index]);
        PDupwindNthAnti1gt11 = PDupwindNthAntifdOrder81(&gt11[index]);
        PDupwindNthSymm1gt11 = PDupwindNthSymmfdOrder81(&gt11[index]);
        PDupwindNthAnti2gt11 = PDupwindNthAntifdOrder82(&gt11[index]);
        PDupwindNthSymm2gt11 = PDupwindNthSymmfdOrder82(&gt11[index]);
        PDupwindNthAnti3gt11 = PDupwindNthAntifdOrder83(&gt11[index]);
        PDupwindNthSymm3gt11 = PDupwindNthSymmfdOrder83(&gt11[index]);
        PDupwindNthAnti1gt12 = PDupwindNthAntifdOrder81(&gt12[index]);
        PDupwindNthSymm1gt12 = PDupwindNthSymmfdOrder81(&gt12[index]);
        PDupwindNthAnti2gt12 = PDupwindNthAntifdOrder82(&gt12[index]);
        PDupwindNthSymm2gt12 = PDupwindNthSymmfdOrder82(&gt12[index]);
        PDupwindNthAnti3gt12 = PDupwindNthAntifdOrder83(&gt12[index]);
        PDupwindNthSymm3gt12 = PDupwindNthSymmfdOrder83(&gt12[index]);
        PDupwindNthAnti1gt13 = PDupwindNthAntifdOrder81(&gt13[index]);
        PDupwindNthSymm1gt13 = PDupwindNthSymmfdOrder81(&gt13[index]);
        PDupwindNthAnti2gt13 = PDupwindNthAntifdOrder82(&gt13[index]);
        PDupwindNthSymm2gt13 = PDupwindNthSymmfdOrder82(&gt13[index]);
        PDupwindNthAnti3gt13 = PDupwindNthAntifdOrder83(&gt13[index]);
        PDupwindNthSymm3gt13 = PDupwindNthSymmfdOrder83(&gt13[index]);
        PDupwindNthAnti1gt22 = PDupwindNthAntifdOrder81(&gt22[index]);
        PDupwindNthSymm1gt22 = PDupwindNthSymmfdOrder81(&gt22[index]);
        PDupwindNthAnti2gt22 = PDupwindNthAntifdOrder82(&gt22[index]);
        PDupwindNthSymm2gt22 = PDupwindNthSymmfdOrder82(&gt22[index]);
        PDupwindNthAnti3gt22 = PDupwindNthAntifdOrder83(&gt22[index]);
        PDupwindNthSymm3gt22 = PDupwindNthSymmfdOrder83(&gt22[index]);
        PDupwindNthAnti1gt23 = PDupwindNthAntifdOrder81(&gt23[index]);
        PDupwindNthSymm1gt23 = PDupwindNthSymmfdOrder81(&gt23[index]);
        PDupwindNthAnti2gt23 = PDupwindNthAntifdOrder82(&gt23[index]);
        PDupwindNthSymm2gt23 = PDupwindNthSymmfdOrder82(&gt23[index]);
        PDupwindNthAnti3gt23 = PDupwindNthAntifdOrder83(&gt23[index]);
        PDupwindNthSymm3gt23 = PDupwindNthSymmfdOrder83(&gt23[index]);
        PDupwindNthAnti1gt33 = PDupwindNthAntifdOrder81(&gt33[index]);
        PDupwindNthSymm1gt33 = PDupwindNthSymmfdOrder81(&gt33[index]);
        PDupwindNthAnti2gt33 = PDupwindNthAntifdOrder82(&gt33[index]);
        PDupwindNthSymm2gt33 = PDupwindNthSymmfdOrder82(&gt33[index]);
        PDupwindNthAnti3gt33 = PDupwindNthAntifdOrder83(&gt33[index]);
        PDupwindNthSymm3gt33 = PDupwindNthSymmfdOrder83(&gt33[index]);
        PDupwindNthAnti1phi = PDupwindNthAntifdOrder81(&phi[index]);
        PDupwindNthSymm1phi = PDupwindNthSymmfdOrder81(&phi[index]);
        PDupwindNthAnti2phi = PDupwindNthAntifdOrder82(&phi[index]);
        PDupwindNthSymm2phi = PDupwindNthSymmfdOrder82(&phi[index]);
        PDupwindNthAnti3phi = PDupwindNthAntifdOrder83(&phi[index]);
        PDupwindNthSymm3phi = PDupwindNthSymmfdOrder83(&phi[index]);
        PDupwindNthAnti1trK = PDupwindNthAntifdOrder81(&trK[index]);
        PDupwindNthSymm1trK = PDupwindNthSymmfdOrder81(&trK[index]);
        PDupwindNthAnti2trK = PDupwindNthAntifdOrder82(&trK[index]);
        PDupwindNthSymm2trK = PDupwindNthSymmfdOrder82(&trK[index]);
        PDupwindNthAnti3trK = PDupwindNthAntifdOrder83(&trK[index]);
        PDupwindNthSymm3trK = PDupwindNthSymmfdOrder83(&trK[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder81(&Xt1[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder81(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder82(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder82(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder83(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder83(&Xt1[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder81(&Xt2[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder81(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder82(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder82(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder83(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder83(&Xt2[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder81(&Xt3[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder81(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder82(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder82(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder83(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder83(&Xt3[index]);
        break;
    }
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    CCTK_REAL_VEC JacPDupwindNthAnti1A;
    CCTK_REAL_VEC JacPDupwindNthAnti1alpha;
    CCTK_REAL_VEC JacPDupwindNthAnti1At11;
    CCTK_REAL_VEC JacPDupwindNthAnti1At12;
    CCTK_REAL_VEC JacPDupwindNthAnti1At13;
    CCTK_REAL_VEC JacPDupwindNthAnti1At22;
    CCTK_REAL_VEC JacPDupwindNthAnti1At23;
    CCTK_REAL_VEC JacPDupwindNthAnti1At33;
    CCTK_REAL_VEC JacPDupwindNthAnti1B1;
    CCTK_REAL_VEC JacPDupwindNthAnti1B2;
    CCTK_REAL_VEC JacPDupwindNthAnti1B3;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta1;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta2;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta3;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt11;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt12;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt13;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt22;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt23;
    CCTK_REAL_VEC JacPDupwindNthAnti1gt33;
    CCTK_REAL_VEC JacPDupwindNthAnti1phi;
    CCTK_REAL_VEC JacPDupwindNthAnti1trK;
    CCTK_REAL_VEC JacPDupwindNthAnti1Xt1;
    CCTK_REAL_VEC JacPDupwindNthAnti1Xt2;
    CCTK_REAL_VEC JacPDupwindNthAnti1Xt3;
    CCTK_REAL_VEC JacPDupwindNthAnti2A;
    CCTK_REAL_VEC JacPDupwindNthAnti2alpha;
    CCTK_REAL_VEC JacPDupwindNthAnti2At11;
    CCTK_REAL_VEC JacPDupwindNthAnti2At12;
    CCTK_REAL_VEC JacPDupwindNthAnti2At13;
    CCTK_REAL_VEC JacPDupwindNthAnti2At22;
    CCTK_REAL_VEC JacPDupwindNthAnti2At23;
    CCTK_REAL_VEC JacPDupwindNthAnti2At33;
    CCTK_REAL_VEC JacPDupwindNthAnti2B1;
    CCTK_REAL_VEC JacPDupwindNthAnti2B2;
    CCTK_REAL_VEC JacPDupwindNthAnti2B3;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta1;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta2;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta3;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt11;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt12;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt13;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt22;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt23;
    CCTK_REAL_VEC JacPDupwindNthAnti2gt33;
    CCTK_REAL_VEC JacPDupwindNthAnti2phi;
    CCTK_REAL_VEC JacPDupwindNthAnti2trK;
    CCTK_REAL_VEC JacPDupwindNthAnti2Xt1;
    CCTK_REAL_VEC JacPDupwindNthAnti2Xt2;
    CCTK_REAL_VEC JacPDupwindNthAnti2Xt3;
    CCTK_REAL_VEC JacPDupwindNthAnti3A;
    CCTK_REAL_VEC JacPDupwindNthAnti3alpha;
    CCTK_REAL_VEC JacPDupwindNthAnti3At11;
    CCTK_REAL_VEC JacPDupwindNthAnti3At12;
    CCTK_REAL_VEC JacPDupwindNthAnti3At13;
    CCTK_REAL_VEC JacPDupwindNthAnti3At22;
    CCTK_REAL_VEC JacPDupwindNthAnti3At23;
    CCTK_REAL_VEC JacPDupwindNthAnti3At33;
    CCTK_REAL_VEC JacPDupwindNthAnti3B1;
    CCTK_REAL_VEC JacPDupwindNthAnti3B2;
    CCTK_REAL_VEC JacPDupwindNthAnti3B3;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta1;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta2;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta3;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt11;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt12;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt13;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt22;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt23;
    CCTK_REAL_VEC JacPDupwindNthAnti3gt33;
    CCTK_REAL_VEC JacPDupwindNthAnti3phi;
    CCTK_REAL_VEC JacPDupwindNthAnti3trK;
    CCTK_REAL_VEC JacPDupwindNthAnti3Xt1;
    CCTK_REAL_VEC JacPDupwindNthAnti3Xt2;
    CCTK_REAL_VEC JacPDupwindNthAnti3Xt3;
    CCTK_REAL_VEC JacPDupwindNthSymm1A;
    CCTK_REAL_VEC JacPDupwindNthSymm1alpha;
    CCTK_REAL_VEC JacPDupwindNthSymm1At11;
    CCTK_REAL_VEC JacPDupwindNthSymm1At12;
    CCTK_REAL_VEC JacPDupwindNthSymm1At13;
    CCTK_REAL_VEC JacPDupwindNthSymm1At22;
    CCTK_REAL_VEC JacPDupwindNthSymm1At23;
    CCTK_REAL_VEC JacPDupwindNthSymm1At33;
    CCTK_REAL_VEC JacPDupwindNthSymm1B1;
    CCTK_REAL_VEC JacPDupwindNthSymm1B2;
    CCTK_REAL_VEC JacPDupwindNthSymm1B3;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta1;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta2;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta3;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt11;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt12;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt13;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt22;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt23;
    CCTK_REAL_VEC JacPDupwindNthSymm1gt33;
    CCTK_REAL_VEC JacPDupwindNthSymm1phi;
    CCTK_REAL_VEC JacPDupwindNthSymm1trK;
    CCTK_REAL_VEC JacPDupwindNthSymm1Xt1;
    CCTK_REAL_VEC JacPDupwindNthSymm1Xt2;
    CCTK_REAL_VEC JacPDupwindNthSymm1Xt3;
    CCTK_REAL_VEC JacPDupwindNthSymm2A;
    CCTK_REAL_VEC JacPDupwindNthSymm2alpha;
    CCTK_REAL_VEC JacPDupwindNthSymm2At11;
    CCTK_REAL_VEC JacPDupwindNthSymm2At12;
    CCTK_REAL_VEC JacPDupwindNthSymm2At13;
    CCTK_REAL_VEC JacPDupwindNthSymm2At22;
    CCTK_REAL_VEC JacPDupwindNthSymm2At23;
    CCTK_REAL_VEC JacPDupwindNthSymm2At33;
    CCTK_REAL_VEC JacPDupwindNthSymm2B1;
    CCTK_REAL_VEC JacPDupwindNthSymm2B2;
    CCTK_REAL_VEC JacPDupwindNthSymm2B3;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta1;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta2;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta3;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt11;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt12;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt13;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt22;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt23;
    CCTK_REAL_VEC JacPDupwindNthSymm2gt33;
    CCTK_REAL_VEC JacPDupwindNthSymm2phi;
    CCTK_REAL_VEC JacPDupwindNthSymm2trK;
    CCTK_REAL_VEC JacPDupwindNthSymm2Xt1;
    CCTK_REAL_VEC JacPDupwindNthSymm2Xt2;
    CCTK_REAL_VEC JacPDupwindNthSymm2Xt3;
    CCTK_REAL_VEC JacPDupwindNthSymm3A;
    CCTK_REAL_VEC JacPDupwindNthSymm3alpha;
    CCTK_REAL_VEC JacPDupwindNthSymm3At11;
    CCTK_REAL_VEC JacPDupwindNthSymm3At12;
    CCTK_REAL_VEC JacPDupwindNthSymm3At13;
    CCTK_REAL_VEC JacPDupwindNthSymm3At22;
    CCTK_REAL_VEC JacPDupwindNthSymm3At23;
    CCTK_REAL_VEC JacPDupwindNthSymm3At33;
    CCTK_REAL_VEC JacPDupwindNthSymm3B1;
    CCTK_REAL_VEC JacPDupwindNthSymm3B2;
    CCTK_REAL_VEC JacPDupwindNthSymm3B3;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta1;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta2;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta3;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt11;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt12;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt13;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt22;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt23;
    CCTK_REAL_VEC JacPDupwindNthSymm3gt33;
    CCTK_REAL_VEC JacPDupwindNthSymm3phi;
    CCTK_REAL_VEC JacPDupwindNthSymm3trK;
    CCTK_REAL_VEC JacPDupwindNthSymm3Xt1;
    CCTK_REAL_VEC JacPDupwindNthSymm3Xt2;
    CCTK_REAL_VEC JacPDupwindNthSymm3Xt3;
    
    if (use_jacobian)
    {
      JacPDupwindNthAnti1A = 
        kmadd(J11L,PDupwindNthAnti1A,kmadd(J21L,PDupwindNthAnti2A,kmul(J31L,PDupwindNthAnti3A)));
      
      JacPDupwindNthAnti1alpha = 
        kmadd(J11L,PDupwindNthAnti1alpha,kmadd(J21L,PDupwindNthAnti2alpha,kmul(J31L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti1At11 = 
        kmadd(J11L,PDupwindNthAnti1At11,kmadd(J21L,PDupwindNthAnti2At11,kmul(J31L,PDupwindNthAnti3At11)));
      
      JacPDupwindNthAnti1At12 = 
        kmadd(J11L,PDupwindNthAnti1At12,kmadd(J21L,PDupwindNthAnti2At12,kmul(J31L,PDupwindNthAnti3At12)));
      
      JacPDupwindNthAnti1At13 = 
        kmadd(J11L,PDupwindNthAnti1At13,kmadd(J21L,PDupwindNthAnti2At13,kmul(J31L,PDupwindNthAnti3At13)));
      
      JacPDupwindNthAnti1At22 = 
        kmadd(J11L,PDupwindNthAnti1At22,kmadd(J21L,PDupwindNthAnti2At22,kmul(J31L,PDupwindNthAnti3At22)));
      
      JacPDupwindNthAnti1At23 = 
        kmadd(J11L,PDupwindNthAnti1At23,kmadd(J21L,PDupwindNthAnti2At23,kmul(J31L,PDupwindNthAnti3At23)));
      
      JacPDupwindNthAnti1At33 = 
        kmadd(J11L,PDupwindNthAnti1At33,kmadd(J21L,PDupwindNthAnti2At33,kmul(J31L,PDupwindNthAnti3At33)));
      
      JacPDupwindNthAnti1B1 = 
        kmadd(J11L,PDupwindNthAnti1B1,kmadd(J21L,PDupwindNthAnti2B1,kmul(J31L,PDupwindNthAnti3B1)));
      
      JacPDupwindNthAnti1B2 = 
        kmadd(J11L,PDupwindNthAnti1B2,kmadd(J21L,PDupwindNthAnti2B2,kmul(J31L,PDupwindNthAnti3B2)));
      
      JacPDupwindNthAnti1B3 = 
        kmadd(J11L,PDupwindNthAnti1B3,kmadd(J21L,PDupwindNthAnti2B3,kmul(J31L,PDupwindNthAnti3B3)));
      
      JacPDupwindNthAnti1beta1 = 
        kmadd(J11L,PDupwindNthAnti1beta1,kmadd(J21L,PDupwindNthAnti2beta1,kmul(J31L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti1beta2 = 
        kmadd(J11L,PDupwindNthAnti1beta2,kmadd(J21L,PDupwindNthAnti2beta2,kmul(J31L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti1beta3 = 
        kmadd(J11L,PDupwindNthAnti1beta3,kmadd(J21L,PDupwindNthAnti2beta3,kmul(J31L,PDupwindNthAnti3beta3)));
      
      JacPDupwindNthAnti1gt11 = 
        kmadd(J11L,PDupwindNthAnti1gt11,kmadd(J21L,PDupwindNthAnti2gt11,kmul(J31L,PDupwindNthAnti3gt11)));
      
      JacPDupwindNthAnti1gt12 = 
        kmadd(J11L,PDupwindNthAnti1gt12,kmadd(J21L,PDupwindNthAnti2gt12,kmul(J31L,PDupwindNthAnti3gt12)));
      
      JacPDupwindNthAnti1gt13 = 
        kmadd(J11L,PDupwindNthAnti1gt13,kmadd(J21L,PDupwindNthAnti2gt13,kmul(J31L,PDupwindNthAnti3gt13)));
      
      JacPDupwindNthAnti1gt22 = 
        kmadd(J11L,PDupwindNthAnti1gt22,kmadd(J21L,PDupwindNthAnti2gt22,kmul(J31L,PDupwindNthAnti3gt22)));
      
      JacPDupwindNthAnti1gt23 = 
        kmadd(J11L,PDupwindNthAnti1gt23,kmadd(J21L,PDupwindNthAnti2gt23,kmul(J31L,PDupwindNthAnti3gt23)));
      
      JacPDupwindNthAnti1gt33 = 
        kmadd(J11L,PDupwindNthAnti1gt33,kmadd(J21L,PDupwindNthAnti2gt33,kmul(J31L,PDupwindNthAnti3gt33)));
      
      JacPDupwindNthAnti1phi = 
        kmadd(J11L,PDupwindNthAnti1phi,kmadd(J21L,PDupwindNthAnti2phi,kmul(J31L,PDupwindNthAnti3phi)));
      
      JacPDupwindNthAnti1trK = 
        kmadd(J11L,PDupwindNthAnti1trK,kmadd(J21L,PDupwindNthAnti2trK,kmul(J31L,PDupwindNthAnti3trK)));
      
      JacPDupwindNthAnti1Xt1 = 
        kmadd(J11L,PDupwindNthAnti1Xt1,kmadd(J21L,PDupwindNthAnti2Xt1,kmul(J31L,PDupwindNthAnti3Xt1)));
      
      JacPDupwindNthAnti1Xt2 = 
        kmadd(J11L,PDupwindNthAnti1Xt2,kmadd(J21L,PDupwindNthAnti2Xt2,kmul(J31L,PDupwindNthAnti3Xt2)));
      
      JacPDupwindNthAnti1Xt3 = 
        kmadd(J11L,PDupwindNthAnti1Xt3,kmadd(J21L,PDupwindNthAnti2Xt3,kmul(J31L,PDupwindNthAnti3Xt3)));
      
      JacPDupwindNthSymm1A = 
        kmadd(J11L,PDupwindNthSymm1A,kmadd(J21L,PDupwindNthSymm2A,kmul(J31L,PDupwindNthSymm3A)));
      
      JacPDupwindNthSymm1alpha = 
        kmadd(J11L,PDupwindNthSymm1alpha,kmadd(J21L,PDupwindNthSymm2alpha,kmul(J31L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm1At11 = 
        kmadd(J11L,PDupwindNthSymm1At11,kmadd(J21L,PDupwindNthSymm2At11,kmul(J31L,PDupwindNthSymm3At11)));
      
      JacPDupwindNthSymm1At12 = 
        kmadd(J11L,PDupwindNthSymm1At12,kmadd(J21L,PDupwindNthSymm2At12,kmul(J31L,PDupwindNthSymm3At12)));
      
      JacPDupwindNthSymm1At13 = 
        kmadd(J11L,PDupwindNthSymm1At13,kmadd(J21L,PDupwindNthSymm2At13,kmul(J31L,PDupwindNthSymm3At13)));
      
      JacPDupwindNthSymm1At22 = 
        kmadd(J11L,PDupwindNthSymm1At22,kmadd(J21L,PDupwindNthSymm2At22,kmul(J31L,PDupwindNthSymm3At22)));
      
      JacPDupwindNthSymm1At23 = 
        kmadd(J11L,PDupwindNthSymm1At23,kmadd(J21L,PDupwindNthSymm2At23,kmul(J31L,PDupwindNthSymm3At23)));
      
      JacPDupwindNthSymm1At33 = 
        kmadd(J11L,PDupwindNthSymm1At33,kmadd(J21L,PDupwindNthSymm2At33,kmul(J31L,PDupwindNthSymm3At33)));
      
      JacPDupwindNthSymm1B1 = 
        kmadd(J11L,PDupwindNthSymm1B1,kmadd(J21L,PDupwindNthSymm2B1,kmul(J31L,PDupwindNthSymm3B1)));
      
      JacPDupwindNthSymm1B2 = 
        kmadd(J11L,PDupwindNthSymm1B2,kmadd(J21L,PDupwindNthSymm2B2,kmul(J31L,PDupwindNthSymm3B2)));
      
      JacPDupwindNthSymm1B3 = 
        kmadd(J11L,PDupwindNthSymm1B3,kmadd(J21L,PDupwindNthSymm2B3,kmul(J31L,PDupwindNthSymm3B3)));
      
      JacPDupwindNthSymm1beta1 = 
        kmadd(J11L,PDupwindNthSymm1beta1,kmadd(J21L,PDupwindNthSymm2beta1,kmul(J31L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm1beta2 = 
        kmadd(J11L,PDupwindNthSymm1beta2,kmadd(J21L,PDupwindNthSymm2beta2,kmul(J31L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm1beta3 = 
        kmadd(J11L,PDupwindNthSymm1beta3,kmadd(J21L,PDupwindNthSymm2beta3,kmul(J31L,PDupwindNthSymm3beta3)));
      
      JacPDupwindNthSymm1gt11 = 
        kmadd(J11L,PDupwindNthSymm1gt11,kmadd(J21L,PDupwindNthSymm2gt11,kmul(J31L,PDupwindNthSymm3gt11)));
      
      JacPDupwindNthSymm1gt12 = 
        kmadd(J11L,PDupwindNthSymm1gt12,kmadd(J21L,PDupwindNthSymm2gt12,kmul(J31L,PDupwindNthSymm3gt12)));
      
      JacPDupwindNthSymm1gt13 = 
        kmadd(J11L,PDupwindNthSymm1gt13,kmadd(J21L,PDupwindNthSymm2gt13,kmul(J31L,PDupwindNthSymm3gt13)));
      
      JacPDupwindNthSymm1gt22 = 
        kmadd(J11L,PDupwindNthSymm1gt22,kmadd(J21L,PDupwindNthSymm2gt22,kmul(J31L,PDupwindNthSymm3gt22)));
      
      JacPDupwindNthSymm1gt23 = 
        kmadd(J11L,PDupwindNthSymm1gt23,kmadd(J21L,PDupwindNthSymm2gt23,kmul(J31L,PDupwindNthSymm3gt23)));
      
      JacPDupwindNthSymm1gt33 = 
        kmadd(J11L,PDupwindNthSymm1gt33,kmadd(J21L,PDupwindNthSymm2gt33,kmul(J31L,PDupwindNthSymm3gt33)));
      
      JacPDupwindNthSymm1phi = 
        kmadd(J11L,PDupwindNthSymm1phi,kmadd(J21L,PDupwindNthSymm2phi,kmul(J31L,PDupwindNthSymm3phi)));
      
      JacPDupwindNthSymm1trK = 
        kmadd(J11L,PDupwindNthSymm1trK,kmadd(J21L,PDupwindNthSymm2trK,kmul(J31L,PDupwindNthSymm3trK)));
      
      JacPDupwindNthSymm1Xt1 = 
        kmadd(J11L,PDupwindNthSymm1Xt1,kmadd(J21L,PDupwindNthSymm2Xt1,kmul(J31L,PDupwindNthSymm3Xt1)));
      
      JacPDupwindNthSymm1Xt2 = 
        kmadd(J11L,PDupwindNthSymm1Xt2,kmadd(J21L,PDupwindNthSymm2Xt2,kmul(J31L,PDupwindNthSymm3Xt2)));
      
      JacPDupwindNthSymm1Xt3 = 
        kmadd(J11L,PDupwindNthSymm1Xt3,kmadd(J21L,PDupwindNthSymm2Xt3,kmul(J31L,PDupwindNthSymm3Xt3)));
      
      JacPDupwindNthAnti2A = 
        kmadd(J12L,PDupwindNthAnti1A,kmadd(J22L,PDupwindNthAnti2A,kmul(J32L,PDupwindNthAnti3A)));
      
      JacPDupwindNthAnti2alpha = 
        kmadd(J12L,PDupwindNthAnti1alpha,kmadd(J22L,PDupwindNthAnti2alpha,kmul(J32L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti2At11 = 
        kmadd(J12L,PDupwindNthAnti1At11,kmadd(J22L,PDupwindNthAnti2At11,kmul(J32L,PDupwindNthAnti3At11)));
      
      JacPDupwindNthAnti2At12 = 
        kmadd(J12L,PDupwindNthAnti1At12,kmadd(J22L,PDupwindNthAnti2At12,kmul(J32L,PDupwindNthAnti3At12)));
      
      JacPDupwindNthAnti2At13 = 
        kmadd(J12L,PDupwindNthAnti1At13,kmadd(J22L,PDupwindNthAnti2At13,kmul(J32L,PDupwindNthAnti3At13)));
      
      JacPDupwindNthAnti2At22 = 
        kmadd(J12L,PDupwindNthAnti1At22,kmadd(J22L,PDupwindNthAnti2At22,kmul(J32L,PDupwindNthAnti3At22)));
      
      JacPDupwindNthAnti2At23 = 
        kmadd(J12L,PDupwindNthAnti1At23,kmadd(J22L,PDupwindNthAnti2At23,kmul(J32L,PDupwindNthAnti3At23)));
      
      JacPDupwindNthAnti2At33 = 
        kmadd(J12L,PDupwindNthAnti1At33,kmadd(J22L,PDupwindNthAnti2At33,kmul(J32L,PDupwindNthAnti3At33)));
      
      JacPDupwindNthAnti2B1 = 
        kmadd(J12L,PDupwindNthAnti1B1,kmadd(J22L,PDupwindNthAnti2B1,kmul(J32L,PDupwindNthAnti3B1)));
      
      JacPDupwindNthAnti2B2 = 
        kmadd(J12L,PDupwindNthAnti1B2,kmadd(J22L,PDupwindNthAnti2B2,kmul(J32L,PDupwindNthAnti3B2)));
      
      JacPDupwindNthAnti2B3 = 
        kmadd(J12L,PDupwindNthAnti1B3,kmadd(J22L,PDupwindNthAnti2B3,kmul(J32L,PDupwindNthAnti3B3)));
      
      JacPDupwindNthAnti2beta1 = 
        kmadd(J12L,PDupwindNthAnti1beta1,kmadd(J22L,PDupwindNthAnti2beta1,kmul(J32L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti2beta2 = 
        kmadd(J12L,PDupwindNthAnti1beta2,kmadd(J22L,PDupwindNthAnti2beta2,kmul(J32L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti2beta3 = 
        kmadd(J12L,PDupwindNthAnti1beta3,kmadd(J22L,PDupwindNthAnti2beta3,kmul(J32L,PDupwindNthAnti3beta3)));
      
      JacPDupwindNthAnti2gt11 = 
        kmadd(J12L,PDupwindNthAnti1gt11,kmadd(J22L,PDupwindNthAnti2gt11,kmul(J32L,PDupwindNthAnti3gt11)));
      
      JacPDupwindNthAnti2gt12 = 
        kmadd(J12L,PDupwindNthAnti1gt12,kmadd(J22L,PDupwindNthAnti2gt12,kmul(J32L,PDupwindNthAnti3gt12)));
      
      JacPDupwindNthAnti2gt13 = 
        kmadd(J12L,PDupwindNthAnti1gt13,kmadd(J22L,PDupwindNthAnti2gt13,kmul(J32L,PDupwindNthAnti3gt13)));
      
      JacPDupwindNthAnti2gt22 = 
        kmadd(J12L,PDupwindNthAnti1gt22,kmadd(J22L,PDupwindNthAnti2gt22,kmul(J32L,PDupwindNthAnti3gt22)));
      
      JacPDupwindNthAnti2gt23 = 
        kmadd(J12L,PDupwindNthAnti1gt23,kmadd(J22L,PDupwindNthAnti2gt23,kmul(J32L,PDupwindNthAnti3gt23)));
      
      JacPDupwindNthAnti2gt33 = 
        kmadd(J12L,PDupwindNthAnti1gt33,kmadd(J22L,PDupwindNthAnti2gt33,kmul(J32L,PDupwindNthAnti3gt33)));
      
      JacPDupwindNthAnti2phi = 
        kmadd(J12L,PDupwindNthAnti1phi,kmadd(J22L,PDupwindNthAnti2phi,kmul(J32L,PDupwindNthAnti3phi)));
      
      JacPDupwindNthAnti2trK = 
        kmadd(J12L,PDupwindNthAnti1trK,kmadd(J22L,PDupwindNthAnti2trK,kmul(J32L,PDupwindNthAnti3trK)));
      
      JacPDupwindNthAnti2Xt1 = 
        kmadd(J12L,PDupwindNthAnti1Xt1,kmadd(J22L,PDupwindNthAnti2Xt1,kmul(J32L,PDupwindNthAnti3Xt1)));
      
      JacPDupwindNthAnti2Xt2 = 
        kmadd(J12L,PDupwindNthAnti1Xt2,kmadd(J22L,PDupwindNthAnti2Xt2,kmul(J32L,PDupwindNthAnti3Xt2)));
      
      JacPDupwindNthAnti2Xt3 = 
        kmadd(J12L,PDupwindNthAnti1Xt3,kmadd(J22L,PDupwindNthAnti2Xt3,kmul(J32L,PDupwindNthAnti3Xt3)));
      
      JacPDupwindNthSymm2A = 
        kmadd(J12L,PDupwindNthSymm1A,kmadd(J22L,PDupwindNthSymm2A,kmul(J32L,PDupwindNthSymm3A)));
      
      JacPDupwindNthSymm2alpha = 
        kmadd(J12L,PDupwindNthSymm1alpha,kmadd(J22L,PDupwindNthSymm2alpha,kmul(J32L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm2At11 = 
        kmadd(J12L,PDupwindNthSymm1At11,kmadd(J22L,PDupwindNthSymm2At11,kmul(J32L,PDupwindNthSymm3At11)));
      
      JacPDupwindNthSymm2At12 = 
        kmadd(J12L,PDupwindNthSymm1At12,kmadd(J22L,PDupwindNthSymm2At12,kmul(J32L,PDupwindNthSymm3At12)));
      
      JacPDupwindNthSymm2At13 = 
        kmadd(J12L,PDupwindNthSymm1At13,kmadd(J22L,PDupwindNthSymm2At13,kmul(J32L,PDupwindNthSymm3At13)));
      
      JacPDupwindNthSymm2At22 = 
        kmadd(J12L,PDupwindNthSymm1At22,kmadd(J22L,PDupwindNthSymm2At22,kmul(J32L,PDupwindNthSymm3At22)));
      
      JacPDupwindNthSymm2At23 = 
        kmadd(J12L,PDupwindNthSymm1At23,kmadd(J22L,PDupwindNthSymm2At23,kmul(J32L,PDupwindNthSymm3At23)));
      
      JacPDupwindNthSymm2At33 = 
        kmadd(J12L,PDupwindNthSymm1At33,kmadd(J22L,PDupwindNthSymm2At33,kmul(J32L,PDupwindNthSymm3At33)));
      
      JacPDupwindNthSymm2B1 = 
        kmadd(J12L,PDupwindNthSymm1B1,kmadd(J22L,PDupwindNthSymm2B1,kmul(J32L,PDupwindNthSymm3B1)));
      
      JacPDupwindNthSymm2B2 = 
        kmadd(J12L,PDupwindNthSymm1B2,kmadd(J22L,PDupwindNthSymm2B2,kmul(J32L,PDupwindNthSymm3B2)));
      
      JacPDupwindNthSymm2B3 = 
        kmadd(J12L,PDupwindNthSymm1B3,kmadd(J22L,PDupwindNthSymm2B3,kmul(J32L,PDupwindNthSymm3B3)));
      
      JacPDupwindNthSymm2beta1 = 
        kmadd(J12L,PDupwindNthSymm1beta1,kmadd(J22L,PDupwindNthSymm2beta1,kmul(J32L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm2beta2 = 
        kmadd(J12L,PDupwindNthSymm1beta2,kmadd(J22L,PDupwindNthSymm2beta2,kmul(J32L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm2beta3 = 
        kmadd(J12L,PDupwindNthSymm1beta3,kmadd(J22L,PDupwindNthSymm2beta3,kmul(J32L,PDupwindNthSymm3beta3)));
      
      JacPDupwindNthSymm2gt11 = 
        kmadd(J12L,PDupwindNthSymm1gt11,kmadd(J22L,PDupwindNthSymm2gt11,kmul(J32L,PDupwindNthSymm3gt11)));
      
      JacPDupwindNthSymm2gt12 = 
        kmadd(J12L,PDupwindNthSymm1gt12,kmadd(J22L,PDupwindNthSymm2gt12,kmul(J32L,PDupwindNthSymm3gt12)));
      
      JacPDupwindNthSymm2gt13 = 
        kmadd(J12L,PDupwindNthSymm1gt13,kmadd(J22L,PDupwindNthSymm2gt13,kmul(J32L,PDupwindNthSymm3gt13)));
      
      JacPDupwindNthSymm2gt22 = 
        kmadd(J12L,PDupwindNthSymm1gt22,kmadd(J22L,PDupwindNthSymm2gt22,kmul(J32L,PDupwindNthSymm3gt22)));
      
      JacPDupwindNthSymm2gt23 = 
        kmadd(J12L,PDupwindNthSymm1gt23,kmadd(J22L,PDupwindNthSymm2gt23,kmul(J32L,PDupwindNthSymm3gt23)));
      
      JacPDupwindNthSymm2gt33 = 
        kmadd(J12L,PDupwindNthSymm1gt33,kmadd(J22L,PDupwindNthSymm2gt33,kmul(J32L,PDupwindNthSymm3gt33)));
      
      JacPDupwindNthSymm2phi = 
        kmadd(J12L,PDupwindNthSymm1phi,kmadd(J22L,PDupwindNthSymm2phi,kmul(J32L,PDupwindNthSymm3phi)));
      
      JacPDupwindNthSymm2trK = 
        kmadd(J12L,PDupwindNthSymm1trK,kmadd(J22L,PDupwindNthSymm2trK,kmul(J32L,PDupwindNthSymm3trK)));
      
      JacPDupwindNthSymm2Xt1 = 
        kmadd(J12L,PDupwindNthSymm1Xt1,kmadd(J22L,PDupwindNthSymm2Xt1,kmul(J32L,PDupwindNthSymm3Xt1)));
      
      JacPDupwindNthSymm2Xt2 = 
        kmadd(J12L,PDupwindNthSymm1Xt2,kmadd(J22L,PDupwindNthSymm2Xt2,kmul(J32L,PDupwindNthSymm3Xt2)));
      
      JacPDupwindNthSymm2Xt3 = 
        kmadd(J12L,PDupwindNthSymm1Xt3,kmadd(J22L,PDupwindNthSymm2Xt3,kmul(J32L,PDupwindNthSymm3Xt3)));
      
      JacPDupwindNthAnti3A = 
        kmadd(J13L,PDupwindNthAnti1A,kmadd(J23L,PDupwindNthAnti2A,kmul(J33L,PDupwindNthAnti3A)));
      
      JacPDupwindNthAnti3alpha = 
        kmadd(J13L,PDupwindNthAnti1alpha,kmadd(J23L,PDupwindNthAnti2alpha,kmul(J33L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti3At11 = 
        kmadd(J13L,PDupwindNthAnti1At11,kmadd(J23L,PDupwindNthAnti2At11,kmul(J33L,PDupwindNthAnti3At11)));
      
      JacPDupwindNthAnti3At12 = 
        kmadd(J13L,PDupwindNthAnti1At12,kmadd(J23L,PDupwindNthAnti2At12,kmul(J33L,PDupwindNthAnti3At12)));
      
      JacPDupwindNthAnti3At13 = 
        kmadd(J13L,PDupwindNthAnti1At13,kmadd(J23L,PDupwindNthAnti2At13,kmul(J33L,PDupwindNthAnti3At13)));
      
      JacPDupwindNthAnti3At22 = 
        kmadd(J13L,PDupwindNthAnti1At22,kmadd(J23L,PDupwindNthAnti2At22,kmul(J33L,PDupwindNthAnti3At22)));
      
      JacPDupwindNthAnti3At23 = 
        kmadd(J13L,PDupwindNthAnti1At23,kmadd(J23L,PDupwindNthAnti2At23,kmul(J33L,PDupwindNthAnti3At23)));
      
      JacPDupwindNthAnti3At33 = 
        kmadd(J13L,PDupwindNthAnti1At33,kmadd(J23L,PDupwindNthAnti2At33,kmul(J33L,PDupwindNthAnti3At33)));
      
      JacPDupwindNthAnti3B1 = 
        kmadd(J13L,PDupwindNthAnti1B1,kmadd(J23L,PDupwindNthAnti2B1,kmul(J33L,PDupwindNthAnti3B1)));
      
      JacPDupwindNthAnti3B2 = 
        kmadd(J13L,PDupwindNthAnti1B2,kmadd(J23L,PDupwindNthAnti2B2,kmul(J33L,PDupwindNthAnti3B2)));
      
      JacPDupwindNthAnti3B3 = 
        kmadd(J13L,PDupwindNthAnti1B3,kmadd(J23L,PDupwindNthAnti2B3,kmul(J33L,PDupwindNthAnti3B3)));
      
      JacPDupwindNthAnti3beta1 = 
        kmadd(J13L,PDupwindNthAnti1beta1,kmadd(J23L,PDupwindNthAnti2beta1,kmul(J33L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti3beta2 = 
        kmadd(J13L,PDupwindNthAnti1beta2,kmadd(J23L,PDupwindNthAnti2beta2,kmul(J33L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti3beta3 = 
        kmadd(J13L,PDupwindNthAnti1beta3,kmadd(J23L,PDupwindNthAnti2beta3,kmul(J33L,PDupwindNthAnti3beta3)));
      
      JacPDupwindNthAnti3gt11 = 
        kmadd(J13L,PDupwindNthAnti1gt11,kmadd(J23L,PDupwindNthAnti2gt11,kmul(J33L,PDupwindNthAnti3gt11)));
      
      JacPDupwindNthAnti3gt12 = 
        kmadd(J13L,PDupwindNthAnti1gt12,kmadd(J23L,PDupwindNthAnti2gt12,kmul(J33L,PDupwindNthAnti3gt12)));
      
      JacPDupwindNthAnti3gt13 = 
        kmadd(J13L,PDupwindNthAnti1gt13,kmadd(J23L,PDupwindNthAnti2gt13,kmul(J33L,PDupwindNthAnti3gt13)));
      
      JacPDupwindNthAnti3gt22 = 
        kmadd(J13L,PDupwindNthAnti1gt22,kmadd(J23L,PDupwindNthAnti2gt22,kmul(J33L,PDupwindNthAnti3gt22)));
      
      JacPDupwindNthAnti3gt23 = 
        kmadd(J13L,PDupwindNthAnti1gt23,kmadd(J23L,PDupwindNthAnti2gt23,kmul(J33L,PDupwindNthAnti3gt23)));
      
      JacPDupwindNthAnti3gt33 = 
        kmadd(J13L,PDupwindNthAnti1gt33,kmadd(J23L,PDupwindNthAnti2gt33,kmul(J33L,PDupwindNthAnti3gt33)));
      
      JacPDupwindNthAnti3phi = 
        kmadd(J13L,PDupwindNthAnti1phi,kmadd(J23L,PDupwindNthAnti2phi,kmul(J33L,PDupwindNthAnti3phi)));
      
      JacPDupwindNthAnti3trK = 
        kmadd(J13L,PDupwindNthAnti1trK,kmadd(J23L,PDupwindNthAnti2trK,kmul(J33L,PDupwindNthAnti3trK)));
      
      JacPDupwindNthAnti3Xt1 = 
        kmadd(J13L,PDupwindNthAnti1Xt1,kmadd(J23L,PDupwindNthAnti2Xt1,kmul(J33L,PDupwindNthAnti3Xt1)));
      
      JacPDupwindNthAnti3Xt2 = 
        kmadd(J13L,PDupwindNthAnti1Xt2,kmadd(J23L,PDupwindNthAnti2Xt2,kmul(J33L,PDupwindNthAnti3Xt2)));
      
      JacPDupwindNthAnti3Xt3 = 
        kmadd(J13L,PDupwindNthAnti1Xt3,kmadd(J23L,PDupwindNthAnti2Xt3,kmul(J33L,PDupwindNthAnti3Xt3)));
      
      JacPDupwindNthSymm3A = 
        kmadd(J13L,PDupwindNthSymm1A,kmadd(J23L,PDupwindNthSymm2A,kmul(J33L,PDupwindNthSymm3A)));
      
      JacPDupwindNthSymm3alpha = 
        kmadd(J13L,PDupwindNthSymm1alpha,kmadd(J23L,PDupwindNthSymm2alpha,kmul(J33L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm3At11 = 
        kmadd(J13L,PDupwindNthSymm1At11,kmadd(J23L,PDupwindNthSymm2At11,kmul(J33L,PDupwindNthSymm3At11)));
      
      JacPDupwindNthSymm3At12 = 
        kmadd(J13L,PDupwindNthSymm1At12,kmadd(J23L,PDupwindNthSymm2At12,kmul(J33L,PDupwindNthSymm3At12)));
      
      JacPDupwindNthSymm3At13 = 
        kmadd(J13L,PDupwindNthSymm1At13,kmadd(J23L,PDupwindNthSymm2At13,kmul(J33L,PDupwindNthSymm3At13)));
      
      JacPDupwindNthSymm3At22 = 
        kmadd(J13L,PDupwindNthSymm1At22,kmadd(J23L,PDupwindNthSymm2At22,kmul(J33L,PDupwindNthSymm3At22)));
      
      JacPDupwindNthSymm3At23 = 
        kmadd(J13L,PDupwindNthSymm1At23,kmadd(J23L,PDupwindNthSymm2At23,kmul(J33L,PDupwindNthSymm3At23)));
      
      JacPDupwindNthSymm3At33 = 
        kmadd(J13L,PDupwindNthSymm1At33,kmadd(J23L,PDupwindNthSymm2At33,kmul(J33L,PDupwindNthSymm3At33)));
      
      JacPDupwindNthSymm3B1 = 
        kmadd(J13L,PDupwindNthSymm1B1,kmadd(J23L,PDupwindNthSymm2B1,kmul(J33L,PDupwindNthSymm3B1)));
      
      JacPDupwindNthSymm3B2 = 
        kmadd(J13L,PDupwindNthSymm1B2,kmadd(J23L,PDupwindNthSymm2B2,kmul(J33L,PDupwindNthSymm3B2)));
      
      JacPDupwindNthSymm3B3 = 
        kmadd(J13L,PDupwindNthSymm1B3,kmadd(J23L,PDupwindNthSymm2B3,kmul(J33L,PDupwindNthSymm3B3)));
      
      JacPDupwindNthSymm3beta1 = 
        kmadd(J13L,PDupwindNthSymm1beta1,kmadd(J23L,PDupwindNthSymm2beta1,kmul(J33L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm3beta2 = 
        kmadd(J13L,PDupwindNthSymm1beta2,kmadd(J23L,PDupwindNthSymm2beta2,kmul(J33L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm3beta3 = 
        kmadd(J13L,PDupwindNthSymm1beta3,kmadd(J23L,PDupwindNthSymm2beta3,kmul(J33L,PDupwindNthSymm3beta3)));
      
      JacPDupwindNthSymm3gt11 = 
        kmadd(J13L,PDupwindNthSymm1gt11,kmadd(J23L,PDupwindNthSymm2gt11,kmul(J33L,PDupwindNthSymm3gt11)));
      
      JacPDupwindNthSymm3gt12 = 
        kmadd(J13L,PDupwindNthSymm1gt12,kmadd(J23L,PDupwindNthSymm2gt12,kmul(J33L,PDupwindNthSymm3gt12)));
      
      JacPDupwindNthSymm3gt13 = 
        kmadd(J13L,PDupwindNthSymm1gt13,kmadd(J23L,PDupwindNthSymm2gt13,kmul(J33L,PDupwindNthSymm3gt13)));
      
      JacPDupwindNthSymm3gt22 = 
        kmadd(J13L,PDupwindNthSymm1gt22,kmadd(J23L,PDupwindNthSymm2gt22,kmul(J33L,PDupwindNthSymm3gt22)));
      
      JacPDupwindNthSymm3gt23 = 
        kmadd(J13L,PDupwindNthSymm1gt23,kmadd(J23L,PDupwindNthSymm2gt23,kmul(J33L,PDupwindNthSymm3gt23)));
      
      JacPDupwindNthSymm3gt33 = 
        kmadd(J13L,PDupwindNthSymm1gt33,kmadd(J23L,PDupwindNthSymm2gt33,kmul(J33L,PDupwindNthSymm3gt33)));
      
      JacPDupwindNthSymm3phi = 
        kmadd(J13L,PDupwindNthSymm1phi,kmadd(J23L,PDupwindNthSymm2phi,kmul(J33L,PDupwindNthSymm3phi)));
      
      JacPDupwindNthSymm3trK = 
        kmadd(J13L,PDupwindNthSymm1trK,kmadd(J23L,PDupwindNthSymm2trK,kmul(J33L,PDupwindNthSymm3trK)));
      
      JacPDupwindNthSymm3Xt1 = 
        kmadd(J13L,PDupwindNthSymm1Xt1,kmadd(J23L,PDupwindNthSymm2Xt1,kmul(J33L,PDupwindNthSymm3Xt1)));
      
      JacPDupwindNthSymm3Xt2 = 
        kmadd(J13L,PDupwindNthSymm1Xt2,kmadd(J23L,PDupwindNthSymm2Xt2,kmul(J33L,PDupwindNthSymm3Xt2)));
      
      JacPDupwindNthSymm3Xt3 = 
        kmadd(J13L,PDupwindNthSymm1Xt3,kmadd(J23L,PDupwindNthSymm2Xt3,kmul(J33L,PDupwindNthSymm3Xt3)));
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
    
    phirhsL = 
      kadd(phirhsL,kmadd(beta1L,JacPDupwindNthAnti1phi,kmadd(beta2L,JacPDupwindNthAnti2phi,kmadd(beta3L,JacPDupwindNthAnti3phi,kmadd(JacPDupwindNthSymm1phi,kfabs(beta1L),kmadd(JacPDupwindNthSymm2phi,kfabs(beta2L),kmul(JacPDupwindNthSymm3phi,kfabs(beta3L))))))));
    
    gt11rhsL = 
      kadd(gt11rhsL,kmadd(beta1L,JacPDupwindNthAnti1gt11,kmadd(beta2L,JacPDupwindNthAnti2gt11,kmadd(beta3L,JacPDupwindNthAnti3gt11,kmadd(JacPDupwindNthSymm1gt11,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt11,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt11,kfabs(beta3L))))))));
    
    gt12rhsL = 
      kadd(gt12rhsL,kmadd(beta1L,JacPDupwindNthAnti1gt12,kmadd(beta2L,JacPDupwindNthAnti2gt12,kmadd(beta3L,JacPDupwindNthAnti3gt12,kmadd(JacPDupwindNthSymm1gt12,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt12,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt12,kfabs(beta3L))))))));
    
    gt13rhsL = 
      kadd(gt13rhsL,kmadd(beta1L,JacPDupwindNthAnti1gt13,kmadd(beta2L,JacPDupwindNthAnti2gt13,kmadd(beta3L,JacPDupwindNthAnti3gt13,kmadd(JacPDupwindNthSymm1gt13,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt13,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt13,kfabs(beta3L))))))));
    
    gt22rhsL = 
      kadd(gt22rhsL,kmadd(beta1L,JacPDupwindNthAnti1gt22,kmadd(beta2L,JacPDupwindNthAnti2gt22,kmadd(beta3L,JacPDupwindNthAnti3gt22,kmadd(JacPDupwindNthSymm1gt22,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt22,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt22,kfabs(beta3L))))))));
    
    gt23rhsL = 
      kadd(gt23rhsL,kmadd(beta1L,JacPDupwindNthAnti1gt23,kmadd(beta2L,JacPDupwindNthAnti2gt23,kmadd(beta3L,JacPDupwindNthAnti3gt23,kmadd(JacPDupwindNthSymm1gt23,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt23,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt23,kfabs(beta3L))))))));
    
    gt33rhsL = 
      kadd(gt33rhsL,kmadd(beta1L,JacPDupwindNthAnti1gt33,kmadd(beta2L,JacPDupwindNthAnti2gt33,kmadd(beta3L,JacPDupwindNthAnti3gt33,kmadd(JacPDupwindNthSymm1gt33,kfabs(beta1L),kmadd(JacPDupwindNthSymm2gt33,kfabs(beta2L),kmul(JacPDupwindNthSymm3gt33,kfabs(beta3L))))))));
    
    Xt1rhsL = 
      kadd(Xt1rhsL,kmadd(beta1L,JacPDupwindNthAnti1Xt1,kmadd(beta2L,JacPDupwindNthAnti2Xt1,kmadd(beta3L,JacPDupwindNthAnti3Xt1,kmadd(JacPDupwindNthSymm1Xt1,kfabs(beta1L),kmadd(JacPDupwindNthSymm2Xt1,kfabs(beta2L),kmul(JacPDupwindNthSymm3Xt1,kfabs(beta3L))))))));
    
    Xt2rhsL = 
      kadd(Xt2rhsL,kmadd(beta1L,JacPDupwindNthAnti1Xt2,kmadd(beta2L,JacPDupwindNthAnti2Xt2,kmadd(beta3L,JacPDupwindNthAnti3Xt2,kmadd(JacPDupwindNthSymm1Xt2,kfabs(beta1L),kmadd(JacPDupwindNthSymm2Xt2,kfabs(beta2L),kmul(JacPDupwindNthSymm3Xt2,kfabs(beta3L))))))));
    
    Xt3rhsL = 
      kadd(Xt3rhsL,kmadd(beta1L,JacPDupwindNthAnti1Xt3,kmadd(beta2L,JacPDupwindNthAnti2Xt3,kmadd(beta3L,JacPDupwindNthAnti3Xt3,kmadd(JacPDupwindNthSymm1Xt3,kfabs(beta1L),kmadd(JacPDupwindNthSymm2Xt3,kfabs(beta2L),kmul(JacPDupwindNthSymm3Xt3,kfabs(beta3L))))))));
    
    trKrhsL = 
      kadd(trKrhsL,kmadd(beta1L,JacPDupwindNthAnti1trK,kmadd(beta2L,JacPDupwindNthAnti2trK,kmadd(beta3L,JacPDupwindNthAnti3trK,kmadd(JacPDupwindNthSymm1trK,kfabs(beta1L),kmadd(JacPDupwindNthSymm2trK,kfabs(beta2L),kmul(JacPDupwindNthSymm3trK,kfabs(beta3L))))))));
    
    At11rhsL = 
      kadd(At11rhsL,kmadd(beta1L,JacPDupwindNthAnti1At11,kmadd(beta2L,JacPDupwindNthAnti2At11,kmadd(beta3L,JacPDupwindNthAnti3At11,kmadd(JacPDupwindNthSymm1At11,kfabs(beta1L),kmadd(JacPDupwindNthSymm2At11,kfabs(beta2L),kmul(JacPDupwindNthSymm3At11,kfabs(beta3L))))))));
    
    At12rhsL = 
      kadd(At12rhsL,kmadd(beta1L,JacPDupwindNthAnti1At12,kmadd(beta2L,JacPDupwindNthAnti2At12,kmadd(beta3L,JacPDupwindNthAnti3At12,kmadd(JacPDupwindNthSymm1At12,kfabs(beta1L),kmadd(JacPDupwindNthSymm2At12,kfabs(beta2L),kmul(JacPDupwindNthSymm3At12,kfabs(beta3L))))))));
    
    At13rhsL = 
      kadd(At13rhsL,kmadd(beta1L,JacPDupwindNthAnti1At13,kmadd(beta2L,JacPDupwindNthAnti2At13,kmadd(beta3L,JacPDupwindNthAnti3At13,kmadd(JacPDupwindNthSymm1At13,kfabs(beta1L),kmadd(JacPDupwindNthSymm2At13,kfabs(beta2L),kmul(JacPDupwindNthSymm3At13,kfabs(beta3L))))))));
    
    At22rhsL = 
      kadd(At22rhsL,kmadd(beta1L,JacPDupwindNthAnti1At22,kmadd(beta2L,JacPDupwindNthAnti2At22,kmadd(beta3L,JacPDupwindNthAnti3At22,kmadd(JacPDupwindNthSymm1At22,kfabs(beta1L),kmadd(JacPDupwindNthSymm2At22,kfabs(beta2L),kmul(JacPDupwindNthSymm3At22,kfabs(beta3L))))))));
    
    At23rhsL = 
      kadd(At23rhsL,kmadd(beta1L,JacPDupwindNthAnti1At23,kmadd(beta2L,JacPDupwindNthAnti2At23,kmadd(beta3L,JacPDupwindNthAnti3At23,kmadd(JacPDupwindNthSymm1At23,kfabs(beta1L),kmadd(JacPDupwindNthSymm2At23,kfabs(beta2L),kmul(JacPDupwindNthSymm3At23,kfabs(beta3L))))))));
    
    At33rhsL = 
      kadd(At33rhsL,kmadd(beta1L,JacPDupwindNthAnti1At33,kmadd(beta2L,JacPDupwindNthAnti2At33,kmadd(beta3L,JacPDupwindNthAnti3At33,kmadd(JacPDupwindNthSymm1At33,kfabs(beta1L),kmadd(JacPDupwindNthSymm2At33,kfabs(beta2L),kmul(JacPDupwindNthSymm3At33,kfabs(beta3L))))))));
    
    alpharhsL = 
      kmadd(kmadd(beta1L,JacPDupwindNthAnti1alpha,kmadd(beta2L,JacPDupwindNthAnti2alpha,kmadd(beta3L,JacPDupwindNthAnti3alpha,kmadd(JacPDupwindNthSymm1alpha,kfabs(beta1L),kmadd(JacPDupwindNthSymm2alpha,kfabs(beta2L),kmul(JacPDupwindNthSymm3alpha,kfabs(beta3L))))))),ToReal(LapseAdvectionCoeff),alpharhsL);
    
    ArhsL = 
      kmadd(ToReal(LapseACoeff),kmsub(kmadd(beta1L,JacPDupwindNthAnti1A,kmadd(beta2L,JacPDupwindNthAnti2A,kmadd(beta3L,JacPDupwindNthAnti3A,kmadd(JacPDupwindNthSymm1A,kfabs(beta1L),kmadd(JacPDupwindNthSymm2A,kfabs(beta2L),kmul(JacPDupwindNthSymm3A,kfabs(beta3L))))))),ToReal(LapseAdvectionCoeff),kmul(kmadd(beta1L,JacPDupwindNthAnti1trK,kmadd(beta2L,JacPDupwindNthAnti2trK,kmadd(beta3L,JacPDupwindNthAnti3trK,kmadd(JacPDupwindNthSymm1trK,kfabs(beta1L),kmadd(JacPDupwindNthSymm2trK,kfabs(beta2L),kmul(JacPDupwindNthSymm3trK,kfabs(beta3L))))))),kadd(ToReal(-1),ToReal(LapseAdvectionCoeff)))),ArhsL);
    
    beta1rhsL = 
      kmadd(kmadd(beta1L,JacPDupwindNthAnti1beta1,kmadd(beta2L,JacPDupwindNthAnti2beta1,kmadd(beta3L,JacPDupwindNthAnti3beta1,kmadd(JacPDupwindNthSymm1beta1,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta1,kfabs(beta2L),kmul(JacPDupwindNthSymm3beta1,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),beta1rhsL);
    
    beta2rhsL = 
      kmadd(kmadd(beta1L,JacPDupwindNthAnti1beta2,kmadd(beta2L,JacPDupwindNthAnti2beta2,kmadd(beta3L,JacPDupwindNthAnti3beta2,kmadd(JacPDupwindNthSymm1beta2,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta2,kfabs(beta2L),kmul(JacPDupwindNthSymm3beta2,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),beta2rhsL);
    
    beta3rhsL = 
      kmadd(kmadd(beta1L,JacPDupwindNthAnti1beta3,kmadd(beta2L,JacPDupwindNthAnti2beta3,kmadd(beta3L,JacPDupwindNthAnti3beta3,kmadd(JacPDupwindNthSymm1beta3,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta3,kfabs(beta2L),kmul(JacPDupwindNthSymm3beta3,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),beta3rhsL);
    
    B1rhsL = 
      kmadd(kmsub(kmadd(beta1L,JacPDupwindNthAnti1B1,kmadd(beta2L,JacPDupwindNthAnti2B1,kmadd(beta3L,JacPDupwindNthAnti3B1,kmadd(JacPDupwindNthSymm1B1,kfabs(beta1L),kmadd(JacPDupwindNthSymm2B1,kfabs(beta2L),kmul(JacPDupwindNthSymm3B1,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),kmul(kmadd(beta1L,JacPDupwindNthAnti1Xt1,kmadd(beta2L,JacPDupwindNthAnti2Xt1,kmadd(beta3L,JacPDupwindNthAnti3Xt1,kmadd(JacPDupwindNthSymm1Xt1,kfabs(beta1L),kmadd(JacPDupwindNthSymm2Xt1,kfabs(beta2L),kmul(JacPDupwindNthSymm3Xt1,kfabs(beta3L))))))),kadd(ToReal(-1),ToReal(ShiftAdvectionCoeff)))),ToReal(ShiftBCoeff),B1rhsL);
    
    B2rhsL = 
      kmadd(kmsub(kmadd(beta1L,JacPDupwindNthAnti1B2,kmadd(beta2L,JacPDupwindNthAnti2B2,kmadd(beta3L,JacPDupwindNthAnti3B2,kmadd(JacPDupwindNthSymm1B2,kfabs(beta1L),kmadd(JacPDupwindNthSymm2B2,kfabs(beta2L),kmul(JacPDupwindNthSymm3B2,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),kmul(kmadd(beta1L,JacPDupwindNthAnti1Xt2,kmadd(beta2L,JacPDupwindNthAnti2Xt2,kmadd(beta3L,JacPDupwindNthAnti3Xt2,kmadd(JacPDupwindNthSymm1Xt2,kfabs(beta1L),kmadd(JacPDupwindNthSymm2Xt2,kfabs(beta2L),kmul(JacPDupwindNthSymm3Xt2,kfabs(beta3L))))))),kadd(ToReal(-1),ToReal(ShiftAdvectionCoeff)))),ToReal(ShiftBCoeff),B2rhsL);
    
    B3rhsL = 
      kmadd(kmsub(kmadd(beta1L,JacPDupwindNthAnti1B3,kmadd(beta2L,JacPDupwindNthAnti2B3,kmadd(beta3L,JacPDupwindNthAnti3B3,kmadd(JacPDupwindNthSymm1B3,kfabs(beta1L),kmadd(JacPDupwindNthSymm2B3,kfabs(beta2L),kmul(JacPDupwindNthSymm3B3,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),kmul(kmadd(beta1L,JacPDupwindNthAnti1Xt3,kmadd(beta2L,JacPDupwindNthAnti2Xt3,kmadd(beta3L,JacPDupwindNthAnti3Xt3,kmadd(JacPDupwindNthSymm1Xt3,kfabs(beta1L),kmadd(JacPDupwindNthSymm2Xt3,kfabs(beta2L),kmul(JacPDupwindNthSymm3Xt3,kfabs(beta3L))))))),kadd(ToReal(-1),ToReal(ShiftAdvectionCoeff)))),ToReal(ShiftBCoeff),B3rhsL);
    
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
  LC_ENDLOOP3VEC (ML_BSSN_Advect);
}

extern "C" void ML_BSSN_Advect(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_Advect_Body");
  }
  
  if (cctk_iteration % ML_BSSN_Advect_calc_every != ML_BSSN_Advect_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_BSSN::ML_curv","ML_BSSN::ML_curvrhs","ML_BSSN::ML_dtlapse","ML_BSSN::ML_dtlapserhs","ML_BSSN::ML_dtshift","ML_BSSN::ML_dtshiftrhs","ML_BSSN::ML_Gamma","ML_BSSN::ML_Gammarhs","ML_BSSN::ML_lapse","ML_BSSN::ML_lapserhs","ML_BSSN::ML_log_confac","ML_BSSN::ML_log_confacrhs","ML_BSSN::ML_metric","ML_BSSN::ML_metricrhs","ML_BSSN::ML_shift","ML_BSSN::ML_shiftrhs","ML_BSSN::ML_trace_curv","ML_BSSN::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Advect", 18, groups);
  
  switch(fdOrder)
  {
    case 2:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Advect", 2, 2, 2);
      break;
    
    case 4:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Advect", 3, 3, 3);
      break;
    
    case 6:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Advect", 4, 4, 4);
      break;
    
    case 8:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Advect", 5, 5, 5);
      break;
  }
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_Advect_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Advect_Body");
  }
}
