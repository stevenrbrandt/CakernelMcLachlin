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

extern "C" void ML_BSSN_UPW_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
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

static void ML_BSSN_UPW_RHS1_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  LC_LOOP3VEC (ML_BSSN_UPW_RHS1,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC AL = vec_load(A[index]);
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC At11L = vec_load(At11[index]);
    CCTK_REAL_VEC At12L = vec_load(At12[index]);
    CCTK_REAL_VEC At13L = vec_load(At13[index]);
    CCTK_REAL_VEC At22L = vec_load(At22[index]);
    CCTK_REAL_VEC At23L = vec_load(At23[index]);
    CCTK_REAL_VEC At33L = vec_load(At33[index]);
    CCTK_REAL_VEC B1L = vec_load(B1[index]);
    CCTK_REAL_VEC B2L = vec_load(B2[index]);
    CCTK_REAL_VEC B3L = vec_load(B3[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC phiL = vec_load(phi[index]);
    CCTK_REAL_VEC rL = vec_load(r[index]);
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    CCTK_REAL_VEC Xt1L = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L = vec_load(Xt3[index]);
    
    CCTK_REAL_VEC eTttL, eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL;
    
    if (*stress_energy_state)
    {
      eTttL = vec_load(eTtt[index]);
      eTtxL = vec_load(eTtx[index]);
      eTtyL = vec_load(eTty[index]);
      eTtzL = vec_load(eTtz[index]);
      eTxxL = vec_load(eTxx[index]);
      eTxyL = vec_load(eTxy[index]);
      eTxzL = vec_load(eTxz[index]);
      eTyyL = vec_load(eTyy[index]);
      eTyzL = vec_load(eTyz[index]);
      eTzzL = vec_load(eTzz[index]);
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
    
    CCTK_REAL_VEC dJ111L, dJ112L, dJ113L, dJ122L, dJ123L, dJ133L, dJ211L, dJ212L, dJ213L, dJ222L, dJ223L, dJ233L, dJ311L, dJ312L, dJ313L, dJ322L, dJ323L, dJ333L, J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L;
    
    if (use_jacobian)
    {
      dJ111L = vec_load(dJ111[index]);
      dJ112L = vec_load(dJ112[index]);
      dJ113L = vec_load(dJ113[index]);
      dJ122L = vec_load(dJ122[index]);
      dJ123L = vec_load(dJ123[index]);
      dJ133L = vec_load(dJ133[index]);
      dJ211L = vec_load(dJ211[index]);
      dJ212L = vec_load(dJ212[index]);
      dJ213L = vec_load(dJ213[index]);
      dJ222L = vec_load(dJ222[index]);
      dJ223L = vec_load(dJ223[index]);
      dJ233L = vec_load(dJ233[index]);
      dJ311L = vec_load(dJ311[index]);
      dJ312L = vec_load(dJ312[index]);
      dJ313L = vec_load(dJ313[index]);
      dJ322L = vec_load(dJ322[index]);
      dJ323L = vec_load(dJ323[index]);
      dJ333L = vec_load(dJ333[index]);
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
    CCTK_REAL_VEC PDstandardNth1alpha;
    CCTK_REAL_VEC PDstandardNth2alpha;
    CCTK_REAL_VEC PDstandardNth3alpha;
    CCTK_REAL_VEC PDstandardNth11alpha;
    CCTK_REAL_VEC PDstandardNth22alpha;
    CCTK_REAL_VEC PDstandardNth33alpha;
    CCTK_REAL_VEC PDstandardNth12alpha;
    CCTK_REAL_VEC PDstandardNth13alpha;
    CCTK_REAL_VEC PDstandardNth23alpha;
    CCTK_REAL_VEC PDstandardNth1beta1;
    CCTK_REAL_VEC PDstandardNth2beta1;
    CCTK_REAL_VEC PDstandardNth3beta1;
    CCTK_REAL_VEC PDstandardNth11beta1;
    CCTK_REAL_VEC PDstandardNth22beta1;
    CCTK_REAL_VEC PDstandardNth33beta1;
    CCTK_REAL_VEC PDstandardNth12beta1;
    CCTK_REAL_VEC PDstandardNth13beta1;
    CCTK_REAL_VEC PDstandardNth23beta1;
    CCTK_REAL_VEC PDstandardNth1beta2;
    CCTK_REAL_VEC PDstandardNth2beta2;
    CCTK_REAL_VEC PDstandardNth3beta2;
    CCTK_REAL_VEC PDstandardNth11beta2;
    CCTK_REAL_VEC PDstandardNth22beta2;
    CCTK_REAL_VEC PDstandardNth33beta2;
    CCTK_REAL_VEC PDstandardNth12beta2;
    CCTK_REAL_VEC PDstandardNth13beta2;
    CCTK_REAL_VEC PDstandardNth23beta2;
    CCTK_REAL_VEC PDstandardNth1beta3;
    CCTK_REAL_VEC PDstandardNth2beta3;
    CCTK_REAL_VEC PDstandardNth3beta3;
    CCTK_REAL_VEC PDstandardNth11beta3;
    CCTK_REAL_VEC PDstandardNth22beta3;
    CCTK_REAL_VEC PDstandardNth33beta3;
    CCTK_REAL_VEC PDstandardNth12beta3;
    CCTK_REAL_VEC PDstandardNth13beta3;
    CCTK_REAL_VEC PDstandardNth23beta3;
    CCTK_REAL_VEC PDstandardNth1gt11;
    CCTK_REAL_VEC PDstandardNth2gt11;
    CCTK_REAL_VEC PDstandardNth3gt11;
    CCTK_REAL_VEC PDstandardNth1gt12;
    CCTK_REAL_VEC PDstandardNth2gt12;
    CCTK_REAL_VEC PDstandardNth3gt12;
    CCTK_REAL_VEC PDstandardNth1gt13;
    CCTK_REAL_VEC PDstandardNth2gt13;
    CCTK_REAL_VEC PDstandardNth3gt13;
    CCTK_REAL_VEC PDstandardNth1gt22;
    CCTK_REAL_VEC PDstandardNth2gt22;
    CCTK_REAL_VEC PDstandardNth3gt22;
    CCTK_REAL_VEC PDstandardNth1gt23;
    CCTK_REAL_VEC PDstandardNth2gt23;
    CCTK_REAL_VEC PDstandardNth3gt23;
    CCTK_REAL_VEC PDstandardNth1gt33;
    CCTK_REAL_VEC PDstandardNth2gt33;
    CCTK_REAL_VEC PDstandardNth3gt33;
    CCTK_REAL_VEC PDstandardNth1phi;
    CCTK_REAL_VEC PDstandardNth2phi;
    CCTK_REAL_VEC PDstandardNth3phi;
    CCTK_REAL_VEC PDstandardNth1trK;
    CCTK_REAL_VEC PDstandardNth2trK;
    CCTK_REAL_VEC PDstandardNth3trK;
    
    switch(fdOrder)
    {
      case 2:
        PDstandardNth1alpha = PDstandardNthfdOrder21(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder22(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder23(&alpha[index]);
        PDstandardNth11alpha = PDstandardNthfdOrder211(&alpha[index]);
        PDstandardNth22alpha = PDstandardNthfdOrder222(&alpha[index]);
        PDstandardNth33alpha = PDstandardNthfdOrder233(&alpha[index]);
        PDstandardNth12alpha = PDstandardNthfdOrder212(&alpha[index]);
        PDstandardNth13alpha = PDstandardNthfdOrder213(&alpha[index]);
        PDstandardNth23alpha = PDstandardNthfdOrder223(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder21(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder22(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder23(&beta1[index]);
        PDstandardNth11beta1 = PDstandardNthfdOrder211(&beta1[index]);
        PDstandardNth22beta1 = PDstandardNthfdOrder222(&beta1[index]);
        PDstandardNth33beta1 = PDstandardNthfdOrder233(&beta1[index]);
        PDstandardNth12beta1 = PDstandardNthfdOrder212(&beta1[index]);
        PDstandardNth13beta1 = PDstandardNthfdOrder213(&beta1[index]);
        PDstandardNth23beta1 = PDstandardNthfdOrder223(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder21(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder22(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder23(&beta2[index]);
        PDstandardNth11beta2 = PDstandardNthfdOrder211(&beta2[index]);
        PDstandardNth22beta2 = PDstandardNthfdOrder222(&beta2[index]);
        PDstandardNth33beta2 = PDstandardNthfdOrder233(&beta2[index]);
        PDstandardNth12beta2 = PDstandardNthfdOrder212(&beta2[index]);
        PDstandardNth13beta2 = PDstandardNthfdOrder213(&beta2[index]);
        PDstandardNth23beta2 = PDstandardNthfdOrder223(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder21(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder22(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder23(&beta3[index]);
        PDstandardNth11beta3 = PDstandardNthfdOrder211(&beta3[index]);
        PDstandardNth22beta3 = PDstandardNthfdOrder222(&beta3[index]);
        PDstandardNth33beta3 = PDstandardNthfdOrder233(&beta3[index]);
        PDstandardNth12beta3 = PDstandardNthfdOrder212(&beta3[index]);
        PDstandardNth13beta3 = PDstandardNthfdOrder213(&beta3[index]);
        PDstandardNth23beta3 = PDstandardNthfdOrder223(&beta3[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder21(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder22(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder23(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder21(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder22(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder23(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder21(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder22(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder23(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder21(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder22(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder23(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder21(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder22(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder23(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder21(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder22(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder23(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder21(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder22(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder23(&phi[index]);
        PDstandardNth1trK = PDstandardNthfdOrder21(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder22(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder23(&trK[index]);
        break;
      
      case 4:
        PDstandardNth1alpha = PDstandardNthfdOrder41(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder42(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder43(&alpha[index]);
        PDstandardNth11alpha = PDstandardNthfdOrder411(&alpha[index]);
        PDstandardNth22alpha = PDstandardNthfdOrder422(&alpha[index]);
        PDstandardNth33alpha = PDstandardNthfdOrder433(&alpha[index]);
        PDstandardNth12alpha = PDstandardNthfdOrder412(&alpha[index]);
        PDstandardNth13alpha = PDstandardNthfdOrder413(&alpha[index]);
        PDstandardNth23alpha = PDstandardNthfdOrder423(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder41(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder42(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder43(&beta1[index]);
        PDstandardNth11beta1 = PDstandardNthfdOrder411(&beta1[index]);
        PDstandardNth22beta1 = PDstandardNthfdOrder422(&beta1[index]);
        PDstandardNth33beta1 = PDstandardNthfdOrder433(&beta1[index]);
        PDstandardNth12beta1 = PDstandardNthfdOrder412(&beta1[index]);
        PDstandardNth13beta1 = PDstandardNthfdOrder413(&beta1[index]);
        PDstandardNth23beta1 = PDstandardNthfdOrder423(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder41(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder42(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder43(&beta2[index]);
        PDstandardNth11beta2 = PDstandardNthfdOrder411(&beta2[index]);
        PDstandardNth22beta2 = PDstandardNthfdOrder422(&beta2[index]);
        PDstandardNth33beta2 = PDstandardNthfdOrder433(&beta2[index]);
        PDstandardNth12beta2 = PDstandardNthfdOrder412(&beta2[index]);
        PDstandardNth13beta2 = PDstandardNthfdOrder413(&beta2[index]);
        PDstandardNth23beta2 = PDstandardNthfdOrder423(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder41(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder42(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder43(&beta3[index]);
        PDstandardNth11beta3 = PDstandardNthfdOrder411(&beta3[index]);
        PDstandardNth22beta3 = PDstandardNthfdOrder422(&beta3[index]);
        PDstandardNth33beta3 = PDstandardNthfdOrder433(&beta3[index]);
        PDstandardNth12beta3 = PDstandardNthfdOrder412(&beta3[index]);
        PDstandardNth13beta3 = PDstandardNthfdOrder413(&beta3[index]);
        PDstandardNth23beta3 = PDstandardNthfdOrder423(&beta3[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder41(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder42(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder43(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder41(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder42(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder43(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder41(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder42(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder43(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder41(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder42(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder43(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder41(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder42(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder43(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder41(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder42(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder43(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder41(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder42(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder43(&phi[index]);
        PDstandardNth1trK = PDstandardNthfdOrder41(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder42(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder43(&trK[index]);
        break;
      
      case 6:
        PDstandardNth1alpha = PDstandardNthfdOrder61(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder62(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder63(&alpha[index]);
        PDstandardNth11alpha = PDstandardNthfdOrder611(&alpha[index]);
        PDstandardNth22alpha = PDstandardNthfdOrder622(&alpha[index]);
        PDstandardNth33alpha = PDstandardNthfdOrder633(&alpha[index]);
        PDstandardNth12alpha = PDstandardNthfdOrder612(&alpha[index]);
        PDstandardNth13alpha = PDstandardNthfdOrder613(&alpha[index]);
        PDstandardNth23alpha = PDstandardNthfdOrder623(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder61(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder62(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder63(&beta1[index]);
        PDstandardNth11beta1 = PDstandardNthfdOrder611(&beta1[index]);
        PDstandardNth22beta1 = PDstandardNthfdOrder622(&beta1[index]);
        PDstandardNth33beta1 = PDstandardNthfdOrder633(&beta1[index]);
        PDstandardNth12beta1 = PDstandardNthfdOrder612(&beta1[index]);
        PDstandardNth13beta1 = PDstandardNthfdOrder613(&beta1[index]);
        PDstandardNth23beta1 = PDstandardNthfdOrder623(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder61(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder62(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder63(&beta2[index]);
        PDstandardNth11beta2 = PDstandardNthfdOrder611(&beta2[index]);
        PDstandardNth22beta2 = PDstandardNthfdOrder622(&beta2[index]);
        PDstandardNth33beta2 = PDstandardNthfdOrder633(&beta2[index]);
        PDstandardNth12beta2 = PDstandardNthfdOrder612(&beta2[index]);
        PDstandardNth13beta2 = PDstandardNthfdOrder613(&beta2[index]);
        PDstandardNth23beta2 = PDstandardNthfdOrder623(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder61(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder62(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder63(&beta3[index]);
        PDstandardNth11beta3 = PDstandardNthfdOrder611(&beta3[index]);
        PDstandardNth22beta3 = PDstandardNthfdOrder622(&beta3[index]);
        PDstandardNth33beta3 = PDstandardNthfdOrder633(&beta3[index]);
        PDstandardNth12beta3 = PDstandardNthfdOrder612(&beta3[index]);
        PDstandardNth13beta3 = PDstandardNthfdOrder613(&beta3[index]);
        PDstandardNth23beta3 = PDstandardNthfdOrder623(&beta3[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder61(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder62(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder63(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder61(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder62(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder63(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder61(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder62(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder63(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder61(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder62(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder63(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder61(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder62(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder63(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder61(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder62(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder63(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder61(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder62(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder63(&phi[index]);
        PDstandardNth1trK = PDstandardNthfdOrder61(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder62(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder63(&trK[index]);
        break;
      
      case 8:
        PDstandardNth1alpha = PDstandardNthfdOrder81(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder82(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder83(&alpha[index]);
        PDstandardNth11alpha = PDstandardNthfdOrder811(&alpha[index]);
        PDstandardNth22alpha = PDstandardNthfdOrder822(&alpha[index]);
        PDstandardNth33alpha = PDstandardNthfdOrder833(&alpha[index]);
        PDstandardNth12alpha = PDstandardNthfdOrder812(&alpha[index]);
        PDstandardNth13alpha = PDstandardNthfdOrder813(&alpha[index]);
        PDstandardNth23alpha = PDstandardNthfdOrder823(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder81(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder82(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder83(&beta1[index]);
        PDstandardNth11beta1 = PDstandardNthfdOrder811(&beta1[index]);
        PDstandardNth22beta1 = PDstandardNthfdOrder822(&beta1[index]);
        PDstandardNth33beta1 = PDstandardNthfdOrder833(&beta1[index]);
        PDstandardNth12beta1 = PDstandardNthfdOrder812(&beta1[index]);
        PDstandardNth13beta1 = PDstandardNthfdOrder813(&beta1[index]);
        PDstandardNth23beta1 = PDstandardNthfdOrder823(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder81(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder82(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder83(&beta2[index]);
        PDstandardNth11beta2 = PDstandardNthfdOrder811(&beta2[index]);
        PDstandardNth22beta2 = PDstandardNthfdOrder822(&beta2[index]);
        PDstandardNth33beta2 = PDstandardNthfdOrder833(&beta2[index]);
        PDstandardNth12beta2 = PDstandardNthfdOrder812(&beta2[index]);
        PDstandardNth13beta2 = PDstandardNthfdOrder813(&beta2[index]);
        PDstandardNth23beta2 = PDstandardNthfdOrder823(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder81(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder82(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder83(&beta3[index]);
        PDstandardNth11beta3 = PDstandardNthfdOrder811(&beta3[index]);
        PDstandardNth22beta3 = PDstandardNthfdOrder822(&beta3[index]);
        PDstandardNth33beta3 = PDstandardNthfdOrder833(&beta3[index]);
        PDstandardNth12beta3 = PDstandardNthfdOrder812(&beta3[index]);
        PDstandardNth13beta3 = PDstandardNthfdOrder813(&beta3[index]);
        PDstandardNth23beta3 = PDstandardNthfdOrder823(&beta3[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder81(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder82(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder83(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder81(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder82(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder83(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder81(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder82(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder83(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder81(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder82(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder83(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder81(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder82(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder83(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder81(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder82(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder83(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder81(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder82(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder83(&phi[index]);
        PDstandardNth1trK = PDstandardNthfdOrder81(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder82(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder83(&trK[index]);
        break;
    }
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    CCTK_REAL_VEC JacPDstandardNth11alpha;
    CCTK_REAL_VEC JacPDstandardNth11beta1;
    CCTK_REAL_VEC JacPDstandardNth11beta2;
    CCTK_REAL_VEC JacPDstandardNth11beta3;
    CCTK_REAL_VEC JacPDstandardNth12alpha;
    CCTK_REAL_VEC JacPDstandardNth12beta1;
    CCTK_REAL_VEC JacPDstandardNth12beta2;
    CCTK_REAL_VEC JacPDstandardNth12beta3;
    CCTK_REAL_VEC JacPDstandardNth13alpha;
    CCTK_REAL_VEC JacPDstandardNth13beta1;
    CCTK_REAL_VEC JacPDstandardNth13beta2;
    CCTK_REAL_VEC JacPDstandardNth13beta3;
    CCTK_REAL_VEC JacPDstandardNth1alpha;
    CCTK_REAL_VEC JacPDstandardNth1beta1;
    CCTK_REAL_VEC JacPDstandardNth1beta2;
    CCTK_REAL_VEC JacPDstandardNth1beta3;
    CCTK_REAL_VEC JacPDstandardNth1gt11;
    CCTK_REAL_VEC JacPDstandardNth1gt12;
    CCTK_REAL_VEC JacPDstandardNth1gt13;
    CCTK_REAL_VEC JacPDstandardNth1gt22;
    CCTK_REAL_VEC JacPDstandardNth1gt23;
    CCTK_REAL_VEC JacPDstandardNth1gt33;
    CCTK_REAL_VEC JacPDstandardNth1phi;
    CCTK_REAL_VEC JacPDstandardNth1trK;
    CCTK_REAL_VEC JacPDstandardNth21alpha;
    CCTK_REAL_VEC JacPDstandardNth21beta1;
    CCTK_REAL_VEC JacPDstandardNth21beta2;
    CCTK_REAL_VEC JacPDstandardNth21beta3;
    CCTK_REAL_VEC JacPDstandardNth22alpha;
    CCTK_REAL_VEC JacPDstandardNth22beta1;
    CCTK_REAL_VEC JacPDstandardNth22beta2;
    CCTK_REAL_VEC JacPDstandardNth22beta3;
    CCTK_REAL_VEC JacPDstandardNth23alpha;
    CCTK_REAL_VEC JacPDstandardNth23beta1;
    CCTK_REAL_VEC JacPDstandardNth23beta2;
    CCTK_REAL_VEC JacPDstandardNth23beta3;
    CCTK_REAL_VEC JacPDstandardNth2alpha;
    CCTK_REAL_VEC JacPDstandardNth2beta1;
    CCTK_REAL_VEC JacPDstandardNth2beta2;
    CCTK_REAL_VEC JacPDstandardNth2beta3;
    CCTK_REAL_VEC JacPDstandardNth2gt11;
    CCTK_REAL_VEC JacPDstandardNth2gt12;
    CCTK_REAL_VEC JacPDstandardNth2gt13;
    CCTK_REAL_VEC JacPDstandardNth2gt22;
    CCTK_REAL_VEC JacPDstandardNth2gt23;
    CCTK_REAL_VEC JacPDstandardNth2gt33;
    CCTK_REAL_VEC JacPDstandardNth2phi;
    CCTK_REAL_VEC JacPDstandardNth2trK;
    CCTK_REAL_VEC JacPDstandardNth31alpha;
    CCTK_REAL_VEC JacPDstandardNth31beta1;
    CCTK_REAL_VEC JacPDstandardNth31beta2;
    CCTK_REAL_VEC JacPDstandardNth31beta3;
    CCTK_REAL_VEC JacPDstandardNth32alpha;
    CCTK_REAL_VEC JacPDstandardNth32beta1;
    CCTK_REAL_VEC JacPDstandardNth32beta2;
    CCTK_REAL_VEC JacPDstandardNth32beta3;
    CCTK_REAL_VEC JacPDstandardNth33alpha;
    CCTK_REAL_VEC JacPDstandardNth33beta1;
    CCTK_REAL_VEC JacPDstandardNth33beta2;
    CCTK_REAL_VEC JacPDstandardNth33beta3;
    CCTK_REAL_VEC JacPDstandardNth3alpha;
    CCTK_REAL_VEC JacPDstandardNth3beta1;
    CCTK_REAL_VEC JacPDstandardNth3beta2;
    CCTK_REAL_VEC JacPDstandardNth3beta3;
    CCTK_REAL_VEC JacPDstandardNth3gt11;
    CCTK_REAL_VEC JacPDstandardNth3gt12;
    CCTK_REAL_VEC JacPDstandardNth3gt13;
    CCTK_REAL_VEC JacPDstandardNth3gt22;
    CCTK_REAL_VEC JacPDstandardNth3gt23;
    CCTK_REAL_VEC JacPDstandardNth3gt33;
    CCTK_REAL_VEC JacPDstandardNth3phi;
    CCTK_REAL_VEC JacPDstandardNth3trK;
    
    if (use_jacobian)
    {
      JacPDstandardNth1alpha = 
        kmadd(J11L,PDstandardNth1alpha,kmadd(J21L,PDstandardNth2alpha,kmul(J31L,PDstandardNth3alpha)));
      
      JacPDstandardNth1beta1 = 
        kmadd(J11L,PDstandardNth1beta1,kmadd(J21L,PDstandardNth2beta1,kmul(J31L,PDstandardNth3beta1)));
      
      JacPDstandardNth1beta2 = 
        kmadd(J11L,PDstandardNth1beta2,kmadd(J21L,PDstandardNth2beta2,kmul(J31L,PDstandardNth3beta2)));
      
      JacPDstandardNth1beta3 = 
        kmadd(J11L,PDstandardNth1beta3,kmadd(J21L,PDstandardNth2beta3,kmul(J31L,PDstandardNth3beta3)));
      
      JacPDstandardNth1gt11 = 
        kmadd(J11L,PDstandardNth1gt11,kmadd(J21L,PDstandardNth2gt11,kmul(J31L,PDstandardNth3gt11)));
      
      JacPDstandardNth1gt12 = 
        kmadd(J11L,PDstandardNth1gt12,kmadd(J21L,PDstandardNth2gt12,kmul(J31L,PDstandardNth3gt12)));
      
      JacPDstandardNth1gt13 = 
        kmadd(J11L,PDstandardNth1gt13,kmadd(J21L,PDstandardNth2gt13,kmul(J31L,PDstandardNth3gt13)));
      
      JacPDstandardNth1gt22 = 
        kmadd(J11L,PDstandardNth1gt22,kmadd(J21L,PDstandardNth2gt22,kmul(J31L,PDstandardNth3gt22)));
      
      JacPDstandardNth1gt23 = 
        kmadd(J11L,PDstandardNth1gt23,kmadd(J21L,PDstandardNth2gt23,kmul(J31L,PDstandardNth3gt23)));
      
      JacPDstandardNth1gt33 = 
        kmadd(J11L,PDstandardNth1gt33,kmadd(J21L,PDstandardNth2gt33,kmul(J31L,PDstandardNth3gt33)));
      
      JacPDstandardNth1phi = 
        kmadd(J11L,PDstandardNth1phi,kmadd(J21L,PDstandardNth2phi,kmul(J31L,PDstandardNth3phi)));
      
      JacPDstandardNth1trK = 
        kmadd(J11L,PDstandardNth1trK,kmadd(J21L,PDstandardNth2trK,kmul(J31L,PDstandardNth3trK)));
      
      JacPDstandardNth2alpha = 
        kmadd(J12L,PDstandardNth1alpha,kmadd(J22L,PDstandardNth2alpha,kmul(J32L,PDstandardNth3alpha)));
      
      JacPDstandardNth2beta1 = 
        kmadd(J12L,PDstandardNth1beta1,kmadd(J22L,PDstandardNth2beta1,kmul(J32L,PDstandardNth3beta1)));
      
      JacPDstandardNth2beta2 = 
        kmadd(J12L,PDstandardNth1beta2,kmadd(J22L,PDstandardNth2beta2,kmul(J32L,PDstandardNth3beta2)));
      
      JacPDstandardNth2beta3 = 
        kmadd(J12L,PDstandardNth1beta3,kmadd(J22L,PDstandardNth2beta3,kmul(J32L,PDstandardNth3beta3)));
      
      JacPDstandardNth2gt11 = 
        kmadd(J12L,PDstandardNth1gt11,kmadd(J22L,PDstandardNth2gt11,kmul(J32L,PDstandardNth3gt11)));
      
      JacPDstandardNth2gt12 = 
        kmadd(J12L,PDstandardNth1gt12,kmadd(J22L,PDstandardNth2gt12,kmul(J32L,PDstandardNth3gt12)));
      
      JacPDstandardNth2gt13 = 
        kmadd(J12L,PDstandardNth1gt13,kmadd(J22L,PDstandardNth2gt13,kmul(J32L,PDstandardNth3gt13)));
      
      JacPDstandardNth2gt22 = 
        kmadd(J12L,PDstandardNth1gt22,kmadd(J22L,PDstandardNth2gt22,kmul(J32L,PDstandardNth3gt22)));
      
      JacPDstandardNth2gt23 = 
        kmadd(J12L,PDstandardNth1gt23,kmadd(J22L,PDstandardNth2gt23,kmul(J32L,PDstandardNth3gt23)));
      
      JacPDstandardNth2gt33 = 
        kmadd(J12L,PDstandardNth1gt33,kmadd(J22L,PDstandardNth2gt33,kmul(J32L,PDstandardNth3gt33)));
      
      JacPDstandardNth2phi = 
        kmadd(J12L,PDstandardNth1phi,kmadd(J22L,PDstandardNth2phi,kmul(J32L,PDstandardNth3phi)));
      
      JacPDstandardNth2trK = 
        kmadd(J12L,PDstandardNth1trK,kmadd(J22L,PDstandardNth2trK,kmul(J32L,PDstandardNth3trK)));
      
      JacPDstandardNth3alpha = 
        kmadd(J13L,PDstandardNth1alpha,kmadd(J23L,PDstandardNth2alpha,kmul(J33L,PDstandardNth3alpha)));
      
      JacPDstandardNth3beta1 = 
        kmadd(J13L,PDstandardNth1beta1,kmadd(J23L,PDstandardNth2beta1,kmul(J33L,PDstandardNth3beta1)));
      
      JacPDstandardNth3beta2 = 
        kmadd(J13L,PDstandardNth1beta2,kmadd(J23L,PDstandardNth2beta2,kmul(J33L,PDstandardNth3beta2)));
      
      JacPDstandardNth3beta3 = 
        kmadd(J13L,PDstandardNth1beta3,kmadd(J23L,PDstandardNth2beta3,kmul(J33L,PDstandardNth3beta3)));
      
      JacPDstandardNth3gt11 = 
        kmadd(J13L,PDstandardNth1gt11,kmadd(J23L,PDstandardNth2gt11,kmul(J33L,PDstandardNth3gt11)));
      
      JacPDstandardNth3gt12 = 
        kmadd(J13L,PDstandardNth1gt12,kmadd(J23L,PDstandardNth2gt12,kmul(J33L,PDstandardNth3gt12)));
      
      JacPDstandardNth3gt13 = 
        kmadd(J13L,PDstandardNth1gt13,kmadd(J23L,PDstandardNth2gt13,kmul(J33L,PDstandardNth3gt13)));
      
      JacPDstandardNth3gt22 = 
        kmadd(J13L,PDstandardNth1gt22,kmadd(J23L,PDstandardNth2gt22,kmul(J33L,PDstandardNth3gt22)));
      
      JacPDstandardNth3gt23 = 
        kmadd(J13L,PDstandardNth1gt23,kmadd(J23L,PDstandardNth2gt23,kmul(J33L,PDstandardNth3gt23)));
      
      JacPDstandardNth3gt33 = 
        kmadd(J13L,PDstandardNth1gt33,kmadd(J23L,PDstandardNth2gt33,kmul(J33L,PDstandardNth3gt33)));
      
      JacPDstandardNth3phi = 
        kmadd(J13L,PDstandardNth1phi,kmadd(J23L,PDstandardNth2phi,kmul(J33L,PDstandardNth3phi)));
      
      JacPDstandardNth3trK = 
        kmadd(J13L,PDstandardNth1trK,kmadd(J23L,PDstandardNth2trK,kmul(J33L,PDstandardNth3trK)));
      
      JacPDstandardNth11alpha = 
        kmadd(dJ111L,PDstandardNth1alpha,kmadd(dJ211L,PDstandardNth2alpha,kmadd(dJ311L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,SQR(J11L),kmadd(PDstandardNth22alpha,SQR(J21L),kmadd(PDstandardNth33alpha,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha)),kmul(J21L,kmul(J31L,PDstandardNth23alpha))),ToReal(2))))))));
      
      JacPDstandardNth11beta1 = 
        kmadd(dJ111L,PDstandardNth1beta1,kmadd(dJ211L,PDstandardNth2beta1,kmadd(dJ311L,PDstandardNth3beta1,kmadd(PDstandardNth11beta1,SQR(J11L),kmadd(PDstandardNth22beta1,SQR(J21L),kmadd(PDstandardNth33beta1,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12beta1,kmul(J31L,PDstandardNth13beta1)),kmul(J21L,kmul(J31L,PDstandardNth23beta1))),ToReal(2))))))));
      
      JacPDstandardNth11beta2 = 
        kmadd(dJ111L,PDstandardNth1beta2,kmadd(dJ211L,PDstandardNth2beta2,kmadd(dJ311L,PDstandardNth3beta2,kmadd(PDstandardNth11beta2,SQR(J11L),kmadd(PDstandardNth22beta2,SQR(J21L),kmadd(PDstandardNth33beta2,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12beta2,kmul(J31L,PDstandardNth13beta2)),kmul(J21L,kmul(J31L,PDstandardNth23beta2))),ToReal(2))))))));
      
      JacPDstandardNth11beta3 = 
        kmadd(dJ111L,PDstandardNth1beta3,kmadd(dJ211L,PDstandardNth2beta3,kmadd(dJ311L,PDstandardNth3beta3,kmadd(PDstandardNth11beta3,SQR(J11L),kmadd(PDstandardNth22beta3,SQR(J21L),kmadd(PDstandardNth33beta3,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12beta3,kmul(J31L,PDstandardNth13beta3)),kmul(J21L,kmul(J31L,PDstandardNth23beta3))),ToReal(2))))))));
      
      JacPDstandardNth22alpha = 
        kmadd(dJ122L,PDstandardNth1alpha,kmadd(dJ222L,PDstandardNth2alpha,kmadd(dJ322L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,SQR(J12L),kmadd(PDstandardNth22alpha,SQR(J22L),kmadd(PDstandardNth33alpha,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha)),kmul(J22L,kmul(J32L,PDstandardNth23alpha))),ToReal(2))))))));
      
      JacPDstandardNth22beta1 = 
        kmadd(dJ122L,PDstandardNth1beta1,kmadd(dJ222L,PDstandardNth2beta1,kmadd(dJ322L,PDstandardNth3beta1,kmadd(PDstandardNth11beta1,SQR(J12L),kmadd(PDstandardNth22beta1,SQR(J22L),kmadd(PDstandardNth33beta1,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12beta1,kmul(J32L,PDstandardNth13beta1)),kmul(J22L,kmul(J32L,PDstandardNth23beta1))),ToReal(2))))))));
      
      JacPDstandardNth22beta2 = 
        kmadd(dJ122L,PDstandardNth1beta2,kmadd(dJ222L,PDstandardNth2beta2,kmadd(dJ322L,PDstandardNth3beta2,kmadd(PDstandardNth11beta2,SQR(J12L),kmadd(PDstandardNth22beta2,SQR(J22L),kmadd(PDstandardNth33beta2,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12beta2,kmul(J32L,PDstandardNth13beta2)),kmul(J22L,kmul(J32L,PDstandardNth23beta2))),ToReal(2))))))));
      
      JacPDstandardNth22beta3 = 
        kmadd(dJ122L,PDstandardNth1beta3,kmadd(dJ222L,PDstandardNth2beta3,kmadd(dJ322L,PDstandardNth3beta3,kmadd(PDstandardNth11beta3,SQR(J12L),kmadd(PDstandardNth22beta3,SQR(J22L),kmadd(PDstandardNth33beta3,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12beta3,kmul(J32L,PDstandardNth13beta3)),kmul(J22L,kmul(J32L,PDstandardNth23beta3))),ToReal(2))))))));
      
      JacPDstandardNth33alpha = 
        kmadd(dJ133L,PDstandardNth1alpha,kmadd(dJ233L,PDstandardNth2alpha,kmadd(dJ333L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,SQR(J13L),kmadd(PDstandardNth22alpha,SQR(J23L),kmadd(PDstandardNth33alpha,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmul(J23L,kmul(J33L,PDstandardNth23alpha))),ToReal(2))))))));
      
      JacPDstandardNth33beta1 = 
        kmadd(dJ133L,PDstandardNth1beta1,kmadd(dJ233L,PDstandardNth2beta1,kmadd(dJ333L,PDstandardNth3beta1,kmadd(PDstandardNth11beta1,SQR(J13L),kmadd(PDstandardNth22beta1,SQR(J23L),kmadd(PDstandardNth33beta1,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12beta1,kmul(J33L,PDstandardNth13beta1)),kmul(J23L,kmul(J33L,PDstandardNth23beta1))),ToReal(2))))))));
      
      JacPDstandardNth33beta2 = 
        kmadd(dJ133L,PDstandardNth1beta2,kmadd(dJ233L,PDstandardNth2beta2,kmadd(dJ333L,PDstandardNth3beta2,kmadd(PDstandardNth11beta2,SQR(J13L),kmadd(PDstandardNth22beta2,SQR(J23L),kmadd(PDstandardNth33beta2,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12beta2,kmul(J33L,PDstandardNth13beta2)),kmul(J23L,kmul(J33L,PDstandardNth23beta2))),ToReal(2))))))));
      
      JacPDstandardNth33beta3 = 
        kmadd(dJ133L,PDstandardNth1beta3,kmadd(dJ233L,PDstandardNth2beta3,kmadd(dJ333L,PDstandardNth3beta3,kmadd(PDstandardNth11beta3,SQR(J13L),kmadd(PDstandardNth22beta3,SQR(J23L),kmadd(PDstandardNth33beta3,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12beta3,kmul(J33L,PDstandardNth13beta3)),kmul(J23L,kmul(J33L,PDstandardNth23beta3))),ToReal(2))))))));
      
      JacPDstandardNth12alpha = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11alpha,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha))),kmadd(J11L,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha)),kmadd(dJ112L,PDstandardNth1alpha,kmadd(J22L,kmadd(J21L,PDstandardNth22alpha,kmul(J31L,PDstandardNth23alpha)),kmadd(dJ212L,PDstandardNth2alpha,kmadd(J32L,kmadd(J21L,PDstandardNth23alpha,kmul(J31L,PDstandardNth33alpha)),kmul(dJ312L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth12beta1 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11beta1,kmadd(J21L,PDstandardNth12beta1,kmul(J31L,PDstandardNth13beta1))),kmadd(J11L,kmadd(J22L,PDstandardNth12beta1,kmul(J32L,PDstandardNth13beta1)),kmadd(dJ112L,PDstandardNth1beta1,kmadd(J22L,kmadd(J21L,PDstandardNth22beta1,kmul(J31L,PDstandardNth23beta1)),kmadd(dJ212L,PDstandardNth2beta1,kmadd(J32L,kmadd(J21L,PDstandardNth23beta1,kmul(J31L,PDstandardNth33beta1)),kmul(dJ312L,PDstandardNth3beta1)))))));
      
      JacPDstandardNth12beta2 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11beta2,kmadd(J21L,PDstandardNth12beta2,kmul(J31L,PDstandardNth13beta2))),kmadd(J11L,kmadd(J22L,PDstandardNth12beta2,kmul(J32L,PDstandardNth13beta2)),kmadd(dJ112L,PDstandardNth1beta2,kmadd(J22L,kmadd(J21L,PDstandardNth22beta2,kmul(J31L,PDstandardNth23beta2)),kmadd(dJ212L,PDstandardNth2beta2,kmadd(J32L,kmadd(J21L,PDstandardNth23beta2,kmul(J31L,PDstandardNth33beta2)),kmul(dJ312L,PDstandardNth3beta2)))))));
      
      JacPDstandardNth12beta3 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11beta3,kmadd(J21L,PDstandardNth12beta3,kmul(J31L,PDstandardNth13beta3))),kmadd(J11L,kmadd(J22L,PDstandardNth12beta3,kmul(J32L,PDstandardNth13beta3)),kmadd(dJ112L,PDstandardNth1beta3,kmadd(J22L,kmadd(J21L,PDstandardNth22beta3,kmul(J31L,PDstandardNth23beta3)),kmadd(dJ212L,PDstandardNth2beta3,kmadd(J32L,kmadd(J21L,PDstandardNth23beta3,kmul(J31L,PDstandardNth33beta3)),kmul(dJ312L,PDstandardNth3beta3)))))));
      
      JacPDstandardNth13alpha = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11alpha,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha))),kmadd(J11L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmadd(dJ113L,PDstandardNth1alpha,kmadd(J23L,kmadd(J21L,PDstandardNth22alpha,kmul(J31L,PDstandardNth23alpha)),kmadd(dJ213L,PDstandardNth2alpha,kmadd(J33L,kmadd(J21L,PDstandardNth23alpha,kmul(J31L,PDstandardNth33alpha)),kmul(dJ313L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth13beta1 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11beta1,kmadd(J21L,PDstandardNth12beta1,kmul(J31L,PDstandardNth13beta1))),kmadd(J11L,kmadd(J23L,PDstandardNth12beta1,kmul(J33L,PDstandardNth13beta1)),kmadd(dJ113L,PDstandardNth1beta1,kmadd(J23L,kmadd(J21L,PDstandardNth22beta1,kmul(J31L,PDstandardNth23beta1)),kmadd(dJ213L,PDstandardNth2beta1,kmadd(J33L,kmadd(J21L,PDstandardNth23beta1,kmul(J31L,PDstandardNth33beta1)),kmul(dJ313L,PDstandardNth3beta1)))))));
      
      JacPDstandardNth13beta2 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11beta2,kmadd(J21L,PDstandardNth12beta2,kmul(J31L,PDstandardNth13beta2))),kmadd(J11L,kmadd(J23L,PDstandardNth12beta2,kmul(J33L,PDstandardNth13beta2)),kmadd(dJ113L,PDstandardNth1beta2,kmadd(J23L,kmadd(J21L,PDstandardNth22beta2,kmul(J31L,PDstandardNth23beta2)),kmadd(dJ213L,PDstandardNth2beta2,kmadd(J33L,kmadd(J21L,PDstandardNth23beta2,kmul(J31L,PDstandardNth33beta2)),kmul(dJ313L,PDstandardNth3beta2)))))));
      
      JacPDstandardNth13beta3 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11beta3,kmadd(J21L,PDstandardNth12beta3,kmul(J31L,PDstandardNth13beta3))),kmadd(J11L,kmadd(J23L,PDstandardNth12beta3,kmul(J33L,PDstandardNth13beta3)),kmadd(dJ113L,PDstandardNth1beta3,kmadd(J23L,kmadd(J21L,PDstandardNth22beta3,kmul(J31L,PDstandardNth23beta3)),kmadd(dJ213L,PDstandardNth2beta3,kmadd(J33L,kmadd(J21L,PDstandardNth23beta3,kmul(J31L,PDstandardNth33beta3)),kmul(dJ313L,PDstandardNth3beta3)))))));
      
      JacPDstandardNth21alpha = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11alpha,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha))),kmadd(J11L,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha)),kmadd(dJ112L,PDstandardNth1alpha,kmadd(J22L,kmadd(J21L,PDstandardNth22alpha,kmul(J31L,PDstandardNth23alpha)),kmadd(dJ212L,PDstandardNth2alpha,kmadd(J32L,kmadd(J21L,PDstandardNth23alpha,kmul(J31L,PDstandardNth33alpha)),kmul(dJ312L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth21beta1 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11beta1,kmadd(J21L,PDstandardNth12beta1,kmul(J31L,PDstandardNth13beta1))),kmadd(J11L,kmadd(J22L,PDstandardNth12beta1,kmul(J32L,PDstandardNth13beta1)),kmadd(dJ112L,PDstandardNth1beta1,kmadd(J22L,kmadd(J21L,PDstandardNth22beta1,kmul(J31L,PDstandardNth23beta1)),kmadd(dJ212L,PDstandardNth2beta1,kmadd(J32L,kmadd(J21L,PDstandardNth23beta1,kmul(J31L,PDstandardNth33beta1)),kmul(dJ312L,PDstandardNth3beta1)))))));
      
      JacPDstandardNth21beta2 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11beta2,kmadd(J21L,PDstandardNth12beta2,kmul(J31L,PDstandardNth13beta2))),kmadd(J11L,kmadd(J22L,PDstandardNth12beta2,kmul(J32L,PDstandardNth13beta2)),kmadd(dJ112L,PDstandardNth1beta2,kmadd(J22L,kmadd(J21L,PDstandardNth22beta2,kmul(J31L,PDstandardNth23beta2)),kmadd(dJ212L,PDstandardNth2beta2,kmadd(J32L,kmadd(J21L,PDstandardNth23beta2,kmul(J31L,PDstandardNth33beta2)),kmul(dJ312L,PDstandardNth3beta2)))))));
      
      JacPDstandardNth21beta3 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11beta3,kmadd(J21L,PDstandardNth12beta3,kmul(J31L,PDstandardNth13beta3))),kmadd(J11L,kmadd(J22L,PDstandardNth12beta3,kmul(J32L,PDstandardNth13beta3)),kmadd(dJ112L,PDstandardNth1beta3,kmadd(J22L,kmadd(J21L,PDstandardNth22beta3,kmul(J31L,PDstandardNth23beta3)),kmadd(dJ212L,PDstandardNth2beta3,kmadd(J32L,kmadd(J21L,PDstandardNth23beta3,kmul(J31L,PDstandardNth33beta3)),kmul(dJ312L,PDstandardNth3beta3)))))));
      
      JacPDstandardNth23alpha = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11alpha,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha))),kmadd(J12L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmadd(dJ123L,PDstandardNth1alpha,kmadd(J23L,kmadd(J22L,PDstandardNth22alpha,kmul(J32L,PDstandardNth23alpha)),kmadd(dJ223L,PDstandardNth2alpha,kmadd(J33L,kmadd(J22L,PDstandardNth23alpha,kmul(J32L,PDstandardNth33alpha)),kmul(dJ323L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth23beta1 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11beta1,kmadd(J22L,PDstandardNth12beta1,kmul(J32L,PDstandardNth13beta1))),kmadd(J12L,kmadd(J23L,PDstandardNth12beta1,kmul(J33L,PDstandardNth13beta1)),kmadd(dJ123L,PDstandardNth1beta1,kmadd(J23L,kmadd(J22L,PDstandardNth22beta1,kmul(J32L,PDstandardNth23beta1)),kmadd(dJ223L,PDstandardNth2beta1,kmadd(J33L,kmadd(J22L,PDstandardNth23beta1,kmul(J32L,PDstandardNth33beta1)),kmul(dJ323L,PDstandardNth3beta1)))))));
      
      JacPDstandardNth23beta2 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11beta2,kmadd(J22L,PDstandardNth12beta2,kmul(J32L,PDstandardNth13beta2))),kmadd(J12L,kmadd(J23L,PDstandardNth12beta2,kmul(J33L,PDstandardNth13beta2)),kmadd(dJ123L,PDstandardNth1beta2,kmadd(J23L,kmadd(J22L,PDstandardNth22beta2,kmul(J32L,PDstandardNth23beta2)),kmadd(dJ223L,PDstandardNth2beta2,kmadd(J33L,kmadd(J22L,PDstandardNth23beta2,kmul(J32L,PDstandardNth33beta2)),kmul(dJ323L,PDstandardNth3beta2)))))));
      
      JacPDstandardNth23beta3 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11beta3,kmadd(J22L,PDstandardNth12beta3,kmul(J32L,PDstandardNth13beta3))),kmadd(J12L,kmadd(J23L,PDstandardNth12beta3,kmul(J33L,PDstandardNth13beta3)),kmadd(dJ123L,PDstandardNth1beta3,kmadd(J23L,kmadd(J22L,PDstandardNth22beta3,kmul(J32L,PDstandardNth23beta3)),kmadd(dJ223L,PDstandardNth2beta3,kmadd(J33L,kmadd(J22L,PDstandardNth23beta3,kmul(J32L,PDstandardNth33beta3)),kmul(dJ323L,PDstandardNth3beta3)))))));
      
      JacPDstandardNth31alpha = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11alpha,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha))),kmadd(J11L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmadd(dJ113L,PDstandardNth1alpha,kmadd(J23L,kmadd(J21L,PDstandardNth22alpha,kmul(J31L,PDstandardNth23alpha)),kmadd(dJ213L,PDstandardNth2alpha,kmadd(J33L,kmadd(J21L,PDstandardNth23alpha,kmul(J31L,PDstandardNth33alpha)),kmul(dJ313L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth31beta1 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11beta1,kmadd(J21L,PDstandardNth12beta1,kmul(J31L,PDstandardNth13beta1))),kmadd(J11L,kmadd(J23L,PDstandardNth12beta1,kmul(J33L,PDstandardNth13beta1)),kmadd(dJ113L,PDstandardNth1beta1,kmadd(J23L,kmadd(J21L,PDstandardNth22beta1,kmul(J31L,PDstandardNth23beta1)),kmadd(dJ213L,PDstandardNth2beta1,kmadd(J33L,kmadd(J21L,PDstandardNth23beta1,kmul(J31L,PDstandardNth33beta1)),kmul(dJ313L,PDstandardNth3beta1)))))));
      
      JacPDstandardNth31beta2 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11beta2,kmadd(J21L,PDstandardNth12beta2,kmul(J31L,PDstandardNth13beta2))),kmadd(J11L,kmadd(J23L,PDstandardNth12beta2,kmul(J33L,PDstandardNth13beta2)),kmadd(dJ113L,PDstandardNth1beta2,kmadd(J23L,kmadd(J21L,PDstandardNth22beta2,kmul(J31L,PDstandardNth23beta2)),kmadd(dJ213L,PDstandardNth2beta2,kmadd(J33L,kmadd(J21L,PDstandardNth23beta2,kmul(J31L,PDstandardNth33beta2)),kmul(dJ313L,PDstandardNth3beta2)))))));
      
      JacPDstandardNth31beta3 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11beta3,kmadd(J21L,PDstandardNth12beta3,kmul(J31L,PDstandardNth13beta3))),kmadd(J11L,kmadd(J23L,PDstandardNth12beta3,kmul(J33L,PDstandardNth13beta3)),kmadd(dJ113L,PDstandardNth1beta3,kmadd(J23L,kmadd(J21L,PDstandardNth22beta3,kmul(J31L,PDstandardNth23beta3)),kmadd(dJ213L,PDstandardNth2beta3,kmadd(J33L,kmadd(J21L,PDstandardNth23beta3,kmul(J31L,PDstandardNth33beta3)),kmul(dJ313L,PDstandardNth3beta3)))))));
      
      JacPDstandardNth32alpha = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11alpha,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha))),kmadd(J12L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmadd(dJ123L,PDstandardNth1alpha,kmadd(J23L,kmadd(J22L,PDstandardNth22alpha,kmul(J32L,PDstandardNth23alpha)),kmadd(dJ223L,PDstandardNth2alpha,kmadd(J33L,kmadd(J22L,PDstandardNth23alpha,kmul(J32L,PDstandardNth33alpha)),kmul(dJ323L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth32beta1 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11beta1,kmadd(J22L,PDstandardNth12beta1,kmul(J32L,PDstandardNth13beta1))),kmadd(J12L,kmadd(J23L,PDstandardNth12beta1,kmul(J33L,PDstandardNth13beta1)),kmadd(dJ123L,PDstandardNth1beta1,kmadd(J23L,kmadd(J22L,PDstandardNth22beta1,kmul(J32L,PDstandardNth23beta1)),kmadd(dJ223L,PDstandardNth2beta1,kmadd(J33L,kmadd(J22L,PDstandardNth23beta1,kmul(J32L,PDstandardNth33beta1)),kmul(dJ323L,PDstandardNth3beta1)))))));
      
      JacPDstandardNth32beta2 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11beta2,kmadd(J22L,PDstandardNth12beta2,kmul(J32L,PDstandardNth13beta2))),kmadd(J12L,kmadd(J23L,PDstandardNth12beta2,kmul(J33L,PDstandardNth13beta2)),kmadd(dJ123L,PDstandardNth1beta2,kmadd(J23L,kmadd(J22L,PDstandardNth22beta2,kmul(J32L,PDstandardNth23beta2)),kmadd(dJ223L,PDstandardNth2beta2,kmadd(J33L,kmadd(J22L,PDstandardNth23beta2,kmul(J32L,PDstandardNth33beta2)),kmul(dJ323L,PDstandardNth3beta2)))))));
      
      JacPDstandardNth32beta3 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11beta3,kmadd(J22L,PDstandardNth12beta3,kmul(J32L,PDstandardNth13beta3))),kmadd(J12L,kmadd(J23L,PDstandardNth12beta3,kmul(J33L,PDstandardNth13beta3)),kmadd(dJ123L,PDstandardNth1beta3,kmadd(J23L,kmadd(J22L,PDstandardNth22beta3,kmul(J32L,PDstandardNth23beta3)),kmadd(dJ223L,PDstandardNth2beta3,kmadd(J33L,kmadd(J22L,PDstandardNth23beta3,kmul(J32L,PDstandardNth33beta3)),kmul(dJ323L,PDstandardNth3beta3)))))));
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
    
    CCTK_REAL_VEC detgt = ToReal(1);
    
    CCTK_REAL_VEC gtu11 = kmul(INV(detgt),kmsub(gt22L,gt33L,SQR(gt23L)));
    
    CCTK_REAL_VEC gtu12 = 
      kmul(INV(detgt),kmsub(gt13L,gt23L,kmul(gt12L,gt33L)));
    
    CCTK_REAL_VEC gtu13 = 
      kmul(INV(detgt),kmsub(gt12L,gt23L,kmul(gt13L,gt22L)));
    
    CCTK_REAL_VEC gtu22 = kmul(INV(detgt),kmsub(gt11L,gt33L,SQR(gt13L)));
    
    CCTK_REAL_VEC gtu23 = 
      kmul(INV(detgt),kmsub(gt12L,gt13L,kmul(gt11L,gt23L)));
    
    CCTK_REAL_VEC gtu33 = kmul(INV(detgt),kmsub(gt11L,gt22L,SQR(gt12L)));
    
    CCTK_REAL_VEC Gtl111 = kmul(JacPDstandardNth1gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl112 = kmul(JacPDstandardNth2gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl113 = kmul(JacPDstandardNth3gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl122 = 
      kmadd(JacPDstandardNth1gt22,ToReal(-0.5),JacPDstandardNth2gt12);
    
    CCTK_REAL_VEC Gtl123 = 
      kmul(kadd(JacPDstandardNth2gt13,ksub(JacPDstandardNth3gt12,JacPDstandardNth1gt23)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl133 = 
      kmadd(JacPDstandardNth1gt33,ToReal(-0.5),JacPDstandardNth3gt13);
    
    CCTK_REAL_VEC Gtl211 = 
      kmadd(JacPDstandardNth2gt11,ToReal(-0.5),JacPDstandardNth1gt12);
    
    CCTK_REAL_VEC Gtl212 = kmul(JacPDstandardNth1gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl213 = 
      kmul(kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth3gt12,JacPDstandardNth2gt13)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl222 = kmul(JacPDstandardNth2gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl223 = kmul(JacPDstandardNth3gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl233 = 
      kmadd(JacPDstandardNth2gt33,ToReal(-0.5),JacPDstandardNth3gt23);
    
    CCTK_REAL_VEC Gtl311 = 
      kmadd(JacPDstandardNth3gt11,ToReal(-0.5),JacPDstandardNth1gt13);
    
    CCTK_REAL_VEC Gtl312 = 
      kmul(kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth2gt13,JacPDstandardNth3gt12)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl313 = kmul(JacPDstandardNth1gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl322 = 
      kmadd(JacPDstandardNth3gt22,ToReal(-0.5),JacPDstandardNth2gt23);
    
    CCTK_REAL_VEC Gtl323 = kmul(JacPDstandardNth2gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl333 = kmul(JacPDstandardNth3gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gt111 = 
      kmadd(Gtl111,gtu11,kmadd(Gtl211,gtu12,kmul(Gtl311,gtu13)));
    
    CCTK_REAL_VEC Gt211 = 
      kmadd(Gtl111,gtu12,kmadd(Gtl211,gtu22,kmul(Gtl311,gtu23)));
    
    CCTK_REAL_VEC Gt311 = 
      kmadd(Gtl111,gtu13,kmadd(Gtl211,gtu23,kmul(Gtl311,gtu33)));
    
    CCTK_REAL_VEC Gt112 = 
      kmadd(Gtl112,gtu11,kmadd(Gtl212,gtu12,kmul(Gtl312,gtu13)));
    
    CCTK_REAL_VEC Gt212 = 
      kmadd(Gtl112,gtu12,kmadd(Gtl212,gtu22,kmul(Gtl312,gtu23)));
    
    CCTK_REAL_VEC Gt312 = 
      kmadd(Gtl112,gtu13,kmadd(Gtl212,gtu23,kmul(Gtl312,gtu33)));
    
    CCTK_REAL_VEC Gt113 = 
      kmadd(Gtl113,gtu11,kmadd(Gtl213,gtu12,kmul(Gtl313,gtu13)));
    
    CCTK_REAL_VEC Gt213 = 
      kmadd(Gtl113,gtu12,kmadd(Gtl213,gtu22,kmul(Gtl313,gtu23)));
    
    CCTK_REAL_VEC Gt313 = 
      kmadd(Gtl113,gtu13,kmadd(Gtl213,gtu23,kmul(Gtl313,gtu33)));
    
    CCTK_REAL_VEC Gt122 = 
      kmadd(Gtl122,gtu11,kmadd(Gtl222,gtu12,kmul(Gtl322,gtu13)));
    
    CCTK_REAL_VEC Gt222 = 
      kmadd(Gtl122,gtu12,kmadd(Gtl222,gtu22,kmul(Gtl322,gtu23)));
    
    CCTK_REAL_VEC Gt322 = 
      kmadd(Gtl122,gtu13,kmadd(Gtl222,gtu23,kmul(Gtl322,gtu33)));
    
    CCTK_REAL_VEC Gt123 = 
      kmadd(Gtl123,gtu11,kmadd(Gtl223,gtu12,kmul(Gtl323,gtu13)));
    
    CCTK_REAL_VEC Gt223 = 
      kmadd(Gtl123,gtu12,kmadd(Gtl223,gtu22,kmul(Gtl323,gtu23)));
    
    CCTK_REAL_VEC Gt323 = 
      kmadd(Gtl123,gtu13,kmadd(Gtl223,gtu23,kmul(Gtl323,gtu33)));
    
    CCTK_REAL_VEC Gt133 = 
      kmadd(Gtl133,gtu11,kmadd(Gtl233,gtu12,kmul(Gtl333,gtu13)));
    
    CCTK_REAL_VEC Gt233 = 
      kmadd(Gtl133,gtu12,kmadd(Gtl233,gtu22,kmul(Gtl333,gtu23)));
    
    CCTK_REAL_VEC Gt333 = 
      kmadd(Gtl133,gtu13,kmadd(Gtl233,gtu23,kmul(Gtl333,gtu33)));
    
    CCTK_REAL_VEC Xtn1 = 
      kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmadd(Gt133,gtu33,kmul(kmadd(Gt112,gtu12,kmadd(Gt113,gtu13,kmul(Gt123,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Xtn2 = 
      kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmadd(Gt233,gtu33,kmul(kmadd(Gt212,gtu12,kmadd(Gt213,gtu13,kmul(Gt223,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Xtn3 = 
      kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmadd(Gt333,gtu33,kmul(kmadd(Gt312,gtu12,kmadd(Gt313,gtu13,kmul(Gt323,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC fac1 = 
      IfThen(conformalMethod,kmul(INV(phiL),ToReal(-0.5)),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 = kmul(fac1,JacPDstandardNth1phi);
    
    CCTK_REAL_VEC cdphi2 = kmul(fac1,JacPDstandardNth2phi);
    
    CCTK_REAL_VEC cdphi3 = kmul(fac1,JacPDstandardNth3phi);
    
    CCTK_REAL_VEC Atm11 = 
      kmadd(At11L,gtu11,kmadd(At12L,gtu12,kmul(At13L,gtu13)));
    
    CCTK_REAL_VEC Atm21 = 
      kmadd(At11L,gtu12,kmadd(At12L,gtu22,kmul(At13L,gtu23)));
    
    CCTK_REAL_VEC Atm31 = 
      kmadd(At11L,gtu13,kmadd(At12L,gtu23,kmul(At13L,gtu33)));
    
    CCTK_REAL_VEC Atm12 = 
      kmadd(At12L,gtu11,kmadd(At22L,gtu12,kmul(At23L,gtu13)));
    
    CCTK_REAL_VEC Atm22 = 
      kmadd(At12L,gtu12,kmadd(At22L,gtu22,kmul(At23L,gtu23)));
    
    CCTK_REAL_VEC Atm32 = 
      kmadd(At12L,gtu13,kmadd(At22L,gtu23,kmul(At23L,gtu33)));
    
    CCTK_REAL_VEC Atm13 = 
      kmadd(At13L,gtu11,kmadd(At23L,gtu12,kmul(At33L,gtu13)));
    
    CCTK_REAL_VEC Atm23 = 
      kmadd(At13L,gtu12,kmadd(At23L,gtu22,kmul(At33L,gtu23)));
    
    CCTK_REAL_VEC Atm33 = 
      kmadd(At13L,gtu13,kmadd(At23L,gtu23,kmul(At33L,gtu33)));
    
    CCTK_REAL_VEC Atu11 = 
      kmadd(Atm11,gtu11,kmadd(Atm12,gtu12,kmul(Atm13,gtu13)));
    
    CCTK_REAL_VEC Atu12 = 
      kmadd(Atm11,gtu12,kmadd(Atm12,gtu22,kmul(Atm13,gtu23)));
    
    CCTK_REAL_VEC Atu13 = 
      kmadd(Atm11,gtu13,kmadd(Atm12,gtu23,kmul(Atm13,gtu33)));
    
    CCTK_REAL_VEC Atu22 = 
      kmadd(Atm21,gtu12,kmadd(Atm22,gtu22,kmul(Atm23,gtu23)));
    
    CCTK_REAL_VEC Atu23 = 
      kmadd(Atm21,gtu13,kmadd(Atm22,gtu23,kmul(Atm23,gtu33)));
    
    CCTK_REAL_VEC Atu33 = 
      kmadd(Atm31,gtu13,kmadd(Atm32,gtu23,kmul(Atm33,gtu33)));
    
    CCTK_REAL_VEC e4phi = 
      IfThen(conformalMethod,INV(SQR(phiL)),kexp(kmul(phiL,ToReal(4))));
    
    CCTK_REAL_VEC em4phi = INV(e4phi);
    
    CCTK_REAL_VEC rho = 
      kmul(INV(SQR(alphaL)),kadd(eTttL,kmadd(eTxxL,SQR(beta1L),kmadd(eTyyL,SQR(beta2L),kmadd(eTzzL,SQR(beta3L),kmadd(kmadd(beta2L,eTtyL,kmul(beta3L,eTtzL)),ToReal(-2),kmul(kmadd(beta2L,kmul(beta3L,eTyzL),kmul(beta1L,kmadd(beta2L,eTxyL,kmsub(beta3L,eTxzL,eTtxL)))),ToReal(2))))))));
    
    CCTK_REAL_VEC S1 = 
      kmul(INV(alphaL),kmadd(beta1L,eTxxL,kmadd(beta2L,eTxyL,kmsub(beta3L,eTxzL,eTtxL))));
    
    CCTK_REAL_VEC S2 = 
      kmul(INV(alphaL),kmadd(beta1L,eTxyL,kmadd(beta2L,eTyyL,kmsub(beta3L,eTyzL,eTtyL))));
    
    CCTK_REAL_VEC S3 = 
      kmul(INV(alphaL),kmadd(beta1L,eTxzL,kmadd(beta2L,eTyzL,kmsub(beta3L,eTzzL,eTtzL))));
    
    CCTK_REAL_VEC trS = 
      kmul(em4phi,kmadd(eTxxL,gtu11,kmadd(eTyyL,gtu22,kmadd(eTzzL,gtu33,kmul(kmadd(eTxyL,gtu12,kmadd(eTxzL,gtu13,kmul(eTyzL,gtu23))),ToReal(2))))));
    
    CCTK_REAL_VEC phirhsL = 
      IfThen(conformalMethod,kmul(phiL,kmadd(kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),ToReal(-0.333333333333333333333333333333),kmul(alphaL,kmul(trKL,ToReal(0.333333333333333333333333333333))))),kmadd(alphaL,kmul(trKL,ToReal(-0.166666666666666666666666666667)),kmul(kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),ToReal(0.166666666666666666666666666667))));
    
    CCTK_REAL_VEC gt11rhsL = 
      kmul(ToReal(-0.666666666666666666666666666667),kmadd(kmadd(gt12L,JacPDstandardNth1beta2,kmul(gt13L,JacPDstandardNth1beta3)),ToReal(-3),kmadd(gt11L,kadd(JacPDstandardNth2beta2,kmadd(JacPDstandardNth1beta1,ToReal(-2),JacPDstandardNth3beta3)),kmul(alphaL,kmul(At11L,ToReal(3))))));
    
    CCTK_REAL_VEC gt12rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At12L,ToReal(-6)),kmadd(gt12L,kadd(JacPDstandardNth1beta1,kmadd(JacPDstandardNth3beta3,ToReal(-2),JacPDstandardNth2beta2)),kmul(kmadd(gt22L,JacPDstandardNth1beta2,kmadd(gt23L,JacPDstandardNth1beta3,kmadd(gt11L,JacPDstandardNth2beta1,kmul(gt13L,JacPDstandardNth2beta3)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt13rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At13L,ToReal(-6)),kmadd(gt13L,kadd(JacPDstandardNth1beta1,kmadd(JacPDstandardNth2beta2,ToReal(-2),JacPDstandardNth3beta3)),kmul(kmadd(gt23L,JacPDstandardNth1beta2,kmadd(gt33L,JacPDstandardNth1beta3,kmadd(gt11L,JacPDstandardNth3beta1,kmul(gt12L,JacPDstandardNth3beta2)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt22rhsL = 
      kmul(ToReal(-0.666666666666666666666666666667),kmadd(kmadd(gt12L,JacPDstandardNth2beta1,kmul(gt23L,JacPDstandardNth2beta3)),ToReal(-3),kmadd(gt22L,kadd(JacPDstandardNth1beta1,kmadd(JacPDstandardNth2beta2,ToReal(-2),JacPDstandardNth3beta3)),kmul(alphaL,kmul(At22L,ToReal(3))))));
    
    CCTK_REAL_VEC gt23rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At23L,ToReal(-6)),kmadd(gt23L,kadd(JacPDstandardNth2beta2,kmadd(JacPDstandardNth1beta1,ToReal(-2),JacPDstandardNth3beta3)),kmul(kmadd(gt13L,JacPDstandardNth2beta1,kmadd(gt33L,JacPDstandardNth2beta3,kmadd(gt12L,JacPDstandardNth3beta1,kmul(gt22L,JacPDstandardNth3beta2)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt33rhsL = 
      kmul(ToReal(-0.666666666666666666666666666667),kmadd(kmadd(gt13L,JacPDstandardNth3beta1,kmul(gt23L,JacPDstandardNth3beta2)),ToReal(-3),kmadd(gt33L,kadd(JacPDstandardNth1beta1,kmadd(JacPDstandardNth3beta3,ToReal(-2),JacPDstandardNth2beta2)),kmul(alphaL,kmul(At33L,ToReal(3))))));
    
    CCTK_REAL_VEC dotXt1 = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(kmadd(Atu11,JacPDstandardNth1alpha,kmadd(Atu12,JacPDstandardNth2alpha,kmul(Atu13,JacPDstandardNth3alpha))),ToReal(-6),kmadd(kmadd(JacPDstandardNth2beta1,Xtn2,kmul(JacPDstandardNth3beta1,Xtn3)),ToReal(-3),kmadd(Xtn1,kmsub(kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3),ToReal(2),JacPDstandardNth1beta1),kmadd(kmadd(gtu12,JacPDstandardNth12beta1,kmadd(gtu13,JacPDstandardNth13beta1,kmadd(gtu22,JacPDstandardNth22beta1,kmadd(gtu23,kadd(JacPDstandardNth23beta1,JacPDstandardNth32beta1),kmul(gtu33,JacPDstandardNth33beta1))))),ToReal(3),kmadd(gtu11,kadd(JacPDstandardNth12beta2,kmadd(JacPDstandardNth11beta1,ToReal(4),JacPDstandardNth13beta3)),kmadd(gtu12,kadd(JacPDstandardNth22beta2,kmadd(JacPDstandardNth21beta1,ToReal(4),JacPDstandardNth23beta3)),kmadd(gtu13,kadd(JacPDstandardNth32beta2,kmadd(JacPDstandardNth31beta1,ToReal(4),JacPDstandardNth33beta3)),kmul(alphaL,kmadd(kmadd(gtu11,S1,kmadd(gtu12,S2,kmul(gtu13,S3))),ToReal(-150.7964473723100754462068823974161384415),kmadd(kmadd(gtu11,JacPDstandardNth1trK,kmadd(gtu12,JacPDstandardNth2trK,kmul(gtu13,JacPDstandardNth3trK))),ToReal(-4),kmadd(ToReal(6),kmadd(Atu22,Gt122,kmadd(Atu33,Gt133,kmul(Atu11,kmadd(cdphi1,ToReal(6),Gt111)))),kmadd(kmadd(Atu23,Gt123,kmul(Atu12,kmadd(cdphi2,ToReal(3),Gt112))),ToReal(12),kmul(Atu13,kmadd(Gt113,ToReal(12),kmul(cdphi3,ToReal(36)))))))))))))))));
    
    CCTK_REAL_VEC dotXt2 = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(kmadd(Atu12,JacPDstandardNth1alpha,kmadd(Atu22,JacPDstandardNth2alpha,kmul(Atu23,JacPDstandardNth3alpha))),ToReal(-6),kmadd(kmadd(JacPDstandardNth1beta2,Xtn1,kmul(JacPDstandardNth3beta2,Xtn3)),ToReal(-3),kmadd(Xtn2,kmsub(kadd(JacPDstandardNth1beta1,JacPDstandardNth3beta3),ToReal(2),JacPDstandardNth2beta2),kmadd(kmadd(gtu11,JacPDstandardNth11beta2,kmadd(gtu23,JacPDstandardNth23beta2,kmadd(gtu13,kadd(JacPDstandardNth13beta2,JacPDstandardNth31beta2),kmul(gtu33,JacPDstandardNth33beta2)))),ToReal(3),kmadd(gtu12,kadd(JacPDstandardNth11beta1,kadd(JacPDstandardNth13beta3,kmadd(JacPDstandardNth21beta2,ToReal(3),kmul(JacPDstandardNth12beta2,ToReal(4))))),kmadd(gtu22,kadd(JacPDstandardNth21beta1,kmadd(JacPDstandardNth22beta2,ToReal(4),JacPDstandardNth23beta3)),kmadd(gtu23,kadd(JacPDstandardNth31beta1,kmadd(JacPDstandardNth32beta2,ToReal(4),JacPDstandardNth33beta3)),kmul(alphaL,kmadd(kmadd(gtu12,S1,kmadd(gtu22,S2,kmul(gtu23,S3))),ToReal(-150.7964473723100754462068823974161384415),kmadd(kmadd(gtu12,JacPDstandardNth1trK,kmadd(gtu22,JacPDstandardNth2trK,kmul(gtu23,JacPDstandardNth3trK))),ToReal(-4),kmadd(ToReal(6),kmadd(Atu11,Gt211,kmadd(Atu33,Gt233,kmul(Atu22,kmadd(cdphi2,ToReal(6),Gt222)))),kmadd(kmadd(Atu13,Gt213,kmul(Atu12,kmadd(cdphi1,ToReal(3),Gt212))),ToReal(12),kmul(Atu23,kmadd(Gt223,ToReal(12),kmul(cdphi3,ToReal(36)))))))))))))))));
    
    CCTK_REAL_VEC dotXt3 = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(kmadd(Atu13,JacPDstandardNth1alpha,kmadd(Atu23,JacPDstandardNth2alpha,kmul(Atu33,JacPDstandardNth3alpha))),ToReal(-6),kmadd(kmadd(JacPDstandardNth1beta3,Xtn1,kmul(JacPDstandardNth2beta3,Xtn2)),ToReal(-3),kmadd(Xtn3,kmsub(kadd(JacPDstandardNth1beta1,JacPDstandardNth2beta2),ToReal(2),JacPDstandardNth3beta3),kmadd(kmadd(gtu11,JacPDstandardNth11beta3,kmadd(gtu12,kadd(JacPDstandardNth12beta3,JacPDstandardNth21beta3),kmadd(gtu22,JacPDstandardNth22beta3,kmul(gtu23,JacPDstandardNth32beta3)))),ToReal(3),kmadd(gtu13,kadd(JacPDstandardNth11beta1,kadd(JacPDstandardNth12beta2,kmadd(JacPDstandardNth31beta3,ToReal(3),kmul(JacPDstandardNth13beta3,ToReal(4))))),kmadd(gtu23,kadd(JacPDstandardNth21beta1,kmadd(JacPDstandardNth23beta3,ToReal(4),JacPDstandardNth22beta2)),kmadd(gtu33,kadd(JacPDstandardNth31beta1,kmadd(JacPDstandardNth33beta3,ToReal(4),JacPDstandardNth32beta2)),kmul(alphaL,kmadd(kmadd(gtu13,S1,kmadd(gtu23,S2,kmul(gtu33,S3))),ToReal(-150.7964473723100754462068823974161384415),kmadd(kmadd(gtu13,JacPDstandardNth1trK,kmadd(gtu23,JacPDstandardNth2trK,kmul(gtu33,JacPDstandardNth3trK))),ToReal(-4),kmadd(kmadd(Atu11,Gt311,kmul(Atu22,Gt322)),ToReal(6),kmadd(kmadd(Atu12,Gt312,kmadd(Atu13,kmadd(cdphi1,ToReal(3),Gt313),kmul(Atu23,kmadd(cdphi2,ToReal(3),Gt323)))),ToReal(12),kmul(Atu33,kmadd(Gt333,ToReal(6),kmul(cdphi3,ToReal(36)))))))))))))))));
    
    CCTK_REAL_VEC Xt1rhsL = dotXt1;
    
    CCTK_REAL_VEC Xt2rhsL = dotXt2;
    
    CCTK_REAL_VEC Xt3rhsL = dotXt3;
    
    CCTK_REAL_VEC dottrK = 
      kmsub(alphaL,kadd(SQR(Atm11),kadd(SQR(Atm22),kadd(SQR(Atm33),kmadd(SQR(trKL),ToReal(0.333333333333333333333333333333),kmadd(kmadd(Atm12,Atm21,kmadd(Atm13,Atm31,kmul(Atm23,Atm32))),ToReal(2),kmul(kadd(rho,trS),ToReal(12.56637061435917295385057353311801153679))))))),kmul(em4phi,knmsub(JacPDstandardNth1alpha,Xtn1,knmsub(JacPDstandardNth2alpha,Xtn2,knmsub(JacPDstandardNth3alpha,Xtn3,kmadd(gtu11,kmadd(cdphi1,kmul(JacPDstandardNth1alpha,ToReal(2)),JacPDstandardNth11alpha),kmadd(gtu12,kadd(JacPDstandardNth12alpha,kadd(JacPDstandardNth21alpha,kmadd(cdphi2,kmul(JacPDstandardNth1alpha,ToReal(2)),kmul(cdphi1,kmul(JacPDstandardNth2alpha,ToReal(2)))))),kmadd(gtu22,kmadd(cdphi2,kmul(JacPDstandardNth2alpha,ToReal(2)),JacPDstandardNth22alpha),kmadd(gtu13,kadd(JacPDstandardNth13alpha,kadd(JacPDstandardNth31alpha,kmadd(cdphi3,kmul(JacPDstandardNth1alpha,ToReal(2)),kmul(cdphi1,kmul(JacPDstandardNth3alpha,ToReal(2)))))),kmadd(gtu23,kadd(JacPDstandardNth23alpha,kadd(JacPDstandardNth32alpha,kmadd(cdphi3,kmul(JacPDstandardNth2alpha,ToReal(2)),kmul(cdphi2,kmul(JacPDstandardNth3alpha,ToReal(2)))))),kmul(gtu33,kmadd(cdphi3,kmul(JacPDstandardNth3alpha,ToReal(2)),JacPDstandardNth33alpha))))))))))));
    
    CCTK_REAL_VEC trKrhsL = dottrK;
    
    CCTK_REAL_VEC alpharhsL = 
      kneg(kmul(kpow(alphaL,harmonicN),kmul(ToReal(harmonicF),kmadd(ksub(AL,trKL),ToReal(LapseACoeff),trKL))));
    
    CCTK_REAL_VEC ArhsL = 
      kmul(knmsub(AL,ToReal(AlphaDriver),dottrK),ToReal(LapseACoeff));
    
    CCTK_REAL_VEC eta = 
      kfmin(ToReal(1),kmul(INV(rL),ToReal(SpatialBetaDriverRadius)));
    
    CCTK_REAL_VEC theta = 
      kfmin(ToReal(1),kexp(knmsub(rL,INV(ToReal(SpatialShiftGammaCoeffRadius)),ToReal(1))));
    
    CCTK_REAL_VEC beta1rhsL;
    CCTK_REAL_VEC beta2rhsL;
    CCTK_REAL_VEC beta3rhsL;
    
    if (harmonicShift)
    {
      beta1rhsL = 
        kmul(alphaL,kmul(phiL,kmul(ToReal(0.5),kmadd(kmadd(gtu11,JacPDstandardNth1alpha,kmadd(gtu12,JacPDstandardNth2alpha,kmul(gtu13,JacPDstandardNth3alpha))),kmul(phiL,ToReal(-2)),kmul(alphaL,kmadd(phiL,kmadd(JacPDstandardNth1gt11,SQR(gtu11),kmul(JacPDstandardNth1gt22,kmul(SQR(gtu12),ToReal(2)))),kmadd(gtu13,kmadd(JacPDstandardNth3phi,ToReal(2),kmul(phiL,kmadd(gtu33,JacPDstandardNth3gt33,kmsub(kmadd(gtu13,JacPDstandardNth1gt33,kmadd(gtu22,JacPDstandardNth2gt23,kmul(gtu23,JacPDstandardNth2gt33))),ToReal(2),kmul(gtu22,JacPDstandardNth3gt22))))),kmadd(gtu11,kmadd(JacPDstandardNth1phi,ToReal(2),kmul(phiL,kmadd(gtu12,JacPDstandardNth2gt11,kmadd(gtu13,JacPDstandardNth3gt11,kmadd(gtu23,kmul(JacPDstandardNth1gt23,ToReal(-2)),knmsub(gtu22,JacPDstandardNth1gt22,kmadd(kmadd(gtu12,JacPDstandardNth1gt12,kmadd(gtu13,JacPDstandardNth1gt13,kmul(gtu22,JacPDstandardNth2gt12))),ToReal(2),kmadd(gtu23,kmul(JacPDstandardNth2gt13,ToReal(2)),kmadd(gtu23,kmul(JacPDstandardNth3gt12,ToReal(2)),kmul(gtu33,kmsub(JacPDstandardNth3gt13,ToReal(2),JacPDstandardNth1gt33))))))))))),kmul(gtu12,kmadd(JacPDstandardNth2phi,ToReal(2),kmul(phiL,kmadd(gtu22,JacPDstandardNth2gt22,kmadd(gtu23,kmul(JacPDstandardNth3gt22,ToReal(2)),kmadd(gtu33,kmsub(JacPDstandardNth3gt23,ToReal(2),JacPDstandardNth2gt33),kmul(gtu13,kmul(JacPDstandardNth1gt23,ToReal(4)))))))))))))))));
      
      beta2rhsL = 
        kmul(alphaL,kmul(phiL,kmul(ToReal(0.5),kmadd(kmadd(gtu12,JacPDstandardNth1alpha,kmadd(gtu22,JacPDstandardNth2alpha,kmul(gtu23,JacPDstandardNth3alpha))),kmul(phiL,ToReal(-2)),kmul(alphaL,kmadd(phiL,kmadd(JacPDstandardNth2gt22,SQR(gtu22),kmul(JacPDstandardNth2gt11,kmul(SQR(gtu12),ToReal(2)))),kmadd(gtu23,kmadd(JacPDstandardNth3phi,ToReal(2),kmul(phiL,kmadd(gtu33,JacPDstandardNth3gt33,kmsub(kmadd(gtu11,JacPDstandardNth1gt13,kmadd(gtu13,JacPDstandardNth1gt33,kmul(gtu23,JacPDstandardNth2gt33))),ToReal(2),kmul(gtu11,JacPDstandardNth3gt11))))),kmadd(gtu22,kmadd(JacPDstandardNth2phi,ToReal(2),kmul(phiL,kmadd(gtu23,JacPDstandardNth3gt22,kmadd(kmadd(gtu23,JacPDstandardNth2gt23,kmul(gtu13,kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth3gt12,JacPDstandardNth2gt13)))),ToReal(2),kmadd(gtu11,kmsub(JacPDstandardNth1gt12,ToReal(2),JacPDstandardNth2gt11),kmul(gtu33,kmsub(JacPDstandardNth3gt23,ToReal(2),JacPDstandardNth2gt33))))))),kmul(gtu12,kmadd(JacPDstandardNth1phi,ToReal(2),kmul(phiL,kmadd(gtu11,JacPDstandardNth1gt11,kmadd(gtu13,kmul(JacPDstandardNth3gt11,ToReal(2)),kmadd(gtu22,kmadd(JacPDstandardNth2gt12,ToReal(2),JacPDstandardNth1gt22),kmadd(gtu33,kmsub(JacPDstandardNth3gt13,ToReal(2),JacPDstandardNth1gt33),kmul(gtu23,kmul(JacPDstandardNth2gt13,ToReal(4))))))))))))))))));
      
      beta3rhsL = 
        kmul(alphaL,kmul(phiL,kmul(ToReal(0.5),kmadd(kmadd(gtu13,JacPDstandardNth1alpha,kmadd(gtu23,JacPDstandardNth2alpha,kmul(gtu33,JacPDstandardNth3alpha))),kmul(phiL,ToReal(-2)),kmul(alphaL,kmadd(phiL,kmul(kmadd(JacPDstandardNth3gt11,SQR(gtu13),kmul(JacPDstandardNth3gt22,SQR(gtu23))),ToReal(2)),kmadd(gtu23,kmadd(JacPDstandardNth2phi,ToReal(2),kmul(phiL,kmadd(gtu22,JacPDstandardNth2gt22,kmadd(gtu33,JacPDstandardNth2gt33,kmsub(kmadd(gtu11,JacPDstandardNth1gt12,kmadd(gtu12,JacPDstandardNth1gt22,kmul(gtu33,JacPDstandardNth3gt23))),ToReal(2),kmul(gtu11,JacPDstandardNth2gt11)))))),kmadd(gtu33,kmadd(JacPDstandardNth3phi,ToReal(2),kmul(phiL,kmadd(gtu33,JacPDstandardNth3gt33,knmsub(gtu22,JacPDstandardNth3gt22,kmadd(kmadd(gtu22,JacPDstandardNth2gt23,kmul(gtu12,kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth2gt13,JacPDstandardNth3gt12)))),ToReal(2),kmul(gtu11,kmsub(JacPDstandardNth1gt13,ToReal(2),JacPDstandardNth3gt11))))))),kmul(gtu13,kmadd(JacPDstandardNth1phi,ToReal(2),kmul(phiL,kmadd(gtu11,JacPDstandardNth1gt11,kmadd(gtu12,kmul(JacPDstandardNth2gt11,ToReal(2)),kmadd(gtu22,kmsub(JacPDstandardNth2gt12,ToReal(2),JacPDstandardNth1gt22),kmadd(gtu33,kmadd(JacPDstandardNth3gt13,ToReal(2),JacPDstandardNth1gt33),kmul(gtu23,kmul(JacPDstandardNth3gt12,ToReal(4))))))))))))))))));
    }
    else
    {
      beta1rhsL = 
        kmul(theta,kmul(kadd(Xt1L,kmadd(beta1L,kmul(eta,ToReal(BetaDriver*(-1 + 
        ShiftBCoeff))),kmul(ksub(B1L,Xt1L),ToReal(ShiftBCoeff)))),ToReal(ShiftGammaCoeff)));
      
      beta2rhsL = 
        kmul(theta,kmul(kadd(Xt2L,kmadd(beta2L,kmul(eta,ToReal(BetaDriver*(-1 + 
        ShiftBCoeff))),kmul(ksub(B2L,Xt2L),ToReal(ShiftBCoeff)))),ToReal(ShiftGammaCoeff)));
      
      beta3rhsL = 
        kmul(theta,kmul(kadd(Xt3L,kmadd(beta3L,kmul(eta,ToReal(BetaDriver*(-1 + 
        ShiftBCoeff))),kmul(ksub(B3L,Xt3L),ToReal(ShiftBCoeff)))),ToReal(ShiftGammaCoeff)));
    }
    
    CCTK_REAL_VEC B1rhsL = 
      kmul(knmsub(B1L,kmul(eta,ToReal(BetaDriver)),dotXt1),ToReal(ShiftBCoeff));
    
    CCTK_REAL_VEC B2rhsL = 
      kmul(knmsub(B2L,kmul(eta,ToReal(BetaDriver)),dotXt2),ToReal(ShiftBCoeff));
    
    CCTK_REAL_VEC B3rhsL = 
      kmul(knmsub(B3L,kmul(eta,ToReal(BetaDriver)),dotXt3),ToReal(ShiftBCoeff));
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(alpharhs[index],alpharhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Arhs[index],ArhsL,elt_count_lo,elt_count_hi);
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
  LC_ENDLOOP3VEC (ML_BSSN_UPW_RHS1);
}

extern "C" void ML_BSSN_UPW_RHS1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_UPW_RHS1_Body");
  }
  
  if (cctk_iteration % ML_BSSN_UPW_RHS1_calc_every != ML_BSSN_UPW_RHS1_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"grid::coordinates","Grid::coordinates","ML_BSSN_UPW::ML_curv","ML_BSSN_UPW::ML_dtlapse","ML_BSSN_UPW::ML_dtlapserhs","ML_BSSN_UPW::ML_dtshift","ML_BSSN_UPW::ML_dtshiftrhs","ML_BSSN_UPW::ML_Gamma","ML_BSSN_UPW::ML_Gammarhs","ML_BSSN_UPW::ML_lapse","ML_BSSN_UPW::ML_lapserhs","ML_BSSN_UPW::ML_log_confac","ML_BSSN_UPW::ML_log_confacrhs","ML_BSSN_UPW::ML_metric","ML_BSSN_UPW::ML_metricrhs","ML_BSSN_UPW::ML_shift","ML_BSSN_UPW::ML_shiftrhs","ML_BSSN_UPW::ML_trace_curv","ML_BSSN_UPW::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_UPW_RHS1", 19, groups);
  
  switch(fdOrder)
  {
    case 2:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_RHS1", 1, 1, 1);
      break;
    
    case 4:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_RHS1", 2, 2, 2);
      break;
    
    case 6:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_RHS1", 3, 3, 3);
      break;
    
    case 8:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_RHS1", 4, 4, 4);
      break;
  }
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_UPW_RHS1_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_UPW_RHS1_Body");
  }
}
