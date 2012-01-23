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

extern "C" void ML_BSSN_UPW_RHS2_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_curvrhs.");
  return;
}

static void ML_BSSN_UPW_RHS2_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  LC_LOOP3VEC (ML_BSSN_UPW_RHS2,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC At11L = vec_load(At11[index]);
    CCTK_REAL_VEC At12L = vec_load(At12[index]);
    CCTK_REAL_VEC At13L = vec_load(At13[index]);
    CCTK_REAL_VEC At22L = vec_load(At22[index]);
    CCTK_REAL_VEC At23L = vec_load(At23[index]);
    CCTK_REAL_VEC At33L = vec_load(At33[index]);
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
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    CCTK_REAL_VEC Xt1L = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L = vec_load(Xt3[index]);
    
    CCTK_REAL_VEC eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL;
    
    if (*stress_energy_state)
    {
      eTxxL = vec_load(eTxx[index]);
      eTxyL = vec_load(eTxy[index]);
      eTxzL = vec_load(eTxz[index]);
      eTyyL = vec_load(eTyy[index]);
      eTyzL = vec_load(eTyz[index]);
      eTzzL = vec_load(eTzz[index]);
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
    CCTK_REAL_VEC PDstandardNth1beta2;
    CCTK_REAL_VEC PDstandardNth2beta2;
    CCTK_REAL_VEC PDstandardNth3beta2;
    CCTK_REAL_VEC PDstandardNth1beta3;
    CCTK_REAL_VEC PDstandardNth2beta3;
    CCTK_REAL_VEC PDstandardNth3beta3;
    CCTK_REAL_VEC PDstandardNth1gt11;
    CCTK_REAL_VEC PDstandardNth2gt11;
    CCTK_REAL_VEC PDstandardNth3gt11;
    CCTK_REAL_VEC PDstandardNth11gt11;
    CCTK_REAL_VEC PDstandardNth22gt11;
    CCTK_REAL_VEC PDstandardNth33gt11;
    CCTK_REAL_VEC PDstandardNth12gt11;
    CCTK_REAL_VEC PDstandardNth13gt11;
    CCTK_REAL_VEC PDstandardNth23gt11;
    CCTK_REAL_VEC PDstandardNth1gt12;
    CCTK_REAL_VEC PDstandardNth2gt12;
    CCTK_REAL_VEC PDstandardNth3gt12;
    CCTK_REAL_VEC PDstandardNth11gt12;
    CCTK_REAL_VEC PDstandardNth22gt12;
    CCTK_REAL_VEC PDstandardNth33gt12;
    CCTK_REAL_VEC PDstandardNth12gt12;
    CCTK_REAL_VEC PDstandardNth13gt12;
    CCTK_REAL_VEC PDstandardNth23gt12;
    CCTK_REAL_VEC PDstandardNth1gt13;
    CCTK_REAL_VEC PDstandardNth2gt13;
    CCTK_REAL_VEC PDstandardNth3gt13;
    CCTK_REAL_VEC PDstandardNth11gt13;
    CCTK_REAL_VEC PDstandardNth22gt13;
    CCTK_REAL_VEC PDstandardNth33gt13;
    CCTK_REAL_VEC PDstandardNth12gt13;
    CCTK_REAL_VEC PDstandardNth13gt13;
    CCTK_REAL_VEC PDstandardNth23gt13;
    CCTK_REAL_VEC PDstandardNth1gt22;
    CCTK_REAL_VEC PDstandardNth2gt22;
    CCTK_REAL_VEC PDstandardNth3gt22;
    CCTK_REAL_VEC PDstandardNth11gt22;
    CCTK_REAL_VEC PDstandardNth22gt22;
    CCTK_REAL_VEC PDstandardNth33gt22;
    CCTK_REAL_VEC PDstandardNth12gt22;
    CCTK_REAL_VEC PDstandardNth13gt22;
    CCTK_REAL_VEC PDstandardNth23gt22;
    CCTK_REAL_VEC PDstandardNth1gt23;
    CCTK_REAL_VEC PDstandardNth2gt23;
    CCTK_REAL_VEC PDstandardNth3gt23;
    CCTK_REAL_VEC PDstandardNth11gt23;
    CCTK_REAL_VEC PDstandardNth22gt23;
    CCTK_REAL_VEC PDstandardNth33gt23;
    CCTK_REAL_VEC PDstandardNth12gt23;
    CCTK_REAL_VEC PDstandardNth13gt23;
    CCTK_REAL_VEC PDstandardNth23gt23;
    CCTK_REAL_VEC PDstandardNth1gt33;
    CCTK_REAL_VEC PDstandardNth2gt33;
    CCTK_REAL_VEC PDstandardNth3gt33;
    CCTK_REAL_VEC PDstandardNth11gt33;
    CCTK_REAL_VEC PDstandardNth22gt33;
    CCTK_REAL_VEC PDstandardNth33gt33;
    CCTK_REAL_VEC PDstandardNth12gt33;
    CCTK_REAL_VEC PDstandardNth13gt33;
    CCTK_REAL_VEC PDstandardNth23gt33;
    CCTK_REAL_VEC PDstandardNth1phi;
    CCTK_REAL_VEC PDstandardNth2phi;
    CCTK_REAL_VEC PDstandardNth3phi;
    CCTK_REAL_VEC PDstandardNth11phi;
    CCTK_REAL_VEC PDstandardNth22phi;
    CCTK_REAL_VEC PDstandardNth33phi;
    CCTK_REAL_VEC PDstandardNth12phi;
    CCTK_REAL_VEC PDstandardNth13phi;
    CCTK_REAL_VEC PDstandardNth23phi;
    CCTK_REAL_VEC PDstandardNth1Xt1;
    CCTK_REAL_VEC PDstandardNth2Xt1;
    CCTK_REAL_VEC PDstandardNth3Xt1;
    CCTK_REAL_VEC PDstandardNth1Xt2;
    CCTK_REAL_VEC PDstandardNth2Xt2;
    CCTK_REAL_VEC PDstandardNth3Xt2;
    CCTK_REAL_VEC PDstandardNth1Xt3;
    CCTK_REAL_VEC PDstandardNth2Xt3;
    CCTK_REAL_VEC PDstandardNth3Xt3;
    
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
        PDstandardNth1beta2 = PDstandardNthfdOrder21(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder22(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder23(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder21(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder22(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder23(&beta3[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder21(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder22(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder23(&gt11[index]);
        PDstandardNth11gt11 = PDstandardNthfdOrder211(&gt11[index]);
        PDstandardNth22gt11 = PDstandardNthfdOrder222(&gt11[index]);
        PDstandardNth33gt11 = PDstandardNthfdOrder233(&gt11[index]);
        PDstandardNth12gt11 = PDstandardNthfdOrder212(&gt11[index]);
        PDstandardNth13gt11 = PDstandardNthfdOrder213(&gt11[index]);
        PDstandardNth23gt11 = PDstandardNthfdOrder223(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder21(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder22(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder23(&gt12[index]);
        PDstandardNth11gt12 = PDstandardNthfdOrder211(&gt12[index]);
        PDstandardNth22gt12 = PDstandardNthfdOrder222(&gt12[index]);
        PDstandardNth33gt12 = PDstandardNthfdOrder233(&gt12[index]);
        PDstandardNth12gt12 = PDstandardNthfdOrder212(&gt12[index]);
        PDstandardNth13gt12 = PDstandardNthfdOrder213(&gt12[index]);
        PDstandardNth23gt12 = PDstandardNthfdOrder223(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder21(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder22(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder23(&gt13[index]);
        PDstandardNth11gt13 = PDstandardNthfdOrder211(&gt13[index]);
        PDstandardNth22gt13 = PDstandardNthfdOrder222(&gt13[index]);
        PDstandardNth33gt13 = PDstandardNthfdOrder233(&gt13[index]);
        PDstandardNth12gt13 = PDstandardNthfdOrder212(&gt13[index]);
        PDstandardNth13gt13 = PDstandardNthfdOrder213(&gt13[index]);
        PDstandardNth23gt13 = PDstandardNthfdOrder223(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder21(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder22(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder23(&gt22[index]);
        PDstandardNth11gt22 = PDstandardNthfdOrder211(&gt22[index]);
        PDstandardNth22gt22 = PDstandardNthfdOrder222(&gt22[index]);
        PDstandardNth33gt22 = PDstandardNthfdOrder233(&gt22[index]);
        PDstandardNth12gt22 = PDstandardNthfdOrder212(&gt22[index]);
        PDstandardNth13gt22 = PDstandardNthfdOrder213(&gt22[index]);
        PDstandardNth23gt22 = PDstandardNthfdOrder223(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder21(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder22(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder23(&gt23[index]);
        PDstandardNth11gt23 = PDstandardNthfdOrder211(&gt23[index]);
        PDstandardNth22gt23 = PDstandardNthfdOrder222(&gt23[index]);
        PDstandardNth33gt23 = PDstandardNthfdOrder233(&gt23[index]);
        PDstandardNth12gt23 = PDstandardNthfdOrder212(&gt23[index]);
        PDstandardNth13gt23 = PDstandardNthfdOrder213(&gt23[index]);
        PDstandardNth23gt23 = PDstandardNthfdOrder223(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder21(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder22(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder23(&gt33[index]);
        PDstandardNth11gt33 = PDstandardNthfdOrder211(&gt33[index]);
        PDstandardNth22gt33 = PDstandardNthfdOrder222(&gt33[index]);
        PDstandardNth33gt33 = PDstandardNthfdOrder233(&gt33[index]);
        PDstandardNth12gt33 = PDstandardNthfdOrder212(&gt33[index]);
        PDstandardNth13gt33 = PDstandardNthfdOrder213(&gt33[index]);
        PDstandardNth23gt33 = PDstandardNthfdOrder223(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder21(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder22(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder23(&phi[index]);
        PDstandardNth11phi = PDstandardNthfdOrder211(&phi[index]);
        PDstandardNth22phi = PDstandardNthfdOrder222(&phi[index]);
        PDstandardNth33phi = PDstandardNthfdOrder233(&phi[index]);
        PDstandardNth12phi = PDstandardNthfdOrder212(&phi[index]);
        PDstandardNth13phi = PDstandardNthfdOrder213(&phi[index]);
        PDstandardNth23phi = PDstandardNthfdOrder223(&phi[index]);
        PDstandardNth1Xt1 = PDstandardNthfdOrder21(&Xt1[index]);
        PDstandardNth2Xt1 = PDstandardNthfdOrder22(&Xt1[index]);
        PDstandardNth3Xt1 = PDstandardNthfdOrder23(&Xt1[index]);
        PDstandardNth1Xt2 = PDstandardNthfdOrder21(&Xt2[index]);
        PDstandardNth2Xt2 = PDstandardNthfdOrder22(&Xt2[index]);
        PDstandardNth3Xt2 = PDstandardNthfdOrder23(&Xt2[index]);
        PDstandardNth1Xt3 = PDstandardNthfdOrder21(&Xt3[index]);
        PDstandardNth2Xt3 = PDstandardNthfdOrder22(&Xt3[index]);
        PDstandardNth3Xt3 = PDstandardNthfdOrder23(&Xt3[index]);
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
        PDstandardNth1beta2 = PDstandardNthfdOrder41(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder42(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder43(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder41(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder42(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder43(&beta3[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder41(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder42(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder43(&gt11[index]);
        PDstandardNth11gt11 = PDstandardNthfdOrder411(&gt11[index]);
        PDstandardNth22gt11 = PDstandardNthfdOrder422(&gt11[index]);
        PDstandardNth33gt11 = PDstandardNthfdOrder433(&gt11[index]);
        PDstandardNth12gt11 = PDstandardNthfdOrder412(&gt11[index]);
        PDstandardNth13gt11 = PDstandardNthfdOrder413(&gt11[index]);
        PDstandardNth23gt11 = PDstandardNthfdOrder423(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder41(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder42(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder43(&gt12[index]);
        PDstandardNth11gt12 = PDstandardNthfdOrder411(&gt12[index]);
        PDstandardNth22gt12 = PDstandardNthfdOrder422(&gt12[index]);
        PDstandardNth33gt12 = PDstandardNthfdOrder433(&gt12[index]);
        PDstandardNth12gt12 = PDstandardNthfdOrder412(&gt12[index]);
        PDstandardNth13gt12 = PDstandardNthfdOrder413(&gt12[index]);
        PDstandardNth23gt12 = PDstandardNthfdOrder423(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder41(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder42(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder43(&gt13[index]);
        PDstandardNth11gt13 = PDstandardNthfdOrder411(&gt13[index]);
        PDstandardNth22gt13 = PDstandardNthfdOrder422(&gt13[index]);
        PDstandardNth33gt13 = PDstandardNthfdOrder433(&gt13[index]);
        PDstandardNth12gt13 = PDstandardNthfdOrder412(&gt13[index]);
        PDstandardNth13gt13 = PDstandardNthfdOrder413(&gt13[index]);
        PDstandardNth23gt13 = PDstandardNthfdOrder423(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder41(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder42(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder43(&gt22[index]);
        PDstandardNth11gt22 = PDstandardNthfdOrder411(&gt22[index]);
        PDstandardNth22gt22 = PDstandardNthfdOrder422(&gt22[index]);
        PDstandardNth33gt22 = PDstandardNthfdOrder433(&gt22[index]);
        PDstandardNth12gt22 = PDstandardNthfdOrder412(&gt22[index]);
        PDstandardNth13gt22 = PDstandardNthfdOrder413(&gt22[index]);
        PDstandardNth23gt22 = PDstandardNthfdOrder423(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder41(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder42(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder43(&gt23[index]);
        PDstandardNth11gt23 = PDstandardNthfdOrder411(&gt23[index]);
        PDstandardNth22gt23 = PDstandardNthfdOrder422(&gt23[index]);
        PDstandardNth33gt23 = PDstandardNthfdOrder433(&gt23[index]);
        PDstandardNth12gt23 = PDstandardNthfdOrder412(&gt23[index]);
        PDstandardNth13gt23 = PDstandardNthfdOrder413(&gt23[index]);
        PDstandardNth23gt23 = PDstandardNthfdOrder423(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder41(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder42(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder43(&gt33[index]);
        PDstandardNth11gt33 = PDstandardNthfdOrder411(&gt33[index]);
        PDstandardNth22gt33 = PDstandardNthfdOrder422(&gt33[index]);
        PDstandardNth33gt33 = PDstandardNthfdOrder433(&gt33[index]);
        PDstandardNth12gt33 = PDstandardNthfdOrder412(&gt33[index]);
        PDstandardNth13gt33 = PDstandardNthfdOrder413(&gt33[index]);
        PDstandardNth23gt33 = PDstandardNthfdOrder423(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder41(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder42(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder43(&phi[index]);
        PDstandardNth11phi = PDstandardNthfdOrder411(&phi[index]);
        PDstandardNth22phi = PDstandardNthfdOrder422(&phi[index]);
        PDstandardNth33phi = PDstandardNthfdOrder433(&phi[index]);
        PDstandardNth12phi = PDstandardNthfdOrder412(&phi[index]);
        PDstandardNth13phi = PDstandardNthfdOrder413(&phi[index]);
        PDstandardNth23phi = PDstandardNthfdOrder423(&phi[index]);
        PDstandardNth1Xt1 = PDstandardNthfdOrder41(&Xt1[index]);
        PDstandardNth2Xt1 = PDstandardNthfdOrder42(&Xt1[index]);
        PDstandardNth3Xt1 = PDstandardNthfdOrder43(&Xt1[index]);
        PDstandardNth1Xt2 = PDstandardNthfdOrder41(&Xt2[index]);
        PDstandardNth2Xt2 = PDstandardNthfdOrder42(&Xt2[index]);
        PDstandardNth3Xt2 = PDstandardNthfdOrder43(&Xt2[index]);
        PDstandardNth1Xt3 = PDstandardNthfdOrder41(&Xt3[index]);
        PDstandardNth2Xt3 = PDstandardNthfdOrder42(&Xt3[index]);
        PDstandardNth3Xt3 = PDstandardNthfdOrder43(&Xt3[index]);
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
        PDstandardNth1beta2 = PDstandardNthfdOrder61(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder62(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder63(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder61(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder62(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder63(&beta3[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder61(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder62(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder63(&gt11[index]);
        PDstandardNth11gt11 = PDstandardNthfdOrder611(&gt11[index]);
        PDstandardNth22gt11 = PDstandardNthfdOrder622(&gt11[index]);
        PDstandardNth33gt11 = PDstandardNthfdOrder633(&gt11[index]);
        PDstandardNth12gt11 = PDstandardNthfdOrder612(&gt11[index]);
        PDstandardNth13gt11 = PDstandardNthfdOrder613(&gt11[index]);
        PDstandardNth23gt11 = PDstandardNthfdOrder623(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder61(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder62(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder63(&gt12[index]);
        PDstandardNth11gt12 = PDstandardNthfdOrder611(&gt12[index]);
        PDstandardNth22gt12 = PDstandardNthfdOrder622(&gt12[index]);
        PDstandardNth33gt12 = PDstandardNthfdOrder633(&gt12[index]);
        PDstandardNth12gt12 = PDstandardNthfdOrder612(&gt12[index]);
        PDstandardNth13gt12 = PDstandardNthfdOrder613(&gt12[index]);
        PDstandardNth23gt12 = PDstandardNthfdOrder623(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder61(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder62(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder63(&gt13[index]);
        PDstandardNth11gt13 = PDstandardNthfdOrder611(&gt13[index]);
        PDstandardNth22gt13 = PDstandardNthfdOrder622(&gt13[index]);
        PDstandardNth33gt13 = PDstandardNthfdOrder633(&gt13[index]);
        PDstandardNth12gt13 = PDstandardNthfdOrder612(&gt13[index]);
        PDstandardNth13gt13 = PDstandardNthfdOrder613(&gt13[index]);
        PDstandardNth23gt13 = PDstandardNthfdOrder623(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder61(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder62(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder63(&gt22[index]);
        PDstandardNth11gt22 = PDstandardNthfdOrder611(&gt22[index]);
        PDstandardNth22gt22 = PDstandardNthfdOrder622(&gt22[index]);
        PDstandardNth33gt22 = PDstandardNthfdOrder633(&gt22[index]);
        PDstandardNth12gt22 = PDstandardNthfdOrder612(&gt22[index]);
        PDstandardNth13gt22 = PDstandardNthfdOrder613(&gt22[index]);
        PDstandardNth23gt22 = PDstandardNthfdOrder623(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder61(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder62(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder63(&gt23[index]);
        PDstandardNth11gt23 = PDstandardNthfdOrder611(&gt23[index]);
        PDstandardNth22gt23 = PDstandardNthfdOrder622(&gt23[index]);
        PDstandardNth33gt23 = PDstandardNthfdOrder633(&gt23[index]);
        PDstandardNth12gt23 = PDstandardNthfdOrder612(&gt23[index]);
        PDstandardNth13gt23 = PDstandardNthfdOrder613(&gt23[index]);
        PDstandardNth23gt23 = PDstandardNthfdOrder623(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder61(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder62(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder63(&gt33[index]);
        PDstandardNth11gt33 = PDstandardNthfdOrder611(&gt33[index]);
        PDstandardNth22gt33 = PDstandardNthfdOrder622(&gt33[index]);
        PDstandardNth33gt33 = PDstandardNthfdOrder633(&gt33[index]);
        PDstandardNth12gt33 = PDstandardNthfdOrder612(&gt33[index]);
        PDstandardNth13gt33 = PDstandardNthfdOrder613(&gt33[index]);
        PDstandardNth23gt33 = PDstandardNthfdOrder623(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder61(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder62(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder63(&phi[index]);
        PDstandardNth11phi = PDstandardNthfdOrder611(&phi[index]);
        PDstandardNth22phi = PDstandardNthfdOrder622(&phi[index]);
        PDstandardNth33phi = PDstandardNthfdOrder633(&phi[index]);
        PDstandardNth12phi = PDstandardNthfdOrder612(&phi[index]);
        PDstandardNth13phi = PDstandardNthfdOrder613(&phi[index]);
        PDstandardNth23phi = PDstandardNthfdOrder623(&phi[index]);
        PDstandardNth1Xt1 = PDstandardNthfdOrder61(&Xt1[index]);
        PDstandardNth2Xt1 = PDstandardNthfdOrder62(&Xt1[index]);
        PDstandardNth3Xt1 = PDstandardNthfdOrder63(&Xt1[index]);
        PDstandardNth1Xt2 = PDstandardNthfdOrder61(&Xt2[index]);
        PDstandardNth2Xt2 = PDstandardNthfdOrder62(&Xt2[index]);
        PDstandardNth3Xt2 = PDstandardNthfdOrder63(&Xt2[index]);
        PDstandardNth1Xt3 = PDstandardNthfdOrder61(&Xt3[index]);
        PDstandardNth2Xt3 = PDstandardNthfdOrder62(&Xt3[index]);
        PDstandardNth3Xt3 = PDstandardNthfdOrder63(&Xt3[index]);
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
        PDstandardNth1beta2 = PDstandardNthfdOrder81(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder82(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder83(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder81(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder82(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder83(&beta3[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder81(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder82(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder83(&gt11[index]);
        PDstandardNth11gt11 = PDstandardNthfdOrder811(&gt11[index]);
        PDstandardNth22gt11 = PDstandardNthfdOrder822(&gt11[index]);
        PDstandardNth33gt11 = PDstandardNthfdOrder833(&gt11[index]);
        PDstandardNth12gt11 = PDstandardNthfdOrder812(&gt11[index]);
        PDstandardNth13gt11 = PDstandardNthfdOrder813(&gt11[index]);
        PDstandardNth23gt11 = PDstandardNthfdOrder823(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder81(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder82(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder83(&gt12[index]);
        PDstandardNth11gt12 = PDstandardNthfdOrder811(&gt12[index]);
        PDstandardNth22gt12 = PDstandardNthfdOrder822(&gt12[index]);
        PDstandardNth33gt12 = PDstandardNthfdOrder833(&gt12[index]);
        PDstandardNth12gt12 = PDstandardNthfdOrder812(&gt12[index]);
        PDstandardNth13gt12 = PDstandardNthfdOrder813(&gt12[index]);
        PDstandardNth23gt12 = PDstandardNthfdOrder823(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder81(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder82(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder83(&gt13[index]);
        PDstandardNth11gt13 = PDstandardNthfdOrder811(&gt13[index]);
        PDstandardNth22gt13 = PDstandardNthfdOrder822(&gt13[index]);
        PDstandardNth33gt13 = PDstandardNthfdOrder833(&gt13[index]);
        PDstandardNth12gt13 = PDstandardNthfdOrder812(&gt13[index]);
        PDstandardNth13gt13 = PDstandardNthfdOrder813(&gt13[index]);
        PDstandardNth23gt13 = PDstandardNthfdOrder823(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder81(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder82(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder83(&gt22[index]);
        PDstandardNth11gt22 = PDstandardNthfdOrder811(&gt22[index]);
        PDstandardNth22gt22 = PDstandardNthfdOrder822(&gt22[index]);
        PDstandardNth33gt22 = PDstandardNthfdOrder833(&gt22[index]);
        PDstandardNth12gt22 = PDstandardNthfdOrder812(&gt22[index]);
        PDstandardNth13gt22 = PDstandardNthfdOrder813(&gt22[index]);
        PDstandardNth23gt22 = PDstandardNthfdOrder823(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder81(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder82(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder83(&gt23[index]);
        PDstandardNth11gt23 = PDstandardNthfdOrder811(&gt23[index]);
        PDstandardNth22gt23 = PDstandardNthfdOrder822(&gt23[index]);
        PDstandardNth33gt23 = PDstandardNthfdOrder833(&gt23[index]);
        PDstandardNth12gt23 = PDstandardNthfdOrder812(&gt23[index]);
        PDstandardNth13gt23 = PDstandardNthfdOrder813(&gt23[index]);
        PDstandardNth23gt23 = PDstandardNthfdOrder823(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder81(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder82(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder83(&gt33[index]);
        PDstandardNth11gt33 = PDstandardNthfdOrder811(&gt33[index]);
        PDstandardNth22gt33 = PDstandardNthfdOrder822(&gt33[index]);
        PDstandardNth33gt33 = PDstandardNthfdOrder833(&gt33[index]);
        PDstandardNth12gt33 = PDstandardNthfdOrder812(&gt33[index]);
        PDstandardNth13gt33 = PDstandardNthfdOrder813(&gt33[index]);
        PDstandardNth23gt33 = PDstandardNthfdOrder823(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder81(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder82(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder83(&phi[index]);
        PDstandardNth11phi = PDstandardNthfdOrder811(&phi[index]);
        PDstandardNth22phi = PDstandardNthfdOrder822(&phi[index]);
        PDstandardNth33phi = PDstandardNthfdOrder833(&phi[index]);
        PDstandardNth12phi = PDstandardNthfdOrder812(&phi[index]);
        PDstandardNth13phi = PDstandardNthfdOrder813(&phi[index]);
        PDstandardNth23phi = PDstandardNthfdOrder823(&phi[index]);
        PDstandardNth1Xt1 = PDstandardNthfdOrder81(&Xt1[index]);
        PDstandardNth2Xt1 = PDstandardNthfdOrder82(&Xt1[index]);
        PDstandardNth3Xt1 = PDstandardNthfdOrder83(&Xt1[index]);
        PDstandardNth1Xt2 = PDstandardNthfdOrder81(&Xt2[index]);
        PDstandardNth2Xt2 = PDstandardNthfdOrder82(&Xt2[index]);
        PDstandardNth3Xt2 = PDstandardNthfdOrder83(&Xt2[index]);
        PDstandardNth1Xt3 = PDstandardNthfdOrder81(&Xt3[index]);
        PDstandardNth2Xt3 = PDstandardNthfdOrder82(&Xt3[index]);
        PDstandardNth3Xt3 = PDstandardNthfdOrder83(&Xt3[index]);
        break;
    }
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    CCTK_REAL_VEC JacPDstandardNth11alpha;
    CCTK_REAL_VEC JacPDstandardNth11gt11;
    CCTK_REAL_VEC JacPDstandardNth11gt12;
    CCTK_REAL_VEC JacPDstandardNth11gt13;
    CCTK_REAL_VEC JacPDstandardNth11gt22;
    CCTK_REAL_VEC JacPDstandardNth11gt23;
    CCTK_REAL_VEC JacPDstandardNth11gt33;
    CCTK_REAL_VEC JacPDstandardNth11phi;
    CCTK_REAL_VEC JacPDstandardNth12alpha;
    CCTK_REAL_VEC JacPDstandardNth12gt11;
    CCTK_REAL_VEC JacPDstandardNth12gt12;
    CCTK_REAL_VEC JacPDstandardNth12gt13;
    CCTK_REAL_VEC JacPDstandardNth12gt22;
    CCTK_REAL_VEC JacPDstandardNth12gt23;
    CCTK_REAL_VEC JacPDstandardNth12gt33;
    CCTK_REAL_VEC JacPDstandardNth12phi;
    CCTK_REAL_VEC JacPDstandardNth13alpha;
    CCTK_REAL_VEC JacPDstandardNth13gt11;
    CCTK_REAL_VEC JacPDstandardNth13gt12;
    CCTK_REAL_VEC JacPDstandardNth13gt13;
    CCTK_REAL_VEC JacPDstandardNth13gt22;
    CCTK_REAL_VEC JacPDstandardNth13gt23;
    CCTK_REAL_VEC JacPDstandardNth13gt33;
    CCTK_REAL_VEC JacPDstandardNth13phi;
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
    CCTK_REAL_VEC JacPDstandardNth1Xt1;
    CCTK_REAL_VEC JacPDstandardNth1Xt2;
    CCTK_REAL_VEC JacPDstandardNth1Xt3;
    CCTK_REAL_VEC JacPDstandardNth21gt11;
    CCTK_REAL_VEC JacPDstandardNth21gt12;
    CCTK_REAL_VEC JacPDstandardNth21gt13;
    CCTK_REAL_VEC JacPDstandardNth21gt22;
    CCTK_REAL_VEC JacPDstandardNth21gt23;
    CCTK_REAL_VEC JacPDstandardNth21gt33;
    CCTK_REAL_VEC JacPDstandardNth22alpha;
    CCTK_REAL_VEC JacPDstandardNth22gt11;
    CCTK_REAL_VEC JacPDstandardNth22gt12;
    CCTK_REAL_VEC JacPDstandardNth22gt13;
    CCTK_REAL_VEC JacPDstandardNth22gt22;
    CCTK_REAL_VEC JacPDstandardNth22gt23;
    CCTK_REAL_VEC JacPDstandardNth22gt33;
    CCTK_REAL_VEC JacPDstandardNth22phi;
    CCTK_REAL_VEC JacPDstandardNth23alpha;
    CCTK_REAL_VEC JacPDstandardNth23gt11;
    CCTK_REAL_VEC JacPDstandardNth23gt12;
    CCTK_REAL_VEC JacPDstandardNth23gt13;
    CCTK_REAL_VEC JacPDstandardNth23gt22;
    CCTK_REAL_VEC JacPDstandardNth23gt23;
    CCTK_REAL_VEC JacPDstandardNth23gt33;
    CCTK_REAL_VEC JacPDstandardNth23phi;
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
    CCTK_REAL_VEC JacPDstandardNth2Xt1;
    CCTK_REAL_VEC JacPDstandardNth2Xt2;
    CCTK_REAL_VEC JacPDstandardNth2Xt3;
    CCTK_REAL_VEC JacPDstandardNth31gt11;
    CCTK_REAL_VEC JacPDstandardNth31gt12;
    CCTK_REAL_VEC JacPDstandardNth31gt13;
    CCTK_REAL_VEC JacPDstandardNth31gt22;
    CCTK_REAL_VEC JacPDstandardNth31gt23;
    CCTK_REAL_VEC JacPDstandardNth31gt33;
    CCTK_REAL_VEC JacPDstandardNth32gt11;
    CCTK_REAL_VEC JacPDstandardNth32gt12;
    CCTK_REAL_VEC JacPDstandardNth32gt13;
    CCTK_REAL_VEC JacPDstandardNth32gt22;
    CCTK_REAL_VEC JacPDstandardNth32gt23;
    CCTK_REAL_VEC JacPDstandardNth32gt33;
    CCTK_REAL_VEC JacPDstandardNth33alpha;
    CCTK_REAL_VEC JacPDstandardNth33gt11;
    CCTK_REAL_VEC JacPDstandardNth33gt12;
    CCTK_REAL_VEC JacPDstandardNth33gt13;
    CCTK_REAL_VEC JacPDstandardNth33gt22;
    CCTK_REAL_VEC JacPDstandardNth33gt23;
    CCTK_REAL_VEC JacPDstandardNth33gt33;
    CCTK_REAL_VEC JacPDstandardNth33phi;
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
    CCTK_REAL_VEC JacPDstandardNth3Xt1;
    CCTK_REAL_VEC JacPDstandardNth3Xt2;
    CCTK_REAL_VEC JacPDstandardNth3Xt3;
    
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
      
      JacPDstandardNth1Xt1 = 
        kmadd(J11L,PDstandardNth1Xt1,kmadd(J21L,PDstandardNth2Xt1,kmul(J31L,PDstandardNth3Xt1)));
      
      JacPDstandardNth1Xt2 = 
        kmadd(J11L,PDstandardNth1Xt2,kmadd(J21L,PDstandardNth2Xt2,kmul(J31L,PDstandardNth3Xt2)));
      
      JacPDstandardNth1Xt3 = 
        kmadd(J11L,PDstandardNth1Xt3,kmadd(J21L,PDstandardNth2Xt3,kmul(J31L,PDstandardNth3Xt3)));
      
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
      
      JacPDstandardNth2Xt1 = 
        kmadd(J12L,PDstandardNth1Xt1,kmadd(J22L,PDstandardNth2Xt1,kmul(J32L,PDstandardNth3Xt1)));
      
      JacPDstandardNth2Xt2 = 
        kmadd(J12L,PDstandardNth1Xt2,kmadd(J22L,PDstandardNth2Xt2,kmul(J32L,PDstandardNth3Xt2)));
      
      JacPDstandardNth2Xt3 = 
        kmadd(J12L,PDstandardNth1Xt3,kmadd(J22L,PDstandardNth2Xt3,kmul(J32L,PDstandardNth3Xt3)));
      
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
      
      JacPDstandardNth3Xt1 = 
        kmadd(J13L,PDstandardNth1Xt1,kmadd(J23L,PDstandardNth2Xt1,kmul(J33L,PDstandardNth3Xt1)));
      
      JacPDstandardNth3Xt2 = 
        kmadd(J13L,PDstandardNth1Xt2,kmadd(J23L,PDstandardNth2Xt2,kmul(J33L,PDstandardNth3Xt2)));
      
      JacPDstandardNth3Xt3 = 
        kmadd(J13L,PDstandardNth1Xt3,kmadd(J23L,PDstandardNth2Xt3,kmul(J33L,PDstandardNth3Xt3)));
      
      JacPDstandardNth11alpha = 
        kmadd(dJ111L,PDstandardNth1alpha,kmadd(dJ211L,PDstandardNth2alpha,kmadd(dJ311L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,SQR(J11L),kmadd(PDstandardNth22alpha,SQR(J21L),kmadd(PDstandardNth33alpha,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha)),kmul(J21L,kmul(J31L,PDstandardNth23alpha))),ToReal(2))))))));
      
      JacPDstandardNth11gt11 = 
        kmadd(dJ111L,PDstandardNth1gt11,kmadd(dJ211L,PDstandardNth2gt11,kmadd(dJ311L,PDstandardNth3gt11,kmadd(PDstandardNth11gt11,SQR(J11L),kmadd(PDstandardNth22gt11,SQR(J21L),kmadd(PDstandardNth33gt11,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12gt11,kmul(J31L,PDstandardNth13gt11)),kmul(J21L,kmul(J31L,PDstandardNth23gt11))),ToReal(2))))))));
      
      JacPDstandardNth11gt12 = 
        kmadd(dJ111L,PDstandardNth1gt12,kmadd(dJ211L,PDstandardNth2gt12,kmadd(dJ311L,PDstandardNth3gt12,kmadd(PDstandardNth11gt12,SQR(J11L),kmadd(PDstandardNth22gt12,SQR(J21L),kmadd(PDstandardNth33gt12,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12gt12,kmul(J31L,PDstandardNth13gt12)),kmul(J21L,kmul(J31L,PDstandardNth23gt12))),ToReal(2))))))));
      
      JacPDstandardNth11gt13 = 
        kmadd(dJ111L,PDstandardNth1gt13,kmadd(dJ211L,PDstandardNth2gt13,kmadd(dJ311L,PDstandardNth3gt13,kmadd(PDstandardNth11gt13,SQR(J11L),kmadd(PDstandardNth22gt13,SQR(J21L),kmadd(PDstandardNth33gt13,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12gt13,kmul(J31L,PDstandardNth13gt13)),kmul(J21L,kmul(J31L,PDstandardNth23gt13))),ToReal(2))))))));
      
      JacPDstandardNth11gt22 = 
        kmadd(dJ111L,PDstandardNth1gt22,kmadd(dJ211L,PDstandardNth2gt22,kmadd(dJ311L,PDstandardNth3gt22,kmadd(PDstandardNth11gt22,SQR(J11L),kmadd(PDstandardNth22gt22,SQR(J21L),kmadd(PDstandardNth33gt22,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12gt22,kmul(J31L,PDstandardNth13gt22)),kmul(J21L,kmul(J31L,PDstandardNth23gt22))),ToReal(2))))))));
      
      JacPDstandardNth11gt23 = 
        kmadd(dJ111L,PDstandardNth1gt23,kmadd(dJ211L,PDstandardNth2gt23,kmadd(dJ311L,PDstandardNth3gt23,kmadd(PDstandardNth11gt23,SQR(J11L),kmadd(PDstandardNth22gt23,SQR(J21L),kmadd(PDstandardNth33gt23,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12gt23,kmul(J31L,PDstandardNth13gt23)),kmul(J21L,kmul(J31L,PDstandardNth23gt23))),ToReal(2))))))));
      
      JacPDstandardNth11gt33 = 
        kmadd(dJ111L,PDstandardNth1gt33,kmadd(dJ211L,PDstandardNth2gt33,kmadd(dJ311L,PDstandardNth3gt33,kmadd(PDstandardNth11gt33,SQR(J11L),kmadd(PDstandardNth22gt33,SQR(J21L),kmadd(PDstandardNth33gt33,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12gt33,kmul(J31L,PDstandardNth13gt33)),kmul(J21L,kmul(J31L,PDstandardNth23gt33))),ToReal(2))))))));
      
      JacPDstandardNth11phi = 
        kmadd(dJ111L,PDstandardNth1phi,kmadd(dJ211L,PDstandardNth2phi,kmadd(dJ311L,PDstandardNth3phi,kmadd(PDstandardNth11phi,SQR(J11L),kmadd(PDstandardNth22phi,SQR(J21L),kmadd(PDstandardNth33phi,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12phi,kmul(J31L,PDstandardNth13phi)),kmul(J21L,kmul(J31L,PDstandardNth23phi))),ToReal(2))))))));
      
      JacPDstandardNth22alpha = 
        kmadd(dJ122L,PDstandardNth1alpha,kmadd(dJ222L,PDstandardNth2alpha,kmadd(dJ322L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,SQR(J12L),kmadd(PDstandardNth22alpha,SQR(J22L),kmadd(PDstandardNth33alpha,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha)),kmul(J22L,kmul(J32L,PDstandardNth23alpha))),ToReal(2))))))));
      
      JacPDstandardNth22gt11 = 
        kmadd(dJ122L,PDstandardNth1gt11,kmadd(dJ222L,PDstandardNth2gt11,kmadd(dJ322L,PDstandardNth3gt11,kmadd(PDstandardNth11gt11,SQR(J12L),kmadd(PDstandardNth22gt11,SQR(J22L),kmadd(PDstandardNth33gt11,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12gt11,kmul(J32L,PDstandardNth13gt11)),kmul(J22L,kmul(J32L,PDstandardNth23gt11))),ToReal(2))))))));
      
      JacPDstandardNth22gt12 = 
        kmadd(dJ122L,PDstandardNth1gt12,kmadd(dJ222L,PDstandardNth2gt12,kmadd(dJ322L,PDstandardNth3gt12,kmadd(PDstandardNth11gt12,SQR(J12L),kmadd(PDstandardNth22gt12,SQR(J22L),kmadd(PDstandardNth33gt12,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12gt12,kmul(J32L,PDstandardNth13gt12)),kmul(J22L,kmul(J32L,PDstandardNth23gt12))),ToReal(2))))))));
      
      JacPDstandardNth22gt13 = 
        kmadd(dJ122L,PDstandardNth1gt13,kmadd(dJ222L,PDstandardNth2gt13,kmadd(dJ322L,PDstandardNth3gt13,kmadd(PDstandardNth11gt13,SQR(J12L),kmadd(PDstandardNth22gt13,SQR(J22L),kmadd(PDstandardNth33gt13,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12gt13,kmul(J32L,PDstandardNth13gt13)),kmul(J22L,kmul(J32L,PDstandardNth23gt13))),ToReal(2))))))));
      
      JacPDstandardNth22gt22 = 
        kmadd(dJ122L,PDstandardNth1gt22,kmadd(dJ222L,PDstandardNth2gt22,kmadd(dJ322L,PDstandardNth3gt22,kmadd(PDstandardNth11gt22,SQR(J12L),kmadd(PDstandardNth22gt22,SQR(J22L),kmadd(PDstandardNth33gt22,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12gt22,kmul(J32L,PDstandardNth13gt22)),kmul(J22L,kmul(J32L,PDstandardNth23gt22))),ToReal(2))))))));
      
      JacPDstandardNth22gt23 = 
        kmadd(dJ122L,PDstandardNth1gt23,kmadd(dJ222L,PDstandardNth2gt23,kmadd(dJ322L,PDstandardNth3gt23,kmadd(PDstandardNth11gt23,SQR(J12L),kmadd(PDstandardNth22gt23,SQR(J22L),kmadd(PDstandardNth33gt23,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12gt23,kmul(J32L,PDstandardNth13gt23)),kmul(J22L,kmul(J32L,PDstandardNth23gt23))),ToReal(2))))))));
      
      JacPDstandardNth22gt33 = 
        kmadd(dJ122L,PDstandardNth1gt33,kmadd(dJ222L,PDstandardNth2gt33,kmadd(dJ322L,PDstandardNth3gt33,kmadd(PDstandardNth11gt33,SQR(J12L),kmadd(PDstandardNth22gt33,SQR(J22L),kmadd(PDstandardNth33gt33,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12gt33,kmul(J32L,PDstandardNth13gt33)),kmul(J22L,kmul(J32L,PDstandardNth23gt33))),ToReal(2))))))));
      
      JacPDstandardNth22phi = 
        kmadd(dJ122L,PDstandardNth1phi,kmadd(dJ222L,PDstandardNth2phi,kmadd(dJ322L,PDstandardNth3phi,kmadd(PDstandardNth11phi,SQR(J12L),kmadd(PDstandardNth22phi,SQR(J22L),kmadd(PDstandardNth33phi,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12phi,kmul(J32L,PDstandardNth13phi)),kmul(J22L,kmul(J32L,PDstandardNth23phi))),ToReal(2))))))));
      
      JacPDstandardNth33alpha = 
        kmadd(dJ133L,PDstandardNth1alpha,kmadd(dJ233L,PDstandardNth2alpha,kmadd(dJ333L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,SQR(J13L),kmadd(PDstandardNth22alpha,SQR(J23L),kmadd(PDstandardNth33alpha,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmul(J23L,kmul(J33L,PDstandardNth23alpha))),ToReal(2))))))));
      
      JacPDstandardNth33gt11 = 
        kmadd(dJ133L,PDstandardNth1gt11,kmadd(dJ233L,PDstandardNth2gt11,kmadd(dJ333L,PDstandardNth3gt11,kmadd(PDstandardNth11gt11,SQR(J13L),kmadd(PDstandardNth22gt11,SQR(J23L),kmadd(PDstandardNth33gt11,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12gt11,kmul(J33L,PDstandardNth13gt11)),kmul(J23L,kmul(J33L,PDstandardNth23gt11))),ToReal(2))))))));
      
      JacPDstandardNth33gt12 = 
        kmadd(dJ133L,PDstandardNth1gt12,kmadd(dJ233L,PDstandardNth2gt12,kmadd(dJ333L,PDstandardNth3gt12,kmadd(PDstandardNth11gt12,SQR(J13L),kmadd(PDstandardNth22gt12,SQR(J23L),kmadd(PDstandardNth33gt12,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12gt12,kmul(J33L,PDstandardNth13gt12)),kmul(J23L,kmul(J33L,PDstandardNth23gt12))),ToReal(2))))))));
      
      JacPDstandardNth33gt13 = 
        kmadd(dJ133L,PDstandardNth1gt13,kmadd(dJ233L,PDstandardNth2gt13,kmadd(dJ333L,PDstandardNth3gt13,kmadd(PDstandardNth11gt13,SQR(J13L),kmadd(PDstandardNth22gt13,SQR(J23L),kmadd(PDstandardNth33gt13,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12gt13,kmul(J33L,PDstandardNth13gt13)),kmul(J23L,kmul(J33L,PDstandardNth23gt13))),ToReal(2))))))));
      
      JacPDstandardNth33gt22 = 
        kmadd(dJ133L,PDstandardNth1gt22,kmadd(dJ233L,PDstandardNth2gt22,kmadd(dJ333L,PDstandardNth3gt22,kmadd(PDstandardNth11gt22,SQR(J13L),kmadd(PDstandardNth22gt22,SQR(J23L),kmadd(PDstandardNth33gt22,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12gt22,kmul(J33L,PDstandardNth13gt22)),kmul(J23L,kmul(J33L,PDstandardNth23gt22))),ToReal(2))))))));
      
      JacPDstandardNth33gt23 = 
        kmadd(dJ133L,PDstandardNth1gt23,kmadd(dJ233L,PDstandardNth2gt23,kmadd(dJ333L,PDstandardNth3gt23,kmadd(PDstandardNth11gt23,SQR(J13L),kmadd(PDstandardNth22gt23,SQR(J23L),kmadd(PDstandardNth33gt23,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12gt23,kmul(J33L,PDstandardNth13gt23)),kmul(J23L,kmul(J33L,PDstandardNth23gt23))),ToReal(2))))))));
      
      JacPDstandardNth33gt33 = 
        kmadd(dJ133L,PDstandardNth1gt33,kmadd(dJ233L,PDstandardNth2gt33,kmadd(dJ333L,PDstandardNth3gt33,kmadd(PDstandardNth11gt33,SQR(J13L),kmadd(PDstandardNth22gt33,SQR(J23L),kmadd(PDstandardNth33gt33,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12gt33,kmul(J33L,PDstandardNth13gt33)),kmul(J23L,kmul(J33L,PDstandardNth23gt33))),ToReal(2))))))));
      
      JacPDstandardNth33phi = 
        kmadd(dJ133L,PDstandardNth1phi,kmadd(dJ233L,PDstandardNth2phi,kmadd(dJ333L,PDstandardNth3phi,kmadd(PDstandardNth11phi,SQR(J13L),kmadd(PDstandardNth22phi,SQR(J23L),kmadd(PDstandardNth33phi,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12phi,kmul(J33L,PDstandardNth13phi)),kmul(J23L,kmul(J33L,PDstandardNth23phi))),ToReal(2))))))));
      
      JacPDstandardNth12alpha = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11alpha,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha))),kmadd(J11L,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha)),kmadd(dJ112L,PDstandardNth1alpha,kmadd(J22L,kmadd(J21L,PDstandardNth22alpha,kmul(J31L,PDstandardNth23alpha)),kmadd(dJ212L,PDstandardNth2alpha,kmadd(J32L,kmadd(J21L,PDstandardNth23alpha,kmul(J31L,PDstandardNth33alpha)),kmul(dJ312L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth12gt11 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt11,kmadd(J21L,PDstandardNth12gt11,kmul(J31L,PDstandardNth13gt11))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt11,kmul(J32L,PDstandardNth13gt11)),kmadd(dJ112L,PDstandardNth1gt11,kmadd(J22L,kmadd(J21L,PDstandardNth22gt11,kmul(J31L,PDstandardNth23gt11)),kmadd(dJ212L,PDstandardNth2gt11,kmadd(J32L,kmadd(J21L,PDstandardNth23gt11,kmul(J31L,PDstandardNth33gt11)),kmul(dJ312L,PDstandardNth3gt11)))))));
      
      JacPDstandardNth12gt12 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt12,kmadd(J21L,PDstandardNth12gt12,kmul(J31L,PDstandardNth13gt12))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt12,kmul(J32L,PDstandardNth13gt12)),kmadd(dJ112L,PDstandardNth1gt12,kmadd(J22L,kmadd(J21L,PDstandardNth22gt12,kmul(J31L,PDstandardNth23gt12)),kmadd(dJ212L,PDstandardNth2gt12,kmadd(J32L,kmadd(J21L,PDstandardNth23gt12,kmul(J31L,PDstandardNth33gt12)),kmul(dJ312L,PDstandardNth3gt12)))))));
      
      JacPDstandardNth12gt13 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt13,kmadd(J21L,PDstandardNth12gt13,kmul(J31L,PDstandardNth13gt13))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt13,kmul(J32L,PDstandardNth13gt13)),kmadd(dJ112L,PDstandardNth1gt13,kmadd(J22L,kmadd(J21L,PDstandardNth22gt13,kmul(J31L,PDstandardNth23gt13)),kmadd(dJ212L,PDstandardNth2gt13,kmadd(J32L,kmadd(J21L,PDstandardNth23gt13,kmul(J31L,PDstandardNth33gt13)),kmul(dJ312L,PDstandardNth3gt13)))))));
      
      JacPDstandardNth12gt22 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt22,kmadd(J21L,PDstandardNth12gt22,kmul(J31L,PDstandardNth13gt22))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt22,kmul(J32L,PDstandardNth13gt22)),kmadd(dJ112L,PDstandardNth1gt22,kmadd(J22L,kmadd(J21L,PDstandardNth22gt22,kmul(J31L,PDstandardNth23gt22)),kmadd(dJ212L,PDstandardNth2gt22,kmadd(J32L,kmadd(J21L,PDstandardNth23gt22,kmul(J31L,PDstandardNth33gt22)),kmul(dJ312L,PDstandardNth3gt22)))))));
      
      JacPDstandardNth12gt23 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt23,kmadd(J21L,PDstandardNth12gt23,kmul(J31L,PDstandardNth13gt23))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt23,kmul(J32L,PDstandardNth13gt23)),kmadd(dJ112L,PDstandardNth1gt23,kmadd(J22L,kmadd(J21L,PDstandardNth22gt23,kmul(J31L,PDstandardNth23gt23)),kmadd(dJ212L,PDstandardNth2gt23,kmadd(J32L,kmadd(J21L,PDstandardNth23gt23,kmul(J31L,PDstandardNth33gt23)),kmul(dJ312L,PDstandardNth3gt23)))))));
      
      JacPDstandardNth12gt33 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt33,kmadd(J21L,PDstandardNth12gt33,kmul(J31L,PDstandardNth13gt33))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt33,kmul(J32L,PDstandardNth13gt33)),kmadd(dJ112L,PDstandardNth1gt33,kmadd(J22L,kmadd(J21L,PDstandardNth22gt33,kmul(J31L,PDstandardNth23gt33)),kmadd(dJ212L,PDstandardNth2gt33,kmadd(J32L,kmadd(J21L,PDstandardNth23gt33,kmul(J31L,PDstandardNth33gt33)),kmul(dJ312L,PDstandardNth3gt33)))))));
      
      JacPDstandardNth12phi = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11phi,kmadd(J21L,PDstandardNth12phi,kmul(J31L,PDstandardNth13phi))),kmadd(J11L,kmadd(J22L,PDstandardNth12phi,kmul(J32L,PDstandardNth13phi)),kmadd(dJ112L,PDstandardNth1phi,kmadd(J22L,kmadd(J21L,PDstandardNth22phi,kmul(J31L,PDstandardNth23phi)),kmadd(dJ212L,PDstandardNth2phi,kmadd(J32L,kmadd(J21L,PDstandardNth23phi,kmul(J31L,PDstandardNth33phi)),kmul(dJ312L,PDstandardNth3phi)))))));
      
      JacPDstandardNth13alpha = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11alpha,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha))),kmadd(J11L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmadd(dJ113L,PDstandardNth1alpha,kmadd(J23L,kmadd(J21L,PDstandardNth22alpha,kmul(J31L,PDstandardNth23alpha)),kmadd(dJ213L,PDstandardNth2alpha,kmadd(J33L,kmadd(J21L,PDstandardNth23alpha,kmul(J31L,PDstandardNth33alpha)),kmul(dJ313L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth13gt11 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt11,kmadd(J21L,PDstandardNth12gt11,kmul(J31L,PDstandardNth13gt11))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt11,kmul(J33L,PDstandardNth13gt11)),kmadd(dJ113L,PDstandardNth1gt11,kmadd(J23L,kmadd(J21L,PDstandardNth22gt11,kmul(J31L,PDstandardNth23gt11)),kmadd(dJ213L,PDstandardNth2gt11,kmadd(J33L,kmadd(J21L,PDstandardNth23gt11,kmul(J31L,PDstandardNth33gt11)),kmul(dJ313L,PDstandardNth3gt11)))))));
      
      JacPDstandardNth13gt12 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt12,kmadd(J21L,PDstandardNth12gt12,kmul(J31L,PDstandardNth13gt12))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt12,kmul(J33L,PDstandardNth13gt12)),kmadd(dJ113L,PDstandardNth1gt12,kmadd(J23L,kmadd(J21L,PDstandardNth22gt12,kmul(J31L,PDstandardNth23gt12)),kmadd(dJ213L,PDstandardNth2gt12,kmadd(J33L,kmadd(J21L,PDstandardNth23gt12,kmul(J31L,PDstandardNth33gt12)),kmul(dJ313L,PDstandardNth3gt12)))))));
      
      JacPDstandardNth13gt13 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt13,kmadd(J21L,PDstandardNth12gt13,kmul(J31L,PDstandardNth13gt13))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt13,kmul(J33L,PDstandardNth13gt13)),kmadd(dJ113L,PDstandardNth1gt13,kmadd(J23L,kmadd(J21L,PDstandardNth22gt13,kmul(J31L,PDstandardNth23gt13)),kmadd(dJ213L,PDstandardNth2gt13,kmadd(J33L,kmadd(J21L,PDstandardNth23gt13,kmul(J31L,PDstandardNth33gt13)),kmul(dJ313L,PDstandardNth3gt13)))))));
      
      JacPDstandardNth13gt22 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt22,kmadd(J21L,PDstandardNth12gt22,kmul(J31L,PDstandardNth13gt22))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt22,kmul(J33L,PDstandardNth13gt22)),kmadd(dJ113L,PDstandardNth1gt22,kmadd(J23L,kmadd(J21L,PDstandardNth22gt22,kmul(J31L,PDstandardNth23gt22)),kmadd(dJ213L,PDstandardNth2gt22,kmadd(J33L,kmadd(J21L,PDstandardNth23gt22,kmul(J31L,PDstandardNth33gt22)),kmul(dJ313L,PDstandardNth3gt22)))))));
      
      JacPDstandardNth13gt23 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt23,kmadd(J21L,PDstandardNth12gt23,kmul(J31L,PDstandardNth13gt23))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt23,kmul(J33L,PDstandardNth13gt23)),kmadd(dJ113L,PDstandardNth1gt23,kmadd(J23L,kmadd(J21L,PDstandardNth22gt23,kmul(J31L,PDstandardNth23gt23)),kmadd(dJ213L,PDstandardNth2gt23,kmadd(J33L,kmadd(J21L,PDstandardNth23gt23,kmul(J31L,PDstandardNth33gt23)),kmul(dJ313L,PDstandardNth3gt23)))))));
      
      JacPDstandardNth13gt33 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt33,kmadd(J21L,PDstandardNth12gt33,kmul(J31L,PDstandardNth13gt33))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt33,kmul(J33L,PDstandardNth13gt33)),kmadd(dJ113L,PDstandardNth1gt33,kmadd(J23L,kmadd(J21L,PDstandardNth22gt33,kmul(J31L,PDstandardNth23gt33)),kmadd(dJ213L,PDstandardNth2gt33,kmadd(J33L,kmadd(J21L,PDstandardNth23gt33,kmul(J31L,PDstandardNth33gt33)),kmul(dJ313L,PDstandardNth3gt33)))))));
      
      JacPDstandardNth13phi = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11phi,kmadd(J21L,PDstandardNth12phi,kmul(J31L,PDstandardNth13phi))),kmadd(J11L,kmadd(J23L,PDstandardNth12phi,kmul(J33L,PDstandardNth13phi)),kmadd(dJ113L,PDstandardNth1phi,kmadd(J23L,kmadd(J21L,PDstandardNth22phi,kmul(J31L,PDstandardNth23phi)),kmadd(dJ213L,PDstandardNth2phi,kmadd(J33L,kmadd(J21L,PDstandardNth23phi,kmul(J31L,PDstandardNth33phi)),kmul(dJ313L,PDstandardNth3phi)))))));
      
      JacPDstandardNth21gt11 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt11,kmadd(J21L,PDstandardNth12gt11,kmul(J31L,PDstandardNth13gt11))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt11,kmul(J32L,PDstandardNth13gt11)),kmadd(dJ112L,PDstandardNth1gt11,kmadd(J22L,kmadd(J21L,PDstandardNth22gt11,kmul(J31L,PDstandardNth23gt11)),kmadd(dJ212L,PDstandardNth2gt11,kmadd(J32L,kmadd(J21L,PDstandardNth23gt11,kmul(J31L,PDstandardNth33gt11)),kmul(dJ312L,PDstandardNth3gt11)))))));
      
      JacPDstandardNth21gt12 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt12,kmadd(J21L,PDstandardNth12gt12,kmul(J31L,PDstandardNth13gt12))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt12,kmul(J32L,PDstandardNth13gt12)),kmadd(dJ112L,PDstandardNth1gt12,kmadd(J22L,kmadd(J21L,PDstandardNth22gt12,kmul(J31L,PDstandardNth23gt12)),kmadd(dJ212L,PDstandardNth2gt12,kmadd(J32L,kmadd(J21L,PDstandardNth23gt12,kmul(J31L,PDstandardNth33gt12)),kmul(dJ312L,PDstandardNth3gt12)))))));
      
      JacPDstandardNth21gt13 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt13,kmadd(J21L,PDstandardNth12gt13,kmul(J31L,PDstandardNth13gt13))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt13,kmul(J32L,PDstandardNth13gt13)),kmadd(dJ112L,PDstandardNth1gt13,kmadd(J22L,kmadd(J21L,PDstandardNth22gt13,kmul(J31L,PDstandardNth23gt13)),kmadd(dJ212L,PDstandardNth2gt13,kmadd(J32L,kmadd(J21L,PDstandardNth23gt13,kmul(J31L,PDstandardNth33gt13)),kmul(dJ312L,PDstandardNth3gt13)))))));
      
      JacPDstandardNth21gt22 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt22,kmadd(J21L,PDstandardNth12gt22,kmul(J31L,PDstandardNth13gt22))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt22,kmul(J32L,PDstandardNth13gt22)),kmadd(dJ112L,PDstandardNth1gt22,kmadd(J22L,kmadd(J21L,PDstandardNth22gt22,kmul(J31L,PDstandardNth23gt22)),kmadd(dJ212L,PDstandardNth2gt22,kmadd(J32L,kmadd(J21L,PDstandardNth23gt22,kmul(J31L,PDstandardNth33gt22)),kmul(dJ312L,PDstandardNth3gt22)))))));
      
      JacPDstandardNth21gt23 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt23,kmadd(J21L,PDstandardNth12gt23,kmul(J31L,PDstandardNth13gt23))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt23,kmul(J32L,PDstandardNth13gt23)),kmadd(dJ112L,PDstandardNth1gt23,kmadd(J22L,kmadd(J21L,PDstandardNth22gt23,kmul(J31L,PDstandardNth23gt23)),kmadd(dJ212L,PDstandardNth2gt23,kmadd(J32L,kmadd(J21L,PDstandardNth23gt23,kmul(J31L,PDstandardNth33gt23)),kmul(dJ312L,PDstandardNth3gt23)))))));
      
      JacPDstandardNth21gt33 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11gt33,kmadd(J21L,PDstandardNth12gt33,kmul(J31L,PDstandardNth13gt33))),kmadd(J11L,kmadd(J22L,PDstandardNth12gt33,kmul(J32L,PDstandardNth13gt33)),kmadd(dJ112L,PDstandardNth1gt33,kmadd(J22L,kmadd(J21L,PDstandardNth22gt33,kmul(J31L,PDstandardNth23gt33)),kmadd(dJ212L,PDstandardNth2gt33,kmadd(J32L,kmadd(J21L,PDstandardNth23gt33,kmul(J31L,PDstandardNth33gt33)),kmul(dJ312L,PDstandardNth3gt33)))))));
      
      JacPDstandardNth23alpha = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11alpha,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha))),kmadd(J12L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmadd(dJ123L,PDstandardNth1alpha,kmadd(J23L,kmadd(J22L,PDstandardNth22alpha,kmul(J32L,PDstandardNth23alpha)),kmadd(dJ223L,PDstandardNth2alpha,kmadd(J33L,kmadd(J22L,PDstandardNth23alpha,kmul(J32L,PDstandardNth33alpha)),kmul(dJ323L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth23gt11 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt11,kmadd(J22L,PDstandardNth12gt11,kmul(J32L,PDstandardNth13gt11))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt11,kmul(J33L,PDstandardNth13gt11)),kmadd(dJ123L,PDstandardNth1gt11,kmadd(J23L,kmadd(J22L,PDstandardNth22gt11,kmul(J32L,PDstandardNth23gt11)),kmadd(dJ223L,PDstandardNth2gt11,kmadd(J33L,kmadd(J22L,PDstandardNth23gt11,kmul(J32L,PDstandardNth33gt11)),kmul(dJ323L,PDstandardNth3gt11)))))));
      
      JacPDstandardNth23gt12 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt12,kmadd(J22L,PDstandardNth12gt12,kmul(J32L,PDstandardNth13gt12))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt12,kmul(J33L,PDstandardNth13gt12)),kmadd(dJ123L,PDstandardNth1gt12,kmadd(J23L,kmadd(J22L,PDstandardNth22gt12,kmul(J32L,PDstandardNth23gt12)),kmadd(dJ223L,PDstandardNth2gt12,kmadd(J33L,kmadd(J22L,PDstandardNth23gt12,kmul(J32L,PDstandardNth33gt12)),kmul(dJ323L,PDstandardNth3gt12)))))));
      
      JacPDstandardNth23gt13 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt13,kmadd(J22L,PDstandardNth12gt13,kmul(J32L,PDstandardNth13gt13))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt13,kmul(J33L,PDstandardNth13gt13)),kmadd(dJ123L,PDstandardNth1gt13,kmadd(J23L,kmadd(J22L,PDstandardNth22gt13,kmul(J32L,PDstandardNth23gt13)),kmadd(dJ223L,PDstandardNth2gt13,kmadd(J33L,kmadd(J22L,PDstandardNth23gt13,kmul(J32L,PDstandardNth33gt13)),kmul(dJ323L,PDstandardNth3gt13)))))));
      
      JacPDstandardNth23gt22 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt22,kmadd(J22L,PDstandardNth12gt22,kmul(J32L,PDstandardNth13gt22))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt22,kmul(J33L,PDstandardNth13gt22)),kmadd(dJ123L,PDstandardNth1gt22,kmadd(J23L,kmadd(J22L,PDstandardNth22gt22,kmul(J32L,PDstandardNth23gt22)),kmadd(dJ223L,PDstandardNth2gt22,kmadd(J33L,kmadd(J22L,PDstandardNth23gt22,kmul(J32L,PDstandardNth33gt22)),kmul(dJ323L,PDstandardNth3gt22)))))));
      
      JacPDstandardNth23gt23 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt23,kmadd(J22L,PDstandardNth12gt23,kmul(J32L,PDstandardNth13gt23))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt23,kmul(J33L,PDstandardNth13gt23)),kmadd(dJ123L,PDstandardNth1gt23,kmadd(J23L,kmadd(J22L,PDstandardNth22gt23,kmul(J32L,PDstandardNth23gt23)),kmadd(dJ223L,PDstandardNth2gt23,kmadd(J33L,kmadd(J22L,PDstandardNth23gt23,kmul(J32L,PDstandardNth33gt23)),kmul(dJ323L,PDstandardNth3gt23)))))));
      
      JacPDstandardNth23gt33 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt33,kmadd(J22L,PDstandardNth12gt33,kmul(J32L,PDstandardNth13gt33))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt33,kmul(J33L,PDstandardNth13gt33)),kmadd(dJ123L,PDstandardNth1gt33,kmadd(J23L,kmadd(J22L,PDstandardNth22gt33,kmul(J32L,PDstandardNth23gt33)),kmadd(dJ223L,PDstandardNth2gt33,kmadd(J33L,kmadd(J22L,PDstandardNth23gt33,kmul(J32L,PDstandardNth33gt33)),kmul(dJ323L,PDstandardNth3gt33)))))));
      
      JacPDstandardNth23phi = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11phi,kmadd(J22L,PDstandardNth12phi,kmul(J32L,PDstandardNth13phi))),kmadd(J12L,kmadd(J23L,PDstandardNth12phi,kmul(J33L,PDstandardNth13phi)),kmadd(dJ123L,PDstandardNth1phi,kmadd(J23L,kmadd(J22L,PDstandardNth22phi,kmul(J32L,PDstandardNth23phi)),kmadd(dJ223L,PDstandardNth2phi,kmadd(J33L,kmadd(J22L,PDstandardNth23phi,kmul(J32L,PDstandardNth33phi)),kmul(dJ323L,PDstandardNth3phi)))))));
      
      JacPDstandardNth31gt11 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt11,kmadd(J21L,PDstandardNth12gt11,kmul(J31L,PDstandardNth13gt11))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt11,kmul(J33L,PDstandardNth13gt11)),kmadd(dJ113L,PDstandardNth1gt11,kmadd(J23L,kmadd(J21L,PDstandardNth22gt11,kmul(J31L,PDstandardNth23gt11)),kmadd(dJ213L,PDstandardNth2gt11,kmadd(J33L,kmadd(J21L,PDstandardNth23gt11,kmul(J31L,PDstandardNth33gt11)),kmul(dJ313L,PDstandardNth3gt11)))))));
      
      JacPDstandardNth31gt12 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt12,kmadd(J21L,PDstandardNth12gt12,kmul(J31L,PDstandardNth13gt12))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt12,kmul(J33L,PDstandardNth13gt12)),kmadd(dJ113L,PDstandardNth1gt12,kmadd(J23L,kmadd(J21L,PDstandardNth22gt12,kmul(J31L,PDstandardNth23gt12)),kmadd(dJ213L,PDstandardNth2gt12,kmadd(J33L,kmadd(J21L,PDstandardNth23gt12,kmul(J31L,PDstandardNth33gt12)),kmul(dJ313L,PDstandardNth3gt12)))))));
      
      JacPDstandardNth31gt13 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt13,kmadd(J21L,PDstandardNth12gt13,kmul(J31L,PDstandardNth13gt13))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt13,kmul(J33L,PDstandardNth13gt13)),kmadd(dJ113L,PDstandardNth1gt13,kmadd(J23L,kmadd(J21L,PDstandardNth22gt13,kmul(J31L,PDstandardNth23gt13)),kmadd(dJ213L,PDstandardNth2gt13,kmadd(J33L,kmadd(J21L,PDstandardNth23gt13,kmul(J31L,PDstandardNth33gt13)),kmul(dJ313L,PDstandardNth3gt13)))))));
      
      JacPDstandardNth31gt22 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt22,kmadd(J21L,PDstandardNth12gt22,kmul(J31L,PDstandardNth13gt22))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt22,kmul(J33L,PDstandardNth13gt22)),kmadd(dJ113L,PDstandardNth1gt22,kmadd(J23L,kmadd(J21L,PDstandardNth22gt22,kmul(J31L,PDstandardNth23gt22)),kmadd(dJ213L,PDstandardNth2gt22,kmadd(J33L,kmadd(J21L,PDstandardNth23gt22,kmul(J31L,PDstandardNth33gt22)),kmul(dJ313L,PDstandardNth3gt22)))))));
      
      JacPDstandardNth31gt23 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt23,kmadd(J21L,PDstandardNth12gt23,kmul(J31L,PDstandardNth13gt23))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt23,kmul(J33L,PDstandardNth13gt23)),kmadd(dJ113L,PDstandardNth1gt23,kmadd(J23L,kmadd(J21L,PDstandardNth22gt23,kmul(J31L,PDstandardNth23gt23)),kmadd(dJ213L,PDstandardNth2gt23,kmadd(J33L,kmadd(J21L,PDstandardNth23gt23,kmul(J31L,PDstandardNth33gt23)),kmul(dJ313L,PDstandardNth3gt23)))))));
      
      JacPDstandardNth31gt33 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11gt33,kmadd(J21L,PDstandardNth12gt33,kmul(J31L,PDstandardNth13gt33))),kmadd(J11L,kmadd(J23L,PDstandardNth12gt33,kmul(J33L,PDstandardNth13gt33)),kmadd(dJ113L,PDstandardNth1gt33,kmadd(J23L,kmadd(J21L,PDstandardNth22gt33,kmul(J31L,PDstandardNth23gt33)),kmadd(dJ213L,PDstandardNth2gt33,kmadd(J33L,kmadd(J21L,PDstandardNth23gt33,kmul(J31L,PDstandardNth33gt33)),kmul(dJ313L,PDstandardNth3gt33)))))));
      
      JacPDstandardNth32gt11 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt11,kmadd(J22L,PDstandardNth12gt11,kmul(J32L,PDstandardNth13gt11))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt11,kmul(J33L,PDstandardNth13gt11)),kmadd(dJ123L,PDstandardNth1gt11,kmadd(J23L,kmadd(J22L,PDstandardNth22gt11,kmul(J32L,PDstandardNth23gt11)),kmadd(dJ223L,PDstandardNth2gt11,kmadd(J33L,kmadd(J22L,PDstandardNth23gt11,kmul(J32L,PDstandardNth33gt11)),kmul(dJ323L,PDstandardNth3gt11)))))));
      
      JacPDstandardNth32gt12 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt12,kmadd(J22L,PDstandardNth12gt12,kmul(J32L,PDstandardNth13gt12))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt12,kmul(J33L,PDstandardNth13gt12)),kmadd(dJ123L,PDstandardNth1gt12,kmadd(J23L,kmadd(J22L,PDstandardNth22gt12,kmul(J32L,PDstandardNth23gt12)),kmadd(dJ223L,PDstandardNth2gt12,kmadd(J33L,kmadd(J22L,PDstandardNth23gt12,kmul(J32L,PDstandardNth33gt12)),kmul(dJ323L,PDstandardNth3gt12)))))));
      
      JacPDstandardNth32gt13 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt13,kmadd(J22L,PDstandardNth12gt13,kmul(J32L,PDstandardNth13gt13))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt13,kmul(J33L,PDstandardNth13gt13)),kmadd(dJ123L,PDstandardNth1gt13,kmadd(J23L,kmadd(J22L,PDstandardNth22gt13,kmul(J32L,PDstandardNth23gt13)),kmadd(dJ223L,PDstandardNth2gt13,kmadd(J33L,kmadd(J22L,PDstandardNth23gt13,kmul(J32L,PDstandardNth33gt13)),kmul(dJ323L,PDstandardNth3gt13)))))));
      
      JacPDstandardNth32gt22 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt22,kmadd(J22L,PDstandardNth12gt22,kmul(J32L,PDstandardNth13gt22))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt22,kmul(J33L,PDstandardNth13gt22)),kmadd(dJ123L,PDstandardNth1gt22,kmadd(J23L,kmadd(J22L,PDstandardNth22gt22,kmul(J32L,PDstandardNth23gt22)),kmadd(dJ223L,PDstandardNth2gt22,kmadd(J33L,kmadd(J22L,PDstandardNth23gt22,kmul(J32L,PDstandardNth33gt22)),kmul(dJ323L,PDstandardNth3gt22)))))));
      
      JacPDstandardNth32gt23 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt23,kmadd(J22L,PDstandardNth12gt23,kmul(J32L,PDstandardNth13gt23))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt23,kmul(J33L,PDstandardNth13gt23)),kmadd(dJ123L,PDstandardNth1gt23,kmadd(J23L,kmadd(J22L,PDstandardNth22gt23,kmul(J32L,PDstandardNth23gt23)),kmadd(dJ223L,PDstandardNth2gt23,kmadd(J33L,kmadd(J22L,PDstandardNth23gt23,kmul(J32L,PDstandardNth33gt23)),kmul(dJ323L,PDstandardNth3gt23)))))));
      
      JacPDstandardNth32gt33 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11gt33,kmadd(J22L,PDstandardNth12gt33,kmul(J32L,PDstandardNth13gt33))),kmadd(J12L,kmadd(J23L,PDstandardNth12gt33,kmul(J33L,PDstandardNth13gt33)),kmadd(dJ123L,PDstandardNth1gt33,kmadd(J23L,kmadd(J22L,PDstandardNth22gt33,kmul(J32L,PDstandardNth23gt33)),kmadd(dJ223L,PDstandardNth2gt33,kmadd(J33L,kmadd(J22L,PDstandardNth23gt33,kmul(J32L,PDstandardNth33gt33)),kmul(dJ323L,PDstandardNth3gt33)))))));
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
      
      JacPDstandardNth1Xt1 = PDstandardNth1Xt1;
      
      JacPDstandardNth1Xt2 = PDstandardNth1Xt2;
      
      JacPDstandardNth1Xt3 = PDstandardNth1Xt3;
      
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
      
      JacPDstandardNth2Xt1 = PDstandardNth2Xt1;
      
      JacPDstandardNth2Xt2 = PDstandardNth2Xt2;
      
      JacPDstandardNth2Xt3 = PDstandardNth2Xt3;
      
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
      
      JacPDstandardNth3Xt1 = PDstandardNth3Xt1;
      
      JacPDstandardNth3Xt2 = PDstandardNth3Xt2;
      
      JacPDstandardNth3Xt3 = PDstandardNth3Xt3;
      
      JacPDstandardNth11alpha = PDstandardNth11alpha;
      
      JacPDstandardNth11gt11 = PDstandardNth11gt11;
      
      JacPDstandardNth11gt12 = PDstandardNth11gt12;
      
      JacPDstandardNth11gt13 = PDstandardNth11gt13;
      
      JacPDstandardNth11gt22 = PDstandardNth11gt22;
      
      JacPDstandardNth11gt23 = PDstandardNth11gt23;
      
      JacPDstandardNth11gt33 = PDstandardNth11gt33;
      
      JacPDstandardNth11phi = PDstandardNth11phi;
      
      JacPDstandardNth22alpha = PDstandardNth22alpha;
      
      JacPDstandardNth22gt11 = PDstandardNth22gt11;
      
      JacPDstandardNth22gt12 = PDstandardNth22gt12;
      
      JacPDstandardNth22gt13 = PDstandardNth22gt13;
      
      JacPDstandardNth22gt22 = PDstandardNth22gt22;
      
      JacPDstandardNth22gt23 = PDstandardNth22gt23;
      
      JacPDstandardNth22gt33 = PDstandardNth22gt33;
      
      JacPDstandardNth22phi = PDstandardNth22phi;
      
      JacPDstandardNth33alpha = PDstandardNth33alpha;
      
      JacPDstandardNth33gt11 = PDstandardNth33gt11;
      
      JacPDstandardNth33gt12 = PDstandardNth33gt12;
      
      JacPDstandardNth33gt13 = PDstandardNth33gt13;
      
      JacPDstandardNth33gt22 = PDstandardNth33gt22;
      
      JacPDstandardNth33gt23 = PDstandardNth33gt23;
      
      JacPDstandardNth33gt33 = PDstandardNth33gt33;
      
      JacPDstandardNth33phi = PDstandardNth33phi;
      
      JacPDstandardNth12alpha = PDstandardNth12alpha;
      
      JacPDstandardNth12gt11 = PDstandardNth12gt11;
      
      JacPDstandardNth12gt12 = PDstandardNth12gt12;
      
      JacPDstandardNth12gt13 = PDstandardNth12gt13;
      
      JacPDstandardNth12gt22 = PDstandardNth12gt22;
      
      JacPDstandardNth12gt23 = PDstandardNth12gt23;
      
      JacPDstandardNth12gt33 = PDstandardNth12gt33;
      
      JacPDstandardNth12phi = PDstandardNth12phi;
      
      JacPDstandardNth13alpha = PDstandardNth13alpha;
      
      JacPDstandardNth13gt11 = PDstandardNth13gt11;
      
      JacPDstandardNth13gt12 = PDstandardNth13gt12;
      
      JacPDstandardNth13gt13 = PDstandardNth13gt13;
      
      JacPDstandardNth13gt22 = PDstandardNth13gt22;
      
      JacPDstandardNth13gt23 = PDstandardNth13gt23;
      
      JacPDstandardNth13gt33 = PDstandardNth13gt33;
      
      JacPDstandardNth13phi = PDstandardNth13phi;
      
      JacPDstandardNth21gt11 = PDstandardNth12gt11;
      
      JacPDstandardNth21gt12 = PDstandardNth12gt12;
      
      JacPDstandardNth21gt13 = PDstandardNth12gt13;
      
      JacPDstandardNth21gt22 = PDstandardNth12gt22;
      
      JacPDstandardNth21gt23 = PDstandardNth12gt23;
      
      JacPDstandardNth21gt33 = PDstandardNth12gt33;
      
      JacPDstandardNth23alpha = PDstandardNth23alpha;
      
      JacPDstandardNth23gt11 = PDstandardNth23gt11;
      
      JacPDstandardNth23gt12 = PDstandardNth23gt12;
      
      JacPDstandardNth23gt13 = PDstandardNth23gt13;
      
      JacPDstandardNth23gt22 = PDstandardNth23gt22;
      
      JacPDstandardNth23gt23 = PDstandardNth23gt23;
      
      JacPDstandardNth23gt33 = PDstandardNth23gt33;
      
      JacPDstandardNth23phi = PDstandardNth23phi;
      
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
    
    CCTK_REAL_VEC detgt = ToReal(1);
    
    CCTK_REAL_VEC gtu11 = 
      kmul(INV(detgt),kmsub(gt22L,gt33L,SQR(gt23L)));
    
    CCTK_REAL_VEC gtu12 = 
      kmul(INV(detgt),kmsub(gt13L,gt23L,kmul(gt12L,gt33L)));
    
    CCTK_REAL_VEC gtu13 = 
      kmul(INV(detgt),kmsub(gt12L,gt23L,kmul(gt13L,gt22L)));
    
    CCTK_REAL_VEC gtu22 = 
      kmul(INV(detgt),kmsub(gt11L,gt33L,SQR(gt13L)));
    
    CCTK_REAL_VEC gtu23 = 
      kmul(INV(detgt),kmsub(gt12L,gt13L,kmul(gt11L,gt23L)));
    
    CCTK_REAL_VEC gtu33 = 
      kmul(INV(detgt),kmsub(gt11L,gt22L,SQR(gt12L)));
    
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
    
    CCTK_REAL_VEC Gtlu111 = 
      kmadd(Gtl111,gtu11,kmadd(Gtl112,gtu12,kmul(Gtl113,gtu13)));
    
    CCTK_REAL_VEC Gtlu112 = 
      kmadd(Gtl111,gtu12,kmadd(Gtl112,gtu22,kmul(Gtl113,gtu23)));
    
    CCTK_REAL_VEC Gtlu113 = 
      kmadd(Gtl111,gtu13,kmadd(Gtl112,gtu23,kmul(Gtl113,gtu33)));
    
    CCTK_REAL_VEC Gtlu121 = 
      kmadd(Gtl112,gtu11,kmadd(Gtl122,gtu12,kmul(Gtl123,gtu13)));
    
    CCTK_REAL_VEC Gtlu122 = 
      kmadd(Gtl112,gtu12,kmadd(Gtl122,gtu22,kmul(Gtl123,gtu23)));
    
    CCTK_REAL_VEC Gtlu123 = 
      kmadd(Gtl112,gtu13,kmadd(Gtl122,gtu23,kmul(Gtl123,gtu33)));
    
    CCTK_REAL_VEC Gtlu131 = 
      kmadd(Gtl113,gtu11,kmadd(Gtl123,gtu12,kmul(Gtl133,gtu13)));
    
    CCTK_REAL_VEC Gtlu132 = 
      kmadd(Gtl113,gtu12,kmadd(Gtl123,gtu22,kmul(Gtl133,gtu23)));
    
    CCTK_REAL_VEC Gtlu133 = 
      kmadd(Gtl113,gtu13,kmadd(Gtl123,gtu23,kmul(Gtl133,gtu33)));
    
    CCTK_REAL_VEC Gtlu211 = 
      kmadd(Gtl211,gtu11,kmadd(Gtl212,gtu12,kmul(Gtl213,gtu13)));
    
    CCTK_REAL_VEC Gtlu212 = 
      kmadd(Gtl211,gtu12,kmadd(Gtl212,gtu22,kmul(Gtl213,gtu23)));
    
    CCTK_REAL_VEC Gtlu213 = 
      kmadd(Gtl211,gtu13,kmadd(Gtl212,gtu23,kmul(Gtl213,gtu33)));
    
    CCTK_REAL_VEC Gtlu221 = 
      kmadd(Gtl212,gtu11,kmadd(Gtl222,gtu12,kmul(Gtl223,gtu13)));
    
    CCTK_REAL_VEC Gtlu222 = 
      kmadd(Gtl212,gtu12,kmadd(Gtl222,gtu22,kmul(Gtl223,gtu23)));
    
    CCTK_REAL_VEC Gtlu223 = 
      kmadd(Gtl212,gtu13,kmadd(Gtl222,gtu23,kmul(Gtl223,gtu33)));
    
    CCTK_REAL_VEC Gtlu231 = 
      kmadd(Gtl213,gtu11,kmadd(Gtl223,gtu12,kmul(Gtl233,gtu13)));
    
    CCTK_REAL_VEC Gtlu232 = 
      kmadd(Gtl213,gtu12,kmadd(Gtl223,gtu22,kmul(Gtl233,gtu23)));
    
    CCTK_REAL_VEC Gtlu233 = 
      kmadd(Gtl213,gtu13,kmadd(Gtl223,gtu23,kmul(Gtl233,gtu33)));
    
    CCTK_REAL_VEC Gtlu311 = 
      kmadd(Gtl311,gtu11,kmadd(Gtl312,gtu12,kmul(Gtl313,gtu13)));
    
    CCTK_REAL_VEC Gtlu312 = 
      kmadd(Gtl311,gtu12,kmadd(Gtl312,gtu22,kmul(Gtl313,gtu23)));
    
    CCTK_REAL_VEC Gtlu313 = 
      kmadd(Gtl311,gtu13,kmadd(Gtl312,gtu23,kmul(Gtl313,gtu33)));
    
    CCTK_REAL_VEC Gtlu321 = 
      kmadd(Gtl312,gtu11,kmadd(Gtl322,gtu12,kmul(Gtl323,gtu13)));
    
    CCTK_REAL_VEC Gtlu322 = 
      kmadd(Gtl312,gtu12,kmadd(Gtl322,gtu22,kmul(Gtl323,gtu23)));
    
    CCTK_REAL_VEC Gtlu323 = 
      kmadd(Gtl312,gtu13,kmadd(Gtl322,gtu23,kmul(Gtl323,gtu33)));
    
    CCTK_REAL_VEC Gtlu331 = 
      kmadd(Gtl313,gtu11,kmadd(Gtl323,gtu12,kmul(Gtl333,gtu13)));
    
    CCTK_REAL_VEC Gtlu332 = 
      kmadd(Gtl313,gtu12,kmadd(Gtl323,gtu22,kmul(Gtl333,gtu23)));
    
    CCTK_REAL_VEC Gtlu333 = 
      kmadd(Gtl313,gtu13,kmadd(Gtl323,gtu23,kmul(Gtl333,gtu33)));
    
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
    
    CCTK_REAL_VEC Rt11 = 
      kmul(ToReal(0.5),knmsub(gtu11,JacPDstandardNth11gt11,knmsub(gtu22,JacPDstandardNth22gt11,knmsub(gtu33,JacPDstandardNth33gt11,knmsub(gtu12,kadd(JacPDstandardNth21gt11,JacPDstandardNth12gt11),knmsub(gtu13,kadd(JacPDstandardNth31gt11,JacPDstandardNth13gt11),knmsub(gtu23,kadd(JacPDstandardNth32gt11,JacPDstandardNth23gt11),kmadd(kmadd(Gt211,Gtlu211,kmadd(Gt212,Gtlu212,kmadd(Gt213,Gtlu213,kmadd(Gt311,Gtlu311,kmadd(Gt312,Gtlu312,kmadd(Gt313,Gtlu313,kmul(gt11L,JacPDstandardNth1Xt1))))))),ToReal(2),kmadd(gt12L,kmul(JacPDstandardNth1Xt2,ToReal(2)),kmadd(gt13L,kmul(JacPDstandardNth1Xt3,ToReal(2)),kmadd(Gtl111,kmul(Xtn1,ToReal(2)),kmadd(Gtl112,kmul(Xtn2,ToReal(2)),kmadd(Gtl113,kmul(Xtn3,ToReal(2)),kmadd(kmadd(Gt211,Gtlu121,kmadd(Gt212,Gtlu122,kmadd(Gt213,Gtlu123,kmadd(Gt311,Gtlu131,kmadd(Gt312,Gtlu132,kmul(Gt313,Gtlu133)))))),ToReal(4),kmul(kmadd(Gt111,Gtlu111,kmadd(Gt112,Gtlu112,kmul(Gt113,Gtlu113))),ToReal(6))))))))))))))));
    
    CCTK_REAL_VEC Rt12 = 
      kmul(ToReal(0.5),kmadd(gt12L,JacPDstandardNth1Xt1,kmadd(gt22L,JacPDstandardNth1Xt2,kmadd(gt23L,JacPDstandardNth1Xt3,kmadd(gt11L,JacPDstandardNth2Xt1,kmadd(gt12L,JacPDstandardNth2Xt2,kmadd(gt13L,JacPDstandardNth2Xt3,kmadd(Gtl112,Xtn1,kmadd(Gtl211,Xtn1,kmadd(Gtl122,Xtn2,kmadd(Gtl212,Xtn2,kmadd(Gtl123,Xtn3,kmadd(Gtl213,Xtn3,knmsub(gtu11,JacPDstandardNth11gt12,knmsub(gtu22,JacPDstandardNth22gt12,knmsub(gtu33,JacPDstandardNth33gt12,knmsub(gtu12,kadd(JacPDstandardNth21gt12,JacPDstandardNth12gt12),knmsub(gtu13,kadd(JacPDstandardNth31gt12,JacPDstandardNth13gt12),knmsub(gtu23,kadd(JacPDstandardNth32gt12,JacPDstandardNth23gt12),kmadd(kmadd(Gt122,Gtlu112,kmadd(Gt123,Gtlu113,kmadd(Gt111,Gtlu121,kmadd(Gt212,Gtlu121,kmadd(Gt222,Gtlu122,kmadd(Gt113,Gtlu123,kmadd(Gt223,Gtlu123,kmadd(Gt312,Gtlu131,kmadd(Gt322,Gtlu132,kmadd(Gt323,Gtlu133,kmadd(Gt111,Gtlu211,kmadd(Gt112,kadd(Gtlu111,kadd(Gtlu122,Gtlu212)),kmadd(Gt113,Gtlu213,kmadd(Gt311,Gtlu231,kmadd(Gt312,Gtlu232,kmadd(Gt313,Gtlu233,kmadd(Gt311,Gtlu321,kmadd(Gt312,Gtlu322,kmul(Gt313,Gtlu323))))))))))))))))))),ToReal(2),kmul(kmadd(Gt211,Gtlu221,kmadd(Gt212,Gtlu222,kmul(Gt213,Gtlu223))),ToReal(4))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt13 = 
      kmul(ToReal(0.5),kmadd(gt13L,JacPDstandardNth1Xt1,kmadd(gt23L,JacPDstandardNth1Xt2,kmadd(gt33L,JacPDstandardNth1Xt3,kmadd(gt11L,JacPDstandardNth3Xt1,kmadd(gt12L,JacPDstandardNth3Xt2,kmadd(gt13L,JacPDstandardNth3Xt3,kmadd(Gtl113,Xtn1,kmadd(Gtl311,Xtn1,kmadd(Gtl123,Xtn2,kmadd(Gtl312,Xtn2,kmadd(Gtl133,Xtn3,kmadd(Gtl313,Xtn3,knmsub(gtu11,JacPDstandardNth11gt13,knmsub(gtu22,JacPDstandardNth22gt13,knmsub(gtu33,JacPDstandardNth33gt13,knmsub(gtu12,kadd(JacPDstandardNth21gt13,JacPDstandardNth12gt13),knmsub(gtu13,kadd(JacPDstandardNth31gt13,JacPDstandardNth13gt13),knmsub(gtu23,kadd(JacPDstandardNth32gt13,JacPDstandardNth23gt13),kmadd(kmadd(Gt123,Gtlu112,kmadd(Gt133,Gtlu113,kmadd(Gt213,Gtlu121,kmadd(Gt223,Gtlu122,kmadd(Gt233,Gtlu123,kmadd(Gt111,Gtlu131,kmadd(Gt313,Gtlu131,kmadd(Gt112,Gtlu132,kmadd(Gt323,Gtlu132,kmadd(Gt333,Gtlu133,kmadd(Gt211,Gtlu231,kmadd(Gt212,Gtlu232,kmadd(Gt213,Gtlu233,kmadd(Gt111,Gtlu311,kmadd(Gt112,Gtlu312,kmadd(Gt113,kadd(Gtlu111,kadd(Gtlu133,Gtlu313)),kmadd(Gt211,Gtlu321,kmadd(Gt212,Gtlu322,kmul(Gt213,Gtlu323))))))))))))))))))),ToReal(2),kmul(kmadd(Gt311,Gtlu331,kmadd(Gt312,Gtlu332,kmul(Gt313,Gtlu333))),ToReal(4))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt22 = 
      kmul(ToReal(0.5),knmsub(gtu11,JacPDstandardNth11gt22,knmsub(gtu22,JacPDstandardNth22gt22,knmsub(gtu33,JacPDstandardNth33gt22,knmsub(gtu12,kadd(JacPDstandardNth21gt22,JacPDstandardNth12gt22),knmsub(gtu13,kadd(JacPDstandardNth31gt22,JacPDstandardNth13gt22),knmsub(gtu23,kadd(JacPDstandardNth32gt22,JacPDstandardNth23gt22),kmadd(gt22L,kmul(JacPDstandardNth2Xt2,ToReal(2)),kmadd(gt23L,kmul(JacPDstandardNth2Xt3,ToReal(2)),kmadd(Gtl212,kmul(Xtn1,ToReal(2)),kmadd(Gtl222,kmul(Xtn2,ToReal(2)),kmadd(Gtl223,kmul(Xtn3,ToReal(2)),kmadd(ToReal(2),kmadd(Gt123,Gtlu123,kmadd(Gt312,Gtlu321,kmadd(Gt322,Gtlu322,kmadd(Gt323,Gtlu323,kmadd(gt12L,JacPDstandardNth2Xt1,kmadd(Gt112,kmadd(Gtlu211,ToReal(2),Gtlu121),kmul(Gt122,kmadd(Gtlu212,ToReal(2),Gtlu122)))))))),kmadd(kmadd(Gt123,Gtlu213,kmadd(Gt312,Gtlu231,kmadd(Gt322,Gtlu232,kmul(Gt323,Gtlu233)))),ToReal(4),kmul(kmadd(Gt212,Gtlu221,kmadd(Gt222,Gtlu222,kmul(Gt223,Gtlu223))),ToReal(6))))))))))))))));
    
    CCTK_REAL_VEC Rt23 = 
      kmul(ToReal(0.5),kmadd(gt13L,JacPDstandardNth2Xt1,kmadd(gt23L,JacPDstandardNth2Xt2,kmadd(gt33L,JacPDstandardNth2Xt3,kmadd(gt12L,JacPDstandardNth3Xt1,kmadd(gt22L,JacPDstandardNth3Xt2,kmadd(gt23L,JacPDstandardNth3Xt3,kmadd(Gtl213,Xtn1,kmadd(Gtl312,Xtn1,kmadd(Gtl223,Xtn2,kmadd(Gtl322,Xtn2,kmadd(Gtl233,Xtn3,kmadd(Gtl323,Xtn3,knmsub(gtu11,JacPDstandardNth11gt23,knmsub(gtu22,JacPDstandardNth22gt23,knmsub(gtu33,JacPDstandardNth33gt23,knmsub(gtu12,kadd(JacPDstandardNth21gt23,JacPDstandardNth12gt23),knmsub(gtu13,kadd(JacPDstandardNth31gt23,JacPDstandardNth13gt23),knmsub(gtu23,kadd(JacPDstandardNth32gt23,JacPDstandardNth23gt23),kmadd(kmadd(Gt123,Gtlu133,kmadd(Gt113,Gtlu211,kmadd(Gt123,Gtlu212,kmadd(Gt133,Gtlu213,kmadd(Gt213,Gtlu221,kmadd(Gt223,Gtlu222,kmadd(Gt233,Gtlu223,kmadd(Gt212,Gtlu231,kmadd(Gt313,Gtlu231,kmadd(Gt222,Gtlu232,kmadd(Gt323,Gtlu232,kmadd(Gt223,Gtlu233,kmadd(Gt333,Gtlu233,kmadd(Gt112,kadd(Gtlu131,Gtlu311),kmadd(Gt122,kadd(Gtlu132,Gtlu312),kmadd(Gt123,Gtlu313,kmadd(Gt212,Gtlu321,kmadd(Gt222,Gtlu322,kmul(Gt223,Gtlu323))))))))))))))))))),ToReal(2),kmul(kmadd(Gt312,Gtlu331,kmadd(Gt322,Gtlu332,kmul(Gt323,Gtlu333))),ToReal(4))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt33 = 
      kmul(ToReal(0.5),knmsub(gtu11,JacPDstandardNth11gt33,knmsub(gtu22,JacPDstandardNth22gt33,knmsub(gtu33,JacPDstandardNth33gt33,knmsub(gtu12,kadd(JacPDstandardNth21gt33,JacPDstandardNth12gt33),knmsub(gtu13,kadd(JacPDstandardNth31gt33,JacPDstandardNth13gt33),knmsub(gtu23,kadd(JacPDstandardNth32gt33,JacPDstandardNth23gt33),kmadd(gt23L,kmul(JacPDstandardNth3Xt2,ToReal(2)),kmadd(gt33L,kmul(JacPDstandardNth3Xt3,ToReal(2)),kmadd(Gtl313,kmul(Xtn1,ToReal(2)),kmadd(Gtl323,kmul(Xtn2,ToReal(2)),kmadd(Gtl333,kmul(Xtn3,ToReal(2)),kmadd(ToReal(2),kmadd(Gt133,Gtlu133,kmadd(Gt213,Gtlu231,kmadd(Gt223,Gtlu232,kmadd(Gt233,Gtlu233,kmadd(gt13L,JacPDstandardNth3Xt1,kmadd(Gt113,kmadd(Gtlu311,ToReal(2),Gtlu131),kmul(Gt123,kmadd(Gtlu312,ToReal(2),Gtlu132)))))))),kmadd(kmadd(Gt133,Gtlu313,kmadd(Gt213,Gtlu321,kmadd(Gt223,Gtlu322,kmul(Gt233,Gtlu323)))),ToReal(4),kmul(kmadd(Gt313,Gtlu331,kmadd(Gt323,Gtlu332,kmul(Gt333,Gtlu333))),ToReal(6))))))))))))))));
    
    CCTK_REAL_VEC fac1 = 
      IfThen(conformalMethod,kmul(INV(phiL),ToReal(-0.5)),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 = kmul(fac1,JacPDstandardNth1phi);
    
    CCTK_REAL_VEC cdphi2 = kmul(fac1,JacPDstandardNth2phi);
    
    CCTK_REAL_VEC cdphi3 = kmul(fac1,JacPDstandardNth3phi);
    
    CCTK_REAL_VEC fac2 = 
      IfThen(conformalMethod,kmul(INV(SQR(phiL)),ToReal(0.5)),ToReal(0));
    
    CCTK_REAL_VEC cdphi211 = 
      kmadd(fac2,SQR(JacPDstandardNth1phi),kmul(fac1,ksub(JacPDstandardNth11phi,kmadd(Gt111,JacPDstandardNth1phi,kmadd(Gt311,JacPDstandardNth3phi,kmul(Gt211,JacPDstandardNth2phi))))));
    
    CCTK_REAL_VEC cdphi212 = 
      kmadd(fac2,kmul(JacPDstandardNth1phi,JacPDstandardNth2phi),kmul(fac1,ksub(JacPDstandardNth12phi,kmadd(Gt112,JacPDstandardNth1phi,kmadd(Gt312,JacPDstandardNth3phi,kmul(Gt212,JacPDstandardNth2phi))))));
    
    CCTK_REAL_VEC cdphi213 = 
      kmadd(fac2,kmul(JacPDstandardNth1phi,JacPDstandardNth3phi),kmul(fac1,ksub(JacPDstandardNth13phi,kmadd(Gt113,JacPDstandardNth1phi,kmadd(Gt313,JacPDstandardNth3phi,kmul(Gt213,JacPDstandardNth2phi))))));
    
    CCTK_REAL_VEC cdphi222 = 
      kmsub(fac2,SQR(JacPDstandardNth2phi),kmul(fac1,kmadd(Gt122,JacPDstandardNth1phi,kmadd(Gt222,JacPDstandardNth2phi,kmsub(Gt322,JacPDstandardNth3phi,JacPDstandardNth22phi)))));
    
    CCTK_REAL_VEC cdphi223 = 
      kmsub(fac2,kmul(JacPDstandardNth2phi,JacPDstandardNth3phi),kmul(fac1,kmadd(Gt123,JacPDstandardNth1phi,kmadd(Gt223,JacPDstandardNth2phi,kmsub(Gt323,JacPDstandardNth3phi,JacPDstandardNth23phi)))));
    
    CCTK_REAL_VEC cdphi233 = 
      kmsub(fac2,SQR(JacPDstandardNth3phi),kmul(fac1,kmadd(Gt133,JacPDstandardNth1phi,kmadd(Gt233,JacPDstandardNth2phi,kmsub(Gt333,JacPDstandardNth3phi,JacPDstandardNth33phi)))));
    
    CCTK_REAL_VEC Rphi11 = 
      kmul(ToReal(-2),kadd(cdphi211,kmadd(SQR(cdphi1),kmul(kmadd(gt11L,gtu11,ToReal(-1)),ToReal(2)),kmul(gt11L,kmadd(cdphi211,gtu11,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu33,SQR(cdphi3))))),ToReal(2),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmul(kmadd(cdphi1,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13)),kmul(cdphi2,kmul(cdphi3,gtu23))),ToReal(4))))))))));
    
    CCTK_REAL_VEC Rphi12 = 
      kmul(ToReal(-2),kadd(cdphi212,kmadd(gt12L,kmadd(cdphi211,gtu11,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu11,SQR(cdphi1))))),ToReal(2),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmadd(gtu33,kmadd(SQR(cdphi3),ToReal(2),cdphi233),kmul(cdphi2,kmul(cdphi3,kmul(gtu23,ToReal(4)))))))),kmul(cdphi1,kmadd(gt12L,kmul(cdphi3,kmul(gtu13,ToReal(4))),kmul(cdphi2,kmadd(gt12L,kmul(gtu12,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi13 = 
      kmul(ToReal(-2),kadd(cdphi213,kmadd(gt13L,kmadd(cdphi211,gtu11,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu11,SQR(cdphi1))))),ToReal(2),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmadd(gtu33,kmadd(SQR(cdphi3),ToReal(2),cdphi233),kmul(cdphi2,kmul(cdphi3,kmul(gtu23,ToReal(4)))))))),kmul(cdphi1,kmadd(gt13L,kmul(cdphi2,kmul(gtu12,ToReal(4))),kmul(cdphi3,kmadd(gt13L,kmul(gtu13,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi22 = 
      kmul(ToReal(-2),kadd(cdphi222,kmadd(SQR(cdphi2),kmul(kmadd(gt22L,gtu22,ToReal(-1)),ToReal(2)),kmul(gt22L,kmadd(cdphi222,gtu22,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu33,SQR(cdphi3))))),ToReal(2),kmadd(gtu11,kmadd(SQR(cdphi1),ToReal(2),cdphi211),kmul(kmadd(cdphi1,kmul(cdphi3,gtu13),kmul(cdphi2,kmadd(cdphi1,gtu12,kmul(cdphi3,gtu23)))),ToReal(4))))))))));
    
    CCTK_REAL_VEC Rphi23 = 
      kmul(ToReal(-2),kadd(cdphi223,kmadd(gt23L,kmadd(cdphi222,gtu22,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu22,SQR(cdphi2))))),ToReal(2),kmadd(gtu11,kmadd(SQR(cdphi1),ToReal(2),cdphi211),kmadd(gtu33,kmadd(SQR(cdphi3),ToReal(2),cdphi233),kmul(cdphi1,kmul(cdphi3,kmul(gtu13,ToReal(4)))))))),kmul(cdphi2,kmadd(gt23L,kmul(cdphi1,kmul(gtu12,ToReal(4))),kmul(cdphi3,kmadd(gt23L,kmul(gtu23,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi33 = 
      kmul(ToReal(-2),kadd(cdphi233,kmadd(SQR(cdphi3),kmul(kmadd(gt33L,gtu33,ToReal(-1)),ToReal(2)),kmul(gt33L,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi213,gtu13,kmul(cdphi223,gtu23)),ToReal(2),kmadd(gtu11,kmadd(SQR(cdphi1),ToReal(2),cdphi211),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmadd(cdphi3,kmul(kmadd(cdphi1,gtu13,kmul(cdphi2,gtu23)),ToReal(4)),kmul(gtu12,kmadd(cdphi212,ToReal(2),kmul(cdphi1,kmul(cdphi2,ToReal(4))))))))))))));
    
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
    
    CCTK_REAL_VEC e4phi = 
      IfThen(conformalMethod,INV(SQR(phiL)),kexp(kmul(phiL,ToReal(4))));
    
    CCTK_REAL_VEC em4phi = INV(e4phi);
    
    CCTK_REAL_VEC g11 = kmul(gt11L,e4phi);
    
    CCTK_REAL_VEC g12 = kmul(gt12L,e4phi);
    
    CCTK_REAL_VEC g13 = kmul(gt13L,e4phi);
    
    CCTK_REAL_VEC g22 = kmul(gt22L,e4phi);
    
    CCTK_REAL_VEC g23 = kmul(gt23L,e4phi);
    
    CCTK_REAL_VEC g33 = kmul(gt33L,e4phi);
    
    CCTK_REAL_VEC gu11 = kmul(em4phi,gtu11);
    
    CCTK_REAL_VEC gu12 = kmul(em4phi,gtu12);
    
    CCTK_REAL_VEC gu13 = kmul(em4phi,gtu13);
    
    CCTK_REAL_VEC gu22 = kmul(em4phi,gtu22);
    
    CCTK_REAL_VEC gu23 = kmul(em4phi,gtu23);
    
    CCTK_REAL_VEC gu33 = kmul(em4phi,gtu33);
    
    CCTK_REAL_VEC R11 = kadd(Rphi11,Rt11);
    
    CCTK_REAL_VEC R12 = kadd(Rphi12,Rt12);
    
    CCTK_REAL_VEC R13 = kadd(Rphi13,Rt13);
    
    CCTK_REAL_VEC R22 = kadd(Rphi22,Rt22);
    
    CCTK_REAL_VEC R23 = kadd(Rphi23,Rt23);
    
    CCTK_REAL_VEC R33 = kadd(Rphi33,Rt33);
    
    CCTK_REAL_VEC trS = 
      kmul(em4phi,kmadd(eTxxL,gtu11,kmadd(eTyyL,gtu22,kmadd(eTzzL,gtu33,kmul(kmadd(eTxyL,gtu12,kmadd(eTxzL,gtu13,kmul(eTyzL,gtu23))),ToReal(2))))));
    
    CCTK_REAL_VEC Ats11 = 
      kmadd(Gt211,JacPDstandardNth2alpha,kmadd(Gt311,JacPDstandardNth3alpha,kmadd(alphaL,R11,kmsub(JacPDstandardNth1alpha,kmadd(cdphi1,ToReal(4),Gt111),JacPDstandardNth11alpha))));
    
    CCTK_REAL_VEC Ats12 = 
      kmadd(Gt312,JacPDstandardNth3alpha,kmadd(alphaL,R12,ksub(kmadd(JacPDstandardNth2alpha,kmadd(cdphi1,ToReal(2),Gt212),kmul(JacPDstandardNth1alpha,kmadd(cdphi2,ToReal(2),Gt112))),JacPDstandardNth12alpha)));
    
    CCTK_REAL_VEC Ats13 = 
      kmadd(Gt213,JacPDstandardNth2alpha,kmadd(alphaL,R13,ksub(kmadd(JacPDstandardNth3alpha,kmadd(cdphi1,ToReal(2),Gt313),kmul(JacPDstandardNth1alpha,kmadd(cdphi3,ToReal(2),Gt113))),JacPDstandardNth13alpha)));
    
    CCTK_REAL_VEC Ats22 = 
      kmadd(Gt122,JacPDstandardNth1alpha,kmadd(Gt322,JacPDstandardNth3alpha,kmadd(alphaL,R22,kmsub(JacPDstandardNth2alpha,kmadd(cdphi2,ToReal(4),Gt222),JacPDstandardNth22alpha))));
    
    CCTK_REAL_VEC Ats23 = 
      kmadd(Gt123,JacPDstandardNth1alpha,kmadd(alphaL,R23,ksub(kmadd(JacPDstandardNth3alpha,kmadd(cdphi2,ToReal(2),Gt323),kmul(JacPDstandardNth2alpha,kmadd(cdphi3,ToReal(2),Gt223))),JacPDstandardNth23alpha)));
    
    CCTK_REAL_VEC Ats33 = 
      kmadd(Gt133,JacPDstandardNth1alpha,kmadd(Gt233,JacPDstandardNth2alpha,kmadd(alphaL,R33,kmsub(JacPDstandardNth3alpha,kmadd(cdphi3,ToReal(4),Gt333),JacPDstandardNth33alpha))));
    
    CCTK_REAL_VEC trAts = 
      kmadd(Ats11,gu11,kmadd(Ats22,gu22,kmadd(Ats33,gu33,kmul(kmadd(Ats12,gu12,kmadd(Ats13,gu13,kmul(Ats23,gu23))),ToReal(2)))));
    
    CCTK_REAL_VEC At11rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(em4phi,kmsub(Ats11,ToReal(3),kmul(g11,trAts)),kmadd(At11L,kmadd(kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3),ToReal(-2),kmul(JacPDstandardNth1beta1,ToReal(4))),kmsub(kmadd(At12L,JacPDstandardNth1beta2,kmul(At13L,JacPDstandardNth1beta3)),ToReal(6),kmul(alphaL,kmadd(kmadd(At12L,Atm21,kmul(At13L,Atm31)),ToReal(6),kmadd(At11L,kmadd(trKL,ToReal(-3),kmul(Atm11,ToReal(6))),kmul(em4phi,kmul(Pi,kmadd(g11,kmul(trS,ToReal(-8)),kmul(eTxxL,ToReal(24))))))))))));
    
    CCTK_REAL_VEC At12rhsL = 
      kmadd(ToReal(0.333333333333333333333333333333),kmadd(At12L,kadd(JacPDstandardNth1beta1,kmadd(JacPDstandardNth3beta3,ToReal(-2),JacPDstandardNth2beta2)),kmsub(kmadd(Ats12,em4phi,kmadd(At22L,JacPDstandardNth1beta2,kmadd(At23L,JacPDstandardNth1beta3,kmadd(At11L,JacPDstandardNth2beta1,kmul(At13L,JacPDstandardNth2beta3))))),ToReal(3),kmul(em4phi,kmul(g12,trAts)))),kmul(alphaL,kmadd(kmadd(At11L,Atm12,kmul(At13L,Atm32)),ToReal(-2),kmadd(At12L,kmadd(Atm22,ToReal(-2),trKL),kmul(em4phi,kmadd(eTxyL,kmul(Pi,ToReal(-8)),kmul(g12,kmul(trS,ToReal(8.37758040957278196923371568875)))))))));
    
    CCTK_REAL_VEC At13rhsL = 
      kmadd(ToReal(0.333333333333333333333333333333),kmadd(At13L,kadd(JacPDstandardNth1beta1,kmadd(JacPDstandardNth2beta2,ToReal(-2),JacPDstandardNth3beta3)),kmsub(kmadd(Ats13,em4phi,kmadd(At23L,JacPDstandardNth1beta2,kmadd(At33L,JacPDstandardNth1beta3,kmadd(At11L,JacPDstandardNth3beta1,kmul(At12L,JacPDstandardNth3beta2))))),ToReal(3),kmul(em4phi,kmul(g13,trAts)))),kmul(alphaL,kmadd(kmadd(At11L,Atm13,kmul(At12L,Atm23)),ToReal(-2),kmadd(At13L,kmadd(Atm33,ToReal(-2),trKL),kmul(em4phi,kmadd(eTxzL,kmul(Pi,ToReal(-8)),kmul(g13,kmul(trS,ToReal(8.37758040957278196923371568875)))))))));
    
    CCTK_REAL_VEC At22rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(em4phi,kmsub(Ats22,ToReal(3),kmul(g22,trAts)),kmadd(At22L,kmadd(kadd(JacPDstandardNth1beta1,JacPDstandardNth3beta3),ToReal(-2),kmul(JacPDstandardNth2beta2,ToReal(4))),kmsub(kmadd(At12L,JacPDstandardNth2beta1,kmul(At23L,JacPDstandardNth2beta3)),ToReal(6),kmul(alphaL,kmadd(kmadd(At12L,Atm12,kmul(At23L,Atm32)),ToReal(6),kmadd(At22L,kmadd(trKL,ToReal(-3),kmul(Atm22,ToReal(6))),kmul(em4phi,kmul(Pi,kmadd(g22,kmul(trS,ToReal(-8)),kmul(eTyyL,ToReal(24))))))))))));
    
    CCTK_REAL_VEC At23rhsL = 
      kmadd(ToReal(0.333333333333333333333333333333),kmadd(At23L,kadd(JacPDstandardNth2beta2,kmadd(JacPDstandardNth1beta1,ToReal(-2),JacPDstandardNth3beta3)),kmsub(kmadd(Ats23,em4phi,kmadd(At13L,JacPDstandardNth2beta1,kmadd(At33L,JacPDstandardNth2beta3,kmadd(At12L,JacPDstandardNth3beta1,kmul(At22L,JacPDstandardNth3beta2))))),ToReal(3),kmul(em4phi,kmul(g23,trAts)))),kmul(alphaL,kmadd(kmadd(At12L,Atm13,kmul(At22L,Atm23)),ToReal(-2),kmadd(At23L,kmadd(Atm33,ToReal(-2),trKL),kmul(em4phi,kmadd(eTyzL,kmul(Pi,ToReal(-8)),kmul(g23,kmul(trS,ToReal(8.37758040957278196923371568875)))))))));
    
    CCTK_REAL_VEC At33rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(em4phi,kmsub(Ats33,ToReal(3),kmul(g33,trAts)),kmadd(At33L,kmadd(kadd(JacPDstandardNth1beta1,JacPDstandardNth2beta2),ToReal(-2),kmul(JacPDstandardNth3beta3,ToReal(4))),kmsub(kmadd(At13L,JacPDstandardNth3beta1,kmul(At23L,JacPDstandardNth3beta2)),ToReal(6),kmul(alphaL,kmadd(kmadd(At13L,Atm13,kmul(At23L,Atm23)),ToReal(6),kmadd(At33L,kmadd(trKL,ToReal(-3),kmul(Atm33,ToReal(6))),kmul(em4phi,kmul(Pi,kmadd(g33,kmul(trS,ToReal(-8)),kmul(eTzzL,ToReal(24))))))))))));
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(At11rhs[index],At11rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At12rhs[index],At12rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At13rhs[index],At13rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At22rhs[index],At22rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At23rhs[index],At23rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At33rhs[index],At33rhsL,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(At11rhs[index],At11rhsL,elt_count);
      vec_store_nta_partial_hi(At12rhs[index],At12rhsL,elt_count);
      vec_store_nta_partial_hi(At13rhs[index],At13rhsL,elt_count);
      vec_store_nta_partial_hi(At22rhs[index],At22rhsL,elt_count);
      vec_store_nta_partial_hi(At23rhs[index],At23rhsL,elt_count);
      vec_store_nta_partial_hi(At33rhs[index],At33rhsL,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(At11rhs[index],At11rhsL,elt_count);
      vec_store_nta_partial_lo(At12rhs[index],At12rhsL,elt_count);
      vec_store_nta_partial_lo(At13rhs[index],At13rhsL,elt_count);
      vec_store_nta_partial_lo(At22rhs[index],At22rhsL,elt_count);
      vec_store_nta_partial_lo(At23rhs[index],At23rhsL,elt_count);
      vec_store_nta_partial_lo(At33rhs[index],At33rhsL,elt_count);
      break;
    }
    vec_store_nta(At11rhs[index],At11rhsL);
    vec_store_nta(At12rhs[index],At12rhsL);
    vec_store_nta(At13rhs[index],At13rhsL);
    vec_store_nta(At22rhs[index],At22rhsL);
    vec_store_nta(At23rhs[index],At23rhsL);
    vec_store_nta(At33rhs[index],At33rhsL);
  }
  LC_ENDLOOP3VEC (ML_BSSN_UPW_RHS2);
}

extern "C" void ML_BSSN_UPW_RHS2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_UPW_RHS2_Body");
  }
  
  if (cctk_iteration % ML_BSSN_UPW_RHS2_calc_every != ML_BSSN_UPW_RHS2_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_BSSN_UPW::ML_curv","ML_BSSN_UPW::ML_curvrhs","ML_BSSN_UPW::ML_Gamma","ML_BSSN_UPW::ML_lapse","ML_BSSN_UPW::ML_log_confac","ML_BSSN_UPW::ML_metric","ML_BSSN_UPW::ML_shift","ML_BSSN_UPW::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_UPW_RHS2", 8, groups);
  
  switch(fdOrder)
  {
    case 2:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_RHS2", 1, 1, 1);
      break;
    
    case 4:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_RHS2", 2, 2, 2);
      break;
    
    case 6:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_RHS2", 3, 3, 3);
      break;
    
    case 8:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_RHS2", 4, 4, 4);
      break;
  }
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_UPW_RHS2_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_UPW_RHS2_Body");
  }
}
