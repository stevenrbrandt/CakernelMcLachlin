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

extern "C" void ML_BSSN_O8_convertToADMBaseDtLapseShift_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ADMBase::dtlapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ADMBase::dtshift.");
  return;
}

static void ML_BSSN_O8_convertToADMBaseDtLapseShift_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O8_convertToADMBaseDtLapseShift_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O8_convertToADMBaseDtLapseShift_calc_every != ML_BSSN_O8_convertToADMBaseDtLapseShift_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::dtlapse","ADMBase::dtshift","grid::coordinates","Grid::coordinates","ML_BSSN_O8::ML_dtlapse","ML_BSSN_O8::ML_dtshift","ML_BSSN_O8::ML_Gamma","ML_BSSN_O8::ML_lapse","ML_BSSN_O8::ML_shift","ML_BSSN_O8::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_O8_convertToADMBaseDtLapseShift", 10, groups);
  
  switch(fdOrder)
  {
    case 2:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_O8_convertToADMBaseDtLapseShift", 2, 2, 2);
      break;
    
    case 4:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_O8_convertToADMBaseDtLapseShift", 3, 3, 3);
      break;
    
    case 6:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_O8_convertToADMBaseDtLapseShift", 4, 4, 4);
      break;
    
    case 8:
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_O8_convertToADMBaseDtLapseShift", 5, 5, 5);
      break;
  }
  
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
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC (ML_BSSN_O8_convertToADMBaseDtLapseShift,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC AL = vec_load(A[index]);
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC B1L = vec_load(B1[index]);
    CCTK_REAL_VEC B2L = vec_load(B2[index]);
    CCTK_REAL_VEC B3L = vec_load(B3[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC rL = vec_load(r[index]);
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    CCTK_REAL_VEC Xt1L = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L = vec_load(Xt3[index]);
    
    
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
    CCTK_REAL_VEC PDupwindNthAnti1alpha;
    CCTK_REAL_VEC PDupwindNthSymm1alpha;
    CCTK_REAL_VEC PDupwindNthAnti2alpha;
    CCTK_REAL_VEC PDupwindNthSymm2alpha;
    CCTK_REAL_VEC PDupwindNthAnti3alpha;
    CCTK_REAL_VEC PDupwindNthSymm3alpha;
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
    
    switch(fdOrder)
    {
      case 2:
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder21(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder21(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder22(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder22(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder23(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder23(&alpha[index]);
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
        break;
      
      case 4:
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder41(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder41(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder42(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder42(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder43(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder43(&alpha[index]);
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
        break;
      
      case 6:
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder61(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder61(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder62(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder62(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder63(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder63(&alpha[index]);
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
        break;
      
      case 8:
        PDupwindNthAnti1alpha = PDupwindNthAntifdOrder81(&alpha[index]);
        PDupwindNthSymm1alpha = PDupwindNthSymmfdOrder81(&alpha[index]);
        PDupwindNthAnti2alpha = PDupwindNthAntifdOrder82(&alpha[index]);
        PDupwindNthSymm2alpha = PDupwindNthSymmfdOrder82(&alpha[index]);
        PDupwindNthAnti3alpha = PDupwindNthAntifdOrder83(&alpha[index]);
        PDupwindNthSymm3alpha = PDupwindNthSymmfdOrder83(&alpha[index]);
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
        break;
    }
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC JacPDupwindNthAnti1alpha;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta1;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta2;
    CCTK_REAL_VEC JacPDupwindNthAnti1beta3;
    CCTK_REAL_VEC JacPDupwindNthAnti2alpha;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta1;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta2;
    CCTK_REAL_VEC JacPDupwindNthAnti2beta3;
    CCTK_REAL_VEC JacPDupwindNthAnti3alpha;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta1;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta2;
    CCTK_REAL_VEC JacPDupwindNthAnti3beta3;
    CCTK_REAL_VEC JacPDupwindNthSymm1alpha;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta1;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta2;
    CCTK_REAL_VEC JacPDupwindNthSymm1beta3;
    CCTK_REAL_VEC JacPDupwindNthSymm2alpha;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta1;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta2;
    CCTK_REAL_VEC JacPDupwindNthSymm2beta3;
    CCTK_REAL_VEC JacPDupwindNthSymm3alpha;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta1;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta2;
    CCTK_REAL_VEC JacPDupwindNthSymm3beta3;
    
    if (use_jacobian)
    {
      JacPDupwindNthAnti1alpha = 
        kmadd(J11L,PDupwindNthAnti1alpha,kmadd(J21L,PDupwindNthAnti2alpha,kmul(J31L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti1beta1 = 
        kmadd(J11L,PDupwindNthAnti1beta1,kmadd(J21L,PDupwindNthAnti2beta1,kmul(J31L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti1beta2 = 
        kmadd(J11L,PDupwindNthAnti1beta2,kmadd(J21L,PDupwindNthAnti2beta2,kmul(J31L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti1beta3 = 
        kmadd(J11L,PDupwindNthAnti1beta3,kmadd(J21L,PDupwindNthAnti2beta3,kmul(J31L,PDupwindNthAnti3beta3)));
      
      JacPDupwindNthSymm1alpha = 
        kmadd(J11L,PDupwindNthSymm1alpha,kmadd(J21L,PDupwindNthSymm2alpha,kmul(J31L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm1beta1 = 
        kmadd(J11L,PDupwindNthSymm1beta1,kmadd(J21L,PDupwindNthSymm2beta1,kmul(J31L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm1beta2 = 
        kmadd(J11L,PDupwindNthSymm1beta2,kmadd(J21L,PDupwindNthSymm2beta2,kmul(J31L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm1beta3 = 
        kmadd(J11L,PDupwindNthSymm1beta3,kmadd(J21L,PDupwindNthSymm2beta3,kmul(J31L,PDupwindNthSymm3beta3)));
      
      JacPDupwindNthAnti2alpha = 
        kmadd(J12L,PDupwindNthAnti1alpha,kmadd(J22L,PDupwindNthAnti2alpha,kmul(J32L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti2beta1 = 
        kmadd(J12L,PDupwindNthAnti1beta1,kmadd(J22L,PDupwindNthAnti2beta1,kmul(J32L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti2beta2 = 
        kmadd(J12L,PDupwindNthAnti1beta2,kmadd(J22L,PDupwindNthAnti2beta2,kmul(J32L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti2beta3 = 
        kmadd(J12L,PDupwindNthAnti1beta3,kmadd(J22L,PDupwindNthAnti2beta3,kmul(J32L,PDupwindNthAnti3beta3)));
      
      JacPDupwindNthSymm2alpha = 
        kmadd(J12L,PDupwindNthSymm1alpha,kmadd(J22L,PDupwindNthSymm2alpha,kmul(J32L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm2beta1 = 
        kmadd(J12L,PDupwindNthSymm1beta1,kmadd(J22L,PDupwindNthSymm2beta1,kmul(J32L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm2beta2 = 
        kmadd(J12L,PDupwindNthSymm1beta2,kmadd(J22L,PDupwindNthSymm2beta2,kmul(J32L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm2beta3 = 
        kmadd(J12L,PDupwindNthSymm1beta3,kmadd(J22L,PDupwindNthSymm2beta3,kmul(J32L,PDupwindNthSymm3beta3)));
      
      JacPDupwindNthAnti3alpha = 
        kmadd(J13L,PDupwindNthAnti1alpha,kmadd(J23L,PDupwindNthAnti2alpha,kmul(J33L,PDupwindNthAnti3alpha)));
      
      JacPDupwindNthAnti3beta1 = 
        kmadd(J13L,PDupwindNthAnti1beta1,kmadd(J23L,PDupwindNthAnti2beta1,kmul(J33L,PDupwindNthAnti3beta1)));
      
      JacPDupwindNthAnti3beta2 = 
        kmadd(J13L,PDupwindNthAnti1beta2,kmadd(J23L,PDupwindNthAnti2beta2,kmul(J33L,PDupwindNthAnti3beta2)));
      
      JacPDupwindNthAnti3beta3 = 
        kmadd(J13L,PDupwindNthAnti1beta3,kmadd(J23L,PDupwindNthAnti2beta3,kmul(J33L,PDupwindNthAnti3beta3)));
      
      JacPDupwindNthSymm3alpha = 
        kmadd(J13L,PDupwindNthSymm1alpha,kmadd(J23L,PDupwindNthSymm2alpha,kmul(J33L,PDupwindNthSymm3alpha)));
      
      JacPDupwindNthSymm3beta1 = 
        kmadd(J13L,PDupwindNthSymm1beta1,kmadd(J23L,PDupwindNthSymm2beta1,kmul(J33L,PDupwindNthSymm3beta1)));
      
      JacPDupwindNthSymm3beta2 = 
        kmadd(J13L,PDupwindNthSymm1beta2,kmadd(J23L,PDupwindNthSymm2beta2,kmul(J33L,PDupwindNthSymm3beta2)));
      
      JacPDupwindNthSymm3beta3 = 
        kmadd(J13L,PDupwindNthSymm1beta3,kmadd(J23L,PDupwindNthSymm2beta3,kmul(J33L,PDupwindNthSymm3beta3)));
    }
    else
    {
      JacPDupwindNthAnti1alpha = PDupwindNthAnti1alpha;
      
      JacPDupwindNthAnti1beta1 = PDupwindNthAnti1beta1;
      
      JacPDupwindNthAnti1beta2 = PDupwindNthAnti1beta2;
      
      JacPDupwindNthAnti1beta3 = PDupwindNthAnti1beta3;
      
      JacPDupwindNthSymm1alpha = PDupwindNthSymm1alpha;
      
      JacPDupwindNthSymm1beta1 = PDupwindNthSymm1beta1;
      
      JacPDupwindNthSymm1beta2 = PDupwindNthSymm1beta2;
      
      JacPDupwindNthSymm1beta3 = PDupwindNthSymm1beta3;
      
      JacPDupwindNthAnti2alpha = PDupwindNthAnti2alpha;
      
      JacPDupwindNthAnti2beta1 = PDupwindNthAnti2beta1;
      
      JacPDupwindNthAnti2beta2 = PDupwindNthAnti2beta2;
      
      JacPDupwindNthAnti2beta3 = PDupwindNthAnti2beta3;
      
      JacPDupwindNthSymm2alpha = PDupwindNthSymm2alpha;
      
      JacPDupwindNthSymm2beta1 = PDupwindNthSymm2beta1;
      
      JacPDupwindNthSymm2beta2 = PDupwindNthSymm2beta2;
      
      JacPDupwindNthSymm2beta3 = PDupwindNthSymm2beta3;
      
      JacPDupwindNthAnti3alpha = PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti3beta1 = PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti3beta2 = PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti3beta3 = PDupwindNthAnti3beta3;
      
      JacPDupwindNthSymm3alpha = PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm3beta1 = PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm3beta2 = PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm3beta3 = PDupwindNthSymm3beta3;
    }
    
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    CCTK_REAL_VEC eta = 
      kfmin(ToReal(1),kmul(INV(rL),ToReal(SpatialBetaDriverRadius)));
    
    CCTK_REAL_VEC theta = 
      kfmin(ToReal(1),kexp(knmsub(rL,INV(ToReal(SpatialShiftGammaCoeffRadius)),ToReal(1))));
    
    CCTK_REAL_VEC dtalpL = 
      kmsub(kmadd(beta1L,JacPDupwindNthAnti1alpha,kmadd(beta2L,JacPDupwindNthAnti2alpha,kmadd(beta3L,JacPDupwindNthAnti3alpha,kmadd(JacPDupwindNthSymm1alpha,kfabs(beta1L),kmadd(JacPDupwindNthSymm2alpha,kfabs(beta2L),kmul(JacPDupwindNthSymm3alpha,kfabs(beta3L))))))),ToReal(LapseAdvectionCoeff),kmul(kpow(alphaL,harmonicN),kmul(ToReal(harmonicF),kmadd(ksub(AL,trKL),ToReal(LapseACoeff),trKL))));
    
    CCTK_REAL_VEC dtbetaxL = 
      kmadd(kmadd(beta1L,JacPDupwindNthAnti1beta1,kmadd(beta2L,JacPDupwindNthAnti2beta1,kmadd(beta3L,JacPDupwindNthAnti3beta1,kmadd(JacPDupwindNthSymm1beta1,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta1,kfabs(beta2L),kmul(JacPDupwindNthSymm3beta1,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),kmul(theta,kmul(kadd(Xt1L,kmadd(beta1L,kmul(eta,ToReal(BetaDriver*(-1 
      + 
      ShiftBCoeff))),kmul(ksub(B1L,Xt1L),ToReal(ShiftBCoeff)))),ToReal(ShiftGammaCoeff))));
    
    CCTK_REAL_VEC dtbetayL = 
      kmadd(kmadd(beta1L,JacPDupwindNthAnti1beta2,kmadd(beta2L,JacPDupwindNthAnti2beta2,kmadd(beta3L,JacPDupwindNthAnti3beta2,kmadd(JacPDupwindNthSymm1beta2,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta2,kfabs(beta2L),kmul(JacPDupwindNthSymm3beta2,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),kmul(theta,kmul(kadd(Xt2L,kmadd(beta2L,kmul(eta,ToReal(BetaDriver*(-1 
      + 
      ShiftBCoeff))),kmul(ksub(B2L,Xt2L),ToReal(ShiftBCoeff)))),ToReal(ShiftGammaCoeff))));
    
    CCTK_REAL_VEC dtbetazL = 
      kmadd(kmadd(beta1L,JacPDupwindNthAnti1beta3,kmadd(beta2L,JacPDupwindNthAnti2beta3,kmadd(beta3L,JacPDupwindNthAnti3beta3,kmadd(JacPDupwindNthSymm1beta3,kfabs(beta1L),kmadd(JacPDupwindNthSymm2beta3,kfabs(beta2L),kmul(JacPDupwindNthSymm3beta3,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),kmul(theta,kmul(kadd(Xt3L,kmadd(beta3L,kmul(eta,ToReal(BetaDriver*(-1 
      + 
      ShiftBCoeff))),kmul(ksub(B3L,Xt3L),ToReal(ShiftBCoeff)))),ToReal(ShiftGammaCoeff))));
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(dtalp[index],dtalpL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(dtbetax[index],dtbetaxL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(dtbetay[index],dtbetayL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(dtbetaz[index],dtbetazL,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(dtalp[index],dtalpL,elt_count);
      vec_store_nta_partial_hi(dtbetax[index],dtbetaxL,elt_count);
      vec_store_nta_partial_hi(dtbetay[index],dtbetayL,elt_count);
      vec_store_nta_partial_hi(dtbetaz[index],dtbetazL,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(dtalp[index],dtalpL,elt_count);
      vec_store_nta_partial_lo(dtbetax[index],dtbetaxL,elt_count);
      vec_store_nta_partial_lo(dtbetay[index],dtbetayL,elt_count);
      vec_store_nta_partial_lo(dtbetaz[index],dtbetazL,elt_count);
      break;
    }
    
    /* Copy local copies back to grid functions */
    vec_store_nta(dtalp[index],dtalpL);
    vec_store_nta(dtbetax[index],dtbetaxL);
    vec_store_nta(dtbetay[index],dtbetayL);
    vec_store_nta(dtbetaz[index],dtbetazL);
  }
  LC_ENDLOOP3VEC (ML_BSSN_O8_convertToADMBaseDtLapseShift);
}

extern "C" void ML_BSSN_O8_convertToADMBaseDtLapseShift(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_O8_convertToADMBaseDtLapseShift_Body);
}
