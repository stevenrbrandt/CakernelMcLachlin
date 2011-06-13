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

extern "C" void ML_BSSN_MP_O8_convertToADMBaseDtLapseShift_SelectBCs(CCTK_ARGUMENTS)
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

static void ML_BSSN_MP_O8_convertToADMBaseDtLapseShift_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_O8_convertToADMBaseDtLapseShift_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_O8_convertToADMBaseDtLapseShift_calc_every != ML_BSSN_MP_O8_convertToADMBaseDtLapseShift_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::dtlapse","ADMBase::dtshift","grid::coordinates","Grid::coordinates","ML_BSSN_MP_O8::ML_dtlapse","ML_BSSN_MP_O8::ML_dtshift","ML_BSSN_MP_O8::ML_Gamma","ML_BSSN_MP_O8::ML_lapse","ML_BSSN_MP_O8::ML_shift","ML_BSSN_MP_O8::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_O8_convertToADMBaseDtLapseShift", 10, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_MP_O8_convertToADMBaseDtLapseShift", 5, 5, 5);
  
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
  LC_LOOP3 (ML_BSSN_MP_O8_convertToADMBaseDtLapseShift,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL AL = A[index];
    CCTK_REAL alphaL = alpha[index];
    CCTK_REAL B1L = B1[index];
    CCTK_REAL B2L = B2[index];
    CCTK_REAL B3L = B3[index];
    CCTK_REAL beta1L = beta1[index];
    CCTK_REAL beta2L = beta2[index];
    CCTK_REAL beta3L = beta3[index];
    CCTK_REAL rL = r[index];
    CCTK_REAL trKL = trK[index];
    CCTK_REAL Xt1L = Xt1[index];
    CCTK_REAL Xt2L = Xt2[index];
    CCTK_REAL Xt3L = Xt3[index];
    
    
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
    CCTK_REAL const PDupwindNthAnti1alpha = PDupwindNthAnti1(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm1alpha = PDupwindNthSymm1(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti2alpha = PDupwindNthAnti2(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm2alpha = PDupwindNthSymm2(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti3alpha = PDupwindNthAnti3(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm3alpha = PDupwindNthSymm3(&alpha[index]);
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
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL JacPDupwindNthAnti1alpha;
    CCTK_REAL JacPDupwindNthAnti1beta1;
    CCTK_REAL JacPDupwindNthAnti1beta2;
    CCTK_REAL JacPDupwindNthAnti1beta3;
    CCTK_REAL JacPDupwindNthAnti2alpha;
    CCTK_REAL JacPDupwindNthAnti2beta1;
    CCTK_REAL JacPDupwindNthAnti2beta2;
    CCTK_REAL JacPDupwindNthAnti2beta3;
    CCTK_REAL JacPDupwindNthAnti3alpha;
    CCTK_REAL JacPDupwindNthAnti3beta1;
    CCTK_REAL JacPDupwindNthAnti3beta2;
    CCTK_REAL JacPDupwindNthAnti3beta3;
    CCTK_REAL JacPDupwindNthSymm1alpha;
    CCTK_REAL JacPDupwindNthSymm1beta1;
    CCTK_REAL JacPDupwindNthSymm1beta2;
    CCTK_REAL JacPDupwindNthSymm1beta3;
    CCTK_REAL JacPDupwindNthSymm2alpha;
    CCTK_REAL JacPDupwindNthSymm2beta1;
    CCTK_REAL JacPDupwindNthSymm2beta2;
    CCTK_REAL JacPDupwindNthSymm2beta3;
    CCTK_REAL JacPDupwindNthSymm3alpha;
    CCTK_REAL JacPDupwindNthSymm3beta1;
    CCTK_REAL JacPDupwindNthSymm3beta2;
    CCTK_REAL JacPDupwindNthSymm3beta3;
    
    if (use_jacobian)
    {
      JacPDupwindNthAnti1alpha = J11L*PDupwindNthAnti1alpha + 
        J21L*PDupwindNthAnti2alpha + J31L*PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti1beta1 = J11L*PDupwindNthAnti1beta1 + 
        J21L*PDupwindNthAnti2beta1 + J31L*PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti1beta2 = J11L*PDupwindNthAnti1beta2 + 
        J21L*PDupwindNthAnti2beta2 + J31L*PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti1beta3 = J11L*PDupwindNthAnti1beta3 + 
        J21L*PDupwindNthAnti2beta3 + J31L*PDupwindNthAnti3beta3;
      
      JacPDupwindNthSymm1alpha = J11L*PDupwindNthSymm1alpha + 
        J21L*PDupwindNthSymm2alpha + J31L*PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm1beta1 = J11L*PDupwindNthSymm1beta1 + 
        J21L*PDupwindNthSymm2beta1 + J31L*PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm1beta2 = J11L*PDupwindNthSymm1beta2 + 
        J21L*PDupwindNthSymm2beta2 + J31L*PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm1beta3 = J11L*PDupwindNthSymm1beta3 + 
        J21L*PDupwindNthSymm2beta3 + J31L*PDupwindNthSymm3beta3;
      
      JacPDupwindNthAnti2alpha = J12L*PDupwindNthAnti1alpha + 
        J22L*PDupwindNthAnti2alpha + J32L*PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti2beta1 = J12L*PDupwindNthAnti1beta1 + 
        J22L*PDupwindNthAnti2beta1 + J32L*PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti2beta2 = J12L*PDupwindNthAnti1beta2 + 
        J22L*PDupwindNthAnti2beta2 + J32L*PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti2beta3 = J12L*PDupwindNthAnti1beta3 + 
        J22L*PDupwindNthAnti2beta3 + J32L*PDupwindNthAnti3beta3;
      
      JacPDupwindNthSymm2alpha = J12L*PDupwindNthSymm1alpha + 
        J22L*PDupwindNthSymm2alpha + J32L*PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm2beta1 = J12L*PDupwindNthSymm1beta1 + 
        J22L*PDupwindNthSymm2beta1 + J32L*PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm2beta2 = J12L*PDupwindNthSymm1beta2 + 
        J22L*PDupwindNthSymm2beta2 + J32L*PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm2beta3 = J12L*PDupwindNthSymm1beta3 + 
        J22L*PDupwindNthSymm2beta3 + J32L*PDupwindNthSymm3beta3;
      
      JacPDupwindNthAnti3alpha = J13L*PDupwindNthAnti1alpha + 
        J23L*PDupwindNthAnti2alpha + J33L*PDupwindNthAnti3alpha;
      
      JacPDupwindNthAnti3beta1 = J13L*PDupwindNthAnti1beta1 + 
        J23L*PDupwindNthAnti2beta1 + J33L*PDupwindNthAnti3beta1;
      
      JacPDupwindNthAnti3beta2 = J13L*PDupwindNthAnti1beta2 + 
        J23L*PDupwindNthAnti2beta2 + J33L*PDupwindNthAnti3beta2;
      
      JacPDupwindNthAnti3beta3 = J13L*PDupwindNthAnti1beta3 + 
        J23L*PDupwindNthAnti2beta3 + J33L*PDupwindNthAnti3beta3;
      
      JacPDupwindNthSymm3alpha = J13L*PDupwindNthSymm1alpha + 
        J23L*PDupwindNthSymm2alpha + J33L*PDupwindNthSymm3alpha;
      
      JacPDupwindNthSymm3beta1 = J13L*PDupwindNthSymm1beta1 + 
        J23L*PDupwindNthSymm2beta1 + J33L*PDupwindNthSymm3beta1;
      
      JacPDupwindNthSymm3beta2 = J13L*PDupwindNthSymm1beta2 + 
        J23L*PDupwindNthSymm2beta2 + J33L*PDupwindNthSymm3beta2;
      
      JacPDupwindNthSymm3beta3 = J13L*PDupwindNthSymm1beta3 + 
        J23L*PDupwindNthSymm2beta3 + J33L*PDupwindNthSymm3beta3;
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
    
    CCTK_REAL eta = fmin(1,INV(rL)*ToReal(SpatialBetaDriverRadius));
    
    CCTK_REAL theta = fmin(1,exp(1 - 
      rL*INV(ToReal(SpatialShiftGammaCoeffRadius))));
    
    CCTK_REAL dtalpL = 
      -(pow(alphaL,ToReal(harmonicN))*ToReal(harmonicF)*(trKL + (AL - 
      trKL)*ToReal(LapseACoeff))) + (beta1L*JacPDupwindNthAnti1alpha + 
      beta2L*JacPDupwindNthAnti2alpha + beta3L*JacPDupwindNthAnti3alpha + 
      JacPDupwindNthSymm1alpha*Abs(beta1L) + 
      JacPDupwindNthSymm2alpha*Abs(beta2L) + 
      JacPDupwindNthSymm3alpha*Abs(beta3L))*ToReal(LapseAdvectionCoeff);
    
    CCTK_REAL dtbetaxL = (beta1L*JacPDupwindNthAnti1beta1 + 
      beta2L*JacPDupwindNthAnti2beta1 + beta3L*JacPDupwindNthAnti3beta1 + 
      JacPDupwindNthSymm1beta1*Abs(beta1L) + 
      JacPDupwindNthSymm2beta1*Abs(beta2L) + 
      JacPDupwindNthSymm3beta1*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      theta*(Xt1L + beta1L*eta*ToReal(BetaDriver)*(-1 + ToReal(ShiftBCoeff)) 
      + (B1L - Xt1L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    CCTK_REAL dtbetayL = (beta1L*JacPDupwindNthAnti1beta2 + 
      beta2L*JacPDupwindNthAnti2beta2 + beta3L*JacPDupwindNthAnti3beta2 + 
      JacPDupwindNthSymm1beta2*Abs(beta1L) + 
      JacPDupwindNthSymm2beta2*Abs(beta2L) + 
      JacPDupwindNthSymm3beta2*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      theta*(Xt2L + beta2L*eta*ToReal(BetaDriver)*(-1 + ToReal(ShiftBCoeff)) 
      + (B2L - Xt2L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    CCTK_REAL dtbetazL = (beta1L*JacPDupwindNthAnti1beta3 + 
      beta2L*JacPDupwindNthAnti2beta3 + beta3L*JacPDupwindNthAnti3beta3 + 
      JacPDupwindNthSymm1beta3*Abs(beta1L) + 
      JacPDupwindNthSymm2beta3*Abs(beta2L) + 
      JacPDupwindNthSymm3beta3*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      theta*(Xt3L + beta3L*eta*ToReal(BetaDriver)*(-1 + ToReal(ShiftBCoeff)) 
      + (B3L - Xt3L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    /* Copy local copies back to grid functions */
    dtalp[index] = dtalpL;
    dtbetax[index] = dtbetaxL;
    dtbetay[index] = dtbetayL;
    dtbetaz[index] = dtbetazL;
  }
  LC_ENDLOOP3 (ML_BSSN_MP_O8_convertToADMBaseDtLapseShift);
}

extern "C" void ML_BSSN_MP_O8_convertToADMBaseDtLapseShift(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_O8_convertToADMBaseDtLapseShift_Body);
}
