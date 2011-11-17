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

extern "C" void ML_ADM_RHS_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_curvrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_shiftrhs.");
  return;
}

static void ML_ADM_RHS_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_REAL_VEC const p1o12dx = kmul(INV(dx),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o12dy = kmul(INV(dy),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o12dz = kmul(INV(dz),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o144dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o144dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o144dydz = kmul(INV(dy),kmul(INV(dz),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o180dx2 = kmul(INV(SQR(dx)),ToReal(0.00555555555555555555555555555556));
  CCTK_REAL_VEC const p1o180dy2 = kmul(INV(SQR(dy)),ToReal(0.00555555555555555555555555555556));
  CCTK_REAL_VEC const p1o180dz2 = kmul(INV(SQR(dz)),ToReal(0.00555555555555555555555555555556));
  CCTK_REAL_VEC const p1o2dx = kmul(INV(dx),ToReal(0.5));
  CCTK_REAL_VEC const p1o2dy = kmul(INV(dy),ToReal(0.5));
  CCTK_REAL_VEC const p1o2dz = kmul(INV(dz),ToReal(0.5));
  CCTK_REAL_VEC const p1o3600dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(0.000277777777777777777777777777778)));
  CCTK_REAL_VEC const p1o3600dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(0.000277777777777777777777777777778)));
  CCTK_REAL_VEC const p1o3600dydz = kmul(INV(dy),kmul(INV(dz),ToReal(0.000277777777777777777777777777778)));
  CCTK_REAL_VEC const p1o4dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(0.25)));
  CCTK_REAL_VEC const p1o4dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(0.25)));
  CCTK_REAL_VEC const p1o4dydz = kmul(INV(dy),kmul(INV(dz),ToReal(0.25)));
  CCTK_REAL_VEC const p1o5040dx2 = kmul(INV(SQR(dx)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dy2 = kmul(INV(SQR(dy)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dz2 = kmul(INV(SQR(dz)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o60dx = kmul(INV(dx),ToReal(0.0166666666666666666666666666667));
  CCTK_REAL_VEC const p1o60dy = kmul(INV(dy),ToReal(0.0166666666666666666666666666667));
  CCTK_REAL_VEC const p1o60dz = kmul(INV(dz),ToReal(0.0166666666666666666666666666667));
  CCTK_REAL_VEC const p1o705600dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o705600dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o705600dydz = kmul(INV(dy),kmul(INV(dz),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o840dx = kmul(INV(dx),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dy = kmul(INV(dy),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dz = kmul(INV(dz),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1odx2 = INV(SQR(dx));
  CCTK_REAL_VEC const p1ody2 = INV(SQR(dy));
  CCTK_REAL_VEC const p1odz2 = INV(SQR(dz));
  CCTK_REAL_VEC const pm1o12dx2 = kmul(INV(SQR(dx)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dy2 = kmul(INV(SQR(dy)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dz2 = kmul(INV(SQR(dz)),ToReal(-0.0833333333333333333333333333333));
  
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
  LC_LOOP3VEC (ML_ADM_RHS,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC g11L = vec_load(g11[index]);
    CCTK_REAL_VEC g12L = vec_load(g12[index]);
    CCTK_REAL_VEC g13L = vec_load(g13[index]);
    CCTK_REAL_VEC g22L = vec_load(g22[index]);
    CCTK_REAL_VEC g23L = vec_load(g23[index]);
    CCTK_REAL_VEC g33L = vec_load(g33[index]);
    CCTK_REAL_VEC K11L = vec_load(K11[index]);
    CCTK_REAL_VEC K12L = vec_load(K12[index]);
    CCTK_REAL_VEC K13L = vec_load(K13[index]);
    CCTK_REAL_VEC K22L = vec_load(K22[index]);
    CCTK_REAL_VEC K23L = vec_load(K23[index]);
    CCTK_REAL_VEC K33L = vec_load(K33[index]);
    
    
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
    CCTK_REAL_VEC PDstandardNth1g11;
    CCTK_REAL_VEC PDstandardNth2g11;
    CCTK_REAL_VEC PDstandardNth3g11;
    CCTK_REAL_VEC PDstandardNth11g11;
    CCTK_REAL_VEC PDstandardNth22g11;
    CCTK_REAL_VEC PDstandardNth33g11;
    CCTK_REAL_VEC PDstandardNth12g11;
    CCTK_REAL_VEC PDstandardNth13g11;
    CCTK_REAL_VEC PDstandardNth23g11;
    CCTK_REAL_VEC PDstandardNth1g12;
    CCTK_REAL_VEC PDstandardNth2g12;
    CCTK_REAL_VEC PDstandardNth3g12;
    CCTK_REAL_VEC PDstandardNth11g12;
    CCTK_REAL_VEC PDstandardNth22g12;
    CCTK_REAL_VEC PDstandardNth33g12;
    CCTK_REAL_VEC PDstandardNth12g12;
    CCTK_REAL_VEC PDstandardNth13g12;
    CCTK_REAL_VEC PDstandardNth23g12;
    CCTK_REAL_VEC PDstandardNth1g13;
    CCTK_REAL_VEC PDstandardNth2g13;
    CCTK_REAL_VEC PDstandardNth3g13;
    CCTK_REAL_VEC PDstandardNth11g13;
    CCTK_REAL_VEC PDstandardNth22g13;
    CCTK_REAL_VEC PDstandardNth33g13;
    CCTK_REAL_VEC PDstandardNth12g13;
    CCTK_REAL_VEC PDstandardNth13g13;
    CCTK_REAL_VEC PDstandardNth23g13;
    CCTK_REAL_VEC PDstandardNth1g22;
    CCTK_REAL_VEC PDstandardNth2g22;
    CCTK_REAL_VEC PDstandardNth3g22;
    CCTK_REAL_VEC PDstandardNth11g22;
    CCTK_REAL_VEC PDstandardNth22g22;
    CCTK_REAL_VEC PDstandardNth33g22;
    CCTK_REAL_VEC PDstandardNth12g22;
    CCTK_REAL_VEC PDstandardNth13g22;
    CCTK_REAL_VEC PDstandardNth23g22;
    CCTK_REAL_VEC PDstandardNth1g23;
    CCTK_REAL_VEC PDstandardNth2g23;
    CCTK_REAL_VEC PDstandardNth3g23;
    CCTK_REAL_VEC PDstandardNth11g23;
    CCTK_REAL_VEC PDstandardNth22g23;
    CCTK_REAL_VEC PDstandardNth33g23;
    CCTK_REAL_VEC PDstandardNth12g23;
    CCTK_REAL_VEC PDstandardNth13g23;
    CCTK_REAL_VEC PDstandardNth23g23;
    CCTK_REAL_VEC PDstandardNth1g33;
    CCTK_REAL_VEC PDstandardNth2g33;
    CCTK_REAL_VEC PDstandardNth3g33;
    CCTK_REAL_VEC PDstandardNth11g33;
    CCTK_REAL_VEC PDstandardNth22g33;
    CCTK_REAL_VEC PDstandardNth33g33;
    CCTK_REAL_VEC PDstandardNth12g33;
    CCTK_REAL_VEC PDstandardNth13g33;
    CCTK_REAL_VEC PDstandardNth23g33;
    CCTK_REAL_VEC PDstandardNth1K11;
    CCTK_REAL_VEC PDstandardNth2K11;
    CCTK_REAL_VEC PDstandardNth3K11;
    CCTK_REAL_VEC PDstandardNth1K12;
    CCTK_REAL_VEC PDstandardNth2K12;
    CCTK_REAL_VEC PDstandardNth3K12;
    CCTK_REAL_VEC PDstandardNth1K13;
    CCTK_REAL_VEC PDstandardNth2K13;
    CCTK_REAL_VEC PDstandardNth3K13;
    CCTK_REAL_VEC PDstandardNth1K22;
    CCTK_REAL_VEC PDstandardNth2K22;
    CCTK_REAL_VEC PDstandardNth3K22;
    CCTK_REAL_VEC PDstandardNth1K23;
    CCTK_REAL_VEC PDstandardNth2K23;
    CCTK_REAL_VEC PDstandardNth3K23;
    CCTK_REAL_VEC PDstandardNth1K33;
    CCTK_REAL_VEC PDstandardNth2K33;
    CCTK_REAL_VEC PDstandardNth3K33;
    
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
        PDstandardNth1g11 = PDstandardNthfdOrder21(&g11[index]);
        PDstandardNth2g11 = PDstandardNthfdOrder22(&g11[index]);
        PDstandardNth3g11 = PDstandardNthfdOrder23(&g11[index]);
        PDstandardNth11g11 = PDstandardNthfdOrder211(&g11[index]);
        PDstandardNth22g11 = PDstandardNthfdOrder222(&g11[index]);
        PDstandardNth33g11 = PDstandardNthfdOrder233(&g11[index]);
        PDstandardNth12g11 = PDstandardNthfdOrder212(&g11[index]);
        PDstandardNth13g11 = PDstandardNthfdOrder213(&g11[index]);
        PDstandardNth23g11 = PDstandardNthfdOrder223(&g11[index]);
        PDstandardNth1g12 = PDstandardNthfdOrder21(&g12[index]);
        PDstandardNth2g12 = PDstandardNthfdOrder22(&g12[index]);
        PDstandardNth3g12 = PDstandardNthfdOrder23(&g12[index]);
        PDstandardNth11g12 = PDstandardNthfdOrder211(&g12[index]);
        PDstandardNth22g12 = PDstandardNthfdOrder222(&g12[index]);
        PDstandardNth33g12 = PDstandardNthfdOrder233(&g12[index]);
        PDstandardNth12g12 = PDstandardNthfdOrder212(&g12[index]);
        PDstandardNth13g12 = PDstandardNthfdOrder213(&g12[index]);
        PDstandardNth23g12 = PDstandardNthfdOrder223(&g12[index]);
        PDstandardNth1g13 = PDstandardNthfdOrder21(&g13[index]);
        PDstandardNth2g13 = PDstandardNthfdOrder22(&g13[index]);
        PDstandardNth3g13 = PDstandardNthfdOrder23(&g13[index]);
        PDstandardNth11g13 = PDstandardNthfdOrder211(&g13[index]);
        PDstandardNth22g13 = PDstandardNthfdOrder222(&g13[index]);
        PDstandardNth33g13 = PDstandardNthfdOrder233(&g13[index]);
        PDstandardNth12g13 = PDstandardNthfdOrder212(&g13[index]);
        PDstandardNth13g13 = PDstandardNthfdOrder213(&g13[index]);
        PDstandardNth23g13 = PDstandardNthfdOrder223(&g13[index]);
        PDstandardNth1g22 = PDstandardNthfdOrder21(&g22[index]);
        PDstandardNth2g22 = PDstandardNthfdOrder22(&g22[index]);
        PDstandardNth3g22 = PDstandardNthfdOrder23(&g22[index]);
        PDstandardNth11g22 = PDstandardNthfdOrder211(&g22[index]);
        PDstandardNth22g22 = PDstandardNthfdOrder222(&g22[index]);
        PDstandardNth33g22 = PDstandardNthfdOrder233(&g22[index]);
        PDstandardNth12g22 = PDstandardNthfdOrder212(&g22[index]);
        PDstandardNth13g22 = PDstandardNthfdOrder213(&g22[index]);
        PDstandardNth23g22 = PDstandardNthfdOrder223(&g22[index]);
        PDstandardNth1g23 = PDstandardNthfdOrder21(&g23[index]);
        PDstandardNth2g23 = PDstandardNthfdOrder22(&g23[index]);
        PDstandardNth3g23 = PDstandardNthfdOrder23(&g23[index]);
        PDstandardNth11g23 = PDstandardNthfdOrder211(&g23[index]);
        PDstandardNth22g23 = PDstandardNthfdOrder222(&g23[index]);
        PDstandardNth33g23 = PDstandardNthfdOrder233(&g23[index]);
        PDstandardNth12g23 = PDstandardNthfdOrder212(&g23[index]);
        PDstandardNth13g23 = PDstandardNthfdOrder213(&g23[index]);
        PDstandardNth23g23 = PDstandardNthfdOrder223(&g23[index]);
        PDstandardNth1g33 = PDstandardNthfdOrder21(&g33[index]);
        PDstandardNth2g33 = PDstandardNthfdOrder22(&g33[index]);
        PDstandardNth3g33 = PDstandardNthfdOrder23(&g33[index]);
        PDstandardNth11g33 = PDstandardNthfdOrder211(&g33[index]);
        PDstandardNth22g33 = PDstandardNthfdOrder222(&g33[index]);
        PDstandardNth33g33 = PDstandardNthfdOrder233(&g33[index]);
        PDstandardNth12g33 = PDstandardNthfdOrder212(&g33[index]);
        PDstandardNth13g33 = PDstandardNthfdOrder213(&g33[index]);
        PDstandardNth23g33 = PDstandardNthfdOrder223(&g33[index]);
        PDstandardNth1K11 = PDstandardNthfdOrder21(&K11[index]);
        PDstandardNth2K11 = PDstandardNthfdOrder22(&K11[index]);
        PDstandardNth3K11 = PDstandardNthfdOrder23(&K11[index]);
        PDstandardNth1K12 = PDstandardNthfdOrder21(&K12[index]);
        PDstandardNth2K12 = PDstandardNthfdOrder22(&K12[index]);
        PDstandardNth3K12 = PDstandardNthfdOrder23(&K12[index]);
        PDstandardNth1K13 = PDstandardNthfdOrder21(&K13[index]);
        PDstandardNth2K13 = PDstandardNthfdOrder22(&K13[index]);
        PDstandardNth3K13 = PDstandardNthfdOrder23(&K13[index]);
        PDstandardNth1K22 = PDstandardNthfdOrder21(&K22[index]);
        PDstandardNth2K22 = PDstandardNthfdOrder22(&K22[index]);
        PDstandardNth3K22 = PDstandardNthfdOrder23(&K22[index]);
        PDstandardNth1K23 = PDstandardNthfdOrder21(&K23[index]);
        PDstandardNth2K23 = PDstandardNthfdOrder22(&K23[index]);
        PDstandardNth3K23 = PDstandardNthfdOrder23(&K23[index]);
        PDstandardNth1K33 = PDstandardNthfdOrder21(&K33[index]);
        PDstandardNth2K33 = PDstandardNthfdOrder22(&K33[index]);
        PDstandardNth3K33 = PDstandardNthfdOrder23(&K33[index]);
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
        PDstandardNth1g11 = PDstandardNthfdOrder41(&g11[index]);
        PDstandardNth2g11 = PDstandardNthfdOrder42(&g11[index]);
        PDstandardNth3g11 = PDstandardNthfdOrder43(&g11[index]);
        PDstandardNth11g11 = PDstandardNthfdOrder411(&g11[index]);
        PDstandardNth22g11 = PDstandardNthfdOrder422(&g11[index]);
        PDstandardNth33g11 = PDstandardNthfdOrder433(&g11[index]);
        PDstandardNth12g11 = PDstandardNthfdOrder412(&g11[index]);
        PDstandardNth13g11 = PDstandardNthfdOrder413(&g11[index]);
        PDstandardNth23g11 = PDstandardNthfdOrder423(&g11[index]);
        PDstandardNth1g12 = PDstandardNthfdOrder41(&g12[index]);
        PDstandardNth2g12 = PDstandardNthfdOrder42(&g12[index]);
        PDstandardNth3g12 = PDstandardNthfdOrder43(&g12[index]);
        PDstandardNth11g12 = PDstandardNthfdOrder411(&g12[index]);
        PDstandardNth22g12 = PDstandardNthfdOrder422(&g12[index]);
        PDstandardNth33g12 = PDstandardNthfdOrder433(&g12[index]);
        PDstandardNth12g12 = PDstandardNthfdOrder412(&g12[index]);
        PDstandardNth13g12 = PDstandardNthfdOrder413(&g12[index]);
        PDstandardNth23g12 = PDstandardNthfdOrder423(&g12[index]);
        PDstandardNth1g13 = PDstandardNthfdOrder41(&g13[index]);
        PDstandardNth2g13 = PDstandardNthfdOrder42(&g13[index]);
        PDstandardNth3g13 = PDstandardNthfdOrder43(&g13[index]);
        PDstandardNth11g13 = PDstandardNthfdOrder411(&g13[index]);
        PDstandardNth22g13 = PDstandardNthfdOrder422(&g13[index]);
        PDstandardNth33g13 = PDstandardNthfdOrder433(&g13[index]);
        PDstandardNth12g13 = PDstandardNthfdOrder412(&g13[index]);
        PDstandardNth13g13 = PDstandardNthfdOrder413(&g13[index]);
        PDstandardNth23g13 = PDstandardNthfdOrder423(&g13[index]);
        PDstandardNth1g22 = PDstandardNthfdOrder41(&g22[index]);
        PDstandardNth2g22 = PDstandardNthfdOrder42(&g22[index]);
        PDstandardNth3g22 = PDstandardNthfdOrder43(&g22[index]);
        PDstandardNth11g22 = PDstandardNthfdOrder411(&g22[index]);
        PDstandardNth22g22 = PDstandardNthfdOrder422(&g22[index]);
        PDstandardNth33g22 = PDstandardNthfdOrder433(&g22[index]);
        PDstandardNth12g22 = PDstandardNthfdOrder412(&g22[index]);
        PDstandardNth13g22 = PDstandardNthfdOrder413(&g22[index]);
        PDstandardNth23g22 = PDstandardNthfdOrder423(&g22[index]);
        PDstandardNth1g23 = PDstandardNthfdOrder41(&g23[index]);
        PDstandardNth2g23 = PDstandardNthfdOrder42(&g23[index]);
        PDstandardNth3g23 = PDstandardNthfdOrder43(&g23[index]);
        PDstandardNth11g23 = PDstandardNthfdOrder411(&g23[index]);
        PDstandardNth22g23 = PDstandardNthfdOrder422(&g23[index]);
        PDstandardNth33g23 = PDstandardNthfdOrder433(&g23[index]);
        PDstandardNth12g23 = PDstandardNthfdOrder412(&g23[index]);
        PDstandardNth13g23 = PDstandardNthfdOrder413(&g23[index]);
        PDstandardNth23g23 = PDstandardNthfdOrder423(&g23[index]);
        PDstandardNth1g33 = PDstandardNthfdOrder41(&g33[index]);
        PDstandardNth2g33 = PDstandardNthfdOrder42(&g33[index]);
        PDstandardNth3g33 = PDstandardNthfdOrder43(&g33[index]);
        PDstandardNth11g33 = PDstandardNthfdOrder411(&g33[index]);
        PDstandardNth22g33 = PDstandardNthfdOrder422(&g33[index]);
        PDstandardNth33g33 = PDstandardNthfdOrder433(&g33[index]);
        PDstandardNth12g33 = PDstandardNthfdOrder412(&g33[index]);
        PDstandardNth13g33 = PDstandardNthfdOrder413(&g33[index]);
        PDstandardNth23g33 = PDstandardNthfdOrder423(&g33[index]);
        PDstandardNth1K11 = PDstandardNthfdOrder41(&K11[index]);
        PDstandardNth2K11 = PDstandardNthfdOrder42(&K11[index]);
        PDstandardNth3K11 = PDstandardNthfdOrder43(&K11[index]);
        PDstandardNth1K12 = PDstandardNthfdOrder41(&K12[index]);
        PDstandardNth2K12 = PDstandardNthfdOrder42(&K12[index]);
        PDstandardNth3K12 = PDstandardNthfdOrder43(&K12[index]);
        PDstandardNth1K13 = PDstandardNthfdOrder41(&K13[index]);
        PDstandardNth2K13 = PDstandardNthfdOrder42(&K13[index]);
        PDstandardNth3K13 = PDstandardNthfdOrder43(&K13[index]);
        PDstandardNth1K22 = PDstandardNthfdOrder41(&K22[index]);
        PDstandardNth2K22 = PDstandardNthfdOrder42(&K22[index]);
        PDstandardNth3K22 = PDstandardNthfdOrder43(&K22[index]);
        PDstandardNth1K23 = PDstandardNthfdOrder41(&K23[index]);
        PDstandardNth2K23 = PDstandardNthfdOrder42(&K23[index]);
        PDstandardNth3K23 = PDstandardNthfdOrder43(&K23[index]);
        PDstandardNth1K33 = PDstandardNthfdOrder41(&K33[index]);
        PDstandardNth2K33 = PDstandardNthfdOrder42(&K33[index]);
        PDstandardNth3K33 = PDstandardNthfdOrder43(&K33[index]);
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
        PDstandardNth1g11 = PDstandardNthfdOrder61(&g11[index]);
        PDstandardNth2g11 = PDstandardNthfdOrder62(&g11[index]);
        PDstandardNth3g11 = PDstandardNthfdOrder63(&g11[index]);
        PDstandardNth11g11 = PDstandardNthfdOrder611(&g11[index]);
        PDstandardNth22g11 = PDstandardNthfdOrder622(&g11[index]);
        PDstandardNth33g11 = PDstandardNthfdOrder633(&g11[index]);
        PDstandardNth12g11 = PDstandardNthfdOrder612(&g11[index]);
        PDstandardNth13g11 = PDstandardNthfdOrder613(&g11[index]);
        PDstandardNth23g11 = PDstandardNthfdOrder623(&g11[index]);
        PDstandardNth1g12 = PDstandardNthfdOrder61(&g12[index]);
        PDstandardNth2g12 = PDstandardNthfdOrder62(&g12[index]);
        PDstandardNth3g12 = PDstandardNthfdOrder63(&g12[index]);
        PDstandardNth11g12 = PDstandardNthfdOrder611(&g12[index]);
        PDstandardNth22g12 = PDstandardNthfdOrder622(&g12[index]);
        PDstandardNth33g12 = PDstandardNthfdOrder633(&g12[index]);
        PDstandardNth12g12 = PDstandardNthfdOrder612(&g12[index]);
        PDstandardNth13g12 = PDstandardNthfdOrder613(&g12[index]);
        PDstandardNth23g12 = PDstandardNthfdOrder623(&g12[index]);
        PDstandardNth1g13 = PDstandardNthfdOrder61(&g13[index]);
        PDstandardNth2g13 = PDstandardNthfdOrder62(&g13[index]);
        PDstandardNth3g13 = PDstandardNthfdOrder63(&g13[index]);
        PDstandardNth11g13 = PDstandardNthfdOrder611(&g13[index]);
        PDstandardNth22g13 = PDstandardNthfdOrder622(&g13[index]);
        PDstandardNth33g13 = PDstandardNthfdOrder633(&g13[index]);
        PDstandardNth12g13 = PDstandardNthfdOrder612(&g13[index]);
        PDstandardNth13g13 = PDstandardNthfdOrder613(&g13[index]);
        PDstandardNth23g13 = PDstandardNthfdOrder623(&g13[index]);
        PDstandardNth1g22 = PDstandardNthfdOrder61(&g22[index]);
        PDstandardNth2g22 = PDstandardNthfdOrder62(&g22[index]);
        PDstandardNth3g22 = PDstandardNthfdOrder63(&g22[index]);
        PDstandardNth11g22 = PDstandardNthfdOrder611(&g22[index]);
        PDstandardNth22g22 = PDstandardNthfdOrder622(&g22[index]);
        PDstandardNth33g22 = PDstandardNthfdOrder633(&g22[index]);
        PDstandardNth12g22 = PDstandardNthfdOrder612(&g22[index]);
        PDstandardNth13g22 = PDstandardNthfdOrder613(&g22[index]);
        PDstandardNth23g22 = PDstandardNthfdOrder623(&g22[index]);
        PDstandardNth1g23 = PDstandardNthfdOrder61(&g23[index]);
        PDstandardNth2g23 = PDstandardNthfdOrder62(&g23[index]);
        PDstandardNth3g23 = PDstandardNthfdOrder63(&g23[index]);
        PDstandardNth11g23 = PDstandardNthfdOrder611(&g23[index]);
        PDstandardNth22g23 = PDstandardNthfdOrder622(&g23[index]);
        PDstandardNth33g23 = PDstandardNthfdOrder633(&g23[index]);
        PDstandardNth12g23 = PDstandardNthfdOrder612(&g23[index]);
        PDstandardNth13g23 = PDstandardNthfdOrder613(&g23[index]);
        PDstandardNth23g23 = PDstandardNthfdOrder623(&g23[index]);
        PDstandardNth1g33 = PDstandardNthfdOrder61(&g33[index]);
        PDstandardNth2g33 = PDstandardNthfdOrder62(&g33[index]);
        PDstandardNth3g33 = PDstandardNthfdOrder63(&g33[index]);
        PDstandardNth11g33 = PDstandardNthfdOrder611(&g33[index]);
        PDstandardNth22g33 = PDstandardNthfdOrder622(&g33[index]);
        PDstandardNth33g33 = PDstandardNthfdOrder633(&g33[index]);
        PDstandardNth12g33 = PDstandardNthfdOrder612(&g33[index]);
        PDstandardNth13g33 = PDstandardNthfdOrder613(&g33[index]);
        PDstandardNth23g33 = PDstandardNthfdOrder623(&g33[index]);
        PDstandardNth1K11 = PDstandardNthfdOrder61(&K11[index]);
        PDstandardNth2K11 = PDstandardNthfdOrder62(&K11[index]);
        PDstandardNth3K11 = PDstandardNthfdOrder63(&K11[index]);
        PDstandardNth1K12 = PDstandardNthfdOrder61(&K12[index]);
        PDstandardNth2K12 = PDstandardNthfdOrder62(&K12[index]);
        PDstandardNth3K12 = PDstandardNthfdOrder63(&K12[index]);
        PDstandardNth1K13 = PDstandardNthfdOrder61(&K13[index]);
        PDstandardNth2K13 = PDstandardNthfdOrder62(&K13[index]);
        PDstandardNth3K13 = PDstandardNthfdOrder63(&K13[index]);
        PDstandardNth1K22 = PDstandardNthfdOrder61(&K22[index]);
        PDstandardNth2K22 = PDstandardNthfdOrder62(&K22[index]);
        PDstandardNth3K22 = PDstandardNthfdOrder63(&K22[index]);
        PDstandardNth1K23 = PDstandardNthfdOrder61(&K23[index]);
        PDstandardNth2K23 = PDstandardNthfdOrder62(&K23[index]);
        PDstandardNth3K23 = PDstandardNthfdOrder63(&K23[index]);
        PDstandardNth1K33 = PDstandardNthfdOrder61(&K33[index]);
        PDstandardNth2K33 = PDstandardNthfdOrder62(&K33[index]);
        PDstandardNth3K33 = PDstandardNthfdOrder63(&K33[index]);
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
        PDstandardNth1g11 = PDstandardNthfdOrder81(&g11[index]);
        PDstandardNth2g11 = PDstandardNthfdOrder82(&g11[index]);
        PDstandardNth3g11 = PDstandardNthfdOrder83(&g11[index]);
        PDstandardNth11g11 = PDstandardNthfdOrder811(&g11[index]);
        PDstandardNth22g11 = PDstandardNthfdOrder822(&g11[index]);
        PDstandardNth33g11 = PDstandardNthfdOrder833(&g11[index]);
        PDstandardNth12g11 = PDstandardNthfdOrder812(&g11[index]);
        PDstandardNth13g11 = PDstandardNthfdOrder813(&g11[index]);
        PDstandardNth23g11 = PDstandardNthfdOrder823(&g11[index]);
        PDstandardNth1g12 = PDstandardNthfdOrder81(&g12[index]);
        PDstandardNth2g12 = PDstandardNthfdOrder82(&g12[index]);
        PDstandardNth3g12 = PDstandardNthfdOrder83(&g12[index]);
        PDstandardNth11g12 = PDstandardNthfdOrder811(&g12[index]);
        PDstandardNth22g12 = PDstandardNthfdOrder822(&g12[index]);
        PDstandardNth33g12 = PDstandardNthfdOrder833(&g12[index]);
        PDstandardNth12g12 = PDstandardNthfdOrder812(&g12[index]);
        PDstandardNth13g12 = PDstandardNthfdOrder813(&g12[index]);
        PDstandardNth23g12 = PDstandardNthfdOrder823(&g12[index]);
        PDstandardNth1g13 = PDstandardNthfdOrder81(&g13[index]);
        PDstandardNth2g13 = PDstandardNthfdOrder82(&g13[index]);
        PDstandardNth3g13 = PDstandardNthfdOrder83(&g13[index]);
        PDstandardNth11g13 = PDstandardNthfdOrder811(&g13[index]);
        PDstandardNth22g13 = PDstandardNthfdOrder822(&g13[index]);
        PDstandardNth33g13 = PDstandardNthfdOrder833(&g13[index]);
        PDstandardNth12g13 = PDstandardNthfdOrder812(&g13[index]);
        PDstandardNth13g13 = PDstandardNthfdOrder813(&g13[index]);
        PDstandardNth23g13 = PDstandardNthfdOrder823(&g13[index]);
        PDstandardNth1g22 = PDstandardNthfdOrder81(&g22[index]);
        PDstandardNth2g22 = PDstandardNthfdOrder82(&g22[index]);
        PDstandardNth3g22 = PDstandardNthfdOrder83(&g22[index]);
        PDstandardNth11g22 = PDstandardNthfdOrder811(&g22[index]);
        PDstandardNth22g22 = PDstandardNthfdOrder822(&g22[index]);
        PDstandardNth33g22 = PDstandardNthfdOrder833(&g22[index]);
        PDstandardNth12g22 = PDstandardNthfdOrder812(&g22[index]);
        PDstandardNth13g22 = PDstandardNthfdOrder813(&g22[index]);
        PDstandardNth23g22 = PDstandardNthfdOrder823(&g22[index]);
        PDstandardNth1g23 = PDstandardNthfdOrder81(&g23[index]);
        PDstandardNth2g23 = PDstandardNthfdOrder82(&g23[index]);
        PDstandardNth3g23 = PDstandardNthfdOrder83(&g23[index]);
        PDstandardNth11g23 = PDstandardNthfdOrder811(&g23[index]);
        PDstandardNth22g23 = PDstandardNthfdOrder822(&g23[index]);
        PDstandardNth33g23 = PDstandardNthfdOrder833(&g23[index]);
        PDstandardNth12g23 = PDstandardNthfdOrder812(&g23[index]);
        PDstandardNth13g23 = PDstandardNthfdOrder813(&g23[index]);
        PDstandardNth23g23 = PDstandardNthfdOrder823(&g23[index]);
        PDstandardNth1g33 = PDstandardNthfdOrder81(&g33[index]);
        PDstandardNth2g33 = PDstandardNthfdOrder82(&g33[index]);
        PDstandardNth3g33 = PDstandardNthfdOrder83(&g33[index]);
        PDstandardNth11g33 = PDstandardNthfdOrder811(&g33[index]);
        PDstandardNth22g33 = PDstandardNthfdOrder822(&g33[index]);
        PDstandardNth33g33 = PDstandardNthfdOrder833(&g33[index]);
        PDstandardNth12g33 = PDstandardNthfdOrder812(&g33[index]);
        PDstandardNth13g33 = PDstandardNthfdOrder813(&g33[index]);
        PDstandardNth23g33 = PDstandardNthfdOrder823(&g33[index]);
        PDstandardNth1K11 = PDstandardNthfdOrder81(&K11[index]);
        PDstandardNth2K11 = PDstandardNthfdOrder82(&K11[index]);
        PDstandardNth3K11 = PDstandardNthfdOrder83(&K11[index]);
        PDstandardNth1K12 = PDstandardNthfdOrder81(&K12[index]);
        PDstandardNth2K12 = PDstandardNthfdOrder82(&K12[index]);
        PDstandardNth3K12 = PDstandardNthfdOrder83(&K12[index]);
        PDstandardNth1K13 = PDstandardNthfdOrder81(&K13[index]);
        PDstandardNth2K13 = PDstandardNthfdOrder82(&K13[index]);
        PDstandardNth3K13 = PDstandardNthfdOrder83(&K13[index]);
        PDstandardNth1K22 = PDstandardNthfdOrder81(&K22[index]);
        PDstandardNth2K22 = PDstandardNthfdOrder82(&K22[index]);
        PDstandardNth3K22 = PDstandardNthfdOrder83(&K22[index]);
        PDstandardNth1K23 = PDstandardNthfdOrder81(&K23[index]);
        PDstandardNth2K23 = PDstandardNthfdOrder82(&K23[index]);
        PDstandardNth3K23 = PDstandardNthfdOrder83(&K23[index]);
        PDstandardNth1K33 = PDstandardNthfdOrder81(&K33[index]);
        PDstandardNth2K33 = PDstandardNthfdOrder82(&K33[index]);
        PDstandardNth3K33 = PDstandardNthfdOrder83(&K33[index]);
        break;
    }
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC JacPDstandardNth11alpha;
    CCTK_REAL_VEC JacPDstandardNth11g22;
    CCTK_REAL_VEC JacPDstandardNth11g23;
    CCTK_REAL_VEC JacPDstandardNth11g33;
    CCTK_REAL_VEC JacPDstandardNth12alpha;
    CCTK_REAL_VEC JacPDstandardNth12g11;
    CCTK_REAL_VEC JacPDstandardNth12g12;
    CCTK_REAL_VEC JacPDstandardNth12g13;
    CCTK_REAL_VEC JacPDstandardNth12g22;
    CCTK_REAL_VEC JacPDstandardNth12g23;
    CCTK_REAL_VEC JacPDstandardNth12g33;
    CCTK_REAL_VEC JacPDstandardNth13alpha;
    CCTK_REAL_VEC JacPDstandardNth13g11;
    CCTK_REAL_VEC JacPDstandardNth13g12;
    CCTK_REAL_VEC JacPDstandardNth13g13;
    CCTK_REAL_VEC JacPDstandardNth13g22;
    CCTK_REAL_VEC JacPDstandardNth13g23;
    CCTK_REAL_VEC JacPDstandardNth13g33;
    CCTK_REAL_VEC JacPDstandardNth1alpha;
    CCTK_REAL_VEC JacPDstandardNth1beta1;
    CCTK_REAL_VEC JacPDstandardNth1beta2;
    CCTK_REAL_VEC JacPDstandardNth1beta3;
    CCTK_REAL_VEC JacPDstandardNth1g11;
    CCTK_REAL_VEC JacPDstandardNth1g12;
    CCTK_REAL_VEC JacPDstandardNth1g13;
    CCTK_REAL_VEC JacPDstandardNth1g22;
    CCTK_REAL_VEC JacPDstandardNth1g23;
    CCTK_REAL_VEC JacPDstandardNth1g33;
    CCTK_REAL_VEC JacPDstandardNth1K11;
    CCTK_REAL_VEC JacPDstandardNth1K12;
    CCTK_REAL_VEC JacPDstandardNth1K13;
    CCTK_REAL_VEC JacPDstandardNth1K22;
    CCTK_REAL_VEC JacPDstandardNth1K23;
    CCTK_REAL_VEC JacPDstandardNth1K33;
    CCTK_REAL_VEC JacPDstandardNth21g11;
    CCTK_REAL_VEC JacPDstandardNth21g12;
    CCTK_REAL_VEC JacPDstandardNth21g13;
    CCTK_REAL_VEC JacPDstandardNth21g22;
    CCTK_REAL_VEC JacPDstandardNth21g23;
    CCTK_REAL_VEC JacPDstandardNth21g33;
    CCTK_REAL_VEC JacPDstandardNth22alpha;
    CCTK_REAL_VEC JacPDstandardNth22g11;
    CCTK_REAL_VEC JacPDstandardNth22g13;
    CCTK_REAL_VEC JacPDstandardNth22g33;
    CCTK_REAL_VEC JacPDstandardNth23alpha;
    CCTK_REAL_VEC JacPDstandardNth23g11;
    CCTK_REAL_VEC JacPDstandardNth23g12;
    CCTK_REAL_VEC JacPDstandardNth23g13;
    CCTK_REAL_VEC JacPDstandardNth23g22;
    CCTK_REAL_VEC JacPDstandardNth23g23;
    CCTK_REAL_VEC JacPDstandardNth23g33;
    CCTK_REAL_VEC JacPDstandardNth2alpha;
    CCTK_REAL_VEC JacPDstandardNth2beta1;
    CCTK_REAL_VEC JacPDstandardNth2beta2;
    CCTK_REAL_VEC JacPDstandardNth2beta3;
    CCTK_REAL_VEC JacPDstandardNth2g11;
    CCTK_REAL_VEC JacPDstandardNth2g12;
    CCTK_REAL_VEC JacPDstandardNth2g13;
    CCTK_REAL_VEC JacPDstandardNth2g22;
    CCTK_REAL_VEC JacPDstandardNth2g23;
    CCTK_REAL_VEC JacPDstandardNth2g33;
    CCTK_REAL_VEC JacPDstandardNth2K11;
    CCTK_REAL_VEC JacPDstandardNth2K12;
    CCTK_REAL_VEC JacPDstandardNth2K13;
    CCTK_REAL_VEC JacPDstandardNth2K22;
    CCTK_REAL_VEC JacPDstandardNth2K23;
    CCTK_REAL_VEC JacPDstandardNth2K33;
    CCTK_REAL_VEC JacPDstandardNth31g11;
    CCTK_REAL_VEC JacPDstandardNth31g12;
    CCTK_REAL_VEC JacPDstandardNth31g13;
    CCTK_REAL_VEC JacPDstandardNth31g22;
    CCTK_REAL_VEC JacPDstandardNth31g23;
    CCTK_REAL_VEC JacPDstandardNth31g33;
    CCTK_REAL_VEC JacPDstandardNth32g11;
    CCTK_REAL_VEC JacPDstandardNth32g12;
    CCTK_REAL_VEC JacPDstandardNth32g13;
    CCTK_REAL_VEC JacPDstandardNth32g22;
    CCTK_REAL_VEC JacPDstandardNth32g23;
    CCTK_REAL_VEC JacPDstandardNth32g33;
    CCTK_REAL_VEC JacPDstandardNth33alpha;
    CCTK_REAL_VEC JacPDstandardNth33g11;
    CCTK_REAL_VEC JacPDstandardNth33g12;
    CCTK_REAL_VEC JacPDstandardNth33g22;
    CCTK_REAL_VEC JacPDstandardNth3alpha;
    CCTK_REAL_VEC JacPDstandardNth3beta1;
    CCTK_REAL_VEC JacPDstandardNth3beta2;
    CCTK_REAL_VEC JacPDstandardNth3beta3;
    CCTK_REAL_VEC JacPDstandardNth3g11;
    CCTK_REAL_VEC JacPDstandardNth3g12;
    CCTK_REAL_VEC JacPDstandardNth3g13;
    CCTK_REAL_VEC JacPDstandardNth3g22;
    CCTK_REAL_VEC JacPDstandardNth3g23;
    CCTK_REAL_VEC JacPDstandardNth3g33;
    CCTK_REAL_VEC JacPDstandardNth3K11;
    CCTK_REAL_VEC JacPDstandardNth3K12;
    CCTK_REAL_VEC JacPDstandardNth3K13;
    CCTK_REAL_VEC JacPDstandardNth3K22;
    CCTK_REAL_VEC JacPDstandardNth3K23;
    CCTK_REAL_VEC JacPDstandardNth3K33;
    
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
      
      JacPDstandardNth1g11 = 
        kmadd(J11L,PDstandardNth1g11,kmadd(J21L,PDstandardNth2g11,kmul(J31L,PDstandardNth3g11)));
      
      JacPDstandardNth1g12 = 
        kmadd(J11L,PDstandardNth1g12,kmadd(J21L,PDstandardNth2g12,kmul(J31L,PDstandardNth3g12)));
      
      JacPDstandardNth1g13 = 
        kmadd(J11L,PDstandardNth1g13,kmadd(J21L,PDstandardNth2g13,kmul(J31L,PDstandardNth3g13)));
      
      JacPDstandardNth1g22 = 
        kmadd(J11L,PDstandardNth1g22,kmadd(J21L,PDstandardNth2g22,kmul(J31L,PDstandardNth3g22)));
      
      JacPDstandardNth1g23 = 
        kmadd(J11L,PDstandardNth1g23,kmadd(J21L,PDstandardNth2g23,kmul(J31L,PDstandardNth3g23)));
      
      JacPDstandardNth1g33 = 
        kmadd(J11L,PDstandardNth1g33,kmadd(J21L,PDstandardNth2g33,kmul(J31L,PDstandardNth3g33)));
      
      JacPDstandardNth1K11 = 
        kmadd(J11L,PDstandardNth1K11,kmadd(J21L,PDstandardNth2K11,kmul(J31L,PDstandardNth3K11)));
      
      JacPDstandardNth1K12 = 
        kmadd(J11L,PDstandardNth1K12,kmadd(J21L,PDstandardNth2K12,kmul(J31L,PDstandardNth3K12)));
      
      JacPDstandardNth1K13 = 
        kmadd(J11L,PDstandardNth1K13,kmadd(J21L,PDstandardNth2K13,kmul(J31L,PDstandardNth3K13)));
      
      JacPDstandardNth1K22 = 
        kmadd(J11L,PDstandardNth1K22,kmadd(J21L,PDstandardNth2K22,kmul(J31L,PDstandardNth3K22)));
      
      JacPDstandardNth1K23 = 
        kmadd(J11L,PDstandardNth1K23,kmadd(J21L,PDstandardNth2K23,kmul(J31L,PDstandardNth3K23)));
      
      JacPDstandardNth1K33 = 
        kmadd(J11L,PDstandardNth1K33,kmadd(J21L,PDstandardNth2K33,kmul(J31L,PDstandardNth3K33)));
      
      JacPDstandardNth2alpha = 
        kmadd(J12L,PDstandardNth1alpha,kmadd(J22L,PDstandardNth2alpha,kmul(J32L,PDstandardNth3alpha)));
      
      JacPDstandardNth2beta1 = 
        kmadd(J12L,PDstandardNth1beta1,kmadd(J22L,PDstandardNth2beta1,kmul(J32L,PDstandardNth3beta1)));
      
      JacPDstandardNth2beta2 = 
        kmadd(J12L,PDstandardNth1beta2,kmadd(J22L,PDstandardNth2beta2,kmul(J32L,PDstandardNth3beta2)));
      
      JacPDstandardNth2beta3 = 
        kmadd(J12L,PDstandardNth1beta3,kmadd(J22L,PDstandardNth2beta3,kmul(J32L,PDstandardNth3beta3)));
      
      JacPDstandardNth2g11 = 
        kmadd(J12L,PDstandardNth1g11,kmadd(J22L,PDstandardNth2g11,kmul(J32L,PDstandardNth3g11)));
      
      JacPDstandardNth2g12 = 
        kmadd(J12L,PDstandardNth1g12,kmadd(J22L,PDstandardNth2g12,kmul(J32L,PDstandardNth3g12)));
      
      JacPDstandardNth2g13 = 
        kmadd(J12L,PDstandardNth1g13,kmadd(J22L,PDstandardNth2g13,kmul(J32L,PDstandardNth3g13)));
      
      JacPDstandardNth2g22 = 
        kmadd(J12L,PDstandardNth1g22,kmadd(J22L,PDstandardNth2g22,kmul(J32L,PDstandardNth3g22)));
      
      JacPDstandardNth2g23 = 
        kmadd(J12L,PDstandardNth1g23,kmadd(J22L,PDstandardNth2g23,kmul(J32L,PDstandardNth3g23)));
      
      JacPDstandardNth2g33 = 
        kmadd(J12L,PDstandardNth1g33,kmadd(J22L,PDstandardNth2g33,kmul(J32L,PDstandardNth3g33)));
      
      JacPDstandardNth2K11 = 
        kmadd(J12L,PDstandardNth1K11,kmadd(J22L,PDstandardNth2K11,kmul(J32L,PDstandardNth3K11)));
      
      JacPDstandardNth2K12 = 
        kmadd(J12L,PDstandardNth1K12,kmadd(J22L,PDstandardNth2K12,kmul(J32L,PDstandardNth3K12)));
      
      JacPDstandardNth2K13 = 
        kmadd(J12L,PDstandardNth1K13,kmadd(J22L,PDstandardNth2K13,kmul(J32L,PDstandardNth3K13)));
      
      JacPDstandardNth2K22 = 
        kmadd(J12L,PDstandardNth1K22,kmadd(J22L,PDstandardNth2K22,kmul(J32L,PDstandardNth3K22)));
      
      JacPDstandardNth2K23 = 
        kmadd(J12L,PDstandardNth1K23,kmadd(J22L,PDstandardNth2K23,kmul(J32L,PDstandardNth3K23)));
      
      JacPDstandardNth2K33 = 
        kmadd(J12L,PDstandardNth1K33,kmadd(J22L,PDstandardNth2K33,kmul(J32L,PDstandardNth3K33)));
      
      JacPDstandardNth3alpha = 
        kmadd(J13L,PDstandardNth1alpha,kmadd(J23L,PDstandardNth2alpha,kmul(J33L,PDstandardNth3alpha)));
      
      JacPDstandardNth3beta1 = 
        kmadd(J13L,PDstandardNth1beta1,kmadd(J23L,PDstandardNth2beta1,kmul(J33L,PDstandardNth3beta1)));
      
      JacPDstandardNth3beta2 = 
        kmadd(J13L,PDstandardNth1beta2,kmadd(J23L,PDstandardNth2beta2,kmul(J33L,PDstandardNth3beta2)));
      
      JacPDstandardNth3beta3 = 
        kmadd(J13L,PDstandardNth1beta3,kmadd(J23L,PDstandardNth2beta3,kmul(J33L,PDstandardNth3beta3)));
      
      JacPDstandardNth3g11 = 
        kmadd(J13L,PDstandardNth1g11,kmadd(J23L,PDstandardNth2g11,kmul(J33L,PDstandardNth3g11)));
      
      JacPDstandardNth3g12 = 
        kmadd(J13L,PDstandardNth1g12,kmadd(J23L,PDstandardNth2g12,kmul(J33L,PDstandardNth3g12)));
      
      JacPDstandardNth3g13 = 
        kmadd(J13L,PDstandardNth1g13,kmadd(J23L,PDstandardNth2g13,kmul(J33L,PDstandardNth3g13)));
      
      JacPDstandardNth3g22 = 
        kmadd(J13L,PDstandardNth1g22,kmadd(J23L,PDstandardNth2g22,kmul(J33L,PDstandardNth3g22)));
      
      JacPDstandardNth3g23 = 
        kmadd(J13L,PDstandardNth1g23,kmadd(J23L,PDstandardNth2g23,kmul(J33L,PDstandardNth3g23)));
      
      JacPDstandardNth3g33 = 
        kmadd(J13L,PDstandardNth1g33,kmadd(J23L,PDstandardNth2g33,kmul(J33L,PDstandardNth3g33)));
      
      JacPDstandardNth3K11 = 
        kmadd(J13L,PDstandardNth1K11,kmadd(J23L,PDstandardNth2K11,kmul(J33L,PDstandardNth3K11)));
      
      JacPDstandardNth3K12 = 
        kmadd(J13L,PDstandardNth1K12,kmadd(J23L,PDstandardNth2K12,kmul(J33L,PDstandardNth3K12)));
      
      JacPDstandardNth3K13 = 
        kmadd(J13L,PDstandardNth1K13,kmadd(J23L,PDstandardNth2K13,kmul(J33L,PDstandardNth3K13)));
      
      JacPDstandardNth3K22 = 
        kmadd(J13L,PDstandardNth1K22,kmadd(J23L,PDstandardNth2K22,kmul(J33L,PDstandardNth3K22)));
      
      JacPDstandardNth3K23 = 
        kmadd(J13L,PDstandardNth1K23,kmadd(J23L,PDstandardNth2K23,kmul(J33L,PDstandardNth3K23)));
      
      JacPDstandardNth3K33 = 
        kmadd(J13L,PDstandardNth1K33,kmadd(J23L,PDstandardNth2K33,kmul(J33L,PDstandardNth3K33)));
      
      JacPDstandardNth11alpha = 
        kmadd(dJ111L,PDstandardNth1alpha,kmadd(dJ211L,PDstandardNth2alpha,kmadd(dJ311L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,SQR(J11L),kmadd(PDstandardNth22alpha,SQR(J21L),kmadd(PDstandardNth33alpha,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha)),kmul(J21L,kmul(J31L,PDstandardNth23alpha))),ToReal(2))))))));
      
      JacPDstandardNth11g22 = 
        kmadd(dJ111L,PDstandardNth1g22,kmadd(dJ211L,PDstandardNth2g22,kmadd(dJ311L,PDstandardNth3g22,kmadd(PDstandardNth11g22,SQR(J11L),kmadd(PDstandardNth22g22,SQR(J21L),kmadd(PDstandardNth33g22,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12g22,kmul(J31L,PDstandardNth13g22)),kmul(J21L,kmul(J31L,PDstandardNth23g22))),ToReal(2))))))));
      
      JacPDstandardNth11g23 = 
        kmadd(dJ111L,PDstandardNth1g23,kmadd(dJ211L,PDstandardNth2g23,kmadd(dJ311L,PDstandardNth3g23,kmadd(PDstandardNth11g23,SQR(J11L),kmadd(PDstandardNth22g23,SQR(J21L),kmadd(PDstandardNth33g23,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12g23,kmul(J31L,PDstandardNth13g23)),kmul(J21L,kmul(J31L,PDstandardNth23g23))),ToReal(2))))))));
      
      JacPDstandardNth11g33 = 
        kmadd(dJ111L,PDstandardNth1g33,kmadd(dJ211L,PDstandardNth2g33,kmadd(dJ311L,PDstandardNth3g33,kmadd(PDstandardNth11g33,SQR(J11L),kmadd(PDstandardNth22g33,SQR(J21L),kmadd(PDstandardNth33g33,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12g33,kmul(J31L,PDstandardNth13g33)),kmul(J21L,kmul(J31L,PDstandardNth23g33))),ToReal(2))))))));
      
      JacPDstandardNth22alpha = 
        kmadd(dJ122L,PDstandardNth1alpha,kmadd(dJ222L,PDstandardNth2alpha,kmadd(dJ322L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,SQR(J12L),kmadd(PDstandardNth22alpha,SQR(J22L),kmadd(PDstandardNth33alpha,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha)),kmul(J22L,kmul(J32L,PDstandardNth23alpha))),ToReal(2))))))));
      
      JacPDstandardNth22g11 = 
        kmadd(dJ122L,PDstandardNth1g11,kmadd(dJ222L,PDstandardNth2g11,kmadd(dJ322L,PDstandardNth3g11,kmadd(PDstandardNth11g11,SQR(J12L),kmadd(PDstandardNth22g11,SQR(J22L),kmadd(PDstandardNth33g11,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12g11,kmul(J32L,PDstandardNth13g11)),kmul(J22L,kmul(J32L,PDstandardNth23g11))),ToReal(2))))))));
      
      JacPDstandardNth22g13 = 
        kmadd(dJ122L,PDstandardNth1g13,kmadd(dJ222L,PDstandardNth2g13,kmadd(dJ322L,PDstandardNth3g13,kmadd(PDstandardNth11g13,SQR(J12L),kmadd(PDstandardNth22g13,SQR(J22L),kmadd(PDstandardNth33g13,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12g13,kmul(J32L,PDstandardNth13g13)),kmul(J22L,kmul(J32L,PDstandardNth23g13))),ToReal(2))))))));
      
      JacPDstandardNth22g33 = 
        kmadd(dJ122L,PDstandardNth1g33,kmadd(dJ222L,PDstandardNth2g33,kmadd(dJ322L,PDstandardNth3g33,kmadd(PDstandardNth11g33,SQR(J12L),kmadd(PDstandardNth22g33,SQR(J22L),kmadd(PDstandardNth33g33,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12g33,kmul(J32L,PDstandardNth13g33)),kmul(J22L,kmul(J32L,PDstandardNth23g33))),ToReal(2))))))));
      
      JacPDstandardNth33alpha = 
        kmadd(dJ133L,PDstandardNth1alpha,kmadd(dJ233L,PDstandardNth2alpha,kmadd(dJ333L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,SQR(J13L),kmadd(PDstandardNth22alpha,SQR(J23L),kmadd(PDstandardNth33alpha,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmul(J23L,kmul(J33L,PDstandardNth23alpha))),ToReal(2))))))));
      
      JacPDstandardNth33g11 = 
        kmadd(dJ133L,PDstandardNth1g11,kmadd(dJ233L,PDstandardNth2g11,kmadd(dJ333L,PDstandardNth3g11,kmadd(PDstandardNth11g11,SQR(J13L),kmadd(PDstandardNth22g11,SQR(J23L),kmadd(PDstandardNth33g11,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12g11,kmul(J33L,PDstandardNth13g11)),kmul(J23L,kmul(J33L,PDstandardNth23g11))),ToReal(2))))))));
      
      JacPDstandardNth33g12 = 
        kmadd(dJ133L,PDstandardNth1g12,kmadd(dJ233L,PDstandardNth2g12,kmadd(dJ333L,PDstandardNth3g12,kmadd(PDstandardNth11g12,SQR(J13L),kmadd(PDstandardNth22g12,SQR(J23L),kmadd(PDstandardNth33g12,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12g12,kmul(J33L,PDstandardNth13g12)),kmul(J23L,kmul(J33L,PDstandardNth23g12))),ToReal(2))))))));
      
      JacPDstandardNth33g22 = 
        kmadd(dJ133L,PDstandardNth1g22,kmadd(dJ233L,PDstandardNth2g22,kmadd(dJ333L,PDstandardNth3g22,kmadd(PDstandardNth11g22,SQR(J13L),kmadd(PDstandardNth22g22,SQR(J23L),kmadd(PDstandardNth33g22,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12g22,kmul(J33L,PDstandardNth13g22)),kmul(J23L,kmul(J33L,PDstandardNth23g22))),ToReal(2))))))));
      
      JacPDstandardNth12alpha = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11alpha,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha))),kmadd(J11L,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha)),kmadd(dJ112L,PDstandardNth1alpha,kmadd(J22L,kmadd(J21L,PDstandardNth22alpha,kmul(J31L,PDstandardNth23alpha)),kmadd(dJ212L,PDstandardNth2alpha,kmadd(J32L,kmadd(J21L,PDstandardNth23alpha,kmul(J31L,PDstandardNth33alpha)),kmul(dJ312L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth12g11 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g11,kmadd(J21L,PDstandardNth12g11,kmul(J31L,PDstandardNth13g11))),kmadd(J11L,kmadd(J22L,PDstandardNth12g11,kmul(J32L,PDstandardNth13g11)),kmadd(dJ112L,PDstandardNth1g11,kmadd(J22L,kmadd(J21L,PDstandardNth22g11,kmul(J31L,PDstandardNth23g11)),kmadd(dJ212L,PDstandardNth2g11,kmadd(J32L,kmadd(J21L,PDstandardNth23g11,kmul(J31L,PDstandardNth33g11)),kmul(dJ312L,PDstandardNth3g11)))))));
      
      JacPDstandardNth12g12 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g12,kmadd(J21L,PDstandardNth12g12,kmul(J31L,PDstandardNth13g12))),kmadd(J11L,kmadd(J22L,PDstandardNth12g12,kmul(J32L,PDstandardNth13g12)),kmadd(dJ112L,PDstandardNth1g12,kmadd(J22L,kmadd(J21L,PDstandardNth22g12,kmul(J31L,PDstandardNth23g12)),kmadd(dJ212L,PDstandardNth2g12,kmadd(J32L,kmadd(J21L,PDstandardNth23g12,kmul(J31L,PDstandardNth33g12)),kmul(dJ312L,PDstandardNth3g12)))))));
      
      JacPDstandardNth12g13 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g13,kmadd(J21L,PDstandardNth12g13,kmul(J31L,PDstandardNth13g13))),kmadd(J11L,kmadd(J22L,PDstandardNth12g13,kmul(J32L,PDstandardNth13g13)),kmadd(dJ112L,PDstandardNth1g13,kmadd(J22L,kmadd(J21L,PDstandardNth22g13,kmul(J31L,PDstandardNth23g13)),kmadd(dJ212L,PDstandardNth2g13,kmadd(J32L,kmadd(J21L,PDstandardNth23g13,kmul(J31L,PDstandardNth33g13)),kmul(dJ312L,PDstandardNth3g13)))))));
      
      JacPDstandardNth12g22 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g22,kmadd(J21L,PDstandardNth12g22,kmul(J31L,PDstandardNth13g22))),kmadd(J11L,kmadd(J22L,PDstandardNth12g22,kmul(J32L,PDstandardNth13g22)),kmadd(dJ112L,PDstandardNth1g22,kmadd(J22L,kmadd(J21L,PDstandardNth22g22,kmul(J31L,PDstandardNth23g22)),kmadd(dJ212L,PDstandardNth2g22,kmadd(J32L,kmadd(J21L,PDstandardNth23g22,kmul(J31L,PDstandardNth33g22)),kmul(dJ312L,PDstandardNth3g22)))))));
      
      JacPDstandardNth12g23 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g23,kmadd(J21L,PDstandardNth12g23,kmul(J31L,PDstandardNth13g23))),kmadd(J11L,kmadd(J22L,PDstandardNth12g23,kmul(J32L,PDstandardNth13g23)),kmadd(dJ112L,PDstandardNth1g23,kmadd(J22L,kmadd(J21L,PDstandardNth22g23,kmul(J31L,PDstandardNth23g23)),kmadd(dJ212L,PDstandardNth2g23,kmadd(J32L,kmadd(J21L,PDstandardNth23g23,kmul(J31L,PDstandardNth33g23)),kmul(dJ312L,PDstandardNth3g23)))))));
      
      JacPDstandardNth12g33 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g33,kmadd(J21L,PDstandardNth12g33,kmul(J31L,PDstandardNth13g33))),kmadd(J11L,kmadd(J22L,PDstandardNth12g33,kmul(J32L,PDstandardNth13g33)),kmadd(dJ112L,PDstandardNth1g33,kmadd(J22L,kmadd(J21L,PDstandardNth22g33,kmul(J31L,PDstandardNth23g33)),kmadd(dJ212L,PDstandardNth2g33,kmadd(J32L,kmadd(J21L,PDstandardNth23g33,kmul(J31L,PDstandardNth33g33)),kmul(dJ312L,PDstandardNth3g33)))))));
      
      JacPDstandardNth13alpha = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11alpha,kmadd(J21L,PDstandardNth12alpha,kmul(J31L,PDstandardNth13alpha))),kmadd(J11L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmadd(dJ113L,PDstandardNth1alpha,kmadd(J23L,kmadd(J21L,PDstandardNth22alpha,kmul(J31L,PDstandardNth23alpha)),kmadd(dJ213L,PDstandardNth2alpha,kmadd(J33L,kmadd(J21L,PDstandardNth23alpha,kmul(J31L,PDstandardNth33alpha)),kmul(dJ313L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth13g11 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g11,kmadd(J21L,PDstandardNth12g11,kmul(J31L,PDstandardNth13g11))),kmadd(J11L,kmadd(J23L,PDstandardNth12g11,kmul(J33L,PDstandardNth13g11)),kmadd(dJ113L,PDstandardNth1g11,kmadd(J23L,kmadd(J21L,PDstandardNth22g11,kmul(J31L,PDstandardNth23g11)),kmadd(dJ213L,PDstandardNth2g11,kmadd(J33L,kmadd(J21L,PDstandardNth23g11,kmul(J31L,PDstandardNth33g11)),kmul(dJ313L,PDstandardNth3g11)))))));
      
      JacPDstandardNth13g12 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g12,kmadd(J21L,PDstandardNth12g12,kmul(J31L,PDstandardNth13g12))),kmadd(J11L,kmadd(J23L,PDstandardNth12g12,kmul(J33L,PDstandardNth13g12)),kmadd(dJ113L,PDstandardNth1g12,kmadd(J23L,kmadd(J21L,PDstandardNth22g12,kmul(J31L,PDstandardNth23g12)),kmadd(dJ213L,PDstandardNth2g12,kmadd(J33L,kmadd(J21L,PDstandardNth23g12,kmul(J31L,PDstandardNth33g12)),kmul(dJ313L,PDstandardNth3g12)))))));
      
      JacPDstandardNth13g13 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g13,kmadd(J21L,PDstandardNth12g13,kmul(J31L,PDstandardNth13g13))),kmadd(J11L,kmadd(J23L,PDstandardNth12g13,kmul(J33L,PDstandardNth13g13)),kmadd(dJ113L,PDstandardNth1g13,kmadd(J23L,kmadd(J21L,PDstandardNth22g13,kmul(J31L,PDstandardNth23g13)),kmadd(dJ213L,PDstandardNth2g13,kmadd(J33L,kmadd(J21L,PDstandardNth23g13,kmul(J31L,PDstandardNth33g13)),kmul(dJ313L,PDstandardNth3g13)))))));
      
      JacPDstandardNth13g22 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g22,kmadd(J21L,PDstandardNth12g22,kmul(J31L,PDstandardNth13g22))),kmadd(J11L,kmadd(J23L,PDstandardNth12g22,kmul(J33L,PDstandardNth13g22)),kmadd(dJ113L,PDstandardNth1g22,kmadd(J23L,kmadd(J21L,PDstandardNth22g22,kmul(J31L,PDstandardNth23g22)),kmadd(dJ213L,PDstandardNth2g22,kmadd(J33L,kmadd(J21L,PDstandardNth23g22,kmul(J31L,PDstandardNth33g22)),kmul(dJ313L,PDstandardNth3g22)))))));
      
      JacPDstandardNth13g23 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g23,kmadd(J21L,PDstandardNth12g23,kmul(J31L,PDstandardNth13g23))),kmadd(J11L,kmadd(J23L,PDstandardNth12g23,kmul(J33L,PDstandardNth13g23)),kmadd(dJ113L,PDstandardNth1g23,kmadd(J23L,kmadd(J21L,PDstandardNth22g23,kmul(J31L,PDstandardNth23g23)),kmadd(dJ213L,PDstandardNth2g23,kmadd(J33L,kmadd(J21L,PDstandardNth23g23,kmul(J31L,PDstandardNth33g23)),kmul(dJ313L,PDstandardNth3g23)))))));
      
      JacPDstandardNth13g33 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g33,kmadd(J21L,PDstandardNth12g33,kmul(J31L,PDstandardNth13g33))),kmadd(J11L,kmadd(J23L,PDstandardNth12g33,kmul(J33L,PDstandardNth13g33)),kmadd(dJ113L,PDstandardNth1g33,kmadd(J23L,kmadd(J21L,PDstandardNth22g33,kmul(J31L,PDstandardNth23g33)),kmadd(dJ213L,PDstandardNth2g33,kmadd(J33L,kmadd(J21L,PDstandardNth23g33,kmul(J31L,PDstandardNth33g33)),kmul(dJ313L,PDstandardNth3g33)))))));
      
      JacPDstandardNth21g11 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g11,kmadd(J21L,PDstandardNth12g11,kmul(J31L,PDstandardNth13g11))),kmadd(J11L,kmadd(J22L,PDstandardNth12g11,kmul(J32L,PDstandardNth13g11)),kmadd(dJ112L,PDstandardNth1g11,kmadd(J22L,kmadd(J21L,PDstandardNth22g11,kmul(J31L,PDstandardNth23g11)),kmadd(dJ212L,PDstandardNth2g11,kmadd(J32L,kmadd(J21L,PDstandardNth23g11,kmul(J31L,PDstandardNth33g11)),kmul(dJ312L,PDstandardNth3g11)))))));
      
      JacPDstandardNth21g12 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g12,kmadd(J21L,PDstandardNth12g12,kmul(J31L,PDstandardNth13g12))),kmadd(J11L,kmadd(J22L,PDstandardNth12g12,kmul(J32L,PDstandardNth13g12)),kmadd(dJ112L,PDstandardNth1g12,kmadd(J22L,kmadd(J21L,PDstandardNth22g12,kmul(J31L,PDstandardNth23g12)),kmadd(dJ212L,PDstandardNth2g12,kmadd(J32L,kmadd(J21L,PDstandardNth23g12,kmul(J31L,PDstandardNth33g12)),kmul(dJ312L,PDstandardNth3g12)))))));
      
      JacPDstandardNth21g13 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g13,kmadd(J21L,PDstandardNth12g13,kmul(J31L,PDstandardNth13g13))),kmadd(J11L,kmadd(J22L,PDstandardNth12g13,kmul(J32L,PDstandardNth13g13)),kmadd(dJ112L,PDstandardNth1g13,kmadd(J22L,kmadd(J21L,PDstandardNth22g13,kmul(J31L,PDstandardNth23g13)),kmadd(dJ212L,PDstandardNth2g13,kmadd(J32L,kmadd(J21L,PDstandardNth23g13,kmul(J31L,PDstandardNth33g13)),kmul(dJ312L,PDstandardNth3g13)))))));
      
      JacPDstandardNth21g22 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g22,kmadd(J21L,PDstandardNth12g22,kmul(J31L,PDstandardNth13g22))),kmadd(J11L,kmadd(J22L,PDstandardNth12g22,kmul(J32L,PDstandardNth13g22)),kmadd(dJ112L,PDstandardNth1g22,kmadd(J22L,kmadd(J21L,PDstandardNth22g22,kmul(J31L,PDstandardNth23g22)),kmadd(dJ212L,PDstandardNth2g22,kmadd(J32L,kmadd(J21L,PDstandardNth23g22,kmul(J31L,PDstandardNth33g22)),kmul(dJ312L,PDstandardNth3g22)))))));
      
      JacPDstandardNth21g23 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g23,kmadd(J21L,PDstandardNth12g23,kmul(J31L,PDstandardNth13g23))),kmadd(J11L,kmadd(J22L,PDstandardNth12g23,kmul(J32L,PDstandardNth13g23)),kmadd(dJ112L,PDstandardNth1g23,kmadd(J22L,kmadd(J21L,PDstandardNth22g23,kmul(J31L,PDstandardNth23g23)),kmadd(dJ212L,PDstandardNth2g23,kmadd(J32L,kmadd(J21L,PDstandardNth23g23,kmul(J31L,PDstandardNth33g23)),kmul(dJ312L,PDstandardNth3g23)))))));
      
      JacPDstandardNth21g33 = 
        kmadd(J12L,kmadd(J11L,PDstandardNth11g33,kmadd(J21L,PDstandardNth12g33,kmul(J31L,PDstandardNth13g33))),kmadd(J11L,kmadd(J22L,PDstandardNth12g33,kmul(J32L,PDstandardNth13g33)),kmadd(dJ112L,PDstandardNth1g33,kmadd(J22L,kmadd(J21L,PDstandardNth22g33,kmul(J31L,PDstandardNth23g33)),kmadd(dJ212L,PDstandardNth2g33,kmadd(J32L,kmadd(J21L,PDstandardNth23g33,kmul(J31L,PDstandardNth33g33)),kmul(dJ312L,PDstandardNth3g33)))))));
      
      JacPDstandardNth23alpha = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11alpha,kmadd(J22L,PDstandardNth12alpha,kmul(J32L,PDstandardNth13alpha))),kmadd(J12L,kmadd(J23L,PDstandardNth12alpha,kmul(J33L,PDstandardNth13alpha)),kmadd(dJ123L,PDstandardNth1alpha,kmadd(J23L,kmadd(J22L,PDstandardNth22alpha,kmul(J32L,PDstandardNth23alpha)),kmadd(dJ223L,PDstandardNth2alpha,kmadd(J33L,kmadd(J22L,PDstandardNth23alpha,kmul(J32L,PDstandardNth33alpha)),kmul(dJ323L,PDstandardNth3alpha)))))));
      
      JacPDstandardNth23g11 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g11,kmadd(J22L,PDstandardNth12g11,kmul(J32L,PDstandardNth13g11))),kmadd(J12L,kmadd(J23L,PDstandardNth12g11,kmul(J33L,PDstandardNth13g11)),kmadd(dJ123L,PDstandardNth1g11,kmadd(J23L,kmadd(J22L,PDstandardNth22g11,kmul(J32L,PDstandardNth23g11)),kmadd(dJ223L,PDstandardNth2g11,kmadd(J33L,kmadd(J22L,PDstandardNth23g11,kmul(J32L,PDstandardNth33g11)),kmul(dJ323L,PDstandardNth3g11)))))));
      
      JacPDstandardNth23g12 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g12,kmadd(J22L,PDstandardNth12g12,kmul(J32L,PDstandardNth13g12))),kmadd(J12L,kmadd(J23L,PDstandardNth12g12,kmul(J33L,PDstandardNth13g12)),kmadd(dJ123L,PDstandardNth1g12,kmadd(J23L,kmadd(J22L,PDstandardNth22g12,kmul(J32L,PDstandardNth23g12)),kmadd(dJ223L,PDstandardNth2g12,kmadd(J33L,kmadd(J22L,PDstandardNth23g12,kmul(J32L,PDstandardNth33g12)),kmul(dJ323L,PDstandardNth3g12)))))));
      
      JacPDstandardNth23g13 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g13,kmadd(J22L,PDstandardNth12g13,kmul(J32L,PDstandardNth13g13))),kmadd(J12L,kmadd(J23L,PDstandardNth12g13,kmul(J33L,PDstandardNth13g13)),kmadd(dJ123L,PDstandardNth1g13,kmadd(J23L,kmadd(J22L,PDstandardNth22g13,kmul(J32L,PDstandardNth23g13)),kmadd(dJ223L,PDstandardNth2g13,kmadd(J33L,kmadd(J22L,PDstandardNth23g13,kmul(J32L,PDstandardNth33g13)),kmul(dJ323L,PDstandardNth3g13)))))));
      
      JacPDstandardNth23g22 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g22,kmadd(J22L,PDstandardNth12g22,kmul(J32L,PDstandardNth13g22))),kmadd(J12L,kmadd(J23L,PDstandardNth12g22,kmul(J33L,PDstandardNth13g22)),kmadd(dJ123L,PDstandardNth1g22,kmadd(J23L,kmadd(J22L,PDstandardNth22g22,kmul(J32L,PDstandardNth23g22)),kmadd(dJ223L,PDstandardNth2g22,kmadd(J33L,kmadd(J22L,PDstandardNth23g22,kmul(J32L,PDstandardNth33g22)),kmul(dJ323L,PDstandardNth3g22)))))));
      
      JacPDstandardNth23g23 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g23,kmadd(J22L,PDstandardNth12g23,kmul(J32L,PDstandardNth13g23))),kmadd(J12L,kmadd(J23L,PDstandardNth12g23,kmul(J33L,PDstandardNth13g23)),kmadd(dJ123L,PDstandardNth1g23,kmadd(J23L,kmadd(J22L,PDstandardNth22g23,kmul(J32L,PDstandardNth23g23)),kmadd(dJ223L,PDstandardNth2g23,kmadd(J33L,kmadd(J22L,PDstandardNth23g23,kmul(J32L,PDstandardNth33g23)),kmul(dJ323L,PDstandardNth3g23)))))));
      
      JacPDstandardNth23g33 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g33,kmadd(J22L,PDstandardNth12g33,kmul(J32L,PDstandardNth13g33))),kmadd(J12L,kmadd(J23L,PDstandardNth12g33,kmul(J33L,PDstandardNth13g33)),kmadd(dJ123L,PDstandardNth1g33,kmadd(J23L,kmadd(J22L,PDstandardNth22g33,kmul(J32L,PDstandardNth23g33)),kmadd(dJ223L,PDstandardNth2g33,kmadd(J33L,kmadd(J22L,PDstandardNth23g33,kmul(J32L,PDstandardNth33g33)),kmul(dJ323L,PDstandardNth3g33)))))));
      
      JacPDstandardNth31g11 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g11,kmadd(J21L,PDstandardNth12g11,kmul(J31L,PDstandardNth13g11))),kmadd(J11L,kmadd(J23L,PDstandardNth12g11,kmul(J33L,PDstandardNth13g11)),kmadd(dJ113L,PDstandardNth1g11,kmadd(J23L,kmadd(J21L,PDstandardNth22g11,kmul(J31L,PDstandardNth23g11)),kmadd(dJ213L,PDstandardNth2g11,kmadd(J33L,kmadd(J21L,PDstandardNth23g11,kmul(J31L,PDstandardNth33g11)),kmul(dJ313L,PDstandardNth3g11)))))));
      
      JacPDstandardNth31g12 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g12,kmadd(J21L,PDstandardNth12g12,kmul(J31L,PDstandardNth13g12))),kmadd(J11L,kmadd(J23L,PDstandardNth12g12,kmul(J33L,PDstandardNth13g12)),kmadd(dJ113L,PDstandardNth1g12,kmadd(J23L,kmadd(J21L,PDstandardNth22g12,kmul(J31L,PDstandardNth23g12)),kmadd(dJ213L,PDstandardNth2g12,kmadd(J33L,kmadd(J21L,PDstandardNth23g12,kmul(J31L,PDstandardNth33g12)),kmul(dJ313L,PDstandardNth3g12)))))));
      
      JacPDstandardNth31g13 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g13,kmadd(J21L,PDstandardNth12g13,kmul(J31L,PDstandardNth13g13))),kmadd(J11L,kmadd(J23L,PDstandardNth12g13,kmul(J33L,PDstandardNth13g13)),kmadd(dJ113L,PDstandardNth1g13,kmadd(J23L,kmadd(J21L,PDstandardNth22g13,kmul(J31L,PDstandardNth23g13)),kmadd(dJ213L,PDstandardNth2g13,kmadd(J33L,kmadd(J21L,PDstandardNth23g13,kmul(J31L,PDstandardNth33g13)),kmul(dJ313L,PDstandardNth3g13)))))));
      
      JacPDstandardNth31g22 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g22,kmadd(J21L,PDstandardNth12g22,kmul(J31L,PDstandardNth13g22))),kmadd(J11L,kmadd(J23L,PDstandardNth12g22,kmul(J33L,PDstandardNth13g22)),kmadd(dJ113L,PDstandardNth1g22,kmadd(J23L,kmadd(J21L,PDstandardNth22g22,kmul(J31L,PDstandardNth23g22)),kmadd(dJ213L,PDstandardNth2g22,kmadd(J33L,kmadd(J21L,PDstandardNth23g22,kmul(J31L,PDstandardNth33g22)),kmul(dJ313L,PDstandardNth3g22)))))));
      
      JacPDstandardNth31g23 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g23,kmadd(J21L,PDstandardNth12g23,kmul(J31L,PDstandardNth13g23))),kmadd(J11L,kmadd(J23L,PDstandardNth12g23,kmul(J33L,PDstandardNth13g23)),kmadd(dJ113L,PDstandardNth1g23,kmadd(J23L,kmadd(J21L,PDstandardNth22g23,kmul(J31L,PDstandardNth23g23)),kmadd(dJ213L,PDstandardNth2g23,kmadd(J33L,kmadd(J21L,PDstandardNth23g23,kmul(J31L,PDstandardNth33g23)),kmul(dJ313L,PDstandardNth3g23)))))));
      
      JacPDstandardNth31g33 = 
        kmadd(J13L,kmadd(J11L,PDstandardNth11g33,kmadd(J21L,PDstandardNth12g33,kmul(J31L,PDstandardNth13g33))),kmadd(J11L,kmadd(J23L,PDstandardNth12g33,kmul(J33L,PDstandardNth13g33)),kmadd(dJ113L,PDstandardNth1g33,kmadd(J23L,kmadd(J21L,PDstandardNth22g33,kmul(J31L,PDstandardNth23g33)),kmadd(dJ213L,PDstandardNth2g33,kmadd(J33L,kmadd(J21L,PDstandardNth23g33,kmul(J31L,PDstandardNth33g33)),kmul(dJ313L,PDstandardNth3g33)))))));
      
      JacPDstandardNth32g11 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g11,kmadd(J22L,PDstandardNth12g11,kmul(J32L,PDstandardNth13g11))),kmadd(J12L,kmadd(J23L,PDstandardNth12g11,kmul(J33L,PDstandardNth13g11)),kmadd(dJ123L,PDstandardNth1g11,kmadd(J23L,kmadd(J22L,PDstandardNth22g11,kmul(J32L,PDstandardNth23g11)),kmadd(dJ223L,PDstandardNth2g11,kmadd(J33L,kmadd(J22L,PDstandardNth23g11,kmul(J32L,PDstandardNth33g11)),kmul(dJ323L,PDstandardNth3g11)))))));
      
      JacPDstandardNth32g12 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g12,kmadd(J22L,PDstandardNth12g12,kmul(J32L,PDstandardNth13g12))),kmadd(J12L,kmadd(J23L,PDstandardNth12g12,kmul(J33L,PDstandardNth13g12)),kmadd(dJ123L,PDstandardNth1g12,kmadd(J23L,kmadd(J22L,PDstandardNth22g12,kmul(J32L,PDstandardNth23g12)),kmadd(dJ223L,PDstandardNth2g12,kmadd(J33L,kmadd(J22L,PDstandardNth23g12,kmul(J32L,PDstandardNth33g12)),kmul(dJ323L,PDstandardNth3g12)))))));
      
      JacPDstandardNth32g13 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g13,kmadd(J22L,PDstandardNth12g13,kmul(J32L,PDstandardNth13g13))),kmadd(J12L,kmadd(J23L,PDstandardNth12g13,kmul(J33L,PDstandardNth13g13)),kmadd(dJ123L,PDstandardNth1g13,kmadd(J23L,kmadd(J22L,PDstandardNth22g13,kmul(J32L,PDstandardNth23g13)),kmadd(dJ223L,PDstandardNth2g13,kmadd(J33L,kmadd(J22L,PDstandardNth23g13,kmul(J32L,PDstandardNth33g13)),kmul(dJ323L,PDstandardNth3g13)))))));
      
      JacPDstandardNth32g22 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g22,kmadd(J22L,PDstandardNth12g22,kmul(J32L,PDstandardNth13g22))),kmadd(J12L,kmadd(J23L,PDstandardNth12g22,kmul(J33L,PDstandardNth13g22)),kmadd(dJ123L,PDstandardNth1g22,kmadd(J23L,kmadd(J22L,PDstandardNth22g22,kmul(J32L,PDstandardNth23g22)),kmadd(dJ223L,PDstandardNth2g22,kmadd(J33L,kmadd(J22L,PDstandardNth23g22,kmul(J32L,PDstandardNth33g22)),kmul(dJ323L,PDstandardNth3g22)))))));
      
      JacPDstandardNth32g23 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g23,kmadd(J22L,PDstandardNth12g23,kmul(J32L,PDstandardNth13g23))),kmadd(J12L,kmadd(J23L,PDstandardNth12g23,kmul(J33L,PDstandardNth13g23)),kmadd(dJ123L,PDstandardNth1g23,kmadd(J23L,kmadd(J22L,PDstandardNth22g23,kmul(J32L,PDstandardNth23g23)),kmadd(dJ223L,PDstandardNth2g23,kmadd(J33L,kmadd(J22L,PDstandardNth23g23,kmul(J32L,PDstandardNth33g23)),kmul(dJ323L,PDstandardNth3g23)))))));
      
      JacPDstandardNth32g33 = 
        kmadd(J13L,kmadd(J12L,PDstandardNth11g33,kmadd(J22L,PDstandardNth12g33,kmul(J32L,PDstandardNth13g33))),kmadd(J12L,kmadd(J23L,PDstandardNth12g33,kmul(J33L,PDstandardNth13g33)),kmadd(dJ123L,PDstandardNth1g33,kmadd(J23L,kmadd(J22L,PDstandardNth22g33,kmul(J32L,PDstandardNth23g33)),kmadd(dJ223L,PDstandardNth2g33,kmadd(J33L,kmadd(J22L,PDstandardNth23g33,kmul(J32L,PDstandardNth33g33)),kmul(dJ323L,PDstandardNth3g33)))))));
    }
    else
    {
      JacPDstandardNth1alpha = PDstandardNth1alpha;
      
      JacPDstandardNth1beta1 = PDstandardNth1beta1;
      
      JacPDstandardNth1beta2 = PDstandardNth1beta2;
      
      JacPDstandardNth1beta3 = PDstandardNth1beta3;
      
      JacPDstandardNth1g11 = PDstandardNth1g11;
      
      JacPDstandardNth1g12 = PDstandardNth1g12;
      
      JacPDstandardNth1g13 = PDstandardNth1g13;
      
      JacPDstandardNth1g22 = PDstandardNth1g22;
      
      JacPDstandardNth1g23 = PDstandardNth1g23;
      
      JacPDstandardNth1g33 = PDstandardNth1g33;
      
      JacPDstandardNth1K11 = PDstandardNth1K11;
      
      JacPDstandardNth1K12 = PDstandardNth1K12;
      
      JacPDstandardNth1K13 = PDstandardNth1K13;
      
      JacPDstandardNth1K22 = PDstandardNth1K22;
      
      JacPDstandardNth1K23 = PDstandardNth1K23;
      
      JacPDstandardNth1K33 = PDstandardNth1K33;
      
      JacPDstandardNth2alpha = PDstandardNth2alpha;
      
      JacPDstandardNth2beta1 = PDstandardNth2beta1;
      
      JacPDstandardNth2beta2 = PDstandardNth2beta2;
      
      JacPDstandardNth2beta3 = PDstandardNth2beta3;
      
      JacPDstandardNth2g11 = PDstandardNth2g11;
      
      JacPDstandardNth2g12 = PDstandardNth2g12;
      
      JacPDstandardNth2g13 = PDstandardNth2g13;
      
      JacPDstandardNth2g22 = PDstandardNth2g22;
      
      JacPDstandardNth2g23 = PDstandardNth2g23;
      
      JacPDstandardNth2g33 = PDstandardNth2g33;
      
      JacPDstandardNth2K11 = PDstandardNth2K11;
      
      JacPDstandardNth2K12 = PDstandardNth2K12;
      
      JacPDstandardNth2K13 = PDstandardNth2K13;
      
      JacPDstandardNth2K22 = PDstandardNth2K22;
      
      JacPDstandardNth2K23 = PDstandardNth2K23;
      
      JacPDstandardNth2K33 = PDstandardNth2K33;
      
      JacPDstandardNth3alpha = PDstandardNth3alpha;
      
      JacPDstandardNth3beta1 = PDstandardNth3beta1;
      
      JacPDstandardNth3beta2 = PDstandardNth3beta2;
      
      JacPDstandardNth3beta3 = PDstandardNth3beta3;
      
      JacPDstandardNth3g11 = PDstandardNth3g11;
      
      JacPDstandardNth3g12 = PDstandardNth3g12;
      
      JacPDstandardNth3g13 = PDstandardNth3g13;
      
      JacPDstandardNth3g22 = PDstandardNth3g22;
      
      JacPDstandardNth3g23 = PDstandardNth3g23;
      
      JacPDstandardNth3g33 = PDstandardNth3g33;
      
      JacPDstandardNth3K11 = PDstandardNth3K11;
      
      JacPDstandardNth3K12 = PDstandardNth3K12;
      
      JacPDstandardNth3K13 = PDstandardNth3K13;
      
      JacPDstandardNth3K22 = PDstandardNth3K22;
      
      JacPDstandardNth3K23 = PDstandardNth3K23;
      
      JacPDstandardNth3K33 = PDstandardNth3K33;
      
      JacPDstandardNth11alpha = PDstandardNth11alpha;
      
      JacPDstandardNth11g22 = PDstandardNth11g22;
      
      JacPDstandardNth11g23 = PDstandardNth11g23;
      
      JacPDstandardNth11g33 = PDstandardNth11g33;
      
      JacPDstandardNth22alpha = PDstandardNth22alpha;
      
      JacPDstandardNth22g11 = PDstandardNth22g11;
      
      JacPDstandardNth22g13 = PDstandardNth22g13;
      
      JacPDstandardNth22g33 = PDstandardNth22g33;
      
      JacPDstandardNth33alpha = PDstandardNth33alpha;
      
      JacPDstandardNth33g11 = PDstandardNth33g11;
      
      JacPDstandardNth33g12 = PDstandardNth33g12;
      
      JacPDstandardNth33g22 = PDstandardNth33g22;
      
      JacPDstandardNth12alpha = PDstandardNth12alpha;
      
      JacPDstandardNth12g11 = PDstandardNth12g11;
      
      JacPDstandardNth12g12 = PDstandardNth12g12;
      
      JacPDstandardNth12g13 = PDstandardNth12g13;
      
      JacPDstandardNth12g22 = PDstandardNth12g22;
      
      JacPDstandardNth12g23 = PDstandardNth12g23;
      
      JacPDstandardNth12g33 = PDstandardNth12g33;
      
      JacPDstandardNth13alpha = PDstandardNth13alpha;
      
      JacPDstandardNth13g11 = PDstandardNth13g11;
      
      JacPDstandardNth13g12 = PDstandardNth13g12;
      
      JacPDstandardNth13g13 = PDstandardNth13g13;
      
      JacPDstandardNth13g22 = PDstandardNth13g22;
      
      JacPDstandardNth13g23 = PDstandardNth13g23;
      
      JacPDstandardNth13g33 = PDstandardNth13g33;
      
      JacPDstandardNth21g11 = PDstandardNth12g11;
      
      JacPDstandardNth21g12 = PDstandardNth12g12;
      
      JacPDstandardNth21g13 = PDstandardNth12g13;
      
      JacPDstandardNth21g22 = PDstandardNth12g22;
      
      JacPDstandardNth21g23 = PDstandardNth12g23;
      
      JacPDstandardNth21g33 = PDstandardNth12g33;
      
      JacPDstandardNth23alpha = PDstandardNth23alpha;
      
      JacPDstandardNth23g11 = PDstandardNth23g11;
      
      JacPDstandardNth23g12 = PDstandardNth23g12;
      
      JacPDstandardNth23g13 = PDstandardNth23g13;
      
      JacPDstandardNth23g22 = PDstandardNth23g22;
      
      JacPDstandardNth23g23 = PDstandardNth23g23;
      
      JacPDstandardNth23g33 = PDstandardNth23g33;
      
      JacPDstandardNth31g11 = PDstandardNth13g11;
      
      JacPDstandardNth31g12 = PDstandardNth13g12;
      
      JacPDstandardNth31g13 = PDstandardNth13g13;
      
      JacPDstandardNth31g22 = PDstandardNth13g22;
      
      JacPDstandardNth31g23 = PDstandardNth13g23;
      
      JacPDstandardNth31g33 = PDstandardNth13g33;
      
      JacPDstandardNth32g11 = PDstandardNth23g11;
      
      JacPDstandardNth32g12 = PDstandardNth23g12;
      
      JacPDstandardNth32g13 = PDstandardNth23g13;
      
      JacPDstandardNth32g22 = PDstandardNth23g22;
      
      JacPDstandardNth32g23 = PDstandardNth23g23;
      
      JacPDstandardNth32g33 = PDstandardNth23g33;
    }
    
    CCTK_REAL_VEC detg = 
      knmsub(g22L,SQR(g13L),knmsub(g11L,SQR(g23L),kmadd(g33L,kmsub(g11L,g22L,SQR(g12L)),kmul(g12L,kmul(g13L,kmul(g23L,ToReal(2)))))));
    
    CCTK_REAL_VEC gu11 = kmul(INV(detg),kmsub(g22L,g33L,SQR(g23L)));
    
    CCTK_REAL_VEC gu12 = kmul(INV(detg),kmsub(g13L,g23L,kmul(g12L,g33L)));
    
    CCTK_REAL_VEC gu13 = kmul(INV(detg),kmsub(g12L,g23L,kmul(g13L,g22L)));
    
    CCTK_REAL_VEC gu21 = kmul(INV(detg),kmsub(g13L,g23L,kmul(g12L,g33L)));
    
    CCTK_REAL_VEC gu22 = kmul(INV(detg),kmsub(g11L,g33L,SQR(g13L)));
    
    CCTK_REAL_VEC gu23 = kmul(INV(detg),kmsub(g12L,g13L,kmul(g11L,g23L)));
    
    CCTK_REAL_VEC gu31 = kmul(INV(detg),kmsub(g12L,g23L,kmul(g13L,g22L)));
    
    CCTK_REAL_VEC gu32 = kmul(INV(detg),kmsub(g12L,g13L,kmul(g11L,g23L)));
    
    CCTK_REAL_VEC gu33 = kmul(INV(detg),kmsub(g11L,g22L,SQR(g12L)));
    
    CCTK_REAL_VEC G111 = 
      kmul(ToReal(0.5),kmadd(gu11,JacPDstandardNth1g11,knmsub(gu12,JacPDstandardNth2g11,kmsub(kmadd(gu12,JacPDstandardNth1g12,kmul(gu13,JacPDstandardNth1g13)),ToReal(2),kmul(gu13,JacPDstandardNth3g11)))));
    
    CCTK_REAL_VEC G211 = 
      kmul(ToReal(0.5),kmadd(gu21,JacPDstandardNth1g11,knmsub(gu22,JacPDstandardNth2g11,kmsub(kmadd(gu22,JacPDstandardNth1g12,kmul(gu23,JacPDstandardNth1g13)),ToReal(2),kmul(gu23,JacPDstandardNth3g11)))));
    
    CCTK_REAL_VEC G311 = 
      kmul(ToReal(0.5),kmadd(gu31,JacPDstandardNth1g11,knmsub(gu32,JacPDstandardNth2g11,kmsub(kmadd(gu32,JacPDstandardNth1g12,kmul(gu33,JacPDstandardNth1g13)),ToReal(2),kmul(gu33,JacPDstandardNth3g11)))));
    
    CCTK_REAL_VEC G112 = 
      kmul(kmadd(gu12,JacPDstandardNth1g22,kmadd(gu11,JacPDstandardNth2g11,kmul(gu13,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth2g13,JacPDstandardNth3g12))))),ToReal(0.5));
    
    CCTK_REAL_VEC G212 = 
      kmul(kmadd(gu22,JacPDstandardNth1g22,kmadd(gu21,JacPDstandardNth2g11,kmul(gu23,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth2g13,JacPDstandardNth3g12))))),ToReal(0.5));
    
    CCTK_REAL_VEC G312 = 
      kmul(kmadd(gu32,JacPDstandardNth1g22,kmadd(gu31,JacPDstandardNth2g11,kmul(gu33,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth2g13,JacPDstandardNth3g12))))),ToReal(0.5));
    
    CCTK_REAL_VEC G113 = 
      kmul(kmadd(gu13,JacPDstandardNth1g33,kmadd(gu11,JacPDstandardNth3g11,kmul(gu12,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth3g12,JacPDstandardNth2g13))))),ToReal(0.5));
    
    CCTK_REAL_VEC G213 = 
      kmul(kmadd(gu23,JacPDstandardNth1g33,kmadd(gu21,JacPDstandardNth3g11,kmul(gu22,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth3g12,JacPDstandardNth2g13))))),ToReal(0.5));
    
    CCTK_REAL_VEC G313 = 
      kmul(kmadd(gu33,JacPDstandardNth1g33,kmadd(gu31,JacPDstandardNth3g11,kmul(gu32,kadd(JacPDstandardNth1g23,ksub(JacPDstandardNth3g12,JacPDstandardNth2g13))))),ToReal(0.5));
    
    CCTK_REAL_VEC G122 = 
      kmul(ToReal(0.5),kmadd(gu12,JacPDstandardNth2g22,kmadd(gu11,kmsub(JacPDstandardNth2g12,ToReal(2),JacPDstandardNth1g22),kmul(gu13,kmsub(JacPDstandardNth2g23,ToReal(2),JacPDstandardNth3g22)))));
    
    CCTK_REAL_VEC G222 = 
      kmul(ToReal(0.5),kmadd(gu22,JacPDstandardNth2g22,kmadd(gu21,kmsub(JacPDstandardNth2g12,ToReal(2),JacPDstandardNth1g22),kmul(gu23,kmsub(JacPDstandardNth2g23,ToReal(2),JacPDstandardNth3g22)))));
    
    CCTK_REAL_VEC G322 = 
      kmul(ToReal(0.5),kmadd(gu32,JacPDstandardNth2g22,kmadd(gu31,kmsub(JacPDstandardNth2g12,ToReal(2),JacPDstandardNth1g22),kmul(gu33,kmsub(JacPDstandardNth2g23,ToReal(2),JacPDstandardNth3g22)))));
    
    CCTK_REAL_VEC G123 = 
      kmul(kmadd(gu13,JacPDstandardNth2g33,kmadd(gu12,JacPDstandardNth3g22,kmul(gu11,kadd(JacPDstandardNth2g13,ksub(JacPDstandardNth3g12,JacPDstandardNth1g23))))),ToReal(0.5));
    
    CCTK_REAL_VEC G223 = 
      kmul(kmadd(gu23,JacPDstandardNth2g33,kmadd(gu22,JacPDstandardNth3g22,kmul(gu21,kadd(JacPDstandardNth2g13,ksub(JacPDstandardNth3g12,JacPDstandardNth1g23))))),ToReal(0.5));
    
    CCTK_REAL_VEC G323 = 
      kmul(kmadd(gu33,JacPDstandardNth2g33,kmadd(gu32,JacPDstandardNth3g22,kmul(gu31,kadd(JacPDstandardNth2g13,ksub(JacPDstandardNth3g12,JacPDstandardNth1g23))))),ToReal(0.5));
    
    CCTK_REAL_VEC G133 = 
      kmul(ToReal(0.5),kmadd(gu13,JacPDstandardNth3g33,kmadd(gu11,kmsub(JacPDstandardNth3g13,ToReal(2),JacPDstandardNth1g33),kmul(gu12,kmsub(JacPDstandardNth3g23,ToReal(2),JacPDstandardNth2g33)))));
    
    CCTK_REAL_VEC G233 = 
      kmul(ToReal(0.5),kmadd(gu23,JacPDstandardNth3g33,kmadd(gu21,kmsub(JacPDstandardNth3g13,ToReal(2),JacPDstandardNth1g33),kmul(gu22,kmsub(JacPDstandardNth3g23,ToReal(2),JacPDstandardNth2g33)))));
    
    CCTK_REAL_VEC G333 = 
      kmul(ToReal(0.5),kmadd(gu33,JacPDstandardNth3g33,kmadd(gu31,kmsub(JacPDstandardNth3g13,ToReal(2),JacPDstandardNth1g33),kmul(gu32,kmsub(JacPDstandardNth3g23,ToReal(2),JacPDstandardNth2g33)))));
    
    CCTK_REAL_VEC R11 = 
      kmul(ToReal(0.5),kmadd(gu23,JacPDstandardNth31g12,kmadd(gu32,JacPDstandardNth31g12,kmadd(G111,kmul(kadd(G212,G313),ToReal(-2)),knmsub(JacPDstandardNth11g23,kadd(gu32,gu23),kmadd(gu12,ksub(JacPDstandardNth21g11,JacPDstandardNth12g11),kmadd(gu13,ksub(JacPDstandardNth31g11,JacPDstandardNth13g11),kmadd(gu23,ksub(JacPDstandardNth21g13,JacPDstandardNth23g11),kmadd(gu32,ksub(JacPDstandardNth21g13,JacPDstandardNth32g11),kmadd(SQR(G212),ToReal(2),kmadd(SQR(G313),ToReal(2),kmadd(G211,kmadd(G222,ToReal(-2),kmadd(G323,ToReal(-2),kmul(G112,ToReal(2)))),kmadd(G311,kmadd(G223,ToReal(-2),kmadd(G333,ToReal(-2),kmul(G113,ToReal(2)))),kmadd(gu22,ksub(kmsub(JacPDstandardNth21g12,ToReal(2),JacPDstandardNth22g11),JacPDstandardNth11g22),kmadd(gu33,ksub(kmsub(JacPDstandardNth31g13,ToReal(2),JacPDstandardNth33g11),JacPDstandardNth11g33),kmul(G213,kmul(G312,ToReal(4))))))))))))))))));
    
    CCTK_REAL_VEC R12 = 
      kmul(ToReal(0.5),kmadd(gu12,kadd(JacPDstandardNth11g22,JacPDstandardNth22g11),kmadd(gu32,JacPDstandardNth31g22,kmadd(gu13,JacPDstandardNth32g11,kmadd(gu23,JacPDstandardNth32g12,kmadd(gu33,JacPDstandardNth32g13,kmadd(kmadd(G112,kadd(G212,G313),kmadd(G212,G323,kmadd(G312,G333,kmul(gu12,JacPDstandardNth12g12)))),ToReal(-2),knmsub(JacPDstandardNth12g23,kadd(gu32,gu23),kmadd(gu22,ksub(JacPDstandardNth21g22,JacPDstandardNth12g22),kmadd(gu13,ksub(JacPDstandardNth11g23,kadd(JacPDstandardNth13g12,JacPDstandardNth12g13)),kmadd(gu23,ksub(JacPDstandardNth21g23,JacPDstandardNth23g12),kmadd(gu32,ksub(JacPDstandardNth22g13,JacPDstandardNth32g12),kmadd(gu33,ksub(JacPDstandardNth31g23,kadd(JacPDstandardNth33g12,JacPDstandardNth12g33)),kmul(kmadd(G122,G211,kmadd(G123,G311,kmadd(G213,G322,kmul(G313,G323)))),ToReal(2)))))))))))))));
    
    CCTK_REAL_VEC R13 = 
      kmul(ToReal(0.5),kmadd(gu22,JacPDstandardNth23g12,kmadd(gu32,JacPDstandardNth31g23,kmadd(gu13,kadd(JacPDstandardNth11g33,JacPDstandardNth33g11),kmadd(gu23,JacPDstandardNth33g12,kmadd(kmadd(G213,G222,kmadd(G223,G313,kmadd(G113,kadd(G212,G313),kmul(gu13,JacPDstandardNth13g13)))),ToReal(-2),knmsub(JacPDstandardNth13g23,kadd(gu32,gu23),kmadd(gu12,kadd(JacPDstandardNth11g23,ksub(JacPDstandardNth23g11,kadd(JacPDstandardNth13g12,JacPDstandardNth12g13))),kmadd(gu33,ksub(JacPDstandardNth31g33,JacPDstandardNth13g33),kmadd(gu22,ksub(JacPDstandardNth21g23,kadd(JacPDstandardNth22g13,JacPDstandardNth13g22)),kmadd(gu23,ksub(JacPDstandardNth21g33,JacPDstandardNth23g13),kmadd(gu32,ksub(JacPDstandardNth23g13,JacPDstandardNth32g13),kmul(kmadd(G123,G211,kmadd(G212,G223,kmadd(G133,G311,kmul(G233,G312)))),ToReal(2))))))))))))));
    
    CCTK_REAL_VEC R22 = 
      kmul(ToReal(0.5),kmadd(kadd(gu13,gu31),JacPDstandardNth12g23,kmadd(gu13,JacPDstandardNth32g12,kmadd(gu31,JacPDstandardNth32g12,kmadd(kmadd(G112,G222,kmadd(G113,G322,kmul(G222,G323))),ToReal(-2),kmadd(gu21,ksub(JacPDstandardNth12g22,JacPDstandardNth21g22),knmsub(gu13,kadd(JacPDstandardNth22g13,JacPDstandardNth13g22),kmadd(gu23,ksub(JacPDstandardNth32g22,JacPDstandardNth23g22),knmsub(gu31,kadd(JacPDstandardNth31g22,JacPDstandardNth22g13),kmadd(SQR(G112),ToReal(2),kmadd(SQR(G323),ToReal(2),kmadd(G122,kmadd(G111,ToReal(-2),kmadd(G313,ToReal(-2),kmul(G212,ToReal(2)))),kmadd(G322,kmadd(G333,ToReal(-2),kmul(G223,ToReal(2))),kmadd(gu11,ksub(kmsub(JacPDstandardNth12g12,ToReal(2),JacPDstandardNth22g11),JacPDstandardNth11g22),kmadd(gu33,ksub(kmsub(JacPDstandardNth32g23,ToReal(2),JacPDstandardNth33g22),JacPDstandardNth22g33),kmul(G123,kmul(G312,ToReal(4))))))))))))))))));
    
    CCTK_REAL_VEC R23 = 
      kmul(ToReal(0.5),kmadd(gu23,kadd(JacPDstandardNth22g33,JacPDstandardNth33g22),kmadd(kmadd(G111,G123,kmadd(kadd(G113,G223),G323,kmul(gu23,JacPDstandardNth23g23))),ToReal(-2),kmadd(gu11,kadd(JacPDstandardNth12g13,ksub(JacPDstandardNth13g12,kadd(JacPDstandardNth23g11,JacPDstandardNth11g23))),kmadd(gu21,kadd(JacPDstandardNth13g22,ksub(JacPDstandardNth22g13,kadd(JacPDstandardNth23g12,JacPDstandardNth21g23))),kmadd(gu13,kadd(JacPDstandardNth12g33,ksub(JacPDstandardNth33g12,kadd(JacPDstandardNth23g13,JacPDstandardNth13g23))),kmadd(gu33,ksub(JacPDstandardNth32g33,JacPDstandardNth23g33),kmadd(gu31,kadd(JacPDstandardNth13g23,ksub(JacPDstandardNth32g13,kadd(JacPDstandardNth31g23,JacPDstandardNth23g13))),kmul(kmadd(G122,G213,kmadd(G133,G312,kmadd(G233,G322,kmul(G112,ksub(G113,G223))))),ToReal(2))))))))));
    
    CCTK_REAL_VEC R33 = 
      kmul(ToReal(0.5),kmadd(gu12,JacPDstandardNth23g13,kmadd(gu21,JacPDstandardNth23g13,kmadd(kmadd(G111,G133,kmadd(G133,G212,kmadd(G112,G233,kmadd(G222,G233,kmadd(G113,G333,kmul(G223,G333)))))),ToReal(-2),kmadd(gu31,ksub(JacPDstandardNth13g33,JacPDstandardNth31g33),kmadd(gu32,ksub(JacPDstandardNth23g33,JacPDstandardNth32g33),knmsub(gu11,kadd(JacPDstandardNth33g11,JacPDstandardNth11g33),kmadd(gu12,ksub(JacPDstandardNth13g23,kadd(JacPDstandardNth33g12,JacPDstandardNth12g33)),kmadd(gu21,ksub(JacPDstandardNth13g23,kadd(JacPDstandardNth33g12,JacPDstandardNth21g33)),kmadd(kmadd(G133,G313,kmadd(G233,G323,kmul(gu11,JacPDstandardNth13g13))),ToReal(2),kmadd(SQR(G113),ToReal(2),kmadd(SQR(G223),ToReal(2),kmadd(gu22,ksub(kmsub(JacPDstandardNth23g23,ToReal(2),JacPDstandardNth33g22),JacPDstandardNth22g33),kmul(G123,kmul(G213,ToReal(4))))))))))))))));
    
    CCTK_REAL_VEC Km11 = 
      kmadd(gu11,K11L,kmadd(gu12,K12L,kmul(gu13,K13L)));
    
    CCTK_REAL_VEC Km21 = 
      kmadd(gu21,K11L,kmadd(gu22,K12L,kmul(gu23,K13L)));
    
    CCTK_REAL_VEC Km31 = 
      kmadd(gu31,K11L,kmadd(gu32,K12L,kmul(gu33,K13L)));
    
    CCTK_REAL_VEC Km12 = 
      kmadd(gu11,K12L,kmadd(gu12,K22L,kmul(gu13,K23L)));
    
    CCTK_REAL_VEC Km22 = 
      kmadd(gu21,K12L,kmadd(gu22,K22L,kmul(gu23,K23L)));
    
    CCTK_REAL_VEC Km32 = 
      kmadd(gu31,K12L,kmadd(gu32,K22L,kmul(gu33,K23L)));
    
    CCTK_REAL_VEC Km13 = 
      kmadd(gu11,K13L,kmadd(gu12,K23L,kmul(gu13,K33L)));
    
    CCTK_REAL_VEC Km23 = 
      kmadd(gu21,K13L,kmadd(gu22,K23L,kmul(gu23,K33L)));
    
    CCTK_REAL_VEC Km33 = 
      kmadd(gu31,K13L,kmadd(gu32,K23L,kmul(gu33,K33L)));
    
    CCTK_REAL_VEC trK = kadd(Km11,kadd(Km22,Km33));
    
    CCTK_REAL_VEC g11rhsL = 
      kmadd(beta1L,JacPDstandardNth1g11,kmadd(beta2L,JacPDstandardNth2g11,kmadd(beta3L,JacPDstandardNth3g11,kmadd(alphaL,kmul(K11L,ToReal(-2)),kmul(kmadd(g11L,JacPDstandardNth1beta1,kmadd(g12L,JacPDstandardNth1beta2,kmul(g13L,JacPDstandardNth1beta3))),ToReal(2))))));
    
    CCTK_REAL_VEC g12rhsL = 
      kmadd(g22L,JacPDstandardNth1beta2,kmadd(g23L,JacPDstandardNth1beta3,kmadd(beta1L,JacPDstandardNth1g12,kmadd(g11L,JacPDstandardNth2beta1,kmadd(g12L,kadd(JacPDstandardNth1beta1,JacPDstandardNth2beta2),kmadd(g13L,JacPDstandardNth2beta3,kmadd(beta2L,JacPDstandardNth2g12,kmadd(beta3L,JacPDstandardNth3g12,kmul(alphaL,kmul(K12L,ToReal(-2)))))))))));
    
    CCTK_REAL_VEC g13rhsL = 
      kmadd(g23L,JacPDstandardNth1beta2,kmadd(g33L,JacPDstandardNth1beta3,kmadd(beta1L,JacPDstandardNth1g13,kmadd(beta2L,JacPDstandardNth2g13,kmadd(g11L,JacPDstandardNth3beta1,kmadd(g12L,JacPDstandardNth3beta2,kmadd(g13L,kadd(JacPDstandardNth1beta1,JacPDstandardNth3beta3),kmadd(beta3L,JacPDstandardNth3g13,kmul(alphaL,kmul(K13L,ToReal(-2)))))))))));
    
    CCTK_REAL_VEC g22rhsL = 
      kmadd(beta1L,JacPDstandardNth1g22,kmadd(beta2L,JacPDstandardNth2g22,kmadd(beta3L,JacPDstandardNth3g22,kmadd(alphaL,kmul(K22L,ToReal(-2)),kmul(kmadd(g12L,JacPDstandardNth2beta1,kmadd(g22L,JacPDstandardNth2beta2,kmul(g23L,JacPDstandardNth2beta3))),ToReal(2))))));
    
    CCTK_REAL_VEC g23rhsL = 
      kmadd(beta1L,JacPDstandardNth1g23,kmadd(g13L,JacPDstandardNth2beta1,kmadd(g33L,JacPDstandardNth2beta3,kmadd(beta2L,JacPDstandardNth2g23,kmadd(g12L,JacPDstandardNth3beta1,kmadd(g22L,JacPDstandardNth3beta2,kmadd(g23L,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3),kmadd(beta3L,JacPDstandardNth3g23,kmul(alphaL,kmul(K23L,ToReal(-2)))))))))));
    
    CCTK_REAL_VEC g33rhsL = 
      kmadd(beta1L,JacPDstandardNth1g33,kmadd(beta2L,JacPDstandardNth2g33,kmadd(beta3L,JacPDstandardNth3g33,kmadd(alphaL,kmul(K33L,ToReal(-2)),kmul(kmadd(g13L,JacPDstandardNth3beta1,kmadd(g23L,JacPDstandardNth3beta2,kmul(g33L,JacPDstandardNth3beta3))),ToReal(2))))));
    
    CCTK_REAL_VEC K11rhsL = 
      kmadd(G111,JacPDstandardNth1alpha,kmadd(beta1L,JacPDstandardNth1K11,kmadd(G211,JacPDstandardNth2alpha,kmadd(beta2L,JacPDstandardNth2K11,kmadd(G311,JacPDstandardNth3alpha,kmadd(beta3L,JacPDstandardNth3K11,kmadd(alphaL,kadd(R11,kmadd(kmadd(K12L,Km21,kmul(K13L,Km31)),ToReal(-2),kmul(K11L,kmadd(Km11,ToReal(-2),trK)))),kmsub(kmadd(JacPDstandardNth1beta1,K11L,kmadd(JacPDstandardNth1beta2,K12L,kmul(JacPDstandardNth1beta3,K13L))),ToReal(2),JacPDstandardNth11alpha))))))));
    
    CCTK_REAL_VEC K12rhsL = 
      kmadd(G112,JacPDstandardNth1alpha,kmadd(beta1L,JacPDstandardNth1K12,kmadd(G212,JacPDstandardNth2alpha,kmadd(beta2L,JacPDstandardNth2K12,kmadd(G312,JacPDstandardNth3alpha,kmadd(beta3L,JacPDstandardNth3K12,kmadd(JacPDstandardNth2beta1,K11L,kmadd(kadd(JacPDstandardNth1beta1,JacPDstandardNth2beta2),K12L,kmadd(JacPDstandardNth2beta3,K13L,kmadd(JacPDstandardNth1beta2,K22L,kmadd(JacPDstandardNth1beta3,K23L,kmsub(alphaL,kadd(R12,kmadd(kmadd(K11L,Km12,kmul(K13L,Km32)),ToReal(-2),kmul(K12L,kmadd(Km22,ToReal(-2),trK)))),JacPDstandardNth12alpha))))))))))));
    
    CCTK_REAL_VEC K13rhsL = 
      kmadd(G113,JacPDstandardNth1alpha,kmadd(beta1L,JacPDstandardNth1K13,kmadd(G213,JacPDstandardNth2alpha,kmadd(beta2L,JacPDstandardNth2K13,kmadd(G313,JacPDstandardNth3alpha,kmadd(beta3L,JacPDstandardNth3K13,kmadd(JacPDstandardNth3beta1,K11L,kmadd(JacPDstandardNth3beta2,K12L,kmadd(kadd(JacPDstandardNth1beta1,JacPDstandardNth3beta3),K13L,kmadd(JacPDstandardNth1beta2,K23L,kmadd(JacPDstandardNth1beta3,K33L,kmsub(alphaL,kadd(R13,kmadd(K13L,trK,kmul(kmadd(K11L,Km13,kmadd(K12L,Km23,kmul(K13L,Km33))),ToReal(-2)))),JacPDstandardNth13alpha))))))))))));
    
    CCTK_REAL_VEC K22rhsL = 
      kmadd(G122,JacPDstandardNth1alpha,kmadd(beta1L,JacPDstandardNth1K22,kmadd(G222,JacPDstandardNth2alpha,kmadd(beta2L,JacPDstandardNth2K22,kmadd(G322,JacPDstandardNth3alpha,kmadd(beta3L,JacPDstandardNth3K22,kmadd(alphaL,kadd(R22,kmadd(kmadd(K12L,Km12,kmul(K23L,Km32)),ToReal(-2),kmul(K22L,kmadd(Km22,ToReal(-2),trK)))),kmsub(kmadd(JacPDstandardNth2beta1,K12L,kmadd(JacPDstandardNth2beta2,K22L,kmul(JacPDstandardNth2beta3,K23L))),ToReal(2),JacPDstandardNth22alpha))))))));
    
    CCTK_REAL_VEC K23rhsL = 
      kmadd(G123,JacPDstandardNth1alpha,kmadd(beta1L,JacPDstandardNth1K23,kmadd(G223,JacPDstandardNth2alpha,kmadd(beta2L,JacPDstandardNth2K23,kmadd(G323,JacPDstandardNth3alpha,kmadd(beta3L,JacPDstandardNth3K23,kmadd(JacPDstandardNth3beta1,K12L,kmadd(JacPDstandardNth2beta1,K13L,kmadd(JacPDstandardNth3beta2,K22L,kmadd(kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3),K23L,kmadd(JacPDstandardNth2beta3,K33L,kmsub(alphaL,kadd(R23,kmadd(K23L,trK,kmul(kmadd(K12L,Km13,kmadd(K22L,Km23,kmul(K23L,Km33))),ToReal(-2)))),JacPDstandardNth23alpha))))))))))));
    
    CCTK_REAL_VEC K33rhsL = 
      kmadd(G133,JacPDstandardNth1alpha,kmadd(beta1L,JacPDstandardNth1K33,kmadd(G233,JacPDstandardNth2alpha,kmadd(beta2L,JacPDstandardNth2K33,kmadd(G333,JacPDstandardNth3alpha,kmadd(beta3L,JacPDstandardNth3K33,kmadd(alphaL,kadd(R33,kmadd(K33L,trK,kmul(kmadd(K13L,Km13,kmadd(K23L,Km23,kmul(K33L,Km33))),ToReal(-2)))),kmsub(kmadd(JacPDstandardNth3beta1,K13L,kmadd(JacPDstandardNth3beta2,K23L,kmul(JacPDstandardNth3beta3,K33L))),ToReal(2),JacPDstandardNth33alpha))))))));
    
    CCTK_REAL_VEC alpharhsL = ToReal(0);
    
    CCTK_REAL_VEC beta1rhsL = ToReal(0);
    
    CCTK_REAL_VEC beta2rhsL = ToReal(0);
    
    CCTK_REAL_VEC beta3rhsL = ToReal(0);
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(alpharhs[index],alpharhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta1rhs[index],beta1rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta2rhs[index],beta2rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta3rhs[index],beta3rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g11rhs[index],g11rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g12rhs[index],g12rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g13rhs[index],g13rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g22rhs[index],g22rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g23rhs[index],g23rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g33rhs[index],g33rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K11rhs[index],K11rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K12rhs[index],K12rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K13rhs[index],K13rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K22rhs[index],K22rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K23rhs[index],K23rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K33rhs[index],K33rhsL,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(alpharhs[index],alpharhsL,elt_count);
      vec_store_nta_partial_hi(beta1rhs[index],beta1rhsL,elt_count);
      vec_store_nta_partial_hi(beta2rhs[index],beta2rhsL,elt_count);
      vec_store_nta_partial_hi(beta3rhs[index],beta3rhsL,elt_count);
      vec_store_nta_partial_hi(g11rhs[index],g11rhsL,elt_count);
      vec_store_nta_partial_hi(g12rhs[index],g12rhsL,elt_count);
      vec_store_nta_partial_hi(g13rhs[index],g13rhsL,elt_count);
      vec_store_nta_partial_hi(g22rhs[index],g22rhsL,elt_count);
      vec_store_nta_partial_hi(g23rhs[index],g23rhsL,elt_count);
      vec_store_nta_partial_hi(g33rhs[index],g33rhsL,elt_count);
      vec_store_nta_partial_hi(K11rhs[index],K11rhsL,elt_count);
      vec_store_nta_partial_hi(K12rhs[index],K12rhsL,elt_count);
      vec_store_nta_partial_hi(K13rhs[index],K13rhsL,elt_count);
      vec_store_nta_partial_hi(K22rhs[index],K22rhsL,elt_count);
      vec_store_nta_partial_hi(K23rhs[index],K23rhsL,elt_count);
      vec_store_nta_partial_hi(K33rhs[index],K33rhsL,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(alpharhs[index],alpharhsL,elt_count);
      vec_store_nta_partial_lo(beta1rhs[index],beta1rhsL,elt_count);
      vec_store_nta_partial_lo(beta2rhs[index],beta2rhsL,elt_count);
      vec_store_nta_partial_lo(beta3rhs[index],beta3rhsL,elt_count);
      vec_store_nta_partial_lo(g11rhs[index],g11rhsL,elt_count);
      vec_store_nta_partial_lo(g12rhs[index],g12rhsL,elt_count);
      vec_store_nta_partial_lo(g13rhs[index],g13rhsL,elt_count);
      vec_store_nta_partial_lo(g22rhs[index],g22rhsL,elt_count);
      vec_store_nta_partial_lo(g23rhs[index],g23rhsL,elt_count);
      vec_store_nta_partial_lo(g33rhs[index],g33rhsL,elt_count);
      vec_store_nta_partial_lo(K11rhs[index],K11rhsL,elt_count);
      vec_store_nta_partial_lo(K12rhs[index],K12rhsL,elt_count);
      vec_store_nta_partial_lo(K13rhs[index],K13rhsL,elt_count);
      vec_store_nta_partial_lo(K22rhs[index],K22rhsL,elt_count);
      vec_store_nta_partial_lo(K23rhs[index],K23rhsL,elt_count);
      vec_store_nta_partial_lo(K33rhs[index],K33rhsL,elt_count);
      break;
    }
    vec_store_nta(alpharhs[index],alpharhsL);
    vec_store_nta(beta1rhs[index],beta1rhsL);
    vec_store_nta(beta2rhs[index],beta2rhsL);
    vec_store_nta(beta3rhs[index],beta3rhsL);
    vec_store_nta(g11rhs[index],g11rhsL);
    vec_store_nta(g12rhs[index],g12rhsL);
    vec_store_nta(g13rhs[index],g13rhsL);
    vec_store_nta(g22rhs[index],g22rhsL);
    vec_store_nta(g23rhs[index],g23rhsL);
    vec_store_nta(g33rhs[index],g33rhsL);
    vec_store_nta(K11rhs[index],K11rhsL);
    vec_store_nta(K12rhs[index],K12rhsL);
    vec_store_nta(K13rhs[index],K13rhsL);
    vec_store_nta(K22rhs[index],K22rhsL);
    vec_store_nta(K23rhs[index],K23rhsL);
    vec_store_nta(K33rhs[index],K33rhsL);
  }
  LC_ENDLOOP3VEC (ML_ADM_RHS);
}

extern "C" void ML_ADM_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADM_RHS_Body");
  }
  
  if (cctk_iteration % ML_ADM_RHS_calc_every != ML_ADM_RHS_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_ADM::ML_curv","ML_ADM::ML_curvrhs","ML_ADM::ML_lapse","ML_ADM::ML_lapserhs","ML_ADM::ML_metric","ML_ADM::ML_metricrhs","ML_ADM::ML_shift","ML_ADM::ML_shiftrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADM_RHS", 8, groups);
  
  switch(fdOrder)
  {
    case 2:
      GenericFD_EnsureStencilFits(cctkGH, "ML_ADM_RHS", 1, 1, 1);
      break;
    
    case 4:
      GenericFD_EnsureStencilFits(cctkGH, "ML_ADM_RHS", 2, 2, 2);
      break;
    
    case 6:
      GenericFD_EnsureStencilFits(cctkGH, "ML_ADM_RHS", 3, 3, 3);
      break;
    
    case 8:
      GenericFD_EnsureStencilFits(cctkGH, "ML_ADM_RHS", 4, 4, 4);
      break;
  }
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADM_RHS_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_ADM_RHS_Body");
  }
}
