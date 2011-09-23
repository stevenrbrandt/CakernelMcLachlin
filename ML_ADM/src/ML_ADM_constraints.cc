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

extern "C" void ML_ADM_constraints_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_Ham","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_Ham.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_mom","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_mom.");
  return;
}

static void ML_ADM_constraints_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADM_constraints_Body");
  }
  
  if (cctk_iteration % ML_ADM_constraints_calc_every != ML_ADM_constraints_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_ADM::ML_curv","ML_ADM::ML_Ham","ML_ADM::ML_metric","ML_ADM::ML_mom"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADM_constraints", 4, groups);
  
  switch(fdOrder)
  {
    case 2:
      GenericFD_EnsureStencilFits(cctkGH, "ML_ADM_constraints", 1, 1, 1);
      break;
    
    case 4:
      GenericFD_EnsureStencilFits(cctkGH, "ML_ADM_constraints", 2, 2, 2);
      break;
    
    case 6:
      GenericFD_EnsureStencilFits(cctkGH, "ML_ADM_constraints", 3, 3, 3);
      break;
    
    case 8:
      GenericFD_EnsureStencilFits(cctkGH, "ML_ADM_constraints", 4, 4, 4);
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
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC (ML_ADM_constraints,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
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
    CCTK_REAL_VEC JacPDstandardNth11g22;
    CCTK_REAL_VEC JacPDstandardNth11g23;
    CCTK_REAL_VEC JacPDstandardNth11g33;
    CCTK_REAL_VEC JacPDstandardNth12g11;
    CCTK_REAL_VEC JacPDstandardNth12g12;
    CCTK_REAL_VEC JacPDstandardNth12g13;
    CCTK_REAL_VEC JacPDstandardNth12g22;
    CCTK_REAL_VEC JacPDstandardNth12g23;
    CCTK_REAL_VEC JacPDstandardNth12g33;
    CCTK_REAL_VEC JacPDstandardNth13g11;
    CCTK_REAL_VEC JacPDstandardNth13g12;
    CCTK_REAL_VEC JacPDstandardNth13g13;
    CCTK_REAL_VEC JacPDstandardNth13g22;
    CCTK_REAL_VEC JacPDstandardNth13g23;
    CCTK_REAL_VEC JacPDstandardNth13g33;
    CCTK_REAL_VEC JacPDstandardNth1g11;
    CCTK_REAL_VEC JacPDstandardNth1g12;
    CCTK_REAL_VEC JacPDstandardNth1g13;
    CCTK_REAL_VEC JacPDstandardNth1g22;
    CCTK_REAL_VEC JacPDstandardNth1g23;
    CCTK_REAL_VEC JacPDstandardNth1g33;
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
    CCTK_REAL_VEC JacPDstandardNth22g11;
    CCTK_REAL_VEC JacPDstandardNth22g13;
    CCTK_REAL_VEC JacPDstandardNth22g33;
    CCTK_REAL_VEC JacPDstandardNth23g11;
    CCTK_REAL_VEC JacPDstandardNth23g12;
    CCTK_REAL_VEC JacPDstandardNth23g13;
    CCTK_REAL_VEC JacPDstandardNth23g22;
    CCTK_REAL_VEC JacPDstandardNth23g23;
    CCTK_REAL_VEC JacPDstandardNth23g33;
    CCTK_REAL_VEC JacPDstandardNth2g11;
    CCTK_REAL_VEC JacPDstandardNth2g12;
    CCTK_REAL_VEC JacPDstandardNth2g13;
    CCTK_REAL_VEC JacPDstandardNth2g22;
    CCTK_REAL_VEC JacPDstandardNth2g23;
    CCTK_REAL_VEC JacPDstandardNth2g33;
    CCTK_REAL_VEC JacPDstandardNth2K11;
    CCTK_REAL_VEC JacPDstandardNth2K12;
    CCTK_REAL_VEC JacPDstandardNth2K13;
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
    CCTK_REAL_VEC JacPDstandardNth33g11;
    CCTK_REAL_VEC JacPDstandardNth33g12;
    CCTK_REAL_VEC JacPDstandardNth33g22;
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
    
    if (use_jacobian)
    {
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
      
      JacPDstandardNth2K23 = 
        kmadd(J12L,PDstandardNth1K23,kmadd(J22L,PDstandardNth2K23,kmul(J32L,PDstandardNth3K23)));
      
      JacPDstandardNth2K33 = 
        kmadd(J12L,PDstandardNth1K33,kmadd(J22L,PDstandardNth2K33,kmul(J32L,PDstandardNth3K33)));
      
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
      
      JacPDstandardNth11g22 = 
        kmadd(dJ111L,PDstandardNth1g22,kmadd(dJ211L,PDstandardNth2g22,kmadd(dJ311L,PDstandardNth3g22,kmadd(PDstandardNth11g22,SQR(J11L),kmadd(PDstandardNth22g22,SQR(J21L),kmadd(PDstandardNth33g22,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12g22,kmul(J31L,PDstandardNth13g22)),kmul(J21L,kmul(J31L,PDstandardNth23g22))),ToReal(2))))))));
      
      JacPDstandardNth11g23 = 
        kmadd(dJ111L,PDstandardNth1g23,kmadd(dJ211L,PDstandardNth2g23,kmadd(dJ311L,PDstandardNth3g23,kmadd(PDstandardNth11g23,SQR(J11L),kmadd(PDstandardNth22g23,SQR(J21L),kmadd(PDstandardNth33g23,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12g23,kmul(J31L,PDstandardNth13g23)),kmul(J21L,kmul(J31L,PDstandardNth23g23))),ToReal(2))))))));
      
      JacPDstandardNth11g33 = 
        kmadd(dJ111L,PDstandardNth1g33,kmadd(dJ211L,PDstandardNth2g33,kmadd(dJ311L,PDstandardNth3g33,kmadd(PDstandardNth11g33,SQR(J11L),kmadd(PDstandardNth22g33,SQR(J21L),kmadd(PDstandardNth33g33,SQR(J31L),kmul(kmadd(J11L,kmadd(J21L,PDstandardNth12g33,kmul(J31L,PDstandardNth13g33)),kmul(J21L,kmul(J31L,PDstandardNth23g33))),ToReal(2))))))));
      
      JacPDstandardNth22g11 = 
        kmadd(dJ122L,PDstandardNth1g11,kmadd(dJ222L,PDstandardNth2g11,kmadd(dJ322L,PDstandardNth3g11,kmadd(PDstandardNth11g11,SQR(J12L),kmadd(PDstandardNth22g11,SQR(J22L),kmadd(PDstandardNth33g11,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12g11,kmul(J32L,PDstandardNth13g11)),kmul(J22L,kmul(J32L,PDstandardNth23g11))),ToReal(2))))))));
      
      JacPDstandardNth22g13 = 
        kmadd(dJ122L,PDstandardNth1g13,kmadd(dJ222L,PDstandardNth2g13,kmadd(dJ322L,PDstandardNth3g13,kmadd(PDstandardNth11g13,SQR(J12L),kmadd(PDstandardNth22g13,SQR(J22L),kmadd(PDstandardNth33g13,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12g13,kmul(J32L,PDstandardNth13g13)),kmul(J22L,kmul(J32L,PDstandardNth23g13))),ToReal(2))))))));
      
      JacPDstandardNth22g33 = 
        kmadd(dJ122L,PDstandardNth1g33,kmadd(dJ222L,PDstandardNth2g33,kmadd(dJ322L,PDstandardNth3g33,kmadd(PDstandardNth11g33,SQR(J12L),kmadd(PDstandardNth22g33,SQR(J22L),kmadd(PDstandardNth33g33,SQR(J32L),kmul(kmadd(J12L,kmadd(J22L,PDstandardNth12g33,kmul(J32L,PDstandardNth13g33)),kmul(J22L,kmul(J32L,PDstandardNth23g33))),ToReal(2))))))));
      
      JacPDstandardNth33g11 = 
        kmadd(dJ133L,PDstandardNth1g11,kmadd(dJ233L,PDstandardNth2g11,kmadd(dJ333L,PDstandardNth3g11,kmadd(PDstandardNth11g11,SQR(J13L),kmadd(PDstandardNth22g11,SQR(J23L),kmadd(PDstandardNth33g11,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12g11,kmul(J33L,PDstandardNth13g11)),kmul(J23L,kmul(J33L,PDstandardNth23g11))),ToReal(2))))))));
      
      JacPDstandardNth33g12 = 
        kmadd(dJ133L,PDstandardNth1g12,kmadd(dJ233L,PDstandardNth2g12,kmadd(dJ333L,PDstandardNth3g12,kmadd(PDstandardNth11g12,SQR(J13L),kmadd(PDstandardNth22g12,SQR(J23L),kmadd(PDstandardNth33g12,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12g12,kmul(J33L,PDstandardNth13g12)),kmul(J23L,kmul(J33L,PDstandardNth23g12))),ToReal(2))))))));
      
      JacPDstandardNth33g22 = 
        kmadd(dJ133L,PDstandardNth1g22,kmadd(dJ233L,PDstandardNth2g22,kmadd(dJ333L,PDstandardNth3g22,kmadd(PDstandardNth11g22,SQR(J13L),kmadd(PDstandardNth22g22,SQR(J23L),kmadd(PDstandardNth33g22,SQR(J33L),kmul(kmadd(J13L,kmadd(J23L,PDstandardNth12g22,kmul(J33L,PDstandardNth13g22)),kmul(J23L,kmul(J33L,PDstandardNth23g22))),ToReal(2))))))));
      
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
      JacPDstandardNth1g11 = PDstandardNth1g11;
      
      JacPDstandardNth1g12 = PDstandardNth1g12;
      
      JacPDstandardNth1g13 = PDstandardNth1g13;
      
      JacPDstandardNth1g22 = PDstandardNth1g22;
      
      JacPDstandardNth1g23 = PDstandardNth1g23;
      
      JacPDstandardNth1g33 = PDstandardNth1g33;
      
      JacPDstandardNth1K12 = PDstandardNth1K12;
      
      JacPDstandardNth1K13 = PDstandardNth1K13;
      
      JacPDstandardNth1K22 = PDstandardNth1K22;
      
      JacPDstandardNth1K23 = PDstandardNth1K23;
      
      JacPDstandardNth1K33 = PDstandardNth1K33;
      
      JacPDstandardNth2g11 = PDstandardNth2g11;
      
      JacPDstandardNth2g12 = PDstandardNth2g12;
      
      JacPDstandardNth2g13 = PDstandardNth2g13;
      
      JacPDstandardNth2g22 = PDstandardNth2g22;
      
      JacPDstandardNth2g23 = PDstandardNth2g23;
      
      JacPDstandardNth2g33 = PDstandardNth2g33;
      
      JacPDstandardNth2K11 = PDstandardNth2K11;
      
      JacPDstandardNth2K12 = PDstandardNth2K12;
      
      JacPDstandardNth2K13 = PDstandardNth2K13;
      
      JacPDstandardNth2K23 = PDstandardNth2K23;
      
      JacPDstandardNth2K33 = PDstandardNth2K33;
      
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
      
      JacPDstandardNth11g22 = PDstandardNth11g22;
      
      JacPDstandardNth11g23 = PDstandardNth11g23;
      
      JacPDstandardNth11g33 = PDstandardNth11g33;
      
      JacPDstandardNth22g11 = PDstandardNth22g11;
      
      JacPDstandardNth22g13 = PDstandardNth22g13;
      
      JacPDstandardNth22g33 = PDstandardNth22g33;
      
      JacPDstandardNth33g11 = PDstandardNth33g11;
      
      JacPDstandardNth33g12 = PDstandardNth33g12;
      
      JacPDstandardNth33g22 = PDstandardNth33g22;
      
      JacPDstandardNth12g11 = PDstandardNth12g11;
      
      JacPDstandardNth12g12 = PDstandardNth12g12;
      
      JacPDstandardNth12g13 = PDstandardNth12g13;
      
      JacPDstandardNth12g22 = PDstandardNth12g22;
      
      JacPDstandardNth12g23 = PDstandardNth12g23;
      
      JacPDstandardNth12g33 = PDstandardNth12g33;
      
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
      kadd(SQR(G111),kadd(SQR(G112),kadd(SQR(G113),kadd(SQR(G211),kadd(SQR(G212),kadd(SQR(G213),kadd(SQR(G311),kadd(SQR(G312),kadd(SQR(G313),knmsub(G111,kadd(G111,kadd(G122,G133)),knmsub(G211,kadd(G211,kadd(G222,G233)),knmsub(G311,kadd(G311,kadd(G322,G333)),kmadd(gu22,kmul(kadd(JacPDstandardNth11g22,kmadd(JacPDstandardNth21g12,ToReal(-2),JacPDstandardNth22g11)),ToReal(-0.5)),kmadd(gu33,kmul(kadd(JacPDstandardNth11g33,kmadd(JacPDstandardNth31g13,ToReal(-2),JacPDstandardNth33g11)),ToReal(-0.5)),kmadd(gu12,kmul(ksub(JacPDstandardNth21g11,JacPDstandardNth12g11),ToReal(0.5)),kmadd(gu13,kmul(ksub(JacPDstandardNth31g11,JacPDstandardNth13g11),ToReal(0.5)),kmadd(gu23,kmul(kadd(JacPDstandardNth21g13,ksub(JacPDstandardNth31g12,kadd(JacPDstandardNth23g11,JacPDstandardNth11g23))),ToReal(0.5)),kmul(gu32,kmul(kadd(JacPDstandardNth21g13,ksub(JacPDstandardNth31g12,kadd(JacPDstandardNth32g11,JacPDstandardNth11g23))),ToReal(0.5))))))))))))))))))));
    
    CCTK_REAL_VEC R12 = 
      kmul(ToReal(0.5),kmadd(gu12,kadd(JacPDstandardNth11g22,JacPDstandardNth22g11),kmadd(gu32,JacPDstandardNth31g22,kmadd(gu13,JacPDstandardNth32g11,kmadd(gu23,JacPDstandardNth32g12,kmadd(gu33,JacPDstandardNth32g13,kmadd(kmadd(G112,G133,kmadd(G212,G233,kmadd(G312,G333,kmul(gu12,JacPDstandardNth12g12)))),ToReal(-2),knmsub(JacPDstandardNth12g23,kadd(gu32,gu23),kmadd(gu22,ksub(JacPDstandardNth21g22,JacPDstandardNth12g22),kmadd(gu13,ksub(JacPDstandardNth11g23,kadd(JacPDstandardNth13g12,JacPDstandardNth12g13)),kmadd(gu23,ksub(JacPDstandardNth21g23,JacPDstandardNth23g12),kmadd(gu32,ksub(JacPDstandardNth22g13,JacPDstandardNth32g12),kmadd(gu33,ksub(JacPDstandardNth31g23,kadd(JacPDstandardNth33g12,JacPDstandardNth12g33)),kmul(kmadd(G113,G123,kmadd(G213,G223,kmul(G313,G323))),ToReal(2)))))))))))))));
    
    CCTK_REAL_VEC R13 = 
      kmul(ToReal(0.5),kmadd(gu22,JacPDstandardNth23g12,kmadd(gu32,JacPDstandardNth31g23,kmadd(gu13,kadd(JacPDstandardNth11g33,JacPDstandardNth33g11),kmadd(gu23,JacPDstandardNth33g12,kmadd(kmadd(G113,G122,kmadd(G213,G222,kmadd(G313,G322,kmul(gu13,JacPDstandardNth13g13)))),ToReal(-2),knmsub(JacPDstandardNth13g23,kadd(gu32,gu23),kmadd(gu12,kadd(JacPDstandardNth11g23,ksub(JacPDstandardNth23g11,kadd(JacPDstandardNth13g12,JacPDstandardNth12g13))),kmadd(gu33,ksub(JacPDstandardNth31g33,JacPDstandardNth13g33),kmadd(gu22,ksub(JacPDstandardNth21g23,kadd(JacPDstandardNth22g13,JacPDstandardNth13g22)),kmadd(gu23,ksub(JacPDstandardNth21g33,JacPDstandardNth23g13),kmadd(gu32,ksub(JacPDstandardNth23g13,JacPDstandardNth32g13),kmul(kmadd(G112,G123,kmadd(G212,G223,kmul(G312,G323))),ToReal(2))))))))))))));
    
    CCTK_REAL_VEC R22 = 
      kadd(SQR(G112),kadd(SQR(G122),kadd(SQR(G123),kadd(SQR(G212),kadd(SQR(G222),kadd(SQR(G223),kadd(SQR(G312),kadd(SQR(G322),kadd(SQR(G323),knmsub(G122,kadd(G111,kadd(G122,G133)),knmsub(G222,kadd(G211,kadd(G222,G233)),knmsub(G322,kadd(G311,kadd(G322,G333)),kmadd(gu11,kmul(kadd(JacPDstandardNth11g22,kmadd(JacPDstandardNth12g12,ToReal(-2),JacPDstandardNth22g11)),ToReal(-0.5)),kmadd(gu33,kmul(kadd(JacPDstandardNth22g33,kmadd(JacPDstandardNth32g23,ToReal(-2),JacPDstandardNth33g22)),ToReal(-0.5)),kmadd(gu21,kmul(ksub(JacPDstandardNth12g22,JacPDstandardNth21g22),ToReal(0.5)),kmadd(gu13,kmul(kadd(JacPDstandardNth12g23,ksub(JacPDstandardNth32g12,kadd(JacPDstandardNth22g13,JacPDstandardNth13g22))),ToReal(0.5)),kmadd(gu23,kmul(ksub(JacPDstandardNth32g22,JacPDstandardNth23g22),ToReal(0.5)),kmul(gu31,kmul(kadd(JacPDstandardNth12g23,ksub(JacPDstandardNth32g12,kadd(JacPDstandardNth31g22,JacPDstandardNth22g13))),ToReal(0.5))))))))))))))))))));
    
    CCTK_REAL_VEC R23 = 
      kmul(ToReal(0.5),kmadd(gu23,kadd(JacPDstandardNth22g33,JacPDstandardNth33g22),kmadd(kmadd(G111,G123,kmadd(G211,G223,kmadd(G311,G323,kmul(gu23,JacPDstandardNth23g23)))),ToReal(-2),kmadd(gu11,kadd(JacPDstandardNth12g13,ksub(JacPDstandardNth13g12,kadd(JacPDstandardNth23g11,JacPDstandardNth11g23))),kmadd(gu21,kadd(JacPDstandardNth13g22,ksub(JacPDstandardNth22g13,kadd(JacPDstandardNth23g12,JacPDstandardNth21g23))),kmadd(gu13,kadd(JacPDstandardNth12g33,ksub(JacPDstandardNth33g12,kadd(JacPDstandardNth23g13,JacPDstandardNth13g23))),kmadd(gu33,ksub(JacPDstandardNth32g33,JacPDstandardNth23g33),kmadd(gu31,kadd(JacPDstandardNth13g23,ksub(JacPDstandardNth32g13,kadd(JacPDstandardNth31g23,JacPDstandardNth23g13))),kmul(kmadd(G112,G113,kmadd(G212,G213,kmul(G312,G313))),ToReal(2))))))))));
    
    CCTK_REAL_VEC R33 = 
      kadd(SQR(G113),kadd(SQR(G123),kadd(SQR(G133),kadd(SQR(G213),kadd(SQR(G223),kadd(SQR(G233),kadd(SQR(G313),kadd(SQR(G323),kadd(SQR(G333),knmsub(G133,kadd(G111,kadd(G122,G133)),knmsub(G233,kadd(G211,kadd(G222,G233)),knmsub(G333,kadd(G311,kadd(G322,G333)),kmadd(gu11,kmul(kadd(JacPDstandardNth11g33,kmadd(JacPDstandardNth13g13,ToReal(-2),JacPDstandardNth33g11)),ToReal(-0.5)),kmadd(gu22,kmul(kadd(JacPDstandardNth22g33,kmadd(JacPDstandardNth23g23,ToReal(-2),JacPDstandardNth33g22)),ToReal(-0.5)),kmadd(gu31,kmul(ksub(JacPDstandardNth13g33,JacPDstandardNth31g33),ToReal(0.5)),kmadd(gu32,kmul(ksub(JacPDstandardNth23g33,JacPDstandardNth32g33),ToReal(0.5)),kmadd(gu12,kmul(kadd(JacPDstandardNth13g23,ksub(JacPDstandardNth23g13,kadd(JacPDstandardNth33g12,JacPDstandardNth12g33))),ToReal(0.5)),kmul(gu21,kmul(kadd(JacPDstandardNth13g23,ksub(JacPDstandardNth23g13,kadd(JacPDstandardNth33g12,JacPDstandardNth21g33))),ToReal(0.5))))))))))))))))))));
    
    CCTK_REAL_VEC trR = 
      kmadd(gu11,R11,kmadd(kadd(gu12,gu21),R12,kmadd(kadd(gu13,gu31),R13,kmadd(gu22,R22,kmadd(kadd(gu23,gu32),R23,kmul(gu33,R33))))));
    
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
    
    CCTK_REAL_VEC HL = 
      kadd(trR,kadd(SQR(trK),kmsub(kmadd(Km12,Km21,kmadd(Km13,Km31,kmul(Km23,Km32))),ToReal(-2),kadd(kadd(SQR(Km33),SQR(Km22)),SQR(Km11)))));
    
    CCTK_REAL_VEC M1L = 
      kmadd(gu21,kadd(JacPDstandardNth2K11,kmadd(G211,K22L,kmadd(G311,K23L,ksub(knmsub(G112,K11L,kmsub(K12L,ksub(G111,G212),kmul(G312,K13L))),JacPDstandardNth1K12)))),kmadd(gu22,kadd(JacPDstandardNth2K12,kmadd(G212,K22L,kmadd(G312,K23L,ksub(knmsub(G122,K11L,kmsub(K12L,ksub(G112,G222),kmul(G322,K13L))),JacPDstandardNth1K22)))),kmadd(gu23,kadd(JacPDstandardNth2K13,kmadd(G213,K22L,kmadd(G313,K23L,ksub(knmsub(G123,K11L,kmsub(K12L,ksub(G113,G223),kmul(G323,K13L))),JacPDstandardNth1K23)))),kmadd(gu31,kadd(JacPDstandardNth3K11,kmadd(G211,K23L,kmadd(G311,K33L,ksub(knmsub(G113,K11L,kmsub(K13L,ksub(G111,G313),kmul(G213,K12L))),JacPDstandardNth1K13)))),kmadd(gu32,kadd(JacPDstandardNth3K12,kmadd(G212,K23L,kmadd(G312,K33L,ksub(knmsub(G123,K11L,kmsub(K13L,ksub(G112,G323),kmul(G223,K12L))),JacPDstandardNth1K23)))),kmul(gu33,kadd(JacPDstandardNth3K13,kmadd(G213,K23L,kmadd(G313,K33L,ksub(knmsub(G133,K11L,kmsub(K13L,ksub(G113,G333),kmul(G233,K12L))),JacPDstandardNth1K33))))))))));
    
    CCTK_REAL_VEC M2L = 
      kmadd(gu11,kadd(JacPDstandardNth1K12,kmadd(G112,K11L,kmadd(G312,K13L,ksub(knmsub(G211,K22L,kmsub(K12L,ksub(G212,G111),kmul(G311,K23L))),JacPDstandardNth2K11)))),kmadd(gu12,kadd(JacPDstandardNth1K22,kmadd(G122,K11L,kmadd(G322,K13L,ksub(knmsub(G212,K22L,kmsub(K12L,ksub(G222,G112),kmul(G312,K23L))),JacPDstandardNth2K12)))),kmadd(gu13,kadd(JacPDstandardNth1K23,kmadd(G123,K11L,kmadd(G323,K13L,ksub(knmsub(G213,K22L,kmsub(K12L,ksub(G223,G113),kmul(G313,K23L))),JacPDstandardNth2K13)))),kmadd(gu31,kadd(JacPDstandardNth3K12,kmadd(G112,K13L,kmadd(G312,K33L,ksub(knmsub(G113,K12L,kmsub(K23L,ksub(G212,G313),kmul(G213,K22L))),JacPDstandardNth2K13)))),kmadd(gu32,kadd(JacPDstandardNth3K22,kmadd(G122,K13L,kmadd(G322,K33L,ksub(knmsub(G123,K12L,kmsub(K23L,ksub(G222,G323),kmul(G223,K22L))),JacPDstandardNth2K23)))),kmul(gu33,kadd(JacPDstandardNth3K23,kmadd(G123,K13L,kmadd(G323,K33L,ksub(knmsub(G133,K12L,kmsub(K23L,ksub(G223,G333),kmul(G233,K22L))),JacPDstandardNth2K33))))))))));
    
    CCTK_REAL_VEC M3L = 
      kmadd(gu11,kadd(JacPDstandardNth1K13,kmadd(G113,K11L,kmadd(G213,K12L,ksub(knmsub(G211,K23L,kmsub(K13L,ksub(G313,G111),kmul(G311,K33L))),JacPDstandardNth3K11)))),kmadd(gu12,kadd(JacPDstandardNth1K23,kmadd(G123,K11L,kmadd(G223,K12L,ksub(knmsub(G212,K23L,kmsub(K13L,ksub(G323,G112),kmul(G312,K33L))),JacPDstandardNth3K12)))),kmadd(gu13,kadd(JacPDstandardNth1K33,kmadd(G133,K11L,kmadd(G233,K12L,ksub(knmsub(G213,K23L,kmsub(K13L,ksub(G333,G113),kmul(G313,K33L))),JacPDstandardNth3K13)))),kmadd(gu21,kadd(JacPDstandardNth2K13,kmadd(G113,K12L,kmadd(G213,K22L,ksub(knmsub(G112,K13L,kmsub(K23L,ksub(G313,G212),kmul(G312,K33L))),JacPDstandardNth3K12)))),kmadd(gu22,kadd(JacPDstandardNth2K23,kmadd(G123,K12L,kmadd(G223,K22L,ksub(knmsub(G122,K13L,kmsub(K23L,ksub(G323,G222),kmul(G322,K33L))),JacPDstandardNth3K22)))),kmul(gu23,kadd(JacPDstandardNth2K33,kmadd(G133,K12L,kmadd(G233,K22L,ksub(knmsub(G123,K13L,kmsub(K23L,ksub(G333,G223),kmul(G323,K33L))),JacPDstandardNth3K23))))))))));
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(H[index],HL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(M1[index],M1L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(M2[index],M2L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(M3[index],M3L,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(H[index],HL,elt_count);
      vec_store_nta_partial_hi(M1[index],M1L,elt_count);
      vec_store_nta_partial_hi(M2[index],M2L,elt_count);
      vec_store_nta_partial_hi(M3[index],M3L,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(H[index],HL,elt_count);
      vec_store_nta_partial_lo(M1[index],M1L,elt_count);
      vec_store_nta_partial_lo(M2[index],M2L,elt_count);
      vec_store_nta_partial_lo(M3[index],M3L,elt_count);
      break;
    }
    
    /* Copy local copies back to grid functions */
    vec_store_nta(H[index],HL);
    vec_store_nta(M1[index],M1L);
    vec_store_nta(M2[index],M2L);
    vec_store_nta(M3[index],M3L);
  }
  LC_ENDLOOP3VEC (ML_ADM_constraints);
}

extern "C" void ML_ADM_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADM_constraints_Body);
}
