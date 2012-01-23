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

extern "C" void ML_ADM_boundary_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_curv","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_curv.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_lapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_lapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_metric","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_metric.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_shift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_shift.");
  return;
}

static void ML_ADM_boundary_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  LC_LOOP3VEC (ML_ADM_boundary,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    switch(fdOrder)
    {
      case 2:
        break;
      
      case 4:
        break;
      
      case 6:
        break;
      
      case 8:
        break;
    }
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC g11L = ToReal(1);
    
    CCTK_REAL_VEC g12L = ToReal(0);
    
    CCTK_REAL_VEC g13L = ToReal(0);
    
    CCTK_REAL_VEC g22L = ToReal(1);
    
    CCTK_REAL_VEC g23L = ToReal(0);
    
    CCTK_REAL_VEC g33L = ToReal(1);
    
    CCTK_REAL_VEC K11L = ToReal(0);
    
    CCTK_REAL_VEC K12L = ToReal(0);
    
    CCTK_REAL_VEC K13L = ToReal(0);
    
    CCTK_REAL_VEC K22L = ToReal(0);
    
    CCTK_REAL_VEC K23L = ToReal(0);
    
    CCTK_REAL_VEC K33L = ToReal(0);
    
    CCTK_REAL_VEC alphaL = ToReal(1);
    
    CCTK_REAL_VEC beta1L = ToReal(0);
    
    CCTK_REAL_VEC beta2L = ToReal(0);
    
    CCTK_REAL_VEC beta3L = ToReal(0);
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(alpha[index],alphaL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta1[index],beta1L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta2[index],beta2L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta3[index],beta3L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g11[index],g11L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g12[index],g12L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g13[index],g13L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g22[index],g22L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g23[index],g23L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(g33[index],g33L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K11[index],K11L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K12[index],K12L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K13[index],K13L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K22[index],K22L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K23[index],K23L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(K33[index],K33L,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(alpha[index],alphaL,elt_count);
      vec_store_nta_partial_hi(beta1[index],beta1L,elt_count);
      vec_store_nta_partial_hi(beta2[index],beta2L,elt_count);
      vec_store_nta_partial_hi(beta3[index],beta3L,elt_count);
      vec_store_nta_partial_hi(g11[index],g11L,elt_count);
      vec_store_nta_partial_hi(g12[index],g12L,elt_count);
      vec_store_nta_partial_hi(g13[index],g13L,elt_count);
      vec_store_nta_partial_hi(g22[index],g22L,elt_count);
      vec_store_nta_partial_hi(g23[index],g23L,elt_count);
      vec_store_nta_partial_hi(g33[index],g33L,elt_count);
      vec_store_nta_partial_hi(K11[index],K11L,elt_count);
      vec_store_nta_partial_hi(K12[index],K12L,elt_count);
      vec_store_nta_partial_hi(K13[index],K13L,elt_count);
      vec_store_nta_partial_hi(K22[index],K22L,elt_count);
      vec_store_nta_partial_hi(K23[index],K23L,elt_count);
      vec_store_nta_partial_hi(K33[index],K33L,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(alpha[index],alphaL,elt_count);
      vec_store_nta_partial_lo(beta1[index],beta1L,elt_count);
      vec_store_nta_partial_lo(beta2[index],beta2L,elt_count);
      vec_store_nta_partial_lo(beta3[index],beta3L,elt_count);
      vec_store_nta_partial_lo(g11[index],g11L,elt_count);
      vec_store_nta_partial_lo(g12[index],g12L,elt_count);
      vec_store_nta_partial_lo(g13[index],g13L,elt_count);
      vec_store_nta_partial_lo(g22[index],g22L,elt_count);
      vec_store_nta_partial_lo(g23[index],g23L,elt_count);
      vec_store_nta_partial_lo(g33[index],g33L,elt_count);
      vec_store_nta_partial_lo(K11[index],K11L,elt_count);
      vec_store_nta_partial_lo(K12[index],K12L,elt_count);
      vec_store_nta_partial_lo(K13[index],K13L,elt_count);
      vec_store_nta_partial_lo(K22[index],K22L,elt_count);
      vec_store_nta_partial_lo(K23[index],K23L,elt_count);
      vec_store_nta_partial_lo(K33[index],K33L,elt_count);
      break;
    }
    vec_store_nta(alpha[index],alphaL);
    vec_store_nta(beta1[index],beta1L);
    vec_store_nta(beta2[index],beta2L);
    vec_store_nta(beta3[index],beta3L);
    vec_store_nta(g11[index],g11L);
    vec_store_nta(g12[index],g12L);
    vec_store_nta(g13[index],g13L);
    vec_store_nta(g22[index],g22L);
    vec_store_nta(g23[index],g23L);
    vec_store_nta(g33[index],g33L);
    vec_store_nta(K11[index],K11L);
    vec_store_nta(K12[index],K12L);
    vec_store_nta(K13[index],K13L);
    vec_store_nta(K22[index],K22L);
    vec_store_nta(K23[index],K23L);
    vec_store_nta(K33[index],K33L);
  }
  LC_ENDLOOP3VEC (ML_ADM_boundary);
}

extern "C" void ML_ADM_boundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADM_boundary_Body");
  }
  
  if (cctk_iteration % ML_ADM_boundary_calc_every != ML_ADM_boundary_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_ADM::ML_curv","ML_ADM::ML_lapse","ML_ADM::ML_metric","ML_ADM::ML_shift"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADM_boundary", 4, groups);
  
  switch(fdOrder)
  {
    case 2:
      break;
    
    case 4:
      break;
    
    case 6:
      break;
    
    case 8:
      break;
  }
  
  GenericFD_LoopOverBoundaryWithGhosts(cctkGH, &ML_ADM_boundary_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_ADM_boundary_Body");
  }
}
