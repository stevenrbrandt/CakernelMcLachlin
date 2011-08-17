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

extern "C" void ML_BSSN_UPW_boundary_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_curv","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_curv.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_dtlapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_dtshift.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_Gamma.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_lapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_lapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_log_confac","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_log_confac.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_metric","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_metric.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_shift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_shift.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_trace_curv","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_trace_curv.");
  return;
}

static void ML_BSSN_UPW_boundary_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_UPW_boundary_Body");
  }
  
  if (cctk_iteration % ML_BSSN_UPW_boundary_calc_every != ML_BSSN_UPW_boundary_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_BSSN_UPW::ML_curv","ML_BSSN_UPW::ML_dtlapse","ML_BSSN_UPW::ML_dtshift","ML_BSSN_UPW::ML_Gamma","ML_BSSN_UPW::ML_lapse","ML_BSSN_UPW::ML_log_confac","ML_BSSN_UPW::ML_metric","ML_BSSN_UPW::ML_shift","ML_BSSN_UPW::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_UPW_boundary", 9, groups);
  
  
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
  CCTK_REAL_VEC const p1o24dx = kmul(INV(dx),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o24dy = kmul(INV(dy),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o24dz = kmul(INV(dz),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o64dx = kmul(INV(dx),ToReal(0.015625));
  CCTK_REAL_VEC const p1o64dy = kmul(INV(dy),ToReal(0.015625));
  CCTK_REAL_VEC const p1o64dz = kmul(INV(dz),ToReal(0.015625));
  CCTK_REAL_VEC const p1odx = INV(dx);
  CCTK_REAL_VEC const p1ody = INV(dy);
  CCTK_REAL_VEC const p1odz = INV(dz);
  CCTK_REAL_VEC const pm1o12dx2 = kmul(INV(SQR(dx)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dy2 = kmul(INV(SQR(dy)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dz2 = kmul(INV(SQR(dz)),ToReal(-0.0833333333333333333333333333333));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC (ML_BSSN_UPW_boundary,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC phiL = IfThen(conformalMethod,ToReal(1),ToReal(0));
    
    CCTK_REAL_VEC gt11L = ToReal(1);
    
    CCTK_REAL_VEC gt12L = ToReal(0);
    
    CCTK_REAL_VEC gt13L = ToReal(0);
    
    CCTK_REAL_VEC gt22L = ToReal(1);
    
    CCTK_REAL_VEC gt23L = ToReal(0);
    
    CCTK_REAL_VEC gt33L = ToReal(1);
    
    CCTK_REAL_VEC trKL = ToReal(0);
    
    CCTK_REAL_VEC At11L = ToReal(0);
    
    CCTK_REAL_VEC At12L = ToReal(0);
    
    CCTK_REAL_VEC At13L = ToReal(0);
    
    CCTK_REAL_VEC At22L = ToReal(0);
    
    CCTK_REAL_VEC At23L = ToReal(0);
    
    CCTK_REAL_VEC At33L = ToReal(0);
    
    CCTK_REAL_VEC Xt1L = ToReal(0);
    
    CCTK_REAL_VEC Xt2L = ToReal(0);
    
    CCTK_REAL_VEC Xt3L = ToReal(0);
    
    CCTK_REAL_VEC alphaL = ToReal(1);
    
    CCTK_REAL_VEC AL = ToReal(0);
    
    CCTK_REAL_VEC beta1L = ToReal(0);
    
    CCTK_REAL_VEC beta2L = ToReal(0);
    
    CCTK_REAL_VEC beta3L = ToReal(0);
    
    CCTK_REAL_VEC B1L = ToReal(0);
    
    CCTK_REAL_VEC B2L = ToReal(0);
    
    CCTK_REAL_VEC B3L = ToReal(0);
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(A[index],AL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(alpha[index],alphaL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At11[index],At11L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At12[index],At12L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At13[index],At13L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At22[index],At22L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At23[index],At23L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(At33[index],At33L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B1[index],B1L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B2[index],B2L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B3[index],B3L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta1[index],beta1L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta2[index],beta2L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta3[index],beta3L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt11[index],gt11L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt12[index],gt12L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt13[index],gt13L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt22[index],gt22L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt23[index],gt23L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt33[index],gt33L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(phi[index],phiL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(trK[index],trKL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt1[index],Xt1L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt2[index],Xt2L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt3[index],Xt3L,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(A[index],AL,elt_count);
      vec_store_nta_partial_hi(alpha[index],alphaL,elt_count);
      vec_store_nta_partial_hi(At11[index],At11L,elt_count);
      vec_store_nta_partial_hi(At12[index],At12L,elt_count);
      vec_store_nta_partial_hi(At13[index],At13L,elt_count);
      vec_store_nta_partial_hi(At22[index],At22L,elt_count);
      vec_store_nta_partial_hi(At23[index],At23L,elt_count);
      vec_store_nta_partial_hi(At33[index],At33L,elt_count);
      vec_store_nta_partial_hi(B1[index],B1L,elt_count);
      vec_store_nta_partial_hi(B2[index],B2L,elt_count);
      vec_store_nta_partial_hi(B3[index],B3L,elt_count);
      vec_store_nta_partial_hi(beta1[index],beta1L,elt_count);
      vec_store_nta_partial_hi(beta2[index],beta2L,elt_count);
      vec_store_nta_partial_hi(beta3[index],beta3L,elt_count);
      vec_store_nta_partial_hi(gt11[index],gt11L,elt_count);
      vec_store_nta_partial_hi(gt12[index],gt12L,elt_count);
      vec_store_nta_partial_hi(gt13[index],gt13L,elt_count);
      vec_store_nta_partial_hi(gt22[index],gt22L,elt_count);
      vec_store_nta_partial_hi(gt23[index],gt23L,elt_count);
      vec_store_nta_partial_hi(gt33[index],gt33L,elt_count);
      vec_store_nta_partial_hi(phi[index],phiL,elt_count);
      vec_store_nta_partial_hi(trK[index],trKL,elt_count);
      vec_store_nta_partial_hi(Xt1[index],Xt1L,elt_count);
      vec_store_nta_partial_hi(Xt2[index],Xt2L,elt_count);
      vec_store_nta_partial_hi(Xt3[index],Xt3L,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(A[index],AL,elt_count);
      vec_store_nta_partial_lo(alpha[index],alphaL,elt_count);
      vec_store_nta_partial_lo(At11[index],At11L,elt_count);
      vec_store_nta_partial_lo(At12[index],At12L,elt_count);
      vec_store_nta_partial_lo(At13[index],At13L,elt_count);
      vec_store_nta_partial_lo(At22[index],At22L,elt_count);
      vec_store_nta_partial_lo(At23[index],At23L,elt_count);
      vec_store_nta_partial_lo(At33[index],At33L,elt_count);
      vec_store_nta_partial_lo(B1[index],B1L,elt_count);
      vec_store_nta_partial_lo(B2[index],B2L,elt_count);
      vec_store_nta_partial_lo(B3[index],B3L,elt_count);
      vec_store_nta_partial_lo(beta1[index],beta1L,elt_count);
      vec_store_nta_partial_lo(beta2[index],beta2L,elt_count);
      vec_store_nta_partial_lo(beta3[index],beta3L,elt_count);
      vec_store_nta_partial_lo(gt11[index],gt11L,elt_count);
      vec_store_nta_partial_lo(gt12[index],gt12L,elt_count);
      vec_store_nta_partial_lo(gt13[index],gt13L,elt_count);
      vec_store_nta_partial_lo(gt22[index],gt22L,elt_count);
      vec_store_nta_partial_lo(gt23[index],gt23L,elt_count);
      vec_store_nta_partial_lo(gt33[index],gt33L,elt_count);
      vec_store_nta_partial_lo(phi[index],phiL,elt_count);
      vec_store_nta_partial_lo(trK[index],trKL,elt_count);
      vec_store_nta_partial_lo(Xt1[index],Xt1L,elt_count);
      vec_store_nta_partial_lo(Xt2[index],Xt2L,elt_count);
      vec_store_nta_partial_lo(Xt3[index],Xt3L,elt_count);
      break;
    }
    
    /* Copy local copies back to grid functions */
    vec_store_nta(A[index],AL);
    vec_store_nta(alpha[index],alphaL);
    vec_store_nta(At11[index],At11L);
    vec_store_nta(At12[index],At12L);
    vec_store_nta(At13[index],At13L);
    vec_store_nta(At22[index],At22L);
    vec_store_nta(At23[index],At23L);
    vec_store_nta(At33[index],At33L);
    vec_store_nta(B1[index],B1L);
    vec_store_nta(B2[index],B2L);
    vec_store_nta(B3[index],B3L);
    vec_store_nta(beta1[index],beta1L);
    vec_store_nta(beta2[index],beta2L);
    vec_store_nta(beta3[index],beta3L);
    vec_store_nta(gt11[index],gt11L);
    vec_store_nta(gt12[index],gt12L);
    vec_store_nta(gt13[index],gt13L);
    vec_store_nta(gt22[index],gt22L);
    vec_store_nta(gt23[index],gt23L);
    vec_store_nta(gt33[index],gt33L);
    vec_store_nta(phi[index],phiL);
    vec_store_nta(trK[index],trKL);
    vec_store_nta(Xt1[index],Xt1L);
    vec_store_nta(Xt2[index],Xt2L);
    vec_store_nta(Xt3[index],Xt3L);
  }
  LC_ENDLOOP3VEC (ML_BSSN_UPW_boundary);
}

extern "C" void ML_BSSN_UPW_boundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverBoundaryWithGhosts(cctkGH, &ML_BSSN_UPW_boundary_Body);
}
