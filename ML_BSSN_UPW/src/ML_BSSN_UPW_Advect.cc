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

static void ML_BSSN_UPW_Advect_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
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
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_UPW_Advect", 3, 3, 3);
  
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
  LC_LOOP3VEC (ML_BSSN_UPW_Advect,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
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
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    phirhsL = 
      kmadd(PDupwindNth1(&phi[index]),beta1L,kmadd(PDupwindNth2(&phi[index]),beta2L,kmadd(PDupwindNth3(&phi[index]),beta3L,phirhsL)));
    
    gt11rhsL = 
      kmadd(PDupwindNth1(&gt11[index]),beta1L,kmadd(PDupwindNth2(&gt11[index]),beta2L,kmadd(PDupwindNth3(&gt11[index]),beta3L,gt11rhsL)));
    
    gt12rhsL = 
      kmadd(PDupwindNth1(&gt12[index]),beta1L,kmadd(PDupwindNth2(&gt12[index]),beta2L,kmadd(PDupwindNth3(&gt12[index]),beta3L,gt12rhsL)));
    
    gt13rhsL = 
      kmadd(PDupwindNth1(&gt13[index]),beta1L,kmadd(PDupwindNth2(&gt13[index]),beta2L,kmadd(PDupwindNth3(&gt13[index]),beta3L,gt13rhsL)));
    
    gt22rhsL = 
      kmadd(PDupwindNth1(&gt22[index]),beta1L,kmadd(PDupwindNth2(&gt22[index]),beta2L,kmadd(PDupwindNth3(&gt22[index]),beta3L,gt22rhsL)));
    
    gt23rhsL = 
      kmadd(PDupwindNth1(&gt23[index]),beta1L,kmadd(PDupwindNth2(&gt23[index]),beta2L,kmadd(PDupwindNth3(&gt23[index]),beta3L,gt23rhsL)));
    
    gt33rhsL = 
      kmadd(PDupwindNth1(&gt33[index]),beta1L,kmadd(PDupwindNth2(&gt33[index]),beta2L,kmadd(PDupwindNth3(&gt33[index]),beta3L,gt33rhsL)));
    
    Xt1rhsL = 
      kmadd(PDupwindNth1(&Xt1[index]),beta1L,kmadd(PDupwindNth2(&Xt1[index]),beta2L,kmadd(PDupwindNth3(&Xt1[index]),beta3L,Xt1rhsL)));
    
    Xt2rhsL = 
      kmadd(PDupwindNth1(&Xt2[index]),beta1L,kmadd(PDupwindNth2(&Xt2[index]),beta2L,kmadd(PDupwindNth3(&Xt2[index]),beta3L,Xt2rhsL)));
    
    Xt3rhsL = 
      kmadd(PDupwindNth1(&Xt3[index]),beta1L,kmadd(PDupwindNth2(&Xt3[index]),beta2L,kmadd(PDupwindNth3(&Xt3[index]),beta3L,Xt3rhsL)));
    
    trKrhsL = 
      kmadd(PDupwindNth1(&trK[index]),beta1L,kmadd(PDupwindNth2(&trK[index]),beta2L,kmadd(PDupwindNth3(&trK[index]),beta3L,trKrhsL)));
    
    At11rhsL = 
      kadd(At11rhsL,kmadd(PDupwindNth1(&At11[index]),beta1L,kmadd(PDupwindNth2(&At11[index]),beta2L,kmul(PDupwindNth3(&At11[index]),beta3L))));
    
    At12rhsL = 
      kadd(At12rhsL,kmadd(PDupwindNth1(&At12[index]),beta1L,kmadd(PDupwindNth2(&At12[index]),beta2L,kmul(PDupwindNth3(&At12[index]),beta3L))));
    
    At13rhsL = 
      kadd(At13rhsL,kmadd(PDupwindNth1(&At13[index]),beta1L,kmadd(PDupwindNth2(&At13[index]),beta2L,kmul(PDupwindNth3(&At13[index]),beta3L))));
    
    At22rhsL = 
      kadd(At22rhsL,kmadd(PDupwindNth1(&At22[index]),beta1L,kmadd(PDupwindNth2(&At22[index]),beta2L,kmul(PDupwindNth3(&At22[index]),beta3L))));
    
    At23rhsL = 
      kadd(At23rhsL,kmadd(PDupwindNth1(&At23[index]),beta1L,kmadd(PDupwindNth2(&At23[index]),beta2L,kmul(PDupwindNth3(&At23[index]),beta3L))));
    
    At33rhsL = 
      kadd(At33rhsL,kmadd(PDupwindNth1(&At33[index]),beta1L,kmadd(PDupwindNth2(&At33[index]),beta2L,kmul(PDupwindNth3(&At33[index]),beta3L))));
    
    alpharhsL = 
      kmadd(kmadd(PDupwindNth1(&alpha[index]),beta1L,kmadd(PDupwindNth2(&alpha[index]),beta2L,kmul(PDupwindNth3(&alpha[index]),beta3L))),ToReal(LapseAdvectionCoeff),alpharhsL);
    
    ArhsL = 
      kmadd(kmadd(PDupwindNth1(&A[index]),beta1L,kmadd(PDupwindNth2(&A[index]),beta2L,kmul(PDupwindNth3(&A[index]),beta3L))),ToReal(LapseAdvectionCoeff),ArhsL);
    
    beta1rhsL = 
      kmadd(kmadd(PDupwindNth1(&beta1[index]),beta1L,kmadd(PDupwindNth2(&beta1[index]),beta2L,kmul(PDupwindNth3(&beta1[index]),beta3L))),ToReal(ShiftAdvectionCoeff),beta1rhsL);
    
    beta2rhsL = 
      kmadd(kmadd(PDupwindNth1(&beta2[index]),beta1L,kmadd(PDupwindNth2(&beta2[index]),beta2L,kmul(PDupwindNth3(&beta2[index]),beta3L))),ToReal(ShiftAdvectionCoeff),beta2rhsL);
    
    beta3rhsL = 
      kmadd(kmadd(PDupwindNth1(&beta3[index]),beta1L,kmadd(PDupwindNth2(&beta3[index]),beta2L,kmul(PDupwindNth3(&beta3[index]),beta3L))),ToReal(ShiftAdvectionCoeff),beta3rhsL);
    
    B1rhsL = 
      kadd(B1rhsL,kmadd(beta1L,kmadd(ksub(PDupwindNth1(&B1[index]),PDupwindNth1(&Xt1[index])),ToReal(ShiftAdvectionCoeff),kmul(PDupwindNth1(&Xt1[index]),ToReal(ShiftBCoeff))),kmadd(beta2L,kmadd(ksub(PDupwindNth2(&B1[index]),PDupwindNth2(&Xt1[index])),ToReal(ShiftAdvectionCoeff),kmul(PDupwindNth2(&Xt1[index]),ToReal(ShiftBCoeff))),kmul(beta3L,kmadd(ksub(PDupwindNth3(&B1[index]),PDupwindNth3(&Xt1[index])),ToReal(ShiftAdvectionCoeff),kmul(PDupwindNth3(&Xt1[index]),ToReal(ShiftBCoeff)))))));
    
    B2rhsL = 
      kadd(B2rhsL,kmadd(beta1L,kmadd(ksub(PDupwindNth1(&B2[index]),PDupwindNth1(&Xt2[index])),ToReal(ShiftAdvectionCoeff),kmul(PDupwindNth1(&Xt2[index]),ToReal(ShiftBCoeff))),kmadd(beta2L,kmadd(ksub(PDupwindNth2(&B2[index]),PDupwindNth2(&Xt2[index])),ToReal(ShiftAdvectionCoeff),kmul(PDupwindNth2(&Xt2[index]),ToReal(ShiftBCoeff))),kmul(beta3L,kmadd(ksub(PDupwindNth3(&B2[index]),PDupwindNth3(&Xt2[index])),ToReal(ShiftAdvectionCoeff),kmul(PDupwindNth3(&Xt2[index]),ToReal(ShiftBCoeff)))))));
    
    B3rhsL = 
      kadd(B3rhsL,kmadd(beta1L,kmadd(ksub(PDupwindNth1(&B3[index]),PDupwindNth1(&Xt3[index])),ToReal(ShiftAdvectionCoeff),kmul(PDupwindNth1(&Xt3[index]),ToReal(ShiftBCoeff))),kmadd(beta2L,kmadd(ksub(PDupwindNth2(&B3[index]),PDupwindNth2(&Xt3[index])),ToReal(ShiftAdvectionCoeff),kmul(PDupwindNth2(&Xt3[index]),ToReal(ShiftBCoeff))),kmul(beta3L,kmadd(ksub(PDupwindNth3(&B3[index]),PDupwindNth3(&Xt3[index])),ToReal(ShiftAdvectionCoeff),kmul(PDupwindNth3(&Xt3[index]),ToReal(ShiftBCoeff)))))));
    
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
    
    /* Copy local copies back to grid functions */
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
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_UPW_Advect_Body);
}
