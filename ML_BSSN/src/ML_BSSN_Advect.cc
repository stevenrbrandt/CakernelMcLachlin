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

static void ML_BSSN_Advect_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
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
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Advect", 3, 3, 3);
  
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
  CCTK_REAL const p1o12dx = 0.0833333333333333333333333333333*INV(dx);
  CCTK_REAL const p1o12dy = 0.0833333333333333333333333333333*INV(dy);
  CCTK_REAL const p1o12dz = 0.0833333333333333333333333333333*INV(dz);
  CCTK_REAL const p1o144dxdy = 0.00694444444444444444444444444444*INV(dx)*INV(dy);
  CCTK_REAL const p1o144dxdz = 0.00694444444444444444444444444444*INV(dx)*INV(dz);
  CCTK_REAL const p1o144dydz = 0.00694444444444444444444444444444*INV(dy)*INV(dz);
  CCTK_REAL const p1o24dx = 0.0416666666666666666666666666667*INV(dx);
  CCTK_REAL const p1o24dy = 0.0416666666666666666666666666667*INV(dy);
  CCTK_REAL const p1o24dz = 0.0416666666666666666666666666667*INV(dz);
  CCTK_REAL const p1o64dx = 0.015625*INV(dx);
  CCTK_REAL const p1o64dy = 0.015625*INV(dy);
  CCTK_REAL const p1o64dz = 0.015625*INV(dz);
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o12dx2 = -0.0833333333333333333333333333333*INV(SQR(dx));
  CCTK_REAL const pm1o12dy2 = -0.0833333333333333333333333333333*INV(SQR(dy));
  CCTK_REAL const pm1o12dz2 = -0.0833333333333333333333333333333*INV(SQR(dz));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_Advect,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL AL = A[index];
    CCTK_REAL alphaL = alpha[index];
    CCTK_REAL alpharhsL = alpharhs[index];
    CCTK_REAL ArhsL = Arhs[index];
    CCTK_REAL At11L = At11[index];
    CCTK_REAL At11rhsL = At11rhs[index];
    CCTK_REAL At12L = At12[index];
    CCTK_REAL At12rhsL = At12rhs[index];
    CCTK_REAL At13L = At13[index];
    CCTK_REAL At13rhsL = At13rhs[index];
    CCTK_REAL At22L = At22[index];
    CCTK_REAL At22rhsL = At22rhs[index];
    CCTK_REAL At23L = At23[index];
    CCTK_REAL At23rhsL = At23rhs[index];
    CCTK_REAL At33L = At33[index];
    CCTK_REAL At33rhsL = At33rhs[index];
    CCTK_REAL B1L = B1[index];
    CCTK_REAL B1rhsL = B1rhs[index];
    CCTK_REAL B2L = B2[index];
    CCTK_REAL B2rhsL = B2rhs[index];
    CCTK_REAL B3L = B3[index];
    CCTK_REAL B3rhsL = B3rhs[index];
    CCTK_REAL beta1L = beta1[index];
    CCTK_REAL beta1rhsL = beta1rhs[index];
    CCTK_REAL beta2L = beta2[index];
    CCTK_REAL beta2rhsL = beta2rhs[index];
    CCTK_REAL beta3L = beta3[index];
    CCTK_REAL beta3rhsL = beta3rhs[index];
    CCTK_REAL gt11L = gt11[index];
    CCTK_REAL gt11rhsL = gt11rhs[index];
    CCTK_REAL gt12L = gt12[index];
    CCTK_REAL gt12rhsL = gt12rhs[index];
    CCTK_REAL gt13L = gt13[index];
    CCTK_REAL gt13rhsL = gt13rhs[index];
    CCTK_REAL gt22L = gt22[index];
    CCTK_REAL gt22rhsL = gt22rhs[index];
    CCTK_REAL gt23L = gt23[index];
    CCTK_REAL gt23rhsL = gt23rhs[index];
    CCTK_REAL gt33L = gt33[index];
    CCTK_REAL gt33rhsL = gt33rhs[index];
    CCTK_REAL phiL = phi[index];
    CCTK_REAL phirhsL = phirhs[index];
    CCTK_REAL trKL = trK[index];
    CCTK_REAL trKrhsL = trKrhs[index];
    CCTK_REAL Xt1L = Xt1[index];
    CCTK_REAL Xt1rhsL = Xt1rhs[index];
    CCTK_REAL Xt2L = Xt2[index];
    CCTK_REAL Xt2rhsL = Xt2rhs[index];
    CCTK_REAL Xt3L = Xt3[index];
    CCTK_REAL Xt3rhsL = Xt3rhs[index];
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDupwindNthAnti1A = PDupwindNthAnti1(&A[index]);
    CCTK_REAL const PDupwindNthSymm1A = PDupwindNthSymm1(&A[index]);
    CCTK_REAL const PDupwindNthAnti2A = PDupwindNthAnti2(&A[index]);
    CCTK_REAL const PDupwindNthSymm2A = PDupwindNthSymm2(&A[index]);
    CCTK_REAL const PDupwindNthAnti3A = PDupwindNthAnti3(&A[index]);
    CCTK_REAL const PDupwindNthSymm3A = PDupwindNthSymm3(&A[index]);
    CCTK_REAL const PDupwindNthAnti1alpha = PDupwindNthAnti1(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm1alpha = PDupwindNthSymm1(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti2alpha = PDupwindNthAnti2(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm2alpha = PDupwindNthSymm2(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti3alpha = PDupwindNthAnti3(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm3alpha = PDupwindNthSymm3(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti1At11 = PDupwindNthAnti1(&At11[index]);
    CCTK_REAL const PDupwindNthSymm1At11 = PDupwindNthSymm1(&At11[index]);
    CCTK_REAL const PDupwindNthAnti2At11 = PDupwindNthAnti2(&At11[index]);
    CCTK_REAL const PDupwindNthSymm2At11 = PDupwindNthSymm2(&At11[index]);
    CCTK_REAL const PDupwindNthAnti3At11 = PDupwindNthAnti3(&At11[index]);
    CCTK_REAL const PDupwindNthSymm3At11 = PDupwindNthSymm3(&At11[index]);
    CCTK_REAL const PDupwindNthAnti1At12 = PDupwindNthAnti1(&At12[index]);
    CCTK_REAL const PDupwindNthSymm1At12 = PDupwindNthSymm1(&At12[index]);
    CCTK_REAL const PDupwindNthAnti2At12 = PDupwindNthAnti2(&At12[index]);
    CCTK_REAL const PDupwindNthSymm2At12 = PDupwindNthSymm2(&At12[index]);
    CCTK_REAL const PDupwindNthAnti3At12 = PDupwindNthAnti3(&At12[index]);
    CCTK_REAL const PDupwindNthSymm3At12 = PDupwindNthSymm3(&At12[index]);
    CCTK_REAL const PDupwindNthAnti1At13 = PDupwindNthAnti1(&At13[index]);
    CCTK_REAL const PDupwindNthSymm1At13 = PDupwindNthSymm1(&At13[index]);
    CCTK_REAL const PDupwindNthAnti2At13 = PDupwindNthAnti2(&At13[index]);
    CCTK_REAL const PDupwindNthSymm2At13 = PDupwindNthSymm2(&At13[index]);
    CCTK_REAL const PDupwindNthAnti3At13 = PDupwindNthAnti3(&At13[index]);
    CCTK_REAL const PDupwindNthSymm3At13 = PDupwindNthSymm3(&At13[index]);
    CCTK_REAL const PDupwindNthAnti1At22 = PDupwindNthAnti1(&At22[index]);
    CCTK_REAL const PDupwindNthSymm1At22 = PDupwindNthSymm1(&At22[index]);
    CCTK_REAL const PDupwindNthAnti2At22 = PDupwindNthAnti2(&At22[index]);
    CCTK_REAL const PDupwindNthSymm2At22 = PDupwindNthSymm2(&At22[index]);
    CCTK_REAL const PDupwindNthAnti3At22 = PDupwindNthAnti3(&At22[index]);
    CCTK_REAL const PDupwindNthSymm3At22 = PDupwindNthSymm3(&At22[index]);
    CCTK_REAL const PDupwindNthAnti1At23 = PDupwindNthAnti1(&At23[index]);
    CCTK_REAL const PDupwindNthSymm1At23 = PDupwindNthSymm1(&At23[index]);
    CCTK_REAL const PDupwindNthAnti2At23 = PDupwindNthAnti2(&At23[index]);
    CCTK_REAL const PDupwindNthSymm2At23 = PDupwindNthSymm2(&At23[index]);
    CCTK_REAL const PDupwindNthAnti3At23 = PDupwindNthAnti3(&At23[index]);
    CCTK_REAL const PDupwindNthSymm3At23 = PDupwindNthSymm3(&At23[index]);
    CCTK_REAL const PDupwindNthAnti1At33 = PDupwindNthAnti1(&At33[index]);
    CCTK_REAL const PDupwindNthSymm1At33 = PDupwindNthSymm1(&At33[index]);
    CCTK_REAL const PDupwindNthAnti2At33 = PDupwindNthAnti2(&At33[index]);
    CCTK_REAL const PDupwindNthSymm2At33 = PDupwindNthSymm2(&At33[index]);
    CCTK_REAL const PDupwindNthAnti3At33 = PDupwindNthAnti3(&At33[index]);
    CCTK_REAL const PDupwindNthSymm3At33 = PDupwindNthSymm3(&At33[index]);
    CCTK_REAL const PDupwindNthAnti1B1 = PDupwindNthAnti1(&B1[index]);
    CCTK_REAL const PDupwindNthSymm1B1 = PDupwindNthSymm1(&B1[index]);
    CCTK_REAL const PDupwindNthAnti2B1 = PDupwindNthAnti2(&B1[index]);
    CCTK_REAL const PDupwindNthSymm2B1 = PDupwindNthSymm2(&B1[index]);
    CCTK_REAL const PDupwindNthAnti3B1 = PDupwindNthAnti3(&B1[index]);
    CCTK_REAL const PDupwindNthSymm3B1 = PDupwindNthSymm3(&B1[index]);
    CCTK_REAL const PDupwindNthAnti1B2 = PDupwindNthAnti1(&B2[index]);
    CCTK_REAL const PDupwindNthSymm1B2 = PDupwindNthSymm1(&B2[index]);
    CCTK_REAL const PDupwindNthAnti2B2 = PDupwindNthAnti2(&B2[index]);
    CCTK_REAL const PDupwindNthSymm2B2 = PDupwindNthSymm2(&B2[index]);
    CCTK_REAL const PDupwindNthAnti3B2 = PDupwindNthAnti3(&B2[index]);
    CCTK_REAL const PDupwindNthSymm3B2 = PDupwindNthSymm3(&B2[index]);
    CCTK_REAL const PDupwindNthAnti1B3 = PDupwindNthAnti1(&B3[index]);
    CCTK_REAL const PDupwindNthSymm1B3 = PDupwindNthSymm1(&B3[index]);
    CCTK_REAL const PDupwindNthAnti2B3 = PDupwindNthAnti2(&B3[index]);
    CCTK_REAL const PDupwindNthSymm2B3 = PDupwindNthSymm2(&B3[index]);
    CCTK_REAL const PDupwindNthAnti3B3 = PDupwindNthAnti3(&B3[index]);
    CCTK_REAL const PDupwindNthSymm3B3 = PDupwindNthSymm3(&B3[index]);
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
    CCTK_REAL const PDupwindNthAnti1gt11 = PDupwindNthAnti1(&gt11[index]);
    CCTK_REAL const PDupwindNthSymm1gt11 = PDupwindNthSymm1(&gt11[index]);
    CCTK_REAL const PDupwindNthAnti2gt11 = PDupwindNthAnti2(&gt11[index]);
    CCTK_REAL const PDupwindNthSymm2gt11 = PDupwindNthSymm2(&gt11[index]);
    CCTK_REAL const PDupwindNthAnti3gt11 = PDupwindNthAnti3(&gt11[index]);
    CCTK_REAL const PDupwindNthSymm3gt11 = PDupwindNthSymm3(&gt11[index]);
    CCTK_REAL const PDupwindNthAnti1gt12 = PDupwindNthAnti1(&gt12[index]);
    CCTK_REAL const PDupwindNthSymm1gt12 = PDupwindNthSymm1(&gt12[index]);
    CCTK_REAL const PDupwindNthAnti2gt12 = PDupwindNthAnti2(&gt12[index]);
    CCTK_REAL const PDupwindNthSymm2gt12 = PDupwindNthSymm2(&gt12[index]);
    CCTK_REAL const PDupwindNthAnti3gt12 = PDupwindNthAnti3(&gt12[index]);
    CCTK_REAL const PDupwindNthSymm3gt12 = PDupwindNthSymm3(&gt12[index]);
    CCTK_REAL const PDupwindNthAnti1gt13 = PDupwindNthAnti1(&gt13[index]);
    CCTK_REAL const PDupwindNthSymm1gt13 = PDupwindNthSymm1(&gt13[index]);
    CCTK_REAL const PDupwindNthAnti2gt13 = PDupwindNthAnti2(&gt13[index]);
    CCTK_REAL const PDupwindNthSymm2gt13 = PDupwindNthSymm2(&gt13[index]);
    CCTK_REAL const PDupwindNthAnti3gt13 = PDupwindNthAnti3(&gt13[index]);
    CCTK_REAL const PDupwindNthSymm3gt13 = PDupwindNthSymm3(&gt13[index]);
    CCTK_REAL const PDupwindNthAnti1gt22 = PDupwindNthAnti1(&gt22[index]);
    CCTK_REAL const PDupwindNthSymm1gt22 = PDupwindNthSymm1(&gt22[index]);
    CCTK_REAL const PDupwindNthAnti2gt22 = PDupwindNthAnti2(&gt22[index]);
    CCTK_REAL const PDupwindNthSymm2gt22 = PDupwindNthSymm2(&gt22[index]);
    CCTK_REAL const PDupwindNthAnti3gt22 = PDupwindNthAnti3(&gt22[index]);
    CCTK_REAL const PDupwindNthSymm3gt22 = PDupwindNthSymm3(&gt22[index]);
    CCTK_REAL const PDupwindNthAnti1gt23 = PDupwindNthAnti1(&gt23[index]);
    CCTK_REAL const PDupwindNthSymm1gt23 = PDupwindNthSymm1(&gt23[index]);
    CCTK_REAL const PDupwindNthAnti2gt23 = PDupwindNthAnti2(&gt23[index]);
    CCTK_REAL const PDupwindNthSymm2gt23 = PDupwindNthSymm2(&gt23[index]);
    CCTK_REAL const PDupwindNthAnti3gt23 = PDupwindNthAnti3(&gt23[index]);
    CCTK_REAL const PDupwindNthSymm3gt23 = PDupwindNthSymm3(&gt23[index]);
    CCTK_REAL const PDupwindNthAnti1gt33 = PDupwindNthAnti1(&gt33[index]);
    CCTK_REAL const PDupwindNthSymm1gt33 = PDupwindNthSymm1(&gt33[index]);
    CCTK_REAL const PDupwindNthAnti2gt33 = PDupwindNthAnti2(&gt33[index]);
    CCTK_REAL const PDupwindNthSymm2gt33 = PDupwindNthSymm2(&gt33[index]);
    CCTK_REAL const PDupwindNthAnti3gt33 = PDupwindNthAnti3(&gt33[index]);
    CCTK_REAL const PDupwindNthSymm3gt33 = PDupwindNthSymm3(&gt33[index]);
    CCTK_REAL const PDupwindNthAnti1phi = PDupwindNthAnti1(&phi[index]);
    CCTK_REAL const PDupwindNthSymm1phi = PDupwindNthSymm1(&phi[index]);
    CCTK_REAL const PDupwindNthAnti2phi = PDupwindNthAnti2(&phi[index]);
    CCTK_REAL const PDupwindNthSymm2phi = PDupwindNthSymm2(&phi[index]);
    CCTK_REAL const PDupwindNthAnti3phi = PDupwindNthAnti3(&phi[index]);
    CCTK_REAL const PDupwindNthSymm3phi = PDupwindNthSymm3(&phi[index]);
    CCTK_REAL const PDupwindNthAnti1trK = PDupwindNthAnti1(&trK[index]);
    CCTK_REAL const PDupwindNthSymm1trK = PDupwindNthSymm1(&trK[index]);
    CCTK_REAL const PDupwindNthAnti2trK = PDupwindNthAnti2(&trK[index]);
    CCTK_REAL const PDupwindNthSymm2trK = PDupwindNthSymm2(&trK[index]);
    CCTK_REAL const PDupwindNthAnti3trK = PDupwindNthAnti3(&trK[index]);
    CCTK_REAL const PDupwindNthSymm3trK = PDupwindNthSymm3(&trK[index]);
    CCTK_REAL const PDupwindNthAnti1Xt1 = PDupwindNthAnti1(&Xt1[index]);
    CCTK_REAL const PDupwindNthSymm1Xt1 = PDupwindNthSymm1(&Xt1[index]);
    CCTK_REAL const PDupwindNthAnti2Xt1 = PDupwindNthAnti2(&Xt1[index]);
    CCTK_REAL const PDupwindNthSymm2Xt1 = PDupwindNthSymm2(&Xt1[index]);
    CCTK_REAL const PDupwindNthAnti3Xt1 = PDupwindNthAnti3(&Xt1[index]);
    CCTK_REAL const PDupwindNthSymm3Xt1 = PDupwindNthSymm3(&Xt1[index]);
    CCTK_REAL const PDupwindNthAnti1Xt2 = PDupwindNthAnti1(&Xt2[index]);
    CCTK_REAL const PDupwindNthSymm1Xt2 = PDupwindNthSymm1(&Xt2[index]);
    CCTK_REAL const PDupwindNthAnti2Xt2 = PDupwindNthAnti2(&Xt2[index]);
    CCTK_REAL const PDupwindNthSymm2Xt2 = PDupwindNthSymm2(&Xt2[index]);
    CCTK_REAL const PDupwindNthAnti3Xt2 = PDupwindNthAnti3(&Xt2[index]);
    CCTK_REAL const PDupwindNthSymm3Xt2 = PDupwindNthSymm3(&Xt2[index]);
    CCTK_REAL const PDupwindNthAnti1Xt3 = PDupwindNthAnti1(&Xt3[index]);
    CCTK_REAL const PDupwindNthSymm1Xt3 = PDupwindNthSymm1(&Xt3[index]);
    CCTK_REAL const PDupwindNthAnti2Xt3 = PDupwindNthAnti2(&Xt3[index]);
    CCTK_REAL const PDupwindNthSymm2Xt3 = PDupwindNthSymm2(&Xt3[index]);
    CCTK_REAL const PDupwindNthAnti3Xt3 = PDupwindNthAnti3(&Xt3[index]);
    CCTK_REAL const PDupwindNthSymm3Xt3 = PDupwindNthSymm3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    phirhsL = beta1L*PDupwindNthAnti1phi + beta2L*PDupwindNthAnti2phi + 
      beta3L*PDupwindNthAnti3phi + phirhsL + PDupwindNthSymm1phi*Abs(beta1L) 
      + PDupwindNthSymm2phi*Abs(beta2L) + PDupwindNthSymm3phi*Abs(beta3L);
    
    gt11rhsL = gt11rhsL + beta1L*PDupwindNthAnti1gt11 + 
      beta2L*PDupwindNthAnti2gt11 + beta3L*PDupwindNthAnti3gt11 + 
      PDupwindNthSymm1gt11*Abs(beta1L) + PDupwindNthSymm2gt11*Abs(beta2L) + 
      PDupwindNthSymm3gt11*Abs(beta3L);
    
    gt12rhsL = gt12rhsL + beta1L*PDupwindNthAnti1gt12 + 
      beta2L*PDupwindNthAnti2gt12 + beta3L*PDupwindNthAnti3gt12 + 
      PDupwindNthSymm1gt12*Abs(beta1L) + PDupwindNthSymm2gt12*Abs(beta2L) + 
      PDupwindNthSymm3gt12*Abs(beta3L);
    
    gt13rhsL = gt13rhsL + beta1L*PDupwindNthAnti1gt13 + 
      beta2L*PDupwindNthAnti2gt13 + beta3L*PDupwindNthAnti3gt13 + 
      PDupwindNthSymm1gt13*Abs(beta1L) + PDupwindNthSymm2gt13*Abs(beta2L) + 
      PDupwindNthSymm3gt13*Abs(beta3L);
    
    gt22rhsL = gt22rhsL + beta1L*PDupwindNthAnti1gt22 + 
      beta2L*PDupwindNthAnti2gt22 + beta3L*PDupwindNthAnti3gt22 + 
      PDupwindNthSymm1gt22*Abs(beta1L) + PDupwindNthSymm2gt22*Abs(beta2L) + 
      PDupwindNthSymm3gt22*Abs(beta3L);
    
    gt23rhsL = gt23rhsL + beta1L*PDupwindNthAnti1gt23 + 
      beta2L*PDupwindNthAnti2gt23 + beta3L*PDupwindNthAnti3gt23 + 
      PDupwindNthSymm1gt23*Abs(beta1L) + PDupwindNthSymm2gt23*Abs(beta2L) + 
      PDupwindNthSymm3gt23*Abs(beta3L);
    
    gt33rhsL = gt33rhsL + beta1L*PDupwindNthAnti1gt33 + 
      beta2L*PDupwindNthAnti2gt33 + beta3L*PDupwindNthAnti3gt33 + 
      PDupwindNthSymm1gt33*Abs(beta1L) + PDupwindNthSymm2gt33*Abs(beta2L) + 
      PDupwindNthSymm3gt33*Abs(beta3L);
    
    Xt1rhsL = beta1L*PDupwindNthAnti1Xt1 + beta2L*PDupwindNthAnti2Xt1 + 
      beta3L*PDupwindNthAnti3Xt1 + Xt1rhsL + PDupwindNthSymm1Xt1*Abs(beta1L) 
      + PDupwindNthSymm2Xt1*Abs(beta2L) + PDupwindNthSymm3Xt1*Abs(beta3L);
    
    Xt2rhsL = beta1L*PDupwindNthAnti1Xt2 + beta2L*PDupwindNthAnti2Xt2 + 
      beta3L*PDupwindNthAnti3Xt2 + Xt2rhsL + PDupwindNthSymm1Xt2*Abs(beta1L) 
      + PDupwindNthSymm2Xt2*Abs(beta2L) + PDupwindNthSymm3Xt2*Abs(beta3L);
    
    Xt3rhsL = beta1L*PDupwindNthAnti1Xt3 + beta2L*PDupwindNthAnti2Xt3 + 
      beta3L*PDupwindNthAnti3Xt3 + Xt3rhsL + PDupwindNthSymm1Xt3*Abs(beta1L) 
      + PDupwindNthSymm2Xt3*Abs(beta2L) + PDupwindNthSymm3Xt3*Abs(beta3L);
    
    trKrhsL = beta1L*PDupwindNthAnti1trK + beta2L*PDupwindNthAnti2trK + 
      beta3L*PDupwindNthAnti3trK + trKrhsL + PDupwindNthSymm1trK*Abs(beta1L) 
      + PDupwindNthSymm2trK*Abs(beta2L) + PDupwindNthSymm3trK*Abs(beta3L);
    
    At11rhsL = At11rhsL + beta1L*PDupwindNthAnti1At11 + 
      beta2L*PDupwindNthAnti2At11 + beta3L*PDupwindNthAnti3At11 + 
      PDupwindNthSymm1At11*Abs(beta1L) + PDupwindNthSymm2At11*Abs(beta2L) + 
      PDupwindNthSymm3At11*Abs(beta3L);
    
    At12rhsL = At12rhsL + beta1L*PDupwindNthAnti1At12 + 
      beta2L*PDupwindNthAnti2At12 + beta3L*PDupwindNthAnti3At12 + 
      PDupwindNthSymm1At12*Abs(beta1L) + PDupwindNthSymm2At12*Abs(beta2L) + 
      PDupwindNthSymm3At12*Abs(beta3L);
    
    At13rhsL = At13rhsL + beta1L*PDupwindNthAnti1At13 + 
      beta2L*PDupwindNthAnti2At13 + beta3L*PDupwindNthAnti3At13 + 
      PDupwindNthSymm1At13*Abs(beta1L) + PDupwindNthSymm2At13*Abs(beta2L) + 
      PDupwindNthSymm3At13*Abs(beta3L);
    
    At22rhsL = At22rhsL + beta1L*PDupwindNthAnti1At22 + 
      beta2L*PDupwindNthAnti2At22 + beta3L*PDupwindNthAnti3At22 + 
      PDupwindNthSymm1At22*Abs(beta1L) + PDupwindNthSymm2At22*Abs(beta2L) + 
      PDupwindNthSymm3At22*Abs(beta3L);
    
    At23rhsL = At23rhsL + beta1L*PDupwindNthAnti1At23 + 
      beta2L*PDupwindNthAnti2At23 + beta3L*PDupwindNthAnti3At23 + 
      PDupwindNthSymm1At23*Abs(beta1L) + PDupwindNthSymm2At23*Abs(beta2L) + 
      PDupwindNthSymm3At23*Abs(beta3L);
    
    At33rhsL = At33rhsL + beta1L*PDupwindNthAnti1At33 + 
      beta2L*PDupwindNthAnti2At33 + beta3L*PDupwindNthAnti3At33 + 
      PDupwindNthSymm1At33*Abs(beta1L) + PDupwindNthSymm2At33*Abs(beta2L) + 
      PDupwindNthSymm3At33*Abs(beta3L);
    
    alpharhsL = alpharhsL + (beta1L*PDupwindNthAnti1alpha + 
      beta2L*PDupwindNthAnti2alpha + beta3L*PDupwindNthAnti3alpha + 
      PDupwindNthSymm1alpha*Abs(beta1L) + PDupwindNthSymm2alpha*Abs(beta2L) + 
      PDupwindNthSymm3alpha*Abs(beta3L))*ToReal(LapseAdvectionCoeff);
    
    ArhsL = ArhsL + (beta1L*PDupwindNthAnti1A + beta2L*PDupwindNthAnti2A + 
      beta3L*PDupwindNthAnti3A + PDupwindNthSymm1A*Abs(beta1L) + 
      PDupwindNthSymm2A*Abs(beta2L) + 
      PDupwindNthSymm3A*Abs(beta3L))*ToReal(LapseAdvectionCoeff);
    
    beta1rhsL = beta1rhsL + (beta1L*PDupwindNthAnti1beta1 + 
      beta2L*PDupwindNthAnti2beta1 + beta3L*PDupwindNthAnti3beta1 + 
      PDupwindNthSymm1beta1*Abs(beta1L) + PDupwindNthSymm2beta1*Abs(beta2L) + 
      PDupwindNthSymm3beta1*Abs(beta3L))*ToReal(ShiftAdvectionCoeff);
    
    beta2rhsL = beta2rhsL + (beta1L*PDupwindNthAnti1beta2 + 
      beta2L*PDupwindNthAnti2beta2 + beta3L*PDupwindNthAnti3beta2 + 
      PDupwindNthSymm1beta2*Abs(beta1L) + PDupwindNthSymm2beta2*Abs(beta2L) + 
      PDupwindNthSymm3beta2*Abs(beta3L))*ToReal(ShiftAdvectionCoeff);
    
    beta3rhsL = beta3rhsL + (beta1L*PDupwindNthAnti1beta3 + 
      beta2L*PDupwindNthAnti2beta3 + beta3L*PDupwindNthAnti3beta3 + 
      PDupwindNthSymm1beta3*Abs(beta1L) + PDupwindNthSymm2beta3*Abs(beta2L) + 
      PDupwindNthSymm3beta3*Abs(beta3L))*ToReal(ShiftAdvectionCoeff);
    
    B1rhsL = B1rhsL + (beta1L*(PDupwindNthAnti1B1 - PDupwindNthAnti1Xt1) + 
      beta2L*(PDupwindNthAnti2B1 - PDupwindNthAnti2Xt1) + 
      beta3L*(PDupwindNthAnti3B1 - PDupwindNthAnti3Xt1) + (PDupwindNthSymm1B1 
      - PDupwindNthSymm1Xt1)*Abs(beta1L) + (PDupwindNthSymm2B1 - 
      PDupwindNthSymm2Xt1)*Abs(beta2L) + (PDupwindNthSymm3B1 - 
      PDupwindNthSymm3Xt1)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      (beta1L*PDupwindNthAnti1Xt1 + beta2L*PDupwindNthAnti2Xt1 + 
      beta3L*PDupwindNthAnti3Xt1 + PDupwindNthSymm1Xt1*Abs(beta1L) + 
      PDupwindNthSymm2Xt1*Abs(beta2L) + 
      PDupwindNthSymm3Xt1*Abs(beta3L))*ToReal(ShiftBCoeff);
    
    B2rhsL = B2rhsL + (beta1L*(PDupwindNthAnti1B2 - PDupwindNthAnti1Xt2) + 
      beta2L*(PDupwindNthAnti2B2 - PDupwindNthAnti2Xt2) + 
      beta3L*(PDupwindNthAnti3B2 - PDupwindNthAnti3Xt2) + (PDupwindNthSymm1B2 
      - PDupwindNthSymm1Xt2)*Abs(beta1L) + (PDupwindNthSymm2B2 - 
      PDupwindNthSymm2Xt2)*Abs(beta2L) + (PDupwindNthSymm3B2 - 
      PDupwindNthSymm3Xt2)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      (beta1L*PDupwindNthAnti1Xt2 + beta2L*PDupwindNthAnti2Xt2 + 
      beta3L*PDupwindNthAnti3Xt2 + PDupwindNthSymm1Xt2*Abs(beta1L) + 
      PDupwindNthSymm2Xt2*Abs(beta2L) + 
      PDupwindNthSymm3Xt2*Abs(beta3L))*ToReal(ShiftBCoeff);
    
    B3rhsL = B3rhsL + (beta1L*(PDupwindNthAnti1B3 - PDupwindNthAnti1Xt3) + 
      beta2L*(PDupwindNthAnti2B3 - PDupwindNthAnti2Xt3) + 
      beta3L*(PDupwindNthAnti3B3 - PDupwindNthAnti3Xt3) + (PDupwindNthSymm1B3 
      - PDupwindNthSymm1Xt3)*Abs(beta1L) + (PDupwindNthSymm2B3 - 
      PDupwindNthSymm2Xt3)*Abs(beta2L) + (PDupwindNthSymm3B3 - 
      PDupwindNthSymm3Xt3)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      (beta1L*PDupwindNthAnti1Xt3 + beta2L*PDupwindNthAnti2Xt3 + 
      beta3L*PDupwindNthAnti3Xt3 + PDupwindNthSymm1Xt3*Abs(beta1L) + 
      PDupwindNthSymm2Xt3*Abs(beta2L) + 
      PDupwindNthSymm3Xt3*Abs(beta3L))*ToReal(ShiftBCoeff);
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    Arhs[index] = ArhsL;
    At11rhs[index] = At11rhsL;
    At12rhs[index] = At12rhsL;
    At13rhs[index] = At13rhsL;
    At22rhs[index] = At22rhsL;
    At23rhs[index] = At23rhsL;
    At33rhs[index] = At33rhsL;
    B1rhs[index] = B1rhsL;
    B2rhs[index] = B2rhsL;
    B3rhs[index] = B3rhsL;
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    gt11rhs[index] = gt11rhsL;
    gt12rhs[index] = gt12rhsL;
    gt13rhs[index] = gt13rhsL;
    gt22rhs[index] = gt22rhsL;
    gt23rhs[index] = gt23rhsL;
    gt33rhs[index] = gt33rhsL;
    phirhs[index] = phirhsL;
    trKrhs[index] = trKrhsL;
    Xt1rhs[index] = Xt1rhsL;
    Xt2rhs[index] = Xt2rhsL;
    Xt3rhs[index] = Xt3rhsL;
  }
  LC_ENDLOOP3 (ML_BSSN_Advect);
}

extern "C" void ML_BSSN_Advect(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_Advect_Body);
}
