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

extern "C" void ML_BSSN_MP_Dissipation_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_curvrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_trace_curvrhs.");
  return;
}

static void ML_BSSN_MP_Dissipation_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_Dissipation_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_Dissipation_calc_every != ML_BSSN_MP_Dissipation_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"Coordinates::jacobian","ML_BSSN_MP::ML_curv","ML_BSSN_MP::ML_curvrhs","ML_BSSN_MP::ML_dtlapse","ML_BSSN_MP::ML_dtlapserhs","ML_BSSN_MP::ML_dtshift","ML_BSSN_MP::ML_dtshiftrhs","ML_BSSN_MP::ML_Gamma","ML_BSSN_MP::ML_Gammarhs","ML_BSSN_MP::ML_lapse","ML_BSSN_MP::ML_lapserhs","ML_BSSN_MP::ML_log_confac","ML_BSSN_MP::ML_log_confacrhs","ML_BSSN_MP::ML_metric","ML_BSSN_MP::ML_metricrhs","ML_BSSN_MP::ML_shift","ML_BSSN_MP::ML_shiftrhs","ML_BSSN_MP::ML_trace_curv","ML_BSSN_MP::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_Dissipation", 19, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_MP_Dissipation", 3, 3, 3);
  
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
  LC_LOOP3 (ML_BSSN_MP_Dissipation,
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
    CCTK_REAL J11L = J11[index];
    CCTK_REAL J12L = J12[index];
    CCTK_REAL J13L = J13[index];
    CCTK_REAL J21L = J21[index];
    CCTK_REAL J22L = J22[index];
    CCTK_REAL J23L = J23[index];
    CCTK_REAL J31L = J31[index];
    CCTK_REAL J32L = J32[index];
    CCTK_REAL J33L = J33[index];
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
    CCTK_REAL const PDdissipationNth1A = PDdissipationNth1(&A[index]);
    CCTK_REAL const PDdissipationNth2A = PDdissipationNth2(&A[index]);
    CCTK_REAL const PDdissipationNth3A = PDdissipationNth3(&A[index]);
    CCTK_REAL const PDdissipationNth1alpha = PDdissipationNth1(&alpha[index]);
    CCTK_REAL const PDdissipationNth2alpha = PDdissipationNth2(&alpha[index]);
    CCTK_REAL const PDdissipationNth3alpha = PDdissipationNth3(&alpha[index]);
    CCTK_REAL const PDdissipationNth1At11 = PDdissipationNth1(&At11[index]);
    CCTK_REAL const PDdissipationNth2At11 = PDdissipationNth2(&At11[index]);
    CCTK_REAL const PDdissipationNth3At11 = PDdissipationNth3(&At11[index]);
    CCTK_REAL const PDdissipationNth1At12 = PDdissipationNth1(&At12[index]);
    CCTK_REAL const PDdissipationNth2At12 = PDdissipationNth2(&At12[index]);
    CCTK_REAL const PDdissipationNth3At12 = PDdissipationNth3(&At12[index]);
    CCTK_REAL const PDdissipationNth1At13 = PDdissipationNth1(&At13[index]);
    CCTK_REAL const PDdissipationNth2At13 = PDdissipationNth2(&At13[index]);
    CCTK_REAL const PDdissipationNth3At13 = PDdissipationNth3(&At13[index]);
    CCTK_REAL const PDdissipationNth1At22 = PDdissipationNth1(&At22[index]);
    CCTK_REAL const PDdissipationNth2At22 = PDdissipationNth2(&At22[index]);
    CCTK_REAL const PDdissipationNth3At22 = PDdissipationNth3(&At22[index]);
    CCTK_REAL const PDdissipationNth1At23 = PDdissipationNth1(&At23[index]);
    CCTK_REAL const PDdissipationNth2At23 = PDdissipationNth2(&At23[index]);
    CCTK_REAL const PDdissipationNth3At23 = PDdissipationNth3(&At23[index]);
    CCTK_REAL const PDdissipationNth1At33 = PDdissipationNth1(&At33[index]);
    CCTK_REAL const PDdissipationNth2At33 = PDdissipationNth2(&At33[index]);
    CCTK_REAL const PDdissipationNth3At33 = PDdissipationNth3(&At33[index]);
    CCTK_REAL const PDdissipationNth1B1 = PDdissipationNth1(&B1[index]);
    CCTK_REAL const PDdissipationNth2B1 = PDdissipationNth2(&B1[index]);
    CCTK_REAL const PDdissipationNth3B1 = PDdissipationNth3(&B1[index]);
    CCTK_REAL const PDdissipationNth1B2 = PDdissipationNth1(&B2[index]);
    CCTK_REAL const PDdissipationNth2B2 = PDdissipationNth2(&B2[index]);
    CCTK_REAL const PDdissipationNth3B2 = PDdissipationNth3(&B2[index]);
    CCTK_REAL const PDdissipationNth1B3 = PDdissipationNth1(&B3[index]);
    CCTK_REAL const PDdissipationNth2B3 = PDdissipationNth2(&B3[index]);
    CCTK_REAL const PDdissipationNth3B3 = PDdissipationNth3(&B3[index]);
    CCTK_REAL const PDdissipationNth1beta1 = PDdissipationNth1(&beta1[index]);
    CCTK_REAL const PDdissipationNth2beta1 = PDdissipationNth2(&beta1[index]);
    CCTK_REAL const PDdissipationNth3beta1 = PDdissipationNth3(&beta1[index]);
    CCTK_REAL const PDdissipationNth1beta2 = PDdissipationNth1(&beta2[index]);
    CCTK_REAL const PDdissipationNth2beta2 = PDdissipationNth2(&beta2[index]);
    CCTK_REAL const PDdissipationNth3beta2 = PDdissipationNth3(&beta2[index]);
    CCTK_REAL const PDdissipationNth1beta3 = PDdissipationNth1(&beta3[index]);
    CCTK_REAL const PDdissipationNth2beta3 = PDdissipationNth2(&beta3[index]);
    CCTK_REAL const PDdissipationNth3beta3 = PDdissipationNth3(&beta3[index]);
    CCTK_REAL const PDdissipationNth1gt11 = PDdissipationNth1(&gt11[index]);
    CCTK_REAL const PDdissipationNth2gt11 = PDdissipationNth2(&gt11[index]);
    CCTK_REAL const PDdissipationNth3gt11 = PDdissipationNth3(&gt11[index]);
    CCTK_REAL const PDdissipationNth1gt12 = PDdissipationNth1(&gt12[index]);
    CCTK_REAL const PDdissipationNth2gt12 = PDdissipationNth2(&gt12[index]);
    CCTK_REAL const PDdissipationNth3gt12 = PDdissipationNth3(&gt12[index]);
    CCTK_REAL const PDdissipationNth1gt13 = PDdissipationNth1(&gt13[index]);
    CCTK_REAL const PDdissipationNth2gt13 = PDdissipationNth2(&gt13[index]);
    CCTK_REAL const PDdissipationNth3gt13 = PDdissipationNth3(&gt13[index]);
    CCTK_REAL const PDdissipationNth1gt22 = PDdissipationNth1(&gt22[index]);
    CCTK_REAL const PDdissipationNth2gt22 = PDdissipationNth2(&gt22[index]);
    CCTK_REAL const PDdissipationNth3gt22 = PDdissipationNth3(&gt22[index]);
    CCTK_REAL const PDdissipationNth1gt23 = PDdissipationNth1(&gt23[index]);
    CCTK_REAL const PDdissipationNth2gt23 = PDdissipationNth2(&gt23[index]);
    CCTK_REAL const PDdissipationNth3gt23 = PDdissipationNth3(&gt23[index]);
    CCTK_REAL const PDdissipationNth1gt33 = PDdissipationNth1(&gt33[index]);
    CCTK_REAL const PDdissipationNth2gt33 = PDdissipationNth2(&gt33[index]);
    CCTK_REAL const PDdissipationNth3gt33 = PDdissipationNth3(&gt33[index]);
    CCTK_REAL const PDdissipationNth1phi = PDdissipationNth1(&phi[index]);
    CCTK_REAL const PDdissipationNth2phi = PDdissipationNth2(&phi[index]);
    CCTK_REAL const PDdissipationNth3phi = PDdissipationNth3(&phi[index]);
    CCTK_REAL const PDdissipationNth1trK = PDdissipationNth1(&trK[index]);
    CCTK_REAL const PDdissipationNth2trK = PDdissipationNth2(&trK[index]);
    CCTK_REAL const PDdissipationNth3trK = PDdissipationNth3(&trK[index]);
    CCTK_REAL const PDdissipationNth1Xt1 = PDdissipationNth1(&Xt1[index]);
    CCTK_REAL const PDdissipationNth2Xt1 = PDdissipationNth2(&Xt1[index]);
    CCTK_REAL const PDdissipationNth3Xt1 = PDdissipationNth3(&Xt1[index]);
    CCTK_REAL const PDdissipationNth1Xt2 = PDdissipationNth1(&Xt2[index]);
    CCTK_REAL const PDdissipationNth2Xt2 = PDdissipationNth2(&Xt2[index]);
    CCTK_REAL const PDdissipationNth3Xt2 = PDdissipationNth3(&Xt2[index]);
    CCTK_REAL const PDdissipationNth1Xt3 = PDdissipationNth1(&Xt3[index]);
    CCTK_REAL const PDdissipationNth2Xt3 = PDdissipationNth2(&Xt3[index]);
    CCTK_REAL const PDdissipationNth3Xt3 = PDdissipationNth3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL epsdiss1 = ToReal(EpsDiss);
    
    CCTK_REAL epsdiss2 = ToReal(EpsDiss);
    
    CCTK_REAL epsdiss3 = ToReal(EpsDiss);
    
    phirhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1phi + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2phi + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3phi + phirhsL;
    
    gt11rhsL = gt11rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1gt11 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2gt11 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3gt11;
    
    gt12rhsL = gt12rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1gt12 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2gt12 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3gt12;
    
    gt13rhsL = gt13rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1gt13 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2gt13 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3gt13;
    
    gt22rhsL = gt22rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1gt22 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2gt22 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3gt22;
    
    gt23rhsL = gt23rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1gt23 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2gt23 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3gt23;
    
    gt33rhsL = gt33rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1gt33 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2gt33 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3gt33;
    
    Xt1rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1Xt1 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2Xt1 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3Xt1 + Xt1rhsL;
    
    Xt2rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1Xt2 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2Xt2 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3Xt2 + Xt2rhsL;
    
    Xt3rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1Xt3 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2Xt3 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3Xt3 + Xt3rhsL;
    
    trKrhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1trK + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2trK + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3trK + trKrhsL;
    
    At11rhsL = At11rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At11 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At11 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At11;
    
    At12rhsL = At12rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At12 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At12 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At12;
    
    At13rhsL = At13rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At13 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At13 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At13;
    
    At22rhsL = At22rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At22 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At22 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At22;
    
    At23rhsL = At23rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At23 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At23 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At23;
    
    At33rhsL = At33rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At33 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At33 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At33;
    
    alpharhsL = alpharhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1alpha + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2alpha + (epsdiss1*J31L + 
      epsdiss2*J32L + epsdiss3*J33L)*PDdissipationNth3alpha;
    
    ArhsL = ArhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1A + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2A + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3A;
    
    beta1rhsL = beta1rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1beta1 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2beta1 + (epsdiss1*J31L + 
      epsdiss2*J32L + epsdiss3*J33L)*PDdissipationNth3beta1;
    
    beta2rhsL = beta2rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1beta2 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2beta2 + (epsdiss1*J31L + 
      epsdiss2*J32L + epsdiss3*J33L)*PDdissipationNth3beta2;
    
    beta3rhsL = beta3rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1beta3 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2beta3 + (epsdiss1*J31L + 
      epsdiss2*J32L + epsdiss3*J33L)*PDdissipationNth3beta3;
    
    B1rhsL = B1rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1B1 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2B1 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3B1;
    
    B2rhsL = B2rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1B2 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2B2 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3B2;
    
    B3rhsL = B3rhsL + (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1B3 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2B3 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3B3;
    
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
  LC_ENDLOOP3 (ML_BSSN_MP_Dissipation);
}

extern "C" void ML_BSSN_MP_Dissipation(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_Dissipation_Body);
}
