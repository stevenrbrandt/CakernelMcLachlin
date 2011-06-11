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
  LC_LOOP3 (ML_BSSN_UPW_Advect,
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
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    phirhsL = PDupwindNth1(&phi[index])*beta1L + 
      PDupwindNth2(&phi[index])*beta2L + PDupwindNth3(&phi[index])*beta3L 
      + phirhsL;
    
    gt11rhsL = PDupwindNth1(&gt11[index])*beta1L + 
      PDupwindNth2(&gt11[index])*beta2L + 
      PDupwindNth3(&gt11[index])*beta3L + gt11rhsL;
    
    gt12rhsL = PDupwindNth1(&gt12[index])*beta1L + 
      PDupwindNth2(&gt12[index])*beta2L + 
      PDupwindNth3(&gt12[index])*beta3L + gt12rhsL;
    
    gt13rhsL = PDupwindNth1(&gt13[index])*beta1L + 
      PDupwindNth2(&gt13[index])*beta2L + 
      PDupwindNth3(&gt13[index])*beta3L + gt13rhsL;
    
    gt22rhsL = PDupwindNth1(&gt22[index])*beta1L + 
      PDupwindNth2(&gt22[index])*beta2L + 
      PDupwindNth3(&gt22[index])*beta3L + gt22rhsL;
    
    gt23rhsL = PDupwindNth1(&gt23[index])*beta1L + 
      PDupwindNth2(&gt23[index])*beta2L + 
      PDupwindNth3(&gt23[index])*beta3L + gt23rhsL;
    
    gt33rhsL = PDupwindNth1(&gt33[index])*beta1L + 
      PDupwindNth2(&gt33[index])*beta2L + 
      PDupwindNth3(&gt33[index])*beta3L + gt33rhsL;
    
    Xt1rhsL = PDupwindNth1(&Xt1[index])*beta1L + 
      PDupwindNth2(&Xt1[index])*beta2L + PDupwindNth3(&Xt1[index])*beta3L 
      + Xt1rhsL;
    
    Xt2rhsL = PDupwindNth1(&Xt2[index])*beta1L + 
      PDupwindNth2(&Xt2[index])*beta2L + PDupwindNth3(&Xt2[index])*beta3L 
      + Xt2rhsL;
    
    Xt3rhsL = PDupwindNth1(&Xt3[index])*beta1L + 
      PDupwindNth2(&Xt3[index])*beta2L + PDupwindNth3(&Xt3[index])*beta3L 
      + Xt3rhsL;
    
    trKrhsL = PDupwindNth1(&trK[index])*beta1L + 
      PDupwindNth2(&trK[index])*beta2L + PDupwindNth3(&trK[index])*beta3L 
      + trKrhsL;
    
    At11rhsL = At11rhsL + PDupwindNth1(&At11[index])*beta1L + 
      PDupwindNth2(&At11[index])*beta2L + 
      PDupwindNth3(&At11[index])*beta3L;
    
    At12rhsL = At12rhsL + PDupwindNth1(&At12[index])*beta1L + 
      PDupwindNth2(&At12[index])*beta2L + 
      PDupwindNth3(&At12[index])*beta3L;
    
    At13rhsL = At13rhsL + PDupwindNth1(&At13[index])*beta1L + 
      PDupwindNth2(&At13[index])*beta2L + 
      PDupwindNth3(&At13[index])*beta3L;
    
    At22rhsL = At22rhsL + PDupwindNth1(&At22[index])*beta1L + 
      PDupwindNth2(&At22[index])*beta2L + 
      PDupwindNth3(&At22[index])*beta3L;
    
    At23rhsL = At23rhsL + PDupwindNth1(&At23[index])*beta1L + 
      PDupwindNth2(&At23[index])*beta2L + 
      PDupwindNth3(&At23[index])*beta3L;
    
    At33rhsL = At33rhsL + PDupwindNth1(&At33[index])*beta1L + 
      PDupwindNth2(&At33[index])*beta2L + 
      PDupwindNth3(&At33[index])*beta3L;
    
    alpharhsL = alpharhsL + (PDupwindNth1(&alpha[index])*beta1L + 
      PDupwindNth2(&alpha[index])*beta2L + 
      PDupwindNth3(&alpha[index])*beta3L)*ToReal(LapseAdvectionCoeff);
    
    ArhsL = ArhsL + (PDupwindNth1(&A[index])*beta1L + 
      PDupwindNth2(&A[index])*beta2L + 
      PDupwindNth3(&A[index])*beta3L)*ToReal(LapseAdvectionCoeff);
    
    beta1rhsL = beta1rhsL + (PDupwindNth1(&beta1[index])*beta1L + 
      PDupwindNth2(&beta1[index])*beta2L + 
      PDupwindNth3(&beta1[index])*beta3L)*ToReal(ShiftAdvectionCoeff);
    
    beta2rhsL = beta2rhsL + (PDupwindNth1(&beta2[index])*beta1L + 
      PDupwindNth2(&beta2[index])*beta2L + 
      PDupwindNth3(&beta2[index])*beta3L)*ToReal(ShiftAdvectionCoeff);
    
    beta3rhsL = beta3rhsL + (PDupwindNth1(&beta3[index])*beta1L + 
      PDupwindNth2(&beta3[index])*beta2L + 
      PDupwindNth3(&beta3[index])*beta3L)*ToReal(ShiftAdvectionCoeff);
    
    B1rhsL = B1rhsL + beta1L*((PDupwindNth1(&B1[index]) - 
      PDupwindNth1(&Xt1[index]))*ToReal(ShiftAdvectionCoeff) + 
      PDupwindNth1(&Xt1[index])*ToReal(ShiftBCoeff)) + 
      beta2L*((PDupwindNth2(&B1[index]) - 
      PDupwindNth2(&Xt1[index]))*ToReal(ShiftAdvectionCoeff) + 
      PDupwindNth2(&Xt1[index])*ToReal(ShiftBCoeff)) + 
      beta3L*((PDupwindNth3(&B1[index]) - 
      PDupwindNth3(&Xt1[index]))*ToReal(ShiftAdvectionCoeff) + 
      PDupwindNth3(&Xt1[index])*ToReal(ShiftBCoeff));
    
    B2rhsL = B2rhsL + beta1L*((PDupwindNth1(&B2[index]) - 
      PDupwindNth1(&Xt2[index]))*ToReal(ShiftAdvectionCoeff) + 
      PDupwindNth1(&Xt2[index])*ToReal(ShiftBCoeff)) + 
      beta2L*((PDupwindNth2(&B2[index]) - 
      PDupwindNth2(&Xt2[index]))*ToReal(ShiftAdvectionCoeff) + 
      PDupwindNth2(&Xt2[index])*ToReal(ShiftBCoeff)) + 
      beta3L*((PDupwindNth3(&B2[index]) - 
      PDupwindNth3(&Xt2[index]))*ToReal(ShiftAdvectionCoeff) + 
      PDupwindNth3(&Xt2[index])*ToReal(ShiftBCoeff));
    
    B3rhsL = B3rhsL + beta1L*((PDupwindNth1(&B3[index]) - 
      PDupwindNth1(&Xt3[index]))*ToReal(ShiftAdvectionCoeff) + 
      PDupwindNth1(&Xt3[index])*ToReal(ShiftBCoeff)) + 
      beta2L*((PDupwindNth2(&B3[index]) - 
      PDupwindNth2(&Xt3[index]))*ToReal(ShiftAdvectionCoeff) + 
      PDupwindNth2(&Xt3[index])*ToReal(ShiftBCoeff)) + 
      beta3L*((PDupwindNth3(&B3[index]) - 
      PDupwindNth3(&Xt3[index]))*ToReal(ShiftAdvectionCoeff) + 
      PDupwindNth3(&Xt3[index])*ToReal(ShiftBCoeff));
    
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
  LC_ENDLOOP3 (ML_BSSN_UPW_Advect);
}

extern "C" void ML_BSSN_UPW_Advect(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_UPW_Advect_Body);
}
