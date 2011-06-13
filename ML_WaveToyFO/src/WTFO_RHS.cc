/*  File produced by Kranc */

#define KRANC_C

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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

extern "C" void WTFO_RHS_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_WaveToyFO::WT_rhorhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_WaveToyFO::WT_rhorhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_WaveToyFO::WT_urhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_WaveToyFO::WT_urhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_WaveToyFO::WT_vrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_WaveToyFO::WT_vrhs.");
  return;
}

static void WTFO_RHS_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WTFO_RHS_Body");
  }
  
  if (cctk_iteration % WTFO_RHS_calc_every != WTFO_RHS_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_WaveToyFO::WT_rho","ML_WaveToyFO::WT_rhorhs","ML_WaveToyFO::WT_urhs","ML_WaveToyFO::WT_v","ML_WaveToyFO::WT_vrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "WTFO_RHS", 5, groups);
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  ptrdiff_t const di = 1;
  ptrdiff_t const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  CCTK_REAL const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL const dz = ToReal(CCTK_DELTA_SPACE(2));
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
  CCTK_REAL const pm1o12dx2 = -0.0833333333333333333333333333333*INV(SQR(dx));
  CCTK_REAL const pm1o12dy2 = -0.0833333333333333333333333333333*INV(SQR(dy));
  CCTK_REAL const pm1o12dz2 = -0.0833333333333333333333333333333*INV(SQR(dz));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (WTFO_RHS,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    CCTK_REAL rhoL = rho[index];
    CCTK_REAL v1L = v1[index];
    CCTK_REAL v2L = v2[index];
    CCTK_REAL v3L = v3[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDstandardNth1rho = PDstandardNth1(&rho[index]);
    CCTK_REAL const PDstandardNth2rho = PDstandardNth2(&rho[index]);
    CCTK_REAL const PDstandardNth3rho = PDstandardNth3(&rho[index]);
    CCTK_REAL const PDstandardNth1v1 = PDstandardNth1(&v1[index]);
    CCTK_REAL const PDstandardNth2v2 = PDstandardNth2(&v2[index]);
    CCTK_REAL const PDstandardNth3v3 = PDstandardNth3(&v3[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL urhsL = rhoL;
    
    CCTK_REAL rhorhsL = PDstandardNth1v1 + PDstandardNth2v2 + 
      PDstandardNth3v3;
    
    CCTK_REAL v1rhsL = PDstandardNth1rho;
    
    CCTK_REAL v2rhsL = PDstandardNth2rho;
    
    CCTK_REAL v3rhsL = PDstandardNth3rho;
    
    
    /* Copy local copies back to grid functions */
    rhorhs[index] = rhorhsL;
    urhs[index] = urhsL;
    v1rhs[index] = v1rhsL;
    v2rhs[index] = v2rhsL;
    v3rhs[index] = v3rhsL;
  }
  LC_ENDLOOP3 (WTFO_RHS);
}

extern "C" void WTFO_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &WTFO_RHS_Body);
}
