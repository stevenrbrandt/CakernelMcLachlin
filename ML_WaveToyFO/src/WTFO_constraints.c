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
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

void WTFO_constraints_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_WaveToyFO::WT_w","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_WaveToyFO::WT_w.");
  return;
}

void WTFO_constraints_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WTFO_constraints_Body");
  }
  
  if (cctk_iteration % WTFO_constraints_calc_every != WTFO_constraints_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_WaveToyFO::WT_v","ML_WaveToyFO::WT_w"};
  GenericFD_AssertGroupStorage(cctkGH, "WTFO_constraints", 2, groups);
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  CCTK_REAL const dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL const dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL const dz = CCTK_DELTA_SPACE(2);
  int const di = 1;
  int const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  int const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  CCTK_REAL const dxi = 1.0 / dx;
  CCTK_REAL const dyi = 1.0 / dy;
  CCTK_REAL const dzi = 1.0 / dz;
  CCTK_REAL const khalf = 0.5;
  CCTK_REAL const kthird = 1/3.0;
  CCTK_REAL const ktwothird = 2.0/3.0;
  CCTK_REAL const kfourthird = 4.0/3.0;
  CCTK_REAL const keightthird = 8.0/3.0;
  CCTK_REAL const hdxi = 0.5 * dxi;
  CCTK_REAL const hdyi = 0.5 * dyi;
  CCTK_REAL const hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  CCTK_REAL const p1o12dx = INV(dx)/12.;
  CCTK_REAL const p1o12dy = INV(dy)/12.;
  CCTK_REAL const p1o12dz = INV(dz)/12.;
  CCTK_REAL const p1o144dxdy = (INV(dx)*INV(dy))/144.;
  CCTK_REAL const p1o144dxdz = (INV(dx)*INV(dz))/144.;
  CCTK_REAL const p1o144dydz = (INV(dy)*INV(dz))/144.;
  CCTK_REAL const pm1o12dx2 = -pow(dx,-2)/12.;
  CCTK_REAL const pm1o12dy2 = -pow(dy,-2)/12.;
  CCTK_REAL const pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (WTFO_constraints,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL PDstandardNth2v1 = INITVALUE;
    // CCTK_REAL PDstandardNth3v1 = INITVALUE;
    // CCTK_REAL PDstandardNth1v2 = INITVALUE;
    // CCTK_REAL PDstandardNth3v2 = INITVALUE;
    // CCTK_REAL PDstandardNth1v3 = INITVALUE;
    // CCTK_REAL PDstandardNth2v3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL  v1L = v1[index];
    CCTK_REAL  v2L = v2[index];
    CCTK_REAL  v3L = v3[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDstandardNth2v1 = PDstandardNth2(v1, i, j, k);
    CCTK_REAL const PDstandardNth3v1 = PDstandardNth3(v1, i, j, k);
    CCTK_REAL const PDstandardNth1v2 = PDstandardNth1(v2, i, j, k);
    CCTK_REAL const PDstandardNth3v2 = PDstandardNth3(v2, i, j, k);
    CCTK_REAL const PDstandardNth1v3 = PDstandardNth1(v3, i, j, k);
    CCTK_REAL const PDstandardNth2v3 = PDstandardNth2(v3, i, j, k);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL w1L = -PDstandardNth2v3 + PDstandardNth3v2;
    
    CCTK_REAL w2L = PDstandardNth1v3 - PDstandardNth3v1;
    
    CCTK_REAL w3L = -PDstandardNth1v2 + PDstandardNth2v1;
    
    
    /* Copy local copies back to grid functions */
    w1[index] = w1L;
    w2[index] = w2L;
    w3[index] = w3L;
  }
  LC_ENDLOOP3 (WTFO_constraints);
}

void WTFO_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &WTFO_constraints_Body);
}
