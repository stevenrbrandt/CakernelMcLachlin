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

void hydro_RHS_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering hydro_RHS_Body");
  }
  
  if (cctk_iteration % hydro_RHS_calc_every != hydro_RHS_calc_offset)
  {
    return;
  }
  
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
  CCTK_REAL const p1o2dx = khalf*INV(dx);
  CCTK_REAL const p1o2dy = khalf*INV(dy);
  CCTK_REAL const p1o2dz = khalf*INV(dz);
  CCTK_REAL const p1o4dxdy = (INV(dx)*INV(dy))/4.;
  CCTK_REAL const p1o4dxdz = (INV(dx)*INV(dz))/4.;
  CCTK_REAL const p1o4dydz = (INV(dy)*INV(dz))/4.;
  CCTK_REAL const p1odx2 = pow(dx,-2);
  CCTK_REAL const p1ody2 = pow(dy,-2);
  CCTK_REAL const p1odz2 = pow(dz,-2);
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (hydro_RHS,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL PDstandardNth1eneflux1 = INITVALUE;
    // CCTK_REAL PDstandardNth2eneflux2 = INITVALUE;
    // CCTK_REAL PDstandardNth3eneflux3 = INITVALUE;
    // CCTK_REAL PDstandardNth1massflux1 = INITVALUE;
    // CCTK_REAL PDstandardNth2massflux2 = INITVALUE;
    // CCTK_REAL PDstandardNth3massflux3 = INITVALUE;
    // CCTK_REAL PDstandardNth1momflux11 = INITVALUE;
    // CCTK_REAL PDstandardNth2momflux12 = INITVALUE;
    // CCTK_REAL PDstandardNth3momflux13 = INITVALUE;
    // CCTK_REAL PDstandardNth1momflux21 = INITVALUE;
    // CCTK_REAL PDstandardNth2momflux22 = INITVALUE;
    // CCTK_REAL PDstandardNth3momflux23 = INITVALUE;
    // CCTK_REAL PDstandardNth1momflux31 = INITVALUE;
    // CCTK_REAL PDstandardNth2momflux32 = INITVALUE;
    // CCTK_REAL PDstandardNth3momflux33 = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL  eneflux1L = eneflux1[index];
    CCTK_REAL  eneflux2L = eneflux2[index];
    CCTK_REAL  eneflux3L = eneflux3[index];
    CCTK_REAL  massflux1L = massflux1[index];
    CCTK_REAL  massflux2L = massflux2[index];
    CCTK_REAL  massflux3L = massflux3[index];
    CCTK_REAL  momflux11L = momflux11[index];
    CCTK_REAL  momflux12L = momflux12[index];
    CCTK_REAL  momflux13L = momflux13[index];
    CCTK_REAL  momflux21L = momflux21[index];
    CCTK_REAL  momflux22L = momflux22[index];
    CCTK_REAL  momflux23L = momflux23[index];
    CCTK_REAL  momflux31L = momflux31[index];
    CCTK_REAL  momflux32L = momflux32[index];
    CCTK_REAL  momflux33L = momflux33[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDstandardNth1eneflux1 = PDstandardNth1(eneflux1, i, j, k);
    CCTK_REAL const PDstandardNth2eneflux2 = PDstandardNth2(eneflux2, i, j, k);
    CCTK_REAL const PDstandardNth3eneflux3 = PDstandardNth3(eneflux3, i, j, k);
    CCTK_REAL const PDstandardNth1massflux1 = PDstandardNth1(massflux1, i, j, k);
    CCTK_REAL const PDstandardNth2massflux2 = PDstandardNth2(massflux2, i, j, k);
    CCTK_REAL const PDstandardNth3massflux3 = PDstandardNth3(massflux3, i, j, k);
    CCTK_REAL const PDstandardNth1momflux11 = PDstandardNth1(momflux11, i, j, k);
    CCTK_REAL const PDstandardNth2momflux12 = PDstandardNth2(momflux12, i, j, k);
    CCTK_REAL const PDstandardNth3momflux13 = PDstandardNth3(momflux13, i, j, k);
    CCTK_REAL const PDstandardNth1momflux21 = PDstandardNth1(momflux21, i, j, k);
    CCTK_REAL const PDstandardNth2momflux22 = PDstandardNth2(momflux22, i, j, k);
    CCTK_REAL const PDstandardNth3momflux23 = PDstandardNth3(momflux23, i, j, k);
    CCTK_REAL const PDstandardNth1momflux31 = PDstandardNth1(momflux31, i, j, k);
    CCTK_REAL const PDstandardNth2momflux32 = PDstandardNth2(momflux32, i, j, k);
    CCTK_REAL const PDstandardNth3momflux33 = PDstandardNth3(momflux33, i, j, k);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL massrhsL = -PDstandardNth1massflux1 - 
      PDstandardNth2massflux2 - PDstandardNth3massflux3;
    
    CCTK_REAL mom1rhsL = -PDstandardNth1momflux11 - 
      PDstandardNth2momflux12 - PDstandardNth3momflux13;
    
    CCTK_REAL mom2rhsL = -PDstandardNth1momflux21 - 
      PDstandardNth2momflux22 - PDstandardNth3momflux23;
    
    CCTK_REAL mom3rhsL = -PDstandardNth1momflux31 - 
      PDstandardNth2momflux32 - PDstandardNth3momflux33;
    
    CCTK_REAL enerhsL = -PDstandardNth1eneflux1 - PDstandardNth2eneflux2 - 
      PDstandardNth3eneflux3;
    
    
    /* Copy local copies back to grid functions */
    enerhs[index] = enerhsL;
    massrhs[index] = massrhsL;
    mom1rhs[index] = mom1rhsL;
    mom2rhs[index] = mom2rhsL;
    mom3rhs[index] = mom3rhsL;
  }
  LC_ENDLOOP3 (hydro_RHS);
}

void hydro_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &hydro_RHS_Body);
}
