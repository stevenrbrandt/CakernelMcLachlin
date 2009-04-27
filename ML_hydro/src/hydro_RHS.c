/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (64-bit) (May 21, 2008) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

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

void hydro_RHS_Body(cGH const * const cctkGH, CCTK_INT const dir, CCTK_INT const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], CCTK_INT const min[3], CCTK_INT const max[3], CCTK_INT const n_subblock_gfs, CCTK_REAL * const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  CCTK_REAL dx = INITVALUE, dy = INITVALUE, dz = INITVALUE;
  CCTK_REAL dxi = INITVALUE, dyi = INITVALUE, dzi = INITVALUE;
  CCTK_REAL khalf = INITVALUE, kthird = INITVALUE, ktwothird = INITVALUE, kfourthird = INITVALUE, keightthird = INITVALUE;
  CCTK_REAL hdxi = INITVALUE, hdyi = INITVALUE, hdzi = INITVALUE;
  
  
  /* Declare predefined quantities */
  CCTK_REAL p1o2dx = INITVALUE;
  CCTK_REAL p1o2dy = INITVALUE;
  CCTK_REAL p1o2dz = INITVALUE;
  CCTK_REAL p1o4dxdy = INITVALUE;
  CCTK_REAL p1o4dxdz = INITVALUE;
  CCTK_REAL p1o4dydz = INITVALUE;
  CCTK_REAL p1odx2 = INITVALUE;
  CCTK_REAL p1ody2 = INITVALUE;
  CCTK_REAL p1odz2 = INITVALUE;
  
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
  dx = CCTK_DELTA_SPACE(0);
  dy = CCTK_DELTA_SPACE(1);
  dz = CCTK_DELTA_SPACE(2);
  dxi = 1.0 / dx;
  dyi = 1.0 / dy;
  dzi = 1.0 / dz;
  khalf = 0.5;
  kthird = 1/3.0;
  ktwothird = 2.0/3.0;
  kfourthird = 4.0/3.0;
  keightthird = 8.0/3.0;
  hdxi = 0.5 * dxi;
  hdyi = 0.5 * dyi;
  hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  p1o2dx = khalf*INV(dx);
  p1o2dy = khalf*INV(dy);
  p1o2dz = khalf*INV(dz);
  p1o4dxdy = (INV(dx)*INV(dy))/4.;
  p1o4dxdz = (INV(dx)*INV(dz))/4.;
  p1o4dydz = (INV(dy)*INV(dz))/4.;
  p1odx2 = pow(dx,-2);
  p1ody2 = pow(dy,-2);
  p1odz2 = pow(dz,-2);
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (hydro_RHS,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lssh[CCTK_LSSH_IDX(0,0)],cctk_lssh[CCTK_LSSH_IDX(0,1)],cctk_lssh[CCTK_LSSH_IDX(0,2)])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    
    /* Declare local copies of grid functions */
    CCTK_REAL eneflux1L = INITVALUE, eneflux2L = INITVALUE, eneflux3L = INITVALUE;
    CCTK_REAL enerhsL = INITVALUE;
    CCTK_REAL massflux1L = INITVALUE, massflux2L = INITVALUE, massflux3L = INITVALUE;
    CCTK_REAL massrhsL = INITVALUE;
    CCTK_REAL mom1rhsL = INITVALUE, mom2rhsL = INITVALUE, mom3rhsL = INITVALUE;
    CCTK_REAL momflux11L = INITVALUE, momflux12L = INITVALUE, momflux13L = INITVALUE, momflux21L = INITVALUE, momflux22L = INITVALUE, momflux23L = INITVALUE;
    CCTK_REAL momflux31L = INITVALUE, momflux32L = INITVALUE, momflux33L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandardNth1eneflux1 = INITVALUE;
    CCTK_REAL PDstandardNth2eneflux2 = INITVALUE;
    CCTK_REAL PDstandardNth3eneflux3 = INITVALUE;
    CCTK_REAL PDstandardNth1massflux1 = INITVALUE;
    CCTK_REAL PDstandardNth2massflux2 = INITVALUE;
    CCTK_REAL PDstandardNth3massflux3 = INITVALUE;
    CCTK_REAL PDstandardNth1momflux11 = INITVALUE;
    CCTK_REAL PDstandardNth2momflux12 = INITVALUE;
    CCTK_REAL PDstandardNth3momflux13 = INITVALUE;
    CCTK_REAL PDstandardNth1momflux21 = INITVALUE;
    CCTK_REAL PDstandardNth2momflux22 = INITVALUE;
    CCTK_REAL PDstandardNth3momflux23 = INITVALUE;
    CCTK_REAL PDstandardNth1momflux31 = INITVALUE;
    CCTK_REAL PDstandardNth2momflux32 = INITVALUE;
    CCTK_REAL PDstandardNth3momflux33 = INITVALUE;
    
    /* Assign local copies of grid functions */
    eneflux1L = eneflux1[index];
    eneflux2L = eneflux2[index];
    eneflux3L = eneflux3[index];
    massflux1L = massflux1[index];
    massflux2L = massflux2[index];
    massflux3L = massflux3[index];
    momflux11L = momflux11[index];
    momflux12L = momflux12[index];
    momflux13L = momflux13[index];
    momflux21L = momflux21[index];
    momflux22L = momflux22[index];
    momflux23L = momflux23[index];
    momflux31L = momflux31[index];
    momflux32L = momflux32[index];
    momflux33L = momflux33[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandardNth1eneflux1 = PDstandardNth1(eneflux1, i, j, k);
    PDstandardNth2eneflux2 = PDstandardNth2(eneflux2, i, j, k);
    PDstandardNth3eneflux3 = PDstandardNth3(eneflux3, i, j, k);
    PDstandardNth1massflux1 = PDstandardNth1(massflux1, i, j, k);
    PDstandardNth2massflux2 = PDstandardNth2(massflux2, i, j, k);
    PDstandardNth3massflux3 = PDstandardNth3(massflux3, i, j, k);
    PDstandardNth1momflux11 = PDstandardNth1(momflux11, i, j, k);
    PDstandardNth2momflux12 = PDstandardNth2(momflux12, i, j, k);
    PDstandardNth3momflux13 = PDstandardNth3(momflux13, i, j, k);
    PDstandardNth1momflux21 = PDstandardNth1(momflux21, i, j, k);
    PDstandardNth2momflux22 = PDstandardNth2(momflux22, i, j, k);
    PDstandardNth3momflux23 = PDstandardNth3(momflux23, i, j, k);
    PDstandardNth1momflux31 = PDstandardNth1(momflux31, i, j, k);
    PDstandardNth2momflux32 = PDstandardNth2(momflux32, i, j, k);
    PDstandardNth3momflux33 = PDstandardNth3(momflux33, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    massrhsL  =  -PDstandardNth1massflux1 - PDstandardNth2massflux2 - 
        PDstandardNth3massflux3;
    
    mom1rhsL  =  -PDstandardNth1momflux11 - PDstandardNth2momflux12 - 
        PDstandardNth3momflux13;
    
    mom2rhsL  =  -PDstandardNth1momflux21 - PDstandardNth2momflux22 - 
        PDstandardNth3momflux23;
    
    mom3rhsL  =  -PDstandardNth1momflux31 - PDstandardNth2momflux32 - 
        PDstandardNth3momflux33;
    
    enerhsL  =  -PDstandardNth1eneflux1 - PDstandardNth2eneflux2 - 
        PDstandardNth3eneflux3;
    
    
    /* Copy local copies back to grid functions */
    enerhs[index] = enerhsL;
    massrhs[index] = massrhsL;
    mom1rhs[index] = mom1rhsL;
    mom2rhs[index] = mom2rhsL;
    mom3rhs[index] = mom3rhsL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (hydro_RHS);
}

void hydro_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &hydro_RHS_Body);
}
