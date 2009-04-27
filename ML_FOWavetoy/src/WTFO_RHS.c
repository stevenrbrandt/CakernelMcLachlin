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

void WTFO_RHS_Body(cGH const * const cctkGH, CCTK_INT const dir, CCTK_INT const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], CCTK_INT const min[3], CCTK_INT const max[3], CCTK_INT const n_subblock_gfs, CCTK_REAL * const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  CCTK_REAL dx = INITVALUE, dy = INITVALUE, dz = INITVALUE;
  CCTK_REAL dxi = INITVALUE, dyi = INITVALUE, dzi = INITVALUE;
  CCTK_REAL khalf = INITVALUE, kthird = INITVALUE, ktwothird = INITVALUE, kfourthird = INITVALUE, keightthird = INITVALUE;
  CCTK_REAL hdxi = INITVALUE, hdyi = INITVALUE, hdzi = INITVALUE;
  
  
  /* Declare predefined quantities */
  CCTK_REAL p1o12dx = INITVALUE;
  CCTK_REAL p1o12dy = INITVALUE;
  CCTK_REAL p1o12dz = INITVALUE;
  CCTK_REAL p1o144dx2dy = INITVALUE;
  CCTK_REAL p1o144dx2dz = INITVALUE;
  CCTK_REAL p1o144dxdy = INITVALUE;
  CCTK_REAL p1o144dxdy2 = INITVALUE;
  CCTK_REAL p1o144dxdz = INITVALUE;
  CCTK_REAL p1o144dxdz2 = INITVALUE;
  CCTK_REAL p1o144dy2dz = INITVALUE;
  CCTK_REAL p1o144dydz = INITVALUE;
  CCTK_REAL p1o144dydz2 = INITVALUE;
  CCTK_REAL p1o1728dxdydz = INITVALUE;
  CCTK_REAL p1o2dx3 = INITVALUE;
  CCTK_REAL p1o2dy3 = INITVALUE;
  CCTK_REAL p1o2dz3 = INITVALUE;
  CCTK_REAL pm1o12dx2 = INITVALUE;
  CCTK_REAL pm1o12dy2 = INITVALUE;
  CCTK_REAL pm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WTFO_RHS_Body");
  }
  
  if (cctk_iteration % WTFO_RHS_calc_every != WTFO_RHS_calc_offset)
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
  p1o12dx = INV(dx)/12.;
  p1o12dy = INV(dy)/12.;
  p1o12dz = INV(dz)/12.;
  p1o144dx2dy = (INV(dy)*pow(dx,-2))/144.;
  p1o144dx2dz = (INV(dz)*pow(dx,-2))/144.;
  p1o144dxdy = (INV(dx)*INV(dy))/144.;
  p1o144dxdy2 = (INV(dx)*pow(dy,-2))/144.;
  p1o144dxdz = (INV(dx)*INV(dz))/144.;
  p1o144dxdz2 = (INV(dx)*pow(dz,-2))/144.;
  p1o144dy2dz = (INV(dz)*pow(dy,-2))/144.;
  p1o144dydz = (INV(dy)*INV(dz))/144.;
  p1o144dydz2 = (INV(dy)*pow(dz,-2))/144.;
  p1o1728dxdydz = (INV(dx)*INV(dy)*INV(dz))/1728.;
  p1o2dx3 = khalf*pow(dx,-3);
  p1o2dy3 = khalf*pow(dy,-3);
  p1o2dz3 = khalf*pow(dz,-3);
  pm1o12dx2 = -pow(dx,-2)/12.;
  pm1o12dy2 = -pow(dy,-2)/12.;
  pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (WTFO_RHS,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lssh[CCTK_LSSH_IDX(0,0)],cctk_lssh[CCTK_LSSH_IDX(0,1)],cctk_lssh[CCTK_LSSH_IDX(0,2)])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    
    /* Declare local copies of grid functions */
    CCTK_REAL rhoL = INITVALUE, rhorhsL = INITVALUE;
    CCTK_REAL urhsL = INITVALUE;
    CCTK_REAL v1L = INITVALUE, v1rhsL = INITVALUE, v2L = INITVALUE, v2rhsL = INITVALUE, v3L = INITVALUE, v3rhsL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandardNth1rho = INITVALUE;
    CCTK_REAL PDstandardNth2rho = INITVALUE;
    CCTK_REAL PDstandardNth3rho = INITVALUE;
    CCTK_REAL PDstandardNth1v1 = INITVALUE;
    CCTK_REAL PDstandardNth2v2 = INITVALUE;
    CCTK_REAL PDstandardNth3v3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    rhoL = rho[index];
    v1L = v1[index];
    v2L = v2[index];
    v3L = v3[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandardNth1rho = PDstandardNth1(rho, i, j, k);
    PDstandardNth2rho = PDstandardNth2(rho, i, j, k);
    PDstandardNth3rho = PDstandardNth3(rho, i, j, k);
    PDstandardNth1v1 = PDstandardNth1(v1, i, j, k);
    PDstandardNth2v2 = PDstandardNth2(v2, i, j, k);
    PDstandardNth3v3 = PDstandardNth3(v3, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    urhsL  =  rhoL;
    
    rhorhsL  =  PDstandardNth1v1 + PDstandardNth2v2 + PDstandardNth3v3;
    
    v1rhsL  =  PDstandardNth1rho;
    
    v2rhsL  =  PDstandardNth2rho;
    
    v3rhsL  =  PDstandardNth3rho;
    
    
    /* Copy local copies back to grid functions */
    rhorhs[index] = rhorhsL;
    urhs[index] = urhsL;
    v1rhs[index] = v1rhsL;
    v2rhs[index] = v2rhsL;
    v3rhs[index] = v3rhsL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (WTFO_RHS);
}

void WTFO_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &WTFO_RHS_Body);
}
