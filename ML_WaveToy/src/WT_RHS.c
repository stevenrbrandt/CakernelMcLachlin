/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (32-bit) (April 20, 2007) */

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

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

void WT_RHS_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  
  /* Declare the variables used for looping over grid points */
  CCTK_INT i = INITVALUE, j = INITVALUE, k = INITVALUE;
  CCTK_INT index = INITVALUE;
  CCTK_INT subblock_index = INITVALUE;
  
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WT_RHS_Body");
  }
  
  if (cctk_iteration % WT_RHS_calc_every != WT_RHS_calc_offset)
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
  for (k = min[2]; k < max[2]; k++)
  {
    for (j = min[1]; j < max[1]; j++)
    {
      for (i = min[0]; i < max[0]; i++)
      {
         index  =  CCTK_GFINDEX3D(cctkGH,i,j,k) ;
         subblock_index  =  i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2])) ;
        
        /* Declare shorthands */
        
        /* Declare local copies of grid functions */
        CCTK_REAL rhoL = INITVALUE, rhorhsL = INITVALUE;
        CCTK_REAL uL = INITVALUE, urhsL = INITVALUE;
        /* Declare precomputed derivatives*/
        
        /* Declare derivatives */
        CCTK_REAL PDstandardNth11u = INITVALUE;
        CCTK_REAL PDstandardNth22u = INITVALUE;
        CCTK_REAL PDstandardNth33u = INITVALUE;
        
        /* Assign local copies of grid functions */
        rhoL = rho[index];
        uL = u[index];
        
        /* Assign local copies of subblock grid functions */
        
        /* Include user supplied include files */
        
        /* Precompute derivatives (new style) */
        PDstandardNth11u = PDstandardNth11(u, i, j, k);
        PDstandardNth22u = PDstandardNth22(u, i, j, k);
        PDstandardNth33u = PDstandardNth33(u, i, j, k);
        
        /* Precompute derivatives (old style) */
        
        /* Calculate temporaries and grid functions */
        urhsL  =  rhoL;
        
        rhorhsL  =  PDstandardNth11u + PDstandardNth22u + PDstandardNth33u;
        
        
        /* Copy local copies back to grid functions */
        rhorhs[index] = rhorhsL;
        urhs[index] = urhsL;
        
        /* Copy local copies back to subblock grid functions */
      }
    }
  }
}

void WT_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &WT_RHS_Body);
}
