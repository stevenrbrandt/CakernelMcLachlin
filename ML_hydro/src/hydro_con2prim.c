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

void hydro_con2prim_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  /* Declare predefined quantities */
  // CCTK_REAL p1o2dx = INITVALUE;
  // CCTK_REAL p1o2dy = INITVALUE;
  // CCTK_REAL p1o2dz = INITVALUE;
  // CCTK_REAL p1o4dxdy = INITVALUE;
  // CCTK_REAL p1o4dxdz = INITVALUE;
  // CCTK_REAL p1o4dydz = INITVALUE;
  // CCTK_REAL p1odx2 = INITVALUE;
  // CCTK_REAL p1ody2 = INITVALUE;
  // CCTK_REAL p1odz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering hydro_con2prim_Body");
  }
  
  if (cctk_iteration % hydro_con2prim_calc_every != hydro_con2prim_calc_offset)
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
  LC_LOOP3 (hydro_con2prim,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    // int subblock_index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int const subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    
    /* Declare local copies of grid functions */
    // CCTK_REAL eneL = INITVALUE;
    // CCTK_REAL epsL = INITVALUE;
    // CCTK_REAL massL = INITVALUE;
    // CCTK_REAL mom1L = INITVALUE, mom2L = INITVALUE, mom3L = INITVALUE;
    // CCTK_REAL pressL = INITVALUE;
    // CCTK_REAL rhoL = INITVALUE;
    // CCTK_REAL vel1L = INITVALUE, vel2L = INITVALUE, vel3L = INITVALUE;
    // CCTK_REAL volL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    // CCTK_REAL PDstandardNth1vel1 = INITVALUE;
    // CCTK_REAL PDstandardNth2vel2 = INITVALUE;
    // CCTK_REAL PDstandardNth3vel3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL const eneL = ene[index];
    CCTK_REAL const massL = mass[index];
    CCTK_REAL const mom1L = mom1[index];
    CCTK_REAL const mom2L = mom2[index];
    CCTK_REAL const mom3L = mom3[index];
    CCTK_REAL const volL = vol[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    CCTK_REAL const PDstandardNth1vel1 = PDstandardNth1(vel1, i, j, k);
    CCTK_REAL const PDstandardNth2vel2 = PDstandardNth2(vel2, i, j, k);
    CCTK_REAL const PDstandardNth3vel3 = PDstandardNth3(vel3, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL const rhoL  =  massL*INV(volL);
    
    CCTK_REAL const vel1L  =  mom1L*INV(massL);
    
    CCTK_REAL const vel2L  =  mom2L*INV(massL);
    
    CCTK_REAL const vel3L  =  mom3L*INV(massL);
    
    CCTK_REAL const epsL  =  khalf*INV(massL)*(2*eneL - massL*(SQR(vel1L) + SQR(vel2L) + SQR(vel3L)));
    
    CCTK_REAL const pressL  =  alpha*(PDstandardNth1vel1 + PDstandardNth2vel2 + PDstandardNth3vel3) + epsL*Gamma*rhoL;
    
    
    /* Copy local copies back to grid functions */
    eps[index] = epsL;
    press[index] = pressL;
    rho[index] = rhoL;
    vel1[index] = vel1L;
    vel2[index] = vel2L;
    vel3[index] = vel3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (hydro_con2prim);
}

void hydro_con2prim(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &hydro_con2prim_Body);
}
