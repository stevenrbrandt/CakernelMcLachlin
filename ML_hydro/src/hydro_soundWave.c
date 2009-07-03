/*  File produced by user eschnett */
/*  Produced with Mathematica Version 7.0 for Mac OS X x86 (64-bit) (February 19, 2009) */

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

void hydro_soundWave_Body(cGH const * const cctkGH, CCTK_INT const dir, CCTK_INT const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], CCTK_INT const min[3], CCTK_INT const max[3], CCTK_INT const n_subblock_gfs, CCTK_REAL * const subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering hydro_soundWave_Body");
  }
  
  if (cctk_iteration % hydro_soundWave_calc_every != hydro_soundWave_calc_offset)
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
  LC_LOOP3 (hydro_soundWave,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lssh[CCTK_LSSH_IDX(0,0)],cctk_lssh[CCTK_LSSH_IDX(0,1)],cctk_lssh[CCTK_LSSH_IDX(0,2)])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    
    /* Declare local copies of grid functions */
    CCTK_REAL epsL = INITVALUE;
    CCTK_REAL rhoL = INITVALUE;
    CCTK_REAL vel1L = INITVALUE, vel2L = INITVALUE, vel3L = INITVALUE;
    CCTK_REAL xL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    xL = x[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    rhoL  =  1.;
    
    vel1L  =  A*Sin(2*Pi*xL*INV(L));
    
    vel2L  =  A*Sin(2*Pi*xL*INV(L));
    
    vel3L  =  A*Sin(2*Pi*xL*INV(L));
    
    epsL  =  1.;
    
    
    /* Copy local copies back to grid functions */
    eps[index] = epsL;
    rho[index] = rhoL;
    vel1[index] = vel1L;
    vel2[index] = vel2L;
    vel3[index] = vel3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (hydro_soundWave);
}

void hydro_soundWave(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &hydro_soundWave_Body);
}
