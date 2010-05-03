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
#include "Vectors.hh"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

static void hydro_con2prim_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
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
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    CCTK_REAL_VEC  eneL = vec_load(ene[index]);
    CCTK_REAL_VEC  epsL = vec_load(eps[index]);
    CCTK_REAL_VEC  massL = vec_load(mass[index]);
    CCTK_REAL_VEC  mom1L = vec_load(mom1[index]);
    CCTK_REAL_VEC  mom2L = vec_load(mom2[index]);
    CCTK_REAL_VEC  mom3L = vec_load(mom3[index]);
    CCTK_REAL_VEC  rhoL = vec_load(rho[index]);
    CCTK_REAL_VEC  vel1L = vec_load(vel1[index]);
    CCTK_REAL_VEC  vel2L = vec_load(vel2[index]);
    CCTK_REAL_VEC  vel3L = vec_load(vel3[index]);
    CCTK_REAL_VEC  volL = vec_load(vol[index]);
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    rhoL = massL*INV(volL);
    
    vel1L = mom1L*INV(massL);
    
    vel2L = mom2L*INV(massL);
    
    vel3L = mom3L*INV(massL);
    
    epsL = khalf*INV(massL)*(2*eneL - massL*(SQR(vel1L) + SQR(vel2L) + 
      SQR(vel3L)));
    
    CCTK_REAL_VEC pressL = epsL*Gamma*rhoL;
    
    
    /* Copy local copies back to grid functions */
    vec_store_nta(eps[index],epsL);
    vec_store_nta(press[index],pressL);
    vec_store_nta(rho[index],rhoL);
    vec_store_nta(vel1[index],vel1L);
    vec_store_nta(vel2[index],vel2L);
    vec_store_nta(vel3[index],vel3L);
  i += CCTK_REAL_VEC_SIZE-1;
  }
  LC_ENDLOOP3 (hydro_con2prim);
}

extern "C" void hydro_con2prim(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &hydro_con2prim_Body);
}
