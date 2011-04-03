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
  
  const char *groups[] = {"ML_hydro::ene_group","ML_hydro::eps_group","ML_hydro::mass_group","ML_hydro::mom_group","ML_hydro::press_group","ML_hydro::rho_group","ML_hydro::vel_group","ML_hydro::vol_group"};
  GenericFD_AssertGroupStorage(cctkGH, "hydro_con2prim", 8, groups);
  
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
  CCTK_REAL const p1o2dx = 0.5*INV(dx);
  CCTK_REAL const p1o2dy = 0.5*INV(dy);
  CCTK_REAL const p1o2dz = 0.5*INV(dz);
  CCTK_REAL const p1o4dxdy = 0.25*INV(dx)*INV(dy);
  CCTK_REAL const p1o4dxdz = 0.25*INV(dx)*INV(dz);
  CCTK_REAL const p1o4dydz = 0.25*INV(dy)*INV(dz);
  CCTK_REAL const p1odx2 = INV(SQR(dx));
  CCTK_REAL const p1ody2 = INV(SQR(dy));
  CCTK_REAL const p1odz2 = INV(SQR(dz));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (hydro_con2prim,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    CCTK_REAL eneL = ene[index];
    CCTK_REAL epsL = eps[index];
    CCTK_REAL massL = mass[index];
    CCTK_REAL mom1L = mom1[index];
    CCTK_REAL mom2L = mom2[index];
    CCTK_REAL mom3L = mom3[index];
    CCTK_REAL rhoL = rho[index];
    CCTK_REAL vel1L = vel1[index];
    CCTK_REAL vel2L = vel2[index];
    CCTK_REAL vel3L = vel3[index];
    CCTK_REAL volL = vol[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    rhoL = massL*INV(volL);
    
    vel1L = mom1L*INV(massL);
    
    vel2L = mom2L*INV(massL);
    
    vel3L = mom3L*INV(massL);
    
    epsL = 0.5*INV(massL)*(2*eneL - massL*(SQR(vel1L) + SQR(vel2L) + 
      SQR(vel3L)));
    
    CCTK_REAL pressL = epsL*rhoL*ToReal(Gamma);
    
    
    /* Copy local copies back to grid functions */
    eps[index] = epsL;
    press[index] = pressL;
    rho[index] = rhoL;
    vel1[index] = vel1L;
    vel2[index] = vel2L;
    vel3[index] = vel3L;
  }
  LC_ENDLOOP3 (hydro_con2prim);
}

extern "C" void hydro_con2prim(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &hydro_con2prim_Body);
}
