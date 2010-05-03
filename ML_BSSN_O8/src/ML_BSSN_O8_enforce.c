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

void ML_BSSN_O8_enforce_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O8_enforce_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O8_enforce_calc_every != ML_BSSN_O8_enforce_calc_offset)
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
  CCTK_REAL const p1o1024dx = INV(dx)/1024.;
  CCTK_REAL const p1o1024dy = INV(dy)/1024.;
  CCTK_REAL const p1o1024dz = INV(dz)/1024.;
  CCTK_REAL const p1o1680dx = INV(dx)/1680.;
  CCTK_REAL const p1o1680dy = INV(dy)/1680.;
  CCTK_REAL const p1o1680dz = INV(dz)/1680.;
  CCTK_REAL const p1o5040dx2 = pow(dx,-2)/5040.;
  CCTK_REAL const p1o5040dy2 = pow(dy,-2)/5040.;
  CCTK_REAL const p1o5040dz2 = pow(dz,-2)/5040.;
  CCTK_REAL const p1o560dx = INV(dx)/560.;
  CCTK_REAL const p1o560dy = INV(dy)/560.;
  CCTK_REAL const p1o560dz = INV(dz)/560.;
  CCTK_REAL const p1o705600dxdy = (INV(dx)*INV(dy))/705600.;
  CCTK_REAL const p1o705600dxdz = (INV(dx)*INV(dz))/705600.;
  CCTK_REAL const p1o705600dydz = (INV(dy)*INV(dz))/705600.;
  CCTK_REAL const p1o840dx = INV(dx)/840.;
  CCTK_REAL const p1o840dy = INV(dy)/840.;
  CCTK_REAL const p1o840dz = INV(dz)/840.;
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o840dx = -INV(dx)/840.;
  CCTK_REAL const pm1o840dy = -INV(dy)/840.;
  CCTK_REAL const pm1o840dz = -INV(dz)/840.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_O8_enforce,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    CCTK_REAL  alphaL = alpha[index];
    CCTK_REAL  At11L = At11[index];
    CCTK_REAL  At12L = At12[index];
    CCTK_REAL  At13L = At13[index];
    CCTK_REAL  At22L = At22[index];
    CCTK_REAL  At23L = At23[index];
    CCTK_REAL  At33L = At33[index];
    CCTK_REAL  gt11L = gt11[index];
    CCTK_REAL  gt12L = gt12[index];
    CCTK_REAL  gt13L = gt13[index];
    CCTK_REAL  gt22L = gt22[index];
    CCTK_REAL  gt23L = gt23[index];
    CCTK_REAL  gt33L = gt33[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu21 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu31 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu32 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL trAt = At11L*gtu11 + At22L*gtu22 + 2*(At12L*gtu21 + 
      At13L*gtu31 + At23L*gtu32) + At33L*gtu33;
    
    At11L = At11L - gt11L*kthird*trAt;
    
    At12L = At12L - gt12L*kthird*trAt;
    
    At13L = At13L - gt13L*kthird*trAt;
    
    At22L = At22L - gt22L*kthird*trAt;
    
    At23L = At23L - gt23L*kthird*trAt;
    
    At33L = At33L - gt33L*kthird*trAt;
    
    alphaL = fmax(alphaL,MinimumLapse);
    
    
    /* Copy local copies back to grid functions */
    alpha[index] = alphaL;
    At11[index] = At11L;
    At12[index] = At12L;
    At13[index] = At13L;
    At22[index] = At22L;
    At23[index] = At23L;
    At33[index] = At33L;
  }
  LC_ENDLOOP3 (ML_BSSN_O8_enforce);
}

void ML_BSSN_O8_enforce(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &ML_BSSN_O8_enforce_Body);
}
