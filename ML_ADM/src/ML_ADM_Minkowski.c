/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (32-bit) (April 20, 2007) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#define KRANC_C

#include <math.h>
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

void ML_ADM_Minkowski_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  
  /* Declare finite differencing variables */
  CCTK_REAL dx = INITVALUE, dy = INITVALUE, dz = INITVALUE;
  CCTK_REAL dxi = INITVALUE, dyi = INITVALUE, dzi = INITVALUE;
  CCTK_REAL khalf = INITVALUE, kthird = INITVALUE, ktwothird = INITVALUE, kfourthird = INITVALUE, keightthird = INITVALUE;
  CCTK_REAL hdxi = INITVALUE, hdyi = INITVALUE, hdzi = INITVALUE;
  
  
  /* Declare predefined quantities */
  CCTK_REAL qp1o12dx = INITVALUE;
  CCTK_REAL qp1o12dy = INITVALUE;
  CCTK_REAL qp1o12dz = INITVALUE;
  CCTK_REAL qp1o144dxdy = INITVALUE;
  CCTK_REAL qp1o144dxdz = INITVALUE;
  CCTK_REAL qp1o144dydz = INITVALUE;
  CCTK_REAL qpm1o12dx2 = INITVALUE;
  CCTK_REAL qpm1o12dy2 = INITVALUE;
  CCTK_REAL qpm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADM_Minkowski_Body");
  }
  
  if (cctk_iteration % ML_ADM_Minkowski_calc_every != ML_ADM_Minkowski_calc_offset)
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
  qp1o12dx = INV(dx)/12.;
  qp1o12dy = INV(dy)/12.;
  qp1o12dz = INV(dz)/12.;
  qp1o144dxdy = (INV(dx)*INV(dy))/144.;
  qp1o144dxdz = (INV(dx)*INV(dz))/144.;
  qp1o144dydz = (INV(dy)*INV(dz))/144.;
  qpm1o12dx2 = -pow(dx,-2)/12.;
  qpm1o12dy2 = -pow(dy,-2)/12.;
  qpm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  LC_LOOP3 (ML_ADM_Minkowski,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    
    /* Declare local copies of grid functions */
    CCTK_REAL alphaL = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    CCTK_REAL g11L = INITVALUE, g12L = INITVALUE, g13L = INITVALUE, g22L = INITVALUE, g23L = INITVALUE, g33L = INITVALUE;
    CCTK_REAL K11L = INITVALUE, K12L = INITVALUE, K13L = INITVALUE, K22L = INITVALUE, K23L = INITVALUE, K33L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    g11L  =  1;
    
    g12L  =  0;
    
    g13L  =  0;
    
    g22L  =  1;
    
    g23L  =  0;
    
    g33L  =  1;
    
    K11L  =  0;
    
    K12L  =  0;
    
    K13L  =  0;
    
    K22L  =  0;
    
    K23L  =  0;
    
    K33L  =  0;
    
    alphaL  =  1;
    
    beta1L  =  0;
    
    beta2L  =  0;
    
    beta3L  =  0;
    
    
    /* Copy local copies back to grid functions */
    alpha[index] = alphaL;
    beta1[index] = beta1L;
    beta2[index] = beta2L;
    beta3[index] = beta3L;
    g11[index] = g11L;
    g12[index] = g12L;
    g13[index] = g13L;
    g22[index] = g22L;
    g23[index] = g23L;
    g33[index] = g33L;
    K11[index] = K11L;
    K12[index] = K12L;
    K13[index] = K13L;
    K22[index] = K22L;
    K23[index] = K23L;
    K33[index] = K33L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_ADM_Minkowski);
}

void ML_ADM_Minkowski(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverEverything(cctkGH, &ML_ADM_Minkowski_Body);
}
