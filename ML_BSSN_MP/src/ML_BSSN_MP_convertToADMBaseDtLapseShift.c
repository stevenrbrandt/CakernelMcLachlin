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

void ML_BSSN_MP_convertToADMBaseDtLapseShift_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  /* Declare predefined quantities */
  // CCTK_REAL p1o12dx = INITVALUE;
  // CCTK_REAL p1o12dy = INITVALUE;
  // CCTK_REAL p1o12dz = INITVALUE;
  // CCTK_REAL p1o144dxdy = INITVALUE;
  // CCTK_REAL p1o144dxdz = INITVALUE;
  // CCTK_REAL p1o144dydz = INITVALUE;
  // CCTK_REAL p1odx = INITVALUE;
  // CCTK_REAL p1ody = INITVALUE;
  // CCTK_REAL p1odz = INITVALUE;
  // CCTK_REAL pm1o12dx2 = INITVALUE;
  // CCTK_REAL pm1o12dy2 = INITVALUE;
  // CCTK_REAL pm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_convertToADMBaseDtLapseShift_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_convertToADMBaseDtLapseShift_calc_every != ML_BSSN_MP_convertToADMBaseDtLapseShift_calc_offset)
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
  CCTK_REAL const p1o12dx = INV(dx)/12.;
  CCTK_REAL const p1o12dy = INV(dy)/12.;
  CCTK_REAL const p1o12dz = INV(dz)/12.;
  CCTK_REAL const p1o144dxdy = (INV(dx)*INV(dy))/144.;
  CCTK_REAL const p1o144dxdz = (INV(dx)*INV(dz))/144.;
  CCTK_REAL const p1o144dydz = (INV(dy)*INV(dz))/144.;
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o12dx2 = -pow(dx,-2)/12.;
  CCTK_REAL const pm1o12dy2 = -pow(dy,-2)/12.;
  CCTK_REAL const pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_MP_convertToADMBaseDtLapseShift,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    // int subblock_index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int const subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    // CCTK_REAL dir1 = INITVALUE, dir2 = INITVALUE, dir3 = INITVALUE;
    
    /* Declare local copies of grid functions */
    // CCTK_REAL AL = INITVALUE;
    // CCTK_REAL alphaL = INITVALUE;
    // CCTK_REAL B1L = INITVALUE, B2L = INITVALUE, B3L = INITVALUE;
    // CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    // CCTK_REAL dtalpL = INITVALUE;
    // CCTK_REAL dtbetaxL = INITVALUE;
    // CCTK_REAL dtbetayL = INITVALUE;
    // CCTK_REAL dtbetazL = INITVALUE;
    // CCTK_REAL J11L = INITVALUE, J12L = INITVALUE, J13L = INITVALUE, J21L = INITVALUE, J22L = INITVALUE, J23L = INITVALUE;
    // CCTK_REAL J31L = INITVALUE, J32L = INITVALUE, J33L = INITVALUE;
    // CCTK_REAL trKL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    CCTK_REAL const AL = A[index];
    CCTK_REAL const alphaL = alpha[index];
    CCTK_REAL const B1L = B1[index];
    CCTK_REAL const B2L = B2[index];
    CCTK_REAL const B3L = B3[index];
    CCTK_REAL const beta1L = beta1[index];
    CCTK_REAL const beta2L = beta2[index];
    CCTK_REAL const beta3L = beta3[index];
    CCTK_REAL const J11L = J11[index];
    CCTK_REAL const J12L = J12[index];
    CCTK_REAL const J13L = J13[index];
    CCTK_REAL const J21L = J21[index];
    CCTK_REAL const J22L = J22[index];
    CCTK_REAL const J23L = J23[index];
    CCTK_REAL const J31L = J31[index];
    CCTK_REAL const J32L = J32[index];
    CCTK_REAL const J33L = J33[index];
    CCTK_REAL const trKL = trK[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    int const dir1  =  Sign(beta1L);
    
    int const dir2  =  Sign(beta2L);
    
    int const dir3  =  Sign(beta3L);
    
    CCTK_REAL const dtalpL  =  (PDupwindNth1(alpha, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
           PDupwindNth2(alpha, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
           PDupwindNth3(alpha, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L))*LapseAdvectionCoeff + 
        harmonicF*(AL*(-1 + LapseAdvectionCoeff) - LapseAdvectionCoeff*trKL)*pow(alphaL,harmonicN);
    
    CCTK_REAL const dtbetaxL  =  (PDupwindNth1(beta1, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
           PDupwindNth2(beta1, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
           PDupwindNth3(beta1, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L))*ShiftAdvectionCoeff + 
        B1L*ShiftGammaCoeff;
    
    CCTK_REAL const dtbetayL  =  (PDupwindNth1(beta2, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
           PDupwindNth2(beta2, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
           PDupwindNth3(beta2, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L))*ShiftAdvectionCoeff + 
        B2L*ShiftGammaCoeff;
    
    CCTK_REAL const dtbetazL  =  (PDupwindNth1(beta3, i, j, k)*(beta1L*J11L + beta2L*J12L + beta3L*J13L) + 
           PDupwindNth2(beta3, i, j, k)*(beta1L*J21L + beta2L*J22L + beta3L*J23L) + 
           PDupwindNth3(beta3, i, j, k)*(beta1L*J31L + beta2L*J32L + beta3L*J33L))*ShiftAdvectionCoeff + 
        B3L*ShiftGammaCoeff;
    
    
    /* Copy local copies back to grid functions */
    dtalp[index] = dtalpL;
    dtbetax[index] = dtbetaxL;
    dtbetay[index] = dtbetayL;
    dtbetaz[index] = dtbetazL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_BSSN_MP_convertToADMBaseDtLapseShift);
}

void ML_BSSN_MP_convertToADMBaseDtLapseShift(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_convertToADMBaseDtLapseShift_Body);
}
