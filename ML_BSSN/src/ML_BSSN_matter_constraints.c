/*  File produced by user diener */
/*  Produced with Mathematica Version 6.0 for Linux x86 (32-bit) (April 20, 2007) */

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

void ML_BSSN_matter_constraints_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_matter_constraints_Body");
  }
  
  if (cctk_iteration % ML_BSSN_matter_constraints_calc_every != ML_BSSN_matter_constraints_calc_offset)
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
  LC_LOOP3 (ML_BSSN_matter_constraints,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL rho = INITVALUE;
    CCTK_REAL S1 = INITVALUE, S2 = INITVALUE, S3 = INITVALUE;
    CCTK_REAL T00 = INITVALUE, T01 = INITVALUE, T02 = INITVALUE, T03 = INITVALUE, T11 = INITVALUE, T12 = INITVALUE;
    CCTK_REAL T13 = INITVALUE, T22 = INITVALUE, T23 = INITVALUE, T33 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL alphaL = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    CCTK_REAL eTttL = INITVALUE;
    CCTK_REAL eTtxL = INITVALUE;
    CCTK_REAL eTtyL = INITVALUE;
    CCTK_REAL eTtzL = INITVALUE;
    CCTK_REAL eTxxL = INITVALUE;
    CCTK_REAL eTxyL = INITVALUE;
    CCTK_REAL eTxzL = INITVALUE;
    CCTK_REAL eTyyL = INITVALUE;
    CCTK_REAL eTyzL = INITVALUE;
    CCTK_REAL eTzzL = INITVALUE;
    CCTK_REAL HL = INITVALUE;
    CCTK_REAL M1L = INITVALUE, M2L = INITVALUE, M3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    alphaL = alpha[index];
    beta1L = beta1[index];
    beta2L = beta2[index];
    beta3L = beta3[index];
    eTttL = eTtt[index];
    eTtxL = eTtx[index];
    eTtyL = eTty[index];
    eTtzL = eTtz[index];
    eTxxL = eTxx[index];
    eTxyL = eTxy[index];
    eTxzL = eTxz[index];
    eTyyL = eTyy[index];
    eTyzL = eTyz[index];
    eTzzL = eTzz[index];
    HL = H[index];
    M1L = M1[index];
    M2L = M2[index];
    M3L = M3[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    T00  =  eTttL;
    
    T01  =  eTtxL;
    
    T02  =  eTtyL;
    
    T03  =  eTtzL;
    
    T11  =  eTxxL;
    
    T12  =  eTxyL;
    
    T13  =  eTxzL;
    
    T22  =  eTyyL;
    
    T23  =  eTyzL;
    
    T33  =  eTzzL;
    
    rho  =  pow(alphaL,-2)*(T00 - 2*(beta2L*T02 + beta3L*T03) + 
          2*(beta1L*(-T01 + beta2L*T12 + beta3L*T13) + beta2L*beta3L*T23) + T11*SQR(beta1L) + T22*SQR(beta2L) + 
          T33*SQR(beta3L));
    
    S1  =  (-T01 + beta1L*T11 + beta2L*T12 + beta3L*T13)*INV(alphaL);
    
    S2  =  (-T02 + beta1L*T12 + beta2L*T22 + beta3L*T23)*INV(alphaL);
    
    S3  =  (-T03 + beta1L*T13 + beta2L*T23 + beta3L*T33)*INV(alphaL);
    
    HL  =  HL - 50.26548245743669181540229413247204614715*rho;
    
    M1L  =  M1L - 25.13274122871834590770114706623602307358*S1;
    
    M2L  =  M2L - 25.13274122871834590770114706623602307358*S2;
    
    M3L  =  M3L - 25.13274122871834590770114706623602307358*S3;
    
    
    /* Copy local copies back to grid functions */
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_BSSN_matter_constraints);
}

void ML_BSSN_matter_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_matter_constraints_Body);
}
