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

void ML_ADM_RHS_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  
  /* Declare finite differencing variables */
  CCTK_REAL dx = INITVALUE, dy = INITVALUE, dz = INITVALUE;
  CCTK_REAL dxi = INITVALUE, dyi = INITVALUE, dzi = INITVALUE;
  CCTK_REAL khalf = INITVALUE, kthird = INITVALUE, ktwothird = INITVALUE, kfourthird = INITVALUE, keightthird = INITVALUE;
  CCTK_REAL hdxi = INITVALUE, hdyi = INITVALUE, hdzi = INITVALUE;
  
  
  /* Declare predefined quantities */
  CCTK_REAL p1o12dx = INITVALUE;
  CCTK_REAL p1o12dy = INITVALUE;
  CCTK_REAL p1o12dz = INITVALUE;
  CCTK_REAL p1o144dxdy = INITVALUE;
  CCTK_REAL p1o144dxdz = INITVALUE;
  CCTK_REAL p1o144dydz = INITVALUE;
  CCTK_REAL pm1o12dx2 = INITVALUE;
  CCTK_REAL pm1o12dy2 = INITVALUE;
  CCTK_REAL pm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADM_RHS_Body");
  }
  
  if (cctk_iteration % ML_ADM_RHS_calc_every != ML_ADM_RHS_calc_offset)
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
  p1o144dxdy = (INV(dx)*INV(dy))/144.;
  p1o144dxdz = (INV(dx)*INV(dz))/144.;
  p1o144dydz = (INV(dy)*INV(dz))/144.;
  pm1o12dx2 = -pow(dx,-2)/12.;
  pm1o12dy2 = -pow(dy,-2)/12.;
  pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  LC_LOOP3 (somename,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL detg = INITVALUE;
    CCTK_REAL G111 = INITVALUE, G112 = INITVALUE, G113 = INITVALUE, G122 = INITVALUE, G123 = INITVALUE, G133 = INITVALUE;
    CCTK_REAL G211 = INITVALUE, G212 = INITVALUE, G213 = INITVALUE, G222 = INITVALUE, G223 = INITVALUE, G233 = INITVALUE;
    CCTK_REAL G311 = INITVALUE, G312 = INITVALUE, G313 = INITVALUE, G322 = INITVALUE, G323 = INITVALUE, G333 = INITVALUE;
    CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    CCTK_REAL Km11 = INITVALUE, Km12 = INITVALUE, Km13 = INITVALUE, Km21 = INITVALUE, Km22 = INITVALUE, Km23 = INITVALUE;
    CCTK_REAL Km31 = INITVALUE, Km32 = INITVALUE, Km33 = INITVALUE;
    CCTK_REAL R11 = INITVALUE, R12 = INITVALUE, R13 = INITVALUE, R22 = INITVALUE, R23 = INITVALUE, R33 = INITVALUE;
    CCTK_REAL trK = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL alphaL = INITVALUE, alpharhsL = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta1rhsL = INITVALUE, beta2L = INITVALUE, beta2rhsL = INITVALUE, beta3L = INITVALUE, beta3rhsL = INITVALUE;
    CCTK_REAL g11L = INITVALUE, g11rhsL = INITVALUE, g12L = INITVALUE, g12rhsL = INITVALUE, g13L = INITVALUE, g13rhsL = INITVALUE;
    CCTK_REAL g22L = INITVALUE, g22rhsL = INITVALUE, g23L = INITVALUE, g23rhsL = INITVALUE, g33L = INITVALUE, g33rhsL = INITVALUE;
    CCTK_REAL K11L = INITVALUE, K11rhsL = INITVALUE, K12L = INITVALUE, K12rhsL = INITVALUE, K13L = INITVALUE, K13rhsL = INITVALUE;
    CCTK_REAL K22L = INITVALUE, K22rhsL = INITVALUE, K23L = INITVALUE, K23rhsL = INITVALUE, K33L = INITVALUE, K33rhsL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandard4th1alpha = INITVALUE;
    CCTK_REAL PDstandard4th2alpha = INITVALUE;
    CCTK_REAL PDstandard4th3alpha = INITVALUE;
    CCTK_REAL PDstandard4th11alpha = INITVALUE;
    CCTK_REAL PDstandard4th22alpha = INITVALUE;
    CCTK_REAL PDstandard4th33alpha = INITVALUE;
    CCTK_REAL PDstandard4th12alpha = INITVALUE;
    CCTK_REAL PDstandard4th13alpha = INITVALUE;
    CCTK_REAL PDstandard4th23alpha = INITVALUE;
    CCTK_REAL PDstandard4th1beta1 = INITVALUE;
    CCTK_REAL PDstandard4th2beta1 = INITVALUE;
    CCTK_REAL PDstandard4th3beta1 = INITVALUE;
    CCTK_REAL PDstandard4th1beta2 = INITVALUE;
    CCTK_REAL PDstandard4th2beta2 = INITVALUE;
    CCTK_REAL PDstandard4th3beta2 = INITVALUE;
    CCTK_REAL PDstandard4th1beta3 = INITVALUE;
    CCTK_REAL PDstandard4th2beta3 = INITVALUE;
    CCTK_REAL PDstandard4th3beta3 = INITVALUE;
    CCTK_REAL PDstandard4th1g11 = INITVALUE;
    CCTK_REAL PDstandard4th2g11 = INITVALUE;
    CCTK_REAL PDstandard4th3g11 = INITVALUE;
    CCTK_REAL PDstandard4th22g11 = INITVALUE;
    CCTK_REAL PDstandard4th33g11 = INITVALUE;
    CCTK_REAL PDstandard4th12g11 = INITVALUE;
    CCTK_REAL PDstandard4th13g11 = INITVALUE;
    CCTK_REAL PDstandard4th21g11 = INITVALUE;
    CCTK_REAL PDstandard4th23g11 = INITVALUE;
    CCTK_REAL PDstandard4th31g11 = INITVALUE;
    CCTK_REAL PDstandard4th32g11 = INITVALUE;
    CCTK_REAL PDstandard4th1g12 = INITVALUE;
    CCTK_REAL PDstandard4th2g12 = INITVALUE;
    CCTK_REAL PDstandard4th3g12 = INITVALUE;
    CCTK_REAL PDstandard4th33g12 = INITVALUE;
    CCTK_REAL PDstandard4th12g12 = INITVALUE;
    CCTK_REAL PDstandard4th13g12 = INITVALUE;
    CCTK_REAL PDstandard4th21g12 = INITVALUE;
    CCTK_REAL PDstandard4th23g12 = INITVALUE;
    CCTK_REAL PDstandard4th31g12 = INITVALUE;
    CCTK_REAL PDstandard4th32g12 = INITVALUE;
    CCTK_REAL PDstandard4th1g13 = INITVALUE;
    CCTK_REAL PDstandard4th2g13 = INITVALUE;
    CCTK_REAL PDstandard4th3g13 = INITVALUE;
    CCTK_REAL PDstandard4th22g13 = INITVALUE;
    CCTK_REAL PDstandard4th12g13 = INITVALUE;
    CCTK_REAL PDstandard4th13g13 = INITVALUE;
    CCTK_REAL PDstandard4th21g13 = INITVALUE;
    CCTK_REAL PDstandard4th23g13 = INITVALUE;
    CCTK_REAL PDstandard4th31g13 = INITVALUE;
    CCTK_REAL PDstandard4th32g13 = INITVALUE;
    CCTK_REAL PDstandard4th1g22 = INITVALUE;
    CCTK_REAL PDstandard4th2g22 = INITVALUE;
    CCTK_REAL PDstandard4th3g22 = INITVALUE;
    CCTK_REAL PDstandard4th11g22 = INITVALUE;
    CCTK_REAL PDstandard4th33g22 = INITVALUE;
    CCTK_REAL PDstandard4th12g22 = INITVALUE;
    CCTK_REAL PDstandard4th13g22 = INITVALUE;
    CCTK_REAL PDstandard4th21g22 = INITVALUE;
    CCTK_REAL PDstandard4th23g22 = INITVALUE;
    CCTK_REAL PDstandard4th31g22 = INITVALUE;
    CCTK_REAL PDstandard4th32g22 = INITVALUE;
    CCTK_REAL PDstandard4th1g23 = INITVALUE;
    CCTK_REAL PDstandard4th2g23 = INITVALUE;
    CCTK_REAL PDstandard4th3g23 = INITVALUE;
    CCTK_REAL PDstandard4th11g23 = INITVALUE;
    CCTK_REAL PDstandard4th12g23 = INITVALUE;
    CCTK_REAL PDstandard4th13g23 = INITVALUE;
    CCTK_REAL PDstandard4th21g23 = INITVALUE;
    CCTK_REAL PDstandard4th23g23 = INITVALUE;
    CCTK_REAL PDstandard4th31g23 = INITVALUE;
    CCTK_REAL PDstandard4th32g23 = INITVALUE;
    CCTK_REAL PDstandard4th1g33 = INITVALUE;
    CCTK_REAL PDstandard4th2g33 = INITVALUE;
    CCTK_REAL PDstandard4th3g33 = INITVALUE;
    CCTK_REAL PDstandard4th11g33 = INITVALUE;
    CCTK_REAL PDstandard4th22g33 = INITVALUE;
    CCTK_REAL PDstandard4th12g33 = INITVALUE;
    CCTK_REAL PDstandard4th13g33 = INITVALUE;
    CCTK_REAL PDstandard4th21g33 = INITVALUE;
    CCTK_REAL PDstandard4th23g33 = INITVALUE;
    CCTK_REAL PDstandard4th31g33 = INITVALUE;
    CCTK_REAL PDstandard4th32g33 = INITVALUE;
    CCTK_REAL PDstandard4th1K11 = INITVALUE;
    CCTK_REAL PDstandard4th2K11 = INITVALUE;
    CCTK_REAL PDstandard4th3K11 = INITVALUE;
    CCTK_REAL PDstandard4th1K12 = INITVALUE;
    CCTK_REAL PDstandard4th2K12 = INITVALUE;
    CCTK_REAL PDstandard4th3K12 = INITVALUE;
    CCTK_REAL PDstandard4th1K13 = INITVALUE;
    CCTK_REAL PDstandard4th2K13 = INITVALUE;
    CCTK_REAL PDstandard4th3K13 = INITVALUE;
    CCTK_REAL PDstandard4th1K22 = INITVALUE;
    CCTK_REAL PDstandard4th2K22 = INITVALUE;
    CCTK_REAL PDstandard4th3K22 = INITVALUE;
    CCTK_REAL PDstandard4th1K23 = INITVALUE;
    CCTK_REAL PDstandard4th2K23 = INITVALUE;
    CCTK_REAL PDstandard4th3K23 = INITVALUE;
    CCTK_REAL PDstandard4th1K33 = INITVALUE;
    CCTK_REAL PDstandard4th2K33 = INITVALUE;
    CCTK_REAL PDstandard4th3K33 = INITVALUE;
    
    /* Assign local copies of grid functions */
    alphaL = alpha[index];
    beta1L = beta1[index];
    beta2L = beta2[index];
    beta3L = beta3[index];
    g11L = g11[index];
    g12L = g12[index];
    g13L = g13[index];
    g22L = g22[index];
    g23L = g23[index];
    g33L = g33[index];
    K11L = K11[index];
    K12L = K12[index];
    K13L = K13[index];
    K22L = K22[index];
    K23L = K23[index];
    K33L = K33[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandard4th1alpha = PDstandard4th1(alpha, i, j, k);
    PDstandard4th2alpha = PDstandard4th2(alpha, i, j, k);
    PDstandard4th3alpha = PDstandard4th3(alpha, i, j, k);
    PDstandard4th11alpha = PDstandard4th11(alpha, i, j, k);
    PDstandard4th22alpha = PDstandard4th22(alpha, i, j, k);
    PDstandard4th33alpha = PDstandard4th33(alpha, i, j, k);
    PDstandard4th12alpha = PDstandard4th12(alpha, i, j, k);
    PDstandard4th13alpha = PDstandard4th13(alpha, i, j, k);
    PDstandard4th23alpha = PDstandard4th23(alpha, i, j, k);
    PDstandard4th1beta1 = PDstandard4th1(beta1, i, j, k);
    PDstandard4th2beta1 = PDstandard4th2(beta1, i, j, k);
    PDstandard4th3beta1 = PDstandard4th3(beta1, i, j, k);
    PDstandard4th1beta2 = PDstandard4th1(beta2, i, j, k);
    PDstandard4th2beta2 = PDstandard4th2(beta2, i, j, k);
    PDstandard4th3beta2 = PDstandard4th3(beta2, i, j, k);
    PDstandard4th1beta3 = PDstandard4th1(beta3, i, j, k);
    PDstandard4th2beta3 = PDstandard4th2(beta3, i, j, k);
    PDstandard4th3beta3 = PDstandard4th3(beta3, i, j, k);
    PDstandard4th1g11 = PDstandard4th1(g11, i, j, k);
    PDstandard4th2g11 = PDstandard4th2(g11, i, j, k);
    PDstandard4th3g11 = PDstandard4th3(g11, i, j, k);
    PDstandard4th22g11 = PDstandard4th22(g11, i, j, k);
    PDstandard4th33g11 = PDstandard4th33(g11, i, j, k);
    PDstandard4th23g11 = PDstandard4th23(g11, i, j, k);
    PDstandard4th1g12 = PDstandard4th1(g12, i, j, k);
    PDstandard4th2g12 = PDstandard4th2(g12, i, j, k);
    PDstandard4th3g12 = PDstandard4th3(g12, i, j, k);
    PDstandard4th33g12 = PDstandard4th33(g12, i, j, k);
    PDstandard4th12g12 = PDstandard4th12(g12, i, j, k);
    PDstandard4th13g12 = PDstandard4th13(g12, i, j, k);
    PDstandard4th23g12 = PDstandard4th23(g12, i, j, k);
    PDstandard4th1g13 = PDstandard4th1(g13, i, j, k);
    PDstandard4th2g13 = PDstandard4th2(g13, i, j, k);
    PDstandard4th3g13 = PDstandard4th3(g13, i, j, k);
    PDstandard4th22g13 = PDstandard4th22(g13, i, j, k);
    PDstandard4th12g13 = PDstandard4th12(g13, i, j, k);
    PDstandard4th13g13 = PDstandard4th13(g13, i, j, k);
    PDstandard4th23g13 = PDstandard4th23(g13, i, j, k);
    PDstandard4th1g22 = PDstandard4th1(g22, i, j, k);
    PDstandard4th2g22 = PDstandard4th2(g22, i, j, k);
    PDstandard4th3g22 = PDstandard4th3(g22, i, j, k);
    PDstandard4th11g22 = PDstandard4th11(g22, i, j, k);
    PDstandard4th33g22 = PDstandard4th33(g22, i, j, k);
    PDstandard4th13g22 = PDstandard4th13(g22, i, j, k);
    PDstandard4th1g23 = PDstandard4th1(g23, i, j, k);
    PDstandard4th2g23 = PDstandard4th2(g23, i, j, k);
    PDstandard4th3g23 = PDstandard4th3(g23, i, j, k);
    PDstandard4th11g23 = PDstandard4th11(g23, i, j, k);
    PDstandard4th12g23 = PDstandard4th12(g23, i, j, k);
    PDstandard4th13g23 = PDstandard4th13(g23, i, j, k);
    PDstandard4th23g23 = PDstandard4th23(g23, i, j, k);
    PDstandard4th1g33 = PDstandard4th1(g33, i, j, k);
    PDstandard4th2g33 = PDstandard4th2(g33, i, j, k);
    PDstandard4th3g33 = PDstandard4th3(g33, i, j, k);
    PDstandard4th11g33 = PDstandard4th11(g33, i, j, k);
    PDstandard4th22g33 = PDstandard4th22(g33, i, j, k);
    PDstandard4th12g33 = PDstandard4th12(g33, i, j, k);
    PDstandard4th1K11 = PDstandard4th1(K11, i, j, k);
    PDstandard4th2K11 = PDstandard4th2(K11, i, j, k);
    PDstandard4th3K11 = PDstandard4th3(K11, i, j, k);
    PDstandard4th1K12 = PDstandard4th1(K12, i, j, k);
    PDstandard4th2K12 = PDstandard4th2(K12, i, j, k);
    PDstandard4th3K12 = PDstandard4th3(K12, i, j, k);
    PDstandard4th1K13 = PDstandard4th1(K13, i, j, k);
    PDstandard4th2K13 = PDstandard4th2(K13, i, j, k);
    PDstandard4th3K13 = PDstandard4th3(K13, i, j, k);
    PDstandard4th1K22 = PDstandard4th1(K22, i, j, k);
    PDstandard4th2K22 = PDstandard4th2(K22, i, j, k);
    PDstandard4th3K22 = PDstandard4th3(K22, i, j, k);
    PDstandard4th1K23 = PDstandard4th1(K23, i, j, k);
    PDstandard4th2K23 = PDstandard4th2(K23, i, j, k);
    PDstandard4th3K23 = PDstandard4th3(K23, i, j, k);
    PDstandard4th1K33 = PDstandard4th1(K33, i, j, k);
    PDstandard4th2K33 = PDstandard4th2(K33, i, j, k);
    PDstandard4th3K33 = PDstandard4th3(K33, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    detg  =  2*g12L*g13L*g23L + g33L*(g11L*g22L - SQR(g12L)) - g22L*SQR(g13L) - g11L*SQR(g23L);
    
    gu11  =  INV(detg)*(g22L*g33L - SQR(g23L));
    
    gu21  =  (g13L*g23L - g12L*g33L)*INV(detg);
    
    gu31  =  (-(g13L*g22L) + g12L*g23L)*INV(detg);
    
    gu22  =  INV(detg)*(g11L*g33L - SQR(g13L));
    
    gu32  =  (g12L*g13L - g11L*g23L)*INV(detg);
    
    gu33  =  INV(detg)*(g11L*g22L - SQR(g12L));
    
    G111  =  khalf*(gu11*PDstandard4th1g11 + 2*(gu21*PDstandard4th1g12 + gu31*PDstandard4th1g13) - gu21*PDstandard4th2g11 - 
          gu31*PDstandard4th3g11);
    
    G211  =  khalf*(gu21*PDstandard4th1g11 + 2*(gu22*PDstandard4th1g12 + gu32*PDstandard4th1g13) - gu22*PDstandard4th2g11 - 
          gu32*PDstandard4th3g11);
    
    G311  =  khalf*(gu31*PDstandard4th1g11 + 2*(gu32*PDstandard4th1g12 + gu33*PDstandard4th1g13) - gu32*PDstandard4th2g11 - 
          gu33*PDstandard4th3g11);
    
    G112  =  khalf*(gu21*PDstandard4th1g22 + gu11*PDstandard4th2g11 + 
          gu31*(PDstandard4th1g23 + PDstandard4th2g13 - PDstandard4th3g12));
    
    G212  =  khalf*(gu22*PDstandard4th1g22 + gu21*PDstandard4th2g11 + 
          gu32*(PDstandard4th1g23 + PDstandard4th2g13 - PDstandard4th3g12));
    
    G312  =  khalf*(gu32*PDstandard4th1g22 + gu31*PDstandard4th2g11 + 
          gu33*(PDstandard4th1g23 + PDstandard4th2g13 - PDstandard4th3g12));
    
    G113  =  khalf*(gu31*PDstandard4th1g33 + gu11*PDstandard4th3g11 + 
          gu21*(PDstandard4th1g23 - PDstandard4th2g13 + PDstandard4th3g12));
    
    G213  =  khalf*(gu32*PDstandard4th1g33 + gu21*PDstandard4th3g11 + 
          gu22*(PDstandard4th1g23 - PDstandard4th2g13 + PDstandard4th3g12));
    
    G313  =  khalf*(gu33*PDstandard4th1g33 + gu31*PDstandard4th3g11 + 
          gu32*(PDstandard4th1g23 - PDstandard4th2g13 + PDstandard4th3g12));
    
    G122  =  khalf*(gu11*(-PDstandard4th1g22 + 2*PDstandard4th2g12) + gu21*PDstandard4th2g22 + 
          gu31*(2*PDstandard4th2g23 - PDstandard4th3g22));
    
    G222  =  khalf*(gu21*(-PDstandard4th1g22 + 2*PDstandard4th2g12) + gu22*PDstandard4th2g22 + 
          gu32*(2*PDstandard4th2g23 - PDstandard4th3g22));
    
    G322  =  khalf*(gu31*(-PDstandard4th1g22 + 2*PDstandard4th2g12) + gu32*PDstandard4th2g22 + 
          gu33*(2*PDstandard4th2g23 - PDstandard4th3g22));
    
    G123  =  khalf*(gu31*PDstandard4th2g33 + gu11*(-PDstandard4th1g23 + PDstandard4th2g13 + PDstandard4th3g12) + 
          gu21*PDstandard4th3g22);
    
    G223  =  khalf*(gu32*PDstandard4th2g33 + gu21*(-PDstandard4th1g23 + PDstandard4th2g13 + PDstandard4th3g12) + 
          gu22*PDstandard4th3g22);
    
    G323  =  khalf*(gu33*PDstandard4th2g33 + gu31*(-PDstandard4th1g23 + PDstandard4th2g13 + PDstandard4th3g12) + 
          gu32*PDstandard4th3g22);
    
    G133  =  khalf*(-(gu11*PDstandard4th1g33) - gu21*PDstandard4th2g33 + 2*gu11*PDstandard4th3g13 + 
          2*gu21*PDstandard4th3g23 + gu31*PDstandard4th3g33);
    
    G233  =  khalf*(-(gu21*PDstandard4th1g33) - gu22*PDstandard4th2g33 + 2*gu21*PDstandard4th3g13 + 
          2*gu22*PDstandard4th3g23 + gu32*PDstandard4th3g33);
    
    G333  =  khalf*(-(gu31*PDstandard4th1g33) - gu32*PDstandard4th2g33 + 2*gu31*PDstandard4th3g13 + 
          2*gu32*PDstandard4th3g23 + gu33*PDstandard4th3g33);
    
    R11  =  khalf*(4*G213*G312 + G211*(2*G112 - 2*G222 - 2*G323) + G311*(2*G113 - 2*G333) - gu22*PDstandard4th11g22 - 
          2*(G111*G212 + G223*G311 + G111*G313 + gu32*PDstandard4th11g23) - gu33*PDstandard4th11g33 + 
          2*gu22*PDstandard4th12g12 + 2*gu32*PDstandard4th12g13 + 2*gu32*PDstandard4th13g12 + 2*gu33*PDstandard4th13g13 - 
          gu22*PDstandard4th22g11 - 2*gu32*PDstandard4th23g11 - gu33*PDstandard4th33g11 + 2*SQR(G212) + 2*SQR(G313));
    
    R12  =  khalf*(2*(G122*G211 + G123*G311 + G213*G322 + G313*G323) - 
          2*(G112*G212 + G112*G313 + G212*G323 + G312*G333 + gu21*PDstandard4th12g12) - gu32*PDstandard4th12g23 - 
          gu33*PDstandard4th12g33 + gu31*(PDstandard4th11g23 - PDstandard4th12g13 - PDstandard4th13g12) + 
          gu32*PDstandard4th13g22 + gu33*PDstandard4th13g23 + gu21*(PDstandard4th11g22 + PDstandard4th22g11) + 
          gu32*PDstandard4th22g13 + gu31*PDstandard4th23g11 - gu32*PDstandard4th23g12 + gu33*PDstandard4th23g13 - 
          gu33*PDstandard4th33g12);
    
    R13  =  khalf*(2*(G123*G211 + G212*G223 + G133*G311 + G233*G312) - 
          2*(G213*G222 + G223*G313 + G113*(G212 + G313) + gu31*PDstandard4th13g13) + 
          gu21*(PDstandard4th11g23 - PDstandard4th12g13 - PDstandard4th13g12 + PDstandard4th23g11) + 
          gu22*(PDstandard4th12g23 - PDstandard4th13g22 - PDstandard4th22g13 + PDstandard4th23g12) + 
          gu31*(PDstandard4th11g33 + PDstandard4th33g11) + 
          gu32*(PDstandard4th12g33 - PDstandard4th13g23 - PDstandard4th23g13 + PDstandard4th33g12));
    
    R22  =  khalf*(4*G123*G312 + G122*(2*G212 - 2*(G111 + G313)) + G322*(2*G223 - 2*G333) - 
          2*(G113*G322 + G222*(G112 + G323) + gu31*PDstandard4th13g22) + 
          gu11*(-PDstandard4th11g22 + 2*PDstandard4th12g12 - PDstandard4th22g11) + 
          gu31*(-2*PDstandard4th22g13 + 2*(PDstandard4th12g23 + PDstandard4th23g12)) + 
          gu33*(-PDstandard4th22g33 + 2*PDstandard4th23g23 - PDstandard4th33g22) + 2*(SQR(G112) + SQR(G323)));
    
    R23  =  khalf*(2*(G112*G113 + G122*G213 + G133*G312 + G233*G322) + 
          gu11*(-PDstandard4th11g23 + PDstandard4th12g13 + PDstandard4th13g12 - PDstandard4th23g11) + 
          gu21*(-PDstandard4th12g23 + PDstandard4th13g22 + PDstandard4th22g13 - PDstandard4th23g12) - 
          2*(G111*G123 + G113*G323 + G223*(G112 + G323) + gu32*PDstandard4th23g23) + 
          gu31*(PDstandard4th12g33 - PDstandard4th13g23 - PDstandard4th23g13 + PDstandard4th33g12) + 
          gu32*(PDstandard4th22g33 + PDstandard4th33g22));
    
    R33  =  khalf*(4*G123*G213 - gu11*PDstandard4th11g33 - 
          2*(G111*G133 + G133*G212 + G112*G233 + G222*G233 + G113*G333 + G223*G333 + gu21*PDstandard4th12g33) + 
          2*(G133*G313 + G233*G323 + gu11*PDstandard4th13g13) + 2*gu21*PDstandard4th13g23 - gu22*PDstandard4th22g33 + 
          2*gu21*PDstandard4th23g13 + 2*gu22*PDstandard4th23g23 - gu11*PDstandard4th33g11 - 2*gu21*PDstandard4th33g12 - 
          gu22*PDstandard4th33g22 + 2*SQR(G113) + 2*SQR(G223));
    
    Km11  =  gu11*K11L + gu21*K12L + gu31*K13L;
    
    Km21  =  gu21*K11L + gu22*K12L + gu32*K13L;
    
    Km31  =  gu31*K11L + gu32*K12L + gu33*K13L;
    
    Km12  =  gu11*K12L + gu21*K22L + gu31*K23L;
    
    Km22  =  gu21*K12L + gu22*K22L + gu32*K23L;
    
    Km32  =  gu31*K12L + gu32*K22L + gu33*K23L;
    
    Km13  =  gu11*K13L + gu21*K23L + gu31*K33L;
    
    Km23  =  gu21*K13L + gu22*K23L + gu32*K33L;
    
    Km33  =  gu31*K13L + gu32*K23L + gu33*K33L;
    
    trK  =  Km11 + Km22 + Km33;
    
    g11rhsL  =  -2*alphaL*K11L + 2*(g11L*PDstandard4th1beta1 + g12L*PDstandard4th1beta2 + g13L*PDstandard4th1beta3) + 
        beta1L*PDstandard4th1g11 + beta2L*PDstandard4th2g11 + beta3L*PDstandard4th3g11;
    
    g12rhsL  =  -2*alphaL*K12L + g22L*PDstandard4th1beta2 + g23L*PDstandard4th1beta3 + beta1L*PDstandard4th1g12 + 
        g11L*PDstandard4th2beta1 + g12L*(PDstandard4th1beta1 + PDstandard4th2beta2) + g13L*PDstandard4th2beta3 + 
        beta2L*PDstandard4th2g12 + beta3L*PDstandard4th3g12;
    
    g13rhsL  =  -2*alphaL*K13L + g23L*PDstandard4th1beta2 + g33L*PDstandard4th1beta3 + beta1L*PDstandard4th1g13 + 
        beta2L*PDstandard4th2g13 + g11L*PDstandard4th3beta1 + g12L*PDstandard4th3beta2 + 
        g13L*(PDstandard4th1beta1 + PDstandard4th3beta3) + beta3L*PDstandard4th3g13;
    
    g22rhsL  =  -2*alphaL*K22L + beta1L*PDstandard4th1g22 + 
        2*(g12L*PDstandard4th2beta1 + g22L*PDstandard4th2beta2 + g23L*PDstandard4th2beta3) + beta2L*PDstandard4th2g22 + 
        beta3L*PDstandard4th3g22;
    
    g23rhsL  =  -2*alphaL*K23L + beta1L*PDstandard4th1g23 + g13L*PDstandard4th2beta1 + g33L*PDstandard4th2beta3 + 
        beta2L*PDstandard4th2g23 + g12L*PDstandard4th3beta1 + g22L*PDstandard4th3beta2 + 
        g23L*(PDstandard4th2beta2 + PDstandard4th3beta3) + beta3L*PDstandard4th3g23;
    
    g33rhsL  =  -2*alphaL*K33L + beta1L*PDstandard4th1g33 + beta2L*PDstandard4th2g33 + 
        2*(g13L*PDstandard4th3beta1 + g23L*PDstandard4th3beta2 + g33L*PDstandard4th3beta3) + beta3L*PDstandard4th3g33;
    
    K11rhsL  =  -PDstandard4th11alpha + G111*PDstandard4th1alpha + 
        2*(K11L*PDstandard4th1beta1 + K12L*PDstandard4th1beta2 + K13L*PDstandard4th1beta3) + beta1L*PDstandard4th1K11 + 
        G211*PDstandard4th2alpha + beta2L*PDstandard4th2K11 + G311*PDstandard4th3alpha + beta3L*PDstandard4th3K11 + 
        alphaL*(-2*(K11L*Km11 + K12L*Km21 + K13L*Km31) + R11 + K11L*trK);
    
    K12rhsL  =  -PDstandard4th12alpha + G112*PDstandard4th1alpha + K22L*PDstandard4th1beta2 + K23L*PDstandard4th1beta3 + 
        beta1L*PDstandard4th1K12 + G212*PDstandard4th2alpha + K11L*PDstandard4th2beta1 + 
        K12L*(PDstandard4th1beta1 + PDstandard4th2beta2) + K13L*PDstandard4th2beta3 + beta2L*PDstandard4th2K12 + 
        G312*PDstandard4th3alpha + beta3L*PDstandard4th3K12 + 
        alphaL*(-2*(K11L*Km12 + K12L*Km22 + K13L*Km32) + R12 + K12L*trK);
    
    K13rhsL  =  -PDstandard4th13alpha + G113*PDstandard4th1alpha + K23L*PDstandard4th1beta2 + K33L*PDstandard4th1beta3 + 
        beta1L*PDstandard4th1K13 + G213*PDstandard4th2alpha + beta2L*PDstandard4th2K13 + G313*PDstandard4th3alpha + 
        K11L*PDstandard4th3beta1 + K12L*PDstandard4th3beta2 + K13L*(PDstandard4th1beta1 + PDstandard4th3beta3) + 
        beta3L*PDstandard4th3K13 + alphaL*(-2*(K11L*Km13 + K12L*Km23 + K13L*Km33) + R13 + K13L*trK);
    
    K22rhsL  =  G122*PDstandard4th1alpha + beta1L*PDstandard4th1K22 - PDstandard4th22alpha + G222*PDstandard4th2alpha + 
        2*(K12L*PDstandard4th2beta1 + K22L*PDstandard4th2beta2 + K23L*PDstandard4th2beta3) + beta2L*PDstandard4th2K22 + 
        G322*PDstandard4th3alpha + beta3L*PDstandard4th3K22 + 
        alphaL*(-2*(K12L*Km12 + K22L*Km22 + K23L*Km32) + R22 + K22L*trK);
    
    K23rhsL  =  G123*PDstandard4th1alpha + beta1L*PDstandard4th1K23 - PDstandard4th23alpha + G223*PDstandard4th2alpha + 
        K13L*PDstandard4th2beta1 + K33L*PDstandard4th2beta3 + beta2L*PDstandard4th2K23 + G323*PDstandard4th3alpha + 
        K12L*PDstandard4th3beta1 + K22L*PDstandard4th3beta2 + K23L*(PDstandard4th2beta2 + PDstandard4th3beta3) + 
        beta3L*PDstandard4th3K23 + alphaL*(-2*(K12L*Km13 + K22L*Km23 + K23L*Km33) + R23 + K23L*trK);
    
    K33rhsL  =  G133*PDstandard4th1alpha + beta1L*PDstandard4th1K33 + G233*PDstandard4th2alpha + beta2L*PDstandard4th2K33 - 
        PDstandard4th33alpha + G333*PDstandard4th3alpha + 
        2*(K13L*PDstandard4th3beta1 + K23L*PDstandard4th3beta2 + K33L*PDstandard4th3beta3) + beta3L*PDstandard4th3K33 + 
        alphaL*(-2*(K13L*Km13 + K23L*Km23 + K33L*Km33) + R33 + K33L*trK);
    
    alpharhsL  =  0;
    
    beta1rhsL  =  0;
    
    beta2rhsL  =  0;
    
    beta3rhsL  =  0;
    
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    g11rhs[index] = g11rhsL;
    g12rhs[index] = g12rhsL;
    g13rhs[index] = g13rhsL;
    g22rhs[index] = g22rhsL;
    g23rhs[index] = g23rhsL;
    g33rhs[index] = g33rhsL;
    K11rhs[index] = K11rhsL;
    K12rhs[index] = K12rhsL;
    K13rhs[index] = K13rhsL;
    K22rhs[index] = K22rhsL;
    K23rhs[index] = K23rhsL;
    K33rhs[index] = K33rhsL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (somename);
}

void ML_ADM_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADM_RHS_Body);
}
