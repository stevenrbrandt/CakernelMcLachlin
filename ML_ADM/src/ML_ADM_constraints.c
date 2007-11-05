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

void ML_ADM_constraints_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADM_constraints_Body");
  }
  
  if (cctk_iteration % ML_ADM_constraints_calc_every != ML_ADM_constraints_calc_offset)
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
    CCTK_REAL trR = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL g11L = INITVALUE, g12L = INITVALUE, g13L = INITVALUE, g22L = INITVALUE, g23L = INITVALUE, g33L = INITVALUE;
    CCTK_REAL HL = INITVALUE;
    CCTK_REAL K11L = INITVALUE, K12L = INITVALUE, K13L = INITVALUE, K22L = INITVALUE, K23L = INITVALUE, K33L = INITVALUE;
    CCTK_REAL M1L = INITVALUE, M2L = INITVALUE, M3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
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
    CCTK_REAL PDstandard4th2K11 = INITVALUE;
    CCTK_REAL PDstandard4th3K11 = INITVALUE;
    CCTK_REAL PDstandard4th1K12 = INITVALUE;
    CCTK_REAL PDstandard4th2K12 = INITVALUE;
    CCTK_REAL PDstandard4th3K12 = INITVALUE;
    CCTK_REAL PDstandard4th1K13 = INITVALUE;
    CCTK_REAL PDstandard4th2K13 = INITVALUE;
    CCTK_REAL PDstandard4th3K13 = INITVALUE;
    CCTK_REAL PDstandard4th1K22 = INITVALUE;
    CCTK_REAL PDstandard4th3K22 = INITVALUE;
    CCTK_REAL PDstandard4th1K23 = INITVALUE;
    CCTK_REAL PDstandard4th2K23 = INITVALUE;
    CCTK_REAL PDstandard4th3K23 = INITVALUE;
    CCTK_REAL PDstandard4th1K33 = INITVALUE;
    CCTK_REAL PDstandard4th2K33 = INITVALUE;
    
    /* Assign local copies of grid functions */
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
    PDstandard4th2K11 = PDstandard4th2(K11, i, j, k);
    PDstandard4th3K11 = PDstandard4th3(K11, i, j, k);
    PDstandard4th1K12 = PDstandard4th1(K12, i, j, k);
    PDstandard4th2K12 = PDstandard4th2(K12, i, j, k);
    PDstandard4th3K12 = PDstandard4th3(K12, i, j, k);
    PDstandard4th1K13 = PDstandard4th1(K13, i, j, k);
    PDstandard4th2K13 = PDstandard4th2(K13, i, j, k);
    PDstandard4th3K13 = PDstandard4th3(K13, i, j, k);
    PDstandard4th1K22 = PDstandard4th1(K22, i, j, k);
    PDstandard4th3K22 = PDstandard4th3(K22, i, j, k);
    PDstandard4th1K23 = PDstandard4th1(K23, i, j, k);
    PDstandard4th2K23 = PDstandard4th2(K23, i, j, k);
    PDstandard4th3K23 = PDstandard4th3(K23, i, j, k);
    PDstandard4th1K33 = PDstandard4th1(K33, i, j, k);
    PDstandard4th2K33 = PDstandard4th2(K33, i, j, k);
    
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
    
    R11  =  khalf*(-(gu22*PDstandard4th11g22) - 2*(G111*G122 + G111*G133 + G211*G222 + G211*G233 + G311*G322 + G311*G333 + 
             gu32*PDstandard4th11g23) - gu33*PDstandard4th11g33 + 2*gu22*PDstandard4th12g12 + 2*gu32*PDstandard4th12g13 + 
          2*gu32*PDstandard4th13g12 + 2*gu33*PDstandard4th13g13 - gu22*PDstandard4th22g11 - 2*gu32*PDstandard4th23g11 - 
          gu33*PDstandard4th33g11 + 2*SQR(G112) + 2*SQR(G113) + 2*SQR(G212) + 2*SQR(G213) + 2*SQR(G312) + 2*SQR(G313));
    
    R12  =  khalf*(2*(G113*G123 + G213*G223 + G313*G323) - 2*(G112*G133 + G212*G233 + G312*G333 + gu21*PDstandard4th12g12) - 
          gu32*PDstandard4th12g23 - gu33*PDstandard4th12g33 + 
          gu31*(PDstandard4th11g23 - PDstandard4th12g13 - PDstandard4th13g12) + gu32*PDstandard4th13g22 + 
          gu33*PDstandard4th13g23 + gu21*(PDstandard4th11g22 + PDstandard4th22g11) + gu32*PDstandard4th22g13 + 
          gu31*PDstandard4th23g11 - gu32*PDstandard4th23g12 + gu33*PDstandard4th23g13 - gu33*PDstandard4th33g12);
    
    R13  =  khalf*(2*(G112*G123 + G212*G223 + G312*G323) - 2*(G113*G122 + G213*G222 + G313*G322 + gu31*PDstandard4th13g13) + 
          gu21*(PDstandard4th11g23 - PDstandard4th12g13 - PDstandard4th13g12 + PDstandard4th23g11) + 
          gu22*(PDstandard4th12g23 - PDstandard4th13g22 - PDstandard4th22g13 + PDstandard4th23g12) + 
          gu31*(PDstandard4th11g33 + PDstandard4th33g11) + 
          gu32*(PDstandard4th12g33 - PDstandard4th13g23 - PDstandard4th23g13 + PDstandard4th33g12));
    
    R22  =  khalf*(-2*(G122*(G111 + G133) + G222*(G211 + G233) + G322*(G311 + G333) + gu31*PDstandard4th13g22) + 
          gu11*(-PDstandard4th11g22 + 2*PDstandard4th12g12 - PDstandard4th22g11) + 
          gu31*(-2*PDstandard4th22g13 + 2*(PDstandard4th12g23 + PDstandard4th23g12)) + 
          gu33*(-PDstandard4th22g33 + 2*PDstandard4th23g23 - PDstandard4th33g22) + 
          2*(SQR(G112) + SQR(G123) + SQR(G212) + SQR(G223) + SQR(G312) + SQR(G323)));
    
    R23  =  khalf*(2*(G112*G113 + G212*G213 + G312*G313) + 
          gu11*(-PDstandard4th11g23 + PDstandard4th12g13 + PDstandard4th13g12 - PDstandard4th23g11) + 
          gu21*(-PDstandard4th12g23 + PDstandard4th13g22 + PDstandard4th22g13 - PDstandard4th23g12) - 
          2*(G111*G123 + G211*G223 + G311*G323 + gu32*PDstandard4th23g23) + 
          gu31*(PDstandard4th12g33 - PDstandard4th13g23 - PDstandard4th23g13 + PDstandard4th33g12) + 
          gu32*(PDstandard4th22g33 + PDstandard4th33g22));
    
    R33  =  khalf*(gu11*(-PDstandard4th11g33 + 2*PDstandard4th13g13 - PDstandard4th33g11) - 
          2*((G111 + G122)*G133 + (G211 + G222)*G233 + (G311 + G322)*G333 + 
             gu21*(PDstandard4th12g33 + PDstandard4th33g12)) + 
          gu22*(-PDstandard4th22g33 + 2*PDstandard4th23g23 - PDstandard4th33g22) + 
          2*(gu21*(PDstandard4th13g23 + PDstandard4th23g13) + SQR(G113) + SQR(G123) + SQR(G213) + SQR(G223) + SQR(G313) + 
             SQR(G323)));
    
    trR  =  gu11*R11 + gu22*R22 + 2*(gu21*R12 + gu31*R13 + gu32*R23) + gu33*R33;
    
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
    
    HL  =  -2*(Km12*Km21 + Km13*Km31 + Km23*Km32) + trR - SQR(Km11) - SQR(Km22) - SQR(Km33) + SQR(trK);
    
    M1L  =  gu21*(-(G112*K11L) + G111*K12L - G212*K12L - G312*K13L + G211*K22L + G311*K23L - PDstandard4th1K12 + 
           PDstandard4th2K11) + gu22*(-(G122*K11L) + G112*K12L - G222*K12L - G322*K13L + G212*K22L + G312*K23L - 
           PDstandard4th1K22 + PDstandard4th2K12) + gu31*
         (-(G113*K11L) - G213*K12L + G111*K13L - G313*K13L + G211*K23L + G311*K33L - PDstandard4th1K13 + PDstandard4th3K11)\
         + gu32*(G113*K12L + G112*K13L + G213*K22L + (G212 + G313)*K23L + G312*K33L - 
           2*(G123*K11L + G223*K12L + G323*K13L + PDstandard4th1K23) + PDstandard4th2K13 + PDstandard4th3K12) + 
        gu33*(-(G133*K11L) - G233*K12L + G113*K13L - G333*K13L + G213*K23L + G313*K33L - PDstandard4th1K33 + 
           PDstandard4th3K13);
    
    M2L  =  gu11*(G112*K11L + (-G111 + G212)*K12L + G312*K13L - G211*K22L - G311*K23L + PDstandard4th1K12 - 
           PDstandard4th2K11) + gu21*(G122*K11L + (-G112 + G222)*K12L + G322*K13L - G212*K22L - G312*K23L + 
           PDstandard4th1K22 - PDstandard4th2K12) + gu31*
         (G123*K11L + (-2*G113 + G223)*K12L + (G112 + G323)*K13L + G212*K23L + G312*K33L + PDstandard4th1K23 - 
           2*(G213*K22L + G313*K23L + PDstandard4th2K13) + PDstandard4th3K12) + 
        gu32*(-(G123*K12L) + G122*K13L - G223*K22L + G222*K23L - G323*K23L + G322*K33L - PDstandard4th2K23 + 
           PDstandard4th3K22) + gu33*(-(G133*K12L) + G123*K13L - G233*K22L + G223*K23L - G333*K23L + G323*K33L - 
           PDstandard4th2K33 + PDstandard4th3K23);
    
    M3L  =  gu11*(G113*K11L + G213*K12L + (-G111 + G313)*K13L - G211*K23L - G311*K33L + PDstandard4th1K13 - 
           PDstandard4th3K11) + gu21*(G123*K11L + (G113 + G223)*K12L + (-2*G112 + G323)*K13L + G213*K22L + 
           (-2*G212 + G313)*K23L + PDstandard4th1K23 + PDstandard4th2K13 - 2*(G312*K33L + PDstandard4th3K12)) + 
        gu31*(G133*K11L + G233*K12L + (-G113 + G333)*K13L - G213*K23L - G313*K33L + PDstandard4th1K33 - 
           PDstandard4th3K13) + gu22*(G123*K12L - G122*K13L + G223*K22L - G222*K23L + G323*K23L - G322*K33L + 
           PDstandard4th2K23 - PDstandard4th3K22) + gu32*
         (G133*K12L - G123*K13L + G233*K22L - G223*K23L + G333*K23L - G323*K33L + PDstandard4th2K33 - PDstandard4th3K23);
    
    
    /* Copy local copies back to grid functions */
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (somename);
}

void ML_ADM_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADM_constraints_Body);
}
