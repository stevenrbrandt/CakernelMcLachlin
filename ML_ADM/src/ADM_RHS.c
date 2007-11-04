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

void ADM_RHS_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ADM_RHS_Body");
  }
  
  if (cctk_iteration % ADM_RHS_calc_every != ADM_RHS_calc_offset)
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
            max[0]-min[0],max[1]-min[1],max[2]-min[2])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL detg = INITVALUE;
    CCTK_REAL G111 = INITVALUE, G121 = INITVALUE, G122 = INITVALUE, G131 = INITVALUE, G132 = INITVALUE, G133 = INITVALUE;
    CCTK_REAL G211 = INITVALUE, G221 = INITVALUE, G222 = INITVALUE, G231 = INITVALUE, G232 = INITVALUE, G233 = INITVALUE;
    CCTK_REAL G311 = INITVALUE, G321 = INITVALUE, G322 = INITVALUE, G331 = INITVALUE, G332 = INITVALUE, G333 = INITVALUE;
    CCTK_REAL gu11 = INITVALUE, gu12 = INITVALUE, gu13 = INITVALUE, gu22 = INITVALUE, gu23 = INITVALUE, gu33 = INITVALUE;
    CCTK_REAL Km11 = INITVALUE, Km12 = INITVALUE, Km13 = INITVALUE, Km21 = INITVALUE, Km22 = INITVALUE, Km23 = INITVALUE;
    CCTK_REAL Km31 = INITVALUE, Km32 = INITVALUE, Km33 = INITVALUE;
    CCTK_REAL R11 = INITVALUE, R21 = INITVALUE, R22 = INITVALUE, R31 = INITVALUE, R32 = INITVALUE, R33 = INITVALUE;
    CCTK_REAL trK = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL alphaL = INITVALUE, alpharhsL = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta1rhsL = INITVALUE, beta2L = INITVALUE, beta2rhsL = INITVALUE, beta3L = INITVALUE, beta3rhsL = INITVALUE;
    CCTK_REAL g11L = INITVALUE, g11rhsL = INITVALUE, g21L = INITVALUE, g21rhsL = INITVALUE, g22L = INITVALUE, g22rhsL = INITVALUE;
    CCTK_REAL g31L = INITVALUE, g31rhsL = INITVALUE, g32L = INITVALUE, g32rhsL = INITVALUE, g33L = INITVALUE, g33rhsL = INITVALUE;
    CCTK_REAL K11L = INITVALUE, K11rhsL = INITVALUE, K21L = INITVALUE, K21rhsL = INITVALUE, K22L = INITVALUE, K22rhsL = INITVALUE;
    CCTK_REAL K31L = INITVALUE, K31rhsL = INITVALUE, K32L = INITVALUE, K32rhsL = INITVALUE, K33L = INITVALUE, K33rhsL = INITVALUE;
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
    CCTK_REAL PDstandard4th1g21 = INITVALUE;
    CCTK_REAL PDstandard4th2g21 = INITVALUE;
    CCTK_REAL PDstandard4th3g21 = INITVALUE;
    CCTK_REAL PDstandard4th33g21 = INITVALUE;
    CCTK_REAL PDstandard4th12g21 = INITVALUE;
    CCTK_REAL PDstandard4th13g21 = INITVALUE;
    CCTK_REAL PDstandard4th21g21 = INITVALUE;
    CCTK_REAL PDstandard4th23g21 = INITVALUE;
    CCTK_REAL PDstandard4th31g21 = INITVALUE;
    CCTK_REAL PDstandard4th32g21 = INITVALUE;
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
    CCTK_REAL PDstandard4th1g31 = INITVALUE;
    CCTK_REAL PDstandard4th2g31 = INITVALUE;
    CCTK_REAL PDstandard4th3g31 = INITVALUE;
    CCTK_REAL PDstandard4th22g31 = INITVALUE;
    CCTK_REAL PDstandard4th12g31 = INITVALUE;
    CCTK_REAL PDstandard4th13g31 = INITVALUE;
    CCTK_REAL PDstandard4th21g31 = INITVALUE;
    CCTK_REAL PDstandard4th23g31 = INITVALUE;
    CCTK_REAL PDstandard4th31g31 = INITVALUE;
    CCTK_REAL PDstandard4th32g31 = INITVALUE;
    CCTK_REAL PDstandard4th1g32 = INITVALUE;
    CCTK_REAL PDstandard4th2g32 = INITVALUE;
    CCTK_REAL PDstandard4th3g32 = INITVALUE;
    CCTK_REAL PDstandard4th11g32 = INITVALUE;
    CCTK_REAL PDstandard4th12g32 = INITVALUE;
    CCTK_REAL PDstandard4th13g32 = INITVALUE;
    CCTK_REAL PDstandard4th21g32 = INITVALUE;
    CCTK_REAL PDstandard4th23g32 = INITVALUE;
    CCTK_REAL PDstandard4th31g32 = INITVALUE;
    CCTK_REAL PDstandard4th32g32 = INITVALUE;
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
    CCTK_REAL PDstandard4th1K21 = INITVALUE;
    CCTK_REAL PDstandard4th2K21 = INITVALUE;
    CCTK_REAL PDstandard4th3K21 = INITVALUE;
    CCTK_REAL PDstandard4th1K22 = INITVALUE;
    CCTK_REAL PDstandard4th2K22 = INITVALUE;
    CCTK_REAL PDstandard4th3K22 = INITVALUE;
    CCTK_REAL PDstandard4th1K31 = INITVALUE;
    CCTK_REAL PDstandard4th2K31 = INITVALUE;
    CCTK_REAL PDstandard4th3K31 = INITVALUE;
    CCTK_REAL PDstandard4th1K32 = INITVALUE;
    CCTK_REAL PDstandard4th2K32 = INITVALUE;
    CCTK_REAL PDstandard4th3K32 = INITVALUE;
    CCTK_REAL PDstandard4th1K33 = INITVALUE;
    CCTK_REAL PDstandard4th2K33 = INITVALUE;
    CCTK_REAL PDstandard4th3K33 = INITVALUE;
    
    /* Assign local copies of grid functions */
    alphaL = alpha[index];
    beta1L = beta1[index];
    beta2L = beta2[index];
    beta3L = beta3[index];
    g11L = g11[index];
    g21L = g21[index];
    g22L = g22[index];
    g31L = g31[index];
    g32L = g32[index];
    g33L = g33[index];
    K11L = K11[index];
    K21L = K21[index];
    K22L = K22[index];
    K31L = K31[index];
    K32L = K32[index];
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
    PDstandard4th1g21 = PDstandard4th1(g21, i, j, k);
    PDstandard4th2g21 = PDstandard4th2(g21, i, j, k);
    PDstandard4th3g21 = PDstandard4th3(g21, i, j, k);
    PDstandard4th33g21 = PDstandard4th33(g21, i, j, k);
    PDstandard4th12g21 = PDstandard4th12(g21, i, j, k);
    PDstandard4th13g21 = PDstandard4th13(g21, i, j, k);
    PDstandard4th23g21 = PDstandard4th23(g21, i, j, k);
    PDstandard4th1g22 = PDstandard4th1(g22, i, j, k);
    PDstandard4th2g22 = PDstandard4th2(g22, i, j, k);
    PDstandard4th3g22 = PDstandard4th3(g22, i, j, k);
    PDstandard4th11g22 = PDstandard4th11(g22, i, j, k);
    PDstandard4th33g22 = PDstandard4th33(g22, i, j, k);
    PDstandard4th13g22 = PDstandard4th13(g22, i, j, k);
    PDstandard4th1g31 = PDstandard4th1(g31, i, j, k);
    PDstandard4th2g31 = PDstandard4th2(g31, i, j, k);
    PDstandard4th3g31 = PDstandard4th3(g31, i, j, k);
    PDstandard4th22g31 = PDstandard4th22(g31, i, j, k);
    PDstandard4th12g31 = PDstandard4th12(g31, i, j, k);
    PDstandard4th13g31 = PDstandard4th13(g31, i, j, k);
    PDstandard4th23g31 = PDstandard4th23(g31, i, j, k);
    PDstandard4th1g32 = PDstandard4th1(g32, i, j, k);
    PDstandard4th2g32 = PDstandard4th2(g32, i, j, k);
    PDstandard4th3g32 = PDstandard4th3(g32, i, j, k);
    PDstandard4th11g32 = PDstandard4th11(g32, i, j, k);
    PDstandard4th12g32 = PDstandard4th12(g32, i, j, k);
    PDstandard4th13g32 = PDstandard4th13(g32, i, j, k);
    PDstandard4th23g32 = PDstandard4th23(g32, i, j, k);
    PDstandard4th1g33 = PDstandard4th1(g33, i, j, k);
    PDstandard4th2g33 = PDstandard4th2(g33, i, j, k);
    PDstandard4th3g33 = PDstandard4th3(g33, i, j, k);
    PDstandard4th11g33 = PDstandard4th11(g33, i, j, k);
    PDstandard4th22g33 = PDstandard4th22(g33, i, j, k);
    PDstandard4th12g33 = PDstandard4th12(g33, i, j, k);
    PDstandard4th1K11 = PDstandard4th1(K11, i, j, k);
    PDstandard4th2K11 = PDstandard4th2(K11, i, j, k);
    PDstandard4th3K11 = PDstandard4th3(K11, i, j, k);
    PDstandard4th1K21 = PDstandard4th1(K21, i, j, k);
    PDstandard4th2K21 = PDstandard4th2(K21, i, j, k);
    PDstandard4th3K21 = PDstandard4th3(K21, i, j, k);
    PDstandard4th1K22 = PDstandard4th1(K22, i, j, k);
    PDstandard4th2K22 = PDstandard4th2(K22, i, j, k);
    PDstandard4th3K22 = PDstandard4th3(K22, i, j, k);
    PDstandard4th1K31 = PDstandard4th1(K31, i, j, k);
    PDstandard4th2K31 = PDstandard4th2(K31, i, j, k);
    PDstandard4th3K31 = PDstandard4th3(K31, i, j, k);
    PDstandard4th1K32 = PDstandard4th1(K32, i, j, k);
    PDstandard4th2K32 = PDstandard4th2(K32, i, j, k);
    PDstandard4th3K32 = PDstandard4th3(K32, i, j, k);
    PDstandard4th1K33 = PDstandard4th1(K33, i, j, k);
    PDstandard4th2K33 = PDstandard4th2(K33, i, j, k);
    PDstandard4th3K33 = PDstandard4th3(K33, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    detg  =  2*g21L*g31L*g32L + g33L*(g11L*g22L - SQR(g21L)) - g22L*SQR(g31L) - g11L*SQR(g32L);
    
    gu11  =  INV(detg)*(g22L*g33L - SQR(g32L));
    
    gu12  =  (g31L*g32L - g21L*g33L)*INV(detg);
    
    gu13  =  (-(g22L*g31L) + g21L*g32L)*INV(detg);
    
    gu22  =  INV(detg)*(g11L*g33L - SQR(g31L));
    
    gu23  =  (g21L*g31L - g11L*g32L)*INV(detg);
    
    gu33  =  INV(detg)*(g11L*g22L - SQR(g21L));
    
    G111  =  khalf*(gu11*PDstandard4th1g11 + 2*(gu12*PDstandard4th1g21 + gu13*PDstandard4th1g31) - gu12*PDstandard4th2g11 - 
          gu13*PDstandard4th3g11);
    
    G211  =  khalf*(gu12*PDstandard4th1g11 + 2*(gu22*PDstandard4th1g21 + gu23*PDstandard4th1g31) - gu22*PDstandard4th2g11 - 
          gu23*PDstandard4th3g11);
    
    G311  =  khalf*(gu13*PDstandard4th1g11 + 2*(gu23*PDstandard4th1g21 + gu33*PDstandard4th1g31) - gu23*PDstandard4th2g11 - 
          gu33*PDstandard4th3g11);
    
    G121  =  khalf*(gu12*PDstandard4th1g22 + gu11*PDstandard4th2g11 + 
          gu13*(PDstandard4th1g32 + PDstandard4th2g31 - PDstandard4th3g21));
    
    G221  =  khalf*(gu22*PDstandard4th1g22 + gu12*PDstandard4th2g11 + 
          gu23*(PDstandard4th1g32 + PDstandard4th2g31 - PDstandard4th3g21));
    
    G321  =  khalf*(gu23*PDstandard4th1g22 + gu13*PDstandard4th2g11 + 
          gu33*(PDstandard4th1g32 + PDstandard4th2g31 - PDstandard4th3g21));
    
    G131  =  khalf*(gu13*PDstandard4th1g33 + gu11*PDstandard4th3g11 + 
          gu12*(PDstandard4th1g32 - PDstandard4th2g31 + PDstandard4th3g21));
    
    G231  =  khalf*(gu23*PDstandard4th1g33 + gu12*PDstandard4th3g11 + 
          gu22*(PDstandard4th1g32 - PDstandard4th2g31 + PDstandard4th3g21));
    
    G331  =  khalf*(gu33*PDstandard4th1g33 + gu13*PDstandard4th3g11 + 
          gu23*(PDstandard4th1g32 - PDstandard4th2g31 + PDstandard4th3g21));
    
    G122  =  khalf*(gu11*(-PDstandard4th1g22 + 2*PDstandard4th2g21) + gu12*PDstandard4th2g22 + 
          gu13*(2*PDstandard4th2g32 - PDstandard4th3g22));
    
    G222  =  khalf*(gu12*(-PDstandard4th1g22 + 2*PDstandard4th2g21) + gu22*PDstandard4th2g22 + 
          gu23*(2*PDstandard4th2g32 - PDstandard4th3g22));
    
    G322  =  khalf*(gu13*(-PDstandard4th1g22 + 2*PDstandard4th2g21) + gu23*PDstandard4th2g22 + 
          gu33*(2*PDstandard4th2g32 - PDstandard4th3g22));
    
    G132  =  khalf*(gu13*PDstandard4th2g33 + gu11*(-PDstandard4th1g32 + PDstandard4th2g31 + PDstandard4th3g21) + 
          gu12*PDstandard4th3g22);
    
    G232  =  khalf*(gu23*PDstandard4th2g33 + gu12*(-PDstandard4th1g32 + PDstandard4th2g31 + PDstandard4th3g21) + 
          gu22*PDstandard4th3g22);
    
    G332  =  khalf*(gu33*PDstandard4th2g33 + gu13*(-PDstandard4th1g32 + PDstandard4th2g31 + PDstandard4th3g21) + 
          gu23*PDstandard4th3g22);
    
    G133  =  khalf*(-(gu11*PDstandard4th1g33) - gu12*PDstandard4th2g33 + 2*gu11*PDstandard4th3g31 + 
          2*gu12*PDstandard4th3g32 + gu13*PDstandard4th3g33);
    
    G233  =  khalf*(-(gu12*PDstandard4th1g33) - gu22*PDstandard4th2g33 + 2*gu12*PDstandard4th3g31 + 
          2*gu22*PDstandard4th3g32 + gu23*PDstandard4th3g33);
    
    G333  =  khalf*(-(gu13*PDstandard4th1g33) - gu23*PDstandard4th2g33 + 2*gu13*PDstandard4th3g31 + 
          2*gu23*PDstandard4th3g32 + gu33*PDstandard4th3g33);
    
    R11  =  khalf*(4*G231*G321 + G211*(2*G121 - 2*G222 - 2*G332) + G311*(2*G131 - 2*G333) - gu22*PDstandard4th11g22 - 
          2*(G111*G221 + G232*G311 + G111*G331 + gu23*PDstandard4th11g32) - gu33*PDstandard4th11g33 + 
          2*gu22*PDstandard4th12g21 + 2*gu23*PDstandard4th12g31 + 2*gu23*PDstandard4th13g21 + 2*gu33*PDstandard4th13g31 - 
          gu22*PDstandard4th22g11 - 2*gu23*PDstandard4th23g11 - gu33*PDstandard4th33g11 + 2*SQR(G221) + 2*SQR(G331));
    
    R21  =  khalf*(2*(G122*G211 + G132*G311 + G231*G322 + G331*G332) - 
          2*(G121*G221 + G121*G331 + G221*G332 + G321*G333 + gu12*PDstandard4th12g21) - gu23*PDstandard4th12g32 - 
          gu33*PDstandard4th12g33 + gu13*(PDstandard4th11g32 - PDstandard4th12g31 - PDstandard4th13g21) + 
          gu23*PDstandard4th13g22 + gu33*PDstandard4th13g32 + gu12*(PDstandard4th11g22 + PDstandard4th22g11) + 
          gu23*PDstandard4th22g31 + gu13*PDstandard4th23g11 - gu23*PDstandard4th23g21 + gu33*PDstandard4th23g31 - 
          gu33*PDstandard4th33g21);
    
    R31  =  khalf*(2*(G132*G211 + G221*G232 + G133*G311 + G233*G321) - 
          2*(G222*G231 + G232*G331 + G131*(G221 + G331) + gu13*PDstandard4th13g31) + 
          gu12*(PDstandard4th11g32 - PDstandard4th12g31 - PDstandard4th13g21 + PDstandard4th23g11) + 
          gu22*(PDstandard4th12g32 - PDstandard4th13g22 - PDstandard4th22g31 + PDstandard4th23g21) + 
          gu13*(PDstandard4th11g33 + PDstandard4th33g11) + 
          gu23*(PDstandard4th12g33 - PDstandard4th13g32 - PDstandard4th23g31 + PDstandard4th33g21));
    
    R22  =  khalf*(4*G132*G321 + G122*(2*G221 - 2*(G111 + G331)) + G322*(2*G232 - 2*G333) - 
          2*(G131*G322 + G222*(G121 + G332) + gu13*PDstandard4th13g22) + 
          gu11*(-PDstandard4th11g22 + 2*PDstandard4th12g21 - PDstandard4th22g11) + 
          gu13*(-2*PDstandard4th22g31 + 2*(PDstandard4th12g32 + PDstandard4th23g21)) + 
          gu33*(-PDstandard4th22g33 + 2*PDstandard4th23g32 - PDstandard4th33g22) + 2*(SQR(G121) + SQR(G332)));
    
    R32  =  khalf*(2*(G121*G131 + G122*G231 + G133*G321 + G233*G322) + 
          gu11*(-PDstandard4th11g32 + PDstandard4th12g31 + PDstandard4th13g21 - PDstandard4th23g11) + 
          gu12*(-PDstandard4th12g32 + PDstandard4th13g22 + PDstandard4th22g31 - PDstandard4th23g21) - 
          2*(G111*G132 + G131*G332 + G232*(G121 + G332) + gu23*PDstandard4th23g32) + 
          gu13*(PDstandard4th12g33 - PDstandard4th13g32 - PDstandard4th23g31 + PDstandard4th33g21) + 
          gu23*(PDstandard4th22g33 + PDstandard4th33g22));
    
    R33  =  khalf*(4*G132*G231 - gu11*PDstandard4th11g33 - 
          2*(G111*G133 + G133*G221 + G121*G233 + G222*G233 + G131*G333 + G232*G333 + gu12*PDstandard4th12g33) + 
          2*(G133*G331 + G233*G332 + gu11*PDstandard4th13g31) + 2*gu12*PDstandard4th13g32 - gu22*PDstandard4th22g33 + 
          2*gu12*PDstandard4th23g31 + 2*gu22*PDstandard4th23g32 - gu11*PDstandard4th33g11 - 2*gu12*PDstandard4th33g21 - 
          gu22*PDstandard4th33g22 + 2*SQR(G131) + 2*SQR(G232));
    
    Km11  =  gu11*K11L + gu12*K21L + gu13*K31L;
    
    Km21  =  gu12*K11L + gu22*K21L + gu23*K31L;
    
    Km31  =  gu13*K11L + gu23*K21L + gu33*K31L;
    
    Km12  =  gu11*K21L + gu12*K22L + gu13*K32L;
    
    Km22  =  gu12*K21L + gu22*K22L + gu23*K32L;
    
    Km32  =  gu13*K21L + gu23*K22L + gu33*K32L;
    
    Km13  =  gu11*K31L + gu12*K32L + gu13*K33L;
    
    Km23  =  gu12*K31L + gu22*K32L + gu23*K33L;
    
    Km33  =  gu13*K31L + gu23*K32L + gu33*K33L;
    
    trK  =  Km11 + Km22 + Km33;
    
    g11rhsL  =  -2*alphaL*K11L + 2*(g11L*PDstandard4th1beta1 + g21L*PDstandard4th1beta2 + g31L*PDstandard4th1beta3) + 
        beta1L*PDstandard4th1g11 + beta2L*PDstandard4th2g11 + beta3L*PDstandard4th3g11;
    
    g21rhsL  =  -2*alphaL*K21L + g22L*PDstandard4th1beta2 + g32L*PDstandard4th1beta3 + beta1L*PDstandard4th1g21 + 
        g11L*PDstandard4th2beta1 + g21L*(PDstandard4th1beta1 + PDstandard4th2beta2) + g31L*PDstandard4th2beta3 + 
        beta2L*PDstandard4th2g21 + beta3L*PDstandard4th3g21;
    
    g31rhsL  =  -2*alphaL*K31L + g32L*PDstandard4th1beta2 + g33L*PDstandard4th1beta3 + beta1L*PDstandard4th1g31 + 
        beta2L*PDstandard4th2g31 + g11L*PDstandard4th3beta1 + g21L*PDstandard4th3beta2 + 
        g31L*(PDstandard4th1beta1 + PDstandard4th3beta3) + beta3L*PDstandard4th3g31;
    
    g22rhsL  =  -2*alphaL*K22L + beta1L*PDstandard4th1g22 + 
        2*(g21L*PDstandard4th2beta1 + g22L*PDstandard4th2beta2 + g32L*PDstandard4th2beta3) + beta2L*PDstandard4th2g22 + 
        beta3L*PDstandard4th3g22;
    
    g32rhsL  =  -2*alphaL*K32L + beta1L*PDstandard4th1g32 + g31L*PDstandard4th2beta1 + g33L*PDstandard4th2beta3 + 
        beta2L*PDstandard4th2g32 + g21L*PDstandard4th3beta1 + g22L*PDstandard4th3beta2 + 
        g32L*(PDstandard4th2beta2 + PDstandard4th3beta3) + beta3L*PDstandard4th3g32;
    
    g33rhsL  =  -2*alphaL*K33L + beta1L*PDstandard4th1g33 + beta2L*PDstandard4th2g33 + 
        2*(g31L*PDstandard4th3beta1 + g32L*PDstandard4th3beta2 + g33L*PDstandard4th3beta3) + beta3L*PDstandard4th3g33;
    
    K11rhsL  =  -PDstandard4th11alpha + G111*PDstandard4th1alpha + 
        2*(K11L*PDstandard4th1beta1 + K21L*PDstandard4th1beta2 + K31L*PDstandard4th1beta3) + beta1L*PDstandard4th1K11 + 
        G211*PDstandard4th2alpha + beta2L*PDstandard4th2K11 + G311*PDstandard4th3alpha + beta3L*PDstandard4th3K11 + 
        alphaL*(-2*(K11L*Km11 + K21L*Km21 + K31L*Km31) + R11 + K11L*trK);
    
    K21rhsL  =  -PDstandard4th12alpha + G121*PDstandard4th1alpha + K22L*PDstandard4th1beta2 + K32L*PDstandard4th1beta3 + 
        beta1L*PDstandard4th1K21 + G221*PDstandard4th2alpha + K11L*PDstandard4th2beta1 + 
        K21L*(PDstandard4th1beta1 + PDstandard4th2beta2) + K31L*PDstandard4th2beta3 + beta2L*PDstandard4th2K21 + 
        G321*PDstandard4th3alpha + beta3L*PDstandard4th3K21 + 
        alphaL*(-2*(K11L*Km12 + K21L*Km22 + K31L*Km32) + R21 + K21L*trK);
    
    K31rhsL  =  -PDstandard4th13alpha + G131*PDstandard4th1alpha + K32L*PDstandard4th1beta2 + K33L*PDstandard4th1beta3 + 
        beta1L*PDstandard4th1K31 + G231*PDstandard4th2alpha + beta2L*PDstandard4th2K31 + G331*PDstandard4th3alpha + 
        K11L*PDstandard4th3beta1 + K21L*PDstandard4th3beta2 + K31L*(PDstandard4th1beta1 + PDstandard4th3beta3) + 
        beta3L*PDstandard4th3K31 + alphaL*(-2*(K11L*Km13 + K21L*Km23 + K31L*Km33) + R31 + K31L*trK);
    
    K22rhsL  =  G122*PDstandard4th1alpha + beta1L*PDstandard4th1K22 - PDstandard4th22alpha + G222*PDstandard4th2alpha + 
        2*(K21L*PDstandard4th2beta1 + K22L*PDstandard4th2beta2 + K32L*PDstandard4th2beta3) + beta2L*PDstandard4th2K22 + 
        G322*PDstandard4th3alpha + beta3L*PDstandard4th3K22 + 
        alphaL*(-2*(K21L*Km12 + K22L*Km22 + K32L*Km32) + R22 + K22L*trK);
    
    K32rhsL  =  G132*PDstandard4th1alpha + beta1L*PDstandard4th1K32 - PDstandard4th23alpha + G232*PDstandard4th2alpha + 
        K31L*PDstandard4th2beta1 + K33L*PDstandard4th2beta3 + beta2L*PDstandard4th2K32 + G332*PDstandard4th3alpha + 
        K21L*PDstandard4th3beta1 + K22L*PDstandard4th3beta2 + K32L*(PDstandard4th2beta2 + PDstandard4th3beta3) + 
        beta3L*PDstandard4th3K32 + alphaL*(-2*(K21L*Km13 + K22L*Km23 + K32L*Km33) + R32 + K32L*trK);
    
    K33rhsL  =  G133*PDstandard4th1alpha + beta1L*PDstandard4th1K33 + G233*PDstandard4th2alpha + beta2L*PDstandard4th2K33 - 
        PDstandard4th33alpha + G333*PDstandard4th3alpha + 
        2*(K31L*PDstandard4th3beta1 + K32L*PDstandard4th3beta2 + K33L*PDstandard4th3beta3) + beta3L*PDstandard4th3K33 + 
        alphaL*(-2*(K31L*Km13 + K32L*Km23 + K33L*Km33) + R33 + K33L*trK);
    
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
    g21rhs[index] = g21rhsL;
    g22rhs[index] = g22rhsL;
    g31rhs[index] = g31rhsL;
    g32rhs[index] = g32rhsL;
    g33rhs[index] = g33rhsL;
    K11rhs[index] = K11rhsL;
    K21rhs[index] = K21rhsL;
    K22rhs[index] = K22rhsL;
    K31rhs[index] = K31rhsL;
    K32rhs[index] = K32rhsL;
    K33rhs[index] = K33rhsL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (somename);
}

void ADM_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverEverything(cctkGH, &ADM_RHS_Body);
}
