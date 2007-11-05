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

void ML_BSSN_RHS_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_RHS_Body");
  }
  
  if (cctk_iteration % ML_BSSN_RHS_calc_every != ML_BSSN_RHS_calc_offset)
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
    CCTK_REAL Atm11 = INITVALUE, Atm12 = INITVALUE, Atm13 = INITVALUE, Atm21 = INITVALUE, Atm22 = INITVALUE, Atm23 = INITVALUE;
    CCTK_REAL Atm31 = INITVALUE, Atm32 = INITVALUE, Atm33 = INITVALUE;
    CCTK_REAL Atu11 = INITVALUE, Atu21 = INITVALUE, Atu22 = INITVALUE, Atu31 = INITVALUE, Atu32 = INITVALUE, Atu33 = INITVALUE;
    CCTK_REAL ddetg1 = INITVALUE, ddetg2 = INITVALUE, ddetg3 = INITVALUE;
    CCTK_REAL ddetgt1 = INITVALUE, ddetgt2 = INITVALUE, ddetgt3 = INITVALUE;
    CCTK_REAL ddgtu1111 = INITVALUE, ddgtu1112 = INITVALUE, ddgtu1113 = INITVALUE, ddgtu2111 = INITVALUE, ddgtu2112 = INITVALUE, ddgtu2113 = INITVALUE;
    CCTK_REAL ddgtu2122 = INITVALUE, ddgtu2123 = INITVALUE, ddgtu2212 = INITVALUE, ddgtu2222 = INITVALUE, ddgtu2223 = INITVALUE, ddgtu3111 = INITVALUE;
    CCTK_REAL ddgtu3112 = INITVALUE, ddgtu3113 = INITVALUE, ddgtu3123 = INITVALUE, ddgtu3133 = INITVALUE, ddgtu3212 = INITVALUE, ddgtu3213 = INITVALUE;
    CCTK_REAL ddgtu3222 = INITVALUE, ddgtu3223 = INITVALUE, ddgtu3233 = INITVALUE, ddgtu3313 = INITVALUE, ddgtu3323 = INITVALUE, ddgtu3333 = INITVALUE;
    CCTK_REAL detg = INITVALUE;
    CCTK_REAL detgt = INITVALUE;
    CCTK_REAL dgtu111 = INITVALUE, dgtu112 = INITVALUE, dgtu113 = INITVALUE, dgtu211 = INITVALUE, dgtu212 = INITVALUE, dgtu213 = INITVALUE;
    CCTK_REAL dgtu221 = INITVALUE, dgtu222 = INITVALUE, dgtu223 = INITVALUE, dgtu311 = INITVALUE, dgtu312 = INITVALUE, dgtu313 = INITVALUE;
    CCTK_REAL dgtu321 = INITVALUE, dgtu322 = INITVALUE, dgtu323 = INITVALUE, dgtu331 = INITVALUE, dgtu332 = INITVALUE, dgtu333 = INITVALUE;
    CCTK_REAL e4phi = INITVALUE;
    CCTK_REAL em4phi = INITVALUE;
    CCTK_REAL g11 = INITVALUE;
    CCTK_REAL G111 = INITVALUE, G112 = INITVALUE, G113 = INITVALUE;
    CCTK_REAL g12 = INITVALUE;
    CCTK_REAL G122 = INITVALUE, G123 = INITVALUE;
    CCTK_REAL g13 = INITVALUE;
    CCTK_REAL G133 = INITVALUE, G211 = INITVALUE, G212 = INITVALUE, G213 = INITVALUE;
    CCTK_REAL g22 = INITVALUE;
    CCTK_REAL G222 = INITVALUE, G223 = INITVALUE;
    CCTK_REAL g23 = INITVALUE;
    CCTK_REAL G233 = INITVALUE, G311 = INITVALUE, G312 = INITVALUE, G313 = INITVALUE, G322 = INITVALUE, G323 = INITVALUE;
    CCTK_REAL g33 = INITVALUE;
    CCTK_REAL G333 = INITVALUE;
    CCTK_REAL Gt111 = INITVALUE, Gt112 = INITVALUE, Gt113 = INITVALUE, Gt122 = INITVALUE, Gt123 = INITVALUE, Gt133 = INITVALUE;
    CCTK_REAL Gt211 = INITVALUE, Gt212 = INITVALUE, Gt213 = INITVALUE, Gt222 = INITVALUE, Gt223 = INITVALUE, Gt233 = INITVALUE;
    CCTK_REAL Gt311 = INITVALUE, Gt312 = INITVALUE, Gt313 = INITVALUE, Gt322 = INITVALUE, Gt323 = INITVALUE, Gt333 = INITVALUE;
    CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    CCTK_REAL R11 = INITVALUE, R12 = INITVALUE, R13 = INITVALUE, R22 = INITVALUE, R23 = INITVALUE, R33 = INITVALUE;
    CCTK_REAL Rphi11 = INITVALUE, Rphi12 = INITVALUE, Rphi13 = INITVALUE, Rphi22 = INITVALUE, Rphi23 = INITVALUE, Rphi33 = INITVALUE;
    CCTK_REAL Rt11 = INITVALUE, Rt12 = INITVALUE, Rt13 = INITVALUE, Rt22 = INITVALUE, Rt23 = INITVALUE, Rt33 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL alphaL = INITVALUE, alpharhsL = INITVALUE;
    CCTK_REAL At11L = INITVALUE, At11rhsL = INITVALUE, At12L = INITVALUE, At12rhsL = INITVALUE, At13L = INITVALUE, At13rhsL = INITVALUE;
    CCTK_REAL At22L = INITVALUE, At22rhsL = INITVALUE, At23L = INITVALUE, At23rhsL = INITVALUE, At33L = INITVALUE, At33rhsL = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta1rhsL = INITVALUE, beta2L = INITVALUE, beta2rhsL = INITVALUE, beta3L = INITVALUE, beta3rhsL = INITVALUE;
    CCTK_REAL dtalphaL = INITVALUE, dtalpharhsL = INITVALUE;
    CCTK_REAL dtbeta1L = INITVALUE, dtbeta1rhsL = INITVALUE, dtbeta2L = INITVALUE, dtbeta2rhsL = INITVALUE, dtbeta3L = INITVALUE, dtbeta3rhsL = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt11rhsL = INITVALUE, gt12L = INITVALUE, gt12rhsL = INITVALUE, gt13L = INITVALUE, gt13rhsL = INITVALUE;
    CCTK_REAL gt22L = INITVALUE, gt22rhsL = INITVALUE, gt23L = INITVALUE, gt23rhsL = INITVALUE, gt33L = INITVALUE, gt33rhsL = INITVALUE;
    CCTK_REAL phiL = INITVALUE, phirhsL = INITVALUE;
    CCTK_REAL trKL = INITVALUE, trKrhsL = INITVALUE;
    CCTK_REAL Xt1L = INITVALUE, Xt1rhsL = INITVALUE, Xt2L = INITVALUE, Xt2rhsL = INITVALUE, Xt3L = INITVALUE, Xt3rhsL = INITVALUE;
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
    CCTK_REAL PDstandard4th21alpha = INITVALUE;
    CCTK_REAL PDstandard4th23alpha = INITVALUE;
    CCTK_REAL PDstandard4th31alpha = INITVALUE;
    CCTK_REAL PDstandard4th32alpha = INITVALUE;
    CCTK_REAL PDstandard4th1At11 = INITVALUE;
    CCTK_REAL PDstandard4th2At11 = INITVALUE;
    CCTK_REAL PDstandard4th3At11 = INITVALUE;
    CCTK_REAL PDstandard4th1At12 = INITVALUE;
    CCTK_REAL PDstandard4th2At12 = INITVALUE;
    CCTK_REAL PDstandard4th3At12 = INITVALUE;
    CCTK_REAL PDstandard4th1At13 = INITVALUE;
    CCTK_REAL PDstandard4th2At13 = INITVALUE;
    CCTK_REAL PDstandard4th3At13 = INITVALUE;
    CCTK_REAL PDstandard4th1At22 = INITVALUE;
    CCTK_REAL PDstandard4th2At22 = INITVALUE;
    CCTK_REAL PDstandard4th3At22 = INITVALUE;
    CCTK_REAL PDstandard4th1At23 = INITVALUE;
    CCTK_REAL PDstandard4th2At23 = INITVALUE;
    CCTK_REAL PDstandard4th3At23 = INITVALUE;
    CCTK_REAL PDstandard4th1At33 = INITVALUE;
    CCTK_REAL PDstandard4th2At33 = INITVALUE;
    CCTK_REAL PDstandard4th3At33 = INITVALUE;
    CCTK_REAL PDstandard4th1beta1 = INITVALUE;
    CCTK_REAL PDstandard4th2beta1 = INITVALUE;
    CCTK_REAL PDstandard4th3beta1 = INITVALUE;
    CCTK_REAL PDstandard4th11beta1 = INITVALUE;
    CCTK_REAL PDstandard4th22beta1 = INITVALUE;
    CCTK_REAL PDstandard4th33beta1 = INITVALUE;
    CCTK_REAL PDstandard4th12beta1 = INITVALUE;
    CCTK_REAL PDstandard4th13beta1 = INITVALUE;
    CCTK_REAL PDstandard4th21beta1 = INITVALUE;
    CCTK_REAL PDstandard4th23beta1 = INITVALUE;
    CCTK_REAL PDstandard4th31beta1 = INITVALUE;
    CCTK_REAL PDstandard4th32beta1 = INITVALUE;
    CCTK_REAL PDstandard4th1beta2 = INITVALUE;
    CCTK_REAL PDstandard4th2beta2 = INITVALUE;
    CCTK_REAL PDstandard4th3beta2 = INITVALUE;
    CCTK_REAL PDstandard4th11beta2 = INITVALUE;
    CCTK_REAL PDstandard4th22beta2 = INITVALUE;
    CCTK_REAL PDstandard4th33beta2 = INITVALUE;
    CCTK_REAL PDstandard4th12beta2 = INITVALUE;
    CCTK_REAL PDstandard4th13beta2 = INITVALUE;
    CCTK_REAL PDstandard4th21beta2 = INITVALUE;
    CCTK_REAL PDstandard4th23beta2 = INITVALUE;
    CCTK_REAL PDstandard4th31beta2 = INITVALUE;
    CCTK_REAL PDstandard4th32beta2 = INITVALUE;
    CCTK_REAL PDstandard4th1beta3 = INITVALUE;
    CCTK_REAL PDstandard4th2beta3 = INITVALUE;
    CCTK_REAL PDstandard4th3beta3 = INITVALUE;
    CCTK_REAL PDstandard4th11beta3 = INITVALUE;
    CCTK_REAL PDstandard4th22beta3 = INITVALUE;
    CCTK_REAL PDstandard4th33beta3 = INITVALUE;
    CCTK_REAL PDstandard4th12beta3 = INITVALUE;
    CCTK_REAL PDstandard4th13beta3 = INITVALUE;
    CCTK_REAL PDstandard4th21beta3 = INITVALUE;
    CCTK_REAL PDstandard4th23beta3 = INITVALUE;
    CCTK_REAL PDstandard4th31beta3 = INITVALUE;
    CCTK_REAL PDstandard4th32beta3 = INITVALUE;
    CCTK_REAL PDstandard4th1gt11 = INITVALUE;
    CCTK_REAL PDstandard4th2gt11 = INITVALUE;
    CCTK_REAL PDstandard4th3gt11 = INITVALUE;
    CCTK_REAL PDstandard4th11gt11 = INITVALUE;
    CCTK_REAL PDstandard4th22gt11 = INITVALUE;
    CCTK_REAL PDstandard4th33gt11 = INITVALUE;
    CCTK_REAL PDstandard4th12gt11 = INITVALUE;
    CCTK_REAL PDstandard4th13gt11 = INITVALUE;
    CCTK_REAL PDstandard4th21gt11 = INITVALUE;
    CCTK_REAL PDstandard4th23gt11 = INITVALUE;
    CCTK_REAL PDstandard4th31gt11 = INITVALUE;
    CCTK_REAL PDstandard4th32gt11 = INITVALUE;
    CCTK_REAL PDstandard4th1gt12 = INITVALUE;
    CCTK_REAL PDstandard4th2gt12 = INITVALUE;
    CCTK_REAL PDstandard4th3gt12 = INITVALUE;
    CCTK_REAL PDstandard4th11gt12 = INITVALUE;
    CCTK_REAL PDstandard4th22gt12 = INITVALUE;
    CCTK_REAL PDstandard4th33gt12 = INITVALUE;
    CCTK_REAL PDstandard4th12gt12 = INITVALUE;
    CCTK_REAL PDstandard4th13gt12 = INITVALUE;
    CCTK_REAL PDstandard4th21gt12 = INITVALUE;
    CCTK_REAL PDstandard4th23gt12 = INITVALUE;
    CCTK_REAL PDstandard4th31gt12 = INITVALUE;
    CCTK_REAL PDstandard4th32gt12 = INITVALUE;
    CCTK_REAL PDstandard4th1gt13 = INITVALUE;
    CCTK_REAL PDstandard4th2gt13 = INITVALUE;
    CCTK_REAL PDstandard4th3gt13 = INITVALUE;
    CCTK_REAL PDstandard4th11gt13 = INITVALUE;
    CCTK_REAL PDstandard4th22gt13 = INITVALUE;
    CCTK_REAL PDstandard4th33gt13 = INITVALUE;
    CCTK_REAL PDstandard4th12gt13 = INITVALUE;
    CCTK_REAL PDstandard4th13gt13 = INITVALUE;
    CCTK_REAL PDstandard4th21gt13 = INITVALUE;
    CCTK_REAL PDstandard4th23gt13 = INITVALUE;
    CCTK_REAL PDstandard4th31gt13 = INITVALUE;
    CCTK_REAL PDstandard4th32gt13 = INITVALUE;
    CCTK_REAL PDstandard4th1gt22 = INITVALUE;
    CCTK_REAL PDstandard4th2gt22 = INITVALUE;
    CCTK_REAL PDstandard4th3gt22 = INITVALUE;
    CCTK_REAL PDstandard4th11gt22 = INITVALUE;
    CCTK_REAL PDstandard4th22gt22 = INITVALUE;
    CCTK_REAL PDstandard4th33gt22 = INITVALUE;
    CCTK_REAL PDstandard4th12gt22 = INITVALUE;
    CCTK_REAL PDstandard4th13gt22 = INITVALUE;
    CCTK_REAL PDstandard4th21gt22 = INITVALUE;
    CCTK_REAL PDstandard4th23gt22 = INITVALUE;
    CCTK_REAL PDstandard4th31gt22 = INITVALUE;
    CCTK_REAL PDstandard4th32gt22 = INITVALUE;
    CCTK_REAL PDstandard4th1gt23 = INITVALUE;
    CCTK_REAL PDstandard4th2gt23 = INITVALUE;
    CCTK_REAL PDstandard4th3gt23 = INITVALUE;
    CCTK_REAL PDstandard4th11gt23 = INITVALUE;
    CCTK_REAL PDstandard4th22gt23 = INITVALUE;
    CCTK_REAL PDstandard4th33gt23 = INITVALUE;
    CCTK_REAL PDstandard4th12gt23 = INITVALUE;
    CCTK_REAL PDstandard4th13gt23 = INITVALUE;
    CCTK_REAL PDstandard4th21gt23 = INITVALUE;
    CCTK_REAL PDstandard4th23gt23 = INITVALUE;
    CCTK_REAL PDstandard4th31gt23 = INITVALUE;
    CCTK_REAL PDstandard4th32gt23 = INITVALUE;
    CCTK_REAL PDstandard4th1gt33 = INITVALUE;
    CCTK_REAL PDstandard4th2gt33 = INITVALUE;
    CCTK_REAL PDstandard4th3gt33 = INITVALUE;
    CCTK_REAL PDstandard4th11gt33 = INITVALUE;
    CCTK_REAL PDstandard4th22gt33 = INITVALUE;
    CCTK_REAL PDstandard4th33gt33 = INITVALUE;
    CCTK_REAL PDstandard4th12gt33 = INITVALUE;
    CCTK_REAL PDstandard4th13gt33 = INITVALUE;
    CCTK_REAL PDstandard4th21gt33 = INITVALUE;
    CCTK_REAL PDstandard4th23gt33 = INITVALUE;
    CCTK_REAL PDstandard4th31gt33 = INITVALUE;
    CCTK_REAL PDstandard4th32gt33 = INITVALUE;
    CCTK_REAL PDstandard4th1phi = INITVALUE;
    CCTK_REAL PDstandard4th2phi = INITVALUE;
    CCTK_REAL PDstandard4th3phi = INITVALUE;
    CCTK_REAL PDstandard4th11phi = INITVALUE;
    CCTK_REAL PDstandard4th22phi = INITVALUE;
    CCTK_REAL PDstandard4th33phi = INITVALUE;
    CCTK_REAL PDstandard4th12phi = INITVALUE;
    CCTK_REAL PDstandard4th13phi = INITVALUE;
    CCTK_REAL PDstandard4th21phi = INITVALUE;
    CCTK_REAL PDstandard4th23phi = INITVALUE;
    CCTK_REAL PDstandard4th31phi = INITVALUE;
    CCTK_REAL PDstandard4th32phi = INITVALUE;
    CCTK_REAL PDstandard4th1trK = INITVALUE;
    CCTK_REAL PDstandard4th2trK = INITVALUE;
    CCTK_REAL PDstandard4th3trK = INITVALUE;
    CCTK_REAL PDstandard4th1Xt1 = INITVALUE;
    CCTK_REAL PDstandard4th2Xt1 = INITVALUE;
    CCTK_REAL PDstandard4th3Xt1 = INITVALUE;
    CCTK_REAL PDstandard4th1Xt2 = INITVALUE;
    CCTK_REAL PDstandard4th2Xt2 = INITVALUE;
    CCTK_REAL PDstandard4th3Xt2 = INITVALUE;
    CCTK_REAL PDstandard4th1Xt3 = INITVALUE;
    CCTK_REAL PDstandard4th2Xt3 = INITVALUE;
    CCTK_REAL PDstandard4th3Xt3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    alphaL = alpha[index];
    At11L = At11[index];
    At12L = At12[index];
    At13L = At13[index];
    At22L = At22[index];
    At23L = At23[index];
    At33L = At33[index];
    beta1L = beta1[index];
    beta2L = beta2[index];
    beta3L = beta3[index];
    dtalphaL = dtalpha[index];
    dtbeta1L = dtbeta1[index];
    dtbeta2L = dtbeta2[index];
    dtbeta3L = dtbeta3[index];
    gt11L = gt11[index];
    gt12L = gt12[index];
    gt13L = gt13[index];
    gt22L = gt22[index];
    gt23L = gt23[index];
    gt33L = gt33[index];
    phiL = phi[index];
    trKL = trK[index];
    trKrhsL = trKrhs[index];
    Xt1L = Xt1[index];
    Xt1rhsL = Xt1rhs[index];
    Xt2L = Xt2[index];
    Xt2rhsL = Xt2rhs[index];
    Xt3L = Xt3[index];
    Xt3rhsL = Xt3rhs[index];
    
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
    PDstandard4th1At11 = PDstandard4th1(At11, i, j, k);
    PDstandard4th2At11 = PDstandard4th2(At11, i, j, k);
    PDstandard4th3At11 = PDstandard4th3(At11, i, j, k);
    PDstandard4th1At12 = PDstandard4th1(At12, i, j, k);
    PDstandard4th2At12 = PDstandard4th2(At12, i, j, k);
    PDstandard4th3At12 = PDstandard4th3(At12, i, j, k);
    PDstandard4th1At13 = PDstandard4th1(At13, i, j, k);
    PDstandard4th2At13 = PDstandard4th2(At13, i, j, k);
    PDstandard4th3At13 = PDstandard4th3(At13, i, j, k);
    PDstandard4th1At22 = PDstandard4th1(At22, i, j, k);
    PDstandard4th2At22 = PDstandard4th2(At22, i, j, k);
    PDstandard4th3At22 = PDstandard4th3(At22, i, j, k);
    PDstandard4th1At23 = PDstandard4th1(At23, i, j, k);
    PDstandard4th2At23 = PDstandard4th2(At23, i, j, k);
    PDstandard4th3At23 = PDstandard4th3(At23, i, j, k);
    PDstandard4th1At33 = PDstandard4th1(At33, i, j, k);
    PDstandard4th2At33 = PDstandard4th2(At33, i, j, k);
    PDstandard4th3At33 = PDstandard4th3(At33, i, j, k);
    PDstandard4th1beta1 = PDstandard4th1(beta1, i, j, k);
    PDstandard4th2beta1 = PDstandard4th2(beta1, i, j, k);
    PDstandard4th3beta1 = PDstandard4th3(beta1, i, j, k);
    PDstandard4th11beta1 = PDstandard4th11(beta1, i, j, k);
    PDstandard4th22beta1 = PDstandard4th22(beta1, i, j, k);
    PDstandard4th33beta1 = PDstandard4th33(beta1, i, j, k);
    PDstandard4th12beta1 = PDstandard4th12(beta1, i, j, k);
    PDstandard4th13beta1 = PDstandard4th13(beta1, i, j, k);
    PDstandard4th23beta1 = PDstandard4th23(beta1, i, j, k);
    PDstandard4th1beta2 = PDstandard4th1(beta2, i, j, k);
    PDstandard4th2beta2 = PDstandard4th2(beta2, i, j, k);
    PDstandard4th3beta2 = PDstandard4th3(beta2, i, j, k);
    PDstandard4th11beta2 = PDstandard4th11(beta2, i, j, k);
    PDstandard4th22beta2 = PDstandard4th22(beta2, i, j, k);
    PDstandard4th33beta2 = PDstandard4th33(beta2, i, j, k);
    PDstandard4th12beta2 = PDstandard4th12(beta2, i, j, k);
    PDstandard4th13beta2 = PDstandard4th13(beta2, i, j, k);
    PDstandard4th23beta2 = PDstandard4th23(beta2, i, j, k);
    PDstandard4th1beta3 = PDstandard4th1(beta3, i, j, k);
    PDstandard4th2beta3 = PDstandard4th2(beta3, i, j, k);
    PDstandard4th3beta3 = PDstandard4th3(beta3, i, j, k);
    PDstandard4th11beta3 = PDstandard4th11(beta3, i, j, k);
    PDstandard4th22beta3 = PDstandard4th22(beta3, i, j, k);
    PDstandard4th33beta3 = PDstandard4th33(beta3, i, j, k);
    PDstandard4th12beta3 = PDstandard4th12(beta3, i, j, k);
    PDstandard4th13beta3 = PDstandard4th13(beta3, i, j, k);
    PDstandard4th23beta3 = PDstandard4th23(beta3, i, j, k);
    PDstandard4th1gt11 = PDstandard4th1(gt11, i, j, k);
    PDstandard4th2gt11 = PDstandard4th2(gt11, i, j, k);
    PDstandard4th3gt11 = PDstandard4th3(gt11, i, j, k);
    PDstandard4th11gt11 = PDstandard4th11(gt11, i, j, k);
    PDstandard4th22gt11 = PDstandard4th22(gt11, i, j, k);
    PDstandard4th33gt11 = PDstandard4th33(gt11, i, j, k);
    PDstandard4th12gt11 = PDstandard4th12(gt11, i, j, k);
    PDstandard4th13gt11 = PDstandard4th13(gt11, i, j, k);
    PDstandard4th23gt11 = PDstandard4th23(gt11, i, j, k);
    PDstandard4th1gt12 = PDstandard4th1(gt12, i, j, k);
    PDstandard4th2gt12 = PDstandard4th2(gt12, i, j, k);
    PDstandard4th3gt12 = PDstandard4th3(gt12, i, j, k);
    PDstandard4th11gt12 = PDstandard4th11(gt12, i, j, k);
    PDstandard4th22gt12 = PDstandard4th22(gt12, i, j, k);
    PDstandard4th33gt12 = PDstandard4th33(gt12, i, j, k);
    PDstandard4th12gt12 = PDstandard4th12(gt12, i, j, k);
    PDstandard4th13gt12 = PDstandard4th13(gt12, i, j, k);
    PDstandard4th23gt12 = PDstandard4th23(gt12, i, j, k);
    PDstandard4th1gt13 = PDstandard4th1(gt13, i, j, k);
    PDstandard4th2gt13 = PDstandard4th2(gt13, i, j, k);
    PDstandard4th3gt13 = PDstandard4th3(gt13, i, j, k);
    PDstandard4th11gt13 = PDstandard4th11(gt13, i, j, k);
    PDstandard4th22gt13 = PDstandard4th22(gt13, i, j, k);
    PDstandard4th33gt13 = PDstandard4th33(gt13, i, j, k);
    PDstandard4th12gt13 = PDstandard4th12(gt13, i, j, k);
    PDstandard4th13gt13 = PDstandard4th13(gt13, i, j, k);
    PDstandard4th23gt13 = PDstandard4th23(gt13, i, j, k);
    PDstandard4th1gt22 = PDstandard4th1(gt22, i, j, k);
    PDstandard4th2gt22 = PDstandard4th2(gt22, i, j, k);
    PDstandard4th3gt22 = PDstandard4th3(gt22, i, j, k);
    PDstandard4th11gt22 = PDstandard4th11(gt22, i, j, k);
    PDstandard4th22gt22 = PDstandard4th22(gt22, i, j, k);
    PDstandard4th33gt22 = PDstandard4th33(gt22, i, j, k);
    PDstandard4th12gt22 = PDstandard4th12(gt22, i, j, k);
    PDstandard4th13gt22 = PDstandard4th13(gt22, i, j, k);
    PDstandard4th23gt22 = PDstandard4th23(gt22, i, j, k);
    PDstandard4th1gt23 = PDstandard4th1(gt23, i, j, k);
    PDstandard4th2gt23 = PDstandard4th2(gt23, i, j, k);
    PDstandard4th3gt23 = PDstandard4th3(gt23, i, j, k);
    PDstandard4th11gt23 = PDstandard4th11(gt23, i, j, k);
    PDstandard4th22gt23 = PDstandard4th22(gt23, i, j, k);
    PDstandard4th33gt23 = PDstandard4th33(gt23, i, j, k);
    PDstandard4th12gt23 = PDstandard4th12(gt23, i, j, k);
    PDstandard4th13gt23 = PDstandard4th13(gt23, i, j, k);
    PDstandard4th23gt23 = PDstandard4th23(gt23, i, j, k);
    PDstandard4th1gt33 = PDstandard4th1(gt33, i, j, k);
    PDstandard4th2gt33 = PDstandard4th2(gt33, i, j, k);
    PDstandard4th3gt33 = PDstandard4th3(gt33, i, j, k);
    PDstandard4th11gt33 = PDstandard4th11(gt33, i, j, k);
    PDstandard4th22gt33 = PDstandard4th22(gt33, i, j, k);
    PDstandard4th33gt33 = PDstandard4th33(gt33, i, j, k);
    PDstandard4th12gt33 = PDstandard4th12(gt33, i, j, k);
    PDstandard4th13gt33 = PDstandard4th13(gt33, i, j, k);
    PDstandard4th23gt33 = PDstandard4th23(gt33, i, j, k);
    PDstandard4th1phi = PDstandard4th1(phi, i, j, k);
    PDstandard4th2phi = PDstandard4th2(phi, i, j, k);
    PDstandard4th3phi = PDstandard4th3(phi, i, j, k);
    PDstandard4th11phi = PDstandard4th11(phi, i, j, k);
    PDstandard4th22phi = PDstandard4th22(phi, i, j, k);
    PDstandard4th33phi = PDstandard4th33(phi, i, j, k);
    PDstandard4th12phi = PDstandard4th12(phi, i, j, k);
    PDstandard4th13phi = PDstandard4th13(phi, i, j, k);
    PDstandard4th23phi = PDstandard4th23(phi, i, j, k);
    PDstandard4th1trK = PDstandard4th1(trK, i, j, k);
    PDstandard4th2trK = PDstandard4th2(trK, i, j, k);
    PDstandard4th3trK = PDstandard4th3(trK, i, j, k);
    PDstandard4th1Xt1 = PDstandard4th1(Xt1, i, j, k);
    PDstandard4th2Xt1 = PDstandard4th2(Xt1, i, j, k);
    PDstandard4th3Xt1 = PDstandard4th3(Xt1, i, j, k);
    PDstandard4th1Xt2 = PDstandard4th1(Xt2, i, j, k);
    PDstandard4th2Xt2 = PDstandard4th2(Xt2, i, j, k);
    PDstandard4th3Xt2 = PDstandard4th3(Xt2, i, j, k);
    PDstandard4th1Xt3 = PDstandard4th1(Xt3, i, j, k);
    PDstandard4th2Xt3 = PDstandard4th2(Xt3, i, j, k);
    PDstandard4th3Xt3 = PDstandard4th3(Xt3, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    detgt  =  1;
    
    ddetgt1  =  0;
    
    ddetgt2  =  0;
    
    ddetgt3  =  0;
    
    gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    gtu21  =  (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    gtu31  =  (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    gtu32  =  (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    dgtu111  =  -2*(gtu11*gtu21*PDstandard4th1gt12 + gtu11*gtu31*PDstandard4th1gt13 + gtu21*gtu31*PDstandard4th1gt23) - 
        PDstandard4th1gt11*SQR(gtu11) - PDstandard4th1gt22*SQR(gtu21) - PDstandard4th1gt33*SQR(gtu31);
    
    dgtu211  =  -(gtu11*gtu21*PDstandard4th1gt11) - gtu11*gtu22*PDstandard4th1gt12 - gtu21*gtu31*PDstandard4th1gt13 - 
        gtu11*gtu32*PDstandard4th1gt13 - gtu21*gtu22*PDstandard4th1gt22 - gtu22*gtu31*PDstandard4th1gt23 - 
        gtu21*gtu32*PDstandard4th1gt23 - gtu31*gtu32*PDstandard4th1gt33 - PDstandard4th1gt12*SQR(gtu21);
    
    dgtu311  =  -(gtu11*gtu31*PDstandard4th1gt11) - gtu21*gtu31*PDstandard4th1gt12 - gtu11*gtu32*PDstandard4th1gt12 - 
        gtu11*gtu33*PDstandard4th1gt13 - gtu21*gtu32*PDstandard4th1gt22 - gtu31*gtu32*PDstandard4th1gt23 - 
        gtu21*gtu33*PDstandard4th1gt23 - gtu31*gtu33*PDstandard4th1gt33 - PDstandard4th1gt13*SQR(gtu31);
    
    dgtu221  =  -2*(gtu21*gtu22*PDstandard4th1gt12 + gtu21*gtu32*PDstandard4th1gt13 + gtu22*gtu32*PDstandard4th1gt23) - 
        PDstandard4th1gt11*SQR(gtu21) - PDstandard4th1gt22*SQR(gtu22) - PDstandard4th1gt33*SQR(gtu32);
    
    dgtu321  =  -(gtu21*gtu31*PDstandard4th1gt11) - gtu22*gtu31*PDstandard4th1gt12 - gtu21*gtu32*PDstandard4th1gt12 - 
        gtu31*gtu32*PDstandard4th1gt13 - gtu21*gtu33*PDstandard4th1gt13 - gtu22*gtu32*PDstandard4th1gt22 - 
        gtu22*gtu33*PDstandard4th1gt23 - gtu32*gtu33*PDstandard4th1gt33 - PDstandard4th1gt23*SQR(gtu32);
    
    dgtu331  =  -2*(gtu31*gtu32*PDstandard4th1gt12 + gtu31*gtu33*PDstandard4th1gt13 + gtu32*gtu33*PDstandard4th1gt23) - 
        PDstandard4th1gt11*SQR(gtu31) - PDstandard4th1gt22*SQR(gtu32) - PDstandard4th1gt33*SQR(gtu33);
    
    dgtu112  =  -2*(gtu11*gtu21*PDstandard4th2gt12 + gtu11*gtu31*PDstandard4th2gt13 + gtu21*gtu31*PDstandard4th2gt23) - 
        PDstandard4th2gt11*SQR(gtu11) - PDstandard4th2gt22*SQR(gtu21) - PDstandard4th2gt33*SQR(gtu31);
    
    dgtu212  =  -(gtu11*gtu21*PDstandard4th2gt11) - gtu11*gtu22*PDstandard4th2gt12 - gtu21*gtu31*PDstandard4th2gt13 - 
        gtu11*gtu32*PDstandard4th2gt13 - gtu21*gtu22*PDstandard4th2gt22 - gtu22*gtu31*PDstandard4th2gt23 - 
        gtu21*gtu32*PDstandard4th2gt23 - gtu31*gtu32*PDstandard4th2gt33 - PDstandard4th2gt12*SQR(gtu21);
    
    dgtu312  =  -(gtu11*gtu31*PDstandard4th2gt11) - gtu21*gtu31*PDstandard4th2gt12 - gtu11*gtu32*PDstandard4th2gt12 - 
        gtu11*gtu33*PDstandard4th2gt13 - gtu21*gtu32*PDstandard4th2gt22 - gtu31*gtu32*PDstandard4th2gt23 - 
        gtu21*gtu33*PDstandard4th2gt23 - gtu31*gtu33*PDstandard4th2gt33 - PDstandard4th2gt13*SQR(gtu31);
    
    dgtu222  =  -2*(gtu21*gtu22*PDstandard4th2gt12 + gtu21*gtu32*PDstandard4th2gt13 + gtu22*gtu32*PDstandard4th2gt23) - 
        PDstandard4th2gt11*SQR(gtu21) - PDstandard4th2gt22*SQR(gtu22) - PDstandard4th2gt33*SQR(gtu32);
    
    dgtu322  =  -(gtu21*gtu31*PDstandard4th2gt11) - gtu22*gtu31*PDstandard4th2gt12 - gtu21*gtu32*PDstandard4th2gt12 - 
        gtu31*gtu32*PDstandard4th2gt13 - gtu21*gtu33*PDstandard4th2gt13 - gtu22*gtu32*PDstandard4th2gt22 - 
        gtu22*gtu33*PDstandard4th2gt23 - gtu32*gtu33*PDstandard4th2gt33 - PDstandard4th2gt23*SQR(gtu32);
    
    dgtu332  =  -2*(gtu31*gtu32*PDstandard4th2gt12 + gtu31*gtu33*PDstandard4th2gt13 + gtu32*gtu33*PDstandard4th2gt23) - 
        PDstandard4th2gt11*SQR(gtu31) - PDstandard4th2gt22*SQR(gtu32) - PDstandard4th2gt33*SQR(gtu33);
    
    dgtu113  =  -2*(gtu11*gtu21*PDstandard4th3gt12 + gtu11*gtu31*PDstandard4th3gt13 + gtu21*gtu31*PDstandard4th3gt23) - 
        PDstandard4th3gt11*SQR(gtu11) - PDstandard4th3gt22*SQR(gtu21) - PDstandard4th3gt33*SQR(gtu31);
    
    dgtu213  =  -(gtu11*gtu21*PDstandard4th3gt11) - gtu11*gtu22*PDstandard4th3gt12 - gtu21*gtu31*PDstandard4th3gt13 - 
        gtu11*gtu32*PDstandard4th3gt13 - gtu21*gtu22*PDstandard4th3gt22 - gtu22*gtu31*PDstandard4th3gt23 - 
        gtu21*gtu32*PDstandard4th3gt23 - gtu31*gtu32*PDstandard4th3gt33 - PDstandard4th3gt12*SQR(gtu21);
    
    dgtu313  =  -(gtu11*gtu31*PDstandard4th3gt11) - gtu21*gtu31*PDstandard4th3gt12 - gtu11*gtu32*PDstandard4th3gt12 - 
        gtu11*gtu33*PDstandard4th3gt13 - gtu21*gtu32*PDstandard4th3gt22 - gtu31*gtu32*PDstandard4th3gt23 - 
        gtu21*gtu33*PDstandard4th3gt23 - gtu31*gtu33*PDstandard4th3gt33 - PDstandard4th3gt13*SQR(gtu31);
    
    dgtu223  =  -2*(gtu21*gtu22*PDstandard4th3gt12 + gtu21*gtu32*PDstandard4th3gt13 + gtu22*gtu32*PDstandard4th3gt23) - 
        PDstandard4th3gt11*SQR(gtu21) - PDstandard4th3gt22*SQR(gtu22) - PDstandard4th3gt33*SQR(gtu32);
    
    dgtu323  =  -(gtu21*gtu31*PDstandard4th3gt11) - gtu22*gtu31*PDstandard4th3gt12 - gtu21*gtu32*PDstandard4th3gt12 - 
        gtu31*gtu32*PDstandard4th3gt13 - gtu21*gtu33*PDstandard4th3gt13 - gtu22*gtu32*PDstandard4th3gt22 - 
        gtu22*gtu33*PDstandard4th3gt23 - gtu32*gtu33*PDstandard4th3gt33 - PDstandard4th3gt23*SQR(gtu32);
    
    dgtu333  =  -2*(gtu31*gtu32*PDstandard4th3gt12 + gtu31*gtu33*PDstandard4th3gt13 + gtu32*gtu33*PDstandard4th3gt23) - 
        PDstandard4th3gt11*SQR(gtu31) - PDstandard4th3gt22*SQR(gtu32) - PDstandard4th3gt33*SQR(gtu33);
    
    ddgtu1111  =  -2*(gtu11*gtu21*PDstandard4th11gt12 + gtu11*gtu31*PDstandard4th11gt13 + gtu21*gtu31*PDstandard4th11gt23 + 
           dgtu111*gtu11*PDstandard4th1gt11 + dgtu211*gtu11*PDstandard4th1gt12 + dgtu111*gtu21*PDstandard4th1gt12 + 
           dgtu311*gtu11*PDstandard4th1gt13 + dgtu111*gtu31*PDstandard4th1gt13 + dgtu211*gtu21*PDstandard4th1gt22 + 
           dgtu311*gtu21*PDstandard4th1gt23 + dgtu211*gtu31*PDstandard4th1gt23 + dgtu311*gtu31*PDstandard4th1gt33) - 
        PDstandard4th11gt11*SQR(gtu11) - PDstandard4th11gt22*SQR(gtu21) - PDstandard4th11gt33*SQR(gtu31);
    
    ddgtu2111  =  -(gtu11*gtu21*PDstandard4th11gt11) - gtu11*gtu22*PDstandard4th11gt12 - gtu21*gtu31*PDstandard4th11gt13 - 
        gtu11*gtu32*PDstandard4th11gt13 - gtu21*gtu22*PDstandard4th11gt22 - gtu22*gtu31*PDstandard4th11gt23 - 
        gtu21*gtu32*PDstandard4th11gt23 - gtu31*gtu32*PDstandard4th11gt33 - dgtu211*gtu11*PDstandard4th1gt11 - 
        dgtu111*gtu21*PDstandard4th1gt11 - dgtu221*gtu11*PDstandard4th1gt12 - 2*dgtu211*gtu21*PDstandard4th1gt12 - 
        dgtu111*gtu22*PDstandard4th1gt12 - dgtu321*gtu11*PDstandard4th1gt13 - dgtu311*gtu21*PDstandard4th1gt13 - 
        dgtu211*gtu31*PDstandard4th1gt13 - dgtu111*gtu32*PDstandard4th1gt13 - dgtu221*gtu21*PDstandard4th1gt22 - 
        dgtu211*gtu22*PDstandard4th1gt22 - dgtu321*gtu21*PDstandard4th1gt23 - dgtu311*gtu22*PDstandard4th1gt23 - 
        dgtu221*gtu31*PDstandard4th1gt23 - dgtu211*gtu32*PDstandard4th1gt23 - dgtu321*gtu31*PDstandard4th1gt33 - 
        dgtu311*gtu32*PDstandard4th1gt33 - PDstandard4th11gt12*SQR(gtu21);
    
    ddgtu3111  =  -(gtu11*gtu31*PDstandard4th11gt11) - gtu21*gtu31*PDstandard4th11gt12 - gtu11*gtu32*PDstandard4th11gt12 - 
        gtu11*gtu33*PDstandard4th11gt13 - gtu21*gtu32*PDstandard4th11gt22 - gtu31*gtu32*PDstandard4th11gt23 - 
        gtu21*gtu33*PDstandard4th11gt23 - gtu31*gtu33*PDstandard4th11gt33 - dgtu311*gtu11*PDstandard4th1gt11 - 
        dgtu111*gtu31*PDstandard4th1gt11 - dgtu321*gtu11*PDstandard4th1gt12 - dgtu311*gtu21*PDstandard4th1gt12 - 
        dgtu211*gtu31*PDstandard4th1gt12 - dgtu111*gtu32*PDstandard4th1gt12 - dgtu331*gtu11*PDstandard4th1gt13 - 
        2*dgtu311*gtu31*PDstandard4th1gt13 - dgtu111*gtu33*PDstandard4th1gt13 - dgtu321*gtu21*PDstandard4th1gt22 - 
        dgtu211*gtu32*PDstandard4th1gt22 - dgtu331*gtu21*PDstandard4th1gt23 - dgtu321*gtu31*PDstandard4th1gt23 - 
        dgtu311*gtu32*PDstandard4th1gt23 - dgtu211*gtu33*PDstandard4th1gt23 - dgtu331*gtu31*PDstandard4th1gt33 - 
        dgtu311*gtu33*PDstandard4th1gt33 - PDstandard4th11gt13*SQR(gtu31);
    
    ddgtu1112  =  -2*(gtu11*gtu21*PDstandard4th12gt12 + gtu11*gtu31*PDstandard4th12gt13 + gtu21*gtu31*PDstandard4th12gt23 + 
           dgtu112*gtu11*PDstandard4th1gt11 + dgtu212*gtu11*PDstandard4th1gt12 + dgtu112*gtu21*PDstandard4th1gt12 + 
           dgtu312*gtu11*PDstandard4th1gt13 + dgtu112*gtu31*PDstandard4th1gt13 + dgtu212*gtu21*PDstandard4th1gt22 + 
           dgtu312*gtu21*PDstandard4th1gt23 + dgtu212*gtu31*PDstandard4th1gt23 + dgtu312*gtu31*PDstandard4th1gt33) - 
        PDstandard4th12gt11*SQR(gtu11) - PDstandard4th12gt22*SQR(gtu21) - PDstandard4th12gt33*SQR(gtu31);
    
    ddgtu2112  =  -(gtu11*gtu21*PDstandard4th12gt11) - gtu11*gtu22*PDstandard4th12gt12 - gtu21*gtu31*PDstandard4th12gt13 - 
        gtu11*gtu32*PDstandard4th12gt13 - gtu21*gtu22*PDstandard4th12gt22 - gtu22*gtu31*PDstandard4th12gt23 - 
        gtu21*gtu32*PDstandard4th12gt23 - gtu31*gtu32*PDstandard4th12gt33 - dgtu212*gtu11*PDstandard4th1gt11 - 
        dgtu112*gtu21*PDstandard4th1gt11 - dgtu222*gtu11*PDstandard4th1gt12 - 2*dgtu212*gtu21*PDstandard4th1gt12 - 
        dgtu112*gtu22*PDstandard4th1gt12 - dgtu322*gtu11*PDstandard4th1gt13 - dgtu312*gtu21*PDstandard4th1gt13 - 
        dgtu212*gtu31*PDstandard4th1gt13 - dgtu112*gtu32*PDstandard4th1gt13 - dgtu222*gtu21*PDstandard4th1gt22 - 
        dgtu212*gtu22*PDstandard4th1gt22 - dgtu322*gtu21*PDstandard4th1gt23 - dgtu312*gtu22*PDstandard4th1gt23 - 
        dgtu222*gtu31*PDstandard4th1gt23 - dgtu212*gtu32*PDstandard4th1gt23 - dgtu322*gtu31*PDstandard4th1gt33 - 
        dgtu312*gtu32*PDstandard4th1gt33 - PDstandard4th12gt12*SQR(gtu21);
    
    ddgtu3112  =  -(gtu11*gtu31*PDstandard4th12gt11) - gtu21*gtu31*PDstandard4th12gt12 - gtu11*gtu32*PDstandard4th12gt12 - 
        gtu11*gtu33*PDstandard4th12gt13 - gtu21*gtu32*PDstandard4th12gt22 - gtu31*gtu32*PDstandard4th12gt23 - 
        gtu21*gtu33*PDstandard4th12gt23 - gtu31*gtu33*PDstandard4th12gt33 - dgtu312*gtu11*PDstandard4th1gt11 - 
        dgtu112*gtu31*PDstandard4th1gt11 - dgtu322*gtu11*PDstandard4th1gt12 - dgtu312*gtu21*PDstandard4th1gt12 - 
        dgtu212*gtu31*PDstandard4th1gt12 - dgtu112*gtu32*PDstandard4th1gt12 - dgtu332*gtu11*PDstandard4th1gt13 - 
        2*dgtu312*gtu31*PDstandard4th1gt13 - dgtu112*gtu33*PDstandard4th1gt13 - dgtu322*gtu21*PDstandard4th1gt22 - 
        dgtu212*gtu32*PDstandard4th1gt22 - dgtu332*gtu21*PDstandard4th1gt23 - dgtu322*gtu31*PDstandard4th1gt23 - 
        dgtu312*gtu32*PDstandard4th1gt23 - dgtu212*gtu33*PDstandard4th1gt23 - dgtu332*gtu31*PDstandard4th1gt33 - 
        dgtu312*gtu33*PDstandard4th1gt33 - PDstandard4th12gt13*SQR(gtu31);
    
    ddgtu2212  =  -2*(gtu21*gtu22*PDstandard4th12gt12 + gtu21*gtu32*PDstandard4th12gt13 + gtu22*gtu32*PDstandard4th12gt23 + 
           dgtu212*gtu21*PDstandard4th1gt11 + dgtu222*gtu21*PDstandard4th1gt12 + dgtu212*gtu22*PDstandard4th1gt12 + 
           dgtu322*gtu21*PDstandard4th1gt13 + dgtu212*gtu32*PDstandard4th1gt13 + dgtu222*gtu22*PDstandard4th1gt22 + 
           dgtu322*gtu22*PDstandard4th1gt23 + dgtu222*gtu32*PDstandard4th1gt23 + dgtu322*gtu32*PDstandard4th1gt33) - 
        PDstandard4th12gt11*SQR(gtu21) - PDstandard4th12gt22*SQR(gtu22) - PDstandard4th12gt33*SQR(gtu32);
    
    ddgtu3212  =  -(gtu21*gtu31*PDstandard4th12gt11) - gtu22*gtu31*PDstandard4th12gt12 - gtu21*gtu32*PDstandard4th12gt12 - 
        gtu31*gtu32*PDstandard4th12gt13 - gtu21*gtu33*PDstandard4th12gt13 - gtu22*gtu32*PDstandard4th12gt22 - 
        gtu22*gtu33*PDstandard4th12gt23 - gtu32*gtu33*PDstandard4th12gt33 - dgtu312*gtu21*PDstandard4th1gt11 - 
        dgtu212*gtu31*PDstandard4th1gt11 - dgtu322*gtu21*PDstandard4th1gt12 - dgtu312*gtu22*PDstandard4th1gt12 - 
        dgtu222*gtu31*PDstandard4th1gt12 - dgtu212*gtu32*PDstandard4th1gt12 - dgtu332*gtu21*PDstandard4th1gt13 - 
        dgtu322*gtu31*PDstandard4th1gt13 - dgtu312*gtu32*PDstandard4th1gt13 - dgtu212*gtu33*PDstandard4th1gt13 - 
        dgtu322*gtu22*PDstandard4th1gt22 - dgtu222*gtu32*PDstandard4th1gt22 - dgtu332*gtu22*PDstandard4th1gt23 - 
        2*dgtu322*gtu32*PDstandard4th1gt23 - dgtu222*gtu33*PDstandard4th1gt23 - dgtu332*gtu32*PDstandard4th1gt33 - 
        dgtu322*gtu33*PDstandard4th1gt33 - PDstandard4th12gt23*SQR(gtu32);
    
    ddgtu1113  =  -2*(gtu11*gtu21*PDstandard4th13gt12 + gtu11*gtu31*PDstandard4th13gt13 + gtu21*gtu31*PDstandard4th13gt23 + 
           dgtu113*gtu11*PDstandard4th1gt11 + dgtu213*gtu11*PDstandard4th1gt12 + dgtu113*gtu21*PDstandard4th1gt12 + 
           dgtu313*gtu11*PDstandard4th1gt13 + dgtu113*gtu31*PDstandard4th1gt13 + dgtu213*gtu21*PDstandard4th1gt22 + 
           dgtu313*gtu21*PDstandard4th1gt23 + dgtu213*gtu31*PDstandard4th1gt23 + dgtu313*gtu31*PDstandard4th1gt33) - 
        PDstandard4th13gt11*SQR(gtu11) - PDstandard4th13gt22*SQR(gtu21) - PDstandard4th13gt33*SQR(gtu31);
    
    ddgtu2113  =  -(gtu11*gtu21*PDstandard4th13gt11) - gtu11*gtu22*PDstandard4th13gt12 - gtu21*gtu31*PDstandard4th13gt13 - 
        gtu11*gtu32*PDstandard4th13gt13 - gtu21*gtu22*PDstandard4th13gt22 - gtu22*gtu31*PDstandard4th13gt23 - 
        gtu21*gtu32*PDstandard4th13gt23 - gtu31*gtu32*PDstandard4th13gt33 - dgtu213*gtu11*PDstandard4th1gt11 - 
        dgtu113*gtu21*PDstandard4th1gt11 - dgtu223*gtu11*PDstandard4th1gt12 - 2*dgtu213*gtu21*PDstandard4th1gt12 - 
        dgtu113*gtu22*PDstandard4th1gt12 - dgtu323*gtu11*PDstandard4th1gt13 - dgtu313*gtu21*PDstandard4th1gt13 - 
        dgtu213*gtu31*PDstandard4th1gt13 - dgtu113*gtu32*PDstandard4th1gt13 - dgtu223*gtu21*PDstandard4th1gt22 - 
        dgtu213*gtu22*PDstandard4th1gt22 - dgtu323*gtu21*PDstandard4th1gt23 - dgtu313*gtu22*PDstandard4th1gt23 - 
        dgtu223*gtu31*PDstandard4th1gt23 - dgtu213*gtu32*PDstandard4th1gt23 - dgtu323*gtu31*PDstandard4th1gt33 - 
        dgtu313*gtu32*PDstandard4th1gt33 - PDstandard4th13gt12*SQR(gtu21);
    
    ddgtu3113  =  -(gtu11*gtu31*PDstandard4th13gt11) - gtu21*gtu31*PDstandard4th13gt12 - gtu11*gtu32*PDstandard4th13gt12 - 
        gtu11*gtu33*PDstandard4th13gt13 - gtu21*gtu32*PDstandard4th13gt22 - gtu31*gtu32*PDstandard4th13gt23 - 
        gtu21*gtu33*PDstandard4th13gt23 - gtu31*gtu33*PDstandard4th13gt33 - dgtu313*gtu11*PDstandard4th1gt11 - 
        dgtu113*gtu31*PDstandard4th1gt11 - dgtu323*gtu11*PDstandard4th1gt12 - dgtu313*gtu21*PDstandard4th1gt12 - 
        dgtu213*gtu31*PDstandard4th1gt12 - dgtu113*gtu32*PDstandard4th1gt12 - dgtu333*gtu11*PDstandard4th1gt13 - 
        2*dgtu313*gtu31*PDstandard4th1gt13 - dgtu113*gtu33*PDstandard4th1gt13 - dgtu323*gtu21*PDstandard4th1gt22 - 
        dgtu213*gtu32*PDstandard4th1gt22 - dgtu333*gtu21*PDstandard4th1gt23 - dgtu323*gtu31*PDstandard4th1gt23 - 
        dgtu313*gtu32*PDstandard4th1gt23 - dgtu213*gtu33*PDstandard4th1gt23 - dgtu333*gtu31*PDstandard4th1gt33 - 
        dgtu313*gtu33*PDstandard4th1gt33 - PDstandard4th13gt13*SQR(gtu31);
    
    ddgtu3213  =  -(gtu21*gtu31*PDstandard4th13gt11) - gtu22*gtu31*PDstandard4th13gt12 - gtu21*gtu32*PDstandard4th13gt12 - 
        gtu31*gtu32*PDstandard4th13gt13 - gtu21*gtu33*PDstandard4th13gt13 - gtu22*gtu32*PDstandard4th13gt22 - 
        gtu22*gtu33*PDstandard4th13gt23 - gtu32*gtu33*PDstandard4th13gt33 - dgtu313*gtu21*PDstandard4th1gt11 - 
        dgtu213*gtu31*PDstandard4th1gt11 - dgtu323*gtu21*PDstandard4th1gt12 - dgtu313*gtu22*PDstandard4th1gt12 - 
        dgtu223*gtu31*PDstandard4th1gt12 - dgtu213*gtu32*PDstandard4th1gt12 - dgtu333*gtu21*PDstandard4th1gt13 - 
        dgtu323*gtu31*PDstandard4th1gt13 - dgtu313*gtu32*PDstandard4th1gt13 - dgtu213*gtu33*PDstandard4th1gt13 - 
        dgtu323*gtu22*PDstandard4th1gt22 - dgtu223*gtu32*PDstandard4th1gt22 - dgtu333*gtu22*PDstandard4th1gt23 - 
        2*dgtu323*gtu32*PDstandard4th1gt23 - dgtu223*gtu33*PDstandard4th1gt23 - dgtu333*gtu32*PDstandard4th1gt33 - 
        dgtu323*gtu33*PDstandard4th1gt33 - PDstandard4th13gt23*SQR(gtu32);
    
    ddgtu3313  =  -2*(gtu31*gtu32*PDstandard4th13gt12 + gtu31*gtu33*PDstandard4th13gt13 + gtu32*gtu33*PDstandard4th13gt23 + 
           dgtu313*gtu31*PDstandard4th1gt11 + dgtu323*gtu31*PDstandard4th1gt12 + dgtu313*gtu32*PDstandard4th1gt12 + 
           dgtu333*gtu31*PDstandard4th1gt13 + dgtu313*gtu33*PDstandard4th1gt13 + dgtu323*gtu32*PDstandard4th1gt22 + 
           dgtu333*gtu32*PDstandard4th1gt23 + dgtu323*gtu33*PDstandard4th1gt23 + dgtu333*gtu33*PDstandard4th1gt33) - 
        PDstandard4th13gt11*SQR(gtu31) - PDstandard4th13gt22*SQR(gtu32) - PDstandard4th13gt33*SQR(gtu33);
    
    ddgtu2122  =  -(gtu11*gtu21*PDstandard4th22gt11) - gtu11*gtu22*PDstandard4th22gt12 - gtu21*gtu31*PDstandard4th22gt13 - 
        gtu11*gtu32*PDstandard4th22gt13 - gtu21*gtu22*PDstandard4th22gt22 - gtu22*gtu31*PDstandard4th22gt23 - 
        gtu21*gtu32*PDstandard4th22gt23 - gtu31*gtu32*PDstandard4th22gt33 - dgtu212*gtu11*PDstandard4th2gt11 - 
        dgtu112*gtu21*PDstandard4th2gt11 - dgtu222*gtu11*PDstandard4th2gt12 - 2*dgtu212*gtu21*PDstandard4th2gt12 - 
        dgtu112*gtu22*PDstandard4th2gt12 - dgtu322*gtu11*PDstandard4th2gt13 - dgtu312*gtu21*PDstandard4th2gt13 - 
        dgtu212*gtu31*PDstandard4th2gt13 - dgtu112*gtu32*PDstandard4th2gt13 - dgtu222*gtu21*PDstandard4th2gt22 - 
        dgtu212*gtu22*PDstandard4th2gt22 - dgtu322*gtu21*PDstandard4th2gt23 - dgtu312*gtu22*PDstandard4th2gt23 - 
        dgtu222*gtu31*PDstandard4th2gt23 - dgtu212*gtu32*PDstandard4th2gt23 - dgtu322*gtu31*PDstandard4th2gt33 - 
        dgtu312*gtu32*PDstandard4th2gt33 - PDstandard4th22gt12*SQR(gtu21);
    
    ddgtu2222  =  -2*(gtu21*gtu22*PDstandard4th22gt12 + gtu21*gtu32*PDstandard4th22gt13 + gtu22*gtu32*PDstandard4th22gt23 + 
           dgtu212*gtu21*PDstandard4th2gt11 + dgtu222*gtu21*PDstandard4th2gt12 + dgtu212*gtu22*PDstandard4th2gt12 + 
           dgtu322*gtu21*PDstandard4th2gt13 + dgtu212*gtu32*PDstandard4th2gt13 + dgtu222*gtu22*PDstandard4th2gt22 + 
           dgtu322*gtu22*PDstandard4th2gt23 + dgtu222*gtu32*PDstandard4th2gt23 + dgtu322*gtu32*PDstandard4th2gt33) - 
        PDstandard4th22gt11*SQR(gtu21) - PDstandard4th22gt22*SQR(gtu22) - PDstandard4th22gt33*SQR(gtu32);
    
    ddgtu3222  =  -(gtu21*gtu31*PDstandard4th22gt11) - gtu22*gtu31*PDstandard4th22gt12 - gtu21*gtu32*PDstandard4th22gt12 - 
        gtu31*gtu32*PDstandard4th22gt13 - gtu21*gtu33*PDstandard4th22gt13 - gtu22*gtu32*PDstandard4th22gt22 - 
        gtu22*gtu33*PDstandard4th22gt23 - gtu32*gtu33*PDstandard4th22gt33 - dgtu312*gtu21*PDstandard4th2gt11 - 
        dgtu212*gtu31*PDstandard4th2gt11 - dgtu322*gtu21*PDstandard4th2gt12 - dgtu312*gtu22*PDstandard4th2gt12 - 
        dgtu222*gtu31*PDstandard4th2gt12 - dgtu212*gtu32*PDstandard4th2gt12 - dgtu332*gtu21*PDstandard4th2gt13 - 
        dgtu322*gtu31*PDstandard4th2gt13 - dgtu312*gtu32*PDstandard4th2gt13 - dgtu212*gtu33*PDstandard4th2gt13 - 
        dgtu322*gtu22*PDstandard4th2gt22 - dgtu222*gtu32*PDstandard4th2gt22 - dgtu332*gtu22*PDstandard4th2gt23 - 
        2*dgtu322*gtu32*PDstandard4th2gt23 - dgtu222*gtu33*PDstandard4th2gt23 - dgtu332*gtu32*PDstandard4th2gt33 - 
        dgtu322*gtu33*PDstandard4th2gt33 - PDstandard4th22gt23*SQR(gtu32);
    
    ddgtu2123  =  -(gtu11*gtu21*PDstandard4th23gt11) - gtu11*gtu22*PDstandard4th23gt12 - gtu21*gtu31*PDstandard4th23gt13 - 
        gtu11*gtu32*PDstandard4th23gt13 - gtu21*gtu22*PDstandard4th23gt22 - gtu22*gtu31*PDstandard4th23gt23 - 
        gtu21*gtu32*PDstandard4th23gt23 - gtu31*gtu32*PDstandard4th23gt33 - dgtu213*gtu11*PDstandard4th2gt11 - 
        dgtu113*gtu21*PDstandard4th2gt11 - dgtu223*gtu11*PDstandard4th2gt12 - 2*dgtu213*gtu21*PDstandard4th2gt12 - 
        dgtu113*gtu22*PDstandard4th2gt12 - dgtu323*gtu11*PDstandard4th2gt13 - dgtu313*gtu21*PDstandard4th2gt13 - 
        dgtu213*gtu31*PDstandard4th2gt13 - dgtu113*gtu32*PDstandard4th2gt13 - dgtu223*gtu21*PDstandard4th2gt22 - 
        dgtu213*gtu22*PDstandard4th2gt22 - dgtu323*gtu21*PDstandard4th2gt23 - dgtu313*gtu22*PDstandard4th2gt23 - 
        dgtu223*gtu31*PDstandard4th2gt23 - dgtu213*gtu32*PDstandard4th2gt23 - dgtu323*gtu31*PDstandard4th2gt33 - 
        dgtu313*gtu32*PDstandard4th2gt33 - PDstandard4th23gt12*SQR(gtu21);
    
    ddgtu3123  =  -(gtu11*gtu31*PDstandard4th23gt11) - gtu21*gtu31*PDstandard4th23gt12 - gtu11*gtu32*PDstandard4th23gt12 - 
        gtu11*gtu33*PDstandard4th23gt13 - gtu21*gtu32*PDstandard4th23gt22 - gtu31*gtu32*PDstandard4th23gt23 - 
        gtu21*gtu33*PDstandard4th23gt23 - gtu31*gtu33*PDstandard4th23gt33 - dgtu313*gtu11*PDstandard4th2gt11 - 
        dgtu113*gtu31*PDstandard4th2gt11 - dgtu323*gtu11*PDstandard4th2gt12 - dgtu313*gtu21*PDstandard4th2gt12 - 
        dgtu213*gtu31*PDstandard4th2gt12 - dgtu113*gtu32*PDstandard4th2gt12 - dgtu333*gtu11*PDstandard4th2gt13 - 
        2*dgtu313*gtu31*PDstandard4th2gt13 - dgtu113*gtu33*PDstandard4th2gt13 - dgtu323*gtu21*PDstandard4th2gt22 - 
        dgtu213*gtu32*PDstandard4th2gt22 - dgtu333*gtu21*PDstandard4th2gt23 - dgtu323*gtu31*PDstandard4th2gt23 - 
        dgtu313*gtu32*PDstandard4th2gt23 - dgtu213*gtu33*PDstandard4th2gt23 - dgtu333*gtu31*PDstandard4th2gt33 - 
        dgtu313*gtu33*PDstandard4th2gt33 - PDstandard4th23gt13*SQR(gtu31);
    
    ddgtu2223  =  -2*(gtu21*gtu22*PDstandard4th23gt12 + gtu21*gtu32*PDstandard4th23gt13 + gtu22*gtu32*PDstandard4th23gt23 + 
           dgtu213*gtu21*PDstandard4th2gt11 + dgtu223*gtu21*PDstandard4th2gt12 + dgtu213*gtu22*PDstandard4th2gt12 + 
           dgtu323*gtu21*PDstandard4th2gt13 + dgtu213*gtu32*PDstandard4th2gt13 + dgtu223*gtu22*PDstandard4th2gt22 + 
           dgtu323*gtu22*PDstandard4th2gt23 + dgtu223*gtu32*PDstandard4th2gt23 + dgtu323*gtu32*PDstandard4th2gt33) - 
        PDstandard4th23gt11*SQR(gtu21) - PDstandard4th23gt22*SQR(gtu22) - PDstandard4th23gt33*SQR(gtu32);
    
    ddgtu3223  =  -(gtu21*gtu31*PDstandard4th23gt11) - gtu22*gtu31*PDstandard4th23gt12 - gtu21*gtu32*PDstandard4th23gt12 - 
        gtu31*gtu32*PDstandard4th23gt13 - gtu21*gtu33*PDstandard4th23gt13 - gtu22*gtu32*PDstandard4th23gt22 - 
        gtu22*gtu33*PDstandard4th23gt23 - gtu32*gtu33*PDstandard4th23gt33 - dgtu313*gtu21*PDstandard4th2gt11 - 
        dgtu213*gtu31*PDstandard4th2gt11 - dgtu323*gtu21*PDstandard4th2gt12 - dgtu313*gtu22*PDstandard4th2gt12 - 
        dgtu223*gtu31*PDstandard4th2gt12 - dgtu213*gtu32*PDstandard4th2gt12 - dgtu333*gtu21*PDstandard4th2gt13 - 
        dgtu323*gtu31*PDstandard4th2gt13 - dgtu313*gtu32*PDstandard4th2gt13 - dgtu213*gtu33*PDstandard4th2gt13 - 
        dgtu323*gtu22*PDstandard4th2gt22 - dgtu223*gtu32*PDstandard4th2gt22 - dgtu333*gtu22*PDstandard4th2gt23 - 
        2*dgtu323*gtu32*PDstandard4th2gt23 - dgtu223*gtu33*PDstandard4th2gt23 - dgtu333*gtu32*PDstandard4th2gt33 - 
        dgtu323*gtu33*PDstandard4th2gt33 - PDstandard4th23gt23*SQR(gtu32);
    
    ddgtu3323  =  -2*(gtu31*gtu32*PDstandard4th23gt12 + gtu31*gtu33*PDstandard4th23gt13 + gtu32*gtu33*PDstandard4th23gt23 + 
           dgtu313*gtu31*PDstandard4th2gt11 + dgtu323*gtu31*PDstandard4th2gt12 + dgtu313*gtu32*PDstandard4th2gt12 + 
           dgtu333*gtu31*PDstandard4th2gt13 + dgtu313*gtu33*PDstandard4th2gt13 + dgtu323*gtu32*PDstandard4th2gt22 + 
           dgtu333*gtu32*PDstandard4th2gt23 + dgtu323*gtu33*PDstandard4th2gt23 + dgtu333*gtu33*PDstandard4th2gt33) - 
        PDstandard4th23gt11*SQR(gtu31) - PDstandard4th23gt22*SQR(gtu32) - PDstandard4th23gt33*SQR(gtu33);
    
    ddgtu3133  =  -(gtu11*gtu31*PDstandard4th33gt11) - gtu21*gtu31*PDstandard4th33gt12 - gtu11*gtu32*PDstandard4th33gt12 - 
        gtu11*gtu33*PDstandard4th33gt13 - gtu21*gtu32*PDstandard4th33gt22 - gtu31*gtu32*PDstandard4th33gt23 - 
        gtu21*gtu33*PDstandard4th33gt23 - gtu31*gtu33*PDstandard4th33gt33 - dgtu313*gtu11*PDstandard4th3gt11 - 
        dgtu113*gtu31*PDstandard4th3gt11 - dgtu323*gtu11*PDstandard4th3gt12 - dgtu313*gtu21*PDstandard4th3gt12 - 
        dgtu213*gtu31*PDstandard4th3gt12 - dgtu113*gtu32*PDstandard4th3gt12 - dgtu333*gtu11*PDstandard4th3gt13 - 
        2*dgtu313*gtu31*PDstandard4th3gt13 - dgtu113*gtu33*PDstandard4th3gt13 - dgtu323*gtu21*PDstandard4th3gt22 - 
        dgtu213*gtu32*PDstandard4th3gt22 - dgtu333*gtu21*PDstandard4th3gt23 - dgtu323*gtu31*PDstandard4th3gt23 - 
        dgtu313*gtu32*PDstandard4th3gt23 - dgtu213*gtu33*PDstandard4th3gt23 - dgtu333*gtu31*PDstandard4th3gt33 - 
        dgtu313*gtu33*PDstandard4th3gt33 - PDstandard4th33gt13*SQR(gtu31);
    
    ddgtu3233  =  -(gtu21*gtu31*PDstandard4th33gt11) - gtu22*gtu31*PDstandard4th33gt12 - gtu21*gtu32*PDstandard4th33gt12 - 
        gtu31*gtu32*PDstandard4th33gt13 - gtu21*gtu33*PDstandard4th33gt13 - gtu22*gtu32*PDstandard4th33gt22 - 
        gtu22*gtu33*PDstandard4th33gt23 - gtu32*gtu33*PDstandard4th33gt33 - dgtu313*gtu21*PDstandard4th3gt11 - 
        dgtu213*gtu31*PDstandard4th3gt11 - dgtu323*gtu21*PDstandard4th3gt12 - dgtu313*gtu22*PDstandard4th3gt12 - 
        dgtu223*gtu31*PDstandard4th3gt12 - dgtu213*gtu32*PDstandard4th3gt12 - dgtu333*gtu21*PDstandard4th3gt13 - 
        dgtu323*gtu31*PDstandard4th3gt13 - dgtu313*gtu32*PDstandard4th3gt13 - dgtu213*gtu33*PDstandard4th3gt13 - 
        dgtu323*gtu22*PDstandard4th3gt22 - dgtu223*gtu32*PDstandard4th3gt22 - dgtu333*gtu22*PDstandard4th3gt23 - 
        2*dgtu323*gtu32*PDstandard4th3gt23 - dgtu223*gtu33*PDstandard4th3gt23 - dgtu333*gtu32*PDstandard4th3gt33 - 
        dgtu323*gtu33*PDstandard4th3gt33 - PDstandard4th33gt23*SQR(gtu32);
    
    ddgtu3333  =  -2*(gtu31*gtu32*PDstandard4th33gt12 + gtu31*gtu33*PDstandard4th33gt13 + gtu32*gtu33*PDstandard4th33gt23 + 
           dgtu313*gtu31*PDstandard4th3gt11 + dgtu323*gtu31*PDstandard4th3gt12 + dgtu313*gtu32*PDstandard4th3gt12 + 
           dgtu333*gtu31*PDstandard4th3gt13 + dgtu313*gtu33*PDstandard4th3gt13 + dgtu323*gtu32*PDstandard4th3gt22 + 
           dgtu333*gtu32*PDstandard4th3gt23 + dgtu323*gtu33*PDstandard4th3gt23 + dgtu333*gtu33*PDstandard4th3gt33) - 
        PDstandard4th33gt11*SQR(gtu31) - PDstandard4th33gt22*SQR(gtu32) - PDstandard4th33gt33*SQR(gtu33);
    
    Gt111  =  khalf*(gtu11*PDstandard4th1gt11 + 2*(gtu21*PDstandard4th1gt12 + gtu31*PDstandard4th1gt13) - 
          gtu21*PDstandard4th2gt11 - gtu31*PDstandard4th3gt11);
    
    Gt211  =  khalf*(gtu21*PDstandard4th1gt11 + 2*(gtu22*PDstandard4th1gt12 + gtu32*PDstandard4th1gt13) - 
          gtu22*PDstandard4th2gt11 - gtu32*PDstandard4th3gt11);
    
    Gt311  =  khalf*(gtu31*PDstandard4th1gt11 + 2*(gtu32*PDstandard4th1gt12 + gtu33*PDstandard4th1gt13) - 
          gtu32*PDstandard4th2gt11 - gtu33*PDstandard4th3gt11);
    
    Gt112  =  khalf*(gtu21*PDstandard4th1gt22 + gtu11*PDstandard4th2gt11 + 
          gtu31*(PDstandard4th1gt23 + PDstandard4th2gt13 - PDstandard4th3gt12));
    
    Gt212  =  khalf*(gtu22*PDstandard4th1gt22 + gtu21*PDstandard4th2gt11 + 
          gtu32*(PDstandard4th1gt23 + PDstandard4th2gt13 - PDstandard4th3gt12));
    
    Gt312  =  khalf*(gtu32*PDstandard4th1gt22 + gtu31*PDstandard4th2gt11 + 
          gtu33*(PDstandard4th1gt23 + PDstandard4th2gt13 - PDstandard4th3gt12));
    
    Gt113  =  khalf*(gtu31*PDstandard4th1gt33 + gtu11*PDstandard4th3gt11 + 
          gtu21*(PDstandard4th1gt23 - PDstandard4th2gt13 + PDstandard4th3gt12));
    
    Gt213  =  khalf*(gtu32*PDstandard4th1gt33 + gtu21*PDstandard4th3gt11 + 
          gtu22*(PDstandard4th1gt23 - PDstandard4th2gt13 + PDstandard4th3gt12));
    
    Gt313  =  khalf*(gtu33*PDstandard4th1gt33 + gtu31*PDstandard4th3gt11 + 
          gtu32*(PDstandard4th1gt23 - PDstandard4th2gt13 + PDstandard4th3gt12));
    
    Gt122  =  khalf*(gtu11*(-PDstandard4th1gt22 + 2*PDstandard4th2gt12) + gtu21*PDstandard4th2gt22 + 
          gtu31*(2*PDstandard4th2gt23 - PDstandard4th3gt22));
    
    Gt222  =  khalf*(gtu21*(-PDstandard4th1gt22 + 2*PDstandard4th2gt12) + gtu22*PDstandard4th2gt22 + 
          gtu32*(2*PDstandard4th2gt23 - PDstandard4th3gt22));
    
    Gt322  =  khalf*(gtu31*(-PDstandard4th1gt22 + 2*PDstandard4th2gt12) + gtu32*PDstandard4th2gt22 + 
          gtu33*(2*PDstandard4th2gt23 - PDstandard4th3gt22));
    
    Gt123  =  khalf*(gtu31*PDstandard4th2gt33 + gtu11*(-PDstandard4th1gt23 + PDstandard4th2gt13 + PDstandard4th3gt12) + 
          gtu21*PDstandard4th3gt22);
    
    Gt223  =  khalf*(gtu32*PDstandard4th2gt33 + gtu21*(-PDstandard4th1gt23 + PDstandard4th2gt13 + PDstandard4th3gt12) + 
          gtu22*PDstandard4th3gt22);
    
    Gt323  =  khalf*(gtu33*PDstandard4th2gt33 + gtu31*(-PDstandard4th1gt23 + PDstandard4th2gt13 + PDstandard4th3gt12) + 
          gtu32*PDstandard4th3gt22);
    
    Gt133  =  khalf*(-(gtu11*PDstandard4th1gt33) - gtu21*PDstandard4th2gt33 + 2*gtu11*PDstandard4th3gt13 + 
          2*gtu21*PDstandard4th3gt23 + gtu31*PDstandard4th3gt33);
    
    Gt233  =  khalf*(-(gtu21*PDstandard4th1gt33) - gtu22*PDstandard4th2gt33 + 2*gtu21*PDstandard4th3gt13 + 
          2*gtu22*PDstandard4th3gt23 + gtu32*PDstandard4th3gt33);
    
    Gt333  =  khalf*(-(gtu31*PDstandard4th1gt33) - gtu32*PDstandard4th2gt33 + 2*gtu31*PDstandard4th3gt13 + 
          2*gtu32*PDstandard4th3gt23 + gtu33*PDstandard4th3gt33);
    
    Rt11  =  -(gtu11*khalf*PDstandard4th11gt11) + gtu21*
         (2*Gt211*Gt212*gt22L + 6*Gt112*gt13L*Gt311 + 4*Gt113*gt11L*Gt312 + 4*gt13L*Gt312*Gt313 + 4*gt13L*Gt211*Gt322 + 
           4*gt13L*Gt311*Gt323 + 2*Gt311*Gt312*gt33L - PDstandard4th12gt11) - gtu31*PDstandard4th13gt11 + 
        2*gt11L*PDstandard4th1Xt1 + gt12L*(6*Gt111*Gt212*gtu21 + 4*Gt211*Gt222*gtu21 + 4*Gt212*Gt222*gtu22 + 
           6*Gt113*Gt211*gtu31 + 6*Gt113*Gt212*gtu32 + 6*Gt113*Gt213*gtu33 + 2*PDstandard4th1Xt2) + 
        gt13L*(6*Gt111*Gt312*gtu21 + 4*Gt212*Gt312*gtu21 + 6*Gt112*Gt312*gtu22 + 6*Gt113*Gt311*gtu31 + 
           6*Gt113*Gt312*gtu32 + 6*Gt113*Gt313*gtu33 + 2*PDstandard4th1Xt3) - gtu22*khalf*PDstandard4th22gt11 + 
        gtu32*(2*Gt212*Gt213*gt22L + 4*gt12L*Gt233*Gt312 + 2*Gt213*gt23L*Gt312 + 6*Gt112*gt13L*Gt313 + 
           4*gt13L*Gt213*Gt322 + 4*gt13L*Gt212*Gt323 + 4*gt13L*Gt313*Gt323 + 4*gt13L*Gt312*Gt333 + 2*Gt312*Gt313*gt33L - 
           PDstandard4th23gt11) - gtu33*khalf*PDstandard4th33gt11 + 
        Gt111*(10*Gt113*gt11L*gtu31 + 6*gt12L*Gt213*gtu31 + 2*gt11L*Xt1L) + 
        Gt211*(4*Gt112*gt11L*gtu11 + 6*Gt111*gt12L*gtu11 + 2*gt23L*Gt311*gtu11 + 4*gt13L*Gt312*gtu11 + 
           6*Gt112*gt12L*gtu21 + 4*gt11L*Gt123*gtu31 + 4*gt12L*Gt223*gtu31 + 2*Gt213*gt22L*gtu31 + 2*gt12L*Xt1L) + 
        Gt311*(6*Gt111*gt13L*gtu11 + 4*gt12L*Gt213*gtu11 + 4*gt13L*Gt313*gtu11 + 4*gt11L*Gt123*gtu21 + 
           4*gt12L*Gt223*gtu21 + 2*Gt212*gt23L*gtu21 + 4*gt11L*Gt133*gtu31 + 4*gt12L*Gt233*gtu31 + 2*Gt213*gt23L*gtu31 + 
           2*Gt313*gt33L*gtu31 + 2*gt13L*Xt1L) + 2*gt12L*Gt212*Xt2L + 
        Gt112*(10*Gt111*gt11L*gtu21 + 4*gt11L*Gt212*gtu21 + 6*gt12L*Gt212*gtu22 + 4*gt11L*Gt213*gtu31 + 
           10*Gt113*gt11L*gtu32 + 2*gt11L*Xt2L) + Gt312*
         (4*gt12L*Gt213*gtu21 + 2*Gt211*gt23L*gtu21 + 4*gt11L*Gt123*gtu22 + 4*gt12L*Gt223*gtu22 + 2*Gt212*gt23L*gtu22 + 
           4*gt13L*Gt213*gtu31 + 4*gt11L*Gt133*gtu32 + 2*gt13L*Xt2L) + 2*Gt113*gt11L*Xt3L + 
        Gt213*(4*gt11L*Gt122*gtu32 + 6*Gt112*gt12L*gtu32 + 4*gt11L*Gt123*gtu33 + 2*gt12L*Xt3L) + 
        Gt313*(4*Gt113*gt11L*gtu31 + 6*Gt111*gt13L*gtu31 + 4*gt12L*Gt213*gtu31 + 2*Gt211*gt23L*gtu31 + 
           4*gt11L*Gt123*gtu32 + 4*gt12L*Gt223*gtu32 + 2*Gt212*gt23L*gtu32 + 4*gt11L*Gt133*gtu33 + 4*gt12L*Gt233*gtu33 + 
           2*Gt213*gt23L*gtu33 + 2*gt13L*Xt3L) + 5*gt11L*gtu11*SQR(Gt111) + 5*gt11L*gtu22*SQR(Gt112) + 
        5*gt11L*gtu33*SQR(Gt113) + gt22L*gtu11*SQR(Gt211) + gt22L*gtu22*SQR(Gt212) + 
        4*(gt12L*Gt211*Gt212*gtu11 + Gt113*gt11L*Gt311*gtu11 + gt11L*Gt122*Gt211*gtu21 + gt11L*Gt122*Gt212*gtu22 + 
           gt13L*Gt212*Gt322*gtu22 + gt13L*Gt312*Gt323*gtu22 + gt12L*Gt212*Gt213*gtu31 + gt13L*Gt211*Gt323*gtu31 + 
           gt13L*Gt311*Gt333*gtu31 + gt11L*Gt123*Gt212*gtu32 + gt12L*Gt213*Gt222*gtu32 + gt12L*Gt212*Gt223*gtu32 + 
           gt12L*Gt213*Gt223*gtu33 + gt13L*Gt213*Gt323*gtu33 + gt13L*Gt313*Gt333*gtu33 + gt12L*gtu21*SQR(Gt212)) + 
        gt22L*gtu33*SQR(Gt213) + gt33L*gtu11*SQR(Gt311) + gt33L*gtu22*SQR(Gt312) + 4*gt13L*gtu31*SQR(Gt313) + 
        gt33L*gtu33*SQR(Gt313);
    
    Rt12  =  gt22L*PDstandard4th1Xt2 + gt23L*PDstandard4th1Xt3 + gt11L*PDstandard4th2Xt1 + 
        gt12L*(PDstandard4th1Xt1 + PDstandard4th2Xt2) + gt13L*PDstandard4th2Xt3 + 
        khalf*(-(gtu11*PDstandard4th11gt12) - 2*gtu21*PDstandard4th12gt12 - 2*gtu31*PDstandard4th13gt12 - 
           gtu22*PDstandard4th22gt12 - 2*gtu32*PDstandard4th23gt12 - gtu33*PDstandard4th33gt12) + 
        (Gt112*gt11L + gt12L*(Gt111 + Gt212) + Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)*Xt1L + 
        (gt11L*Gt122 + gt12L*(Gt112 + Gt222) + Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322)*Xt2L + 
        (gt11L*Gt123 + gt12L*(Gt113 + Gt223) + Gt213*gt22L + gt23L*Gt313 + gt13L*Gt323)*Xt3L + 
        gtu11*(gt11L*(Gt112*(3*Gt111 + 2*Gt212) + 2*Gt113*Gt312) + 
           Gt111*(gt12L*Gt212 + 2*(Gt211*gt22L + gt23L*Gt311) + gt13L*Gt312) + 
           Gt211*(5*Gt112*gt12L + 3*(Gt212*gt22L + gt23L*Gt312)) + 
           Gt311*(3*Gt112*gt13L + Gt212*gt23L + 2*(Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + Gt312*gt33L) + 
           2*(Gt312*(gt12L*Gt213 + gt13L*(Gt212 + Gt313)) + gt12L*(SQR(Gt111) + SQR(Gt212)))) + 
        gtu22*(Gt112*(3*gt11L*Gt122 + gt12L*Gt222 + 2*(Gt212*gt22L + gt23L*Gt312) + gt13L*Gt322) + 
           Gt222*(3*Gt212*gt22L + gt23L*Gt312 + 2*(gt11L*Gt122 + gt13L*Gt322)) + 
           Gt312*(3*Gt122*gt13L + 2*(Gt223*gt22L + gt23L*Gt323)) + 
           Gt322*(3*Gt212*gt23L + 2*(gt11L*Gt123 + gt13L*Gt323) + Gt312*gt33L) + 
           gt12L*(5*Gt122*Gt212 + 2*(Gt123*Gt312 + Gt223*Gt322 + SQR(Gt112) + SQR(Gt222)))) + 
        gtu33*(Gt113*(3*gt11L*Gt123 + gt12L*Gt223 + 2*(Gt213*gt22L + gt23L*Gt313) + gt13L*Gt323) + 
           Gt223*(3*Gt213*gt22L + gt23L*Gt313 + 2*(gt11L*Gt123 + gt13L*Gt323)) + 
           Gt313*(3*Gt123*gt13L + 2*(gt22L*Gt233 + gt23L*Gt333)) + 
           Gt323*(3*Gt213*gt23L + 2*(gt11L*Gt133 + gt13L*Gt333) + Gt313*gt33L) + 
           gt12L*(5*Gt123*Gt213 + 2*(Gt133*Gt313 + Gt233*Gt323 + SQR(Gt113) + SQR(Gt223)))) + 
        gtu21*(Gt111*(3*gt11L*Gt122 + gt12L*(2*Gt112 + Gt222) + gt13L*Gt322) + 
           2*(Gt112*(gt11L*Gt222 + Gt211*gt22L) + Gt223*gt22L*Gt311 + Gt111*(Gt112*gt12L + Gt212*gt22L + gt23L*Gt312) + 
              (Gt113*gt11L + gt13L*Gt313)*Gt322 + gt13L*(Gt222*Gt312 + Gt212*Gt322)) + 
           gt12L*(5*Gt122*Gt211 + Gt212*(4*Gt112 + 2*Gt222) + 
              2*(Gt212*(Gt112 + Gt222) + Gt123*Gt311 + (Gt113 + Gt223)*Gt312 + Gt213*Gt322)) + 
           Gt312*(2*(Gt213*gt22L + gt23L*Gt313) + gt13L*(4*Gt112 + 2*Gt323)) + 
           gt23L*(4*Gt212*Gt312 + Gt311*(Gt222 + 2*(Gt112 + Gt323))) + 
           gt11L*(2*(Gt122*Gt212 + Gt123*Gt312) + 3*SQR(Gt112)) + 
           3*(Gt122*gt13L*Gt311 + Gt211*(Gt222*gt22L + gt23L*Gt322) + gt22L*SQR(Gt212)) + gt33L*(Gt311*Gt322 + SQR(Gt312)))\
         + gtu31*(Gt111*(3*gt11L*Gt123 + gt12L*(2*Gt113 + Gt223) + gt13L*Gt323) + 
           gt11L*(2*(Gt123*Gt212 + Gt133*Gt312) + Gt113*(3*Gt112 + 2*Gt323)) + 
           3*((Gt212*Gt213 + Gt211*Gt223)*gt22L + Gt123*gt13L*Gt311 + gt23L*(Gt213*Gt312 + Gt211*Gt323)) + 
           gt12L*(5*Gt123*Gt211 + 3*Gt112*Gt213 + Gt212*(Gt113 + 2*Gt223) + 
              2*(Gt212*Gt223 + Gt133*Gt311 + Gt233*Gt312 + Gt113*Gt313 + Gt213*(Gt112 + Gt323))) + 
           gt13L*(3*Gt112*Gt313 + 2*Gt212*Gt323 + Gt312*(Gt113 + 2*Gt333)) + 
           gt23L*(Gt212*Gt313 + Gt311*(Gt223 + 2*(Gt113 + Gt333))) + (Gt312*Gt313 + Gt311*Gt323)*gt33L + 
           2*(Gt223*(Gt112*gt11L + gt13L*Gt312) + gt22L*(Gt113*Gt211 + Gt233*Gt311 + Gt213*Gt313) + 
              Gt111*(Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + gt13L*Gt313*Gt323 + gt23L*SQR(Gt313))) + 
        gtu32*(3*((Gt213*Gt222 + Gt212*Gt223)*gt22L + gt13L*(Gt123*Gt312 + Gt122*Gt313) + Gt213*gt23L*Gt322) + 
           Gt112*(3*gt11L*Gt123 + gt12L*(2*Gt113 + Gt223) + gt13L*Gt323) + 
           gt12L*(5*(Gt123*Gt212 + Gt122*Gt213) + Gt222*(Gt113 + 4*Gt223) + 
              2*(Gt133*Gt312 + Gt123*Gt313 + Gt233*Gt322 + Gt223*Gt323)) + 
           gt11L*(3*Gt113*Gt122 + 2*(Gt133*Gt322 + Gt123*(Gt222 + Gt323))) + 
           gt13L*(2*Gt222*Gt323 + Gt322*(Gt113 + 2*Gt333)) + 
           gt23L*(Gt222*Gt313 + 3*Gt212*Gt323 + Gt312*(Gt223 + 2*(Gt113 + Gt333))) + (Gt313*Gt322 + Gt312*Gt323)*gt33L + 
           2*(gt22L*(Gt113*Gt212 + Gt233*Gt312) + Gt112*(Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + 
              Gt223*(gt11L*Gt122 + gt22L*Gt313 + gt13L*Gt322) + gt23L*Gt313*Gt323 + gt13L*SQR(Gt323)));
    
    Rt13  =  gt23L*PDstandard4th1Xt2 + gt33L*PDstandard4th1Xt3 + 
        khalf*(-(gtu11*PDstandard4th11gt13) - 2*gtu21*PDstandard4th12gt13 - 2*gtu31*PDstandard4th13gt13 - 
           gtu22*PDstandard4th22gt13 - 2*gtu32*PDstandard4th23gt13 - gtu33*PDstandard4th33gt13) + gt11L*PDstandard4th3Xt1 + 
        gt12L*PDstandard4th3Xt2 + gt13L*(PDstandard4th1Xt1 + PDstandard4th3Xt3) + 
        (Gt113*gt11L + gt12L*Gt213 + Gt211*gt23L + gt13L*(Gt111 + Gt313) + Gt311*gt33L)*Xt1L + 
        (gt11L*Gt123 + gt12L*Gt223 + Gt212*gt23L + gt13L*(Gt112 + Gt323) + Gt312*gt33L)*Xt2L + 
        (gt11L*Gt133 + gt12L*Gt233 + Gt213*gt23L + gt13L*(Gt113 + Gt333) + Gt313*gt33L)*Xt3L + 
        gtu21*((Gt212*Gt213 + Gt211*Gt223)*gt22L + Gt123*(Gt111*gt11L + gt12L*Gt211 + 3*gt13L*Gt311) + 
           gt12L*(3*Gt113*Gt212 + Gt112*Gt213 + Gt223*(Gt111 + 2*Gt313)) + 
           gt11L*(3*Gt112*Gt113 + 2*(Gt122*Gt213 + Gt123*Gt313)) + gt13L*(Gt112*Gt313 + Gt111*Gt323) + 
           gt23L*(Gt213*Gt312 + Gt212*Gt313 + Gt211*Gt323) + Gt312*Gt313*gt33L + 
           3*(Gt113*gt13L*Gt312 + Gt311*(Gt223*gt23L + Gt323*gt33L)) + 
           2*(Gt213*(gt12L*Gt222 + gt13L*Gt322) + gt11L*(Gt111*Gt123 + Gt112*Gt223 + Gt113*Gt323) + 
              gt12L*(Gt123*Gt211 + Gt212*Gt223 + Gt213*Gt323) + 
              gt13L*(Gt112*(Gt111 + Gt212) + Gt123*Gt311 + (Gt113 + Gt223)*Gt312 + 2*Gt313*Gt323) + 
              (Gt112*Gt311 + Gt312*(Gt212 + Gt313))*gt33L + Gt111*(Gt112*gt13L + Gt212*gt23L + Gt312*gt33L) + 
              Gt211*(Gt122*gt13L + Gt222*gt23L + Gt322*gt33L) + gt23L*(Gt112*Gt211 + Gt213*Gt312 + SQR(Gt212)))) + 
        gtu32*(gt11L*(Gt123*(3*Gt113 + 2*Gt223) + 2*Gt133*Gt323) + Gt223*(gt23L*(2*Gt212 + 3*Gt313) + 2*gt13L*Gt323) + 
           Gt233*(Gt212*gt22L + 3*gt23L*Gt312 + 2*(gt11L*Gt122 + gt12L*Gt323)) + 
           Gt123*(3*gt12L*Gt213 + gt13L*(2*Gt212 + 5*Gt313) + 2*gt11L*Gt333) + 
           gt13L*(5*Gt133*Gt312 + 2*Gt233*Gt322 + 4*Gt323*Gt333) + 3*(Gt313*Gt323 + Gt312*Gt333)*gt33L + 
           Gt212*(gt23L*Gt333 + 2*Gt323*gt33L) + Gt113*(gt12L*Gt223 + gt13L*Gt323 + 2*(Gt212*gt23L + Gt312*gt33L)) + 
           Gt112*(3*gt11L*Gt133 + gt12L*Gt233 + gt13L*(2*Gt113 + Gt333) + 2*(Gt113*gt13L + Gt213*gt23L + Gt313*gt33L)) + 
           Gt213*(Gt223*gt22L + gt23L*Gt323 + 2*(Gt122*gt13L + Gt222*gt23L + Gt322*gt33L)) + 
           gt12L*(3*Gt133*Gt212 + 2*(Gt222*Gt233 + Gt223*Gt333 + SQR(Gt223)))) + 
        gtu31*(Gt133*(Gt111*gt11L + gt12L*Gt211 + 3*gt13L*Gt311) + Gt233*(Gt211*gt22L + 3*gt23L*Gt311) + 
           gt12L*(4*Gt113*Gt213 + Gt233*(Gt111 + 2*Gt313)) + gt13L*(4*Gt113*Gt313 + 2*Gt213*Gt323 + Gt111*Gt333) + 
           Gt211*(gt23L*Gt333 + 2*Gt323*gt33L) + gt11L*(2*(Gt123*Gt213 + Gt133*Gt313) + 3*SQR(Gt113)) + gt22L*SQR(Gt213) + 
           gt33L*(3*Gt311*Gt333 + SQR(Gt313)) + 2*(Gt211*(Gt123*gt13L + Gt223*gt23L) + Gt213*(gt12L*Gt223 + gt23L*Gt313) + 
              gt23L*(Gt113*Gt211 + Gt213*(Gt212 + Gt313)) + gt11L*(Gt111*Gt133 + Gt112*Gt233 + Gt113*Gt333) + 
              gt12L*(Gt133*Gt211 + Gt212*Gt233 + Gt213*Gt333) + 
              gt13L*(Gt112*Gt213 + Gt133*Gt311 + Gt233*Gt312 + Gt113*(Gt111 + Gt313) + 2*Gt313*Gt333) + 
              Gt111*(Gt113*gt13L + Gt213*gt23L + Gt313*gt33L) + gt33L*(Gt113*Gt311 + Gt213*Gt312 + SQR(Gt313)))) + 
        gtu11*(Gt211*(3*Gt113*gt12L + 2*Gt112*gt13L + Gt213*gt22L + gt23L*Gt313) + 
           gt11L*(2*Gt112*Gt213 + Gt113*(3*Gt111 + 2*Gt313)) + 
           Gt111*(gt12L*Gt213 + gt13L*Gt313 + 2*(Gt211*gt23L + Gt311*gt33L)) + 
           Gt311*(5*Gt113*gt13L + 3*(Gt213*gt23L + Gt313*gt33L)) + 
           2*(Gt212*(gt12L*Gt213 + Gt211*gt23L) + Gt213*(gt13L*Gt312 + gt12L*Gt313) + Gt211*Gt312*gt33L + 
              gt13L*(SQR(Gt111) + SQR(Gt313)))) + gtu22*
         (Gt212*(3*Gt123*gt12L + 2*Gt122*gt13L + Gt223*gt22L + gt23L*Gt323) + 
           gt11L*(2*Gt122*Gt223 + Gt123*(3*Gt112 + 2*Gt323)) + 
           Gt112*(gt12L*Gt223 + gt13L*Gt323 + 2*(Gt212*gt23L + Gt312*gt33L)) + 
           Gt312*(5*Gt123*gt13L + 3*(Gt223*gt23L + Gt323*gt33L)) + 
           2*(Gt222*(gt12L*Gt223 + Gt212*gt23L) + Gt223*(gt13L*Gt322 + gt12L*Gt323) + Gt212*Gt322*gt33L + 
              gt13L*(SQR(Gt112) + SQR(Gt323)))) + gtu33*
         (gt12L*(3*Gt133*Gt213 + 2*Gt233*(Gt223 + Gt333)) + Gt213*(gt22L*Gt233 + gt23L*(2*Gt223 + Gt333)) + 
           Gt113*(3*gt11L*Gt133 + gt12L*Gt233 + gt13L*Gt333 + 2*(Gt213*gt23L + Gt313*gt33L)) + 
           Gt313*(5*Gt133*gt13L + 3*(Gt233*gt23L + Gt333*gt33L)) + 
           2*(Gt123*(gt13L*Gt213 + gt11L*Gt233) + gt11L*Gt133*Gt333 + Gt323*(gt13L*Gt233 + Gt213*gt33L) + 
              gt13L*(SQR(Gt113) + SQR(Gt333))));
    
    Rt22  =  6*(Gt122*gt12L*Gt212*gtu21 + Gt112*gt12L*Gt222*gtu21 + Gt122*gt12L*Gt222*gtu22 + Gt123*gt12L*Gt212*gtu31 + 
           Gt123*gt12L*Gt222*gtu32 + Gt123*gt12L*Gt223*gtu33) - gtu11*khalf*PDstandard4th11gt22 + 
        gtu21*(10*Gt212*Gt222*gt22L + 4*Gt122*gt23L*Gt311 + 2*Gt122*gt13L*Gt312 + 6*Gt222*gt23L*Gt312 + 
           4*Gt113*gt12L*Gt322 + 4*gt23L*Gt312*Gt323 + 2*Gt312*Gt322*gt33L - PDstandard4th12gt22) + 
        gtu31*(10*Gt212*Gt223*gt22L + 2*Gt123*gt13L*Gt312 + 4*Gt112*gt23L*Gt313 + 4*Gt113*gt12L*Gt323 + 
           4*gt23L*Gt312*Gt333 + 2*Gt312*Gt323*gt33L - PDstandard4th13gt22) - gtu22*khalf*PDstandard4th22gt22 + 
        gtu32*(6*Gt223*gt23L*Gt322 + 2*Gt122*gt13L*Gt323 + 4*gt23L*Gt322*Gt333 + 2*Gt322*Gt323*gt33L - 
           PDstandard4th23gt22) + gt12L*(4*Gt111*Gt123*gtu31 + 6*Gt112*Gt223*gtu31 + 4*Gt113*Gt122*gtu32 + 
           4*Gt113*Gt123*gtu33 + 2*PDstandard4th2Xt1) + 
        gt22L*(4*Gt122*Gt213*gtu32 + 10*Gt222*Gt223*gtu32 + 4*Gt123*Gt213*gtu33 + 2*PDstandard4th2Xt2) + 
        gt23L*(6*Gt212*Gt322*gtu21 + 4*Gt313*Gt322*gtu21 + 6*Gt222*Gt322*gtu22 + 4*Gt123*Gt311*gtu31 + 
           6*Gt212*Gt323*gtu31 + 4*Gt313*Gt323*gtu31 + 4*Gt122*Gt313*gtu32 + 6*Gt223*Gt323*gtu33 + 2*PDstandard4th2Xt3) - 
        gtu33*khalf*PDstandard4th33gt22 + 2*Gt212*gt22L*Xt1L + 
        Gt112*(4*Gt111*gt12L*gtu11 + 6*gt12L*Gt212*gtu11 + 2*gt11L*Gt122*gtu21 + 4*Gt122*gt12L*gtu22 + 
           2*gt11L*Gt123*gtu31 + 4*Gt123*gt12L*gtu32 + 2*gt12L*Xt1L) + 
        Gt312*(4*Gt113*gt12L*gtu11 + 2*Gt112*gt13L*gtu11 + 4*Gt213*gt22L*gtu11 + 6*Gt212*gt23L*gtu11 + 
           4*gt23L*Gt313*gtu11 + 4*Gt123*gt12L*gtu21 + 4*Gt122*gt23L*gtu22 + 4*gt12L*Gt133*gtu31 + 4*gt22L*Gt233*gtu31 + 
           6*Gt223*gt23L*gtu31 + 4*Gt123*gt23L*gtu32 + 2*gt23L*Xt1L) + 2*Gt222*gt22L*Xt2L + 
        Gt122*(2*gt11L*Gt123*gtu32 + 6*gt12L*Gt223*gtu32 + 2*gt12L*Xt2L) + 
        Gt322*(2*Gt112*gt13L*gtu21 + 4*Gt213*gt22L*gtu21 + 4*Gt123*gt12L*gtu22 + 2*Gt122*gt13L*gtu22 + 
           4*Gt223*gt22L*gtu22 + 4*gt23L*Gt323*gtu22 + 4*gt12L*Gt133*gtu32 + 2*Gt123*gt13L*gtu32 + 4*gt22L*Gt233*gtu32 + 
           2*gt23L*Xt2L) + 2*Gt123*gt12L*Xt3L + 2*Gt223*gt22L*Xt3L + 
        Gt323*(2*Gt112*gt13L*gtu31 + 4*Gt213*gt22L*gtu31 + 4*Gt123*gt12L*gtu32 + 4*Gt223*gt22L*gtu32 + 
           6*Gt222*gt23L*gtu32 + 4*gt12L*Gt133*gtu33 + 2*Gt123*gt13L*gtu33 + 4*gt22L*Gt233*gtu33 + 4*gt23L*Gt333*gtu33 + 
           2*gt23L*Xt3L) + gt11L*gtu11*SQR(Gt112) + 4*
         (Gt112*Gt211*gt22L*gtu11 + Gt112*gt23L*Gt311*gtu11 + Gt111*Gt122*gt12L*gtu21 + Gt122*Gt211*gt22L*gtu21 + 
           Gt112*Gt212*gt22L*gtu21 + Gt223*gt22L*Gt312*gtu21 + Gt112*gt23L*Gt312*gtu21 + Gt122*Gt212*gt22L*gtu22 + 
           Gt112*Gt113*gt12L*gtu31 + Gt123*Gt211*gt22L*gtu31 + Gt112*Gt213*gt22L*gtu31 + Gt123*Gt212*gt22L*gtu32 + 
           Gt123*gt23L*Gt313*gtu33 + gt12L*gtu21*SQR(Gt112)) + gt11L*gtu22*SQR(Gt122) + gt11L*gtu33*SQR(Gt123) + 
        5*gt22L*gtu11*SQR(Gt212) + 5*gt22L*gtu22*SQR(Gt222) + 5*gt22L*gtu33*SQR(Gt223) + gt33L*gtu11*SQR(Gt312) + 
        gt33L*gtu22*SQR(Gt322) + 4*gt23L*gtu32*SQR(Gt323) + gt33L*gtu33*SQR(Gt323);
    
    Rt23  =  gt13L*PDstandard4th2Xt1 + gt33L*PDstandard4th2Xt3 + 
        khalf*(-(gtu11*PDstandard4th11gt23) - 2*gtu21*PDstandard4th12gt23 - 2*gtu31*PDstandard4th13gt23 - 
           gtu22*PDstandard4th22gt23 - 2*gtu32*PDstandard4th23gt23 - gtu33*PDstandard4th33gt23) + gt12L*PDstandard4th3Xt1 + 
        gt22L*PDstandard4th3Xt2 + gt23L*(PDstandard4th2Xt2 + PDstandard4th3Xt3) + 
        (Gt113*gt12L + Gt112*gt13L + Gt213*gt22L + gt23L*(Gt212 + Gt313) + Gt312*gt33L)*Xt1L + 
        (Gt123*gt12L + Gt122*gt13L + Gt223*gt22L + gt23L*(Gt222 + Gt323) + Gt322*gt33L)*Xt2L + 
        (gt12L*Gt133 + Gt123*gt13L + gt22L*Gt233 + gt23L*(Gt223 + Gt333) + Gt323*gt33L)*Xt3L + 
        gtu21*(gt11L*(Gt113*Gt122 + Gt112*Gt123) + gt12L*
            (Gt113*Gt222 + 3*(Gt122*Gt213 + Gt112*Gt223) + Gt123*(Gt212 + 2*Gt313)) + 
           gt23L*(2*Gt123*Gt311 + Gt222*(4*Gt212 + Gt313) + 5*(Gt223*Gt312 + Gt213*Gt322) + (Gt212 + 4*Gt313)*Gt323) + 
           gt22L*(Gt212*(2*Gt113 + 3*Gt223) + 2*Gt223*Gt313 + Gt213*(3*Gt222 + 2*Gt323)) + 
           3*(Gt313*Gt322 + Gt312*Gt323)*gt33L + 2*(Gt111*(Gt123*gt12L + Gt122*gt13L) + gt13L*(Gt122*Gt212 + Gt112*Gt222) + 
              Gt211*(Gt123*gt22L + Gt122*gt23L) + Gt113*(gt23L*Gt312 + gt12L*Gt323) + 
              (Gt122*Gt311 + Gt222*Gt312 + Gt212*Gt322)*gt33L + Gt112*(Gt113*gt12L + Gt212*gt23L + Gt312*gt33L)) + 
           gt13L*(Gt122*Gt313 + 3*(Gt123*Gt312 + Gt113*Gt322) + Gt112*Gt323 + 2*SQR(Gt112))) + 
        gtu31*(Gt133*(Gt112*gt11L + gt12L*Gt212 + 2*(Gt211*gt22L + gt23L*Gt311)) + 
           gt23L*(5*Gt233*Gt312 + Gt223*(4*Gt212 + Gt313) + Gt213*(2*Gt112 + 5*Gt323) + Gt212*Gt333) + 
           Gt113*(gt11L*Gt123 + 2*(Gt112*gt13L + Gt213*gt22L + gt23L*Gt313) + gt13L*Gt323 + gt12L*(Gt223 + 2*Gt333)) + 
           Gt123*(3*gt12L*Gt213 + gt13L*(2*Gt212 + Gt313) + 2*Gt311*gt33L) + 
           Gt333*(Gt112*gt13L + 4*gt23L*Gt313 + 3*Gt312*gt33L) + 
           3*(Gt112*gt12L*Gt233 + gt22L*(Gt213*Gt223 + Gt212*Gt233) + Gt133*gt13L*Gt312 + Gt313*Gt323*gt33L) + 
           2*(Gt111*(gt12L*Gt133 + Gt123*gt13L) + Gt123*Gt211*gt23L + gt13L*(Gt112*Gt223 + Gt113*Gt323) + 
              Gt213*gt22L*Gt333 + (Gt223*Gt312 + Gt212*Gt323)*gt33L + Gt313*(gt12L*Gt133 + gt22L*Gt233 + Gt112*gt33L) + 
              gt12L*SQR(Gt113))) + gtu11*(Gt113*(Gt112*gt11L + gt12L*Gt212 + 2*(Gt211*gt22L + gt23L*Gt311)) + 
           (Gt112*gt13L + Gt212*gt23L)*Gt313 + Gt213*(5*gt23L*Gt312 + 2*gt22L*Gt313) + 
           3*(Gt213*(Gt112*gt12L + Gt212*gt22L) + Gt312*(Gt113*gt13L + Gt313*gt33L)) + 
           2*(Gt111*(Gt113*gt12L + Gt112*gt13L) + Gt113*gt12L*Gt313 + Gt212*Gt312*gt33L + 
              Gt112*(gt13L*Gt212 + Gt211*gt23L + Gt311*gt33L) + gt23L*(SQR(Gt212) + SQR(Gt313)))) + 
        gtu32*(Gt133*(gt11L*Gt122 + gt12L*Gt222 + 3*gt13L*Gt322) + gt12L*(4*Gt123*Gt223 + Gt122*Gt233 + 2*Gt133*Gt323) + 
           Gt233*(3*gt23L*Gt322 + 2*gt22L*Gt323) + gt23L*Gt323*(4*Gt223 + 2*Gt333) + Gt122*(gt13L*Gt333 + 2*Gt313*gt33L) + 
           Gt222*(gt22L*Gt233 + gt23L*(2*Gt223 + Gt333) + 2*Gt323*gt33L) + gt11L*SQR(Gt123) + 
           3*(Gt322*Gt333*gt33L + gt22L*SQR(Gt223)) + gt33L*SQR(Gt323) + 
           2*(Gt113*(Gt123*gt12L + Gt122*gt13L) + Gt133*(Gt112*gt12L + Gt212*gt22L + gt23L*Gt312) + 
              Gt233*(Gt122*gt12L + Gt222*gt22L + gt23L*Gt322) + gt13L*(Gt122*Gt223 + Gt123*(Gt112 + Gt323)) + 
              Gt223*gt22L*Gt333 + Gt123*(Gt213*gt22L + gt23L*Gt313 + gt13L*(Gt222 + Gt323) + gt12L*Gt333) + 
              gt23L*(Gt123*Gt212 + Gt122*Gt213 + Gt223*(Gt222 + Gt323) + Gt323*Gt333) + 
              gt33L*(Gt123*Gt312 + Gt223*Gt322 + SQR(Gt323)))) + 
        gtu22*(Gt123*(gt11L*Gt122 + gt12L*Gt222 + 2*(Gt212*gt22L + gt23L*Gt312)) + (Gt122*gt13L + Gt222*gt23L)*Gt323 + 
           Gt223*(5*gt23L*Gt322 + 2*gt22L*Gt323) + 3*
            (Gt223*(Gt122*gt12L + Gt222*gt22L) + Gt322*(Gt123*gt13L + Gt323*gt33L)) + 
           2*(Gt112*(Gt123*gt12L + Gt122*gt13L) + Gt123*gt12L*Gt323 + Gt222*Gt322*gt33L + 
              Gt122*(gt13L*Gt222 + Gt212*gt23L + Gt312*gt33L) + gt23L*(SQR(Gt222) + SQR(Gt323)))) + 
        gtu33*(Gt133*(gt11L*Gt123 + gt12L*Gt223 + 2*(Gt213*gt22L + gt23L*Gt313)) + (Gt123*gt13L + Gt223*gt23L)*Gt333 + 
           Gt233*(5*gt23L*Gt323 + 2*gt22L*Gt333) + 3*
            ((Gt123*gt12L + Gt223*gt22L)*Gt233 + Gt323*(Gt133*gt13L + Gt333*gt33L)) + 
           2*(Gt113*(gt12L*Gt133 + Gt123*gt13L) + gt12L*Gt133*Gt333 + Gt223*Gt323*gt33L + 
              Gt123*(gt13L*Gt223 + Gt213*gt23L + Gt313*gt33L) + gt23L*(SQR(Gt223) + SQR(Gt333))));
    
    Rt33  =  6*(Gt133*gt13L*Gt313*gtu31 + Gt233*gt23L*Gt313*gtu31 + Gt113*gt13L*Gt333*gtu31 + Gt213*gt23L*Gt333*gtu31 + 
           Gt133*gt13L*Gt323*gtu32 + Gt133*gt13L*Gt333*gtu33) + 
        gtu21*(4*Gt122*gt13L*Gt213 + 2*Gt113*gt12L*Gt223 + 4*Gt112*gt13L*Gt223 + 2*Gt213*Gt223*gt22L + 
           4*Gt123*Gt211*gt23L + 4*Gt212*Gt223*gt23L + 6*Gt123*gt13L*Gt313 + 6*Gt223*gt23L*Gt313 + 6*Gt113*gt13L*Gt323 + 
           6*Gt213*gt23L*Gt323 + 4*Gt123*Gt311*gt33L - PDstandard4th12gt33) - gtu31*PDstandard4th13gt33 + 
        gtu22*(4*Gt222*Gt223*gt23L + 6*Gt123*gt13L*Gt323 + 6*Gt223*gt23L*Gt323 + 4*Gt223*Gt322*gt33L - 
           khalf*PDstandard4th22gt33) + gtu32*(2*Gt223*gt22L*Gt233 + 4*Gt133*Gt212*gt23L + 6*Gt233*gt23L*Gt323 + 
           6*Gt123*gt13L*Gt333 + 6*Gt223*gt23L*Gt333 + 4*Gt123*Gt313*gt33L + 10*Gt323*Gt333*gt33L - PDstandard4th23gt33) - 
        gtu33*khalf*PDstandard4th33gt33 + 2*(gt12L*Gt133*Gt213*gtu31 + Gt113*gt12L*Gt233*gtu31 + gt12L*Gt133*Gt223*gtu32 + 
           gt13L*PDstandard4th3Xt1) + 2*gt23L*PDstandard4th3Xt2 + 
        gt33L*(4*Gt213*Gt322*gtu21 + 10*Gt313*Gt323*gtu21 + 4*Gt123*Gt312*gtu22 + 4*Gt133*Gt311*gtu31 + 
           4*Gt213*Gt323*gtu31 + 10*Gt313*Gt333*gtu31 + 4*Gt133*Gt312*gtu32 + 4*Gt133*Gt313*gtu33 + 2*PDstandard4th3Xt3) + 
        2*Gt213*gt23L*Xt1L + 2*Gt313*gt33L*Xt1L + Gt113*
         (4*Gt111*gt13L*gtu11 + 2*gt12L*Gt213*gtu11 + 2*gt11L*Gt123*gtu21 + 2*gt11L*Gt133*gtu31 + 4*Gt123*gt13L*gtu32 + 
           4*Gt133*gt13L*gtu33 + 2*gt13L*Xt1L) + 2*Gt223*gt23L*Xt2L + 2*Gt323*gt33L*Xt2L + 
        Gt123*(4*Gt111*gt13L*gtu21 + 2*gt12L*Gt213*gtu21 + 4*Gt112*gt13L*gtu22 + 2*gt12L*Gt223*gtu22 + 
           4*Gt212*gt23L*gtu22 + 4*gt13L*Gt213*gtu31 + 2*gt11L*Gt133*gtu32 + 4*gt13L*Gt233*gtu33 + 2*gt13L*Xt2L) + 
        2*Gt133*gt13L*Xt3L + 2*Gt333*gt33L*Xt3L + Gt233*
         (4*Gt112*gt13L*gtu31 + 2*Gt213*gt22L*gtu31 + 2*Gt123*gt12L*gtu32 + 2*gt12L*Gt133*gtu33 + 4*Gt223*gt23L*gtu33 + 
           6*gt23L*Gt333*gtu33 + 4*Gt323*gt33L*gtu33 + 2*gt23L*Xt3L) + 
        gtu11*(4*Gt212*Gt213*gt23L + 6*Gt113*gt13L*Gt313 + 6*Gt213*gt23L*Gt313 + 4*Gt113*Gt311*gt33L + 
           4*Gt213*Gt312*gt33L - khalf*PDstandard4th11gt33 + gt11L*SQR(Gt113)) + 
        4*(Gt112*gt13L*Gt213*gtu11 + Gt113*Gt211*gt23L*gtu11 + Gt112*Gt113*gt13L*gtu21 + Gt113*Gt212*gt23L*gtu21 + 
           Gt213*Gt222*gt23L*gtu21 + Gt113*Gt312*gt33L*gtu21 + Gt223*Gt312*gt33L*gtu21 + Gt122*gt13L*Gt223*gtu22 + 
           Gt111*Gt133*gt13L*gtu31 + Gt133*Gt211*gt23L*gtu31 + Gt113*Gt213*gt23L*gtu31 + Gt213*Gt223*gt23L*gtu31 + 
           Gt212*Gt233*gt23L*gtu31 + Gt233*Gt312*gt33L*gtu31 + Gt113*Gt313*gt33L*gtu31 + Gt112*Gt133*gt13L*gtu32 + 
           Gt123*gt13L*Gt223*gtu32 + Gt122*gt13L*Gt233*gtu32 + Gt123*Gt213*gt23L*gtu32 + Gt222*Gt233*gt23L*gtu32 + 
           Gt233*Gt322*gt33L*gtu32 + Gt223*Gt323*gt33L*gtu32 + Gt133*Gt213*gt23L*gtu33 + gt13L*gtu31*SQR(Gt113)) + 
        gt11L*gtu22*SQR(Gt123) + gt11L*gtu33*SQR(Gt133) + gt22L*gtu11*SQR(Gt213) + gt22L*gtu22*SQR(Gt223) + 
        4*gt23L*gtu32*SQR(Gt223) + gt22L*gtu33*SQR(Gt233) + 5*gt33L*gtu11*SQR(Gt313) + 5*gt33L*gtu22*SQR(Gt323) + 
        5*gt33L*gtu33*SQR(Gt333);
    
    Rphi11  =  2*(-PDstandard4th11phi - gt11L*gtu11*PDstandard4th11phi - 2*gt11L*gtu21*PDstandard4th12phi - 
          2*gt11L*gtu31*PDstandard4th13phi - gt11L*gtu22*PDstandard4th22phi - 2*gt11L*gtu32*PDstandard4th23phi - 
          gt11L*gtu33*PDstandard4th33phi + Gt311*PDstandard4th3phi + gt11L*Gt311*gtu11*PDstandard4th3phi + 
          2*gt11L*Gt312*gtu21*PDstandard4th3phi + gt11L*Gt322*gtu22*PDstandard4th3phi + 
          2*gt11L*Gt313*gtu31*PDstandard4th3phi + 2*gt11L*Gt323*gtu32*PDstandard4th3phi + 
          gt11L*Gt333*gtu33*PDstandard4th3phi + PDstandard4th1phi*
           (Gt111 + Gt111*gt11L*gtu11 + 2*Gt112*gt11L*gtu21 + gt11L*Gt122*gtu22 + 2*Gt113*gt11L*gtu31 + 
             2*gt11L*Gt123*gtu32 + gt11L*Gt133*gtu33 - 4*gt11L*gtu21*PDstandard4th2phi - 4*gt11L*gtu31*PDstandard4th3phi) + 
          PDstandard4th2phi*(Gt211 + gt11L*Gt211*gtu11 + 
             gt11L*(2*Gt212*gtu21 + Gt222*gtu22 + 2*Gt213*gtu31 + 2*Gt223*gtu32 + Gt233*gtu33) - 
             4*gt11L*gtu32*PDstandard4th3phi) + (2 - 2*gt11L*gtu11)*SQR(PDstandard4th1phi) - 
          2*gt11L*gtu22*SQR(PDstandard4th2phi) - 2*gt11L*gtu33*SQR(PDstandard4th3phi));
    
    Rphi12  =  2*(-(gt12L*gtu11*PDstandard4th11phi) - PDstandard4th12phi - 2*gt12L*gtu21*PDstandard4th12phi - 
          2*gt12L*gtu31*PDstandard4th13phi - gt12L*gtu22*PDstandard4th22phi - 2*gt12L*gtu32*PDstandard4th23phi - 
          gt12L*gtu33*PDstandard4th33phi + Gt312*PDstandard4th3phi + gt12L*Gt311*gtu11*PDstandard4th3phi + 
          2*gt12L*Gt312*gtu21*PDstandard4th3phi + gt12L*Gt322*gtu22*PDstandard4th3phi + 
          2*gt12L*Gt313*gtu31*PDstandard4th3phi + 2*gt12L*Gt323*gtu32*PDstandard4th3phi + 
          gt12L*Gt333*gtu33*PDstandard4th3phi + PDstandard4th1phi*
           (Gt112 + Gt111*gt12L*gtu11 + 2*Gt112*gt12L*gtu21 + Gt122*gt12L*gtu22 + 2*Gt113*gt12L*gtu31 + 
             2*Gt123*gt12L*gtu32 + gt12L*Gt133*gtu33 + (2 - 4*gt12L*gtu21)*PDstandard4th2phi - 
             4*gt12L*gtu31*PDstandard4th3phi) + PDstandard4th2phi*
           (Gt212 + 2*gt12L*Gt212*gtu21 + gt12L*(Gt211*gtu11 + Gt222*gtu22 + 2*Gt213*gtu31 + 2*Gt223*gtu32 + Gt233*gtu33) - 
             4*gt12L*gtu32*PDstandard4th3phi) - 2*gt12L*gtu11*SQR(PDstandard4th1phi) - 
          2*gt12L*gtu22*SQR(PDstandard4th2phi) - 2*gt12L*gtu33*SQR(PDstandard4th3phi));
    
    Rphi13  =  2*(-PDstandard4th13phi + gt13L*(-(gtu11*PDstandard4th11phi) - 2*gtu21*PDstandard4th12phi - 
             2*gtu31*PDstandard4th13phi) - gt13L*gtu22*PDstandard4th22phi - 2*gt13L*gtu32*PDstandard4th23phi - 
          gt13L*gtu33*PDstandard4th33phi + Gt313*PDstandard4th3phi + gt13L*Gt311*gtu11*PDstandard4th3phi + 
          2*gt13L*Gt312*gtu21*PDstandard4th3phi + gt13L*Gt322*gtu22*PDstandard4th3phi + 
          2*gt13L*Gt313*gtu31*PDstandard4th3phi + 2*gt13L*Gt323*gtu32*PDstandard4th3phi + 
          gt13L*Gt333*gtu33*PDstandard4th3phi + PDstandard4th1phi*
           (Gt113 + Gt111*gt13L*gtu11 + 2*Gt112*gt13L*gtu21 + Gt122*gt13L*gtu22 + 2*Gt113*gt13L*gtu31 + 
             2*Gt123*gt13L*gtu32 + Gt133*gt13L*gtu33 - 4*gt13L*gtu21*PDstandard4th2phi + 
             (2 - 4*gt13L*gtu31)*PDstandard4th3phi) + 
          PDstandard4th2phi*(Gt213 + 2*gt13L*Gt213*gtu31 + 
             gt13L*(Gt211*gtu11 + 2*Gt212*gtu21 + Gt222*gtu22 + 2*Gt223*gtu32 + Gt233*gtu33) - 
             4*gt13L*gtu32*PDstandard4th3phi) - 2*gt13L*gtu11*SQR(PDstandard4th1phi) - 
          2*gt13L*gtu22*SQR(PDstandard4th2phi) - 2*gt13L*gtu33*SQR(PDstandard4th3phi));
    
    Rphi22  =  2*(-PDstandard4th22phi + gt22L*(-(gtu11*PDstandard4th11phi) - 2*gtu21*PDstandard4th12phi - 
             2*gtu31*PDstandard4th13phi - gtu22*PDstandard4th22phi) - 2*gt22L*gtu32*PDstandard4th23phi - 
          gt22L*gtu33*PDstandard4th33phi + Gt322*PDstandard4th3phi + gt22L*Gt311*gtu11*PDstandard4th3phi + 
          2*gt22L*Gt312*gtu21*PDstandard4th3phi + gt22L*Gt322*gtu22*PDstandard4th3phi + 
          2*gt22L*Gt313*gtu31*PDstandard4th3phi + 2*gt22L*Gt323*gtu32*PDstandard4th3phi + 
          gt22L*Gt333*gtu33*PDstandard4th3phi + PDstandard4th1phi*
           (Gt122 + Gt111*gt22L*gtu11 + 2*Gt112*gt22L*gtu21 + Gt122*gt22L*gtu22 + 2*Gt113*gt22L*gtu31 + 
             2*Gt123*gt22L*gtu32 + Gt133*gt22L*gtu33 - 4*gt22L*gtu21*PDstandard4th2phi - 4*gt22L*gtu31*PDstandard4th3phi) + 
          PDstandard4th2phi*(Gt222 + Gt222*gt22L*gtu22 + 
             gt22L*(Gt211*gtu11 + 2*Gt212*gtu21 + 2*Gt213*gtu31 + 2*Gt223*gtu32 + Gt233*gtu33) - 
             4*gt22L*gtu32*PDstandard4th3phi) - 2*gt22L*gtu11*SQR(PDstandard4th1phi) + 
          (2 - 2*gt22L*gtu22)*SQR(PDstandard4th2phi) - 2*gt22L*gtu33*SQR(PDstandard4th3phi));
    
    Rphi23  =  2*(-PDstandard4th23phi + gt23L*(-(gtu11*PDstandard4th11phi) - 2*gtu21*PDstandard4th12phi - 
             2*gtu31*PDstandard4th13phi - gtu22*PDstandard4th22phi - 2*gtu32*PDstandard4th23phi) - 
          gt23L*gtu33*PDstandard4th33phi + Gt323*PDstandard4th3phi + gt23L*Gt311*gtu11*PDstandard4th3phi + 
          2*gt23L*Gt312*gtu21*PDstandard4th3phi + gt23L*Gt322*gtu22*PDstandard4th3phi + 
          2*gt23L*Gt313*gtu31*PDstandard4th3phi + 2*gt23L*Gt323*gtu32*PDstandard4th3phi + 
          gt23L*Gt333*gtu33*PDstandard4th3phi + PDstandard4th1phi*
           (Gt123 + Gt111*gt23L*gtu11 + 2*Gt112*gt23L*gtu21 + Gt122*gt23L*gtu22 + 2*Gt113*gt23L*gtu31 + 
             2*Gt123*gt23L*gtu32 + Gt133*gt23L*gtu33 - 4*gt23L*gtu21*PDstandard4th2phi - 4*gt23L*gtu31*PDstandard4th3phi) + 
          PDstandard4th2phi*(Gt223 + 2*Gt223*gt23L*gtu32 + 
             gt23L*(Gt211*gtu11 + 2*Gt212*gtu21 + Gt222*gtu22 + 2*Gt213*gtu31 + Gt233*gtu33) + 
             (2 - 4*gt23L*gtu32)*PDstandard4th3phi) - 2*gt23L*gtu11*SQR(PDstandard4th1phi) - 
          2*gt23L*gtu22*SQR(PDstandard4th2phi) - 2*gt23L*gtu33*SQR(PDstandard4th3phi));
    
    Rphi33  =  2*(-PDstandard4th33phi + (Gt333 + gt33L*
              (Gt322*gtu22 + 2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33))*PDstandard4th3phi + 
          PDstandard4th2phi*(Gt233 + gt33L*(Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + 
                Gt233*gtu33 - 4*gtu32*PDstandard4th3phi)) + 
          PDstandard4th1phi*(Gt133 + gt33L*(Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + 
                Gt133*gtu33 - 4*(gtu21*PDstandard4th2phi + gtu31*PDstandard4th3phi))) + 2*SQR(PDstandard4th3phi) + 
          gt33L*(-(gtu11*PDstandard4th11phi) - 2*gtu21*PDstandard4th12phi - 2*gtu31*PDstandard4th13phi - 
             gtu22*PDstandard4th22phi - 2*gtu32*PDstandard4th23phi - gtu33*PDstandard4th33phi + 
             Gt311*gtu11*PDstandard4th3phi - 2*gtu11*SQR(PDstandard4th1phi) - 2*gtu22*SQR(PDstandard4th2phi) - 
             2*gtu33*SQR(PDstandard4th3phi)));
    
    Atm11  =  At11L*gtu11 + At12L*gtu21 + At13L*gtu31;
    
    Atm21  =  At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    Atm31  =  At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    Atm12  =  At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    Atm22  =  At12L*gtu21 + At22L*gtu22 + At23L*gtu32;
    
    Atm32  =  At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    Atm13  =  At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    Atm23  =  At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    Atm33  =  At13L*gtu31 + At23L*gtu32 + At33L*gtu33;
    
    Atu11  =  Atm11*gtu11 + Atm12*gtu21 + Atm13*gtu31;
    
    Atu21  =  Atm11*gtu21 + Atm12*gtu22 + Atm13*gtu32;
    
    Atu31  =  Atm11*gtu31 + Atm12*gtu32 + Atm13*gtu33;
    
    Atu22  =  Atm21*gtu21 + Atm22*gtu22 + Atm23*gtu32;
    
    Atu32  =  Atm21*gtu31 + Atm22*gtu32 + Atm23*gtu33;
    
    Atu33  =  Atm31*gtu31 + Atm32*gtu32 + Atm33*gtu33;
    
    e4phi  =  exp(4*phiL);
    
    em4phi  =  INV(e4phi);
    
    g11  =  e4phi*gt11L;
    
    g12  =  e4phi*gt12L;
    
    g13  =  e4phi*gt13L;
    
    g22  =  e4phi*gt22L;
    
    g23  =  e4phi*gt23L;
    
    g33  =  e4phi*gt33L;
    
    detg  =  2*g12*g13*g23 + g33*(g11*g22 - SQR(g12)) - g22*SQR(g13) - g11*SQR(g23);
    
    gu11  =  em4phi*gtu11;
    
    gu21  =  em4phi*gtu21;
    
    gu31  =  em4phi*gtu31;
    
    gu22  =  em4phi*gtu22;
    
    gu32  =  em4phi*gtu32;
    
    gu33  =  em4phi*gtu33;
    
    ddetg1  =  e4phi*(ddetgt1 + 4*detgt*PDstandard4th1phi);
    
    ddetg2  =  e4phi*(ddetgt2 + 4*detgt*PDstandard4th2phi);
    
    ddetg3  =  e4phi*(ddetgt3 + 4*detgt*PDstandard4th3phi);
    
    G111  =  -((-6*detg*Gt111 + ddetg1*(-6 + g11*gu11) + g11*(ddetg2*gu21 + ddetg3*gu31))*INV(detg))/6.;
    
    G211  =  ((6*detg*Gt211 - g11*(ddetg1*gu21 + ddetg2*gu22 + ddetg3*gu32))*INV(detg))/6.;
    
    G311  =  ((6*detg*Gt311 - g11*(ddetg1*gu31 + ddetg2*gu32 + ddetg3*gu33))*INV(detg))/6.;
    
    G112  =  -((-6*detg*Gt112 + ddetg2*(-3 + g12*gu21) + g12*(ddetg1*gu11 + ddetg3*gu31))*INV(detg))/6.;
    
    G212  =  -((-6*detg*Gt212 + ddetg1*(-3 + g12*gu21) + g12*(ddetg2*gu22 + ddetg3*gu32))*INV(detg))/6.;
    
    G312  =  ((6*detg*Gt312 - g12*(ddetg1*gu31 + ddetg2*gu32 + ddetg3*gu33))*INV(detg))/6.;
    
    G113  =  -((-6*detg*Gt113 + g13*(ddetg1*gu11 + ddetg2*gu21) + ddetg3*(-3 + g13*gu31))*INV(detg))/6.;
    
    G213  =  ((6*detg*Gt213 - g13*(ddetg1*gu21 + ddetg2*gu22 + ddetg3*gu32))*INV(detg))/6.;
    
    G313  =  -((-6*detg*Gt313 + ddetg1*(-3 + g13*gu31) + g13*(ddetg2*gu32 + ddetg3*gu33))*INV(detg))/6.;
    
    G122  =  ((6*detg*Gt122 - g22*(ddetg1*gu11 + ddetg2*gu21 + ddetg3*gu31))*INV(detg))/6.;
    
    G222  =  -((-6*detg*Gt222 + ddetg2*(-6 + g22*gu22) + g22*(ddetg1*gu21 + ddetg3*gu32))*INV(detg))/6.;
    
    G322  =  ((6*detg*Gt322 - g22*(ddetg1*gu31 + ddetg2*gu32 + ddetg3*gu33))*INV(detg))/6.;
    
    G123  =  ((6*detg*Gt123 - g23*(ddetg1*gu11 + ddetg2*gu21 + ddetg3*gu31))*INV(detg))/6.;
    
    G223  =  -((-6*detg*Gt223 + g23*(ddetg1*gu21 + ddetg2*gu22) + ddetg3*(-3 + g23*gu32))*INV(detg))/6.;
    
    G323  =  -((-6*detg*Gt323 + ddetg2*(-3 + g23*gu32) + g23*(ddetg1*gu31 + ddetg3*gu33))*INV(detg))/6.;
    
    G133  =  ((6*detg*Gt133 - g33*(ddetg1*gu11 + ddetg2*gu21 + ddetg3*gu31))*INV(detg))/6.;
    
    G233  =  ((6*detg*Gt233 - g33*(ddetg1*gu21 + ddetg2*gu22 + ddetg3*gu32))*INV(detg))/6.;
    
    G333  =  -((-6*detg*Gt333 + g33*(ddetg1*gu31 + ddetg2*gu32) + ddetg3*(-6 + g33*gu33))*INV(detg))/6.;
    
    R11  =  Rphi11 + Rt11;
    
    R12  =  Rphi12 + Rt12;
    
    R13  =  Rphi13 + Rt13;
    
    R22  =  Rphi22 + Rt22;
    
    R23  =  Rphi23 + Rt23;
    
    R33  =  Rphi33 + Rt33;
    
    phirhsL  =  beta1L*PDstandard4th1phi + beta2L*PDstandard4th2phi + beta3L*PDstandard4th3phi + 
        ((PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3)*phiL)/6. - (alphaL*trKL)/6.;
    
    gt11rhsL  =  -2*alphaL*At11L + 2*(gt11L*PDstandard4th1beta1 + gt12L*PDstandard4th1beta2 + gt13L*PDstandard4th1beta3) + 
        beta1L*PDstandard4th1gt11 + beta2L*PDstandard4th2gt11 - 
        gt11L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) + beta3L*PDstandard4th3gt11;
    
    gt12rhsL  =  -2*alphaL*At12L + gt22L*PDstandard4th1beta2 + gt23L*PDstandard4th1beta3 + beta1L*PDstandard4th1gt12 + 
        gt11L*PDstandard4th2beta1 + gt13L*PDstandard4th2beta3 + beta2L*PDstandard4th2gt12 + 
        gt12L*(kthird*(PDstandard4th1beta1 + PDstandard4th2beta2) - ktwothird*PDstandard4th3beta3) + 
        beta3L*PDstandard4th3gt12;
    
    gt13rhsL  =  -2*alphaL*At13L + gt23L*PDstandard4th1beta2 + gt33L*PDstandard4th1beta3 + beta1L*PDstandard4th1gt13 + 
        beta2L*PDstandard4th2gt13 + gt11L*PDstandard4th3beta1 + gt12L*PDstandard4th3beta2 + 
        gt13L*(-(ktwothird*PDstandard4th2beta2) + kthird*(PDstandard4th1beta1 + PDstandard4th3beta3)) + 
        beta3L*PDstandard4th3gt13;
    
    gt22rhsL  =  -2*alphaL*At22L + beta1L*PDstandard4th1gt22 + 
        2*(gt12L*PDstandard4th2beta1 + gt22L*PDstandard4th2beta2 + gt23L*PDstandard4th2beta3) + beta2L*PDstandard4th2gt22 - 
        gt22L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) + beta3L*PDstandard4th3gt22;
    
    gt23rhsL  =  -2*alphaL*At23L + beta1L*PDstandard4th1gt23 + gt13L*PDstandard4th2beta1 + gt33L*PDstandard4th2beta3 + 
        beta2L*PDstandard4th2gt23 + gt12L*PDstandard4th3beta1 + gt22L*PDstandard4th3beta2 + 
        gt23L*(-(ktwothird*PDstandard4th1beta1) + kthird*(PDstandard4th2beta2 + PDstandard4th3beta3)) + 
        beta3L*PDstandard4th3gt23;
    
    gt33rhsL  =  -2*alphaL*At33L + beta1L*PDstandard4th1gt33 + beta2L*PDstandard4th2gt33 - 
        gt33L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) + 
        2*(gt13L*PDstandard4th3beta1 + gt23L*PDstandard4th3beta2 + gt33L*PDstandard4th3beta3) + beta3L*PDstandard4th3gt33;
    
    Xt1rhsL  =  kthird*(-3*(beta1L*(ddgtu1111 + ddgtu2112 + ddgtu3113) + beta2L*(ddgtu1112 + ddgtu2122 + ddgtu3123) + 
             beta3L*(ddgtu1113 + ddgtu2123 + ddgtu3133)) + 16*(gtu21*PDstandard4th12beta1 + gtu31*PDstandard4th13beta1) + 
          (7*dgtu111 + 4*(dgtu212 + dgtu313))*PDstandard4th1beta1 + 
          gtu11*(10*PDstandard4th11beta1 + 4*(PDstandard4th12beta2 + PDstandard4th13beta3) - 4*alphaL*PDstandard4th1trK) + 
          12*gtu32*PDstandard4th23beta1 + (9*dgtu211 + 6*dgtu323)*PDstandard4th2beta1 + 
          (dgtu212 - 2*(dgtu111 + dgtu313))*PDstandard4th2beta2 + 
          4*(gtu21*(PDstandard4th22beta2 + PDstandard4th23beta3) + gtu31*(PDstandard4th23beta2 + PDstandard4th33beta3)) - 
          6*(Atu11*PDstandard4th1alpha + Atu21*PDstandard4th2alpha + Atu31*PDstandard4th3alpha) + 
          (9*dgtu311 + 6*dgtu333)*PDstandard4th3beta1 + 
          6*(gtu22*PDstandard4th22beta1 + dgtu222*PDstandard4th2beta1 + gtu33*PDstandard4th33beta1 + 
             dgtu322*PDstandard4th3beta1) + 3*(dgtu112*PDstandard4th1beta2 + dgtu113*PDstandard4th1beta3 + 
             dgtu213*PDstandard4th2beta3 + dgtu312*PDstandard4th3beta2) + 
          (-2*(dgtu111 + dgtu212) + dgtu313)*PDstandard4th3beta3 + 
          alphaL*(12*(Atu21*Gt112 + Atu31*Gt113 + Atu32*Gt123) + 6*(Atu11*Gt111 + Atu22*Gt122 + Atu33*Gt133) + 
             36*(Atu11*PDstandard4th1phi + Atu21*PDstandard4th2phi + Atu31*PDstandard4th3phi) - 
             4*(gtu21*PDstandard4th2trK + gtu31*PDstandard4th3trK)));
    
    Xt2rhsL  =  kthird*(-3*(beta1L*(ddgtu2111 + ddgtu2212 + ddgtu3213) + beta2L*(ddgtu2112 + ddgtu2222 + ddgtu3223) + 
             beta3L*(ddgtu2113 + ddgtu2223 + ddgtu3233)) + 12*gtu31*PDstandard4th13beta2 + 
          (dgtu211 - 2*(dgtu222 + dgtu323))*PDstandard4th1beta1 + (9*dgtu212 + 6*dgtu313)*PDstandard4th1beta2 + 
          gtu21*(16*PDstandard4th12beta2 + 4*PDstandard4th13beta3 - 4*alphaL*PDstandard4th1trK) + 
          (7*dgtu222 + 4*dgtu323)*PDstandard4th2beta2 + 
          4*(gtu21*PDstandard4th11beta1 + gtu32*PDstandard4th13beta1 + 
             gtu22*(PDstandard4th12beta1 + PDstandard4th23beta3) + dgtu211*PDstandard4th2beta2) + 
          gtu22*(10*PDstandard4th22beta2 - 4*alphaL*PDstandard4th2trK) - 
          6*(Atu21*PDstandard4th1alpha + Atu22*PDstandard4th2alpha + Atu32*PDstandard4th3alpha) + 
          3*(dgtu213*PDstandard4th1beta3 + dgtu221*PDstandard4th2beta1 + dgtu223*PDstandard4th2beta3 + 
             dgtu321*PDstandard4th3beta1) + (9*dgtu322 + 6*dgtu333)*PDstandard4th3beta2 + 
          6*(gtu11*PDstandard4th11beta2 + dgtu111*PDstandard4th1beta2 + gtu33*PDstandard4th33beta2 + 
             dgtu311*PDstandard4th3beta2) + (-2*(dgtu211 + dgtu222) + dgtu323)*PDstandard4th3beta3 + 
          alphaL*(12*(Atu21*Gt212 + Atu31*Gt213 + Atu32*Gt223) + 6*(Atu11*Gt211 + Atu22*Gt222 + Atu33*Gt233) + 
             36*(Atu21*PDstandard4th1phi + Atu22*PDstandard4th2phi + Atu32*PDstandard4th3phi)) + 
          gtu32*(16*PDstandard4th23beta2 + 4*PDstandard4th33beta3 - 4*alphaL*PDstandard4th3trK));
    
    Xt3rhsL  =  kthird*(-3*(beta1L*(ddgtu3111 + ddgtu3212 + ddgtu3313) + beta2L*(ddgtu3112 + ddgtu3222 + ddgtu3323) + 
             beta3L*(ddgtu3113 + ddgtu3223 + ddgtu3333)) + 12*gtu21*PDstandard4th12beta3 + 
          (dgtu311 - 2*(dgtu322 + dgtu333))*PDstandard4th1beta1 + (6*dgtu212 + 9*dgtu313)*PDstandard4th1beta3 + 
          gtu31*(16*PDstandard4th13beta3 - 4*alphaL*PDstandard4th1trK) + 
          (dgtu322 - 2*(dgtu311 + dgtu333))*PDstandard4th2beta2 + (6*dgtu222 + 9*dgtu323)*PDstandard4th2beta3 + 
          6*(gtu11*PDstandard4th11beta3 + dgtu111*PDstandard4th1beta3 + gtu22*PDstandard4th22beta3 + 
             dgtu211*PDstandard4th2beta3) + gtu32*(16*PDstandard4th23beta3 - 4*alphaL*PDstandard4th2trK) - 
          6*(Atu31*PDstandard4th1alpha + Atu32*PDstandard4th2alpha + Atu33*PDstandard4th3alpha) + 
          3*(dgtu312*PDstandard4th1beta2 + dgtu321*PDstandard4th2beta1 + dgtu331*PDstandard4th3beta1 + 
             dgtu332*PDstandard4th3beta2) + (4*dgtu322 + 7*dgtu333)*PDstandard4th3beta3 + 
          4*(gtu31*(PDstandard4th11beta1 + PDstandard4th12beta2) + gtu32*(PDstandard4th12beta1 + PDstandard4th22beta2) + 
             gtu33*(PDstandard4th13beta1 + PDstandard4th23beta2) + dgtu311*PDstandard4th3beta3) + 
          alphaL*(12*(Atu21*Gt312 + Atu31*Gt313 + Atu32*Gt323) + 6*(Atu11*Gt311 + Atu22*Gt322 + Atu33*Gt333) + 
             36*(Atu31*PDstandard4th1phi + Atu32*PDstandard4th2phi + Atu33*PDstandard4th3phi)) + 
          gtu33*(10*PDstandard4th33beta3 - 4*alphaL*PDstandard4th3trK));
    
    trKrhsL  =  -(gu11*PDstandard4th11alpha) - 2*gu21*PDstandard4th12alpha - 2*gu31*PDstandard4th13alpha + 
        (G111*gu11 + 2*G112*gu21 + G122*gu22 + 2*G113*gu31 + 2*G123*gu32 + G133*gu33)*PDstandard4th1alpha + 
        beta1L*PDstandard4th1trK - gu22*PDstandard4th22alpha - 2*gu32*PDstandard4th23alpha + 
        (G211*gu11 + 2*G212*gu21 + G222*gu22 + 2*G213*gu31 + 2*G223*gu32 + G233*gu33)*PDstandard4th2alpha + 
        beta2L*PDstandard4th2trK - gu33*PDstandard4th33alpha + G311*gu11*PDstandard4th3alpha + 
        G322*gu22*PDstandard4th3alpha + 2*G313*gu31*PDstandard4th3alpha + 2*G323*gu32*PDstandard4th3alpha + 
        G333*gu33*PDstandard4th3alpha + 2*(alphaL*Atm12*Atm21 + alphaL*Atm13*Atm31 + alphaL*Atm23*Atm32 + 
           G312*gu21*PDstandard4th3alpha) + beta3L*PDstandard4th3trK + alphaL*SQR(Atm11) + alphaL*SQR(Atm22) + 
        alphaL*SQR(Atm33) + alphaL*kthird*SQR(trKL);
    
    At11rhsL  =  beta1L*PDstandard4th1At11 + 2*(At11L*PDstandard4th1beta1 + At12L*PDstandard4th1beta2 + 
           At13L*PDstandard4th1beta3) + beta2L*PDstandard4th2At11 + beta3L*PDstandard4th3At11 - 
        At11L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) - 
        em4phi*kthird*(3*PDstandard4th11alpha - g11*gu11*PDstandard4th11alpha - 2*g11*gu21*PDstandard4th12alpha - 
           2*g11*gu31*PDstandard4th13alpha + (G111*(-3 + g11*gu11) + 
              g11*(2*G112*gu21 + G122*gu22 + 2*G113*gu31 + 2*G123*gu32 + G133*gu33))*PDstandard4th1alpha - 
           g11*gu22*PDstandard4th22alpha - 2*g11*gu32*PDstandard4th23alpha + 
           (G211*(-3 + g11*gu11) + g11*(2*G212*gu21 + G222*gu22 + 2*G213*gu31 + 2*G223*gu32 + G233*gu33))*
            PDstandard4th2alpha - g11*gu33*PDstandard4th33alpha - 3*G311*PDstandard4th3alpha + 
           g11*G311*gu11*PDstandard4th3alpha + 2*g11*G312*gu21*PDstandard4th3alpha + g11*G322*gu22*PDstandard4th3alpha + 
           2*g11*G313*gu31*PDstandard4th3alpha + 2*g11*G323*gu32*PDstandard4th3alpha + g11*G333*gu33*PDstandard4th3alpha - 
           3*alphaL*R11 + alphaL*g11*gu11*R11 + 2*alphaL*g11*gu21*R12 + 2*alphaL*g11*gu31*R13 + alphaL*g11*gu22*R22 + 
           2*alphaL*g11*gu32*R23 + alphaL*g11*gu33*R33) + alphaL*(-2*(At11L*Atm11 + At12L*Atm21 + At13L*Atm31) + At11L*trKL);
    
    At12rhsL  =  beta1L*PDstandard4th1At12 + At22L*PDstandard4th1beta2 + At23L*PDstandard4th1beta3 + 
        beta2L*PDstandard4th2At12 + At11L*PDstandard4th2beta1 + At13L*PDstandard4th2beta3 + beta3L*PDstandard4th3At12 + 
        At12L*(PDstandard4th1beta1 + PDstandard4th2beta2 - 
           ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3)) - 
        em4phi*kthird*(3*PDstandard4th12alpha + (G112*(-3 + 2*g12*gu21) + 
              g12*(G111*gu11 + G122*gu22 + 2*(G113*gu31 + G123*gu32) + G133*gu33))*PDstandard4th1alpha + 
           (G212*(-3 + 2*g12*gu21) + g12*(G211*gu11 + G222*gu22 + 2*(G213*gu31 + G223*gu32) + G233*gu33))*
            PDstandard4th2alpha + G312*(-3 + 2*g12*gu21)*PDstandard4th3alpha + 
           alphaL*(-3*R12 + 2*g12*(gu31*R13 + gu32*R23)) + 
           g12*(-(gu11*PDstandard4th11alpha) - 2*gu21*PDstandard4th12alpha - 2*gu31*PDstandard4th13alpha - 
              gu22*PDstandard4th22alpha - 2*gu32*PDstandard4th23alpha - gu33*PDstandard4th33alpha + 
              G311*gu11*PDstandard4th3alpha + G322*gu22*PDstandard4th3alpha + 2*G313*gu31*PDstandard4th3alpha + 
              2*G323*gu32*PDstandard4th3alpha + G333*gu33*PDstandard4th3alpha + alphaL*gu11*R11 + 2*alphaL*gu21*R12 + 
              alphaL*gu22*R22 + alphaL*gu33*R33)) + alphaL*(-2*(At11L*Atm12 + At12L*Atm22 + At13L*Atm32) + At12L*trKL);
    
    At13rhsL  =  beta1L*PDstandard4th1At13 + At23L*PDstandard4th1beta2 + At33L*PDstandard4th1beta3 + 
        beta2L*PDstandard4th2At13 + beta3L*PDstandard4th3At13 + At11L*PDstandard4th3beta1 + At12L*PDstandard4th3beta2 + 
        At13L*(PDstandard4th1beta1 + PDstandard4th3beta3 - 
           ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3)) - 
        em4phi*kthird*(3*PDstandard4th13alpha + (G113*(-3 + 2*g13*gu31) + 
              g13*(G111*gu11 + G122*gu22 + 2*(G112*gu21 + G123*gu32) + G133*gu33))*PDstandard4th1alpha + 
           (G213*(-3 + 2*g13*gu31) + g13*(G211*gu11 + G222*gu22 + 2*(G212*gu21 + G223*gu32) + G233*gu33))*
            PDstandard4th2alpha + (G313*(-3 + 2*g13*gu31) + g13*(G322*gu22 + 2*(G312*gu21 + G323*gu32) + G333*gu33))*
            PDstandard4th3alpha + g13*(-(gu11*PDstandard4th11alpha) - 2*gu21*PDstandard4th12alpha - 
              2*gu31*PDstandard4th13alpha - gu22*PDstandard4th22alpha - 2*gu32*PDstandard4th23alpha - 
              gu33*PDstandard4th33alpha + G311*gu11*PDstandard4th3alpha + alphaL*gu11*R11 + 2*alphaL*gu21*R12 + 
              2*alphaL*gu31*R13 + 2*alphaL*gu32*R23) + alphaL*(-3*R13 + g13*(gu22*R22 + gu33*R33))) + 
        alphaL*(-2*(At11L*Atm13 + At12L*Atm23 + At13L*Atm33) + At13L*trKL);
    
    At22rhsL  =  beta1L*PDstandard4th1At22 + beta2L*PDstandard4th2At22 + 
        2*(At12L*PDstandard4th2beta1 + At22L*PDstandard4th2beta2 + At23L*PDstandard4th2beta3) + beta3L*PDstandard4th3At22 - 
        At22L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) - 
        em4phi*kthird*(-(g22*gu11*PDstandard4th11alpha) - 2*g22*gu21*PDstandard4th12alpha - 
           2*g22*gu31*PDstandard4th13alpha + (G122*(-3 + g22*gu22) + 
              g22*(G111*gu11 + 2*G112*gu21 + 2*G113*gu31 + 2*G123*gu32 + G133*gu33))*PDstandard4th1alpha + 
           3*PDstandard4th22alpha - g22*gu22*PDstandard4th22alpha - 2*g22*gu32*PDstandard4th23alpha + 
           (G222*(-3 + g22*gu22) + g22*(G211*gu11 + 2*G212*gu21 + 2*G213*gu31 + 2*G223*gu32 + G233*gu33))*
            PDstandard4th2alpha - g22*gu33*PDstandard4th33alpha - 3*G322*PDstandard4th3alpha + 
           g22*G311*gu11*PDstandard4th3alpha + 2*g22*G312*gu21*PDstandard4th3alpha + g22*G322*gu22*PDstandard4th3alpha + 
           2*g22*G313*gu31*PDstandard4th3alpha + 2*g22*G323*gu32*PDstandard4th3alpha + g22*G333*gu33*PDstandard4th3alpha + 
           alphaL*g22*gu11*R11 + 2*alphaL*g22*gu21*R12 + 2*alphaL*g22*gu31*R13 - 3*alphaL*R22 + alphaL*g22*gu22*R22 + 
           2*alphaL*g22*gu32*R23 + alphaL*g22*gu33*R33) + alphaL*(-2*(At12L*Atm12 + At22L*Atm22 + At23L*Atm32) + At22L*trKL);
    
    At23rhsL  =  beta1L*PDstandard4th1At23 + beta2L*PDstandard4th2At23 + At13L*PDstandard4th2beta1 + 
        At33L*PDstandard4th2beta3 + beta3L*PDstandard4th3At23 + At12L*PDstandard4th3beta1 + At22L*PDstandard4th3beta2 + 
        At23L*(PDstandard4th2beta2 + PDstandard4th3beta3 - 
           ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3)) - 
        em4phi*kthird*((G123*(-3 + 2*g23*gu32) + g23*(G111*gu11 + G122*gu22 + 2*(G112*gu21 + G113*gu31) + G133*gu33))*
            PDstandard4th1alpha + 3*PDstandard4th23alpha + 
           (G223*(-3 + 2*g23*gu32) + g23*(G211*gu11 + G222*gu22 + 2*(G212*gu21 + G213*gu31) + G233*gu33))*
            PDstandard4th2alpha + (G323*(-3 + 2*g23*gu32) + g23*(G322*gu22 + 2*(G312*gu21 + G313*gu31) + G333*gu33))*
            PDstandard4th3alpha + g23*(-(gu11*PDstandard4th11alpha) - 2*gu21*PDstandard4th12alpha - 
              2*gu31*PDstandard4th13alpha - gu22*PDstandard4th22alpha - 2*gu32*PDstandard4th23alpha - 
              gu33*PDstandard4th33alpha + G311*gu11*PDstandard4th3alpha + alphaL*gu11*R11 + 2*alphaL*gu21*R12 + 
              2*alphaL*gu31*R13 + alphaL*gu22*R22 + 2*alphaL*gu32*R23) + alphaL*(-3*R23 + g23*gu33*R33)) + 
        alphaL*(-2*(At12L*Atm13 + At22L*Atm23 + At23L*Atm33) + At23L*trKL);
    
    At33rhsL  =  beta1L*PDstandard4th1At33 + beta2L*PDstandard4th2At33 + beta3L*PDstandard4th3At33 - 
        At33L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) + 
        2*(At13L*PDstandard4th3beta1 + At23L*PDstandard4th3beta2 + At33L*PDstandard4th3beta3) - 
        em4phi*kthird*(-(g33*gu11*PDstandard4th11alpha) - 2*g33*gu21*PDstandard4th12alpha - 
           2*g33*gu31*PDstandard4th13alpha + (g33*(G111*gu11 + 2*G112*gu21 + G122*gu22 + 2*G113*gu31 + 2*G123*gu32) + 
              G133*(-3 + g33*gu33))*PDstandard4th1alpha - g33*gu22*PDstandard4th22alpha - 2*g33*gu32*PDstandard4th23alpha + 
           (g33*(G211*gu11 + 2*G212*gu21 + G222*gu22 + 2*G213*gu31 + 2*G223*gu32) + G233*(-3 + g33*gu33))*
            PDstandard4th2alpha + 3*PDstandard4th33alpha - g33*gu33*PDstandard4th33alpha - 3*G333*PDstandard4th3alpha + 
           G311*g33*gu11*PDstandard4th3alpha + 2*G312*g33*gu21*PDstandard4th3alpha + G322*g33*gu22*PDstandard4th3alpha + 
           2*G313*g33*gu31*PDstandard4th3alpha + 2*G323*g33*gu32*PDstandard4th3alpha + g33*G333*gu33*PDstandard4th3alpha + 
           alphaL*g33*gu11*R11 + 2*alphaL*g33*gu21*R12 + 2*alphaL*g33*gu31*R13 + alphaL*g33*gu22*R22 + 
           2*alphaL*g33*gu32*R23 - 3*alphaL*R33 + alphaL*g33*gu33*R33) + 
        alphaL*(-2*(At13L*Atm13 + At23L*Atm23 + At33L*Atm33) + At33L*trKL);
    
    alpharhsL  =  beta1L*PDstandard4th1alpha + beta2L*PDstandard4th2alpha + beta3L*PDstandard4th3alpha - 
        dtalphaL*harmonicF*pow(alphaL,harmonicN);
    
    dtalpharhsL  =  -(AlphaDriver*dtalphaL) + trKrhsL;
    
    beta1rhsL  =  dtbeta1L*ShiftGammaCoeff*pow(alphaL,ShiftAlphaPower);
    
    beta2rhsL  =  dtbeta2L*ShiftGammaCoeff*pow(alphaL,ShiftAlphaPower);
    
    beta3rhsL  =  dtbeta3L*ShiftGammaCoeff*pow(alphaL,ShiftAlphaPower);
    
    dtbeta1rhsL  =  -(BetaDriver*dtbeta1L) + Xt1rhsL;
    
    dtbeta2rhsL  =  -(BetaDriver*dtbeta2L) + Xt2rhsL;
    
    dtbeta3rhsL  =  -(BetaDriver*dtbeta3L) + Xt3rhsL;
    
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    At11rhs[index] = At11rhsL;
    At12rhs[index] = At12rhsL;
    At13rhs[index] = At13rhsL;
    At22rhs[index] = At22rhsL;
    At23rhs[index] = At23rhsL;
    At33rhs[index] = At33rhsL;
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    dtalpharhs[index] = dtalpharhsL;
    dtbeta1rhs[index] = dtbeta1rhsL;
    dtbeta2rhs[index] = dtbeta2rhsL;
    dtbeta3rhs[index] = dtbeta3rhsL;
    gt11rhs[index] = gt11rhsL;
    gt12rhs[index] = gt12rhsL;
    gt13rhs[index] = gt13rhsL;
    gt22rhs[index] = gt22rhsL;
    gt23rhs[index] = gt23rhsL;
    gt33rhs[index] = gt33rhsL;
    phirhs[index] = phirhsL;
    trKrhs[index] = trKrhsL;
    Xt1rhs[index] = Xt1rhsL;
    Xt2rhs[index] = Xt2rhsL;
    Xt3rhs[index] = Xt3rhsL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (somename);
}

void ML_BSSN_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_RHS_Body);
}
