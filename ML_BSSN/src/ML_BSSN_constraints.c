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

void ML_BSSN_constraints_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_constraints_Body");
  }
  
  if (cctk_iteration % ML_BSSN_constraints_calc_every != ML_BSSN_constraints_calc_offset)
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
    CCTK_REAL ddetg1 = INITVALUE, ddetg2 = INITVALUE, ddetg3 = INITVALUE;
    CCTK_REAL ddetgt1 = INITVALUE, ddetgt2 = INITVALUE, ddetgt3 = INITVALUE;
    CCTK_REAL detg = INITVALUE;
    CCTK_REAL detgt = INITVALUE;
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
    CCTK_REAL gK112 = INITVALUE, gK113 = INITVALUE, gK121 = INITVALUE, gK122 = INITVALUE, gK123 = INITVALUE, gK131 = INITVALUE;
    CCTK_REAL gK132 = INITVALUE, gK133 = INITVALUE, gK221 = INITVALUE, gK223 = INITVALUE, gK231 = INITVALUE, gK232 = INITVALUE;
    CCTK_REAL gK233 = INITVALUE, gK331 = INITVALUE, gK332 = INITVALUE;
    CCTK_REAL Gt111 = INITVALUE, Gt112 = INITVALUE, Gt113 = INITVALUE, Gt122 = INITVALUE, Gt123 = INITVALUE, Gt133 = INITVALUE;
    CCTK_REAL Gt211 = INITVALUE, Gt212 = INITVALUE, Gt213 = INITVALUE, Gt222 = INITVALUE, Gt223 = INITVALUE, Gt233 = INITVALUE;
    CCTK_REAL Gt311 = INITVALUE, Gt312 = INITVALUE, Gt313 = INITVALUE, Gt322 = INITVALUE, Gt323 = INITVALUE, Gt333 = INITVALUE;
    CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    CCTK_REAL R11 = INITVALUE, R12 = INITVALUE, R13 = INITVALUE, R22 = INITVALUE, R23 = INITVALUE, R33 = INITVALUE;
    CCTK_REAL Rphi11 = INITVALUE, Rphi12 = INITVALUE, Rphi13 = INITVALUE, Rphi22 = INITVALUE, Rphi23 = INITVALUE, Rphi33 = INITVALUE;
    CCTK_REAL Rt11 = INITVALUE, Rt12 = INITVALUE, Rt13 = INITVALUE, Rt22 = INITVALUE, Rt23 = INITVALUE, Rt33 = INITVALUE;
    CCTK_REAL trR = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL At11L = INITVALUE, At12L = INITVALUE, At13L = INITVALUE, At22L = INITVALUE, At23L = INITVALUE, At33L = INITVALUE;
    CCTK_REAL cAL = INITVALUE;
    CCTK_REAL cSL = INITVALUE;
    CCTK_REAL cXt1L = INITVALUE, cXt2L = INITVALUE, cXt3L = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt12L = INITVALUE, gt13L = INITVALUE, gt22L = INITVALUE, gt23L = INITVALUE, gt33L = INITVALUE;
    CCTK_REAL HL = INITVALUE;
    CCTK_REAL M1L = INITVALUE, M2L = INITVALUE, M3L = INITVALUE;
    CCTK_REAL phiL = INITVALUE;
    CCTK_REAL trKL = INITVALUE;
    CCTK_REAL Xt1L = INITVALUE, Xt2L = INITVALUE, Xt3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandard4th2At11 = INITVALUE;
    CCTK_REAL PDstandard4th3At11 = INITVALUE;
    CCTK_REAL PDstandard4th1At12 = INITVALUE;
    CCTK_REAL PDstandard4th2At12 = INITVALUE;
    CCTK_REAL PDstandard4th3At12 = INITVALUE;
    CCTK_REAL PDstandard4th1At13 = INITVALUE;
    CCTK_REAL PDstandard4th2At13 = INITVALUE;
    CCTK_REAL PDstandard4th3At13 = INITVALUE;
    CCTK_REAL PDstandard4th1At22 = INITVALUE;
    CCTK_REAL PDstandard4th3At22 = INITVALUE;
    CCTK_REAL PDstandard4th1At23 = INITVALUE;
    CCTK_REAL PDstandard4th2At23 = INITVALUE;
    CCTK_REAL PDstandard4th3At23 = INITVALUE;
    CCTK_REAL PDstandard4th1At33 = INITVALUE;
    CCTK_REAL PDstandard4th2At33 = INITVALUE;
    CCTK_REAL PDstandard4th1gt11 = INITVALUE;
    CCTK_REAL PDstandard4th2gt11 = INITVALUE;
    CCTK_REAL PDstandard4th3gt11 = INITVALUE;
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
    CCTK_REAL PDstandard4th22gt13 = INITVALUE;
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
    
    /* Assign local copies of grid functions */
    At11L = At11[index];
    At12L = At12[index];
    At13L = At13[index];
    At22L = At22[index];
    At23L = At23[index];
    At33L = At33[index];
    gt11L = gt11[index];
    gt12L = gt12[index];
    gt13L = gt13[index];
    gt22L = gt22[index];
    gt23L = gt23[index];
    gt33L = gt33[index];
    phiL = phi[index];
    trKL = trK[index];
    Xt1L = Xt1[index];
    Xt2L = Xt2[index];
    Xt3L = Xt3[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandard4th2At11 = PDstandard4th2(At11, i, j, k);
    PDstandard4th3At11 = PDstandard4th3(At11, i, j, k);
    PDstandard4th1At12 = PDstandard4th1(At12, i, j, k);
    PDstandard4th2At12 = PDstandard4th2(At12, i, j, k);
    PDstandard4th3At12 = PDstandard4th3(At12, i, j, k);
    PDstandard4th1At13 = PDstandard4th1(At13, i, j, k);
    PDstandard4th2At13 = PDstandard4th2(At13, i, j, k);
    PDstandard4th3At13 = PDstandard4th3(At13, i, j, k);
    PDstandard4th1At22 = PDstandard4th1(At22, i, j, k);
    PDstandard4th3At22 = PDstandard4th3(At22, i, j, k);
    PDstandard4th1At23 = PDstandard4th1(At23, i, j, k);
    PDstandard4th2At23 = PDstandard4th2(At23, i, j, k);
    PDstandard4th3At23 = PDstandard4th3(At23, i, j, k);
    PDstandard4th1At33 = PDstandard4th1(At33, i, j, k);
    PDstandard4th2At33 = PDstandard4th2(At33, i, j, k);
    PDstandard4th1gt11 = PDstandard4th1(gt11, i, j, k);
    PDstandard4th2gt11 = PDstandard4th2(gt11, i, j, k);
    PDstandard4th3gt11 = PDstandard4th3(gt11, i, j, k);
    PDstandard4th22gt11 = PDstandard4th22(gt11, i, j, k);
    PDstandard4th33gt11 = PDstandard4th33(gt11, i, j, k);
    PDstandard4th23gt11 = PDstandard4th23(gt11, i, j, k);
    PDstandard4th1gt12 = PDstandard4th1(gt12, i, j, k);
    PDstandard4th2gt12 = PDstandard4th2(gt12, i, j, k);
    PDstandard4th3gt12 = PDstandard4th3(gt12, i, j, k);
    PDstandard4th33gt12 = PDstandard4th33(gt12, i, j, k);
    PDstandard4th12gt12 = PDstandard4th12(gt12, i, j, k);
    PDstandard4th13gt12 = PDstandard4th13(gt12, i, j, k);
    PDstandard4th23gt12 = PDstandard4th23(gt12, i, j, k);
    PDstandard4th1gt13 = PDstandard4th1(gt13, i, j, k);
    PDstandard4th2gt13 = PDstandard4th2(gt13, i, j, k);
    PDstandard4th3gt13 = PDstandard4th3(gt13, i, j, k);
    PDstandard4th22gt13 = PDstandard4th22(gt13, i, j, k);
    PDstandard4th12gt13 = PDstandard4th12(gt13, i, j, k);
    PDstandard4th13gt13 = PDstandard4th13(gt13, i, j, k);
    PDstandard4th23gt13 = PDstandard4th23(gt13, i, j, k);
    PDstandard4th1gt22 = PDstandard4th1(gt22, i, j, k);
    PDstandard4th2gt22 = PDstandard4th2(gt22, i, j, k);
    PDstandard4th3gt22 = PDstandard4th3(gt22, i, j, k);
    PDstandard4th11gt22 = PDstandard4th11(gt22, i, j, k);
    PDstandard4th33gt22 = PDstandard4th33(gt22, i, j, k);
    PDstandard4th13gt22 = PDstandard4th13(gt22, i, j, k);
    PDstandard4th1gt23 = PDstandard4th1(gt23, i, j, k);
    PDstandard4th2gt23 = PDstandard4th2(gt23, i, j, k);
    PDstandard4th3gt23 = PDstandard4th3(gt23, i, j, k);
    PDstandard4th11gt23 = PDstandard4th11(gt23, i, j, k);
    PDstandard4th12gt23 = PDstandard4th12(gt23, i, j, k);
    PDstandard4th13gt23 = PDstandard4th13(gt23, i, j, k);
    PDstandard4th23gt23 = PDstandard4th23(gt23, i, j, k);
    PDstandard4th1gt33 = PDstandard4th1(gt33, i, j, k);
    PDstandard4th2gt33 = PDstandard4th2(gt33, i, j, k);
    PDstandard4th3gt33 = PDstandard4th3(gt33, i, j, k);
    PDstandard4th11gt33 = PDstandard4th11(gt33, i, j, k);
    PDstandard4th22gt33 = PDstandard4th22(gt33, i, j, k);
    PDstandard4th12gt33 = PDstandard4th12(gt33, i, j, k);
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
    
    Rt11  =  khalf*(-(gtu22*PDstandard4th11gt22) - 2*(Gt111*Gt122 + Gt111*Gt133 + Gt211*Gt222 + Gt211*Gt233 + Gt311*Gt322 + 
             Gt311*Gt333 + gtu32*PDstandard4th11gt23) - gtu33*PDstandard4th11gt33 + 2*gtu22*PDstandard4th12gt12 + 
          2*gtu32*PDstandard4th12gt13 + 2*gtu32*PDstandard4th13gt12 + 2*gtu33*PDstandard4th13gt13 - 
          gtu22*PDstandard4th22gt11 - 2*gtu32*PDstandard4th23gt11 - gtu33*PDstandard4th33gt11 + 2*SQR(Gt112) + 
          2*SQR(Gt113) + 2*SQR(Gt212) + 2*SQR(Gt213) + 2*SQR(Gt312) + 2*SQR(Gt313));
    
    Rt12  =  khalf*(2*(Gt113*Gt123 + Gt213*Gt223 + Gt313*Gt323) - 
          2*(Gt112*Gt133 + Gt212*Gt233 + Gt312*Gt333 + gtu21*PDstandard4th12gt12) - gtu32*PDstandard4th12gt23 - 
          gtu33*PDstandard4th12gt33 + gtu31*(PDstandard4th11gt23 - PDstandard4th12gt13 - PDstandard4th13gt12) + 
          gtu32*PDstandard4th13gt22 + gtu33*PDstandard4th13gt23 + gtu21*(PDstandard4th11gt22 + PDstandard4th22gt11) + 
          gtu32*PDstandard4th22gt13 + gtu31*PDstandard4th23gt11 - gtu32*PDstandard4th23gt12 + gtu33*PDstandard4th23gt13 - 
          gtu33*PDstandard4th33gt12);
    
    Rt13  =  khalf*(2*(Gt112*Gt123 + Gt212*Gt223 + Gt312*Gt323) - 
          2*(Gt113*Gt122 + Gt213*Gt222 + Gt313*Gt322 + gtu31*PDstandard4th13gt13) + 
          gtu21*(PDstandard4th11gt23 - PDstandard4th12gt13 - PDstandard4th13gt12 + PDstandard4th23gt11) + 
          gtu22*(PDstandard4th12gt23 - PDstandard4th13gt22 - PDstandard4th22gt13 + PDstandard4th23gt12) + 
          gtu31*(PDstandard4th11gt33 + PDstandard4th33gt11) + 
          gtu32*(PDstandard4th12gt33 - PDstandard4th13gt23 - PDstandard4th23gt13 + PDstandard4th33gt12));
    
    Rt22  =  khalf*(-2*(Gt122*(Gt111 + Gt133) + Gt222*(Gt211 + Gt233) + Gt322*(Gt311 + Gt333) + gtu31*PDstandard4th13gt22) + 
          gtu11*(-PDstandard4th11gt22 + 2*PDstandard4th12gt12 - PDstandard4th22gt11) + 
          gtu31*(-2*PDstandard4th22gt13 + 2*(PDstandard4th12gt23 + PDstandard4th23gt12)) + 
          gtu33*(-PDstandard4th22gt33 + 2*PDstandard4th23gt23 - PDstandard4th33gt22) + 
          2*(SQR(Gt112) + SQR(Gt123) + SQR(Gt212) + SQR(Gt223) + SQR(Gt312) + SQR(Gt323)));
    
    Rt23  =  khalf*(2*(Gt112*Gt113 + Gt212*Gt213 + Gt312*Gt313) + 
          gtu11*(-PDstandard4th11gt23 + PDstandard4th12gt13 + PDstandard4th13gt12 - PDstandard4th23gt11) + 
          gtu21*(-PDstandard4th12gt23 + PDstandard4th13gt22 + PDstandard4th22gt13 - PDstandard4th23gt12) - 
          2*(Gt111*Gt123 + Gt211*Gt223 + Gt311*Gt323 + gtu32*PDstandard4th23gt23) + 
          gtu31*(PDstandard4th12gt33 - PDstandard4th13gt23 - PDstandard4th23gt13 + PDstandard4th33gt12) + 
          gtu32*(PDstandard4th22gt33 + PDstandard4th33gt22));
    
    Rt33  =  khalf*(gtu11*(-PDstandard4th11gt33 + 2*PDstandard4th13gt13 - PDstandard4th33gt11) - 
          2*((Gt111 + Gt122)*Gt133 + (Gt211 + Gt222)*Gt233 + (Gt311 + Gt322)*Gt333 + 
             gtu21*(PDstandard4th12gt33 + PDstandard4th33gt12)) + 
          gtu22*(-PDstandard4th22gt33 + 2*PDstandard4th23gt23 - PDstandard4th33gt22) + 
          2*(gtu21*(PDstandard4th13gt23 + PDstandard4th23gt13) + SQR(Gt113) + SQR(Gt123) + SQR(Gt213) + SQR(Gt223) + 
             SQR(Gt313) + SQR(Gt323)));
    
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
    
    e4phi  =  exp(4*phiL);
    
    em4phi  =  INV(e4phi);
    
    g11  =  e4phi*gt11L;
    
    g12  =  e4phi*gt12L;
    
    g13  =  e4phi*gt13L;
    
    g22  =  e4phi*gt22L;
    
    g23  =  e4phi*gt23L;
    
    g33  =  e4phi*gt33L;
    
    detg  =  CUB(e4phi);
    
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
    
    trR  =  gu11*R11 + gu22*R22 + 2*(gu21*R12 + gu31*R13 + gu32*R23) + gu33*R33;
    
    Atm11  =  At11L*gtu11 + At12L*gtu21 + At13L*gtu31;
    
    Atm21  =  At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    Atm31  =  At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    Atm12  =  At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    Atm22  =  At12L*gtu21 + At22L*gtu22 + At23L*gtu32;
    
    Atm32  =  At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    Atm13  =  At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    Atm23  =  At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    Atm33  =  At13L*gtu31 + At23L*gtu32 + At33L*gtu33;
    
    HL  =  -2*(Atm12*Atm21 + Atm13*Atm31 + Atm23*Atm32) + trR - SQR(Atm11) - SQR(Atm22) - SQR(Atm33) + ktwothird*SQR(trKL);
    
    gK112  =  e4phi*(-2*(At11L*G112 + At12L*G212 + At13L*G312) + PDstandard4th2At11 + 4*At11L*PDstandard4th2phi) + 
        g11*kthird*PDstandard4th2trK;
    
    gK113  =  e4phi*(-2*(At11L*G113 + At12L*G213 + At13L*G313) + PDstandard4th3At11 + 4*At11L*PDstandard4th3phi) + 
        g11*kthird*PDstandard4th3trK;
    
    gK121  =  e4phi*(-(At11L*G112) - At22L*G211 - At12L*(G111 + G212) - At23L*G311 - At13L*G312 + PDstandard4th1At12 + 
           4*At12L*PDstandard4th1phi) + g12*kthird*PDstandard4th1trK;
    
    gK122  =  e4phi*(-(At11L*G122) - At22L*G212 - At12L*(G112 + G222) - At23L*G312 - At13L*G322 + PDstandard4th2At12 + 
           4*At12L*PDstandard4th2phi) + g12*kthird*PDstandard4th2trK;
    
    gK123  =  e4phi*(-(At11L*G123) - At22L*G213 - At12L*(G113 + G223) - At23L*G313 - At13L*G323 + PDstandard4th3At12 + 
           4*At12L*PDstandard4th3phi) + g12*kthird*PDstandard4th3trK;
    
    gK131  =  e4phi*(-(At11L*G113) - At23L*G211 - At12L*G213 - At33L*G311 - At13L*(G111 + G313) + PDstandard4th1At13 + 
           4*At13L*PDstandard4th1phi) + g13*kthird*PDstandard4th1trK;
    
    gK132  =  e4phi*(-(At11L*G123) - At23L*G212 - At12L*G223 - At33L*G312 - At13L*(G112 + G323) + PDstandard4th2At13 + 
           4*At13L*PDstandard4th2phi) + g13*kthird*PDstandard4th2trK;
    
    gK133  =  e4phi*(-(At11L*G133) - At23L*G213 - At12L*G233 - At33L*G313 - At13L*(G113 + G333) + PDstandard4th3At13 + 
           4*At13L*PDstandard4th3phi) + g13*kthird*PDstandard4th3trK;
    
    gK221  =  e4phi*(-2*(At12L*G112 + At22L*G212 + At23L*G312) + PDstandard4th1At22 + 4*At22L*PDstandard4th1phi) + 
        g22*kthird*PDstandard4th1trK;
    
    gK223  =  e4phi*(-2*(At12L*G123 + At22L*G223 + At23L*G323) + PDstandard4th3At22 + 4*At22L*PDstandard4th3phi) + 
        g22*kthird*PDstandard4th3trK;
    
    gK231  =  e4phi*(-(At13L*G112) - At12L*G113 - At23L*G212 - At22L*G213 - At33L*G312 - At23L*G313 + PDstandard4th1At23 + 
           4*At23L*PDstandard4th1phi) + g23*kthird*PDstandard4th1trK;
    
    gK232  =  e4phi*(-(At13L*G122) - At12L*G123 - At23L*G222 - At22L*G223 - At33L*G322 - At23L*G323 + PDstandard4th2At23 + 
           4*At23L*PDstandard4th2phi) + g23*kthird*PDstandard4th2trK;
    
    gK233  =  e4phi*(-(At13L*G123) - At12L*G133 - At23L*G223 - At22L*G233 - At33L*G323 - At23L*G333 + PDstandard4th3At23 + 
           4*At23L*PDstandard4th3phi) + g23*kthird*PDstandard4th3trK;
    
    gK331  =  e4phi*(-2*(At13L*G113 + At23L*G213 + At33L*G313) + PDstandard4th1At33 + 4*At33L*PDstandard4th1phi) + 
        g33*kthird*PDstandard4th1trK;
    
    gK332  =  e4phi*(-2*(At13L*G123 + At23L*G223 + At33L*G323) + PDstandard4th2At33 + 4*At33L*PDstandard4th2phi) + 
        g33*kthird*PDstandard4th2trK;
    
    M1L  =  (gK112 - gK121)*gu21 + (gK122 - gK221)*gu22 + (gK113 - gK131)*gu31 + (gK123 + gK132 - 2*gK231)*gu32 + 
        (gK133 - gK331)*gu33;
    
    M2L  =  (-gK112 + gK121)*gu11 + (-gK122 + gK221)*gu21 + (gK123 - 2*gK132 + gK231)*gu31 + (gK223 - gK232)*gu32 + 
        (gK233 - gK332)*gu33;
    
    M3L  =  (-gK113 + gK131)*gu11 + (-2*gK123 + gK132 + gK231)*gu21 + (-gK223 + gK232)*gu22 + (-gK133 + gK331)*gu31 + 
        (-gK233 + gK332)*gu32;
    
    cSL  =  trR;
    
    cXt1L  =  Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33 - Xt1L;
    
    cXt2L  =  Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33 - Xt2L;
    
    cXt3L  =  Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33 - Xt3L;
    
    cAL  =  At11L*gtu11 + At22L*gtu22 + 2*(At12L*gtu21 + At13L*gtu31 + At23L*gtu32) + At33L*gtu33;
    
    
    /* Copy local copies back to grid functions */
    cA[index] = cAL;
    cS[index] = cSL;
    cXt1[index] = cXt1L;
    cXt2[index] = cXt2L;
    cXt3[index] = cXt3L;
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (somename);
}

void ML_BSSN_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_constraints_Body);
}
