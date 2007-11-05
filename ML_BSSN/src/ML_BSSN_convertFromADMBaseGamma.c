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

void ML_BSSN_convertFromADMBaseGamma_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_convertFromADMBaseGamma_Body");
  }
  
  if (cctk_iteration % ML_BSSN_convertFromADMBaseGamma_calc_every != ML_BSSN_convertFromADMBaseGamma_calc_offset)
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
    CCTK_REAL detgt = INITVALUE;
    CCTK_REAL Gt111 = INITVALUE, Gt112 = INITVALUE, Gt113 = INITVALUE, Gt122 = INITVALUE, Gt123 = INITVALUE, Gt133 = INITVALUE;
    CCTK_REAL Gt211 = INITVALUE, Gt212 = INITVALUE, Gt213 = INITVALUE, Gt222 = INITVALUE, Gt223 = INITVALUE, Gt233 = INITVALUE;
    CCTK_REAL Gt311 = INITVALUE, Gt312 = INITVALUE, Gt313 = INITVALUE, Gt322 = INITVALUE, Gt323 = INITVALUE, Gt333 = INITVALUE;
    CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL gt11L = INITVALUE, gt12L = INITVALUE, gt13L = INITVALUE, gt22L = INITVALUE, gt23L = INITVALUE, gt33L = INITVALUE;
    CCTK_REAL Xt1L = INITVALUE, Xt2L = INITVALUE, Xt3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandard4th1gt11 = INITVALUE;
    CCTK_REAL PDstandard4th2gt11 = INITVALUE;
    CCTK_REAL PDstandard4th3gt11 = INITVALUE;
    CCTK_REAL PDstandard4th1gt12 = INITVALUE;
    CCTK_REAL PDstandard4th2gt12 = INITVALUE;
    CCTK_REAL PDstandard4th3gt12 = INITVALUE;
    CCTK_REAL PDstandard4th1gt13 = INITVALUE;
    CCTK_REAL PDstandard4th2gt13 = INITVALUE;
    CCTK_REAL PDstandard4th3gt13 = INITVALUE;
    CCTK_REAL PDstandard4th1gt22 = INITVALUE;
    CCTK_REAL PDstandard4th2gt22 = INITVALUE;
    CCTK_REAL PDstandard4th3gt22 = INITVALUE;
    CCTK_REAL PDstandard4th1gt23 = INITVALUE;
    CCTK_REAL PDstandard4th2gt23 = INITVALUE;
    CCTK_REAL PDstandard4th3gt23 = INITVALUE;
    CCTK_REAL PDstandard4th1gt33 = INITVALUE;
    CCTK_REAL PDstandard4th2gt33 = INITVALUE;
    CCTK_REAL PDstandard4th3gt33 = INITVALUE;
    
    /* Assign local copies of grid functions */
    gt11L = gt11[index];
    gt12L = gt12[index];
    gt13L = gt13[index];
    gt22L = gt22[index];
    gt23L = gt23[index];
    gt33L = gt33[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandard4th1gt11 = PDstandard4th1(gt11, i, j, k);
    PDstandard4th2gt11 = PDstandard4th2(gt11, i, j, k);
    PDstandard4th3gt11 = PDstandard4th3(gt11, i, j, k);
    PDstandard4th1gt12 = PDstandard4th1(gt12, i, j, k);
    PDstandard4th2gt12 = PDstandard4th2(gt12, i, j, k);
    PDstandard4th3gt12 = PDstandard4th3(gt12, i, j, k);
    PDstandard4th1gt13 = PDstandard4th1(gt13, i, j, k);
    PDstandard4th2gt13 = PDstandard4th2(gt13, i, j, k);
    PDstandard4th3gt13 = PDstandard4th3(gt13, i, j, k);
    PDstandard4th1gt22 = PDstandard4th1(gt22, i, j, k);
    PDstandard4th2gt22 = PDstandard4th2(gt22, i, j, k);
    PDstandard4th3gt22 = PDstandard4th3(gt22, i, j, k);
    PDstandard4th1gt23 = PDstandard4th1(gt23, i, j, k);
    PDstandard4th2gt23 = PDstandard4th2(gt23, i, j, k);
    PDstandard4th3gt23 = PDstandard4th3(gt23, i, j, k);
    PDstandard4th1gt33 = PDstandard4th1(gt33, i, j, k);
    PDstandard4th2gt33 = PDstandard4th2(gt33, i, j, k);
    PDstandard4th3gt33 = PDstandard4th3(gt33, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    detgt  =  1;
    
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
    
    Xt1L  =  Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    Xt2L  =  Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    Xt3L  =  Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    
    /* Copy local copies back to grid functions */
    Xt1[index] = Xt1L;
    Xt2[index] = Xt2L;
    Xt3[index] = Xt3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (somename);
}

void ML_BSSN_convertFromADMBaseGamma(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_convertFromADMBaseGamma_Body);
}
