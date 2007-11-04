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

void BSSN_convertFromADMBaseGamma_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering BSSN_convertFromADMBaseGamma_Body");
  }
  
  if (cctk_iteration % BSSN_convertFromADMBaseGamma_calc_every != BSSN_convertFromADMBaseGamma_calc_offset)
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
    CCTK_REAL detgt = INITVALUE;
    CCTK_REAL Gt111 = INITVALUE, Gt121 = INITVALUE, Gt122 = INITVALUE, Gt131 = INITVALUE, Gt132 = INITVALUE, Gt133 = INITVALUE;
    CCTK_REAL Gt211 = INITVALUE, Gt221 = INITVALUE, Gt222 = INITVALUE, Gt231 = INITVALUE, Gt232 = INITVALUE, Gt233 = INITVALUE;
    CCTK_REAL Gt311 = INITVALUE, Gt321 = INITVALUE, Gt322 = INITVALUE, Gt331 = INITVALUE, Gt332 = INITVALUE, Gt333 = INITVALUE;
    CCTK_REAL gtu11 = INITVALUE, gtu12 = INITVALUE, gtu13 = INITVALUE, gtu22 = INITVALUE, gtu23 = INITVALUE, gtu33 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL gt11L = INITVALUE, gt21L = INITVALUE, gt22L = INITVALUE, gt31L = INITVALUE, gt32L = INITVALUE, gt33L = INITVALUE;
    CCTK_REAL Xt1L = INITVALUE, Xt2L = INITVALUE, Xt3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandard4th1gt11 = INITVALUE;
    CCTK_REAL PDstandard4th2gt11 = INITVALUE;
    CCTK_REAL PDstandard4th3gt11 = INITVALUE;
    CCTK_REAL PDstandard4th1gt21 = INITVALUE;
    CCTK_REAL PDstandard4th2gt21 = INITVALUE;
    CCTK_REAL PDstandard4th3gt21 = INITVALUE;
    CCTK_REAL PDstandard4th1gt22 = INITVALUE;
    CCTK_REAL PDstandard4th2gt22 = INITVALUE;
    CCTK_REAL PDstandard4th3gt22 = INITVALUE;
    CCTK_REAL PDstandard4th1gt31 = INITVALUE;
    CCTK_REAL PDstandard4th2gt31 = INITVALUE;
    CCTK_REAL PDstandard4th3gt31 = INITVALUE;
    CCTK_REAL PDstandard4th1gt32 = INITVALUE;
    CCTK_REAL PDstandard4th2gt32 = INITVALUE;
    CCTK_REAL PDstandard4th3gt32 = INITVALUE;
    CCTK_REAL PDstandard4th1gt33 = INITVALUE;
    CCTK_REAL PDstandard4th2gt33 = INITVALUE;
    CCTK_REAL PDstandard4th3gt33 = INITVALUE;
    
    /* Assign local copies of grid functions */
    gt11L = gt11[index];
    gt21L = gt21[index];
    gt22L = gt22[index];
    gt31L = gt31[index];
    gt32L = gt32[index];
    gt33L = gt33[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandard4th1gt11 = PDstandard4th1(gt11, i, j, k);
    PDstandard4th2gt11 = PDstandard4th2(gt11, i, j, k);
    PDstandard4th3gt11 = PDstandard4th3(gt11, i, j, k);
    PDstandard4th1gt21 = PDstandard4th1(gt21, i, j, k);
    PDstandard4th2gt21 = PDstandard4th2(gt21, i, j, k);
    PDstandard4th3gt21 = PDstandard4th3(gt21, i, j, k);
    PDstandard4th1gt22 = PDstandard4th1(gt22, i, j, k);
    PDstandard4th2gt22 = PDstandard4th2(gt22, i, j, k);
    PDstandard4th3gt22 = PDstandard4th3(gt22, i, j, k);
    PDstandard4th1gt31 = PDstandard4th1(gt31, i, j, k);
    PDstandard4th2gt31 = PDstandard4th2(gt31, i, j, k);
    PDstandard4th3gt31 = PDstandard4th3(gt31, i, j, k);
    PDstandard4th1gt32 = PDstandard4th1(gt32, i, j, k);
    PDstandard4th2gt32 = PDstandard4th2(gt32, i, j, k);
    PDstandard4th3gt32 = PDstandard4th3(gt32, i, j, k);
    PDstandard4th1gt33 = PDstandard4th1(gt33, i, j, k);
    PDstandard4th2gt33 = PDstandard4th2(gt33, i, j, k);
    PDstandard4th3gt33 = PDstandard4th3(gt33, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    detgt  =  2*gt21L*gt31L*gt32L + gt33L*(gt11L*gt22L - SQR(gt21L)) - gt22L*SQR(gt31L) - gt11L*SQR(gt32L);
    
    gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt32L));
    
    gtu12  =  (gt31L*gt32L - gt21L*gt33L)*INV(detgt);
    
    gtu13  =  (-(gt22L*gt31L) + gt21L*gt32L)*INV(detgt);
    
    gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt31L));
    
    gtu23  =  (gt21L*gt31L - gt11L*gt32L)*INV(detgt);
    
    gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt21L));
    
    Gt111  =  khalf*(gtu11*PDstandard4th1gt11 + 2*(gtu12*PDstandard4th1gt21 + gtu13*PDstandard4th1gt31) - 
          gtu12*PDstandard4th2gt11 - gtu13*PDstandard4th3gt11);
    
    Gt211  =  khalf*(gtu12*PDstandard4th1gt11 + 2*(gtu22*PDstandard4th1gt21 + gtu23*PDstandard4th1gt31) - 
          gtu22*PDstandard4th2gt11 - gtu23*PDstandard4th3gt11);
    
    Gt311  =  khalf*(gtu13*PDstandard4th1gt11 + 2*(gtu23*PDstandard4th1gt21 + gtu33*PDstandard4th1gt31) - 
          gtu23*PDstandard4th2gt11 - gtu33*PDstandard4th3gt11);
    
    Gt121  =  khalf*(gtu12*PDstandard4th1gt22 + gtu11*PDstandard4th2gt11 + 
          gtu13*(PDstandard4th1gt32 + PDstandard4th2gt31 - PDstandard4th3gt21));
    
    Gt221  =  khalf*(gtu22*PDstandard4th1gt22 + gtu12*PDstandard4th2gt11 + 
          gtu23*(PDstandard4th1gt32 + PDstandard4th2gt31 - PDstandard4th3gt21));
    
    Gt321  =  khalf*(gtu23*PDstandard4th1gt22 + gtu13*PDstandard4th2gt11 + 
          gtu33*(PDstandard4th1gt32 + PDstandard4th2gt31 - PDstandard4th3gt21));
    
    Gt131  =  khalf*(gtu13*PDstandard4th1gt33 + gtu11*PDstandard4th3gt11 + 
          gtu12*(PDstandard4th1gt32 - PDstandard4th2gt31 + PDstandard4th3gt21));
    
    Gt231  =  khalf*(gtu23*PDstandard4th1gt33 + gtu12*PDstandard4th3gt11 + 
          gtu22*(PDstandard4th1gt32 - PDstandard4th2gt31 + PDstandard4th3gt21));
    
    Gt331  =  khalf*(gtu33*PDstandard4th1gt33 + gtu13*PDstandard4th3gt11 + 
          gtu23*(PDstandard4th1gt32 - PDstandard4th2gt31 + PDstandard4th3gt21));
    
    Gt122  =  khalf*(gtu11*(-PDstandard4th1gt22 + 2*PDstandard4th2gt21) + gtu12*PDstandard4th2gt22 + 
          gtu13*(2*PDstandard4th2gt32 - PDstandard4th3gt22));
    
    Gt222  =  khalf*(gtu12*(-PDstandard4th1gt22 + 2*PDstandard4th2gt21) + gtu22*PDstandard4th2gt22 + 
          gtu23*(2*PDstandard4th2gt32 - PDstandard4th3gt22));
    
    Gt322  =  khalf*(gtu13*(-PDstandard4th1gt22 + 2*PDstandard4th2gt21) + gtu23*PDstandard4th2gt22 + 
          gtu33*(2*PDstandard4th2gt32 - PDstandard4th3gt22));
    
    Gt132  =  khalf*(gtu13*PDstandard4th2gt33 + gtu11*(-PDstandard4th1gt32 + PDstandard4th2gt31 + PDstandard4th3gt21) + 
          gtu12*PDstandard4th3gt22);
    
    Gt232  =  khalf*(gtu23*PDstandard4th2gt33 + gtu12*(-PDstandard4th1gt32 + PDstandard4th2gt31 + PDstandard4th3gt21) + 
          gtu22*PDstandard4th3gt22);
    
    Gt332  =  khalf*(gtu33*PDstandard4th2gt33 + gtu13*(-PDstandard4th1gt32 + PDstandard4th2gt31 + PDstandard4th3gt21) + 
          gtu23*PDstandard4th3gt22);
    
    Gt133  =  khalf*(-(gtu11*PDstandard4th1gt33) - gtu12*PDstandard4th2gt33 + 2*gtu11*PDstandard4th3gt31 + 
          2*gtu12*PDstandard4th3gt32 + gtu13*PDstandard4th3gt33);
    
    Gt233  =  khalf*(-(gtu12*PDstandard4th1gt33) - gtu22*PDstandard4th2gt33 + 2*gtu12*PDstandard4th3gt31 + 
          2*gtu22*PDstandard4th3gt32 + gtu23*PDstandard4th3gt33);
    
    Gt333  =  khalf*(-(gtu13*PDstandard4th1gt33) - gtu23*PDstandard4th2gt33 + 2*gtu13*PDstandard4th3gt31 + 
          2*gtu23*PDstandard4th3gt32 + gtu33*PDstandard4th3gt33);
    
    Xt1L  =  Gt111*gtu11 + Gt122*gtu22 + 2*(Gt121*gtu12 + Gt131*gtu13 + Gt132*gtu23) + Gt133*gtu33;
    
    Xt2L  =  Gt211*gtu11 + Gt222*gtu22 + 2*(Gt221*gtu12 + Gt231*gtu13 + Gt232*gtu23) + Gt233*gtu33;
    
    Xt3L  =  Gt311*gtu11 + Gt322*gtu22 + 2*(Gt321*gtu12 + Gt331*gtu13 + Gt332*gtu23) + Gt333*gtu33;
    
    
    /* Copy local copies back to grid functions */
    Xt1[index] = Xt1L;
    Xt2[index] = Xt2L;
    Xt3[index] = Xt3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (somename);
}

void BSSN_convertFromADMBaseGamma(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverEverything(cctkGH, &BSSN_convertFromADMBaseGamma_Body);
}
