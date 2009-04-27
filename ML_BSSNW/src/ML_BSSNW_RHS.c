/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (64-bit) (May 21, 2008) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

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

void ML_BSSNW_RHS_Body(cGH const * const cctkGH, CCTK_INT const dir, CCTK_INT const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], CCTK_INT const min[3], CCTK_INT const max[3], CCTK_INT const n_subblock_gfs, CCTK_REAL * const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
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
  CCTK_REAL pm1o12dx = INITVALUE;
  CCTK_REAL pm1o12dx2 = INITVALUE;
  CCTK_REAL pm1o12dy = INITVALUE;
  CCTK_REAL pm1o12dy2 = INITVALUE;
  CCTK_REAL pm1o12dz = INITVALUE;
  CCTK_REAL pm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSNW_RHS_Body");
  }
  
  if (cctk_iteration % ML_BSSNW_RHS_calc_every != ML_BSSNW_RHS_calc_offset)
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
  pm1o12dx = -INV(dx)/12.;
  pm1o12dx2 = -pow(dx,-2)/12.;
  pm1o12dy = -INV(dy)/12.;
  pm1o12dy2 = -pow(dy,-2)/12.;
  pm1o12dz = -INV(dz)/12.;
  pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSNW_RHS,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lssh[CCTK_LSSH_IDX(0,0)],cctk_lssh[CCTK_LSSH_IDX(0,1)],cctk_lssh[CCTK_LSSH_IDX(0,2)])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL Atm11 = INITVALUE, Atm12 = INITVALUE, Atm13 = INITVALUE, Atm21 = INITVALUE, Atm22 = INITVALUE, Atm23 = INITVALUE;
    CCTK_REAL Atm31 = INITVALUE, Atm32 = INITVALUE, Atm33 = INITVALUE;
    CCTK_REAL Ats11 = INITVALUE, Ats12 = INITVALUE, Ats13 = INITVALUE, Ats22 = INITVALUE, Ats23 = INITVALUE, Ats33 = INITVALUE;
    CCTK_REAL Atu11 = INITVALUE, Atu21 = INITVALUE, Atu22 = INITVALUE, Atu31 = INITVALUE, Atu32 = INITVALUE, Atu33 = INITVALUE;
    CCTK_REAL betam1 = INITVALUE, betam2 = INITVALUE, betam3 = INITVALUE;
    CCTK_REAL betap1 = INITVALUE, betap2 = INITVALUE, betap3 = INITVALUE;
    CCTK_REAL cdphi211 = INITVALUE, cdphi212 = INITVALUE, cdphi213 = INITVALUE, cdphi222 = INITVALUE, cdphi223 = INITVALUE, cdphi233 = INITVALUE;
    CCTK_REAL detgt = INITVALUE;
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
    CCTK_REAL invW = INITVALUE, invW2 = INITVALUE;
    CCTK_REAL pdphi1 = INITVALUE, pdphi2 = INITVALUE, pdphi3 = INITVALUE;
    CCTK_REAL R11 = INITVALUE, R12 = INITVALUE, R13 = INITVALUE, R22 = INITVALUE, R23 = INITVALUE, R33 = INITVALUE;
    CCTK_REAL Rphi11 = INITVALUE, Rphi12 = INITVALUE, Rphi13 = INITVALUE, Rphi22 = INITVALUE, Rphi23 = INITVALUE, Rphi33 = INITVALUE;
    CCTK_REAL Rt11 = INITVALUE, Rt12 = INITVALUE, Rt13 = INITVALUE, Rt22 = INITVALUE, Rt23 = INITVALUE, Rt33 = INITVALUE;
    CCTK_REAL trAts = INITVALUE;
    CCTK_REAL W2 = INITVALUE;
    CCTK_REAL Xtn1 = INITVALUE, Xtn2 = INITVALUE, Xtn3 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL AL = INITVALUE;
    CCTK_REAL alphaL = INITVALUE, alpharhsL = INITVALUE;
    CCTK_REAL ArhsL = INITVALUE;
    CCTK_REAL At11L = INITVALUE, At11rhsL = INITVALUE, At12L = INITVALUE, At12rhsL = INITVALUE, At13L = INITVALUE, At13rhsL = INITVALUE;
    CCTK_REAL At22L = INITVALUE, At22rhsL = INITVALUE, At23L = INITVALUE, At23rhsL = INITVALUE, At33L = INITVALUE, At33rhsL = INITVALUE;
    CCTK_REAL B1L = INITVALUE, B1rhsL = INITVALUE, B2L = INITVALUE, B2rhsL = INITVALUE, B3L = INITVALUE, B3rhsL = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta1rhsL = INITVALUE, beta2L = INITVALUE, beta2rhsL = INITVALUE, beta3L = INITVALUE, beta3rhsL = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt11rhsL = INITVALUE, gt12L = INITVALUE, gt12rhsL = INITVALUE, gt13L = INITVALUE, gt13rhsL = INITVALUE;
    CCTK_REAL gt22L = INITVALUE, gt22rhsL = INITVALUE, gt23L = INITVALUE, gt23rhsL = INITVALUE, gt33L = INITVALUE, gt33rhsL = INITVALUE;
    CCTK_REAL trKL = INITVALUE, trKrhsL = INITVALUE;
    CCTK_REAL WL = INITVALUE, WrhsL = INITVALUE;
    CCTK_REAL Xt1L = INITVALUE, Xt1rhsL = INITVALUE, Xt2L = INITVALUE, Xt2rhsL = INITVALUE, Xt3L = INITVALUE, Xt3rhsL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandardNth1alpha = INITVALUE;
    CCTK_REAL PDstandardNth2alpha = INITVALUE;
    CCTK_REAL PDstandardNth3alpha = INITVALUE;
    CCTK_REAL PDstandardNth11alpha = INITVALUE;
    CCTK_REAL PDstandardNth22alpha = INITVALUE;
    CCTK_REAL PDstandardNth33alpha = INITVALUE;
    CCTK_REAL PDstandardNth12alpha = INITVALUE;
    CCTK_REAL PDstandardNth13alpha = INITVALUE;
    CCTK_REAL PDstandardNth23alpha = INITVALUE;
    CCTK_REAL PDupwindpNth1alpha = INITVALUE;
    CCTK_REAL PDupwindpNth2alpha = INITVALUE;
    CCTK_REAL PDupwindpNth3alpha = INITVALUE;
    CCTK_REAL PDupwindmNth1alpha = INITVALUE;
    CCTK_REAL PDupwindmNth2alpha = INITVALUE;
    CCTK_REAL PDupwindmNth3alpha = INITVALUE;
    CCTK_REAL PDupwindpNth1At11 = INITVALUE;
    CCTK_REAL PDupwindpNth2At11 = INITVALUE;
    CCTK_REAL PDupwindpNth3At11 = INITVALUE;
    CCTK_REAL PDupwindmNth1At11 = INITVALUE;
    CCTK_REAL PDupwindmNth2At11 = INITVALUE;
    CCTK_REAL PDupwindmNth3At11 = INITVALUE;
    CCTK_REAL PDupwindpNth1At12 = INITVALUE;
    CCTK_REAL PDupwindpNth2At12 = INITVALUE;
    CCTK_REAL PDupwindpNth3At12 = INITVALUE;
    CCTK_REAL PDupwindmNth1At12 = INITVALUE;
    CCTK_REAL PDupwindmNth2At12 = INITVALUE;
    CCTK_REAL PDupwindmNth3At12 = INITVALUE;
    CCTK_REAL PDupwindpNth1At13 = INITVALUE;
    CCTK_REAL PDupwindpNth2At13 = INITVALUE;
    CCTK_REAL PDupwindpNth3At13 = INITVALUE;
    CCTK_REAL PDupwindmNth1At13 = INITVALUE;
    CCTK_REAL PDupwindmNth2At13 = INITVALUE;
    CCTK_REAL PDupwindmNth3At13 = INITVALUE;
    CCTK_REAL PDupwindpNth1At22 = INITVALUE;
    CCTK_REAL PDupwindpNth2At22 = INITVALUE;
    CCTK_REAL PDupwindpNth3At22 = INITVALUE;
    CCTK_REAL PDupwindmNth1At22 = INITVALUE;
    CCTK_REAL PDupwindmNth2At22 = INITVALUE;
    CCTK_REAL PDupwindmNth3At22 = INITVALUE;
    CCTK_REAL PDupwindpNth1At23 = INITVALUE;
    CCTK_REAL PDupwindpNth2At23 = INITVALUE;
    CCTK_REAL PDupwindpNth3At23 = INITVALUE;
    CCTK_REAL PDupwindmNth1At23 = INITVALUE;
    CCTK_REAL PDupwindmNth2At23 = INITVALUE;
    CCTK_REAL PDupwindmNth3At23 = INITVALUE;
    CCTK_REAL PDupwindpNth1At33 = INITVALUE;
    CCTK_REAL PDupwindpNth2At33 = INITVALUE;
    CCTK_REAL PDupwindpNth3At33 = INITVALUE;
    CCTK_REAL PDupwindmNth1At33 = INITVALUE;
    CCTK_REAL PDupwindmNth2At33 = INITVALUE;
    CCTK_REAL PDupwindmNth3At33 = INITVALUE;
    CCTK_REAL PDupwindpNth1B1 = INITVALUE;
    CCTK_REAL PDupwindpNth2B1 = INITVALUE;
    CCTK_REAL PDupwindpNth3B1 = INITVALUE;
    CCTK_REAL PDupwindmNth1B1 = INITVALUE;
    CCTK_REAL PDupwindmNth2B1 = INITVALUE;
    CCTK_REAL PDupwindmNth3B1 = INITVALUE;
    CCTK_REAL PDupwindpNth1B2 = INITVALUE;
    CCTK_REAL PDupwindpNth2B2 = INITVALUE;
    CCTK_REAL PDupwindpNth3B2 = INITVALUE;
    CCTK_REAL PDupwindmNth1B2 = INITVALUE;
    CCTK_REAL PDupwindmNth2B2 = INITVALUE;
    CCTK_REAL PDupwindmNth3B2 = INITVALUE;
    CCTK_REAL PDupwindpNth1B3 = INITVALUE;
    CCTK_REAL PDupwindpNth2B3 = INITVALUE;
    CCTK_REAL PDupwindpNth3B3 = INITVALUE;
    CCTK_REAL PDupwindmNth1B3 = INITVALUE;
    CCTK_REAL PDupwindmNth2B3 = INITVALUE;
    CCTK_REAL PDupwindmNth3B3 = INITVALUE;
    CCTK_REAL PDstandardNth1beta1 = INITVALUE;
    CCTK_REAL PDstandardNth2beta1 = INITVALUE;
    CCTK_REAL PDstandardNth3beta1 = INITVALUE;
    CCTK_REAL PDstandardNth11beta1 = INITVALUE;
    CCTK_REAL PDstandardNth22beta1 = INITVALUE;
    CCTK_REAL PDstandardNth33beta1 = INITVALUE;
    CCTK_REAL PDstandardNth12beta1 = INITVALUE;
    CCTK_REAL PDstandardNth13beta1 = INITVALUE;
    CCTK_REAL PDstandardNth23beta1 = INITVALUE;
    CCTK_REAL PDupwindpNth1beta1 = INITVALUE;
    CCTK_REAL PDupwindpNth2beta1 = INITVALUE;
    CCTK_REAL PDupwindpNth3beta1 = INITVALUE;
    CCTK_REAL PDupwindmNth1beta1 = INITVALUE;
    CCTK_REAL PDupwindmNth2beta1 = INITVALUE;
    CCTK_REAL PDupwindmNth3beta1 = INITVALUE;
    CCTK_REAL PDstandardNth1beta2 = INITVALUE;
    CCTK_REAL PDstandardNth2beta2 = INITVALUE;
    CCTK_REAL PDstandardNth3beta2 = INITVALUE;
    CCTK_REAL PDstandardNth11beta2 = INITVALUE;
    CCTK_REAL PDstandardNth22beta2 = INITVALUE;
    CCTK_REAL PDstandardNth33beta2 = INITVALUE;
    CCTK_REAL PDstandardNth12beta2 = INITVALUE;
    CCTK_REAL PDstandardNth13beta2 = INITVALUE;
    CCTK_REAL PDstandardNth23beta2 = INITVALUE;
    CCTK_REAL PDupwindpNth1beta2 = INITVALUE;
    CCTK_REAL PDupwindpNth2beta2 = INITVALUE;
    CCTK_REAL PDupwindpNth3beta2 = INITVALUE;
    CCTK_REAL PDupwindmNth1beta2 = INITVALUE;
    CCTK_REAL PDupwindmNth2beta2 = INITVALUE;
    CCTK_REAL PDupwindmNth3beta2 = INITVALUE;
    CCTK_REAL PDstandardNth1beta3 = INITVALUE;
    CCTK_REAL PDstandardNth2beta3 = INITVALUE;
    CCTK_REAL PDstandardNth3beta3 = INITVALUE;
    CCTK_REAL PDstandardNth11beta3 = INITVALUE;
    CCTK_REAL PDstandardNth22beta3 = INITVALUE;
    CCTK_REAL PDstandardNth33beta3 = INITVALUE;
    CCTK_REAL PDstandardNth12beta3 = INITVALUE;
    CCTK_REAL PDstandardNth13beta3 = INITVALUE;
    CCTK_REAL PDstandardNth23beta3 = INITVALUE;
    CCTK_REAL PDupwindpNth1beta3 = INITVALUE;
    CCTK_REAL PDupwindpNth2beta3 = INITVALUE;
    CCTK_REAL PDupwindpNth3beta3 = INITVALUE;
    CCTK_REAL PDupwindmNth1beta3 = INITVALUE;
    CCTK_REAL PDupwindmNth2beta3 = INITVALUE;
    CCTK_REAL PDupwindmNth3beta3 = INITVALUE;
    CCTK_REAL PDstandardNth1gt11 = INITVALUE;
    CCTK_REAL PDstandardNth2gt11 = INITVALUE;
    CCTK_REAL PDstandardNth3gt11 = INITVALUE;
    CCTK_REAL PDstandardNth11gt11 = INITVALUE;
    CCTK_REAL PDstandardNth22gt11 = INITVALUE;
    CCTK_REAL PDstandardNth33gt11 = INITVALUE;
    CCTK_REAL PDstandardNth12gt11 = INITVALUE;
    CCTK_REAL PDstandardNth13gt11 = INITVALUE;
    CCTK_REAL PDstandardNth23gt11 = INITVALUE;
    CCTK_REAL PDupwindpNth1gt11 = INITVALUE;
    CCTK_REAL PDupwindpNth2gt11 = INITVALUE;
    CCTK_REAL PDupwindpNth3gt11 = INITVALUE;
    CCTK_REAL PDupwindmNth1gt11 = INITVALUE;
    CCTK_REAL PDupwindmNth2gt11 = INITVALUE;
    CCTK_REAL PDupwindmNth3gt11 = INITVALUE;
    CCTK_REAL PDstandardNth1gt12 = INITVALUE;
    CCTK_REAL PDstandardNth2gt12 = INITVALUE;
    CCTK_REAL PDstandardNth3gt12 = INITVALUE;
    CCTK_REAL PDstandardNth11gt12 = INITVALUE;
    CCTK_REAL PDstandardNth22gt12 = INITVALUE;
    CCTK_REAL PDstandardNth33gt12 = INITVALUE;
    CCTK_REAL PDstandardNth12gt12 = INITVALUE;
    CCTK_REAL PDstandardNth13gt12 = INITVALUE;
    CCTK_REAL PDstandardNth23gt12 = INITVALUE;
    CCTK_REAL PDupwindpNth1gt12 = INITVALUE;
    CCTK_REAL PDupwindpNth2gt12 = INITVALUE;
    CCTK_REAL PDupwindpNth3gt12 = INITVALUE;
    CCTK_REAL PDupwindmNth1gt12 = INITVALUE;
    CCTK_REAL PDupwindmNth2gt12 = INITVALUE;
    CCTK_REAL PDupwindmNth3gt12 = INITVALUE;
    CCTK_REAL PDstandardNth1gt13 = INITVALUE;
    CCTK_REAL PDstandardNth2gt13 = INITVALUE;
    CCTK_REAL PDstandardNth3gt13 = INITVALUE;
    CCTK_REAL PDstandardNth11gt13 = INITVALUE;
    CCTK_REAL PDstandardNth22gt13 = INITVALUE;
    CCTK_REAL PDstandardNth33gt13 = INITVALUE;
    CCTK_REAL PDstandardNth12gt13 = INITVALUE;
    CCTK_REAL PDstandardNth13gt13 = INITVALUE;
    CCTK_REAL PDstandardNth23gt13 = INITVALUE;
    CCTK_REAL PDupwindpNth1gt13 = INITVALUE;
    CCTK_REAL PDupwindpNth2gt13 = INITVALUE;
    CCTK_REAL PDupwindpNth3gt13 = INITVALUE;
    CCTK_REAL PDupwindmNth1gt13 = INITVALUE;
    CCTK_REAL PDupwindmNth2gt13 = INITVALUE;
    CCTK_REAL PDupwindmNth3gt13 = INITVALUE;
    CCTK_REAL PDstandardNth1gt22 = INITVALUE;
    CCTK_REAL PDstandardNth2gt22 = INITVALUE;
    CCTK_REAL PDstandardNth3gt22 = INITVALUE;
    CCTK_REAL PDstandardNth11gt22 = INITVALUE;
    CCTK_REAL PDstandardNth22gt22 = INITVALUE;
    CCTK_REAL PDstandardNth33gt22 = INITVALUE;
    CCTK_REAL PDstandardNth12gt22 = INITVALUE;
    CCTK_REAL PDstandardNth13gt22 = INITVALUE;
    CCTK_REAL PDstandardNth23gt22 = INITVALUE;
    CCTK_REAL PDupwindpNth1gt22 = INITVALUE;
    CCTK_REAL PDupwindpNth2gt22 = INITVALUE;
    CCTK_REAL PDupwindpNth3gt22 = INITVALUE;
    CCTK_REAL PDupwindmNth1gt22 = INITVALUE;
    CCTK_REAL PDupwindmNth2gt22 = INITVALUE;
    CCTK_REAL PDupwindmNth3gt22 = INITVALUE;
    CCTK_REAL PDstandardNth1gt23 = INITVALUE;
    CCTK_REAL PDstandardNth2gt23 = INITVALUE;
    CCTK_REAL PDstandardNth3gt23 = INITVALUE;
    CCTK_REAL PDstandardNth11gt23 = INITVALUE;
    CCTK_REAL PDstandardNth22gt23 = INITVALUE;
    CCTK_REAL PDstandardNth33gt23 = INITVALUE;
    CCTK_REAL PDstandardNth12gt23 = INITVALUE;
    CCTK_REAL PDstandardNth13gt23 = INITVALUE;
    CCTK_REAL PDstandardNth23gt23 = INITVALUE;
    CCTK_REAL PDupwindpNth1gt23 = INITVALUE;
    CCTK_REAL PDupwindpNth2gt23 = INITVALUE;
    CCTK_REAL PDupwindpNth3gt23 = INITVALUE;
    CCTK_REAL PDupwindmNth1gt23 = INITVALUE;
    CCTK_REAL PDupwindmNth2gt23 = INITVALUE;
    CCTK_REAL PDupwindmNth3gt23 = INITVALUE;
    CCTK_REAL PDstandardNth1gt33 = INITVALUE;
    CCTK_REAL PDstandardNth2gt33 = INITVALUE;
    CCTK_REAL PDstandardNth3gt33 = INITVALUE;
    CCTK_REAL PDstandardNth11gt33 = INITVALUE;
    CCTK_REAL PDstandardNth22gt33 = INITVALUE;
    CCTK_REAL PDstandardNth33gt33 = INITVALUE;
    CCTK_REAL PDstandardNth12gt33 = INITVALUE;
    CCTK_REAL PDstandardNth13gt33 = INITVALUE;
    CCTK_REAL PDstandardNth23gt33 = INITVALUE;
    CCTK_REAL PDupwindpNth1gt33 = INITVALUE;
    CCTK_REAL PDupwindpNth2gt33 = INITVALUE;
    CCTK_REAL PDupwindpNth3gt33 = INITVALUE;
    CCTK_REAL PDupwindmNth1gt33 = INITVALUE;
    CCTK_REAL PDupwindmNth2gt33 = INITVALUE;
    CCTK_REAL PDupwindmNth3gt33 = INITVALUE;
    CCTK_REAL PDstandardNth1trK = INITVALUE;
    CCTK_REAL PDstandardNth2trK = INITVALUE;
    CCTK_REAL PDstandardNth3trK = INITVALUE;
    CCTK_REAL PDupwindpNth1trK = INITVALUE;
    CCTK_REAL PDupwindpNth2trK = INITVALUE;
    CCTK_REAL PDupwindpNth3trK = INITVALUE;
    CCTK_REAL PDupwindmNth1trK = INITVALUE;
    CCTK_REAL PDupwindmNth2trK = INITVALUE;
    CCTK_REAL PDupwindmNth3trK = INITVALUE;
    CCTK_REAL PDstandardNth1W = INITVALUE;
    CCTK_REAL PDstandardNth2W = INITVALUE;
    CCTK_REAL PDstandardNth3W = INITVALUE;
    CCTK_REAL PDstandardNth11W = INITVALUE;
    CCTK_REAL PDstandardNth22W = INITVALUE;
    CCTK_REAL PDstandardNth33W = INITVALUE;
    CCTK_REAL PDstandardNth12W = INITVALUE;
    CCTK_REAL PDstandardNth13W = INITVALUE;
    CCTK_REAL PDstandardNth23W = INITVALUE;
    CCTK_REAL PDupwindpNth1W = INITVALUE;
    CCTK_REAL PDupwindpNth2W = INITVALUE;
    CCTK_REAL PDupwindpNth3W = INITVALUE;
    CCTK_REAL PDupwindmNth1W = INITVALUE;
    CCTK_REAL PDupwindmNth2W = INITVALUE;
    CCTK_REAL PDupwindmNth3W = INITVALUE;
    CCTK_REAL PDstandardNth1Xt1 = INITVALUE;
    CCTK_REAL PDstandardNth2Xt1 = INITVALUE;
    CCTK_REAL PDstandardNth3Xt1 = INITVALUE;
    CCTK_REAL PDupwindpNth1Xt1 = INITVALUE;
    CCTK_REAL PDupwindpNth2Xt1 = INITVALUE;
    CCTK_REAL PDupwindpNth3Xt1 = INITVALUE;
    CCTK_REAL PDupwindmNth1Xt1 = INITVALUE;
    CCTK_REAL PDupwindmNth2Xt1 = INITVALUE;
    CCTK_REAL PDupwindmNth3Xt1 = INITVALUE;
    CCTK_REAL PDstandardNth1Xt2 = INITVALUE;
    CCTK_REAL PDstandardNth2Xt2 = INITVALUE;
    CCTK_REAL PDstandardNth3Xt2 = INITVALUE;
    CCTK_REAL PDupwindpNth1Xt2 = INITVALUE;
    CCTK_REAL PDupwindpNth2Xt2 = INITVALUE;
    CCTK_REAL PDupwindpNth3Xt2 = INITVALUE;
    CCTK_REAL PDupwindmNth1Xt2 = INITVALUE;
    CCTK_REAL PDupwindmNth2Xt2 = INITVALUE;
    CCTK_REAL PDupwindmNth3Xt2 = INITVALUE;
    CCTK_REAL PDstandardNth1Xt3 = INITVALUE;
    CCTK_REAL PDstandardNth2Xt3 = INITVALUE;
    CCTK_REAL PDstandardNth3Xt3 = INITVALUE;
    CCTK_REAL PDupwindpNth1Xt3 = INITVALUE;
    CCTK_REAL PDupwindpNth2Xt3 = INITVALUE;
    CCTK_REAL PDupwindpNth3Xt3 = INITVALUE;
    CCTK_REAL PDupwindmNth1Xt3 = INITVALUE;
    CCTK_REAL PDupwindmNth2Xt3 = INITVALUE;
    CCTK_REAL PDupwindmNth3Xt3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    AL = A[index];
    alphaL = alpha[index];
    At11L = At11[index];
    At12L = At12[index];
    At13L = At13[index];
    At22L = At22[index];
    At23L = At23[index];
    At33L = At33[index];
    B1L = B1[index];
    B2L = B2[index];
    B3L = B3[index];
    beta1L = beta1[index];
    beta2L = beta2[index];
    beta3L = beta3[index];
    gt11L = gt11[index];
    gt12L = gt12[index];
    gt13L = gt13[index];
    gt22L = gt22[index];
    gt23L = gt23[index];
    gt33L = gt33[index];
    trKL = trK[index];
    trKrhsL = trKrhs[index];
    WL = W[index];
    Xt1L = Xt1[index];
    Xt1rhsL = Xt1rhs[index];
    Xt2L = Xt2[index];
    Xt2rhsL = Xt2rhs[index];
    Xt3L = Xt3[index];
    Xt3rhsL = Xt3rhs[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandardNth1alpha = PDstandardNth1(alpha, i, j, k);
    PDstandardNth2alpha = PDstandardNth2(alpha, i, j, k);
    PDstandardNth3alpha = PDstandardNth3(alpha, i, j, k);
    PDstandardNth11alpha = PDstandardNth11(alpha, i, j, k);
    PDstandardNth22alpha = PDstandardNth22(alpha, i, j, k);
    PDstandardNth33alpha = PDstandardNth33(alpha, i, j, k);
    PDstandardNth12alpha = PDstandardNth12(alpha, i, j, k);
    PDstandardNth13alpha = PDstandardNth13(alpha, i, j, k);
    PDstandardNth23alpha = PDstandardNth23(alpha, i, j, k);
    PDupwindpNth1alpha = PDupwindpNth1(alpha, i, j, k);
    PDupwindpNth2alpha = PDupwindpNth2(alpha, i, j, k);
    PDupwindpNth3alpha = PDupwindpNth3(alpha, i, j, k);
    PDupwindmNth1alpha = PDupwindmNth1(alpha, i, j, k);
    PDupwindmNth2alpha = PDupwindmNth2(alpha, i, j, k);
    PDupwindmNth3alpha = PDupwindmNth3(alpha, i, j, k);
    PDupwindpNth1At11 = PDupwindpNth1(At11, i, j, k);
    PDupwindpNth2At11 = PDupwindpNth2(At11, i, j, k);
    PDupwindpNth3At11 = PDupwindpNth3(At11, i, j, k);
    PDupwindmNth1At11 = PDupwindmNth1(At11, i, j, k);
    PDupwindmNth2At11 = PDupwindmNth2(At11, i, j, k);
    PDupwindmNth3At11 = PDupwindmNth3(At11, i, j, k);
    PDupwindpNth1At12 = PDupwindpNth1(At12, i, j, k);
    PDupwindpNth2At12 = PDupwindpNth2(At12, i, j, k);
    PDupwindpNth3At12 = PDupwindpNth3(At12, i, j, k);
    PDupwindmNth1At12 = PDupwindmNth1(At12, i, j, k);
    PDupwindmNth2At12 = PDupwindmNth2(At12, i, j, k);
    PDupwindmNth3At12 = PDupwindmNth3(At12, i, j, k);
    PDupwindpNth1At13 = PDupwindpNth1(At13, i, j, k);
    PDupwindpNth2At13 = PDupwindpNth2(At13, i, j, k);
    PDupwindpNth3At13 = PDupwindpNth3(At13, i, j, k);
    PDupwindmNth1At13 = PDupwindmNth1(At13, i, j, k);
    PDupwindmNth2At13 = PDupwindmNth2(At13, i, j, k);
    PDupwindmNth3At13 = PDupwindmNth3(At13, i, j, k);
    PDupwindpNth1At22 = PDupwindpNth1(At22, i, j, k);
    PDupwindpNth2At22 = PDupwindpNth2(At22, i, j, k);
    PDupwindpNth3At22 = PDupwindpNth3(At22, i, j, k);
    PDupwindmNth1At22 = PDupwindmNth1(At22, i, j, k);
    PDupwindmNth2At22 = PDupwindmNth2(At22, i, j, k);
    PDupwindmNth3At22 = PDupwindmNth3(At22, i, j, k);
    PDupwindpNth1At23 = PDupwindpNth1(At23, i, j, k);
    PDupwindpNth2At23 = PDupwindpNth2(At23, i, j, k);
    PDupwindpNth3At23 = PDupwindpNth3(At23, i, j, k);
    PDupwindmNth1At23 = PDupwindmNth1(At23, i, j, k);
    PDupwindmNth2At23 = PDupwindmNth2(At23, i, j, k);
    PDupwindmNth3At23 = PDupwindmNth3(At23, i, j, k);
    PDupwindpNth1At33 = PDupwindpNth1(At33, i, j, k);
    PDupwindpNth2At33 = PDupwindpNth2(At33, i, j, k);
    PDupwindpNth3At33 = PDupwindpNth3(At33, i, j, k);
    PDupwindmNth1At33 = PDupwindmNth1(At33, i, j, k);
    PDupwindmNth2At33 = PDupwindmNth2(At33, i, j, k);
    PDupwindmNth3At33 = PDupwindmNth3(At33, i, j, k);
    PDupwindpNth1B1 = PDupwindpNth1(B1, i, j, k);
    PDupwindpNth2B1 = PDupwindpNth2(B1, i, j, k);
    PDupwindpNth3B1 = PDupwindpNth3(B1, i, j, k);
    PDupwindmNth1B1 = PDupwindmNth1(B1, i, j, k);
    PDupwindmNth2B1 = PDupwindmNth2(B1, i, j, k);
    PDupwindmNth3B1 = PDupwindmNth3(B1, i, j, k);
    PDupwindpNth1B2 = PDupwindpNth1(B2, i, j, k);
    PDupwindpNth2B2 = PDupwindpNth2(B2, i, j, k);
    PDupwindpNth3B2 = PDupwindpNth3(B2, i, j, k);
    PDupwindmNth1B2 = PDupwindmNth1(B2, i, j, k);
    PDupwindmNth2B2 = PDupwindmNth2(B2, i, j, k);
    PDupwindmNth3B2 = PDupwindmNth3(B2, i, j, k);
    PDupwindpNth1B3 = PDupwindpNth1(B3, i, j, k);
    PDupwindpNth2B3 = PDupwindpNth2(B3, i, j, k);
    PDupwindpNth3B3 = PDupwindpNth3(B3, i, j, k);
    PDupwindmNth1B3 = PDupwindmNth1(B3, i, j, k);
    PDupwindmNth2B3 = PDupwindmNth2(B3, i, j, k);
    PDupwindmNth3B3 = PDupwindmNth3(B3, i, j, k);
    PDstandardNth1beta1 = PDstandardNth1(beta1, i, j, k);
    PDstandardNth2beta1 = PDstandardNth2(beta1, i, j, k);
    PDstandardNth3beta1 = PDstandardNth3(beta1, i, j, k);
    PDstandardNth11beta1 = PDstandardNth11(beta1, i, j, k);
    PDstandardNth22beta1 = PDstandardNth22(beta1, i, j, k);
    PDstandardNth33beta1 = PDstandardNth33(beta1, i, j, k);
    PDstandardNth12beta1 = PDstandardNth12(beta1, i, j, k);
    PDstandardNth13beta1 = PDstandardNth13(beta1, i, j, k);
    PDstandardNth23beta1 = PDstandardNth23(beta1, i, j, k);
    PDupwindpNth1beta1 = PDupwindpNth1(beta1, i, j, k);
    PDupwindpNth2beta1 = PDupwindpNth2(beta1, i, j, k);
    PDupwindpNth3beta1 = PDupwindpNth3(beta1, i, j, k);
    PDupwindmNth1beta1 = PDupwindmNth1(beta1, i, j, k);
    PDupwindmNth2beta1 = PDupwindmNth2(beta1, i, j, k);
    PDupwindmNth3beta1 = PDupwindmNth3(beta1, i, j, k);
    PDstandardNth1beta2 = PDstandardNth1(beta2, i, j, k);
    PDstandardNth2beta2 = PDstandardNth2(beta2, i, j, k);
    PDstandardNth3beta2 = PDstandardNth3(beta2, i, j, k);
    PDstandardNth11beta2 = PDstandardNth11(beta2, i, j, k);
    PDstandardNth22beta2 = PDstandardNth22(beta2, i, j, k);
    PDstandardNth33beta2 = PDstandardNth33(beta2, i, j, k);
    PDstandardNth12beta2 = PDstandardNth12(beta2, i, j, k);
    PDstandardNth13beta2 = PDstandardNth13(beta2, i, j, k);
    PDstandardNth23beta2 = PDstandardNth23(beta2, i, j, k);
    PDupwindpNth1beta2 = PDupwindpNth1(beta2, i, j, k);
    PDupwindpNth2beta2 = PDupwindpNth2(beta2, i, j, k);
    PDupwindpNth3beta2 = PDupwindpNth3(beta2, i, j, k);
    PDupwindmNth1beta2 = PDupwindmNth1(beta2, i, j, k);
    PDupwindmNth2beta2 = PDupwindmNth2(beta2, i, j, k);
    PDupwindmNth3beta2 = PDupwindmNth3(beta2, i, j, k);
    PDstandardNth1beta3 = PDstandardNth1(beta3, i, j, k);
    PDstandardNth2beta3 = PDstandardNth2(beta3, i, j, k);
    PDstandardNth3beta3 = PDstandardNth3(beta3, i, j, k);
    PDstandardNth11beta3 = PDstandardNth11(beta3, i, j, k);
    PDstandardNth22beta3 = PDstandardNth22(beta3, i, j, k);
    PDstandardNth33beta3 = PDstandardNth33(beta3, i, j, k);
    PDstandardNth12beta3 = PDstandardNth12(beta3, i, j, k);
    PDstandardNth13beta3 = PDstandardNth13(beta3, i, j, k);
    PDstandardNth23beta3 = PDstandardNth23(beta3, i, j, k);
    PDupwindpNth1beta3 = PDupwindpNth1(beta3, i, j, k);
    PDupwindpNth2beta3 = PDupwindpNth2(beta3, i, j, k);
    PDupwindpNth3beta3 = PDupwindpNth3(beta3, i, j, k);
    PDupwindmNth1beta3 = PDupwindmNth1(beta3, i, j, k);
    PDupwindmNth2beta3 = PDupwindmNth2(beta3, i, j, k);
    PDupwindmNth3beta3 = PDupwindmNth3(beta3, i, j, k);
    PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    PDstandardNth11gt11 = PDstandardNth11(gt11, i, j, k);
    PDstandardNth22gt11 = PDstandardNth22(gt11, i, j, k);
    PDstandardNth33gt11 = PDstandardNth33(gt11, i, j, k);
    PDstandardNth12gt11 = PDstandardNth12(gt11, i, j, k);
    PDstandardNth13gt11 = PDstandardNth13(gt11, i, j, k);
    PDstandardNth23gt11 = PDstandardNth23(gt11, i, j, k);
    PDupwindpNth1gt11 = PDupwindpNth1(gt11, i, j, k);
    PDupwindpNth2gt11 = PDupwindpNth2(gt11, i, j, k);
    PDupwindpNth3gt11 = PDupwindpNth3(gt11, i, j, k);
    PDupwindmNth1gt11 = PDupwindmNth1(gt11, i, j, k);
    PDupwindmNth2gt11 = PDupwindmNth2(gt11, i, j, k);
    PDupwindmNth3gt11 = PDupwindmNth3(gt11, i, j, k);
    PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    PDstandardNth11gt12 = PDstandardNth11(gt12, i, j, k);
    PDstandardNth22gt12 = PDstandardNth22(gt12, i, j, k);
    PDstandardNth33gt12 = PDstandardNth33(gt12, i, j, k);
    PDstandardNth12gt12 = PDstandardNth12(gt12, i, j, k);
    PDstandardNth13gt12 = PDstandardNth13(gt12, i, j, k);
    PDstandardNth23gt12 = PDstandardNth23(gt12, i, j, k);
    PDupwindpNth1gt12 = PDupwindpNth1(gt12, i, j, k);
    PDupwindpNth2gt12 = PDupwindpNth2(gt12, i, j, k);
    PDupwindpNth3gt12 = PDupwindpNth3(gt12, i, j, k);
    PDupwindmNth1gt12 = PDupwindmNth1(gt12, i, j, k);
    PDupwindmNth2gt12 = PDupwindmNth2(gt12, i, j, k);
    PDupwindmNth3gt12 = PDupwindmNth3(gt12, i, j, k);
    PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    PDstandardNth11gt13 = PDstandardNth11(gt13, i, j, k);
    PDstandardNth22gt13 = PDstandardNth22(gt13, i, j, k);
    PDstandardNth33gt13 = PDstandardNth33(gt13, i, j, k);
    PDstandardNth12gt13 = PDstandardNth12(gt13, i, j, k);
    PDstandardNth13gt13 = PDstandardNth13(gt13, i, j, k);
    PDstandardNth23gt13 = PDstandardNth23(gt13, i, j, k);
    PDupwindpNth1gt13 = PDupwindpNth1(gt13, i, j, k);
    PDupwindpNth2gt13 = PDupwindpNth2(gt13, i, j, k);
    PDupwindpNth3gt13 = PDupwindpNth3(gt13, i, j, k);
    PDupwindmNth1gt13 = PDupwindmNth1(gt13, i, j, k);
    PDupwindmNth2gt13 = PDupwindmNth2(gt13, i, j, k);
    PDupwindmNth3gt13 = PDupwindmNth3(gt13, i, j, k);
    PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    PDstandardNth11gt22 = PDstandardNth11(gt22, i, j, k);
    PDstandardNth22gt22 = PDstandardNth22(gt22, i, j, k);
    PDstandardNth33gt22 = PDstandardNth33(gt22, i, j, k);
    PDstandardNth12gt22 = PDstandardNth12(gt22, i, j, k);
    PDstandardNth13gt22 = PDstandardNth13(gt22, i, j, k);
    PDstandardNth23gt22 = PDstandardNth23(gt22, i, j, k);
    PDupwindpNth1gt22 = PDupwindpNth1(gt22, i, j, k);
    PDupwindpNth2gt22 = PDupwindpNth2(gt22, i, j, k);
    PDupwindpNth3gt22 = PDupwindpNth3(gt22, i, j, k);
    PDupwindmNth1gt22 = PDupwindmNth1(gt22, i, j, k);
    PDupwindmNth2gt22 = PDupwindmNth2(gt22, i, j, k);
    PDupwindmNth3gt22 = PDupwindmNth3(gt22, i, j, k);
    PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    PDstandardNth11gt23 = PDstandardNth11(gt23, i, j, k);
    PDstandardNth22gt23 = PDstandardNth22(gt23, i, j, k);
    PDstandardNth33gt23 = PDstandardNth33(gt23, i, j, k);
    PDstandardNth12gt23 = PDstandardNth12(gt23, i, j, k);
    PDstandardNth13gt23 = PDstandardNth13(gt23, i, j, k);
    PDstandardNth23gt23 = PDstandardNth23(gt23, i, j, k);
    PDupwindpNth1gt23 = PDupwindpNth1(gt23, i, j, k);
    PDupwindpNth2gt23 = PDupwindpNth2(gt23, i, j, k);
    PDupwindpNth3gt23 = PDupwindpNth3(gt23, i, j, k);
    PDupwindmNth1gt23 = PDupwindmNth1(gt23, i, j, k);
    PDupwindmNth2gt23 = PDupwindmNth2(gt23, i, j, k);
    PDupwindmNth3gt23 = PDupwindmNth3(gt23, i, j, k);
    PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    PDstandardNth11gt33 = PDstandardNth11(gt33, i, j, k);
    PDstandardNth22gt33 = PDstandardNth22(gt33, i, j, k);
    PDstandardNth33gt33 = PDstandardNth33(gt33, i, j, k);
    PDstandardNth12gt33 = PDstandardNth12(gt33, i, j, k);
    PDstandardNth13gt33 = PDstandardNth13(gt33, i, j, k);
    PDstandardNth23gt33 = PDstandardNth23(gt33, i, j, k);
    PDupwindpNth1gt33 = PDupwindpNth1(gt33, i, j, k);
    PDupwindpNth2gt33 = PDupwindpNth2(gt33, i, j, k);
    PDupwindpNth3gt33 = PDupwindpNth3(gt33, i, j, k);
    PDupwindmNth1gt33 = PDupwindmNth1(gt33, i, j, k);
    PDupwindmNth2gt33 = PDupwindmNth2(gt33, i, j, k);
    PDupwindmNth3gt33 = PDupwindmNth3(gt33, i, j, k);
    PDstandardNth1trK = PDstandardNth1(trK, i, j, k);
    PDstandardNth2trK = PDstandardNth2(trK, i, j, k);
    PDstandardNth3trK = PDstandardNth3(trK, i, j, k);
    PDupwindpNth1trK = PDupwindpNth1(trK, i, j, k);
    PDupwindpNth2trK = PDupwindpNth2(trK, i, j, k);
    PDupwindpNth3trK = PDupwindpNth3(trK, i, j, k);
    PDupwindmNth1trK = PDupwindmNth1(trK, i, j, k);
    PDupwindmNth2trK = PDupwindmNth2(trK, i, j, k);
    PDupwindmNth3trK = PDupwindmNth3(trK, i, j, k);
    PDstandardNth1W = PDstandardNth1(W, i, j, k);
    PDstandardNth2W = PDstandardNth2(W, i, j, k);
    PDstandardNth3W = PDstandardNth3(W, i, j, k);
    PDstandardNth11W = PDstandardNth11(W, i, j, k);
    PDstandardNth22W = PDstandardNth22(W, i, j, k);
    PDstandardNth33W = PDstandardNth33(W, i, j, k);
    PDstandardNth12W = PDstandardNth12(W, i, j, k);
    PDstandardNth13W = PDstandardNth13(W, i, j, k);
    PDstandardNth23W = PDstandardNth23(W, i, j, k);
    PDupwindpNth1W = PDupwindpNth1(W, i, j, k);
    PDupwindpNth2W = PDupwindpNth2(W, i, j, k);
    PDupwindpNth3W = PDupwindpNth3(W, i, j, k);
    PDupwindmNth1W = PDupwindmNth1(W, i, j, k);
    PDupwindmNth2W = PDupwindmNth2(W, i, j, k);
    PDupwindmNth3W = PDupwindmNth3(W, i, j, k);
    PDstandardNth1Xt1 = PDstandardNth1(Xt1, i, j, k);
    PDstandardNth2Xt1 = PDstandardNth2(Xt1, i, j, k);
    PDstandardNth3Xt1 = PDstandardNth3(Xt1, i, j, k);
    PDupwindpNth1Xt1 = PDupwindpNth1(Xt1, i, j, k);
    PDupwindpNth2Xt1 = PDupwindpNth2(Xt1, i, j, k);
    PDupwindpNth3Xt1 = PDupwindpNth3(Xt1, i, j, k);
    PDupwindmNth1Xt1 = PDupwindmNth1(Xt1, i, j, k);
    PDupwindmNth2Xt1 = PDupwindmNth2(Xt1, i, j, k);
    PDupwindmNth3Xt1 = PDupwindmNth3(Xt1, i, j, k);
    PDstandardNth1Xt2 = PDstandardNth1(Xt2, i, j, k);
    PDstandardNth2Xt2 = PDstandardNth2(Xt2, i, j, k);
    PDstandardNth3Xt2 = PDstandardNth3(Xt2, i, j, k);
    PDupwindpNth1Xt2 = PDupwindpNth1(Xt2, i, j, k);
    PDupwindpNth2Xt2 = PDupwindpNth2(Xt2, i, j, k);
    PDupwindpNth3Xt2 = PDupwindpNth3(Xt2, i, j, k);
    PDupwindmNth1Xt2 = PDupwindmNth1(Xt2, i, j, k);
    PDupwindmNth2Xt2 = PDupwindmNth2(Xt2, i, j, k);
    PDupwindmNth3Xt2 = PDupwindmNth3(Xt2, i, j, k);
    PDstandardNth1Xt3 = PDstandardNth1(Xt3, i, j, k);
    PDstandardNth2Xt3 = PDstandardNth2(Xt3, i, j, k);
    PDstandardNth3Xt3 = PDstandardNth3(Xt3, i, j, k);
    PDupwindpNth1Xt3 = PDupwindpNth1(Xt3, i, j, k);
    PDupwindpNth2Xt3 = PDupwindpNth2(Xt3, i, j, k);
    PDupwindpNth3Xt3 = PDupwindpNth3(Xt3, i, j, k);
    PDupwindmNth1Xt3 = PDupwindmNth1(Xt3, i, j, k);
    PDupwindmNth2Xt3 = PDupwindmNth2(Xt3, i, j, k);
    PDupwindmNth3Xt3 = PDupwindmNth3(Xt3, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL const T1000001  =  Abs(beta1L);
    
    CCTK_REAL const T1000002  =  Abs(beta2L);
    
    CCTK_REAL const T1000003  =  Abs(beta3L);
    
    CCTK_REAL const T1000005  =  -PDstandardNth3gt12;
    
    CCTK_REAL const T1000007  =  -PDstandardNth2gt13;
    
    CCTK_REAL const T1000009  =  -PDstandardNth1gt23;
    
    CCTK_REAL const T1000171  =  
       PDstandardNth1beta1 + PDstandardNth2beta2 + PDstandardNth3beta3;
    
    CCTK_REAL const T1000172  =  
       PDstandardNth11beta1 + PDstandardNth12beta2 + PDstandardNth13beta3;
    
    CCTK_REAL const T1000173  =  
       PDstandardNth12beta1 + PDstandardNth22beta2 + PDstandardNth23beta3;
    
    CCTK_REAL const T1000174  =  
       PDstandardNth13beta1 + PDstandardNth23beta2 + PDstandardNth33beta3;
    
    CCTK_REAL const T1000175  =  -1 + LapseAdvectionCoeff;
    
    detgt  =  1;
    
    invW  =  INV(WL);
    
    W2  =  SQR(WL);
    
    CCTK_REAL const T1000004  =  INV(detgt);
    
    CCTK_REAL const T1000006  =  PDstandardNth1gt23 + PDstandardNth2gt13 + T1000005;
    
    CCTK_REAL const T1000008  =  PDstandardNth1gt23 + PDstandardNth3gt12 + T1000007;
    
    CCTK_REAL const T1000010  =  PDstandardNth2gt13 + PDstandardNth3gt12 + T1000009;
    
    CCTK_REAL const T1000149  =  invW*PDstandardNth3W;
    
    betam1  =  khalf*(beta1L - T1000001);
    
    betam2  =  khalf*(beta2L - T1000002);
    
    betam3  =  khalf*(beta3L - T1000003);
    
    betap1  =  khalf*(beta1L + T1000001);
    
    betap2  =  khalf*(beta2L + T1000002);
    
    betap3  =  khalf*(beta3L + T1000003);
    
    invW2  =  SQR(invW);
    
    pdphi1  =  -(invW*khalf*PDstandardNth1W);
    
    pdphi2  =  -(invW*khalf*PDstandardNth2W);
    
    pdphi3  =  -(invW*khalf*PDstandardNth3W);
    
    CCTK_REAL const T1000150  =  SQR(pdphi1);
    
    CCTK_REAL const T1000151  =  SQR(pdphi2);
    
    CCTK_REAL const T1000152  =  SQR(pdphi3);
    
    CCTK_REAL const T1000168  =  -pdphi1;
    
    CCTK_REAL const T1000169  =  -pdphi3;
    
    CCTK_REAL const T1000170  =  -pdphi2;
    
    gtu11  =  T1000004*(gt22L*gt33L - SQR(gt23L));
    
    gtu21  =  (gt13L*gt23L - gt12L*gt33L)*T1000004;
    
    gtu31  =  (-(gt13L*gt22L) + gt12L*gt23L)*T1000004;
    
    gtu22  =  T1000004*(gt11L*gt33L - SQR(gt13L));
    
    gtu32  =  (gt12L*gt13L - gt11L*gt23L)*T1000004;
    
    gtu33  =  T1000004*(gt11L*gt22L - SQR(gt12L));
    
    g11  =  gt11L*invW2;
    
    g12  =  gt12L*invW2;
    
    g13  =  gt13L*invW2;
    
    g22  =  gt22L*invW2;
    
    g23  =  gt23L*invW2;
    
    g33  =  gt33L*invW2;
    
    WrhsL  =  betam1*PDupwindmNth1W + betam2*PDupwindmNth2W + 
        betam3*PDupwindmNth3W + betap1*PDupwindpNth1W + betap2*PDupwindpNth2W + 
        betap3*PDupwindpNth3W + kthird*(-T1000171 + alphaL*trKL)*WL;
    
    gt11rhsL  =  -2*alphaL*At11L + 2*(gt11L*PDstandardNth1beta1 + 
           gt12L*PDstandardNth1beta2 + gt13L*PDstandardNth1beta3) + 
        betam1*PDupwindmNth1gt11 + betam2*PDupwindmNth2gt11 + 
        betam3*PDupwindmNth3gt11 + betap1*PDupwindpNth1gt11 + 
        betap2*PDupwindpNth2gt11 + betap3*PDupwindpNth3gt11 - 
        gt11L*ktwothird*T1000171;
    
    gt12rhsL  =  -2*alphaL*At12L + gt22L*PDstandardNth1beta2 + 
        gt23L*PDstandardNth1beta3 + gt11L*PDstandardNth2beta1 + 
        gt13L*PDstandardNth2beta3 + gt12L*
         (kthird*(PDstandardNth1beta1 + PDstandardNth2beta2) - 
           ktwothird*PDstandardNth3beta3) + betam1*PDupwindmNth1gt12 + 
        betam2*PDupwindmNth2gt12 + betam3*PDupwindmNth3gt12 + 
        betap1*PDupwindpNth1gt12 + betap2*PDupwindpNth2gt12 + 
        betap3*PDupwindpNth3gt12;
    
    gt13rhsL  =  -2*alphaL*At13L + gt23L*PDstandardNth1beta2 + 
        gt33L*PDstandardNth1beta3 + gt11L*PDstandardNth3beta1 + 
        gt12L*PDstandardNth3beta2 + gt13L*
         (-(ktwothird*PDstandardNth2beta2) + 
           kthird*(PDstandardNth1beta1 + PDstandardNth3beta3)) + 
        betam1*PDupwindmNth1gt13 + betam2*PDupwindmNth2gt13 + 
        betam3*PDupwindmNth3gt13 + betap1*PDupwindpNth1gt13 + 
        betap2*PDupwindpNth2gt13 + betap3*PDupwindpNth3gt13;
    
    gt22rhsL  =  -2*alphaL*At22L + 2*(gt12L*PDstandardNth2beta1 + 
           gt22L*PDstandardNth2beta2 + gt23L*PDstandardNth2beta3) + 
        betam1*PDupwindmNth1gt22 + betam2*PDupwindmNth2gt22 + 
        betam3*PDupwindmNth3gt22 + betap1*PDupwindpNth1gt22 + 
        betap2*PDupwindpNth2gt22 + betap3*PDupwindpNth3gt22 - 
        gt22L*ktwothird*T1000171;
    
    gt23rhsL  =  -2*alphaL*At23L + gt13L*PDstandardNth2beta1 + 
        gt33L*PDstandardNth2beta3 + gt12L*PDstandardNth3beta1 + 
        gt22L*PDstandardNth3beta2 + gt23L*
         (-(ktwothird*PDstandardNth1beta1) + 
           kthird*(PDstandardNth2beta2 + PDstandardNth3beta3)) + 
        betam1*PDupwindmNth1gt23 + betam2*PDupwindmNth2gt23 + 
        betam3*PDupwindmNth3gt23 + betap1*PDupwindpNth1gt23 + 
        betap2*PDupwindpNth2gt23 + betap3*PDupwindpNth3gt23;
    
    gt33rhsL  =  -2*alphaL*At33L + 2*(gt13L*PDstandardNth3beta1 + 
           gt23L*PDstandardNth3beta2 + gt33L*PDstandardNth3beta3) + 
        betam1*PDupwindmNth1gt33 + betam2*PDupwindmNth2gt33 + 
        betam3*PDupwindmNth3gt33 + betap1*PDupwindpNth1gt33 + 
        betap2*PDupwindpNth2gt33 + betap3*PDupwindpNth3gt33 - 
        gt33L*ktwothird*T1000171;
    
    alpharhsL  =  LapseAdvectionCoeff*
         (betam1*PDupwindmNth1alpha + betam2*PDupwindmNth2alpha + 
           betam3*PDupwindmNth3alpha + betap1*PDupwindpNth1alpha + 
           betap2*PDupwindpNth2alpha + betap3*PDupwindpNth3alpha) + 
        harmonicF*(AL*T1000175 - LapseAdvectionCoeff*trKL)*pow(alphaL,harmonicN);
    
    beta1rhsL  =  (betam1*PDupwindmNth1beta1 + betam2*PDupwindmNth2beta1 + 
           betam3*PDupwindmNth3beta1 + betap1*PDupwindpNth1beta1 + 
           betap2*PDupwindpNth2beta1 + betap3*PDupwindpNth3beta1)*
         ShiftAdvectionCoeff + B1L*ShiftGammaCoeff;
    
    beta2rhsL  =  (betam1*PDupwindmNth1beta2 + betam2*PDupwindmNth2beta2 + 
           betam3*PDupwindmNth3beta2 + betap1*PDupwindpNth1beta2 + 
           betap2*PDupwindpNth2beta2 + betap3*PDupwindpNth3beta2)*
         ShiftAdvectionCoeff + B2L*ShiftGammaCoeff;
    
    beta3rhsL  =  (betam1*PDupwindmNth1beta3 + betam2*PDupwindmNth2beta3 + 
           betam3*PDupwindmNth3beta3 + betap1*PDupwindpNth1beta3 + 
           betap2*PDupwindpNth2beta3 + betap3*PDupwindpNth3beta3)*
         ShiftAdvectionCoeff + B3L*ShiftGammaCoeff;
    
    CCTK_REAL const T1000153  =  At12L*gtu21;
    
    CCTK_REAL const T1000154  =  At13L*gtu31;
    
    CCTK_REAL const T1000155  =  At23L*gtu32;
    
    CCTK_REAL const T1000156  =  gtu31*pdphi1;
    
    CCTK_REAL const T1000157  =  gtu32*pdphi2;
    
    CCTK_REAL const T1000158  =  gtu33*pdphi3;
    
    CCTK_REAL const T1000160  =  gtu11*pdphi1;
    
    CCTK_REAL const T1000161  =  gtu21*pdphi2;
    
    CCTK_REAL const T1000162  =  gtu31*pdphi3;
    
    CCTK_REAL const T1000164  =  gtu21*pdphi1;
    
    CCTK_REAL const T1000165  =  gtu22*pdphi2;
    
    CCTK_REAL const T1000166  =  gtu32*pdphi3;
    
    Gt111  =  khalf*(gtu11*PDstandardNth1gt11 + 
          2*(gtu21*PDstandardNth1gt12 + gtu31*PDstandardNth1gt13) - 
          gtu21*PDstandardNth2gt11 - gtu31*PDstandardNth3gt11);
    
    Gt211  =  khalf*(gtu21*PDstandardNth1gt11 + 
          2*(gtu22*PDstandardNth1gt12 + gtu32*PDstandardNth1gt13) - 
          gtu22*PDstandardNth2gt11 - gtu32*PDstandardNth3gt11);
    
    Gt311  =  khalf*(gtu31*PDstandardNth1gt11 + 
          2*(gtu32*PDstandardNth1gt12 + gtu33*PDstandardNth1gt13) - 
          gtu32*PDstandardNth2gt11 - gtu33*PDstandardNth3gt11);
    
    Gt112  =  khalf*(gtu21*PDstandardNth1gt22 + gtu11*PDstandardNth2gt11 + 
          gtu31*T1000006);
    
    Gt212  =  khalf*(gtu22*PDstandardNth1gt22 + gtu21*PDstandardNth2gt11 + 
          gtu32*T1000006);
    
    Gt312  =  khalf*(gtu32*PDstandardNth1gt22 + gtu31*PDstandardNth2gt11 + 
          gtu33*T1000006);
    
    Gt113  =  khalf*(gtu31*PDstandardNth1gt33 + gtu11*PDstandardNth3gt11 + 
          gtu21*T1000008);
    
    Gt213  =  khalf*(gtu32*PDstandardNth1gt33 + gtu21*PDstandardNth3gt11 + 
          gtu22*T1000008);
    
    Gt313  =  khalf*(gtu33*PDstandardNth1gt33 + gtu31*PDstandardNth3gt11 + 
          gtu32*T1000008);
    
    Gt122  =  khalf*(gtu11*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + 
          gtu21*PDstandardNth2gt22 + 
          gtu31*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    Gt222  =  khalf*(gtu21*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + 
          gtu22*PDstandardNth2gt22 + 
          gtu32*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    Gt322  =  khalf*(gtu31*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + 
          gtu32*PDstandardNth2gt22 + 
          gtu33*(2*PDstandardNth2gt23 - PDstandardNth3gt22));
    
    Gt123  =  khalf*(gtu31*PDstandardNth2gt33 + gtu21*PDstandardNth3gt22 + 
          gtu11*T1000010);
    
    Gt223  =  khalf*(gtu32*PDstandardNth2gt33 + gtu22*PDstandardNth3gt22 + 
          gtu21*T1000010);
    
    Gt323  =  khalf*(gtu33*PDstandardNth2gt33 + gtu32*PDstandardNth3gt22 + 
          gtu31*T1000010);
    
    Gt133  =  khalf*(-(gtu11*PDstandardNth1gt33) - gtu21*PDstandardNth2gt33 + 
          2*gtu11*PDstandardNth3gt13 + 2*gtu21*PDstandardNth3gt23 + 
          gtu31*PDstandardNth3gt33);
    
    Gt233  =  khalf*(-(gtu21*PDstandardNth1gt33) - gtu22*PDstandardNth2gt33 + 
          2*gtu21*PDstandardNth3gt13 + 2*gtu22*PDstandardNth3gt23 + 
          gtu32*PDstandardNth3gt33);
    
    Gt333  =  khalf*(-(gtu31*PDstandardNth1gt33) - gtu32*PDstandardNth2gt33 + 
          2*gtu31*PDstandardNth3gt13 + 2*gtu32*PDstandardNth3gt23 + 
          gtu33*PDstandardNth3gt33);
    
    Atm21  =  At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    Atm31  =  At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    Atm12  =  At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    Atm32  =  At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    Atm13  =  At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    Atm23  =  At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    gu11  =  gtu11*W2;
    
    gu21  =  gtu21*W2;
    
    gu31  =  gtu31*W2;
    
    gu22  =  gtu22*W2;
    
    gu32  =  gtu32*W2;
    
    gu33  =  gtu33*W2;
    
    CCTK_REAL const T1000011  =  SQR(Gt212);
    
    CCTK_REAL const T1000012  =  SQR(Gt313);
    
    CCTK_REAL const T1000013  =  SQR(Gt111);
    
    CCTK_REAL const T1000014  =  SQR(Gt112);
    
    CCTK_REAL const T1000015  =  SQR(Gt312);
    
    CCTK_REAL const T1000017  =  Gt211*Gt222*gt22L;
    
    CCTK_REAL const T1000018  =  Gt211*gt23L*Gt322;
    
    CCTK_REAL const T1000019  =  Gt112*gt12L;
    
    CCTK_REAL const T1000020  =  Gt112*Gt212;
    
    CCTK_REAL const T1000021  =  Gt113*Gt312;
    
    CCTK_REAL const T1000022  =  2*gt11L*Gt122;
    
    CCTK_REAL const T1000023  =  gt12L*Gt222;
    
    CCTK_REAL const T1000024  =  Gt212*gt22L;
    
    CCTK_REAL const T1000025  =  gt23L*Gt312;
    
    CCTK_REAL const T1000026  =  gt13L*Gt322;
    
    CCTK_REAL const T1000027  =  Gt123*Gt312;
    
    CCTK_REAL const T1000028  =  Gt212*Gt213*gt22L;
    
    CCTK_REAL const T1000029  =  Gt211*Gt223*gt22L;
    
    CCTK_REAL const T1000030  =  Gt213*gt23L*Gt312;
    
    CCTK_REAL const T1000031  =  Gt211*gt23L*Gt323;
    
    CCTK_REAL const T1000032  =  Gt113*gt12L;
    
    CCTK_REAL const T1000033  =  Gt213*gt22L;
    
    CCTK_REAL const T1000034  =  gt23L*Gt313;
    
    CCTK_REAL const T1000036  =  Gt213*Gt222*gt22L;
    
    CCTK_REAL const T1000037  =  Gt212*Gt223*gt22L;
    
    CCTK_REAL const T1000038  =  Gt213*gt23L*Gt322;
    
    CCTK_REAL const T1000039  =  Gt212*gt23L*Gt323;
    
    CCTK_REAL const T1000040  =  2*gt11L*Gt123;
    
    CCTK_REAL const T1000041  =  gt12L*Gt223;
    
    CCTK_REAL const T1000042  =  gt13L*Gt323;
    
    CCTK_REAL const T1000044  =  SQR(Gt113);
    
    CCTK_REAL const T1000045  =  Gt211*gt22L;
    
    CCTK_REAL const T1000046  =  gt23L*Gt311;
    
    CCTK_REAL const T1000047  =  gt12L*Gt212;
    
    CCTK_REAL const T1000048  =  gt13L*Gt312;
    
    CCTK_REAL const T1000050  =  Gt112*Gt213;
    
    CCTK_REAL const T1000051  =  Gt113*Gt313;
    
    CCTK_REAL const T1000052  =  gt12L*Gt212*Gt223;
    
    CCTK_REAL const T1000053  =  2*Gt123*gt13L*Gt311;
    
    CCTK_REAL const T1000054  =  Gt223*gt23L*Gt311;
    
    CCTK_REAL const T1000055  =  Gt113*gt13L*Gt312;
    
    CCTK_REAL const T1000056  =  gt13L*Gt223*Gt312;
    
    CCTK_REAL const T1000057  =  gt13L*Gt313*Gt323;
    
    CCTK_REAL const T1000058  =  Gt312*Gt313*gt33L;
    
    CCTK_REAL const T1000059  =  Gt311*Gt323*gt33L;
    
    CCTK_REAL const T1000060  =  Gt212*gt23L*Gt313;
    
    CCTK_REAL const T1000061  =  2*Gt112*Gt113;
    
    CCTK_REAL const T1000062  =  Gt122*Gt213;
    
    CCTK_REAL const T1000063  =  Gt112*gt13L;
    
    CCTK_REAL const T1000064  =  gt12L*Gt222*Gt223;
    
    CCTK_REAL const T1000065  =  gt13L*Gt223*Gt322;
    
    CCTK_REAL const T1000066  =  SQR(Gt323);
    
    CCTK_REAL const T1000068  =  Gt212*gt23L;
    
    CCTK_REAL const T1000069  =  Gt312*gt33L;
    
    CCTK_REAL const T1000070  =  SQR(Gt213);
    
    CCTK_REAL const T1000071  =  Gt233*gt23L*Gt311;
    
    CCTK_REAL const T1000072  =  Gt133*Gt313;
    
    CCTK_REAL const T1000073  =  gt13L*Gt313*Gt333;
    
    CCTK_REAL const T1000075  =  Gt311*Gt333*gt33L;
    
    CCTK_REAL const T1000076  =  Gt113*gt13L;
    
    CCTK_REAL const T1000077  =  Gt223*gt23L*Gt313;
    
    CCTK_REAL const T1000078  =  gt13L*Gt323*Gt333;
    
    CCTK_REAL const T1000079  =  2*gt11L*Gt133;
    
    CCTK_REAL const T1000080  =  gt12L*Gt233;
    
    CCTK_REAL const T1000081  =  gt13L*Gt333;
    
    CCTK_REAL const T1000083  =  Gt313*Gt323*gt33L;
    
    CCTK_REAL const T1000084  =  SQR(Gt223);
    
    CCTK_REAL const T1000085  =  Gt233*gt23L*Gt312;
    
    CCTK_REAL const T1000086  =  2*Gt123*gt13L*Gt313;
    
    CCTK_REAL const T1000087  =  gt13L*Gt223*Gt323;
    
    CCTK_REAL const T1000088  =  Gt312*Gt333*gt33L;
    
    CCTK_REAL const T1000089  =  Gt213*gt23L;
    
    CCTK_REAL const T1000090  =  Gt313*gt33L;
    
    CCTK_REAL const T1000092  =  gt12L*Gt213;
    
    CCTK_REAL const T1000093  =  gt13L*Gt313;
    
    CCTK_REAL const T1000094  =  Gt211*gt23L;
    
    CCTK_REAL const T1000095  =  Gt311*gt33L;
    
    CCTK_REAL const T1000096  =  gt11L*Gt123;
    
    CCTK_REAL const T1000099  =  SQR(Gt222);
    
    CCTK_REAL const T1000100  =  Gt113*gt12L*Gt212;
    
    CCTK_REAL const T1000101  =  Gt112*gt13L*Gt212;
    
    CCTK_REAL const T1000102  =  2*Gt112*gt12L*Gt213;
    
    CCTK_REAL const T1000103  =  Gt113*Gt211*gt22L;
    
    CCTK_REAL const T1000104  =  Gt112*Gt211*gt23L;
    
    CCTK_REAL const T1000106  =  Gt113*gt23L*Gt311;
    
    CCTK_REAL const T1000107  =  2*Gt113*gt13L*Gt312;
    
    CCTK_REAL const T1000108  =  Gt112*gt13L*Gt313;
    
    CCTK_REAL const T1000109  =  Gt213*gt22L*Gt313;
    
    CCTK_REAL const T1000111  =  Gt112*Gt311*gt33L;
    
    CCTK_REAL const T1000112  =  Gt212*Gt312*gt33L;
    
    CCTK_REAL const T1000113  =  Gt123*gt12L*Gt212;
    
    CCTK_REAL const T1000114  =  Gt212*Gt222*gt23L;
    
    CCTK_REAL const T1000115  =  2*Gt223*gt23L*Gt312;
    
    CCTK_REAL const T1000116  =  Gt113*gt13L*Gt322;
    
    CCTK_REAL const T1000117  =  gt23L*Gt313*Gt323;
    
    CCTK_REAL const T1000118  =  Gt313*Gt322*gt33L;
    
    CCTK_REAL const T1000119  =  Gt312*Gt323*gt33L;
    
    CCTK_REAL const T1000121  =  Gt122*gt13L*Gt212;
    
    CCTK_REAL const T1000122  =  2*Gt122*gt12L*Gt213;
    
    CCTK_REAL const T1000123  =  Gt113*gt12L*Gt222;
    
    CCTK_REAL const T1000124  =  Gt113*Gt212*gt22L;
    
    CCTK_REAL const T1000125  =  Gt123*gt13L*Gt312;
    
    CCTK_REAL const T1000126  =  Gt113*gt23L*Gt312;
    
    CCTK_REAL const T1000127  =  Gt223*gt23L*Gt312;
    
    CCTK_REAL const T1000128  =  Gt223*gt22L*Gt313;
    
    CCTK_REAL const T1000129  =  Gt222*gt23L*Gt313;
    
    CCTK_REAL const T1000130  =  Gt212*Gt322*gt33L;
    
    CCTK_REAL const T1000131  =  Gt212*Gt223*gt23L;
    
    CCTK_REAL const T1000132  =  Gt133*gt13L*Gt312;
    
    CCTK_REAL const T1000133  =  Gt113*gt13L*Gt323;
    
    CCTK_REAL const T1000134  =  Gt213*gt23L*Gt323;
    
    CCTK_REAL const T1000135  =  Gt212*gt23L*Gt333;
    
    CCTK_REAL const T1000136  =  gt23L*Gt313*Gt333;
    
    CCTK_REAL const T1000138  =  2*Gt213*Gt223*gt22L;
    
    CCTK_REAL const T1000139  =  gt22L*Gt233*Gt313;
    
    CCTK_REAL const T1000140  =  2*Gt213*gt23L*Gt323;
    
    CCTK_REAL const T1000141  =  Gt212*Gt323*gt33L;
    
    CCTK_REAL const T1000142  =  SQR(Gt123);
    
    CCTK_REAL const T1000143  =  Gt222*Gt223*gt23L;
    
    CCTK_REAL const T1000144  =  Gt133*gt13L*Gt322;
    
    CCTK_REAL const T1000145  =  gt23L*Gt323*Gt333;
    
    CCTK_REAL const T1000147  =  Gt322*Gt333*gt33L;
    
    CCTK_REAL const T1000148  =  SQR(Gt333);
    
    CCTK_REAL const T1000159  =  T1000156 + T1000157 + T1000158;
    
    CCTK_REAL const T1000163  =  T1000160 + T1000161 + T1000162;
    
    CCTK_REAL const T1000167  =  T1000164 + T1000165 + T1000166;
    
    Xtn1  =  Gt111*gtu11 + Gt122*gtu22 + 
        2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    Xtn2  =  Gt211*gtu11 + Gt222*gtu22 + 
        2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    Xtn3  =  Gt311*gtu11 + Gt322*gtu22 + 
        2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    cdphi211  =  invW*khalf*(-PDstandardNth11W + Gt111*PDstandardNth1W + 
          Gt211*PDstandardNth2W + Gt311*PDstandardNth3W + invW*SQR(PDstandardNth1W));
    
    cdphi212  =  invW*khalf*(-PDstandardNth12W + Gt212*PDstandardNth2W + 
          PDstandardNth1W*(Gt112 + invW*PDstandardNth2W) + Gt312*PDstandardNth3W);
    
    cdphi213  =  invW*khalf*(-PDstandardNth13W + Gt213*PDstandardNth2W + 
          Gt313*PDstandardNth3W + PDstandardNth1W*(Gt113 + T1000149));
    
    cdphi222  =  invW*khalf*(Gt122*PDstandardNth1W - PDstandardNth22W + 
          Gt222*PDstandardNth2W + Gt322*PDstandardNth3W + invW*SQR(PDstandardNth2W));
    
    cdphi223  =  invW*khalf*(Gt123*PDstandardNth1W - PDstandardNth23W + 
          Gt323*PDstandardNth3W + PDstandardNth2W*(Gt223 + T1000149));
    
    cdphi233  =  invW*khalf*(Gt133*PDstandardNth1W + Gt233*PDstandardNth2W - 
          PDstandardNth33W + Gt333*PDstandardNth3W + invW*SQR(PDstandardNth3W));
    
    Atm11  =  At11L*gtu11 + T1000153 + T1000154;
    
    Atm22  =  At22L*gtu22 + T1000153 + T1000155;
    
    Atm33  =  At33L*gtu33 + T1000154 + T1000155;
    
    G111  =  Gt111 - 2*((-2 + gt11L*gtu11)*pdphi1 + 
           gt11L*(gtu21*pdphi2 + gtu31*pdphi3));
    
    G212  =  Gt212 - 2*(gt12L*(gtu21*pdphi1 + gtu22*pdphi2 + gtu32*pdphi3) + 
           T1000168);
    
    G313  =  Gt313 - 2*(gt13L*(gtu31*pdphi1 + gtu32*pdphi2 + gtu33*pdphi3) + 
           T1000168);
    
    CCTK_REAL const T1000016  =  gt22L*T1000011;
    
    CCTK_REAL const T1000035  =  T1000032 + T1000033 + T1000034;
    
    CCTK_REAL const T1000043  =  T1000032 + T1000040 + T1000041 + T1000042;
    
    CCTK_REAL const T1000049  =  T1000019 + T1000024 + T1000025;
    
    CCTK_REAL const T1000067  =  gt13L*T1000066;
    
    CCTK_REAL const T1000074  =  gt33L*T1000012;
    
    CCTK_REAL const T1000082  =  T1000076 + T1000079 + T1000080 + T1000081;
    
    CCTK_REAL const T1000091  =  T1000076 + T1000089 + T1000090;
    
    CCTK_REAL const T1000097  =  T1000041 + T1000042 + T1000096;
    
    CCTK_REAL const T1000098  =  T1000063 + T1000068 + T1000069;
    
    CCTK_REAL const T1000105  =  gt23L*T1000011;
    
    CCTK_REAL const T1000110  =  gt23L*T1000012;
    
    CCTK_REAL const T1000120  =  gt13L*T1000014;
    
    CCTK_REAL const T1000137  =  gt12L*T1000044;
    
    CCTK_REAL const T1000146  =  gt33L*T1000066;
    
    Rt11  =  -(gtu11*khalf*PDstandardNth11gt11) + 
        gtu21*(2*Gt211*Gt212*gt22L + 4*Gt112*gt13L*Gt311 + 2*Gt113*gt11L*Gt312 + 
           2*gt13L*Gt312*Gt313 + 2*gt13L*Gt211*Gt322 + 2*gt13L*Gt311*Gt323 + 
           2*Gt311*Gt312*gt33L - PDstandardNth12gt11) - gtu31*PDstandardNth13gt11 + 
        gt11L*PDstandardNth1Xt1 + gt12L*
         (4*Gt111*Gt212*gtu21 + 2*Gt211*Gt222*gtu21 + 2*Gt212*Gt222*gtu22 + 
           4*Gt113*Gt211*gtu31 + 4*Gt113*Gt212*gtu32 + 4*Gt113*Gt213*gtu33 + 
           PDstandardNth1Xt2) + gt13L*
         (4*Gt111*Gt312*gtu21 + 2*Gt212*Gt312*gtu21 + 4*Gt112*Gt312*gtu22 + 
           4*Gt113*Gt311*gtu31 + 4*Gt113*Gt312*gtu32 + 4*Gt113*Gt313*gtu33 + 
           PDstandardNth1Xt3) - gtu22*khalf*PDstandardNth22gt11 - 
        gtu32*PDstandardNth23gt11 - gtu33*khalf*PDstandardNth33gt11 + 
        gt22L*gtu22*T1000011 + 2*(gt12L*Gt211*Gt212*gtu11 + 
           Gt113*gt11L*Gt311*gtu11 + Gt211*gt23L*Gt311*gtu11 + 
           gt13L*Gt211*Gt312*gtu11 + Gt112*gt11L*Gt212*gtu21 + 
           gt12L*Gt223*Gt311*gtu21 + Gt212*gt23L*Gt311*gtu21 + 
           gt12L*Gt213*Gt312*gtu21 + Gt211*gt23L*Gt312*gtu21 + 
           gt11L*Gt122*Gt212*gtu22 + gt11L*Gt123*Gt312*gtu22 + 
           gt12L*Gt223*Gt312*gtu22 + Gt212*gt23L*Gt312*gtu22 + 
           gt13L*Gt212*Gt322*gtu22 + gt13L*Gt312*Gt323*gtu22 + 
           gt12L*Gt212*Gt213*gtu31 + gt12L*Gt211*Gt223*gtu31 + 
           Gt211*Gt213*gt22L*gtu31 + gt12L*Gt233*Gt311*gtu31 + 
           Gt213*gt23L*Gt311*gtu31 + gt13L*Gt213*Gt312*gtu31 + 
           Gt113*gt11L*Gt313*gtu31 + Gt211*gt23L*Gt313*gtu31 + 
           gt13L*Gt211*Gt323*gtu31 + gt13L*Gt311*Gt333*gtu31 + 
           Gt311*Gt313*gt33L*gtu31 + gt11L*Gt123*Gt212*gtu32 + 
           gt12L*Gt213*Gt222*gtu32 + gt12L*Gt212*Gt223*gtu32 + 
           Gt212*Gt213*gt22L*gtu32 + gt11L*Gt133*Gt312*gtu32 + 
           gt12L*Gt233*Gt312*gtu32 + Gt213*gt23L*Gt312*gtu32 + 
           Gt212*gt23L*Gt313*gtu32 + gt13L*Gt213*Gt322*gtu32 + 
           gt13L*Gt212*Gt323*gtu32 + gt13L*Gt313*Gt323*gtu32 + 
           gt13L*Gt312*Gt333*gtu32 + Gt312*Gt313*gt33L*gtu32 + 
           gt12L*Gt213*Gt223*gtu33 + gt12L*Gt233*Gt313*gtu33 + 
           Gt213*gt23L*Gt313*gtu33 + gt13L*Gt213*Gt323*gtu33 + 
           gt13L*Gt313*Gt333*gtu33 + gt12L*gtu21*T1000011) + 
        2*gt13L*gtu31*T1000012 + gt33L*gtu33*T1000012 + 3*gt11L*gtu11*T1000013 + 
        3*gt11L*gtu22*T1000014 + gt33L*gtu22*T1000015 + 3*gt11L*gtu33*T1000044 + 
        gt22L*gtu33*T1000070 + Gt111*
         (6*Gt113*gt11L*gtu31 + 4*gt12L*Gt213*gtu31 + gt11L*Xtn1) + 
        Gt211*(2*Gt112*gt11L*gtu11 + 4*Gt111*gt12L*gtu11 + 2*gt11L*Gt122*gtu21 + 
           2*gt11L*Gt123*gtu31 + gt12L*Xtn1) + 
        Gt311*(4*Gt111*gt13L*gtu11 + 2*gt12L*Gt213*gtu11 + 2*gt13L*Gt313*gtu11 + 
           2*gt11L*Gt123*gtu21 + 2*gt11L*Gt133*gtu31 + gt13L*Xtn1) + 
        gt12L*Gt212*Xtn2 + gt13L*Gt312*Xtn2 + 
        Gt112*(6*Gt111*gt11L*gtu21 + 4*gt12L*Gt211*gtu21 + 4*gt12L*Gt212*gtu22 + 
           2*gt11L*Gt213*gtu31 + 6*Gt113*gt11L*gtu32 + gt11L*Xtn2) + 
        Gt113*gt11L*Xtn3 + Gt213*(2*gt11L*Gt122*gtu32 + 4*Gt112*gt12L*gtu32 + 
           2*gt11L*Gt123*gtu33 + gt12L*Xtn3) + 
        Gt313*(4*Gt111*gt13L*gtu31 + 2*gt12L*Gt213*gtu31 + 2*gt11L*Gt123*gtu32 + 
           4*Gt112*gt13L*gtu32 + 2*gt12L*Gt223*gtu32 + 2*gt11L*Gt133*gtu33 + 
           gt13L*Xtn3) + gt22L*gtu11*SQR(Gt211) + gt33L*gtu11*SQR(Gt311);
    
    Rt22  =  4*(Gt122*gt12L*Gt212*gtu21 + Gt112*gt12L*Gt222*gtu21 + 
           Gt122*gt12L*Gt222*gtu22 + Gt123*gt12L*Gt212*gtu31 + 
           Gt123*gt12L*Gt222*gtu32 + Gt123*gt12L*Gt223*gtu33) - 
        gtu11*khalf*PDstandardNth11gt22 + 
        gtu21*(6*Gt212*Gt222*gt22L + 2*Gt122*gt23L*Gt311 + 2*Gt122*gt13L*Gt312 + 
           4*Gt222*gt23L*Gt312 + 2*Gt113*gt12L*Gt322 + 2*gt23L*Gt312*Gt323 + 
           2*Gt312*Gt322*gt33L - PDstandardNth12gt22) + 
        gtu31*(6*Gt212*Gt223*gt22L + 2*Gt123*gt13L*Gt312 + 2*Gt112*gt23L*Gt313 + 
           2*Gt113*gt12L*Gt323 + 2*gt23L*Gt312*Gt333 + 2*Gt312*Gt323*gt33L - 
           PDstandardNth13gt22) - gtu22*khalf*PDstandardNth22gt22 + 
        gtu32*(4*Gt122*gt12L*Gt223 + 2*Gt123*Gt212*gt22L + 2*gt12L*Gt133*Gt322 + 
           4*Gt223*gt23L*Gt322 + 2*Gt123*gt12L*Gt323 + 4*Gt222*gt23L*Gt323 + 
           2*gt23L*Gt322*Gt333 + 2*Gt322*Gt323*gt33L - PDstandardNth23gt22) + 
        gt12L*(2*Gt111*Gt123*gtu31 + 4*Gt112*Gt223*gtu31 + 2*Gt113*Gt122*gtu32 + 
           2*Gt113*Gt123*gtu33 + PDstandardNth2Xt1) + 
        gt22L*(2*Gt122*Gt213*gtu32 + 6*Gt222*Gt223*gtu32 + 2*Gt123*Gt213*gtu33 + 
           PDstandardNth2Xt2) + gt23L*
         (4*Gt212*Gt322*gtu21 + 2*Gt313*Gt322*gtu21 + 4*Gt222*Gt322*gtu22 + 
           2*Gt123*Gt311*gtu31 + 4*Gt212*Gt323*gtu31 + 2*Gt313*Gt323*gtu31 + 
           2*Gt122*Gt313*gtu32 + 2*Gt123*Gt313*gtu33 + 4*Gt223*Gt323*gtu33 + 
           2*Gt323*Gt333*gtu33 + PDstandardNth2Xt3) - 
        gtu33*khalf*PDstandardNth33gt22 + 3*gt22L*gtu11*T1000011 + 
        gt11L*gtu11*T1000014 + 2*(Gt112*Gt211*gt22L*gtu11 + 
           Gt112*gt23L*Gt311*gtu11 + Gt113*gt12L*Gt312*gtu11 + 
           Gt112*gt13L*Gt312*gtu11 + Gt111*Gt122*gt12L*gtu21 + 
           Gt122*Gt211*gt22L*gtu21 + Gt112*Gt212*gt22L*gtu21 + 
           Gt223*gt22L*Gt312*gtu21 + Gt112*gt23L*Gt312*gtu21 + 
           Gt112*gt13L*Gt322*gtu21 + Gt213*gt22L*Gt322*gtu21 + 
           Gt122*Gt212*gt22L*gtu22 + Gt123*gt12L*Gt322*gtu22 + 
           Gt122*gt13L*Gt322*gtu22 + Gt223*gt22L*Gt322*gtu22 + 
           gt23L*Gt322*Gt323*gtu22 + Gt112*Gt113*gt12L*gtu31 + 
           Gt123*Gt211*gt22L*gtu31 + Gt112*Gt213*gt22L*gtu31 + 
           Gt112*gt13L*Gt323*gtu31 + Gt213*gt22L*Gt323*gtu31 + 
           gt11L*Gt122*Gt123*gtu32 + Gt123*gt13L*Gt322*gtu32 + 
           gt22L*Gt233*Gt322*gtu32 + Gt122*gt13L*Gt323*gtu32 + 
           Gt223*gt22L*Gt323*gtu32 + gt12L*Gt133*Gt323*gtu33 + 
           Gt123*gt13L*Gt323*gtu33 + gt22L*Gt233*Gt323*gtu33 + gt12L*gtu21*T1000014)
          + gt33L*gtu11*T1000015 + 2*gt23L*gtu32*T1000066 + gt33L*gtu33*T1000066 + 
        3*gt22L*gtu33*T1000084 + 3*gt22L*gtu22*T1000099 + gt11L*gtu33*T1000142 + 
        Gt212*gt22L*Xtn1 + Gt112*(2*Gt111*gt12L*gtu11 + 4*gt12L*Gt212*gtu11 + 
           2*gt11L*Gt122*gtu21 + 2*Gt122*gt12L*gtu22 + 2*gt11L*Gt123*gtu31 + 
           2*Gt123*gt12L*gtu32 + gt12L*Xtn1) + 
        Gt312*(2*Gt213*gt22L*gtu11 + 4*Gt212*gt23L*gtu11 + 2*gt23L*Gt313*gtu11 + 
           2*Gt123*gt12L*gtu21 + 2*Gt122*gt23L*gtu22 + 2*gt12L*Gt133*gtu31 + 
           2*gt22L*Gt233*gtu31 + 4*Gt223*gt23L*gtu31 + 2*Gt123*gt23L*gtu32 + 
           gt23L*Xtn1) + Gt122*gt12L*Xtn2 + Gt222*gt22L*Xtn2 + gt23L*Gt322*Xtn2 + 
        Gt123*gt12L*Xtn3 + Gt223*gt22L*Xtn3 + gt23L*Gt323*Xtn3 + 
        gt11L*gtu22*SQR(Gt122) + gt33L*gtu22*SQR(Gt322);
    
    Rt33  =  4*(Gt123*gt13L*Gt323*gtu22 + Gt223*gt23L*Gt323*gtu22 + 
           Gt133*gt13L*Gt313*gtu31 + Gt233*gt23L*Gt313*gtu31 + 
           Gt113*gt13L*Gt333*gtu31 + Gt133*gt13L*Gt323*gtu32 + 
           Gt233*gt23L*Gt323*gtu32 + Gt123*gt13L*Gt333*gtu32 + 
           Gt133*gt13L*Gt333*gtu33) + 
        gtu21*(2*Gt212*Gt223*gt23L + 4*Gt123*gt13L*Gt313 + 4*Gt223*gt23L*Gt313 + 
           4*Gt113*gt13L*Gt323 + 4*Gt213*gt23L*Gt323 + 2*Gt123*Gt311*gt33L - 
           PDstandardNth12gt33) + gtu31*
         (4*Gt213*gt23L*Gt333 + 2*Gt233*Gt312*gt33L + 6*Gt313*Gt333*gt33L - 
           PDstandardNth13gt33) - gtu22*khalf*PDstandardNth22gt33 + 
        gtu32*(4*Gt223*gt23L*Gt333 + 2*Gt123*Gt313*gt33L + 6*Gt323*Gt333*gt33L - 
           PDstandardNth23gt33) - gtu33*khalf*PDstandardNth33gt33 + 
        gt13L*PDstandardNth3Xt1 + gt23L*PDstandardNth3Xt2 + 
        gt33L*(2*Gt213*Gt322*gtu21 + 6*Gt313*Gt323*gtu21 + 2*Gt123*Gt312*gtu22 + 
           2*Gt133*Gt311*gtu31 + 2*Gt133*Gt312*gtu32 + 2*Gt133*Gt313*gtu33 + 
           PDstandardNth3Xt3) + gtu11*
         (2*Gt212*Gt213*gt23L + 4*Gt113*gt13L*Gt313 + 4*Gt213*gt23L*Gt313 + 
           2*Gt113*Gt311*gt33L + 2*Gt213*Gt312*gt33L - khalf*PDstandardNth11gt33 + 
           3*gt33L*T1000012) + gt11L*gtu11*T1000044 + 
        2*(Gt111*Gt113*gt13L*gtu11 + Gt113*gt12L*Gt213*gtu11 + 
           Gt112*gt13L*Gt213*gtu11 + Gt113*Gt211*gt23L*gtu11 + 
           Gt113*gt11L*Gt123*gtu21 + Gt112*Gt113*gt13L*gtu21 + 
           Gt111*Gt123*gt13L*gtu21 + Gt123*gt12L*Gt213*gtu21 + 
           Gt122*gt13L*Gt213*gtu21 + Gt113*gt12L*Gt223*gtu21 + 
           Gt112*gt13L*Gt223*gtu21 + Gt213*Gt223*gt22L*gtu21 + 
           Gt123*Gt211*gt23L*gtu21 + Gt113*Gt212*gt23L*gtu21 + 
           Gt213*Gt222*gt23L*gtu21 + Gt113*Gt312*gt33L*gtu21 + 
           Gt223*Gt312*gt33L*gtu21 + Gt112*Gt123*gt13L*gtu22 + 
           Gt123*gt12L*Gt223*gtu22 + Gt122*gt13L*Gt223*gtu22 + 
           Gt123*Gt212*gt23L*gtu22 + Gt222*Gt223*gt23L*gtu22 + 
           Gt223*Gt322*gt33L*gtu22 + Gt113*gt11L*Gt133*gtu31 + 
           Gt111*Gt133*gt13L*gtu31 + gt12L*Gt133*Gt213*gtu31 + 
           Gt123*gt13L*Gt213*gtu31 + Gt113*gt12L*Gt233*gtu31 + 
           Gt112*gt13L*Gt233*gtu31 + Gt213*gt22L*Gt233*gtu31 + 
           Gt133*Gt211*gt23L*gtu31 + Gt113*Gt213*gt23L*gtu31 + 
           Gt213*Gt223*gt23L*gtu31 + Gt212*Gt233*gt23L*gtu31 + 
           Gt113*Gt313*gt33L*gtu31 + Gt213*Gt323*gt33L*gtu31 + 
           gt11L*Gt123*Gt133*gtu32 + Gt113*Gt123*gt13L*gtu32 + 
           Gt112*Gt133*gt13L*gtu32 + gt12L*Gt133*Gt223*gtu32 + 
           Gt123*gt13L*Gt223*gtu32 + Gt123*gt12L*Gt233*gtu32 + 
           Gt122*gt13L*Gt233*gtu32 + Gt223*gt22L*Gt233*gtu32 + 
           Gt133*Gt212*gt23L*gtu32 + Gt123*Gt213*gt23L*gtu32 + 
           Gt222*Gt233*gt23L*gtu32 + Gt233*Gt322*gt33L*gtu32 + 
           Gt223*Gt323*gt33L*gtu32 + Gt113*Gt133*gt13L*gtu33 + 
           gt12L*Gt133*Gt233*gtu33 + Gt123*gt13L*Gt233*gtu33 + 
           Gt133*Gt213*gt23L*gtu33 + Gt223*Gt233*gt23L*gtu33 + gt13L*gtu31*T1000044)
          + 3*gt33L*gtu22*T1000066 + gt22L*gtu11*T1000070 + gt22L*gtu22*T1000084 + 
        2*gt23L*gtu32*T1000084 + gt11L*gtu22*T1000142 + 3*gt33L*gtu33*T1000148 + 
        Gt113*gt13L*Xtn1 + Gt213*gt23L*Xtn1 + Gt313*gt33L*Xtn1 + Gt123*gt13L*Xtn2 + 
        Gt223*gt23L*Xtn2 + Gt323*gt33L*Xtn2 + Gt133*gt13L*Xtn3 + Gt333*gt33L*Xtn3 + 
        Gt233*(4*gt23L*Gt333*gtu33 + 2*Gt323*gt33L*gtu33 + gt23L*Xtn3) + 
        gt11L*gtu33*SQR(Gt133) + gt22L*gtu33*SQR(Gt233);
    
    Rphi11  =  -2*(cdphi211 - 2*T1000150 + 
          gt11L*(cdphi233*gtu33 + 4*(gtu32*pdphi2*pdphi3 + 
                pdphi1*(gtu21*pdphi2 + gtu31*pdphi3)) + 
             gtu11*(cdphi211 + 2*T1000150) + gtu22*(cdphi222 + 2*T1000151) + 
             2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + gtu33*T1000152)));
    
    Rphi12  =  -2*(cdphi212 + pdphi1*(-2*pdphi2 + 4*gt12L*gtu31*pdphi3) + 
          gt12L*(cdphi233*gtu33 + 4*pdphi2*(gtu21*pdphi1 + gtu32*pdphi3) + 
             gtu11*(cdphi211 + 2*T1000150) + gtu22*(cdphi222 + 2*T1000151) + 
             2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + gtu33*T1000152)));
    
    Rphi13  =  -2*(cdphi213 + (-2*pdphi1 + 4*gt13L*gtu32*pdphi2)*pdphi3 + 
          gt13L*(cdphi233*gtu33 + 4*pdphi1*(gtu21*pdphi2 + gtu31*pdphi3) + 
             gtu11*(cdphi211 + 2*T1000150) + gtu22*(cdphi222 + 2*T1000151) + 
             2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + gtu33*T1000152)));
    
    Rphi22  =  -2*(cdphi222 - 2*T1000151 + 
          gt22L*(cdphi233*gtu33 + 4*(gtu32*pdphi2*pdphi3 + 
                pdphi1*(gtu21*pdphi2 + gtu31*pdphi3)) + 
             gtu11*(cdphi211 + 2*T1000150) + gtu22*(cdphi222 + 2*T1000151) + 
             2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + gtu33*T1000152)));
    
    Rphi23  =  -2*(cdphi223 + (-2 + 4*gt23L*gtu32)*pdphi2*pdphi3 + 
          gt23L*(cdphi233*gtu33 + 4*pdphi1*(gtu21*pdphi2 + gtu31*pdphi3) + 
             gtu11*(cdphi211 + 2*T1000150) + gtu22*(cdphi222 + 2*T1000151) + 
             2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + gtu33*T1000152)));
    
    Rphi33  =  -2*(cdphi233 - 2*T1000152 + 
          gt33L*(cdphi233*gtu33 + 4*(gtu32*pdphi2*pdphi3 + 
                pdphi1*(gtu21*pdphi2 + gtu31*pdphi3)) + 
             gtu11*(cdphi211 + 2*T1000150) + gtu22*(cdphi222 + 2*T1000151) + 
             2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + gtu33*T1000152)));
    
    Atu11  =  Atm11*gtu11 + Atm12*gtu21 + Atm13*gtu31;
    
    Atu21  =  Atm11*gtu21 + Atm12*gtu22 + Atm13*gtu32;
    
    Atu31  =  Atm11*gtu31 + Atm12*gtu32 + Atm13*gtu33;
    
    Atu22  =  Atm21*gtu21 + Atm22*gtu22 + Atm23*gtu32;
    
    Atu32  =  Atm21*gtu31 + Atm22*gtu32 + Atm23*gtu33;
    
    Atu33  =  Atm31*gtu31 + Atm32*gtu32 + Atm33*gtu33;
    
    G211  =  Gt211 - 2*gt11L*T1000167;
    
    G311  =  Gt311 - 2*gt11L*T1000159;
    
    G112  =  Gt112 - 2*(gt12L*T1000163 + T1000170);
    
    G312  =  Gt312 - 2*gt12L*T1000159;
    
    G113  =  Gt113 - 2*(gt13L*T1000163 + T1000169);
    
    G213  =  Gt213 - 2*gt13L*T1000167;
    
    G122  =  Gt122 - 2*gt22L*T1000163;
    
    G222  =  Gt222 - 2*(-2*pdphi2 + gt22L*T1000167);
    
    G322  =  Gt322 - 2*gt22L*T1000159;
    
    G123  =  Gt123 - 2*gt23L*T1000163;
    
    G223  =  Gt223 - 2*(gt23L*T1000167 + T1000169);
    
    G323  =  Gt323 - 2*(gt23L*T1000159 + T1000170);
    
    G133  =  Gt133 - 2*gt33L*T1000163;
    
    G233  =  Gt233 - 2*gt33L*T1000167;
    
    G333  =  Gt333 - 2*(-2*pdphi3 + gt33L*T1000159);
    
    Rt12  =  khalf*(-(gtu11*PDstandardNth11gt12) - 2*gtu21*PDstandardNth12gt12 - 
          2*gtu31*PDstandardNth13gt12 + gt12L*PDstandardNth1Xt1 + 
          gt22L*PDstandardNth1Xt2 + gt23L*PDstandardNth1Xt3 - 
          gtu22*PDstandardNth22gt12 - 2*gtu32*PDstandardNth23gt12 + 
          gt11L*PDstandardNth2Xt1 + gt12L*PDstandardNth2Xt2 + 
          gt13L*PDstandardNth2Xt3 - gtu33*PDstandardNth33gt12 + 
          2*gtu21*(Gt112*gt11L*Gt222 + Gt112*Gt211*gt22L + 2*Gt122*gt13L*Gt311 + 
             Gt112*gt23L*Gt311 + Gt222*gt23L*Gt311 + gt13L*Gt222*Gt312 + 
             Gt213*gt22L*Gt312 + Gt212*gt23L*Gt312 + gt23L*Gt312*Gt313 + 
             Gt113*gt11L*Gt322 + gt13L*Gt313*Gt322 + Gt311*Gt322*gt33L + T1000016 + 
             T1000017 + T1000018 + gt12L*
              (2*Gt122*Gt211 + Gt212*Gt222 + Gt213*Gt322 + T1000020 + T1000021) + 
             Gt111*(T1000019 + T1000022 + T1000023 + T1000026)) + 
          2*gtu11*(3*Gt112*gt12L*Gt211 + 2*Gt211*Gt212*gt22L + Gt113*gt12L*Gt311 + 
             2*Gt112*gt13L*Gt311 + Gt213*gt22L*Gt311 + Gt212*gt23L*Gt311 + 
             gt13L*Gt212*Gt312 + gt12L*Gt213*Gt312 + 2*Gt211*gt23L*Gt312 + 
             gt23L*Gt311*Gt313 + gt13L*Gt312*Gt313 + Gt311*Gt312*gt33L + 
             gt12L*T1000011 + gt12L*T1000013 + 
             gt11L*(2*Gt111*Gt112 + T1000020 + T1000021) + 
             Gt111*(T1000045 + T1000046 + T1000047 + T1000048)) + 
          2*gtu21*(Gt122*gt12L*Gt211 + 3*Gt112*gt12L*Gt212 + gt12L*Gt212*Gt222 + 
             Gt123*gt12L*Gt311 + Gt223*gt22L*Gt311 + 3*Gt112*gt13L*Gt312 + 
             gt12L*Gt223*Gt312 + 2*Gt212*gt23L*Gt312 + gt13L*Gt212*Gt322 + 
             gt23L*Gt311*Gt323 + gt13L*Gt312*Gt323 + gt33L*T1000015 + T1000016 + 
             T1000017 + T1000018 + gt11L*(Gt122*Gt212 + 2*T1000014 + T1000027) + 
             Gt111*T1000049) + 2*gtu22*
           (gt11L*Gt122*Gt222 + 2*Gt212*Gt222*gt22L + 2*Gt122*gt13L*Gt312 + 
             Gt223*gt22L*Gt312 + Gt222*gt23L*Gt312 + gt11L*Gt123*Gt322 + 
             gt13L*Gt222*Gt322 + 2*Gt212*gt23L*Gt322 + gt23L*Gt312*Gt323 + 
             gt13L*Gt322*Gt323 + Gt312*Gt322*gt33L + gt12L*T1000014 + 
             Gt112*(T1000022 + T1000023 + T1000024 + T1000025 + T1000026) + 
             gt12L*(3*Gt122*Gt212 + Gt223*Gt322 + T1000027 + T1000099)) + 
          2*gtu31*(Gt123*gt12L*Gt211 + gt12L*Gt133*Gt311 + gt22L*Gt233*Gt311 + 
             gt12L*Gt233*Gt312 + 2*Gt112*gt13L*Gt313 + gt13L*Gt212*Gt323 + 
             gt23L*Gt311*Gt333 + gt13L*Gt312*Gt333 + T1000028 + T1000029 + 
             T1000030 + T1000031 + Gt111*T1000035 + T1000052 + T1000055 + 
             T1000058 + T1000060 + gt11L*(Gt123*Gt212 + Gt133*Gt312 + T1000061) + 
             T1000100 + T1000102) + 2*gtu31*
           (Gt112*gt11L*Gt223 + Gt113*gt11L*Gt323 + T1000028 + T1000029 + 
             T1000030 + T1000031 + Gt111*T1000043 + 
             gt12L*(2*Gt123*Gt211 + Gt212*Gt223 + Gt213*Gt323 + T1000050 + 
                T1000051) + T1000053 + T1000054 + T1000056 + T1000057 + T1000059 + 
             T1000103 + T1000106 + T1000109 + T1000110) + 
          2*gtu32*(gt11L*Gt122*Gt223 + 2*Gt123*gt13L*Gt312 + gt11L*Gt123*Gt323 + 
             T1000036 + T1000037 + T1000038 + T1000039 + Gt112*T1000043 + 
             gt12L*(Gt123*(2*Gt212 + Gt313) + Gt223*(Gt222 + Gt323) + T1000062) + 
             T1000065 + T1000067 + T1000117 + T1000119 + T1000124 + T1000126 + 
             T1000127 + T1000128) + 2*gtu32*
           (gt12L*Gt133*Gt312 + gt22L*Gt233*Gt312 + 2*Gt122*gt13L*Gt313 + 
             gt12L*Gt233*Gt322 + gt11L*
              (2*Gt113*Gt122 + Gt123*Gt222 + Gt133*Gt322) + gt13L*Gt222*Gt323 + 
             gt23L*Gt312*Gt333 + gt13L*Gt322*Gt333 + Gt112*T1000035 + T1000036 + 
             T1000037 + T1000038 + T1000039 + T1000064 + T1000113 + T1000116 + 
             T1000118 + T1000122 + T1000123 + T1000129) + 
          2*gtu33*(gt11L*Gt123*Gt223 + gt11L*Gt133*Gt323 + 
             Gt113*(T1000033 + T1000034 + T1000040 + T1000041 + T1000042) + 
             T1000077 + T1000078 + T1000083 + 
             gt12L*(3*Gt123*Gt213 + Gt233*Gt323 + T1000072 + T1000084) + T1000086 + 
             T1000087 + T1000136 + T1000137 + T1000138 + T1000139 + T1000140) + 
          (Gt111*gt12L + T1000045 + T1000046)*Xtn1 + 
          (Gt112*gt11L + T1000047 + T1000048)*Xtn1 + 
          (gt11L*Gt122 + T1000023 + T1000026)*Xtn2 + T1000049*Xtn2 + 
          T1000035*Xtn3 + T1000097*Xtn3);
    
    Rt13  =  khalf*(-(gtu11*PDstandardNth11gt13) - 2*gtu21*PDstandardNth12gt13 - 
          2*gtu31*PDstandardNth13gt13 + gt13L*PDstandardNth1Xt1 + 
          gt23L*PDstandardNth1Xt2 + gt33L*PDstandardNth1Xt3 - 
          gtu22*PDstandardNth22gt13 - 2*gtu32*PDstandardNth23gt13 - 
          gtu33*PDstandardNth33gt13 + gt11L*PDstandardNth3Xt1 + 
          gt12L*PDstandardNth3Xt2 + gt13L*PDstandardNth3Xt3 + 
          2*gtu31*(Gt112*gt13L*Gt213 + Gt112*gt11L*Gt233 + Gt211*gt22L*Gt233 + 
             Gt113*Gt211*gt23L + Gt212*Gt213*gt23L + 2*Gt133*gt13L*Gt311 + 
             gt13L*Gt233*Gt312 + Gt113*gt13L*Gt313 + Gt213*gt23L*Gt313 + 
             Gt113*gt11L*Gt333 + Gt211*gt23L*Gt333 + 
             gt12L*(2*Gt133*Gt211 + Gt212*Gt233 + Gt213*Gt333) + 
             Gt113*Gt311*gt33L + Gt213*Gt312*gt33L + T1000071 + T1000073 + 
             T1000074 + T1000075 + Gt111*T1000082) + 
          2*gtu31*(Gt123*gt13L*Gt211 + 3*Gt113*gt12L*Gt213 + gt12L*Gt213*Gt223 + 
             Gt211*Gt223*gt23L + Gt133*gt13L*Gt311 + 3*Gt113*gt13L*Gt313 + 
             gt12L*Gt233*Gt313 + 2*Gt213*gt23L*Gt313 + gt13L*Gt213*Gt323 + 
             Gt211*Gt323*gt33L + gt22L*T1000070 + T1000071 + 
             gt11L*(Gt123*Gt213 + 2*T1000044 + T1000072) + T1000073 + T1000074 + 
             T1000075 + Gt111*T1000091) + 
          2*gtu11*(2*Gt113*gt12L*Gt211 + Gt112*gt13L*Gt211 + gt12L*Gt212*Gt213 + 
             Gt211*Gt213*gt22L + Gt211*Gt212*gt23L + 3*Gt113*gt13L*Gt311 + 
             2*Gt213*gt23L*Gt311 + gt13L*Gt213*Gt312 + gt12L*Gt213*Gt313 + 
             Gt211*gt23L*Gt313 + Gt211*Gt312*gt33L + 2*Gt311*Gt313*gt33L + 
             gt13L*T1000012 + gt13L*T1000013 + 
             gt11L*(2*Gt111*Gt113 + T1000050 + T1000051) + 
             Gt111*(T1000092 + T1000093 + T1000094 + T1000095)) + 
          2*gtu21*(Gt122*gt13L*Gt211 + 2*Gt113*gt12L*Gt212 + Gt112*gt12L*Gt213 + 
             gt12L*Gt213*Gt222 + Gt211*Gt222*gt23L + Gt123*gt13L*Gt311 + 
             gt12L*Gt223*Gt313 + gt13L*Gt213*Gt322 + Gt211*Gt322*gt33L + T1000028 + 
             T1000030 + T1000054 + T1000057 + T1000058 + T1000059 + T1000060 + 
             gt11L*(Gt123*Gt313 + T1000061 + T1000062) + Gt111*T1000098 + 
             T1000107 + T1000108) + 2*gtu21*
           (2*Gt123*gt12L*Gt211 + gt12L*Gt213*Gt323 + 
             gt11L*(2*Gt111*Gt123 + Gt112*Gt223 + Gt113*Gt323) + T1000029 + 
             T1000030 + T1000031 + T1000052 + T1000053 + T1000054 + T1000055 + 
             T1000056 + T1000057 + T1000058 + T1000059 + 
             Gt111*(T1000041 + T1000042 + T1000063) + T1000101 + T1000104 + 
             T1000105 + T1000111 + T1000112) + 
          2*gtu22*(2*Gt123*gt12L*Gt212 + 3*Gt123*gt13L*Gt312 + gt12L*Gt223*Gt323 + 
             gt11L*(2*Gt112*Gt123 + Gt122*Gt223 + Gt123*Gt323) + 
             2*Gt312*Gt323*gt33L + T1000037 + T1000039 + T1000064 + T1000065 + 
             T1000067 + Gt112*(T1000041 + T1000042 + T1000068 + T1000069) + 
             T1000114 + T1000115 + T1000120 + T1000121 + T1000130) + 
          2*gtu32*(Gt122*gt13L*Gt213 + gt11L*Gt122*Gt233 + Gt212*gt22L*Gt233 + 
             Gt113*Gt212*gt23L + Gt213*Gt222*gt23L + 2*Gt133*gt13L*Gt312 + 
             Gt123*gt13L*Gt313 + gt13L*Gt233*Gt322 + gt11L*Gt123*Gt333 + 
             gt12L*(2*Gt133*Gt212 + Gt222*Gt233 + Gt223*Gt333) + 
             Gt113*Gt312*gt33L + Gt213*Gt322*gt33L + T1000077 + T1000078 + 
             Gt112*T1000082 + T1000083 + T1000085 + T1000088 + T1000135) + 
          2*gtu32*(Gt123*gt13L*Gt212 + 2*Gt123*gt12L*Gt213 + Gt113*gt12L*Gt223 + 
             Gt213*Gt223*gt22L + gt12L*Gt233*Gt323 + 
             gt11L*(2*Gt113*Gt123 + Gt123*Gt223 + Gt133*Gt323) + T1000077 + 
             T1000078 + T1000083 + gt12L*T1000084 + T1000085 + T1000086 + 
             T1000087 + T1000088 + Gt112*T1000091 + T1000131 + T1000132 + 
             T1000133 + T1000134 + T1000141) + 
          2*gtu33*(2*gt12L*Gt133*Gt213 + Gt123*gt13L*Gt213 + gt11L*Gt123*Gt233 + 
             gt12L*Gt223*Gt233 + Gt213*gt22L*Gt233 + Gt213*Gt223*gt23L + 
             3*Gt133*gt13L*Gt313 + 2*Gt233*gt23L*Gt313 + gt13L*Gt233*Gt323 + 
             gt11L*Gt133*Gt333 + gt12L*Gt233*Gt333 + Gt213*gt23L*Gt333 + 
             Gt213*Gt323*gt33L + 2*Gt313*Gt333*gt33L + gt13L*T1000044 + 
             Gt113*(T1000079 + T1000080 + T1000081 + T1000089 + T1000090) + 
             gt13L*T1000148) + (Gt113*gt11L + T1000092 + T1000093)*Xtn1 + 
          (Gt111*gt13L + T1000094 + T1000095)*Xtn1 + T1000097*Xtn2 + 
          T1000098*Xtn2 + (gt11L*Gt133 + T1000080 + T1000081)*Xtn3 + T1000091*Xtn3);
    
    Rt23  =  khalf*(-(gtu11*PDstandardNth11gt23) - 2*gtu21*PDstandardNth12gt23 - 
          2*gtu31*PDstandardNth13gt23 - gtu22*PDstandardNth22gt23 - 
          2*gtu32*PDstandardNth23gt23 + gt13L*PDstandardNth2Xt1 + 
          gt23L*PDstandardNth2Xt2 + gt33L*PDstandardNth2Xt3 - 
          gtu33*PDstandardNth33gt23 + gt12L*PDstandardNth3Xt1 + 
          gt22L*PDstandardNth3Xt2 + gt23L*PDstandardNth3Xt3 + 
          2*gtu22*(gt11L*Gt122*Gt123 + Gt112*Gt123*gt12L + Gt112*Gt122*gt13L + 
             Gt123*gt12L*Gt222 + Gt122*gt13L*Gt222 + 2*Gt122*gt12L*Gt223 + 
             Gt123*Gt212*gt22L + 2*Gt222*Gt223*gt22L + Gt122*Gt212*gt23L + 
             Gt123*gt23L*Gt312 + 2*Gt123*gt13L*Gt322 + 3*Gt223*gt23L*Gt322 + 
             Gt123*gt12L*Gt323 + Gt122*gt13L*Gt323 + Gt223*gt22L*Gt323 + 
             Gt222*gt23L*Gt323 + Gt122*Gt312*gt33L + Gt222*Gt322*gt33L + 
             2*Gt322*Gt323*gt33L + gt23L*T1000066 + gt23L*T1000099) + 
          2*gtu11*(Gt112*Gt113*gt11L + Gt111*Gt113*gt12L + Gt111*Gt112*gt13L + 
             2*Gt212*Gt213*gt22L + 3*Gt213*gt23L*Gt312 + Gt113*gt12L*Gt313 + 
             2*Gt312*Gt313*gt33L + T1000060 + T1000100 + T1000101 + T1000102 + 
             T1000103 + T1000104 + T1000105 + T1000106 + T1000107 + T1000108 + 
             T1000109 + T1000110 + T1000111 + T1000112) + 
          2*gtu21*(Gt112*gt11L*Gt123 + Gt111*Gt123*gt12L + Gt111*Gt122*gt13L + 
             Gt112*gt13L*Gt222 + 2*Gt112*gt12L*Gt223 + Gt123*Gt211*gt22L + 
             2*Gt212*Gt223*gt22L + Gt122*Gt211*gt23L + Gt123*gt23L*Gt311 + 
             Gt113*gt12L*Gt323 + Gt112*gt13L*Gt323 + Gt213*gt22L*Gt323 + 
             Gt122*Gt311*gt33L + Gt222*Gt312*gt33L + T1000038 + T1000039 + 
             T1000113 + T1000114 + T1000115 + T1000116 + T1000117 + T1000118 + 
             T1000119 + T1000125) + 2*gtu21*
           (Gt113*gt11L*Gt122 + 2*Gt213*Gt222*gt22L + Gt123*gt12L*Gt313 + 
             Gt122*gt13L*Gt313 + 2*Gt213*gt23L*Gt322 + 
             Gt112*(T1000032 + T1000068 + T1000069) + T1000114 + T1000116 + 
             T1000117 + T1000118 + T1000119 + T1000120 + T1000121 + T1000122 + 
             T1000123 + T1000124 + T1000125 + T1000126 + T1000127 + T1000128 + 
             T1000129 + T1000130) + 2*gtu31*
           (Gt112*gt11L*Gt133 + Gt111*gt12L*Gt133 + Gt111*Gt123*gt13L + 
             gt12L*Gt133*Gt212 + Gt112*gt13L*Gt223 + Gt133*Gt211*gt22L + 
             2*Gt112*gt12L*Gt233 + 2*Gt212*gt22L*Gt233 + Gt123*Gt211*gt23L + 
             Gt133*gt23L*Gt311 + 2*Gt233*gt23L*Gt312 + Gt113*gt12L*Gt333 + 
             Gt112*gt13L*Gt333 + Gt213*gt22L*Gt333 + Gt123*Gt311*gt33L + 
             Gt223*Gt312*gt33L + T1000083 + T1000088 + T1000131 + T1000132 + 
             T1000133 + T1000134 + T1000135 + T1000136) + 
          2*gtu31*(Gt112*Gt213*gt23L + gt12L*Gt133*Gt313 + 
             Gt123*(2*gt12L*Gt213 + gt13L*(Gt212 + Gt313)) + Gt112*Gt313*gt33L + 
             T1000077 + T1000083 + T1000085 + T1000088 + 
             Gt113*(T1000033 + T1000034 + T1000041 + T1000042 + T1000063 + 
                T1000096) + T1000131 + T1000132 + T1000136 + T1000137 + T1000138 + 
             T1000139 + T1000140 + T1000141) + 
          2*gtu32*(gt11L*Gt122*Gt133 + Gt112*gt12L*Gt133 + Gt112*Gt123*gt13L + 
             gt12L*Gt133*Gt222 + Gt122*gt13L*Gt223 + Gt133*Gt212*gt22L + 
             2*Gt122*gt12L*Gt233 + 2*Gt222*gt22L*Gt233 + Gt123*Gt212*gt23L + 
             Gt133*gt23L*Gt312 + 2*Gt233*gt23L*Gt322 + Gt123*gt13L*Gt323 + 
             Gt223*gt23L*Gt323 + Gt123*gt12L*Gt333 + Gt122*gt13L*Gt333 + 
             Gt223*gt22L*Gt333 + Gt222*gt23L*Gt333 + Gt123*Gt312*gt33L + 
             Gt223*Gt322*gt33L + T1000143 + T1000144 + T1000145 + T1000146 + 
             T1000147) + 2*gtu32*(Gt113*Gt123*gt12L + Gt113*Gt122*gt13L + 
             Gt123*gt13L*Gt222 + 3*Gt123*gt12L*Gt223 + Gt123*Gt213*gt22L + 
             Gt122*Gt213*gt23L + Gt123*gt23L*Gt313 + Gt233*gt23L*Gt322 + 
             gt12L*Gt133*Gt323 + 2*Gt123*gt13L*Gt323 + gt22L*Gt233*Gt323 + 
             3*Gt223*gt23L*Gt323 + Gt122*Gt313*gt33L + Gt222*Gt323*gt33L + 
             2*gt22L*T1000084 + gt11L*T1000142 + T1000143 + T1000144 + T1000145 + 
             T1000146 + T1000147) + 2*gtu33*
           (gt11L*Gt123*Gt133 + Gt113*gt12L*Gt133 + Gt113*Gt123*gt13L + 
             gt12L*Gt133*Gt223 + Gt123*gt13L*Gt223 + Gt133*Gt213*gt22L + 
             2*Gt123*gt12L*Gt233 + 2*Gt223*gt22L*Gt233 + Gt123*Gt213*gt23L + 
             Gt133*gt23L*Gt313 + 2*Gt133*gt13L*Gt323 + 3*Gt233*gt23L*Gt323 + 
             gt12L*Gt133*Gt333 + Gt123*gt13L*Gt333 + gt22L*Gt233*Gt333 + 
             Gt223*gt23L*Gt333 + Gt123*Gt313*gt33L + Gt223*Gt323*gt33L + 
             2*Gt323*Gt333*gt33L + gt23L*T1000084 + gt23L*T1000148) + 
          T1000035*Xtn1 + T1000098*Xtn1 + 
          (Gt123*gt12L + Gt223*gt22L + gt23L*Gt323)*Xtn2 + 
          (Gt122*gt13L + Gt222*gt23L + Gt322*gt33L)*Xtn2 + 
          (gt12L*Gt133 + gt22L*Gt233 + gt23L*Gt333)*Xtn3 + 
          (Gt123*gt13L + Gt223*gt23L + Gt323*gt33L)*Xtn3);
    
    R11  =  Rphi11 + Rt11;
    
    R22  =  Rphi22 + Rt22;
    
    R33  =  Rphi33 + Rt33;
    
    Xt1rhsL  =  gtu11*PDstandardNth11beta1 + gtu22*PDstandardNth22beta1 + 
        gtu33*PDstandardNth33beta1 - 
        2*(Atu11*PDstandardNth1alpha + Atu21*PDstandardNth2alpha + 
           Atu31*PDstandardNth3alpha) + 
        2*(gtu21*PDstandardNth12beta1 + gtu31*PDstandardNth13beta1 + 
           gtu32*PDstandardNth23beta1 + 
           alphaL*(Atu11*Gt111 + 2*Atu21*Gt112 + 2*Atu31*Gt113 + Atu22*Gt122 + 
              2*Atu32*Gt123 + Atu33*Gt133 + 
              6*(Atu11*pdphi1 + Atu21*pdphi2 + Atu31*pdphi3) - 
              ktwothird*(gtu11*PDstandardNth1trK + gtu21*PDstandardNth2trK + 
                 gtu31*PDstandardNth3trK))) + betam1*PDupwindmNth1Xt1 + 
        betam2*PDupwindmNth2Xt1 + betam3*PDupwindmNth3Xt1 + 
        betap1*PDupwindpNth1Xt1 + betap2*PDupwindpNth2Xt1 + 
        betap3*PDupwindpNth3Xt1 + kthird*
         (gtu11*T1000172 + gtu21*T1000173 + gtu31*T1000174) + 
        (-PDstandardNth1beta1 + ktwothird*T1000171)*Xtn1 - 
        PDstandardNth2beta1*Xtn2 - PDstandardNth3beta1*Xtn3;
    
    Xt2rhsL  =  gtu11*PDstandardNth11beta2 + gtu22*PDstandardNth22beta2 + 
        gtu33*PDstandardNth33beta2 - 
        2*(Atu21*PDstandardNth1alpha + Atu22*PDstandardNth2alpha + 
           Atu32*PDstandardNth3alpha) + 
        2*(gtu21*PDstandardNth12beta2 + gtu31*PDstandardNth13beta2 + 
           gtu32*PDstandardNth23beta2 + 
           alphaL*(Atu11*Gt211 + 2*Atu21*Gt212 + 2*Atu31*Gt213 + Atu22*Gt222 + 
              2*Atu32*Gt223 + Atu33*Gt233 + 
              6*(Atu21*pdphi1 + Atu22*pdphi2 + Atu32*pdphi3) - 
              ktwothird*(gtu21*PDstandardNth1trK + gtu22*PDstandardNth2trK + 
                 gtu32*PDstandardNth3trK))) + betam1*PDupwindmNth1Xt2 + 
        betam2*PDupwindmNth2Xt2 + betam3*PDupwindmNth3Xt2 + 
        betap1*PDupwindpNth1Xt2 + betap2*PDupwindpNth2Xt2 + 
        betap3*PDupwindpNth3Xt2 + kthird*
         (gtu21*T1000172 + gtu22*T1000173 + gtu32*T1000174) - 
        PDstandardNth1beta2*Xtn1 - PDstandardNth2beta2*Xtn2 + 
        ktwothird*T1000171*Xtn2 - PDstandardNth3beta2*Xtn3;
    
    Xt3rhsL  =  gtu11*PDstandardNth11beta3 + gtu22*PDstandardNth22beta3 + 
        gtu33*PDstandardNth33beta3 - 
        2*(Atu31*PDstandardNth1alpha + Atu32*PDstandardNth2alpha + 
           Atu33*PDstandardNth3alpha) + 
        2*(gtu21*PDstandardNth12beta3 + gtu31*PDstandardNth13beta3 + 
           gtu32*PDstandardNth23beta3 + 
           alphaL*(Atu11*Gt311 + 2*Atu21*Gt312 + 2*Atu31*Gt313 + Atu22*Gt322 + 
              2*Atu32*Gt323 + Atu33*Gt333 + 
              6*(Atu31*pdphi1 + Atu32*pdphi2 + Atu33*pdphi3) - 
              ktwothird*(gtu31*PDstandardNth1trK + gtu32*PDstandardNth2trK + 
                 gtu33*PDstandardNth3trK))) + betam1*PDupwindmNth1Xt3 + 
        betam2*PDupwindmNth2Xt3 + betam3*PDupwindmNth3Xt3 + 
        betap1*PDupwindpNth1Xt3 + betap2*PDupwindpNth2Xt3 + 
        betap3*PDupwindpNth3Xt3 + kthird*
         (gtu31*T1000172 + gtu32*T1000173 + gtu33*T1000174) - 
        PDstandardNth1beta3*Xtn1 - PDstandardNth2beta3*Xtn2 - 
        PDstandardNth3beta3*Xtn3 + ktwothird*T1000171*Xtn3;
    
    trKrhsL  =  -(gu11*PDstandardNth11alpha) - 2*gu21*PDstandardNth12alpha - 
        2*gu31*PDstandardNth13alpha + 
        (G111*gu11 + 2*G112*gu21 + G122*gu22 + 2*G113*gu31 + 2*G123*gu32 + 
           G133*gu33)*PDstandardNth1alpha - gu22*PDstandardNth22alpha - 
        2*gu32*PDstandardNth23alpha + 
        (G211*gu11 + 2*G212*gu21 + G222*gu22 + 2*G213*gu31 + 2*G223*gu32 + 
           G233*gu33)*PDstandardNth2alpha - gu33*PDstandardNth33alpha + 
        G311*gu11*PDstandardNth3alpha + G322*gu22*PDstandardNth3alpha + 
        2*G313*gu31*PDstandardNth3alpha + 2*G323*gu32*PDstandardNth3alpha + 
        G333*gu33*PDstandardNth3alpha + 
        2*(alphaL*Atm12*Atm21 + alphaL*Atm13*Atm31 + alphaL*Atm23*Atm32 + 
           G312*gu21*PDstandardNth3alpha) + betam1*PDupwindmNth1trK + 
        betam2*PDupwindmNth2trK + betam3*PDupwindmNth3trK + 
        betap1*PDupwindpNth1trK + betap2*PDupwindpNth2trK + 
        betap3*PDupwindpNth3trK + alphaL*SQR(Atm11) + alphaL*SQR(Atm22) + 
        alphaL*SQR(Atm33) + alphaL*kthird*SQR(trKL);
    
    R12  =  Rphi12 + Rt12;
    
    R13  =  Rphi13 + Rt13;
    
    R23  =  Rphi23 + Rt23;
    
    Ats11  =  -PDstandardNth11alpha + G111*PDstandardNth1alpha + 
        G211*PDstandardNth2alpha + G311*PDstandardNth3alpha + alphaL*R11;
    
    Ats22  =  G122*PDstandardNth1alpha - PDstandardNth22alpha + 
        G222*PDstandardNth2alpha + G322*PDstandardNth3alpha + alphaL*R22;
    
    Ats33  =  G133*PDstandardNth1alpha + G233*PDstandardNth2alpha - 
        PDstandardNth33alpha + G333*PDstandardNth3alpha + alphaL*R33;
    
    ArhsL  =  T1000175*(AL*AlphaDriver - trKrhsL);
    
    B1rhsL  =  -(B1L*BetaDriver) + (betam1*(PDupwindmNth1B1 - PDupwindmNth1Xt1) + 
           betam2*(PDupwindmNth2B1 - PDupwindmNth2Xt1) + 
           betam3*(PDupwindmNth3B1 - PDupwindmNth3Xt1) + 
           betap1*(PDupwindpNth1B1 - PDupwindpNth1Xt1) + 
           betap2*(PDupwindpNth2B1 - PDupwindpNth2Xt1) + 
           betap3*(PDupwindpNth3B1 - PDupwindpNth3Xt1))*ShiftAdvectionCoeff + 
        Xt1rhsL;
    
    B2rhsL  =  -(B2L*BetaDriver) + (betam1*(PDupwindmNth1B2 - PDupwindmNth1Xt2) + 
           betam2*(PDupwindmNth2B2 - PDupwindmNth2Xt2) + 
           betam3*(PDupwindmNth3B2 - PDupwindmNth3Xt2) + 
           betap1*(PDupwindpNth1B2 - PDupwindpNth1Xt2) + 
           betap2*(PDupwindpNth2B2 - PDupwindpNth2Xt2) + 
           betap3*(PDupwindpNth3B2 - PDupwindpNth3Xt2))*ShiftAdvectionCoeff + 
        Xt2rhsL;
    
    B3rhsL  =  -(B3L*BetaDriver) + (betam1*(PDupwindmNth1B3 - PDupwindmNth1Xt3) + 
           betam2*(PDupwindmNth2B3 - PDupwindmNth2Xt3) + 
           betam3*(PDupwindmNth3B3 - PDupwindmNth3Xt3) + 
           betap1*(PDupwindpNth1B3 - PDupwindpNth1Xt3) + 
           betap2*(PDupwindpNth2B3 - PDupwindpNth2Xt3) + 
           betap3*(PDupwindpNth3B3 - PDupwindpNth3Xt3))*ShiftAdvectionCoeff + 
        Xt3rhsL;
    
    Ats12  =  -PDstandardNth12alpha + G112*PDstandardNth1alpha + 
        G212*PDstandardNth2alpha + G312*PDstandardNth3alpha + alphaL*R12;
    
    Ats13  =  -PDstandardNth13alpha + G113*PDstandardNth1alpha + 
        G213*PDstandardNth2alpha + G313*PDstandardNth3alpha + alphaL*R13;
    
    Ats23  =  G123*PDstandardNth1alpha - PDstandardNth23alpha + 
        G223*PDstandardNth2alpha + G323*PDstandardNth3alpha + alphaL*R23;
    
    trAts  =  Ats11*gu11 + Ats22*gu22 + 2*(Ats12*gu21 + Ats13*gu31 + Ats23*gu32) + 
        Ats33*gu33;
    
    At11rhsL  =  2*(At11L*PDstandardNth1beta1 + At12L*PDstandardNth1beta2 + 
           At13L*PDstandardNth1beta3) + betam1*PDupwindmNth1At11 + 
        betam2*PDupwindmNth2At11 + betam3*PDupwindmNth3At11 + 
        betap1*PDupwindpNth1At11 + betap2*PDupwindpNth2At11 + 
        betap3*PDupwindpNth3At11 - At11L*ktwothird*T1000171 + 
        alphaL*(-2*(At11L*Atm11 + At12L*Atm21 + At13L*Atm31) + At11L*trKL) + 
        (Ats11 - g11*kthird*trAts)*W2;
    
    At12rhsL  =  At22L*PDstandardNth1beta2 + At23L*PDstandardNth1beta3 + 
        At11L*PDstandardNth2beta1 + At13L*PDstandardNth2beta3 + 
        betam1*PDupwindmNth1At12 + betam2*PDupwindmNth2At12 + 
        betam3*PDupwindmNth3At12 + betap1*PDupwindpNth1At12 + 
        betap2*PDupwindpNth2At12 + betap3*PDupwindpNth3At12 + 
        At12L*(PDstandardNth1beta1 + PDstandardNth2beta2 - ktwothird*T1000171) + 
        alphaL*(-2*(At11L*Atm12 + At12L*Atm22 + At13L*Atm32) + At12L*trKL) + 
        (Ats12 - g12*kthird*trAts)*W2;
    
    At13rhsL  =  At23L*PDstandardNth1beta2 + At33L*PDstandardNth1beta3 + 
        At11L*PDstandardNth3beta1 + At12L*PDstandardNth3beta2 + 
        betam1*PDupwindmNth1At13 + betam2*PDupwindmNth2At13 + 
        betam3*PDupwindmNth3At13 + betap1*PDupwindpNth1At13 + 
        betap2*PDupwindpNth2At13 + betap3*PDupwindpNth3At13 + 
        At13L*(PDstandardNth1beta1 + PDstandardNth3beta3 - ktwothird*T1000171) + 
        alphaL*(-2*(At11L*Atm13 + At12L*Atm23 + At13L*Atm33) + At13L*trKL) + 
        (Ats13 - g13*kthird*trAts)*W2;
    
    At22rhsL  =  2*(At12L*PDstandardNth2beta1 + At22L*PDstandardNth2beta2 + 
           At23L*PDstandardNth2beta3) + betam1*PDupwindmNth1At22 + 
        betam2*PDupwindmNth2At22 + betam3*PDupwindmNth3At22 + 
        betap1*PDupwindpNth1At22 + betap2*PDupwindpNth2At22 + 
        betap3*PDupwindpNth3At22 - At22L*ktwothird*T1000171 + 
        alphaL*(-2*(At12L*Atm12 + At22L*Atm22 + At23L*Atm32) + At22L*trKL) + 
        (Ats22 - g22*kthird*trAts)*W2;
    
    At23rhsL  =  At13L*PDstandardNth2beta1 + At33L*PDstandardNth2beta3 + 
        At12L*PDstandardNth3beta1 + At22L*PDstandardNth3beta2 + 
        betam1*PDupwindmNth1At23 + betam2*PDupwindmNth2At23 + 
        betam3*PDupwindmNth3At23 + betap1*PDupwindpNth1At23 + 
        betap2*PDupwindpNth2At23 + betap3*PDupwindpNth3At23 + 
        At23L*(PDstandardNth2beta2 + PDstandardNth3beta3 - ktwothird*T1000171) + 
        alphaL*(-2*(At12L*Atm13 + At22L*Atm23 + At23L*Atm33) + At23L*trKL) + 
        (Ats23 - g23*kthird*trAts)*W2;
    
    At33rhsL  =  2*(At13L*PDstandardNth3beta1 + At23L*PDstandardNth3beta2 + 
           At33L*PDstandardNth3beta3) + betam1*PDupwindmNth1At33 + 
        betam2*PDupwindmNth2At33 + betam3*PDupwindmNth3At33 + 
        betap1*PDupwindpNth1At33 + betap2*PDupwindpNth2At33 + 
        betap3*PDupwindpNth3At33 - At33L*ktwothird*T1000171 + 
        alphaL*(-2*(At13L*Atm13 + At23L*Atm23 + At33L*Atm33) + At33L*trKL) + 
        (Ats33 - g33*kthird*trAts)*W2;
    
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    Arhs[index] = ArhsL;
    At11rhs[index] = At11rhsL;
    At12rhs[index] = At12rhsL;
    At13rhs[index] = At13rhsL;
    At22rhs[index] = At22rhsL;
    At23rhs[index] = At23rhsL;
    At33rhs[index] = At33rhsL;
    B1rhs[index] = B1rhsL;
    B2rhs[index] = B2rhsL;
    B3rhs[index] = B3rhsL;
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    gt11rhs[index] = gt11rhsL;
    gt12rhs[index] = gt12rhsL;
    gt13rhs[index] = gt13rhsL;
    gt22rhs[index] = gt22rhsL;
    gt23rhs[index] = gt23rhsL;
    gt33rhs[index] = gt33rhsL;
    trKrhs[index] = trKrhsL;
    Wrhs[index] = WrhsL;
    Xt1rhs[index] = Xt1rhsL;
    Xt2rhs[index] = Xt2rhsL;
    Xt3rhs[index] = Xt3rhsL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_BSSNW_RHS);
}

void ML_BSSNW_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSNW_RHS_Body);
}
