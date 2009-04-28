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

void ML_BSSN_MP_constraints_Body(cGH const * const cctkGH, CCTK_INT const dir, CCTK_INT const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], CCTK_INT const min[3], CCTK_INT const max[3], CCTK_INT const n_subblock_gfs, CCTK_REAL * const subblock_gfs[])
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
  CCTK_REAL p1odx = INITVALUE;
  CCTK_REAL p1ody = INITVALUE;
  CCTK_REAL p1odz = INITVALUE;
  CCTK_REAL pm1o12dx2 = INITVALUE;
  CCTK_REAL pm1o12dy2 = INITVALUE;
  CCTK_REAL pm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_constraints_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_constraints_calc_every != ML_BSSN_MP_constraints_calc_offset)
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
  p1odx = INV(dx);
  p1ody = INV(dy);
  p1odz = INV(dz);
  pm1o12dx2 = -pow(dx,-2)/12.;
  pm1o12dy2 = -pow(dy,-2)/12.;
  pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_MP_constraints,
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
    CCTK_REAL detgt = INITVALUE;
    CCTK_REAL e4phi = INITVALUE;
    CCTK_REAL em4phi = INITVALUE;
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
    CCTK_REAL dJ111L = INITVALUE, dJ112L = INITVALUE, dJ113L = INITVALUE, dJ122L = INITVALUE, dJ123L = INITVALUE, dJ133L = INITVALUE;
    CCTK_REAL dJ211L = INITVALUE, dJ212L = INITVALUE, dJ213L = INITVALUE, dJ222L = INITVALUE, dJ223L = INITVALUE, dJ233L = INITVALUE;
    CCTK_REAL dJ311L = INITVALUE, dJ312L = INITVALUE, dJ313L = INITVALUE, dJ322L = INITVALUE, dJ323L = INITVALUE, dJ333L = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt12L = INITVALUE, gt13L = INITVALUE, gt22L = INITVALUE, gt23L = INITVALUE, gt33L = INITVALUE;
    CCTK_REAL HL = INITVALUE;
    CCTK_REAL J11L = INITVALUE, J12L = INITVALUE, J13L = INITVALUE, J21L = INITVALUE, J22L = INITVALUE, J23L = INITVALUE;
    CCTK_REAL J31L = INITVALUE, J32L = INITVALUE, J33L = INITVALUE;
    CCTK_REAL M1L = INITVALUE, M2L = INITVALUE, M3L = INITVALUE;
    CCTK_REAL phiL = INITVALUE;
    CCTK_REAL trKL = INITVALUE;
    CCTK_REAL Xt1L = INITVALUE, Xt2L = INITVALUE, Xt3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandardNth1At11 = INITVALUE;
    CCTK_REAL PDstandardNth2At11 = INITVALUE;
    CCTK_REAL PDstandardNth3At11 = INITVALUE;
    CCTK_REAL PDstandardNth1At12 = INITVALUE;
    CCTK_REAL PDstandardNth2At12 = INITVALUE;
    CCTK_REAL PDstandardNth3At12 = INITVALUE;
    CCTK_REAL PDstandardNth1At13 = INITVALUE;
    CCTK_REAL PDstandardNth2At13 = INITVALUE;
    CCTK_REAL PDstandardNth3At13 = INITVALUE;
    CCTK_REAL PDstandardNth1At22 = INITVALUE;
    CCTK_REAL PDstandardNth2At22 = INITVALUE;
    CCTK_REAL PDstandardNth3At22 = INITVALUE;
    CCTK_REAL PDstandardNth1At23 = INITVALUE;
    CCTK_REAL PDstandardNth2At23 = INITVALUE;
    CCTK_REAL PDstandardNth3At23 = INITVALUE;
    CCTK_REAL PDstandardNth1At33 = INITVALUE;
    CCTK_REAL PDstandardNth2At33 = INITVALUE;
    CCTK_REAL PDstandardNth3At33 = INITVALUE;
    CCTK_REAL PDstandardNth1gt11 = INITVALUE;
    CCTK_REAL PDstandardNth2gt11 = INITVALUE;
    CCTK_REAL PDstandardNth3gt11 = INITVALUE;
    CCTK_REAL PDstandardNth11gt11 = INITVALUE;
    CCTK_REAL PDstandardNth22gt11 = INITVALUE;
    CCTK_REAL PDstandardNth33gt11 = INITVALUE;
    CCTK_REAL PDstandardNth12gt11 = INITVALUE;
    CCTK_REAL PDstandardNth13gt11 = INITVALUE;
    CCTK_REAL PDstandardNth23gt11 = INITVALUE;
    CCTK_REAL PDstandardNth1gt12 = INITVALUE;
    CCTK_REAL PDstandardNth2gt12 = INITVALUE;
    CCTK_REAL PDstandardNth3gt12 = INITVALUE;
    CCTK_REAL PDstandardNth11gt12 = INITVALUE;
    CCTK_REAL PDstandardNth22gt12 = INITVALUE;
    CCTK_REAL PDstandardNth33gt12 = INITVALUE;
    CCTK_REAL PDstandardNth12gt12 = INITVALUE;
    CCTK_REAL PDstandardNth13gt12 = INITVALUE;
    CCTK_REAL PDstandardNth23gt12 = INITVALUE;
    CCTK_REAL PDstandardNth1gt13 = INITVALUE;
    CCTK_REAL PDstandardNth2gt13 = INITVALUE;
    CCTK_REAL PDstandardNth3gt13 = INITVALUE;
    CCTK_REAL PDstandardNth11gt13 = INITVALUE;
    CCTK_REAL PDstandardNth22gt13 = INITVALUE;
    CCTK_REAL PDstandardNth33gt13 = INITVALUE;
    CCTK_REAL PDstandardNth12gt13 = INITVALUE;
    CCTK_REAL PDstandardNth13gt13 = INITVALUE;
    CCTK_REAL PDstandardNth23gt13 = INITVALUE;
    CCTK_REAL PDstandardNth1gt22 = INITVALUE;
    CCTK_REAL PDstandardNth2gt22 = INITVALUE;
    CCTK_REAL PDstandardNth3gt22 = INITVALUE;
    CCTK_REAL PDstandardNth11gt22 = INITVALUE;
    CCTK_REAL PDstandardNth22gt22 = INITVALUE;
    CCTK_REAL PDstandardNth33gt22 = INITVALUE;
    CCTK_REAL PDstandardNth12gt22 = INITVALUE;
    CCTK_REAL PDstandardNth13gt22 = INITVALUE;
    CCTK_REAL PDstandardNth23gt22 = INITVALUE;
    CCTK_REAL PDstandardNth1gt23 = INITVALUE;
    CCTK_REAL PDstandardNth2gt23 = INITVALUE;
    CCTK_REAL PDstandardNth3gt23 = INITVALUE;
    CCTK_REAL PDstandardNth11gt23 = INITVALUE;
    CCTK_REAL PDstandardNth22gt23 = INITVALUE;
    CCTK_REAL PDstandardNth33gt23 = INITVALUE;
    CCTK_REAL PDstandardNth12gt23 = INITVALUE;
    CCTK_REAL PDstandardNth13gt23 = INITVALUE;
    CCTK_REAL PDstandardNth23gt23 = INITVALUE;
    CCTK_REAL PDstandardNth1gt33 = INITVALUE;
    CCTK_REAL PDstandardNth2gt33 = INITVALUE;
    CCTK_REAL PDstandardNth3gt33 = INITVALUE;
    CCTK_REAL PDstandardNth11gt33 = INITVALUE;
    CCTK_REAL PDstandardNth22gt33 = INITVALUE;
    CCTK_REAL PDstandardNth33gt33 = INITVALUE;
    CCTK_REAL PDstandardNth12gt33 = INITVALUE;
    CCTK_REAL PDstandardNth13gt33 = INITVALUE;
    CCTK_REAL PDstandardNth23gt33 = INITVALUE;
    CCTK_REAL PDstandardNth1phi = INITVALUE;
    CCTK_REAL PDstandardNth2phi = INITVALUE;
    CCTK_REAL PDstandardNth3phi = INITVALUE;
    CCTK_REAL PDstandardNth11phi = INITVALUE;
    CCTK_REAL PDstandardNth22phi = INITVALUE;
    CCTK_REAL PDstandardNth33phi = INITVALUE;
    CCTK_REAL PDstandardNth12phi = INITVALUE;
    CCTK_REAL PDstandardNth13phi = INITVALUE;
    CCTK_REAL PDstandardNth23phi = INITVALUE;
    CCTK_REAL PDstandardNth1trK = INITVALUE;
    CCTK_REAL PDstandardNth2trK = INITVALUE;
    CCTK_REAL PDstandardNth3trK = INITVALUE;
    CCTK_REAL PDstandardNth1Xt1 = INITVALUE;
    CCTK_REAL PDstandardNth2Xt1 = INITVALUE;
    CCTK_REAL PDstandardNth3Xt1 = INITVALUE;
    CCTK_REAL PDstandardNth1Xt2 = INITVALUE;
    CCTK_REAL PDstandardNth2Xt2 = INITVALUE;
    CCTK_REAL PDstandardNth3Xt2 = INITVALUE;
    CCTK_REAL PDstandardNth1Xt3 = INITVALUE;
    CCTK_REAL PDstandardNth2Xt3 = INITVALUE;
    CCTK_REAL PDstandardNth3Xt3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    At11L = At11[index];
    At12L = At12[index];
    At13L = At13[index];
    At22L = At22[index];
    At23L = At23[index];
    At33L = At33[index];
    dJ111L = dJ111[index];
    dJ112L = dJ112[index];
    dJ113L = dJ113[index];
    dJ122L = dJ122[index];
    dJ123L = dJ123[index];
    dJ133L = dJ133[index];
    dJ211L = dJ211[index];
    dJ212L = dJ212[index];
    dJ213L = dJ213[index];
    dJ222L = dJ222[index];
    dJ223L = dJ223[index];
    dJ233L = dJ233[index];
    dJ311L = dJ311[index];
    dJ312L = dJ312[index];
    dJ313L = dJ313[index];
    dJ322L = dJ322[index];
    dJ323L = dJ323[index];
    dJ333L = dJ333[index];
    gt11L = gt11[index];
    gt12L = gt12[index];
    gt13L = gt13[index];
    gt22L = gt22[index];
    gt23L = gt23[index];
    gt33L = gt33[index];
    J11L = J11[index];
    J12L = J12[index];
    J13L = J13[index];
    J21L = J21[index];
    J22L = J22[index];
    J23L = J23[index];
    J31L = J31[index];
    J32L = J32[index];
    J33L = J33[index];
    phiL = phi[index];
    trKL = trK[index];
    Xt1L = Xt1[index];
    Xt2L = Xt2[index];
    Xt3L = Xt3[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandardNth1At11 = PDstandardNth1(At11, i, j, k);
    PDstandardNth2At11 = PDstandardNth2(At11, i, j, k);
    PDstandardNth3At11 = PDstandardNth3(At11, i, j, k);
    PDstandardNth1At12 = PDstandardNth1(At12, i, j, k);
    PDstandardNth2At12 = PDstandardNth2(At12, i, j, k);
    PDstandardNth3At12 = PDstandardNth3(At12, i, j, k);
    PDstandardNth1At13 = PDstandardNth1(At13, i, j, k);
    PDstandardNth2At13 = PDstandardNth2(At13, i, j, k);
    PDstandardNth3At13 = PDstandardNth3(At13, i, j, k);
    PDstandardNth1At22 = PDstandardNth1(At22, i, j, k);
    PDstandardNth2At22 = PDstandardNth2(At22, i, j, k);
    PDstandardNth3At22 = PDstandardNth3(At22, i, j, k);
    PDstandardNth1At23 = PDstandardNth1(At23, i, j, k);
    PDstandardNth2At23 = PDstandardNth2(At23, i, j, k);
    PDstandardNth3At23 = PDstandardNth3(At23, i, j, k);
    PDstandardNth1At33 = PDstandardNth1(At33, i, j, k);
    PDstandardNth2At33 = PDstandardNth2(At33, i, j, k);
    PDstandardNth3At33 = PDstandardNth3(At33, i, j, k);
    PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    PDstandardNth11gt11 = PDstandardNth11(gt11, i, j, k);
    PDstandardNth22gt11 = PDstandardNth22(gt11, i, j, k);
    PDstandardNth33gt11 = PDstandardNth33(gt11, i, j, k);
    PDstandardNth12gt11 = PDstandardNth12(gt11, i, j, k);
    PDstandardNth13gt11 = PDstandardNth13(gt11, i, j, k);
    PDstandardNth23gt11 = PDstandardNth23(gt11, i, j, k);
    PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    PDstandardNth11gt12 = PDstandardNth11(gt12, i, j, k);
    PDstandardNth22gt12 = PDstandardNth22(gt12, i, j, k);
    PDstandardNth33gt12 = PDstandardNth33(gt12, i, j, k);
    PDstandardNth12gt12 = PDstandardNth12(gt12, i, j, k);
    PDstandardNth13gt12 = PDstandardNth13(gt12, i, j, k);
    PDstandardNth23gt12 = PDstandardNth23(gt12, i, j, k);
    PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    PDstandardNth11gt13 = PDstandardNth11(gt13, i, j, k);
    PDstandardNth22gt13 = PDstandardNth22(gt13, i, j, k);
    PDstandardNth33gt13 = PDstandardNth33(gt13, i, j, k);
    PDstandardNth12gt13 = PDstandardNth12(gt13, i, j, k);
    PDstandardNth13gt13 = PDstandardNth13(gt13, i, j, k);
    PDstandardNth23gt13 = PDstandardNth23(gt13, i, j, k);
    PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    PDstandardNth11gt22 = PDstandardNth11(gt22, i, j, k);
    PDstandardNth22gt22 = PDstandardNth22(gt22, i, j, k);
    PDstandardNth33gt22 = PDstandardNth33(gt22, i, j, k);
    PDstandardNth12gt22 = PDstandardNth12(gt22, i, j, k);
    PDstandardNth13gt22 = PDstandardNth13(gt22, i, j, k);
    PDstandardNth23gt22 = PDstandardNth23(gt22, i, j, k);
    PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    PDstandardNth11gt23 = PDstandardNth11(gt23, i, j, k);
    PDstandardNth22gt23 = PDstandardNth22(gt23, i, j, k);
    PDstandardNth33gt23 = PDstandardNth33(gt23, i, j, k);
    PDstandardNth12gt23 = PDstandardNth12(gt23, i, j, k);
    PDstandardNth13gt23 = PDstandardNth13(gt23, i, j, k);
    PDstandardNth23gt23 = PDstandardNth23(gt23, i, j, k);
    PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    PDstandardNth11gt33 = PDstandardNth11(gt33, i, j, k);
    PDstandardNth22gt33 = PDstandardNth22(gt33, i, j, k);
    PDstandardNth33gt33 = PDstandardNth33(gt33, i, j, k);
    PDstandardNth12gt33 = PDstandardNth12(gt33, i, j, k);
    PDstandardNth13gt33 = PDstandardNth13(gt33, i, j, k);
    PDstandardNth23gt33 = PDstandardNth23(gt33, i, j, k);
    PDstandardNth1phi = PDstandardNth1(phi, i, j, k);
    PDstandardNth2phi = PDstandardNth2(phi, i, j, k);
    PDstandardNth3phi = PDstandardNth3(phi, i, j, k);
    PDstandardNth11phi = PDstandardNth11(phi, i, j, k);
    PDstandardNth22phi = PDstandardNth22(phi, i, j, k);
    PDstandardNth33phi = PDstandardNth33(phi, i, j, k);
    PDstandardNth12phi = PDstandardNth12(phi, i, j, k);
    PDstandardNth13phi = PDstandardNth13(phi, i, j, k);
    PDstandardNth23phi = PDstandardNth23(phi, i, j, k);
    PDstandardNth1trK = PDstandardNth1(trK, i, j, k);
    PDstandardNth2trK = PDstandardNth2(trK, i, j, k);
    PDstandardNth3trK = PDstandardNth3(trK, i, j, k);
    PDstandardNth1Xt1 = PDstandardNth1(Xt1, i, j, k);
    PDstandardNth2Xt1 = PDstandardNth2(Xt1, i, j, k);
    PDstandardNth3Xt1 = PDstandardNth3(Xt1, i, j, k);
    PDstandardNth1Xt2 = PDstandardNth1(Xt2, i, j, k);
    PDstandardNth2Xt2 = PDstandardNth2(Xt2, i, j, k);
    PDstandardNth3Xt2 = PDstandardNth3(Xt2, i, j, k);
    PDstandardNth1Xt3 = PDstandardNth1(Xt3, i, j, k);
    PDstandardNth2Xt3 = PDstandardNth2(Xt3, i, j, k);
    PDstandardNth3Xt3 = PDstandardNth3(Xt3, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    detgt  =  1;
    
    gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    gtu21  =  (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    gtu31  =  (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    gtu32  =  (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    Gt111  =  khalf*((gtu11*J11L - gtu21*J12L - gtu31*J13L)*PDstandardNth1gt11 + 
          (gtu11*J21L - gtu21*J22L - gtu31*J23L)*PDstandardNth2gt11 + 
          (gtu11*J31L - gtu21*J32L - gtu31*J33L)*PDstandardNth3gt11 + 
          2*(J11L*(gtu21*PDstandardNth1gt12 + gtu31*PDstandardNth1gt13) + 
             J21L*(gtu21*PDstandardNth2gt12 + gtu31*PDstandardNth2gt13) + 
             J31L*(gtu21*PDstandardNth3gt12 + gtu31*PDstandardNth3gt13)));
    
    Gt211  =  khalf*((gtu21*J11L - gtu22*J12L - gtu32*J13L)*PDstandardNth1gt11 + 
          (gtu21*J21L - gtu22*J22L - gtu32*J23L)*PDstandardNth2gt11 + 
          (gtu21*J31L - gtu22*J32L - gtu32*J33L)*PDstandardNth3gt11 + 
          2*(J11L*(gtu22*PDstandardNth1gt12 + gtu32*PDstandardNth1gt13) + 
             J21L*(gtu22*PDstandardNth2gt12 + gtu32*PDstandardNth2gt13) + 
             J31L*(gtu22*PDstandardNth3gt12 + gtu32*PDstandardNth3gt13)));
    
    Gt311  =  khalf*((gtu31*J11L - gtu32*J12L - gtu33*J13L)*PDstandardNth1gt11 + 
          (gtu31*J21L - gtu32*J22L - gtu33*J23L)*PDstandardNth2gt11 + 
          (gtu31*J31L - gtu32*J32L - gtu33*J33L)*PDstandardNth3gt11 + 
          2*(J11L*(gtu32*PDstandardNth1gt12 + gtu33*PDstandardNth1gt13) + 
             J21L*(gtu32*PDstandardNth2gt12 + gtu33*PDstandardNth2gt13) + 
             J31L*(gtu32*PDstandardNth3gt12 + gtu33*PDstandardNth3gt13)));
    
    Gt112  =  khalf*(gtu11*(J12L*PDstandardNth1gt11 + J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
          gtu21*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22) + 
          gtu31*(-(J13L*PDstandardNth1gt12) + J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23));
    
    Gt212  =  khalf*(gtu21*(J12L*PDstandardNth1gt11 + J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
          gtu22*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22) + 
          gtu32*(-(J13L*PDstandardNth1gt12) + J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23));
    
    Gt312  =  khalf*(gtu31*(J12L*PDstandardNth1gt11 + J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
          gtu32*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22) + 
          gtu33*(-(J13L*PDstandardNth1gt12) + J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23));
    
    Gt113  =  khalf*(gtu11*(J13L*PDstandardNth1gt11 + J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
          gtu21*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23) + gtu31*(J11L*PDstandardNth1gt33 + J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    Gt213  =  khalf*(gtu21*(J13L*PDstandardNth1gt11 + J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
          gtu22*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23) + gtu32*(J11L*PDstandardNth1gt33 + J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    Gt313  =  khalf*(gtu31*(J13L*PDstandardNth1gt11 + J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
          gtu32*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
             J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
             J31L*PDstandardNth3gt23) + gtu33*(J11L*PDstandardNth1gt33 + J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    Gt122  =  khalf*(gtu11*(-(J11L*PDstandardNth1gt22) + 2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
             J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - J31L*PDstandardNth3gt22) + 
          gtu21*(J12L*PDstandardNth1gt22 + J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
          gtu31*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22 - 
             2*(J12L*PDstandardNth1gt23 + J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    Gt222  =  khalf*(gtu21*(-(J11L*PDstandardNth1gt22) + 2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
             J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - J31L*PDstandardNth3gt22) + 
          gtu22*(J12L*PDstandardNth1gt22 + J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
          gtu32*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22 - 
             2*(J12L*PDstandardNth1gt23 + J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    Gt322  =  khalf*(gtu31*(-(J11L*PDstandardNth1gt22) + 2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
             J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - J31L*PDstandardNth3gt22) + 
          gtu32*(J12L*PDstandardNth1gt22 + J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
          gtu33*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22 - 
             2*(J12L*PDstandardNth1gt23 + J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    Gt123  =  khalf*(gtu21*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
          gtu11*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt23) + gtu31*(J12L*PDstandardNth1gt33 + J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    Gt223  =  khalf*(gtu22*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
          gtu21*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt23) + gtu32*(J12L*PDstandardNth1gt33 + J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    Gt323  =  khalf*(gtu32*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
          gtu31*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
             J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
             J31L*PDstandardNth3gt23) + gtu33*(J12L*PDstandardNth1gt33 + J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    Gt133  =  khalf*(gtu11*(-(J11L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
             J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - J31L*PDstandardNth3gt33) + 
          gtu21*(-(J12L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
             J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - J32L*PDstandardNth3gt33) + 
          gtu31*(J13L*PDstandardNth1gt33 + J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    Gt233  =  khalf*(gtu21*(-(J11L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
             J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - J31L*PDstandardNth3gt33) + 
          gtu22*(-(J12L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
             J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - J32L*PDstandardNth3gt33) + 
          gtu32*(J13L*PDstandardNth1gt33 + J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    Gt333  =  khalf*(gtu31*(-(J11L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
             J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - J31L*PDstandardNth3gt33) + 
          gtu32*(-(J12L*PDstandardNth1gt33) + 2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
             J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - J32L*PDstandardNth3gt33) + 
          gtu33*(J13L*PDstandardNth1gt33 + J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    Rt11  =  (Gt113*(gt13L*Gt312 + 3*(gt12L*Gt212 + gt13L*Gt312)) + 
           gt12L*(Gt213*(4*Gt112 + 2*Gt222) + Gt212*(Gt113 + 2*Gt223) + 2*(Gt233*Gt312 + Gt223*Gt313)) + 
           gt11L*(6*Gt112*Gt113 + 2*(Gt122*Gt213 + Gt133*Gt312 + Gt123*(Gt212 + Gt313))) + 
           gt13L*(2*Gt213*Gt322 + Gt313*(4*Gt112 + 2*Gt323)) + 
           2*(Gt213*(Gt212*gt22L + gt23L*Gt312) + Gt212*(gt23L*Gt313 + gt13L*Gt323) + Gt312*(gt13L*Gt333 + Gt313*gt33L)))*
         gtu32 + J11L*(gt11L*PDstandardNth1Xt1 + gt12L*PDstandardNth1Xt2 + gt13L*PDstandardNth1Xt3) + 
        J21L*(gt11L*PDstandardNth2Xt1 + gt12L*PDstandardNth2Xt2 + gt13L*PDstandardNth2Xt3) + 
        J31L*(gt11L*PDstandardNth3Xt1 + gt12L*PDstandardNth3Xt2 + gt13L*PDstandardNth3Xt3) + 
        (Gt111*gt11L + gt12L*Gt211 + gt13L*Gt311)*Xt1L + (Gt112*gt11L + gt12L*Gt212 + gt13L*Gt312)*Xt2L + 
        (Gt113*gt11L + gt12L*Gt213 + gt13L*Gt313)*Xt3L + 
        gtu21*(Gt112*(gt13L*Gt311 + 3*(gt12L*Gt211 + gt13L*Gt311)) + 
           gt11L*(Gt112*(6*Gt111 + 2*Gt212) + 2*(Gt122*Gt211 + Gt123*Gt311 + Gt113*Gt312)) + 
           2*(Gt212*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + gt13L*(Gt211*Gt322 + Gt311*Gt323)) + 
           Gt312*(gt13L*(4*Gt111 + 2*Gt313) + 2*(Gt211*gt23L + Gt311*gt33L)) + 
           gt12L*(4*Gt111*Gt212 + Gt211*(Gt112 + 2*Gt222) + 2*(Gt223*Gt311 + Gt213*Gt312 + SQR(Gt212)))) + 
        gtu11*(4*Gt111*(gt12L*Gt211 + gt13L*Gt311) + 
           2*(gt12L*(Gt211*Gt212 + Gt213*Gt311) + Gt211*(gt23L*Gt311 + gt13L*Gt312) + gt13L*Gt311*Gt313) + 
           gt11L*(2*(Gt112*Gt211 + Gt113*Gt311) + 3*SQR(Gt111)) + gt22L*SQR(Gt211) + gt33L*SQR(Gt311)) + 
        gtu22*(4*Gt112*(gt12L*Gt212 + gt13L*Gt312) + 
           2*(gt12L*(Gt212*Gt222 + Gt223*Gt312) + Gt212*(gt23L*Gt312 + gt13L*Gt322) + gt13L*Gt312*Gt323) + 
           gt11L*(2*(Gt122*Gt212 + Gt123*Gt312) + 3*SQR(Gt112)) + gt22L*SQR(Gt212) + gt33L*SQR(Gt312)) + 
        gtu33*(4*Gt113*(gt12L*Gt213 + gt13L*Gt313) + 
           2*(gt12L*(Gt213*Gt223 + Gt233*Gt313) + Gt213*(gt23L*Gt313 + gt13L*Gt323) + gt13L*Gt313*Gt333) + 
           gt11L*(2*(Gt123*Gt213 + Gt133*Gt313) + 3*SQR(Gt113)) + gt22L*SQR(Gt213) + gt33L*SQR(Gt313)) + 
        gtu31*(Gt113*(gt13L*Gt311 + 3*(gt12L*Gt211 + gt13L*Gt311)) + 
           gt11L*(2*(Gt123*Gt211 + Gt112*Gt213 + Gt133*Gt311) + Gt113*(6*Gt111 + 2*Gt313)) + 
           gt12L*(Gt211*(Gt113 + 2*Gt223) + 2*Gt233*Gt311 + Gt213*(4*Gt111 + 2*(Gt212 + Gt313))) + 
           2*(Gt213*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + gt13L*Gt311*Gt333 + Gt313*(Gt211*gt23L + Gt311*gt33L)) + 
           gt13L*(4*Gt111*Gt313 + 2*(Gt211*Gt323 + SQR(Gt313)))) + 
        khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt11 + J12L*J21L*PDstandardNth12gt11 + J11L*J22L*PDstandardNth12gt11 + 
                 J12L*J31L*PDstandardNth13gt11 + J11L*J32L*PDstandardNth13gt11 + dJ112L*PDstandardNth1gt11 + 
                 J21L*J22L*PDstandardNth22gt11 + J22L*J31L*PDstandardNth23gt11 + J21L*J32L*PDstandardNth23gt11 + 
                 dJ212L*PDstandardNth2gt11 + J31L*J32L*PDstandardNth33gt11 + dJ312L*PDstandardNth3gt11) + 
              gtu31*(J11L*J13L*PDstandardNth11gt11 + J13L*J21L*PDstandardNth12gt11 + J11L*J23L*PDstandardNth12gt11 + 
                 J13L*J31L*PDstandardNth13gt11 + J11L*J33L*PDstandardNth13gt11 + dJ113L*PDstandardNth1gt11 + 
                 J21L*J23L*PDstandardNth22gt11 + J23L*J31L*PDstandardNth23gt11 + J21L*J33L*PDstandardNth23gt11 + 
                 dJ213L*PDstandardNth2gt11 + J31L*J33L*PDstandardNth33gt11 + dJ313L*PDstandardNth3gt11) + 
              gtu32*(J12L*J13L*PDstandardNth11gt11 + J13L*J22L*PDstandardNth12gt11 + J12L*J23L*PDstandardNth12gt11 + 
                 J13L*J32L*PDstandardNth13gt11 + J12L*J33L*PDstandardNth13gt11 + dJ123L*PDstandardNth1gt11 + 
                 J22L*J23L*PDstandardNth22gt11 + J23L*J32L*PDstandardNth23gt11 + J22L*J33L*PDstandardNth23gt11 + 
                 dJ223L*PDstandardNth2gt11 + J32L*J33L*PDstandardNth33gt11 + dJ323L*PDstandardNth3gt11)) - 
           gtu11*(2*J11L*J21L*PDstandardNth12gt11 + 2*J11L*J31L*PDstandardNth13gt11 + dJ111L*PDstandardNth1gt11 + 
              2*J21L*J31L*PDstandardNth23gt11 + dJ211L*PDstandardNth2gt11 + dJ311L*PDstandardNth3gt11 + 
              PDstandardNth11gt11*SQR(J11L) + PDstandardNth22gt11*SQR(J21L) + PDstandardNth33gt11*SQR(J31L)) - 
           gtu22*(2*J12L*J22L*PDstandardNth12gt11 + 2*J12L*J32L*PDstandardNth13gt11 + dJ122L*PDstandardNth1gt11 + 
              2*J22L*J32L*PDstandardNth23gt11 + dJ222L*PDstandardNth2gt11 + dJ322L*PDstandardNth3gt11 + 
              PDstandardNth11gt11*SQR(J12L) + PDstandardNth22gt11*SQR(J22L) + PDstandardNth33gt11*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt11 + 2*J13L*J33L*PDstandardNth13gt11 + dJ133L*PDstandardNth1gt11 + 
              2*J23L*J33L*PDstandardNth23gt11 + dJ233L*PDstandardNth2gt11 + dJ333L*PDstandardNth3gt11 + 
              PDstandardNth11gt11*SQR(J13L) + PDstandardNth22gt11*SQR(J23L) + PDstandardNth33gt11*SQR(J33L)));
    
    Rt12  =  khalf*((gt12L*J11L + gt11L*J12L)*PDstandardNth1Xt1 + (gt22L*J11L + gt12L*J12L)*PDstandardNth1Xt2 + 
          (gt23L*J11L + gt13L*J12L)*PDstandardNth1Xt3 + (gt12L*J21L + gt11L*J22L)*PDstandardNth2Xt1 + 
          (gt22L*J21L + gt12L*J22L)*PDstandardNth2Xt2 + (gt23L*J21L + gt13L*J22L)*PDstandardNth2Xt3 - 
          2*(gtu21*(J11L*J12L*PDstandardNth11gt12 + J12L*J21L*PDstandardNth12gt12 + J11L*J22L*PDstandardNth12gt12 + 
                J12L*J31L*PDstandardNth13gt12 + J11L*J32L*PDstandardNth13gt12 + dJ112L*PDstandardNth1gt12 + 
                J21L*J22L*PDstandardNth22gt12 + J22L*J31L*PDstandardNth23gt12 + J21L*J32L*PDstandardNth23gt12 + 
                dJ212L*PDstandardNth2gt12 + J31L*J32L*PDstandardNth33gt12 + dJ312L*PDstandardNth3gt12) + 
             gtu31*(J11L*J13L*PDstandardNth11gt12 + J13L*J21L*PDstandardNth12gt12 + J11L*J23L*PDstandardNth12gt12 + 
                J13L*J31L*PDstandardNth13gt12 + J11L*J33L*PDstandardNth13gt12 + dJ113L*PDstandardNth1gt12 + 
                J21L*J23L*PDstandardNth22gt12 + J23L*J31L*PDstandardNth23gt12 + J21L*J33L*PDstandardNth23gt12 + 
                dJ213L*PDstandardNth2gt12 + J31L*J33L*PDstandardNth33gt12 + dJ313L*PDstandardNth3gt12) + 
             gtu32*(J12L*J13L*PDstandardNth11gt12 + J13L*J22L*PDstandardNth12gt12 + J12L*J23L*PDstandardNth12gt12 + 
                J13L*J32L*PDstandardNth13gt12 + J12L*J33L*PDstandardNth13gt12 + dJ123L*PDstandardNth1gt12 + 
                J22L*J23L*PDstandardNth22gt12 + J23L*J32L*PDstandardNth23gt12 + J22L*J33L*PDstandardNth23gt12 + 
                dJ223L*PDstandardNth2gt12 + J32L*J33L*PDstandardNth33gt12 + dJ323L*PDstandardNth3gt12)) + 
          (gt12L*J31L + gt11L*J32L)*PDstandardNth3Xt1 + (gt22L*J31L + gt12L*J32L)*PDstandardNth3Xt2 + 
          (gt23L*J31L + gt13L*J32L)*PDstandardNth3Xt3 + 
          (Gt112*gt11L + Gt111*gt12L + gt12L*Gt212 + Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)*Xt1L + 
          (gt11L*Gt122 + Gt112*gt12L + gt12L*Gt222 + Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322)*Xt2L + 
          (gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + Gt213*gt22L + gt23L*Gt313 + gt13L*Gt323)*Xt3L + 
          2*((Gt123*gt12L*Gt211 + Gt113*gt12L*Gt212 + 2*Gt112*gt12L*Gt213 + gt12L*Gt212*Gt223 + Gt212*Gt213*gt22L + 
                Gt211*Gt223*gt22L + gt12L*Gt133*Gt311 + gt22L*Gt233*Gt311 + Gt113*gt13L*Gt312 + gt12L*Gt233*Gt312 + 
                Gt213*gt23L*Gt312 + gt11L*(2*Gt112*Gt113 + Gt123*Gt212 + Gt133*Gt312) + 2*Gt112*gt13L*Gt313 + 
                Gt212*gt23L*Gt313 + Gt111*(Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + gt13L*Gt212*Gt323 + 
                Gt211*gt23L*Gt323 + gt23L*Gt311*Gt333 + gt13L*Gt312*Gt333 + Gt312*Gt313*gt33L)*gtu31 + 
             (Gt123*gt12L*Gt212 + 2*Gt122*gt12L*Gt213 + Gt113*gt12L*Gt222 + gt12L*Gt222*Gt223 + Gt213*Gt222*gt22L + 
                Gt212*Gt223*gt22L + gt12L*Gt133*Gt312 + gt22L*Gt233*Gt312 + 2*Gt122*gt13L*Gt313 + Gt222*gt23L*Gt313 + 
                Gt112*(Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + Gt113*gt13L*Gt322 + gt12L*Gt233*Gt322 + 
                Gt213*gt23L*Gt322 + gt11L*(2*Gt113*Gt122 + Gt123*Gt222 + Gt133*Gt322) + gt13L*Gt222*Gt323 + 
                Gt212*gt23L*Gt323 + gt23L*Gt312*Gt333 + gt13L*Gt322*Gt333 + Gt313*Gt322*gt33L)*gtu32 + 
             gtu11*(3*Gt112*gt12L*Gt211 + 2*Gt211*Gt212*gt22L + Gt113*gt12L*Gt311 + 2*Gt112*gt13L*Gt311 + 
                Gt213*gt22L*Gt311 + Gt212*gt23L*Gt311 + gt13L*Gt212*Gt312 + gt12L*Gt213*Gt312 + 2*Gt211*gt23L*Gt312 + 
                gt11L*(2*Gt111*Gt112 + Gt112*Gt212 + Gt113*Gt312) + 
                Gt111*(gt12L*Gt212 + Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + gt23L*Gt311*Gt313 + gt13L*Gt312*Gt313 + 
                Gt311*Gt312*gt33L + gt12L*SQR(Gt111) + gt12L*SQR(Gt212)) + 
             gtu21*(Gt112*gt11L*Gt222 + Gt112*Gt211*gt22L + Gt211*Gt222*gt22L + 2*Gt122*gt13L*Gt311 + Gt112*gt23L*Gt311 + 
                Gt222*gt23L*Gt311 + gt13L*Gt222*Gt312 + Gt213*gt22L*Gt312 + Gt212*gt23L*Gt312 + gt23L*Gt312*Gt313 + 
                Gt113*gt11L*Gt322 + Gt211*gt23L*Gt322 + gt13L*Gt313*Gt322 + 
                Gt111*(2*gt11L*Gt122 + Gt112*gt12L + gt12L*Gt222 + gt13L*Gt322) + 
                gt12L*(2*Gt122*Gt211 + Gt112*Gt212 + Gt212*Gt222 + Gt113*Gt312 + Gt213*Gt322) + Gt311*Gt322*gt33L + 
                gt22L*SQR(Gt212)) + gtu22*(gt11L*Gt122*Gt222 + 2*Gt212*Gt222*gt22L + 2*Gt122*gt13L*Gt312 + 
                Gt223*gt22L*Gt312 + Gt222*gt23L*Gt312 + gt11L*Gt123*Gt322 + gt13L*Gt222*Gt322 + 2*Gt212*gt23L*Gt322 + 
                Gt112*(2*gt11L*Gt122 + gt12L*Gt222 + Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322) + gt23L*Gt312*Gt323 + 
                gt13L*Gt322*Gt323 + Gt312*Gt322*gt33L + gt12L*SQR(Gt112) + 
                gt12L*(3*Gt122*Gt212 + Gt123*Gt312 + Gt223*Gt322 + SQR(Gt222))) + 
             gtu33*(gt11L*Gt123*Gt223 + 2*Gt213*Gt223*gt22L + 2*Gt123*gt13L*Gt313 + gt22L*Gt233*Gt313 + Gt223*gt23L*Gt313 + 
                gt11L*Gt133*Gt323 + gt13L*Gt223*Gt323 + 2*Gt213*gt23L*Gt323 + 
                Gt113*(2*gt11L*Gt123 + gt12L*Gt223 + Gt213*gt22L + gt23L*Gt313 + gt13L*Gt323) + gt23L*Gt313*Gt333 + 
                gt13L*Gt323*Gt333 + Gt313*Gt323*gt33L + gt12L*SQR(Gt113) + 
                gt12L*(3*Gt123*Gt213 + Gt133*Gt313 + Gt233*Gt323 + SQR(Gt223))) + 
             gtu21*(Gt122*gt12L*Gt211 + 3*Gt112*gt12L*Gt212 + gt12L*Gt212*Gt222 + Gt211*Gt222*gt22L + Gt123*gt12L*Gt311 + 
                Gt223*gt22L*Gt311 + 3*Gt112*gt13L*Gt312 + gt12L*Gt223*Gt312 + 2*Gt212*gt23L*Gt312 + 
                Gt111*(Gt112*gt12L + Gt212*gt22L + gt23L*Gt312) + gt13L*Gt212*Gt322 + Gt211*gt23L*Gt322 + 
                gt23L*Gt311*Gt323 + gt13L*Gt312*Gt323 + gt11L*(Gt122*Gt212 + Gt123*Gt312 + 2*SQR(Gt112)) + 
                gt22L*SQR(Gt212) + gt33L*SQR(Gt312)) + 
             gtu31*(Gt112*gt11L*Gt223 + Gt113*Gt211*gt22L + Gt212*Gt213*gt22L + Gt211*Gt223*gt22L + 2*Gt123*gt13L*Gt311 + 
                Gt113*gt23L*Gt311 + Gt223*gt23L*Gt311 + gt13L*Gt223*Gt312 + Gt213*gt23L*Gt312 + Gt213*gt22L*Gt313 + 
                Gt113*gt11L*Gt323 + Gt211*gt23L*Gt323 + gt13L*Gt313*Gt323 + 
                Gt111*(2*gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + gt13L*Gt323) + 
                gt12L*(2*Gt123*Gt211 + Gt112*Gt213 + Gt212*Gt223 + Gt113*Gt313 + Gt213*Gt323) + Gt311*Gt323*gt33L + 
                gt23L*SQR(Gt313)) + gtu32*(gt11L*Gt122*Gt223 + Gt113*Gt212*gt22L + Gt213*Gt222*gt22L + Gt212*Gt223*gt22L + 
                2*Gt123*gt13L*Gt312 + Gt113*gt23L*Gt312 + Gt223*gt23L*Gt312 + Gt223*gt22L*Gt313 + gt13L*Gt223*Gt322 + 
                Gt213*gt23L*Gt322 + gt11L*Gt123*Gt323 + Gt212*gt23L*Gt323 + gt23L*Gt313*Gt323 + 
                Gt112*(2*gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + gt13L*Gt323) + 
                gt12L*(Gt122*Gt213 + Gt123*(2*Gt212 + Gt313) + Gt223*(Gt222 + Gt323)) + Gt312*Gt323*gt33L + gt13L*SQR(Gt323)
                )) - gtu11*(2*J11L*J21L*PDstandardNth12gt12 + 2*J11L*J31L*PDstandardNth13gt12 + dJ111L*PDstandardNth1gt12 + 
             2*J21L*J31L*PDstandardNth23gt12 + dJ211L*PDstandardNth2gt12 + dJ311L*PDstandardNth3gt12 + 
             PDstandardNth11gt12*SQR(J11L) + PDstandardNth22gt12*SQR(J21L) + PDstandardNth33gt12*SQR(J31L)) - 
          gtu22*(2*J12L*J22L*PDstandardNth12gt12 + 2*J12L*J32L*PDstandardNth13gt12 + dJ122L*PDstandardNth1gt12 + 
             2*J22L*J32L*PDstandardNth23gt12 + dJ222L*PDstandardNth2gt12 + dJ322L*PDstandardNth3gt12 + 
             PDstandardNth11gt12*SQR(J12L) + PDstandardNth22gt12*SQR(J22L) + PDstandardNth33gt12*SQR(J32L)) - 
          gtu33*(2*J13L*J23L*PDstandardNth12gt12 + 2*J13L*J33L*PDstandardNth13gt12 + dJ133L*PDstandardNth1gt12 + 
             2*J23L*J33L*PDstandardNth23gt12 + dJ233L*PDstandardNth2gt12 + dJ333L*PDstandardNth3gt12 + 
             PDstandardNth11gt12*SQR(J13L) + PDstandardNth22gt12*SQR(J23L) + PDstandardNth33gt12*SQR(J33L)));
    
    Rt13  =  khalf*((gt13L*J11L + gt11L*J13L)*PDstandardNth1Xt1 + (gt23L*J11L + gt12L*J13L)*PDstandardNth1Xt2 + 
          (gt33L*J11L + gt13L*J13L)*PDstandardNth1Xt3 + (gt13L*J21L + gt11L*J23L)*PDstandardNth2Xt1 + 
          (gt23L*J21L + gt12L*J23L)*PDstandardNth2Xt2 + (gt33L*J21L + gt13L*J23L)*PDstandardNth2Xt3 - 
          2*(gtu21*(J11L*J12L*PDstandardNth11gt13 + J12L*J21L*PDstandardNth12gt13 + J11L*J22L*PDstandardNth12gt13 + 
                J12L*J31L*PDstandardNth13gt13 + J11L*J32L*PDstandardNth13gt13 + dJ112L*PDstandardNth1gt13 + 
                J21L*J22L*PDstandardNth22gt13 + J22L*J31L*PDstandardNth23gt13 + J21L*J32L*PDstandardNth23gt13 + 
                dJ212L*PDstandardNth2gt13 + J31L*J32L*PDstandardNth33gt13 + dJ312L*PDstandardNth3gt13) + 
             gtu31*(J11L*J13L*PDstandardNth11gt13 + J13L*J21L*PDstandardNth12gt13 + J11L*J23L*PDstandardNth12gt13 + 
                J13L*J31L*PDstandardNth13gt13 + J11L*J33L*PDstandardNth13gt13 + dJ113L*PDstandardNth1gt13 + 
                J21L*J23L*PDstandardNth22gt13 + J23L*J31L*PDstandardNth23gt13 + J21L*J33L*PDstandardNth23gt13 + 
                dJ213L*PDstandardNth2gt13 + J31L*J33L*PDstandardNth33gt13 + dJ313L*PDstandardNth3gt13) + 
             gtu32*(J12L*J13L*PDstandardNth11gt13 + J13L*J22L*PDstandardNth12gt13 + J12L*J23L*PDstandardNth12gt13 + 
                J13L*J32L*PDstandardNth13gt13 + J12L*J33L*PDstandardNth13gt13 + dJ123L*PDstandardNth1gt13 + 
                J22L*J23L*PDstandardNth22gt13 + J23L*J32L*PDstandardNth23gt13 + J22L*J33L*PDstandardNth23gt13 + 
                dJ223L*PDstandardNth2gt13 + J32L*J33L*PDstandardNth33gt13 + dJ323L*PDstandardNth3gt13)) + 
          (gt13L*J31L + gt11L*J33L)*PDstandardNth3Xt1 + (gt23L*J31L + gt12L*J33L)*PDstandardNth3Xt2 + 
          (gt33L*J31L + gt13L*J33L)*PDstandardNth3Xt3 + 
          (Gt113*gt11L + Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L + gt13L*Gt313 + Gt311*gt33L)*Xt1L + 
          (gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + gt13L*Gt323 + Gt312*gt33L)*Xt2L + 
          (gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + Gt213*gt23L + gt13L*Gt333 + Gt313*gt33L)*Xt3L + 
          2*((Gt122*gt13L*Gt211 + 2*Gt113*gt12L*Gt212 + Gt112*gt12L*Gt213 + gt12L*Gt213*Gt222 + Gt212*Gt213*gt22L + 
                Gt211*Gt222*gt23L + Gt123*gt13L*Gt311 + Gt223*gt23L*Gt311 + 2*Gt113*gt13L*Gt312 + Gt213*gt23L*Gt312 + 
                Gt112*gt13L*Gt313 + gt12L*Gt223*Gt313 + Gt212*gt23L*Gt313 + 
                gt11L*(2*Gt112*Gt113 + Gt122*Gt213 + Gt123*Gt313) + gt13L*Gt213*Gt322 + gt13L*Gt313*Gt323 + 
                Gt312*Gt313*gt33L + Gt211*Gt322*gt33L + Gt311*Gt323*gt33L + Gt111*(Gt112*gt13L + Gt212*gt23L + Gt312*gt33L))
               *gtu21 + (Gt122*gt13L*Gt213 + gt11L*Gt122*Gt233 + Gt212*gt22L*Gt233 + Gt113*Gt212*gt23L + 
                Gt213*Gt222*gt23L + 2*Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + Gt123*gt13L*Gt313 + Gt223*gt23L*Gt313 + 
                gt13L*Gt233*Gt322 + gt11L*Gt123*Gt333 + Gt212*gt23L*Gt333 + gt13L*Gt323*Gt333 + 
                Gt112*(2*gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + gt13L*Gt333) + 
                gt12L*(2*Gt133*Gt212 + Gt222*Gt233 + Gt223*Gt333) + Gt113*Gt312*gt33L + Gt213*Gt322*gt33L + 
                Gt313*Gt323*gt33L + Gt312*Gt333*gt33L)*gtu32 + 
             gtu21*(2*Gt123*gt12L*Gt211 + Gt112*gt13L*Gt212 + gt12L*Gt212*Gt223 + Gt211*Gt223*gt22L + Gt112*Gt211*gt23L + 
                2*Gt123*gt13L*Gt311 + Gt223*gt23L*Gt311 + Gt113*gt13L*Gt312 + gt13L*Gt223*Gt312 + Gt213*gt23L*Gt312 + 
                gt12L*Gt213*Gt323 + Gt211*gt23L*Gt323 + gt13L*Gt313*Gt323 + 
                gt11L*(2*Gt111*Gt123 + Gt112*Gt223 + Gt113*Gt323) + Gt111*(Gt112*gt13L + gt12L*Gt223 + gt13L*Gt323) + 
                Gt112*Gt311*gt33L + Gt212*Gt312*gt33L + Gt312*Gt313*gt33L + Gt311*Gt323*gt33L + gt23L*SQR(Gt212)) + 
             gtu32*(Gt123*gt13L*Gt212 + 2*Gt123*gt12L*Gt213 + Gt113*gt12L*Gt223 + Gt213*Gt223*gt22L + Gt212*Gt223*gt23L + 
                Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + 2*Gt123*gt13L*Gt313 + Gt223*gt23L*Gt313 + Gt113*gt13L*Gt323 + 
                gt13L*Gt223*Gt323 + gt12L*Gt233*Gt323 + Gt213*gt23L*Gt323 + 
                gt11L*(2*Gt113*Gt123 + Gt123*Gt223 + Gt133*Gt323) + gt13L*Gt323*Gt333 + Gt212*Gt323*gt33L + 
                Gt313*Gt323*gt33L + Gt312*Gt333*gt33L + Gt112*(Gt113*gt13L + Gt213*gt23L + Gt313*gt33L) + gt12L*SQR(Gt223))\
              + gtu11*(2*Gt113*gt12L*Gt211 + Gt112*gt13L*Gt211 + gt12L*Gt212*Gt213 + Gt211*Gt213*gt22L + 
                Gt211*Gt212*gt23L + 3*Gt113*gt13L*Gt311 + 2*Gt213*gt23L*Gt311 + gt13L*Gt213*Gt312 + gt12L*Gt213*Gt313 + 
                Gt211*gt23L*Gt313 + gt11L*(2*Gt111*Gt113 + Gt112*Gt213 + Gt113*Gt313) + Gt211*Gt312*gt33L + 
                2*Gt311*Gt313*gt33L + Gt111*(gt12L*Gt213 + Gt211*gt23L + gt13L*Gt313 + Gt311*gt33L) + gt13L*SQR(Gt111) + 
                gt13L*SQR(Gt313)) + gtu31*(Gt112*gt13L*Gt213 + Gt112*gt11L*Gt233 + Gt211*gt22L*Gt233 + Gt113*Gt211*gt23L + 
                Gt212*Gt213*gt23L + 2*Gt133*gt13L*Gt311 + Gt233*gt23L*Gt311 + gt13L*Gt233*Gt312 + Gt113*gt13L*Gt313 + 
                Gt213*gt23L*Gt313 + Gt113*gt11L*Gt333 + Gt211*gt23L*Gt333 + gt13L*Gt313*Gt333 + 
                Gt111*(2*gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + gt13L*Gt333) + 
                gt12L*(2*Gt133*Gt211 + Gt212*Gt233 + Gt213*Gt333) + Gt113*Gt311*gt33L + Gt213*Gt312*gt33L + 
                Gt311*Gt333*gt33L + gt33L*SQR(Gt313)) + 
             gtu31*(Gt123*gt13L*Gt211 + 3*Gt113*gt12L*Gt213 + gt12L*Gt213*Gt223 + Gt211*Gt223*gt23L + Gt133*gt13L*Gt311 + 
                Gt233*gt23L*Gt311 + 3*Gt113*gt13L*Gt313 + gt12L*Gt233*Gt313 + 2*Gt213*gt23L*Gt313 + gt13L*Gt213*Gt323 + 
                gt13L*Gt313*Gt333 + Gt211*Gt323*gt33L + Gt311*Gt333*gt33L + 
                Gt111*(Gt113*gt13L + Gt213*gt23L + Gt313*gt33L) + gt11L*(Gt123*Gt213 + Gt133*Gt313 + 2*SQR(Gt113)) + 
                gt22L*SQR(Gt213) + gt33L*SQR(Gt313)) + 
             gtu22*(2*Gt123*gt12L*Gt212 + Gt122*gt13L*Gt212 + gt12L*Gt222*Gt223 + Gt212*Gt223*gt22L + Gt212*Gt222*gt23L + 
                3*Gt123*gt13L*Gt312 + 2*Gt223*gt23L*Gt312 + gt13L*Gt223*Gt322 + gt12L*Gt223*Gt323 + Gt212*gt23L*Gt323 + 
                gt11L*(2*Gt112*Gt123 + Gt122*Gt223 + Gt123*Gt323) + Gt212*Gt322*gt33L + 2*Gt312*Gt323*gt33L + 
                Gt112*(gt12L*Gt223 + Gt212*gt23L + gt13L*Gt323 + Gt312*gt33L) + gt13L*SQR(Gt112) + gt13L*SQR(Gt323)) + 
             gtu33*(2*gt12L*Gt133*Gt213 + Gt123*gt13L*Gt213 + gt11L*Gt123*Gt233 + gt12L*Gt223*Gt233 + Gt213*gt22L*Gt233 + 
                Gt213*Gt223*gt23L + 3*Gt133*gt13L*Gt313 + 2*Gt233*gt23L*Gt313 + gt13L*Gt233*Gt323 + gt11L*Gt133*Gt333 + 
                gt12L*Gt233*Gt333 + Gt213*gt23L*Gt333 + Gt213*Gt323*gt33L + 2*Gt313*Gt333*gt33L + 
                Gt113*(2*gt11L*Gt133 + gt12L*Gt233 + Gt213*gt23L + gt13L*Gt333 + Gt313*gt33L) + gt13L*SQR(Gt113) + 
                gt13L*SQR(Gt333))) - gtu11*(2*J11L*J21L*PDstandardNth12gt13 + 2*J11L*J31L*PDstandardNth13gt13 + 
             dJ111L*PDstandardNth1gt13 + 2*J21L*J31L*PDstandardNth23gt13 + dJ211L*PDstandardNth2gt13 + 
             dJ311L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J11L) + PDstandardNth22gt13*SQR(J21L) + 
             PDstandardNth33gt13*SQR(J31L)) - gtu22*(2*J12L*J22L*PDstandardNth12gt13 + 2*J12L*J32L*PDstandardNth13gt13 + 
             dJ122L*PDstandardNth1gt13 + 2*J22L*J32L*PDstandardNth23gt13 + dJ222L*PDstandardNth2gt13 + 
             dJ322L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J12L) + PDstandardNth22gt13*SQR(J22L) + 
             PDstandardNth33gt13*SQR(J32L)) - gtu33*(2*J13L*J23L*PDstandardNth12gt13 + 2*J13L*J33L*PDstandardNth13gt13 + 
             dJ133L*PDstandardNth1gt13 + 2*J23L*J33L*PDstandardNth23gt13 + dJ233L*PDstandardNth2gt13 + 
             dJ333L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J13L) + PDstandardNth22gt13*SQR(J23L) + 
             PDstandardNth33gt13*SQR(J33L)));
    
    Rt22  =  (Gt223*(3*Gt112*gt12L + 6*Gt212*gt22L + 4*gt23L*Gt312) + 
           Gt123*(Gt112*gt11L + gt12L*(2*Gt111 + 4*Gt212) + 2*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)) + 
           Gt112*(gt11L*Gt123 + gt12L*(2*Gt113 + Gt223) + 2*(Gt213*gt22L + gt23L*Gt313) + gt13L*Gt323) + 
           2*(Gt113*gt12L*Gt323 + Gt312*(gt12L*Gt133 + gt22L*Gt233 + gt23L*Gt333)) + 
           Gt323*(Gt112*gt13L + 4*Gt212*gt23L + 2*(Gt213*gt22L + gt23L*Gt313 + Gt312*gt33L)))*gtu31 + 
        J12L*(gt12L*PDstandardNth1Xt1 + gt22L*PDstandardNth1Xt2 + gt23L*PDstandardNth1Xt3) + 
        J22L*(gt12L*PDstandardNth2Xt1 + gt22L*PDstandardNth2Xt2 + gt23L*PDstandardNth2Xt3) + 
        J32L*(gt12L*PDstandardNth3Xt1 + gt22L*PDstandardNth3Xt2 + gt23L*PDstandardNth3Xt3) + 
        (Gt112*gt12L + Gt212*gt22L + gt23L*Gt312)*Xt1L + (Gt122*gt12L + Gt222*gt22L + gt23L*Gt322)*Xt2L + 
        (Gt123*gt12L + Gt223*gt22L + gt23L*Gt323)*Xt3L + 
        gtu21*(Gt222*(3*Gt112*gt12L + 6*Gt212*gt22L + 4*gt23L*Gt312) + 
           Gt122*(Gt112*gt11L + gt12L*(2*Gt111 + 4*Gt212) + 2*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)) + 
           Gt112*(gt11L*Gt122 + gt12L*Gt222 + 2*(Gt212*gt22L + gt23L*Gt312) + gt13L*Gt322) + 
           Gt322*(Gt112*gt13L + 4*Gt212*gt23L + 2*(Gt213*gt22L + gt23L*Gt313 + Gt312*gt33L)) + 
           2*(Gt312*(Gt123*gt12L + Gt223*gt22L + gt23L*Gt323) + gt12L*(Gt113*Gt322 + SQR(Gt112)))) + 
        gtu11*(Gt112*(gt12L*(2*Gt111 + 4*Gt212) + 2*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)) + 
           Gt312*(2*(Gt113*gt12L + Gt213*gt22L) + gt23L*(4*Gt212 + 2*Gt313)) + gt11L*SQR(Gt112) + 3*gt22L*SQR(Gt212) + 
           gt33L*SQR(Gt312)) + gtu22*(Gt122*(gt12L*(2*Gt112 + 4*Gt222) + 2*(Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322)) + 
           Gt322*(2*(Gt123*gt12L + Gt223*gt22L) + gt23L*(4*Gt222 + 2*Gt323)) + gt11L*SQR(Gt122) + 3*gt22L*SQR(Gt222) + 
           gt33L*SQR(Gt322)) + gtu33*(Gt123*(gt12L*(2*Gt113 + 4*Gt223) + 2*(Gt213*gt22L + gt23L*Gt313 + gt13L*Gt323)) + 
           Gt323*(2*(gt12L*Gt133 + gt22L*Gt233) + gt23L*(4*Gt223 + 2*Gt333)) + gt11L*SQR(Gt123) + 3*gt22L*SQR(Gt223) + 
           gt33L*SQR(Gt323)) + gtu32*(gt22L*(2*(Gt122*Gt213 + Gt233*Gt322) + Gt223*(6*Gt222 + 2*Gt323)) + 
           4*(gt12L*(Gt123*Gt222 + Gt122*Gt223) + gt23L*(Gt223*Gt322 + Gt222*Gt323)) + 
           2*(Gt123*(Gt112*gt12L + Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322) + gt12L*(Gt133*Gt322 + Gt123*Gt323) + 
              Gt122*(gt11L*Gt123 + Gt113*gt12L + gt23L*Gt313 + gt13L*Gt323) + Gt322*(gt23L*Gt333 + Gt323*gt33L) + 
              gt23L*SQR(Gt323))) + khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt22 + J12L*J21L*PDstandardNth12gt22 + 
                 J11L*J22L*PDstandardNth12gt22 + J12L*J31L*PDstandardNth13gt22 + J11L*J32L*PDstandardNth13gt22 + 
                 dJ112L*PDstandardNth1gt22 + J21L*J22L*PDstandardNth22gt22 + J22L*J31L*PDstandardNth23gt22 + 
                 J21L*J32L*PDstandardNth23gt22 + dJ212L*PDstandardNth2gt22 + J31L*J32L*PDstandardNth33gt22 + 
                 dJ312L*PDstandardNth3gt22) + gtu31*(J11L*J13L*PDstandardNth11gt22 + J13L*J21L*PDstandardNth12gt22 + 
                 J11L*J23L*PDstandardNth12gt22 + J13L*J31L*PDstandardNth13gt22 + J11L*J33L*PDstandardNth13gt22 + 
                 dJ113L*PDstandardNth1gt22 + J21L*J23L*PDstandardNth22gt22 + J23L*J31L*PDstandardNth23gt22 + 
                 J21L*J33L*PDstandardNth23gt22 + dJ213L*PDstandardNth2gt22 + J31L*J33L*PDstandardNth33gt22 + 
                 dJ313L*PDstandardNth3gt22) + gtu32*(J12L*J13L*PDstandardNth11gt22 + J13L*J22L*PDstandardNth12gt22 + 
                 J12L*J23L*PDstandardNth12gt22 + J13L*J32L*PDstandardNth13gt22 + J12L*J33L*PDstandardNth13gt22 + 
                 dJ123L*PDstandardNth1gt22 + J22L*J23L*PDstandardNth22gt22 + J23L*J32L*PDstandardNth23gt22 + 
                 J22L*J33L*PDstandardNth23gt22 + dJ223L*PDstandardNth2gt22 + J32L*J33L*PDstandardNth33gt22 + 
                 dJ323L*PDstandardNth3gt22)) - gtu11*
            (2*J11L*J21L*PDstandardNth12gt22 + 2*J11L*J31L*PDstandardNth13gt22 + dJ111L*PDstandardNth1gt22 + 
              2*J21L*J31L*PDstandardNth23gt22 + dJ211L*PDstandardNth2gt22 + dJ311L*PDstandardNth3gt22 + 
              PDstandardNth11gt22*SQR(J11L) + PDstandardNth22gt22*SQR(J21L) + PDstandardNth33gt22*SQR(J31L)) - 
           gtu22*(2*J12L*J22L*PDstandardNth12gt22 + 2*J12L*J32L*PDstandardNth13gt22 + dJ122L*PDstandardNth1gt22 + 
              2*J22L*J32L*PDstandardNth23gt22 + dJ222L*PDstandardNth2gt22 + dJ322L*PDstandardNth3gt22 + 
              PDstandardNth11gt22*SQR(J12L) + PDstandardNth22gt22*SQR(J22L) + PDstandardNth33gt22*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt22 + 2*J13L*J33L*PDstandardNth13gt22 + dJ133L*PDstandardNth1gt22 + 
              2*J23L*J33L*PDstandardNth23gt22 + dJ233L*PDstandardNth2gt22 + dJ333L*PDstandardNth3gt22 + 
              PDstandardNth11gt22*SQR(J13L) + PDstandardNth22gt22*SQR(J23L) + PDstandardNth33gt22*SQR(J33L)));
    
    Rt23  =  khalf*((gt13L*J12L + gt12L*J13L)*PDstandardNth1Xt1 + (gt23L*J12L + gt22L*J13L)*PDstandardNth1Xt2 + 
          (gt33L*J12L + gt23L*J13L)*PDstandardNth1Xt3 + (gt13L*J22L + gt12L*J23L)*PDstandardNth2Xt1 + 
          (gt23L*J22L + gt22L*J23L)*PDstandardNth2Xt2 + (gt33L*J22L + gt23L*J23L)*PDstandardNth2Xt3 - 
          2*(gtu21*(J11L*J12L*PDstandardNth11gt23 + J12L*J21L*PDstandardNth12gt23 + J11L*J22L*PDstandardNth12gt23 + 
                J12L*J31L*PDstandardNth13gt23 + J11L*J32L*PDstandardNth13gt23 + dJ112L*PDstandardNth1gt23 + 
                J21L*J22L*PDstandardNth22gt23 + J22L*J31L*PDstandardNth23gt23 + J21L*J32L*PDstandardNth23gt23 + 
                dJ212L*PDstandardNth2gt23 + J31L*J32L*PDstandardNth33gt23 + dJ312L*PDstandardNth3gt23) + 
             gtu31*(J11L*J13L*PDstandardNth11gt23 + J13L*J21L*PDstandardNth12gt23 + J11L*J23L*PDstandardNth12gt23 + 
                J13L*J31L*PDstandardNth13gt23 + J11L*J33L*PDstandardNth13gt23 + dJ113L*PDstandardNth1gt23 + 
                J21L*J23L*PDstandardNth22gt23 + J23L*J31L*PDstandardNth23gt23 + J21L*J33L*PDstandardNth23gt23 + 
                dJ213L*PDstandardNth2gt23 + J31L*J33L*PDstandardNth33gt23 + dJ313L*PDstandardNth3gt23) + 
             gtu32*(J12L*J13L*PDstandardNth11gt23 + J13L*J22L*PDstandardNth12gt23 + J12L*J23L*PDstandardNth12gt23 + 
                J13L*J32L*PDstandardNth13gt23 + J12L*J33L*PDstandardNth13gt23 + dJ123L*PDstandardNth1gt23 + 
                J22L*J23L*PDstandardNth22gt23 + J23L*J32L*PDstandardNth23gt23 + J22L*J33L*PDstandardNth23gt23 + 
                dJ223L*PDstandardNth2gt23 + J32L*J33L*PDstandardNth33gt23 + dJ323L*PDstandardNth3gt23)) + 
          (gt13L*J32L + gt12L*J33L)*PDstandardNth3Xt1 + (gt23L*J32L + gt22L*J33L)*PDstandardNth3Xt2 + 
          (gt33L*J32L + gt23L*J33L)*PDstandardNth3Xt3 + 
          (Gt113*gt12L + Gt112*gt13L + Gt213*gt22L + Gt212*gt23L + gt23L*Gt313 + Gt312*gt33L)*Xt1L + 
          (Gt123*gt12L + Gt122*gt13L + Gt223*gt22L + Gt222*gt23L + gt23L*Gt323 + Gt322*gt33L)*Xt2L + 
          (gt12L*Gt133 + Gt123*gt13L + gt22L*Gt233 + Gt223*gt23L + gt23L*Gt333 + Gt323*gt33L)*Xt3L + 
          2*((Gt112*gt11L*Gt123 + Gt111*Gt123*gt12L + Gt111*Gt122*gt13L + Gt123*gt12L*Gt212 + Gt112*gt13L*Gt222 + 
                2*Gt112*gt12L*Gt223 + Gt123*Gt211*gt22L + 2*Gt212*Gt223*gt22L + Gt122*Gt211*gt23L + Gt212*Gt222*gt23L + 
                Gt123*gt23L*Gt311 + Gt123*gt13L*Gt312 + 2*Gt223*gt23L*Gt312 + Gt113*gt13L*Gt322 + Gt213*gt23L*Gt322 + 
                Gt113*gt12L*Gt323 + Gt112*gt13L*Gt323 + Gt213*gt22L*Gt323 + Gt212*gt23L*Gt323 + gt23L*Gt313*Gt323 + 
                Gt122*Gt311*gt33L + Gt222*Gt312*gt33L + Gt313*Gt322*gt33L + Gt312*Gt323*gt33L)*gtu21 + 
             (Gt112*gt11L*Gt133 + Gt111*gt12L*Gt133 + Gt111*Gt123*gt13L + gt12L*Gt133*Gt212 + Gt112*gt13L*Gt223 + 
                Gt133*Gt211*gt22L + 2*Gt112*gt12L*Gt233 + 2*Gt212*gt22L*Gt233 + Gt123*Gt211*gt23L + Gt212*Gt223*gt23L + 
                Gt133*gt23L*Gt311 + Gt133*gt13L*Gt312 + 2*Gt233*gt23L*Gt312 + Gt113*gt13L*Gt323 + Gt213*gt23L*Gt323 + 
                Gt113*gt12L*Gt333 + Gt112*gt13L*Gt333 + Gt213*gt22L*Gt333 + Gt212*gt23L*Gt333 + gt23L*Gt313*Gt333 + 
                Gt123*Gt311*gt33L + Gt223*Gt312*gt33L + Gt313*Gt323*gt33L + Gt312*Gt333*gt33L)*gtu31 + 
             gtu21*(Gt113*gt11L*Gt122 + Gt122*gt13L*Gt212 + 2*Gt122*gt12L*Gt213 + Gt113*gt12L*Gt222 + Gt113*Gt212*gt22L + 
                2*Gt213*Gt222*gt22L + Gt212*Gt222*gt23L + Gt123*gt13L*Gt312 + Gt113*gt23L*Gt312 + Gt223*gt23L*Gt312 + 
                Gt123*gt12L*Gt313 + Gt122*gt13L*Gt313 + Gt223*gt22L*Gt313 + Gt222*gt23L*Gt313 + Gt113*gt13L*Gt322 + 
                2*Gt213*gt23L*Gt322 + gt23L*Gt313*Gt323 + Gt212*Gt322*gt33L + Gt313*Gt322*gt33L + Gt312*Gt323*gt33L + 
                Gt112*(Gt113*gt12L + Gt212*gt23L + Gt312*gt33L) + gt13L*SQR(Gt112)) + 
             gtu31*(2*Gt213*Gt223*gt22L + Gt112*Gt213*gt23L + Gt212*Gt223*gt23L + Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + 
                gt12L*Gt133*Gt313 + gt22L*Gt233*Gt313 + Gt223*gt23L*Gt313 + Gt123*(2*gt12L*Gt213 + gt13L*(Gt212 + Gt313)) + 
                2*Gt213*gt23L*Gt323 + Gt113*(gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt213*gt22L + gt23L*Gt313 + 
                   gt13L*Gt323) + gt23L*Gt313*Gt333 + Gt112*Gt313*gt33L + Gt212*Gt323*gt33L + Gt313*Gt323*gt33L + 
                Gt312*Gt333*gt33L + gt12L*SQR(Gt113)) + 
             gtu11*(Gt112*Gt113*gt11L + Gt111*Gt113*gt12L + Gt111*Gt112*gt13L + Gt113*gt12L*Gt212 + Gt112*gt13L*Gt212 + 
                2*Gt112*gt12L*Gt213 + Gt113*Gt211*gt22L + 2*Gt212*Gt213*gt22L + Gt112*Gt211*gt23L + Gt113*gt23L*Gt311 + 
                2*Gt113*gt13L*Gt312 + 3*Gt213*gt23L*Gt312 + Gt113*gt12L*Gt313 + Gt112*gt13L*Gt313 + Gt213*gt22L*Gt313 + 
                Gt212*gt23L*Gt313 + Gt112*Gt311*gt33L + Gt212*Gt312*gt33L + 2*Gt312*Gt313*gt33L + gt23L*SQR(Gt212) + 
                gt23L*SQR(Gt313)) + gtu22*(gt11L*Gt122*Gt123 + Gt112*Gt123*gt12L + Gt112*Gt122*gt13L + Gt123*gt12L*Gt222 + 
                Gt122*gt13L*Gt222 + 2*Gt122*gt12L*Gt223 + Gt123*Gt212*gt22L + 2*Gt222*Gt223*gt22L + Gt122*Gt212*gt23L + 
                Gt123*gt23L*Gt312 + 2*Gt123*gt13L*Gt322 + 3*Gt223*gt23L*Gt322 + Gt123*gt12L*Gt323 + Gt122*gt13L*Gt323 + 
                Gt223*gt22L*Gt323 + Gt222*gt23L*Gt323 + Gt122*Gt312*gt33L + Gt222*Gt322*gt33L + 2*Gt322*Gt323*gt33L + 
                gt23L*SQR(Gt222) + gt23L*SQR(Gt323)) + 
             gtu32*(gt11L*Gt122*Gt133 + Gt112*gt12L*Gt133 + Gt112*Gt123*gt13L + gt12L*Gt133*Gt222 + Gt122*gt13L*Gt223 + 
                Gt133*Gt212*gt22L + 2*Gt122*gt12L*Gt233 + 2*Gt222*gt22L*Gt233 + Gt123*Gt212*gt23L + Gt222*Gt223*gt23L + 
                Gt133*gt23L*Gt312 + Gt133*gt13L*Gt322 + 2*Gt233*gt23L*Gt322 + Gt123*gt13L*Gt323 + Gt223*gt23L*Gt323 + 
                Gt123*gt12L*Gt333 + Gt122*gt13L*Gt333 + Gt223*gt22L*Gt333 + Gt222*gt23L*Gt333 + gt23L*Gt323*Gt333 + 
                Gt123*Gt312*gt33L + Gt223*Gt322*gt33L + Gt322*Gt333*gt33L + gt33L*SQR(Gt323)) + 
             gtu32*(Gt113*Gt123*gt12L + Gt113*Gt122*gt13L + Gt123*gt13L*Gt222 + 3*Gt123*gt12L*Gt223 + Gt123*Gt213*gt22L + 
                Gt122*Gt213*gt23L + Gt222*Gt223*gt23L + Gt123*gt23L*Gt313 + Gt133*gt13L*Gt322 + Gt233*gt23L*Gt322 + 
                gt12L*Gt133*Gt323 + 2*Gt123*gt13L*Gt323 + gt22L*Gt233*Gt323 + 3*Gt223*gt23L*Gt323 + gt23L*Gt323*Gt333 + 
                Gt122*Gt313*gt33L + Gt222*Gt323*gt33L + Gt322*Gt333*gt33L + gt11L*SQR(Gt123) + 2*gt22L*SQR(Gt223) + 
                gt33L*SQR(Gt323)) + gtu33*(gt11L*Gt123*Gt133 + Gt113*gt12L*Gt133 + Gt113*Gt123*gt13L + gt12L*Gt133*Gt223 + 
                Gt123*gt13L*Gt223 + Gt133*Gt213*gt22L + 2*Gt123*gt12L*Gt233 + 2*Gt223*gt22L*Gt233 + Gt123*Gt213*gt23L + 
                Gt133*gt23L*Gt313 + 2*Gt133*gt13L*Gt323 + 3*Gt233*gt23L*Gt323 + gt12L*Gt133*Gt333 + Gt123*gt13L*Gt333 + 
                gt22L*Gt233*Gt333 + Gt223*gt23L*Gt333 + Gt123*Gt313*gt33L + Gt223*Gt323*gt33L + 2*Gt323*Gt333*gt33L + 
                gt23L*SQR(Gt223) + gt23L*SQR(Gt333))) - 
          gtu11*(2*J11L*J21L*PDstandardNth12gt23 + 2*J11L*J31L*PDstandardNth13gt23 + dJ111L*PDstandardNth1gt23 + 
             2*J21L*J31L*PDstandardNth23gt23 + dJ211L*PDstandardNth2gt23 + dJ311L*PDstandardNth3gt23 + 
             PDstandardNth11gt23*SQR(J11L) + PDstandardNth22gt23*SQR(J21L) + PDstandardNth33gt23*SQR(J31L)) - 
          gtu22*(2*J12L*J22L*PDstandardNth12gt23 + 2*J12L*J32L*PDstandardNth13gt23 + dJ122L*PDstandardNth1gt23 + 
             2*J22L*J32L*PDstandardNth23gt23 + dJ222L*PDstandardNth2gt23 + dJ322L*PDstandardNth3gt23 + 
             PDstandardNth11gt23*SQR(J12L) + PDstandardNth22gt23*SQR(J22L) + PDstandardNth33gt23*SQR(J32L)) - 
          gtu33*(2*J13L*J23L*PDstandardNth12gt23 + 2*J13L*J33L*PDstandardNth13gt23 + dJ133L*PDstandardNth1gt23 + 
             2*J23L*J33L*PDstandardNth23gt23 + dJ233L*PDstandardNth2gt23 + dJ333L*PDstandardNth3gt23 + 
             PDstandardNth11gt23*SQR(J13L) + PDstandardNth22gt23*SQR(J23L) + PDstandardNth33gt23*SQR(J33L)));
    
    Rt33  =  (4*((Gt123*gt13L + Gt223*gt23L)*Gt313 + (Gt113*gt13L + Gt213*gt23L)*Gt323) + 
           (2*Gt213*Gt322 + 6*Gt313*Gt323)*gt33L + 2*
            (gt13L*(Gt122*Gt213 + Gt112*Gt223) + Gt213*(Gt223*gt22L + Gt222*gt23L) + 
              Gt123*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L + Gt311*gt33L) + Gt223*(Gt212*gt23L + Gt312*gt33L) + 
              Gt113*(gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + Gt312*gt33L)))*gtu21 + 
        J13L*(gt13L*PDstandardNth1Xt1 + gt23L*PDstandardNth1Xt2 + gt33L*PDstandardNth1Xt3) + 
        J23L*(gt13L*PDstandardNth2Xt1 + gt23L*PDstandardNth2Xt2 + gt33L*PDstandardNth2Xt3) + 
        J33L*(gt13L*PDstandardNth3Xt1 + gt23L*PDstandardNth3Xt2 + gt33L*PDstandardNth3Xt3) + 
        (Gt113*gt13L + Gt213*gt23L + Gt313*gt33L)*Xt1L + (Gt123*gt13L + Gt223*gt23L + Gt323*gt33L)*Xt2L + 
        (Gt133*gt13L + Gt233*gt23L + Gt333*gt33L)*Xt3L + 
        gtu31*(Gt133*(Gt113*gt11L + 2*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L) + 4*gt13L*Gt313) + 
           Gt333*(3*Gt113*gt13L + 4*Gt213*gt23L + 6*Gt313*gt33L) + 
           Gt233*(Gt113*gt12L + 4*gt23L*Gt313 + 2*(Gt213*gt22L + Gt212*gt23L + Gt312*gt33L)) + 
           Gt113*(gt11L*Gt133 + gt12L*Gt233 + gt13L*Gt333 + 2*(Gt213*gt23L + Gt313*gt33L)) + 
           2*(Gt133*Gt311*gt33L + Gt213*(Gt223*gt23L + Gt323*gt33L) + gt13L*(Gt123*Gt213 + Gt112*Gt233 + SQR(Gt113)))) + 
        gtu32*(4*((Gt133*gt13L + Gt233*gt23L)*Gt323 + (Gt123*gt13L + Gt223*gt23L)*Gt333) + 
           Gt323*(2*Gt223 + 6*Gt333)*gt33L + 2*(Gt133*(Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + Gt312*gt33L) + 
              Gt123*(gt11L*Gt133 + gt13L*(Gt113 + Gt223) + gt12L*Gt233 + Gt213*gt23L + Gt313*gt33L) + 
              Gt233*(Gt122*gt13L + Gt223*gt22L + Gt222*gt23L + Gt322*gt33L) + gt23L*SQR(Gt223))) + 
        gtu11*(4*(Gt113*gt13L + Gt213*gt23L)*Gt313 + 
           2*(Gt113*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L + Gt311*gt33L) + 
              Gt213*(Gt112*gt13L + Gt212*gt23L + Gt312*gt33L)) + gt11L*SQR(Gt113) + gt22L*SQR(Gt213) + 3*gt33L*SQR(Gt313))\
         + gtu22*(4*(Gt123*gt13L + Gt223*gt23L)*Gt323 + 
           2*(Gt123*(Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + Gt312*gt33L) + 
              Gt223*(Gt122*gt13L + Gt222*gt23L + Gt322*gt33L)) + gt11L*SQR(Gt123) + gt22L*SQR(Gt223) + 3*gt33L*SQR(Gt323))\
         + gtu33*(4*(Gt133*gt13L + Gt233*gt23L)*Gt333 + 
           2*(Gt133*(Gt113*gt13L + gt12L*Gt233 + Gt213*gt23L + Gt313*gt33L) + 
              Gt233*(Gt123*gt13L + Gt223*gt23L + Gt323*gt33L)) + gt11L*SQR(Gt133) + gt22L*SQR(Gt233) + 3*gt33L*SQR(Gt333))\
         + khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt33 + J12L*J21L*PDstandardNth12gt33 + 
                 J11L*J22L*PDstandardNth12gt33 + J12L*J31L*PDstandardNth13gt33 + J11L*J32L*PDstandardNth13gt33 + 
                 dJ112L*PDstandardNth1gt33 + J21L*J22L*PDstandardNth22gt33 + J22L*J31L*PDstandardNth23gt33 + 
                 J21L*J32L*PDstandardNth23gt33 + dJ212L*PDstandardNth2gt33 + J31L*J32L*PDstandardNth33gt33 + 
                 dJ312L*PDstandardNth3gt33) + gtu31*(J11L*J13L*PDstandardNth11gt33 + J13L*J21L*PDstandardNth12gt33 + 
                 J11L*J23L*PDstandardNth12gt33 + J13L*J31L*PDstandardNth13gt33 + J11L*J33L*PDstandardNth13gt33 + 
                 dJ113L*PDstandardNth1gt33 + J21L*J23L*PDstandardNth22gt33 + J23L*J31L*PDstandardNth23gt33 + 
                 J21L*J33L*PDstandardNth23gt33 + dJ213L*PDstandardNth2gt33 + J31L*J33L*PDstandardNth33gt33 + 
                 dJ313L*PDstandardNth3gt33) + gtu32*(J12L*J13L*PDstandardNth11gt33 + J13L*J22L*PDstandardNth12gt33 + 
                 J12L*J23L*PDstandardNth12gt33 + J13L*J32L*PDstandardNth13gt33 + J12L*J33L*PDstandardNth13gt33 + 
                 dJ123L*PDstandardNth1gt33 + J22L*J23L*PDstandardNth22gt33 + J23L*J32L*PDstandardNth23gt33 + 
                 J22L*J33L*PDstandardNth23gt33 + dJ223L*PDstandardNth2gt33 + J32L*J33L*PDstandardNth33gt33 + 
                 dJ323L*PDstandardNth3gt33)) - gtu11*
            (2*J11L*J21L*PDstandardNth12gt33 + 2*J11L*J31L*PDstandardNth13gt33 + dJ111L*PDstandardNth1gt33 + 
              2*J21L*J31L*PDstandardNth23gt33 + dJ211L*PDstandardNth2gt33 + dJ311L*PDstandardNth3gt33 + 
              PDstandardNth11gt33*SQR(J11L) + PDstandardNth22gt33*SQR(J21L) + PDstandardNth33gt33*SQR(J31L)) - 
           gtu22*(2*J12L*J22L*PDstandardNth12gt33 + 2*J12L*J32L*PDstandardNth13gt33 + dJ122L*PDstandardNth1gt33 + 
              2*J22L*J32L*PDstandardNth23gt33 + dJ222L*PDstandardNth2gt33 + dJ322L*PDstandardNth3gt33 + 
              PDstandardNth11gt33*SQR(J12L) + PDstandardNth22gt33*SQR(J22L) + PDstandardNth33gt33*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt33 + 2*J13L*J33L*PDstandardNth13gt33 + dJ133L*PDstandardNth1gt33 + 
              2*J23L*J33L*PDstandardNth23gt33 + dJ233L*PDstandardNth2gt33 + dJ333L*PDstandardNth3gt33 + 
              PDstandardNth11gt33*SQR(J13L) + PDstandardNth22gt33*SQR(J23L) + PDstandardNth33gt33*SQR(J33L)));
    
    Rphi11  =  -2*((dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L)*PDstandardNth1phi + 
           2*(J11L*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + J21L*J31L*PDstandardNth23phi) + 
           (dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L)*PDstandardNth2phi + 
           (dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L)*PDstandardNth3phi + PDstandardNth11phi*SQR(J11L) + 
           PDstandardNth22phi*SQR(J21L) + PDstandardNth33phi*SQR(J31L) + 
           gt11L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*PDstandardNth1phi + 
              (2*(dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*PDstandardNth2phi + 
              (dJ322L*gtu22 + 2*(dJ313L*gtu31 + dJ323L*gtu32) + (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + (Gt213*gtu31 + Gt223*gtu32)*J32L + (Gt313*gtu31 + Gt323*gtu32)*J33L))
                *PDstandardNth3phi + 2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L)) + 
              PDstandardNth22phi*(gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L)))) - 
        4*gt11L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*PDstandardNth3phi + 
              PDstandardNth1phi*(((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*
                  PDstandardNth2phi + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*
                  PDstandardNth3phi)) + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + gtu11*SQR(J11L) + 
              gtu22*SQR(J12L) + gtu33*SQR(J13L))*SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + gtu11*SQR(J31L) + 
              gtu22*SQR(J32L) + gtu33*SQR(J33L))*SQR(PDstandardNth3phi)) + 
        4*SQR(J11L*PDstandardNth1phi + J21L*PDstandardNth2phi + J31L*PDstandardNth3phi);
    
    Rphi12  =  4*(J11L*PDstandardNth1phi + J21L*PDstandardNth2phi + J31L*PDstandardNth3phi)*
         (J12L*PDstandardNth1phi + J22L*PDstandardNth2phi + J32L*PDstandardNth3phi) - 
        2*(J12L*(J11L*PDstandardNth11phi + J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + 
           J11L*(J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + 
           (dJ112L - Gt112*J11L - Gt212*J12L - Gt312*J13L)*PDstandardNth1phi + 
           J22L*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi) + 
           (dJ212L - Gt112*J21L - Gt212*J22L - Gt312*J23L)*PDstandardNth2phi + 
           J32L*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi) + 
           (dJ312L - Gt112*J31L - Gt212*J32L - Gt312*J33L)*PDstandardNth3phi + 
           gt12L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*PDstandardNth1phi + 
              (2*(dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*PDstandardNth2phi + 
              (dJ322L*gtu22 + 2*(dJ313L*gtu31 + dJ323L*gtu32) + (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + (Gt213*gtu31 + Gt223*gtu32)*J32L + (Gt313*gtu31 + Gt323*gtu32)*J33L))
                *PDstandardNth3phi + 2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L)) + 
              PDstandardNth22phi*(gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L)))) - 
        4*gt12L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*PDstandardNth3phi + 
              PDstandardNth1phi*(((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*
                  PDstandardNth2phi + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*
                  PDstandardNth3phi)) + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + gtu11*SQR(J11L) + 
              gtu22*SQR(J12L) + gtu33*SQR(J13L))*SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + gtu11*SQR(J31L) + 
              gtu22*SQR(J32L) + gtu33*SQR(J33L))*SQR(PDstandardNth3phi));
    
    Rphi13  =  4*(J11L*PDstandardNth1phi + J21L*PDstandardNth2phi + J31L*PDstandardNth3phi)*
         (J13L*PDstandardNth1phi + J23L*PDstandardNth2phi + J33L*PDstandardNth3phi) - 
        2*(J13L*(J11L*PDstandardNth11phi + J21L*PDstandardNth12phi + J31L*PDstandardNth13phi) + 
           J11L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
           (dJ113L - Gt113*J11L - Gt213*J12L - Gt313*J13L)*PDstandardNth1phi + 
           J23L*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi) + 
           (dJ213L - Gt113*J21L - Gt213*J22L - Gt313*J23L)*PDstandardNth2phi + 
           J33L*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi) + 
           (dJ313L - Gt113*J31L - Gt213*J32L - Gt313*J33L)*PDstandardNth3phi + 
           gt13L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*PDstandardNth1phi + 
              (2*(dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*PDstandardNth2phi + 
              (dJ322L*gtu22 + 2*(dJ313L*gtu31 + dJ323L*gtu32) + (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + (Gt213*gtu31 + Gt223*gtu32)*J32L + (Gt313*gtu31 + Gt323*gtu32)*J33L))
                *PDstandardNth3phi + 2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L)) + 
              PDstandardNth22phi*(gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L)))) - 
        4*gt13L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*PDstandardNth3phi + 
              PDstandardNth1phi*(((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*
                  PDstandardNth2phi + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*
                  PDstandardNth3phi)) + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + gtu11*SQR(J11L) + 
              gtu22*SQR(J12L) + gtu33*SQR(J13L))*SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + gtu11*SQR(J31L) + 
              gtu22*SQR(J32L) + gtu33*SQR(J33L))*SQR(PDstandardNth3phi));
    
    Rphi22  =  -2*((dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L)*PDstandardNth1phi + 
           2*(J12L*(J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + J22L*J32L*PDstandardNth23phi) + 
           (dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L)*PDstandardNth2phi + 
           (dJ322L - Gt122*J31L - Gt222*J32L - Gt322*J33L)*PDstandardNth3phi + PDstandardNth11phi*SQR(J12L) + 
           PDstandardNth22phi*SQR(J22L) + PDstandardNth33phi*SQR(J32L) + 
           gt22L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*PDstandardNth1phi + 
              (2*(dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*PDstandardNth2phi + 
              (dJ322L*gtu22 + 2*(dJ313L*gtu31 + dJ323L*gtu32) + (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + (Gt213*gtu31 + Gt223*gtu32)*J32L + (Gt313*gtu31 + Gt323*gtu32)*J33L))
                *PDstandardNth3phi + 2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L)) + 
              PDstandardNth22phi*(gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L)))) - 
        4*gt22L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*PDstandardNth3phi + 
              PDstandardNth1phi*(((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*
                  PDstandardNth2phi + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*
                  PDstandardNth3phi)) + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + gtu11*SQR(J11L) + 
              gtu22*SQR(J12L) + gtu33*SQR(J13L))*SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + gtu11*SQR(J31L) + 
              gtu22*SQR(J32L) + gtu33*SQR(J33L))*SQR(PDstandardNth3phi)) + 
        4*SQR(J12L*PDstandardNth1phi + J22L*PDstandardNth2phi + J32L*PDstandardNth3phi);
    
    Rphi23  =  4*(J12L*PDstandardNth1phi + J22L*PDstandardNth2phi + J32L*PDstandardNth3phi)*
         (J13L*PDstandardNth1phi + J23L*PDstandardNth2phi + J33L*PDstandardNth3phi) - 
        2*(J13L*(J12L*PDstandardNth11phi + J22L*PDstandardNth12phi + J32L*PDstandardNth13phi) + 
           J12L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + 
           (dJ123L - Gt123*J11L - Gt223*J12L - Gt323*J13L)*PDstandardNth1phi + 
           J23L*(J22L*PDstandardNth22phi + J32L*PDstandardNth23phi) + 
           (dJ223L - Gt123*J21L - Gt223*J22L - Gt323*J23L)*PDstandardNth2phi + 
           J33L*(J22L*PDstandardNth23phi + J32L*PDstandardNth33phi) + 
           (dJ323L - Gt123*J31L - Gt223*J32L - Gt323*J33L)*PDstandardNth3phi + 
           gt23L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*PDstandardNth1phi + 
              (2*(dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*PDstandardNth2phi + 
              (dJ322L*gtu22 + 2*(dJ313L*gtu31 + dJ323L*gtu32) + (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + (Gt213*gtu31 + Gt223*gtu32)*J32L + (Gt313*gtu31 + Gt323*gtu32)*J33L))
                *PDstandardNth3phi + 2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L)) + 
              PDstandardNth22phi*(gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L)))) - 
        4*gt23L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*PDstandardNth3phi + 
              PDstandardNth1phi*(((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*
                  PDstandardNth2phi + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*
                  PDstandardNth3phi)) + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + gtu11*SQR(J11L) + 
              gtu22*SQR(J12L) + gtu33*SQR(J13L))*SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + gtu11*SQR(J31L) + 
              gtu22*SQR(J32L) + gtu33*SQR(J33L))*SQR(PDstandardNth3phi));
    
    Rphi33  =  -2*((dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L)*PDstandardNth1phi + 
           2*(J13L*(J23L*PDstandardNth12phi + J33L*PDstandardNth13phi) + J23L*J33L*PDstandardNth23phi) + 
           (dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L)*PDstandardNth2phi + 
           (dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L)*PDstandardNth3phi + PDstandardNth11phi*SQR(J13L) + 
           PDstandardNth22phi*SQR(J23L) + PDstandardNth33phi*SQR(J33L) + 
           gt33L*((2*(dJ112L*gtu21 + dJ113L*gtu31 + dJ123L*gtu32) + gtu11*(dJ111L - Gt111*J11L - Gt211*J12L - Gt311*J13L) + 
                 gtu22*(dJ122L - Gt122*J11L - Gt222*J12L - Gt322*J13L) + 
                 gtu33*(dJ133L - Gt133*J11L - Gt233*J12L - Gt333*J13L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J11L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J12L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J13L))*PDstandardNth1phi + 
              (2*(dJ212L*gtu21 + dJ213L*gtu31 + dJ223L*gtu32) + gtu11*(dJ211L - Gt111*J21L - Gt211*J22L - Gt311*J23L) + 
                 gtu22*(dJ222L - Gt122*J21L - Gt222*J22L - Gt322*J23L) + 
                 gtu33*(dJ233L - Gt133*J21L - Gt233*J22L - Gt333*J23L) - 
                 2*((Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32)*J21L + (Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32)*J22L + 
                    (Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32)*J23L))*PDstandardNth2phi + 
              (dJ322L*gtu22 + 2*(dJ313L*gtu31 + dJ323L*gtu32) + (-2*Gt112*gtu21 - Gt122*gtu22)*J31L + 
                 (-2*Gt212*gtu21 - Gt222*gtu22)*J32L + (-2*Gt312*gtu21 - Gt322*gtu22)*J33L + 
                 gtu11*(dJ311L - Gt111*J31L - Gt211*J32L - Gt311*J33L) + 
                 gtu33*(dJ333L - Gt133*J31L - Gt233*J32L - Gt333*J33L) - 
                 2*((Gt113*gtu31 + Gt123*gtu32)*J31L + (Gt213*gtu31 + Gt223*gtu32)*J32L + (Gt313*gtu31 + Gt323*gtu32)*J33L))
                *PDstandardNth3phi + 2*((gtu21*J11L*J22L + gtu31*(J13L*J21L + J11L*J23L))*PDstandardNth12phi + 
                 (gtu21*J11L*J32L + gtu31*(J13L*J31L + J11L*J33L))*PDstandardNth13phi + 
                 J11L*((gtu21*J12L + gtu31*J13L)*PDstandardNth11phi + 
                    gtu11*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 J12L*(gtu32*J13L*PDstandardNth11phi + gtu21*(J21L*PDstandardNth12phi + J31L*PDstandardNth13phi)) + 
                 (gtu11*J21L*J31L + (gtu22*J22L + gtu32*J23L)*J32L + (gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth23phi + 
                 J22L*((gtu22*J12L + gtu32*J13L)*PDstandardNth12phi + (gtu21*J21L + gtu32*J23L)*PDstandardNth22phi + 
                    gtu21*J31L*PDstandardNth23phi) + 
                 J23L*((gtu32*J12L + gtu33*J13L)*PDstandardNth12phi + 
                    gtu31*(J21L*PDstandardNth22phi + J31L*PDstandardNth23phi)) + 
                 J32L*((gtu22*J12L + gtu32*J13L)*PDstandardNth13phi + gtu32*J33L*PDstandardNth33phi + 
                    gtu21*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + 
                 J33L*((gtu32*J12L + gtu33*J13L)*PDstandardNth13phi + 
                    gtu31*(J21L*PDstandardNth23phi + J31L*PDstandardNth33phi)) + dJ312L*gtu21*PDstandardNth3phi) + 
              PDstandardNth11phi*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + gtu33*SQR(J13L)) + 
              PDstandardNth22phi*(gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
              PDstandardNth33phi*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + gtu33*SQR(J33L)))) - 
        4*gt33L*(2*(((gtu11*J21L + gtu21*J22L + gtu31*J23L)*J31L + (gtu21*J21L + gtu22*J22L + gtu32*J23L)*J32L + 
                 (gtu31*J21L + gtu32*J22L + gtu33*J23L)*J33L)*PDstandardNth2phi*PDstandardNth3phi + 
              PDstandardNth1phi*(((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J21L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J22L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J23L)*
                  PDstandardNth2phi + ((gtu11*J11L + gtu21*J12L + gtu31*J13L)*J31L + 
                    (gtu21*J11L + gtu22*J12L + gtu32*J13L)*J32L + (gtu31*J11L + gtu32*J12L + gtu33*J13L)*J33L)*
                  PDstandardNth3phi)) + (2*(gtu32*J12L*J13L + J11L*(gtu21*J12L + gtu31*J13L)) + gtu11*SQR(J11L) + 
              gtu22*SQR(J12L) + gtu33*SQR(J13L))*SQR(PDstandardNth1phi) + 
           (2*(gtu32*J22L*J23L + J21L*(gtu21*J22L + gtu31*J23L)) + gtu11*SQR(J21L) + gtu22*SQR(J22L) + gtu33*SQR(J23L))*
            SQR(PDstandardNth2phi) + (2*(gtu32*J32L*J33L + J31L*(gtu21*J32L + gtu31*J33L)) + gtu11*SQR(J31L) + 
              gtu22*SQR(J32L) + gtu33*SQR(J33L))*SQR(PDstandardNth3phi)) + 
        4*SQR(J13L*PDstandardNth1phi + J23L*PDstandardNth2phi + J33L*PDstandardNth3phi);
    
    e4phi  =  exp(4*phiL);
    
    em4phi  =  INV(e4phi);
    
    gu11  =  em4phi*gtu11;
    
    gu21  =  em4phi*gtu21;
    
    gu31  =  em4phi*gtu31;
    
    gu22  =  em4phi*gtu22;
    
    gu32  =  em4phi*gtu32;
    
    gu33  =  em4phi*gtu33;
    
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
    
    M1L  =  gtu11*(-2*At11L*Gt111 - 2*At12L*Gt211 - 2*At13L*Gt311 + J11L*PDstandardNth1At11 + 
           6*At11L*J11L*PDstandardNth1phi + J21L*PDstandardNth2At11 + 6*At11L*J21L*PDstandardNth2phi + 
           J31L*PDstandardNth3At11 + 6*At11L*J31L*PDstandardNth3phi) - 
        gtu21*(At12L*Gt111 + At11L*Gt112 + At22L*Gt211 + At12L*Gt212 + At23L*Gt311 + At13L*Gt312 - 
           J11L*PDstandardNth1At12 - 6*At12L*J11L*PDstandardNth1phi - J21L*PDstandardNth2At12 - 
           6*At12L*J21L*PDstandardNth2phi - J31L*PDstandardNth3At12 - 6*At12L*J31L*PDstandardNth3phi) - 
        gtu31*(At13L*Gt111 + At11L*Gt113 + At23L*Gt211 + At12L*Gt213 + At33L*Gt311 + At13L*Gt313 - 
           J11L*PDstandardNth1At13 - 6*At13L*J11L*PDstandardNth1phi - J21L*PDstandardNth2At13 - 
           6*At13L*J21L*PDstandardNth2phi - J31L*PDstandardNth3At13 - 6*At13L*J31L*PDstandardNth3phi) + 
        gtu21*(-2*At11L*Gt112 - 2*At12L*Gt212 - 2*At13L*Gt312 + J12L*PDstandardNth1At11 + 6*At11L*J12L*PDstandardNth1phi + 
           J22L*PDstandardNth2At11 + 6*At11L*J22L*PDstandardNth2phi + J32L*PDstandardNth3At11 + 
           6*At11L*J32L*PDstandardNth3phi) - gtu22*(At12L*Gt112 + At11L*Gt122 + At22L*Gt212 + At12L*Gt222 + At23L*Gt312 + 
           At13L*Gt322 - J12L*PDstandardNth1At12 - 6*At12L*J12L*PDstandardNth1phi - J22L*PDstandardNth2At12 - 
           6*At12L*J22L*PDstandardNth2phi - J32L*PDstandardNth3At12 - 6*At12L*J32L*PDstandardNth3phi) - 
        gtu32*(At13L*Gt112 + At11L*Gt123 + At23L*Gt212 + At12L*Gt223 + At33L*Gt312 + At13L*Gt323 - 
           J12L*PDstandardNth1At13 - 6*At13L*J12L*PDstandardNth1phi - J22L*PDstandardNth2At13 - 
           6*At13L*J22L*PDstandardNth2phi - J32L*PDstandardNth3At13 - 6*At13L*J32L*PDstandardNth3phi) + 
        gtu31*(-2*At11L*Gt113 - 2*At12L*Gt213 - 2*At13L*Gt313 + J13L*PDstandardNth1At11 + 6*At11L*J13L*PDstandardNth1phi + 
           J23L*PDstandardNth2At11 + 6*At11L*J23L*PDstandardNth2phi + J33L*PDstandardNth3At11 + 
           6*At11L*J33L*PDstandardNth3phi) - gtu32*(At12L*Gt113 + At11L*Gt123 + At22L*Gt213 + At12L*Gt223 + At23L*Gt313 + 
           At13L*Gt323 - J13L*PDstandardNth1At12 - 6*At12L*J13L*PDstandardNth1phi - J23L*PDstandardNth2At12 - 
           6*At12L*J23L*PDstandardNth2phi - J33L*PDstandardNth3At12 - 6*At12L*J33L*PDstandardNth3phi) - 
        gtu33*(At13L*Gt113 + At11L*Gt133 + At23L*Gt213 + At12L*Gt233 + At33L*Gt313 + At13L*Gt333 - 
           J13L*PDstandardNth1At13 - 6*At13L*J13L*PDstandardNth1phi - J23L*PDstandardNth2At13 - 
           6*At13L*J23L*PDstandardNth2phi - J33L*PDstandardNth3At13 - 6*At13L*J33L*PDstandardNth3phi) - 
        ktwothird*(J11L*PDstandardNth1trK + J21L*PDstandardNth2trK + J31L*PDstandardNth3trK);
    
    M2L  =  -(gtu11*(At12L*Gt111 + At11L*Gt112 + At22L*Gt211 + At12L*Gt212 + At23L*Gt311 + At13L*Gt312 - 
             J11L*PDstandardNth1At12 - 6*At12L*J11L*PDstandardNth1phi - J21L*PDstandardNth2At12 - 
             6*At12L*J21L*PDstandardNth2phi - J31L*PDstandardNth3At12 - 6*At12L*J31L*PDstandardNth3phi)) + 
        gtu21*(-2*At12L*Gt112 - 2*At22L*Gt212 - 2*At23L*Gt312 + J11L*PDstandardNth1At22 + 6*At22L*J11L*PDstandardNth1phi + 
           J21L*PDstandardNth2At22 + 6*At22L*J21L*PDstandardNth2phi + J31L*PDstandardNth3At22 + 
           6*At22L*J31L*PDstandardNth3phi) - gtu31*(At13L*Gt112 + At12L*Gt113 + At23L*Gt212 + At22L*Gt213 + At33L*Gt312 + 
           At23L*Gt313 - J11L*PDstandardNth1At23 - 6*At23L*J11L*PDstandardNth1phi - J21L*PDstandardNth2At23 - 
           6*At23L*J21L*PDstandardNth2phi - J31L*PDstandardNth3At23 - 6*At23L*J31L*PDstandardNth3phi) - 
        gtu21*(At12L*Gt112 + At11L*Gt122 + At22L*Gt212 + At12L*Gt222 + At23L*Gt312 + At13L*Gt322 - 
           J12L*PDstandardNth1At12 - 6*At12L*J12L*PDstandardNth1phi - J22L*PDstandardNth2At12 - 
           6*At12L*J22L*PDstandardNth2phi - J32L*PDstandardNth3At12 - 6*At12L*J32L*PDstandardNth3phi) + 
        gtu22*(-2*At12L*Gt122 - 2*At22L*Gt222 - 2*At23L*Gt322 + J12L*PDstandardNth1At22 + 6*At22L*J12L*PDstandardNth1phi + 
           J22L*PDstandardNth2At22 + 6*At22L*J22L*PDstandardNth2phi + J32L*PDstandardNth3At22 + 
           6*At22L*J32L*PDstandardNth3phi) - gtu32*(At13L*Gt122 + At12L*Gt123 + At23L*Gt222 + At22L*Gt223 + At33L*Gt322 + 
           At23L*Gt323 - J12L*PDstandardNth1At23 - 6*At23L*J12L*PDstandardNth1phi - J22L*PDstandardNth2At23 - 
           6*At23L*J22L*PDstandardNth2phi - J32L*PDstandardNth3At23 - 6*At23L*J32L*PDstandardNth3phi) - 
        gtu31*(At12L*Gt113 + At11L*Gt123 + At22L*Gt213 + At12L*Gt223 + At23L*Gt313 + At13L*Gt323 - 
           J13L*PDstandardNth1At12 - 6*At12L*J13L*PDstandardNth1phi - J23L*PDstandardNth2At12 - 
           6*At12L*J23L*PDstandardNth2phi - J33L*PDstandardNth3At12 - 6*At12L*J33L*PDstandardNth3phi) + 
        gtu32*(-2*At12L*Gt123 - 2*At22L*Gt223 - 2*At23L*Gt323 + J13L*PDstandardNth1At22 + 6*At22L*J13L*PDstandardNth1phi + 
           J23L*PDstandardNth2At22 + 6*At22L*J23L*PDstandardNth2phi + J33L*PDstandardNth3At22 + 
           6*At22L*J33L*PDstandardNth3phi) - gtu33*(At13L*Gt123 + At12L*Gt133 + At23L*Gt223 + At22L*Gt233 + At33L*Gt323 + 
           At23L*Gt333 - J13L*PDstandardNth1At23 - 6*At23L*J13L*PDstandardNth1phi - J23L*PDstandardNth2At23 - 
           6*At23L*J23L*PDstandardNth2phi - J33L*PDstandardNth3At23 - 6*At23L*J33L*PDstandardNth3phi) - 
        ktwothird*(J12L*PDstandardNth1trK + J22L*PDstandardNth2trK + J32L*PDstandardNth3trK);
    
    M3L  =  -(gtu11*(At13L*Gt111 + At11L*Gt113 + At23L*Gt211 + At12L*Gt213 + At33L*Gt311 + At13L*Gt313 - 
             J11L*PDstandardNth1At13 - 6*At13L*J11L*PDstandardNth1phi - J21L*PDstandardNth2At13 - 
             6*At13L*J21L*PDstandardNth2phi - J31L*PDstandardNth3At13 - 6*At13L*J31L*PDstandardNth3phi)) - 
        gtu21*(At13L*Gt112 + At12L*Gt113 + At23L*Gt212 + At22L*Gt213 + At33L*Gt312 + At23L*Gt313 - 
           J11L*PDstandardNth1At23 - 6*At23L*J11L*PDstandardNth1phi - J21L*PDstandardNth2At23 - 
           6*At23L*J21L*PDstandardNth2phi - J31L*PDstandardNth3At23 - 6*At23L*J31L*PDstandardNth3phi) + 
        gtu31*(-2*At13L*Gt113 - 2*At23L*Gt213 - 2*At33L*Gt313 + J11L*PDstandardNth1At33 + 6*At33L*J11L*PDstandardNth1phi + 
           J21L*PDstandardNth2At33 + 6*At33L*J21L*PDstandardNth2phi + J31L*PDstandardNth3At33 + 
           6*At33L*J31L*PDstandardNth3phi) - gtu21*(At13L*Gt112 + At11L*Gt123 + At23L*Gt212 + At12L*Gt223 + At33L*Gt312 + 
           At13L*Gt323 - J12L*PDstandardNth1At13 - 6*At13L*J12L*PDstandardNth1phi - J22L*PDstandardNth2At13 - 
           6*At13L*J22L*PDstandardNth2phi - J32L*PDstandardNth3At13 - 6*At13L*J32L*PDstandardNth3phi) - 
        gtu22*(At13L*Gt122 + At12L*Gt123 + At23L*Gt222 + At22L*Gt223 + At33L*Gt322 + At23L*Gt323 - 
           J12L*PDstandardNth1At23 - 6*At23L*J12L*PDstandardNth1phi - J22L*PDstandardNth2At23 - 
           6*At23L*J22L*PDstandardNth2phi - J32L*PDstandardNth3At23 - 6*At23L*J32L*PDstandardNth3phi) + 
        gtu32*(-2*At13L*Gt123 - 2*At23L*Gt223 - 2*At33L*Gt323 + J12L*PDstandardNth1At33 + 6*At33L*J12L*PDstandardNth1phi + 
           J22L*PDstandardNth2At33 + 6*At33L*J22L*PDstandardNth2phi + J32L*PDstandardNth3At33 + 
           6*At33L*J32L*PDstandardNth3phi) - gtu31*(At13L*Gt113 + At11L*Gt133 + At23L*Gt213 + At12L*Gt233 + At33L*Gt313 + 
           At13L*Gt333 - J13L*PDstandardNth1At13 - 6*At13L*J13L*PDstandardNth1phi - J23L*PDstandardNth2At13 - 
           6*At13L*J23L*PDstandardNth2phi - J33L*PDstandardNth3At13 - 6*At13L*J33L*PDstandardNth3phi) - 
        gtu32*(At13L*Gt123 + At12L*Gt133 + At23L*Gt223 + At22L*Gt233 + At33L*Gt323 + At23L*Gt333 - 
           J13L*PDstandardNth1At23 - 6*At23L*J13L*PDstandardNth1phi - J23L*PDstandardNth2At23 - 
           6*At23L*J23L*PDstandardNth2phi - J33L*PDstandardNth3At23 - 6*At23L*J33L*PDstandardNth3phi) + 
        gtu33*(-2*At13L*Gt133 - 2*At23L*Gt233 - 2*At33L*Gt333 + J13L*PDstandardNth1At33 + 6*At33L*J13L*PDstandardNth1phi + 
           J23L*PDstandardNth2At33 + 6*At33L*J23L*PDstandardNth2phi + J33L*PDstandardNth3At33 + 
           6*At33L*J33L*PDstandardNth3phi) - ktwothird*
         (J13L*PDstandardNth1trK + J23L*PDstandardNth2trK + J33L*PDstandardNth3trK);
    
    cSL  =  Log(detgt);
    
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
  LC_ENDLOOP3 (ML_BSSN_MP_constraints);
}

void ML_BSSN_MP_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_constraints_Body);
}
