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

void ML_ADM_RHS_Body(cGH const * const cctkGH, CCTK_INT const dir, CCTK_INT const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], CCTK_INT const min[3], CCTK_INT const max[3], CCTK_INT const n_subblock_gfs, CCTK_REAL * const subblock_gfs[])
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
  #pragma omp parallel
  LC_LOOP3 (ML_ADM_RHS,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lssh[CCTK_LSSH_IDX(0,0)],cctk_lssh[CCTK_LSSH_IDX(0,1)],cctk_lssh[CCTK_LSSH_IDX(0,2)])
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
    CCTK_REAL gu11 = INITVALUE, gu12 = INITVALUE, gu13 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu23 = INITVALUE;
    CCTK_REAL gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
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
    CCTK_REAL PDstandardNth1alpha = INITVALUE;
    CCTK_REAL PDstandardNth2alpha = INITVALUE;
    CCTK_REAL PDstandardNth3alpha = INITVALUE;
    CCTK_REAL PDstandardNth11alpha = INITVALUE;
    CCTK_REAL PDstandardNth22alpha = INITVALUE;
    CCTK_REAL PDstandardNth33alpha = INITVALUE;
    CCTK_REAL PDstandardNth12alpha = INITVALUE;
    CCTK_REAL PDstandardNth13alpha = INITVALUE;
    CCTK_REAL PDstandardNth23alpha = INITVALUE;
    CCTK_REAL PDstandardNth1beta1 = INITVALUE;
    CCTK_REAL PDstandardNth2beta1 = INITVALUE;
    CCTK_REAL PDstandardNth3beta1 = INITVALUE;
    CCTK_REAL PDstandardNth1beta2 = INITVALUE;
    CCTK_REAL PDstandardNth2beta2 = INITVALUE;
    CCTK_REAL PDstandardNth3beta2 = INITVALUE;
    CCTK_REAL PDstandardNth1beta3 = INITVALUE;
    CCTK_REAL PDstandardNth2beta3 = INITVALUE;
    CCTK_REAL PDstandardNth3beta3 = INITVALUE;
    CCTK_REAL PDstandardNth1g11 = INITVALUE;
    CCTK_REAL PDstandardNth2g11 = INITVALUE;
    CCTK_REAL PDstandardNth3g11 = INITVALUE;
    CCTK_REAL PDstandardNth22g11 = INITVALUE;
    CCTK_REAL PDstandardNth33g11 = INITVALUE;
    CCTK_REAL PDstandardNth23g11 = INITVALUE;
    CCTK_REAL PDstandardNth1g12 = INITVALUE;
    CCTK_REAL PDstandardNth2g12 = INITVALUE;
    CCTK_REAL PDstandardNth3g12 = INITVALUE;
    CCTK_REAL PDstandardNth33g12 = INITVALUE;
    CCTK_REAL PDstandardNth12g12 = INITVALUE;
    CCTK_REAL PDstandardNth13g12 = INITVALUE;
    CCTK_REAL PDstandardNth23g12 = INITVALUE;
    CCTK_REAL PDstandardNth1g13 = INITVALUE;
    CCTK_REAL PDstandardNth2g13 = INITVALUE;
    CCTK_REAL PDstandardNth3g13 = INITVALUE;
    CCTK_REAL PDstandardNth22g13 = INITVALUE;
    CCTK_REAL PDstandardNth12g13 = INITVALUE;
    CCTK_REAL PDstandardNth13g13 = INITVALUE;
    CCTK_REAL PDstandardNth23g13 = INITVALUE;
    CCTK_REAL PDstandardNth1g22 = INITVALUE;
    CCTK_REAL PDstandardNth2g22 = INITVALUE;
    CCTK_REAL PDstandardNth3g22 = INITVALUE;
    CCTK_REAL PDstandardNth11g22 = INITVALUE;
    CCTK_REAL PDstandardNth33g22 = INITVALUE;
    CCTK_REAL PDstandardNth13g22 = INITVALUE;
    CCTK_REAL PDstandardNth1g23 = INITVALUE;
    CCTK_REAL PDstandardNth2g23 = INITVALUE;
    CCTK_REAL PDstandardNth3g23 = INITVALUE;
    CCTK_REAL PDstandardNth11g23 = INITVALUE;
    CCTK_REAL PDstandardNth12g23 = INITVALUE;
    CCTK_REAL PDstandardNth13g23 = INITVALUE;
    CCTK_REAL PDstandardNth23g23 = INITVALUE;
    CCTK_REAL PDstandardNth1g33 = INITVALUE;
    CCTK_REAL PDstandardNth2g33 = INITVALUE;
    CCTK_REAL PDstandardNth3g33 = INITVALUE;
    CCTK_REAL PDstandardNth11g33 = INITVALUE;
    CCTK_REAL PDstandardNth22g33 = INITVALUE;
    CCTK_REAL PDstandardNth12g33 = INITVALUE;
    CCTK_REAL PDstandardNth1K11 = INITVALUE;
    CCTK_REAL PDstandardNth2K11 = INITVALUE;
    CCTK_REAL PDstandardNth3K11 = INITVALUE;
    CCTK_REAL PDstandardNth1K12 = INITVALUE;
    CCTK_REAL PDstandardNth2K12 = INITVALUE;
    CCTK_REAL PDstandardNth3K12 = INITVALUE;
    CCTK_REAL PDstandardNth1K13 = INITVALUE;
    CCTK_REAL PDstandardNth2K13 = INITVALUE;
    CCTK_REAL PDstandardNth3K13 = INITVALUE;
    CCTK_REAL PDstandardNth1K22 = INITVALUE;
    CCTK_REAL PDstandardNth2K22 = INITVALUE;
    CCTK_REAL PDstandardNth3K22 = INITVALUE;
    CCTK_REAL PDstandardNth1K23 = INITVALUE;
    CCTK_REAL PDstandardNth2K23 = INITVALUE;
    CCTK_REAL PDstandardNth3K23 = INITVALUE;
    CCTK_REAL PDstandardNth1K33 = INITVALUE;
    CCTK_REAL PDstandardNth2K33 = INITVALUE;
    CCTK_REAL PDstandardNth3K33 = INITVALUE;
    
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
    PDstandardNth1alpha = PDstandardNth1(alpha, i, j, k);
    PDstandardNth2alpha = PDstandardNth2(alpha, i, j, k);
    PDstandardNth3alpha = PDstandardNth3(alpha, i, j, k);
    PDstandardNth11alpha = PDstandardNth11(alpha, i, j, k);
    PDstandardNth22alpha = PDstandardNth22(alpha, i, j, k);
    PDstandardNth33alpha = PDstandardNth33(alpha, i, j, k);
    PDstandardNth12alpha = PDstandardNth12(alpha, i, j, k);
    PDstandardNth13alpha = PDstandardNth13(alpha, i, j, k);
    PDstandardNth23alpha = PDstandardNth23(alpha, i, j, k);
    PDstandardNth1beta1 = PDstandardNth1(beta1, i, j, k);
    PDstandardNth2beta1 = PDstandardNth2(beta1, i, j, k);
    PDstandardNth3beta1 = PDstandardNth3(beta1, i, j, k);
    PDstandardNth1beta2 = PDstandardNth1(beta2, i, j, k);
    PDstandardNth2beta2 = PDstandardNth2(beta2, i, j, k);
    PDstandardNth3beta2 = PDstandardNth3(beta2, i, j, k);
    PDstandardNth1beta3 = PDstandardNth1(beta3, i, j, k);
    PDstandardNth2beta3 = PDstandardNth2(beta3, i, j, k);
    PDstandardNth3beta3 = PDstandardNth3(beta3, i, j, k);
    PDstandardNth1g11 = PDstandardNth1(g11, i, j, k);
    PDstandardNth2g11 = PDstandardNth2(g11, i, j, k);
    PDstandardNth3g11 = PDstandardNth3(g11, i, j, k);
    PDstandardNth22g11 = PDstandardNth22(g11, i, j, k);
    PDstandardNth33g11 = PDstandardNth33(g11, i, j, k);
    PDstandardNth23g11 = PDstandardNth23(g11, i, j, k);
    PDstandardNth1g12 = PDstandardNth1(g12, i, j, k);
    PDstandardNth2g12 = PDstandardNth2(g12, i, j, k);
    PDstandardNth3g12 = PDstandardNth3(g12, i, j, k);
    PDstandardNth33g12 = PDstandardNth33(g12, i, j, k);
    PDstandardNth12g12 = PDstandardNth12(g12, i, j, k);
    PDstandardNth13g12 = PDstandardNth13(g12, i, j, k);
    PDstandardNth23g12 = PDstandardNth23(g12, i, j, k);
    PDstandardNth1g13 = PDstandardNth1(g13, i, j, k);
    PDstandardNth2g13 = PDstandardNth2(g13, i, j, k);
    PDstandardNth3g13 = PDstandardNth3(g13, i, j, k);
    PDstandardNth22g13 = PDstandardNth22(g13, i, j, k);
    PDstandardNth12g13 = PDstandardNth12(g13, i, j, k);
    PDstandardNth13g13 = PDstandardNth13(g13, i, j, k);
    PDstandardNth23g13 = PDstandardNth23(g13, i, j, k);
    PDstandardNth1g22 = PDstandardNth1(g22, i, j, k);
    PDstandardNth2g22 = PDstandardNth2(g22, i, j, k);
    PDstandardNth3g22 = PDstandardNth3(g22, i, j, k);
    PDstandardNth11g22 = PDstandardNth11(g22, i, j, k);
    PDstandardNth33g22 = PDstandardNth33(g22, i, j, k);
    PDstandardNth13g22 = PDstandardNth13(g22, i, j, k);
    PDstandardNth1g23 = PDstandardNth1(g23, i, j, k);
    PDstandardNth2g23 = PDstandardNth2(g23, i, j, k);
    PDstandardNth3g23 = PDstandardNth3(g23, i, j, k);
    PDstandardNth11g23 = PDstandardNth11(g23, i, j, k);
    PDstandardNth12g23 = PDstandardNth12(g23, i, j, k);
    PDstandardNth13g23 = PDstandardNth13(g23, i, j, k);
    PDstandardNth23g23 = PDstandardNth23(g23, i, j, k);
    PDstandardNth1g33 = PDstandardNth1(g33, i, j, k);
    PDstandardNth2g33 = PDstandardNth2(g33, i, j, k);
    PDstandardNth3g33 = PDstandardNth3(g33, i, j, k);
    PDstandardNth11g33 = PDstandardNth11(g33, i, j, k);
    PDstandardNth22g33 = PDstandardNth22(g33, i, j, k);
    PDstandardNth12g33 = PDstandardNth12(g33, i, j, k);
    PDstandardNth1K11 = PDstandardNth1(K11, i, j, k);
    PDstandardNth2K11 = PDstandardNth2(K11, i, j, k);
    PDstandardNth3K11 = PDstandardNth3(K11, i, j, k);
    PDstandardNth1K12 = PDstandardNth1(K12, i, j, k);
    PDstandardNth2K12 = PDstandardNth2(K12, i, j, k);
    PDstandardNth3K12 = PDstandardNth3(K12, i, j, k);
    PDstandardNth1K13 = PDstandardNth1(K13, i, j, k);
    PDstandardNth2K13 = PDstandardNth2(K13, i, j, k);
    PDstandardNth3K13 = PDstandardNth3(K13, i, j, k);
    PDstandardNth1K22 = PDstandardNth1(K22, i, j, k);
    PDstandardNth2K22 = PDstandardNth2(K22, i, j, k);
    PDstandardNth3K22 = PDstandardNth3(K22, i, j, k);
    PDstandardNth1K23 = PDstandardNth1(K23, i, j, k);
    PDstandardNth2K23 = PDstandardNth2(K23, i, j, k);
    PDstandardNth3K23 = PDstandardNth3(K23, i, j, k);
    PDstandardNth1K33 = PDstandardNth1(K33, i, j, k);
    PDstandardNth2K33 = PDstandardNth2(K33, i, j, k);
    PDstandardNth3K33 = PDstandardNth3(K33, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL const T1000001  =  SQR(g23L);
    
    CCTK_REAL const T1000003  =  g13L*g23L;
    
    CCTK_REAL const T1000004  =  -(g12L*g33L);
    
    CCTK_REAL const T1000007  =  SQR(g13L);
    
    CCTK_REAL const T1000008  =  -(g13L*g22L);
    
    CCTK_REAL const T1000009  =  g12L*g23L;
    
    CCTK_REAL const T1000012  =  g12L*g13L;
    
    CCTK_REAL const T1000013  =  -(g11L*g23L);
    
    CCTK_REAL const T1000016  =  SQR(g12L);
    
    CCTK_REAL const T1000017  =  -PDstandardNth3g12;
    
    CCTK_REAL const T1000019  =  -PDstandardNth2g13;
    
    CCTK_REAL const T1000021  =  -PDstandardNth1g23;
    
    CCTK_REAL const T1000023  =  -PDstandardNth11g23;
    
    CCTK_REAL const T1000024  =  -PDstandardNth23g11;
    
    CCTK_REAL const T1000029  =  -2*PDstandardNth12g12;
    
    CCTK_REAL const T1000031  =  -PDstandardNth13g22;
    
    CCTK_REAL const T1000032  =  -PDstandardNth22g13;
    
    CCTK_REAL const T1000034  =  -2*PDstandardNth13g13;
    
    CCTK_REAL const T1000036  =  -PDstandardNth12g33;
    
    CCTK_REAL const T1000037  =  -PDstandardNth33g12;
    
    CCTK_REAL const T1000039  =  -2*PDstandardNth23g23;
    
    g11rhsL  =  -2*alphaL*K11L + 2*(g11L*PDstandardNth1beta1 + 
           g12L*PDstandardNth1beta2 + g13L*PDstandardNth1beta3) + 
        beta1L*PDstandardNth1g11 + beta2L*PDstandardNth2g11 + 
        beta3L*PDstandardNth3g11;
    
    g12rhsL  =  -2*alphaL*K12L + g22L*PDstandardNth1beta2 + 
        g23L*PDstandardNth1beta3 + beta1L*PDstandardNth1g12 + 
        g11L*PDstandardNth2beta1 + g12L*
         (PDstandardNth1beta1 + PDstandardNth2beta2) + g13L*PDstandardNth2beta3 + 
        beta2L*PDstandardNth2g12 + beta3L*PDstandardNth3g12;
    
    g13rhsL  =  -2*alphaL*K13L + g23L*PDstandardNth1beta2 + 
        g33L*PDstandardNth1beta3 + beta1L*PDstandardNth1g13 + 
        beta2L*PDstandardNth2g13 + g11L*PDstandardNth3beta1 + 
        g12L*PDstandardNth3beta2 + g13L*
         (PDstandardNth1beta1 + PDstandardNth3beta3) + beta3L*PDstandardNth3g13;
    
    g22rhsL  =  -2*alphaL*K22L + beta1L*PDstandardNth1g22 + 
        2*(g12L*PDstandardNth2beta1 + g22L*PDstandardNth2beta2 + 
           g23L*PDstandardNth2beta3) + beta2L*PDstandardNth2g22 + 
        beta3L*PDstandardNth3g22;
    
    g23rhsL  =  -2*alphaL*K23L + beta1L*PDstandardNth1g23 + 
        g13L*PDstandardNth2beta1 + g33L*PDstandardNth2beta3 + 
        beta2L*PDstandardNth2g23 + g12L*PDstandardNth3beta1 + 
        g22L*PDstandardNth3beta2 + g23L*
         (PDstandardNth2beta2 + PDstandardNth3beta3) + beta3L*PDstandardNth3g23;
    
    g33rhsL  =  -2*alphaL*K33L + beta1L*PDstandardNth1g33 + 
        beta2L*PDstandardNth2g33 + 2*
         (g13L*PDstandardNth3beta1 + g23L*PDstandardNth3beta2 + 
           g33L*PDstandardNth3beta3) + beta3L*PDstandardNth3g33;
    
    alpharhsL  =  0;
    
    beta1rhsL  =  0;
    
    beta2rhsL  =  0;
    
    beta3rhsL  =  0;
    
    CCTK_REAL const T1000005  =  T1000003 + T1000004;
    
    CCTK_REAL const T1000010  =  T1000008 + T1000009;
    
    CCTK_REAL const T1000014  =  T1000012 + T1000013;
    
    CCTK_REAL const T1000018  =  PDstandardNth1g23 + PDstandardNth2g13 + T1000017;
    
    CCTK_REAL const T1000020  =  PDstandardNth1g23 + PDstandardNth3g12 + T1000019;
    
    CCTK_REAL const T1000022  =  PDstandardNth2g13 + PDstandardNth3g12 + T1000021;
    
    CCTK_REAL const T1000025  =  
       PDstandardNth12g13 + PDstandardNth13g12 + T1000023 + T1000024;
    
    CCTK_REAL const T1000030  =  PDstandardNth11g22 + PDstandardNth22g11 + T1000029;
    
    CCTK_REAL const T1000033  =  
       PDstandardNth12g23 + PDstandardNth23g12 + T1000031 + T1000032;
    
    CCTK_REAL const T1000035  =  PDstandardNth11g33 + PDstandardNth33g11 + T1000034;
    
    CCTK_REAL const T1000038  =  
       PDstandardNth13g23 + PDstandardNth23g13 + T1000036 + T1000037;
    
    CCTK_REAL const T1000040  =  PDstandardNth22g33 + PDstandardNth33g22 + T1000039;
    
    detg  =  2*g12L*g13L*g23L + g11L*(g22L*g33L - T1000001) - g22L*T1000007 - 
        g33L*T1000016;
    
    CCTK_REAL const T1000002  =  INV(detg);
    
    CCTK_REAL const T1000006  =  T1000002*T1000005;
    
    CCTK_REAL const T1000011  =  T1000002*T1000010;
    
    CCTK_REAL const T1000015  =  T1000002*T1000014;
    
    gu11  =  (g22L*g33L - T1000001)*T1000002;
    
    gu22  =  T1000002*(g11L*g33L - T1000007);
    
    gu33  =  T1000002*(g11L*g22L - T1000016);
    
    gu12  =  T1000006;
    
    gu13  =  T1000011;
    
    gu21  =  T1000006;
    
    gu23  =  T1000015;
    
    gu31  =  T1000011;
    
    gu32  =  T1000015;
    
    G111  =  khalf*(gu11*PDstandardNth1g11 + 
          2*(gu12*PDstandardNth1g12 + gu13*PDstandardNth1g13) - 
          gu12*PDstandardNth2g11 - gu13*PDstandardNth3g11);
    
    G211  =  khalf*(gu21*PDstandardNth1g11 + 
          2*(gu22*PDstandardNth1g12 + gu23*PDstandardNth1g13) - 
          gu22*PDstandardNth2g11 - gu23*PDstandardNth3g11);
    
    G311  =  khalf*(gu31*PDstandardNth1g11 + 
          2*(gu32*PDstandardNth1g12 + gu33*PDstandardNth1g13) - 
          gu32*PDstandardNth2g11 - gu33*PDstandardNth3g11);
    
    G112  =  khalf*(gu12*PDstandardNth1g22 + gu11*PDstandardNth2g11 + gu13*T1000018);
    
    G212  =  khalf*(gu22*PDstandardNth1g22 + gu21*PDstandardNth2g11 + gu23*T1000018);
    
    G312  =  khalf*(gu32*PDstandardNth1g22 + gu31*PDstandardNth2g11 + gu33*T1000018);
    
    G113  =  khalf*(gu13*PDstandardNth1g33 + gu11*PDstandardNth3g11 + gu12*T1000020);
    
    G213  =  khalf*(gu23*PDstandardNth1g33 + gu21*PDstandardNth3g11 + gu22*T1000020);
    
    G313  =  khalf*(gu33*PDstandardNth1g33 + gu31*PDstandardNth3g11 + gu32*T1000020);
    
    G122  =  khalf*(gu11*(-PDstandardNth1g22 + 2*PDstandardNth2g12) + 
          gu12*PDstandardNth2g22 + gu13*(2*PDstandardNth2g23 - PDstandardNth3g22));
    
    G222  =  khalf*(gu21*(-PDstandardNth1g22 + 2*PDstandardNth2g12) + 
          gu22*PDstandardNth2g22 + gu23*(2*PDstandardNth2g23 - PDstandardNth3g22));
    
    G322  =  khalf*(gu31*(-PDstandardNth1g22 + 2*PDstandardNth2g12) + 
          gu32*PDstandardNth2g22 + gu33*(2*PDstandardNth2g23 - PDstandardNth3g22));
    
    G123  =  khalf*(gu13*PDstandardNth2g33 + gu12*PDstandardNth3g22 + gu11*T1000022);
    
    G223  =  khalf*(gu23*PDstandardNth2g33 + gu22*PDstandardNth3g22 + gu21*T1000022);
    
    G323  =  khalf*(gu33*PDstandardNth2g33 + gu32*PDstandardNth3g22 + gu31*T1000022);
    
    G133  =  khalf*(-(gu11*PDstandardNth1g33) - gu12*PDstandardNth2g33 + 
          2*gu11*PDstandardNth3g13 + 2*gu12*PDstandardNth3g23 + 
          gu13*PDstandardNth3g33);
    
    G233  =  khalf*(-(gu21*PDstandardNth1g33) - gu22*PDstandardNth2g33 + 
          2*gu21*PDstandardNth3g13 + 2*gu22*PDstandardNth3g23 + 
          gu23*PDstandardNth3g33);
    
    G333  =  khalf*(-(gu31*PDstandardNth1g33) - gu32*PDstandardNth2g33 + 
          2*gu31*PDstandardNth3g13 + 2*gu32*PDstandardNth3g23 + 
          gu33*PDstandardNth3g33);
    
    Km11  =  gu11*K11L + gu12*K12L + gu13*K13L;
    
    Km21  =  gu21*K11L + gu22*K12L + gu23*K13L;
    
    Km31  =  gu31*K11L + gu32*K12L + gu33*K13L;
    
    Km12  =  gu11*K12L + gu12*K22L + gu13*K23L;
    
    Km22  =  gu21*K12L + gu22*K22L + gu23*K23L;
    
    Km32  =  gu31*K12L + gu32*K22L + gu33*K23L;
    
    Km13  =  gu11*K13L + gu12*K23L + gu13*K33L;
    
    Km23  =  gu21*K13L + gu22*K23L + gu23*K33L;
    
    Km33  =  gu31*K13L + gu32*K23L + gu33*K33L;
    
    CCTK_REAL const T1000026  =  G111 + G212 + G313;
    
    CCTK_REAL const T1000027  =  G112 + G222 + G323;
    
    CCTK_REAL const T1000028  =  G113 + G223 + G333;
    
    R12  =  khalf*(2*(G122*G211 + G123*G311 + G213*G322 + G313*G323) - 
          2*(G112*G212 + G112*G313 + G212*G323 + G312*G333 + 
             gu12*PDstandardNth12g12) - gu32*PDstandardNth12g23 - 
          gu33*PDstandardNth12g33 + gu13*
           (PDstandardNth11g23 - PDstandardNth12g13 - PDstandardNth13g12) + 
          gu32*PDstandardNth13g22 + gu33*PDstandardNth13g23 + 
          gu12*(PDstandardNth11g22 + PDstandardNth22g11) + 
          gu32*PDstandardNth22g13 + gu13*PDstandardNth23g11 - 
          gu32*PDstandardNth23g12 + gu33*PDstandardNth23g13 - 
          gu33*PDstandardNth33g12);
    
    R13  =  khalf*(2*(G123*G211 + G212*G223 + G133*G311 + G233*G312) - 
          2*(G213*G222 + G223*G313 + G113*(G212 + G313) + 
             gu13*PDstandardNth13g13) + 
          gu12*(PDstandardNth11g23 - PDstandardNth12g13 - PDstandardNth13g12 + 
             PDstandardNth23g11) + gu22*
           (PDstandardNth12g23 - PDstandardNth13g22 - PDstandardNth22g13 + 
             PDstandardNth23g12) + gu13*(PDstandardNth11g33 + PDstandardNth33g11) + 
          gu23*(PDstandardNth12g33 - PDstandardNth13g23 - PDstandardNth23g13 + 
             PDstandardNth33g12));
    
    R23  =  khalf*(2*(G112*G113 + G122*G213 + G133*G312 + G233*G322) + 
          gu11*(-PDstandardNth11g23 + PDstandardNth12g13 + PDstandardNth13g12 - 
             PDstandardNth23g11) + gu21*
           (-PDstandardNth12g23 + PDstandardNth13g22 + PDstandardNth22g13 - 
             PDstandardNth23g12) - 2*
           (G111*G123 + G113*G323 + G223*(G112 + G323) + gu23*PDstandardNth23g23) + 
          gu13*(PDstandardNth12g33 - PDstandardNth13g23 - PDstandardNth23g13 + 
             PDstandardNth33g12) + gu23*(PDstandardNth22g33 + PDstandardNth33g22));
    
    trK  =  Km11 + Km22 + Km33;
    
    R11  =  2*(G112*G211 + G113*G311 + G213*G312) - G111*T1000026 - G211*T1000027 - 
        G311*T1000028 + khalf*(gu23*T1000025 + gu32*T1000025 - gu22*T1000030 - 
           gu33*T1000035) + SQR(G111) + SQR(G212) + SQR(G313);
    
    R22  =  2*(G122*G212 + G123*G312 + G223*G322) - G122*T1000026 - G222*T1000027 - 
        G322*T1000028 + khalf*(-(gu11*T1000030) + gu13*T1000033 + gu31*T1000033 - 
           gu33*T1000040) + SQR(G112) + SQR(G222) + SQR(G323);
    
    R33  =  2*(G123*G213 + G133*G313 + G233*G323) - G133*T1000026 - G233*T1000027 - 
        G333*T1000028 + khalf*(-(gu11*T1000035) + gu12*T1000038 + gu21*T1000038 - 
           gu22*T1000040) + SQR(G113) + SQR(G223) + SQR(G333);
    
    K12rhsL  =  -PDstandardNth12alpha + G112*PDstandardNth1alpha + 
        K22L*PDstandardNth1beta2 + K23L*PDstandardNth1beta3 + 
        beta1L*PDstandardNth1K12 + G212*PDstandardNth2alpha + 
        K11L*PDstandardNth2beta1 + K12L*
         (PDstandardNth1beta1 + PDstandardNth2beta2) + K13L*PDstandardNth2beta3 + 
        beta2L*PDstandardNth2K12 + G312*PDstandardNth3alpha + 
        beta3L*PDstandardNth3K12 + alphaL*
         (-2*(K11L*Km12 + K12L*Km22 + K13L*Km32) + R12 + K12L*trK);
    
    K13rhsL  =  -PDstandardNth13alpha + G113*PDstandardNth1alpha + 
        K23L*PDstandardNth1beta2 + K33L*PDstandardNth1beta3 + 
        beta1L*PDstandardNth1K13 + G213*PDstandardNth2alpha + 
        beta2L*PDstandardNth2K13 + G313*PDstandardNth3alpha + 
        K11L*PDstandardNth3beta1 + K12L*PDstandardNth3beta2 + 
        K13L*(PDstandardNth1beta1 + PDstandardNth3beta3) + 
        beta3L*PDstandardNth3K13 + alphaL*
         (-2*(K11L*Km13 + K12L*Km23 + K13L*Km33) + R13 + K13L*trK);
    
    K23rhsL  =  G123*PDstandardNth1alpha + beta1L*PDstandardNth1K23 - 
        PDstandardNth23alpha + G223*PDstandardNth2alpha + 
        K13L*PDstandardNth2beta1 + K33L*PDstandardNth2beta3 + 
        beta2L*PDstandardNth2K23 + G323*PDstandardNth3alpha + 
        K12L*PDstandardNth3beta1 + K22L*PDstandardNth3beta2 + 
        K23L*(PDstandardNth2beta2 + PDstandardNth3beta3) + 
        beta3L*PDstandardNth3K23 + alphaL*
         (-2*(K12L*Km13 + K22L*Km23 + K23L*Km33) + R23 + K23L*trK);
    
    K11rhsL  =  -PDstandardNth11alpha + G111*PDstandardNth1alpha + 
        2*(K11L*PDstandardNth1beta1 + K12L*PDstandardNth1beta2 + 
           K13L*PDstandardNth1beta3) + beta1L*PDstandardNth1K11 + 
        G211*PDstandardNth2alpha + beta2L*PDstandardNth2K11 + 
        G311*PDstandardNth3alpha + beta3L*PDstandardNth3K11 + 
        alphaL*(-2*(K11L*Km11 + K12L*Km21 + K13L*Km31) + R11 + K11L*trK);
    
    K22rhsL  =  G122*PDstandardNth1alpha + beta1L*PDstandardNth1K22 - 
        PDstandardNth22alpha + G222*PDstandardNth2alpha + 
        2*(K12L*PDstandardNth2beta1 + K22L*PDstandardNth2beta2 + 
           K23L*PDstandardNth2beta3) + beta2L*PDstandardNth2K22 + 
        G322*PDstandardNth3alpha + beta3L*PDstandardNth3K22 + 
        alphaL*(-2*(K12L*Km12 + K22L*Km22 + K23L*Km32) + R22 + K22L*trK);
    
    K33rhsL  =  G133*PDstandardNth1alpha + beta1L*PDstandardNth1K33 + 
        G233*PDstandardNth2alpha + beta2L*PDstandardNth2K33 - 
        PDstandardNth33alpha + G333*PDstandardNth3alpha + 
        2*(K13L*PDstandardNth3beta1 + K23L*PDstandardNth3beta2 + 
           K33L*PDstandardNth3beta3) + beta3L*PDstandardNth3K33 + 
        alphaL*(-2*(K13L*Km13 + K23L*Km23 + K33L*Km33) + R33 + K33L*trK);
    
    
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
  LC_ENDLOOP3 (ML_ADM_RHS);
}

void ML_ADM_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADM_RHS_Body);
}
