/*  File produced by Kranc */

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

void ML_ADM_constraints_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  /* Declare predefined quantities */
  // CCTK_REAL p1o12dx = INITVALUE;
  // CCTK_REAL p1o12dy = INITVALUE;
  // CCTK_REAL p1o12dz = INITVALUE;
  // CCTK_REAL p1o144dxdy = INITVALUE;
  // CCTK_REAL p1o144dxdz = INITVALUE;
  // CCTK_REAL p1o144dydz = INITVALUE;
  // CCTK_REAL pm1o12dx2 = INITVALUE;
  // CCTK_REAL pm1o12dy2 = INITVALUE;
  // CCTK_REAL pm1o12dz2 = INITVALUE;
  
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
  CCTK_REAL const dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL const dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL const dz = CCTK_DELTA_SPACE(2);
  int const di = 1;
  int const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  int const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  CCTK_REAL const dxi = 1.0 / dx;
  CCTK_REAL const dyi = 1.0 / dy;
  CCTK_REAL const dzi = 1.0 / dz;
  CCTK_REAL const khalf = 0.5;
  CCTK_REAL const kthird = 1/3.0;
  CCTK_REAL const ktwothird = 2.0/3.0;
  CCTK_REAL const kfourthird = 4.0/3.0;
  CCTK_REAL const keightthird = 8.0/3.0;
  CCTK_REAL const hdxi = 0.5 * dxi;
  CCTK_REAL const hdyi = 0.5 * dyi;
  CCTK_REAL const hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  CCTK_REAL const p1o12dx = INV(dx)/12.;
  CCTK_REAL const p1o12dy = INV(dy)/12.;
  CCTK_REAL const p1o12dz = INV(dz)/12.;
  CCTK_REAL const p1o144dxdy = (INV(dx)*INV(dy))/144.;
  CCTK_REAL const p1o144dxdz = (INV(dx)*INV(dz))/144.;
  CCTK_REAL const p1o144dydz = (INV(dy)*INV(dz))/144.;
  CCTK_REAL const pm1o12dx2 = -pow(dx,-2)/12.;
  CCTK_REAL const pm1o12dy2 = -pow(dy,-2)/12.;
  CCTK_REAL const pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_ADM_constraints,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    // int subblock_index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int const subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    // CCTK_REAL detg = INITVALUE;
    // CCTK_REAL G111 = INITVALUE, G112 = INITVALUE, G113 = INITVALUE, G122 = INITVALUE, G123 = INITVALUE, G133 = INITVALUE;
    // CCTK_REAL G211 = INITVALUE, G212 = INITVALUE, G213 = INITVALUE, G222 = INITVALUE, G223 = INITVALUE, G233 = INITVALUE;
    // CCTK_REAL G311 = INITVALUE, G312 = INITVALUE, G313 = INITVALUE, G322 = INITVALUE, G323 = INITVALUE, G333 = INITVALUE;
    // CCTK_REAL gu11 = INITVALUE, gu12 = INITVALUE, gu13 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu23 = INITVALUE;
    // CCTK_REAL gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    // CCTK_REAL Km11 = INITVALUE, Km12 = INITVALUE, Km13 = INITVALUE, Km21 = INITVALUE, Km22 = INITVALUE, Km23 = INITVALUE;
    // CCTK_REAL Km31 = INITVALUE, Km32 = INITVALUE, Km33 = INITVALUE;
    // CCTK_REAL R11 = INITVALUE, R12 = INITVALUE, R13 = INITVALUE, R22 = INITVALUE, R23 = INITVALUE, R33 = INITVALUE;
    // CCTK_REAL trK = INITVALUE;
    // CCTK_REAL trR = INITVALUE;
    
    /* Declare local copies of grid functions */
    // CCTK_REAL g11L = INITVALUE, g12L = INITVALUE, g13L = INITVALUE, g22L = INITVALUE, g23L = INITVALUE, g33L = INITVALUE;
    // CCTK_REAL HL = INITVALUE;
    // CCTK_REAL K11L = INITVALUE, K12L = INITVALUE, K13L = INITVALUE, K22L = INITVALUE, K23L = INITVALUE, K33L = INITVALUE;
    // CCTK_REAL M1L = INITVALUE, M2L = INITVALUE, M3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    // CCTK_REAL PDstandardNth1g11 = INITVALUE;
    // CCTK_REAL PDstandardNth2g11 = INITVALUE;
    // CCTK_REAL PDstandardNth3g11 = INITVALUE;
    // CCTK_REAL PDstandardNth22g11 = INITVALUE;
    // CCTK_REAL PDstandardNth33g11 = INITVALUE;
    // CCTK_REAL PDstandardNth23g11 = INITVALUE;
    // CCTK_REAL PDstandardNth1g12 = INITVALUE;
    // CCTK_REAL PDstandardNth2g12 = INITVALUE;
    // CCTK_REAL PDstandardNth3g12 = INITVALUE;
    // CCTK_REAL PDstandardNth33g12 = INITVALUE;
    // CCTK_REAL PDstandardNth12g12 = INITVALUE;
    // CCTK_REAL PDstandardNth13g12 = INITVALUE;
    // CCTK_REAL PDstandardNth23g12 = INITVALUE;
    // CCTK_REAL PDstandardNth1g13 = INITVALUE;
    // CCTK_REAL PDstandardNth2g13 = INITVALUE;
    // CCTK_REAL PDstandardNth3g13 = INITVALUE;
    // CCTK_REAL PDstandardNth22g13 = INITVALUE;
    // CCTK_REAL PDstandardNth12g13 = INITVALUE;
    // CCTK_REAL PDstandardNth13g13 = INITVALUE;
    // CCTK_REAL PDstandardNth23g13 = INITVALUE;
    // CCTK_REAL PDstandardNth1g22 = INITVALUE;
    // CCTK_REAL PDstandardNth2g22 = INITVALUE;
    // CCTK_REAL PDstandardNth3g22 = INITVALUE;
    // CCTK_REAL PDstandardNth11g22 = INITVALUE;
    // CCTK_REAL PDstandardNth33g22 = INITVALUE;
    // CCTK_REAL PDstandardNth13g22 = INITVALUE;
    // CCTK_REAL PDstandardNth1g23 = INITVALUE;
    // CCTK_REAL PDstandardNth2g23 = INITVALUE;
    // CCTK_REAL PDstandardNth3g23 = INITVALUE;
    // CCTK_REAL PDstandardNth11g23 = INITVALUE;
    // CCTK_REAL PDstandardNth12g23 = INITVALUE;
    // CCTK_REAL PDstandardNth13g23 = INITVALUE;
    // CCTK_REAL PDstandardNth23g23 = INITVALUE;
    // CCTK_REAL PDstandardNth1g33 = INITVALUE;
    // CCTK_REAL PDstandardNth2g33 = INITVALUE;
    // CCTK_REAL PDstandardNth3g33 = INITVALUE;
    // CCTK_REAL PDstandardNth11g33 = INITVALUE;
    // CCTK_REAL PDstandardNth22g33 = INITVALUE;
    // CCTK_REAL PDstandardNth12g33 = INITVALUE;
    // CCTK_REAL PDstandardNth2K11 = INITVALUE;
    // CCTK_REAL PDstandardNth3K11 = INITVALUE;
    // CCTK_REAL PDstandardNth1K12 = INITVALUE;
    // CCTK_REAL PDstandardNth2K12 = INITVALUE;
    // CCTK_REAL PDstandardNth3K12 = INITVALUE;
    // CCTK_REAL PDstandardNth1K13 = INITVALUE;
    // CCTK_REAL PDstandardNth2K13 = INITVALUE;
    // CCTK_REAL PDstandardNth3K13 = INITVALUE;
    // CCTK_REAL PDstandardNth1K22 = INITVALUE;
    // CCTK_REAL PDstandardNth3K22 = INITVALUE;
    // CCTK_REAL PDstandardNth1K23 = INITVALUE;
    // CCTK_REAL PDstandardNth2K23 = INITVALUE;
    // CCTK_REAL PDstandardNth3K23 = INITVALUE;
    // CCTK_REAL PDstandardNth1K33 = INITVALUE;
    // CCTK_REAL PDstandardNth2K33 = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL const g11L = g11[index];
    CCTK_REAL const g12L = g12[index];
    CCTK_REAL const g13L = g13[index];
    CCTK_REAL const g22L = g22[index];
    CCTK_REAL const g23L = g23[index];
    CCTK_REAL const g33L = g33[index];
    CCTK_REAL const K11L = K11[index];
    CCTK_REAL const K12L = K12[index];
    CCTK_REAL const K13L = K13[index];
    CCTK_REAL const K22L = K22[index];
    CCTK_REAL const K23L = K23[index];
    CCTK_REAL const K33L = K33[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    CCTK_REAL const PDstandardNth1g11 = PDstandardNth1(g11, i, j, k);
    CCTK_REAL const PDstandardNth2g11 = PDstandardNth2(g11, i, j, k);
    CCTK_REAL const PDstandardNth3g11 = PDstandardNth3(g11, i, j, k);
    CCTK_REAL const PDstandardNth22g11 = PDstandardNth22(g11, i, j, k);
    CCTK_REAL const PDstandardNth33g11 = PDstandardNth33(g11, i, j, k);
    CCTK_REAL const PDstandardNth23g11 = PDstandardNth23(g11, i, j, k);
    CCTK_REAL const PDstandardNth1g12 = PDstandardNth1(g12, i, j, k);
    CCTK_REAL const PDstandardNth2g12 = PDstandardNth2(g12, i, j, k);
    CCTK_REAL const PDstandardNth3g12 = PDstandardNth3(g12, i, j, k);
    CCTK_REAL const PDstandardNth33g12 = PDstandardNth33(g12, i, j, k);
    CCTK_REAL const PDstandardNth12g12 = PDstandardNth12(g12, i, j, k);
    CCTK_REAL const PDstandardNth13g12 = PDstandardNth13(g12, i, j, k);
    CCTK_REAL const PDstandardNth23g12 = PDstandardNth23(g12, i, j, k);
    CCTK_REAL const PDstandardNth1g13 = PDstandardNth1(g13, i, j, k);
    CCTK_REAL const PDstandardNth2g13 = PDstandardNth2(g13, i, j, k);
    CCTK_REAL const PDstandardNth3g13 = PDstandardNth3(g13, i, j, k);
    CCTK_REAL const PDstandardNth22g13 = PDstandardNth22(g13, i, j, k);
    CCTK_REAL const PDstandardNth12g13 = PDstandardNth12(g13, i, j, k);
    CCTK_REAL const PDstandardNth13g13 = PDstandardNth13(g13, i, j, k);
    CCTK_REAL const PDstandardNth23g13 = PDstandardNth23(g13, i, j, k);
    CCTK_REAL const PDstandardNth1g22 = PDstandardNth1(g22, i, j, k);
    CCTK_REAL const PDstandardNth2g22 = PDstandardNth2(g22, i, j, k);
    CCTK_REAL const PDstandardNth3g22 = PDstandardNth3(g22, i, j, k);
    CCTK_REAL const PDstandardNth11g22 = PDstandardNth11(g22, i, j, k);
    CCTK_REAL const PDstandardNth33g22 = PDstandardNth33(g22, i, j, k);
    CCTK_REAL const PDstandardNth13g22 = PDstandardNth13(g22, i, j, k);
    CCTK_REAL const PDstandardNth1g23 = PDstandardNth1(g23, i, j, k);
    CCTK_REAL const PDstandardNth2g23 = PDstandardNth2(g23, i, j, k);
    CCTK_REAL const PDstandardNth3g23 = PDstandardNth3(g23, i, j, k);
    CCTK_REAL const PDstandardNth11g23 = PDstandardNth11(g23, i, j, k);
    CCTK_REAL const PDstandardNth12g23 = PDstandardNth12(g23, i, j, k);
    CCTK_REAL const PDstandardNth13g23 = PDstandardNth13(g23, i, j, k);
    CCTK_REAL const PDstandardNth23g23 = PDstandardNth23(g23, i, j, k);
    CCTK_REAL const PDstandardNth1g33 = PDstandardNth1(g33, i, j, k);
    CCTK_REAL const PDstandardNth2g33 = PDstandardNth2(g33, i, j, k);
    CCTK_REAL const PDstandardNth3g33 = PDstandardNth3(g33, i, j, k);
    CCTK_REAL const PDstandardNth11g33 = PDstandardNth11(g33, i, j, k);
    CCTK_REAL const PDstandardNth22g33 = PDstandardNth22(g33, i, j, k);
    CCTK_REAL const PDstandardNth12g33 = PDstandardNth12(g33, i, j, k);
    CCTK_REAL const PDstandardNth2K11 = PDstandardNth2(K11, i, j, k);
    CCTK_REAL const PDstandardNth3K11 = PDstandardNth3(K11, i, j, k);
    CCTK_REAL const PDstandardNth1K12 = PDstandardNth1(K12, i, j, k);
    CCTK_REAL const PDstandardNth2K12 = PDstandardNth2(K12, i, j, k);
    CCTK_REAL const PDstandardNth3K12 = PDstandardNth3(K12, i, j, k);
    CCTK_REAL const PDstandardNth1K13 = PDstandardNth1(K13, i, j, k);
    CCTK_REAL const PDstandardNth2K13 = PDstandardNth2(K13, i, j, k);
    CCTK_REAL const PDstandardNth3K13 = PDstandardNth3(K13, i, j, k);
    CCTK_REAL const PDstandardNth1K22 = PDstandardNth1(K22, i, j, k);
    CCTK_REAL const PDstandardNth3K22 = PDstandardNth3(K22, i, j, k);
    CCTK_REAL const PDstandardNth1K23 = PDstandardNth1(K23, i, j, k);
    CCTK_REAL const PDstandardNth2K23 = PDstandardNth2(K23, i, j, k);
    CCTK_REAL const PDstandardNth3K23 = PDstandardNth3(K23, i, j, k);
    CCTK_REAL const PDstandardNth1K33 = PDstandardNth1(K33, i, j, k);
    CCTK_REAL const PDstandardNth2K33 = PDstandardNth2(K33, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL const detg  =  2*g12L*g13L*g23L + g33L*(g11L*g22L - SQR(g12L)) - g22L*SQR(g13L) - g11L*SQR(g23L);
    
    CCTK_REAL const gu11  =  INV(detg)*(g22L*g33L - SQR(g23L));
    
    CCTK_REAL const gu12  =  (g13L*g23L - g12L*g33L)*INV(detg);
    
    CCTK_REAL const gu13  =  (-(g13L*g22L) + g12L*g23L)*INV(detg);
    
    CCTK_REAL const gu21  =  (g13L*g23L - g12L*g33L)*INV(detg);
    
    CCTK_REAL const gu22  =  INV(detg)*(g11L*g33L - SQR(g13L));
    
    CCTK_REAL const gu23  =  (g12L*g13L - g11L*g23L)*INV(detg);
    
    CCTK_REAL const gu31  =  (-(g13L*g22L) + g12L*g23L)*INV(detg);
    
    CCTK_REAL const gu32  =  (g12L*g13L - g11L*g23L)*INV(detg);
    
    CCTK_REAL const gu33  =  INV(detg)*(g11L*g22L - SQR(g12L));
    
    CCTK_REAL const G111  =  khalf*(gu11*PDstandardNth1g11 + 2*(gu12*PDstandardNth1g12 + gu13*PDstandardNth1g13) - gu12*PDstandardNth2g11 - 
          gu13*PDstandardNth3g11);
    
    CCTK_REAL const G211  =  khalf*(gu21*PDstandardNth1g11 + 2*(gu22*PDstandardNth1g12 + gu23*PDstandardNth1g13) - gu22*PDstandardNth2g11 - 
          gu23*PDstandardNth3g11);
    
    CCTK_REAL const G311  =  khalf*(gu31*PDstandardNth1g11 + 2*(gu32*PDstandardNth1g12 + gu33*PDstandardNth1g13) - gu32*PDstandardNth2g11 - 
          gu33*PDstandardNth3g11);
    
    CCTK_REAL const G112  =  khalf*(gu12*PDstandardNth1g22 + gu11*PDstandardNth2g11 + 
          gu13*(PDstandardNth1g23 + PDstandardNth2g13 - PDstandardNth3g12));
    
    CCTK_REAL const G212  =  khalf*(gu22*PDstandardNth1g22 + gu21*PDstandardNth2g11 + 
          gu23*(PDstandardNth1g23 + PDstandardNth2g13 - PDstandardNth3g12));
    
    CCTK_REAL const G312  =  khalf*(gu32*PDstandardNth1g22 + gu31*PDstandardNth2g11 + 
          gu33*(PDstandardNth1g23 + PDstandardNth2g13 - PDstandardNth3g12));
    
    CCTK_REAL const G113  =  khalf*(gu13*PDstandardNth1g33 + gu11*PDstandardNth3g11 + 
          gu12*(PDstandardNth1g23 - PDstandardNth2g13 + PDstandardNth3g12));
    
    CCTK_REAL const G213  =  khalf*(gu23*PDstandardNth1g33 + gu21*PDstandardNth3g11 + 
          gu22*(PDstandardNth1g23 - PDstandardNth2g13 + PDstandardNth3g12));
    
    CCTK_REAL const G313  =  khalf*(gu33*PDstandardNth1g33 + gu31*PDstandardNth3g11 + 
          gu32*(PDstandardNth1g23 - PDstandardNth2g13 + PDstandardNth3g12));
    
    CCTK_REAL const G122  =  khalf*(gu11*(-PDstandardNth1g22 + 2*PDstandardNth2g12) + gu12*PDstandardNth2g22 + 
          gu13*(2*PDstandardNth2g23 - PDstandardNth3g22));
    
    CCTK_REAL const G222  =  khalf*(gu21*(-PDstandardNth1g22 + 2*PDstandardNth2g12) + gu22*PDstandardNth2g22 + 
          gu23*(2*PDstandardNth2g23 - PDstandardNth3g22));
    
    CCTK_REAL const G322  =  khalf*(gu31*(-PDstandardNth1g22 + 2*PDstandardNth2g12) + gu32*PDstandardNth2g22 + 
          gu33*(2*PDstandardNth2g23 - PDstandardNth3g22));
    
    CCTK_REAL const G123  =  khalf*(gu13*PDstandardNth2g33 + gu11*(-PDstandardNth1g23 + PDstandardNth2g13 + PDstandardNth3g12) + 
          gu12*PDstandardNth3g22);
    
    CCTK_REAL const G223  =  khalf*(gu23*PDstandardNth2g33 + gu21*(-PDstandardNth1g23 + PDstandardNth2g13 + PDstandardNth3g12) + 
          gu22*PDstandardNth3g22);
    
    CCTK_REAL const G323  =  khalf*(gu33*PDstandardNth2g33 + gu31*(-PDstandardNth1g23 + PDstandardNth2g13 + PDstandardNth3g12) + 
          gu32*PDstandardNth3g22);
    
    CCTK_REAL const G133  =  khalf*(-(gu11*PDstandardNth1g33) - gu12*PDstandardNth2g33 + 2*gu11*PDstandardNth3g13 + 
          2*gu12*PDstandardNth3g23 + gu13*PDstandardNth3g33);
    
    CCTK_REAL const G233  =  khalf*(-(gu21*PDstandardNth1g33) - gu22*PDstandardNth2g33 + 2*gu21*PDstandardNth3g13 + 
          2*gu22*PDstandardNth3g23 + gu23*PDstandardNth3g33);
    
    CCTK_REAL const G333  =  khalf*(-(gu31*PDstandardNth1g33) - gu32*PDstandardNth2g33 + 2*gu31*PDstandardNth3g13 + 
          2*gu32*PDstandardNth3g23 + gu33*PDstandardNth3g33);
    
    CCTK_REAL const R11  =  -(G111*(G111 + G122 + G133)) - G211*(G211 + G222 + G233) - G311*(G311 + G322 + G333) + 
        khalf*(-(gu22*(PDstandardNth11g22 - 2*PDstandardNth12g12 + PDstandardNth22g11)) + 
           gu23*(-PDstandardNth11g23 + PDstandardNth12g13 + PDstandardNth13g12 - PDstandardNth23g11) + 
           gu32*(-PDstandardNth11g23 + PDstandardNth12g13 + PDstandardNth13g12 - PDstandardNth23g11) - 
           gu33*(PDstandardNth11g33 - 2*PDstandardNth13g13 + PDstandardNth33g11)) + SQR(G111) + SQR(G112) + SQR(G113) + 
        SQR(G211) + SQR(G212) + SQR(G213) + SQR(G311) + SQR(G312) + SQR(G313);
    
    CCTK_REAL const R12  =  khalf*(2*(G113*G123 + G213*G223 + G313*G323) - 2*(G112*G133 + G212*G233 + G312*G333 + gu12*PDstandardNth12g12) - 
          gu32*PDstandardNth12g23 - gu33*PDstandardNth12g33 + 
          gu13*(PDstandardNth11g23 - PDstandardNth12g13 - PDstandardNth13g12) + gu32*PDstandardNth13g22 + 
          gu33*PDstandardNth13g23 + gu12*(PDstandardNth11g22 + PDstandardNth22g11) + gu32*PDstandardNth22g13 + 
          gu13*PDstandardNth23g11 - gu32*PDstandardNth23g12 + gu33*PDstandardNth23g13 - gu33*PDstandardNth33g12);
    
    CCTK_REAL const R13  =  khalf*(2*(G112*G123 + G212*G223 + G312*G323) - 2*(G113*G122 + G213*G222 + G313*G322 + gu13*PDstandardNth13g13) + 
          gu12*(PDstandardNth11g23 - PDstandardNth12g13 - PDstandardNth13g12 + PDstandardNth23g11) + 
          gu22*(PDstandardNth12g23 - PDstandardNth13g22 - PDstandardNth22g13 + PDstandardNth23g12) + 
          gu13*(PDstandardNth11g33 + PDstandardNth33g11) + 
          gu23*(PDstandardNth12g33 - PDstandardNth13g23 - PDstandardNth23g13 + PDstandardNth33g12));
    
    CCTK_REAL const R22  =  -(G122*(G111 + G122 + G133)) - G222*(G211 + G222 + G233) - G322*(G311 + G322 + G333) + 
        khalf*(-(gu11*(PDstandardNth11g22 - 2*PDstandardNth12g12 + PDstandardNth22g11)) + 
           gu13*(PDstandardNth12g23 - PDstandardNth13g22 - PDstandardNth22g13 + PDstandardNth23g12) + 
           gu31*(PDstandardNth12g23 - PDstandardNth13g22 - PDstandardNth22g13 + PDstandardNth23g12) - 
           gu33*(PDstandardNth22g33 - 2*PDstandardNth23g23 + PDstandardNth33g22)) + SQR(G112) + SQR(G122) + SQR(G123) + 
        SQR(G212) + SQR(G222) + SQR(G223) + SQR(G312) + SQR(G322) + SQR(G323);
    
    CCTK_REAL const R23  =  khalf*(2*(G112*G113 + G212*G213 + G312*G313) + 
          gu11*(-PDstandardNth11g23 + PDstandardNth12g13 + PDstandardNth13g12 - PDstandardNth23g11) + 
          gu21*(-PDstandardNth12g23 + PDstandardNth13g22 + PDstandardNth22g13 - PDstandardNth23g12) - 
          2*(G111*G123 + G211*G223 + G311*G323 + gu23*PDstandardNth23g23) + 
          gu13*(PDstandardNth12g33 - PDstandardNth13g23 - PDstandardNth23g13 + PDstandardNth33g12) + 
          gu23*(PDstandardNth22g33 + PDstandardNth33g22));
    
    CCTK_REAL const R33  =  -(G133*(G111 + G122 + G133)) - G233*(G211 + G222 + G233) - G333*(G311 + G322 + G333) + 
        khalf*(-(gu11*(PDstandardNth11g33 - 2*PDstandardNth13g13 + PDstandardNth33g11)) + 
           gu12*(-PDstandardNth12g33 + PDstandardNth13g23 + PDstandardNth23g13 - PDstandardNth33g12) + 
           gu21*(-PDstandardNth12g33 + PDstandardNth13g23 + PDstandardNth23g13 - PDstandardNth33g12) - 
           gu22*(PDstandardNth22g33 - 2*PDstandardNth23g23 + PDstandardNth33g22)) + SQR(G113) + SQR(G123) + SQR(G133) + 
        SQR(G213) + SQR(G223) + SQR(G233) + SQR(G313) + SQR(G323) + SQR(G333);
    
    CCTK_REAL const trR  =  gu11*R11 + (gu12 + gu21)*R12 + (gu13 + gu31)*R13 + gu22*R22 + (gu23 + gu32)*R23 + gu33*R33;
    
    CCTK_REAL const Km11  =  gu11*K11L + gu12*K12L + gu13*K13L;
    
    CCTK_REAL const Km21  =  gu21*K11L + gu22*K12L + gu23*K13L;
    
    CCTK_REAL const Km31  =  gu31*K11L + gu32*K12L + gu33*K13L;
    
    CCTK_REAL const Km12  =  gu11*K12L + gu12*K22L + gu13*K23L;
    
    CCTK_REAL const Km22  =  gu21*K12L + gu22*K22L + gu23*K23L;
    
    CCTK_REAL const Km32  =  gu31*K12L + gu32*K22L + gu33*K23L;
    
    CCTK_REAL const Km13  =  gu11*K13L + gu12*K23L + gu13*K33L;
    
    CCTK_REAL const Km23  =  gu21*K13L + gu22*K23L + gu23*K33L;
    
    CCTK_REAL const Km33  =  gu31*K13L + gu32*K23L + gu33*K33L;
    
    CCTK_REAL const trK  =  Km11 + Km22 + Km33;
    
    CCTK_REAL const HL  =  -2*(Km12*Km21 + Km13*Km31 + Km23*Km32) + trR - SQR(Km11) - SQR(Km22) - SQR(Km33) + SQR(trK);
    
    CCTK_REAL const M1L  =  gu21*(-(G112*K11L) + G111*K12L - G212*K12L - G312*K13L + G211*K22L + G311*K23L - PDstandardNth1K12 + 
           PDstandardNth2K11) + gu22*(-(G122*K11L) + G112*K12L - G222*K12L - G322*K13L + G212*K22L + G312*K23L - 
           PDstandardNth1K22 + PDstandardNth2K12) + gu23*
         (-(G123*K11L) + G113*K12L - G223*K12L - G323*K13L + G213*K22L + G313*K23L - PDstandardNth1K23 + PDstandardNth2K13)\
         + gu31*(-(G113*K11L) - G213*K12L + G111*K13L - G313*K13L + G211*K23L + G311*K33L - PDstandardNth1K13 + 
           PDstandardNth3K11) + gu32*(-(G123*K11L) - G223*K12L + G112*K13L - G323*K13L + G212*K23L + G312*K33L - 
           PDstandardNth1K23 + PDstandardNth3K12) + gu33*
         (-(G133*K11L) - G233*K12L + G113*K13L - G333*K13L + G213*K23L + G313*K33L - PDstandardNth1K33 + PDstandardNth3K13);
    
    CCTK_REAL const M2L  =  gu11*(G112*K11L + (-G111 + G212)*K12L + G312*K13L - G211*K22L - G311*K23L + PDstandardNth1K12 - 
           PDstandardNth2K11) + gu12*(G122*K11L + (-G112 + G222)*K12L + G322*K13L - G212*K22L - G312*K23L + 
           PDstandardNth1K22 - PDstandardNth2K12) + gu13*
         (G123*K11L + (-G113 + G223)*K12L + G323*K13L - G213*K22L - G313*K23L + PDstandardNth1K23 - PDstandardNth2K13) + 
        gu31*(-(G113*K12L) + G112*K13L - G213*K22L + G212*K23L - G313*K23L + G312*K33L - PDstandardNth2K13 + 
           PDstandardNth3K12) + gu32*(-(G123*K12L) + G122*K13L - G223*K22L + G222*K23L - G323*K23L + G322*K33L - 
           PDstandardNth2K23 + PDstandardNth3K22) + gu33*
         (-(G133*K12L) + G123*K13L - G233*K22L + G223*K23L - G333*K23L + G323*K33L - PDstandardNth2K33 + PDstandardNth3K23);
    
    CCTK_REAL const M3L  =  gu11*(G113*K11L + G213*K12L + (-G111 + G313)*K13L - G211*K23L - G311*K33L + PDstandardNth1K13 - 
           PDstandardNth3K11) + gu12*(G123*K11L + G223*K12L + (-G112 + G323)*K13L - G212*K23L - G312*K33L + 
           PDstandardNth1K23 - PDstandardNth3K12) + gu21*
         (G113*K12L - G112*K13L + G213*K22L - G212*K23L + G313*K23L - G312*K33L + PDstandardNth2K13 - PDstandardNth3K12) + 
        gu13*(G133*K11L + G233*K12L + (-G113 + G333)*K13L - G213*K23L - G313*K33L + PDstandardNth1K33 - 
           PDstandardNth3K13) + gu22*(G123*K12L - G122*K13L + G223*K22L - G222*K23L + G323*K23L - G322*K33L + 
           PDstandardNth2K23 - PDstandardNth3K22) + gu23*
         (G133*K12L - G123*K13L + G233*K22L - G223*K23L + G333*K23L - G323*K33L + PDstandardNth2K33 - PDstandardNth3K23);
    
    
    /* Copy local copies back to grid functions */
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_ADM_constraints);
}

void ML_ADM_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADM_constraints_Body);
}
