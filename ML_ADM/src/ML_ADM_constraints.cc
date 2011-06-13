/*  File produced by Kranc */

#define KRANC_C

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"
#include "Differencing.h"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define QAD(x) (SQR(SQR(x)))
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))

extern "C" void ML_ADM_constraints_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_Ham","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_Ham.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_mom","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_mom.");
  return;
}

static void ML_ADM_constraints_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADM_constraints_Body");
  }
  
  if (cctk_iteration % ML_ADM_constraints_calc_every != ML_ADM_constraints_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_ADM::ML_curv","ML_ADM::ML_Ham","ML_ADM::ML_metric","ML_ADM::ML_mom"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADM_constraints", 4, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_ADM_constraints", 2, 2, 2);
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  ptrdiff_t const di = 1;
  ptrdiff_t const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const cdi = sizeof(CCTK_REAL) * di;
  ptrdiff_t const cdj = sizeof(CCTK_REAL) * dj;
  ptrdiff_t const cdk = sizeof(CCTK_REAL) * dk;
  CCTK_REAL const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL const dz = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL const dt = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL const dxi = INV(dx);
  CCTK_REAL const dyi = INV(dy);
  CCTK_REAL const dzi = INV(dz);
  CCTK_REAL const khalf = 0.5;
  CCTK_REAL const kthird = 1/3.0;
  CCTK_REAL const ktwothird = 2.0/3.0;
  CCTK_REAL const kfourthird = 4.0/3.0;
  CCTK_REAL const keightthird = 8.0/3.0;
  CCTK_REAL const hdxi = 0.5 * dxi;
  CCTK_REAL const hdyi = 0.5 * dyi;
  CCTK_REAL const hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  CCTK_REAL const p1o12dx = 0.0833333333333333333333333333333*INV(dx);
  CCTK_REAL const p1o12dy = 0.0833333333333333333333333333333*INV(dy);
  CCTK_REAL const p1o12dz = 0.0833333333333333333333333333333*INV(dz);
  CCTK_REAL const p1o144dxdy = 0.00694444444444444444444444444444*INV(dx)*INV(dy);
  CCTK_REAL const p1o144dxdz = 0.00694444444444444444444444444444*INV(dx)*INV(dz);
  CCTK_REAL const p1o144dydz = 0.00694444444444444444444444444444*INV(dy)*INV(dz);
  CCTK_REAL const pm1o12dx2 = -0.0833333333333333333333333333333*INV(SQR(dx));
  CCTK_REAL const pm1o12dy2 = -0.0833333333333333333333333333333*INV(SQR(dy));
  CCTK_REAL const pm1o12dz2 = -0.0833333333333333333333333333333*INV(SQR(dz));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_ADM_constraints,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL g11L = g11[index];
    CCTK_REAL g12L = g12[index];
    CCTK_REAL g13L = g13[index];
    CCTK_REAL g22L = g22[index];
    CCTK_REAL g23L = g23[index];
    CCTK_REAL g33L = g33[index];
    CCTK_REAL K11L = K11[index];
    CCTK_REAL K12L = K12[index];
    CCTK_REAL K13L = K13[index];
    CCTK_REAL K22L = K22[index];
    CCTK_REAL K23L = K23[index];
    CCTK_REAL K33L = K33[index];
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDstandardNth1g11 = PDstandardNth1(&g11[index]);
    CCTK_REAL const PDstandardNth2g11 = PDstandardNth2(&g11[index]);
    CCTK_REAL const PDstandardNth3g11 = PDstandardNth3(&g11[index]);
    CCTK_REAL const PDstandardNth22g11 = PDstandardNth22(&g11[index]);
    CCTK_REAL const PDstandardNth33g11 = PDstandardNth33(&g11[index]);
    CCTK_REAL const PDstandardNth23g11 = PDstandardNth23(&g11[index]);
    CCTK_REAL const PDstandardNth1g12 = PDstandardNth1(&g12[index]);
    CCTK_REAL const PDstandardNth2g12 = PDstandardNth2(&g12[index]);
    CCTK_REAL const PDstandardNth3g12 = PDstandardNth3(&g12[index]);
    CCTK_REAL const PDstandardNth33g12 = PDstandardNth33(&g12[index]);
    CCTK_REAL const PDstandardNth12g12 = PDstandardNth12(&g12[index]);
    CCTK_REAL const PDstandardNth13g12 = PDstandardNth13(&g12[index]);
    CCTK_REAL const PDstandardNth23g12 = PDstandardNth23(&g12[index]);
    CCTK_REAL const PDstandardNth1g13 = PDstandardNth1(&g13[index]);
    CCTK_REAL const PDstandardNth2g13 = PDstandardNth2(&g13[index]);
    CCTK_REAL const PDstandardNth3g13 = PDstandardNth3(&g13[index]);
    CCTK_REAL const PDstandardNth22g13 = PDstandardNth22(&g13[index]);
    CCTK_REAL const PDstandardNth12g13 = PDstandardNth12(&g13[index]);
    CCTK_REAL const PDstandardNth13g13 = PDstandardNth13(&g13[index]);
    CCTK_REAL const PDstandardNth23g13 = PDstandardNth23(&g13[index]);
    CCTK_REAL const PDstandardNth1g22 = PDstandardNth1(&g22[index]);
    CCTK_REAL const PDstandardNth2g22 = PDstandardNth2(&g22[index]);
    CCTK_REAL const PDstandardNth3g22 = PDstandardNth3(&g22[index]);
    CCTK_REAL const PDstandardNth11g22 = PDstandardNth11(&g22[index]);
    CCTK_REAL const PDstandardNth33g22 = PDstandardNth33(&g22[index]);
    CCTK_REAL const PDstandardNth13g22 = PDstandardNth13(&g22[index]);
    CCTK_REAL const PDstandardNth1g23 = PDstandardNth1(&g23[index]);
    CCTK_REAL const PDstandardNth2g23 = PDstandardNth2(&g23[index]);
    CCTK_REAL const PDstandardNth3g23 = PDstandardNth3(&g23[index]);
    CCTK_REAL const PDstandardNth11g23 = PDstandardNth11(&g23[index]);
    CCTK_REAL const PDstandardNth12g23 = PDstandardNth12(&g23[index]);
    CCTK_REAL const PDstandardNth13g23 = PDstandardNth13(&g23[index]);
    CCTK_REAL const PDstandardNth23g23 = PDstandardNth23(&g23[index]);
    CCTK_REAL const PDstandardNth1g33 = PDstandardNth1(&g33[index]);
    CCTK_REAL const PDstandardNth2g33 = PDstandardNth2(&g33[index]);
    CCTK_REAL const PDstandardNth3g33 = PDstandardNth3(&g33[index]);
    CCTK_REAL const PDstandardNth11g33 = PDstandardNth11(&g33[index]);
    CCTK_REAL const PDstandardNth22g33 = PDstandardNth22(&g33[index]);
    CCTK_REAL const PDstandardNth12g33 = PDstandardNth12(&g33[index]);
    CCTK_REAL const PDstandardNth2K11 = PDstandardNth2(&K11[index]);
    CCTK_REAL const PDstandardNth3K11 = PDstandardNth3(&K11[index]);
    CCTK_REAL const PDstandardNth1K12 = PDstandardNth1(&K12[index]);
    CCTK_REAL const PDstandardNth2K12 = PDstandardNth2(&K12[index]);
    CCTK_REAL const PDstandardNth3K12 = PDstandardNth3(&K12[index]);
    CCTK_REAL const PDstandardNth1K13 = PDstandardNth1(&K13[index]);
    CCTK_REAL const PDstandardNth2K13 = PDstandardNth2(&K13[index]);
    CCTK_REAL const PDstandardNth3K13 = PDstandardNth3(&K13[index]);
    CCTK_REAL const PDstandardNth1K22 = PDstandardNth1(&K22[index]);
    CCTK_REAL const PDstandardNth3K22 = PDstandardNth3(&K22[index]);
    CCTK_REAL const PDstandardNth1K23 = PDstandardNth1(&K23[index]);
    CCTK_REAL const PDstandardNth2K23 = PDstandardNth2(&K23[index]);
    CCTK_REAL const PDstandardNth3K23 = PDstandardNth3(&K23[index]);
    CCTK_REAL const PDstandardNth1K33 = PDstandardNth1(&K33[index]);
    CCTK_REAL const PDstandardNth2K33 = PDstandardNth2(&K33[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL detg = 2*g12L*g13L*g23L + g33L*(g11L*g22L - SQR(g12L)) - 
      g22L*SQR(g13L) - g11L*SQR(g23L);
    
    CCTK_REAL gu11 = INV(detg)*(g22L*g33L - SQR(g23L));
    
    CCTK_REAL gu12 = (g13L*g23L - g12L*g33L)*INV(detg);
    
    CCTK_REAL gu13 = (-(g13L*g22L) + g12L*g23L)*INV(detg);
    
    CCTK_REAL gu21 = (g13L*g23L - g12L*g33L)*INV(detg);
    
    CCTK_REAL gu22 = INV(detg)*(g11L*g33L - SQR(g13L));
    
    CCTK_REAL gu23 = (g12L*g13L - g11L*g23L)*INV(detg);
    
    CCTK_REAL gu31 = (-(g13L*g22L) + g12L*g23L)*INV(detg);
    
    CCTK_REAL gu32 = (g12L*g13L - g11L*g23L)*INV(detg);
    
    CCTK_REAL gu33 = INV(detg)*(g11L*g22L - SQR(g12L));
    
    CCTK_REAL G111 = 0.5*(gu11*PDstandardNth1g11 + 
      2*(gu12*PDstandardNth1g12 + gu13*PDstandardNth1g13) - 
      gu12*PDstandardNth2g11 - gu13*PDstandardNth3g11);
    
    CCTK_REAL G211 = 0.5*(gu21*PDstandardNth1g11 + 
      2*(gu22*PDstandardNth1g12 + gu23*PDstandardNth1g13) - 
      gu22*PDstandardNth2g11 - gu23*PDstandardNth3g11);
    
    CCTK_REAL G311 = 0.5*(gu31*PDstandardNth1g11 + 
      2*(gu32*PDstandardNth1g12 + gu33*PDstandardNth1g13) - 
      gu32*PDstandardNth2g11 - gu33*PDstandardNth3g11);
    
    CCTK_REAL G112 = 0.5*(gu12*PDstandardNth1g22 + gu11*PDstandardNth2g11 
      + gu13*(PDstandardNth1g23 + PDstandardNth2g13 - PDstandardNth3g12));
    
    CCTK_REAL G212 = 0.5*(gu22*PDstandardNth1g22 + gu21*PDstandardNth2g11 
      + gu23*(PDstandardNth1g23 + PDstandardNth2g13 - PDstandardNth3g12));
    
    CCTK_REAL G312 = 0.5*(gu32*PDstandardNth1g22 + gu31*PDstandardNth2g11 
      + gu33*(PDstandardNth1g23 + PDstandardNth2g13 - PDstandardNth3g12));
    
    CCTK_REAL G113 = 0.5*(gu13*PDstandardNth1g33 + gu11*PDstandardNth3g11 
      + gu12*(PDstandardNth1g23 - PDstandardNth2g13 + PDstandardNth3g12));
    
    CCTK_REAL G213 = 0.5*(gu23*PDstandardNth1g33 + gu21*PDstandardNth3g11 
      + gu22*(PDstandardNth1g23 - PDstandardNth2g13 + PDstandardNth3g12));
    
    CCTK_REAL G313 = 0.5*(gu33*PDstandardNth1g33 + gu31*PDstandardNth3g11 
      + gu32*(PDstandardNth1g23 - PDstandardNth2g13 + PDstandardNth3g12));
    
    CCTK_REAL G122 = 0.5*(gu11*(-PDstandardNth1g22 + 2*PDstandardNth2g12) 
      + gu12*PDstandardNth2g22 + gu13*(2*PDstandardNth2g23 - 
      PDstandardNth3g22));
    
    CCTK_REAL G222 = 0.5*(gu21*(-PDstandardNth1g22 + 2*PDstandardNth2g12) 
      + gu22*PDstandardNth2g22 + gu23*(2*PDstandardNth2g23 - 
      PDstandardNth3g22));
    
    CCTK_REAL G322 = 0.5*(gu31*(-PDstandardNth1g22 + 2*PDstandardNth2g12) 
      + gu32*PDstandardNth2g22 + gu33*(2*PDstandardNth2g23 - 
      PDstandardNth3g22));
    
    CCTK_REAL G123 = 0.5*(gu13*PDstandardNth2g33 + 
      gu11*(-PDstandardNth1g23 + PDstandardNth2g13 + PDstandardNth3g12) + 
      gu12*PDstandardNth3g22);
    
    CCTK_REAL G223 = 0.5*(gu23*PDstandardNth2g33 + 
      gu21*(-PDstandardNth1g23 + PDstandardNth2g13 + PDstandardNth3g12) + 
      gu22*PDstandardNth3g22);
    
    CCTK_REAL G323 = 0.5*(gu33*PDstandardNth2g33 + 
      gu31*(-PDstandardNth1g23 + PDstandardNth2g13 + PDstandardNth3g12) + 
      gu32*PDstandardNth3g22);
    
    CCTK_REAL G133 = 0.5*(gu11*(-PDstandardNth1g33 + 2*PDstandardNth3g13) 
      + gu12*(-PDstandardNth2g33 + 2*PDstandardNth3g23) + 
      gu13*PDstandardNth3g33);
    
    CCTK_REAL G233 = 0.5*(gu21*(-PDstandardNth1g33 + 2*PDstandardNth3g13) 
      + gu22*(-PDstandardNth2g33 + 2*PDstandardNth3g23) + 
      gu23*PDstandardNth3g33);
    
    CCTK_REAL G333 = 0.5*(gu31*(-PDstandardNth1g33 + 2*PDstandardNth3g13) 
      + gu32*(-PDstandardNth2g33 + 2*PDstandardNth3g23) + 
      gu33*PDstandardNth3g33);
    
    CCTK_REAL R11 = -(G111*(G111 + G122 + G133)) - G211*(G211 + G222 + 
      G233) - G311*(G311 + G322 + G333) - 0.5*gu22*(PDstandardNth11g22 - 
      2*PDstandardNth12g12 + PDstandardNth22g11) + 
      0.5*gu23*(-PDstandardNth11g23 + PDstandardNth12g13 + PDstandardNth13g12 
      - PDstandardNth23g11) + 0.5*gu32*(-PDstandardNth11g23 + 
      PDstandardNth12g13 + PDstandardNth13g12 - PDstandardNth23g11) - 
      0.5*gu33*(PDstandardNth11g33 - 2*PDstandardNth13g13 + 
      PDstandardNth33g11) + SQR(G111) + SQR(G112) + SQR(G113) + SQR(G211) + 
      SQR(G212) + SQR(G213) + SQR(G311) + SQR(G312) + SQR(G313);
    
    CCTK_REAL R12 = 0.5*(2*(G113*G123 + G213*G223 + G313*G323) - 
      2*(G112*G133 + G212*G233 + G312*G333 + gu12*PDstandardNth12g12) + 
      gu13*(PDstandardNth11g23 - PDstandardNth12g13 - PDstandardNth13g12) + 
      gu12*(PDstandardNth11g22 + PDstandardNth22g11) + 
      gu32*PDstandardNth22g13 + gu13*PDstandardNth23g11 + 
      gu32*(-PDstandardNth12g23 + PDstandardNth13g22 - PDstandardNth23g12) + 
      gu33*PDstandardNth23g13 + gu33*(-PDstandardNth12g33 + 
      PDstandardNth13g23 - PDstandardNth33g12));
    
    CCTK_REAL R13 = 0.5*(2*(G112*G123 + G212*G223 + G312*G323) - 
      2*(G113*G122 + G213*G222 + G313*G322 + gu13*PDstandardNth13g13) + 
      gu12*(PDstandardNth11g23 - PDstandardNth12g13 - PDstandardNth13g12 + 
      PDstandardNth23g11) + gu22*(PDstandardNth12g23 - PDstandardNth13g22 - 
      PDstandardNth22g13 + PDstandardNth23g12) + gu13*(PDstandardNth11g33 + 
      PDstandardNth33g11) + gu23*(PDstandardNth12g33 - PDstandardNth13g23 - 
      PDstandardNth23g13 + PDstandardNth33g12));
    
    CCTK_REAL R22 = -(G122*(G111 + G122 + G133)) - G222*(G211 + G222 + 
      G233) - G322*(G311 + G322 + G333) - 0.5*gu11*(PDstandardNth11g22 - 
      2*PDstandardNth12g12 + PDstandardNth22g11) + 
      0.5*gu13*(PDstandardNth12g23 - PDstandardNth13g22 - PDstandardNth22g13 
      + PDstandardNth23g12) + 0.5*gu31*(PDstandardNth12g23 - 
      PDstandardNth13g22 - PDstandardNth22g13 + PDstandardNth23g12) - 
      0.5*gu33*(PDstandardNth22g33 - 2*PDstandardNth23g23 + 
      PDstandardNth33g22) + SQR(G112) + SQR(G122) + SQR(G123) + SQR(G212) + 
      SQR(G222) + SQR(G223) + SQR(G312) + SQR(G322) + SQR(G323);
    
    CCTK_REAL R23 = 0.5*(2*(G112*G113 + G212*G213 + G312*G313) + 
      gu11*(-PDstandardNth11g23 + PDstandardNth12g13 + PDstandardNth13g12 - 
      PDstandardNth23g11) + gu21*(-PDstandardNth12g23 + PDstandardNth13g22 + 
      PDstandardNth22g13 - PDstandardNth23g12) - 2*(G111*G123 + G211*G223 + 
      G311*G323 + gu23*PDstandardNth23g23) + gu13*(PDstandardNth12g33 - 
      PDstandardNth13g23 - PDstandardNth23g13 + PDstandardNth33g12) + 
      gu23*(PDstandardNth22g33 + PDstandardNth33g22));
    
    CCTK_REAL R33 = -(G133*(G111 + G122 + G133)) - G233*(G211 + G222 + 
      G233) - G333*(G311 + G322 + G333) - 0.5*gu11*(PDstandardNth11g33 - 
      2*PDstandardNth13g13 + PDstandardNth33g11) + 
      0.5*gu12*(-PDstandardNth12g33 + PDstandardNth13g23 + PDstandardNth23g13 
      - PDstandardNth33g12) + 0.5*gu21*(-PDstandardNth12g33 + 
      PDstandardNth13g23 + PDstandardNth23g13 - PDstandardNth33g12) - 
      0.5*gu22*(PDstandardNth22g33 - 2*PDstandardNth23g23 + 
      PDstandardNth33g22) + SQR(G113) + SQR(G123) + SQR(G133) + SQR(G213) + 
      SQR(G223) + SQR(G233) + SQR(G313) + SQR(G323) + SQR(G333);
    
    CCTK_REAL trR = gu11*R11 + (gu12 + gu21)*R12 + (gu13 + gu31)*R13 + 
      gu22*R22 + (gu23 + gu32)*R23 + gu33*R33;
    
    CCTK_REAL Km11 = gu11*K11L + gu12*K12L + gu13*K13L;
    
    CCTK_REAL Km21 = gu21*K11L + gu22*K12L + gu23*K13L;
    
    CCTK_REAL Km31 = gu31*K11L + gu32*K12L + gu33*K13L;
    
    CCTK_REAL Km12 = gu11*K12L + gu12*K22L + gu13*K23L;
    
    CCTK_REAL Km22 = gu21*K12L + gu22*K22L + gu23*K23L;
    
    CCTK_REAL Km32 = gu31*K12L + gu32*K22L + gu33*K23L;
    
    CCTK_REAL Km13 = gu11*K13L + gu12*K23L + gu13*K33L;
    
    CCTK_REAL Km23 = gu21*K13L + gu22*K23L + gu23*K33L;
    
    CCTK_REAL Km33 = gu31*K13L + gu32*K23L + gu33*K33L;
    
    CCTK_REAL trK = Km11 + Km22 + Km33;
    
    CCTK_REAL HL = -2*(Km12*Km21 + Km13*Km31 + Km23*Km32) + trR - 
      SQR(Km11) - SQR(Km22) - SQR(Km33) + SQR(trK);
    
    CCTK_REAL M1L = gu21*(-(G112*K11L) + (G111 - G212)*K12L - G312*K13L + 
      G211*K22L + G311*K23L - PDstandardNth1K12 + PDstandardNth2K11) + 
      gu22*(-(G122*K11L) + (G112 - G222)*K12L - G322*K13L + G212*K22L + 
      G312*K23L - PDstandardNth1K22 + PDstandardNth2K12) + gu23*(-(G123*K11L) 
      + (G113 - G223)*K12L - G323*K13L + G213*K22L + G313*K23L - 
      PDstandardNth1K23 + PDstandardNth2K13) + gu31*(-(G113*K11L) - G213*K12L 
      + (G111 - G313)*K13L + G211*K23L + G311*K33L - PDstandardNth1K13 + 
      PDstandardNth3K11) + gu32*(-(G123*K11L) - G223*K12L + (G112 - 
      G323)*K13L + G212*K23L + G312*K33L - PDstandardNth1K23 + 
      PDstandardNth3K12) + gu33*(-(G133*K11L) - G233*K12L + (G113 - 
      G333)*K13L + G213*K23L + G313*K33L - PDstandardNth1K33 + 
      PDstandardNth3K13);
    
    CCTK_REAL M2L = gu11*(G112*K11L + (-G111 + G212)*K12L + G312*K13L - 
      G211*K22L - G311*K23L + PDstandardNth1K12 - PDstandardNth2K11) + 
      gu12*(G122*K11L + (-G112 + G222)*K12L + G322*K13L - G212*K22L - 
      G312*K23L + PDstandardNth1K22 - PDstandardNth2K12) + gu13*(G123*K11L + 
      (-G113 + G223)*K12L + G323*K13L - G213*K22L - G313*K23L + 
      PDstandardNth1K23 - PDstandardNth2K13) + gu31*(-(G113*K12L) + G112*K13L 
      - G213*K22L + (G212 - G313)*K23L + G312*K33L - PDstandardNth2K13 + 
      PDstandardNth3K12) + gu32*(-(G123*K12L) + G122*K13L - G223*K22L + (G222 
      - G323)*K23L + G322*K33L - PDstandardNth2K23 + PDstandardNth3K22) + 
      gu33*(-(G133*K12L) + G123*K13L - G233*K22L + (G223 - G333)*K23L + 
      G323*K33L - PDstandardNth2K33 + PDstandardNth3K23);
    
    CCTK_REAL M3L = gu11*(G113*K11L + G213*K12L + (-G111 + G313)*K13L - 
      G211*K23L - G311*K33L + PDstandardNth1K13 - PDstandardNth3K11) + 
      gu12*(G123*K11L + G223*K12L + (-G112 + G323)*K13L - G212*K23L - 
      G312*K33L + PDstandardNth1K23 - PDstandardNth3K12) + gu21*(G113*K12L - 
      G112*K13L + G213*K22L + (-G212 + G313)*K23L - G312*K33L + 
      PDstandardNth2K13 - PDstandardNth3K12) + gu13*(G133*K11L + G233*K12L + 
      (-G113 + G333)*K13L - G213*K23L - G313*K33L + PDstandardNth1K33 - 
      PDstandardNth3K13) + gu22*(G123*K12L - G122*K13L + G223*K22L + (-G222 + 
      G323)*K23L - G322*K33L + PDstandardNth2K23 - PDstandardNth3K22) + 
      gu23*(G133*K12L - G123*K13L + G233*K22L + (-G223 + G333)*K23L - 
      G323*K33L + PDstandardNth2K33 - PDstandardNth3K23);
    
    /* Copy local copies back to grid functions */
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
  }
  LC_ENDLOOP3 (ML_ADM_constraints);
}

extern "C" void ML_ADM_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADM_constraints_Body);
}
