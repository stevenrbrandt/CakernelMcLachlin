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

void ML_ADM_RHS_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_curvrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADM::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADM::ML_shiftrhs.");
  return;
}

void ML_ADM_RHS_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
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
  LC_LOOP3 (ML_ADM_RHS,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL PDstandardNth1alpha = INITVALUE;
    // CCTK_REAL PDstandardNth2alpha = INITVALUE;
    // CCTK_REAL PDstandardNth3alpha = INITVALUE;
    // CCTK_REAL PDstandardNth11alpha = INITVALUE;
    // CCTK_REAL PDstandardNth22alpha = INITVALUE;
    // CCTK_REAL PDstandardNth33alpha = INITVALUE;
    // CCTK_REAL PDstandardNth12alpha = INITVALUE;
    // CCTK_REAL PDstandardNth13alpha = INITVALUE;
    // CCTK_REAL PDstandardNth23alpha = INITVALUE;
    // CCTK_REAL PDstandardNth1beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth2beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth3beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth1beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth2beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth3beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth1beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth2beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth3beta3 = INITVALUE;
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
    // CCTK_REAL PDstandardNth1K11 = INITVALUE;
    // CCTK_REAL PDstandardNth2K11 = INITVALUE;
    // CCTK_REAL PDstandardNth3K11 = INITVALUE;
    // CCTK_REAL PDstandardNth1K12 = INITVALUE;
    // CCTK_REAL PDstandardNth2K12 = INITVALUE;
    // CCTK_REAL PDstandardNth3K12 = INITVALUE;
    // CCTK_REAL PDstandardNth1K13 = INITVALUE;
    // CCTK_REAL PDstandardNth2K13 = INITVALUE;
    // CCTK_REAL PDstandardNth3K13 = INITVALUE;
    // CCTK_REAL PDstandardNth1K22 = INITVALUE;
    // CCTK_REAL PDstandardNth2K22 = INITVALUE;
    // CCTK_REAL PDstandardNth3K22 = INITVALUE;
    // CCTK_REAL PDstandardNth1K23 = INITVALUE;
    // CCTK_REAL PDstandardNth2K23 = INITVALUE;
    // CCTK_REAL PDstandardNth3K23 = INITVALUE;
    // CCTK_REAL PDstandardNth1K33 = INITVALUE;
    // CCTK_REAL PDstandardNth2K33 = INITVALUE;
    // CCTK_REAL PDstandardNth3K33 = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL  alphaL = alpha[index];
    CCTK_REAL  beta1L = beta1[index];
    CCTK_REAL  beta2L = beta2[index];
    CCTK_REAL  beta3L = beta3[index];
    CCTK_REAL  g11L = g11[index];
    CCTK_REAL  g12L = g12[index];
    CCTK_REAL  g13L = g13[index];
    CCTK_REAL  g22L = g22[index];
    CCTK_REAL  g23L = g23[index];
    CCTK_REAL  g33L = g33[index];
    CCTK_REAL  K11L = K11[index];
    CCTK_REAL  K12L = K12[index];
    CCTK_REAL  K13L = K13[index];
    CCTK_REAL  K22L = K22[index];
    CCTK_REAL  K23L = K23[index];
    CCTK_REAL  K33L = K33[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDstandardNth1alpha = PDstandardNth1(alpha, i, j, k);
    CCTK_REAL const PDstandardNth2alpha = PDstandardNth2(alpha, i, j, k);
    CCTK_REAL const PDstandardNth3alpha = PDstandardNth3(alpha, i, j, k);
    CCTK_REAL const PDstandardNth11alpha = PDstandardNth11(alpha, i, j, k);
    CCTK_REAL const PDstandardNth22alpha = PDstandardNth22(alpha, i, j, k);
    CCTK_REAL const PDstandardNth33alpha = PDstandardNth33(alpha, i, j, k);
    CCTK_REAL const PDstandardNth12alpha = PDstandardNth12(alpha, i, j, k);
    CCTK_REAL const PDstandardNth13alpha = PDstandardNth13(alpha, i, j, k);
    CCTK_REAL const PDstandardNth23alpha = PDstandardNth23(alpha, i, j, k);
    CCTK_REAL const PDstandardNth1beta1 = PDstandardNth1(beta1, i, j, k);
    CCTK_REAL const PDstandardNth2beta1 = PDstandardNth2(beta1, i, j, k);
    CCTK_REAL const PDstandardNth3beta1 = PDstandardNth3(beta1, i, j, k);
    CCTK_REAL const PDstandardNth1beta2 = PDstandardNth1(beta2, i, j, k);
    CCTK_REAL const PDstandardNth2beta2 = PDstandardNth2(beta2, i, j, k);
    CCTK_REAL const PDstandardNth3beta2 = PDstandardNth3(beta2, i, j, k);
    CCTK_REAL const PDstandardNth1beta3 = PDstandardNth1(beta3, i, j, k);
    CCTK_REAL const PDstandardNth2beta3 = PDstandardNth2(beta3, i, j, k);
    CCTK_REAL const PDstandardNth3beta3 = PDstandardNth3(beta3, i, j, k);
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
    CCTK_REAL const PDstandardNth1K11 = PDstandardNth1(K11, i, j, k);
    CCTK_REAL const PDstandardNth2K11 = PDstandardNth2(K11, i, j, k);
    CCTK_REAL const PDstandardNth3K11 = PDstandardNth3(K11, i, j, k);
    CCTK_REAL const PDstandardNth1K12 = PDstandardNth1(K12, i, j, k);
    CCTK_REAL const PDstandardNth2K12 = PDstandardNth2(K12, i, j, k);
    CCTK_REAL const PDstandardNth3K12 = PDstandardNth3(K12, i, j, k);
    CCTK_REAL const PDstandardNth1K13 = PDstandardNth1(K13, i, j, k);
    CCTK_REAL const PDstandardNth2K13 = PDstandardNth2(K13, i, j, k);
    CCTK_REAL const PDstandardNth3K13 = PDstandardNth3(K13, i, j, k);
    CCTK_REAL const PDstandardNth1K22 = PDstandardNth1(K22, i, j, k);
    CCTK_REAL const PDstandardNth2K22 = PDstandardNth2(K22, i, j, k);
    CCTK_REAL const PDstandardNth3K22 = PDstandardNth3(K22, i, j, k);
    CCTK_REAL const PDstandardNth1K23 = PDstandardNth1(K23, i, j, k);
    CCTK_REAL const PDstandardNth2K23 = PDstandardNth2(K23, i, j, k);
    CCTK_REAL const PDstandardNth3K23 = PDstandardNth3(K23, i, j, k);
    CCTK_REAL const PDstandardNth1K33 = PDstandardNth1(K33, i, j, k);
    CCTK_REAL const PDstandardNth2K33 = PDstandardNth2(K33, i, j, k);
    CCTK_REAL const PDstandardNth3K33 = PDstandardNth3(K33, i, j, k);
    
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
    
    CCTK_REAL G111 = khalf*(gu11*PDstandardNth1g11 + 
      2*(gu12*PDstandardNth1g12 + gu13*PDstandardNth1g13) - 
      gu12*PDstandardNth2g11 - gu13*PDstandardNth3g11);
    
    CCTK_REAL G211 = khalf*(gu21*PDstandardNth1g11 + 
      2*(gu22*PDstandardNth1g12 + gu23*PDstandardNth1g13) - 
      gu22*PDstandardNth2g11 - gu23*PDstandardNth3g11);
    
    CCTK_REAL G311 = khalf*(gu31*PDstandardNth1g11 + 
      2*(gu32*PDstandardNth1g12 + gu33*PDstandardNth1g13) - 
      gu32*PDstandardNth2g11 - gu33*PDstandardNth3g11);
    
    CCTK_REAL G112 = khalf*(gu12*PDstandardNth1g22 + 
      gu11*PDstandardNth2g11 + gu13*(PDstandardNth1g23 + PDstandardNth2g13 - 
      PDstandardNth3g12));
    
    CCTK_REAL G212 = khalf*(gu22*PDstandardNth1g22 + 
      gu21*PDstandardNth2g11 + gu23*(PDstandardNth1g23 + PDstandardNth2g13 - 
      PDstandardNth3g12));
    
    CCTK_REAL G312 = khalf*(gu32*PDstandardNth1g22 + 
      gu31*PDstandardNth2g11 + gu33*(PDstandardNth1g23 + PDstandardNth2g13 - 
      PDstandardNth3g12));
    
    CCTK_REAL G113 = khalf*(gu13*PDstandardNth1g33 + 
      gu11*PDstandardNth3g11 + gu12*(PDstandardNth1g23 - PDstandardNth2g13 + 
      PDstandardNth3g12));
    
    CCTK_REAL G213 = khalf*(gu23*PDstandardNth1g33 + 
      gu21*PDstandardNth3g11 + gu22*(PDstandardNth1g23 - PDstandardNth2g13 + 
      PDstandardNth3g12));
    
    CCTK_REAL G313 = khalf*(gu33*PDstandardNth1g33 + 
      gu31*PDstandardNth3g11 + gu32*(PDstandardNth1g23 - PDstandardNth2g13 + 
      PDstandardNth3g12));
    
    CCTK_REAL G122 = khalf*(gu11*(-PDstandardNth1g22 + 
      2*PDstandardNth2g12) + gu12*PDstandardNth2g22 + 
      gu13*(2*PDstandardNth2g23 - PDstandardNth3g22));
    
    CCTK_REAL G222 = khalf*(gu21*(-PDstandardNth1g22 + 
      2*PDstandardNth2g12) + gu22*PDstandardNth2g22 + 
      gu23*(2*PDstandardNth2g23 - PDstandardNth3g22));
    
    CCTK_REAL G322 = khalf*(gu31*(-PDstandardNth1g22 + 
      2*PDstandardNth2g12) + gu32*PDstandardNth2g22 + 
      gu33*(2*PDstandardNth2g23 - PDstandardNth3g22));
    
    CCTK_REAL G123 = khalf*(gu13*PDstandardNth2g33 + 
      gu11*(-PDstandardNth1g23 + PDstandardNth2g13 + PDstandardNth3g12) + 
      gu12*PDstandardNth3g22);
    
    CCTK_REAL G223 = khalf*(gu23*PDstandardNth2g33 + 
      gu21*(-PDstandardNth1g23 + PDstandardNth2g13 + PDstandardNth3g12) + 
      gu22*PDstandardNth3g22);
    
    CCTK_REAL G323 = khalf*(gu33*PDstandardNth2g33 + 
      gu31*(-PDstandardNth1g23 + PDstandardNth2g13 + PDstandardNth3g12) + 
      gu32*PDstandardNth3g22);
    
    CCTK_REAL G133 = khalf*(-(gu11*PDstandardNth1g33) - 
      gu12*PDstandardNth2g33 + 2*gu11*PDstandardNth3g13 + 
      2*gu12*PDstandardNth3g23 + gu13*PDstandardNth3g33);
    
    CCTK_REAL G233 = khalf*(-(gu21*PDstandardNth1g33) - 
      gu22*PDstandardNth2g33 + 2*gu21*PDstandardNth3g13 + 
      2*gu22*PDstandardNth3g23 + gu23*PDstandardNth3g33);
    
    CCTK_REAL G333 = khalf*(-(gu31*PDstandardNth1g33) - 
      gu32*PDstandardNth2g33 + 2*gu31*PDstandardNth3g13 + 
      2*gu32*PDstandardNth3g23 + gu33*PDstandardNth3g33);
    
    CCTK_REAL R11 = 2*(G112*G211 + G113*G311 + G213*G312) - G111*(G111 + 
      G212 + G313) - G211*(G112 + G222 + G323) - G311*(G113 + G223 + G333) + 
      khalf*(-(gu22*(PDstandardNth11g22 - 2*PDstandardNth12g12 + 
      PDstandardNth22g11)) + gu23*(-PDstandardNth11g23 + PDstandardNth12g13 + 
      PDstandardNth13g12 - PDstandardNth23g11) + gu32*(-PDstandardNth11g23 + 
      PDstandardNth12g13 + PDstandardNth13g12 - PDstandardNth23g11) - 
      gu33*(PDstandardNth11g33 - 2*PDstandardNth13g13 + PDstandardNth33g11)) 
      + SQR(G111) + SQR(G212) + SQR(G313);
    
    CCTK_REAL R12 = khalf*(2*(G122*G211 + G123*G311 + G213*G322 + 
      G313*G323) - 2*(G112*G212 + G112*G313 + G212*G323 + G312*G333 + 
      gu12*PDstandardNth12g12) - gu32*PDstandardNth12g23 - 
      gu33*PDstandardNth12g33 + gu13*(PDstandardNth11g23 - PDstandardNth12g13 
      - PDstandardNth13g12) + gu32*PDstandardNth13g22 + 
      gu33*PDstandardNth13g23 + gu12*(PDstandardNth11g22 + 
      PDstandardNth22g11) + gu32*PDstandardNth22g13 + gu13*PDstandardNth23g11 
      - gu32*PDstandardNth23g12 + gu33*PDstandardNth23g13 - 
      gu33*PDstandardNth33g12);
    
    CCTK_REAL R13 = khalf*(2*(G123*G211 + G212*G223 + G133*G311 + 
      G233*G312) - 2*(G213*G222 + G223*G313 + G113*(G212 + G313) + 
      gu13*PDstandardNth13g13) + gu12*(PDstandardNth11g23 - 
      PDstandardNth12g13 - PDstandardNth13g12 + PDstandardNth23g11) + 
      gu22*(PDstandardNth12g23 - PDstandardNth13g22 - PDstandardNth22g13 + 
      PDstandardNth23g12) + gu13*(PDstandardNth11g33 + PDstandardNth33g11) + 
      gu23*(PDstandardNth12g33 - PDstandardNth13g23 - PDstandardNth23g13 + 
      PDstandardNth33g12));
    
    CCTK_REAL R22 = -(G122*(G111 + G212 + G313)) + 2*(G122*G212 + 
      G123*G312 + G223*G322) - G222*(G112 + G222 + G323) - G322*(G113 + G223 
      + G333) + khalf*(-(gu11*(PDstandardNth11g22 - 2*PDstandardNth12g12 + 
      PDstandardNth22g11)) + gu13*(PDstandardNth12g23 - PDstandardNth13g22 - 
      PDstandardNth22g13 + PDstandardNth23g12) + gu31*(PDstandardNth12g23 - 
      PDstandardNth13g22 - PDstandardNth22g13 + PDstandardNth23g12) - 
      gu33*(PDstandardNth22g33 - 2*PDstandardNth23g23 + PDstandardNth33g22)) 
      + SQR(G112) + SQR(G222) + SQR(G323);
    
    CCTK_REAL R23 = khalf*(2*(G112*G113 + G122*G213 + G133*G312 + 
      G233*G322) + gu11*(-PDstandardNth11g23 + PDstandardNth12g13 + 
      PDstandardNth13g12 - PDstandardNth23g11) + gu21*(-PDstandardNth12g23 + 
      PDstandardNth13g22 + PDstandardNth22g13 - PDstandardNth23g12) - 
      2*(G111*G123 + G113*G323 + G223*(G112 + G323) + 
      gu23*PDstandardNth23g23) + gu13*(PDstandardNth12g33 - 
      PDstandardNth13g23 - PDstandardNth23g13 + PDstandardNth33g12) + 
      gu23*(PDstandardNth22g33 + PDstandardNth33g22));
    
    CCTK_REAL R33 = -(G133*(G111 + G212 + G313)) + 2*(G123*G213 + 
      G133*G313) + 2*G233*G323 - G233*(G112 + G222 + G323) - G333*(G113 + 
      G223 + G333) + khalf*(-(gu11*(PDstandardNth11g33 - 2*PDstandardNth13g13 
      + PDstandardNth33g11)) + gu12*(-PDstandardNth12g33 + PDstandardNth13g23 
      + PDstandardNth23g13 - PDstandardNth33g12) + gu21*(-PDstandardNth12g33 
      + PDstandardNth13g23 + PDstandardNth23g13 - PDstandardNth33g12) - 
      gu22*(PDstandardNth22g33 - 2*PDstandardNth23g23 + PDstandardNth33g22)) 
      + SQR(G113) + SQR(G223) + SQR(G333);
    
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
    
    CCTK_REAL g11rhsL = -2*alphaL*K11L + 2*(g11L*PDstandardNth1beta1 + 
      g12L*PDstandardNth1beta2 + g13L*PDstandardNth1beta3) + 
      beta1L*PDstandardNth1g11 + beta2L*PDstandardNth2g11 + 
      beta3L*PDstandardNth3g11;
    
    CCTK_REAL g12rhsL = -2*alphaL*K12L + g22L*PDstandardNth1beta2 + 
      g23L*PDstandardNth1beta3 + beta1L*PDstandardNth1g12 + 
      g11L*PDstandardNth2beta1 + g12L*(PDstandardNth1beta1 + 
      PDstandardNth2beta2) + g13L*PDstandardNth2beta3 + 
      beta2L*PDstandardNth2g12 + beta3L*PDstandardNth3g12;
    
    CCTK_REAL g13rhsL = -2*alphaL*K13L + g23L*PDstandardNth1beta2 + 
      g33L*PDstandardNth1beta3 + beta1L*PDstandardNth1g13 + 
      beta2L*PDstandardNth2g13 + g11L*PDstandardNth3beta1 + 
      g12L*PDstandardNth3beta2 + g13L*(PDstandardNth1beta1 + 
      PDstandardNth3beta3) + beta3L*PDstandardNth3g13;
    
    CCTK_REAL g22rhsL = -2*alphaL*K22L + beta1L*PDstandardNth1g22 + 
      2*(g12L*PDstandardNth2beta1 + g22L*PDstandardNth2beta2 + 
      g23L*PDstandardNth2beta3) + beta2L*PDstandardNth2g22 + 
      beta3L*PDstandardNth3g22;
    
    CCTK_REAL g23rhsL = -2*alphaL*K23L + beta1L*PDstandardNth1g23 + 
      g13L*PDstandardNth2beta1 + g33L*PDstandardNth2beta3 + 
      beta2L*PDstandardNth2g23 + g12L*PDstandardNth3beta1 + 
      g22L*PDstandardNth3beta2 + g23L*(PDstandardNth2beta2 + 
      PDstandardNth3beta3) + beta3L*PDstandardNth3g23;
    
    CCTK_REAL g33rhsL = -2*alphaL*K33L + beta1L*PDstandardNth1g33 + 
      beta2L*PDstandardNth2g33 + 2*(g13L*PDstandardNth3beta1 + 
      g23L*PDstandardNth3beta2 + g33L*PDstandardNth3beta3) + 
      beta3L*PDstandardNth3g33;
    
    CCTK_REAL K11rhsL = -PDstandardNth11alpha + G111*PDstandardNth1alpha + 
      2*(K11L*PDstandardNth1beta1 + K12L*PDstandardNth1beta2 + 
      K13L*PDstandardNth1beta3) + beta1L*PDstandardNth1K11 + 
      G211*PDstandardNth2alpha + beta2L*PDstandardNth2K11 + 
      G311*PDstandardNth3alpha + beta3L*PDstandardNth3K11 + 
      alphaL*(-2*(K11L*Km11 + K12L*Km21 + K13L*Km31) + R11 + K11L*trK);
    
    CCTK_REAL K12rhsL = -PDstandardNth12alpha + G112*PDstandardNth1alpha + 
      K22L*PDstandardNth1beta2 + K23L*PDstandardNth1beta3 + 
      beta1L*PDstandardNth1K12 + G212*PDstandardNth2alpha + 
      K11L*PDstandardNth2beta1 + K12L*(PDstandardNth1beta1 + 
      PDstandardNth2beta2) + K13L*PDstandardNth2beta3 + 
      beta2L*PDstandardNth2K12 + G312*PDstandardNth3alpha + 
      beta3L*PDstandardNth3K12 + alphaL*(-2*(K11L*Km12 + K12L*Km22 + 
      K13L*Km32) + R12 + K12L*trK);
    
    CCTK_REAL K13rhsL = -PDstandardNth13alpha + G113*PDstandardNth1alpha + 
      K23L*PDstandardNth1beta2 + K33L*PDstandardNth1beta3 + 
      beta1L*PDstandardNth1K13 + G213*PDstandardNth2alpha + 
      beta2L*PDstandardNth2K13 + G313*PDstandardNth3alpha + 
      K11L*PDstandardNth3beta1 + K12L*PDstandardNth3beta2 + 
      K13L*(PDstandardNth1beta1 + PDstandardNth3beta3) + 
      beta3L*PDstandardNth3K13 + alphaL*(-2*(K11L*Km13 + K12L*Km23 + 
      K13L*Km33) + R13 + K13L*trK);
    
    CCTK_REAL K22rhsL = G122*PDstandardNth1alpha + 
      beta1L*PDstandardNth1K22 - PDstandardNth22alpha + 
      G222*PDstandardNth2alpha + 2*(K12L*PDstandardNth2beta1 + 
      K22L*PDstandardNth2beta2 + K23L*PDstandardNth2beta3) + 
      beta2L*PDstandardNth2K22 + G322*PDstandardNth3alpha + 
      beta3L*PDstandardNth3K22 + alphaL*(-2*(K12L*Km12 + K22L*Km22 + 
      K23L*Km32) + R22 + K22L*trK);
    
    CCTK_REAL K23rhsL = G123*PDstandardNth1alpha + 
      beta1L*PDstandardNth1K23 - PDstandardNth23alpha + 
      G223*PDstandardNth2alpha + K13L*PDstandardNth2beta1 + 
      K33L*PDstandardNth2beta3 + beta2L*PDstandardNth2K23 + 
      G323*PDstandardNth3alpha + K12L*PDstandardNth3beta1 + 
      K22L*PDstandardNth3beta2 + K23L*(PDstandardNth2beta2 + 
      PDstandardNth3beta3) + beta3L*PDstandardNth3K23 + alphaL*(-2*(K12L*Km13 
      + K22L*Km23 + K23L*Km33) + R23 + K23L*trK);
    
    CCTK_REAL K33rhsL = G133*PDstandardNth1alpha + 
      beta1L*PDstandardNth1K33 + G233*PDstandardNth2alpha + 
      beta2L*PDstandardNth2K33 - PDstandardNth33alpha + 
      G333*PDstandardNth3alpha + 2*(K13L*PDstandardNth3beta1 + 
      K23L*PDstandardNth3beta2 + K33L*PDstandardNth3beta3) + 
      beta3L*PDstandardNth3K33 + alphaL*(-2*(K13L*Km13 + K23L*Km23 + 
      K33L*Km33) + R33 + K33L*trK);
    
    CCTK_REAL alpharhsL = 0;
    
    CCTK_REAL beta1rhsL = 0;
    
    CCTK_REAL beta2rhsL = 0;
    
    CCTK_REAL beta3rhsL = 0;
    
    
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
  }
  LC_ENDLOOP3 (ML_ADM_RHS);
}

void ML_ADM_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADM_RHS_Body);
}
