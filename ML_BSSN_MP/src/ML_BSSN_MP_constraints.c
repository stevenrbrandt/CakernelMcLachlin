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

void ML_BSSN_MP_constraints_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
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
  CCTK_REAL const p1o24dx = INV(dx)/24.;
  CCTK_REAL const p1o24dy = INV(dy)/24.;
  CCTK_REAL const p1o24dz = INV(dz)/24.;
  CCTK_REAL const p1o64dx = INV(dx)/64.;
  CCTK_REAL const p1o64dy = INV(dy)/64.;
  CCTK_REAL const p1o64dz = INV(dz)/64.;
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o12dx2 = -pow(dx,-2)/12.;
  CCTK_REAL const pm1o12dy2 = -pow(dy,-2)/12.;
  CCTK_REAL const pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_MP_constraints,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL PDstandardNth1At11 = INITVALUE;
    // CCTK_REAL PDstandardNth2At11 = INITVALUE;
    // CCTK_REAL PDstandardNth3At11 = INITVALUE;
    // CCTK_REAL PDstandardNth1At12 = INITVALUE;
    // CCTK_REAL PDstandardNth2At12 = INITVALUE;
    // CCTK_REAL PDstandardNth3At12 = INITVALUE;
    // CCTK_REAL PDstandardNth1At13 = INITVALUE;
    // CCTK_REAL PDstandardNth2At13 = INITVALUE;
    // CCTK_REAL PDstandardNth3At13 = INITVALUE;
    // CCTK_REAL PDstandardNth1At22 = INITVALUE;
    // CCTK_REAL PDstandardNth2At22 = INITVALUE;
    // CCTK_REAL PDstandardNth3At22 = INITVALUE;
    // CCTK_REAL PDstandardNth1At23 = INITVALUE;
    // CCTK_REAL PDstandardNth2At23 = INITVALUE;
    // CCTK_REAL PDstandardNth3At23 = INITVALUE;
    // CCTK_REAL PDstandardNth1At33 = INITVALUE;
    // CCTK_REAL PDstandardNth2At33 = INITVALUE;
    // CCTK_REAL PDstandardNth3At33 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth11gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth22gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth33gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth12gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth13gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth23gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth1phi = INITVALUE;
    // CCTK_REAL PDstandardNth2phi = INITVALUE;
    // CCTK_REAL PDstandardNth3phi = INITVALUE;
    // CCTK_REAL PDstandardNth11phi = INITVALUE;
    // CCTK_REAL PDstandardNth22phi = INITVALUE;
    // CCTK_REAL PDstandardNth33phi = INITVALUE;
    // CCTK_REAL PDstandardNth12phi = INITVALUE;
    // CCTK_REAL PDstandardNth13phi = INITVALUE;
    // CCTK_REAL PDstandardNth23phi = INITVALUE;
    // CCTK_REAL PDstandardNth1trK = INITVALUE;
    // CCTK_REAL PDstandardNth2trK = INITVALUE;
    // CCTK_REAL PDstandardNth3trK = INITVALUE;
    // CCTK_REAL PDstandardNth1Xt1 = INITVALUE;
    // CCTK_REAL PDstandardNth2Xt1 = INITVALUE;
    // CCTK_REAL PDstandardNth3Xt1 = INITVALUE;
    // CCTK_REAL PDstandardNth1Xt2 = INITVALUE;
    // CCTK_REAL PDstandardNth2Xt2 = INITVALUE;
    // CCTK_REAL PDstandardNth3Xt2 = INITVALUE;
    // CCTK_REAL PDstandardNth1Xt3 = INITVALUE;
    // CCTK_REAL PDstandardNth2Xt3 = INITVALUE;
    // CCTK_REAL PDstandardNth3Xt3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL  alphaL = alpha[index];
    CCTK_REAL  At11L = At11[index];
    CCTK_REAL  At12L = At12[index];
    CCTK_REAL  At13L = At13[index];
    CCTK_REAL  At22L = At22[index];
    CCTK_REAL  At23L = At23[index];
    CCTK_REAL  At33L = At33[index];
    CCTK_REAL  beta1L = beta1[index];
    CCTK_REAL  beta2L = beta2[index];
    CCTK_REAL  beta3L = beta3[index];
    CCTK_REAL  dJ111L = dJ111[index];
    CCTK_REAL  dJ112L = dJ112[index];
    CCTK_REAL  dJ113L = dJ113[index];
    CCTK_REAL  dJ122L = dJ122[index];
    CCTK_REAL  dJ123L = dJ123[index];
    CCTK_REAL  dJ133L = dJ133[index];
    CCTK_REAL  dJ211L = dJ211[index];
    CCTK_REAL  dJ212L = dJ212[index];
    CCTK_REAL  dJ213L = dJ213[index];
    CCTK_REAL  dJ222L = dJ222[index];
    CCTK_REAL  dJ223L = dJ223[index];
    CCTK_REAL  dJ233L = dJ233[index];
    CCTK_REAL  dJ311L = dJ311[index];
    CCTK_REAL  dJ312L = dJ312[index];
    CCTK_REAL  dJ313L = dJ313[index];
    CCTK_REAL  dJ322L = dJ322[index];
    CCTK_REAL  dJ323L = dJ323[index];
    CCTK_REAL  dJ333L = dJ333[index];
    CCTK_REAL  eTttL = (*stress_energy_state) ? (eTtt[index]) : 0.0;
    CCTK_REAL  eTtxL = (*stress_energy_state) ? (eTtx[index]) : 0.0;
    CCTK_REAL  eTtyL = (*stress_energy_state) ? (eTty[index]) : 0.0;
    CCTK_REAL  eTtzL = (*stress_energy_state) ? (eTtz[index]) : 0.0;
    CCTK_REAL  eTxxL = (*stress_energy_state) ? (eTxx[index]) : 0.0;
    CCTK_REAL  eTxyL = (*stress_energy_state) ? (eTxy[index]) : 0.0;
    CCTK_REAL  eTxzL = (*stress_energy_state) ? (eTxz[index]) : 0.0;
    CCTK_REAL  eTyyL = (*stress_energy_state) ? (eTyy[index]) : 0.0;
    CCTK_REAL  eTyzL = (*stress_energy_state) ? (eTyz[index]) : 0.0;
    CCTK_REAL  eTzzL = (*stress_energy_state) ? (eTzz[index]) : 0.0;
    CCTK_REAL  gt11L = gt11[index];
    CCTK_REAL  gt12L = gt12[index];
    CCTK_REAL  gt13L = gt13[index];
    CCTK_REAL  gt22L = gt22[index];
    CCTK_REAL  gt23L = gt23[index];
    CCTK_REAL  gt33L = gt33[index];
    CCTK_REAL  J11L = J11[index];
    CCTK_REAL  J12L = J12[index];
    CCTK_REAL  J13L = J13[index];
    CCTK_REAL  J21L = J21[index];
    CCTK_REAL  J22L = J22[index];
    CCTK_REAL  J23L = J23[index];
    CCTK_REAL  J31L = J31[index];
    CCTK_REAL  J32L = J32[index];
    CCTK_REAL  J33L = J33[index];
    CCTK_REAL  phiL = phi[index];
    CCTK_REAL  trKL = trK[index];
    CCTK_REAL  Xt1L = Xt1[index];
    CCTK_REAL  Xt2L = Xt2[index];
    CCTK_REAL  Xt3L = Xt3[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDstandardNth1At11 = PDstandardNth1(At11, i, j, k);
    CCTK_REAL const PDstandardNth2At11 = PDstandardNth2(At11, i, j, k);
    CCTK_REAL const PDstandardNth3At11 = PDstandardNth3(At11, i, j, k);
    CCTK_REAL const PDstandardNth1At12 = PDstandardNth1(At12, i, j, k);
    CCTK_REAL const PDstandardNth2At12 = PDstandardNth2(At12, i, j, k);
    CCTK_REAL const PDstandardNth3At12 = PDstandardNth3(At12, i, j, k);
    CCTK_REAL const PDstandardNth1At13 = PDstandardNth1(At13, i, j, k);
    CCTK_REAL const PDstandardNth2At13 = PDstandardNth2(At13, i, j, k);
    CCTK_REAL const PDstandardNth3At13 = PDstandardNth3(At13, i, j, k);
    CCTK_REAL const PDstandardNth1At22 = PDstandardNth1(At22, i, j, k);
    CCTK_REAL const PDstandardNth2At22 = PDstandardNth2(At22, i, j, k);
    CCTK_REAL const PDstandardNth3At22 = PDstandardNth3(At22, i, j, k);
    CCTK_REAL const PDstandardNth1At23 = PDstandardNth1(At23, i, j, k);
    CCTK_REAL const PDstandardNth2At23 = PDstandardNth2(At23, i, j, k);
    CCTK_REAL const PDstandardNth3At23 = PDstandardNth3(At23, i, j, k);
    CCTK_REAL const PDstandardNth1At33 = PDstandardNth1(At33, i, j, k);
    CCTK_REAL const PDstandardNth2At33 = PDstandardNth2(At33, i, j, k);
    CCTK_REAL const PDstandardNth3At33 = PDstandardNth3(At33, i, j, k);
    CCTK_REAL const PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    CCTK_REAL const PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    CCTK_REAL const PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    CCTK_REAL const PDstandardNth11gt11 = PDstandardNth11(gt11, i, j, k);
    CCTK_REAL const PDstandardNth22gt11 = PDstandardNth22(gt11, i, j, k);
    CCTK_REAL const PDstandardNth33gt11 = PDstandardNth33(gt11, i, j, k);
    CCTK_REAL const PDstandardNth12gt11 = PDstandardNth12(gt11, i, j, k);
    CCTK_REAL const PDstandardNth13gt11 = PDstandardNth13(gt11, i, j, k);
    CCTK_REAL const PDstandardNth23gt11 = PDstandardNth23(gt11, i, j, k);
    CCTK_REAL const PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    CCTK_REAL const PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    CCTK_REAL const PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    CCTK_REAL const PDstandardNth11gt12 = PDstandardNth11(gt12, i, j, k);
    CCTK_REAL const PDstandardNth22gt12 = PDstandardNth22(gt12, i, j, k);
    CCTK_REAL const PDstandardNth33gt12 = PDstandardNth33(gt12, i, j, k);
    CCTK_REAL const PDstandardNth12gt12 = PDstandardNth12(gt12, i, j, k);
    CCTK_REAL const PDstandardNth13gt12 = PDstandardNth13(gt12, i, j, k);
    CCTK_REAL const PDstandardNth23gt12 = PDstandardNth23(gt12, i, j, k);
    CCTK_REAL const PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    CCTK_REAL const PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    CCTK_REAL const PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    CCTK_REAL const PDstandardNth11gt13 = PDstandardNth11(gt13, i, j, k);
    CCTK_REAL const PDstandardNth22gt13 = PDstandardNth22(gt13, i, j, k);
    CCTK_REAL const PDstandardNth33gt13 = PDstandardNth33(gt13, i, j, k);
    CCTK_REAL const PDstandardNth12gt13 = PDstandardNth12(gt13, i, j, k);
    CCTK_REAL const PDstandardNth13gt13 = PDstandardNth13(gt13, i, j, k);
    CCTK_REAL const PDstandardNth23gt13 = PDstandardNth23(gt13, i, j, k);
    CCTK_REAL const PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    CCTK_REAL const PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    CCTK_REAL const PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    CCTK_REAL const PDstandardNth11gt22 = PDstandardNth11(gt22, i, j, k);
    CCTK_REAL const PDstandardNth22gt22 = PDstandardNth22(gt22, i, j, k);
    CCTK_REAL const PDstandardNth33gt22 = PDstandardNth33(gt22, i, j, k);
    CCTK_REAL const PDstandardNth12gt22 = PDstandardNth12(gt22, i, j, k);
    CCTK_REAL const PDstandardNth13gt22 = PDstandardNth13(gt22, i, j, k);
    CCTK_REAL const PDstandardNth23gt22 = PDstandardNth23(gt22, i, j, k);
    CCTK_REAL const PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    CCTK_REAL const PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    CCTK_REAL const PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    CCTK_REAL const PDstandardNth11gt23 = PDstandardNth11(gt23, i, j, k);
    CCTK_REAL const PDstandardNth22gt23 = PDstandardNth22(gt23, i, j, k);
    CCTK_REAL const PDstandardNth33gt23 = PDstandardNth33(gt23, i, j, k);
    CCTK_REAL const PDstandardNth12gt23 = PDstandardNth12(gt23, i, j, k);
    CCTK_REAL const PDstandardNth13gt23 = PDstandardNth13(gt23, i, j, k);
    CCTK_REAL const PDstandardNth23gt23 = PDstandardNth23(gt23, i, j, k);
    CCTK_REAL const PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    CCTK_REAL const PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    CCTK_REAL const PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    CCTK_REAL const PDstandardNth11gt33 = PDstandardNth11(gt33, i, j, k);
    CCTK_REAL const PDstandardNth22gt33 = PDstandardNth22(gt33, i, j, k);
    CCTK_REAL const PDstandardNth33gt33 = PDstandardNth33(gt33, i, j, k);
    CCTK_REAL const PDstandardNth12gt33 = PDstandardNth12(gt33, i, j, k);
    CCTK_REAL const PDstandardNth13gt33 = PDstandardNth13(gt33, i, j, k);
    CCTK_REAL const PDstandardNth23gt33 = PDstandardNth23(gt33, i, j, k);
    CCTK_REAL const PDstandardNth1phi = PDstandardNth1(phi, i, j, k);
    CCTK_REAL const PDstandardNth2phi = PDstandardNth2(phi, i, j, k);
    CCTK_REAL const PDstandardNth3phi = PDstandardNth3(phi, i, j, k);
    CCTK_REAL const PDstandardNth11phi = PDstandardNth11(phi, i, j, k);
    CCTK_REAL const PDstandardNth22phi = PDstandardNth22(phi, i, j, k);
    CCTK_REAL const PDstandardNth33phi = PDstandardNth33(phi, i, j, k);
    CCTK_REAL const PDstandardNth12phi = PDstandardNth12(phi, i, j, k);
    CCTK_REAL const PDstandardNth13phi = PDstandardNth13(phi, i, j, k);
    CCTK_REAL const PDstandardNth23phi = PDstandardNth23(phi, i, j, k);
    CCTK_REAL const PDstandardNth1trK = PDstandardNth1(trK, i, j, k);
    CCTK_REAL const PDstandardNth2trK = PDstandardNth2(trK, i, j, k);
    CCTK_REAL const PDstandardNth3trK = PDstandardNth3(trK, i, j, k);
    CCTK_REAL const PDstandardNth1Xt1 = PDstandardNth1(Xt1, i, j, k);
    CCTK_REAL const PDstandardNth2Xt1 = PDstandardNth2(Xt1, i, j, k);
    CCTK_REAL const PDstandardNth3Xt1 = PDstandardNth3(Xt1, i, j, k);
    CCTK_REAL const PDstandardNth1Xt2 = PDstandardNth1(Xt2, i, j, k);
    CCTK_REAL const PDstandardNth2Xt2 = PDstandardNth2(Xt2, i, j, k);
    CCTK_REAL const PDstandardNth3Xt2 = PDstandardNth3(Xt2, i, j, k);
    CCTK_REAL const PDstandardNth1Xt3 = PDstandardNth1(Xt3, i, j, k);
    CCTK_REAL const PDstandardNth2Xt3 = PDstandardNth2(Xt3, i, j, k);
    CCTK_REAL const PDstandardNth3Xt3 = PDstandardNth3(Xt3, i, j, k);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu21 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu31 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu32 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL Gt111 = khalf*((gtu11*J11L - gtu21*J12L - 
      gtu31*J13L)*PDstandardNth1gt11 + (gtu11*J21L - gtu21*J22L - 
      gtu31*J23L)*PDstandardNth2gt11 + (gtu11*J31L - gtu21*J32L - 
      gtu31*J33L)*PDstandardNth3gt11 + 2*(J11L*(gtu21*PDstandardNth1gt12 + 
      gtu31*PDstandardNth1gt13) + J21L*(gtu21*PDstandardNth2gt12 + 
      gtu31*PDstandardNth2gt13) + J31L*(gtu21*PDstandardNth3gt12 + 
      gtu31*PDstandardNth3gt13)));
    
    CCTK_REAL Gt211 = khalf*((gtu21*J11L - gtu22*J12L - 
      gtu32*J13L)*PDstandardNth1gt11 + (gtu21*J21L - gtu22*J22L - 
      gtu32*J23L)*PDstandardNth2gt11 + (gtu21*J31L - gtu22*J32L - 
      gtu32*J33L)*PDstandardNth3gt11 + 2*(J11L*(gtu22*PDstandardNth1gt12 + 
      gtu32*PDstandardNth1gt13) + J21L*(gtu22*PDstandardNth2gt12 + 
      gtu32*PDstandardNth2gt13) + J31L*(gtu22*PDstandardNth3gt12 + 
      gtu32*PDstandardNth3gt13)));
    
    CCTK_REAL Gt311 = khalf*((gtu31*J11L - gtu32*J12L - 
      gtu33*J13L)*PDstandardNth1gt11 + (gtu31*J21L - gtu32*J22L - 
      gtu33*J23L)*PDstandardNth2gt11 + (gtu31*J31L - gtu32*J32L - 
      gtu33*J33L)*PDstandardNth3gt11 + 2*(J11L*(gtu32*PDstandardNth1gt12 + 
      gtu33*PDstandardNth1gt13) + J21L*(gtu32*PDstandardNth2gt12 + 
      gtu33*PDstandardNth2gt13) + J31L*(gtu32*PDstandardNth3gt12 + 
      gtu33*PDstandardNth3gt13)));
    
    CCTK_REAL Gt112 = khalf*(gtu11*(J12L*PDstandardNth1gt11 + 
      J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
      gtu21*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
      J31L*PDstandardNth3gt22) + gtu31*(-(J13L*PDstandardNth1gt12) + 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23));
    
    CCTK_REAL Gt212 = khalf*(gtu21*(J12L*PDstandardNth1gt11 + 
      J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
      gtu22*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
      J31L*PDstandardNth3gt22) + gtu32*(-(J13L*PDstandardNth1gt12) + 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23));
    
    CCTK_REAL Gt312 = khalf*(gtu31*(J12L*PDstandardNth1gt11 + 
      J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
      gtu32*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
      J31L*PDstandardNth3gt22) + gtu33*(-(J13L*PDstandardNth1gt12) + 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23));
    
    CCTK_REAL Gt113 = khalf*(gtu11*(J13L*PDstandardNth1gt11 + 
      J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
      gtu21*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
      J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
      J31L*PDstandardNth3gt23) + gtu31*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL Gt213 = khalf*(gtu21*(J13L*PDstandardNth1gt11 + 
      J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
      gtu22*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
      J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
      J31L*PDstandardNth3gt23) + gtu32*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL Gt313 = khalf*(gtu31*(J13L*PDstandardNth1gt11 + 
      J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
      gtu32*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
      J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
      J31L*PDstandardNth3gt23) + gtu33*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL Gt122 = khalf*(gtu11*(-(J11L*PDstandardNth1gt22) + 
      2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
      J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
      J31L*PDstandardNth3gt22) + gtu21*(J12L*PDstandardNth1gt22 + 
      J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
      gtu31*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
      J33L*PDstandardNth3gt22 - 2*(J12L*PDstandardNth1gt23 + 
      J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL Gt222 = khalf*(gtu21*(-(J11L*PDstandardNth1gt22) + 
      2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
      J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
      J31L*PDstandardNth3gt22) + gtu22*(J12L*PDstandardNth1gt22 + 
      J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
      gtu32*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
      J33L*PDstandardNth3gt22 - 2*(J12L*PDstandardNth1gt23 + 
      J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL Gt322 = khalf*(gtu31*(-(J11L*PDstandardNth1gt22) + 
      2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
      J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
      J31L*PDstandardNth3gt22) + gtu32*(J12L*PDstandardNth1gt22 + 
      J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
      gtu33*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
      J33L*PDstandardNth3gt22 - 2*(J12L*PDstandardNth1gt23 + 
      J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL Gt123 = khalf*(gtu21*(J13L*PDstandardNth1gt22 + 
      J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
      gtu11*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
      J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt23) + gtu31*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL Gt223 = khalf*(gtu22*(J13L*PDstandardNth1gt22 + 
      J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
      gtu21*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
      J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt23) + gtu32*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL Gt323 = khalf*(gtu32*(J13L*PDstandardNth1gt22 + 
      J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
      gtu31*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
      J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt23) + gtu33*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL Gt133 = khalf*(gtu11*(-(J11L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
      J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt33) + gtu21*(-(J12L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
      J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
      J32L*PDstandardNth3gt33) + gtu31*(J13L*PDstandardNth1gt33 + 
      J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL Gt233 = khalf*(gtu21*(-(J11L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
      J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt33) + gtu22*(-(J12L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
      J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
      J32L*PDstandardNth3gt33) + gtu32*(J13L*PDstandardNth1gt33 + 
      J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL Gt333 = khalf*(gtu31*(-(J11L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
      J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt33) + gtu32*(-(J12L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
      J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
      J32L*PDstandardNth3gt33) + gtu33*(J13L*PDstandardNth1gt33 + 
      J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL Rt11 = (Gt113*(gt13L*Gt312 + 3*(gt12L*Gt212 + gt13L*Gt312)) 
      + gt12L*(Gt213*(4*Gt112 + 2*Gt222) + Gt212*(Gt113 + 2*Gt223) + 
      2*(Gt233*Gt312 + Gt223*Gt313)) + gt11L*(6*Gt112*Gt113 + 2*(Gt122*Gt213 
      + Gt133*Gt312 + Gt123*(Gt212 + Gt313))) + gt13L*(2*Gt213*Gt322 + 
      Gt313*(4*Gt112 + 2*Gt323)) + 2*(Gt213*(Gt212*gt22L + gt23L*Gt312) + 
      Gt212*(gt23L*Gt313 + gt13L*Gt323) + Gt312*(gt13L*Gt333 + 
      Gt313*gt33L)))*gtu32 + J11L*(gt11L*PDstandardNth1Xt1 + 
      gt12L*PDstandardNth1Xt2 + gt13L*PDstandardNth1Xt3) + 
      J21L*(gt11L*PDstandardNth2Xt1 + gt12L*PDstandardNth2Xt2 + 
      gt13L*PDstandardNth2Xt3) + J31L*(gt11L*PDstandardNth3Xt1 + 
      gt12L*PDstandardNth3Xt2 + gt13L*PDstandardNth3Xt3) + (Gt111*gt11L + 
      gt12L*Gt211 + gt13L*Gt311)*Xt1L + (Gt112*gt11L + gt12L*Gt212 + 
      gt13L*Gt312)*Xt2L + (Gt113*gt11L + gt12L*Gt213 + gt13L*Gt313)*Xt3L + 
      gtu21*(Gt112*(gt13L*Gt311 + 3*(gt12L*Gt211 + gt13L*Gt311)) + 
      gt11L*(Gt112*(6*Gt111 + 2*Gt212) + 2*(Gt122*Gt211 + Gt123*Gt311 + 
      Gt113*Gt312)) + 2*(Gt212*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + 
      gt13L*(Gt211*Gt322 + Gt311*Gt323)) + Gt312*(gt13L*(4*Gt111 + 2*Gt313) + 
      2*(Gt211*gt23L + Gt311*gt33L)) + gt12L*(4*Gt111*Gt212 + Gt211*(Gt112 + 
      2*Gt222) + 2*(Gt223*Gt311 + Gt213*Gt312 + SQR(Gt212)))) + 
      gtu11*(4*Gt111*(gt12L*Gt211 + gt13L*Gt311) + 2*(gt12L*(Gt211*Gt212 + 
      Gt213*Gt311) + Gt211*(gt23L*Gt311 + gt13L*Gt312) + gt13L*Gt311*Gt313) + 
      gt11L*(2*(Gt112*Gt211 + Gt113*Gt311) + 3*SQR(Gt111)) + gt22L*SQR(Gt211) 
      + gt33L*SQR(Gt311)) + gtu22*(4*Gt112*(gt12L*Gt212 + gt13L*Gt312) + 
      2*(gt12L*(Gt212*Gt222 + Gt223*Gt312) + Gt212*(gt23L*Gt312 + 
      gt13L*Gt322) + gt13L*Gt312*Gt323) + gt11L*(2*(Gt122*Gt212 + 
      Gt123*Gt312) + 3*SQR(Gt112)) + gt22L*SQR(Gt212) + gt33L*SQR(Gt312)) + 
      gtu33*(4*Gt113*(gt12L*Gt213 + gt13L*Gt313) + 2*(gt12L*(Gt213*Gt223 + 
      Gt233*Gt313) + Gt213*(gt23L*Gt313 + gt13L*Gt323) + gt13L*Gt313*Gt333) + 
      gt11L*(2*(Gt123*Gt213 + Gt133*Gt313) + 3*SQR(Gt113)) + gt22L*SQR(Gt213) 
      + gt33L*SQR(Gt313)) + gtu31*(Gt113*(gt13L*Gt311 + 3*(gt12L*Gt211 + 
      gt13L*Gt311)) + gt11L*(2*(Gt123*Gt211 + Gt112*Gt213 + Gt133*Gt311) + 
      Gt113*(6*Gt111 + 2*Gt313)) + gt12L*(Gt211*(Gt113 + 2*Gt223) + 
      2*Gt233*Gt311 + Gt213*(4*Gt111 + 2*(Gt212 + Gt313))) + 
      2*(Gt213*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312) + gt13L*Gt311*Gt333 
      + Gt313*(Gt211*gt23L + Gt311*gt33L)) + gt13L*(4*Gt111*Gt313 + 
      2*(Gt211*Gt323 + SQR(Gt313)))) + 
      khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt11 + 
      J12L*J21L*PDstandardNth12gt11 + J11L*J22L*PDstandardNth12gt11 + 
      J12L*J31L*PDstandardNth13gt11 + J11L*J32L*PDstandardNth13gt11 + 
      dJ112L*PDstandardNth1gt11 + J21L*J22L*PDstandardNth22gt11 + 
      J22L*J31L*PDstandardNth23gt11 + J21L*J32L*PDstandardNth23gt11 + 
      dJ212L*PDstandardNth2gt11 + J31L*J32L*PDstandardNth33gt11 + 
      dJ312L*PDstandardNth3gt11) + gtu31*(J11L*J13L*PDstandardNth11gt11 + 
      J13L*J21L*PDstandardNth12gt11 + J11L*J23L*PDstandardNth12gt11 + 
      J13L*J31L*PDstandardNth13gt11 + J11L*J33L*PDstandardNth13gt11 + 
      dJ113L*PDstandardNth1gt11 + J21L*J23L*PDstandardNth22gt11 + 
      J23L*J31L*PDstandardNth23gt11 + J21L*J33L*PDstandardNth23gt11 + 
      dJ213L*PDstandardNth2gt11 + J31L*J33L*PDstandardNth33gt11 + 
      dJ313L*PDstandardNth3gt11) + gtu32*(J12L*J13L*PDstandardNth11gt11 + 
      J13L*J22L*PDstandardNth12gt11 + J12L*J23L*PDstandardNth12gt11 + 
      J13L*J32L*PDstandardNth13gt11 + J12L*J33L*PDstandardNth13gt11 + 
      dJ123L*PDstandardNth1gt11 + J22L*J23L*PDstandardNth22gt11 + 
      J23L*J32L*PDstandardNth23gt11 + J22L*J33L*PDstandardNth23gt11 + 
      dJ223L*PDstandardNth2gt11 + J32L*J33L*PDstandardNth33gt11 + 
      dJ323L*PDstandardNth3gt11)) - gtu11*(2*J11L*J21L*PDstandardNth12gt11 + 
      2*J11L*J31L*PDstandardNth13gt11 + dJ111L*PDstandardNth1gt11 + 
      2*J21L*J31L*PDstandardNth23gt11 + dJ211L*PDstandardNth2gt11 + 
      dJ311L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J11L) + 
      PDstandardNth22gt11*SQR(J21L) + PDstandardNth33gt11*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt11 + 
      2*J12L*J32L*PDstandardNth13gt11 + dJ122L*PDstandardNth1gt11 + 
      2*J22L*J32L*PDstandardNth23gt11 + dJ222L*PDstandardNth2gt11 + 
      dJ322L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J12L) + 
      PDstandardNth22gt11*SQR(J22L) + PDstandardNth33gt11*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt11 + 
      2*J13L*J33L*PDstandardNth13gt11 + dJ133L*PDstandardNth1gt11 + 
      2*J23L*J33L*PDstandardNth23gt11 + dJ233L*PDstandardNth2gt11 + 
      dJ333L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J13L) + 
      PDstandardNth22gt11*SQR(J23L) + PDstandardNth33gt11*SQR(J33L)));
    
    CCTK_REAL Rt12 = khalf*((gt12L*J11L + gt11L*J12L)*PDstandardNth1Xt1 + 
      (gt22L*J11L + gt12L*J12L)*PDstandardNth1Xt2 + (gt23L*J11L + 
      gt13L*J12L)*PDstandardNth1Xt3 + (gt12L*J21L + 
      gt11L*J22L)*PDstandardNth2Xt1 + (gt22L*J21L + 
      gt12L*J22L)*PDstandardNth2Xt2 + (gt23L*J21L + 
      gt13L*J22L)*PDstandardNth2Xt3 - 2*(gtu21*(J11L*J12L*PDstandardNth11gt12 
      + J12L*J21L*PDstandardNth12gt12 + J11L*J22L*PDstandardNth12gt12 + 
      J12L*J31L*PDstandardNth13gt12 + J11L*J32L*PDstandardNth13gt12 + 
      dJ112L*PDstandardNth1gt12 + J21L*J22L*PDstandardNth22gt12 + 
      J22L*J31L*PDstandardNth23gt12 + J21L*J32L*PDstandardNth23gt12 + 
      dJ212L*PDstandardNth2gt12 + J31L*J32L*PDstandardNth33gt12 + 
      dJ312L*PDstandardNth3gt12) + gtu31*(J11L*J13L*PDstandardNth11gt12 + 
      J13L*J21L*PDstandardNth12gt12 + J11L*J23L*PDstandardNth12gt12 + 
      J13L*J31L*PDstandardNth13gt12 + J11L*J33L*PDstandardNth13gt12 + 
      dJ113L*PDstandardNth1gt12 + J21L*J23L*PDstandardNth22gt12 + 
      J23L*J31L*PDstandardNth23gt12 + J21L*J33L*PDstandardNth23gt12 + 
      dJ213L*PDstandardNth2gt12 + J31L*J33L*PDstandardNth33gt12 + 
      dJ313L*PDstandardNth3gt12) + gtu32*(J12L*J13L*PDstandardNth11gt12 + 
      J13L*J22L*PDstandardNth12gt12 + J12L*J23L*PDstandardNth12gt12 + 
      J13L*J32L*PDstandardNth13gt12 + J12L*J33L*PDstandardNth13gt12 + 
      dJ123L*PDstandardNth1gt12 + J22L*J23L*PDstandardNth22gt12 + 
      J23L*J32L*PDstandardNth23gt12 + J22L*J33L*PDstandardNth23gt12 + 
      dJ223L*PDstandardNth2gt12 + J32L*J33L*PDstandardNth33gt12 + 
      dJ323L*PDstandardNth3gt12)) + (gt12L*J31L + 
      gt11L*J32L)*PDstandardNth3Xt1 + (gt22L*J31L + 
      gt12L*J32L)*PDstandardNth3Xt2 + (gt23L*J31L + 
      gt13L*J32L)*PDstandardNth3Xt3 + (Gt112*gt11L + Gt111*gt12L + 
      gt12L*Gt212 + Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)*Xt1L + 
      (gt11L*Gt122 + Gt112*gt12L + gt12L*Gt222 + Gt212*gt22L + gt23L*Gt312 + 
      gt13L*Gt322)*Xt2L + (gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + 
      Gt213*gt22L + gt23L*Gt313 + gt13L*Gt323)*Xt3L + 2*((Gt123*gt12L*Gt211 + 
      Gt113*gt12L*Gt212 + 2*Gt112*gt12L*Gt213 + gt12L*Gt212*Gt223 + 
      Gt212*Gt213*gt22L + Gt211*Gt223*gt22L + gt12L*Gt133*Gt311 + 
      gt22L*Gt233*Gt311 + Gt113*gt13L*Gt312 + gt12L*Gt233*Gt312 + 
      Gt213*gt23L*Gt312 + gt11L*(2*Gt112*Gt113 + Gt123*Gt212 + Gt133*Gt312) + 
      2*Gt112*gt13L*Gt313 + Gt212*gt23L*Gt313 + Gt111*(Gt113*gt12L + 
      Gt213*gt22L + gt23L*Gt313) + gt13L*Gt212*Gt323 + Gt211*gt23L*Gt323 + 
      gt23L*Gt311*Gt333 + gt13L*Gt312*Gt333 + Gt312*Gt313*gt33L)*gtu31 + 
      (Gt123*gt12L*Gt212 + 2*Gt122*gt12L*Gt213 + Gt113*gt12L*Gt222 + 
      gt12L*Gt222*Gt223 + Gt213*Gt222*gt22L + Gt212*Gt223*gt22L + 
      gt12L*Gt133*Gt312 + gt22L*Gt233*Gt312 + 2*Gt122*gt13L*Gt313 + 
      Gt222*gt23L*Gt313 + Gt112*(Gt113*gt12L + Gt213*gt22L + gt23L*Gt313) + 
      Gt113*gt13L*Gt322 + gt12L*Gt233*Gt322 + Gt213*gt23L*Gt322 + 
      gt11L*(2*Gt113*Gt122 + Gt123*Gt222 + Gt133*Gt322) + gt13L*Gt222*Gt323 + 
      Gt212*gt23L*Gt323 + gt23L*Gt312*Gt333 + gt13L*Gt322*Gt333 + 
      Gt313*Gt322*gt33L)*gtu32 + gtu11*(3*Gt112*gt12L*Gt211 + 
      2*Gt211*Gt212*gt22L + Gt113*gt12L*Gt311 + 2*Gt112*gt13L*Gt311 + 
      Gt213*gt22L*Gt311 + Gt212*gt23L*Gt311 + gt13L*Gt212*Gt312 + 
      gt12L*Gt213*Gt312 + 2*Gt211*gt23L*Gt312 + gt11L*(2*Gt111*Gt112 + 
      Gt112*Gt212 + Gt113*Gt312) + Gt111*(gt12L*Gt212 + Gt211*gt22L + 
      gt23L*Gt311 + gt13L*Gt312) + gt23L*Gt311*Gt313 + gt13L*Gt312*Gt313 + 
      Gt311*Gt312*gt33L + gt12L*SQR(Gt111) + gt12L*SQR(Gt212)) + 
      gtu21*(Gt112*gt11L*Gt222 + Gt112*Gt211*gt22L + Gt211*Gt222*gt22L + 
      2*Gt122*gt13L*Gt311 + Gt112*gt23L*Gt311 + Gt222*gt23L*Gt311 + 
      gt13L*Gt222*Gt312 + Gt213*gt22L*Gt312 + Gt212*gt23L*Gt312 + 
      gt23L*Gt312*Gt313 + Gt113*gt11L*Gt322 + Gt211*gt23L*Gt322 + 
      gt13L*Gt313*Gt322 + Gt111*(2*gt11L*Gt122 + Gt112*gt12L + gt12L*Gt222 + 
      gt13L*Gt322) + gt12L*(2*Gt122*Gt211 + Gt112*Gt212 + Gt212*Gt222 + 
      Gt113*Gt312 + Gt213*Gt322) + Gt311*Gt322*gt33L + gt22L*SQR(Gt212)) + 
      gtu22*(gt11L*Gt122*Gt222 + 2*Gt212*Gt222*gt22L + 2*Gt122*gt13L*Gt312 + 
      Gt223*gt22L*Gt312 + Gt222*gt23L*Gt312 + gt11L*Gt123*Gt322 + 
      gt13L*Gt222*Gt322 + 2*Gt212*gt23L*Gt322 + Gt112*(2*gt11L*Gt122 + 
      gt12L*Gt222 + Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322) + 
      gt23L*Gt312*Gt323 + gt13L*Gt322*Gt323 + Gt312*Gt322*gt33L + 
      gt12L*SQR(Gt112) + gt12L*(3*Gt122*Gt212 + Gt123*Gt312 + Gt223*Gt322 + 
      SQR(Gt222))) + gtu33*(gt11L*Gt123*Gt223 + 2*Gt213*Gt223*gt22L + 
      2*Gt123*gt13L*Gt313 + gt22L*Gt233*Gt313 + Gt223*gt23L*Gt313 + 
      gt11L*Gt133*Gt323 + gt13L*Gt223*Gt323 + 2*Gt213*gt23L*Gt323 + 
      Gt113*(2*gt11L*Gt123 + gt12L*Gt223 + Gt213*gt22L + gt23L*Gt313 + 
      gt13L*Gt323) + gt23L*Gt313*Gt333 + gt13L*Gt323*Gt333 + 
      Gt313*Gt323*gt33L + gt12L*SQR(Gt113) + gt12L*(3*Gt123*Gt213 + 
      Gt133*Gt313 + Gt233*Gt323 + SQR(Gt223))) + gtu21*(Gt122*gt12L*Gt211 + 
      3*Gt112*gt12L*Gt212 + gt12L*Gt212*Gt222 + Gt211*Gt222*gt22L + 
      Gt123*gt12L*Gt311 + Gt223*gt22L*Gt311 + 3*Gt112*gt13L*Gt312 + 
      gt12L*Gt223*Gt312 + 2*Gt212*gt23L*Gt312 + Gt111*(Gt112*gt12L + 
      Gt212*gt22L + gt23L*Gt312) + gt13L*Gt212*Gt322 + Gt211*gt23L*Gt322 + 
      gt23L*Gt311*Gt323 + gt13L*Gt312*Gt323 + gt11L*(Gt122*Gt212 + 
      Gt123*Gt312 + 2*SQR(Gt112)) + gt22L*SQR(Gt212) + gt33L*SQR(Gt312)) + 
      gtu31*(Gt112*gt11L*Gt223 + Gt113*Gt211*gt22L + Gt212*Gt213*gt22L + 
      Gt211*Gt223*gt22L + 2*Gt123*gt13L*Gt311 + Gt113*gt23L*Gt311 + 
      Gt223*gt23L*Gt311 + gt13L*Gt223*Gt312 + Gt213*gt23L*Gt312 + 
      Gt213*gt22L*Gt313 + Gt113*gt11L*Gt323 + Gt211*gt23L*Gt323 + 
      gt13L*Gt313*Gt323 + Gt111*(2*gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + 
      gt13L*Gt323) + gt12L*(2*Gt123*Gt211 + Gt112*Gt213 + Gt212*Gt223 + 
      Gt113*Gt313 + Gt213*Gt323) + Gt311*Gt323*gt33L + gt23L*SQR(Gt313)) + 
      gtu32*(gt11L*Gt122*Gt223 + Gt113*Gt212*gt22L + Gt213*Gt222*gt22L + 
      Gt212*Gt223*gt22L + 2*Gt123*gt13L*Gt312 + Gt113*gt23L*Gt312 + 
      Gt223*gt23L*Gt312 + Gt223*gt22L*Gt313 + gt13L*Gt223*Gt322 + 
      Gt213*gt23L*Gt322 + gt11L*Gt123*Gt323 + Gt212*gt23L*Gt323 + 
      gt23L*Gt313*Gt323 + Gt112*(2*gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + 
      gt13L*Gt323) + gt12L*(Gt122*Gt213 + Gt123*(2*Gt212 + Gt313) + 
      Gt223*(Gt222 + Gt323)) + Gt312*Gt323*gt33L + gt13L*SQR(Gt323))) - 
      gtu11*(2*J11L*J21L*PDstandardNth12gt12 + 
      2*J11L*J31L*PDstandardNth13gt12 + dJ111L*PDstandardNth1gt12 + 
      2*J21L*J31L*PDstandardNth23gt12 + dJ211L*PDstandardNth2gt12 + 
      dJ311L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J11L) + 
      PDstandardNth22gt12*SQR(J21L) + PDstandardNth33gt12*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt12 + 
      2*J12L*J32L*PDstandardNth13gt12 + dJ122L*PDstandardNth1gt12 + 
      2*J22L*J32L*PDstandardNth23gt12 + dJ222L*PDstandardNth2gt12 + 
      dJ322L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J12L) + 
      PDstandardNth22gt12*SQR(J22L) + PDstandardNth33gt12*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt12 + 
      2*J13L*J33L*PDstandardNth13gt12 + dJ133L*PDstandardNth1gt12 + 
      2*J23L*J33L*PDstandardNth23gt12 + dJ233L*PDstandardNth2gt12 + 
      dJ333L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J13L) + 
      PDstandardNth22gt12*SQR(J23L) + PDstandardNth33gt12*SQR(J33L)));
    
    CCTK_REAL Rt13 = khalf*((gt13L*J11L + gt11L*J13L)*PDstandardNth1Xt1 + 
      (gt23L*J11L + gt12L*J13L)*PDstandardNth1Xt2 + (gt33L*J11L + 
      gt13L*J13L)*PDstandardNth1Xt3 + (gt13L*J21L + 
      gt11L*J23L)*PDstandardNth2Xt1 + (gt23L*J21L + 
      gt12L*J23L)*PDstandardNth2Xt2 + (gt33L*J21L + 
      gt13L*J23L)*PDstandardNth2Xt3 - 2*(gtu21*(J11L*J12L*PDstandardNth11gt13 
      + J12L*J21L*PDstandardNth12gt13 + J11L*J22L*PDstandardNth12gt13 + 
      J12L*J31L*PDstandardNth13gt13 + J11L*J32L*PDstandardNth13gt13 + 
      dJ112L*PDstandardNth1gt13 + J21L*J22L*PDstandardNth22gt13 + 
      J22L*J31L*PDstandardNth23gt13 + J21L*J32L*PDstandardNth23gt13 + 
      dJ212L*PDstandardNth2gt13 + J31L*J32L*PDstandardNth33gt13 + 
      dJ312L*PDstandardNth3gt13) + gtu31*(J11L*J13L*PDstandardNth11gt13 + 
      J13L*J21L*PDstandardNth12gt13 + J11L*J23L*PDstandardNth12gt13 + 
      J13L*J31L*PDstandardNth13gt13 + J11L*J33L*PDstandardNth13gt13 + 
      dJ113L*PDstandardNth1gt13 + J21L*J23L*PDstandardNth22gt13 + 
      J23L*J31L*PDstandardNth23gt13 + J21L*J33L*PDstandardNth23gt13 + 
      dJ213L*PDstandardNth2gt13 + J31L*J33L*PDstandardNth33gt13 + 
      dJ313L*PDstandardNth3gt13) + gtu32*(J12L*J13L*PDstandardNth11gt13 + 
      J13L*J22L*PDstandardNth12gt13 + J12L*J23L*PDstandardNth12gt13 + 
      J13L*J32L*PDstandardNth13gt13 + J12L*J33L*PDstandardNth13gt13 + 
      dJ123L*PDstandardNth1gt13 + J22L*J23L*PDstandardNth22gt13 + 
      J23L*J32L*PDstandardNth23gt13 + J22L*J33L*PDstandardNth23gt13 + 
      dJ223L*PDstandardNth2gt13 + J32L*J33L*PDstandardNth33gt13 + 
      dJ323L*PDstandardNth3gt13)) + (gt13L*J31L + 
      gt11L*J33L)*PDstandardNth3Xt1 + (gt23L*J31L + 
      gt12L*J33L)*PDstandardNth3Xt2 + (gt33L*J31L + 
      gt13L*J33L)*PDstandardNth3Xt3 + (Gt113*gt11L + Gt111*gt13L + 
      gt12L*Gt213 + Gt211*gt23L + gt13L*Gt313 + Gt311*gt33L)*Xt1L + 
      (gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + gt13L*Gt323 + 
      Gt312*gt33L)*Xt2L + (gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + 
      Gt213*gt23L + gt13L*Gt333 + Gt313*gt33L)*Xt3L + 2*((Gt122*gt13L*Gt211 + 
      2*Gt113*gt12L*Gt212 + Gt112*gt12L*Gt213 + gt12L*Gt213*Gt222 + 
      Gt212*Gt213*gt22L + Gt211*Gt222*gt23L + Gt123*gt13L*Gt311 + 
      Gt223*gt23L*Gt311 + 2*Gt113*gt13L*Gt312 + Gt213*gt23L*Gt312 + 
      Gt112*gt13L*Gt313 + gt12L*Gt223*Gt313 + Gt212*gt23L*Gt313 + 
      gt11L*(2*Gt112*Gt113 + Gt122*Gt213 + Gt123*Gt313) + gt13L*Gt213*Gt322 + 
      gt13L*Gt313*Gt323 + Gt312*Gt313*gt33L + Gt211*Gt322*gt33L + 
      Gt311*Gt323*gt33L + Gt111*(Gt112*gt13L + Gt212*gt23L + 
      Gt312*gt33L))*gtu21 + (Gt122*gt13L*Gt213 + gt11L*Gt122*Gt233 + 
      Gt212*gt22L*Gt233 + Gt113*Gt212*gt23L + Gt213*Gt222*gt23L + 
      2*Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + Gt123*gt13L*Gt313 + 
      Gt223*gt23L*Gt313 + gt13L*Gt233*Gt322 + gt11L*Gt123*Gt333 + 
      Gt212*gt23L*Gt333 + gt13L*Gt323*Gt333 + Gt112*(2*gt11L*Gt133 + 
      Gt113*gt13L + gt12L*Gt233 + gt13L*Gt333) + gt12L*(2*Gt133*Gt212 + 
      Gt222*Gt233 + Gt223*Gt333) + Gt113*Gt312*gt33L + Gt213*Gt322*gt33L + 
      Gt313*Gt323*gt33L + Gt312*Gt333*gt33L)*gtu32 + 
      gtu21*(2*Gt123*gt12L*Gt211 + Gt112*gt13L*Gt212 + gt12L*Gt212*Gt223 + 
      Gt211*Gt223*gt22L + Gt112*Gt211*gt23L + 2*Gt123*gt13L*Gt311 + 
      Gt223*gt23L*Gt311 + Gt113*gt13L*Gt312 + gt13L*Gt223*Gt312 + 
      Gt213*gt23L*Gt312 + gt12L*Gt213*Gt323 + Gt211*gt23L*Gt323 + 
      gt13L*Gt313*Gt323 + gt11L*(2*Gt111*Gt123 + Gt112*Gt223 + Gt113*Gt323) + 
      Gt111*(Gt112*gt13L + gt12L*Gt223 + gt13L*Gt323) + Gt112*Gt311*gt33L + 
      Gt212*Gt312*gt33L + Gt312*Gt313*gt33L + Gt311*Gt323*gt33L + 
      gt23L*SQR(Gt212)) + gtu32*(Gt123*gt13L*Gt212 + 2*Gt123*gt12L*Gt213 + 
      Gt113*gt12L*Gt223 + Gt213*Gt223*gt22L + Gt212*Gt223*gt23L + 
      Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + 2*Gt123*gt13L*Gt313 + 
      Gt223*gt23L*Gt313 + Gt113*gt13L*Gt323 + gt13L*Gt223*Gt323 + 
      gt12L*Gt233*Gt323 + Gt213*gt23L*Gt323 + gt11L*(2*Gt113*Gt123 + 
      Gt123*Gt223 + Gt133*Gt323) + gt13L*Gt323*Gt333 + Gt212*Gt323*gt33L + 
      Gt313*Gt323*gt33L + Gt312*Gt333*gt33L + Gt112*(Gt113*gt13L + 
      Gt213*gt23L + Gt313*gt33L) + gt12L*SQR(Gt223)) + 
      gtu11*(2*Gt113*gt12L*Gt211 + Gt112*gt13L*Gt211 + gt12L*Gt212*Gt213 + 
      Gt211*Gt213*gt22L + Gt211*Gt212*gt23L + 3*Gt113*gt13L*Gt311 + 
      2*Gt213*gt23L*Gt311 + gt13L*Gt213*Gt312 + gt12L*Gt213*Gt313 + 
      Gt211*gt23L*Gt313 + gt11L*(2*Gt111*Gt113 + Gt112*Gt213 + Gt113*Gt313) + 
      Gt211*Gt312*gt33L + 2*Gt311*Gt313*gt33L + Gt111*(gt12L*Gt213 + 
      Gt211*gt23L + gt13L*Gt313 + Gt311*gt33L) + gt13L*SQR(Gt111) + 
      gt13L*SQR(Gt313)) + gtu31*(Gt112*gt13L*Gt213 + Gt112*gt11L*Gt233 + 
      Gt211*gt22L*Gt233 + Gt113*Gt211*gt23L + Gt212*Gt213*gt23L + 
      2*Gt133*gt13L*Gt311 + Gt233*gt23L*Gt311 + gt13L*Gt233*Gt312 + 
      Gt113*gt13L*Gt313 + Gt213*gt23L*Gt313 + Gt113*gt11L*Gt333 + 
      Gt211*gt23L*Gt333 + gt13L*Gt313*Gt333 + Gt111*(2*gt11L*Gt133 + 
      Gt113*gt13L + gt12L*Gt233 + gt13L*Gt333) + gt12L*(2*Gt133*Gt211 + 
      Gt212*Gt233 + Gt213*Gt333) + Gt113*Gt311*gt33L + Gt213*Gt312*gt33L + 
      Gt311*Gt333*gt33L + gt33L*SQR(Gt313)) + gtu31*(Gt123*gt13L*Gt211 + 
      3*Gt113*gt12L*Gt213 + gt12L*Gt213*Gt223 + Gt211*Gt223*gt23L + 
      Gt133*gt13L*Gt311 + Gt233*gt23L*Gt311 + 3*Gt113*gt13L*Gt313 + 
      gt12L*Gt233*Gt313 + 2*Gt213*gt23L*Gt313 + gt13L*Gt213*Gt323 + 
      gt13L*Gt313*Gt333 + Gt211*Gt323*gt33L + Gt311*Gt333*gt33L + 
      Gt111*(Gt113*gt13L + Gt213*gt23L + Gt313*gt33L) + gt11L*(Gt123*Gt213 + 
      Gt133*Gt313 + 2*SQR(Gt113)) + gt22L*SQR(Gt213) + gt33L*SQR(Gt313)) + 
      gtu22*(2*Gt123*gt12L*Gt212 + Gt122*gt13L*Gt212 + gt12L*Gt222*Gt223 + 
      Gt212*Gt223*gt22L + Gt212*Gt222*gt23L + 3*Gt123*gt13L*Gt312 + 
      2*Gt223*gt23L*Gt312 + gt13L*Gt223*Gt322 + gt12L*Gt223*Gt323 + 
      Gt212*gt23L*Gt323 + gt11L*(2*Gt112*Gt123 + Gt122*Gt223 + Gt123*Gt323) + 
      Gt212*Gt322*gt33L + 2*Gt312*Gt323*gt33L + Gt112*(gt12L*Gt223 + 
      Gt212*gt23L + gt13L*Gt323 + Gt312*gt33L) + gt13L*SQR(Gt112) + 
      gt13L*SQR(Gt323)) + gtu33*(2*gt12L*Gt133*Gt213 + Gt123*gt13L*Gt213 + 
      gt11L*Gt123*Gt233 + gt12L*Gt223*Gt233 + Gt213*gt22L*Gt233 + 
      Gt213*Gt223*gt23L + 3*Gt133*gt13L*Gt313 + 2*Gt233*gt23L*Gt313 + 
      gt13L*Gt233*Gt323 + gt11L*Gt133*Gt333 + gt12L*Gt233*Gt333 + 
      Gt213*gt23L*Gt333 + Gt213*Gt323*gt33L + 2*Gt313*Gt333*gt33L + 
      Gt113*(2*gt11L*Gt133 + gt12L*Gt233 + Gt213*gt23L + gt13L*Gt333 + 
      Gt313*gt33L) + gt13L*SQR(Gt113) + gt13L*SQR(Gt333))) - 
      gtu11*(2*J11L*J21L*PDstandardNth12gt13 + 
      2*J11L*J31L*PDstandardNth13gt13 + dJ111L*PDstandardNth1gt13 + 
      2*J21L*J31L*PDstandardNth23gt13 + dJ211L*PDstandardNth2gt13 + 
      dJ311L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J11L) + 
      PDstandardNth22gt13*SQR(J21L) + PDstandardNth33gt13*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt13 + 
      2*J12L*J32L*PDstandardNth13gt13 + dJ122L*PDstandardNth1gt13 + 
      2*J22L*J32L*PDstandardNth23gt13 + dJ222L*PDstandardNth2gt13 + 
      dJ322L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J12L) + 
      PDstandardNth22gt13*SQR(J22L) + PDstandardNth33gt13*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt13 + 
      2*J13L*J33L*PDstandardNth13gt13 + dJ133L*PDstandardNth1gt13 + 
      2*J23L*J33L*PDstandardNth23gt13 + dJ233L*PDstandardNth2gt13 + 
      dJ333L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J13L) + 
      PDstandardNth22gt13*SQR(J23L) + PDstandardNth33gt13*SQR(J33L)));
    
    CCTK_REAL Rt22 = (Gt223*(3*Gt112*gt12L + 6*Gt212*gt22L + 
      4*gt23L*Gt312) + Gt123*(Gt112*gt11L + gt12L*(2*Gt111 + 4*Gt212) + 
      2*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)) + Gt112*(gt11L*Gt123 + 
      gt12L*(2*Gt113 + Gt223) + 2*(Gt213*gt22L + gt23L*Gt313) + gt13L*Gt323) 
      + 2*(Gt113*gt12L*Gt323 + Gt312*(gt12L*Gt133 + gt22L*Gt233 + 
      gt23L*Gt333)) + Gt323*(Gt112*gt13L + 4*Gt212*gt23L + 2*(Gt213*gt22L + 
      gt23L*Gt313 + Gt312*gt33L)))*gtu31 + J12L*(gt12L*PDstandardNth1Xt1 + 
      gt22L*PDstandardNth1Xt2 + gt23L*PDstandardNth1Xt3) + 
      J22L*(gt12L*PDstandardNth2Xt1 + gt22L*PDstandardNth2Xt2 + 
      gt23L*PDstandardNth2Xt3) + J32L*(gt12L*PDstandardNth3Xt1 + 
      gt22L*PDstandardNth3Xt2 + gt23L*PDstandardNth3Xt3) + (Gt112*gt12L + 
      Gt212*gt22L + gt23L*Gt312)*Xt1L + (Gt122*gt12L + Gt222*gt22L + 
      gt23L*Gt322)*Xt2L + (Gt123*gt12L + Gt223*gt22L + gt23L*Gt323)*Xt3L + 
      gtu21*(Gt222*(3*Gt112*gt12L + 6*Gt212*gt22L + 4*gt23L*Gt312) + 
      Gt122*(Gt112*gt11L + gt12L*(2*Gt111 + 4*Gt212) + 2*(Gt211*gt22L + 
      gt23L*Gt311 + gt13L*Gt312)) + Gt112*(gt11L*Gt122 + gt12L*Gt222 + 
      2*(Gt212*gt22L + gt23L*Gt312) + gt13L*Gt322) + Gt322*(Gt112*gt13L + 
      4*Gt212*gt23L + 2*(Gt213*gt22L + gt23L*Gt313 + Gt312*gt33L)) + 
      2*(Gt312*(Gt123*gt12L + Gt223*gt22L + gt23L*Gt323) + gt12L*(Gt113*Gt322 
      + SQR(Gt112)))) + gtu11*(Gt112*(gt12L*(2*Gt111 + 4*Gt212) + 
      2*(Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)) + Gt312*(2*(Gt113*gt12L + 
      Gt213*gt22L) + gt23L*(4*Gt212 + 2*Gt313)) + gt11L*SQR(Gt112) + 
      3*gt22L*SQR(Gt212) + gt33L*SQR(Gt312)) + gtu22*(Gt122*(gt12L*(2*Gt112 + 
      4*Gt222) + 2*(Gt212*gt22L + gt23L*Gt312 + gt13L*Gt322)) + 
      Gt322*(2*(Gt123*gt12L + Gt223*gt22L) + gt23L*(4*Gt222 + 2*Gt323)) + 
      gt11L*SQR(Gt122) + 3*gt22L*SQR(Gt222) + gt33L*SQR(Gt322)) + 
      gtu33*(Gt123*(gt12L*(2*Gt113 + 4*Gt223) + 2*(Gt213*gt22L + gt23L*Gt313 
      + gt13L*Gt323)) + Gt323*(2*(gt12L*Gt133 + gt22L*Gt233) + gt23L*(4*Gt223 
      + 2*Gt333)) + gt11L*SQR(Gt123) + 3*gt22L*SQR(Gt223) + gt33L*SQR(Gt323)) 
      + gtu32*(gt22L*(2*(Gt122*Gt213 + Gt233*Gt322) + Gt223*(6*Gt222 + 
      2*Gt323)) + 4*(gt12L*(Gt123*Gt222 + Gt122*Gt223) + gt23L*(Gt223*Gt322 + 
      Gt222*Gt323)) + 2*(Gt123*(Gt112*gt12L + Gt212*gt22L + gt23L*Gt312 + 
      gt13L*Gt322) + gt12L*(Gt133*Gt322 + Gt123*Gt323) + Gt122*(gt11L*Gt123 + 
      Gt113*gt12L + gt23L*Gt313 + gt13L*Gt323) + Gt322*(gt23L*Gt333 + 
      Gt323*gt33L) + gt23L*SQR(Gt323))) + 
      khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt22 + 
      J12L*J21L*PDstandardNth12gt22 + J11L*J22L*PDstandardNth12gt22 + 
      J12L*J31L*PDstandardNth13gt22 + J11L*J32L*PDstandardNth13gt22 + 
      dJ112L*PDstandardNth1gt22 + J21L*J22L*PDstandardNth22gt22 + 
      J22L*J31L*PDstandardNth23gt22 + J21L*J32L*PDstandardNth23gt22 + 
      dJ212L*PDstandardNth2gt22 + J31L*J32L*PDstandardNth33gt22 + 
      dJ312L*PDstandardNth3gt22) + gtu31*(J11L*J13L*PDstandardNth11gt22 + 
      J13L*J21L*PDstandardNth12gt22 + J11L*J23L*PDstandardNth12gt22 + 
      J13L*J31L*PDstandardNth13gt22 + J11L*J33L*PDstandardNth13gt22 + 
      dJ113L*PDstandardNth1gt22 + J21L*J23L*PDstandardNth22gt22 + 
      J23L*J31L*PDstandardNth23gt22 + J21L*J33L*PDstandardNth23gt22 + 
      dJ213L*PDstandardNth2gt22 + J31L*J33L*PDstandardNth33gt22 + 
      dJ313L*PDstandardNth3gt22) + gtu32*(J12L*J13L*PDstandardNth11gt22 + 
      J13L*J22L*PDstandardNth12gt22 + J12L*J23L*PDstandardNth12gt22 + 
      J13L*J32L*PDstandardNth13gt22 + J12L*J33L*PDstandardNth13gt22 + 
      dJ123L*PDstandardNth1gt22 + J22L*J23L*PDstandardNth22gt22 + 
      J23L*J32L*PDstandardNth23gt22 + J22L*J33L*PDstandardNth23gt22 + 
      dJ223L*PDstandardNth2gt22 + J32L*J33L*PDstandardNth33gt22 + 
      dJ323L*PDstandardNth3gt22)) - gtu11*(2*J11L*J21L*PDstandardNth12gt22 + 
      2*J11L*J31L*PDstandardNth13gt22 + dJ111L*PDstandardNth1gt22 + 
      2*J21L*J31L*PDstandardNth23gt22 + dJ211L*PDstandardNth2gt22 + 
      dJ311L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J11L) + 
      PDstandardNth22gt22*SQR(J21L) + PDstandardNth33gt22*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt22 + 
      2*J12L*J32L*PDstandardNth13gt22 + dJ122L*PDstandardNth1gt22 + 
      2*J22L*J32L*PDstandardNth23gt22 + dJ222L*PDstandardNth2gt22 + 
      dJ322L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J12L) + 
      PDstandardNth22gt22*SQR(J22L) + PDstandardNth33gt22*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt22 + 
      2*J13L*J33L*PDstandardNth13gt22 + dJ133L*PDstandardNth1gt22 + 
      2*J23L*J33L*PDstandardNth23gt22 + dJ233L*PDstandardNth2gt22 + 
      dJ333L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J13L) + 
      PDstandardNth22gt22*SQR(J23L) + PDstandardNth33gt22*SQR(J33L)));
    
    CCTK_REAL Rt23 = khalf*((gt13L*J12L + gt12L*J13L)*PDstandardNth1Xt1 + 
      (gt23L*J12L + gt22L*J13L)*PDstandardNth1Xt2 + (gt33L*J12L + 
      gt23L*J13L)*PDstandardNth1Xt3 + (gt13L*J22L + 
      gt12L*J23L)*PDstandardNth2Xt1 + (gt23L*J22L + 
      gt22L*J23L)*PDstandardNth2Xt2 + (gt33L*J22L + 
      gt23L*J23L)*PDstandardNth2Xt3 - 2*(gtu21*(J11L*J12L*PDstandardNth11gt23 
      + J12L*J21L*PDstandardNth12gt23 + J11L*J22L*PDstandardNth12gt23 + 
      J12L*J31L*PDstandardNth13gt23 + J11L*J32L*PDstandardNth13gt23 + 
      dJ112L*PDstandardNth1gt23 + J21L*J22L*PDstandardNth22gt23 + 
      J22L*J31L*PDstandardNth23gt23 + J21L*J32L*PDstandardNth23gt23 + 
      dJ212L*PDstandardNth2gt23 + J31L*J32L*PDstandardNth33gt23 + 
      dJ312L*PDstandardNth3gt23) + gtu31*(J11L*J13L*PDstandardNth11gt23 + 
      J13L*J21L*PDstandardNth12gt23 + J11L*J23L*PDstandardNth12gt23 + 
      J13L*J31L*PDstandardNth13gt23 + J11L*J33L*PDstandardNth13gt23 + 
      dJ113L*PDstandardNth1gt23 + J21L*J23L*PDstandardNth22gt23 + 
      J23L*J31L*PDstandardNth23gt23 + J21L*J33L*PDstandardNth23gt23 + 
      dJ213L*PDstandardNth2gt23 + J31L*J33L*PDstandardNth33gt23 + 
      dJ313L*PDstandardNth3gt23) + gtu32*(J12L*J13L*PDstandardNth11gt23 + 
      J13L*J22L*PDstandardNth12gt23 + J12L*J23L*PDstandardNth12gt23 + 
      J13L*J32L*PDstandardNth13gt23 + J12L*J33L*PDstandardNth13gt23 + 
      dJ123L*PDstandardNth1gt23 + J22L*J23L*PDstandardNth22gt23 + 
      J23L*J32L*PDstandardNth23gt23 + J22L*J33L*PDstandardNth23gt23 + 
      dJ223L*PDstandardNth2gt23 + J32L*J33L*PDstandardNth33gt23 + 
      dJ323L*PDstandardNth3gt23)) + (gt13L*J32L + 
      gt12L*J33L)*PDstandardNth3Xt1 + (gt23L*J32L + 
      gt22L*J33L)*PDstandardNth3Xt2 + (gt33L*J32L + 
      gt23L*J33L)*PDstandardNth3Xt3 + (Gt113*gt12L + Gt112*gt13L + 
      Gt213*gt22L + Gt212*gt23L + gt23L*Gt313 + Gt312*gt33L)*Xt1L + 
      (Gt123*gt12L + Gt122*gt13L + Gt223*gt22L + Gt222*gt23L + gt23L*Gt323 + 
      Gt322*gt33L)*Xt2L + (gt12L*Gt133 + Gt123*gt13L + gt22L*Gt233 + 
      Gt223*gt23L + gt23L*Gt333 + Gt323*gt33L)*Xt3L + 2*((Gt112*gt11L*Gt123 + 
      Gt111*Gt123*gt12L + Gt111*Gt122*gt13L + Gt123*gt12L*Gt212 + 
      Gt112*gt13L*Gt222 + 2*Gt112*gt12L*Gt223 + Gt123*Gt211*gt22L + 
      2*Gt212*Gt223*gt22L + Gt122*Gt211*gt23L + Gt212*Gt222*gt23L + 
      Gt123*gt23L*Gt311 + Gt123*gt13L*Gt312 + 2*Gt223*gt23L*Gt312 + 
      Gt113*gt13L*Gt322 + Gt213*gt23L*Gt322 + Gt113*gt12L*Gt323 + 
      Gt112*gt13L*Gt323 + Gt213*gt22L*Gt323 + Gt212*gt23L*Gt323 + 
      gt23L*Gt313*Gt323 + Gt122*Gt311*gt33L + Gt222*Gt312*gt33L + 
      Gt313*Gt322*gt33L + Gt312*Gt323*gt33L)*gtu21 + (Gt112*gt11L*Gt133 + 
      Gt111*gt12L*Gt133 + Gt111*Gt123*gt13L + gt12L*Gt133*Gt212 + 
      Gt112*gt13L*Gt223 + Gt133*Gt211*gt22L + 2*Gt112*gt12L*Gt233 + 
      2*Gt212*gt22L*Gt233 + Gt123*Gt211*gt23L + Gt212*Gt223*gt23L + 
      Gt133*gt23L*Gt311 + Gt133*gt13L*Gt312 + 2*Gt233*gt23L*Gt312 + 
      Gt113*gt13L*Gt323 + Gt213*gt23L*Gt323 + Gt113*gt12L*Gt333 + 
      Gt112*gt13L*Gt333 + Gt213*gt22L*Gt333 + Gt212*gt23L*Gt333 + 
      gt23L*Gt313*Gt333 + Gt123*Gt311*gt33L + Gt223*Gt312*gt33L + 
      Gt313*Gt323*gt33L + Gt312*Gt333*gt33L)*gtu31 + gtu21*(Gt113*gt11L*Gt122 
      + Gt122*gt13L*Gt212 + 2*Gt122*gt12L*Gt213 + Gt113*gt12L*Gt222 + 
      Gt113*Gt212*gt22L + 2*Gt213*Gt222*gt22L + Gt212*Gt222*gt23L + 
      Gt123*gt13L*Gt312 + Gt113*gt23L*Gt312 + Gt223*gt23L*Gt312 + 
      Gt123*gt12L*Gt313 + Gt122*gt13L*Gt313 + Gt223*gt22L*Gt313 + 
      Gt222*gt23L*Gt313 + Gt113*gt13L*Gt322 + 2*Gt213*gt23L*Gt322 + 
      gt23L*Gt313*Gt323 + Gt212*Gt322*gt33L + Gt313*Gt322*gt33L + 
      Gt312*Gt323*gt33L + Gt112*(Gt113*gt12L + Gt212*gt23L + Gt312*gt33L) + 
      gt13L*SQR(Gt112)) + gtu31*(2*Gt213*Gt223*gt22L + Gt112*Gt213*gt23L + 
      Gt212*Gt223*gt23L + Gt133*gt13L*Gt312 + Gt233*gt23L*Gt312 + 
      gt12L*Gt133*Gt313 + gt22L*Gt233*Gt313 + Gt223*gt23L*Gt313 + 
      Gt123*(2*gt12L*Gt213 + gt13L*(Gt212 + Gt313)) + 2*Gt213*gt23L*Gt323 + 
      Gt113*(gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt213*gt22L + 
      gt23L*Gt313 + gt13L*Gt323) + gt23L*Gt313*Gt333 + Gt112*Gt313*gt33L + 
      Gt212*Gt323*gt33L + Gt313*Gt323*gt33L + Gt312*Gt333*gt33L + 
      gt12L*SQR(Gt113)) + gtu11*(Gt112*Gt113*gt11L + Gt111*Gt113*gt12L + 
      Gt111*Gt112*gt13L + Gt113*gt12L*Gt212 + Gt112*gt13L*Gt212 + 
      2*Gt112*gt12L*Gt213 + Gt113*Gt211*gt22L + 2*Gt212*Gt213*gt22L + 
      Gt112*Gt211*gt23L + Gt113*gt23L*Gt311 + 2*Gt113*gt13L*Gt312 + 
      3*Gt213*gt23L*Gt312 + Gt113*gt12L*Gt313 + Gt112*gt13L*Gt313 + 
      Gt213*gt22L*Gt313 + Gt212*gt23L*Gt313 + Gt112*Gt311*gt33L + 
      Gt212*Gt312*gt33L + 2*Gt312*Gt313*gt33L + gt23L*SQR(Gt212) + 
      gt23L*SQR(Gt313)) + gtu22*(gt11L*Gt122*Gt123 + Gt112*Gt123*gt12L + 
      Gt112*Gt122*gt13L + Gt123*gt12L*Gt222 + Gt122*gt13L*Gt222 + 
      2*Gt122*gt12L*Gt223 + Gt123*Gt212*gt22L + 2*Gt222*Gt223*gt22L + 
      Gt122*Gt212*gt23L + Gt123*gt23L*Gt312 + 2*Gt123*gt13L*Gt322 + 
      3*Gt223*gt23L*Gt322 + Gt123*gt12L*Gt323 + Gt122*gt13L*Gt323 + 
      Gt223*gt22L*Gt323 + Gt222*gt23L*Gt323 + Gt122*Gt312*gt33L + 
      Gt222*Gt322*gt33L + 2*Gt322*Gt323*gt33L + gt23L*SQR(Gt222) + 
      gt23L*SQR(Gt323)) + gtu32*(gt11L*Gt122*Gt133 + Gt112*gt12L*Gt133 + 
      Gt112*Gt123*gt13L + gt12L*Gt133*Gt222 + Gt122*gt13L*Gt223 + 
      Gt133*Gt212*gt22L + 2*Gt122*gt12L*Gt233 + 2*Gt222*gt22L*Gt233 + 
      Gt123*Gt212*gt23L + Gt222*Gt223*gt23L + Gt133*gt23L*Gt312 + 
      Gt133*gt13L*Gt322 + 2*Gt233*gt23L*Gt322 + Gt123*gt13L*Gt323 + 
      Gt223*gt23L*Gt323 + Gt123*gt12L*Gt333 + Gt122*gt13L*Gt333 + 
      Gt223*gt22L*Gt333 + Gt222*gt23L*Gt333 + gt23L*Gt323*Gt333 + 
      Gt123*Gt312*gt33L + Gt223*Gt322*gt33L + Gt322*Gt333*gt33L + 
      gt33L*SQR(Gt323)) + gtu32*(Gt113*Gt123*gt12L + Gt113*Gt122*gt13L + 
      Gt123*gt13L*Gt222 + 3*Gt123*gt12L*Gt223 + Gt123*Gt213*gt22L + 
      Gt122*Gt213*gt23L + Gt222*Gt223*gt23L + Gt123*gt23L*Gt313 + 
      Gt133*gt13L*Gt322 + Gt233*gt23L*Gt322 + gt12L*Gt133*Gt323 + 
      2*Gt123*gt13L*Gt323 + gt22L*Gt233*Gt323 + 3*Gt223*gt23L*Gt323 + 
      gt23L*Gt323*Gt333 + Gt122*Gt313*gt33L + Gt222*Gt323*gt33L + 
      Gt322*Gt333*gt33L + gt11L*SQR(Gt123) + 2*gt22L*SQR(Gt223) + 
      gt33L*SQR(Gt323)) + gtu33*(gt11L*Gt123*Gt133 + Gt113*gt12L*Gt133 + 
      Gt113*Gt123*gt13L + gt12L*Gt133*Gt223 + Gt123*gt13L*Gt223 + 
      Gt133*Gt213*gt22L + 2*Gt123*gt12L*Gt233 + 2*Gt223*gt22L*Gt233 + 
      Gt123*Gt213*gt23L + Gt133*gt23L*Gt313 + 2*Gt133*gt13L*Gt323 + 
      3*Gt233*gt23L*Gt323 + gt12L*Gt133*Gt333 + Gt123*gt13L*Gt333 + 
      gt22L*Gt233*Gt333 + Gt223*gt23L*Gt333 + Gt123*Gt313*gt33L + 
      Gt223*Gt323*gt33L + 2*Gt323*Gt333*gt33L + gt23L*SQR(Gt223) + 
      gt23L*SQR(Gt333))) - gtu11*(2*J11L*J21L*PDstandardNth12gt23 + 
      2*J11L*J31L*PDstandardNth13gt23 + dJ111L*PDstandardNth1gt23 + 
      2*J21L*J31L*PDstandardNth23gt23 + dJ211L*PDstandardNth2gt23 + 
      dJ311L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J11L) + 
      PDstandardNth22gt23*SQR(J21L) + PDstandardNth33gt23*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt23 + 
      2*J12L*J32L*PDstandardNth13gt23 + dJ122L*PDstandardNth1gt23 + 
      2*J22L*J32L*PDstandardNth23gt23 + dJ222L*PDstandardNth2gt23 + 
      dJ322L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J12L) + 
      PDstandardNth22gt23*SQR(J22L) + PDstandardNth33gt23*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt23 + 
      2*J13L*J33L*PDstandardNth13gt23 + dJ133L*PDstandardNth1gt23 + 
      2*J23L*J33L*PDstandardNth23gt23 + dJ233L*PDstandardNth2gt23 + 
      dJ333L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J13L) + 
      PDstandardNth22gt23*SQR(J23L) + PDstandardNth33gt23*SQR(J33L)));
    
    CCTK_REAL Rt33 = (4*((Gt123*gt13L + Gt223*gt23L)*Gt313 + (Gt113*gt13L 
      + Gt213*gt23L)*Gt323) + (2*Gt213*Gt322 + 6*Gt313*Gt323)*gt33L + 
      2*(gt13L*(Gt122*Gt213 + Gt112*Gt223) + Gt213*(Gt223*gt22L + 
      Gt222*gt23L) + Gt123*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L + 
      Gt311*gt33L) + Gt223*(Gt212*gt23L + Gt312*gt33L) + Gt113*(gt11L*Gt123 + 
      Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + Gt312*gt33L)))*gtu21 + 
      J13L*(gt13L*PDstandardNth1Xt1 + gt23L*PDstandardNth1Xt2 + 
      gt33L*PDstandardNth1Xt3) + J23L*(gt13L*PDstandardNth2Xt1 + 
      gt23L*PDstandardNth2Xt2 + gt33L*PDstandardNth2Xt3) + 
      J33L*(gt13L*PDstandardNth3Xt1 + gt23L*PDstandardNth3Xt2 + 
      gt33L*PDstandardNth3Xt3) + (Gt113*gt13L + Gt213*gt23L + 
      Gt313*gt33L)*Xt1L + (Gt123*gt13L + Gt223*gt23L + Gt323*gt33L)*Xt2L + 
      (Gt133*gt13L + Gt233*gt23L + Gt333*gt33L)*Xt3L + 
      gtu31*(Gt133*(Gt113*gt11L + 2*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L) 
      + 4*gt13L*Gt313) + Gt333*(3*Gt113*gt13L + 4*Gt213*gt23L + 
      6*Gt313*gt33L) + Gt233*(Gt113*gt12L + 4*gt23L*Gt313 + 2*(Gt213*gt22L + 
      Gt212*gt23L + Gt312*gt33L)) + Gt113*(gt11L*Gt133 + gt12L*Gt233 + 
      gt13L*Gt333 + 2*(Gt213*gt23L + Gt313*gt33L)) + 2*(Gt133*Gt311*gt33L + 
      Gt213*(Gt223*gt23L + Gt323*gt33L) + gt13L*(Gt123*Gt213 + Gt112*Gt233 + 
      SQR(Gt113)))) + gtu32*(4*((Gt133*gt13L + Gt233*gt23L)*Gt323 + 
      (Gt123*gt13L + Gt223*gt23L)*Gt333) + Gt323*(2*Gt223 + 6*Gt333)*gt33L + 
      2*(Gt133*(Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + Gt312*gt33L) + 
      Gt123*(gt11L*Gt133 + gt13L*(Gt113 + Gt223) + gt12L*Gt233 + Gt213*gt23L 
      + Gt313*gt33L) + Gt233*(Gt122*gt13L + Gt223*gt22L + Gt222*gt23L + 
      Gt322*gt33L) + gt23L*SQR(Gt223))) + gtu11*(4*(Gt113*gt13L + 
      Gt213*gt23L)*Gt313 + 2*(Gt113*(Gt111*gt13L + gt12L*Gt213 + Gt211*gt23L 
      + Gt311*gt33L) + Gt213*(Gt112*gt13L + Gt212*gt23L + Gt312*gt33L)) + 
      gt11L*SQR(Gt113) + gt22L*SQR(Gt213) + 3*gt33L*SQR(Gt313)) + 
      gtu22*(4*(Gt123*gt13L + Gt223*gt23L)*Gt323 + 2*(Gt123*(Gt112*gt13L + 
      gt12L*Gt223 + Gt212*gt23L + Gt312*gt33L) + Gt223*(Gt122*gt13L + 
      Gt222*gt23L + Gt322*gt33L)) + gt11L*SQR(Gt123) + gt22L*SQR(Gt223) + 
      3*gt33L*SQR(Gt323)) + gtu33*(4*(Gt133*gt13L + Gt233*gt23L)*Gt333 + 
      2*(Gt133*(Gt113*gt13L + gt12L*Gt233 + Gt213*gt23L + Gt313*gt33L) + 
      Gt233*(Gt123*gt13L + Gt223*gt23L + Gt323*gt33L)) + gt11L*SQR(Gt133) + 
      gt22L*SQR(Gt233) + 3*gt33L*SQR(Gt333)) + 
      khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt33 + 
      J12L*J21L*PDstandardNth12gt33 + J11L*J22L*PDstandardNth12gt33 + 
      J12L*J31L*PDstandardNth13gt33 + J11L*J32L*PDstandardNth13gt33 + 
      dJ112L*PDstandardNth1gt33 + J21L*J22L*PDstandardNth22gt33 + 
      J22L*J31L*PDstandardNth23gt33 + J21L*J32L*PDstandardNth23gt33 + 
      dJ212L*PDstandardNth2gt33 + J31L*J32L*PDstandardNth33gt33 + 
      dJ312L*PDstandardNth3gt33) + gtu31*(J11L*J13L*PDstandardNth11gt33 + 
      J13L*J21L*PDstandardNth12gt33 + J11L*J23L*PDstandardNth12gt33 + 
      J13L*J31L*PDstandardNth13gt33 + J11L*J33L*PDstandardNth13gt33 + 
      dJ113L*PDstandardNth1gt33 + J21L*J23L*PDstandardNth22gt33 + 
      J23L*J31L*PDstandardNth23gt33 + J21L*J33L*PDstandardNth23gt33 + 
      dJ213L*PDstandardNth2gt33 + J31L*J33L*PDstandardNth33gt33 + 
      dJ313L*PDstandardNth3gt33) + gtu32*(J12L*J13L*PDstandardNth11gt33 + 
      J13L*J22L*PDstandardNth12gt33 + J12L*J23L*PDstandardNth12gt33 + 
      J13L*J32L*PDstandardNth13gt33 + J12L*J33L*PDstandardNth13gt33 + 
      dJ123L*PDstandardNth1gt33 + J22L*J23L*PDstandardNth22gt33 + 
      J23L*J32L*PDstandardNth23gt33 + J22L*J33L*PDstandardNth23gt33 + 
      dJ223L*PDstandardNth2gt33 + J32L*J33L*PDstandardNth33gt33 + 
      dJ323L*PDstandardNth3gt33)) - gtu11*(2*J11L*J21L*PDstandardNth12gt33 + 
      2*J11L*J31L*PDstandardNth13gt33 + dJ111L*PDstandardNth1gt33 + 
      2*J21L*J31L*PDstandardNth23gt33 + dJ211L*PDstandardNth2gt33 + 
      dJ311L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J11L) + 
      PDstandardNth22gt33*SQR(J21L) + PDstandardNth33gt33*SQR(J31L)) - 
      gtu22*(2*J12L*J22L*PDstandardNth12gt33 + 
      2*J12L*J32L*PDstandardNth13gt33 + dJ122L*PDstandardNth1gt33 + 
      2*J22L*J32L*PDstandardNth23gt33 + dJ222L*PDstandardNth2gt33 + 
      dJ322L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J12L) + 
      PDstandardNth22gt33*SQR(J22L) + PDstandardNth33gt33*SQR(J32L)) - 
      gtu33*(2*J13L*J23L*PDstandardNth12gt33 + 
      2*J13L*J33L*PDstandardNth13gt33 + dJ133L*PDstandardNth1gt33 + 
      2*J23L*J33L*PDstandardNth23gt33 + dJ233L*PDstandardNth2gt33 + 
      dJ333L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J13L) + 
      PDstandardNth22gt33*SQR(J23L) + PDstandardNth33gt33*SQR(J33L)));
    
    CCTK_REAL fac1 = IfThen(conformalMethod,-(khalf*INV(phiL)),1);
    
    CCTK_REAL cdphi1 = fac1*(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + J31L*PDstandardNth3phi);
    
    CCTK_REAL cdphi2 = fac1*(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + J32L*PDstandardNth3phi);
    
    CCTK_REAL cdphi3 = fac1*(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi);
    
    CCTK_REAL fac2 = IfThen(conformalMethod,khalf*pow(phiL,-2),0);
    
    CCTK_REAL cdphi211 = fac1*((dJ111L - Gt111*J11L - Gt211*J12L - 
      Gt311*J13L)*PDstandardNth1phi + 2*(J11L*(J21L*PDstandardNth12phi + 
      J31L*PDstandardNth13phi) + J21L*J31L*PDstandardNth23phi) + (dJ211L - 
      Gt111*J21L - Gt211*J22L - Gt311*J23L)*PDstandardNth2phi + (dJ311L - 
      Gt111*J31L - Gt211*J32L - Gt311*J33L)*PDstandardNth3phi + 
      PDstandardNth11phi*SQR(J11L) + PDstandardNth22phi*SQR(J21L) + 
      PDstandardNth33phi*SQR(J31L)) + fac2*SQR(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + J31L*PDstandardNth3phi);
    
    CCTK_REAL cdphi212 = fac2*(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + 
      J31L*PDstandardNth3phi)*(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + J32L*PDstandardNth3phi) + 
      fac1*(J12L*(J11L*PDstandardNth11phi + J21L*PDstandardNth12phi + 
      J31L*PDstandardNth13phi) + J11L*(J22L*PDstandardNth12phi + 
      J32L*PDstandardNth13phi) + (dJ112L - Gt112*J11L - Gt212*J12L - 
      Gt312*J13L)*PDstandardNth1phi + J22L*(J21L*PDstandardNth22phi + 
      J31L*PDstandardNth23phi) + (dJ212L - Gt112*J21L - Gt212*J22L - 
      Gt312*J23L)*PDstandardNth2phi + J32L*(J21L*PDstandardNth23phi + 
      J31L*PDstandardNth33phi) + (dJ312L - Gt112*J31L - Gt212*J32L - 
      Gt312*J33L)*PDstandardNth3phi);
    
    CCTK_REAL cdphi213 = fac2*(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + 
      J31L*PDstandardNth3phi)*(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi) + 
      fac1*(J13L*(J11L*PDstandardNth11phi + J21L*PDstandardNth12phi + 
      J31L*PDstandardNth13phi) + J11L*(J23L*PDstandardNth12phi + 
      J33L*PDstandardNth13phi) + (dJ113L - Gt113*J11L - Gt213*J12L - 
      Gt313*J13L)*PDstandardNth1phi + J23L*(J21L*PDstandardNth22phi + 
      J31L*PDstandardNth23phi) + (dJ213L - Gt113*J21L - Gt213*J22L - 
      Gt313*J23L)*PDstandardNth2phi + J33L*(J21L*PDstandardNth23phi + 
      J31L*PDstandardNth33phi) + (dJ313L - Gt113*J31L - Gt213*J32L - 
      Gt313*J33L)*PDstandardNth3phi);
    
    CCTK_REAL cdphi222 = fac1*((dJ122L - Gt122*J11L - Gt222*J12L - 
      Gt322*J13L)*PDstandardNth1phi + 2*(J12L*(J22L*PDstandardNth12phi + 
      J32L*PDstandardNth13phi) + J22L*J32L*PDstandardNth23phi) + (dJ222L - 
      Gt122*J21L - Gt222*J22L - Gt322*J23L)*PDstandardNth2phi + (dJ322L - 
      Gt122*J31L - Gt222*J32L - Gt322*J33L)*PDstandardNth3phi + 
      PDstandardNth11phi*SQR(J12L) + PDstandardNth22phi*SQR(J22L) + 
      PDstandardNth33phi*SQR(J32L)) + fac2*SQR(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + J32L*PDstandardNth3phi);
    
    CCTK_REAL cdphi223 = fac2*(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + 
      J32L*PDstandardNth3phi)*(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi) + 
      fac1*(J13L*(J12L*PDstandardNth11phi + J22L*PDstandardNth12phi + 
      J32L*PDstandardNth13phi) + J12L*(J23L*PDstandardNth12phi + 
      J33L*PDstandardNth13phi) + (dJ123L - Gt123*J11L - Gt223*J12L - 
      Gt323*J13L)*PDstandardNth1phi + J23L*(J22L*PDstandardNth22phi + 
      J32L*PDstandardNth23phi) + (dJ223L - Gt123*J21L - Gt223*J22L - 
      Gt323*J23L)*PDstandardNth2phi + J33L*(J22L*PDstandardNth23phi + 
      J32L*PDstandardNth33phi) + (dJ323L - Gt123*J31L - Gt223*J32L - 
      Gt323*J33L)*PDstandardNth3phi);
    
    CCTK_REAL cdphi233 = fac1*((dJ133L - Gt133*J11L - Gt233*J12L - 
      Gt333*J13L)*PDstandardNth1phi + 2*(J13L*(J23L*PDstandardNth12phi + 
      J33L*PDstandardNth13phi) + J23L*J33L*PDstandardNth23phi) + (dJ233L - 
      Gt133*J21L - Gt233*J22L - Gt333*J23L)*PDstandardNth2phi + (dJ333L - 
      Gt133*J31L - Gt233*J32L - Gt333*J33L)*PDstandardNth3phi + 
      PDstandardNth11phi*SQR(J13L) + PDstandardNth22phi*SQR(J23L) + 
      PDstandardNth33phi*SQR(J33L)) + fac2*SQR(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi);
    
    CCTK_REAL Rphi11 = -2*(cdphi211 + 2*(-1 + gt11L*gtu11)*SQR(cdphi1) + 
      gt11L*(cdphi211*gtu11 + 4*(cdphi1*(cdphi2*gtu21 + cdphi3*gtu31) + 
      cdphi2*cdphi3*gtu32) + cdphi233*gtu33 + gtu22*(cdphi222 + 
      2*SQR(cdphi2)) + 2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + 
      gtu33*SQR(cdphi3))));
    
    CCTK_REAL Rphi12 = -2*(cdphi212 + cdphi1*(cdphi2*(-2 + 4*gt12L*gtu21) 
      + 4*cdphi3*gt12L*gtu31) + gt12L*(cdphi211*gtu11 + 4*cdphi2*cdphi3*gtu32 
      + 2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + 
      gtu11*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2)) + gtu33*(cdphi233 
      + 2*SQR(cdphi3))));
    
    CCTK_REAL Rphi13 = -2*(cdphi213 + cdphi1*(4*cdphi2*gt13L*gtu21 + 
      cdphi3*(-2 + 4*gt13L*gtu31)) + gt13L*(cdphi211*gtu11 + 
      4*cdphi2*cdphi3*gtu32 + 2*(cdphi212*gtu21 + cdphi213*gtu31 + 
      cdphi223*gtu32 + gtu11*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2)) 
      + gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL Rphi22 = -2*(cdphi222 + 2*(-1 + gt22L*gtu22)*SQR(cdphi2) + 
      gt22L*(cdphi222*gtu22 + 4*(cdphi1*cdphi3*gtu31 + cdphi2*(cdphi1*gtu21 + 
      cdphi3*gtu32)) + cdphi233*gtu33 + gtu11*(cdphi211 + 2*SQR(cdphi1)) + 
      2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + 
      gtu33*SQR(cdphi3))));
    
    CCTK_REAL Rphi23 = -2*(cdphi223 + cdphi2*(4*cdphi1*gt23L*gtu21 + 
      cdphi3*(-2 + 4*gt23L*gtu32)) + gt23L*(cdphi222*gtu22 + 
      4*cdphi1*cdphi3*gtu31 + gtu11*(cdphi211 + 2*SQR(cdphi1)) + 
      2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + 
      gtu22*SQR(cdphi2)) + gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL Rphi33 = -2*(cdphi233 + gt33L*((4*cdphi1*cdphi2 + 
      2*cdphi212)*gtu21 + 4*cdphi3*(cdphi1*gtu31 + cdphi2*gtu32) + 
      2*(cdphi213*gtu31 + cdphi223*gtu32) + cdphi233*gtu33 + gtu11*(cdphi211 
      + 2*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2))) + 2*(-1 + 
      gt33L*gtu33)*SQR(cdphi3));
    
    CCTK_REAL e4phi = IfThen(conformalMethod,pow(phiL,-2),exp(4*phiL));
    
    CCTK_REAL em4phi = INV(e4phi);
    
    CCTK_REAL gu11 = em4phi*gtu11;
    
    CCTK_REAL gu21 = em4phi*gtu21;
    
    CCTK_REAL gu31 = em4phi*gtu31;
    
    CCTK_REAL gu22 = em4phi*gtu22;
    
    CCTK_REAL gu32 = em4phi*gtu32;
    
    CCTK_REAL gu33 = em4phi*gtu33;
    
    CCTK_REAL R11 = Rphi11 + Rt11;
    
    CCTK_REAL R12 = Rphi12 + Rt12;
    
    CCTK_REAL R13 = Rphi13 + Rt13;
    
    CCTK_REAL R22 = Rphi22 + Rt22;
    
    CCTK_REAL R23 = Rphi23 + Rt23;
    
    CCTK_REAL R33 = Rphi33 + Rt33;
    
    CCTK_REAL trR = gu11*R11 + gu22*R22 + 2*(gu21*R12 + gu31*R13 + 
      gu32*R23) + gu33*R33;
    
    CCTK_REAL Atm11 = At11L*gtu11 + At12L*gtu21 + At13L*gtu31;
    
    CCTK_REAL Atm21 = At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    CCTK_REAL Atm31 = At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    CCTK_REAL Atm12 = At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    CCTK_REAL Atm22 = At12L*gtu21 + At22L*gtu22 + At23L*gtu32;
    
    CCTK_REAL Atm32 = At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    CCTK_REAL Atm13 = At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    CCTK_REAL Atm23 = At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    CCTK_REAL Atm33 = At13L*gtu31 + At23L*gtu32 + At33L*gtu33;
    
    CCTK_REAL rho = pow(alphaL,-2)*(eTttL - 2*(beta2L*eTtyL + 
      beta3L*eTtzL) + 2*(beta1L*(-eTtxL + beta2L*eTxyL + beta3L*eTxzL) + 
      beta2L*beta3L*eTyzL) + eTxxL*SQR(beta1L) + eTyyL*SQR(beta2L) + 
      eTzzL*SQR(beta3L));
    
    CCTK_REAL S1 = (-eTtxL + beta1L*eTxxL + beta2L*eTxyL + 
      beta3L*eTxzL)*INV(alphaL);
    
    CCTK_REAL S2 = (-eTtyL + beta1L*eTxyL + beta2L*eTyyL + 
      beta3L*eTyzL)*INV(alphaL);
    
    CCTK_REAL S3 = (-eTtzL + beta1L*eTxzL + beta2L*eTyzL + 
      beta3L*eTzzL)*INV(alphaL);
    
    CCTK_REAL HL = -2.*(Atm12*Atm21 + Atm13*Atm31 + Atm23*Atm32) - 
      50.26548245743669181540229413247204614715*rho + trR - 1.*(SQR(Atm11) + 
      SQR(Atm22) + SQR(Atm33)) + 
      0.6666666666666666666666666666666666666667*SQR(trKL);
    
    CCTK_REAL M1L = -(At12L*Gt112*gtu22) - At22L*Gt212*gtu22 - 
      At12L*Gt222*gtu22 - At23L*Gt312*gtu22 - At13L*Gt322*gtu22 + 
      At12L*(6*cdphi1*gtu21 - Gt111*gtu21 + 6*cdphi2*gtu22) + 
      At11L*(6*cdphi1*gtu11 - 2*Gt111*gtu11 + 6*cdphi2*gtu21 - Gt122*gtu22) + 
      6*At13L*cdphi1*gtu31 + 6*At11L*cdphi3*gtu31 - At13L*Gt111*gtu31 - 
      3*At11L*Gt113*gtu31 - At23L*Gt211*gtu31 - 3*At12L*Gt213*gtu31 - 
      At33L*Gt311*gtu31 - 3*At13L*Gt313*gtu31 + 6*At13L*cdphi2*gtu32 + 
      6*At12L*cdphi3*gtu32 - At13L*Gt112*gtu32 - At12L*Gt113*gtu32 - 
      At23L*Gt212*gtu32 - At22L*Gt213*gtu32 - 2*At12L*Gt223*gtu32 - 
      At33L*Gt312*gtu32 - At23L*Gt313*gtu32 - 2*At13L*Gt323*gtu32 - 
      2*(At12L*Gt211*gtu11 + At13L*Gt311*gtu11 + At11L*Gt123*gtu32) + 
      6*At13L*cdphi3*gtu33 - At13L*Gt113*gtu33 - At11L*Gt133*gtu33 - 
      At23L*Gt213*gtu33 - At12L*Gt233*gtu33 - At33L*Gt313*gtu33 - 
      At13L*Gt333*gtu33 + (gtu11*J11L + gtu21*J12L + 
      gtu31*J13L)*PDstandardNth1At11 + gtu22*J12L*PDstandardNth1At12 + 
      gtu32*J13L*PDstandardNth1At12 + gtu21*(-3*At11L*Gt112 - At22L*Gt211 - 
      3*At12L*Gt212 - At23L*Gt311 - 3*At13L*Gt312 + J11L*PDstandardNth1At12) 
      + gtu31*J11L*PDstandardNth1At13 + gtu32*J12L*PDstandardNth1At13 + 
      gtu33*J13L*PDstandardNth1At13 - J11L*ktwothird*PDstandardNth1trK + 
      (gtu11*J21L + gtu21*J22L + gtu31*J23L)*PDstandardNth2At11 + 
      gtu21*J21L*PDstandardNth2At12 + gtu22*J22L*PDstandardNth2At12 + 
      gtu32*J23L*PDstandardNth2At12 + gtu31*J21L*PDstandardNth2At13 + 
      gtu32*J22L*PDstandardNth2At13 + gtu33*J23L*PDstandardNth2At13 - 
      J21L*ktwothird*PDstandardNth2trK + gtu11*J31L*PDstandardNth3At11 + 
      gtu21*J32L*PDstandardNth3At11 + gtu31*J33L*PDstandardNth3At11 + 
      gtu21*J31L*PDstandardNth3At12 + gtu22*J32L*PDstandardNth3At12 + 
      gtu32*J33L*PDstandardNth3At12 + gtu31*J31L*PDstandardNth3At13 + 
      gtu32*J32L*PDstandardNth3At13 + gtu33*J33L*PDstandardNth3At13 - 
      J31L*ktwothird*PDstandardNth3trK - 
      25.13274122871834590770114706623602307358*S1;
    
    CCTK_REAL M2L = -(At11L*Gt112*gtu11) - At22L*Gt211*gtu11 - 
      At12L*Gt212*gtu11 - At23L*Gt311*gtu11 - At13L*Gt312*gtu11 + 
      At12L*(6*cdphi1*gtu11 - Gt111*gtu11) + 6*At22L*cdphi1*gtu21 + 
      6*At12L*cdphi2*gtu21 - 3*At12L*Gt112*gtu21 - At11L*Gt122*gtu21 - 
      3*At22L*Gt212*gtu21 - At12L*Gt222*gtu21 - 3*At23L*Gt312*gtu21 - 
      At13L*Gt322*gtu21 + 6*At22L*cdphi2*gtu22 - 2*At12L*Gt122*gtu22 - 
      2*At22L*Gt222*gtu22 - 2*At23L*Gt322*gtu22 + 6*At23L*cdphi1*gtu31 + 
      6*At12L*cdphi3*gtu31 - At13L*Gt112*gtu31 - 2*At12L*Gt113*gtu31 - 
      At11L*Gt123*gtu31 - At23L*Gt212*gtu31 - 2*At22L*Gt213*gtu31 - 
      At12L*Gt223*gtu31 - At33L*Gt312*gtu31 - 2*At23L*Gt313*gtu31 - 
      At13L*Gt323*gtu31 + 6*At23L*cdphi2*gtu32 + 6*At22L*cdphi3*gtu32 - 
      At13L*Gt122*gtu32 - 3*At12L*Gt123*gtu32 - At23L*Gt222*gtu32 - 
      3*At22L*Gt223*gtu32 - At33L*Gt322*gtu32 - 3*At23L*Gt323*gtu32 + 
      6*At23L*cdphi3*gtu33 - At13L*Gt123*gtu33 - At12L*Gt133*gtu33 - 
      At23L*Gt223*gtu33 - At22L*Gt233*gtu33 - At33L*Gt323*gtu33 - 
      At23L*Gt333*gtu33 + (gtu11*J11L + gtu21*J12L + 
      gtu31*J13L)*PDstandardNth1At12 + gtu21*J11L*PDstandardNth1At22 + 
      gtu22*J12L*PDstandardNth1At22 + gtu32*J13L*PDstandardNth1At22 + 
      gtu31*J11L*PDstandardNth1At23 + gtu32*J12L*PDstandardNth1At23 + 
      gtu33*J13L*PDstandardNth1At23 - J12L*ktwothird*PDstandardNth1trK + 
      (gtu11*J21L + gtu21*J22L + gtu31*J23L)*PDstandardNth2At12 + 
      gtu21*J21L*PDstandardNth2At22 + gtu22*J22L*PDstandardNth2At22 + 
      gtu32*J23L*PDstandardNth2At22 + gtu31*J21L*PDstandardNth2At23 + 
      gtu32*J22L*PDstandardNth2At23 + gtu33*J23L*PDstandardNth2At23 - 
      J22L*ktwothird*PDstandardNth2trK + gtu11*J31L*PDstandardNth3At12 + 
      gtu21*J32L*PDstandardNth3At12 + gtu31*J33L*PDstandardNth3At12 + 
      gtu21*J31L*PDstandardNth3At22 + gtu22*J32L*PDstandardNth3At22 + 
      gtu32*J33L*PDstandardNth3At22 + gtu31*J31L*PDstandardNth3At23 + 
      gtu32*J32L*PDstandardNth3At23 + gtu33*J33L*PDstandardNth3At23 - 
      J32L*ktwothird*PDstandardNth3trK - 
      25.13274122871834590770114706623602307358*S2;
    
    CCTK_REAL M3L = -(At11L*Gt113*gtu11) - At23L*Gt211*gtu11 - 
      At12L*Gt213*gtu11 - At33L*Gt311*gtu11 - At13L*Gt313*gtu11 + 
      At13L*(6*cdphi1*gtu11 - Gt111*gtu11) + 6*At23L*cdphi1*gtu21 + 
      6*At13L*cdphi2*gtu21 - 2*At13L*Gt112*gtu21 - At12L*Gt113*gtu21 - 
      At11L*Gt123*gtu21 - 2*At23L*Gt212*gtu21 - At22L*Gt213*gtu21 - 
      At12L*Gt223*gtu21 - 2*At33L*Gt312*gtu21 - At23L*Gt313*gtu21 - 
      At13L*Gt323*gtu21 + 6*At23L*cdphi2*gtu22 - At13L*Gt122*gtu22 - 
      At12L*Gt123*gtu22 - At23L*Gt222*gtu22 - At22L*Gt223*gtu22 - 
      At33L*Gt322*gtu22 - At23L*Gt323*gtu22 + 6*At33L*cdphi1*gtu31 + 
      6*At13L*cdphi3*gtu31 - 3*At13L*Gt113*gtu31 - At11L*Gt133*gtu31 - 
      3*At23L*Gt213*gtu31 - At12L*Gt233*gtu31 - 3*At33L*Gt313*gtu31 - 
      At13L*Gt333*gtu31 + 6*At33L*cdphi2*gtu32 + 6*At23L*cdphi3*gtu32 - 
      3*At13L*Gt123*gtu32 - At12L*Gt133*gtu32 - 3*At23L*Gt223*gtu32 - 
      At22L*Gt233*gtu32 - 3*At33L*Gt323*gtu32 - At23L*Gt333*gtu32 + 
      6*At33L*cdphi3*gtu33 - 2*At13L*Gt133*gtu33 - 2*At23L*Gt233*gtu33 - 
      2*At33L*Gt333*gtu33 + (gtu11*J11L + gtu21*J12L + 
      gtu31*J13L)*PDstandardNth1At13 + gtu21*J11L*PDstandardNth1At23 + 
      gtu22*J12L*PDstandardNth1At23 + gtu32*J13L*PDstandardNth1At23 + 
      gtu31*J11L*PDstandardNth1At33 + gtu32*J12L*PDstandardNth1At33 + 
      gtu33*J13L*PDstandardNth1At33 - J13L*ktwothird*PDstandardNth1trK + 
      (gtu11*J21L + gtu21*J22L + gtu31*J23L)*PDstandardNth2At13 + 
      gtu21*J21L*PDstandardNth2At23 + gtu22*J22L*PDstandardNth2At23 + 
      gtu32*J23L*PDstandardNth2At23 + gtu31*J21L*PDstandardNth2At33 + 
      gtu32*J22L*PDstandardNth2At33 + gtu33*J23L*PDstandardNth2At33 - 
      J23L*ktwothird*PDstandardNth2trK + gtu11*J31L*PDstandardNth3At13 + 
      gtu21*J32L*PDstandardNth3At13 + gtu31*J33L*PDstandardNth3At13 + 
      gtu21*J31L*PDstandardNth3At23 + gtu22*J32L*PDstandardNth3At23 + 
      gtu32*J33L*PDstandardNth3At23 + gtu31*J31L*PDstandardNth3At33 + 
      gtu32*J32L*PDstandardNth3At33 + gtu33*J33L*PDstandardNth3At33 - 
      J33L*ktwothird*PDstandardNth3trK - 
      25.13274122871834590770114706623602307358*S3;
    
    CCTK_REAL cSL = Log(detgt);
    
    CCTK_REAL cXt1L = Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + 
      Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33 - Xt1L;
    
    CCTK_REAL cXt2L = Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + 
      Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33 - Xt2L;
    
    CCTK_REAL cXt3L = Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + 
      Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33 - Xt3L;
    
    CCTK_REAL cAL = At11L*gtu11 + At22L*gtu22 + 2*(At12L*gtu21 + 
      At13L*gtu31 + At23L*gtu32) + At33L*gtu33;
    
    
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
  }
  LC_ENDLOOP3 (ML_BSSN_MP_constraints);
}

void ML_BSSN_MP_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_constraints_Body);
}
