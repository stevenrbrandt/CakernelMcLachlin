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
#include "Vectors.hh"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

static void ML_BSSN_MP_RHS2_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_RHS2_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_RHS2_calc_every != ML_BSSN_MP_RHS2_calc_offset)
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
  LC_LOOP3 (ML_BSSN_MP_RHS2,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL_VEC PDstandardNth1alpha = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2alpha = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3alpha = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth11alpha = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth22alpha = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth33alpha = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth12alpha = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth13alpha = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth23alpha = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth1At11 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth2At11 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth3At11 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti1At11 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm1At11 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti2At11 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm2At11 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti3At11 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm3At11 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth1At12 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth2At12 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth3At12 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti1At12 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm1At12 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti2At12 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm2At12 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti3At12 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm3At12 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth1At13 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth2At13 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth3At13 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti1At13 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm1At13 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti2At13 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm2At13 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti3At13 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm3At13 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth1At22 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth2At22 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth3At22 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti1At22 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm1At22 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti2At22 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm2At22 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti3At22 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm3At22 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth1At23 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth2At23 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth3At23 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti1At23 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm1At23 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti2At23 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm2At23 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti3At23 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm3At23 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth1At33 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth2At33 = INITVALUE;
    // CCTK_REAL_VEC PDdissipationNth3At33 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti1At33 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm1At33 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti2At33 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm2At33 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthAnti3At33 = INITVALUE;
    // CCTK_REAL_VEC PDupwindNthSymm3At33 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1beta1 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2beta1 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3beta1 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1beta2 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2beta2 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3beta2 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1beta3 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2beta3 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3beta3 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1gt11 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2gt11 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3gt11 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth11gt11 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth22gt11 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth33gt11 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth12gt11 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth13gt11 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth23gt11 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1gt12 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2gt12 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3gt12 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth11gt12 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth22gt12 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth33gt12 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth12gt12 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth13gt12 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth23gt12 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1gt13 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2gt13 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3gt13 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth11gt13 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth22gt13 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth33gt13 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth12gt13 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth13gt13 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth23gt13 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1gt22 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2gt22 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3gt22 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth11gt22 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth22gt22 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth33gt22 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth12gt22 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth13gt22 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth23gt22 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1gt23 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2gt23 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3gt23 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth11gt23 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth22gt23 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth33gt23 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth12gt23 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth13gt23 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth23gt23 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1gt33 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2gt33 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3gt33 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth11gt33 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth22gt33 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth33gt33 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth12gt33 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth13gt33 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth23gt33 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1phi = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2phi = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3phi = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth11phi = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth22phi = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth33phi = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth12phi = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth13phi = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth23phi = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1Xt1 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2Xt1 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3Xt1 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1Xt2 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2Xt2 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3Xt2 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth1Xt3 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth2Xt3 = INITVALUE;
    // CCTK_REAL_VEC PDstandardNth3Xt3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL_VEC  alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC  At11L = vec_load(At11[index]);
    CCTK_REAL_VEC  At12L = vec_load(At12[index]);
    CCTK_REAL_VEC  At13L = vec_load(At13[index]);
    CCTK_REAL_VEC  At22L = vec_load(At22[index]);
    CCTK_REAL_VEC  At23L = vec_load(At23[index]);
    CCTK_REAL_VEC  At33L = vec_load(At33[index]);
    CCTK_REAL_VEC  beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC  beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC  beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC  dJ111L = vec_load(dJ111[index]);
    CCTK_REAL_VEC  dJ112L = vec_load(dJ112[index]);
    CCTK_REAL_VEC  dJ113L = vec_load(dJ113[index]);
    CCTK_REAL_VEC  dJ122L = vec_load(dJ122[index]);
    CCTK_REAL_VEC  dJ123L = vec_load(dJ123[index]);
    CCTK_REAL_VEC  dJ133L = vec_load(dJ133[index]);
    CCTK_REAL_VEC  dJ211L = vec_load(dJ211[index]);
    CCTK_REAL_VEC  dJ212L = vec_load(dJ212[index]);
    CCTK_REAL_VEC  dJ213L = vec_load(dJ213[index]);
    CCTK_REAL_VEC  dJ222L = vec_load(dJ222[index]);
    CCTK_REAL_VEC  dJ223L = vec_load(dJ223[index]);
    CCTK_REAL_VEC  dJ233L = vec_load(dJ233[index]);
    CCTK_REAL_VEC  dJ311L = vec_load(dJ311[index]);
    CCTK_REAL_VEC  dJ312L = vec_load(dJ312[index]);
    CCTK_REAL_VEC  dJ313L = vec_load(dJ313[index]);
    CCTK_REAL_VEC  dJ322L = vec_load(dJ322[index]);
    CCTK_REAL_VEC  dJ323L = vec_load(dJ323[index]);
    CCTK_REAL_VEC  dJ333L = vec_load(dJ333[index]);
    CCTK_REAL_VEC  eTxxL = (*stress_energy_state) ? vec_load(eTxx[index]) : 0.0;
    CCTK_REAL_VEC  eTxyL = (*stress_energy_state) ? vec_load(eTxy[index]) : 0.0;
    CCTK_REAL_VEC  eTxzL = (*stress_energy_state) ? vec_load(eTxz[index]) : 0.0;
    CCTK_REAL_VEC  eTyyL = (*stress_energy_state) ? vec_load(eTyy[index]) : 0.0;
    CCTK_REAL_VEC  eTyzL = (*stress_energy_state) ? vec_load(eTyz[index]) : 0.0;
    CCTK_REAL_VEC  eTzzL = (*stress_energy_state) ? vec_load(eTzz[index]) : 0.0;
    CCTK_REAL_VEC  gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC  gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC  gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC  gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC  gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC  gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC  J11L = vec_load(J11[index]);
    CCTK_REAL_VEC  J12L = vec_load(J12[index]);
    CCTK_REAL_VEC  J13L = vec_load(J13[index]);
    CCTK_REAL_VEC  J21L = vec_load(J21[index]);
    CCTK_REAL_VEC  J22L = vec_load(J22[index]);
    CCTK_REAL_VEC  J23L = vec_load(J23[index]);
    CCTK_REAL_VEC  J31L = vec_load(J31[index]);
    CCTK_REAL_VEC  J32L = vec_load(J32[index]);
    CCTK_REAL_VEC  J33L = vec_load(J33[index]);
    CCTK_REAL_VEC  phiL = vec_load(phi[index]);
    CCTK_REAL_VEC  trKL = vec_load(trK[index]);
    CCTK_REAL_VEC  Xt1L = vec_load(Xt1[index]);
    CCTK_REAL_VEC  Xt2L = vec_load(Xt2[index]);
    CCTK_REAL_VEC  Xt3L = vec_load(Xt3[index]);
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC const PDstandardNth1alpha = PDstandardNth1(alpha, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2alpha = PDstandardNth2(alpha, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3alpha = PDstandardNth3(alpha, i, j, k);
    CCTK_REAL_VEC const PDstandardNth11alpha = PDstandardNth11(alpha, i, j, k);
    CCTK_REAL_VEC const PDstandardNth22alpha = PDstandardNth22(alpha, i, j, k);
    CCTK_REAL_VEC const PDstandardNth33alpha = PDstandardNth33(alpha, i, j, k);
    CCTK_REAL_VEC const PDstandardNth12alpha = PDstandardNth12(alpha, i, j, k);
    CCTK_REAL_VEC const PDstandardNth13alpha = PDstandardNth13(alpha, i, j, k);
    CCTK_REAL_VEC const PDstandardNth23alpha = PDstandardNth23(alpha, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth1At11 = PDdissipationNth1(At11, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth2At11 = PDdissipationNth2(At11, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth3At11 = PDdissipationNth3(At11, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti1At11 = PDupwindNthAnti1(At11, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm1At11 = PDupwindNthSymm1(At11, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti2At11 = PDupwindNthAnti2(At11, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm2At11 = PDupwindNthSymm2(At11, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti3At11 = PDupwindNthAnti3(At11, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm3At11 = PDupwindNthSymm3(At11, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth1At12 = PDdissipationNth1(At12, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth2At12 = PDdissipationNth2(At12, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth3At12 = PDdissipationNth3(At12, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti1At12 = PDupwindNthAnti1(At12, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm1At12 = PDupwindNthSymm1(At12, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti2At12 = PDupwindNthAnti2(At12, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm2At12 = PDupwindNthSymm2(At12, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti3At12 = PDupwindNthAnti3(At12, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm3At12 = PDupwindNthSymm3(At12, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth1At13 = PDdissipationNth1(At13, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth2At13 = PDdissipationNth2(At13, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth3At13 = PDdissipationNth3(At13, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti1At13 = PDupwindNthAnti1(At13, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm1At13 = PDupwindNthSymm1(At13, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti2At13 = PDupwindNthAnti2(At13, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm2At13 = PDupwindNthSymm2(At13, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti3At13 = PDupwindNthAnti3(At13, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm3At13 = PDupwindNthSymm3(At13, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth1At22 = PDdissipationNth1(At22, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth2At22 = PDdissipationNth2(At22, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth3At22 = PDdissipationNth3(At22, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti1At22 = PDupwindNthAnti1(At22, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm1At22 = PDupwindNthSymm1(At22, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti2At22 = PDupwindNthAnti2(At22, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm2At22 = PDupwindNthSymm2(At22, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti3At22 = PDupwindNthAnti3(At22, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm3At22 = PDupwindNthSymm3(At22, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth1At23 = PDdissipationNth1(At23, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth2At23 = PDdissipationNth2(At23, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth3At23 = PDdissipationNth3(At23, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti1At23 = PDupwindNthAnti1(At23, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm1At23 = PDupwindNthSymm1(At23, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti2At23 = PDupwindNthAnti2(At23, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm2At23 = PDupwindNthSymm2(At23, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti3At23 = PDupwindNthAnti3(At23, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm3At23 = PDupwindNthSymm3(At23, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth1At33 = PDdissipationNth1(At33, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth2At33 = PDdissipationNth2(At33, i, j, k);
    CCTK_REAL_VEC const PDdissipationNth3At33 = PDdissipationNth3(At33, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti1At33 = PDupwindNthAnti1(At33, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm1At33 = PDupwindNthSymm1(At33, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti2At33 = PDupwindNthAnti2(At33, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm2At33 = PDupwindNthSymm2(At33, i, j, k);
    CCTK_REAL_VEC const PDupwindNthAnti3At33 = PDupwindNthAnti3(At33, i, j, k);
    CCTK_REAL_VEC const PDupwindNthSymm3At33 = PDupwindNthSymm3(At33, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1beta1 = PDstandardNth1(beta1, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2beta1 = PDstandardNth2(beta1, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3beta1 = PDstandardNth3(beta1, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1beta2 = PDstandardNth1(beta2, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2beta2 = PDstandardNth2(beta2, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3beta2 = PDstandardNth3(beta2, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1beta3 = PDstandardNth1(beta3, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2beta3 = PDstandardNth2(beta3, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3beta3 = PDstandardNth3(beta3, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    CCTK_REAL_VEC const PDstandardNth11gt11 = PDstandardNth11(gt11, i, j, k);
    CCTK_REAL_VEC const PDstandardNth22gt11 = PDstandardNth22(gt11, i, j, k);
    CCTK_REAL_VEC const PDstandardNth33gt11 = PDstandardNth33(gt11, i, j, k);
    CCTK_REAL_VEC const PDstandardNth12gt11 = PDstandardNth12(gt11, i, j, k);
    CCTK_REAL_VEC const PDstandardNth13gt11 = PDstandardNth13(gt11, i, j, k);
    CCTK_REAL_VEC const PDstandardNth23gt11 = PDstandardNth23(gt11, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    CCTK_REAL_VEC const PDstandardNth11gt12 = PDstandardNth11(gt12, i, j, k);
    CCTK_REAL_VEC const PDstandardNth22gt12 = PDstandardNth22(gt12, i, j, k);
    CCTK_REAL_VEC const PDstandardNth33gt12 = PDstandardNth33(gt12, i, j, k);
    CCTK_REAL_VEC const PDstandardNth12gt12 = PDstandardNth12(gt12, i, j, k);
    CCTK_REAL_VEC const PDstandardNth13gt12 = PDstandardNth13(gt12, i, j, k);
    CCTK_REAL_VEC const PDstandardNth23gt12 = PDstandardNth23(gt12, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    CCTK_REAL_VEC const PDstandardNth11gt13 = PDstandardNth11(gt13, i, j, k);
    CCTK_REAL_VEC const PDstandardNth22gt13 = PDstandardNth22(gt13, i, j, k);
    CCTK_REAL_VEC const PDstandardNth33gt13 = PDstandardNth33(gt13, i, j, k);
    CCTK_REAL_VEC const PDstandardNth12gt13 = PDstandardNth12(gt13, i, j, k);
    CCTK_REAL_VEC const PDstandardNth13gt13 = PDstandardNth13(gt13, i, j, k);
    CCTK_REAL_VEC const PDstandardNth23gt13 = PDstandardNth23(gt13, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    CCTK_REAL_VEC const PDstandardNth11gt22 = PDstandardNth11(gt22, i, j, k);
    CCTK_REAL_VEC const PDstandardNth22gt22 = PDstandardNth22(gt22, i, j, k);
    CCTK_REAL_VEC const PDstandardNth33gt22 = PDstandardNth33(gt22, i, j, k);
    CCTK_REAL_VEC const PDstandardNth12gt22 = PDstandardNth12(gt22, i, j, k);
    CCTK_REAL_VEC const PDstandardNth13gt22 = PDstandardNth13(gt22, i, j, k);
    CCTK_REAL_VEC const PDstandardNth23gt22 = PDstandardNth23(gt22, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    CCTK_REAL_VEC const PDstandardNth11gt23 = PDstandardNth11(gt23, i, j, k);
    CCTK_REAL_VEC const PDstandardNth22gt23 = PDstandardNth22(gt23, i, j, k);
    CCTK_REAL_VEC const PDstandardNth33gt23 = PDstandardNth33(gt23, i, j, k);
    CCTK_REAL_VEC const PDstandardNth12gt23 = PDstandardNth12(gt23, i, j, k);
    CCTK_REAL_VEC const PDstandardNth13gt23 = PDstandardNth13(gt23, i, j, k);
    CCTK_REAL_VEC const PDstandardNth23gt23 = PDstandardNth23(gt23, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    CCTK_REAL_VEC const PDstandardNth11gt33 = PDstandardNth11(gt33, i, j, k);
    CCTK_REAL_VEC const PDstandardNth22gt33 = PDstandardNth22(gt33, i, j, k);
    CCTK_REAL_VEC const PDstandardNth33gt33 = PDstandardNth33(gt33, i, j, k);
    CCTK_REAL_VEC const PDstandardNth12gt33 = PDstandardNth12(gt33, i, j, k);
    CCTK_REAL_VEC const PDstandardNth13gt33 = PDstandardNth13(gt33, i, j, k);
    CCTK_REAL_VEC const PDstandardNth23gt33 = PDstandardNth23(gt33, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1phi = PDstandardNth1(phi, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2phi = PDstandardNth2(phi, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3phi = PDstandardNth3(phi, i, j, k);
    CCTK_REAL_VEC const PDstandardNth11phi = PDstandardNth11(phi, i, j, k);
    CCTK_REAL_VEC const PDstandardNth22phi = PDstandardNth22(phi, i, j, k);
    CCTK_REAL_VEC const PDstandardNth33phi = PDstandardNth33(phi, i, j, k);
    CCTK_REAL_VEC const PDstandardNth12phi = PDstandardNth12(phi, i, j, k);
    CCTK_REAL_VEC const PDstandardNth13phi = PDstandardNth13(phi, i, j, k);
    CCTK_REAL_VEC const PDstandardNth23phi = PDstandardNth23(phi, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1Xt1 = PDstandardNth1(Xt1, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2Xt1 = PDstandardNth2(Xt1, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3Xt1 = PDstandardNth3(Xt1, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1Xt2 = PDstandardNth1(Xt2, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2Xt2 = PDstandardNth2(Xt2, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3Xt2 = PDstandardNth3(Xt2, i, j, k);
    CCTK_REAL_VEC const PDstandardNth1Xt3 = PDstandardNth1(Xt3, i, j, k);
    CCTK_REAL_VEC const PDstandardNth2Xt3 = PDstandardNth2(Xt3, i, j, k);
    CCTK_REAL_VEC const PDstandardNth3Xt3 = PDstandardNth3(Xt3, i, j, k);
    
    /* Calculate temporaries and grid functions */
    int dir1 = Sign(beta1L);
    
    int dir2 = Sign(beta2L);
    
    int dir3 = Sign(beta3L);
    
    CCTK_REAL_VEC epsdiss1 = EpsDiss;
    
    CCTK_REAL_VEC epsdiss2 = EpsDiss;
    
    CCTK_REAL_VEC epsdiss3 = EpsDiss;
    
    CCTK_REAL_VEC detgt = 1;
    
    CCTK_REAL_VEC gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL_VEC gtu21 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL_VEC gtu31 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL_VEC gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL_VEC gtu32 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL_VEC gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL_VEC Gt111 = khalf*((gtu11*J11L - gtu21*J12L - 
      gtu31*J13L)*PDstandardNth1gt11 + (gtu11*J21L - gtu21*J22L - 
      gtu31*J23L)*PDstandardNth2gt11 + (gtu11*J31L - gtu21*J32L - 
      gtu31*J33L)*PDstandardNth3gt11 + 2*(J11L*(gtu21*PDstandardNth1gt12 + 
      gtu31*PDstandardNth1gt13) + J21L*(gtu21*PDstandardNth2gt12 + 
      gtu31*PDstandardNth2gt13) + J31L*(gtu21*PDstandardNth3gt12 + 
      gtu31*PDstandardNth3gt13)));
    
    CCTK_REAL_VEC Gt211 = khalf*((gtu21*J11L - gtu22*J12L - 
      gtu32*J13L)*PDstandardNth1gt11 + (gtu21*J21L - gtu22*J22L - 
      gtu32*J23L)*PDstandardNth2gt11 + (gtu21*J31L - gtu22*J32L - 
      gtu32*J33L)*PDstandardNth3gt11 + 2*(J11L*(gtu22*PDstandardNth1gt12 + 
      gtu32*PDstandardNth1gt13) + J21L*(gtu22*PDstandardNth2gt12 + 
      gtu32*PDstandardNth2gt13) + J31L*(gtu22*PDstandardNth3gt12 + 
      gtu32*PDstandardNth3gt13)));
    
    CCTK_REAL_VEC Gt311 = khalf*((gtu31*J11L - gtu32*J12L - 
      gtu33*J13L)*PDstandardNth1gt11 + (gtu31*J21L - gtu32*J22L - 
      gtu33*J23L)*PDstandardNth2gt11 + (gtu31*J31L - gtu32*J32L - 
      gtu33*J33L)*PDstandardNth3gt11 + 2*(J11L*(gtu32*PDstandardNth1gt12 + 
      gtu33*PDstandardNth1gt13) + J21L*(gtu32*PDstandardNth2gt12 + 
      gtu33*PDstandardNth2gt13) + J31L*(gtu32*PDstandardNth3gt12 + 
      gtu33*PDstandardNth3gt13)));
    
    CCTK_REAL_VEC Gt112 = khalf*(gtu11*(J12L*PDstandardNth1gt11 + 
      J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
      gtu21*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
      J31L*PDstandardNth3gt22) + gtu31*(-(J13L*PDstandardNth1gt12) + 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23));
    
    CCTK_REAL_VEC Gt212 = khalf*(gtu21*(J12L*PDstandardNth1gt11 + 
      J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
      gtu22*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
      J31L*PDstandardNth3gt22) + gtu32*(-(J13L*PDstandardNth1gt12) + 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23));
    
    CCTK_REAL_VEC Gt312 = khalf*(gtu31*(J12L*PDstandardNth1gt11 + 
      J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
      gtu32*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
      J31L*PDstandardNth3gt22) + gtu33*(-(J13L*PDstandardNth1gt12) + 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23));
    
    CCTK_REAL_VEC Gt113 = khalf*(gtu11*(J13L*PDstandardNth1gt11 + 
      J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
      gtu21*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
      J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
      J31L*PDstandardNth3gt23) + gtu31*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL_VEC Gt213 = khalf*(gtu21*(J13L*PDstandardNth1gt11 + 
      J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
      gtu22*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
      J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
      J31L*PDstandardNth3gt23) + gtu32*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL_VEC Gt313 = khalf*(gtu31*(J13L*PDstandardNth1gt11 + 
      J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
      gtu32*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
      J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
      J31L*PDstandardNth3gt23) + gtu33*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL_VEC Gt122 = khalf*(gtu11*(-(J11L*PDstandardNth1gt22) + 
      2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
      J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
      J31L*PDstandardNth3gt22) + gtu21*(J12L*PDstandardNth1gt22 + 
      J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
      gtu31*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
      J33L*PDstandardNth3gt22 - 2*(J12L*PDstandardNth1gt23 + 
      J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL_VEC Gt222 = khalf*(gtu21*(-(J11L*PDstandardNth1gt22) + 
      2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
      J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
      J31L*PDstandardNth3gt22) + gtu22*(J12L*PDstandardNth1gt22 + 
      J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
      gtu32*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
      J33L*PDstandardNth3gt22 - 2*(J12L*PDstandardNth1gt23 + 
      J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL_VEC Gt322 = khalf*(gtu31*(-(J11L*PDstandardNth1gt22) + 
      2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
      J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
      J31L*PDstandardNth3gt22) + gtu32*(J12L*PDstandardNth1gt22 + 
      J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
      gtu33*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
      J33L*PDstandardNth3gt22 - 2*(J12L*PDstandardNth1gt23 + 
      J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL_VEC Gt123 = khalf*(gtu21*(J13L*PDstandardNth1gt22 + 
      J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
      gtu11*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
      J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt23) + gtu31*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL_VEC Gt223 = khalf*(gtu22*(J13L*PDstandardNth1gt22 + 
      J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
      gtu21*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
      J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt23) + gtu32*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL_VEC Gt323 = khalf*(gtu32*(J13L*PDstandardNth1gt22 + 
      J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
      gtu31*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
      J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt23) + gtu33*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL_VEC Gt133 = khalf*(gtu11*(-(J11L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
      J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt33) + gtu21*(-(J12L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
      J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
      J32L*PDstandardNth3gt33) + gtu31*(J13L*PDstandardNth1gt33 + 
      J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL_VEC Gt233 = khalf*(gtu21*(-(J11L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
      J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt33) + gtu22*(-(J12L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
      J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
      J32L*PDstandardNth3gt33) + gtu32*(J13L*PDstandardNth1gt33 + 
      J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL_VEC Gt333 = khalf*(gtu31*(-(J11L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
      J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt33) + gtu32*(-(J12L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
      J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
      J32L*PDstandardNth3gt33) + gtu33*(J13L*PDstandardNth1gt33 + 
      J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL_VEC Xtn1 = Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + 
      Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    CCTK_REAL_VEC Xtn2 = Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + 
      Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    CCTK_REAL_VEC Xtn3 = Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + 
      Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    CCTK_REAL_VEC Rt11 = (Gt113*(gt13L*Gt312 + 3*(gt12L*Gt212 + 
      gt13L*Gt312)) + gt12L*(Gt213*(4*Gt112 + 2*Gt222) + Gt212*(Gt113 + 
      2*Gt223) + 2*(Gt233*Gt312 + Gt223*Gt313)) + gt11L*(6*Gt112*Gt113 + 
      2*(Gt122*Gt213 + Gt133*Gt312 + Gt123*(Gt212 + Gt313))) + 
      gt13L*(2*Gt213*Gt322 + Gt313*(4*Gt112 + 2*Gt323)) + 
      2*(Gt213*(Gt212*gt22L + gt23L*Gt312) + Gt212*(gt23L*Gt313 + 
      gt13L*Gt323) + Gt312*(gt13L*Gt333 + Gt313*gt33L)))*gtu32 + 
      J11L*(gt11L*PDstandardNth1Xt1 + gt12L*PDstandardNth1Xt2 + 
      gt13L*PDstandardNth1Xt3) + J21L*(gt11L*PDstandardNth2Xt1 + 
      gt12L*PDstandardNth2Xt2 + gt13L*PDstandardNth2Xt3) + 
      J31L*(gt11L*PDstandardNth3Xt1 + gt12L*PDstandardNth3Xt2 + 
      gt13L*PDstandardNth3Xt3) + (Gt111*gt11L + gt12L*Gt211 + 
      gt13L*Gt311)*Xtn1 + (Gt112*gt11L + gt12L*Gt212 + gt13L*Gt312)*Xtn2 + 
      (Gt113*gt11L + gt12L*Gt213 + gt13L*Gt313)*Xtn3 + 
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
    
    CCTK_REAL_VEC Rt12 = khalf*((gt12L*J11L + 
      gt11L*J12L)*PDstandardNth1Xt1 + (gt22L*J11L + 
      gt12L*J12L)*PDstandardNth1Xt2 + (gt23L*J11L + 
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
      gt12L*Gt212 + Gt211*gt22L + gt23L*Gt311 + gt13L*Gt312)*Xtn1 + 
      (gt11L*Gt122 + Gt112*gt12L + gt12L*Gt222 + Gt212*gt22L + gt23L*Gt312 + 
      gt13L*Gt322)*Xtn2 + (gt11L*Gt123 + Gt113*gt12L + gt12L*Gt223 + 
      Gt213*gt22L + gt23L*Gt313 + gt13L*Gt323)*Xtn3 + 2*((Gt123*gt12L*Gt211 + 
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
    
    CCTK_REAL_VEC Rt13 = khalf*((gt13L*J11L + 
      gt11L*J13L)*PDstandardNth1Xt1 + (gt23L*J11L + 
      gt12L*J13L)*PDstandardNth1Xt2 + (gt33L*J11L + 
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
      gt12L*Gt213 + Gt211*gt23L + gt13L*Gt313 + Gt311*gt33L)*Xtn1 + 
      (gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + gt13L*Gt323 + 
      Gt312*gt33L)*Xtn2 + (gt11L*Gt133 + Gt113*gt13L + gt12L*Gt233 + 
      Gt213*gt23L + gt13L*Gt333 + Gt313*gt33L)*Xtn3 + 2*((Gt122*gt13L*Gt211 + 
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
    
    CCTK_REAL_VEC Rt22 = (Gt223*(3*Gt112*gt12L + 6*Gt212*gt22L + 
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
      Gt212*gt22L + gt23L*Gt312)*Xtn1 + (Gt122*gt12L + Gt222*gt22L + 
      gt23L*Gt322)*Xtn2 + (Gt123*gt12L + Gt223*gt22L + gt23L*Gt323)*Xtn3 + 
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
    
    CCTK_REAL_VEC Rt23 = khalf*((gt13L*J12L + 
      gt12L*J13L)*PDstandardNth1Xt1 + (gt23L*J12L + 
      gt22L*J13L)*PDstandardNth1Xt2 + (gt33L*J12L + 
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
      Gt213*gt22L + Gt212*gt23L + gt23L*Gt313 + Gt312*gt33L)*Xtn1 + 
      (Gt123*gt12L + Gt122*gt13L + Gt223*gt22L + Gt222*gt23L + gt23L*Gt323 + 
      Gt322*gt33L)*Xtn2 + (gt12L*Gt133 + Gt123*gt13L + gt22L*Gt233 + 
      Gt223*gt23L + gt23L*Gt333 + Gt323*gt33L)*Xtn3 + 2*((Gt112*gt11L*Gt123 + 
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
    
    CCTK_REAL_VEC Rt33 = (4*((Gt123*gt13L + Gt223*gt23L)*Gt313 + 
      (Gt113*gt13L + Gt213*gt23L)*Gt323) + (2*Gt213*Gt322 + 
      6*Gt313*Gt323)*gt33L + 2*(gt13L*(Gt122*Gt213 + Gt112*Gt223) + 
      Gt213*(Gt223*gt22L + Gt222*gt23L) + Gt123*(Gt111*gt13L + gt12L*Gt213 + 
      Gt211*gt23L + Gt311*gt33L) + Gt223*(Gt212*gt23L + Gt312*gt33L) + 
      Gt113*(gt11L*Gt123 + Gt112*gt13L + gt12L*Gt223 + Gt212*gt23L + 
      Gt312*gt33L)))*gtu21 + J13L*(gt13L*PDstandardNth1Xt1 + 
      gt23L*PDstandardNth1Xt2 + gt33L*PDstandardNth1Xt3) + 
      J23L*(gt13L*PDstandardNth2Xt1 + gt23L*PDstandardNth2Xt2 + 
      gt33L*PDstandardNth2Xt3) + J33L*(gt13L*PDstandardNth3Xt1 + 
      gt23L*PDstandardNth3Xt2 + gt33L*PDstandardNth3Xt3) + (Gt113*gt13L + 
      Gt213*gt23L + Gt313*gt33L)*Xtn1 + (Gt123*gt13L + Gt223*gt23L + 
      Gt323*gt33L)*Xtn2 + (Gt133*gt13L + Gt233*gt23L + Gt333*gt33L)*Xtn3 + 
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
    
    CCTK_REAL_VEC fac1 = IfThen(conformalMethod,-(khalf*INV(phiL)),1);
    
    CCTK_REAL_VEC cdphi1 = fac1*(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + J31L*PDstandardNth3phi);
    
    CCTK_REAL_VEC cdphi2 = fac1*(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + J32L*PDstandardNth3phi);
    
    CCTK_REAL_VEC cdphi3 = fac1*(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi);
    
    CCTK_REAL_VEC fac2 = IfThen(conformalMethod,khalf*pow(phiL,-2),0);
    
    CCTK_REAL_VEC cdphi211 = fac1*((dJ111L - Gt111*J11L - Gt211*J12L - 
      Gt311*J13L)*PDstandardNth1phi + 2*(J11L*(J21L*PDstandardNth12phi + 
      J31L*PDstandardNth13phi) + J21L*J31L*PDstandardNth23phi) + (dJ211L - 
      Gt111*J21L - Gt211*J22L - Gt311*J23L)*PDstandardNth2phi + (dJ311L - 
      Gt111*J31L - Gt211*J32L - Gt311*J33L)*PDstandardNth3phi + 
      PDstandardNth11phi*SQR(J11L) + PDstandardNth22phi*SQR(J21L) + 
      PDstandardNth33phi*SQR(J31L)) + fac2*SQR(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + J31L*PDstandardNth3phi);
    
    CCTK_REAL_VEC cdphi212 = fac2*(J11L*PDstandardNth1phi + 
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
    
    CCTK_REAL_VEC cdphi213 = fac2*(J11L*PDstandardNth1phi + 
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
    
    CCTK_REAL_VEC cdphi222 = fac1*((dJ122L - Gt122*J11L - Gt222*J12L - 
      Gt322*J13L)*PDstandardNth1phi + 2*(J12L*(J22L*PDstandardNth12phi + 
      J32L*PDstandardNth13phi) + J22L*J32L*PDstandardNth23phi) + (dJ222L - 
      Gt122*J21L - Gt222*J22L - Gt322*J23L)*PDstandardNth2phi + (dJ322L - 
      Gt122*J31L - Gt222*J32L - Gt322*J33L)*PDstandardNth3phi + 
      PDstandardNth11phi*SQR(J12L) + PDstandardNth22phi*SQR(J22L) + 
      PDstandardNth33phi*SQR(J32L)) + fac2*SQR(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + J32L*PDstandardNth3phi);
    
    CCTK_REAL_VEC cdphi223 = fac2*(J12L*PDstandardNth1phi + 
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
    
    CCTK_REAL_VEC cdphi233 = fac1*((dJ133L - Gt133*J11L - Gt233*J12L - 
      Gt333*J13L)*PDstandardNth1phi + 2*(J13L*(J23L*PDstandardNth12phi + 
      J33L*PDstandardNth13phi) + J23L*J33L*PDstandardNth23phi) + (dJ233L - 
      Gt133*J21L - Gt233*J22L - Gt333*J23L)*PDstandardNth2phi + (dJ333L - 
      Gt133*J31L - Gt233*J32L - Gt333*J33L)*PDstandardNth3phi + 
      PDstandardNth11phi*SQR(J13L) + PDstandardNth22phi*SQR(J23L) + 
      PDstandardNth33phi*SQR(J33L)) + fac2*SQR(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi);
    
    CCTK_REAL_VEC Rphi11 = -2*(cdphi211 + 2*(-1 + gt11L*gtu11)*SQR(cdphi1) 
      + gt11L*(cdphi211*gtu11 + 4*(cdphi1*(cdphi2*gtu21 + cdphi3*gtu31) + 
      cdphi2*cdphi3*gtu32) + cdphi233*gtu33 + gtu22*(cdphi222 + 
      2*SQR(cdphi2)) + 2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + 
      gtu33*SQR(cdphi3))));
    
    CCTK_REAL_VEC Rphi12 = -2*(cdphi212 + cdphi1*(cdphi2*(-2 + 
      4*gt12L*gtu21) + 4*cdphi3*gt12L*gtu31) + gt12L*(cdphi211*gtu11 + 
      4*cdphi2*cdphi3*gtu32 + 2*(cdphi212*gtu21 + cdphi213*gtu31 + 
      cdphi223*gtu32 + gtu11*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2)) 
      + gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL_VEC Rphi13 = -2*(cdphi213 + cdphi1*(4*cdphi2*gt13L*gtu21 + 
      cdphi3*(-2 + 4*gt13L*gtu31)) + gt13L*(cdphi211*gtu11 + 
      4*cdphi2*cdphi3*gtu32 + 2*(cdphi212*gtu21 + cdphi213*gtu31 + 
      cdphi223*gtu32 + gtu11*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2)) 
      + gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL_VEC Rphi22 = -2*(cdphi222 + 2*(-1 + gt22L*gtu22)*SQR(cdphi2) 
      + gt22L*(cdphi222*gtu22 + 4*(cdphi1*cdphi3*gtu31 + cdphi2*(cdphi1*gtu21 
      + cdphi3*gtu32)) + cdphi233*gtu33 + gtu11*(cdphi211 + 2*SQR(cdphi1)) + 
      2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + 
      gtu33*SQR(cdphi3))));
    
    CCTK_REAL_VEC Rphi23 = -2*(cdphi223 + cdphi2*(4*cdphi1*gt23L*gtu21 + 
      cdphi3*(-2 + 4*gt23L*gtu32)) + gt23L*(cdphi222*gtu22 + 
      4*cdphi1*cdphi3*gtu31 + gtu11*(cdphi211 + 2*SQR(cdphi1)) + 
      2*(cdphi212*gtu21 + cdphi213*gtu31 + cdphi223*gtu32 + 
      gtu22*SQR(cdphi2)) + gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL_VEC Rphi33 = -2*(cdphi233 + gt33L*((4*cdphi1*cdphi2 + 
      2*cdphi212)*gtu21 + 4*cdphi3*(cdphi1*gtu31 + cdphi2*gtu32) + 
      2*(cdphi213*gtu31 + cdphi223*gtu32) + cdphi233*gtu33 + gtu11*(cdphi211 
      + 2*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2))) + 2*(-1 + 
      gt33L*gtu33)*SQR(cdphi3));
    
    CCTK_REAL_VEC Atm11 = At11L*gtu11 + At12L*gtu21 + At13L*gtu31;
    
    CCTK_REAL_VEC Atm21 = At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    CCTK_REAL_VEC Atm31 = At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    CCTK_REAL_VEC Atm12 = At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    CCTK_REAL_VEC Atm22 = At12L*gtu21 + At22L*gtu22 + At23L*gtu32;
    
    CCTK_REAL_VEC Atm32 = At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    CCTK_REAL_VEC Atm13 = At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    CCTK_REAL_VEC Atm23 = At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    CCTK_REAL_VEC Atm33 = At13L*gtu31 + At23L*gtu32 + At33L*gtu33;
    
    CCTK_REAL_VEC e4phi = 
      IfThen(conformalMethod,pow(phiL,-2),exp(4*phiL));
    
    CCTK_REAL_VEC em4phi = INV(e4phi);
    
    CCTK_REAL_VEC g11 = e4phi*gt11L;
    
    CCTK_REAL_VEC g12 = e4phi*gt12L;
    
    CCTK_REAL_VEC g13 = e4phi*gt13L;
    
    CCTK_REAL_VEC g22 = e4phi*gt22L;
    
    CCTK_REAL_VEC g23 = e4phi*gt23L;
    
    CCTK_REAL_VEC g33 = e4phi*gt33L;
    
    CCTK_REAL_VEC gu11 = em4phi*gtu11;
    
    CCTK_REAL_VEC gu21 = em4phi*gtu21;
    
    CCTK_REAL_VEC gu31 = em4phi*gtu31;
    
    CCTK_REAL_VEC gu22 = em4phi*gtu22;
    
    CCTK_REAL_VEC gu32 = em4phi*gtu32;
    
    CCTK_REAL_VEC gu33 = em4phi*gtu33;
    
    CCTK_REAL_VEC R11 = Rphi11 + Rt11;
    
    CCTK_REAL_VEC R12 = Rphi12 + Rt12;
    
    CCTK_REAL_VEC R13 = Rphi13 + Rt13;
    
    CCTK_REAL_VEC R22 = Rphi22 + Rt22;
    
    CCTK_REAL_VEC R23 = Rphi23 + Rt23;
    
    CCTK_REAL_VEC R33 = Rphi33 + Rt33;
    
    CCTK_REAL_VEC trS = em4phi*(eTxxL*gtu11 + eTyyL*gtu22 + 2*(eTxyL*gtu21 
      + eTxzL*gtu31 + eTyzL*gtu32) + eTzzL*gtu33);
    
    CCTK_REAL_VEC Ats11 = (-dJ111L + 4*cdphi1*J11L + Gt111*J11L + 
      Gt211*J12L + Gt311*J13L)*PDstandardNth1alpha - 
      2*(J11L*J21L*PDstandardNth12alpha + J11L*J31L*PDstandardNth13alpha + 
      J21L*J31L*PDstandardNth23alpha) + (-dJ211L + 4*cdphi1*J21L + Gt111*J21L 
      + Gt211*J22L + Gt311*J23L)*PDstandardNth2alpha + (-dJ311L + 
      4*cdphi1*J31L + Gt111*J31L + Gt211*J32L + 
      Gt311*J33L)*PDstandardNth3alpha + alphaL*R11 - 
      PDstandardNth11alpha*SQR(J11L) - PDstandardNth22alpha*SQR(J21L) - 
      PDstandardNth33alpha*SQR(J31L);
    
    CCTK_REAL_VEC Ats12 = -(J11L*J12L*PDstandardNth11alpha) - 
      J12L*J21L*PDstandardNth12alpha - J11L*J22L*PDstandardNth12alpha - 
      J12L*J31L*PDstandardNth13alpha - J11L*J32L*PDstandardNth13alpha + 
      (-dJ112L + 2*cdphi2*J11L + Gt112*J11L + 2*cdphi1*J12L + Gt212*J12L + 
      Gt312*J13L)*PDstandardNth1alpha - J21L*J22L*PDstandardNth22alpha - 
      J22L*J31L*PDstandardNth23alpha - J21L*J32L*PDstandardNth23alpha + 
      (-dJ212L + 2*cdphi2*J21L + Gt112*J21L + 2*cdphi1*J22L + Gt212*J22L + 
      Gt312*J23L)*PDstandardNth2alpha - J31L*J32L*PDstandardNth33alpha - 
      dJ312L*PDstandardNth3alpha + 2*cdphi2*J31L*PDstandardNth3alpha + 
      Gt112*J31L*PDstandardNth3alpha + 2*cdphi1*J32L*PDstandardNth3alpha + 
      Gt212*J32L*PDstandardNth3alpha + Gt312*J33L*PDstandardNth3alpha + 
      alphaL*R12;
    
    CCTK_REAL_VEC Ats13 = -(J11L*J13L*PDstandardNth11alpha) - 
      J13L*J21L*PDstandardNth12alpha - J11L*J23L*PDstandardNth12alpha - 
      J13L*J31L*PDstandardNth13alpha - J11L*J33L*PDstandardNth13alpha + 
      (-dJ113L + 2*cdphi3*J11L + Gt113*J11L + Gt213*J12L + 2*cdphi1*J13L + 
      Gt313*J13L)*PDstandardNth1alpha - J21L*J23L*PDstandardNth22alpha - 
      J23L*J31L*PDstandardNth23alpha - J21L*J33L*PDstandardNth23alpha + 
      (-dJ213L + 2*cdphi3*J21L + Gt113*J21L + Gt213*J22L + 2*cdphi1*J23L + 
      Gt313*J23L)*PDstandardNth2alpha - J31L*J33L*PDstandardNth33alpha - 
      dJ313L*PDstandardNth3alpha + 2*cdphi3*J31L*PDstandardNth3alpha + 
      Gt113*J31L*PDstandardNth3alpha + Gt213*J32L*PDstandardNth3alpha + 
      2*cdphi1*J33L*PDstandardNth3alpha + Gt313*J33L*PDstandardNth3alpha + 
      alphaL*R13;
    
    CCTK_REAL_VEC Ats22 = (-dJ122L + Gt122*J11L + 4*cdphi2*J12L + 
      Gt222*J12L + Gt322*J13L)*PDstandardNth1alpha - 
      2*(J12L*J22L*PDstandardNth12alpha + J12L*J32L*PDstandardNth13alpha + 
      J22L*J32L*PDstandardNth23alpha) + (-dJ222L + Gt122*J21L + 4*cdphi2*J22L 
      + Gt222*J22L + Gt322*J23L)*PDstandardNth2alpha + (-dJ322L + Gt122*J31L 
      + 4*cdphi2*J32L + Gt222*J32L + Gt322*J33L)*PDstandardNth3alpha + 
      alphaL*R22 - PDstandardNth11alpha*SQR(J12L) - 
      PDstandardNth22alpha*SQR(J22L) - PDstandardNth33alpha*SQR(J32L);
    
    CCTK_REAL_VEC Ats23 = -(J12L*J13L*PDstandardNth11alpha) - 
      J13L*J22L*PDstandardNth12alpha - J12L*J23L*PDstandardNth12alpha - 
      J13L*J32L*PDstandardNth13alpha - J12L*J33L*PDstandardNth13alpha + 
      (-dJ123L + Gt123*J11L + 2*cdphi3*J12L + Gt223*J12L + 2*cdphi2*J13L + 
      Gt323*J13L)*PDstandardNth1alpha - J22L*J23L*PDstandardNth22alpha - 
      J23L*J32L*PDstandardNth23alpha - J22L*J33L*PDstandardNth23alpha + 
      (-dJ223L + Gt123*J21L + 2*cdphi3*J22L + Gt223*J22L + 2*cdphi2*J23L + 
      Gt323*J23L)*PDstandardNth2alpha - J32L*J33L*PDstandardNth33alpha - 
      dJ323L*PDstandardNth3alpha + Gt123*J31L*PDstandardNth3alpha + 
      2*cdphi3*J32L*PDstandardNth3alpha + Gt223*J32L*PDstandardNth3alpha + 
      2*cdphi2*J33L*PDstandardNth3alpha + Gt323*J33L*PDstandardNth3alpha + 
      alphaL*R23;
    
    CCTK_REAL_VEC Ats33 = (-dJ133L + Gt133*J11L + Gt233*J12L + 
      4*cdphi3*J13L + Gt333*J13L)*PDstandardNth1alpha - 
      2*(J13L*J23L*PDstandardNth12alpha + J13L*J33L*PDstandardNth13alpha + 
      J23L*J33L*PDstandardNth23alpha) + (-dJ233L + Gt133*J21L + Gt233*J22L + 
      4*cdphi3*J23L + Gt333*J23L)*PDstandardNth2alpha + (-dJ333L + Gt133*J31L 
      + Gt233*J32L + 4*cdphi3*J33L + Gt333*J33L)*PDstandardNth3alpha + 
      alphaL*R33 - PDstandardNth11alpha*SQR(J13L) - 
      PDstandardNth22alpha*SQR(J23L) - PDstandardNth33alpha*SQR(J33L);
    
    CCTK_REAL_VEC trAts = Ats11*gu11 + Ats22*gu22 + 2*(Ats12*gu21 + 
      Ats13*gu31 + Ats23*gu32) + Ats33*gu33;
    
    CCTK_REAL_VEC At11rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At11 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At11 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At11 - 
      At11L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + 2*(J11L*(At11L*PDstandardNth1beta1 + 
      At12L*PDstandardNth1beta2 + At13L*PDstandardNth1beta3) + 
      J21L*(At11L*PDstandardNth2beta1 + At12L*PDstandardNth2beta2 + 
      At13L*PDstandardNth2beta3) + J31L*(At11L*PDstandardNth3beta1 + 
      At12L*PDstandardNth3beta2 + At13L*PDstandardNth3beta3)) + (beta1L*J11L 
      + beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1At11 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2At11 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3At11 + 
      alphaL*(-2*(At11L*Atm11 + At12L*Atm21 + At13L*Atm31) + At11L*trKL) + 
      em4phi*(Ats11 - g11*kthird*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTxxL - 
      g11*kthird*trS)) + (J11L*PDupwindNthSymm1At11 + 
      J21L*PDupwindNthSymm2At11 + J31L*PDupwindNthSymm3At11)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At11 + J22L*PDupwindNthSymm2At11 + 
      J32L*PDupwindNthSymm3At11)*Abs(beta2L) + (J13L*PDupwindNthSymm1At11 + 
      J23L*PDupwindNthSymm2At11 + J33L*PDupwindNthSymm3At11)*Abs(beta3L);
    
    CCTK_REAL_VEC At12rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At12 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At12 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At12 + (At12L*J11L + 
      At11L*J12L)*PDstandardNth1beta1 + (At22L*J11L + 
      At12L*J12L)*PDstandardNth1beta2 + (At23L*J11L + 
      At13L*J12L)*PDstandardNth1beta3 + (At12L*J21L + 
      At11L*J22L)*PDstandardNth2beta1 + (At22L*J21L + 
      At12L*J22L)*PDstandardNth2beta2 + (At23L*J21L + 
      At13L*J22L)*PDstandardNth2beta3 + (At12L*J31L + 
      At11L*J32L)*PDstandardNth3beta1 + (At22L*J31L + 
      At12L*J32L)*PDstandardNth3beta2 + (At23L*J31L + 
      At13L*J32L)*PDstandardNth3beta3 - 
      At12L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1At12 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2At12 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3At12 + alphaL*(-2*(At11L*Atm12 + 
      At12L*Atm22 + At13L*Atm32) + At12L*trKL) + em4phi*(Ats12 - 
      g12*kthird*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTxyL - 
      g12*kthird*trS)) + (J11L*PDupwindNthSymm1At12 + 
      J21L*PDupwindNthSymm2At12 + J31L*PDupwindNthSymm3At12)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At12 + J22L*PDupwindNthSymm2At12 + 
      J32L*PDupwindNthSymm3At12)*Abs(beta2L) + (J13L*PDupwindNthSymm1At12 + 
      J23L*PDupwindNthSymm2At12 + J33L*PDupwindNthSymm3At12)*Abs(beta3L);
    
    CCTK_REAL_VEC At13rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At13 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At13 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At13 + (At13L*J11L + 
      At11L*J13L)*PDstandardNth1beta1 + (At23L*J11L + 
      At12L*J13L)*PDstandardNth1beta2 + (At33L*J11L + 
      At13L*J13L)*PDstandardNth1beta3 + (At13L*J21L + 
      At11L*J23L)*PDstandardNth2beta1 + (At23L*J21L + 
      At12L*J23L)*PDstandardNth2beta2 + (At33L*J21L + 
      At13L*J23L)*PDstandardNth2beta3 + (At13L*J31L + 
      At11L*J33L)*PDstandardNth3beta1 + (At23L*J31L + 
      At12L*J33L)*PDstandardNth3beta2 + (At33L*J31L + 
      At13L*J33L)*PDstandardNth3beta3 - 
      At13L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1At13 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2At13 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3At13 + alphaL*(-2*(At11L*Atm13 + 
      At12L*Atm23 + At13L*Atm33) + At13L*trKL) + em4phi*(Ats13 - 
      g13*kthird*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTxzL - 
      g13*kthird*trS)) + (J11L*PDupwindNthSymm1At13 + 
      J21L*PDupwindNthSymm2At13 + J31L*PDupwindNthSymm3At13)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At13 + J22L*PDupwindNthSymm2At13 + 
      J32L*PDupwindNthSymm3At13)*Abs(beta2L) + (J13L*PDupwindNthSymm1At13 + 
      J23L*PDupwindNthSymm2At13 + J33L*PDupwindNthSymm3At13)*Abs(beta3L);
    
    CCTK_REAL_VEC At22rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At22 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At22 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At22 - 
      At22L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + 2*(J12L*(At12L*PDstandardNth1beta1 + 
      At22L*PDstandardNth1beta2 + At23L*PDstandardNth1beta3) + 
      J22L*(At12L*PDstandardNth2beta1 + At22L*PDstandardNth2beta2 + 
      At23L*PDstandardNth2beta3) + J32L*(At12L*PDstandardNth3beta1 + 
      At22L*PDstandardNth3beta2 + At23L*PDstandardNth3beta3)) + (beta1L*J11L 
      + beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1At22 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2At22 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3At22 + 
      alphaL*(-2*(At12L*Atm12 + At22L*Atm22 + At23L*Atm32) + At22L*trKL) + 
      em4phi*(Ats22 - g22*kthird*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTyyL - 
      g22*kthird*trS)) + (J11L*PDupwindNthSymm1At22 + 
      J21L*PDupwindNthSymm2At22 + J31L*PDupwindNthSymm3At22)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At22 + J22L*PDupwindNthSymm2At22 + 
      J32L*PDupwindNthSymm3At22)*Abs(beta2L) + (J13L*PDupwindNthSymm1At22 + 
      J23L*PDupwindNthSymm2At22 + J33L*PDupwindNthSymm3At22)*Abs(beta3L);
    
    CCTK_REAL_VEC At23rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At23 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At23 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At23 + (At13L*J12L + 
      At12L*J13L)*PDstandardNth1beta1 + (At23L*J12L + 
      At22L*J13L)*PDstandardNth1beta2 + (At33L*J12L + 
      At23L*J13L)*PDstandardNth1beta3 + (At13L*J22L + 
      At12L*J23L)*PDstandardNth2beta1 + (At23L*J22L + 
      At22L*J23L)*PDstandardNth2beta2 + (At33L*J22L + 
      At23L*J23L)*PDstandardNth2beta3 + (At13L*J32L + 
      At12L*J33L)*PDstandardNth3beta1 + (At23L*J32L + 
      At22L*J33L)*PDstandardNth3beta2 + (At33L*J32L + 
      At23L*J33L)*PDstandardNth3beta3 - 
      At23L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1At23 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2At23 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3At23 + alphaL*(-2*(At12L*Atm13 + 
      At22L*Atm23 + At23L*Atm33) + At23L*trKL) + em4phi*(Ats23 - 
      g23*kthird*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTyzL - 
      g23*kthird*trS)) + (J11L*PDupwindNthSymm1At23 + 
      J21L*PDupwindNthSymm2At23 + J31L*PDupwindNthSymm3At23)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At23 + J22L*PDupwindNthSymm2At23 + 
      J32L*PDupwindNthSymm3At23)*Abs(beta2L) + (J13L*PDupwindNthSymm1At23 + 
      J23L*PDupwindNthSymm2At23 + J33L*PDupwindNthSymm3At23)*Abs(beta3L);
    
    CCTK_REAL_VEC At33rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1At33 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2At33 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3At33 - 
      At33L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + 2*(J13L*(At13L*PDstandardNth1beta1 + 
      At23L*PDstandardNth1beta2 + At33L*PDstandardNth1beta3) + 
      J23L*(At13L*PDstandardNth2beta1 + At23L*PDstandardNth2beta2 + 
      At33L*PDstandardNth2beta3) + J33L*(At13L*PDstandardNth3beta1 + 
      At23L*PDstandardNth3beta2 + At33L*PDstandardNth3beta3)) + (beta1L*J11L 
      + beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1At33 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2At33 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3At33 + 
      alphaL*(-2*(At13L*Atm13 + At23L*Atm23 + At33L*Atm33) + At33L*trKL) + 
      em4phi*(Ats33 - g33*kthird*trAts - 
      25.13274122871834590770114706623602307358*alphaL*(eTzzL - 
      g33*kthird*trS)) + (J11L*PDupwindNthSymm1At33 + 
      J21L*PDupwindNthSymm2At33 + J31L*PDupwindNthSymm3At33)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1At33 + J22L*PDupwindNthSymm2At33 + 
      J32L*PDupwindNthSymm3At33)*Abs(beta2L) + (J13L*PDupwindNthSymm1At33 + 
      J23L*PDupwindNthSymm2At33 + J33L*PDupwindNthSymm3At33)*Abs(beta3L);
    
    
    /* Copy local copies back to grid functions */
    vec_store_nta(At11rhs[index],At11rhsL);
    vec_store_nta(At12rhs[index],At12rhsL);
    vec_store_nta(At13rhs[index],At13rhsL);
    vec_store_nta(At22rhs[index],At22rhsL);
    vec_store_nta(At23rhs[index],At23rhsL);
    vec_store_nta(At33rhs[index],At33rhsL);
  i += CCTK_REAL_VEC_SIZE-1;
  }
  LC_ENDLOOP3 (ML_BSSN_MP_RHS2);
}

extern "C" void ML_BSSN_MP_RHS2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_RHS2_Body);
}
