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

void ML_ADMQuantities_MP_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  // CCTK_REAL p1odx = INITVALUE;
  // CCTK_REAL p1ody = INITVALUE;
  // CCTK_REAL p1odz = INITVALUE;
  // CCTK_REAL pm1o12dx2 = INITVALUE;
  // CCTK_REAL pm1o12dy2 = INITVALUE;
  // CCTK_REAL pm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMQuantities_MP_Body");
  }
  
  if (cctk_iteration % ML_ADMQuantities_MP_calc_every != ML_ADMQuantities_MP_calc_offset)
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
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o12dx2 = -pow(dx,-2)/12.;
  CCTK_REAL const pm1o12dy2 = -pow(dy,-2)/12.;
  CCTK_REAL const pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_ADMQuantities_MP,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    // int subblock_index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int const subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    // CCTK_REAL Atm11 = INITVALUE, Atm12 = INITVALUE, Atm13 = INITVALUE, Atm21 = INITVALUE, Atm22 = INITVALUE, Atm23 = INITVALUE;
    // CCTK_REAL Atm31 = INITVALUE, Atm32 = INITVALUE, Atm33 = INITVALUE;
    // CCTK_REAL detgt = INITVALUE;
    // CCTK_REAL dgtu111 = INITVALUE, dgtu112 = INITVALUE, dgtu113 = INITVALUE, dgtu211 = INITVALUE, dgtu212 = INITVALUE, dgtu213 = INITVALUE;
    // CCTK_REAL dgtu221 = INITVALUE, dgtu222 = INITVALUE, dgtu223 = INITVALUE, dgtu311 = INITVALUE, dgtu312 = INITVALUE, dgtu313 = INITVALUE;
    // CCTK_REAL dgtu321 = INITVALUE, dgtu322 = INITVALUE, dgtu323 = INITVALUE, dgtu331 = INITVALUE, dgtu332 = INITVALUE, dgtu333 = INITVALUE;
    // CCTK_REAL ephi = INITVALUE;
    // CCTK_REAL Gt111 = INITVALUE, Gt112 = INITVALUE, Gt113 = INITVALUE, Gt122 = INITVALUE, Gt123 = INITVALUE, Gt133 = INITVALUE;
    // CCTK_REAL Gt211 = INITVALUE, Gt212 = INITVALUE, Gt213 = INITVALUE, Gt222 = INITVALUE, Gt223 = INITVALUE, Gt233 = INITVALUE;
    // CCTK_REAL Gt311 = INITVALUE, Gt312 = INITVALUE, Gt313 = INITVALUE, Gt322 = INITVALUE, Gt323 = INITVALUE, Gt333 = INITVALUE;
    // CCTK_REAL Gtl111 = INITVALUE, Gtl112 = INITVALUE, Gtl113 = INITVALUE, Gtl122 = INITVALUE, Gtl123 = INITVALUE, Gtl133 = INITVALUE;
    // CCTK_REAL Gtl211 = INITVALUE, Gtl212 = INITVALUE, Gtl213 = INITVALUE, Gtl222 = INITVALUE, Gtl223 = INITVALUE, Gtl233 = INITVALUE;
    // CCTK_REAL Gtl311 = INITVALUE, Gtl312 = INITVALUE, Gtl313 = INITVALUE, Gtl322 = INITVALUE, Gtl323 = INITVALUE, Gtl333 = INITVALUE;
    // CCTK_REAL Gtlu111 = INITVALUE, Gtlu112 = INITVALUE, Gtlu113 = INITVALUE, Gtlu121 = INITVALUE, Gtlu122 = INITVALUE, Gtlu123 = INITVALUE;
    // CCTK_REAL Gtlu131 = INITVALUE, Gtlu132 = INITVALUE, Gtlu133 = INITVALUE, Gtlu211 = INITVALUE, Gtlu212 = INITVALUE, Gtlu213 = INITVALUE;
    // CCTK_REAL Gtlu221 = INITVALUE, Gtlu222 = INITVALUE, Gtlu223 = INITVALUE, Gtlu231 = INITVALUE, Gtlu232 = INITVALUE, Gtlu233 = INITVALUE;
    // CCTK_REAL Gtlu311 = INITVALUE, Gtlu312 = INITVALUE, Gtlu313 = INITVALUE, Gtlu321 = INITVALUE, Gtlu322 = INITVALUE, Gtlu323 = INITVALUE;
    // CCTK_REAL Gtlu331 = INITVALUE, Gtlu332 = INITVALUE, Gtlu333 = INITVALUE;
    // CCTK_REAL gtu11 = INITVALUE, gtu21 = INITVALUE, gtu22 = INITVALUE, gtu31 = INITVALUE, gtu32 = INITVALUE, gtu33 = INITVALUE;
    // CCTK_REAL rho = INITVALUE;
    // CCTK_REAL Rt11 = INITVALUE, Rt12 = INITVALUE, Rt13 = INITVALUE, Rt22 = INITVALUE, Rt23 = INITVALUE, Rt33 = INITVALUE;
    // CCTK_REAL S1 = INITVALUE, S2 = INITVALUE, S3 = INITVALUE;
    // CCTK_REAL trRt = INITVALUE;
    // CCTK_REAL Xtn1 = INITVALUE, Xtn2 = INITVALUE, Xtn3 = INITVALUE;
    
    /* Declare local copies of grid functions */
    // CCTK_REAL alphaL = INITVALUE;
    // CCTK_REAL At11L = INITVALUE, At12L = INITVALUE, At13L = INITVALUE, At22L = INITVALUE, At23L = INITVALUE, At33L = INITVALUE;
    // CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    // CCTK_REAL dJ111L = INITVALUE, dJ112L = INITVALUE, dJ113L = INITVALUE, dJ122L = INITVALUE, dJ123L = INITVALUE, dJ133L = INITVALUE;
    // CCTK_REAL dJ211L = INITVALUE, dJ212L = INITVALUE, dJ213L = INITVALUE, dJ222L = INITVALUE, dJ223L = INITVALUE, dJ233L = INITVALUE;
    // CCTK_REAL dJ311L = INITVALUE, dJ312L = INITVALUE, dJ313L = INITVALUE, dJ322L = INITVALUE, dJ323L = INITVALUE, dJ333L = INITVALUE;
    // CCTK_REAL eTttL = INITVALUE;
    // CCTK_REAL eTtxL = INITVALUE;
    // CCTK_REAL eTtyL = INITVALUE;
    // CCTK_REAL eTtzL = INITVALUE;
    // CCTK_REAL eTxxL = INITVALUE;
    // CCTK_REAL eTxyL = INITVALUE;
    // CCTK_REAL eTxzL = INITVALUE;
    // CCTK_REAL eTyyL = INITVALUE;
    // CCTK_REAL eTyzL = INITVALUE;
    // CCTK_REAL eTzzL = INITVALUE;
    // CCTK_REAL gt11L = INITVALUE, gt12L = INITVALUE, gt13L = INITVALUE, gt22L = INITVALUE, gt23L = INITVALUE, gt33L = INITVALUE;
    // CCTK_REAL J11L = INITVALUE, J12L = INITVALUE, J13L = INITVALUE, J21L = INITVALUE, J22L = INITVALUE, J23L = INITVALUE;
    // CCTK_REAL J31L = INITVALUE, J32L = INITVALUE, J33L = INITVALUE;
    // CCTK_REAL Jadm1L = INITVALUE, Jadm2L = INITVALUE, Jadm3L = INITVALUE;
    // CCTK_REAL MadmL = INITVALUE;
    // CCTK_REAL phiL = INITVALUE;
    // CCTK_REAL trKL = INITVALUE;
    // CCTK_REAL xL = INITVALUE;
    // CCTK_REAL Xt1L = INITVALUE, Xt2L = INITVALUE, Xt3L = INITVALUE;
    // CCTK_REAL yL = INITVALUE;
    // CCTK_REAL zL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
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
    CCTK_REAL const alphaL = alpha[index];
    CCTK_REAL const At11L = At11[index];
    CCTK_REAL const At12L = At12[index];
    CCTK_REAL const At13L = At13[index];
    CCTK_REAL const At22L = At22[index];
    CCTK_REAL const At23L = At23[index];
    CCTK_REAL const At33L = At33[index];
    CCTK_REAL const beta1L = beta1[index];
    CCTK_REAL const beta2L = beta2[index];
    CCTK_REAL const beta3L = beta3[index];
    CCTK_REAL const dJ111L = dJ111[index];
    CCTK_REAL const dJ112L = dJ112[index];
    CCTK_REAL const dJ113L = dJ113[index];
    CCTK_REAL const dJ122L = dJ122[index];
    CCTK_REAL const dJ123L = dJ123[index];
    CCTK_REAL const dJ133L = dJ133[index];
    CCTK_REAL const dJ211L = dJ211[index];
    CCTK_REAL const dJ212L = dJ212[index];
    CCTK_REAL const dJ213L = dJ213[index];
    CCTK_REAL const dJ222L = dJ222[index];
    CCTK_REAL const dJ223L = dJ223[index];
    CCTK_REAL const dJ233L = dJ233[index];
    CCTK_REAL const dJ311L = dJ311[index];
    CCTK_REAL const dJ312L = dJ312[index];
    CCTK_REAL const dJ313L = dJ313[index];
    CCTK_REAL const dJ322L = dJ322[index];
    CCTK_REAL const dJ323L = dJ323[index];
    CCTK_REAL const dJ333L = dJ333[index];
    CCTK_REAL const eTttL = (stress_energy_state) ? (eTtt[index]) : 0.0;
    CCTK_REAL const eTtxL = (stress_energy_state) ? (eTtx[index]) : 0.0;
    CCTK_REAL const eTtyL = (stress_energy_state) ? (eTty[index]) : 0.0;
    CCTK_REAL const eTtzL = (stress_energy_state) ? (eTtz[index]) : 0.0;
    CCTK_REAL const eTxxL = (stress_energy_state) ? (eTxx[index]) : 0.0;
    CCTK_REAL const eTxyL = (stress_energy_state) ? (eTxy[index]) : 0.0;
    CCTK_REAL const eTxzL = (stress_energy_state) ? (eTxz[index]) : 0.0;
    CCTK_REAL const eTyyL = (stress_energy_state) ? (eTyy[index]) : 0.0;
    CCTK_REAL const eTyzL = (stress_energy_state) ? (eTyz[index]) : 0.0;
    CCTK_REAL const eTzzL = (stress_energy_state) ? (eTzz[index]) : 0.0;
    CCTK_REAL const gt11L = gt11[index];
    CCTK_REAL const gt12L = gt12[index];
    CCTK_REAL const gt13L = gt13[index];
    CCTK_REAL const gt22L = gt22[index];
    CCTK_REAL const gt23L = gt23[index];
    CCTK_REAL const gt33L = gt33[index];
    CCTK_REAL const J11L = J11[index];
    CCTK_REAL const J12L = J12[index];
    CCTK_REAL const J13L = J13[index];
    CCTK_REAL const J21L = J21[index];
    CCTK_REAL const J22L = J22[index];
    CCTK_REAL const J23L = J23[index];
    CCTK_REAL const J31L = J31[index];
    CCTK_REAL const J32L = J32[index];
    CCTK_REAL const J33L = J33[index];
    CCTK_REAL const phiL = phi[index];
    CCTK_REAL const trKL = trK[index];
    CCTK_REAL const xL = x[index];
    CCTK_REAL const Xt1L = Xt1[index];
    CCTK_REAL const Xt2L = Xt2[index];
    CCTK_REAL const Xt3L = Xt3[index];
    CCTK_REAL const yL = y[index];
    CCTK_REAL const zL = z[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
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
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL const detgt  =  1;
    
    CCTK_REAL const gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL const gtu21  =  (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL const gtu31  =  (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL const gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL const gtu32  =  (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL const gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL const dgtu111  =  -2*(gtu11*gtu21*J11L*PDstandardNth1gt12 + gtu11*gtu31*J11L*PDstandardNth1gt13 + 
           gtu21*gtu31*J11L*PDstandardNth1gt23 + gtu11*gtu21*J21L*PDstandardNth2gt12 + 
           gtu11*gtu31*J21L*PDstandardNth2gt13 + gtu21*gtu31*J21L*PDstandardNth2gt23 + 
           gtu11*gtu21*J31L*PDstandardNth3gt12 + gtu11*gtu31*J31L*PDstandardNth3gt13 + gtu21*gtu31*J31L*PDstandardNth3gt23)\
         - J11L*PDstandardNth1gt11*SQR(gtu11) - J21L*PDstandardNth2gt11*SQR(gtu11) - J31L*PDstandardNth3gt11*SQR(gtu11) - 
        J11L*PDstandardNth1gt22*SQR(gtu21) - J21L*PDstandardNth2gt22*SQR(gtu21) - J31L*PDstandardNth3gt22*SQR(gtu21) - 
        J11L*PDstandardNth1gt33*SQR(gtu31) - J21L*PDstandardNth2gt33*SQR(gtu31) - J31L*PDstandardNth3gt33*SQR(gtu31);
    
    CCTK_REAL const dgtu211  =  -(gtu11*gtu21*J11L*PDstandardNth1gt11) - gtu11*gtu22*J11L*PDstandardNth1gt12 - 
        gtu21*gtu31*J11L*PDstandardNth1gt13 - gtu11*gtu32*J11L*PDstandardNth1gt13 - gtu21*gtu22*J11L*PDstandardNth1gt22 - 
        gtu22*gtu31*J11L*PDstandardNth1gt23 - gtu21*gtu32*J11L*PDstandardNth1gt23 - gtu31*gtu32*J11L*PDstandardNth1gt33 - 
        gtu11*gtu21*J21L*PDstandardNth2gt11 - gtu11*gtu22*J21L*PDstandardNth2gt12 - gtu21*gtu31*J21L*PDstandardNth2gt13 - 
        gtu11*gtu32*J21L*PDstandardNth2gt13 - gtu21*gtu22*J21L*PDstandardNth2gt22 - gtu22*gtu31*J21L*PDstandardNth2gt23 - 
        gtu21*gtu32*J21L*PDstandardNth2gt23 - gtu31*gtu32*J21L*PDstandardNth2gt33 - gtu11*gtu21*J31L*PDstandardNth3gt11 - 
        gtu11*gtu22*J31L*PDstandardNth3gt12 - gtu21*gtu31*J31L*PDstandardNth3gt13 - gtu11*gtu32*J31L*PDstandardNth3gt13 - 
        gtu21*gtu22*J31L*PDstandardNth3gt22 - gtu22*gtu31*J31L*PDstandardNth3gt23 - gtu21*gtu32*J31L*PDstandardNth3gt23 - 
        gtu31*gtu32*J31L*PDstandardNth3gt33 - J11L*PDstandardNth1gt12*SQR(gtu21) - J21L*PDstandardNth2gt12*SQR(gtu21) - 
        J31L*PDstandardNth3gt12*SQR(gtu21);
    
    CCTK_REAL const dgtu311  =  -(gtu11*gtu31*J11L*PDstandardNth1gt11) - gtu21*gtu31*J11L*PDstandardNth1gt12 - 
        gtu11*gtu32*J11L*PDstandardNth1gt12 - gtu11*gtu33*J11L*PDstandardNth1gt13 - gtu21*gtu32*J11L*PDstandardNth1gt22 - 
        gtu31*gtu32*J11L*PDstandardNth1gt23 - gtu21*gtu33*J11L*PDstandardNth1gt23 - gtu31*gtu33*J11L*PDstandardNth1gt33 - 
        gtu11*gtu31*J21L*PDstandardNth2gt11 - gtu21*gtu31*J21L*PDstandardNth2gt12 - gtu11*gtu32*J21L*PDstandardNth2gt12 - 
        gtu11*gtu33*J21L*PDstandardNth2gt13 - gtu21*gtu32*J21L*PDstandardNth2gt22 - gtu31*gtu32*J21L*PDstandardNth2gt23 - 
        gtu21*gtu33*J21L*PDstandardNth2gt23 - gtu31*gtu33*J21L*PDstandardNth2gt33 - gtu11*gtu31*J31L*PDstandardNth3gt11 - 
        gtu21*gtu31*J31L*PDstandardNth3gt12 - gtu11*gtu32*J31L*PDstandardNth3gt12 - gtu11*gtu33*J31L*PDstandardNth3gt13 - 
        gtu21*gtu32*J31L*PDstandardNth3gt22 - gtu31*gtu32*J31L*PDstandardNth3gt23 - gtu21*gtu33*J31L*PDstandardNth3gt23 - 
        gtu31*gtu33*J31L*PDstandardNth3gt33 - J11L*PDstandardNth1gt13*SQR(gtu31) - J21L*PDstandardNth2gt13*SQR(gtu31) - 
        J31L*PDstandardNth3gt13*SQR(gtu31);
    
    CCTK_REAL const dgtu221  =  -2*(gtu21*gtu22*J11L*PDstandardNth1gt12 + gtu21*gtu32*J11L*PDstandardNth1gt13 + 
           gtu22*gtu32*J11L*PDstandardNth1gt23 + gtu21*gtu22*J21L*PDstandardNth2gt12 + 
           gtu21*gtu32*J21L*PDstandardNth2gt13 + gtu22*gtu32*J21L*PDstandardNth2gt23 + 
           gtu21*gtu22*J31L*PDstandardNth3gt12 + gtu21*gtu32*J31L*PDstandardNth3gt13 + gtu22*gtu32*J31L*PDstandardNth3gt23)\
         - J11L*PDstandardNth1gt11*SQR(gtu21) - J21L*PDstandardNth2gt11*SQR(gtu21) - J31L*PDstandardNth3gt11*SQR(gtu21) - 
        J11L*PDstandardNth1gt22*SQR(gtu22) - J21L*PDstandardNth2gt22*SQR(gtu22) - J31L*PDstandardNth3gt22*SQR(gtu22) - 
        J11L*PDstandardNth1gt33*SQR(gtu32) - J21L*PDstandardNth2gt33*SQR(gtu32) - J31L*PDstandardNth3gt33*SQR(gtu32);
    
    CCTK_REAL const dgtu321  =  -(gtu21*gtu31*J11L*PDstandardNth1gt11) - gtu22*gtu31*J11L*PDstandardNth1gt12 - 
        gtu21*gtu32*J11L*PDstandardNth1gt12 - gtu31*gtu32*J11L*PDstandardNth1gt13 - gtu21*gtu33*J11L*PDstandardNth1gt13 - 
        gtu22*gtu32*J11L*PDstandardNth1gt22 - gtu22*gtu33*J11L*PDstandardNth1gt23 - gtu32*gtu33*J11L*PDstandardNth1gt33 - 
        gtu21*gtu31*J21L*PDstandardNth2gt11 - gtu22*gtu31*J21L*PDstandardNth2gt12 - gtu21*gtu32*J21L*PDstandardNth2gt12 - 
        gtu31*gtu32*J21L*PDstandardNth2gt13 - gtu21*gtu33*J21L*PDstandardNth2gt13 - gtu22*gtu32*J21L*PDstandardNth2gt22 - 
        gtu22*gtu33*J21L*PDstandardNth2gt23 - gtu32*gtu33*J21L*PDstandardNth2gt33 - gtu21*gtu31*J31L*PDstandardNth3gt11 - 
        gtu22*gtu31*J31L*PDstandardNth3gt12 - gtu21*gtu32*J31L*PDstandardNth3gt12 - gtu31*gtu32*J31L*PDstandardNth3gt13 - 
        gtu21*gtu33*J31L*PDstandardNth3gt13 - gtu22*gtu32*J31L*PDstandardNth3gt22 - gtu22*gtu33*J31L*PDstandardNth3gt23 - 
        gtu32*gtu33*J31L*PDstandardNth3gt33 - J11L*PDstandardNth1gt23*SQR(gtu32) - J21L*PDstandardNth2gt23*SQR(gtu32) - 
        J31L*PDstandardNth3gt23*SQR(gtu32);
    
    CCTK_REAL const dgtu331  =  -2*(gtu31*gtu32*J11L*PDstandardNth1gt12 + gtu31*gtu33*J11L*PDstandardNth1gt13 + 
           gtu32*gtu33*J11L*PDstandardNth1gt23 + gtu31*gtu32*J21L*PDstandardNth2gt12 + 
           gtu31*gtu33*J21L*PDstandardNth2gt13 + gtu32*gtu33*J21L*PDstandardNth2gt23 + 
           gtu31*gtu32*J31L*PDstandardNth3gt12 + gtu31*gtu33*J31L*PDstandardNth3gt13 + gtu32*gtu33*J31L*PDstandardNth3gt23)\
         - J11L*PDstandardNth1gt11*SQR(gtu31) - J21L*PDstandardNth2gt11*SQR(gtu31) - J31L*PDstandardNth3gt11*SQR(gtu31) - 
        J11L*PDstandardNth1gt22*SQR(gtu32) - J21L*PDstandardNth2gt22*SQR(gtu32) - J31L*PDstandardNth3gt22*SQR(gtu32) - 
        J11L*PDstandardNth1gt33*SQR(gtu33) - J21L*PDstandardNth2gt33*SQR(gtu33) - J31L*PDstandardNth3gt33*SQR(gtu33);
    
    CCTK_REAL const dgtu112  =  -2*(gtu11*gtu21*J12L*PDstandardNth1gt12 + gtu11*gtu31*J12L*PDstandardNth1gt13 + 
           gtu21*gtu31*J12L*PDstandardNth1gt23 + gtu11*gtu21*J22L*PDstandardNth2gt12 + 
           gtu11*gtu31*J22L*PDstandardNth2gt13 + gtu21*gtu31*J22L*PDstandardNth2gt23 + 
           gtu11*gtu21*J32L*PDstandardNth3gt12 + gtu11*gtu31*J32L*PDstandardNth3gt13 + gtu21*gtu31*J32L*PDstandardNth3gt23)\
         - J12L*PDstandardNth1gt11*SQR(gtu11) - J22L*PDstandardNth2gt11*SQR(gtu11) - J32L*PDstandardNth3gt11*SQR(gtu11) - 
        J12L*PDstandardNth1gt22*SQR(gtu21) - J22L*PDstandardNth2gt22*SQR(gtu21) - J32L*PDstandardNth3gt22*SQR(gtu21) - 
        J12L*PDstandardNth1gt33*SQR(gtu31) - J22L*PDstandardNth2gt33*SQR(gtu31) - J32L*PDstandardNth3gt33*SQR(gtu31);
    
    CCTK_REAL const dgtu212  =  -(gtu11*gtu21*J12L*PDstandardNth1gt11) - gtu11*gtu22*J12L*PDstandardNth1gt12 - 
        gtu21*gtu31*J12L*PDstandardNth1gt13 - gtu11*gtu32*J12L*PDstandardNth1gt13 - gtu21*gtu22*J12L*PDstandardNth1gt22 - 
        gtu22*gtu31*J12L*PDstandardNth1gt23 - gtu21*gtu32*J12L*PDstandardNth1gt23 - gtu31*gtu32*J12L*PDstandardNth1gt33 - 
        gtu11*gtu21*J22L*PDstandardNth2gt11 - gtu11*gtu22*J22L*PDstandardNth2gt12 - gtu21*gtu31*J22L*PDstandardNth2gt13 - 
        gtu11*gtu32*J22L*PDstandardNth2gt13 - gtu21*gtu22*J22L*PDstandardNth2gt22 - gtu22*gtu31*J22L*PDstandardNth2gt23 - 
        gtu21*gtu32*J22L*PDstandardNth2gt23 - gtu31*gtu32*J22L*PDstandardNth2gt33 - gtu11*gtu21*J32L*PDstandardNth3gt11 - 
        gtu11*gtu22*J32L*PDstandardNth3gt12 - gtu21*gtu31*J32L*PDstandardNth3gt13 - gtu11*gtu32*J32L*PDstandardNth3gt13 - 
        gtu21*gtu22*J32L*PDstandardNth3gt22 - gtu22*gtu31*J32L*PDstandardNth3gt23 - gtu21*gtu32*J32L*PDstandardNth3gt23 - 
        gtu31*gtu32*J32L*PDstandardNth3gt33 - J12L*PDstandardNth1gt12*SQR(gtu21) - J22L*PDstandardNth2gt12*SQR(gtu21) - 
        J32L*PDstandardNth3gt12*SQR(gtu21);
    
    CCTK_REAL const dgtu312  =  -(gtu11*gtu31*J12L*PDstandardNth1gt11) - gtu21*gtu31*J12L*PDstandardNth1gt12 - 
        gtu11*gtu32*J12L*PDstandardNth1gt12 - gtu11*gtu33*J12L*PDstandardNth1gt13 - gtu21*gtu32*J12L*PDstandardNth1gt22 - 
        gtu31*gtu32*J12L*PDstandardNth1gt23 - gtu21*gtu33*J12L*PDstandardNth1gt23 - gtu31*gtu33*J12L*PDstandardNth1gt33 - 
        gtu11*gtu31*J22L*PDstandardNth2gt11 - gtu21*gtu31*J22L*PDstandardNth2gt12 - gtu11*gtu32*J22L*PDstandardNth2gt12 - 
        gtu11*gtu33*J22L*PDstandardNth2gt13 - gtu21*gtu32*J22L*PDstandardNth2gt22 - gtu31*gtu32*J22L*PDstandardNth2gt23 - 
        gtu21*gtu33*J22L*PDstandardNth2gt23 - gtu31*gtu33*J22L*PDstandardNth2gt33 - gtu11*gtu31*J32L*PDstandardNth3gt11 - 
        gtu21*gtu31*J32L*PDstandardNth3gt12 - gtu11*gtu32*J32L*PDstandardNth3gt12 - gtu11*gtu33*J32L*PDstandardNth3gt13 - 
        gtu21*gtu32*J32L*PDstandardNth3gt22 - gtu31*gtu32*J32L*PDstandardNth3gt23 - gtu21*gtu33*J32L*PDstandardNth3gt23 - 
        gtu31*gtu33*J32L*PDstandardNth3gt33 - J12L*PDstandardNth1gt13*SQR(gtu31) - J22L*PDstandardNth2gt13*SQR(gtu31) - 
        J32L*PDstandardNth3gt13*SQR(gtu31);
    
    CCTK_REAL const dgtu222  =  -2*(gtu21*gtu22*J12L*PDstandardNth1gt12 + gtu21*gtu32*J12L*PDstandardNth1gt13 + 
           gtu22*gtu32*J12L*PDstandardNth1gt23 + gtu21*gtu22*J22L*PDstandardNth2gt12 + 
           gtu21*gtu32*J22L*PDstandardNth2gt13 + gtu22*gtu32*J22L*PDstandardNth2gt23 + 
           gtu21*gtu22*J32L*PDstandardNth3gt12 + gtu21*gtu32*J32L*PDstandardNth3gt13 + gtu22*gtu32*J32L*PDstandardNth3gt23)\
         - J12L*PDstandardNth1gt11*SQR(gtu21) - J22L*PDstandardNth2gt11*SQR(gtu21) - J32L*PDstandardNth3gt11*SQR(gtu21) - 
        J12L*PDstandardNth1gt22*SQR(gtu22) - J22L*PDstandardNth2gt22*SQR(gtu22) - J32L*PDstandardNth3gt22*SQR(gtu22) - 
        J12L*PDstandardNth1gt33*SQR(gtu32) - J22L*PDstandardNth2gt33*SQR(gtu32) - J32L*PDstandardNth3gt33*SQR(gtu32);
    
    CCTK_REAL const dgtu322  =  -(gtu21*gtu31*J12L*PDstandardNth1gt11) - gtu22*gtu31*J12L*PDstandardNth1gt12 - 
        gtu21*gtu32*J12L*PDstandardNth1gt12 - gtu31*gtu32*J12L*PDstandardNth1gt13 - gtu21*gtu33*J12L*PDstandardNth1gt13 - 
        gtu22*gtu32*J12L*PDstandardNth1gt22 - gtu22*gtu33*J12L*PDstandardNth1gt23 - gtu32*gtu33*J12L*PDstandardNth1gt33 - 
        gtu21*gtu31*J22L*PDstandardNth2gt11 - gtu22*gtu31*J22L*PDstandardNth2gt12 - gtu21*gtu32*J22L*PDstandardNth2gt12 - 
        gtu31*gtu32*J22L*PDstandardNth2gt13 - gtu21*gtu33*J22L*PDstandardNth2gt13 - gtu22*gtu32*J22L*PDstandardNth2gt22 - 
        gtu22*gtu33*J22L*PDstandardNth2gt23 - gtu32*gtu33*J22L*PDstandardNth2gt33 - gtu21*gtu31*J32L*PDstandardNth3gt11 - 
        gtu22*gtu31*J32L*PDstandardNth3gt12 - gtu21*gtu32*J32L*PDstandardNth3gt12 - gtu31*gtu32*J32L*PDstandardNth3gt13 - 
        gtu21*gtu33*J32L*PDstandardNth3gt13 - gtu22*gtu32*J32L*PDstandardNth3gt22 - gtu22*gtu33*J32L*PDstandardNth3gt23 - 
        gtu32*gtu33*J32L*PDstandardNth3gt33 - J12L*PDstandardNth1gt23*SQR(gtu32) - J22L*PDstandardNth2gt23*SQR(gtu32) - 
        J32L*PDstandardNth3gt23*SQR(gtu32);
    
    CCTK_REAL const dgtu332  =  -2*(gtu31*gtu32*J12L*PDstandardNth1gt12 + gtu31*gtu33*J12L*PDstandardNth1gt13 + 
           gtu32*gtu33*J12L*PDstandardNth1gt23 + gtu31*gtu32*J22L*PDstandardNth2gt12 + 
           gtu31*gtu33*J22L*PDstandardNth2gt13 + gtu32*gtu33*J22L*PDstandardNth2gt23 + 
           gtu31*gtu32*J32L*PDstandardNth3gt12 + gtu31*gtu33*J32L*PDstandardNth3gt13 + gtu32*gtu33*J32L*PDstandardNth3gt23)\
         - J12L*PDstandardNth1gt11*SQR(gtu31) - J22L*PDstandardNth2gt11*SQR(gtu31) - J32L*PDstandardNth3gt11*SQR(gtu31) - 
        J12L*PDstandardNth1gt22*SQR(gtu32) - J22L*PDstandardNth2gt22*SQR(gtu32) - J32L*PDstandardNth3gt22*SQR(gtu32) - 
        J12L*PDstandardNth1gt33*SQR(gtu33) - J22L*PDstandardNth2gt33*SQR(gtu33) - J32L*PDstandardNth3gt33*SQR(gtu33);
    
    CCTK_REAL const dgtu113  =  -2*(gtu11*gtu21*J13L*PDstandardNth1gt12 + gtu11*gtu31*J13L*PDstandardNth1gt13 + 
           gtu21*gtu31*J13L*PDstandardNth1gt23 + gtu11*gtu21*J23L*PDstandardNth2gt12 + 
           gtu11*gtu31*J23L*PDstandardNth2gt13 + gtu21*gtu31*J23L*PDstandardNth2gt23 + 
           gtu11*gtu21*J33L*PDstandardNth3gt12 + gtu11*gtu31*J33L*PDstandardNth3gt13 + gtu21*gtu31*J33L*PDstandardNth3gt23)\
         - J13L*PDstandardNth1gt11*SQR(gtu11) - J23L*PDstandardNth2gt11*SQR(gtu11) - J33L*PDstandardNth3gt11*SQR(gtu11) - 
        J13L*PDstandardNth1gt22*SQR(gtu21) - J23L*PDstandardNth2gt22*SQR(gtu21) - J33L*PDstandardNth3gt22*SQR(gtu21) - 
        J13L*PDstandardNth1gt33*SQR(gtu31) - J23L*PDstandardNth2gt33*SQR(gtu31) - J33L*PDstandardNth3gt33*SQR(gtu31);
    
    CCTK_REAL const dgtu213  =  -(gtu11*gtu21*J13L*PDstandardNth1gt11) - gtu11*gtu22*J13L*PDstandardNth1gt12 - 
        gtu21*gtu31*J13L*PDstandardNth1gt13 - gtu11*gtu32*J13L*PDstandardNth1gt13 - gtu21*gtu22*J13L*PDstandardNth1gt22 - 
        gtu22*gtu31*J13L*PDstandardNth1gt23 - gtu21*gtu32*J13L*PDstandardNth1gt23 - gtu31*gtu32*J13L*PDstandardNth1gt33 - 
        gtu11*gtu21*J23L*PDstandardNth2gt11 - gtu11*gtu22*J23L*PDstandardNth2gt12 - gtu21*gtu31*J23L*PDstandardNth2gt13 - 
        gtu11*gtu32*J23L*PDstandardNth2gt13 - gtu21*gtu22*J23L*PDstandardNth2gt22 - gtu22*gtu31*J23L*PDstandardNth2gt23 - 
        gtu21*gtu32*J23L*PDstandardNth2gt23 - gtu31*gtu32*J23L*PDstandardNth2gt33 - gtu11*gtu21*J33L*PDstandardNth3gt11 - 
        gtu11*gtu22*J33L*PDstandardNth3gt12 - gtu21*gtu31*J33L*PDstandardNth3gt13 - gtu11*gtu32*J33L*PDstandardNth3gt13 - 
        gtu21*gtu22*J33L*PDstandardNth3gt22 - gtu22*gtu31*J33L*PDstandardNth3gt23 - gtu21*gtu32*J33L*PDstandardNth3gt23 - 
        gtu31*gtu32*J33L*PDstandardNth3gt33 - J13L*PDstandardNth1gt12*SQR(gtu21) - J23L*PDstandardNth2gt12*SQR(gtu21) - 
        J33L*PDstandardNth3gt12*SQR(gtu21);
    
    CCTK_REAL const dgtu313  =  -(gtu11*gtu31*J13L*PDstandardNth1gt11) - gtu21*gtu31*J13L*PDstandardNth1gt12 - 
        gtu11*gtu32*J13L*PDstandardNth1gt12 - gtu11*gtu33*J13L*PDstandardNth1gt13 - gtu21*gtu32*J13L*PDstandardNth1gt22 - 
        gtu31*gtu32*J13L*PDstandardNth1gt23 - gtu21*gtu33*J13L*PDstandardNth1gt23 - gtu31*gtu33*J13L*PDstandardNth1gt33 - 
        gtu11*gtu31*J23L*PDstandardNth2gt11 - gtu21*gtu31*J23L*PDstandardNth2gt12 - gtu11*gtu32*J23L*PDstandardNth2gt12 - 
        gtu11*gtu33*J23L*PDstandardNth2gt13 - gtu21*gtu32*J23L*PDstandardNth2gt22 - gtu31*gtu32*J23L*PDstandardNth2gt23 - 
        gtu21*gtu33*J23L*PDstandardNth2gt23 - gtu31*gtu33*J23L*PDstandardNth2gt33 - gtu11*gtu31*J33L*PDstandardNth3gt11 - 
        gtu21*gtu31*J33L*PDstandardNth3gt12 - gtu11*gtu32*J33L*PDstandardNth3gt12 - gtu11*gtu33*J33L*PDstandardNth3gt13 - 
        gtu21*gtu32*J33L*PDstandardNth3gt22 - gtu31*gtu32*J33L*PDstandardNth3gt23 - gtu21*gtu33*J33L*PDstandardNth3gt23 - 
        gtu31*gtu33*J33L*PDstandardNth3gt33 - J13L*PDstandardNth1gt13*SQR(gtu31) - J23L*PDstandardNth2gt13*SQR(gtu31) - 
        J33L*PDstandardNth3gt13*SQR(gtu31);
    
    CCTK_REAL const dgtu223  =  -2*(gtu21*gtu22*J13L*PDstandardNth1gt12 + gtu21*gtu32*J13L*PDstandardNth1gt13 + 
           gtu22*gtu32*J13L*PDstandardNth1gt23 + gtu21*gtu22*J23L*PDstandardNth2gt12 + 
           gtu21*gtu32*J23L*PDstandardNth2gt13 + gtu22*gtu32*J23L*PDstandardNth2gt23 + 
           gtu21*gtu22*J33L*PDstandardNth3gt12 + gtu21*gtu32*J33L*PDstandardNth3gt13 + gtu22*gtu32*J33L*PDstandardNth3gt23)\
         - J13L*PDstandardNth1gt11*SQR(gtu21) - J23L*PDstandardNth2gt11*SQR(gtu21) - J33L*PDstandardNth3gt11*SQR(gtu21) - 
        J13L*PDstandardNth1gt22*SQR(gtu22) - J23L*PDstandardNth2gt22*SQR(gtu22) - J33L*PDstandardNth3gt22*SQR(gtu22) - 
        J13L*PDstandardNth1gt33*SQR(gtu32) - J23L*PDstandardNth2gt33*SQR(gtu32) - J33L*PDstandardNth3gt33*SQR(gtu32);
    
    CCTK_REAL const dgtu323  =  -(gtu21*gtu31*J13L*PDstandardNth1gt11) - gtu22*gtu31*J13L*PDstandardNth1gt12 - 
        gtu21*gtu32*J13L*PDstandardNth1gt12 - gtu31*gtu32*J13L*PDstandardNth1gt13 - gtu21*gtu33*J13L*PDstandardNth1gt13 - 
        gtu22*gtu32*J13L*PDstandardNth1gt22 - gtu22*gtu33*J13L*PDstandardNth1gt23 - gtu32*gtu33*J13L*PDstandardNth1gt33 - 
        gtu21*gtu31*J23L*PDstandardNth2gt11 - gtu22*gtu31*J23L*PDstandardNth2gt12 - gtu21*gtu32*J23L*PDstandardNth2gt12 - 
        gtu31*gtu32*J23L*PDstandardNth2gt13 - gtu21*gtu33*J23L*PDstandardNth2gt13 - gtu22*gtu32*J23L*PDstandardNth2gt22 - 
        gtu22*gtu33*J23L*PDstandardNth2gt23 - gtu32*gtu33*J23L*PDstandardNth2gt33 - gtu21*gtu31*J33L*PDstandardNth3gt11 - 
        gtu22*gtu31*J33L*PDstandardNth3gt12 - gtu21*gtu32*J33L*PDstandardNth3gt12 - gtu31*gtu32*J33L*PDstandardNth3gt13 - 
        gtu21*gtu33*J33L*PDstandardNth3gt13 - gtu22*gtu32*J33L*PDstandardNth3gt22 - gtu22*gtu33*J33L*PDstandardNth3gt23 - 
        gtu32*gtu33*J33L*PDstandardNth3gt33 - J13L*PDstandardNth1gt23*SQR(gtu32) - J23L*PDstandardNth2gt23*SQR(gtu32) - 
        J33L*PDstandardNth3gt23*SQR(gtu32);
    
    CCTK_REAL const dgtu333  =  -2*(gtu31*gtu32*J13L*PDstandardNth1gt12 + gtu31*gtu33*J13L*PDstandardNth1gt13 + 
           gtu32*gtu33*J13L*PDstandardNth1gt23 + gtu31*gtu32*J23L*PDstandardNth2gt12 + 
           gtu31*gtu33*J23L*PDstandardNth2gt13 + gtu32*gtu33*J23L*PDstandardNth2gt23 + 
           gtu31*gtu32*J33L*PDstandardNth3gt12 + gtu31*gtu33*J33L*PDstandardNth3gt13 + gtu32*gtu33*J33L*PDstandardNth3gt23)\
         - J13L*PDstandardNth1gt11*SQR(gtu31) - J23L*PDstandardNth2gt11*SQR(gtu31) - J33L*PDstandardNth3gt11*SQR(gtu31) - 
        J13L*PDstandardNth1gt22*SQR(gtu32) - J23L*PDstandardNth2gt22*SQR(gtu32) - J33L*PDstandardNth3gt22*SQR(gtu32) - 
        J13L*PDstandardNth1gt33*SQR(gtu33) - J23L*PDstandardNth2gt33*SQR(gtu33) - J33L*PDstandardNth3gt33*SQR(gtu33);
    
    CCTK_REAL const Gtl111  =  khalf*(J11L*PDstandardNth1gt11 + J21L*PDstandardNth2gt11 + J31L*PDstandardNth3gt11);
    
    CCTK_REAL const Gtl112  =  khalf*(J12L*PDstandardNth1gt11 + J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11);
    
    CCTK_REAL const Gtl113  =  khalf*(J13L*PDstandardNth1gt11 + J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11);
    
    CCTK_REAL const Gtl122  =  J12L*PDstandardNth1gt12 - J11L*khalf*PDstandardNth1gt22 + J22L*PDstandardNth2gt12 - 
        J21L*khalf*PDstandardNth2gt22 + J32L*PDstandardNth3gt12 - J31L*khalf*PDstandardNth3gt22;
    
    CCTK_REAL const Gtl123  =  khalf*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - J11L*PDstandardNth1gt23 + 
          J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 + 
          J32L*PDstandardNth3gt13 - J31L*PDstandardNth3gt23);
    
    CCTK_REAL const Gtl133  =  J13L*PDstandardNth1gt13 - J11L*khalf*PDstandardNth1gt33 + J23L*PDstandardNth2gt13 - 
        J21L*khalf*PDstandardNth2gt33 + J33L*PDstandardNth3gt13 - J31L*khalf*PDstandardNth3gt33;
    
    CCTK_REAL const Gtl211  =  khalf*(-(J12L*PDstandardNth1gt11) + 2*J11L*PDstandardNth1gt12 - J22L*PDstandardNth2gt11 + 
          2*J21L*PDstandardNth2gt12 - J32L*PDstandardNth3gt11 + 2*J31L*PDstandardNth3gt12);
    
    CCTK_REAL const Gtl212  =  khalf*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22);
    
    CCTK_REAL const Gtl213  =  khalf*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 + 
          J23L*PDstandardNth2gt12 - J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 - 
          J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23);
    
    CCTK_REAL const Gtl222  =  khalf*(J12L*PDstandardNth1gt22 + J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22);
    
    CCTK_REAL const Gtl223  =  khalf*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22);
    
    CCTK_REAL const Gtl233  =  J13L*PDstandardNth1gt23 - J12L*khalf*PDstandardNth1gt33 + J23L*PDstandardNth2gt23 - 
        J22L*khalf*PDstandardNth2gt33 + J33L*PDstandardNth3gt23 - J32L*khalf*PDstandardNth3gt33;
    
    CCTK_REAL const Gtl311  =  khalf*(-(J13L*PDstandardNth1gt11) + 2*J11L*PDstandardNth1gt13 - J23L*PDstandardNth2gt11 + 
          2*J21L*PDstandardNth2gt13 - J33L*PDstandardNth3gt11 + 2*J31L*PDstandardNth3gt13);
    
    CCTK_REAL const Gtl312  =  khalf*(-(J13L*PDstandardNth1gt12) + J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
          J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
          J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23);
    
    CCTK_REAL const Gtl313  =  khalf*(J11L*PDstandardNth1gt33 + J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33);
    
    CCTK_REAL const Gtl322  =  khalf*(-(J13L*PDstandardNth1gt22) + 2*J12L*PDstandardNth1gt23 - J23L*PDstandardNth2gt22 + 
          2*J22L*PDstandardNth2gt23 - J33L*PDstandardNth3gt22 + 2*J32L*PDstandardNth3gt23);
    
    CCTK_REAL const Gtl323  =  khalf*(J12L*PDstandardNth1gt33 + J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33);
    
    CCTK_REAL const Gtl333  =  khalf*(J13L*PDstandardNth1gt33 + J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33);
    
    CCTK_REAL const Gtlu111  =  Gtl111*gtu11 + Gtl112*gtu21 + Gtl113*gtu31;
    
    CCTK_REAL const Gtlu112  =  Gtl111*gtu21 + Gtl112*gtu22 + Gtl113*gtu32;
    
    CCTK_REAL const Gtlu113  =  Gtl111*gtu31 + Gtl112*gtu32 + Gtl113*gtu33;
    
    CCTK_REAL const Gtlu121  =  Gtl112*gtu11 + Gtl122*gtu21 + Gtl123*gtu31;
    
    CCTK_REAL const Gtlu122  =  Gtl112*gtu21 + Gtl122*gtu22 + Gtl123*gtu32;
    
    CCTK_REAL const Gtlu123  =  Gtl112*gtu31 + Gtl122*gtu32 + Gtl123*gtu33;
    
    CCTK_REAL const Gtlu131  =  Gtl113*gtu11 + Gtl123*gtu21 + Gtl133*gtu31;
    
    CCTK_REAL const Gtlu132  =  Gtl113*gtu21 + Gtl123*gtu22 + Gtl133*gtu32;
    
    CCTK_REAL const Gtlu133  =  Gtl113*gtu31 + Gtl123*gtu32 + Gtl133*gtu33;
    
    CCTK_REAL const Gtlu211  =  Gtl211*gtu11 + Gtl212*gtu21 + Gtl213*gtu31;
    
    CCTK_REAL const Gtlu212  =  Gtl211*gtu21 + Gtl212*gtu22 + Gtl213*gtu32;
    
    CCTK_REAL const Gtlu213  =  Gtl211*gtu31 + Gtl212*gtu32 + Gtl213*gtu33;
    
    CCTK_REAL const Gtlu221  =  Gtl212*gtu11 + Gtl222*gtu21 + Gtl223*gtu31;
    
    CCTK_REAL const Gtlu222  =  Gtl212*gtu21 + Gtl222*gtu22 + Gtl223*gtu32;
    
    CCTK_REAL const Gtlu223  =  Gtl212*gtu31 + Gtl222*gtu32 + Gtl223*gtu33;
    
    CCTK_REAL const Gtlu231  =  Gtl213*gtu11 + Gtl223*gtu21 + Gtl233*gtu31;
    
    CCTK_REAL const Gtlu232  =  Gtl213*gtu21 + Gtl223*gtu22 + Gtl233*gtu32;
    
    CCTK_REAL const Gtlu233  =  Gtl213*gtu31 + Gtl223*gtu32 + Gtl233*gtu33;
    
    CCTK_REAL const Gtlu311  =  Gtl311*gtu11 + Gtl312*gtu21 + Gtl313*gtu31;
    
    CCTK_REAL const Gtlu312  =  Gtl311*gtu21 + Gtl312*gtu22 + Gtl313*gtu32;
    
    CCTK_REAL const Gtlu313  =  Gtl311*gtu31 + Gtl312*gtu32 + Gtl313*gtu33;
    
    CCTK_REAL const Gtlu321  =  Gtl312*gtu11 + Gtl322*gtu21 + Gtl323*gtu31;
    
    CCTK_REAL const Gtlu322  =  Gtl312*gtu21 + Gtl322*gtu22 + Gtl323*gtu32;
    
    CCTK_REAL const Gtlu323  =  Gtl312*gtu31 + Gtl322*gtu32 + Gtl323*gtu33;
    
    CCTK_REAL const Gtlu331  =  Gtl313*gtu11 + Gtl323*gtu21 + Gtl333*gtu31;
    
    CCTK_REAL const Gtlu332  =  Gtl313*gtu21 + Gtl323*gtu22 + Gtl333*gtu32;
    
    CCTK_REAL const Gtlu333  =  Gtl313*gtu31 + Gtl323*gtu32 + Gtl333*gtu33;
    
    CCTK_REAL const Gt111  =  Gtl111*gtu11 + Gtl211*gtu21 + Gtl311*gtu31;
    
    CCTK_REAL const Gt211  =  Gtl111*gtu21 + Gtl211*gtu22 + Gtl311*gtu32;
    
    CCTK_REAL const Gt311  =  Gtl111*gtu31 + Gtl211*gtu32 + Gtl311*gtu33;
    
    CCTK_REAL const Gt112  =  Gtl112*gtu11 + Gtl212*gtu21 + Gtl312*gtu31;
    
    CCTK_REAL const Gt212  =  Gtl112*gtu21 + Gtl212*gtu22 + Gtl312*gtu32;
    
    CCTK_REAL const Gt312  =  Gtl112*gtu31 + Gtl212*gtu32 + Gtl312*gtu33;
    
    CCTK_REAL const Gt113  =  Gtl113*gtu11 + Gtl213*gtu21 + Gtl313*gtu31;
    
    CCTK_REAL const Gt213  =  Gtl113*gtu21 + Gtl213*gtu22 + Gtl313*gtu32;
    
    CCTK_REAL const Gt313  =  Gtl113*gtu31 + Gtl213*gtu32 + Gtl313*gtu33;
    
    CCTK_REAL const Gt122  =  Gtl122*gtu11 + Gtl222*gtu21 + Gtl322*gtu31;
    
    CCTK_REAL const Gt222  =  Gtl122*gtu21 + Gtl222*gtu22 + Gtl322*gtu32;
    
    CCTK_REAL const Gt322  =  Gtl122*gtu31 + Gtl222*gtu32 + Gtl322*gtu33;
    
    CCTK_REAL const Gt123  =  Gtl123*gtu11 + Gtl223*gtu21 + Gtl323*gtu31;
    
    CCTK_REAL const Gt223  =  Gtl123*gtu21 + Gtl223*gtu22 + Gtl323*gtu32;
    
    CCTK_REAL const Gt323  =  Gtl123*gtu31 + Gtl223*gtu32 + Gtl323*gtu33;
    
    CCTK_REAL const Gt133  =  Gtl133*gtu11 + Gtl233*gtu21 + Gtl333*gtu31;
    
    CCTK_REAL const Gt233  =  Gtl133*gtu21 + Gtl233*gtu22 + Gtl333*gtu32;
    
    CCTK_REAL const Gt333  =  Gtl133*gtu31 + Gtl233*gtu32 + Gtl333*gtu33;
    
    CCTK_REAL const Xtn1  =  Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    CCTK_REAL const Xtn2  =  Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    CCTK_REAL const Xtn3  =  Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    CCTK_REAL const Rt11  =  3*(Gt111*Gtlu111 + Gt112*Gtlu112 + Gt113*Gtlu113) + 
        2*(Gt211*Gtlu121 + Gt212*Gtlu122 + Gt213*Gtlu123 + Gt311*Gtlu131 + Gt312*Gtlu132 + Gt313*Gtlu133) + Gt211*Gtlu211 + 
        Gt212*Gtlu212 + Gt213*Gtlu213 + Gt311*Gtlu311 + Gt312*Gtlu312 + Gt313*Gtlu313 + 
        J11L*(gt11L*PDstandardNth1Xt1 + gt12L*PDstandardNth1Xt2 + gt13L*PDstandardNth1Xt3) + 
        J21L*(gt11L*PDstandardNth2Xt1 + gt12L*PDstandardNth2Xt2 + gt13L*PDstandardNth2Xt3) + 
        J31L*(gt11L*PDstandardNth3Xt1 + gt12L*PDstandardNth3Xt2 + gt13L*PDstandardNth3Xt3) + Gtl111*Xtn1 + Gtl112*Xtn2 + 
        Gtl113*Xtn3 + khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt11 + J12L*J21L*PDstandardNth12gt11 + 
                 J11L*J22L*PDstandardNth12gt11 + J12L*J31L*PDstandardNth13gt11 + J11L*J32L*PDstandardNth13gt11 + 
                 dJ112L*PDstandardNth1gt11 + J21L*J22L*PDstandardNth22gt11 + J22L*J31L*PDstandardNth23gt11 + 
                 J21L*J32L*PDstandardNth23gt11 + dJ212L*PDstandardNth2gt11 + J31L*J32L*PDstandardNth33gt11 + 
                 dJ312L*PDstandardNth3gt11) + gtu31*(J11L*J13L*PDstandardNth11gt11 + J13L*J21L*PDstandardNth12gt11 + 
                 J11L*J23L*PDstandardNth12gt11 + J13L*J31L*PDstandardNth13gt11 + J11L*J33L*PDstandardNth13gt11 + 
                 dJ113L*PDstandardNth1gt11 + J21L*J23L*PDstandardNth22gt11 + J23L*J31L*PDstandardNth23gt11 + 
                 J21L*J33L*PDstandardNth23gt11 + dJ213L*PDstandardNth2gt11 + J31L*J33L*PDstandardNth33gt11 + 
                 dJ313L*PDstandardNth3gt11) + gtu32*(J12L*J13L*PDstandardNth11gt11 + J13L*J22L*PDstandardNth12gt11 + 
                 J12L*J23L*PDstandardNth12gt11 + J13L*J32L*PDstandardNth13gt11 + J12L*J33L*PDstandardNth13gt11 + 
                 dJ123L*PDstandardNth1gt11 + J22L*J23L*PDstandardNth22gt11 + J23L*J32L*PDstandardNth23gt11 + 
                 J22L*J33L*PDstandardNth23gt11 + dJ223L*PDstandardNth2gt11 + J32L*J33L*PDstandardNth33gt11 + 
                 dJ323L*PDstandardNth3gt11)) - gtu11*
            (2*J11L*J21L*PDstandardNth12gt11 + 2*J11L*J31L*PDstandardNth13gt11 + dJ111L*PDstandardNth1gt11 + 
              2*J21L*J31L*PDstandardNth23gt11 + dJ211L*PDstandardNth2gt11 + dJ311L*PDstandardNth3gt11 + 
              PDstandardNth11gt11*SQR(J11L) + PDstandardNth22gt11*SQR(J21L) + PDstandardNth33gt11*SQR(J31L)) - 
           gtu22*(2*J12L*J22L*PDstandardNth12gt11 + 2*J12L*J32L*PDstandardNth13gt11 + dJ122L*PDstandardNth1gt11 + 
              2*J22L*J32L*PDstandardNth23gt11 + dJ222L*PDstandardNth2gt11 + dJ322L*PDstandardNth3gt11 + 
              PDstandardNth11gt11*SQR(J12L) + PDstandardNth22gt11*SQR(J22L) + PDstandardNth33gt11*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt11 + 2*J13L*J33L*PDstandardNth13gt11 + dJ133L*PDstandardNth1gt11 + 
              2*J23L*J33L*PDstandardNth23gt11 + dJ233L*PDstandardNth2gt11 + dJ333L*PDstandardNth3gt11 + 
              PDstandardNth11gt11*SQR(J13L) + PDstandardNth22gt11*SQR(J23L) + PDstandardNth33gt11*SQR(J33L)));
    
    CCTK_REAL const Rt12  =  Gt122*Gtlu112 + Gt123*Gtlu113 + (Gt111 + Gt212)*Gtlu121 + Gt222*Gtlu122 + (Gt113 + Gt223)*Gtlu123 + 
        Gt322*Gtlu132 + Gt323*Gtlu133 + Gt111*Gtlu211 + Gt112*(Gtlu111 + Gtlu122 + Gtlu212) + Gt113*Gtlu213 + 
        2*(Gt211*Gtlu221 + Gt212*Gtlu222 + Gt213*Gtlu223) + Gt311*(Gtlu231 + Gtlu321) + 
        Gt312*(Gtlu131 + Gtlu232 + Gtlu322) + Gt313*(Gtlu233 + Gtlu323) + 
        khalf*((gt12L*J11L + gt11L*J12L)*PDstandardNth1Xt1 + (gt22L*J11L + gt12L*J12L)*PDstandardNth1Xt2 + 
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
           (gt23L*J31L + gt13L*J32L)*PDstandardNth3Xt3 + (Gtl112 + Gtl211)*Xtn1 + (Gtl122 + Gtl212)*Xtn2 + 
           (Gtl123 + Gtl213)*Xtn3 - gtu11*(2*J11L*J21L*PDstandardNth12gt12 + 2*J11L*J31L*PDstandardNth13gt12 + 
              dJ111L*PDstandardNth1gt12 + 2*J21L*J31L*PDstandardNth23gt12 + dJ211L*PDstandardNth2gt12 + 
              dJ311L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J11L) + PDstandardNth22gt12*SQR(J21L) + 
              PDstandardNth33gt12*SQR(J31L)) - gtu22*
            (2*J12L*J22L*PDstandardNth12gt12 + 2*J12L*J32L*PDstandardNth13gt12 + dJ122L*PDstandardNth1gt12 + 
              2*J22L*J32L*PDstandardNth23gt12 + dJ222L*PDstandardNth2gt12 + dJ322L*PDstandardNth3gt12 + 
              PDstandardNth11gt12*SQR(J12L) + PDstandardNth22gt12*SQR(J22L) + PDstandardNth33gt12*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt12 + 2*J13L*J33L*PDstandardNth13gt12 + dJ133L*PDstandardNth1gt12 + 
              2*J23L*J33L*PDstandardNth23gt12 + dJ233L*PDstandardNth2gt12 + dJ333L*PDstandardNth3gt12 + 
              PDstandardNth11gt12*SQR(J13L) + PDstandardNth22gt12*SQR(J23L) + PDstandardNth33gt12*SQR(J33L)));
    
    CCTK_REAL const Rt13  =  Gt123*Gtlu112 + Gt133*Gtlu113 + Gt223*Gtlu122 + Gt233*Gtlu123 + (Gt111 + Gt313)*Gtlu131 + 
        (Gt112 + Gt323)*Gtlu132 + Gt333*Gtlu133 + Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*(Gtlu111 + Gtlu133 + Gtlu313) + 
        Gt211*(Gtlu231 + Gtlu321) + Gt212*(Gtlu232 + Gtlu322) + Gt213*(Gtlu121 + Gtlu233 + Gtlu323) + 
        2*(Gt311*Gtlu331 + Gt312*Gtlu332 + Gt313*Gtlu333) + 
        khalf*((gt13L*J11L + gt11L*J13L)*PDstandardNth1Xt1 + (gt23L*J11L + gt12L*J13L)*PDstandardNth1Xt2 + 
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
           (gt33L*J31L + gt13L*J33L)*PDstandardNth3Xt3 + (Gtl113 + Gtl311)*Xtn1 + (Gtl123 + Gtl312)*Xtn2 + 
           (Gtl133 + Gtl313)*Xtn3 - gtu11*(2*J11L*J21L*PDstandardNth12gt13 + 2*J11L*J31L*PDstandardNth13gt13 + 
              dJ111L*PDstandardNth1gt13 + 2*J21L*J31L*PDstandardNth23gt13 + dJ211L*PDstandardNth2gt13 + 
              dJ311L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J11L) + PDstandardNth22gt13*SQR(J21L) + 
              PDstandardNth33gt13*SQR(J31L)) - gtu22*
            (2*J12L*J22L*PDstandardNth12gt13 + 2*J12L*J32L*PDstandardNth13gt13 + dJ122L*PDstandardNth1gt13 + 
              2*J22L*J32L*PDstandardNth23gt13 + dJ222L*PDstandardNth2gt13 + dJ322L*PDstandardNth3gt13 + 
              PDstandardNth11gt13*SQR(J12L) + PDstandardNth22gt13*SQR(J22L) + PDstandardNth33gt13*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt13 + 2*J13L*J33L*PDstandardNth13gt13 + dJ133L*PDstandardNth1gt13 + 
              2*J23L*J33L*PDstandardNth23gt13 + dJ233L*PDstandardNth2gt13 + dJ333L*PDstandardNth3gt13 + 
              PDstandardNth11gt13*SQR(J13L) + PDstandardNth22gt13*SQR(J23L) + PDstandardNth33gt13*SQR(J33L)));
    
    CCTK_REAL const Rt22  =  Gt112*(Gtlu121 + 2*Gtlu211) + Gt122*(Gtlu122 + 2*Gtlu212) + Gt123*(Gtlu123 + 2*Gtlu213) + 
        3*(Gt212*Gtlu221 + Gt222*Gtlu222 + Gt223*Gtlu223) + 2*(Gt312*Gtlu231 + Gt322*Gtlu232 + Gt323*Gtlu233) + 
        Gt312*Gtlu321 + Gt322*Gtlu322 + Gt323*Gtlu323 + 
        J12L*(gt12L*PDstandardNth1Xt1 + gt22L*PDstandardNth1Xt2 + gt23L*PDstandardNth1Xt3) + 
        J22L*(gt12L*PDstandardNth2Xt1 + gt22L*PDstandardNth2Xt2 + gt23L*PDstandardNth2Xt3) + 
        J32L*(gt12L*PDstandardNth3Xt1 + gt22L*PDstandardNth3Xt2 + gt23L*PDstandardNth3Xt3) + Gtl212*Xtn1 + Gtl222*Xtn2 + 
        Gtl223*Xtn3 + khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt22 + J12L*J21L*PDstandardNth12gt22 + 
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
    
    CCTK_REAL const Rt23  =  Gt113*Gtlu211 + Gt133*Gtlu213 + Gt213*Gtlu221 + Gt233*Gtlu223 + (Gt212 + Gt313)*Gtlu231 + 
        (Gt222 + Gt323)*Gtlu232 + Gt333*Gtlu233 + Gt112*(Gtlu131 + Gtlu311) + Gt122*(Gtlu132 + Gtlu312) + 
        Gt123*(Gtlu133 + Gtlu212 + Gtlu313) + Gt212*Gtlu321 + Gt222*Gtlu322 + Gt223*(Gtlu222 + Gtlu233 + Gtlu323) + 
        2*(Gt312*Gtlu331 + Gt322*Gtlu332 + Gt323*Gtlu333) + 
        khalf*((gt13L*J12L + gt12L*J13L)*PDstandardNth1Xt1 + (gt23L*J12L + gt22L*J13L)*PDstandardNth1Xt2 + 
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
           (gt33L*J32L + gt23L*J33L)*PDstandardNth3Xt3 + (Gtl213 + Gtl312)*Xtn1 + (Gtl223 + Gtl322)*Xtn2 + 
           (Gtl233 + Gtl323)*Xtn3 - gtu11*(2*J11L*J21L*PDstandardNth12gt23 + 2*J11L*J31L*PDstandardNth13gt23 + 
              dJ111L*PDstandardNth1gt23 + 2*J21L*J31L*PDstandardNth23gt23 + dJ211L*PDstandardNth2gt23 + 
              dJ311L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J11L) + PDstandardNth22gt23*SQR(J21L) + 
              PDstandardNth33gt23*SQR(J31L)) - gtu22*
            (2*J12L*J22L*PDstandardNth12gt23 + 2*J12L*J32L*PDstandardNth13gt23 + dJ122L*PDstandardNth1gt23 + 
              2*J22L*J32L*PDstandardNth23gt23 + dJ222L*PDstandardNth2gt23 + dJ322L*PDstandardNth3gt23 + 
              PDstandardNth11gt23*SQR(J12L) + PDstandardNth22gt23*SQR(J22L) + PDstandardNth33gt23*SQR(J32L)) - 
           gtu33*(2*J13L*J23L*PDstandardNth12gt23 + 2*J13L*J33L*PDstandardNth13gt23 + dJ133L*PDstandardNth1gt23 + 
              2*J23L*J33L*PDstandardNth23gt23 + dJ233L*PDstandardNth2gt23 + dJ333L*PDstandardNth3gt23 + 
              PDstandardNth11gt23*SQR(J13L) + PDstandardNth22gt23*SQR(J23L) + PDstandardNth33gt23*SQR(J33L)));
    
    CCTK_REAL const Rt33  =  Gt113*(Gtlu131 + 2*Gtlu311) + Gt123*(Gtlu132 + 2*Gtlu312) + Gt133*(Gtlu133 + 2*Gtlu313) + 
        Gt213*(Gtlu231 + 2*Gtlu321) + Gt223*(Gtlu232 + 2*Gtlu322) + Gt233*(Gtlu233 + 2*Gtlu323) + 
        3*(Gt313*Gtlu331 + Gt323*Gtlu332 + Gt333*Gtlu333) + 
        J13L*(gt13L*PDstandardNth1Xt1 + gt23L*PDstandardNth1Xt2 + gt33L*PDstandardNth1Xt3) + 
        J23L*(gt13L*PDstandardNth2Xt1 + gt23L*PDstandardNth2Xt2 + gt33L*PDstandardNth2Xt3) + 
        J33L*(gt13L*PDstandardNth3Xt1 + gt23L*PDstandardNth3Xt2 + gt33L*PDstandardNth3Xt3) + Gtl313*Xtn1 + Gtl323*Xtn2 + 
        Gtl333*Xtn3 + khalf*(-2*(gtu21*(J11L*J12L*PDstandardNth11gt33 + J12L*J21L*PDstandardNth12gt33 + 
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
    
    CCTK_REAL const trRt  =  gtu11*Rt11 + gtu22*Rt22 + 2*(gtu21*Rt12 + gtu31*Rt13 + gtu32*Rt23) + gtu33*Rt33;
    
    CCTK_REAL const ephi  =  IfThen(conformalMethod,pow(phiL,-khalf),exp(phiL));
    
    CCTK_REAL const Atm11  =  At11L*gtu11 + At12L*gtu21 + At13L*gtu31;
    
    CCTK_REAL const Atm21  =  At11L*gtu21 + At12L*gtu22 + At13L*gtu32;
    
    CCTK_REAL const Atm31  =  At11L*gtu31 + At12L*gtu32 + At13L*gtu33;
    
    CCTK_REAL const Atm12  =  At12L*gtu11 + At22L*gtu21 + At23L*gtu31;
    
    CCTK_REAL const Atm22  =  At12L*gtu21 + At22L*gtu22 + At23L*gtu32;
    
    CCTK_REAL const Atm32  =  At12L*gtu31 + At22L*gtu32 + At23L*gtu33;
    
    CCTK_REAL const Atm13  =  At13L*gtu11 + At23L*gtu21 + At33L*gtu31;
    
    CCTK_REAL const Atm23  =  At13L*gtu21 + At23L*gtu22 + At33L*gtu32;
    
    CCTK_REAL const Atm33  =  At13L*gtu31 + At23L*gtu32 + At33L*gtu33;
    
    CCTK_REAL const rho  =  pow(alphaL,-2)*(eTttL - 2*(beta2L*eTtyL + beta3L*eTtzL) + 
          2*(beta1L*(-eTtxL + beta2L*eTxyL + beta3L*eTxzL) + beta2L*beta3L*eTyzL) + eTxxL*SQR(beta1L) + eTyyL*SQR(beta2L) + 
          eTzzL*SQR(beta3L));
    
    CCTK_REAL const S1  =  (-eTtxL + beta1L*eTxxL + beta2L*eTxyL + beta3L*eTxzL)*INV(alphaL);
    
    CCTK_REAL const S2  =  (-eTtyL + beta1L*eTxyL + beta2L*eTyyL + beta3L*eTyzL)*INV(alphaL);
    
    CCTK_REAL const S3  =  (-eTtzL + beta1L*eTxzL + beta2L*eTyzL + beta3L*eTzzL)*INV(alphaL);
    
    CCTK_REAL const MadmL  =  0.01989436788648691697111047042156429525431*
        (-((Gt111*Gtlu111 + Gt112*Gtlu112 + Gt113*Gtlu113 + Gt211*Gtlu121 + Gt212*Gtlu122 + Gt213*Gtlu123 + Gt311*Gtlu131 + 
               Gt312*Gtlu132 + Gt313*Gtlu133)*gtu11) - 
          (Gt112*Gtlu111 + Gt122*Gtlu112 + Gt123*Gtlu113 + Gt212*Gtlu121 + Gt222*Gtlu122 + Gt223*Gtlu123 + Gt312*Gtlu131 + 
             Gt322*Gtlu132 + Gt323*Gtlu133)*gtu21 - (Gt111*Gtlu211 + Gt112*Gtlu212 + Gt113*Gtlu213 + Gt211*Gtlu221 + 
             Gt212*Gtlu222 + Gt213*Gtlu223 + Gt311*Gtlu231 + Gt312*Gtlu232 + Gt313*Gtlu233)*gtu21 - 
          (Gt112*Gtlu211 + Gt122*Gtlu212 + Gt123*Gtlu213 + Gt212*Gtlu221 + Gt222*Gtlu222 + Gt223*Gtlu223 + Gt312*Gtlu231 + 
             Gt322*Gtlu232 + Gt323*Gtlu233)*gtu22 - (Gt113*Gtlu111 + Gt123*Gtlu112 + Gt133*Gtlu113 + Gt213*Gtlu121 + 
             Gt223*Gtlu122 + Gt233*Gtlu123 + Gt313*Gtlu131 + Gt323*Gtlu132 + Gt333*Gtlu133)*gtu31 - 
          (Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*Gtlu313 + Gt211*Gtlu321 + Gt212*Gtlu322 + Gt213*Gtlu323 + Gt311*Gtlu331 + 
             Gt312*Gtlu332 + Gt313*Gtlu333)*gtu31 - (Gt113*Gtlu211 + Gt123*Gtlu212 + Gt133*Gtlu213 + Gt213*Gtlu221 + 
             Gt223*Gtlu222 + Gt233*Gtlu223 + Gt313*Gtlu231 + Gt323*Gtlu232 + Gt333*Gtlu233)*gtu32 - 
          (Gt112*Gtlu311 + Gt122*Gtlu312 + Gt123*Gtlu313 + Gt212*Gtlu321 + Gt222*Gtlu322 + Gt223*Gtlu323 + Gt312*Gtlu331 + 
             Gt322*Gtlu332 + Gt323*Gtlu333)*gtu32 - (Gt113*Gtlu311 + Gt123*Gtlu312 + Gt133*Gtlu313 + Gt213*Gtlu321 + 
             Gt223*Gtlu322 + Gt233*Gtlu323 + Gt313*Gtlu331 + Gt323*Gtlu332 + Gt333*Gtlu333)*gtu33 + trRt - ephi*trRt + 
          pow(ephi,5)*(2*Atm12*Atm21 + 2.*Atm13*Atm31 + 2.*Atm23*Atm32 + 50.26548245743669181540229413247204614715*rho + 
             SQR(Atm11) + SQR(Atm22) + SQR(Atm33) - 0.6666666666666666666666666666666666666667*SQR(trKL)));
    
    CCTK_REAL const Jadm1L  =  0.01989436788648691697111047042156429525431*
        (2*Atm23 - 2*Atm32 + (-(At11L*dgtu113) - 2*At12L*dgtu213 - At22L*dgtu223 - 2*At13L*dgtu313 - 2*At23L*dgtu323 - 
             At33L*dgtu333 + kfourthird*(J13L*PDstandardNth1trK + J23L*PDstandardNth2trK + J33L*PDstandardNth3trK) + 
             50.26548245743669181540229413247204614715*S3)*yL + 
          (At11L*dgtu112 + At22L*dgtu222 + 2*(At12L*dgtu212 + At13L*dgtu312 + At23L*dgtu322) + At33L*dgtu332 - 
             kfourthird*(J12L*PDstandardNth1trK + J22L*PDstandardNth2trK + J32L*PDstandardNth3trK) - 
             50.26548245743669181540229413247204614715*S2)*zL)*pow(ephi,6);
    
    CCTK_REAL const Jadm2L  =  0.01989436788648691697111047042156429525431*
        (-2*Atm13 + 2*Atm31 + (At11L*dgtu113 + At22L*dgtu223 + 2*(At12L*dgtu213 + At13L*dgtu313 + At23L*dgtu323) + 
             At33L*dgtu333 - kfourthird*(J13L*PDstandardNth1trK + J23L*PDstandardNth2trK + J33L*PDstandardNth3trK) - 
             50.26548245743669181540229413247204614715*S3)*xL + 
          (-(At11L*dgtu111) - 2*At12L*dgtu211 - At22L*dgtu221 - 2*At13L*dgtu311 - 2*At23L*dgtu321 - At33L*dgtu331 + 
             kfourthird*(J11L*PDstandardNth1trK + J21L*PDstandardNth2trK + J31L*PDstandardNth3trK) + 
             50.26548245743669181540229413247204614715*S1)*zL)*pow(ephi,6);
    
    CCTK_REAL const Jadm3L  =  0.01989436788648691697111047042156429525431*
        (2*Atm12 - 2*Atm21 + (-(At11L*dgtu112) - 2*At12L*dgtu212 - At22L*dgtu222 - 2*At13L*dgtu312 - 2*At23L*dgtu322 - 
             At33L*dgtu332 + kfourthird*(J12L*PDstandardNth1trK + J22L*PDstandardNth2trK + J32L*PDstandardNth3trK) + 
             50.26548245743669181540229413247204614715*S2)*xL + 
          (At11L*dgtu111 + At22L*dgtu221 + 2*(At12L*dgtu211 + At13L*dgtu311 + At23L*dgtu321) + At33L*dgtu331 - 
             kfourthird*(J11L*PDstandardNth1trK + J21L*PDstandardNth2trK + J31L*PDstandardNth3trK) - 
             50.26548245743669181540229413247204614715*S1)*yL)*pow(ephi,6);
    
    
    /* Copy local copies back to grid functions */
    Jadm1[index] = Jadm1L;
    Jadm2[index] = Jadm2L;
    Jadm3[index] = Jadm3L;
    Madm[index] = MadmL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_ADMQuantities_MP);
}

void ML_ADMQuantities_MP(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADMQuantities_MP_Body);
}
