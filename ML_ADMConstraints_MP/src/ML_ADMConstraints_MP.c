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

void ML_ADMConstraints_MP_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMConstraints_MP_Body");
  }
  
  if (cctk_iteration % ML_ADMConstraints_MP_calc_every != ML_ADMConstraints_MP_calc_offset)
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
  LC_LOOP3 (ML_ADMConstraints_MP,
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
    // CCTK_REAL gu11 = INITVALUE, gu21 = INITVALUE, gu22 = INITVALUE, gu31 = INITVALUE, gu32 = INITVALUE, gu33 = INITVALUE;
    // CCTK_REAL Km11 = INITVALUE, Km12 = INITVALUE, Km13 = INITVALUE, Km21 = INITVALUE, Km22 = INITVALUE, Km23 = INITVALUE;
    // CCTK_REAL Km31 = INITVALUE, Km32 = INITVALUE, Km33 = INITVALUE;
    // CCTK_REAL R11 = INITVALUE, R12 = INITVALUE, R13 = INITVALUE, R22 = INITVALUE, R23 = INITVALUE, R33 = INITVALUE;
    // CCTK_REAL rho = INITVALUE;
    // CCTK_REAL S1 = INITVALUE, S2 = INITVALUE, S3 = INITVALUE;
    // CCTK_REAL trK = INITVALUE;
    // CCTK_REAL trR = INITVALUE;
    
    /* Declare local copies of grid functions */
    // CCTK_REAL alpL = INITVALUE;
    // CCTK_REAL betaxL = INITVALUE;
    // CCTK_REAL betayL = INITVALUE;
    // CCTK_REAL betazL = INITVALUE;
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
    // CCTK_REAL gxxL = INITVALUE;
    // CCTK_REAL gxyL = INITVALUE;
    // CCTK_REAL gxzL = INITVALUE;
    // CCTK_REAL gyyL = INITVALUE;
    // CCTK_REAL gyzL = INITVALUE;
    // CCTK_REAL gzzL = INITVALUE;
    // CCTK_REAL HL = INITVALUE;
    // CCTK_REAL J11L = INITVALUE, J12L = INITVALUE, J13L = INITVALUE, J21L = INITVALUE, J22L = INITVALUE, J23L = INITVALUE;
    // CCTK_REAL J31L = INITVALUE, J32L = INITVALUE, J33L = INITVALUE;
    // CCTK_REAL kxxL = INITVALUE;
    // CCTK_REAL kxyL = INITVALUE;
    // CCTK_REAL kxzL = INITVALUE;
    // CCTK_REAL kyyL = INITVALUE;
    // CCTK_REAL kyzL = INITVALUE;
    // CCTK_REAL kzzL = INITVALUE;
    // CCTK_REAL M1L = INITVALUE, M2L = INITVALUE, M3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    // CCTK_REAL PDstandardNth1gxx = INITVALUE;
    // CCTK_REAL PDstandardNth2gxx = INITVALUE;
    // CCTK_REAL PDstandardNth3gxx = INITVALUE;
    // CCTK_REAL PDstandardNth11gxx = INITVALUE;
    // CCTK_REAL PDstandardNth22gxx = INITVALUE;
    // CCTK_REAL PDstandardNth33gxx = INITVALUE;
    // CCTK_REAL PDstandardNth12gxx = INITVALUE;
    // CCTK_REAL PDstandardNth13gxx = INITVALUE;
    // CCTK_REAL PDstandardNth23gxx = INITVALUE;
    // CCTK_REAL PDstandardNth1gxy = INITVALUE;
    // CCTK_REAL PDstandardNth2gxy = INITVALUE;
    // CCTK_REAL PDstandardNth3gxy = INITVALUE;
    // CCTK_REAL PDstandardNth11gxy = INITVALUE;
    // CCTK_REAL PDstandardNth22gxy = INITVALUE;
    // CCTK_REAL PDstandardNth33gxy = INITVALUE;
    // CCTK_REAL PDstandardNth12gxy = INITVALUE;
    // CCTK_REAL PDstandardNth13gxy = INITVALUE;
    // CCTK_REAL PDstandardNth23gxy = INITVALUE;
    // CCTK_REAL PDstandardNth1gxz = INITVALUE;
    // CCTK_REAL PDstandardNth2gxz = INITVALUE;
    // CCTK_REAL PDstandardNth3gxz = INITVALUE;
    // CCTK_REAL PDstandardNth11gxz = INITVALUE;
    // CCTK_REAL PDstandardNth22gxz = INITVALUE;
    // CCTK_REAL PDstandardNth33gxz = INITVALUE;
    // CCTK_REAL PDstandardNth12gxz = INITVALUE;
    // CCTK_REAL PDstandardNth13gxz = INITVALUE;
    // CCTK_REAL PDstandardNth23gxz = INITVALUE;
    // CCTK_REAL PDstandardNth1gyy = INITVALUE;
    // CCTK_REAL PDstandardNth2gyy = INITVALUE;
    // CCTK_REAL PDstandardNth3gyy = INITVALUE;
    // CCTK_REAL PDstandardNth11gyy = INITVALUE;
    // CCTK_REAL PDstandardNth22gyy = INITVALUE;
    // CCTK_REAL PDstandardNth33gyy = INITVALUE;
    // CCTK_REAL PDstandardNth12gyy = INITVALUE;
    // CCTK_REAL PDstandardNth13gyy = INITVALUE;
    // CCTK_REAL PDstandardNth23gyy = INITVALUE;
    // CCTK_REAL PDstandardNth1gyz = INITVALUE;
    // CCTK_REAL PDstandardNth2gyz = INITVALUE;
    // CCTK_REAL PDstandardNth3gyz = INITVALUE;
    // CCTK_REAL PDstandardNth11gyz = INITVALUE;
    // CCTK_REAL PDstandardNth22gyz = INITVALUE;
    // CCTK_REAL PDstandardNth33gyz = INITVALUE;
    // CCTK_REAL PDstandardNth12gyz = INITVALUE;
    // CCTK_REAL PDstandardNth13gyz = INITVALUE;
    // CCTK_REAL PDstandardNth23gyz = INITVALUE;
    // CCTK_REAL PDstandardNth1gzz = INITVALUE;
    // CCTK_REAL PDstandardNth2gzz = INITVALUE;
    // CCTK_REAL PDstandardNth3gzz = INITVALUE;
    // CCTK_REAL PDstandardNth11gzz = INITVALUE;
    // CCTK_REAL PDstandardNth22gzz = INITVALUE;
    // CCTK_REAL PDstandardNth33gzz = INITVALUE;
    // CCTK_REAL PDstandardNth12gzz = INITVALUE;
    // CCTK_REAL PDstandardNth13gzz = INITVALUE;
    // CCTK_REAL PDstandardNth23gzz = INITVALUE;
    // CCTK_REAL PDstandardNth1kxx = INITVALUE;
    // CCTK_REAL PDstandardNth2kxx = INITVALUE;
    // CCTK_REAL PDstandardNth3kxx = INITVALUE;
    // CCTK_REAL PDstandardNth1kxy = INITVALUE;
    // CCTK_REAL PDstandardNth2kxy = INITVALUE;
    // CCTK_REAL PDstandardNth3kxy = INITVALUE;
    // CCTK_REAL PDstandardNth1kxz = INITVALUE;
    // CCTK_REAL PDstandardNth2kxz = INITVALUE;
    // CCTK_REAL PDstandardNth3kxz = INITVALUE;
    // CCTK_REAL PDstandardNth1kyy = INITVALUE;
    // CCTK_REAL PDstandardNth2kyy = INITVALUE;
    // CCTK_REAL PDstandardNth3kyy = INITVALUE;
    // CCTK_REAL PDstandardNth1kyz = INITVALUE;
    // CCTK_REAL PDstandardNth2kyz = INITVALUE;
    // CCTK_REAL PDstandardNth3kyz = INITVALUE;
    // CCTK_REAL PDstandardNth1kzz = INITVALUE;
    // CCTK_REAL PDstandardNth2kzz = INITVALUE;
    // CCTK_REAL PDstandardNth3kzz = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL const alpL = alp[index];
    CCTK_REAL const betaxL = betax[index];
    CCTK_REAL const betayL = betay[index];
    CCTK_REAL const betazL = betaz[index];
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
    CCTK_REAL const gxxL = gxx[index];
    CCTK_REAL const gxyL = gxy[index];
    CCTK_REAL const gxzL = gxz[index];
    CCTK_REAL const gyyL = gyy[index];
    CCTK_REAL const gyzL = gyz[index];
    CCTK_REAL const gzzL = gzz[index];
    CCTK_REAL const J11L = J11[index];
    CCTK_REAL const J12L = J12[index];
    CCTK_REAL const J13L = J13[index];
    CCTK_REAL const J21L = J21[index];
    CCTK_REAL const J22L = J22[index];
    CCTK_REAL const J23L = J23[index];
    CCTK_REAL const J31L = J31[index];
    CCTK_REAL const J32L = J32[index];
    CCTK_REAL const J33L = J33[index];
    CCTK_REAL const kxxL = kxx[index];
    CCTK_REAL const kxyL = kxy[index];
    CCTK_REAL const kxzL = kxz[index];
    CCTK_REAL const kyyL = kyy[index];
    CCTK_REAL const kyzL = kyz[index];
    CCTK_REAL const kzzL = kzz[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    CCTK_REAL const PDstandardNth1gxx = PDstandardNth1(gxx, i, j, k);
    CCTK_REAL const PDstandardNth2gxx = PDstandardNth2(gxx, i, j, k);
    CCTK_REAL const PDstandardNth3gxx = PDstandardNth3(gxx, i, j, k);
    CCTK_REAL const PDstandardNth11gxx = PDstandardNth11(gxx, i, j, k);
    CCTK_REAL const PDstandardNth22gxx = PDstandardNth22(gxx, i, j, k);
    CCTK_REAL const PDstandardNth33gxx = PDstandardNth33(gxx, i, j, k);
    CCTK_REAL const PDstandardNth12gxx = PDstandardNth12(gxx, i, j, k);
    CCTK_REAL const PDstandardNth13gxx = PDstandardNth13(gxx, i, j, k);
    CCTK_REAL const PDstandardNth23gxx = PDstandardNth23(gxx, i, j, k);
    CCTK_REAL const PDstandardNth1gxy = PDstandardNth1(gxy, i, j, k);
    CCTK_REAL const PDstandardNth2gxy = PDstandardNth2(gxy, i, j, k);
    CCTK_REAL const PDstandardNth3gxy = PDstandardNth3(gxy, i, j, k);
    CCTK_REAL const PDstandardNth11gxy = PDstandardNth11(gxy, i, j, k);
    CCTK_REAL const PDstandardNth22gxy = PDstandardNth22(gxy, i, j, k);
    CCTK_REAL const PDstandardNth33gxy = PDstandardNth33(gxy, i, j, k);
    CCTK_REAL const PDstandardNth12gxy = PDstandardNth12(gxy, i, j, k);
    CCTK_REAL const PDstandardNth13gxy = PDstandardNth13(gxy, i, j, k);
    CCTK_REAL const PDstandardNth23gxy = PDstandardNth23(gxy, i, j, k);
    CCTK_REAL const PDstandardNth1gxz = PDstandardNth1(gxz, i, j, k);
    CCTK_REAL const PDstandardNth2gxz = PDstandardNth2(gxz, i, j, k);
    CCTK_REAL const PDstandardNth3gxz = PDstandardNth3(gxz, i, j, k);
    CCTK_REAL const PDstandardNth11gxz = PDstandardNth11(gxz, i, j, k);
    CCTK_REAL const PDstandardNth22gxz = PDstandardNth22(gxz, i, j, k);
    CCTK_REAL const PDstandardNth33gxz = PDstandardNth33(gxz, i, j, k);
    CCTK_REAL const PDstandardNth12gxz = PDstandardNth12(gxz, i, j, k);
    CCTK_REAL const PDstandardNth13gxz = PDstandardNth13(gxz, i, j, k);
    CCTK_REAL const PDstandardNth23gxz = PDstandardNth23(gxz, i, j, k);
    CCTK_REAL const PDstandardNth1gyy = PDstandardNth1(gyy, i, j, k);
    CCTK_REAL const PDstandardNth2gyy = PDstandardNth2(gyy, i, j, k);
    CCTK_REAL const PDstandardNth3gyy = PDstandardNth3(gyy, i, j, k);
    CCTK_REAL const PDstandardNth11gyy = PDstandardNth11(gyy, i, j, k);
    CCTK_REAL const PDstandardNth22gyy = PDstandardNth22(gyy, i, j, k);
    CCTK_REAL const PDstandardNth33gyy = PDstandardNth33(gyy, i, j, k);
    CCTK_REAL const PDstandardNth12gyy = PDstandardNth12(gyy, i, j, k);
    CCTK_REAL const PDstandardNth13gyy = PDstandardNth13(gyy, i, j, k);
    CCTK_REAL const PDstandardNth23gyy = PDstandardNth23(gyy, i, j, k);
    CCTK_REAL const PDstandardNth1gyz = PDstandardNth1(gyz, i, j, k);
    CCTK_REAL const PDstandardNth2gyz = PDstandardNth2(gyz, i, j, k);
    CCTK_REAL const PDstandardNth3gyz = PDstandardNth3(gyz, i, j, k);
    CCTK_REAL const PDstandardNth11gyz = PDstandardNth11(gyz, i, j, k);
    CCTK_REAL const PDstandardNth22gyz = PDstandardNth22(gyz, i, j, k);
    CCTK_REAL const PDstandardNth33gyz = PDstandardNth33(gyz, i, j, k);
    CCTK_REAL const PDstandardNth12gyz = PDstandardNth12(gyz, i, j, k);
    CCTK_REAL const PDstandardNth13gyz = PDstandardNth13(gyz, i, j, k);
    CCTK_REAL const PDstandardNth23gyz = PDstandardNth23(gyz, i, j, k);
    CCTK_REAL const PDstandardNth1gzz = PDstandardNth1(gzz, i, j, k);
    CCTK_REAL const PDstandardNth2gzz = PDstandardNth2(gzz, i, j, k);
    CCTK_REAL const PDstandardNth3gzz = PDstandardNth3(gzz, i, j, k);
    CCTK_REAL const PDstandardNth11gzz = PDstandardNth11(gzz, i, j, k);
    CCTK_REAL const PDstandardNth22gzz = PDstandardNth22(gzz, i, j, k);
    CCTK_REAL const PDstandardNth33gzz = PDstandardNth33(gzz, i, j, k);
    CCTK_REAL const PDstandardNth12gzz = PDstandardNth12(gzz, i, j, k);
    CCTK_REAL const PDstandardNth13gzz = PDstandardNth13(gzz, i, j, k);
    CCTK_REAL const PDstandardNth23gzz = PDstandardNth23(gzz, i, j, k);
    CCTK_REAL const PDstandardNth1kxx = PDstandardNth1(kxx, i, j, k);
    CCTK_REAL const PDstandardNth2kxx = PDstandardNth2(kxx, i, j, k);
    CCTK_REAL const PDstandardNth3kxx = PDstandardNth3(kxx, i, j, k);
    CCTK_REAL const PDstandardNth1kxy = PDstandardNth1(kxy, i, j, k);
    CCTK_REAL const PDstandardNth2kxy = PDstandardNth2(kxy, i, j, k);
    CCTK_REAL const PDstandardNth3kxy = PDstandardNth3(kxy, i, j, k);
    CCTK_REAL const PDstandardNth1kxz = PDstandardNth1(kxz, i, j, k);
    CCTK_REAL const PDstandardNth2kxz = PDstandardNth2(kxz, i, j, k);
    CCTK_REAL const PDstandardNth3kxz = PDstandardNth3(kxz, i, j, k);
    CCTK_REAL const PDstandardNth1kyy = PDstandardNth1(kyy, i, j, k);
    CCTK_REAL const PDstandardNth2kyy = PDstandardNth2(kyy, i, j, k);
    CCTK_REAL const PDstandardNth3kyy = PDstandardNth3(kyy, i, j, k);
    CCTK_REAL const PDstandardNth1kyz = PDstandardNth1(kyz, i, j, k);
    CCTK_REAL const PDstandardNth2kyz = PDstandardNth2(kyz, i, j, k);
    CCTK_REAL const PDstandardNth3kyz = PDstandardNth3(kyz, i, j, k);
    CCTK_REAL const PDstandardNth1kzz = PDstandardNth1(kzz, i, j, k);
    CCTK_REAL const PDstandardNth2kzz = PDstandardNth2(kzz, i, j, k);
    CCTK_REAL const PDstandardNth3kzz = PDstandardNth3(kzz, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL const detg  =  2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - SQR(gxyL)) - gyyL*SQR(gxzL) - gxxL*SQR(gyzL);
    
    CCTK_REAL const gu11  =  INV(detg)*(gyyL*gzzL - SQR(gyzL));
    
    CCTK_REAL const gu21  =  (gxzL*gyzL - gxyL*gzzL)*INV(detg);
    
    CCTK_REAL const gu31  =  (-(gxzL*gyyL) + gxyL*gyzL)*INV(detg);
    
    CCTK_REAL const gu22  =  INV(detg)*(gxxL*gzzL - SQR(gxzL));
    
    CCTK_REAL const gu32  =  (gxyL*gxzL - gxxL*gyzL)*INV(detg);
    
    CCTK_REAL const gu33  =  INV(detg)*(gxxL*gyyL - SQR(gxyL));
    
    CCTK_REAL const G111  =  khalf*((gu11*J11L - gu21*J12L - gu31*J13L)*PDstandardNth1gxx + 
          (gu11*J21L - gu21*J22L - gu31*J23L)*PDstandardNth2gxx + (gu11*J31L - gu21*J32L - gu31*J33L)*PDstandardNth3gxx + 
          2*(J11L*(gu21*PDstandardNth1gxy + gu31*PDstandardNth1gxz) + 
             J21L*(gu21*PDstandardNth2gxy + gu31*PDstandardNth2gxz) + J31L*(gu21*PDstandardNth3gxy + gu31*PDstandardNth3gxz)
             ));
    
    CCTK_REAL const G211  =  khalf*((gu21*J11L - gu22*J12L - gu32*J13L)*PDstandardNth1gxx + 
          (gu21*J21L - gu22*J22L - gu32*J23L)*PDstandardNth2gxx + (gu21*J31L - gu22*J32L - gu32*J33L)*PDstandardNth3gxx + 
          2*(J11L*(gu22*PDstandardNth1gxy + gu32*PDstandardNth1gxz) + 
             J21L*(gu22*PDstandardNth2gxy + gu32*PDstandardNth2gxz) + J31L*(gu22*PDstandardNth3gxy + gu32*PDstandardNth3gxz)
             ));
    
    CCTK_REAL const G311  =  khalf*((gu31*J11L - gu32*J12L - gu33*J13L)*PDstandardNth1gxx + 
          (gu31*J21L - gu32*J22L - gu33*J23L)*PDstandardNth2gxx + (gu31*J31L - gu32*J32L - gu33*J33L)*PDstandardNth3gxx + 
          2*(J11L*(gu32*PDstandardNth1gxy + gu33*PDstandardNth1gxz) + 
             J21L*(gu32*PDstandardNth2gxy + gu33*PDstandardNth2gxz) + J31L*(gu32*PDstandardNth3gxy + gu33*PDstandardNth3gxz)
             ));
    
    CCTK_REAL const G112  =  khalf*(gu11*(J12L*PDstandardNth1gxx + J22L*PDstandardNth2gxx + J32L*PDstandardNth3gxx) + 
          gu21*(J11L*PDstandardNth1gyy + J21L*PDstandardNth2gyy + J31L*PDstandardNth3gyy) + 
          gu31*(-(J13L*PDstandardNth1gxy) + J12L*PDstandardNth1gxz + J11L*PDstandardNth1gyz - J23L*PDstandardNth2gxy + 
             J22L*PDstandardNth2gxz + J21L*PDstandardNth2gyz - J33L*PDstandardNth3gxy + J32L*PDstandardNth3gxz + 
             J31L*PDstandardNth3gyz));
    
    CCTK_REAL const G212  =  khalf*(gu21*(J12L*PDstandardNth1gxx + J22L*PDstandardNth2gxx + J32L*PDstandardNth3gxx) + 
          gu22*(J11L*PDstandardNth1gyy + J21L*PDstandardNth2gyy + J31L*PDstandardNth3gyy) + 
          gu32*(-(J13L*PDstandardNth1gxy) + J12L*PDstandardNth1gxz + J11L*PDstandardNth1gyz - J23L*PDstandardNth2gxy + 
             J22L*PDstandardNth2gxz + J21L*PDstandardNth2gyz - J33L*PDstandardNth3gxy + J32L*PDstandardNth3gxz + 
             J31L*PDstandardNth3gyz));
    
    CCTK_REAL const G312  =  khalf*(gu31*(J12L*PDstandardNth1gxx + J22L*PDstandardNth2gxx + J32L*PDstandardNth3gxx) + 
          gu32*(J11L*PDstandardNth1gyy + J21L*PDstandardNth2gyy + J31L*PDstandardNth3gyy) + 
          gu33*(-(J13L*PDstandardNth1gxy) + J12L*PDstandardNth1gxz + J11L*PDstandardNth1gyz - J23L*PDstandardNth2gxy + 
             J22L*PDstandardNth2gxz + J21L*PDstandardNth2gyz - J33L*PDstandardNth3gxy + J32L*PDstandardNth3gxz + 
             J31L*PDstandardNth3gyz));
    
    CCTK_REAL const G113  =  khalf*(gu11*(J13L*PDstandardNth1gxx + J23L*PDstandardNth2gxx + J33L*PDstandardNth3gxx) + 
          gu21*(J13L*PDstandardNth1gxy - J12L*PDstandardNth1gxz + J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy - 
             J22L*PDstandardNth2gxz + J21L*PDstandardNth2gyz + J33L*PDstandardNth3gxy - J32L*PDstandardNth3gxz + 
             J31L*PDstandardNth3gyz) + gu31*(J11L*PDstandardNth1gzz + J21L*PDstandardNth2gzz + J31L*PDstandardNth3gzz));
    
    CCTK_REAL const G213  =  khalf*(gu21*(J13L*PDstandardNth1gxx + J23L*PDstandardNth2gxx + J33L*PDstandardNth3gxx) + 
          gu22*(J13L*PDstandardNth1gxy - J12L*PDstandardNth1gxz + J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy - 
             J22L*PDstandardNth2gxz + J21L*PDstandardNth2gyz + J33L*PDstandardNth3gxy - J32L*PDstandardNth3gxz + 
             J31L*PDstandardNth3gyz) + gu32*(J11L*PDstandardNth1gzz + J21L*PDstandardNth2gzz + J31L*PDstandardNth3gzz));
    
    CCTK_REAL const G313  =  khalf*(gu31*(J13L*PDstandardNth1gxx + J23L*PDstandardNth2gxx + J33L*PDstandardNth3gxx) + 
          gu32*(J13L*PDstandardNth1gxy - J12L*PDstandardNth1gxz + J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy - 
             J22L*PDstandardNth2gxz + J21L*PDstandardNth2gyz + J33L*PDstandardNth3gxy - J32L*PDstandardNth3gxz + 
             J31L*PDstandardNth3gyz) + gu33*(J11L*PDstandardNth1gzz + J21L*PDstandardNth2gzz + J31L*PDstandardNth3gzz));
    
    CCTK_REAL const G122  =  khalf*(gu11*(-(J11L*PDstandardNth1gyy) + 2*(J12L*PDstandardNth1gxy + J22L*PDstandardNth2gxy) - 
             J21L*PDstandardNth2gyy + 2*J32L*PDstandardNth3gxy - J31L*PDstandardNth3gyy) + 
          gu21*(J12L*PDstandardNth1gyy + J22L*PDstandardNth2gyy + J32L*PDstandardNth3gyy) - 
          gu31*(J13L*PDstandardNth1gyy + J23L*PDstandardNth2gyy + J33L*PDstandardNth3gyy - 
             2*(J12L*PDstandardNth1gyz + J22L*PDstandardNth2gyz + J32L*PDstandardNth3gyz)));
    
    CCTK_REAL const G222  =  khalf*(gu21*(-(J11L*PDstandardNth1gyy) + 2*(J12L*PDstandardNth1gxy + J22L*PDstandardNth2gxy) - 
             J21L*PDstandardNth2gyy + 2*J32L*PDstandardNth3gxy - J31L*PDstandardNth3gyy) + 
          gu22*(J12L*PDstandardNth1gyy + J22L*PDstandardNth2gyy + J32L*PDstandardNth3gyy) - 
          gu32*(J13L*PDstandardNth1gyy + J23L*PDstandardNth2gyy + J33L*PDstandardNth3gyy - 
             2*(J12L*PDstandardNth1gyz + J22L*PDstandardNth2gyz + J32L*PDstandardNth3gyz)));
    
    CCTK_REAL const G322  =  khalf*(gu31*(-(J11L*PDstandardNth1gyy) + 2*(J12L*PDstandardNth1gxy + J22L*PDstandardNth2gxy) - 
             J21L*PDstandardNth2gyy + 2*J32L*PDstandardNth3gxy - J31L*PDstandardNth3gyy) + 
          gu32*(J12L*PDstandardNth1gyy + J22L*PDstandardNth2gyy + J32L*PDstandardNth3gyy) - 
          gu33*(J13L*PDstandardNth1gyy + J23L*PDstandardNth2gyy + J33L*PDstandardNth3gyy - 
             2*(J12L*PDstandardNth1gyz + J22L*PDstandardNth2gyz + J32L*PDstandardNth3gyz)));
    
    CCTK_REAL const G123  =  khalf*(gu21*(J13L*PDstandardNth1gyy + J23L*PDstandardNth2gyy + J33L*PDstandardNth3gyy) + 
          gu11*(J13L*PDstandardNth1gxy + J12L*PDstandardNth1gxz - J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy + 
             J22L*PDstandardNth2gxz - J21L*PDstandardNth2gyz + J33L*PDstandardNth3gxy + J32L*PDstandardNth3gxz - 
             J31L*PDstandardNth3gyz) + gu31*(J12L*PDstandardNth1gzz + J22L*PDstandardNth2gzz + J32L*PDstandardNth3gzz));
    
    CCTK_REAL const G223  =  khalf*(gu22*(J13L*PDstandardNth1gyy + J23L*PDstandardNth2gyy + J33L*PDstandardNth3gyy) + 
          gu21*(J13L*PDstandardNth1gxy + J12L*PDstandardNth1gxz - J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy + 
             J22L*PDstandardNth2gxz - J21L*PDstandardNth2gyz + J33L*PDstandardNth3gxy + J32L*PDstandardNth3gxz - 
             J31L*PDstandardNth3gyz) + gu32*(J12L*PDstandardNth1gzz + J22L*PDstandardNth2gzz + J32L*PDstandardNth3gzz));
    
    CCTK_REAL const G323  =  khalf*(gu32*(J13L*PDstandardNth1gyy + J23L*PDstandardNth2gyy + J33L*PDstandardNth3gyy) + 
          gu31*(J13L*PDstandardNth1gxy + J12L*PDstandardNth1gxz - J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy + 
             J22L*PDstandardNth2gxz - J21L*PDstandardNth2gyz + J33L*PDstandardNth3gxy + J32L*PDstandardNth3gxz - 
             J31L*PDstandardNth3gyz) + gu33*(J12L*PDstandardNth1gzz + J22L*PDstandardNth2gzz + J32L*PDstandardNth3gzz));
    
    CCTK_REAL const G133  =  khalf*(gu11*(-(J11L*PDstandardNth1gzz) + 2*(J13L*PDstandardNth1gxz + J23L*PDstandardNth2gxz) - 
             J21L*PDstandardNth2gzz + 2*J33L*PDstandardNth3gxz - J31L*PDstandardNth3gzz) + 
          gu21*(-(J12L*PDstandardNth1gzz) + 2*(J13L*PDstandardNth1gyz + J23L*PDstandardNth2gyz) - J22L*PDstandardNth2gzz + 
             2*J33L*PDstandardNth3gyz - J32L*PDstandardNth3gzz) + 
          gu31*(J13L*PDstandardNth1gzz + J23L*PDstandardNth2gzz + J33L*PDstandardNth3gzz));
    
    CCTK_REAL const G233  =  khalf*(gu21*(-(J11L*PDstandardNth1gzz) + 2*(J13L*PDstandardNth1gxz + J23L*PDstandardNth2gxz) - 
             J21L*PDstandardNth2gzz + 2*J33L*PDstandardNth3gxz - J31L*PDstandardNth3gzz) + 
          gu22*(-(J12L*PDstandardNth1gzz) + 2*(J13L*PDstandardNth1gyz + J23L*PDstandardNth2gyz) - J22L*PDstandardNth2gzz + 
             2*J33L*PDstandardNth3gyz - J32L*PDstandardNth3gzz) + 
          gu32*(J13L*PDstandardNth1gzz + J23L*PDstandardNth2gzz + J33L*PDstandardNth3gzz));
    
    CCTK_REAL const G333  =  khalf*(gu31*(-(J11L*PDstandardNth1gzz) + 2*(J13L*PDstandardNth1gxz + J23L*PDstandardNth2gxz) - 
             J21L*PDstandardNth2gzz + 2*J33L*PDstandardNth3gxz - J31L*PDstandardNth3gzz) + 
          gu32*(-(J12L*PDstandardNth1gzz) + 2*(J13L*PDstandardNth1gyz + J23L*PDstandardNth2gyz) - J22L*PDstandardNth2gzz + 
             2*J33L*PDstandardNth3gyz - J32L*PDstandardNth3gzz) + 
          gu33*(J13L*PDstandardNth1gzz + J23L*PDstandardNth2gzz + J33L*PDstandardNth3gzz));
    
    CCTK_REAL const R11  =  khalf*(-((dJ122L*gu22 + 2*dJ123L*gu32 + dJ133L*gu33)*PDstandardNth1gxx) + 
          gu32*(-4*J11L*J21L*PDstandardNth12gyz - 2*J13L*J32L*PDstandardNth13gxx + 2*J13L*J31L*PDstandardNth13gxy + 
             2*J12L*J31L*PDstandardNth13gxz - 4*J11L*J31L*PDstandardNth13gyz + 2*dJ113L*PDstandardNth1gxy) + 
          2*(gu22*J11L*J12L*PDstandardNth11gxy + gu32*J11L*J13L*PDstandardNth11gxy + gu32*J11L*J12L*PDstandardNth11gxz + 
             gu33*J11L*J13L*PDstandardNth11gxz + gu22*J12L*J21L*PDstandardNth12gxy + gu32*J13L*J21L*PDstandardNth12gxy + 
             gu22*J11L*J22L*PDstandardNth12gxy + gu32*J11L*J23L*PDstandardNth12gxy + gu32*J12L*J21L*PDstandardNth12gxz + 
             gu33*J13L*J21L*PDstandardNth12gxz + gu32*J11L*J22L*PDstandardNth12gxz + gu33*J11L*J23L*PDstandardNth12gxz + 
             gu22*J12L*J31L*PDstandardNth13gxy + gu22*J11L*J32L*PDstandardNth13gxy + gu32*J11L*J33L*PDstandardNth13gxy + 
             gu33*J13L*J31L*PDstandardNth13gxz + gu32*J11L*J32L*PDstandardNth13gxz + gu33*J11L*J33L*PDstandardNth13gxz + 
             dJ112L*gu22*PDstandardNth1gxy) + 2*dJ112L*gu32*PDstandardNth1gxz + 2*dJ113L*gu33*PDstandardNth1gxz - 
          dJ111L*gu22*PDstandardNth1gyy - 2*(G111*G122 + G111*G133 + G211*G222 + G211*G233 + G311*G322 + G311*G333 + 
             gu32*J12L*J13L*PDstandardNth11gxx + gu22*J12L*J22L*PDstandardNth12gxx + gu32*J13L*J22L*PDstandardNth12gxx + 
             gu32*J12L*J23L*PDstandardNth12gxx + gu33*J13L*J23L*PDstandardNth12gxx + gu22*J11L*J21L*PDstandardNth12gyy + 
             gu33*J11L*J21L*PDstandardNth12gzz + gu22*J12L*J32L*PDstandardNth13gxx + gu32*J12L*J33L*PDstandardNth13gxx + 
             gu33*J13L*J33L*PDstandardNth13gxx + gu22*J11L*J31L*PDstandardNth13gyy + gu33*J11L*J31L*PDstandardNth13gzz + 
             dJ111L*gu32*PDstandardNth1gyz) - dJ111L*gu33*PDstandardNth1gzz - 2*gu32*J22L*J23L*PDstandardNth22gxx + 
          2*gu22*J21L*J22L*PDstandardNth22gxy + 2*gu32*J21L*J23L*PDstandardNth22gxy + 2*gu32*J21L*J22L*PDstandardNth22gxz + 
          2*gu33*J21L*J23L*PDstandardNth22gxz - 2*gu22*J22L*J32L*PDstandardNth23gxx - 2*gu32*J23L*J32L*PDstandardNth23gxx - 
          2*gu32*J22L*J33L*PDstandardNth23gxx - 2*gu33*J23L*J33L*PDstandardNth23gxx + 2*gu22*J22L*J31L*PDstandardNth23gxy + 
          2*gu32*J23L*J31L*PDstandardNth23gxy + 2*gu22*J21L*J32L*PDstandardNth23gxy + 2*gu32*J21L*J33L*PDstandardNth23gxy + 
          2*gu32*J22L*J31L*PDstandardNth23gxz + 2*gu33*J23L*J31L*PDstandardNth23gxz + 2*gu32*J21L*J32L*PDstandardNth23gxz + 
          2*gu33*J21L*J33L*PDstandardNth23gxz - 2*gu22*J21L*J31L*PDstandardNth23gyy - 4*gu32*J21L*J31L*PDstandardNth23gyz - 
          2*gu33*J21L*J31L*PDstandardNth23gzz - (dJ222L*gu22 + 2*dJ223L*gu32 + dJ233L*gu33)*PDstandardNth2gxx + 
          2*dJ212L*gu22*PDstandardNth2gxy + 2*dJ213L*gu32*PDstandardNth2gxy + 2*dJ212L*gu32*PDstandardNth2gxz + 
          2*dJ213L*gu33*PDstandardNth2gxz - dJ211L*gu22*PDstandardNth2gyy - 2*dJ211L*gu32*PDstandardNth2gyz - 
          dJ211L*gu33*PDstandardNth2gzz - 2*gu32*J32L*J33L*PDstandardNth33gxx + 2*gu22*J31L*J32L*PDstandardNth33gxy + 
          2*gu32*J31L*J33L*PDstandardNth33gxy + 2*gu32*J31L*J32L*PDstandardNth33gxz + 2*gu33*J31L*J33L*PDstandardNth33gxz - 
          dJ322L*gu22*PDstandardNth3gxx - 2*dJ323L*gu32*PDstandardNth3gxx - dJ333L*gu33*PDstandardNth3gxx + 
          2*dJ312L*gu22*PDstandardNth3gxy + 2*dJ313L*gu32*PDstandardNth3gxy + 2*dJ312L*gu32*PDstandardNth3gxz + 
          2*dJ313L*gu33*PDstandardNth3gxz - dJ311L*gu22*PDstandardNth3gyy - 2*dJ311L*gu32*PDstandardNth3gyz - 
          dJ311L*gu33*PDstandardNth3gzz + 2*SQR(G112) + 2*SQR(G113) + 2*SQR(G212) + 2*SQR(G213) + 2*SQR(G312) + 
          2*SQR(G313) - gu22*PDstandardNth11gyy*SQR(J11L) - 2*gu32*PDstandardNth11gyz*SQR(J11L) - 
          gu33*PDstandardNth11gzz*SQR(J11L) - gu22*PDstandardNth11gxx*SQR(J12L) - gu33*PDstandardNth11gxx*SQR(J13L) - 
          gu22*PDstandardNth22gyy*SQR(J21L) - 2*gu32*PDstandardNth22gyz*SQR(J21L) - gu33*PDstandardNth22gzz*SQR(J21L) - 
          gu22*PDstandardNth22gxx*SQR(J22L) - gu33*PDstandardNth22gxx*SQR(J23L) - gu22*PDstandardNth33gyy*SQR(J31L) - 
          2*gu32*PDstandardNth33gyz*SQR(J31L) - gu33*PDstandardNth33gzz*SQR(J31L) - gu22*PDstandardNth33gxx*SQR(J32L) - 
          gu33*PDstandardNth33gxx*SQR(J33L));
    
    CCTK_REAL const R12  =  G113*G123 - G112*G133 + G211*G212 + G212*G222 + G213*G223 - G212*(G211 + G222 + G233) + G311*G312 + G312*G322 + 
        G313*G323 - G312*(G311 + G322 + G333) + khalf*
         (gu32*(-(J13L*J22L*PDstandardNth12gxy) + J12L*J23L*PDstandardNth12gxy + 
              (J13L*J22L - J12L*J23L)*PDstandardNth12gxy - J13L*J32L*PDstandardNth13gxy + J12L*J33L*PDstandardNth13gxy + 
              (J13L*J32L - J12L*J33L)*PDstandardNth13gxy) + 
           gu22*(-(J12L*J21L*PDstandardNth12gyy) + J11L*J22L*PDstandardNth12gyy + 
              (J12L*J21L - J11L*J22L)*PDstandardNth12gyy - J12L*J31L*PDstandardNth13gyy + J11L*J32L*PDstandardNth13gyy + 
              (J12L*J31L - J11L*J32L)*PDstandardNth13gyy) + 
           gu31*(J12L*J13L*PDstandardNth11gxx - J11L*J13L*PDstandardNth11gxy - J11L*J12L*PDstandardNth11gxz + 
              J13L*J22L*PDstandardNth12gxx + J12L*J23L*PDstandardNth12gxx - J13L*J21L*PDstandardNth12gxy - 
              J11L*J23L*PDstandardNth12gxy - J12L*J21L*PDstandardNth12gxz - J11L*J22L*PDstandardNth12gxz + 
              2*J11L*J21L*PDstandardNth12gyz + J13L*J32L*PDstandardNth13gxx + J12L*J33L*PDstandardNth13gxx - 
              J13L*J31L*PDstandardNth13gxy - J11L*J33L*PDstandardNth13gxy - J12L*J31L*PDstandardNth13gxz - 
              J11L*J32L*PDstandardNth13gxz + 2*J11L*J31L*PDstandardNth13gyz + dJ123L*PDstandardNth1gxx - 
              dJ113L*PDstandardNth1gxy - dJ112L*PDstandardNth1gxz + dJ111L*PDstandardNth1gyz + 
              J22L*J23L*PDstandardNth22gxx - J21L*J23L*PDstandardNth22gxy - J21L*J22L*PDstandardNth22gxz + 
              J23L*J32L*PDstandardNth23gxx + J22L*J33L*PDstandardNth23gxx - J23L*J31L*PDstandardNth23gxy - 
              J21L*J33L*PDstandardNth23gxy - J22L*J31L*PDstandardNth23gxz - J21L*J32L*PDstandardNth23gxz + 
              2*J21L*J31L*PDstandardNth23gyz + dJ223L*PDstandardNth2gxx - dJ213L*PDstandardNth2gxy - 
              dJ212L*PDstandardNth2gxz + dJ211L*PDstandardNth2gyz + J32L*J33L*PDstandardNth33gxx - 
              J31L*J33L*PDstandardNth33gxy - J31L*J32L*PDstandardNth33gxz + dJ323L*PDstandardNth3gxx - 
              dJ313L*PDstandardNth3gxy - dJ312L*PDstandardNth3gxz + dJ311L*PDstandardNth3gyz + 
              PDstandardNth11gyz*SQR(J11L) + PDstandardNth22gyz*SQR(J21L) + PDstandardNth33gyz*SQR(J31L)) + 
           gu21*(-2*J11L*J12L*PDstandardNth11gxy + 2*J12L*J22L*PDstandardNth12gxx - 2*J12L*J21L*PDstandardNth12gxy - 
              2*J11L*J22L*PDstandardNth12gxy + 2*J11L*J21L*PDstandardNth12gyy + 2*J12L*J32L*PDstandardNth13gxx - 
              2*J12L*J31L*PDstandardNth13gxy - 2*J11L*J32L*PDstandardNth13gxy + 2*J11L*J31L*PDstandardNth13gyy + 
              dJ122L*PDstandardNth1gxx - 2*dJ112L*PDstandardNth1gxy + dJ111L*PDstandardNth1gyy - 
              2*J21L*J22L*PDstandardNth22gxy + 2*J22L*J32L*PDstandardNth23gxx - 2*J22L*J31L*PDstandardNth23gxy - 
              2*J21L*J32L*PDstandardNth23gxy + 2*J21L*J31L*PDstandardNth23gyy + dJ222L*PDstandardNth2gxx - 
              2*dJ212L*PDstandardNth2gxy + dJ211L*PDstandardNth2gyy - 2*J31L*J32L*PDstandardNth33gxy + 
              dJ322L*PDstandardNth3gxx - 2*dJ312L*PDstandardNth3gxy + dJ311L*PDstandardNth3gyy + 
              PDstandardNth11gyy*SQR(J11L) + PDstandardNth11gxx*SQR(J12L) + PDstandardNth22gyy*SQR(J21L) + 
              PDstandardNth22gxx*SQR(J22L) + PDstandardNth33gyy*SQR(J31L) + PDstandardNth33gxx*SQR(J32L)) + 
           gu32*(-(J12L*J13L*PDstandardNth11gxy) + J11L*J13L*PDstandardNth11gyy - J11L*J12L*PDstandardNth11gyz - 
              J13L*J22L*PDstandardNth12gxy - J12L*J23L*PDstandardNth12gxy + 2*J12L*J22L*PDstandardNth12gxz + 
              J13L*J21L*PDstandardNth12gyy + J11L*J23L*PDstandardNth12gyy - J12L*J21L*PDstandardNth12gyz - 
              J11L*J22L*PDstandardNth12gyz - J13L*J32L*PDstandardNth13gxy - J12L*J33L*PDstandardNth13gxy + 
              2*J12L*J32L*PDstandardNth13gxz + J13L*J31L*PDstandardNth13gyy + J11L*J33L*PDstandardNth13gyy - 
              J12L*J31L*PDstandardNth13gyz - J11L*J32L*PDstandardNth13gyz - dJ123L*PDstandardNth1gxy + 
              dJ122L*PDstandardNth1gxz + dJ113L*PDstandardNth1gyy - dJ112L*PDstandardNth1gyz - 
              J22L*J23L*PDstandardNth22gxy + J21L*J23L*PDstandardNth22gyy - J21L*J22L*PDstandardNth22gyz - 
              J23L*J32L*PDstandardNth23gxy - J22L*J33L*PDstandardNth23gxy + 2*J22L*J32L*PDstandardNth23gxz + 
              J23L*J31L*PDstandardNth23gyy + J21L*J33L*PDstandardNth23gyy - J22L*J31L*PDstandardNth23gyz - 
              J21L*J32L*PDstandardNth23gyz - dJ223L*PDstandardNth2gxy + dJ222L*PDstandardNth2gxz + 
              dJ213L*PDstandardNth2gyy - dJ212L*PDstandardNth2gyz - J32L*J33L*PDstandardNth33gxy + 
              J31L*J33L*PDstandardNth33gyy - J31L*J32L*PDstandardNth33gyz - dJ323L*PDstandardNth3gxy + 
              dJ322L*PDstandardNth3gxz + dJ313L*PDstandardNth3gyy - dJ312L*PDstandardNth3gyz + 
              PDstandardNth11gxz*SQR(J12L) + PDstandardNth22gxz*SQR(J22L) + PDstandardNth33gxz*SQR(J32L)) + 
           gu33*(J12L*J13L*PDstandardNth11gxz + J11L*J13L*PDstandardNth11gyz - J11L*J12L*PDstandardNth11gzz - 
              2*J13L*J23L*PDstandardNth12gxy + J13L*J22L*PDstandardNth12gxz + J12L*J23L*PDstandardNth12gxz + 
              J13L*J21L*PDstandardNth12gyz + J11L*J23L*PDstandardNth12gyz - J12L*J21L*PDstandardNth12gzz - 
              J11L*J22L*PDstandardNth12gzz - 2*J13L*J33L*PDstandardNth13gxy + J13L*J32L*PDstandardNth13gxz + 
              J12L*J33L*PDstandardNth13gxz + J13L*J31L*PDstandardNth13gyz + J11L*J33L*PDstandardNth13gyz - 
              J12L*J31L*PDstandardNth13gzz - J11L*J32L*PDstandardNth13gzz - dJ133L*PDstandardNth1gxy + 
              dJ123L*PDstandardNth1gxz + dJ113L*PDstandardNth1gyz - dJ112L*PDstandardNth1gzz + 
              J22L*J23L*PDstandardNth22gxz + J21L*J23L*PDstandardNth22gyz - J21L*J22L*PDstandardNth22gzz - 
              2*J23L*J33L*PDstandardNth23gxy + J23L*J32L*PDstandardNth23gxz + J22L*J33L*PDstandardNth23gxz + 
              J23L*J31L*PDstandardNth23gyz + J21L*J33L*PDstandardNth23gyz - J22L*J31L*PDstandardNth23gzz - 
              J21L*J32L*PDstandardNth23gzz - dJ233L*PDstandardNth2gxy + dJ223L*PDstandardNth2gxz + 
              dJ213L*PDstandardNth2gyz - dJ212L*PDstandardNth2gzz + J32L*J33L*PDstandardNth33gxz + 
              J31L*J33L*PDstandardNth33gyz - J31L*J32L*PDstandardNth33gzz - dJ333L*PDstandardNth3gxy + 
              dJ323L*PDstandardNth3gxz + dJ313L*PDstandardNth3gyz - dJ312L*PDstandardNth3gzz - 
              PDstandardNth11gxy*SQR(J13L) - PDstandardNth22gxy*SQR(J23L) - PDstandardNth33gxy*SQR(J33L)));
    
    CCTK_REAL const R13  =  -(G113*G122) + G112*G123 + G211*G213 + G212*G223 + G213*G233 - G213*(G211 + G222 + G233) + G311*G313 + 
        G312*G323 + G313*G333 - G313*(G311 + G322 + G333) + 
        khalf*(gu32*(J13L*J22L*PDstandardNth12gxz - J12L*J23L*PDstandardNth12gxz + 
              (-(J13L*J22L) + J12L*J23L)*PDstandardNth12gxz + J13L*J32L*PDstandardNth13gxz - J12L*J33L*PDstandardNth13gxz + 
              (-(J13L*J32L) + J12L*J33L)*PDstandardNth13gxz) + 
           gu33*(-(J13L*J21L*PDstandardNth12gzz) + J11L*J23L*PDstandardNth12gzz + 
              (J13L*J21L - J11L*J23L)*PDstandardNth12gzz - J13L*J31L*PDstandardNth13gzz + J11L*J33L*PDstandardNth13gzz + 
              (J13L*J31L - J11L*J33L)*PDstandardNth13gzz) + 
           gu21*(J12L*J13L*PDstandardNth11gxx - J11L*J13L*PDstandardNth11gxy - J11L*J12L*PDstandardNth11gxz + 
              J13L*J22L*PDstandardNth12gxx + J12L*J23L*PDstandardNth12gxx - J13L*J21L*PDstandardNth12gxy - 
              J11L*J23L*PDstandardNth12gxy - J12L*J21L*PDstandardNth12gxz - J11L*J22L*PDstandardNth12gxz + 
              2*J11L*J21L*PDstandardNth12gyz + J13L*J32L*PDstandardNth13gxx + J12L*J33L*PDstandardNth13gxx - 
              J13L*J31L*PDstandardNth13gxy - J11L*J33L*PDstandardNth13gxy - J12L*J31L*PDstandardNth13gxz - 
              J11L*J32L*PDstandardNth13gxz + 2*J11L*J31L*PDstandardNth13gyz + dJ123L*PDstandardNth1gxx - 
              dJ113L*PDstandardNth1gxy - dJ112L*PDstandardNth1gxz + dJ111L*PDstandardNth1gyz + 
              J22L*J23L*PDstandardNth22gxx - J21L*J23L*PDstandardNth22gxy - J21L*J22L*PDstandardNth22gxz + 
              J23L*J32L*PDstandardNth23gxx + J22L*J33L*PDstandardNth23gxx - J23L*J31L*PDstandardNth23gxy - 
              J21L*J33L*PDstandardNth23gxy - J22L*J31L*PDstandardNth23gxz - J21L*J32L*PDstandardNth23gxz + 
              2*J21L*J31L*PDstandardNth23gyz + dJ223L*PDstandardNth2gxx - dJ213L*PDstandardNth2gxy - 
              dJ212L*PDstandardNth2gxz + dJ211L*PDstandardNth2gyz + J32L*J33L*PDstandardNth33gxx - 
              J31L*J33L*PDstandardNth33gxy - J31L*J32L*PDstandardNth33gxz + dJ323L*PDstandardNth3gxx - 
              dJ313L*PDstandardNth3gxy - dJ312L*PDstandardNth3gxz + dJ311L*PDstandardNth3gyz + 
              PDstandardNth11gyz*SQR(J11L) + PDstandardNth22gyz*SQR(J21L) + PDstandardNth33gyz*SQR(J31L)) + 
           gu22*(J12L*J13L*PDstandardNth11gxy - J11L*J13L*PDstandardNth11gyy + J11L*J12L*PDstandardNth11gyz + 
              J13L*J22L*PDstandardNth12gxy + J12L*J23L*PDstandardNth12gxy - 2*J12L*J22L*PDstandardNth12gxz - 
              J13L*J21L*PDstandardNth12gyy - J11L*J23L*PDstandardNth12gyy + J12L*J21L*PDstandardNth12gyz + 
              J11L*J22L*PDstandardNth12gyz + J13L*J32L*PDstandardNth13gxy + J12L*J33L*PDstandardNth13gxy - 
              2*J12L*J32L*PDstandardNth13gxz - J13L*J31L*PDstandardNth13gyy - J11L*J33L*PDstandardNth13gyy + 
              J12L*J31L*PDstandardNth13gyz + J11L*J32L*PDstandardNth13gyz + dJ123L*PDstandardNth1gxy - 
              dJ122L*PDstandardNth1gxz - dJ113L*PDstandardNth1gyy + dJ112L*PDstandardNth1gyz + 
              J22L*J23L*PDstandardNth22gxy - J21L*J23L*PDstandardNth22gyy + J21L*J22L*PDstandardNth22gyz + 
              J23L*J32L*PDstandardNth23gxy + J22L*J33L*PDstandardNth23gxy - 2*J22L*J32L*PDstandardNth23gxz - 
              J23L*J31L*PDstandardNth23gyy - J21L*J33L*PDstandardNth23gyy + J22L*J31L*PDstandardNth23gyz + 
              J21L*J32L*PDstandardNth23gyz + dJ223L*PDstandardNth2gxy - dJ222L*PDstandardNth2gxz - 
              dJ213L*PDstandardNth2gyy + dJ212L*PDstandardNth2gyz + J32L*J33L*PDstandardNth33gxy - 
              J31L*J33L*PDstandardNth33gyy + J31L*J32L*PDstandardNth33gyz + dJ323L*PDstandardNth3gxy - 
              dJ322L*PDstandardNth3gxz - dJ313L*PDstandardNth3gyy + dJ312L*PDstandardNth3gyz - 
              PDstandardNth11gxz*SQR(J12L) - PDstandardNth22gxz*SQR(J22L) - PDstandardNth33gxz*SQR(J32L)) + 
           gu31*(-2*J11L*J13L*PDstandardNth11gxz + 2*J13L*J23L*PDstandardNth12gxx - 2*J13L*J21L*PDstandardNth12gxz - 
              2*J11L*J23L*PDstandardNth12gxz + 2*J11L*J21L*PDstandardNth12gzz + 2*J13L*J33L*PDstandardNth13gxx - 
              2*J13L*J31L*PDstandardNth13gxz - 2*J11L*J33L*PDstandardNth13gxz + 2*J11L*J31L*PDstandardNth13gzz + 
              dJ133L*PDstandardNth1gxx - 2*dJ113L*PDstandardNth1gxz + dJ111L*PDstandardNth1gzz - 
              2*J21L*J23L*PDstandardNth22gxz + 2*J23L*J33L*PDstandardNth23gxx - 2*J23L*J31L*PDstandardNth23gxz - 
              2*J21L*J33L*PDstandardNth23gxz + 2*J21L*J31L*PDstandardNth23gzz + dJ233L*PDstandardNth2gxx - 
              2*dJ213L*PDstandardNth2gxz + dJ211L*PDstandardNth2gzz - 2*J31L*J33L*PDstandardNth33gxz + 
              dJ333L*PDstandardNth3gxx - 2*dJ313L*PDstandardNth3gxz + dJ311L*PDstandardNth3gzz + 
              PDstandardNth11gzz*SQR(J11L) + PDstandardNth11gxx*SQR(J13L) + PDstandardNth22gzz*SQR(J21L) + 
              PDstandardNth22gxx*SQR(J23L) + PDstandardNth33gzz*SQR(J31L) + PDstandardNth33gxx*SQR(J33L)) + 
           gu32*(-(J12L*J13L*PDstandardNth11gxz) - J11L*J13L*PDstandardNth11gyz + J11L*J12L*PDstandardNth11gzz + 
              2*J13L*J23L*PDstandardNth12gxy - J13L*J22L*PDstandardNth12gxz - J12L*J23L*PDstandardNth12gxz - 
              J13L*J21L*PDstandardNth12gyz - J11L*J23L*PDstandardNth12gyz + J12L*J21L*PDstandardNth12gzz + 
              J11L*J22L*PDstandardNth12gzz + 2*J13L*J33L*PDstandardNth13gxy - J13L*J32L*PDstandardNth13gxz - 
              J12L*J33L*PDstandardNth13gxz - J13L*J31L*PDstandardNth13gyz - J11L*J33L*PDstandardNth13gyz + 
              J12L*J31L*PDstandardNth13gzz + J11L*J32L*PDstandardNth13gzz + dJ133L*PDstandardNth1gxy - 
              dJ123L*PDstandardNth1gxz - dJ113L*PDstandardNth1gyz + dJ112L*PDstandardNth1gzz - 
              J22L*J23L*PDstandardNth22gxz - J21L*J23L*PDstandardNth22gyz + J21L*J22L*PDstandardNth22gzz + 
              2*J23L*J33L*PDstandardNth23gxy - J23L*J32L*PDstandardNth23gxz - J22L*J33L*PDstandardNth23gxz - 
              J23L*J31L*PDstandardNth23gyz - J21L*J33L*PDstandardNth23gyz + J22L*J31L*PDstandardNth23gzz + 
              J21L*J32L*PDstandardNth23gzz + dJ233L*PDstandardNth2gxy - dJ223L*PDstandardNth2gxz - 
              dJ213L*PDstandardNth2gyz + dJ212L*PDstandardNth2gzz - J32L*J33L*PDstandardNth33gxz - 
              J31L*J33L*PDstandardNth33gyz + J31L*J32L*PDstandardNth33gzz + dJ333L*PDstandardNth3gxy - 
              dJ323L*PDstandardNth3gxz - dJ313L*PDstandardNth3gyz + dJ312L*PDstandardNth3gzz + 
              PDstandardNth11gxy*SQR(J13L) + PDstandardNth22gxy*SQR(J23L) + PDstandardNth33gxy*SQR(J33L)));
    
    CCTK_REAL const R22  =  khalf*(-(dJ122L*gu11*PDstandardNth1gxx) + 
          gu31*(-4*J12L*J22L*PDstandardNth12gxz - 2*J13L*J21L*PDstandardNth12gyy + 2*J11L*J22L*PDstandardNth12gyz + 
             2*J13L*J32L*PDstandardNth13gxy - 4*J12L*J32L*PDstandardNth13gxz - 2*J13L*J31L*PDstandardNth13gyy + 
             2*J11L*J32L*PDstandardNth13gyz + 2*dJ123L*PDstandardNth1gxy) + 
          2*(gu11*J11L*J12L*PDstandardNth11gxy + gu31*J12L*J13L*PDstandardNth11gxy + gu31*J11L*J12L*PDstandardNth11gyz + 
             gu33*J12L*J13L*PDstandardNth11gyz + gu11*J12L*J21L*PDstandardNth12gxy + gu11*J11L*J22L*PDstandardNth12gxy + 
             gu31*J13L*J22L*PDstandardNth12gxy + gu31*J12L*J23L*PDstandardNth12gxy + gu31*J12L*J21L*PDstandardNth12gyz + 
             gu33*J13L*J22L*PDstandardNth12gyz + gu33*J12L*J23L*PDstandardNth12gyz + gu11*J12L*J31L*PDstandardNth13gxy + 
             gu11*J11L*J32L*PDstandardNth13gxy + gu31*J12L*J33L*PDstandardNth13gxy + gu31*J12L*J31L*PDstandardNth13gyz + 
             gu33*J13L*J32L*PDstandardNth13gyz + gu33*J12L*J33L*PDstandardNth13gyz + dJ112L*gu11*PDstandardNth1gxy) - 
          2*(G111*G122 + G122*G133 + G211*G222 + G222*G233 + G311*G322 + G322*G333 + gu31*J11L*J13L*PDstandardNth11gyy + 
             gu11*J12L*J22L*PDstandardNth12gxx + gu11*J11L*J21L*PDstandardNth12gyy + gu31*J11L*J23L*PDstandardNth12gyy + 
             gu33*J13L*J23L*PDstandardNth12gyy + gu33*J12L*J22L*PDstandardNth12gzz + gu11*J12L*J32L*PDstandardNth13gxx + 
             gu11*J11L*J31L*PDstandardNth13gyy + gu31*J11L*J33L*PDstandardNth13gyy + gu33*J13L*J33L*PDstandardNth13gyy + 
             gu33*J12L*J32L*PDstandardNth13gzz + dJ122L*gu31*PDstandardNth1gxz) - dJ111L*gu11*PDstandardNth1gyy - 
          2*dJ113L*gu31*PDstandardNth1gyy - dJ133L*gu33*PDstandardNth1gyy + 2*dJ112L*gu31*PDstandardNth1gyz + 
          2*dJ123L*gu33*PDstandardNth1gyz - dJ122L*gu33*PDstandardNth1gzz + 2*gu11*J21L*J22L*PDstandardNth22gxy + 
          2*gu31*J22L*J23L*PDstandardNth22gxy - 2*gu31*J21L*J23L*PDstandardNth22gyy + 2*gu31*J21L*J22L*PDstandardNth22gyz + 
          2*gu33*J22L*J23L*PDstandardNth22gyz - 2*gu11*J22L*J32L*PDstandardNth23gxx + 2*gu11*J22L*J31L*PDstandardNth23gxy + 
          2*gu11*J21L*J32L*PDstandardNth23gxy + 2*gu31*J23L*J32L*PDstandardNth23gxy + 2*gu31*J22L*J33L*PDstandardNth23gxy - 
          4*gu31*J22L*J32L*PDstandardNth23gxz - 2*gu11*J21L*J31L*PDstandardNth23gyy - 2*gu31*J23L*J31L*PDstandardNth23gyy - 
          2*gu31*J21L*J33L*PDstandardNth23gyy - 2*gu33*J23L*J33L*PDstandardNth23gyy + 2*gu31*J22L*J31L*PDstandardNth23gyz + 
          2*gu31*J21L*J32L*PDstandardNth23gyz + 2*gu33*J23L*J32L*PDstandardNth23gyz + 2*gu33*J22L*J33L*PDstandardNth23gyz - 
          2*gu33*J22L*J32L*PDstandardNth23gzz - dJ222L*gu11*PDstandardNth2gxx + 2*dJ212L*gu11*PDstandardNth2gxy + 
          2*dJ223L*gu31*PDstandardNth2gxy - 2*dJ222L*gu31*PDstandardNth2gxz - dJ211L*gu11*PDstandardNth2gyy - 
          2*dJ213L*gu31*PDstandardNth2gyy - dJ233L*gu33*PDstandardNth2gyy + 2*dJ212L*gu31*PDstandardNth2gyz + 
          2*dJ223L*gu33*PDstandardNth2gyz - dJ222L*gu33*PDstandardNth2gzz + 2*gu11*J31L*J32L*PDstandardNth33gxy + 
          2*gu31*J32L*J33L*PDstandardNth33gxy - 2*gu31*J31L*J33L*PDstandardNth33gyy + 2*gu31*J31L*J32L*PDstandardNth33gyz + 
          2*gu33*J32L*J33L*PDstandardNth33gyz - dJ322L*gu11*PDstandardNth3gxx + 2*dJ312L*gu11*PDstandardNth3gxy + 
          2*dJ323L*gu31*PDstandardNth3gxy - 2*dJ322L*gu31*PDstandardNth3gxz - dJ311L*gu11*PDstandardNth3gyy - 
          2*dJ313L*gu31*PDstandardNth3gyy - dJ333L*gu33*PDstandardNth3gyy + 2*dJ312L*gu31*PDstandardNth3gyz + 
          2*dJ323L*gu33*PDstandardNth3gyz - dJ322L*gu33*PDstandardNth3gzz + 2*SQR(G112) + 2*SQR(G123) + 2*SQR(G212) + 
          2*SQR(G223) + 2*SQR(G312) + 2*SQR(G323) - gu11*PDstandardNth11gyy*SQR(J11L) - gu11*PDstandardNth11gxx*SQR(J12L) - 
          2*gu31*PDstandardNth11gxz*SQR(J12L) - gu33*PDstandardNth11gzz*SQR(J12L) - gu33*PDstandardNth11gyy*SQR(J13L) - 
          gu11*PDstandardNth22gyy*SQR(J21L) - gu11*PDstandardNth22gxx*SQR(J22L) - 2*gu31*PDstandardNth22gxz*SQR(J22L) - 
          gu33*PDstandardNth22gzz*SQR(J22L) - gu33*PDstandardNth22gyy*SQR(J23L) - gu11*PDstandardNth33gyy*SQR(J31L) - 
          gu11*PDstandardNth33gxx*SQR(J32L) - 2*gu31*PDstandardNth33gxz*SQR(J32L) - gu33*PDstandardNth33gzz*SQR(J32L) - 
          gu33*PDstandardNth33gyy*SQR(J33L));
    
    CCTK_REAL const R23  =  G112*G113 - G111*G123 + G212*G213 + G222*G223 + G223*G233 - G223*(G211 + G222 + G233) + G312*G313 + G322*G323 + 
        G323*G333 - G323*(G311 + G322 + G333) + khalf*
         (gu31*(-(J13L*J22L*PDstandardNth12gxz) + J12L*J23L*PDstandardNth12gxz + 
              (J13L*J22L - J12L*J23L)*PDstandardNth12gxz - J13L*J32L*PDstandardNth13gxz + J12L*J33L*PDstandardNth13gxz + 
              (J13L*J32L - J12L*J33L)*PDstandardNth13gxz) + 
           gu33*(-(J13L*J22L*PDstandardNth12gzz) + J12L*J23L*PDstandardNth12gzz + 
              (J13L*J22L - J12L*J23L)*PDstandardNth12gzz - J13L*J32L*PDstandardNth13gzz + J12L*J33L*PDstandardNth13gzz + 
              (J13L*J32L - J12L*J33L)*PDstandardNth13gzz) + 
           gu11*(-(J12L*J13L*PDstandardNth11gxx) + J11L*J13L*PDstandardNth11gxy + J11L*J12L*PDstandardNth11gxz - 
              J13L*J22L*PDstandardNth12gxx - J12L*J23L*PDstandardNth12gxx + J13L*J21L*PDstandardNth12gxy + 
              J11L*J23L*PDstandardNth12gxy + J12L*J21L*PDstandardNth12gxz + J11L*J22L*PDstandardNth12gxz - 
              2*J11L*J21L*PDstandardNth12gyz - J13L*J32L*PDstandardNth13gxx - J12L*J33L*PDstandardNth13gxx + 
              J13L*J31L*PDstandardNth13gxy + J11L*J33L*PDstandardNth13gxy + J12L*J31L*PDstandardNth13gxz + 
              J11L*J32L*PDstandardNth13gxz - 2*J11L*J31L*PDstandardNth13gyz - dJ123L*PDstandardNth1gxx + 
              dJ113L*PDstandardNth1gxy + dJ112L*PDstandardNth1gxz - dJ111L*PDstandardNth1gyz - 
              J22L*J23L*PDstandardNth22gxx + J21L*J23L*PDstandardNth22gxy + J21L*J22L*PDstandardNth22gxz - 
              J23L*J32L*PDstandardNth23gxx - J22L*J33L*PDstandardNth23gxx + J23L*J31L*PDstandardNth23gxy + 
              J21L*J33L*PDstandardNth23gxy + J22L*J31L*PDstandardNth23gxz + J21L*J32L*PDstandardNth23gxz - 
              2*J21L*J31L*PDstandardNth23gyz - dJ223L*PDstandardNth2gxx + dJ213L*PDstandardNth2gxy + 
              dJ212L*PDstandardNth2gxz - dJ211L*PDstandardNth2gyz - J32L*J33L*PDstandardNth33gxx + 
              J31L*J33L*PDstandardNth33gxy + J31L*J32L*PDstandardNth33gxz - dJ323L*PDstandardNth3gxx + 
              dJ313L*PDstandardNth3gxy + dJ312L*PDstandardNth3gxz - dJ311L*PDstandardNth3gyz - 
              PDstandardNth11gyz*SQR(J11L) - PDstandardNth22gyz*SQR(J21L) - PDstandardNth33gyz*SQR(J31L)) + 
           gu21*(-(J12L*J13L*PDstandardNth11gxy) + J11L*J13L*PDstandardNth11gyy - J11L*J12L*PDstandardNth11gyz - 
              J13L*J22L*PDstandardNth12gxy - J12L*J23L*PDstandardNth12gxy + 2*J12L*J22L*PDstandardNth12gxz + 
              J13L*J21L*PDstandardNth12gyy + J11L*J23L*PDstandardNth12gyy - J12L*J21L*PDstandardNth12gyz - 
              J11L*J22L*PDstandardNth12gyz - J13L*J32L*PDstandardNth13gxy - J12L*J33L*PDstandardNth13gxy + 
              2*J12L*J32L*PDstandardNth13gxz + J13L*J31L*PDstandardNth13gyy + J11L*J33L*PDstandardNth13gyy - 
              J12L*J31L*PDstandardNth13gyz - J11L*J32L*PDstandardNth13gyz - dJ123L*PDstandardNth1gxy + 
              dJ122L*PDstandardNth1gxz + dJ113L*PDstandardNth1gyy - dJ112L*PDstandardNth1gyz - 
              J22L*J23L*PDstandardNth22gxy + J21L*J23L*PDstandardNth22gyy - J21L*J22L*PDstandardNth22gyz - 
              J23L*J32L*PDstandardNth23gxy - J22L*J33L*PDstandardNth23gxy + 2*J22L*J32L*PDstandardNth23gxz + 
              J23L*J31L*PDstandardNth23gyy + J21L*J33L*PDstandardNth23gyy - J22L*J31L*PDstandardNth23gyz - 
              J21L*J32L*PDstandardNth23gyz - dJ223L*PDstandardNth2gxy + dJ222L*PDstandardNth2gxz + 
              dJ213L*PDstandardNth2gyy - dJ212L*PDstandardNth2gyz - J32L*J33L*PDstandardNth33gxy + 
              J31L*J33L*PDstandardNth33gyy - J31L*J32L*PDstandardNth33gyz - dJ323L*PDstandardNth3gxy + 
              dJ322L*PDstandardNth3gxz + dJ313L*PDstandardNth3gyy - dJ312L*PDstandardNth3gyz + 
              PDstandardNth11gxz*SQR(J12L) + PDstandardNth22gxz*SQR(J22L) + PDstandardNth33gxz*SQR(J32L)) + 
           gu31*(-(J12L*J13L*PDstandardNth11gxz) - J11L*J13L*PDstandardNth11gyz + J11L*J12L*PDstandardNth11gzz + 
              2*J13L*J23L*PDstandardNth12gxy - J13L*J22L*PDstandardNth12gxz - J12L*J23L*PDstandardNth12gxz - 
              J13L*J21L*PDstandardNth12gyz - J11L*J23L*PDstandardNth12gyz + J12L*J21L*PDstandardNth12gzz + 
              J11L*J22L*PDstandardNth12gzz + 2*J13L*J33L*PDstandardNth13gxy - J13L*J32L*PDstandardNth13gxz - 
              J12L*J33L*PDstandardNth13gxz - J13L*J31L*PDstandardNth13gyz - J11L*J33L*PDstandardNth13gyz + 
              J12L*J31L*PDstandardNth13gzz + J11L*J32L*PDstandardNth13gzz + dJ133L*PDstandardNth1gxy - 
              dJ123L*PDstandardNth1gxz - dJ113L*PDstandardNth1gyz + dJ112L*PDstandardNth1gzz - 
              J22L*J23L*PDstandardNth22gxz - J21L*J23L*PDstandardNth22gyz + J21L*J22L*PDstandardNth22gzz + 
              2*J23L*J33L*PDstandardNth23gxy - J23L*J32L*PDstandardNth23gxz - J22L*J33L*PDstandardNth23gxz - 
              J23L*J31L*PDstandardNth23gyz - J21L*J33L*PDstandardNth23gyz + J22L*J31L*PDstandardNth23gzz + 
              J21L*J32L*PDstandardNth23gzz + dJ233L*PDstandardNth2gxy - dJ223L*PDstandardNth2gxz - 
              dJ213L*PDstandardNth2gyz + dJ212L*PDstandardNth2gzz - J32L*J33L*PDstandardNth33gxz - 
              J31L*J33L*PDstandardNth33gyz + J31L*J32L*PDstandardNth33gzz + dJ333L*PDstandardNth3gxy - 
              dJ323L*PDstandardNth3gxz - dJ313L*PDstandardNth3gyz + dJ312L*PDstandardNth3gzz + 
              PDstandardNth11gxy*SQR(J13L) + PDstandardNth22gxy*SQR(J23L) + PDstandardNth33gxy*SQR(J33L)) + 
           gu32*(-2*J12L*J13L*PDstandardNth11gyz + 2*J13L*J23L*PDstandardNth12gyy - 2*J13L*J22L*PDstandardNth12gyz - 
              2*J12L*J23L*PDstandardNth12gyz + 2*J12L*J22L*PDstandardNth12gzz + 2*J13L*J33L*PDstandardNth13gyy - 
              2*J13L*J32L*PDstandardNth13gyz - 2*J12L*J33L*PDstandardNth13gyz + 2*J12L*J32L*PDstandardNth13gzz + 
              dJ133L*PDstandardNth1gyy - 2*dJ123L*PDstandardNth1gyz + dJ122L*PDstandardNth1gzz - 
              2*J22L*J23L*PDstandardNth22gyz + 2*J23L*J33L*PDstandardNth23gyy - 2*J23L*J32L*PDstandardNth23gyz - 
              2*J22L*J33L*PDstandardNth23gyz + 2*J22L*J32L*PDstandardNth23gzz + dJ233L*PDstandardNth2gyy - 
              2*dJ223L*PDstandardNth2gyz + dJ222L*PDstandardNth2gzz - 2*J32L*J33L*PDstandardNth33gyz + 
              dJ333L*PDstandardNth3gyy - 2*dJ323L*PDstandardNth3gyz + dJ322L*PDstandardNth3gzz + 
              PDstandardNth11gzz*SQR(J12L) + PDstandardNth11gyy*SQR(J13L) + PDstandardNth22gzz*SQR(J22L) + 
              PDstandardNth22gyy*SQR(J23L) + PDstandardNth33gzz*SQR(J32L) + PDstandardNth33gyy*SQR(J33L)));
    
    CCTK_REAL const R33  =  khalf*(-(dJ133L*gu11*PDstandardNth1gxx) + 
          gu21*(-4*J13L*J23L*PDstandardNth12gxy + 2*J13L*J22L*PDstandardNth12gxz + 2*J13L*J21L*PDstandardNth12gyz - 
             2*J12L*J21L*PDstandardNth12gzz - 4*J13L*J33L*PDstandardNth13gxy + 2*J13L*J32L*PDstandardNth13gxz + 
             2*J13L*J31L*PDstandardNth13gyz - 2*J12L*J31L*PDstandardNth13gzz - 2*dJ133L*PDstandardNth1gxy) + 
          2*dJ123L*gu21*PDstandardNth1gxz + 2*(gu11*J11L*J13L*PDstandardNth11gxz + gu21*J12L*J13L*PDstandardNth11gxz + 
             gu21*J11L*J13L*PDstandardNth11gyz + gu22*J12L*J13L*PDstandardNth11gyz + gu11*J13L*J21L*PDstandardNth12gxz + 
             gu11*J11L*J23L*PDstandardNth12gxz + gu21*J12L*J23L*PDstandardNth12gxz + gu22*J13L*J22L*PDstandardNth12gyz + 
             gu21*J11L*J23L*PDstandardNth12gyz + gu22*J12L*J23L*PDstandardNth12gyz + gu11*J13L*J31L*PDstandardNth13gxz + 
             gu11*J11L*J33L*PDstandardNth13gxz + gu21*J12L*J33L*PDstandardNth13gxz + gu22*J13L*J32L*PDstandardNth13gyz + 
             gu21*J11L*J33L*PDstandardNth13gyz + gu22*J12L*J33L*PDstandardNth13gyz + dJ113L*gu11*PDstandardNth1gxz) - 
          dJ133L*gu22*PDstandardNth1gyy + 2*dJ113L*gu21*PDstandardNth1gyz + 2*dJ123L*gu22*PDstandardNth1gyz - 
          dJ111L*gu11*PDstandardNth1gzz - dJ122L*gu22*PDstandardNth1gzz - 
          2*(G111*G133 + G122*G133 + G211*G233 + G222*G233 + G311*G333 + G322*G333 + gu21*J11L*J12L*PDstandardNth11gzz + 
             gu11*J13L*J23L*PDstandardNth12gxx + gu22*J13L*J23L*PDstandardNth12gyy + gu11*J11L*J21L*PDstandardNth12gzz + 
             gu21*J11L*J22L*PDstandardNth12gzz + gu22*J12L*J22L*PDstandardNth12gzz + gu11*J13L*J33L*PDstandardNth13gxx + 
             gu22*J13L*J33L*PDstandardNth13gyy + gu11*J11L*J31L*PDstandardNth13gzz + gu21*J11L*J32L*PDstandardNth13gzz + 
             gu22*J12L*J32L*PDstandardNth13gzz + dJ112L*gu21*PDstandardNth1gzz) + 2*gu11*J21L*J23L*PDstandardNth22gxz + 
          2*gu21*J22L*J23L*PDstandardNth22gxz + 2*gu21*J21L*J23L*PDstandardNth22gyz + 2*gu22*J22L*J23L*PDstandardNth22gyz - 
          2*gu21*J21L*J22L*PDstandardNth22gzz - 2*gu11*J23L*J33L*PDstandardNth23gxx - 4*gu21*J23L*J33L*PDstandardNth23gxy + 
          2*gu11*J23L*J31L*PDstandardNth23gxz + 2*gu21*J23L*J32L*PDstandardNth23gxz + 2*gu11*J21L*J33L*PDstandardNth23gxz + 
          2*gu21*J22L*J33L*PDstandardNth23gxz - 2*gu22*J23L*J33L*PDstandardNth23gyy + 2*gu21*J23L*J31L*PDstandardNth23gyz + 
          2*gu22*J23L*J32L*PDstandardNth23gyz + 2*gu21*J21L*J33L*PDstandardNth23gyz + 2*gu22*J22L*J33L*PDstandardNth23gyz - 
          2*gu11*J21L*J31L*PDstandardNth23gzz - 2*gu21*J22L*J31L*PDstandardNth23gzz - 2*gu21*J21L*J32L*PDstandardNth23gzz - 
          2*gu22*J22L*J32L*PDstandardNth23gzz - dJ233L*gu11*PDstandardNth2gxx - 2*dJ233L*gu21*PDstandardNth2gxy + 
          2*dJ213L*gu11*PDstandardNth2gxz + 2*dJ223L*gu21*PDstandardNth2gxz - dJ233L*gu22*PDstandardNth2gyy + 
          2*dJ213L*gu21*PDstandardNth2gyz + 2*dJ223L*gu22*PDstandardNth2gyz - dJ211L*gu11*PDstandardNth2gzz - 
          2*dJ212L*gu21*PDstandardNth2gzz - dJ222L*gu22*PDstandardNth2gzz + 2*gu11*J31L*J33L*PDstandardNth33gxz + 
          2*gu21*J32L*J33L*PDstandardNth33gxz + 2*gu21*J31L*J33L*PDstandardNth33gyz + 2*gu22*J32L*J33L*PDstandardNth33gyz - 
          2*gu21*J31L*J32L*PDstandardNth33gzz - dJ333L*gu11*PDstandardNth3gxx - 2*dJ333L*gu21*PDstandardNth3gxy + 
          2*dJ313L*gu11*PDstandardNth3gxz + 2*dJ323L*gu21*PDstandardNth3gxz - dJ333L*gu22*PDstandardNth3gyy + 
          2*dJ313L*gu21*PDstandardNth3gyz + 2*dJ323L*gu22*PDstandardNth3gyz - dJ311L*gu11*PDstandardNth3gzz - 
          2*dJ312L*gu21*PDstandardNth3gzz - dJ322L*gu22*PDstandardNth3gzz + 2*SQR(G113) + 2*SQR(G123) + 2*SQR(G213) + 
          2*SQR(G223) + 2*SQR(G313) + 2*SQR(G323) - gu11*PDstandardNth11gzz*SQR(J11L) - gu22*PDstandardNth11gzz*SQR(J12L) - 
          gu11*PDstandardNth11gxx*SQR(J13L) - 2*gu21*PDstandardNth11gxy*SQR(J13L) - gu22*PDstandardNth11gyy*SQR(J13L) - 
          gu11*PDstandardNth22gzz*SQR(J21L) - gu22*PDstandardNth22gzz*SQR(J22L) - gu11*PDstandardNth22gxx*SQR(J23L) - 
          2*gu21*PDstandardNth22gxy*SQR(J23L) - gu22*PDstandardNth22gyy*SQR(J23L) - gu11*PDstandardNth33gzz*SQR(J31L) - 
          gu22*PDstandardNth33gzz*SQR(J32L) - gu11*PDstandardNth33gxx*SQR(J33L) - 2*gu21*PDstandardNth33gxy*SQR(J33L) - 
          gu22*PDstandardNth33gyy*SQR(J33L));
    
    CCTK_REAL const trR  =  gu11*R11 + gu22*R22 + 2*(gu21*R12 + gu31*R13 + gu32*R23) + gu33*R33;
    
    CCTK_REAL const Km11  =  gu11*kxxL + gu21*kxyL + gu31*kxzL;
    
    CCTK_REAL const Km21  =  gu21*kxxL + gu22*kxyL + gu32*kxzL;
    
    CCTK_REAL const Km31  =  gu31*kxxL + gu32*kxyL + gu33*kxzL;
    
    CCTK_REAL const Km12  =  gu11*kxyL + gu21*kyyL + gu31*kyzL;
    
    CCTK_REAL const Km22  =  gu21*kxyL + gu22*kyyL + gu32*kyzL;
    
    CCTK_REAL const Km32  =  gu31*kxyL + gu32*kyyL + gu33*kyzL;
    
    CCTK_REAL const Km13  =  gu11*kxzL + gu21*kyzL + gu31*kzzL;
    
    CCTK_REAL const Km23  =  gu21*kxzL + gu22*kyzL + gu32*kzzL;
    
    CCTK_REAL const Km33  =  gu31*kxzL + gu32*kyzL + gu33*kzzL;
    
    CCTK_REAL const trK  =  Km11 + Km22 + Km33;
    
    CCTK_REAL const rho  =  pow(alpL,-2)*(eTttL - 2*(betayL*eTtyL + betazL*eTtzL) + 
          2*(betaxL*(-eTtxL + betayL*eTxyL + betazL*eTxzL) + betayL*betazL*eTyzL) + eTxxL*SQR(betaxL) + eTyyL*SQR(betayL) + 
          eTzzL*SQR(betazL));
    
    CCTK_REAL const S1  =  (-eTtxL + betaxL*eTxxL + betayL*eTxyL + betazL*eTxzL)*INV(alpL);
    
    CCTK_REAL const S2  =  (-eTtyL + betaxL*eTxyL + betayL*eTyyL + betazL*eTyzL)*INV(alpL);
    
    CCTK_REAL const S3  =  (-eTtzL + betaxL*eTxzL + betayL*eTyzL + betazL*eTzzL)*INV(alpL);
    
    CCTK_REAL const HL  =  -2*(Km12*Km21 + Km13*Km31 + Km23*Km32) - 50.26548245743669181540229413247204614715*rho + trR - SQR(Km11) - 
        SQR(Km22) - SQR(Km33) + SQR(trK);
    
    CCTK_REAL const M1L  =  gu21*(-(G112*kxxL) + G111*kxyL - G212*kxyL - G312*kxzL + G211*kyyL + G311*kyzL + J12L*PDstandardNth1kxx - 
           J11L*PDstandardNth1kxy + J22L*PDstandardNth2kxx - J21L*PDstandardNth2kxy + J32L*PDstandardNth3kxx - 
           J31L*PDstandardNth3kxy) + gu31*(-(G113*kxxL) - G213*kxyL + G111*kxzL - G313*kxzL + G211*kyzL + G311*kzzL + 
           J13L*PDstandardNth1kxx - J11L*PDstandardNth1kxz + J23L*PDstandardNth2kxx - J21L*PDstandardNth2kxz + 
           J33L*PDstandardNth3kxx - J31L*PDstandardNth3kxz) + 
        gu22*(-(G122*kxxL) + G112*kxyL - G222*kxyL - G322*kxzL + G212*kyyL + G312*kyzL + J12L*PDstandardNth1kxy - 
           J11L*PDstandardNth1kyy + J22L*PDstandardNth2kxy - J21L*PDstandardNth2kyy + J32L*PDstandardNth3kxy - 
           J31L*PDstandardNth3kyy) + gu32*(G113*kxyL + G112*kxzL + G213*kyyL + (G212 + G313)*kyzL + G312*kzzL + 
           J13L*PDstandardNth1kxy + J12L*PDstandardNth1kxz + J23L*PDstandardNth2kxy + J22L*PDstandardNth2kxz + 
           J33L*PDstandardNth3kxy + J32L*PDstandardNth3kxz - 
           2*(G123*kxxL + G223*kxyL + G323*kxzL + J11L*PDstandardNth1kyz + J21L*PDstandardNth2kyz + J31L*PDstandardNth3kyz))
          + gu33*(-(G133*kxxL) - G233*kxyL + G113*kxzL - G333*kxzL + G213*kyzL + G313*kzzL + J13L*PDstandardNth1kxz - 
           J11L*PDstandardNth1kzz + J23L*PDstandardNth2kxz - J21L*PDstandardNth2kzz + J33L*PDstandardNth3kxz - 
           J31L*PDstandardNth3kzz) - 25.13274122871834590770114706623602307358*S1;
    
    CCTK_REAL const M2L  =  gu11*(G112*kxxL + (-G111 + G212)*kxyL + G312*kxzL - G211*kyyL - G311*kyzL - J12L*PDstandardNth1kxx + 
           J11L*PDstandardNth1kxy - J22L*PDstandardNth2kxx + J21L*PDstandardNth2kxy - J32L*PDstandardNth3kxx + 
           J31L*PDstandardNth3kxy) + gu21*(G122*kxxL + (-G112 + G222)*kxyL + G322*kxzL - G212*kyyL - G312*kyzL - 
           J12L*PDstandardNth1kxy + J11L*PDstandardNth1kyy - J22L*PDstandardNth2kxy + J21L*PDstandardNth2kyy - 
           J32L*PDstandardNth3kxy + J31L*PDstandardNth3kyy) + 
        gu31*(G123*kxxL + (-2*G113 + G223)*kxyL + (G112 + G323)*kxzL + G212*kyzL + G312*kzzL + J13L*PDstandardNth1kxy + 
           J11L*PDstandardNth1kyz + J23L*PDstandardNth2kxy + J21L*PDstandardNth2kyz + J33L*PDstandardNth3kxy - 
           2*(G213*kyyL + G313*kyzL + J12L*PDstandardNth1kxz + J22L*PDstandardNth2kxz + J32L*PDstandardNth3kxz) + 
           J31L*PDstandardNth3kyz) + gu32*(-(G123*kxyL) + G122*kxzL - G223*kyyL + G222*kyzL - G323*kyzL + G322*kzzL + 
           J13L*PDstandardNth1kyy - J12L*PDstandardNth1kyz + J23L*PDstandardNth2kyy - J22L*PDstandardNth2kyz + 
           J33L*PDstandardNth3kyy - J32L*PDstandardNth3kyz) + 
        gu33*(-(G133*kxyL) + G123*kxzL - G233*kyyL + G223*kyzL - G333*kyzL + G323*kzzL + J13L*PDstandardNth1kyz - 
           J12L*PDstandardNth1kzz + J23L*PDstandardNth2kyz - J22L*PDstandardNth2kzz + J33L*PDstandardNth3kyz - 
           J32L*PDstandardNth3kzz) - 25.13274122871834590770114706623602307358*S2;
    
    CCTK_REAL const M3L  =  gu11*(G113*kxxL + G213*kxyL + (-G111 + G313)*kxzL - G211*kyzL - G311*kzzL - J13L*PDstandardNth1kxx + 
           J11L*PDstandardNth1kxz - J23L*PDstandardNth2kxx + J21L*PDstandardNth2kxz - J33L*PDstandardNth3kxx + 
           J31L*PDstandardNth3kxz) + gu21*(G123*kxxL + (G113 + G223)*kxyL + (-2*G112 + G323)*kxzL + G213*kyyL + 
           (-2*G212 + G313)*kyzL + J12L*PDstandardNth1kxz + J11L*PDstandardNth1kyz + J22L*PDstandardNth2kxz + 
           J21L*PDstandardNth2kyz - 2*(G312*kzzL + J13L*PDstandardNth1kxy + J23L*PDstandardNth2kxy + 
              J33L*PDstandardNth3kxy) + J32L*PDstandardNth3kxz + J31L*PDstandardNth3kyz) + 
        gu22*(G123*kxyL - G122*kxzL + G223*kyyL - G222*kyzL + G323*kyzL - G322*kzzL - J13L*PDstandardNth1kyy + 
           J12L*PDstandardNth1kyz - J23L*PDstandardNth2kyy + J22L*PDstandardNth2kyz - J33L*PDstandardNth3kyy + 
           J32L*PDstandardNth3kyz) + gu31*(G133*kxxL + G233*kxyL + (-G113 + G333)*kxzL - G213*kyzL - G313*kzzL - 
           J13L*PDstandardNth1kxz + J11L*PDstandardNth1kzz - J23L*PDstandardNth2kxz + J21L*PDstandardNth2kzz - 
           J33L*PDstandardNth3kxz + J31L*PDstandardNth3kzz) + 
        gu32*(G133*kxyL - G123*kxzL + G233*kyyL - G223*kyzL + G333*kyzL - G323*kzzL - J13L*PDstandardNth1kyz + 
           J12L*PDstandardNth1kzz - J23L*PDstandardNth2kyz + J22L*PDstandardNth2kzz - J33L*PDstandardNth3kyz + 
           J32L*PDstandardNth3kzz) - 25.13274122871834590770114706623602307358*S3;
    
    
    /* Copy local copies back to grid functions */
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_ADMConstraints_MP);
}

void ML_ADMConstraints_MP(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADMConstraints_MP_Body);
}
