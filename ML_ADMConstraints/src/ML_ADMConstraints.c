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

void ML_ADMConstraints_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMConstraints_Body");
  }
  
  if (cctk_iteration % ML_ADMConstraints_calc_every != ML_ADMConstraints_calc_offset)
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
  LC_LOOP3 (ML_ADMConstraints,
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
    // CCTK_REAL PDstandardNth22gxx = INITVALUE;
    // CCTK_REAL PDstandardNth33gxx = INITVALUE;
    // CCTK_REAL PDstandardNth23gxx = INITVALUE;
    // CCTK_REAL PDstandardNth1gxy = INITVALUE;
    // CCTK_REAL PDstandardNth2gxy = INITVALUE;
    // CCTK_REAL PDstandardNth3gxy = INITVALUE;
    // CCTK_REAL PDstandardNth33gxy = INITVALUE;
    // CCTK_REAL PDstandardNth12gxy = INITVALUE;
    // CCTK_REAL PDstandardNth13gxy = INITVALUE;
    // CCTK_REAL PDstandardNth23gxy = INITVALUE;
    // CCTK_REAL PDstandardNth1gxz = INITVALUE;
    // CCTK_REAL PDstandardNth2gxz = INITVALUE;
    // CCTK_REAL PDstandardNth3gxz = INITVALUE;
    // CCTK_REAL PDstandardNth22gxz = INITVALUE;
    // CCTK_REAL PDstandardNth12gxz = INITVALUE;
    // CCTK_REAL PDstandardNth13gxz = INITVALUE;
    // CCTK_REAL PDstandardNth23gxz = INITVALUE;
    // CCTK_REAL PDstandardNth1gyy = INITVALUE;
    // CCTK_REAL PDstandardNth2gyy = INITVALUE;
    // CCTK_REAL PDstandardNth3gyy = INITVALUE;
    // CCTK_REAL PDstandardNth11gyy = INITVALUE;
    // CCTK_REAL PDstandardNth33gyy = INITVALUE;
    // CCTK_REAL PDstandardNth13gyy = INITVALUE;
    // CCTK_REAL PDstandardNth1gyz = INITVALUE;
    // CCTK_REAL PDstandardNth2gyz = INITVALUE;
    // CCTK_REAL PDstandardNth3gyz = INITVALUE;
    // CCTK_REAL PDstandardNth11gyz = INITVALUE;
    // CCTK_REAL PDstandardNth12gyz = INITVALUE;
    // CCTK_REAL PDstandardNth13gyz = INITVALUE;
    // CCTK_REAL PDstandardNth23gyz = INITVALUE;
    // CCTK_REAL PDstandardNth1gzz = INITVALUE;
    // CCTK_REAL PDstandardNth2gzz = INITVALUE;
    // CCTK_REAL PDstandardNth3gzz = INITVALUE;
    // CCTK_REAL PDstandardNth11gzz = INITVALUE;
    // CCTK_REAL PDstandardNth22gzz = INITVALUE;
    // CCTK_REAL PDstandardNth12gzz = INITVALUE;
    // CCTK_REAL PDstandardNth2kxx = INITVALUE;
    // CCTK_REAL PDstandardNth3kxx = INITVALUE;
    // CCTK_REAL PDstandardNth1kxy = INITVALUE;
    // CCTK_REAL PDstandardNth2kxy = INITVALUE;
    // CCTK_REAL PDstandardNth3kxy = INITVALUE;
    // CCTK_REAL PDstandardNth1kxz = INITVALUE;
    // CCTK_REAL PDstandardNth2kxz = INITVALUE;
    // CCTK_REAL PDstandardNth3kxz = INITVALUE;
    // CCTK_REAL PDstandardNth1kyy = INITVALUE;
    // CCTK_REAL PDstandardNth3kyy = INITVALUE;
    // CCTK_REAL PDstandardNth1kyz = INITVALUE;
    // CCTK_REAL PDstandardNth2kyz = INITVALUE;
    // CCTK_REAL PDstandardNth3kyz = INITVALUE;
    // CCTK_REAL PDstandardNth1kzz = INITVALUE;
    // CCTK_REAL PDstandardNth2kzz = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL const alpL = alp[index];
    CCTK_REAL const betaxL = betax[index];
    CCTK_REAL const betayL = betay[index];
    CCTK_REAL const betazL = betaz[index];
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
    CCTK_REAL const PDstandardNth22gxx = PDstandardNth22(gxx, i, j, k);
    CCTK_REAL const PDstandardNth33gxx = PDstandardNth33(gxx, i, j, k);
    CCTK_REAL const PDstandardNth23gxx = PDstandardNth23(gxx, i, j, k);
    CCTK_REAL const PDstandardNth1gxy = PDstandardNth1(gxy, i, j, k);
    CCTK_REAL const PDstandardNth2gxy = PDstandardNth2(gxy, i, j, k);
    CCTK_REAL const PDstandardNth3gxy = PDstandardNth3(gxy, i, j, k);
    CCTK_REAL const PDstandardNth33gxy = PDstandardNth33(gxy, i, j, k);
    CCTK_REAL const PDstandardNth12gxy = PDstandardNth12(gxy, i, j, k);
    CCTK_REAL const PDstandardNth13gxy = PDstandardNth13(gxy, i, j, k);
    CCTK_REAL const PDstandardNth23gxy = PDstandardNth23(gxy, i, j, k);
    CCTK_REAL const PDstandardNth1gxz = PDstandardNth1(gxz, i, j, k);
    CCTK_REAL const PDstandardNth2gxz = PDstandardNth2(gxz, i, j, k);
    CCTK_REAL const PDstandardNth3gxz = PDstandardNth3(gxz, i, j, k);
    CCTK_REAL const PDstandardNth22gxz = PDstandardNth22(gxz, i, j, k);
    CCTK_REAL const PDstandardNth12gxz = PDstandardNth12(gxz, i, j, k);
    CCTK_REAL const PDstandardNth13gxz = PDstandardNth13(gxz, i, j, k);
    CCTK_REAL const PDstandardNth23gxz = PDstandardNth23(gxz, i, j, k);
    CCTK_REAL const PDstandardNth1gyy = PDstandardNth1(gyy, i, j, k);
    CCTK_REAL const PDstandardNth2gyy = PDstandardNth2(gyy, i, j, k);
    CCTK_REAL const PDstandardNth3gyy = PDstandardNth3(gyy, i, j, k);
    CCTK_REAL const PDstandardNth11gyy = PDstandardNth11(gyy, i, j, k);
    CCTK_REAL const PDstandardNth33gyy = PDstandardNth33(gyy, i, j, k);
    CCTK_REAL const PDstandardNth13gyy = PDstandardNth13(gyy, i, j, k);
    CCTK_REAL const PDstandardNth1gyz = PDstandardNth1(gyz, i, j, k);
    CCTK_REAL const PDstandardNth2gyz = PDstandardNth2(gyz, i, j, k);
    CCTK_REAL const PDstandardNth3gyz = PDstandardNth3(gyz, i, j, k);
    CCTK_REAL const PDstandardNth11gyz = PDstandardNth11(gyz, i, j, k);
    CCTK_REAL const PDstandardNth12gyz = PDstandardNth12(gyz, i, j, k);
    CCTK_REAL const PDstandardNth13gyz = PDstandardNth13(gyz, i, j, k);
    CCTK_REAL const PDstandardNth23gyz = PDstandardNth23(gyz, i, j, k);
    CCTK_REAL const PDstandardNth1gzz = PDstandardNth1(gzz, i, j, k);
    CCTK_REAL const PDstandardNth2gzz = PDstandardNth2(gzz, i, j, k);
    CCTK_REAL const PDstandardNth3gzz = PDstandardNth3(gzz, i, j, k);
    CCTK_REAL const PDstandardNth11gzz = PDstandardNth11(gzz, i, j, k);
    CCTK_REAL const PDstandardNth22gzz = PDstandardNth22(gzz, i, j, k);
    CCTK_REAL const PDstandardNth12gzz = PDstandardNth12(gzz, i, j, k);
    CCTK_REAL const PDstandardNth2kxx = PDstandardNth2(kxx, i, j, k);
    CCTK_REAL const PDstandardNth3kxx = PDstandardNth3(kxx, i, j, k);
    CCTK_REAL const PDstandardNth1kxy = PDstandardNth1(kxy, i, j, k);
    CCTK_REAL const PDstandardNth2kxy = PDstandardNth2(kxy, i, j, k);
    CCTK_REAL const PDstandardNth3kxy = PDstandardNth3(kxy, i, j, k);
    CCTK_REAL const PDstandardNth1kxz = PDstandardNth1(kxz, i, j, k);
    CCTK_REAL const PDstandardNth2kxz = PDstandardNth2(kxz, i, j, k);
    CCTK_REAL const PDstandardNth3kxz = PDstandardNth3(kxz, i, j, k);
    CCTK_REAL const PDstandardNth1kyy = PDstandardNth1(kyy, i, j, k);
    CCTK_REAL const PDstandardNth3kyy = PDstandardNth3(kyy, i, j, k);
    CCTK_REAL const PDstandardNth1kyz = PDstandardNth1(kyz, i, j, k);
    CCTK_REAL const PDstandardNth2kyz = PDstandardNth2(kyz, i, j, k);
    CCTK_REAL const PDstandardNth3kyz = PDstandardNth3(kyz, i, j, k);
    CCTK_REAL const PDstandardNth1kzz = PDstandardNth1(kzz, i, j, k);
    CCTK_REAL const PDstandardNth2kzz = PDstandardNth2(kzz, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL const detg  =  2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - SQR(gxyL)) - gyyL*SQR(gxzL) - gxxL*SQR(gyzL);
    
    CCTK_REAL const gu11  =  INV(detg)*(gyyL*gzzL - SQR(gyzL));
    
    CCTK_REAL const gu21  =  (gxzL*gyzL - gxyL*gzzL)*INV(detg);
    
    CCTK_REAL const gu31  =  (-(gxzL*gyyL) + gxyL*gyzL)*INV(detg);
    
    CCTK_REAL const gu22  =  INV(detg)*(gxxL*gzzL - SQR(gxzL));
    
    CCTK_REAL const gu32  =  (gxyL*gxzL - gxxL*gyzL)*INV(detg);
    
    CCTK_REAL const gu33  =  INV(detg)*(gxxL*gyyL - SQR(gxyL));
    
    CCTK_REAL const G111  =  khalf*(gu11*PDstandardNth1gxx + 2*(gu21*PDstandardNth1gxy + gu31*PDstandardNth1gxz) - gu21*PDstandardNth2gxx - 
          gu31*PDstandardNth3gxx);
    
    CCTK_REAL const G211  =  khalf*(gu21*PDstandardNth1gxx + 2*(gu22*PDstandardNth1gxy + gu32*PDstandardNth1gxz) - gu22*PDstandardNth2gxx - 
          gu32*PDstandardNth3gxx);
    
    CCTK_REAL const G311  =  khalf*(gu31*PDstandardNth1gxx + 2*(gu32*PDstandardNth1gxy + gu33*PDstandardNth1gxz) - gu32*PDstandardNth2gxx - 
          gu33*PDstandardNth3gxx);
    
    CCTK_REAL const G112  =  khalf*(gu21*PDstandardNth1gyy + gu11*PDstandardNth2gxx + 
          gu31*(PDstandardNth1gyz + PDstandardNth2gxz - PDstandardNth3gxy));
    
    CCTK_REAL const G212  =  khalf*(gu22*PDstandardNth1gyy + gu21*PDstandardNth2gxx + 
          gu32*(PDstandardNth1gyz + PDstandardNth2gxz - PDstandardNth3gxy));
    
    CCTK_REAL const G312  =  khalf*(gu32*PDstandardNth1gyy + gu31*PDstandardNth2gxx + 
          gu33*(PDstandardNth1gyz + PDstandardNth2gxz - PDstandardNth3gxy));
    
    CCTK_REAL const G113  =  khalf*(gu31*PDstandardNth1gzz + gu11*PDstandardNth3gxx + 
          gu21*(PDstandardNth1gyz - PDstandardNth2gxz + PDstandardNth3gxy));
    
    CCTK_REAL const G213  =  khalf*(gu32*PDstandardNth1gzz + gu21*PDstandardNth3gxx + 
          gu22*(PDstandardNth1gyz - PDstandardNth2gxz + PDstandardNth3gxy));
    
    CCTK_REAL const G313  =  khalf*(gu33*PDstandardNth1gzz + gu31*PDstandardNth3gxx + 
          gu32*(PDstandardNth1gyz - PDstandardNth2gxz + PDstandardNth3gxy));
    
    CCTK_REAL const G122  =  khalf*(gu11*(-PDstandardNth1gyy + 2*PDstandardNth2gxy) + gu21*PDstandardNth2gyy + 
          gu31*(2*PDstandardNth2gyz - PDstandardNth3gyy));
    
    CCTK_REAL const G222  =  khalf*(gu21*(-PDstandardNth1gyy + 2*PDstandardNth2gxy) + gu22*PDstandardNth2gyy + 
          gu32*(2*PDstandardNth2gyz - PDstandardNth3gyy));
    
    CCTK_REAL const G322  =  khalf*(gu31*(-PDstandardNth1gyy + 2*PDstandardNth2gxy) + gu32*PDstandardNth2gyy + 
          gu33*(2*PDstandardNth2gyz - PDstandardNth3gyy));
    
    CCTK_REAL const G123  =  khalf*(gu31*PDstandardNth2gzz + gu11*(-PDstandardNth1gyz + PDstandardNth2gxz + PDstandardNth3gxy) + 
          gu21*PDstandardNth3gyy);
    
    CCTK_REAL const G223  =  khalf*(gu32*PDstandardNth2gzz + gu21*(-PDstandardNth1gyz + PDstandardNth2gxz + PDstandardNth3gxy) + 
          gu22*PDstandardNth3gyy);
    
    CCTK_REAL const G323  =  khalf*(gu33*PDstandardNth2gzz + gu31*(-PDstandardNth1gyz + PDstandardNth2gxz + PDstandardNth3gxy) + 
          gu32*PDstandardNth3gyy);
    
    CCTK_REAL const G133  =  khalf*(-(gu11*PDstandardNth1gzz) - gu21*PDstandardNth2gzz + 2*gu11*PDstandardNth3gxz + 
          2*gu21*PDstandardNth3gyz + gu31*PDstandardNth3gzz);
    
    CCTK_REAL const G233  =  khalf*(-(gu21*PDstandardNth1gzz) - gu22*PDstandardNth2gzz + 2*gu21*PDstandardNth3gxz + 
          2*gu22*PDstandardNth3gyz + gu32*PDstandardNth3gzz);
    
    CCTK_REAL const G333  =  khalf*(-(gu31*PDstandardNth1gzz) - gu32*PDstandardNth2gzz + 2*gu31*PDstandardNth3gxz + 
          2*gu32*PDstandardNth3gyz + gu33*PDstandardNth3gzz);
    
    CCTK_REAL const R11  =  khalf*(-(gu22*PDstandardNth11gyy) - 2*(G111*G122 + G111*G133 + G211*G222 + G211*G233 + G311*G322 + G311*G333 + 
             gu32*PDstandardNth11gyz) - gu33*PDstandardNth11gzz + 2*gu22*PDstandardNth12gxy + 2*gu32*PDstandardNth12gxz + 
          2*gu32*PDstandardNth13gxy + 2*gu33*PDstandardNth13gxz - gu22*PDstandardNth22gxx - 2*gu32*PDstandardNth23gxx - 
          gu33*PDstandardNth33gxx + 2*SQR(G112) + 2*SQR(G113) + 2*SQR(G212) + 2*SQR(G213) + 2*SQR(G312) + 2*SQR(G313));
    
    CCTK_REAL const R12  =  khalf*(2*(G113*G123 + G213*G223 + G313*G323) - 2*(G112*G133 + G212*G233 + G312*G333 + gu21*PDstandardNth12gxy) - 
          gu32*PDstandardNth12gyz - gu33*PDstandardNth12gzz + 
          gu31*(PDstandardNth11gyz - PDstandardNth12gxz - PDstandardNth13gxy) + gu32*PDstandardNth13gyy + 
          gu33*PDstandardNth13gyz + gu21*(PDstandardNth11gyy + PDstandardNth22gxx) + gu32*PDstandardNth22gxz + 
          gu31*PDstandardNth23gxx - gu32*PDstandardNth23gxy + gu33*PDstandardNth23gxz - gu33*PDstandardNth33gxy);
    
    CCTK_REAL const R13  =  khalf*(2*(G112*G123 + G212*G223 + G312*G323) - 2*(G113*G122 + G213*G222 + G313*G322 + gu31*PDstandardNth13gxz) + 
          gu21*(PDstandardNth11gyz - PDstandardNth12gxz - PDstandardNth13gxy + PDstandardNth23gxx) + 
          gu22*(PDstandardNth12gyz - PDstandardNth13gyy - PDstandardNth22gxz + PDstandardNth23gxy) + 
          gu31*(PDstandardNth11gzz + PDstandardNth33gxx) + 
          gu32*(PDstandardNth12gzz - PDstandardNth13gyz - PDstandardNth23gxz + PDstandardNth33gxy));
    
    CCTK_REAL const R22  =  khalf*(-2*(G122*(G111 + G133) + G222*(G211 + G233) + G322*(G311 + G333) + gu31*PDstandardNth13gyy) + 
          gu11*(-PDstandardNth11gyy + 2*PDstandardNth12gxy - PDstandardNth22gxx) + 
          gu31*(-2*PDstandardNth22gxz + 2*(PDstandardNth12gyz + PDstandardNth23gxy)) + 
          gu33*(-PDstandardNth22gzz + 2*PDstandardNth23gyz - PDstandardNth33gyy) + 
          2*(SQR(G112) + SQR(G123) + SQR(G212) + SQR(G223) + SQR(G312) + SQR(G323)));
    
    CCTK_REAL const R23  =  khalf*(2*(G112*G113 + G212*G213 + G312*G313) + 
          gu11*(-PDstandardNth11gyz + PDstandardNth12gxz + PDstandardNth13gxy - PDstandardNth23gxx) + 
          gu21*(-PDstandardNth12gyz + PDstandardNth13gyy + PDstandardNth22gxz - PDstandardNth23gxy) - 
          2*(G111*G123 + G211*G223 + G311*G323 + gu32*PDstandardNth23gyz) + 
          gu31*(PDstandardNth12gzz - PDstandardNth13gyz - PDstandardNth23gxz + PDstandardNth33gxy) + 
          gu32*(PDstandardNth22gzz + PDstandardNth33gyy));
    
    CCTK_REAL const R33  =  khalf*(gu11*(-PDstandardNth11gzz + 2*PDstandardNth13gxz - PDstandardNth33gxx) - 
          2*((G111 + G122)*G133 + (G211 + G222)*G233 + (G311 + G322)*G333 + 
             gu21*(PDstandardNth12gzz + PDstandardNth33gxy)) + 
          gu22*(-PDstandardNth22gzz + 2*PDstandardNth23gyz - PDstandardNth33gyy) + 
          2*(gu21*(PDstandardNth13gyz + PDstandardNth23gxz) + SQR(G113) + SQR(G123) + SQR(G213) + SQR(G223) + SQR(G313) + 
             SQR(G323)));
    
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
    
    CCTK_REAL const M1L  =  gu21*(-(G112*kxxL) + G111*kxyL - G212*kxyL - G312*kxzL + G211*kyyL + G311*kyzL - PDstandardNth1kxy + 
           PDstandardNth2kxx) + gu22*(-(G122*kxxL) + G112*kxyL - G222*kxyL - G322*kxzL + G212*kyyL + G312*kyzL - 
           PDstandardNth1kyy + PDstandardNth2kxy) + gu31*
         (-(G113*kxxL) - G213*kxyL + G111*kxzL - G313*kxzL + G211*kyzL + G311*kzzL - PDstandardNth1kxz + PDstandardNth3kxx)\
         + gu32*(G113*kxyL + G112*kxzL + G213*kyyL + (G212 + G313)*kyzL + G312*kzzL - 
           2*(G123*kxxL + G223*kxyL + G323*kxzL + PDstandardNth1kyz) + PDstandardNth2kxz + PDstandardNth3kxy) + 
        gu33*(-(G133*kxxL) - G233*kxyL + G113*kxzL - G333*kxzL + G213*kyzL + G313*kzzL - PDstandardNth1kzz + 
           PDstandardNth3kxz) - 25.13274122871834590770114706623602307358*S1;
    
    CCTK_REAL const M2L  =  gu11*(G112*kxxL + (-G111 + G212)*kxyL + G312*kxzL - G211*kyyL - G311*kyzL + PDstandardNth1kxy - 
           PDstandardNth2kxx) + gu21*(G122*kxxL + (-G112 + G222)*kxyL + G322*kxzL - G212*kyyL - G312*kyzL + 
           PDstandardNth1kyy - PDstandardNth2kxy) + gu31*
         (G123*kxxL + (-2*G113 + G223)*kxyL + (G112 + G323)*kxzL + G212*kyzL + G312*kzzL + PDstandardNth1kyz - 
           2*(G213*kyyL + G313*kyzL + PDstandardNth2kxz) + PDstandardNth3kxy) + 
        gu32*(-(G123*kxyL) + G122*kxzL - G223*kyyL + G222*kyzL - G323*kyzL + G322*kzzL - PDstandardNth2kyz + 
           PDstandardNth3kyy) + gu33*(-(G133*kxyL) + G123*kxzL - G233*kyyL + G223*kyzL - G333*kyzL + G323*kzzL - 
           PDstandardNth2kzz + PDstandardNth3kyz) - 25.13274122871834590770114706623602307358*S2;
    
    CCTK_REAL const M3L  =  gu11*(G113*kxxL + G213*kxyL + (-G111 + G313)*kxzL - G211*kyzL - G311*kzzL + PDstandardNth1kxz - 
           PDstandardNth3kxx) + gu21*(G123*kxxL + (G113 + G223)*kxyL + (-2*G112 + G323)*kxzL + G213*kyyL + 
           (-2*G212 + G313)*kyzL + PDstandardNth1kyz + PDstandardNth2kxz - 2*(G312*kzzL + PDstandardNth3kxy)) + 
        gu31*(G133*kxxL + G233*kxyL + (-G113 + G333)*kxzL - G213*kyzL - G313*kzzL + PDstandardNth1kzz - 
           PDstandardNth3kxz) + gu22*(G123*kxyL - G122*kxzL + G223*kyyL - G222*kyzL + G323*kyzL - G322*kzzL + 
           PDstandardNth2kyz - PDstandardNth3kyy) + gu32*
         (G133*kxyL - G123*kxzL + G233*kyyL - G223*kyzL + G333*kyzL - G323*kzzL + PDstandardNth2kzz - PDstandardNth3kyz) - 
        25.13274122871834590770114706623602307358*S3;
    
    
    /* Copy local copies back to grid functions */
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_ADMConstraints);
}

void ML_ADMConstraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADMConstraints_Body);
}
