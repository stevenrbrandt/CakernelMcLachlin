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

extern "C" void ML_ADMConstraints_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMConstraints::ML_Ham","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMConstraints::ML_Ham.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMConstraints::ML_mom","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMConstraints::ML_mom.");
  return;
}

static void ML_ADMConstraints_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
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
  CCTK_REAL const t = ToReal(cctk_time);
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
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_ADMConstraints,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL = alp[index];
    CCTK_REAL betaxL = betax[index];
    CCTK_REAL betayL = betay[index];
    CCTK_REAL betazL = betaz[index];
    CCTK_REAL gxxL = gxx[index];
    CCTK_REAL gxyL = gxy[index];
    CCTK_REAL gxzL = gxz[index];
    CCTK_REAL gyyL = gyy[index];
    CCTK_REAL gyzL = gyz[index];
    CCTK_REAL gzzL = gzz[index];
    CCTK_REAL kxxL = kxx[index];
    CCTK_REAL kxyL = kxy[index];
    CCTK_REAL kxzL = kxz[index];
    CCTK_REAL kyyL = kyy[index];
    CCTK_REAL kyzL = kyz[index];
    CCTK_REAL kzzL = kzz[index];
    
    CCTK_REAL eTttL, eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL;
    
    if (*stress_energy_state)
    {
      eTttL = eTtt[index];
      eTtxL = eTtx[index];
      eTtyL = eTty[index];
      eTtzL = eTtz[index];
      eTxxL = eTxx[index];
      eTxyL = eTxy[index];
      eTxzL = eTxz[index];
      eTyyL = eTyy[index];
      eTyzL = eTyz[index];
      eTzzL = eTzz[index];
    }
    else
    {
      eTttL = ToReal(0.0);
      eTtxL = ToReal(0.0);
      eTtyL = ToReal(0.0);
      eTtzL = ToReal(0.0);
      eTxxL = ToReal(0.0);
      eTxyL = ToReal(0.0);
      eTxzL = ToReal(0.0);
      eTyyL = ToReal(0.0);
      eTyzL = ToReal(0.0);
      eTzzL = ToReal(0.0);
    }
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDstandardNth1gxx = PDstandardNth1(&gxx[index]);
    CCTK_REAL const PDstandardNth2gxx = PDstandardNth2(&gxx[index]);
    CCTK_REAL const PDstandardNth3gxx = PDstandardNth3(&gxx[index]);
    CCTK_REAL const PDstandardNth22gxx = PDstandardNth22(&gxx[index]);
    CCTK_REAL const PDstandardNth33gxx = PDstandardNth33(&gxx[index]);
    CCTK_REAL const PDstandardNth23gxx = PDstandardNth23(&gxx[index]);
    CCTK_REAL const PDstandardNth1gxy = PDstandardNth1(&gxy[index]);
    CCTK_REAL const PDstandardNth2gxy = PDstandardNth2(&gxy[index]);
    CCTK_REAL const PDstandardNth3gxy = PDstandardNth3(&gxy[index]);
    CCTK_REAL const PDstandardNth33gxy = PDstandardNth33(&gxy[index]);
    CCTK_REAL const PDstandardNth12gxy = PDstandardNth12(&gxy[index]);
    CCTK_REAL const PDstandardNth13gxy = PDstandardNth13(&gxy[index]);
    CCTK_REAL const PDstandardNth23gxy = PDstandardNth23(&gxy[index]);
    CCTK_REAL const PDstandardNth1gxz = PDstandardNth1(&gxz[index]);
    CCTK_REAL const PDstandardNth2gxz = PDstandardNth2(&gxz[index]);
    CCTK_REAL const PDstandardNth3gxz = PDstandardNth3(&gxz[index]);
    CCTK_REAL const PDstandardNth22gxz = PDstandardNth22(&gxz[index]);
    CCTK_REAL const PDstandardNth12gxz = PDstandardNth12(&gxz[index]);
    CCTK_REAL const PDstandardNth13gxz = PDstandardNth13(&gxz[index]);
    CCTK_REAL const PDstandardNth23gxz = PDstandardNth23(&gxz[index]);
    CCTK_REAL const PDstandardNth1gyy = PDstandardNth1(&gyy[index]);
    CCTK_REAL const PDstandardNth2gyy = PDstandardNth2(&gyy[index]);
    CCTK_REAL const PDstandardNth3gyy = PDstandardNth3(&gyy[index]);
    CCTK_REAL const PDstandardNth11gyy = PDstandardNth11(&gyy[index]);
    CCTK_REAL const PDstandardNth33gyy = PDstandardNth33(&gyy[index]);
    CCTK_REAL const PDstandardNth13gyy = PDstandardNth13(&gyy[index]);
    CCTK_REAL const PDstandardNth1gyz = PDstandardNth1(&gyz[index]);
    CCTK_REAL const PDstandardNth2gyz = PDstandardNth2(&gyz[index]);
    CCTK_REAL const PDstandardNth3gyz = PDstandardNth3(&gyz[index]);
    CCTK_REAL const PDstandardNth11gyz = PDstandardNth11(&gyz[index]);
    CCTK_REAL const PDstandardNth12gyz = PDstandardNth12(&gyz[index]);
    CCTK_REAL const PDstandardNth13gyz = PDstandardNth13(&gyz[index]);
    CCTK_REAL const PDstandardNth23gyz = PDstandardNth23(&gyz[index]);
    CCTK_REAL const PDstandardNth1gzz = PDstandardNth1(&gzz[index]);
    CCTK_REAL const PDstandardNth2gzz = PDstandardNth2(&gzz[index]);
    CCTK_REAL const PDstandardNth3gzz = PDstandardNth3(&gzz[index]);
    CCTK_REAL const PDstandardNth11gzz = PDstandardNth11(&gzz[index]);
    CCTK_REAL const PDstandardNth22gzz = PDstandardNth22(&gzz[index]);
    CCTK_REAL const PDstandardNth12gzz = PDstandardNth12(&gzz[index]);
    CCTK_REAL const PDstandardNth2kxx = PDstandardNth2(&kxx[index]);
    CCTK_REAL const PDstandardNth3kxx = PDstandardNth3(&kxx[index]);
    CCTK_REAL const PDstandardNth1kxy = PDstandardNth1(&kxy[index]);
    CCTK_REAL const PDstandardNth2kxy = PDstandardNth2(&kxy[index]);
    CCTK_REAL const PDstandardNth3kxy = PDstandardNth3(&kxy[index]);
    CCTK_REAL const PDstandardNth1kxz = PDstandardNth1(&kxz[index]);
    CCTK_REAL const PDstandardNth2kxz = PDstandardNth2(&kxz[index]);
    CCTK_REAL const PDstandardNth3kxz = PDstandardNth3(&kxz[index]);
    CCTK_REAL const PDstandardNth1kyy = PDstandardNth1(&kyy[index]);
    CCTK_REAL const PDstandardNth3kyy = PDstandardNth3(&kyy[index]);
    CCTK_REAL const PDstandardNth1kyz = PDstandardNth1(&kyz[index]);
    CCTK_REAL const PDstandardNth2kyz = PDstandardNth2(&kyz[index]);
    CCTK_REAL const PDstandardNth3kyz = PDstandardNth3(&kyz[index]);
    CCTK_REAL const PDstandardNth1kzz = PDstandardNth1(&kzz[index]);
    CCTK_REAL const PDstandardNth2kzz = PDstandardNth2(&kzz[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL detg = 2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - 
      SQR(gxyL)) - gyyL*SQR(gxzL) - gxxL*SQR(gyzL);
    
    CCTK_REAL gu11 = INV(detg)*(gyyL*gzzL - SQR(gyzL));
    
    CCTK_REAL gu21 = (gxzL*gyzL - gxyL*gzzL)*INV(detg);
    
    CCTK_REAL gu31 = (-(gxzL*gyyL) + gxyL*gyzL)*INV(detg);
    
    CCTK_REAL gu22 = INV(detg)*(gxxL*gzzL - SQR(gxzL));
    
    CCTK_REAL gu32 = (gxyL*gxzL - gxxL*gyzL)*INV(detg);
    
    CCTK_REAL gu33 = INV(detg)*(gxxL*gyyL - SQR(gxyL));
    
    CCTK_REAL G111 = 0.5*(gu11*PDstandardNth1gxx + 
      2*(gu21*PDstandardNth1gxy + gu31*PDstandardNth1gxz) - 
      gu21*PDstandardNth2gxx - gu31*PDstandardNth3gxx);
    
    CCTK_REAL G211 = 0.5*(gu21*PDstandardNth1gxx + 
      2*(gu22*PDstandardNth1gxy + gu32*PDstandardNth1gxz) - 
      gu22*PDstandardNth2gxx - gu32*PDstandardNth3gxx);
    
    CCTK_REAL G311 = 0.5*(gu31*PDstandardNth1gxx + 
      2*(gu32*PDstandardNth1gxy + gu33*PDstandardNth1gxz) - 
      gu32*PDstandardNth2gxx - gu33*PDstandardNth3gxx);
    
    CCTK_REAL G112 = 0.5*(gu21*PDstandardNth1gyy + gu11*PDstandardNth2gxx 
      + gu31*(PDstandardNth1gyz + PDstandardNth2gxz - PDstandardNth3gxy));
    
    CCTK_REAL G212 = 0.5*(gu22*PDstandardNth1gyy + gu21*PDstandardNth2gxx 
      + gu32*(PDstandardNth1gyz + PDstandardNth2gxz - PDstandardNth3gxy));
    
    CCTK_REAL G312 = 0.5*(gu32*PDstandardNth1gyy + gu31*PDstandardNth2gxx 
      + gu33*(PDstandardNth1gyz + PDstandardNth2gxz - PDstandardNth3gxy));
    
    CCTK_REAL G113 = 0.5*(gu31*PDstandardNth1gzz + gu11*PDstandardNth3gxx 
      + gu21*(PDstandardNth1gyz - PDstandardNth2gxz + PDstandardNth3gxy));
    
    CCTK_REAL G213 = 0.5*(gu32*PDstandardNth1gzz + gu21*PDstandardNth3gxx 
      + gu22*(PDstandardNth1gyz - PDstandardNth2gxz + PDstandardNth3gxy));
    
    CCTK_REAL G313 = 0.5*(gu33*PDstandardNth1gzz + gu31*PDstandardNth3gxx 
      + gu32*(PDstandardNth1gyz - PDstandardNth2gxz + PDstandardNth3gxy));
    
    CCTK_REAL G122 = 0.5*(gu11*(-PDstandardNth1gyy + 2*PDstandardNth2gxy) 
      + gu21*PDstandardNth2gyy + gu31*(2*PDstandardNth2gyz - 
      PDstandardNth3gyy));
    
    CCTK_REAL G222 = 0.5*(gu21*(-PDstandardNth1gyy + 2*PDstandardNth2gxy) 
      + gu22*PDstandardNth2gyy + gu32*(2*PDstandardNth2gyz - 
      PDstandardNth3gyy));
    
    CCTK_REAL G322 = 0.5*(gu31*(-PDstandardNth1gyy + 2*PDstandardNth2gxy) 
      + gu32*PDstandardNth2gyy + gu33*(2*PDstandardNth2gyz - 
      PDstandardNth3gyy));
    
    CCTK_REAL G123 = 0.5*(gu31*PDstandardNth2gzz + 
      gu11*(-PDstandardNth1gyz + PDstandardNth2gxz + PDstandardNth3gxy) + 
      gu21*PDstandardNth3gyy);
    
    CCTK_REAL G223 = 0.5*(gu32*PDstandardNth2gzz + 
      gu21*(-PDstandardNth1gyz + PDstandardNth2gxz + PDstandardNth3gxy) + 
      gu22*PDstandardNth3gyy);
    
    CCTK_REAL G323 = 0.5*(gu33*PDstandardNth2gzz + 
      gu31*(-PDstandardNth1gyz + PDstandardNth2gxz + PDstandardNth3gxy) + 
      gu32*PDstandardNth3gyy);
    
    CCTK_REAL G133 = 0.5*(gu11*(-PDstandardNth1gzz + 2*PDstandardNth3gxz) 
      + gu21*(-PDstandardNth2gzz + 2*PDstandardNth3gyz) + 
      gu31*PDstandardNth3gzz);
    
    CCTK_REAL G233 = 0.5*(gu21*(-PDstandardNth1gzz + 2*PDstandardNth3gxz) 
      + gu22*(-PDstandardNth2gzz + 2*PDstandardNth3gyz) + 
      gu32*PDstandardNth3gzz);
    
    CCTK_REAL G333 = 0.5*(gu31*(-PDstandardNth1gzz + 2*PDstandardNth3gxz) 
      + gu32*(-PDstandardNth2gzz + 2*PDstandardNth3gyz) + 
      gu33*PDstandardNth3gzz);
    
    CCTK_REAL R11 = 0.5*(gu32*(4*(gxyL*(-(G123*G211) + G113*G212) + 
      gyyL*(G212*G213 - G211*G223) + gxzL*G113*G312 + gyzL*G212*G313 + 
      gzzL*G312*G313 + G112*(gxxL*G113 + gxyL*G213 + gxzL*G313) - 
      G111*(gxxL*G123 + gxyL*G223 + gxzL*G323) + G311*(-(gxzL*G123) - 
      gyzL*G223 - gzzL*G323) + gyzL*(G213*G312 - G211*G323)) + 
      2*(-PDstandardNth11gyz + PDstandardNth12gxz + PDstandardNth13gxy - 
      PDstandardNth23gxx)) + gu22*(-PDstandardNth11gyy + 2*PDstandardNth12gxy 
      - PDstandardNth22gxx + 2*(2*gyzL*G212*G312 + 2*G112*(gxyL*G212 + 
      gxzL*G312) - G111*(gxxL*G122 + gxyL*G222 + gxzL*G322) + 
      G211*(-(gxyL*G122) - gyyL*G222 - gyzL*G322) + 
      G311*(-(gxzL*G122) - gyzL*G222 - gzzL*G322) + gxxL*SQR(G112) + 
      gyyL*SQR(G212) + gzzL*SQR(G312))) + gu33*(-PDstandardNth11gzz + 
      2*PDstandardNth13gxz - PDstandardNth33gxx + 2*(2*gyzL*G213*G313 + 
      2*G113*(gxyL*G213 + gxzL*G313) - G111*(gxxL*G133 + gxyL*G233 + 
      gxzL*G333) + G211*(-(gxyL*G133) - gyyL*G233 - gyzL*G333) + 
      G311*(-(gxzL*G133) - gyzL*G233 - gzzL*G333) + gxxL*SQR(G113) + 
      gyyL*SQR(G213) + gzzL*SQR(G313))));
    
    CCTK_REAL R12 = 0.5*(gu21*(PDstandardNth11gyy - 2*PDstandardNth12gxy + 
      PDstandardNth22gxx) + gu31*(PDstandardNth11gyz - PDstandardNth12gxz - 
      PDstandardNth13gxy + PDstandardNth23gxx) + gu32*(-PDstandardNth12gyz + 
      PDstandardNth13gyy + PDstandardNth22gxz - PDstandardNth23gxy) + 
      gu33*(-PDstandardNth12gzz + PDstandardNth13gyz + PDstandardNth23gxz - 
      PDstandardNth33gxy) + 2*((gxyL*(G123*G211 - G113*G212) + 
      gyyL*(-(G212*G213) + G211*G223) - G112*(gxxL*G113 + gxyL*G213 + 
      gxzL*G313) + G312*(-(gxzL*G113) - gyzL*G213 - gzzL*G313) + 
      G111*(gxxL*G123 + gxyL*G223 + gxzL*G323) + G311*(gxzL*G123 + 
      gyzL*G223 + gzzL*G323) + gyzL*(-(G212*G313) + G211*G323))*gu31 + 
      (gxyL*(-(G123*G212) + G122*G213) + gyyL*(G213*G222 - G212*G223) + 
      gxzL*G122*G313 + gyzL*G213*G322 + gzzL*G313*G322 + 
      G113*(gxxL*G122 + gxyL*G222 + gxzL*G322) - G112*(gxxL*G123 + 
      gxyL*G223 + gxzL*G323) + G312*(-(gxzL*G123) - gyzL*G223 - 
      gzzL*G323) + gyzL*(G222*G313 - G212*G323))*gu32 + 
      (gxyL*(-(G133*G212) + G123*G213) + gyyL*(G213*G223 - G212*G233) + 
      gxzL*G123*G313 + gyzL*G213*G323 + gzzL*G313*G323 + 
      G113*(gxxL*G123 + gxyL*G223 + gxzL*G323) - G112*(gxxL*G133 + 
      gxyL*G233 + gxzL*G333) + G312*(-(gxzL*G133) - gyzL*G233 - 
      gzzL*G333) + gyzL*(G223*G313 - G212*G333))*gu33 + 
      gu21*(-2*(gyzL*G212*G312 + G112*(gxyL*G212 + gxzL*G312)) + 
      G111*(gxxL*G122 + gxyL*G222 + gxzL*G322) + G211*(gxyL*G122 + 
      gyyL*G222 + gyzL*G322) + G311*(gxzL*G122 + gyzL*G222 + 
      gzzL*G322) - gxxL*SQR(G112) - gyyL*SQR(G212) - 
      gzzL*SQR(G312))));
    
    CCTK_REAL R13 = 0.5*(gu21*(PDstandardNth11gyz - PDstandardNth12gxz - 
      PDstandardNth13gxy + PDstandardNth23gxx) + gu22*(PDstandardNth12gyz - 
      PDstandardNth13gyy - PDstandardNth22gxz + PDstandardNth23gxy) + 
      gu31*(PDstandardNth11gzz - 2*PDstandardNth13gxz + PDstandardNth33gxx) + 
      gu32*(PDstandardNth12gzz - PDstandardNth13gyz - PDstandardNth23gxz + 
      PDstandardNth33gxy) + 2*((gxyL*(G123*G211 - G113*G212) + 
      gyyL*(-(G212*G213) + G211*G223) - G112*(gxxL*G113 + gxyL*G213 + 
      gxzL*G313) + G312*(-(gxzL*G113) - gyzL*G213 - gzzL*G313) + 
      G111*(gxxL*G123 + gxyL*G223 + gxzL*G323) + G311*(gxzL*G123 + 
      gyzL*G223 + gzzL*G323) + gyzL*(-(G212*G313) + G211*G323))*gu21 + 
      (gxyL*(G123*G212 - G122*G213) + gyyL*(-(G213*G222) + G212*G223) - 
      G113*(gxxL*G122 + gxyL*G222 + gxzL*G322) + G313*(-(gxzL*G122) - 
      gyzL*G222 - gzzL*G322) + G112*(gxxL*G123 + gxyL*G223 + 
      gxzL*G323) + G312*(gxzL*G123 + gyzL*G223 + gzzL*G323) + 
      gyzL*(-(G213*G322) + G212*G323))*gu22 + (gxyL*(G133*G212 - 
      G123*G213) + gyyL*(-(G213*G223) + G212*G233) - G113*(gxxL*G123 + 
      gxyL*G223 + gxzL*G323) + G313*(-(gxzL*G123) - gyzL*G223 - 
      gzzL*G323) + G112*(gxxL*G133 + gxyL*G233 + gxzL*G333) + 
      G312*(gxzL*G133 + gyzL*G233 + gzzL*G333) + gyzL*(-(G213*G323) + 
      G212*G333))*gu32 + gu31*(-2*(gyzL*G213*G313 + G113*(gxyL*G213 + 
      gxzL*G313)) + G111*(gxxL*G133 + gxyL*G233 + gxzL*G333) + 
      G211*(gxyL*G133 + gyyL*G233 + gyzL*G333) + G311*(gxzL*G133 + 
      gyzL*G233 + gzzL*G333) - gxxL*SQR(G113) - gyyL*SQR(G213) - 
      gzzL*SQR(G313))));
    
    CCTK_REAL R22 = 0.5*(gu31*(4*(gxyL*(G123*G212 - G122*G213) + 
      gyyL*(-(G213*G222) + G212*G223) + gxzL*(G123*G312 - G122*G313) - 
      G113*(gxxL*G122 + gxyL*G222 + gxzL*G322) + gyzL*(G223*G312 - 
      G222*G313 - G213*G322) + gyzL*G212*G323 + G112*(gxxL*G123 + 
      gxyL*G223 + gxzL*G323) + gzzL*(-(G313*G322) + G312*G323)) + 
      2*(PDstandardNth12gyz - PDstandardNth13gyy - PDstandardNth22gxz + 
      PDstandardNth23gxy)) + gu11*(-PDstandardNth11gyy + 2*PDstandardNth12gxy 
      - PDstandardNth22gxx + 2*(2*gyzL*G212*G312 + 2*G112*(gxyL*G212 + 
      gxzL*G312) - G111*(gxxL*G122 + gxyL*G222 + gxzL*G322) + 
      G211*(-(gxyL*G122) - gyyL*G222 - gyzL*G322) + 
      G311*(-(gxzL*G122) - gyzL*G222 - gzzL*G322) + gxxL*SQR(G112) + 
      gyyL*SQR(G212) + gzzL*SQR(G312))) + gu33*(-PDstandardNth22gzz + 
      2*PDstandardNth23gyz - PDstandardNth33gyy + 2*(2*gyzL*G223*G323 + 
      2*G123*(gxyL*G223 + gxzL*G323) - G122*(gxxL*G133 + gxyL*G233 + 
      gxzL*G333) + G222*(-(gxyL*G133) - gyyL*G233 - gyzL*G333) + 
      G322*(-(gxzL*G133) - gyzL*G233 - gzzL*G333) + gxxL*SQR(G123) + 
      gyyL*SQR(G223) + gzzL*SQR(G323))));
    
    CCTK_REAL R23 = 0.5*(gu11*(-PDstandardNth11gyz + PDstandardNth12gxz + 
      PDstandardNth13gxy - PDstandardNth23gxx) + gu21*(-PDstandardNth12gyz + 
      PDstandardNth13gyy + PDstandardNth22gxz - PDstandardNth23gxy) + 
      gu31*(PDstandardNth12gzz - PDstandardNth13gyz - PDstandardNth23gxz + 
      PDstandardNth33gxy) + gu32*(PDstandardNth22gzz - 2*PDstandardNth23gyz + 
      PDstandardNth33gyy) + 2*((gxyL*(-(G123*G211) + G113*G212) + 
      gyyL*(G212*G213 - G211*G223) + gxzL*G113*G312 + gyzL*G212*G313 + 
      gzzL*G312*G313 + G112*(gxxL*G113 + gxyL*G213 + gxzL*G313) - 
      G111*(gxxL*G123 + gxyL*G223 + gxzL*G323) + G311*(-(gxzL*G123) - 
      gyzL*G223 - gzzL*G323) + gyzL*(G213*G312 - G211*G323))*gu11 + 
      (gxyL*(-(G123*G212) + G122*G213) + gyyL*(G213*G222 - G212*G223) + 
      gxzL*G122*G313 + gyzL*G213*G322 + gzzL*G313*G322 + 
      G113*(gxxL*G122 + gxyL*G222 + gxzL*G322) - G112*(gxxL*G123 + 
      gxyL*G223 + gxzL*G323) + G312*(-(gxzL*G123) - gyzL*G223 - 
      gzzL*G323) + gyzL*(G222*G313 - G212*G323))*gu21 + 
      (gxyL*(G133*G212 - G123*G213) + gyyL*(-(G213*G223) + G212*G233) - 
      G113*(gxxL*G123 + gxyL*G223 + gxzL*G323) + G313*(-(gxzL*G123) - 
      gyzL*G223 - gzzL*G323) + G112*(gxxL*G133 + gxyL*G233 + 
      gxzL*G333) + G312*(gxzL*G133 + gyzL*G233 + gzzL*G333) + 
      gyzL*(-(G213*G323) + G212*G333))*gu31 + gu32*(-2*(gyzL*G223*G323 + 
      G123*(gxyL*G223 + gxzL*G323)) + G122*(gxxL*G133 + gxyL*G233 + 
      gxzL*G333) + G222*(gxyL*G133 + gyyL*G233 + gyzL*G333) + 
      G322*(gxzL*G133 + gyzL*G233 + gzzL*G333) - gxxL*SQR(G123) - 
      gyyL*SQR(G223) - gzzL*SQR(G323))));
    
    CCTK_REAL R33 = 0.5*(gu21*(4*(gxyL*(-(G133*G212) + G123*G213) + 
      gyyL*(G213*G223 - G212*G233) + gxzL*G123*G313 + gyzL*G213*G323 + 
      gzzL*G313*G323 + G113*(gxxL*G123 + gxyL*G223 + gxzL*G323) - 
      G112*(gxxL*G133 + gxyL*G233 + gxzL*G333) + G312*(-(gxzL*G133) - 
      gyzL*G233 - gzzL*G333) + gyzL*(G223*G313 - G212*G333)) + 
      2*(-PDstandardNth12gzz + PDstandardNth13gyz + PDstandardNth23gxz - 
      PDstandardNth33gxy)) + gu11*(-PDstandardNth11gzz + 2*PDstandardNth13gxz 
      - PDstandardNth33gxx + 2*(2*gyzL*G213*G313 + 2*G113*(gxyL*G213 + 
      gxzL*G313) - G111*(gxxL*G133 + gxyL*G233 + gxzL*G333) + 
      G211*(-(gxyL*G133) - gyyL*G233 - gyzL*G333) + 
      G311*(-(gxzL*G133) - gyzL*G233 - gzzL*G333) + gxxL*SQR(G113) + 
      gyyL*SQR(G213) + gzzL*SQR(G313))) + gu22*(-PDstandardNth22gzz + 
      2*PDstandardNth23gyz - PDstandardNth33gyy + 2*(2*gyzL*G223*G323 + 
      2*G123*(gxyL*G223 + gxzL*G323) - G122*(gxxL*G133 + gxyL*G233 + 
      gxzL*G333) + G222*(-(gxyL*G133) - gyyL*G233 - gyzL*G333) + 
      G322*(-(gxzL*G133) - gyzL*G233 - gzzL*G333) + gxxL*SQR(G123) + 
      gyyL*SQR(G223) + gzzL*SQR(G323))));
    
    CCTK_REAL trR = gu11*R11 + gu22*R22 + 2*(gu21*R12 + gu31*R13 + 
      gu32*R23) + gu33*R33;
    
    CCTK_REAL Km11 = kxxL*gu11 + kxyL*gu21 + kxzL*gu31;
    
    CCTK_REAL Km21 = kxxL*gu21 + kxyL*gu22 + kxzL*gu32;
    
    CCTK_REAL Km31 = kxxL*gu31 + kxyL*gu32 + kxzL*gu33;
    
    CCTK_REAL Km12 = kxyL*gu11 + kyyL*gu21 + kyzL*gu31;
    
    CCTK_REAL Km22 = kxyL*gu21 + kyyL*gu22 + kyzL*gu32;
    
    CCTK_REAL Km32 = kxyL*gu31 + kyyL*gu32 + kyzL*gu33;
    
    CCTK_REAL Km13 = kxzL*gu11 + kyzL*gu21 + kzzL*gu31;
    
    CCTK_REAL Km23 = kxzL*gu21 + kyzL*gu22 + kzzL*gu32;
    
    CCTK_REAL Km33 = kxzL*gu31 + kyzL*gu32 + kzzL*gu33;
    
    CCTK_REAL trK = Km11 + Km22 + Km33;
    
    CCTK_REAL rho = INV(SQR(alpL))*(eTttL - 2*(betayL*eTtyL + 
      betazL*eTtzL) + 2*(betaxL*(-eTtxL + betayL*eTxyL + 
      betazL*eTxzL) + betayL*betazL*eTyzL) + eTxxL*SQR(betaxL) 
      + eTyyL*SQR(betayL) + eTzzL*SQR(betazL));
    
    CCTK_REAL S1 = (-eTtxL + betaxL*eTxxL + betayL*eTxyL + 
      betazL*eTxzL)*INV(alpL);
    
    CCTK_REAL S2 = (-eTtyL + betaxL*eTxyL + betayL*eTyyL + 
      betazL*eTyzL)*INV(alpL);
    
    CCTK_REAL S3 = (-eTtzL + betaxL*eTxzL + betayL*eTyzL + 
      betazL*eTzzL)*INV(alpL);
    
    CCTK_REAL HL = -2*(Km12*Km21 + Km13*Km31 + Km23*Km32) - 16*Pi*rho + 
      trR - SQR(Km11) - SQR(Km22) - SQR(Km33) + SQR(trK);
    
    CCTK_REAL M1L = gu21*(-(kxxL*G112) + kyyL*G211 + kxyL*(G111 - 
      G212) + kyzL*G311 - kxzL*G312 - PDstandardNth1kxy + 
      PDstandardNth2kxx) + gu22*(-(kxxL*G122) + kyyL*G212 + kxyL*(G112 
      - G222) + kyzL*G312 - kxzL*G322 - PDstandardNth1kyy + 
      PDstandardNth2kxy) + gu31*(-(kxxL*G113) + kyzL*G211 - kxyL*G213 + 
      kzzL*G311 + kxzL*(G111 - G313) - PDstandardNth1kxz + 
      PDstandardNth3kxx) + gu32*(kyyL*G213 + kxyL*(G113 - 2*G223) + 
      kzzL*G312 + kyzL*(G212 + G313) + kxzL*(G112 - 2*G323) - 
      2*(kxxL*G123 + PDstandardNth1kyz) + PDstandardNth2kxz + 
      PDstandardNth3kxy) + gu33*(-(kxxL*G133) + kyzL*G213 - kxyL*G233 + 
      kzzL*G313 + kxzL*(G113 - G333) - PDstandardNth1kzz + 
      PDstandardNth3kxz) - 8*Pi*S1;
    
    CCTK_REAL M2L = gu11*(kxxL*G112 - kyyL*G211 + kxyL*(-G111 + 
      G212) - kyzL*G311 + kxzL*G312 + PDstandardNth1kxy - 
      PDstandardNth2kxx) + gu21*(kxxL*G122 - kyyL*G212 + kxyL*(-G112 + 
      G222) - kyzL*G312 + kxzL*G322 + PDstandardNth1kyy - 
      PDstandardNth2kxy) + gu31*(kxxL*G123 + kxyL*G223 + kzzL*G312 + 
      kyzL*(G212 - 2*G313) + kxzL*(G112 + G323) + PDstandardNth1kyz - 
      2*(kxyL*G113 + kyyL*G213 + PDstandardNth2kxz) + PDstandardNth3kxy) 
      + gu32*(kxzL*G122 - kxyL*G123 - kyyL*G223 + kzzL*G322 + 
      kyzL*(G222 - G323) - PDstandardNth2kyz + PDstandardNth3kyy) + 
      gu33*(kxzL*G123 - kxyL*G133 - kyyL*G233 + kzzL*G323 + 
      kyzL*(G223 - G333) - PDstandardNth2kzz + PDstandardNth3kyz) - 
      8*Pi*S2;
    
    CCTK_REAL M3L = gu11*(kxxL*G113 - kyzL*G211 + kxyL*G213 - 
      kzzL*G311 + kxzL*(-G111 + G313) + PDstandardNth1kxz - 
      PDstandardNth3kxx) + gu21*(kxxL*G123 + kyyL*G213 + kxyL*(G113 + 
      G223) + kyzL*G313 + kxzL*G323 + PDstandardNth1kyz + 
      PDstandardNth2kxz - 2*(kxzL*G112 + kyzL*G212 + kzzL*G312 + 
      PDstandardNth3kxy)) + gu31*(kxxL*G133 - kyzL*G213 + kxyL*G233 - 
      kzzL*G313 + kxzL*(-G113 + G333) + PDstandardNth1kzz - 
      PDstandardNth3kxz) + gu22*(-(kxzL*G122) + kxyL*G123 + kyyL*G223 - 
      kzzL*G322 + kyzL*(-G222 + G323) + PDstandardNth2kyz - 
      PDstandardNth3kyy) + gu32*(-(kxzL*G123) + kxyL*G133 + kyyL*G233 - 
      kzzL*G323 + kyzL*(-G223 + G333) + PDstandardNth2kzz - 
      PDstandardNth3kyz) - 8*Pi*S3;
    
    /* Copy local copies back to grid functions */
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
  }
  LC_ENDLOOP3 (ML_ADMConstraints);
}

extern "C" void ML_ADMConstraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMConstraints_Body");
  }
  
  if (cctk_iteration % ML_ADMConstraints_calc_every != ML_ADMConstraints_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::curv","ADMBase::lapse","ADMBase::metric","ADMBase::shift","ML_ADMConstraints::ML_Ham","ML_ADMConstraints::ML_mom"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADMConstraints", 6, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_ADMConstraints", 2, 2, 2);
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADMConstraints_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_ADMConstraints_Body");
  }
}
