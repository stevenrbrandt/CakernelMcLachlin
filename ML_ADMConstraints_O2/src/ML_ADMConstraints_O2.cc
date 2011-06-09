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

extern "C" void ML_ADMConstraints_O2_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMConstraints_O2::ML_Ham","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMConstraints_O2::ML_Ham.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMConstraints_O2::ML_mom","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMConstraints_O2::ML_mom.");
  return;
}

static void ML_ADMConstraints_O2_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMConstraints_O2_Body");
  }
  
  if (cctk_iteration % ML_ADMConstraints_O2_calc_every != ML_ADMConstraints_O2_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::curv","ADMBase::lapse","ADMBase::metric","ADMBase::shift","ML_ADMConstraints_O2::ML_Ham","ML_ADMConstraints_O2::ML_mom"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADMConstraints_O2", 6, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_ADMConstraints_O2", 1, 1, 1);
  
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
  CCTK_REAL const p1o2dx = 0.5*INV(dx);
  CCTK_REAL const p1o2dy = 0.5*INV(dy);
  CCTK_REAL const p1o2dz = 0.5*INV(dz);
  CCTK_REAL const p1o4dxdy = 0.25*INV(dx)*INV(dy);
  CCTK_REAL const p1o4dxdz = 0.25*INV(dx)*INV(dz);
  CCTK_REAL const p1o4dydz = 0.25*INV(dy)*INV(dz);
  CCTK_REAL const p1odx2 = INV(SQR(dx));
  CCTK_REAL const p1ody2 = INV(SQR(dy));
  CCTK_REAL const p1odz2 = INV(SQR(dz));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_ADMConstraints_O2,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
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
    CCTK_REAL detg = 2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - SQR(gxyL)) - 
      gyyL*SQR(gxzL) - gxxL*SQR(gyzL);
    
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
    
    CCTK_REAL R11 = 0.5*(-2*(G111*G122 + G111*G133 + G211*G222 + G211*G233 
      + G311*G322 + G311*G333 + gu32*PDstandardNth11gyz) + 
      2*gu32*PDstandardNth12gxz + 2*gu32*PDstandardNth13gxy + 
      gu22*(-PDstandardNth11gyy + 2*PDstandardNth12gxy - PDstandardNth22gxx) 
      - 2*gu32*PDstandardNth23gxx + gu33*(-PDstandardNth11gzz + 
      2*PDstandardNth13gxz - PDstandardNth33gxx) + 2*SQR(G112) + 2*SQR(G113) 
      + 2*SQR(G212) + 2*SQR(G213) + 2*SQR(G312) + 2*SQR(G313));
    
    CCTK_REAL R12 = 0.5*(2*(G113*G123 + G213*G223 + G313*G323) - 
      2*(G112*G133 + G212*G233 + G312*G333 + gu21*PDstandardNth12gxy) + 
      gu31*(PDstandardNth11gyz - PDstandardNth12gxz - PDstandardNth13gxy) + 
      gu21*(PDstandardNth11gyy + PDstandardNth22gxx) + 
      gu32*PDstandardNth22gxz + gu31*PDstandardNth23gxx + 
      gu32*(-PDstandardNth12gyz + PDstandardNth13gyy - PDstandardNth23gxy) + 
      gu33*PDstandardNth23gxz + gu33*(-PDstandardNth12gzz + 
      PDstandardNth13gyz - PDstandardNth33gxy));
    
    CCTK_REAL R13 = 0.5*(2*(G112*G123 + G212*G223 + G312*G323) - 
      2*(G113*G122 + G213*G222 + G313*G322 + gu31*PDstandardNth13gxz) + 
      gu21*(PDstandardNth11gyz - PDstandardNth12gxz - PDstandardNth13gxy + 
      PDstandardNth23gxx) + gu22*(PDstandardNth12gyz - PDstandardNth13gyy - 
      PDstandardNth22gxz + PDstandardNth23gxy) + gu31*(PDstandardNth11gzz + 
      PDstandardNth33gxx) + gu32*(PDstandardNth12gzz - PDstandardNth13gyz - 
      PDstandardNth23gxz + PDstandardNth33gxy));
    
    CCTK_REAL R22 = 0.5*(-2*(G122*(G111 + G133) + G222*(G211 + G233) + 
      G322*(G311 + G333) + gu31*PDstandardNth13gyy) + 
      gu11*(-PDstandardNth11gyy + 2*PDstandardNth12gxy - PDstandardNth22gxx) 
      + gu31*(-2*PDstandardNth22gxz + 2*(PDstandardNth12gyz + 
      PDstandardNth23gxy)) + gu33*(-PDstandardNth22gzz + 2*PDstandardNth23gyz 
      - PDstandardNth33gyy) + 2*(SQR(G112) + SQR(G123) + SQR(G212) + 
      SQR(G223) + SQR(G312) + SQR(G323)));
    
    CCTK_REAL R23 = 0.5*(2*(G112*G113 + G212*G213 + G312*G313) + 
      gu11*(-PDstandardNth11gyz + PDstandardNth12gxz + PDstandardNth13gxy - 
      PDstandardNth23gxx) + gu21*(-PDstandardNth12gyz + PDstandardNth13gyy + 
      PDstandardNth22gxz - PDstandardNth23gxy) - 2*(G111*G123 + G211*G223 + 
      G311*G323 + gu32*PDstandardNth23gyz) + gu31*(PDstandardNth12gzz - 
      PDstandardNth13gyz - PDstandardNth23gxz + PDstandardNth33gxy) + 
      gu32*(PDstandardNth22gzz + PDstandardNth33gyy));
    
    CCTK_REAL R33 = 0.5*(gu11*(-PDstandardNth11gzz + 2*PDstandardNth13gxz 
      - PDstandardNth33gxx) - 2*((G111 + G122)*G133 + (G211 + G222)*G233 + 
      (G311 + G322)*G333 + gu21*(PDstandardNth12gzz + PDstandardNth33gxy)) + 
      gu22*(-PDstandardNth22gzz + 2*PDstandardNth23gyz - PDstandardNth33gyy) 
      + 2*(gu21*(PDstandardNth13gyz + PDstandardNth23gxz) + SQR(G113) + 
      SQR(G123) + SQR(G213) + SQR(G223) + SQR(G313) + SQR(G323)));
    
    CCTK_REAL trR = gu11*R11 + gu22*R22 + 2*(gu21*R12 + gu31*R13 + 
      gu32*R23) + gu33*R33;
    
    CCTK_REAL Km11 = gu11*kxxL + gu21*kxyL + gu31*kxzL;
    
    CCTK_REAL Km21 = gu21*kxxL + gu22*kxyL + gu32*kxzL;
    
    CCTK_REAL Km31 = gu31*kxxL + gu32*kxyL + gu33*kxzL;
    
    CCTK_REAL Km12 = gu11*kxyL + gu21*kyyL + gu31*kyzL;
    
    CCTK_REAL Km22 = gu21*kxyL + gu22*kyyL + gu32*kyzL;
    
    CCTK_REAL Km32 = gu31*kxyL + gu32*kyyL + gu33*kyzL;
    
    CCTK_REAL Km13 = gu11*kxzL + gu21*kyzL + gu31*kzzL;
    
    CCTK_REAL Km23 = gu21*kxzL + gu22*kyzL + gu32*kzzL;
    
    CCTK_REAL Km33 = gu31*kxzL + gu32*kyzL + gu33*kzzL;
    
    CCTK_REAL trK = Km11 + Km22 + Km33;
    
    CCTK_REAL rho = INV(SQR(alpL))*(eTttL - 2*(betayL*eTtyL + 
      betazL*eTtzL) + 2*(betaxL*(-eTtxL + betayL*eTxyL + betazL*eTxzL) + 
      betayL*betazL*eTyzL) + eTxxL*SQR(betaxL) + eTyyL*SQR(betayL) + 
      eTzzL*SQR(betazL));
    
    CCTK_REAL S1 = (-eTtxL + betaxL*eTxxL + betayL*eTxyL + 
      betazL*eTxzL)*INV(alpL);
    
    CCTK_REAL S2 = (-eTtyL + betaxL*eTxyL + betayL*eTyyL + 
      betazL*eTyzL)*INV(alpL);
    
    CCTK_REAL S3 = (-eTtzL + betaxL*eTxzL + betayL*eTyzL + 
      betazL*eTzzL)*INV(alpL);
    
    CCTK_REAL HL = -2*(Km12*Km21 + Km13*Km31 + Km23*Km32) - 
      50.26548245743669181540229413247204614715*rho + trR - SQR(Km11) - 
      SQR(Km22) - SQR(Km33) + SQR(trK);
    
    CCTK_REAL M1L = gu21*(-(G112*kxxL) + (G111 - G212)*kxyL - G312*kxzL + 
      G211*kyyL + G311*kyzL - PDstandardNth1kxy + PDstandardNth2kxx) + 
      gu22*(-(G122*kxxL) + (G112 - G222)*kxyL - G322*kxzL + G212*kyyL + 
      G312*kyzL - PDstandardNth1kyy + PDstandardNth2kxy) + gu31*(-(G113*kxxL) 
      - G213*kxyL + (G111 - G313)*kxzL + G211*kyzL + G311*kzzL - 
      PDstandardNth1kxz + PDstandardNth3kxx) + gu32*(G113*kxyL + G112*kxzL + 
      G213*kyyL + (G212 + G313)*kyzL + G312*kzzL - 2*(G123*kxxL + G223*kxyL + 
      G323*kxzL + PDstandardNth1kyz) + PDstandardNth2kxz + PDstandardNth3kxy) 
      + gu33*(-(G133*kxxL) - G233*kxyL + (G113 - G333)*kxzL + G213*kyzL + 
      G313*kzzL - PDstandardNth1kzz + PDstandardNth3kxz) - 
      25.13274122871834590770114706623602307358*S1;
    
    CCTK_REAL M2L = gu11*(G112*kxxL + (-G111 + G212)*kxyL + G312*kxzL - 
      G211*kyyL - G311*kyzL + PDstandardNth1kxy - PDstandardNth2kxx) + 
      gu21*(G122*kxxL + (-G112 + G222)*kxyL + G322*kxzL - G212*kyyL - 
      G312*kyzL + PDstandardNth1kyy - PDstandardNth2kxy) + gu31*(G123*kxxL + 
      (-2*G113 + G223)*kxyL + (G112 + G323)*kxzL + G212*kyzL + G312*kzzL + 
      PDstandardNth1kyz - 2*(G213*kyyL + G313*kyzL + PDstandardNth2kxz) + 
      PDstandardNth3kxy) + gu32*(-(G123*kxyL) + G122*kxzL - G223*kyyL + (G222 
      - G323)*kyzL + G322*kzzL - PDstandardNth2kyz + PDstandardNth3kyy) + 
      gu33*(-(G133*kxyL) + G123*kxzL - G233*kyyL + (G223 - G333)*kyzL + 
      G323*kzzL - PDstandardNth2kzz + PDstandardNth3kyz) - 
      25.13274122871834590770114706623602307358*S2;
    
    CCTK_REAL M3L = gu11*(G113*kxxL + G213*kxyL + (-G111 + G313)*kxzL - 
      G211*kyzL - G311*kzzL + PDstandardNth1kxz - PDstandardNth3kxx) + 
      gu21*(G123*kxxL + (G113 + G223)*kxyL + (-2*G112 + G323)*kxzL + 
      G213*kyyL + (-2*G212 + G313)*kyzL + PDstandardNth1kyz + 
      PDstandardNth2kxz - 2*(G312*kzzL + PDstandardNth3kxy)) + 
      gu31*(G133*kxxL + G233*kxyL + (-G113 + G333)*kxzL - G213*kyzL - 
      G313*kzzL + PDstandardNth1kzz - PDstandardNth3kxz) + gu22*(G123*kxyL - 
      G122*kxzL + G223*kyyL + (-G222 + G323)*kyzL - G322*kzzL + 
      PDstandardNth2kyz - PDstandardNth3kyy) + gu32*(G133*kxyL - G123*kxzL + 
      G233*kyyL + (-G223 + G333)*kyzL - G323*kzzL + PDstandardNth2kzz - 
      PDstandardNth3kyz) - 25.13274122871834590770114706623602307358*S3;
    
    /* Copy local copies back to grid functions */
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
  }
  LC_ENDLOOP3 (ML_ADMConstraints_O2);
}

extern "C" void ML_ADMConstraints_O2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADMConstraints_O2_Body);
}
