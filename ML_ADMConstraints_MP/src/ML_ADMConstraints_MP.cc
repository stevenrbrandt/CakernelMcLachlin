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

extern "C" void ML_ADMConstraints_MP_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMConstraints_MP::ML_Ham","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMConstraints_MP::ML_Ham.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMConstraints_MP::ML_mom","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_ADMConstraints_MP::ML_mom.");
  return;
}

static void ML_ADMConstraints_MP_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMConstraints_MP_Body");
  }
  
  if (cctk_iteration % ML_ADMConstraints_MP_calc_every != ML_ADMConstraints_MP_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::curv","ADMBase::lapse","ADMBase::metric","ADMBase::shift","Coordinates::jacobian","Coordinates::jacobian2","ML_ADMConstraints_MP::ML_Ham","ML_ADMConstraints_MP::ML_mom"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADMConstraints_MP", 8, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_ADMConstraints_MP", 2, 2, 2);
  
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
  LC_LOOP3 (ML_ADMConstraints_MP,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL = alp[index];
    CCTK_REAL betaxL = betax[index];
    CCTK_REAL betayL = betay[index];
    CCTK_REAL betazL = betaz[index];
    CCTK_REAL dJ111L = dJ111[index];
    CCTK_REAL dJ112L = dJ112[index];
    CCTK_REAL dJ113L = dJ113[index];
    CCTK_REAL dJ122L = dJ122[index];
    CCTK_REAL dJ123L = dJ123[index];
    CCTK_REAL dJ133L = dJ133[index];
    CCTK_REAL dJ211L = dJ211[index];
    CCTK_REAL dJ212L = dJ212[index];
    CCTK_REAL dJ213L = dJ213[index];
    CCTK_REAL dJ222L = dJ222[index];
    CCTK_REAL dJ223L = dJ223[index];
    CCTK_REAL dJ233L = dJ233[index];
    CCTK_REAL dJ311L = dJ311[index];
    CCTK_REAL dJ312L = dJ312[index];
    CCTK_REAL dJ313L = dJ313[index];
    CCTK_REAL dJ322L = dJ322[index];
    CCTK_REAL dJ323L = dJ323[index];
    CCTK_REAL dJ333L = dJ333[index];
    CCTK_REAL gxxL = gxx[index];
    CCTK_REAL gxyL = gxy[index];
    CCTK_REAL gxzL = gxz[index];
    CCTK_REAL gyyL = gyy[index];
    CCTK_REAL gyzL = gyz[index];
    CCTK_REAL gzzL = gzz[index];
    CCTK_REAL J11L = J11[index];
    CCTK_REAL J12L = J12[index];
    CCTK_REAL J13L = J13[index];
    CCTK_REAL J21L = J21[index];
    CCTK_REAL J22L = J22[index];
    CCTK_REAL J23L = J23[index];
    CCTK_REAL J31L = J31[index];
    CCTK_REAL J32L = J32[index];
    CCTK_REAL J33L = J33[index];
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
    CCTK_REAL const PDstandardNth11gxx = PDstandardNth11(&gxx[index]);
    CCTK_REAL const PDstandardNth22gxx = PDstandardNth22(&gxx[index]);
    CCTK_REAL const PDstandardNth33gxx = PDstandardNth33(&gxx[index]);
    CCTK_REAL const PDstandardNth12gxx = PDstandardNth12(&gxx[index]);
    CCTK_REAL const PDstandardNth13gxx = PDstandardNth13(&gxx[index]);
    CCTK_REAL const PDstandardNth23gxx = PDstandardNth23(&gxx[index]);
    CCTK_REAL const PDstandardNth1gxy = PDstandardNth1(&gxy[index]);
    CCTK_REAL const PDstandardNth2gxy = PDstandardNth2(&gxy[index]);
    CCTK_REAL const PDstandardNth3gxy = PDstandardNth3(&gxy[index]);
    CCTK_REAL const PDstandardNth11gxy = PDstandardNth11(&gxy[index]);
    CCTK_REAL const PDstandardNth22gxy = PDstandardNth22(&gxy[index]);
    CCTK_REAL const PDstandardNth33gxy = PDstandardNth33(&gxy[index]);
    CCTK_REAL const PDstandardNth12gxy = PDstandardNth12(&gxy[index]);
    CCTK_REAL const PDstandardNth13gxy = PDstandardNth13(&gxy[index]);
    CCTK_REAL const PDstandardNth23gxy = PDstandardNth23(&gxy[index]);
    CCTK_REAL const PDstandardNth1gxz = PDstandardNth1(&gxz[index]);
    CCTK_REAL const PDstandardNth2gxz = PDstandardNth2(&gxz[index]);
    CCTK_REAL const PDstandardNth3gxz = PDstandardNth3(&gxz[index]);
    CCTK_REAL const PDstandardNth11gxz = PDstandardNth11(&gxz[index]);
    CCTK_REAL const PDstandardNth22gxz = PDstandardNth22(&gxz[index]);
    CCTK_REAL const PDstandardNth33gxz = PDstandardNth33(&gxz[index]);
    CCTK_REAL const PDstandardNth12gxz = PDstandardNth12(&gxz[index]);
    CCTK_REAL const PDstandardNth13gxz = PDstandardNth13(&gxz[index]);
    CCTK_REAL const PDstandardNth23gxz = PDstandardNth23(&gxz[index]);
    CCTK_REAL const PDstandardNth1gyy = PDstandardNth1(&gyy[index]);
    CCTK_REAL const PDstandardNth2gyy = PDstandardNth2(&gyy[index]);
    CCTK_REAL const PDstandardNth3gyy = PDstandardNth3(&gyy[index]);
    CCTK_REAL const PDstandardNth11gyy = PDstandardNth11(&gyy[index]);
    CCTK_REAL const PDstandardNth22gyy = PDstandardNth22(&gyy[index]);
    CCTK_REAL const PDstandardNth33gyy = PDstandardNth33(&gyy[index]);
    CCTK_REAL const PDstandardNth12gyy = PDstandardNth12(&gyy[index]);
    CCTK_REAL const PDstandardNth13gyy = PDstandardNth13(&gyy[index]);
    CCTK_REAL const PDstandardNth23gyy = PDstandardNth23(&gyy[index]);
    CCTK_REAL const PDstandardNth1gyz = PDstandardNth1(&gyz[index]);
    CCTK_REAL const PDstandardNth2gyz = PDstandardNth2(&gyz[index]);
    CCTK_REAL const PDstandardNth3gyz = PDstandardNth3(&gyz[index]);
    CCTK_REAL const PDstandardNth11gyz = PDstandardNth11(&gyz[index]);
    CCTK_REAL const PDstandardNth22gyz = PDstandardNth22(&gyz[index]);
    CCTK_REAL const PDstandardNth33gyz = PDstandardNth33(&gyz[index]);
    CCTK_REAL const PDstandardNth12gyz = PDstandardNth12(&gyz[index]);
    CCTK_REAL const PDstandardNth13gyz = PDstandardNth13(&gyz[index]);
    CCTK_REAL const PDstandardNth23gyz = PDstandardNth23(&gyz[index]);
    CCTK_REAL const PDstandardNth1gzz = PDstandardNth1(&gzz[index]);
    CCTK_REAL const PDstandardNth2gzz = PDstandardNth2(&gzz[index]);
    CCTK_REAL const PDstandardNth3gzz = PDstandardNth3(&gzz[index]);
    CCTK_REAL const PDstandardNth11gzz = PDstandardNth11(&gzz[index]);
    CCTK_REAL const PDstandardNth22gzz = PDstandardNth22(&gzz[index]);
    CCTK_REAL const PDstandardNth33gzz = PDstandardNth33(&gzz[index]);
    CCTK_REAL const PDstandardNth12gzz = PDstandardNth12(&gzz[index]);
    CCTK_REAL const PDstandardNth13gzz = PDstandardNth13(&gzz[index]);
    CCTK_REAL const PDstandardNth23gzz = PDstandardNth23(&gzz[index]);
    CCTK_REAL const PDstandardNth1kxx = PDstandardNth1(&kxx[index]);
    CCTK_REAL const PDstandardNth2kxx = PDstandardNth2(&kxx[index]);
    CCTK_REAL const PDstandardNth3kxx = PDstandardNth3(&kxx[index]);
    CCTK_REAL const PDstandardNth1kxy = PDstandardNth1(&kxy[index]);
    CCTK_REAL const PDstandardNth2kxy = PDstandardNth2(&kxy[index]);
    CCTK_REAL const PDstandardNth3kxy = PDstandardNth3(&kxy[index]);
    CCTK_REAL const PDstandardNth1kxz = PDstandardNth1(&kxz[index]);
    CCTK_REAL const PDstandardNth2kxz = PDstandardNth2(&kxz[index]);
    CCTK_REAL const PDstandardNth3kxz = PDstandardNth3(&kxz[index]);
    CCTK_REAL const PDstandardNth1kyy = PDstandardNth1(&kyy[index]);
    CCTK_REAL const PDstandardNth2kyy = PDstandardNth2(&kyy[index]);
    CCTK_REAL const PDstandardNth3kyy = PDstandardNth3(&kyy[index]);
    CCTK_REAL const PDstandardNth1kyz = PDstandardNth1(&kyz[index]);
    CCTK_REAL const PDstandardNth2kyz = PDstandardNth2(&kyz[index]);
    CCTK_REAL const PDstandardNth3kyz = PDstandardNth3(&kyz[index]);
    CCTK_REAL const PDstandardNth1kzz = PDstandardNth1(&kzz[index]);
    CCTK_REAL const PDstandardNth2kzz = PDstandardNth2(&kzz[index]);
    CCTK_REAL const PDstandardNth3kzz = PDstandardNth3(&kzz[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL detg = 2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - SQR(gxyL)) - 
      gyyL*SQR(gxzL) - gxxL*SQR(gyzL);
    
    CCTK_REAL gu11 = INV(detg)*(gyyL*gzzL - SQR(gyzL));
    
    CCTK_REAL gu21 = (gxzL*gyzL - gxyL*gzzL)*INV(detg);
    
    CCTK_REAL gu31 = (-(gxzL*gyyL) + gxyL*gyzL)*INV(detg);
    
    CCTK_REAL gu22 = INV(detg)*(gxxL*gzzL - SQR(gxzL));
    
    CCTK_REAL gu32 = (gxyL*gxzL - gxxL*gyzL)*INV(detg);
    
    CCTK_REAL gu33 = INV(detg)*(gxxL*gyyL - SQR(gxyL));
    
    CCTK_REAL G111 = 0.5*((gu11*J11L - gu21*J12L - 
      gu31*J13L)*PDstandardNth1gxx + (gu11*J21L - gu21*J22L - 
      gu31*J23L)*PDstandardNth2gxx + (gu11*J31L - gu21*J32L - 
      gu31*J33L)*PDstandardNth3gxx + 2*(J11L*(gu21*PDstandardNth1gxy + 
      gu31*PDstandardNth1gxz) + J21L*(gu21*PDstandardNth2gxy + 
      gu31*PDstandardNth2gxz) + J31L*(gu21*PDstandardNth3gxy + 
      gu31*PDstandardNth3gxz)));
    
    CCTK_REAL G211 = 0.5*((gu21*J11L - gu22*J12L - 
      gu32*J13L)*PDstandardNth1gxx + (gu21*J21L - gu22*J22L - 
      gu32*J23L)*PDstandardNth2gxx + (gu21*J31L - gu22*J32L - 
      gu32*J33L)*PDstandardNth3gxx + 2*(J11L*(gu22*PDstandardNth1gxy + 
      gu32*PDstandardNth1gxz) + J21L*(gu22*PDstandardNth2gxy + 
      gu32*PDstandardNth2gxz) + J31L*(gu22*PDstandardNth3gxy + 
      gu32*PDstandardNth3gxz)));
    
    CCTK_REAL G311 = 0.5*((gu31*J11L - gu32*J12L - 
      gu33*J13L)*PDstandardNth1gxx + (gu31*J21L - gu32*J22L - 
      gu33*J23L)*PDstandardNth2gxx + (gu31*J31L - gu32*J32L - 
      gu33*J33L)*PDstandardNth3gxx + 2*(J11L*(gu32*PDstandardNth1gxy + 
      gu33*PDstandardNth1gxz) + J21L*(gu32*PDstandardNth2gxy + 
      gu33*PDstandardNth2gxz) + J31L*(gu32*PDstandardNth3gxy + 
      gu33*PDstandardNth3gxz)));
    
    CCTK_REAL G112 = 0.5*(gu11*(J12L*PDstandardNth1gxx + 
      J22L*PDstandardNth2gxx + J32L*PDstandardNth3gxx) + 
      gu21*(J11L*PDstandardNth1gyy + J21L*PDstandardNth2gyy + 
      J31L*PDstandardNth3gyy) + gu31*(-(J13L*PDstandardNth1gxy) + 
      J12L*PDstandardNth1gxz + J11L*PDstandardNth1gyz - 
      J23L*PDstandardNth2gxy + J22L*PDstandardNth2gxz + 
      J21L*PDstandardNth2gyz - J33L*PDstandardNth3gxy + 
      J32L*PDstandardNth3gxz + J31L*PDstandardNth3gyz));
    
    CCTK_REAL G212 = 0.5*(gu21*(J12L*PDstandardNth1gxx + 
      J22L*PDstandardNth2gxx + J32L*PDstandardNth3gxx) + 
      gu22*(J11L*PDstandardNth1gyy + J21L*PDstandardNth2gyy + 
      J31L*PDstandardNth3gyy) + gu32*(-(J13L*PDstandardNth1gxy) + 
      J12L*PDstandardNth1gxz + J11L*PDstandardNth1gyz - 
      J23L*PDstandardNth2gxy + J22L*PDstandardNth2gxz + 
      J21L*PDstandardNth2gyz - J33L*PDstandardNth3gxy + 
      J32L*PDstandardNth3gxz + J31L*PDstandardNth3gyz));
    
    CCTK_REAL G312 = 0.5*(gu31*(J12L*PDstandardNth1gxx + 
      J22L*PDstandardNth2gxx + J32L*PDstandardNth3gxx) + 
      gu32*(J11L*PDstandardNth1gyy + J21L*PDstandardNth2gyy + 
      J31L*PDstandardNth3gyy) + gu33*(-(J13L*PDstandardNth1gxy) + 
      J12L*PDstandardNth1gxz + J11L*PDstandardNth1gyz - 
      J23L*PDstandardNth2gxy + J22L*PDstandardNth2gxz + 
      J21L*PDstandardNth2gyz - J33L*PDstandardNth3gxy + 
      J32L*PDstandardNth3gxz + J31L*PDstandardNth3gyz));
    
    CCTK_REAL G113 = 0.5*(gu11*(J13L*PDstandardNth1gxx + 
      J23L*PDstandardNth2gxx + J33L*PDstandardNth3gxx) + 
      gu21*(J13L*PDstandardNth1gxy - J12L*PDstandardNth1gxz + 
      J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy - 
      J22L*PDstandardNth2gxz + J21L*PDstandardNth2gyz + 
      J33L*PDstandardNth3gxy - J32L*PDstandardNth3gxz + 
      J31L*PDstandardNth3gyz) + gu31*(J11L*PDstandardNth1gzz + 
      J21L*PDstandardNth2gzz + J31L*PDstandardNth3gzz));
    
    CCTK_REAL G213 = 0.5*(gu21*(J13L*PDstandardNth1gxx + 
      J23L*PDstandardNth2gxx + J33L*PDstandardNth3gxx) + 
      gu22*(J13L*PDstandardNth1gxy - J12L*PDstandardNth1gxz + 
      J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy - 
      J22L*PDstandardNth2gxz + J21L*PDstandardNth2gyz + 
      J33L*PDstandardNth3gxy - J32L*PDstandardNth3gxz + 
      J31L*PDstandardNth3gyz) + gu32*(J11L*PDstandardNth1gzz + 
      J21L*PDstandardNth2gzz + J31L*PDstandardNth3gzz));
    
    CCTK_REAL G313 = 0.5*(gu31*(J13L*PDstandardNth1gxx + 
      J23L*PDstandardNth2gxx + J33L*PDstandardNth3gxx) + 
      gu32*(J13L*PDstandardNth1gxy - J12L*PDstandardNth1gxz + 
      J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy - 
      J22L*PDstandardNth2gxz + J21L*PDstandardNth2gyz + 
      J33L*PDstandardNth3gxy - J32L*PDstandardNth3gxz + 
      J31L*PDstandardNth3gyz) + gu33*(J11L*PDstandardNth1gzz + 
      J21L*PDstandardNth2gzz + J31L*PDstandardNth3gzz));
    
    CCTK_REAL G122 = 0.5*(gu11*(-(J11L*PDstandardNth1gyy) + 
      2*(J12L*PDstandardNth1gxy + J22L*PDstandardNth2gxy) - 
      J21L*PDstandardNth2gyy + 2*J32L*PDstandardNth3gxy - 
      J31L*PDstandardNth3gyy) + gu21*(J12L*PDstandardNth1gyy + 
      J22L*PDstandardNth2gyy + J32L*PDstandardNth3gyy) - 
      gu31*(J13L*PDstandardNth1gyy + J23L*PDstandardNth2gyy + 
      J33L*PDstandardNth3gyy - 2*(J12L*PDstandardNth1gyz + 
      J22L*PDstandardNth2gyz + J32L*PDstandardNth3gyz)));
    
    CCTK_REAL G222 = 0.5*(gu21*(-(J11L*PDstandardNth1gyy) + 
      2*(J12L*PDstandardNth1gxy + J22L*PDstandardNth2gxy) - 
      J21L*PDstandardNth2gyy + 2*J32L*PDstandardNth3gxy - 
      J31L*PDstandardNth3gyy) + gu22*(J12L*PDstandardNth1gyy + 
      J22L*PDstandardNth2gyy + J32L*PDstandardNth3gyy) - 
      gu32*(J13L*PDstandardNth1gyy + J23L*PDstandardNth2gyy + 
      J33L*PDstandardNth3gyy - 2*(J12L*PDstandardNth1gyz + 
      J22L*PDstandardNth2gyz + J32L*PDstandardNth3gyz)));
    
    CCTK_REAL G322 = 0.5*(gu31*(-(J11L*PDstandardNth1gyy) + 
      2*(J12L*PDstandardNth1gxy + J22L*PDstandardNth2gxy) - 
      J21L*PDstandardNth2gyy + 2*J32L*PDstandardNth3gxy - 
      J31L*PDstandardNth3gyy) + gu32*(J12L*PDstandardNth1gyy + 
      J22L*PDstandardNth2gyy + J32L*PDstandardNth3gyy) - 
      gu33*(J13L*PDstandardNth1gyy + J23L*PDstandardNth2gyy + 
      J33L*PDstandardNth3gyy - 2*(J12L*PDstandardNth1gyz + 
      J22L*PDstandardNth2gyz + J32L*PDstandardNth3gyz)));
    
    CCTK_REAL G123 = 0.5*(gu21*(J13L*PDstandardNth1gyy + 
      J23L*PDstandardNth2gyy + J33L*PDstandardNth3gyy) + 
      gu11*(J13L*PDstandardNth1gxy + J12L*PDstandardNth1gxz - 
      J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy + 
      J22L*PDstandardNth2gxz - J21L*PDstandardNth2gyz + 
      J33L*PDstandardNth3gxy + J32L*PDstandardNth3gxz - 
      J31L*PDstandardNth3gyz) + gu31*(J12L*PDstandardNth1gzz + 
      J22L*PDstandardNth2gzz + J32L*PDstandardNth3gzz));
    
    CCTK_REAL G223 = 0.5*(gu22*(J13L*PDstandardNth1gyy + 
      J23L*PDstandardNth2gyy + J33L*PDstandardNth3gyy) + 
      gu21*(J13L*PDstandardNth1gxy + J12L*PDstandardNth1gxz - 
      J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy + 
      J22L*PDstandardNth2gxz - J21L*PDstandardNth2gyz + 
      J33L*PDstandardNth3gxy + J32L*PDstandardNth3gxz - 
      J31L*PDstandardNth3gyz) + gu32*(J12L*PDstandardNth1gzz + 
      J22L*PDstandardNth2gzz + J32L*PDstandardNth3gzz));
    
    CCTK_REAL G323 = 0.5*(gu32*(J13L*PDstandardNth1gyy + 
      J23L*PDstandardNth2gyy + J33L*PDstandardNth3gyy) + 
      gu31*(J13L*PDstandardNth1gxy + J12L*PDstandardNth1gxz - 
      J11L*PDstandardNth1gyz + J23L*PDstandardNth2gxy + 
      J22L*PDstandardNth2gxz - J21L*PDstandardNth2gyz + 
      J33L*PDstandardNth3gxy + J32L*PDstandardNth3gxz - 
      J31L*PDstandardNth3gyz) + gu33*(J12L*PDstandardNth1gzz + 
      J22L*PDstandardNth2gzz + J32L*PDstandardNth3gzz));
    
    CCTK_REAL G133 = 0.5*(gu11*(-(J11L*PDstandardNth1gzz) + 
      2*(J13L*PDstandardNth1gxz + J23L*PDstandardNth2gxz) - 
      J21L*PDstandardNth2gzz + 2*J33L*PDstandardNth3gxz - 
      J31L*PDstandardNth3gzz) + gu21*(-(J12L*PDstandardNth1gzz) + 
      2*(J13L*PDstandardNth1gyz + J23L*PDstandardNth2gyz) - 
      J22L*PDstandardNth2gzz + 2*J33L*PDstandardNth3gyz - 
      J32L*PDstandardNth3gzz) + gu31*(J13L*PDstandardNth1gzz + 
      J23L*PDstandardNth2gzz + J33L*PDstandardNth3gzz));
    
    CCTK_REAL G233 = 0.5*(gu21*(-(J11L*PDstandardNth1gzz) + 
      2*(J13L*PDstandardNth1gxz + J23L*PDstandardNth2gxz) - 
      J21L*PDstandardNth2gzz + 2*J33L*PDstandardNth3gxz - 
      J31L*PDstandardNth3gzz) + gu22*(-(J12L*PDstandardNth1gzz) + 
      2*(J13L*PDstandardNth1gyz + J23L*PDstandardNth2gyz) - 
      J22L*PDstandardNth2gzz + 2*J33L*PDstandardNth3gyz - 
      J32L*PDstandardNth3gzz) + gu32*(J13L*PDstandardNth1gzz + 
      J23L*PDstandardNth2gzz + J33L*PDstandardNth3gzz));
    
    CCTK_REAL G333 = 0.5*(gu31*(-(J11L*PDstandardNth1gzz) + 
      2*(J13L*PDstandardNth1gxz + J23L*PDstandardNth2gxz) - 
      J21L*PDstandardNth2gzz + 2*J33L*PDstandardNth3gxz - 
      J31L*PDstandardNth3gzz) + gu32*(-(J12L*PDstandardNth1gzz) + 
      2*(J13L*PDstandardNth1gyz + J23L*PDstandardNth2gyz) - 
      J22L*PDstandardNth2gzz + 2*J33L*PDstandardNth3gyz - 
      J32L*PDstandardNth3gzz) + gu33*(J13L*PDstandardNth1gzz + 
      J23L*PDstandardNth2gzz + J33L*PDstandardNth3gzz));
    
    CCTK_REAL R11 = 0.5*(-((dJ122L*gu22 + 2*dJ123L*gu32 + 
      dJ133L*gu33)*PDstandardNth1gxx) + gu32*(-4*J11L*J21L*PDstandardNth12gyz 
      - 2*J13L*J32L*PDstandardNth13gxx + 2*J13L*J31L*PDstandardNth13gxy + 
      2*J12L*J31L*PDstandardNth13gxz - 4*J11L*J31L*PDstandardNth13gyz + 
      2*dJ113L*PDstandardNth1gxy) + 2*(gu22*J11L*J12L*PDstandardNth11gxy + 
      gu32*J11L*J13L*PDstandardNth11gxy + gu32*J11L*J12L*PDstandardNth11gxz + 
      gu33*J11L*J13L*PDstandardNth11gxz + gu22*J12L*J21L*PDstandardNth12gxy + 
      gu32*J13L*J21L*PDstandardNth12gxy + gu22*J11L*J22L*PDstandardNth12gxy + 
      gu32*J11L*J23L*PDstandardNth12gxy + gu32*J12L*J21L*PDstandardNth12gxz + 
      gu33*J13L*J21L*PDstandardNth12gxz + gu32*J11L*J22L*PDstandardNth12gxz + 
      gu33*J11L*J23L*PDstandardNth12gxz + gu22*J12L*J31L*PDstandardNth13gxy + 
      gu22*J11L*J32L*PDstandardNth13gxy + gu32*J11L*J33L*PDstandardNth13gxy + 
      gu33*J13L*J31L*PDstandardNth13gxz + gu32*J11L*J32L*PDstandardNth13gxz + 
      gu33*J11L*J33L*PDstandardNth13gxz + dJ112L*gu22*PDstandardNth1gxy) + 
      2*dJ112L*gu32*PDstandardNth1gxz - 2*(G111*G122 + G111*G133 + G211*G222 
      + G211*G233 + G311*G322 + G311*G333 + gu32*J12L*J13L*PDstandardNth11gxx 
      + gu22*J12L*J22L*PDstandardNth12gxx + gu32*J13L*J22L*PDstandardNth12gxx 
      + gu32*J12L*J23L*PDstandardNth12gxx + gu33*J13L*J23L*PDstandardNth12gxx 
      + gu22*J11L*J21L*PDstandardNth12gyy + gu33*J11L*J21L*PDstandardNth12gzz 
      + gu22*J12L*J32L*PDstandardNth13gxx + gu32*J12L*J33L*PDstandardNth13gxx 
      + gu33*J13L*J33L*PDstandardNth13gxx + gu22*J11L*J31L*PDstandardNth13gyy 
      + gu33*J11L*J31L*PDstandardNth13gzz + dJ111L*gu32*PDstandardNth1gyz) + 
      2*gu32*J21L*J23L*PDstandardNth22gxy + 
      2*gu32*J21L*J22L*PDstandardNth22gxz - 
      2*gu32*J23L*J32L*PDstandardNth23gxx - 
      2*gu32*J22L*J33L*PDstandardNth23gxx + 
      2*gu32*J23L*J31L*PDstandardNth23gxy + 
      2*gu32*J21L*J33L*PDstandardNth23gxy + 
      2*gu32*J22L*J31L*PDstandardNth23gxz + 
      2*gu32*J21L*J32L*PDstandardNth23gxz - 
      4*gu32*J21L*J31L*PDstandardNth23gyz - (dJ222L*gu22 + 2*dJ223L*gu32 + 
      dJ233L*gu33)*PDstandardNth2gxx + 2*dJ213L*gu32*PDstandardNth2gxy + 
      2*dJ212L*gu32*PDstandardNth2gxz - 2*dJ211L*gu32*PDstandardNth2gyz + 
      2*gu22*J31L*J32L*PDstandardNth33gxy + 
      2*gu32*J31L*J33L*PDstandardNth33gxy + 
      2*gu32*J31L*J32L*PDstandardNth33gxz - 2*dJ323L*gu32*PDstandardNth3gxx + 
      gu22*(-2*J22L*J32L*PDstandardNth23gxx - dJ322L*PDstandardNth3gxx) + 
      gu33*(2*J21L*J23L*PDstandardNth22gxz - dJ333L*PDstandardNth3gxx) + 
      2*dJ312L*gu22*PDstandardNth3gxy + 2*dJ313L*gu32*PDstandardNth3gxy + 
      2*dJ312L*gu32*PDstandardNth3gxz + 2*dJ313L*gu33*PDstandardNth3gxz + 
      gu22*(2*J22L*J31L*PDstandardNth23gxy - dJ311L*PDstandardNth3gyy) - 
      2*dJ311L*gu32*PDstandardNth3gyz + gu33*(-2*J23L*J33L*PDstandardNth23gxx 
      - dJ311L*PDstandardNth3gzz) + 2*SQR(G112) + 2*SQR(G113) + 2*SQR(G212) + 
      2*SQR(G213) + 2*SQR(G312) + 2*SQR(G313) - 
      2*gu32*PDstandardNth11gyz*SQR(J11L) + 
      gu22*(2*J21L*J32L*PDstandardNth23gxy - PDstandardNth11gyy*SQR(J11L)) + 
      gu33*(2*J23L*J31L*PDstandardNth23gxz - PDstandardNth11gzz*SQR(J11L)) + 
      gu22*(-2*J21L*J31L*PDstandardNth23gyy - PDstandardNth11gxx*SQR(J12L)) + 
      gu33*(2*J21L*J33L*PDstandardNth23gxz - PDstandardNth11gxx*SQR(J13L)) - 
      2*gu32*PDstandardNth22gyz*SQR(J21L) + gu22*(2*dJ212L*PDstandardNth2gxy 
      - PDstandardNth22gyy*SQR(J21L)) + gu33*(-2*J21L*J31L*PDstandardNth23gzz 
      - PDstandardNth22gzz*SQR(J21L)) + PDstandardNth22gxx*(-2*gu32*J22L*J23L 
      - gu22*SQR(J22L)) + gu33*(2*dJ213L*PDstandardNth2gxz - 
      PDstandardNth22gxx*SQR(J23L)) - 2*gu32*PDstandardNth33gyz*SQR(J31L) + 
      gu22*(-(dJ111L*PDstandardNth1gyy) + 2*J21L*J22L*PDstandardNth22gxy - 
      dJ211L*PDstandardNth2gyy - PDstandardNth33gyy*SQR(J31L)) + 
      gu33*(2*dJ113L*PDstandardNth1gxz - dJ111L*PDstandardNth1gzz - 
      dJ211L*PDstandardNth2gzz - PDstandardNth33gzz*SQR(J31L)) + 
      PDstandardNth33gxx*(-2*gu32*J32L*J33L - gu22*SQR(J32L)) + 
      gu33*(2*J31L*J33L*PDstandardNth33gxz - PDstandardNth33gxx*SQR(J33L)));
    
    CCTK_REAL R12 = 0.5*(-2*(G112*G133 + G212*G233 + G312*G333 + 
      gu21*J11L*J12L*PDstandardNth11gxy) + 2*(G113*G123 + G213*G223 + 
      G313*G323 + gu21*J12L*J22L*PDstandardNth12gxx) + 
      J12L*(-(gu32*J13L*PDstandardNth11gxy) + 
      J11L*(-(gu31*PDstandardNth11gxz) - gu33*PDstandardNth11gzz) - 
      gu31*J21L*PDstandardNth12gxz) + J11L*J22L*(-2*gu21*PDstandardNth12gxy - 
      gu32*PDstandardNth12gyz) + J12L*(gu31*J23L*PDstandardNth12gxx - 
      gu32*J21L*PDstandardNth12gyz) + J12L*J21L*(-2*gu21*PDstandardNth12gxy - 
      gu33*PDstandardNth12gzz) + J11L*(gu33*J13L*PDstandardNth11gyz + 
      gu31*(-(J23L*PDstandardNth12gxy) - J33L*PDstandardNth13gxy)) + 
      J12L*(gu33*J13L*PDstandardNth11gxz + gu32*(-(J11L*PDstandardNth11gyz) - 
      J23L*PDstandardNth12gxy - J33L*PDstandardNth13gxy)) + 
      J12L*(2*gu32*J22L*PDstandardNth12gxz - gu31*J31L*PDstandardNth13gxz) + 
      J12L*(gu33*J23L*PDstandardNth12gxz - gu32*J31L*PDstandardNth13gyz) + 
      J11L*(2*gu21*J21L*PDstandardNth12gyy - gu32*J32L*PDstandardNth13gyz) + 
      J11L*(gu32*J23L*PDstandardNth12gyy - gu33*J32L*PDstandardNth13gzz) + 
      (dJ122L*gu21 + dJ123L*gu31)*PDstandardNth1gxx - 
      2*dJ112L*gu21*PDstandardNth1gxy + dJ111L*gu21*PDstandardNth1gyy + 
      gu32*J21L*J23L*PDstandardNth22gyy + gu33*J21L*J23L*PDstandardNth22gyz + 
      gu32*(J13L*J31L*PDstandardNth13gyy - J21L*J22L*PDstandardNth22gyz) + 
      gu33*(-2*J13L*J23L*PDstandardNth12gxy + 
      J22L*(-(J11L*PDstandardNth12gzz) - J21L*PDstandardNth22gzz)) + 
      2*gu21*J22L*J32L*PDstandardNth23gxx + gu31*J23L*J32L*PDstandardNth23gxx 
      + gu31*J22L*J33L*PDstandardNth23gxx - 
      2*gu21*J22L*J31L*PDstandardNth23gxy - 
      2*gu21*J21L*J32L*PDstandardNth23gxy - 
      2*gu33*J23L*J33L*PDstandardNth23gxy + 
      J33L*(-2*gu33*J13L*PDstandardNth13gxy - gu31*J21L*PDstandardNth23gxy) + 
      J33L*(gu33*J12L*PDstandardNth13gxz - gu32*J22L*PDstandardNth23gxy) + 
      J31L*(-2*gu21*J12L*PDstandardNth13gxy - gu31*J23L*PDstandardNth23gxy) + 
      J32L*(2*gu21*J12L*PDstandardNth13gxx - gu32*J23L*PDstandardNth23gxy) + 
      2*gu32*J22L*J32L*PDstandardNth23gxz + gu33*J23L*J32L*PDstandardNth23gxz 
      + gu33*J22L*J33L*PDstandardNth23gxz + 
      J32L*(-2*gu21*J11L*PDstandardNth13gxy - gu31*J21L*PDstandardNth23gxz) + 
      gu31*(-(J11L*J32L*PDstandardNth13gxz) + J22L*(J13L*PDstandardNth12gxx - 
      J11L*PDstandardNth12gxz - J31L*PDstandardNth23gxz)) + 
      2*gu21*J21L*J31L*PDstandardNth23gyy + gu32*J23L*J31L*PDstandardNth23gyy 
      + gu32*J21L*J33L*PDstandardNth23gyy + 
      2*gu31*J21L*J31L*PDstandardNth23gyz + gu33*J23L*J31L*PDstandardNth23gyz 
      + gu33*J21L*J33L*PDstandardNth23gyz + 
      J32L*(gu33*J13L*PDstandardNth13gxz - gu32*J21L*PDstandardNth23gyz) + 
      J31L*(2*gu21*J11L*PDstandardNth13gyy - gu32*J22L*PDstandardNth23gyz) + 
      J31L*(2*gu31*J11L*PDstandardNth13gyz - gu33*J22L*PDstandardNth23gzz) + 
      gu33*(J13L*J31L*PDstandardNth13gyz - J21L*J32L*PDstandardNth23gzz) + 
      (dJ222L*gu21 + dJ223L*gu31)*PDstandardNth2gxx - 
      2*dJ212L*gu21*PDstandardNth2gxy + gu31*(2*J11L*J21L*PDstandardNth12gyz 
      - dJ113L*PDstandardNth1gxy - dJ213L*PDstandardNth2gxy) + 
      gu32*(J11L*J33L*PDstandardNth13gyy - dJ223L*PDstandardNth2gxy) + 
      gu33*(J11L*J33L*PDstandardNth13gyz - dJ233L*PDstandardNth2gxy) + 
      dJ222L*gu32*PDstandardNth2gxz + dJ223L*gu33*PDstandardNth2gxz + 
      gu31*(J13L*J32L*PDstandardNth13gxx - dJ112L*PDstandardNth1gxz - 
      dJ212L*PDstandardNth2gxz) + dJ211L*gu21*PDstandardNth2gyy + 
      dJ213L*gu32*PDstandardNth2gyy + dJ211L*gu31*PDstandardNth2gyz + 
      dJ213L*gu33*PDstandardNth2gyz + gu32*(J13L*J21L*PDstandardNth12gyy - 
      dJ123L*PDstandardNth1gxy - dJ212L*PDstandardNth2gyz) + 
      gu33*(J13L*J22L*PDstandardNth12gxz - J12L*J31L*PDstandardNth13gzz - 
      dJ212L*PDstandardNth2gzz) + gu31*J32L*J33L*PDstandardNth33gxx - 
      2*gu21*J31L*J32L*PDstandardNth33gxy + gu31*(dJ111L*PDstandardNth1gyz - 
      J31L*J33L*PDstandardNth33gxy) + gu32*(dJ122L*PDstandardNth1gxz - 
      J32L*J33L*PDstandardNth33gxy) + gu33*J32L*J33L*PDstandardNth33gxz + 
      gu31*(J22L*J23L*PDstandardNth22gxx - J31L*J32L*PDstandardNth33gxz) + 
      gu32*J31L*J33L*PDstandardNth33gyy + gu33*J31L*J33L*PDstandardNth33gyz + 
      gu32*(dJ113L*PDstandardNth1gyy - J31L*J32L*PDstandardNth33gyz) + 
      gu33*(J13L*J21L*PDstandardNth12gyz - dJ133L*PDstandardNth1gxy - 
      J31L*J32L*PDstandardNth33gzz) + dJ322L*gu21*PDstandardNth3gxx + 
      dJ323L*gu31*PDstandardNth3gxx - 2*dJ312L*gu21*PDstandardNth3gxy + 
      gu31*(J12L*J33L*PDstandardNth13gxx - J21L*J23L*PDstandardNth22gxy - 
      dJ313L*PDstandardNth3gxy) + gu32*(J13L*(J11L*PDstandardNth11gyy - 
      J22L*PDstandardNth12gxy - J32L*PDstandardNth13gxy) - 
      dJ112L*PDstandardNth1gyz - dJ323L*PDstandardNth3gxy) + 
      gu33*(dJ123L*PDstandardNth1gxz - dJ333L*PDstandardNth3gxy) + 
      dJ322L*gu32*PDstandardNth3gxz + dJ323L*gu33*PDstandardNth3gxz + 
      gu31*(J13L*(J12L*PDstandardNth11gxx - J11L*PDstandardNth11gxy - 
      J21L*PDstandardNth12gxy - J31L*PDstandardNth13gxy) - 
      J21L*J22L*PDstandardNth22gxz - dJ312L*PDstandardNth3gxz) + 
      dJ311L*gu21*PDstandardNth3gyy + dJ313L*gu32*PDstandardNth3gyy + 
      dJ311L*gu31*PDstandardNth3gyz + dJ313L*gu33*PDstandardNth3gyz + 
      gu32*(2*J12L*J32L*PDstandardNth13gxz - J22L*J23L*PDstandardNth22gxy - 
      dJ312L*PDstandardNth3gyz) + gu33*(dJ113L*PDstandardNth1gyz - 
      dJ312L*PDstandardNth3gzz) + gu21*PDstandardNth11gyy*SQR(J11L) + 
      gu31*PDstandardNth11gyz*SQR(J11L) + gu21*PDstandardNth11gxx*SQR(J12L) + 
      gu32*PDstandardNth11gxz*SQR(J12L) + gu33*(J11L*J23L*PDstandardNth12gyz 
      - dJ112L*PDstandardNth1gzz - PDstandardNth11gxy*SQR(J13L)) + 
      gu21*PDstandardNth22gyy*SQR(J21L) + gu31*PDstandardNth22gyz*SQR(J21L) + 
      gu21*PDstandardNth22gxx*SQR(J22L) + gu32*PDstandardNth22gxz*SQR(J22L) + 
      PDstandardNth22gxy*(-2*gu21*J21L*J22L - gu33*SQR(J23L)) + 
      gu21*PDstandardNth33gyy*SQR(J31L) + gu31*PDstandardNth33gyz*SQR(J31L) + 
      gu21*PDstandardNth33gxx*SQR(J32L) + gu32*PDstandardNth33gxz*SQR(J32L) + 
      gu33*(J22L*J23L*PDstandardNth22gxz - PDstandardNth33gxy*SQR(J33L)));
    
    CCTK_REAL R13 = 0.5*(2*(G112*G123 + G212*G223 + G312*G323 + 
      gu31*J13L*J23L*PDstandardNth12gxx) + J12L*J23L*(gu22*PDstandardNth12gxy 
      - gu32*PDstandardNth12gxz) - 2*(G113*G122 + G213*G222 + G313*G322 + 
      gu31*J11L*J13L*PDstandardNth11gxz + gu31*J13L*J21L*PDstandardNth12gxz) 
      + J11L*J23L*(-2*gu31*PDstandardNth12gxz - gu32*PDstandardNth12gyz) + 
      J12L*(-2*gu22*J22L*PDstandardNth12gxz - gu21*J31L*PDstandardNth13gxz) + 
      J12L*(J13L*(gu22*PDstandardNth11gxy - gu32*PDstandardNth11gxz) + 
      J11L*(-(gu21*PDstandardNth11gxz) + gu22*PDstandardNth11gyz) - 
      gu21*J21L*PDstandardNth12gxz - gu32*J33L*PDstandardNth13gxz) + 
      J13L*(-(gu22*J21L*PDstandardNth12gyy) + gu32*(2*J23L*PDstandardNth12gxy 
      - J32L*PDstandardNth13gxz)) + J11L*(gu32*J12L*PDstandardNth11gzz + 
      gu21*(-(J22L*PDstandardNth12gxz) - J32L*PDstandardNth13gxz)) + 
      J11L*(J13L*(-(gu22*PDstandardNth11gyy) - gu32*PDstandardNth11gyz) + 
      J23L*(-(gu21*PDstandardNth12gxy) - gu22*PDstandardNth12gyy) - 
      gu22*J33L*PDstandardNth13gyy) + J11L*(2*gu21*J21L*PDstandardNth12gyz - 
      gu32*J33L*PDstandardNth13gyz) + J13L*(J22L*(gu22*PDstandardNth12gxy - 
      gu32*PDstandardNth12gxz) + gu32*(-(J21L*PDstandardNth12gyz) - 
      J31L*PDstandardNth13gyz)) + (dJ123L*gu21 + 
      dJ133L*gu31)*PDstandardNth1gxx - 2*dJ113L*gu31*PDstandardNth1gxz + 
      dJ111L*gu31*PDstandardNth1gzz + J21L*(2*gu31*J11L*PDstandardNth12gzz - 
      gu21*J23L*PDstandardNth22gxy) + gu32*J21L*J22L*PDstandardNth22gzz + 
      gu21*J23L*J32L*PDstandardNth23gxx + gu21*J22L*J33L*PDstandardNth23gxx + 
      2*gu31*J23L*J33L*PDstandardNth23gxx + gu22*J23L*J32L*PDstandardNth23gxy 
      + gu22*J22L*J33L*PDstandardNth23gxy + 
      2*gu32*J23L*J33L*PDstandardNth23gxy + 
      J33L*(2*gu31*J13L*PDstandardNth13gxx - gu21*J21L*PDstandardNth23gxy) + 
      gu21*(-(J13L*J21L*PDstandardNth12gxy) - J11L*J33L*PDstandardNth13gxy + 
      J23L*(J12L*PDstandardNth12gxx - J31L*PDstandardNth23gxy)) - 
      2*gu31*J23L*J31L*PDstandardNth23gxz - 
      2*gu22*J22L*J32L*PDstandardNth23gxz - 
      2*gu31*J21L*J33L*PDstandardNth23gxz + 
      J32L*(-2*gu22*J12L*PDstandardNth13gxz - gu21*J21L*PDstandardNth23gxz) + 
      J31L*(-2*gu31*J13L*PDstandardNth13gxz - gu21*J22L*PDstandardNth23gxz) + 
      J33L*(-2*gu31*J11L*PDstandardNth13gxz - gu32*J22L*PDstandardNth23gxz) + 
      J32L*(gu22*J11L*PDstandardNth13gyz - gu32*J23L*PDstandardNth23gxz) + 
      gu22*(J12L*J31L*PDstandardNth13gyz - J21L*J33L*PDstandardNth23gyy) + 
      gu22*(J12L*J21L*PDstandardNth12gyz + J31L*(-(J13L*PDstandardNth13gyy) - 
      J23L*PDstandardNth23gyy)) + 2*gu21*J21L*J31L*PDstandardNth23gyz + 
      gu22*J22L*J31L*PDstandardNth23gyz + gu22*J21L*J32L*PDstandardNth23gyz + 
      J31L*(2*gu21*J11L*PDstandardNth13gyz - gu32*J23L*PDstandardNth23gyz) + 
      gu32*(J11L*J32L*PDstandardNth13gzz - J21L*J33L*PDstandardNth23gyz) + 
      2*gu31*J21L*J31L*PDstandardNth23gzz + gu32*J22L*J31L*PDstandardNth23gzz 
      + gu32*J21L*J32L*PDstandardNth23gzz + (dJ223L*gu21 + 
      dJ233L*gu31)*PDstandardNth2gxx + dJ223L*gu22*PDstandardNth2gxy + 
      dJ233L*gu32*PDstandardNth2gxy + gu21*(J13L*J32L*PDstandardNth13gxx - 
      dJ113L*PDstandardNth1gxy - dJ213L*PDstandardNth2gxy) - 
      2*dJ213L*gu31*PDstandardNth2gxz + gu21*(J12L*J33L*PDstandardNth13gxx - 
      dJ112L*PDstandardNth1gxz - dJ212L*PDstandardNth2gxz) + 
      gu22*(dJ123L*PDstandardNth1gxy - dJ222L*PDstandardNth2gxz) + 
      gu32*(dJ133L*PDstandardNth1gxy - dJ223L*PDstandardNth2gxz) + 
      gu22*(J11L*J22L*PDstandardNth12gyz - dJ122L*PDstandardNth1gxz - 
      dJ213L*PDstandardNth2gyy) + dJ211L*gu21*PDstandardNth2gyz + 
      dJ212L*gu22*PDstandardNth2gyz + gu32*(J12L*J21L*PDstandardNth12gzz - 
      dJ123L*PDstandardNth1gxz - dJ213L*PDstandardNth2gyz) + 
      dJ211L*gu31*PDstandardNth2gzz + dJ212L*gu32*PDstandardNth2gzz + 
      gu21*J32L*J33L*PDstandardNth33gxx + gu22*J32L*J33L*PDstandardNth33gxy + 
      J31L*(2*gu31*J11L*PDstandardNth13gzz - gu21*J33L*PDstandardNth33gxy) - 
      2*gu31*J31L*J33L*PDstandardNth33gxz + gu21*(dJ111L*PDstandardNth1gyz - 
      J31L*J32L*PDstandardNth33gxz) + gu32*(J11L*J22L*PDstandardNth12gzz - 
      dJ113L*PDstandardNth1gyz - J32L*J33L*PDstandardNth33gxz) + 
      gu22*(J13L*J32L*PDstandardNth13gxy - dJ113L*PDstandardNth1gyy - 
      J31L*J33L*PDstandardNth33gyy) + gu22*J31L*J32L*PDstandardNth33gyz + 
      gu32*(dJ112L*PDstandardNth1gzz - J31L*J33L*PDstandardNth33gyz) + 
      gu32*J31L*J32L*PDstandardNth33gzz + dJ323L*gu21*PDstandardNth3gxx + 
      dJ333L*gu31*PDstandardNth3gxx + dJ323L*gu22*PDstandardNth3gxy + 
      dJ333L*gu32*PDstandardNth3gxy + gu21*(J22L*J23L*PDstandardNth22gxx - 
      dJ313L*PDstandardNth3gxy) - 2*dJ313L*gu31*PDstandardNth3gxz + 
      gu21*(J13L*(J12L*PDstandardNth11gxx - J11L*PDstandardNth11gxy - 
      J31L*PDstandardNth13gxy) + J22L*(J13L*PDstandardNth12gxx - 
      J21L*PDstandardNth22gxz) - dJ312L*PDstandardNth3gxz) + 
      gu22*(dJ112L*PDstandardNth1gyz - dJ322L*PDstandardNth3gxz) + 
      gu32*(2*J13L*J33L*PDstandardNth13gxy - J22L*J23L*PDstandardNth22gxz - 
      dJ323L*PDstandardNth3gxz) + gu22*(J22L*J23L*PDstandardNth22gxy - 
      dJ313L*PDstandardNth3gyy) + dJ311L*gu21*PDstandardNth3gyz + 
      dJ312L*gu22*PDstandardNth3gyz + gu32*(J12L*J31L*PDstandardNth13gzz - 
      J21L*J23L*PDstandardNth22gyz - dJ313L*PDstandardNth3gyz) + 
      dJ311L*gu31*PDstandardNth3gzz + dJ312L*gu32*PDstandardNth3gzz + 
      gu21*PDstandardNth11gyz*SQR(J11L) + gu31*PDstandardNth11gzz*SQR(J11L) + 
      gu22*(J12L*J33L*PDstandardNth13gxy - J21L*J23L*PDstandardNth22gyy - 
      PDstandardNth11gxz*SQR(J12L)) + gu31*PDstandardNth11gxx*SQR(J13L) + 
      gu32*PDstandardNth11gxy*SQR(J13L) + gu21*PDstandardNth22gyz*SQR(J21L) + 
      gu31*PDstandardNth22gzz*SQR(J21L) + 
      PDstandardNth22gxz*(-2*gu31*J21L*J23L - gu22*SQR(J22L)) + 
      gu31*PDstandardNth22gxx*SQR(J23L) + gu32*PDstandardNth22gxy*SQR(J23L) + 
      gu21*PDstandardNth33gyz*SQR(J31L) + gu31*PDstandardNth33gzz*SQR(J31L) + 
      gu22*(J21L*J22L*PDstandardNth22gyz - PDstandardNth33gxz*SQR(J32L)) + 
      gu31*PDstandardNth33gxx*SQR(J33L) + gu32*PDstandardNth33gxy*SQR(J33L));
    
    CCTK_REAL R22 = 0.5*(gu31*(-4*J12L*J22L*PDstandardNth12gxz - 
      2*J13L*J21L*PDstandardNth12gyy + 2*J11L*J22L*PDstandardNth12gyz + 
      2*J13L*J32L*PDstandardNth13gxy - 4*J12L*J32L*PDstandardNth13gxz - 
      2*J13L*J31L*PDstandardNth13gyy + 2*J11L*J32L*PDstandardNth13gyz + 
      2*dJ123L*PDstandardNth1gxy) + 2*(gu11*J11L*J12L*PDstandardNth11gxy + 
      gu31*J12L*J13L*PDstandardNth11gxy + gu31*J11L*J12L*PDstandardNth11gyz + 
      gu33*J12L*J13L*PDstandardNth11gyz + gu11*J12L*J21L*PDstandardNth12gxy + 
      gu11*J11L*J22L*PDstandardNth12gxy + gu31*J13L*J22L*PDstandardNth12gxy + 
      gu31*J12L*J23L*PDstandardNth12gxy + gu31*J12L*J21L*PDstandardNth12gyz + 
      gu33*J13L*J22L*PDstandardNth12gyz + gu33*J12L*J23L*PDstandardNth12gyz + 
      gu11*J12L*J31L*PDstandardNth13gxy + gu11*J11L*J32L*PDstandardNth13gxy + 
      gu31*J12L*J33L*PDstandardNth13gxy + gu31*J12L*J31L*PDstandardNth13gyz + 
      gu33*J13L*J32L*PDstandardNth13gyz + gu33*J12L*J33L*PDstandardNth13gyz + 
      dJ112L*gu11*PDstandardNth1gxy) - 2*(G111*G122 + G122*G133 + G211*G222 + 
      G222*G233 + G311*G322 + G322*G333 + gu31*J11L*J13L*PDstandardNth11gyy + 
      gu11*J12L*J22L*PDstandardNth12gxx + gu11*J11L*J21L*PDstandardNth12gyy + 
      gu31*J11L*J23L*PDstandardNth12gyy + gu33*J13L*J23L*PDstandardNth12gyy + 
      gu33*J12L*J22L*PDstandardNth12gzz + gu11*J12L*J32L*PDstandardNth13gxx + 
      gu11*J11L*J31L*PDstandardNth13gyy + gu31*J11L*J33L*PDstandardNth13gyy + 
      gu33*J13L*J33L*PDstandardNth13gyy + gu33*J12L*J32L*PDstandardNth13gzz + 
      dJ122L*gu31*PDstandardNth1gxz) + (-2*dJ113L*gu31 - 
      dJ133L*gu33)*PDstandardNth1gyy + 2*dJ112L*gu31*PDstandardNth1gyz + 
      2*gu31*J22L*J23L*PDstandardNth22gxy + 
      2*gu31*J21L*J22L*PDstandardNth22gyz + 
      2*gu31*J23L*J32L*PDstandardNth23gxy + 
      2*gu31*J22L*J33L*PDstandardNth23gxy - 
      4*gu31*J22L*J32L*PDstandardNth23gxz - 
      2*gu31*J23L*J31L*PDstandardNth23gyy - 
      2*gu31*J21L*J33L*PDstandardNth23gyy + 
      2*gu31*J22L*J31L*PDstandardNth23gyz + 
      2*gu31*J21L*J32L*PDstandardNth23gyz + 2*dJ223L*gu31*PDstandardNth2gxy - 
      2*dJ222L*gu31*PDstandardNth2gxz - 2*dJ213L*gu31*PDstandardNth2gyy + 
      2*dJ212L*gu31*PDstandardNth2gyz + 2*gu11*J31L*J32L*PDstandardNth33gxy + 
      2*gu31*J32L*J33L*PDstandardNth33gxy + 
      2*gu31*J31L*J32L*PDstandardNth33gyz + 
      2*gu33*J32L*J33L*PDstandardNth33gyz + 
      gu11*(-2*J22L*J32L*PDstandardNth23gxx - dJ322L*PDstandardNth3gxx) + 
      2*dJ312L*gu11*PDstandardNth3gxy + 2*dJ323L*gu31*PDstandardNth3gxy - 
      2*dJ322L*gu31*PDstandardNth3gxz - 2*dJ313L*gu31*PDstandardNth3gyy + 
      gu11*(2*J22L*J31L*PDstandardNth23gxy - dJ311L*PDstandardNth3gyy) + 
      gu33*(-2*J23L*J33L*PDstandardNth23gyy - dJ333L*PDstandardNth3gyy) + 
      2*dJ312L*gu31*PDstandardNth3gyz + 2*dJ323L*gu33*PDstandardNth3gyz + 
      gu33*(2*J23L*J32L*PDstandardNth23gyz - dJ322L*PDstandardNth3gzz) + 
      2*SQR(G112) + 2*SQR(G123) + 2*SQR(G212) + 2*SQR(G223) + 2*SQR(G312) + 
      2*SQR(G323) + gu11*(2*J21L*J32L*PDstandardNth23gxy - 
      PDstandardNth11gyy*SQR(J11L)) - 2*gu31*PDstandardNth11gxz*SQR(J12L) + 
      gu11*(-2*J21L*J31L*PDstandardNth23gyy - PDstandardNth11gxx*SQR(J12L)) + 
      gu33*(2*J22L*J33L*PDstandardNth23gyz - PDstandardNth11gzz*SQR(J12L)) + 
      gu33*(-2*J22L*J32L*PDstandardNth23gzz - PDstandardNth11gyy*SQR(J13L)) + 
      PDstandardNth22gyy*(-2*gu31*J21L*J23L - gu11*SQR(J21L)) - 
      2*gu31*PDstandardNth22gxz*SQR(J22L) + gu11*(-(dJ122L*PDstandardNth1gxx) 
      - dJ111L*PDstandardNth1gyy - dJ222L*PDstandardNth2gxx - 
      PDstandardNth22gxx*SQR(J22L)) + gu33*(2*dJ123L*PDstandardNth1gyz - 
      dJ122L*PDstandardNth1gzz - dJ233L*PDstandardNth2gyy - 
      PDstandardNth22gzz*SQR(J22L)) + gu33*(2*dJ223L*PDstandardNth2gyz - 
      PDstandardNth22gyy*SQR(J23L)) + gu11*(2*dJ212L*PDstandardNth2gxy - 
      PDstandardNth33gyy*SQR(J31L)) - 2*gu31*PDstandardNth33gxz*SQR(J32L) + 
      gu11*(2*J21L*J22L*PDstandardNth22gxy - dJ211L*PDstandardNth2gyy - 
      PDstandardNth33gxx*SQR(J32L)) + gu33*(2*J22L*J23L*PDstandardNth22gyz - 
      dJ222L*PDstandardNth2gzz - PDstandardNth33gzz*SQR(J32L)) + 
      PDstandardNth33gyy*(-2*gu31*J31L*J33L - gu33*SQR(J33L)));
    
    CCTK_REAL R23 = 0.5*(2*(G112*G113 + G212*G213 + G312*G313 + 
      gu31*J13L*J23L*PDstandardNth12gxy + gu32*J13L*J23L*PDstandardNth12gyy) 
      - 2*(G111*G123 + G211*G223 + G311*G323 + 
      gu32*J12L*J13L*PDstandardNth11gyz + gu11*J11L*J21L*PDstandardNth12gyz) 
      + J12L*(-2*gu32*J23L*PDstandardNth12gyz - gu11*J33L*PDstandardNth13gxx) 
      + J12L*(J13L*(-(gu21*PDstandardNth11gxy) - gu31*PDstandardNth11gxz) + 
      J11L*(gu11*PDstandardNth11gxz - gu21*PDstandardNth11gyz) + 
      gu31*J11L*PDstandardNth11gzz + J23L*(-(gu11*PDstandardNth12gxx) - 
      gu21*PDstandardNth12gxy) + gu21*(2*J22L*PDstandardNth12gxz - 
      J21L*PDstandardNth12gyz - J33L*PDstandardNth13gxy)) + 
      J13L*(-2*gu32*J22L*PDstandardNth12gyz - gu31*J32L*PDstandardNth13gxz) + 
      J11L*J33L*(gu11*PDstandardNth13gxy - gu31*PDstandardNth13gyz) + 
      J12L*(gu31*J21L*PDstandardNth12gzz - gu21*J31L*PDstandardNth13gyz) + 
      J11L*(J13L*(gu21*PDstandardNth11gyy - gu31*PDstandardNth11gyz) + 
      J23L*(gu11*PDstandardNth12gxy - gu31*PDstandardNth12gyz) + 
      gu21*(J23L*PDstandardNth12gyy - J32L*PDstandardNth13gyz)) + 
      dJ133L*gu32*PDstandardNth1gyy - 2*dJ123L*gu32*PDstandardNth1gyz + 
      dJ122L*gu32*PDstandardNth1gzz + J22L*(gu21*(-(J13L*PDstandardNth12gxy) 
      - J11L*PDstandardNth12gyz) + gu11*(J11L*PDstandardNth12gxz - 
      J23L*PDstandardNth22gxx)) + J22L*(2*gu32*J12L*PDstandardNth12gzz - 
      gu21*J23L*PDstandardNth22gxy) + gu11*(J11L*J32L*PDstandardNth13gxz - 
      J22L*J33L*PDstandardNth23gxx) + gu11*J23L*J31L*PDstandardNth23gxy + 
      gu11*J21L*J33L*PDstandardNth23gxy + 2*gu31*J23L*J33L*PDstandardNth23gxy 
      + J33L*(2*gu32*J13L*PDstandardNth13gyy - gu21*J22L*PDstandardNth23gxy) 
      + gu21*(J11L*J33L*PDstandardNth13gyy - J23L*J32L*PDstandardNth23gxy) + 
      gu11*J22L*J31L*PDstandardNth23gxz + gu11*J21L*J32L*PDstandardNth23gxz + 
      2*gu21*J22L*J32L*PDstandardNth23gxz + 
      J33L*(-2*gu32*J12L*PDstandardNth13gyz - gu31*J22L*PDstandardNth23gxz) + 
      J32L*(-2*gu32*J13L*PDstandardNth13gyz - gu31*J23L*PDstandardNth23gxz) + 
      gu21*J23L*J31L*PDstandardNth23gyy + gu21*J21L*J33L*PDstandardNth23gyy + 
      2*gu32*J23L*J33L*PDstandardNth23gyy - 
      2*gu11*J21L*J31L*PDstandardNth23gyz - 
      2*gu32*J23L*J32L*PDstandardNth23gyz - 
      2*gu32*J22L*J33L*PDstandardNth23gyz + 
      J32L*(2*gu32*J12L*PDstandardNth13gzz - gu21*J21L*PDstandardNth23gyz) + 
      J31L*(-2*gu11*J11L*PDstandardNth13gyz - gu21*J22L*PDstandardNth23gyz) + 
      gu31*(J11L*J32L*PDstandardNth13gzz - J23L*J31L*PDstandardNth23gyz) + 
      gu31*(dJ133L*PDstandardNth1gxy - J21L*J33L*PDstandardNth23gyz) + 
      gu31*J22L*J31L*PDstandardNth23gzz + gu31*J21L*J32L*PDstandardNth23gzz + 
      2*gu32*J22L*J32L*PDstandardNth23gzz + 
      gu11*(J13L*J21L*PDstandardNth12gxy + J12L*J21L*PDstandardNth12gxz + 
      J13L*(-(J12L*PDstandardNth11gxx) + J11L*PDstandardNth11gxy - 
      J22L*PDstandardNth12gxx - J32L*PDstandardNth13gxx) - 
      dJ123L*PDstandardNth1gxx - dJ223L*PDstandardNth2gxx) + 
      dJ213L*gu11*PDstandardNth2gxy + dJ233L*gu31*PDstandardNth2gxy + 
      gu21*(J13L*(J21L*PDstandardNth12gyy - J32L*PDstandardNth13gxy) - 
      dJ123L*PDstandardNth1gxy - dJ223L*PDstandardNth2gxy) + 
      dJ212L*gu11*PDstandardNth2gxz + dJ222L*gu21*PDstandardNth2gxz + 
      gu31*(2*J13L*J33L*PDstandardNth13gxy - dJ123L*PDstandardNth1gxz - 
      dJ223L*PDstandardNth2gxz) + dJ213L*gu21*PDstandardNth2gyy + 
      dJ233L*gu32*PDstandardNth2gyy - 2*dJ223L*gu32*PDstandardNth2gyz + 
      gu11*(dJ113L*PDstandardNth1gxy - dJ211L*PDstandardNth2gyz) + 
      gu21*(dJ122L*PDstandardNth1gxz - dJ212L*PDstandardNth2gyz) + 
      gu31*((-(J13L*J22L) - J12L*J23L)*PDstandardNth12gxz - 
      J13L*J21L*PDstandardNth12gyz - J12L*J33L*PDstandardNth13gxz - 
      dJ113L*PDstandardNth1gyz - dJ213L*PDstandardNth2gyz) + 
      dJ212L*gu31*PDstandardNth2gzz + dJ222L*gu32*PDstandardNth2gzz + 
      gu11*(dJ112L*PDstandardNth1gxz - J32L*J33L*PDstandardNth33gxx) + 
      gu11*J31L*J33L*PDstandardNth33gxy + gu21*(dJ113L*PDstandardNth1gyy - 
      J32L*J33L*PDstandardNth33gxy) + gu11*J31L*J32L*PDstandardNth33gxz + 
      gu31*(dJ112L*PDstandardNth1gzz - J32L*J33L*PDstandardNth33gxz) + 
      gu21*J31L*J33L*PDstandardNth33gyy - 2*gu32*J32L*J33L*PDstandardNth33gyz 
      + gu21*(-(dJ112L*PDstandardNth1gyz) + J32L*(2*J12L*PDstandardNth13gxz - 
      J31L*PDstandardNth33gyz)) + gu31*(J22L*(J11L*PDstandardNth12gzz - 
      J23L*PDstandardNth22gxz) + J31L*(-(J13L*PDstandardNth13gyz) - 
      J33L*PDstandardNth33gyz)) + gu31*J31L*J32L*PDstandardNth33gzz + 
      gu11*(J13L*J31L*PDstandardNth13gxy - dJ111L*PDstandardNth1gyz - 
      dJ323L*PDstandardNth3gxx) + dJ313L*gu11*PDstandardNth3gxy + 
      dJ333L*gu31*PDstandardNth3gxy + gu21*(J21L*J23L*PDstandardNth22gyy - 
      dJ323L*PDstandardNth3gxy) + dJ312L*gu11*PDstandardNth3gxz + 
      dJ322L*gu21*PDstandardNth3gxz + gu31*(J12L*J31L*PDstandardNth13gzz - 
      J21L*J23L*PDstandardNth22gyz - dJ323L*PDstandardNth3gxz) + 
      dJ313L*gu21*PDstandardNth3gyy + dJ333L*gu32*PDstandardNth3gyy - 
      2*dJ323L*gu32*PDstandardNth3gyz + gu11*(J21L*J23L*PDstandardNth22gxy - 
      dJ311L*PDstandardNth3gyz) + gu21*(J13L*J31L*PDstandardNth13gyy - 
      J21L*J22L*PDstandardNth22gyz - dJ312L*PDstandardNth3gyz) + 
      gu31*(J21L*J22L*PDstandardNth22gzz - dJ313L*PDstandardNth3gyz) + 
      dJ312L*gu31*PDstandardNth3gzz + dJ322L*gu32*PDstandardNth3gzz + 
      gu11*(J21L*J22L*PDstandardNth22gxz - PDstandardNth11gyz*SQR(J11L)) + 
      gu21*PDstandardNth11gxz*SQR(J12L) + gu32*PDstandardNth11gzz*SQR(J12L) + 
      gu31*PDstandardNth11gxy*SQR(J13L) + gu32*PDstandardNth11gyy*SQR(J13L) + 
      PDstandardNth22gyz*(-2*gu32*J22L*J23L - gu11*SQR(J21L)) + 
      gu21*PDstandardNth22gxz*SQR(J22L) + gu32*PDstandardNth22gzz*SQR(J22L) + 
      gu31*PDstandardNth22gxy*SQR(J23L) + gu32*PDstandardNth22gyy*SQR(J23L) + 
      gu11*(J12L*J31L*PDstandardNth13gxz - J23L*J32L*PDstandardNth23gxx - 
      PDstandardNth33gyz*SQR(J31L)) + gu21*PDstandardNth33gxz*SQR(J32L) + 
      gu32*PDstandardNth33gzz*SQR(J32L) + gu31*PDstandardNth33gxy*SQR(J33L) + 
      gu32*PDstandardNth33gyy*SQR(J33L));
    
    CCTK_REAL R33 = 0.5*(gu21*(-4*J13L*J23L*PDstandardNth12gxy + 
      2*J13L*J22L*PDstandardNth12gxz + 2*J13L*J21L*PDstandardNth12gyz - 
      2*J12L*J21L*PDstandardNth12gzz - 4*J13L*J33L*PDstandardNth13gxy + 
      2*J13L*J32L*PDstandardNth13gxz + 2*J13L*J31L*PDstandardNth13gyz - 
      2*J12L*J31L*PDstandardNth13gzz - 2*dJ133L*PDstandardNth1gxy) + 
      2*dJ123L*gu21*PDstandardNth1gxz + 2*(gu11*J11L*J13L*PDstandardNth11gxz 
      + gu21*J12L*J13L*PDstandardNth11gxz + gu21*J11L*J13L*PDstandardNth11gyz 
      + gu22*J12L*J13L*PDstandardNth11gyz + gu11*J13L*J21L*PDstandardNth12gxz 
      + gu11*J11L*J23L*PDstandardNth12gxz + gu21*J12L*J23L*PDstandardNth12gxz 
      + gu22*J13L*J22L*PDstandardNth12gyz + gu21*J11L*J23L*PDstandardNth12gyz 
      + gu22*J12L*J23L*PDstandardNth12gyz + gu11*J13L*J31L*PDstandardNth13gxz 
      + gu11*J11L*J33L*PDstandardNth13gxz + gu21*J12L*J33L*PDstandardNth13gxz 
      + gu22*J13L*J32L*PDstandardNth13gyz + gu21*J11L*J33L*PDstandardNth13gyz 
      + gu22*J12L*J33L*PDstandardNth13gyz + dJ113L*gu11*PDstandardNth1gxz) + 
      dJ133L*(-(gu11*PDstandardNth1gxx) - gu22*PDstandardNth1gyy) + 
      2*dJ113L*gu21*PDstandardNth1gyz - 2*(G111*G133 + G122*G133 + G211*G233 
      + G222*G233 + G311*G333 + G322*G333 + gu21*J11L*J12L*PDstandardNth11gzz 
      + gu11*J13L*J23L*PDstandardNth12gxx + gu22*J13L*J23L*PDstandardNth12gyy 
      + gu11*J11L*J21L*PDstandardNth12gzz + gu21*J11L*J22L*PDstandardNth12gzz 
      + gu22*J12L*J22L*PDstandardNth12gzz + gu11*J13L*J33L*PDstandardNth13gxx 
      + gu22*J13L*J33L*PDstandardNth13gyy + gu11*J11L*J31L*PDstandardNth13gzz 
      + gu21*J11L*J32L*PDstandardNth13gzz + gu22*J12L*J32L*PDstandardNth13gzz 
      + dJ112L*gu21*PDstandardNth1gzz) + 2*gu21*J22L*J23L*PDstandardNth22gxz 
      + 2*gu21*J21L*J23L*PDstandardNth22gyz - 
      4*gu21*J23L*J33L*PDstandardNth23gxy + 
      2*gu21*J23L*J32L*PDstandardNth23gxz + 
      2*gu21*J22L*J33L*PDstandardNth23gxz + 
      2*gu21*J23L*J31L*PDstandardNth23gyz + 
      2*gu21*J21L*J33L*PDstandardNth23gyz - 
      2*gu21*J22L*J31L*PDstandardNth23gzz - 
      2*gu21*J21L*J32L*PDstandardNth23gzz - 2*dJ233L*gu21*PDstandardNth2gxy + 
      2*dJ223L*gu21*PDstandardNth2gxz + 2*dJ213L*gu21*PDstandardNth2gyz - 
      2*dJ212L*gu21*PDstandardNth2gzz + 2*gu21*J32L*J33L*PDstandardNth33gxz + 
      2*gu21*J31L*J33L*PDstandardNth33gyz - 
      2*gu21*J31L*J32L*PDstandardNth33gzz + 
      gu11*(2*J23L*J31L*PDstandardNth23gxz - dJ333L*PDstandardNth3gxx) - 
      2*dJ333L*gu21*PDstandardNth3gxy + 2*dJ313L*gu11*PDstandardNth3gxz + 
      2*dJ323L*gu21*PDstandardNth3gxz + gu22*(-2*J23L*J33L*PDstandardNth23gyy 
      - dJ333L*PDstandardNth3gyy) + 2*dJ313L*gu21*PDstandardNth3gyz + 
      2*dJ323L*gu22*PDstandardNth3gyz - 2*dJ312L*gu21*PDstandardNth3gzz + 
      gu11*(2*J21L*J33L*PDstandardNth23gxz - dJ311L*PDstandardNth3gzz) + 
      gu22*(2*J23L*J32L*PDstandardNth23gyz - dJ322L*PDstandardNth3gzz) + 
      2*SQR(G113) + 2*SQR(G123) + 2*SQR(G213) + 2*SQR(G223) + 2*SQR(G313) + 
      2*SQR(G323) + gu11*(-2*J21L*J31L*PDstandardNth23gzz - 
      PDstandardNth11gzz*SQR(J11L)) + gu22*(2*J22L*J33L*PDstandardNth23gyz - 
      PDstandardNth11gzz*SQR(J12L)) - 2*gu21*PDstandardNth11gxy*SQR(J13L) + 
      gu11*(-(dJ111L*PDstandardNth1gzz) + 2*J21L*J23L*PDstandardNth22gxz - 
      dJ233L*PDstandardNth2gxx - PDstandardNth11gxx*SQR(J13L)) + 
      gu22*(-2*J22L*J32L*PDstandardNth23gzz - PDstandardNth11gyy*SQR(J13L)) + 
      PDstandardNth22gzz*(-2*gu21*J21L*J22L - gu11*SQR(J21L)) + 
      gu22*(2*dJ123L*PDstandardNth1gyz - dJ122L*PDstandardNth1gzz - 
      dJ233L*PDstandardNth2gyy - PDstandardNth22gzz*SQR(J22L)) - 
      2*gu21*PDstandardNth22gxy*SQR(J23L) + gu11*(2*dJ213L*PDstandardNth2gxz 
      - PDstandardNth22gxx*SQR(J23L)) + gu22*(2*dJ223L*PDstandardNth2gyz - 
      PDstandardNth22gyy*SQR(J23L)) + gu11*(-2*J23L*J33L*PDstandardNth23gxx - 
      dJ211L*PDstandardNth2gzz - PDstandardNth33gzz*SQR(J31L)) + 
      gu22*(2*J22L*J23L*PDstandardNth22gyz - dJ222L*PDstandardNth2gzz - 
      PDstandardNth33gzz*SQR(J32L)) - 2*gu21*PDstandardNth33gxy*SQR(J33L) + 
      gu11*(2*J31L*J33L*PDstandardNth33gxz - PDstandardNth33gxx*SQR(J33L)) + 
      gu22*(2*J32L*J33L*PDstandardNth33gyz - PDstandardNth33gyy*SQR(J33L)));
    
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
      G211*kyyL + G311*kyzL + J12L*PDstandardNth1kxx - J11L*PDstandardNth1kxy 
      + J22L*PDstandardNth2kxx - J21L*PDstandardNth2kxy + 
      J32L*PDstandardNth3kxx - J31L*PDstandardNth3kxy) + gu31*(-(G113*kxxL) - 
      G213*kxyL + (G111 - G313)*kxzL + G211*kyzL + G311*kzzL + 
      J13L*PDstandardNth1kxx - J11L*PDstandardNth1kxz + 
      J23L*PDstandardNth2kxx - J21L*PDstandardNth2kxz + 
      J33L*PDstandardNth3kxx - J31L*PDstandardNth3kxz) + gu22*(-(G122*kxxL) + 
      (G112 - G222)*kxyL - G322*kxzL + G212*kyyL + G312*kyzL + 
      J12L*PDstandardNth1kxy - J11L*PDstandardNth1kyy + 
      J22L*PDstandardNth2kxy - J21L*PDstandardNth2kyy + 
      J32L*PDstandardNth3kxy - J31L*PDstandardNth3kyy) + gu32*(G113*kxyL + 
      G112*kxzL + G213*kyyL + (G212 + G313)*kyzL + G312*kzzL + 
      J13L*PDstandardNth1kxy + J12L*PDstandardNth1kxz + 
      J23L*PDstandardNth2kxy + J22L*PDstandardNth2kxz + 
      J33L*PDstandardNth3kxy + J32L*PDstandardNth3kxz - 2*(G123*kxxL + 
      G223*kxyL + G323*kxzL + J11L*PDstandardNth1kyz + J21L*PDstandardNth2kyz 
      + J31L*PDstandardNth3kyz)) + gu33*(-(G133*kxxL) - G233*kxyL + (G113 - 
      G333)*kxzL + G213*kyzL + G313*kzzL + J13L*PDstandardNth1kxz - 
      J11L*PDstandardNth1kzz + J23L*PDstandardNth2kxz - 
      J21L*PDstandardNth2kzz + J33L*PDstandardNth3kxz - 
      J31L*PDstandardNth3kzz) - 25.13274122871834590770114706623602307358*S1;
    
    CCTK_REAL M2L = gu11*(G112*kxxL + (-G111 + G212)*kxyL + G312*kxzL - 
      G211*kyyL - G311*kyzL - J12L*PDstandardNth1kxx + J11L*PDstandardNth1kxy 
      - J22L*PDstandardNth2kxx + J21L*PDstandardNth2kxy - 
      J32L*PDstandardNth3kxx + J31L*PDstandardNth3kxy) + gu21*(G122*kxxL + 
      (-G112 + G222)*kxyL + G322*kxzL - G212*kyyL - G312*kyzL - 
      J12L*PDstandardNth1kxy + J11L*PDstandardNth1kyy - 
      J22L*PDstandardNth2kxy + J21L*PDstandardNth2kyy - 
      J32L*PDstandardNth3kxy + J31L*PDstandardNth3kyy) + gu31*(G123*kxxL + 
      (-2*G113 + G223)*kxyL + (G112 + G323)*kxzL + G212*kyzL + G312*kzzL + 
      J13L*PDstandardNth1kxy + J11L*PDstandardNth1kyz + 
      J23L*PDstandardNth2kxy + J21L*PDstandardNth2kyz + 
      J33L*PDstandardNth3kxy - 2*(G213*kyyL + G313*kyzL + 
      J12L*PDstandardNth1kxz + J22L*PDstandardNth2kxz + 
      J32L*PDstandardNth3kxz) + J31L*PDstandardNth3kyz) + gu32*(-(G123*kxyL) 
      + G122*kxzL - G223*kyyL + (G222 - G323)*kyzL + G322*kzzL + 
      J13L*PDstandardNth1kyy - J12L*PDstandardNth1kyz + 
      J23L*PDstandardNth2kyy - J22L*PDstandardNth2kyz + 
      J33L*PDstandardNth3kyy - J32L*PDstandardNth3kyz) + gu33*(-(G133*kxyL) + 
      G123*kxzL - G233*kyyL + (G223 - G333)*kyzL + G323*kzzL + 
      J13L*PDstandardNth1kyz - J12L*PDstandardNth1kzz + 
      J23L*PDstandardNth2kyz - J22L*PDstandardNth2kzz + 
      J33L*PDstandardNth3kyz - J32L*PDstandardNth3kzz) - 
      25.13274122871834590770114706623602307358*S2;
    
    CCTK_REAL M3L = gu11*(G113*kxxL + G213*kxyL + (-G111 + G313)*kxzL - 
      G211*kyzL - G311*kzzL - J13L*PDstandardNth1kxx + J11L*PDstandardNth1kxz 
      - J23L*PDstandardNth2kxx + J21L*PDstandardNth2kxz - 
      J33L*PDstandardNth3kxx + J31L*PDstandardNth3kxz) + gu21*(G123*kxxL + 
      (G113 + G223)*kxyL + (-2*G112 + G323)*kxzL + G213*kyyL + (-2*G212 + 
      G313)*kyzL + J12L*PDstandardNth1kxz + J11L*PDstandardNth1kyz + 
      J22L*PDstandardNth2kxz + J21L*PDstandardNth2kyz - 2*(G312*kzzL + 
      J13L*PDstandardNth1kxy + J23L*PDstandardNth2kxy + 
      J33L*PDstandardNth3kxy) + J32L*PDstandardNth3kxz + 
      J31L*PDstandardNth3kyz) + gu22*(G123*kxyL - G122*kxzL + G223*kyyL + 
      (-G222 + G323)*kyzL - G322*kzzL - J13L*PDstandardNth1kyy + 
      J12L*PDstandardNth1kyz - J23L*PDstandardNth2kyy + 
      J22L*PDstandardNth2kyz - J33L*PDstandardNth3kyy + 
      J32L*PDstandardNth3kyz) + gu31*(G133*kxxL + G233*kxyL + (-G113 + 
      G333)*kxzL - G213*kyzL - G313*kzzL - J13L*PDstandardNth1kxz + 
      J11L*PDstandardNth1kzz - J23L*PDstandardNth2kxz + 
      J21L*PDstandardNth2kzz - J33L*PDstandardNth3kxz + 
      J31L*PDstandardNth3kzz) + gu32*(G133*kxyL - G123*kxzL + G233*kyyL + 
      (-G223 + G333)*kyzL - G323*kzzL - J13L*PDstandardNth1kyz + 
      J12L*PDstandardNth1kzz - J23L*PDstandardNth2kyz + 
      J22L*PDstandardNth2kzz - J33L*PDstandardNth3kyz + 
      J32L*PDstandardNth3kzz) - 25.13274122871834590770114706623602307358*S3;
    
    /* Copy local copies back to grid functions */
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
  }
  LC_ENDLOOP3 (ML_ADMConstraints_MP);
}

extern "C" void ML_ADMConstraints_MP(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADMConstraints_MP_Body);
}
