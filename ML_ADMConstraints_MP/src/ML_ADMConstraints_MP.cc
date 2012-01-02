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

static void ML_ADMConstraints_MP_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  
  /* Jacobian variable pointers */
  bool const use_jacobian = (!CCTK_IsFunctionAliased("MultiPatch_GetMap") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)
                       && strlen(jacobian_group) > 0;
  if (use_jacobian && strlen(jacobian_derivative_group) == 0)
  {
    CCTK_WARN (1, "GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names");
  }
  
  CCTK_REAL const *restrict jacobian_ptrs[9];
  if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_group,
                                                9, jacobian_ptrs);
  
  CCTK_REAL const *restrict const J11 = use_jacobian ? jacobian_ptrs[0] : 0;
  CCTK_REAL const *restrict const J12 = use_jacobian ? jacobian_ptrs[1] : 0;
  CCTK_REAL const *restrict const J13 = use_jacobian ? jacobian_ptrs[2] : 0;
  CCTK_REAL const *restrict const J21 = use_jacobian ? jacobian_ptrs[3] : 0;
  CCTK_REAL const *restrict const J22 = use_jacobian ? jacobian_ptrs[4] : 0;
  CCTK_REAL const *restrict const J23 = use_jacobian ? jacobian_ptrs[5] : 0;
  CCTK_REAL const *restrict const J31 = use_jacobian ? jacobian_ptrs[6] : 0;
  CCTK_REAL const *restrict const J32 = use_jacobian ? jacobian_ptrs[7] : 0;
  CCTK_REAL const *restrict const J33 = use_jacobian ? jacobian_ptrs[8] : 0;
  
  CCTK_REAL const *restrict jacobian_derivative_ptrs[18];
  if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_derivative_group,
                                                18, jacobian_derivative_ptrs);
  
  CCTK_REAL const *restrict const dJ111 = use_jacobian ? jacobian_derivative_ptrs[0] : 0;
  CCTK_REAL const *restrict const dJ112 = use_jacobian ? jacobian_derivative_ptrs[1] : 0;
  CCTK_REAL const *restrict const dJ113 = use_jacobian ? jacobian_derivative_ptrs[2] : 0;
  CCTK_REAL const *restrict const dJ122 = use_jacobian ? jacobian_derivative_ptrs[3] : 0;
  CCTK_REAL const *restrict const dJ123 = use_jacobian ? jacobian_derivative_ptrs[4] : 0;
  CCTK_REAL const *restrict const dJ133 = use_jacobian ? jacobian_derivative_ptrs[5] : 0;
  CCTK_REAL const *restrict const dJ211 = use_jacobian ? jacobian_derivative_ptrs[6] : 0;
  CCTK_REAL const *restrict const dJ212 = use_jacobian ? jacobian_derivative_ptrs[7] : 0;
  CCTK_REAL const *restrict const dJ213 = use_jacobian ? jacobian_derivative_ptrs[8] : 0;
  CCTK_REAL const *restrict const dJ222 = use_jacobian ? jacobian_derivative_ptrs[9] : 0;
  CCTK_REAL const *restrict const dJ223 = use_jacobian ? jacobian_derivative_ptrs[10] : 0;
  CCTK_REAL const *restrict const dJ233 = use_jacobian ? jacobian_derivative_ptrs[11] : 0;
  CCTK_REAL const *restrict const dJ311 = use_jacobian ? jacobian_derivative_ptrs[12] : 0;
  CCTK_REAL const *restrict const dJ312 = use_jacobian ? jacobian_derivative_ptrs[13] : 0;
  CCTK_REAL const *restrict const dJ313 = use_jacobian ? jacobian_derivative_ptrs[14] : 0;
  CCTK_REAL const *restrict const dJ322 = use_jacobian ? jacobian_derivative_ptrs[15] : 0;
  CCTK_REAL const *restrict const dJ323 = use_jacobian ? jacobian_derivative_ptrs[16] : 0;
  CCTK_REAL const *restrict const dJ333 = use_jacobian ? jacobian_derivative_ptrs[17] : 0;
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_ADMConstraints_MP,
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
    
    CCTK_REAL dJ111L, dJ112L, dJ113L, dJ122L, dJ123L, dJ133L, dJ211L, dJ212L, dJ213L, dJ222L, dJ223L, dJ233L, dJ311L, dJ312L, dJ313L, dJ322L, dJ323L, dJ333L, J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L;
    
    if (use_jacobian)
    {
      dJ111L = dJ111[index];
      dJ112L = dJ112[index];
      dJ113L = dJ113[index];
      dJ122L = dJ122[index];
      dJ123L = dJ123[index];
      dJ133L = dJ133[index];
      dJ211L = dJ211[index];
      dJ212L = dJ212[index];
      dJ213L = dJ213[index];
      dJ222L = dJ222[index];
      dJ223L = dJ223[index];
      dJ233L = dJ233[index];
      dJ311L = dJ311[index];
      dJ312L = dJ312[index];
      dJ313L = dJ313[index];
      dJ322L = dJ322[index];
      dJ323L = dJ323[index];
      dJ333L = dJ333[index];
      J11L = J11[index];
      J12L = J12[index];
      J13L = J13[index];
      J21L = J21[index];
      J22L = J22[index];
      J23L = J23[index];
      J31L = J31[index];
      J32L = J32[index];
      J33L = J33[index];
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
    CCTK_REAL JacPDstandardNth11gyy;
    CCTK_REAL JacPDstandardNth11gyz;
    CCTK_REAL JacPDstandardNth11gzz;
    CCTK_REAL JacPDstandardNth12gxx;
    CCTK_REAL JacPDstandardNth12gxy;
    CCTK_REAL JacPDstandardNth12gxz;
    CCTK_REAL JacPDstandardNth12gyy;
    CCTK_REAL JacPDstandardNth12gyz;
    CCTK_REAL JacPDstandardNth12gzz;
    CCTK_REAL JacPDstandardNth13gxx;
    CCTK_REAL JacPDstandardNth13gxy;
    CCTK_REAL JacPDstandardNth13gxz;
    CCTK_REAL JacPDstandardNth13gyy;
    CCTK_REAL JacPDstandardNth13gyz;
    CCTK_REAL JacPDstandardNth13gzz;
    CCTK_REAL JacPDstandardNth1gxx;
    CCTK_REAL JacPDstandardNth1gxy;
    CCTK_REAL JacPDstandardNth1gxz;
    CCTK_REAL JacPDstandardNth1gyy;
    CCTK_REAL JacPDstandardNth1gyz;
    CCTK_REAL JacPDstandardNth1gzz;
    CCTK_REAL JacPDstandardNth1kxy;
    CCTK_REAL JacPDstandardNth1kxz;
    CCTK_REAL JacPDstandardNth1kyy;
    CCTK_REAL JacPDstandardNth1kyz;
    CCTK_REAL JacPDstandardNth1kzz;
    CCTK_REAL JacPDstandardNth21gxx;
    CCTK_REAL JacPDstandardNth21gxy;
    CCTK_REAL JacPDstandardNth21gxz;
    CCTK_REAL JacPDstandardNth21gyy;
    CCTK_REAL JacPDstandardNth21gyz;
    CCTK_REAL JacPDstandardNth21gzz;
    CCTK_REAL JacPDstandardNth22gxx;
    CCTK_REAL JacPDstandardNth22gxz;
    CCTK_REAL JacPDstandardNth22gzz;
    CCTK_REAL JacPDstandardNth23gxx;
    CCTK_REAL JacPDstandardNth23gxy;
    CCTK_REAL JacPDstandardNth23gxz;
    CCTK_REAL JacPDstandardNth23gyy;
    CCTK_REAL JacPDstandardNth23gyz;
    CCTK_REAL JacPDstandardNth23gzz;
    CCTK_REAL JacPDstandardNth2gxx;
    CCTK_REAL JacPDstandardNth2gxy;
    CCTK_REAL JacPDstandardNth2gxz;
    CCTK_REAL JacPDstandardNth2gyy;
    CCTK_REAL JacPDstandardNth2gyz;
    CCTK_REAL JacPDstandardNth2gzz;
    CCTK_REAL JacPDstandardNth2kxx;
    CCTK_REAL JacPDstandardNth2kxy;
    CCTK_REAL JacPDstandardNth2kxz;
    CCTK_REAL JacPDstandardNth2kyz;
    CCTK_REAL JacPDstandardNth2kzz;
    CCTK_REAL JacPDstandardNth31gxx;
    CCTK_REAL JacPDstandardNth31gxy;
    CCTK_REAL JacPDstandardNth31gxz;
    CCTK_REAL JacPDstandardNth31gyy;
    CCTK_REAL JacPDstandardNth31gyz;
    CCTK_REAL JacPDstandardNth31gzz;
    CCTK_REAL JacPDstandardNth32gxx;
    CCTK_REAL JacPDstandardNth32gxy;
    CCTK_REAL JacPDstandardNth32gxz;
    CCTK_REAL JacPDstandardNth32gyy;
    CCTK_REAL JacPDstandardNth32gyz;
    CCTK_REAL JacPDstandardNth32gzz;
    CCTK_REAL JacPDstandardNth33gxx;
    CCTK_REAL JacPDstandardNth33gxy;
    CCTK_REAL JacPDstandardNth33gyy;
    CCTK_REAL JacPDstandardNth3gxx;
    CCTK_REAL JacPDstandardNth3gxy;
    CCTK_REAL JacPDstandardNth3gxz;
    CCTK_REAL JacPDstandardNth3gyy;
    CCTK_REAL JacPDstandardNth3gyz;
    CCTK_REAL JacPDstandardNth3gzz;
    CCTK_REAL JacPDstandardNth3kxx;
    CCTK_REAL JacPDstandardNth3kxy;
    CCTK_REAL JacPDstandardNth3kxz;
    CCTK_REAL JacPDstandardNth3kyy;
    CCTK_REAL JacPDstandardNth3kyz;
    
    if (use_jacobian)
    {
      JacPDstandardNth1gxx = J11L*PDstandardNth1gxx + 
        J21L*PDstandardNth2gxx + J31L*PDstandardNth3gxx;
      
      JacPDstandardNth1gxy = J11L*PDstandardNth1gxy + 
        J21L*PDstandardNth2gxy + J31L*PDstandardNth3gxy;
      
      JacPDstandardNth1gxz = J11L*PDstandardNth1gxz + 
        J21L*PDstandardNth2gxz + J31L*PDstandardNth3gxz;
      
      JacPDstandardNth1gyy = J11L*PDstandardNth1gyy + 
        J21L*PDstandardNth2gyy + J31L*PDstandardNth3gyy;
      
      JacPDstandardNth1gyz = J11L*PDstandardNth1gyz + 
        J21L*PDstandardNth2gyz + J31L*PDstandardNth3gyz;
      
      JacPDstandardNth1gzz = J11L*PDstandardNth1gzz + 
        J21L*PDstandardNth2gzz + J31L*PDstandardNth3gzz;
      
      JacPDstandardNth1kxy = J11L*PDstandardNth1kxy + 
        J21L*PDstandardNth2kxy + J31L*PDstandardNth3kxy;
      
      JacPDstandardNth1kxz = J11L*PDstandardNth1kxz + 
        J21L*PDstandardNth2kxz + J31L*PDstandardNth3kxz;
      
      JacPDstandardNth1kyy = J11L*PDstandardNth1kyy + 
        J21L*PDstandardNth2kyy + J31L*PDstandardNth3kyy;
      
      JacPDstandardNth1kyz = J11L*PDstandardNth1kyz + 
        J21L*PDstandardNth2kyz + J31L*PDstandardNth3kyz;
      
      JacPDstandardNth1kzz = J11L*PDstandardNth1kzz + 
        J21L*PDstandardNth2kzz + J31L*PDstandardNth3kzz;
      
      JacPDstandardNth2gxx = J12L*PDstandardNth1gxx + 
        J22L*PDstandardNth2gxx + J32L*PDstandardNth3gxx;
      
      JacPDstandardNth2gxy = J12L*PDstandardNth1gxy + 
        J22L*PDstandardNth2gxy + J32L*PDstandardNth3gxy;
      
      JacPDstandardNth2gxz = J12L*PDstandardNth1gxz + 
        J22L*PDstandardNth2gxz + J32L*PDstandardNth3gxz;
      
      JacPDstandardNth2gyy = J12L*PDstandardNth1gyy + 
        J22L*PDstandardNth2gyy + J32L*PDstandardNth3gyy;
      
      JacPDstandardNth2gyz = J12L*PDstandardNth1gyz + 
        J22L*PDstandardNth2gyz + J32L*PDstandardNth3gyz;
      
      JacPDstandardNth2gzz = J12L*PDstandardNth1gzz + 
        J22L*PDstandardNth2gzz + J32L*PDstandardNth3gzz;
      
      JacPDstandardNth2kxx = J12L*PDstandardNth1kxx + 
        J22L*PDstandardNth2kxx + J32L*PDstandardNth3kxx;
      
      JacPDstandardNth2kxy = J12L*PDstandardNth1kxy + 
        J22L*PDstandardNth2kxy + J32L*PDstandardNth3kxy;
      
      JacPDstandardNth2kxz = J12L*PDstandardNth1kxz + 
        J22L*PDstandardNth2kxz + J32L*PDstandardNth3kxz;
      
      JacPDstandardNth2kyz = J12L*PDstandardNth1kyz + 
        J22L*PDstandardNth2kyz + J32L*PDstandardNth3kyz;
      
      JacPDstandardNth2kzz = J12L*PDstandardNth1kzz + 
        J22L*PDstandardNth2kzz + J32L*PDstandardNth3kzz;
      
      JacPDstandardNth3gxx = J13L*PDstandardNth1gxx + 
        J23L*PDstandardNth2gxx + J33L*PDstandardNth3gxx;
      
      JacPDstandardNth3gxy = J13L*PDstandardNth1gxy + 
        J23L*PDstandardNth2gxy + J33L*PDstandardNth3gxy;
      
      JacPDstandardNth3gxz = J13L*PDstandardNth1gxz + 
        J23L*PDstandardNth2gxz + J33L*PDstandardNth3gxz;
      
      JacPDstandardNth3gyy = J13L*PDstandardNth1gyy + 
        J23L*PDstandardNth2gyy + J33L*PDstandardNth3gyy;
      
      JacPDstandardNth3gyz = J13L*PDstandardNth1gyz + 
        J23L*PDstandardNth2gyz + J33L*PDstandardNth3gyz;
      
      JacPDstandardNth3gzz = J13L*PDstandardNth1gzz + 
        J23L*PDstandardNth2gzz + J33L*PDstandardNth3gzz;
      
      JacPDstandardNth3kxx = J13L*PDstandardNth1kxx + 
        J23L*PDstandardNth2kxx + J33L*PDstandardNth3kxx;
      
      JacPDstandardNth3kxy = J13L*PDstandardNth1kxy + 
        J23L*PDstandardNth2kxy + J33L*PDstandardNth3kxy;
      
      JacPDstandardNth3kxz = J13L*PDstandardNth1kxz + 
        J23L*PDstandardNth2kxz + J33L*PDstandardNth3kxz;
      
      JacPDstandardNth3kyy = J13L*PDstandardNth1kyy + 
        J23L*PDstandardNth2kyy + J33L*PDstandardNth3kyy;
      
      JacPDstandardNth3kyz = J13L*PDstandardNth1kyz + 
        J23L*PDstandardNth2kyz + J33L*PDstandardNth3kyz;
      
      JacPDstandardNth11gyy = dJ111L*PDstandardNth1gyy + 
        2*(J11L*(J21L*PDstandardNth12gyy + J31L*PDstandardNth13gyy) + 
        J21L*J31L*PDstandardNth23gyy) + dJ211L*PDstandardNth2gyy + 
        dJ311L*PDstandardNth3gyy + PDstandardNth11gyy*SQR(J11L) + 
        PDstandardNth22gyy*SQR(J21L) + PDstandardNth33gyy*SQR(J31L);
      
      JacPDstandardNth11gyz = dJ111L*PDstandardNth1gyz + 
        2*(J11L*(J21L*PDstandardNth12gyz + J31L*PDstandardNth13gyz) + 
        J21L*J31L*PDstandardNth23gyz) + dJ211L*PDstandardNth2gyz + 
        dJ311L*PDstandardNth3gyz + PDstandardNth11gyz*SQR(J11L) + 
        PDstandardNth22gyz*SQR(J21L) + PDstandardNth33gyz*SQR(J31L);
      
      JacPDstandardNth11gzz = dJ111L*PDstandardNth1gzz + 
        2*(J11L*(J21L*PDstandardNth12gzz + J31L*PDstandardNth13gzz) + 
        J21L*J31L*PDstandardNth23gzz) + dJ211L*PDstandardNth2gzz + 
        dJ311L*PDstandardNth3gzz + PDstandardNth11gzz*SQR(J11L) + 
        PDstandardNth22gzz*SQR(J21L) + PDstandardNth33gzz*SQR(J31L);
      
      JacPDstandardNth22gxx = dJ122L*PDstandardNth1gxx + 
        2*(J12L*(J22L*PDstandardNth12gxx + J32L*PDstandardNth13gxx) + 
        J22L*J32L*PDstandardNth23gxx) + dJ222L*PDstandardNth2gxx + 
        dJ322L*PDstandardNth3gxx + PDstandardNth11gxx*SQR(J12L) + 
        PDstandardNth22gxx*SQR(J22L) + PDstandardNth33gxx*SQR(J32L);
      
      JacPDstandardNth22gxz = dJ122L*PDstandardNth1gxz + 
        2*(J12L*(J22L*PDstandardNth12gxz + J32L*PDstandardNth13gxz) + 
        J22L*J32L*PDstandardNth23gxz) + dJ222L*PDstandardNth2gxz + 
        dJ322L*PDstandardNth3gxz + PDstandardNth11gxz*SQR(J12L) + 
        PDstandardNth22gxz*SQR(J22L) + PDstandardNth33gxz*SQR(J32L);
      
      JacPDstandardNth22gzz = dJ122L*PDstandardNth1gzz + 
        2*(J12L*(J22L*PDstandardNth12gzz + J32L*PDstandardNth13gzz) + 
        J22L*J32L*PDstandardNth23gzz) + dJ222L*PDstandardNth2gzz + 
        dJ322L*PDstandardNth3gzz + PDstandardNth11gzz*SQR(J12L) + 
        PDstandardNth22gzz*SQR(J22L) + PDstandardNth33gzz*SQR(J32L);
      
      JacPDstandardNth33gxx = dJ133L*PDstandardNth1gxx + 
        2*(J13L*(J23L*PDstandardNth12gxx + J33L*PDstandardNth13gxx) + 
        J23L*J33L*PDstandardNth23gxx) + dJ233L*PDstandardNth2gxx + 
        dJ333L*PDstandardNth3gxx + PDstandardNth11gxx*SQR(J13L) + 
        PDstandardNth22gxx*SQR(J23L) + PDstandardNth33gxx*SQR(J33L);
      
      JacPDstandardNth33gxy = dJ133L*PDstandardNth1gxy + 
        2*(J13L*(J23L*PDstandardNth12gxy + J33L*PDstandardNth13gxy) + 
        J23L*J33L*PDstandardNth23gxy) + dJ233L*PDstandardNth2gxy + 
        dJ333L*PDstandardNth3gxy + PDstandardNth11gxy*SQR(J13L) + 
        PDstandardNth22gxy*SQR(J23L) + PDstandardNth33gxy*SQR(J33L);
      
      JacPDstandardNth33gyy = dJ133L*PDstandardNth1gyy + 
        2*(J13L*(J23L*PDstandardNth12gyy + J33L*PDstandardNth13gyy) + 
        J23L*J33L*PDstandardNth23gyy) + dJ233L*PDstandardNth2gyy + 
        dJ333L*PDstandardNth3gyy + PDstandardNth11gyy*SQR(J13L) + 
        PDstandardNth22gyy*SQR(J23L) + PDstandardNth33gyy*SQR(J33L);
      
      JacPDstandardNth12gxx = J12L*(J11L*PDstandardNth11gxx + 
        J21L*PDstandardNth12gxx + J31L*PDstandardNth13gxx) + 
        J11L*(J22L*PDstandardNth12gxx + J32L*PDstandardNth13gxx) + 
        dJ112L*PDstandardNth1gxx + J22L*(J21L*PDstandardNth22gxx + 
        J31L*PDstandardNth23gxx) + dJ212L*PDstandardNth2gxx + 
        J32L*(J21L*PDstandardNth23gxx + J31L*PDstandardNth33gxx) + 
        dJ312L*PDstandardNth3gxx;
      
      JacPDstandardNth12gxy = J12L*(J11L*PDstandardNth11gxy + 
        J21L*PDstandardNth12gxy + J31L*PDstandardNth13gxy) + 
        J11L*(J22L*PDstandardNth12gxy + J32L*PDstandardNth13gxy) + 
        dJ112L*PDstandardNth1gxy + J22L*(J21L*PDstandardNth22gxy + 
        J31L*PDstandardNth23gxy) + dJ212L*PDstandardNth2gxy + 
        J32L*(J21L*PDstandardNth23gxy + J31L*PDstandardNth33gxy) + 
        dJ312L*PDstandardNth3gxy;
      
      JacPDstandardNth12gxz = J12L*(J11L*PDstandardNth11gxz + 
        J21L*PDstandardNth12gxz + J31L*PDstandardNth13gxz) + 
        J11L*(J22L*PDstandardNth12gxz + J32L*PDstandardNth13gxz) + 
        dJ112L*PDstandardNth1gxz + J22L*(J21L*PDstandardNth22gxz + 
        J31L*PDstandardNth23gxz) + dJ212L*PDstandardNth2gxz + 
        J32L*(J21L*PDstandardNth23gxz + J31L*PDstandardNth33gxz) + 
        dJ312L*PDstandardNth3gxz;
      
      JacPDstandardNth12gyy = J12L*(J11L*PDstandardNth11gyy + 
        J21L*PDstandardNth12gyy + J31L*PDstandardNth13gyy) + 
        J11L*(J22L*PDstandardNth12gyy + J32L*PDstandardNth13gyy) + 
        dJ112L*PDstandardNth1gyy + J22L*(J21L*PDstandardNth22gyy + 
        J31L*PDstandardNth23gyy) + dJ212L*PDstandardNth2gyy + 
        J32L*(J21L*PDstandardNth23gyy + J31L*PDstandardNth33gyy) + 
        dJ312L*PDstandardNth3gyy;
      
      JacPDstandardNth12gyz = J12L*(J11L*PDstandardNth11gyz + 
        J21L*PDstandardNth12gyz + J31L*PDstandardNth13gyz) + 
        J11L*(J22L*PDstandardNth12gyz + J32L*PDstandardNth13gyz) + 
        dJ112L*PDstandardNth1gyz + J22L*(J21L*PDstandardNth22gyz + 
        J31L*PDstandardNth23gyz) + dJ212L*PDstandardNth2gyz + 
        J32L*(J21L*PDstandardNth23gyz + J31L*PDstandardNth33gyz) + 
        dJ312L*PDstandardNth3gyz;
      
      JacPDstandardNth12gzz = J12L*(J11L*PDstandardNth11gzz + 
        J21L*PDstandardNth12gzz + J31L*PDstandardNth13gzz) + 
        J11L*(J22L*PDstandardNth12gzz + J32L*PDstandardNth13gzz) + 
        dJ112L*PDstandardNth1gzz + J22L*(J21L*PDstandardNth22gzz + 
        J31L*PDstandardNth23gzz) + dJ212L*PDstandardNth2gzz + 
        J32L*(J21L*PDstandardNth23gzz + J31L*PDstandardNth33gzz) + 
        dJ312L*PDstandardNth3gzz;
      
      JacPDstandardNth13gxx = J13L*(J11L*PDstandardNth11gxx + 
        J21L*PDstandardNth12gxx + J31L*PDstandardNth13gxx) + 
        J11L*(J23L*PDstandardNth12gxx + J33L*PDstandardNth13gxx) + 
        dJ113L*PDstandardNth1gxx + J23L*(J21L*PDstandardNth22gxx + 
        J31L*PDstandardNth23gxx) + dJ213L*PDstandardNth2gxx + 
        J33L*(J21L*PDstandardNth23gxx + J31L*PDstandardNth33gxx) + 
        dJ313L*PDstandardNth3gxx;
      
      JacPDstandardNth13gxy = J13L*(J11L*PDstandardNth11gxy + 
        J21L*PDstandardNth12gxy + J31L*PDstandardNth13gxy) + 
        J11L*(J23L*PDstandardNth12gxy + J33L*PDstandardNth13gxy) + 
        dJ113L*PDstandardNth1gxy + J23L*(J21L*PDstandardNth22gxy + 
        J31L*PDstandardNth23gxy) + dJ213L*PDstandardNth2gxy + 
        J33L*(J21L*PDstandardNth23gxy + J31L*PDstandardNth33gxy) + 
        dJ313L*PDstandardNth3gxy;
      
      JacPDstandardNth13gxz = J13L*(J11L*PDstandardNth11gxz + 
        J21L*PDstandardNth12gxz + J31L*PDstandardNth13gxz) + 
        J11L*(J23L*PDstandardNth12gxz + J33L*PDstandardNth13gxz) + 
        dJ113L*PDstandardNth1gxz + J23L*(J21L*PDstandardNth22gxz + 
        J31L*PDstandardNth23gxz) + dJ213L*PDstandardNth2gxz + 
        J33L*(J21L*PDstandardNth23gxz + J31L*PDstandardNth33gxz) + 
        dJ313L*PDstandardNth3gxz;
      
      JacPDstandardNth13gyy = J13L*(J11L*PDstandardNth11gyy + 
        J21L*PDstandardNth12gyy + J31L*PDstandardNth13gyy) + 
        J11L*(J23L*PDstandardNth12gyy + J33L*PDstandardNth13gyy) + 
        dJ113L*PDstandardNth1gyy + J23L*(J21L*PDstandardNth22gyy + 
        J31L*PDstandardNth23gyy) + dJ213L*PDstandardNth2gyy + 
        J33L*(J21L*PDstandardNth23gyy + J31L*PDstandardNth33gyy) + 
        dJ313L*PDstandardNth3gyy;
      
      JacPDstandardNth13gyz = J13L*(J11L*PDstandardNth11gyz + 
        J21L*PDstandardNth12gyz + J31L*PDstandardNth13gyz) + 
        J11L*(J23L*PDstandardNth12gyz + J33L*PDstandardNth13gyz) + 
        dJ113L*PDstandardNth1gyz + J23L*(J21L*PDstandardNth22gyz + 
        J31L*PDstandardNth23gyz) + dJ213L*PDstandardNth2gyz + 
        J33L*(J21L*PDstandardNth23gyz + J31L*PDstandardNth33gyz) + 
        dJ313L*PDstandardNth3gyz;
      
      JacPDstandardNth13gzz = J13L*(J11L*PDstandardNth11gzz + 
        J21L*PDstandardNth12gzz + J31L*PDstandardNth13gzz) + 
        J11L*(J23L*PDstandardNth12gzz + J33L*PDstandardNth13gzz) + 
        dJ113L*PDstandardNth1gzz + J23L*(J21L*PDstandardNth22gzz + 
        J31L*PDstandardNth23gzz) + dJ213L*PDstandardNth2gzz + 
        J33L*(J21L*PDstandardNth23gzz + J31L*PDstandardNth33gzz) + 
        dJ313L*PDstandardNth3gzz;
      
      JacPDstandardNth21gxx = J12L*(J11L*PDstandardNth11gxx + 
        J21L*PDstandardNth12gxx + J31L*PDstandardNth13gxx) + 
        J11L*(J22L*PDstandardNth12gxx + J32L*PDstandardNth13gxx) + 
        dJ112L*PDstandardNth1gxx + J22L*(J21L*PDstandardNth22gxx + 
        J31L*PDstandardNth23gxx) + dJ212L*PDstandardNth2gxx + 
        J32L*(J21L*PDstandardNth23gxx + J31L*PDstandardNth33gxx) + 
        dJ312L*PDstandardNth3gxx;
      
      JacPDstandardNth21gxy = J12L*(J11L*PDstandardNth11gxy + 
        J21L*PDstandardNth12gxy + J31L*PDstandardNth13gxy) + 
        J11L*(J22L*PDstandardNth12gxy + J32L*PDstandardNth13gxy) + 
        dJ112L*PDstandardNth1gxy + J22L*(J21L*PDstandardNth22gxy + 
        J31L*PDstandardNth23gxy) + dJ212L*PDstandardNth2gxy + 
        J32L*(J21L*PDstandardNth23gxy + J31L*PDstandardNth33gxy) + 
        dJ312L*PDstandardNth3gxy;
      
      JacPDstandardNth21gxz = J12L*(J11L*PDstandardNth11gxz + 
        J21L*PDstandardNth12gxz + J31L*PDstandardNth13gxz) + 
        J11L*(J22L*PDstandardNth12gxz + J32L*PDstandardNth13gxz) + 
        dJ112L*PDstandardNth1gxz + J22L*(J21L*PDstandardNth22gxz + 
        J31L*PDstandardNth23gxz) + dJ212L*PDstandardNth2gxz + 
        J32L*(J21L*PDstandardNth23gxz + J31L*PDstandardNth33gxz) + 
        dJ312L*PDstandardNth3gxz;
      
      JacPDstandardNth21gyy = J12L*(J11L*PDstandardNth11gyy + 
        J21L*PDstandardNth12gyy + J31L*PDstandardNth13gyy) + 
        J11L*(J22L*PDstandardNth12gyy + J32L*PDstandardNth13gyy) + 
        dJ112L*PDstandardNth1gyy + J22L*(J21L*PDstandardNth22gyy + 
        J31L*PDstandardNth23gyy) + dJ212L*PDstandardNth2gyy + 
        J32L*(J21L*PDstandardNth23gyy + J31L*PDstandardNth33gyy) + 
        dJ312L*PDstandardNth3gyy;
      
      JacPDstandardNth21gyz = J12L*(J11L*PDstandardNth11gyz + 
        J21L*PDstandardNth12gyz + J31L*PDstandardNth13gyz) + 
        J11L*(J22L*PDstandardNth12gyz + J32L*PDstandardNth13gyz) + 
        dJ112L*PDstandardNth1gyz + J22L*(J21L*PDstandardNth22gyz + 
        J31L*PDstandardNth23gyz) + dJ212L*PDstandardNth2gyz + 
        J32L*(J21L*PDstandardNth23gyz + J31L*PDstandardNth33gyz) + 
        dJ312L*PDstandardNth3gyz;
      
      JacPDstandardNth21gzz = J12L*(J11L*PDstandardNth11gzz + 
        J21L*PDstandardNth12gzz + J31L*PDstandardNth13gzz) + 
        J11L*(J22L*PDstandardNth12gzz + J32L*PDstandardNth13gzz) + 
        dJ112L*PDstandardNth1gzz + J22L*(J21L*PDstandardNth22gzz + 
        J31L*PDstandardNth23gzz) + dJ212L*PDstandardNth2gzz + 
        J32L*(J21L*PDstandardNth23gzz + J31L*PDstandardNth33gzz) + 
        dJ312L*PDstandardNth3gzz;
      
      JacPDstandardNth23gxx = J13L*(J12L*PDstandardNth11gxx + 
        J22L*PDstandardNth12gxx + J32L*PDstandardNth13gxx) + 
        J12L*(J23L*PDstandardNth12gxx + J33L*PDstandardNth13gxx) + 
        dJ123L*PDstandardNth1gxx + J23L*(J22L*PDstandardNth22gxx + 
        J32L*PDstandardNth23gxx) + dJ223L*PDstandardNth2gxx + 
        J33L*(J22L*PDstandardNth23gxx + J32L*PDstandardNth33gxx) + 
        dJ323L*PDstandardNth3gxx;
      
      JacPDstandardNth23gxy = J13L*(J12L*PDstandardNth11gxy + 
        J22L*PDstandardNth12gxy + J32L*PDstandardNth13gxy) + 
        J12L*(J23L*PDstandardNth12gxy + J33L*PDstandardNth13gxy) + 
        dJ123L*PDstandardNth1gxy + J23L*(J22L*PDstandardNth22gxy + 
        J32L*PDstandardNth23gxy) + dJ223L*PDstandardNth2gxy + 
        J33L*(J22L*PDstandardNth23gxy + J32L*PDstandardNth33gxy) + 
        dJ323L*PDstandardNth3gxy;
      
      JacPDstandardNth23gxz = J13L*(J12L*PDstandardNth11gxz + 
        J22L*PDstandardNth12gxz + J32L*PDstandardNth13gxz) + 
        J12L*(J23L*PDstandardNth12gxz + J33L*PDstandardNth13gxz) + 
        dJ123L*PDstandardNth1gxz + J23L*(J22L*PDstandardNth22gxz + 
        J32L*PDstandardNth23gxz) + dJ223L*PDstandardNth2gxz + 
        J33L*(J22L*PDstandardNth23gxz + J32L*PDstandardNth33gxz) + 
        dJ323L*PDstandardNth3gxz;
      
      JacPDstandardNth23gyy = J13L*(J12L*PDstandardNth11gyy + 
        J22L*PDstandardNth12gyy + J32L*PDstandardNth13gyy) + 
        J12L*(J23L*PDstandardNth12gyy + J33L*PDstandardNth13gyy) + 
        dJ123L*PDstandardNth1gyy + J23L*(J22L*PDstandardNth22gyy + 
        J32L*PDstandardNth23gyy) + dJ223L*PDstandardNth2gyy + 
        J33L*(J22L*PDstandardNth23gyy + J32L*PDstandardNth33gyy) + 
        dJ323L*PDstandardNth3gyy;
      
      JacPDstandardNth23gyz = J13L*(J12L*PDstandardNth11gyz + 
        J22L*PDstandardNth12gyz + J32L*PDstandardNth13gyz) + 
        J12L*(J23L*PDstandardNth12gyz + J33L*PDstandardNth13gyz) + 
        dJ123L*PDstandardNth1gyz + J23L*(J22L*PDstandardNth22gyz + 
        J32L*PDstandardNth23gyz) + dJ223L*PDstandardNth2gyz + 
        J33L*(J22L*PDstandardNth23gyz + J32L*PDstandardNth33gyz) + 
        dJ323L*PDstandardNth3gyz;
      
      JacPDstandardNth23gzz = J13L*(J12L*PDstandardNth11gzz + 
        J22L*PDstandardNth12gzz + J32L*PDstandardNth13gzz) + 
        J12L*(J23L*PDstandardNth12gzz + J33L*PDstandardNth13gzz) + 
        dJ123L*PDstandardNth1gzz + J23L*(J22L*PDstandardNth22gzz + 
        J32L*PDstandardNth23gzz) + dJ223L*PDstandardNth2gzz + 
        J33L*(J22L*PDstandardNth23gzz + J32L*PDstandardNth33gzz) + 
        dJ323L*PDstandardNth3gzz;
      
      JacPDstandardNth31gxx = J13L*(J11L*PDstandardNth11gxx + 
        J21L*PDstandardNth12gxx + J31L*PDstandardNth13gxx) + 
        J11L*(J23L*PDstandardNth12gxx + J33L*PDstandardNth13gxx) + 
        dJ113L*PDstandardNth1gxx + J23L*(J21L*PDstandardNth22gxx + 
        J31L*PDstandardNth23gxx) + dJ213L*PDstandardNth2gxx + 
        J33L*(J21L*PDstandardNth23gxx + J31L*PDstandardNth33gxx) + 
        dJ313L*PDstandardNth3gxx;
      
      JacPDstandardNth31gxy = J13L*(J11L*PDstandardNth11gxy + 
        J21L*PDstandardNth12gxy + J31L*PDstandardNth13gxy) + 
        J11L*(J23L*PDstandardNth12gxy + J33L*PDstandardNth13gxy) + 
        dJ113L*PDstandardNth1gxy + J23L*(J21L*PDstandardNth22gxy + 
        J31L*PDstandardNth23gxy) + dJ213L*PDstandardNth2gxy + 
        J33L*(J21L*PDstandardNth23gxy + J31L*PDstandardNth33gxy) + 
        dJ313L*PDstandardNth3gxy;
      
      JacPDstandardNth31gxz = J13L*(J11L*PDstandardNth11gxz + 
        J21L*PDstandardNth12gxz + J31L*PDstandardNth13gxz) + 
        J11L*(J23L*PDstandardNth12gxz + J33L*PDstandardNth13gxz) + 
        dJ113L*PDstandardNth1gxz + J23L*(J21L*PDstandardNth22gxz + 
        J31L*PDstandardNth23gxz) + dJ213L*PDstandardNth2gxz + 
        J33L*(J21L*PDstandardNth23gxz + J31L*PDstandardNth33gxz) + 
        dJ313L*PDstandardNth3gxz;
      
      JacPDstandardNth31gyy = J13L*(J11L*PDstandardNth11gyy + 
        J21L*PDstandardNth12gyy + J31L*PDstandardNth13gyy) + 
        J11L*(J23L*PDstandardNth12gyy + J33L*PDstandardNth13gyy) + 
        dJ113L*PDstandardNth1gyy + J23L*(J21L*PDstandardNth22gyy + 
        J31L*PDstandardNth23gyy) + dJ213L*PDstandardNth2gyy + 
        J33L*(J21L*PDstandardNth23gyy + J31L*PDstandardNth33gyy) + 
        dJ313L*PDstandardNth3gyy;
      
      JacPDstandardNth31gyz = J13L*(J11L*PDstandardNth11gyz + 
        J21L*PDstandardNth12gyz + J31L*PDstandardNth13gyz) + 
        J11L*(J23L*PDstandardNth12gyz + J33L*PDstandardNth13gyz) + 
        dJ113L*PDstandardNth1gyz + J23L*(J21L*PDstandardNth22gyz + 
        J31L*PDstandardNth23gyz) + dJ213L*PDstandardNth2gyz + 
        J33L*(J21L*PDstandardNth23gyz + J31L*PDstandardNth33gyz) + 
        dJ313L*PDstandardNth3gyz;
      
      JacPDstandardNth31gzz = J13L*(J11L*PDstandardNth11gzz + 
        J21L*PDstandardNth12gzz + J31L*PDstandardNth13gzz) + 
        J11L*(J23L*PDstandardNth12gzz + J33L*PDstandardNth13gzz) + 
        dJ113L*PDstandardNth1gzz + J23L*(J21L*PDstandardNth22gzz + 
        J31L*PDstandardNth23gzz) + dJ213L*PDstandardNth2gzz + 
        J33L*(J21L*PDstandardNth23gzz + J31L*PDstandardNth33gzz) + 
        dJ313L*PDstandardNth3gzz;
      
      JacPDstandardNth32gxx = J13L*(J12L*PDstandardNth11gxx + 
        J22L*PDstandardNth12gxx + J32L*PDstandardNth13gxx) + 
        J12L*(J23L*PDstandardNth12gxx + J33L*PDstandardNth13gxx) + 
        dJ123L*PDstandardNth1gxx + J23L*(J22L*PDstandardNth22gxx + 
        J32L*PDstandardNth23gxx) + dJ223L*PDstandardNth2gxx + 
        J33L*(J22L*PDstandardNth23gxx + J32L*PDstandardNth33gxx) + 
        dJ323L*PDstandardNth3gxx;
      
      JacPDstandardNth32gxy = J13L*(J12L*PDstandardNth11gxy + 
        J22L*PDstandardNth12gxy + J32L*PDstandardNth13gxy) + 
        J12L*(J23L*PDstandardNth12gxy + J33L*PDstandardNth13gxy) + 
        dJ123L*PDstandardNth1gxy + J23L*(J22L*PDstandardNth22gxy + 
        J32L*PDstandardNth23gxy) + dJ223L*PDstandardNth2gxy + 
        J33L*(J22L*PDstandardNth23gxy + J32L*PDstandardNth33gxy) + 
        dJ323L*PDstandardNth3gxy;
      
      JacPDstandardNth32gxz = J13L*(J12L*PDstandardNth11gxz + 
        J22L*PDstandardNth12gxz + J32L*PDstandardNth13gxz) + 
        J12L*(J23L*PDstandardNth12gxz + J33L*PDstandardNth13gxz) + 
        dJ123L*PDstandardNth1gxz + J23L*(J22L*PDstandardNth22gxz + 
        J32L*PDstandardNth23gxz) + dJ223L*PDstandardNth2gxz + 
        J33L*(J22L*PDstandardNth23gxz + J32L*PDstandardNth33gxz) + 
        dJ323L*PDstandardNth3gxz;
      
      JacPDstandardNth32gyy = J13L*(J12L*PDstandardNth11gyy + 
        J22L*PDstandardNth12gyy + J32L*PDstandardNth13gyy) + 
        J12L*(J23L*PDstandardNth12gyy + J33L*PDstandardNth13gyy) + 
        dJ123L*PDstandardNth1gyy + J23L*(J22L*PDstandardNth22gyy + 
        J32L*PDstandardNth23gyy) + dJ223L*PDstandardNth2gyy + 
        J33L*(J22L*PDstandardNth23gyy + J32L*PDstandardNth33gyy) + 
        dJ323L*PDstandardNth3gyy;
      
      JacPDstandardNth32gyz = J13L*(J12L*PDstandardNth11gyz + 
        J22L*PDstandardNth12gyz + J32L*PDstandardNth13gyz) + 
        J12L*(J23L*PDstandardNth12gyz + J33L*PDstandardNth13gyz) + 
        dJ123L*PDstandardNth1gyz + J23L*(J22L*PDstandardNth22gyz + 
        J32L*PDstandardNth23gyz) + dJ223L*PDstandardNth2gyz + 
        J33L*(J22L*PDstandardNth23gyz + J32L*PDstandardNth33gyz) + 
        dJ323L*PDstandardNth3gyz;
      
      JacPDstandardNth32gzz = J13L*(J12L*PDstandardNth11gzz + 
        J22L*PDstandardNth12gzz + J32L*PDstandardNth13gzz) + 
        J12L*(J23L*PDstandardNth12gzz + J33L*PDstandardNth13gzz) + 
        dJ123L*PDstandardNth1gzz + J23L*(J22L*PDstandardNth22gzz + 
        J32L*PDstandardNth23gzz) + dJ223L*PDstandardNth2gzz + 
        J33L*(J22L*PDstandardNth23gzz + J32L*PDstandardNth33gzz) + 
        dJ323L*PDstandardNth3gzz;
    }
    else
    {
      JacPDstandardNth1gxx = PDstandardNth1gxx;
      
      JacPDstandardNth1gxy = PDstandardNth1gxy;
      
      JacPDstandardNth1gxz = PDstandardNth1gxz;
      
      JacPDstandardNth1gyy = PDstandardNth1gyy;
      
      JacPDstandardNth1gyz = PDstandardNth1gyz;
      
      JacPDstandardNth1gzz = PDstandardNth1gzz;
      
      JacPDstandardNth1kxy = PDstandardNth1kxy;
      
      JacPDstandardNth1kxz = PDstandardNth1kxz;
      
      JacPDstandardNth1kyy = PDstandardNth1kyy;
      
      JacPDstandardNth1kyz = PDstandardNth1kyz;
      
      JacPDstandardNth1kzz = PDstandardNth1kzz;
      
      JacPDstandardNth2gxx = PDstandardNth2gxx;
      
      JacPDstandardNth2gxy = PDstandardNth2gxy;
      
      JacPDstandardNth2gxz = PDstandardNth2gxz;
      
      JacPDstandardNth2gyy = PDstandardNth2gyy;
      
      JacPDstandardNth2gyz = PDstandardNth2gyz;
      
      JacPDstandardNth2gzz = PDstandardNth2gzz;
      
      JacPDstandardNth2kxx = PDstandardNth2kxx;
      
      JacPDstandardNth2kxy = PDstandardNth2kxy;
      
      JacPDstandardNth2kxz = PDstandardNth2kxz;
      
      JacPDstandardNth2kyz = PDstandardNth2kyz;
      
      JacPDstandardNth2kzz = PDstandardNth2kzz;
      
      JacPDstandardNth3gxx = PDstandardNth3gxx;
      
      JacPDstandardNth3gxy = PDstandardNth3gxy;
      
      JacPDstandardNth3gxz = PDstandardNth3gxz;
      
      JacPDstandardNth3gyy = PDstandardNth3gyy;
      
      JacPDstandardNth3gyz = PDstandardNth3gyz;
      
      JacPDstandardNth3gzz = PDstandardNth3gzz;
      
      JacPDstandardNth3kxx = PDstandardNth3kxx;
      
      JacPDstandardNth3kxy = PDstandardNth3kxy;
      
      JacPDstandardNth3kxz = PDstandardNth3kxz;
      
      JacPDstandardNth3kyy = PDstandardNth3kyy;
      
      JacPDstandardNth3kyz = PDstandardNth3kyz;
      
      JacPDstandardNth11gyy = PDstandardNth11gyy;
      
      JacPDstandardNth11gyz = PDstandardNth11gyz;
      
      JacPDstandardNth11gzz = PDstandardNth11gzz;
      
      JacPDstandardNth22gxx = PDstandardNth22gxx;
      
      JacPDstandardNth22gxz = PDstandardNth22gxz;
      
      JacPDstandardNth22gzz = PDstandardNth22gzz;
      
      JacPDstandardNth33gxx = PDstandardNth33gxx;
      
      JacPDstandardNth33gxy = PDstandardNth33gxy;
      
      JacPDstandardNth33gyy = PDstandardNth33gyy;
      
      JacPDstandardNth12gxx = PDstandardNth12gxx;
      
      JacPDstandardNth12gxy = PDstandardNth12gxy;
      
      JacPDstandardNth12gxz = PDstandardNth12gxz;
      
      JacPDstandardNth12gyy = PDstandardNth12gyy;
      
      JacPDstandardNth12gyz = PDstandardNth12gyz;
      
      JacPDstandardNth12gzz = PDstandardNth12gzz;
      
      JacPDstandardNth13gxx = PDstandardNth13gxx;
      
      JacPDstandardNth13gxy = PDstandardNth13gxy;
      
      JacPDstandardNth13gxz = PDstandardNth13gxz;
      
      JacPDstandardNth13gyy = PDstandardNth13gyy;
      
      JacPDstandardNth13gyz = PDstandardNth13gyz;
      
      JacPDstandardNth13gzz = PDstandardNth13gzz;
      
      JacPDstandardNth21gxx = PDstandardNth12gxx;
      
      JacPDstandardNth21gxy = PDstandardNth12gxy;
      
      JacPDstandardNth21gxz = PDstandardNth12gxz;
      
      JacPDstandardNth21gyy = PDstandardNth12gyy;
      
      JacPDstandardNth21gyz = PDstandardNth12gyz;
      
      JacPDstandardNth21gzz = PDstandardNth12gzz;
      
      JacPDstandardNth23gxx = PDstandardNth23gxx;
      
      JacPDstandardNth23gxy = PDstandardNth23gxy;
      
      JacPDstandardNth23gxz = PDstandardNth23gxz;
      
      JacPDstandardNth23gyy = PDstandardNth23gyy;
      
      JacPDstandardNth23gyz = PDstandardNth23gyz;
      
      JacPDstandardNth23gzz = PDstandardNth23gzz;
      
      JacPDstandardNth31gxx = PDstandardNth13gxx;
      
      JacPDstandardNth31gxy = PDstandardNth13gxy;
      
      JacPDstandardNth31gxz = PDstandardNth13gxz;
      
      JacPDstandardNth31gyy = PDstandardNth13gyy;
      
      JacPDstandardNth31gyz = PDstandardNth13gyz;
      
      JacPDstandardNth31gzz = PDstandardNth13gzz;
      
      JacPDstandardNth32gxx = PDstandardNth23gxx;
      
      JacPDstandardNth32gxy = PDstandardNth23gxy;
      
      JacPDstandardNth32gxz = PDstandardNth23gxz;
      
      JacPDstandardNth32gyy = PDstandardNth23gyy;
      
      JacPDstandardNth32gyz = PDstandardNth23gyz;
      
      JacPDstandardNth32gzz = PDstandardNth23gzz;
    }
    
    CCTK_REAL detg = 2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - 
      SQR(gxyL)) - gyyL*SQR(gxzL) - gxxL*SQR(gyzL);
    
    CCTK_REAL gu11 = INV(detg)*(gyyL*gzzL - SQR(gyzL));
    
    CCTK_REAL gu21 = (gxzL*gyzL - gxyL*gzzL)*INV(detg);
    
    CCTK_REAL gu31 = (-(gxzL*gyyL) + gxyL*gyzL)*INV(detg);
    
    CCTK_REAL gu22 = INV(detg)*(gxxL*gzzL - SQR(gxzL));
    
    CCTK_REAL gu32 = (gxyL*gxzL - gxxL*gyzL)*INV(detg);
    
    CCTK_REAL gu33 = INV(detg)*(gxxL*gyyL - SQR(gxyL));
    
    CCTK_REAL G111 = 0.5*(gu11*JacPDstandardNth1gxx + 
      2*(gu21*JacPDstandardNth1gxy + gu31*JacPDstandardNth1gxz) - 
      gu21*JacPDstandardNth2gxx - gu31*JacPDstandardNth3gxx);
    
    CCTK_REAL G211 = 0.5*(gu21*JacPDstandardNth1gxx + 
      2*(gu22*JacPDstandardNth1gxy + gu32*JacPDstandardNth1gxz) - 
      gu22*JacPDstandardNth2gxx - gu32*JacPDstandardNth3gxx);
    
    CCTK_REAL G311 = 0.5*(gu31*JacPDstandardNth1gxx + 
      2*(gu32*JacPDstandardNth1gxy + gu33*JacPDstandardNth1gxz) - 
      gu32*JacPDstandardNth2gxx - gu33*JacPDstandardNth3gxx);
    
    CCTK_REAL G112 = 0.5*(gu21*JacPDstandardNth1gyy + 
      gu11*JacPDstandardNth2gxx + gu31*(JacPDstandardNth1gyz + 
      JacPDstandardNth2gxz - JacPDstandardNth3gxy));
    
    CCTK_REAL G212 = 0.5*(gu22*JacPDstandardNth1gyy + 
      gu21*JacPDstandardNth2gxx + gu32*(JacPDstandardNth1gyz + 
      JacPDstandardNth2gxz - JacPDstandardNth3gxy));
    
    CCTK_REAL G312 = 0.5*(gu32*JacPDstandardNth1gyy + 
      gu31*JacPDstandardNth2gxx + gu33*(JacPDstandardNth1gyz + 
      JacPDstandardNth2gxz - JacPDstandardNth3gxy));
    
    CCTK_REAL G113 = 0.5*(gu31*JacPDstandardNth1gzz + 
      gu11*JacPDstandardNth3gxx + gu21*(JacPDstandardNth1gyz - 
      JacPDstandardNth2gxz + JacPDstandardNth3gxy));
    
    CCTK_REAL G213 = 0.5*(gu32*JacPDstandardNth1gzz + 
      gu21*JacPDstandardNth3gxx + gu22*(JacPDstandardNth1gyz - 
      JacPDstandardNth2gxz + JacPDstandardNth3gxy));
    
    CCTK_REAL G313 = 0.5*(gu33*JacPDstandardNth1gzz + 
      gu31*JacPDstandardNth3gxx + gu32*(JacPDstandardNth1gyz - 
      JacPDstandardNth2gxz + JacPDstandardNth3gxy));
    
    CCTK_REAL G122 = 0.5*(gu11*(-JacPDstandardNth1gyy + 
      2*JacPDstandardNth2gxy) + gu21*JacPDstandardNth2gyy + 
      gu31*(2*JacPDstandardNth2gyz - JacPDstandardNth3gyy));
    
    CCTK_REAL G222 = 0.5*(gu21*(-JacPDstandardNth1gyy + 
      2*JacPDstandardNth2gxy) + gu22*JacPDstandardNth2gyy + 
      gu32*(2*JacPDstandardNth2gyz - JacPDstandardNth3gyy));
    
    CCTK_REAL G322 = 0.5*(gu31*(-JacPDstandardNth1gyy + 
      2*JacPDstandardNth2gxy) + gu32*JacPDstandardNth2gyy + 
      gu33*(2*JacPDstandardNth2gyz - JacPDstandardNth3gyy));
    
    CCTK_REAL G123 = 0.5*(gu31*JacPDstandardNth2gzz + 
      gu11*(-JacPDstandardNth1gyz + JacPDstandardNth2gxz + 
      JacPDstandardNth3gxy) + gu21*JacPDstandardNth3gyy);
    
    CCTK_REAL G223 = 0.5*(gu32*JacPDstandardNth2gzz + 
      gu21*(-JacPDstandardNth1gyz + JacPDstandardNth2gxz + 
      JacPDstandardNth3gxy) + gu22*JacPDstandardNth3gyy);
    
    CCTK_REAL G323 = 0.5*(gu33*JacPDstandardNth2gzz + 
      gu31*(-JacPDstandardNth1gyz + JacPDstandardNth2gxz + 
      JacPDstandardNth3gxy) + gu32*JacPDstandardNth3gyy);
    
    CCTK_REAL G133 = 0.5*(gu11*(-JacPDstandardNth1gzz + 
      2*JacPDstandardNth3gxz) + gu21*(-JacPDstandardNth2gzz + 
      2*JacPDstandardNth3gyz) + gu31*JacPDstandardNth3gzz);
    
    CCTK_REAL G233 = 0.5*(gu21*(-JacPDstandardNth1gzz + 
      2*JacPDstandardNth3gxz) + gu22*(-JacPDstandardNth2gzz + 
      2*JacPDstandardNth3gyz) + gu32*JacPDstandardNth3gzz);
    
    CCTK_REAL G333 = 0.5*(gu31*(-JacPDstandardNth1gzz + 
      2*JacPDstandardNth3gxz) + gu32*(-JacPDstandardNth2gzz + 
      2*JacPDstandardNth3gyz) + gu33*JacPDstandardNth3gzz);
    
    CCTK_REAL R11 = 0.5*(gu21*(-JacPDstandardNth12gxx + 
      JacPDstandardNth21gxx) + gu31*(-JacPDstandardNth13gxx + 
      JacPDstandardNth31gxx) + gu32*(4*(gxyL*(-(G123*G211) + G113*G212) + 
      gyyL*(G212*G213 - G211*G223) + gxzL*G113*G312 + gyzL*G212*G313 + 
      gzzL*G312*G313 + G112*(gxxL*G113 + gxyL*G213 + gxzL*G313) - 
      G111*(gxxL*G123 + gxyL*G223 + gxzL*G323) + G311*(-(gxzL*G123) - 
      gyzL*G223 - gzzL*G323) + gyzL*(G213*G312 - G211*G323)) - 
      JacPDstandardNth11gyz + JacPDstandardNth21gxz - JacPDstandardNth23gxx + 
      JacPDstandardNth31gxy) + gu32*(-JacPDstandardNth11gyz + 
      JacPDstandardNth21gxz + JacPDstandardNth31gxy - JacPDstandardNth32gxx) 
      + gu22*(-JacPDstandardNth11gyy + 2*JacPDstandardNth21gxy - 
      JacPDstandardNth22gxx + 2*(2*gyzL*G212*G312 + 2*G112*(gxyL*G212 + 
      gxzL*G312) - G111*(gxxL*G122 + gxyL*G222 + gxzL*G322) + 
      G211*(-(gxyL*G122) - gyyL*G222 - gyzL*G322) + 
      G311*(-(gxzL*G122) - gyzL*G222 - gzzL*G322) + gxxL*SQR(G112) + 
      gyyL*SQR(G212) + gzzL*SQR(G312))) + gu33*(-JacPDstandardNth11gzz + 
      2*JacPDstandardNth31gxz - JacPDstandardNth33gxx + 2*(2*gyzL*G213*G313 
      + 2*G113*(gxyL*G213 + gxzL*G313) - G111*(gxxL*G133 + gxyL*G233 
      + gxzL*G333) + G211*(-(gxyL*G133) - gyyL*G233 - gyzL*G333) + 
      G311*(-(gxzL*G133) - gyzL*G233 - gzzL*G333) + gxxL*SQR(G113) + 
      gyyL*SQR(G213) + gzzL*SQR(G313))));
    
    CCTK_REAL R12 = 0.5*(gu22*(-JacPDstandardNth12gyy + 
      JacPDstandardNth21gyy) + gu21*(JacPDstandardNth11gyy - 
      2*JacPDstandardNth12gxy + JacPDstandardNth22gxx) + 
      gu32*(-2*JacPDstandardNth12gyz + JacPDstandardNth21gyz + 
      JacPDstandardNth22gxz - JacPDstandardNth23gxy + JacPDstandardNth31gyy) 
      + gu31*(JacPDstandardNth11gyz - JacPDstandardNth12gxz - 
      JacPDstandardNth13gxy + JacPDstandardNth32gxx) + 
      gu33*(-JacPDstandardNth12gzz + JacPDstandardNth31gyz + 
      JacPDstandardNth32gxz - JacPDstandardNth33gxy) + 2*((gxyL*(G123*G211 
      - G113*G212) + gyyL*(-(G212*G213) + G211*G223) - G112*(gxxL*G113 + 
      gxyL*G213 + gxzL*G313) + G312*(-(gxzL*G113) - gyzL*G213 - 
      gzzL*G313) + G111*(gxxL*G123 + gxyL*G223 + gxzL*G323) + 
      G311*(gxzL*G123 + gyzL*G223 + gzzL*G323) + gyzL*(-(G212*G313) + 
      G211*G323))*gu31 + (gxyL*(-(G123*G212) + G122*G213) + 
      gyyL*(G213*G222 - G212*G223) + gxzL*G122*G313 + gyzL*G213*G322 + 
      gzzL*G313*G322 + G113*(gxxL*G122 + gxyL*G222 + gxzL*G322) - 
      G112*(gxxL*G123 + gxyL*G223 + gxzL*G323) + G312*(-(gxzL*G123) - 
      gyzL*G223 - gzzL*G323) + gyzL*(G222*G313 - G212*G323))*gu32 + 
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
    
    CCTK_REAL R13 = 0.5*(gu21*(JacPDstandardNth11gyz - 
      JacPDstandardNth12gxz - JacPDstandardNth13gxy + JacPDstandardNth23gxx) 
      + gu22*(-JacPDstandardNth13gyy + JacPDstandardNth21gyz - 
      JacPDstandardNth22gxz + JacPDstandardNth23gxy) + 
      gu33*(-JacPDstandardNth13gzz + JacPDstandardNth31gzz) + 
      gu31*(JacPDstandardNth11gzz - 2*JacPDstandardNth13gxz + 
      JacPDstandardNth33gxx) + gu32*(-2*JacPDstandardNth13gyz + 
      JacPDstandardNth21gzz + JacPDstandardNth31gyz - JacPDstandardNth32gxz + 
      JacPDstandardNth33gxy) + 2*((gxyL*(G123*G211 - G113*G212) + 
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
    
    CCTK_REAL R22 = 0.5*(gu21*(JacPDstandardNth12gyy - 
      JacPDstandardNth21gyy) + gu31*(4*(gxyL*(G123*G212 - G122*G213) + 
      gyyL*(-(G213*G222) + G212*G223) + gxzL*(G123*G312 - G122*G313) - 
      G113*(gxxL*G122 + gxyL*G222 + gxzL*G322) + gyzL*(G223*G312 - 
      G222*G313 - G213*G322) + gyzL*G212*G323 + G112*(gxxL*G123 + 
      gxyL*G223 + gxzL*G323) + gzzL*(-(G313*G322) + G312*G323)) + 
      JacPDstandardNth12gyz - JacPDstandardNth13gyy - JacPDstandardNth22gxz + 
      JacPDstandardNth32gxy) + gu31*(JacPDstandardNth12gyz - 
      JacPDstandardNth22gxz - JacPDstandardNth31gyy + JacPDstandardNth32gxy) 
      + gu32*(-JacPDstandardNth23gyy + JacPDstandardNth32gyy) + 
      gu11*(-JacPDstandardNth11gyy + 2*JacPDstandardNth12gxy - 
      JacPDstandardNth22gxx + 2*(2*gyzL*G212*G312 + 2*G112*(gxyL*G212 + 
      gxzL*G312) - G111*(gxxL*G122 + gxyL*G222 + gxzL*G322) + 
      G211*(-(gxyL*G122) - gyyL*G222 - gyzL*G322) + 
      G311*(-(gxzL*G122) - gyzL*G222 - gzzL*G322) + gxxL*SQR(G112) + 
      gyyL*SQR(G212) + gzzL*SQR(G312))) + gu33*(-JacPDstandardNth22gzz + 
      2*JacPDstandardNth32gyz - JacPDstandardNth33gyy + 2*(2*gyzL*G223*G323 
      + 2*G123*(gxyL*G223 + gxzL*G323) - G122*(gxxL*G133 + gxyL*G233 
      + gxzL*G333) + G222*(-(gxyL*G133) - gyyL*G233 - gyzL*G333) + 
      G322*(-(gxzL*G133) - gyzL*G233 - gzzL*G333) + gxxL*SQR(G123) + 
      gyyL*SQR(G223) + gzzL*SQR(G323))));
    
    CCTK_REAL R23 = 0.5*(gu11*(-JacPDstandardNth11gyz + 
      JacPDstandardNth12gxz + JacPDstandardNth13gxy - JacPDstandardNth23gxx) 
      + gu21*(JacPDstandardNth13gyy - JacPDstandardNth21gyz + 
      JacPDstandardNth22gxz - JacPDstandardNth23gxy) + 
      gu33*(-JacPDstandardNth23gzz + JacPDstandardNth32gzz) + 
      gu31*(JacPDstandardNth12gzz - 2*JacPDstandardNth23gxz - 
      JacPDstandardNth31gyz + JacPDstandardNth32gxz + JacPDstandardNth33gxy) 
      + gu32*(JacPDstandardNth22gzz - 2*JacPDstandardNth23gyz + 
      JacPDstandardNth33gyy) + 2*((gxyL*(-(G123*G211) + G113*G212) + 
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
    
    CCTK_REAL R33 = 0.5*(gu31*(JacPDstandardNth13gzz - 
      JacPDstandardNth31gzz) + gu32*(JacPDstandardNth23gzz - 
      JacPDstandardNth32gzz) + gu21*(4*(gxyL*(-(G133*G212) + G123*G213) + 
      gyyL*(G213*G223 - G212*G233) + gxzL*G123*G313 + gyzL*G213*G323 + 
      gzzL*G313*G323 + G113*(gxxL*G123 + gxyL*G223 + gxzL*G323) - 
      G112*(gxxL*G133 + gxyL*G233 + gxzL*G333) + G312*(-(gxzL*G133) - 
      gyzL*G233 - gzzL*G333) + gyzL*(G223*G313 - G212*G333)) - 
      JacPDstandardNth12gzz + JacPDstandardNth13gyz + JacPDstandardNth23gxz - 
      JacPDstandardNth33gxy) + gu21*(JacPDstandardNth13gyz - 
      JacPDstandardNth21gzz + JacPDstandardNth23gxz - JacPDstandardNth33gxy) 
      + gu11*(-JacPDstandardNth11gzz + 2*JacPDstandardNth13gxz - 
      JacPDstandardNth33gxx + 2*(2*gyzL*G213*G313 + 2*G113*(gxyL*G213 + 
      gxzL*G313) - G111*(gxxL*G133 + gxyL*G233 + gxzL*G333) + 
      G211*(-(gxyL*G133) - gyyL*G233 - gyzL*G333) + 
      G311*(-(gxzL*G133) - gyzL*G233 - gzzL*G333) + gxxL*SQR(G113) + 
      gyyL*SQR(G213) + gzzL*SQR(G313))) + gu22*(-JacPDstandardNth22gzz + 
      2*JacPDstandardNth23gyz - JacPDstandardNth33gyy + 2*(2*gyzL*G223*G323 
      + 2*G123*(gxyL*G223 + gxzL*G323) - G122*(gxxL*G133 + gxyL*G233 
      + gxzL*G333) + G222*(-(gxyL*G133) - gyyL*G233 - gyzL*G333) + 
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
      G212) + kyzL*G311 - kxzL*G312 - JacPDstandardNth1kxy + 
      JacPDstandardNth2kxx) + gu22*(-(kxxL*G122) + kyyL*G212 + 
      kxyL*(G112 - G222) + kyzL*G312 - kxzL*G322 - JacPDstandardNth1kyy 
      + JacPDstandardNth2kxy) + gu31*(-(kxxL*G113) + kyzL*G211 - 
      kxyL*G213 + kzzL*G311 + kxzL*(G111 - G313) - JacPDstandardNth1kxz 
      + JacPDstandardNth3kxx) + gu32*(kyyL*G213 + kxyL*(G113 - 2*G223) + 
      kzzL*G312 + kyzL*(G212 + G313) + kxzL*(G112 - 2*G323) - 
      2*(kxxL*G123 + JacPDstandardNth1kyz) + JacPDstandardNth2kxz + 
      JacPDstandardNth3kxy) + gu33*(-(kxxL*G133) + kyzL*G213 - 
      kxyL*G233 + kzzL*G313 + kxzL*(G113 - G333) - JacPDstandardNth1kzz 
      + JacPDstandardNth3kxz) - 8*Pi*S1;
    
    CCTK_REAL M2L = gu11*(kxxL*G112 - kyyL*G211 + kxyL*(-G111 + 
      G212) - kyzL*G311 + kxzL*G312 + JacPDstandardNth1kxy - 
      JacPDstandardNth2kxx) + gu21*(kxxL*G122 - kyyL*G212 + kxyL*(-G112 
      + G222) - kyzL*G312 + kxzL*G322 + JacPDstandardNth1kyy - 
      JacPDstandardNth2kxy) + gu31*(kxxL*G123 + kxyL*G223 + kzzL*G312 + 
      kyzL*(G212 - 2*G313) + kxzL*(G112 + G323) + JacPDstandardNth1kyz - 
      2*(kxyL*G113 + kyyL*G213 + JacPDstandardNth2kxz) + 
      JacPDstandardNth3kxy) + gu32*(kxzL*G122 - kxyL*G123 - kyyL*G223 + 
      kzzL*G322 + kyzL*(G222 - G323) - JacPDstandardNth2kyz + 
      JacPDstandardNth3kyy) + gu33*(kxzL*G123 - kxyL*G133 - kyyL*G233 + 
      kzzL*G323 + kyzL*(G223 - G333) - JacPDstandardNth2kzz + 
      JacPDstandardNth3kyz) - 8*Pi*S2;
    
    CCTK_REAL M3L = gu11*(kxxL*G113 - kyzL*G211 + kxyL*G213 - 
      kzzL*G311 + kxzL*(-G111 + G313) + JacPDstandardNth1kxz - 
      JacPDstandardNth3kxx) + gu21*(kxxL*G123 + kyyL*G213 + kxyL*(G113 
      + G223) + kyzL*G313 + kxzL*G323 + JacPDstandardNth1kyz + 
      JacPDstandardNth2kxz - 2*(kxzL*G112 + kyzL*G212 + kzzL*G312 + 
      JacPDstandardNth3kxy)) + gu31*(kxxL*G133 - kyzL*G213 + kxyL*G233 
      - kzzL*G313 + kxzL*(-G113 + G333) + JacPDstandardNth1kzz - 
      JacPDstandardNth3kxz) + gu22*(-(kxzL*G122) + kxyL*G123 + 
      kyyL*G223 - kzzL*G322 + kyzL*(-G222 + G323) + 
      JacPDstandardNth2kyz - JacPDstandardNth3kyy) + gu32*(-(kxzL*G123) + 
      kxyL*G133 + kyyL*G233 - kzzL*G323 + kyzL*(-G223 + G333) + 
      JacPDstandardNth2kzz - JacPDstandardNth3kyz) - 8*Pi*S3;
    
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
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMConstraints_MP_Body");
  }
  
  if (cctk_iteration % ML_ADMConstraints_MP_calc_every != ML_ADMConstraints_MP_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::curv","ADMBase::lapse","ADMBase::metric","ADMBase::shift","ML_ADMConstraints_MP::ML_Ham","ML_ADMConstraints_MP::ML_mom"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_ADMConstraints_MP", 6, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_ADMConstraints_MP", 2, 2, 2);
  
  GenericFD_LoopOverInterior(cctkGH, &ML_ADMConstraints_MP_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_ADMConstraints_MP_Body");
  }
}
