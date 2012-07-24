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
#include "cctk_Loop.h"
#include "loopcontrol.h"
#include "vectors.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define ScalarINV(x) ((CCTK_REAL)1.0 / (x))
#define ScalarSQR(x) ((x) * (x))
#define ScalarCUB(x) ((x) * ScalarSQR(x))
#define ScalarQAD(x) (ScalarSQR(ScalarSQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))
#define QAD(x) (SQR(SQR(x)))

extern "C" void ML_BSSN_Host_RHS2_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_curvrhs.");
  return;
}

static void ML_BSSN_Host_RHS2_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  ptrdiff_t const di = 1;
  ptrdiff_t const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const cdi = sizeof(CCTK_REAL) * di;
  ptrdiff_t const cdj = sizeof(CCTK_REAL) * dj;
  ptrdiff_t const cdk = sizeof(CCTK_REAL) * dk;
  CCTK_REAL_VEC const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL_VEC const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL_VEC const dz = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL_VEC const dt = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL_VEC const t = ToReal(cctk_time);
  CCTK_REAL_VEC const dxi = INV(dx);
  CCTK_REAL_VEC const dyi = INV(dy);
  CCTK_REAL_VEC const dzi = INV(dz);
  CCTK_REAL_VEC const khalf = ToReal(0.5);
  CCTK_REAL_VEC const kthird = ToReal(1.0/3.0);
  CCTK_REAL_VEC const ktwothird = ToReal(2.0/3.0);
  CCTK_REAL_VEC const kfourthird = ToReal(4.0/3.0);
  CCTK_REAL_VEC const keightthird = ToReal(8.0/3.0);
  CCTK_REAL_VEC const hdxi = kmul(ToReal(0.5), dxi);
  CCTK_REAL_VEC const hdyi = kmul(ToReal(0.5), dyi);
  CCTK_REAL_VEC const hdzi = kmul(ToReal(0.5), dzi);
  
  /* Initialize predefined quantities */
  CCTK_REAL_VEC const p1o1024dx = kdiv(ToReal(0.0009765625),dx);
  CCTK_REAL_VEC const p1o1024dy = kdiv(ToReal(0.0009765625),dy);
  CCTK_REAL_VEC const p1o1024dz = kdiv(ToReal(0.0009765625),dz);
  CCTK_REAL_VEC const p1o1680dx = kdiv(ToReal(0.000595238095238095238095238095238),dx);
  CCTK_REAL_VEC const p1o1680dy = kdiv(ToReal(0.000595238095238095238095238095238),dy);
  CCTK_REAL_VEC const p1o1680dz = kdiv(ToReal(0.000595238095238095238095238095238),dz);
  CCTK_REAL_VEC const p1o2dx = kdiv(ToReal(0.5),dx);
  CCTK_REAL_VEC const p1o2dy = kdiv(ToReal(0.5),dy);
  CCTK_REAL_VEC const p1o2dz = kdiv(ToReal(0.5),dz);
  CCTK_REAL_VEC const p1o5040dx2 = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dx,dx));
  CCTK_REAL_VEC const p1o5040dy2 = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dy,dy));
  CCTK_REAL_VEC const p1o5040dz2 = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dz,dz));
  CCTK_REAL_VEC const p1o560dx = kdiv(ToReal(0.00178571428571428571428571428571),dx);
  CCTK_REAL_VEC const p1o560dy = kdiv(ToReal(0.00178571428571428571428571428571),dy);
  CCTK_REAL_VEC const p1o560dz = kdiv(ToReal(0.00178571428571428571428571428571),dz);
  CCTK_REAL_VEC const p1o705600dxdy = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dy,dx));
  CCTK_REAL_VEC const p1o705600dxdz = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dz,dx));
  CCTK_REAL_VEC const p1o705600dydz = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dz,dy));
  CCTK_REAL_VEC const p1o840dx = kdiv(ToReal(0.00119047619047619047619047619048),dx);
  CCTK_REAL_VEC const p1o840dy = kdiv(ToReal(0.00119047619047619047619047619048),dy);
  CCTK_REAL_VEC const p1o840dz = kdiv(ToReal(0.00119047619047619047619047619048),dz);
  CCTK_REAL_VEC const p1odx = kdiv(ToReal(1),dx);
  CCTK_REAL_VEC const p1ody = kdiv(ToReal(1),dy);
  CCTK_REAL_VEC const p1odz = kdiv(ToReal(1),dz);
  CCTK_REAL_VEC const pm1o2dx = kdiv(ToReal(-0.5),dx);
  CCTK_REAL_VEC const pm1o2dy = kdiv(ToReal(-0.5),dy);
  CCTK_REAL_VEC const pm1o2dz = kdiv(ToReal(-0.5),dz);
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC(ML_BSSN_Host_RHS2,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC At11L = vec_load(At11[index]);
    CCTK_REAL_VEC At12L = vec_load(At12[index]);
    CCTK_REAL_VEC At13L = vec_load(At13[index]);
    CCTK_REAL_VEC At22L = vec_load(At22[index]);
    CCTK_REAL_VEC At23L = vec_load(At23[index]);
    CCTK_REAL_VEC At33L = vec_load(At33[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC phiL = vec_load(phi[index]);
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    CCTK_REAL_VEC Xt1L = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L = vec_load(Xt3[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC const PDstandardNth1alpha = PDstandardNth1(&alpha[index]);
    CCTK_REAL_VEC const PDstandardNth2alpha = PDstandardNth2(&alpha[index]);
    CCTK_REAL_VEC const PDstandardNth3alpha = PDstandardNth3(&alpha[index]);
    CCTK_REAL_VEC const PDstandardNth11alpha = PDstandardNth11(&alpha[index]);
    CCTK_REAL_VEC const PDstandardNth22alpha = PDstandardNth22(&alpha[index]);
    CCTK_REAL_VEC const PDstandardNth33alpha = PDstandardNth33(&alpha[index]);
    CCTK_REAL_VEC const PDstandardNth12alpha = PDstandardNth12(&alpha[index]);
    CCTK_REAL_VEC const PDstandardNth13alpha = PDstandardNth13(&alpha[index]);
    CCTK_REAL_VEC const PDstandardNth23alpha = PDstandardNth23(&alpha[index]);
    CCTK_REAL_VEC const PDstandardNth1beta1 = PDstandardNth1(&beta1[index]);
    CCTK_REAL_VEC const PDstandardNth2beta1 = PDstandardNth2(&beta1[index]);
    CCTK_REAL_VEC const PDstandardNth3beta1 = PDstandardNth3(&beta1[index]);
    CCTK_REAL_VEC const PDstandardNth1beta2 = PDstandardNth1(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth2beta2 = PDstandardNth2(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth3beta2 = PDstandardNth3(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth1beta3 = PDstandardNth1(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth2beta3 = PDstandardNth2(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth3beta3 = PDstandardNth3(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth1gt11 = PDstandardNth1(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth2gt11 = PDstandardNth2(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth3gt11 = PDstandardNth3(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth11gt11 = PDstandardNth11(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth22gt11 = PDstandardNth22(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth33gt11 = PDstandardNth33(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth12gt11 = PDstandardNth12(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth13gt11 = PDstandardNth13(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth23gt11 = PDstandardNth23(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth1gt12 = PDstandardNth1(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth2gt12 = PDstandardNth2(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth3gt12 = PDstandardNth3(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth11gt12 = PDstandardNth11(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth22gt12 = PDstandardNth22(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth33gt12 = PDstandardNth33(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth12gt12 = PDstandardNth12(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth13gt12 = PDstandardNth13(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth23gt12 = PDstandardNth23(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth1gt13 = PDstandardNth1(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth2gt13 = PDstandardNth2(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth3gt13 = PDstandardNth3(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth11gt13 = PDstandardNth11(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth22gt13 = PDstandardNth22(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth33gt13 = PDstandardNth33(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth12gt13 = PDstandardNth12(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth13gt13 = PDstandardNth13(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth23gt13 = PDstandardNth23(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth1gt22 = PDstandardNth1(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth2gt22 = PDstandardNth2(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth3gt22 = PDstandardNth3(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth11gt22 = PDstandardNth11(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth22gt22 = PDstandardNth22(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth33gt22 = PDstandardNth33(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth12gt22 = PDstandardNth12(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth13gt22 = PDstandardNth13(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth23gt22 = PDstandardNth23(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth1gt23 = PDstandardNth1(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth2gt23 = PDstandardNth2(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth3gt23 = PDstandardNth3(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth11gt23 = PDstandardNth11(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth22gt23 = PDstandardNth22(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth33gt23 = PDstandardNth33(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth12gt23 = PDstandardNth12(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth13gt23 = PDstandardNth13(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth23gt23 = PDstandardNth23(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth1gt33 = PDstandardNth1(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth2gt33 = PDstandardNth2(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth3gt33 = PDstandardNth3(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth11gt33 = PDstandardNth11(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth22gt33 = PDstandardNth22(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth33gt33 = PDstandardNth33(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth12gt33 = PDstandardNth12(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth13gt33 = PDstandardNth13(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth23gt33 = PDstandardNth23(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth1phi = PDstandardNth1(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth2phi = PDstandardNth2(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth3phi = PDstandardNth3(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth11phi = PDstandardNth11(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth22phi = PDstandardNth22(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth33phi = PDstandardNth33(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth12phi = PDstandardNth12(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth13phi = PDstandardNth13(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth23phi = PDstandardNth23(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth1Xt1 = PDstandardNth1(&Xt1[index]);
    CCTK_REAL_VEC const PDstandardNth2Xt1 = PDstandardNth2(&Xt1[index]);
    CCTK_REAL_VEC const PDstandardNth3Xt1 = PDstandardNth3(&Xt1[index]);
    CCTK_REAL_VEC const PDstandardNth1Xt2 = PDstandardNth1(&Xt2[index]);
    CCTK_REAL_VEC const PDstandardNth2Xt2 = PDstandardNth2(&Xt2[index]);
    CCTK_REAL_VEC const PDstandardNth3Xt2 = PDstandardNth3(&Xt2[index]);
    CCTK_REAL_VEC const PDstandardNth1Xt3 = PDstandardNth1(&Xt3[index]);
    CCTK_REAL_VEC const PDstandardNth2Xt3 = PDstandardNth2(&Xt3[index]);
    CCTK_REAL_VEC const PDstandardNth3Xt3 = PDstandardNth3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    CCTK_REAL_VEC detgt = ToReal(1);
    
    CCTK_REAL_VEC gtu11 = 
      kdiv(kmsub(gt22L,gt33L,kmul(gt23L,gt23L)),detgt);
    
    CCTK_REAL_VEC gtu12 = 
      kdiv(kmsub(gt13L,gt23L,kmul(gt12L,gt33L)),detgt);
    
    CCTK_REAL_VEC gtu13 = 
      kdiv(kmsub(gt12L,gt23L,kmul(gt13L,gt22L)),detgt);
    
    CCTK_REAL_VEC gtu22 = 
      kdiv(kmsub(gt11L,gt33L,kmul(gt13L,gt13L)),detgt);
    
    CCTK_REAL_VEC gtu23 = 
      kdiv(kmsub(gt12L,gt13L,kmul(gt11L,gt23L)),detgt);
    
    CCTK_REAL_VEC gtu33 = 
      kdiv(kmsub(gt11L,gt22L,kmul(gt12L,gt12L)),detgt);
    
    CCTK_REAL_VEC Gtl111 = kmul(PDstandardNth1gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl112 = kmul(PDstandardNth2gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl113 = kmul(PDstandardNth3gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl122 = 
      kmadd(PDstandardNth1gt22,ToReal(-0.5),PDstandardNth2gt12);
    
    CCTK_REAL_VEC Gtl123 = 
      kmul(kadd(PDstandardNth2gt13,ksub(PDstandardNth3gt12,PDstandardNth1gt23)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl133 = 
      kmadd(PDstandardNth1gt33,ToReal(-0.5),PDstandardNth3gt13);
    
    CCTK_REAL_VEC Gtl211 = 
      kmadd(PDstandardNth2gt11,ToReal(-0.5),PDstandardNth1gt12);
    
    CCTK_REAL_VEC Gtl212 = kmul(PDstandardNth1gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl213 = 
      kmul(kadd(PDstandardNth1gt23,ksub(PDstandardNth3gt12,PDstandardNth2gt13)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl222 = kmul(PDstandardNth2gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl223 = kmul(PDstandardNth3gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl233 = 
      kmadd(PDstandardNth2gt33,ToReal(-0.5),PDstandardNth3gt23);
    
    CCTK_REAL_VEC Gtl311 = 
      kmadd(PDstandardNth3gt11,ToReal(-0.5),PDstandardNth1gt13);
    
    CCTK_REAL_VEC Gtl312 = 
      kmul(kadd(PDstandardNth1gt23,ksub(PDstandardNth2gt13,PDstandardNth3gt12)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl313 = kmul(PDstandardNth1gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl322 = 
      kmadd(PDstandardNth3gt22,ToReal(-0.5),PDstandardNth2gt23);
    
    CCTK_REAL_VEC Gtl323 = kmul(PDstandardNth2gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl333 = kmul(PDstandardNth3gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtlu111 = 
      kmadd(Gtl111,gtu11,kmadd(Gtl112,gtu12,kmul(Gtl113,gtu13)));
    
    CCTK_REAL_VEC Gtlu112 = 
      kmadd(Gtl111,gtu12,kmadd(Gtl112,gtu22,kmul(Gtl113,gtu23)));
    
    CCTK_REAL_VEC Gtlu113 = 
      kmadd(Gtl111,gtu13,kmadd(Gtl112,gtu23,kmul(Gtl113,gtu33)));
    
    CCTK_REAL_VEC Gtlu121 = 
      kmadd(Gtl112,gtu11,kmadd(Gtl122,gtu12,kmul(Gtl123,gtu13)));
    
    CCTK_REAL_VEC Gtlu122 = 
      kmadd(Gtl112,gtu12,kmadd(Gtl122,gtu22,kmul(Gtl123,gtu23)));
    
    CCTK_REAL_VEC Gtlu123 = 
      kmadd(Gtl112,gtu13,kmadd(Gtl122,gtu23,kmul(Gtl123,gtu33)));
    
    CCTK_REAL_VEC Gtlu131 = 
      kmadd(Gtl113,gtu11,kmadd(Gtl123,gtu12,kmul(Gtl133,gtu13)));
    
    CCTK_REAL_VEC Gtlu132 = 
      kmadd(Gtl113,gtu12,kmadd(Gtl123,gtu22,kmul(Gtl133,gtu23)));
    
    CCTK_REAL_VEC Gtlu133 = 
      kmadd(Gtl113,gtu13,kmadd(Gtl123,gtu23,kmul(Gtl133,gtu33)));
    
    CCTK_REAL_VEC Gtlu211 = 
      kmadd(Gtl211,gtu11,kmadd(Gtl212,gtu12,kmul(Gtl213,gtu13)));
    
    CCTK_REAL_VEC Gtlu212 = 
      kmadd(Gtl211,gtu12,kmadd(Gtl212,gtu22,kmul(Gtl213,gtu23)));
    
    CCTK_REAL_VEC Gtlu213 = 
      kmadd(Gtl211,gtu13,kmadd(Gtl212,gtu23,kmul(Gtl213,gtu33)));
    
    CCTK_REAL_VEC Gtlu221 = 
      kmadd(Gtl212,gtu11,kmadd(Gtl222,gtu12,kmul(Gtl223,gtu13)));
    
    CCTK_REAL_VEC Gtlu222 = 
      kmadd(Gtl212,gtu12,kmadd(Gtl222,gtu22,kmul(Gtl223,gtu23)));
    
    CCTK_REAL_VEC Gtlu223 = 
      kmadd(Gtl212,gtu13,kmadd(Gtl222,gtu23,kmul(Gtl223,gtu33)));
    
    CCTK_REAL_VEC Gtlu231 = 
      kmadd(Gtl213,gtu11,kmadd(Gtl223,gtu12,kmul(Gtl233,gtu13)));
    
    CCTK_REAL_VEC Gtlu232 = 
      kmadd(Gtl213,gtu12,kmadd(Gtl223,gtu22,kmul(Gtl233,gtu23)));
    
    CCTK_REAL_VEC Gtlu233 = 
      kmadd(Gtl213,gtu13,kmadd(Gtl223,gtu23,kmul(Gtl233,gtu33)));
    
    CCTK_REAL_VEC Gtlu311 = 
      kmadd(Gtl311,gtu11,kmadd(Gtl312,gtu12,kmul(Gtl313,gtu13)));
    
    CCTK_REAL_VEC Gtlu312 = 
      kmadd(Gtl311,gtu12,kmadd(Gtl312,gtu22,kmul(Gtl313,gtu23)));
    
    CCTK_REAL_VEC Gtlu313 = 
      kmadd(Gtl311,gtu13,kmadd(Gtl312,gtu23,kmul(Gtl313,gtu33)));
    
    CCTK_REAL_VEC Gtlu321 = 
      kmadd(Gtl312,gtu11,kmadd(Gtl322,gtu12,kmul(Gtl323,gtu13)));
    
    CCTK_REAL_VEC Gtlu322 = 
      kmadd(Gtl312,gtu12,kmadd(Gtl322,gtu22,kmul(Gtl323,gtu23)));
    
    CCTK_REAL_VEC Gtlu323 = 
      kmadd(Gtl312,gtu13,kmadd(Gtl322,gtu23,kmul(Gtl323,gtu33)));
    
    CCTK_REAL_VEC Gtlu331 = 
      kmadd(Gtl313,gtu11,kmadd(Gtl323,gtu12,kmul(Gtl333,gtu13)));
    
    CCTK_REAL_VEC Gtlu332 = 
      kmadd(Gtl313,gtu12,kmadd(Gtl323,gtu22,kmul(Gtl333,gtu23)));
    
    CCTK_REAL_VEC Gtlu333 = 
      kmadd(Gtl313,gtu13,kmadd(Gtl323,gtu23,kmul(Gtl333,gtu33)));
    
    CCTK_REAL_VEC Gt111 = 
      kmadd(Gtl111,gtu11,kmadd(Gtl211,gtu12,kmul(Gtl311,gtu13)));
    
    CCTK_REAL_VEC Gt211 = 
      kmadd(Gtl111,gtu12,kmadd(Gtl211,gtu22,kmul(Gtl311,gtu23)));
    
    CCTK_REAL_VEC Gt311 = 
      kmadd(Gtl111,gtu13,kmadd(Gtl211,gtu23,kmul(Gtl311,gtu33)));
    
    CCTK_REAL_VEC Gt112 = 
      kmadd(Gtl112,gtu11,kmadd(Gtl212,gtu12,kmul(Gtl312,gtu13)));
    
    CCTK_REAL_VEC Gt212 = 
      kmadd(Gtl112,gtu12,kmadd(Gtl212,gtu22,kmul(Gtl312,gtu23)));
    
    CCTK_REAL_VEC Gt312 = 
      kmadd(Gtl112,gtu13,kmadd(Gtl212,gtu23,kmul(Gtl312,gtu33)));
    
    CCTK_REAL_VEC Gt113 = 
      kmadd(Gtl113,gtu11,kmadd(Gtl213,gtu12,kmul(Gtl313,gtu13)));
    
    CCTK_REAL_VEC Gt213 = 
      kmadd(Gtl113,gtu12,kmadd(Gtl213,gtu22,kmul(Gtl313,gtu23)));
    
    CCTK_REAL_VEC Gt313 = 
      kmadd(Gtl113,gtu13,kmadd(Gtl213,gtu23,kmul(Gtl313,gtu33)));
    
    CCTK_REAL_VEC Gt122 = 
      kmadd(Gtl122,gtu11,kmadd(Gtl222,gtu12,kmul(Gtl322,gtu13)));
    
    CCTK_REAL_VEC Gt222 = 
      kmadd(Gtl122,gtu12,kmadd(Gtl222,gtu22,kmul(Gtl322,gtu23)));
    
    CCTK_REAL_VEC Gt322 = 
      kmadd(Gtl122,gtu13,kmadd(Gtl222,gtu23,kmul(Gtl322,gtu33)));
    
    CCTK_REAL_VEC Gt123 = 
      kmadd(Gtl123,gtu11,kmadd(Gtl223,gtu12,kmul(Gtl323,gtu13)));
    
    CCTK_REAL_VEC Gt223 = 
      kmadd(Gtl123,gtu12,kmadd(Gtl223,gtu22,kmul(Gtl323,gtu23)));
    
    CCTK_REAL_VEC Gt323 = 
      kmadd(Gtl123,gtu13,kmadd(Gtl223,gtu23,kmul(Gtl323,gtu33)));
    
    CCTK_REAL_VEC Gt133 = 
      kmadd(Gtl133,gtu11,kmadd(Gtl233,gtu12,kmul(Gtl333,gtu13)));
    
    CCTK_REAL_VEC Gt233 = 
      kmadd(Gtl133,gtu12,kmadd(Gtl233,gtu22,kmul(Gtl333,gtu23)));
    
    CCTK_REAL_VEC Gt333 = 
      kmadd(Gtl133,gtu13,kmadd(Gtl233,gtu23,kmul(Gtl333,gtu33)));
    
    CCTK_REAL_VEC Xtn1 = 
      kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmadd(Gt133,gtu33,kmul(kmadd(Gt112,gtu12,kmadd(Gt113,gtu13,kmul(Gt123,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Xtn2 = 
      kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmadd(Gt233,gtu33,kmul(kmadd(Gt212,gtu12,kmadd(Gt213,gtu13,kmul(Gt223,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Xtn3 = 
      kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmadd(Gt333,gtu33,kmul(kmadd(Gt312,gtu12,kmadd(Gt313,gtu13,kmul(Gt323,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Rt11 = 
      kmul(ToReal(0.5),kmadd(gtu12,kmul(PDstandardNth12gt11,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt11,ToReal(-2)),kmadd(gtu23,kmul(PDstandardNth23gt11,ToReal(-2)),knmsub(gtu11,PDstandardNth11gt11,knmsub(gtu22,PDstandardNth22gt11,knmsub(gtu33,PDstandardNth33gt11,kmadd(kmadd(Gt211,Gtlu211,kmadd(Gt212,Gtlu212,kmadd(Gt213,Gtlu213,kmadd(Gt311,Gtlu311,kmadd(Gt312,Gtlu312,kmadd(Gt313,Gtlu313,kmul(gt11L,PDstandardNth1Xt1))))))),ToReal(2),kmadd(gt12L,kmul(PDstandardNth1Xt2,ToReal(2)),kmadd(gt13L,kmul(PDstandardNth1Xt3,ToReal(2)),kmadd(Gtl111,kmul(Xtn1,ToReal(2)),kmadd(Gtl112,kmul(Xtn2,ToReal(2)),kmadd(Gtl113,kmul(Xtn3,ToReal(2)),kmadd(kmadd(Gt211,Gtlu121,kmadd(Gt212,Gtlu122,kmadd(Gt213,Gtlu123,kmadd(Gt311,Gtlu131,kmadd(Gt312,Gtlu132,kmul(Gt313,Gtlu133)))))),ToReal(4),kmul(kmadd(Gt111,Gtlu111,kmadd(Gt112,Gtlu112,kmul(Gt113,Gtlu113))),ToReal(6))))))))))))))));
    
    CCTK_REAL_VEC Rt12 = 
      kmul(ToReal(0.5),kmadd(gt12L,PDstandardNth1Xt1,kmadd(gt22L,PDstandardNth1Xt2,kmadd(gt23L,PDstandardNth1Xt3,kmadd(gt11L,PDstandardNth2Xt1,kmadd(gt12L,PDstandardNth2Xt2,kmadd(gt13L,PDstandardNth2Xt3,kmadd(Gtl112,Xtn1,kmadd(Gtl211,Xtn1,kmadd(Gtl122,Xtn2,kmadd(Gtl212,Xtn2,kmadd(Gtl123,Xtn3,kmadd(Gtl213,Xtn3,kmadd(gtu12,kmul(PDstandardNth12gt12,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt12,ToReal(-2)),kmadd(gtu23,kmul(PDstandardNth23gt12,ToReal(-2)),knmsub(gtu11,PDstandardNth11gt12,knmsub(gtu22,PDstandardNth22gt12,knmsub(gtu33,PDstandardNth33gt12,kmadd(kmadd(Gt112,Gtlu111,kmadd(Gt122,Gtlu112,kmadd(Gt123,Gtlu113,kmadd(Gt111,Gtlu121,kmadd(Gt212,Gtlu121,kmadd(Gt112,Gtlu122,kmadd(Gt222,Gtlu122,kmadd(Gt113,Gtlu123,kmadd(Gt223,Gtlu123,kmadd(Gt312,Gtlu131,kmadd(Gt322,Gtlu132,kmadd(Gt323,Gtlu133,kmadd(Gt111,Gtlu211,kmadd(Gt112,Gtlu212,kmadd(Gt113,Gtlu213,kmadd(Gt311,Gtlu231,kmadd(Gt312,Gtlu232,kmadd(Gt313,Gtlu233,kmadd(Gt311,Gtlu321,kmadd(Gt312,Gtlu322,kmul(Gt313,Gtlu323))))))))))))))))))))),ToReal(2),kmul(kmadd(Gt211,Gtlu221,kmadd(Gt212,Gtlu222,kmul(Gt213,Gtlu223))),ToReal(4))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt13 = 
      kmul(ToReal(0.5),kmadd(gt13L,PDstandardNth1Xt1,kmadd(gt23L,PDstandardNth1Xt2,kmadd(gt33L,PDstandardNth1Xt3,kmadd(gt11L,PDstandardNth3Xt1,kmadd(gt12L,PDstandardNth3Xt2,kmadd(gt13L,PDstandardNth3Xt3,kmadd(Gtl113,Xtn1,kmadd(Gtl311,Xtn1,kmadd(Gtl123,Xtn2,kmadd(Gtl312,Xtn2,kmadd(Gtl133,Xtn3,kmadd(Gtl313,Xtn3,kmadd(gtu12,kmul(PDstandardNth12gt13,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt13,ToReal(-2)),kmadd(gtu23,kmul(PDstandardNth23gt13,ToReal(-2)),knmsub(gtu11,PDstandardNth11gt13,knmsub(gtu22,PDstandardNth22gt13,knmsub(gtu33,PDstandardNth33gt13,kmadd(kmadd(Gt113,Gtlu111,kmadd(Gt123,Gtlu112,kmadd(Gt133,Gtlu113,kmadd(Gt213,Gtlu121,kmadd(Gt223,Gtlu122,kmadd(Gt233,Gtlu123,kmadd(Gt111,Gtlu131,kmadd(Gt313,Gtlu131,kmadd(Gt112,Gtlu132,kmadd(Gt323,Gtlu132,kmadd(Gt113,Gtlu133,kmadd(Gt333,Gtlu133,kmadd(Gt211,Gtlu231,kmadd(Gt212,Gtlu232,kmadd(Gt213,Gtlu233,kmadd(Gt111,Gtlu311,kmadd(Gt112,Gtlu312,kmadd(Gt113,Gtlu313,kmadd(Gt211,Gtlu321,kmadd(Gt212,Gtlu322,kmul(Gt213,Gtlu323))))))))))))))))))))),ToReal(2),kmul(kmadd(Gt311,Gtlu331,kmadd(Gt312,Gtlu332,kmul(Gt313,Gtlu333))),ToReal(4))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt22 = 
      kmul(ToReal(0.5),kmadd(gtu12,kmul(PDstandardNth12gt22,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt22,ToReal(-2)),kmadd(gtu23,kmul(PDstandardNth23gt22,ToReal(-2)),knmsub(gtu11,PDstandardNth11gt22,knmsub(gtu22,PDstandardNth22gt22,knmsub(gtu33,PDstandardNth33gt22,kmadd(kmadd(Gt112,Gtlu121,kmadd(Gt122,Gtlu122,kmadd(Gt123,Gtlu123,kmadd(Gt312,Gtlu321,kmadd(Gt322,Gtlu322,kmadd(Gt323,Gtlu323,kmul(gt12L,PDstandardNth2Xt1))))))),ToReal(2),kmadd(gt22L,kmul(PDstandardNth2Xt2,ToReal(2)),kmadd(gt23L,kmul(PDstandardNth2Xt3,ToReal(2)),kmadd(Gtl212,kmul(Xtn1,ToReal(2)),kmadd(Gtl222,kmul(Xtn2,ToReal(2)),kmadd(Gtl223,kmul(Xtn3,ToReal(2)),kmadd(kmadd(Gt112,Gtlu211,kmadd(Gt122,Gtlu212,kmadd(Gt123,Gtlu213,kmadd(Gt312,Gtlu231,kmadd(Gt322,Gtlu232,kmul(Gt323,Gtlu233)))))),ToReal(4),kmul(kmadd(Gt212,Gtlu221,kmadd(Gt222,Gtlu222,kmul(Gt223,Gtlu223))),ToReal(6))))))))))))))));
    
    CCTK_REAL_VEC Rt23 = 
      kmul(ToReal(0.5),kmadd(gt13L,PDstandardNth2Xt1,kmadd(gt23L,PDstandardNth2Xt2,kmadd(gt33L,PDstandardNth2Xt3,kmadd(gt12L,PDstandardNth3Xt1,kmadd(gt22L,PDstandardNth3Xt2,kmadd(gt23L,PDstandardNth3Xt3,kmadd(Gtl213,Xtn1,kmadd(Gtl312,Xtn1,kmadd(Gtl223,Xtn2,kmadd(Gtl322,Xtn2,kmadd(Gtl233,Xtn3,kmadd(Gtl323,Xtn3,kmadd(gtu12,kmul(PDstandardNth12gt23,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt23,ToReal(-2)),kmadd(gtu23,kmul(PDstandardNth23gt23,ToReal(-2)),knmsub(gtu11,PDstandardNth11gt23,knmsub(gtu22,PDstandardNth22gt23,knmsub(gtu33,PDstandardNth33gt23,kmadd(kmadd(Gt112,Gtlu131,kmadd(Gt122,Gtlu132,kmadd(Gt123,Gtlu133,kmadd(Gt113,Gtlu211,kmadd(Gt123,Gtlu212,kmadd(Gt133,Gtlu213,kmadd(Gt213,Gtlu221,kmadd(Gt223,Gtlu222,kmadd(Gt233,Gtlu223,kmadd(Gt212,Gtlu231,kmadd(Gt313,Gtlu231,kmadd(Gt222,Gtlu232,kmadd(Gt323,Gtlu232,kmadd(Gt223,Gtlu233,kmadd(Gt333,Gtlu233,kmadd(Gt112,Gtlu311,kmadd(Gt122,Gtlu312,kmadd(Gt123,Gtlu313,kmadd(Gt212,Gtlu321,kmadd(Gt222,Gtlu322,kmul(Gt223,Gtlu323))))))))))))))))))))),ToReal(2),kmul(kmadd(Gt312,Gtlu331,kmadd(Gt322,Gtlu332,kmul(Gt323,Gtlu333))),ToReal(4))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt33 = 
      kmul(ToReal(0.5),kmadd(gtu12,kmul(PDstandardNth12gt33,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt33,ToReal(-2)),kmadd(gtu23,kmul(PDstandardNth23gt33,ToReal(-2)),knmsub(gtu11,PDstandardNth11gt33,knmsub(gtu22,PDstandardNth22gt33,knmsub(gtu33,PDstandardNth33gt33,kmadd(kmadd(Gt113,Gtlu131,kmadd(Gt123,Gtlu132,kmadd(Gt133,Gtlu133,kmadd(Gt213,Gtlu231,kmadd(Gt223,Gtlu232,kmadd(Gt233,Gtlu233,kmul(gt13L,PDstandardNth3Xt1))))))),ToReal(2),kmadd(gt23L,kmul(PDstandardNth3Xt2,ToReal(2)),kmadd(gt33L,kmul(PDstandardNth3Xt3,ToReal(2)),kmadd(Gtl313,kmul(Xtn1,ToReal(2)),kmadd(Gtl323,kmul(Xtn2,ToReal(2)),kmadd(Gtl333,kmul(Xtn3,ToReal(2)),kmadd(kmadd(Gt113,Gtlu311,kmadd(Gt123,Gtlu312,kmadd(Gt133,Gtlu313,kmadd(Gt213,Gtlu321,kmadd(Gt223,Gtlu322,kmul(Gt233,Gtlu323)))))),ToReal(4),kmul(kmadd(Gt313,Gtlu331,kmadd(Gt323,Gtlu332,kmul(Gt333,Gtlu333))),ToReal(6))))))))))))))));
    
    CCTK_REAL_VEC fac1 = 
      IfThen(conformalMethod,kdiv(ToReal(-0.5),phiL),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 = kmul(fac1,PDstandardNth1phi);
    
    CCTK_REAL_VEC cdphi2 = kmul(fac1,PDstandardNth2phi);
    
    CCTK_REAL_VEC cdphi3 = kmul(fac1,PDstandardNth3phi);
    
    CCTK_REAL_VEC fac2 = 
      IfThen(conformalMethod,kdiv(ToReal(0.5),kmul(phiL,phiL)),ToReal(0));
    
    CCTK_REAL_VEC cdphi211 = 
      kmsub(fac2,kmul(PDstandardNth1phi,PDstandardNth1phi),kmul(fac1,kmadd(Gt111,PDstandardNth1phi,kmadd(Gt211,PDstandardNth2phi,kmsub(Gt311,PDstandardNth3phi,PDstandardNth11phi)))));
    
    CCTK_REAL_VEC cdphi212 = 
      kmsub(fac2,kmul(PDstandardNth1phi,PDstandardNth2phi),kmul(fac1,kmadd(Gt112,PDstandardNth1phi,kmadd(Gt212,PDstandardNth2phi,kmsub(Gt312,PDstandardNth3phi,PDstandardNth12phi)))));
    
    CCTK_REAL_VEC cdphi213 = 
      kmsub(fac2,kmul(PDstandardNth1phi,PDstandardNth3phi),kmul(fac1,kmadd(Gt113,PDstandardNth1phi,kmadd(Gt213,PDstandardNth2phi,kmsub(Gt313,PDstandardNth3phi,PDstandardNth13phi)))));
    
    CCTK_REAL_VEC cdphi222 = 
      kmsub(fac2,kmul(PDstandardNth2phi,PDstandardNth2phi),kmul(fac1,kmadd(Gt122,PDstandardNth1phi,kmadd(Gt222,PDstandardNth2phi,kmsub(Gt322,PDstandardNth3phi,PDstandardNth22phi)))));
    
    CCTK_REAL_VEC cdphi223 = 
      kmsub(fac2,kmul(PDstandardNth2phi,PDstandardNth3phi),kmul(fac1,kmadd(Gt123,PDstandardNth1phi,kmadd(Gt223,PDstandardNth2phi,kmsub(Gt323,PDstandardNth3phi,PDstandardNth23phi)))));
    
    CCTK_REAL_VEC cdphi233 = 
      kmsub(fac2,kmul(PDstandardNth3phi,PDstandardNth3phi),kmul(fac1,kmadd(Gt133,PDstandardNth1phi,kmadd(Gt233,PDstandardNth2phi,kmsub(Gt333,PDstandardNth3phi,PDstandardNth33phi)))));
    
    CCTK_REAL_VEC Rphi11 = 
      kmul(ToReal(-2),kadd(cdphi211,kmadd(kmul(cdphi1,cdphi1),kmul(kmadd(gt11L,gtu11,ToReal(-1)),ToReal(2)),kmul(gt11L,kmadd(cdphi211,gtu11,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu33,kmul(cdphi3,cdphi3))))),ToReal(2),kmadd(gtu22,kmadd(kmul(cdphi2,cdphi2),ToReal(2),cdphi222),kmul(kmadd(cdphi1,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13)),kmul(cdphi2,kmul(cdphi3,gtu23))),ToReal(4))))))))));
    
    CCTK_REAL_VEC Rphi12 = 
      kmul(ToReal(-2),kadd(cdphi212,kmadd(gt12L,kmadd(cdphi211,gtu11,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu11,kmul(cdphi1,cdphi1))))),ToReal(2),kmadd(gtu22,kmadd(kmul(cdphi2,cdphi2),ToReal(2),cdphi222),kmadd(gtu33,kmadd(kmul(cdphi3,cdphi3),ToReal(2),cdphi233),kmul(cdphi2,kmul(cdphi3,kmul(gtu23,ToReal(4)))))))),kmul(cdphi1,kmadd(gt12L,kmul(cdphi3,kmul(gtu13,ToReal(4))),kmul(cdphi2,kmadd(gt12L,kmul(gtu12,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi13 = 
      kmul(ToReal(-2),kadd(cdphi213,kmadd(gt13L,kmadd(cdphi211,gtu11,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu11,kmul(cdphi1,cdphi1))))),ToReal(2),kmadd(gtu22,kmadd(kmul(cdphi2,cdphi2),ToReal(2),cdphi222),kmadd(gtu33,kmadd(kmul(cdphi3,cdphi3),ToReal(2),cdphi233),kmul(cdphi2,kmul(cdphi3,kmul(gtu23,ToReal(4)))))))),kmul(cdphi1,kmadd(gt13L,kmul(cdphi2,kmul(gtu12,ToReal(4))),kmul(cdphi3,kmadd(gt13L,kmul(gtu13,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi22 = 
      kmul(ToReal(-2),kadd(cdphi222,kmadd(kmul(cdphi2,cdphi2),kmul(kmadd(gt22L,gtu22,ToReal(-1)),ToReal(2)),kmul(gt22L,kmadd(cdphi222,gtu22,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu33,kmul(cdphi3,cdphi3))))),ToReal(2),kmadd(gtu11,kmadd(kmul(cdphi1,cdphi1),ToReal(2),cdphi211),kmul(kmadd(cdphi1,kmul(cdphi3,gtu13),kmul(cdphi2,kmadd(cdphi1,gtu12,kmul(cdphi3,gtu23)))),ToReal(4))))))))));
    
    CCTK_REAL_VEC Rphi23 = 
      kmul(ToReal(-2),kadd(cdphi223,kmadd(gt23L,kmadd(cdphi222,gtu22,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu22,kmul(cdphi2,cdphi2))))),ToReal(2),kmadd(gtu11,kmadd(kmul(cdphi1,cdphi1),ToReal(2),cdphi211),kmadd(gtu33,kmadd(kmul(cdphi3,cdphi3),ToReal(2),cdphi233),kmul(cdphi1,kmul(cdphi3,kmul(gtu13,ToReal(4)))))))),kmul(cdphi2,kmadd(gt23L,kmul(cdphi1,kmul(gtu12,ToReal(4))),kmul(cdphi3,kmadd(gt23L,kmul(gtu23,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi33 = 
      kmul(ToReal(-2),kadd(cdphi233,kmadd(kmul(cdphi3,cdphi3),kmul(kmadd(gt33L,gtu33,ToReal(-1)),ToReal(2)),kmul(gt33L,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi213,gtu13,kmul(cdphi223,gtu23)),ToReal(2),kmadd(gtu11,kmadd(kmul(cdphi1,cdphi1),ToReal(2),cdphi211),kmadd(gtu22,kmadd(kmul(cdphi2,cdphi2),ToReal(2),cdphi222),kmadd(cdphi3,kmul(kmadd(cdphi1,gtu13,kmul(cdphi2,gtu23)),ToReal(4)),kmul(gtu12,kmadd(cdphi212,ToReal(2),kmul(cdphi1,kmul(cdphi2,ToReal(4))))))))))))));
    
    CCTK_REAL_VEC Atm11 = 
      kmadd(At11L,gtu11,kmadd(At12L,gtu12,kmul(At13L,gtu13)));
    
    CCTK_REAL_VEC Atm21 = 
      kmadd(At11L,gtu12,kmadd(At12L,gtu22,kmul(At13L,gtu23)));
    
    CCTK_REAL_VEC Atm31 = 
      kmadd(At11L,gtu13,kmadd(At12L,gtu23,kmul(At13L,gtu33)));
    
    CCTK_REAL_VEC Atm12 = 
      kmadd(At12L,gtu11,kmadd(At22L,gtu12,kmul(At23L,gtu13)));
    
    CCTK_REAL_VEC Atm22 = 
      kmadd(At12L,gtu12,kmadd(At22L,gtu22,kmul(At23L,gtu23)));
    
    CCTK_REAL_VEC Atm32 = 
      kmadd(At12L,gtu13,kmadd(At22L,gtu23,kmul(At23L,gtu33)));
    
    CCTK_REAL_VEC Atm13 = 
      kmadd(At13L,gtu11,kmadd(At23L,gtu12,kmul(At33L,gtu13)));
    
    CCTK_REAL_VEC Atm23 = 
      kmadd(At13L,gtu12,kmadd(At23L,gtu22,kmul(At33L,gtu23)));
    
    CCTK_REAL_VEC Atm33 = 
      kmadd(At13L,gtu13,kmadd(At23L,gtu23,kmul(At33L,gtu33)));
    
    CCTK_REAL_VEC e4phi = 
      IfThen(conformalMethod,kdiv(ToReal(1),kmul(phiL,phiL)),kexp(kmul(phiL,ToReal(4))));
    
    CCTK_REAL_VEC em4phi = kdiv(ToReal(1),e4phi);
    
    CCTK_REAL_VEC g11 = kmul(gt11L,e4phi);
    
    CCTK_REAL_VEC g12 = kmul(gt12L,e4phi);
    
    CCTK_REAL_VEC g13 = kmul(gt13L,e4phi);
    
    CCTK_REAL_VEC g22 = kmul(gt22L,e4phi);
    
    CCTK_REAL_VEC g23 = kmul(gt23L,e4phi);
    
    CCTK_REAL_VEC g33 = kmul(gt33L,e4phi);
    
    CCTK_REAL_VEC gu11 = kmul(em4phi,gtu11);
    
    CCTK_REAL_VEC gu12 = kmul(em4phi,gtu12);
    
    CCTK_REAL_VEC gu13 = kmul(em4phi,gtu13);
    
    CCTK_REAL_VEC gu22 = kmul(em4phi,gtu22);
    
    CCTK_REAL_VEC gu23 = kmul(em4phi,gtu23);
    
    CCTK_REAL_VEC gu33 = kmul(em4phi,gtu33);
    
    CCTK_REAL_VEC R11 = kadd(Rphi11,Rt11);
    
    CCTK_REAL_VEC R12 = kadd(Rphi12,Rt12);
    
    CCTK_REAL_VEC R13 = kadd(Rphi13,Rt13);
    
    CCTK_REAL_VEC R22 = kadd(Rphi22,Rt22);
    
    CCTK_REAL_VEC R23 = kadd(Rphi23,Rt23);
    
    CCTK_REAL_VEC R33 = kadd(Rphi33,Rt33);
    
    CCTK_REAL_VEC Ats11 = 
      kmadd(Gt211,PDstandardNth2alpha,kmadd(Gt311,PDstandardNth3alpha,kmadd(alphaL,R11,kmsub(PDstandardNth1alpha,kmadd(cdphi1,ToReal(4),Gt111),PDstandardNth11alpha))));
    
    CCTK_REAL_VEC Ats12 = 
      kmadd(Gt312,PDstandardNth3alpha,kmadd(alphaL,R12,ksub(kmadd(PDstandardNth2alpha,kmadd(cdphi1,ToReal(2),Gt212),kmul(PDstandardNth1alpha,kmadd(cdphi2,ToReal(2),Gt112))),PDstandardNth12alpha)));
    
    CCTK_REAL_VEC Ats13 = 
      kmadd(Gt213,PDstandardNth2alpha,kmadd(alphaL,R13,ksub(kmadd(PDstandardNth3alpha,kmadd(cdphi1,ToReal(2),Gt313),kmul(PDstandardNth1alpha,kmadd(cdphi3,ToReal(2),Gt113))),PDstandardNth13alpha)));
    
    CCTK_REAL_VEC Ats22 = 
      kmadd(Gt122,PDstandardNth1alpha,kmadd(Gt322,PDstandardNth3alpha,kmadd(alphaL,R22,kmsub(PDstandardNth2alpha,kmadd(cdphi2,ToReal(4),Gt222),PDstandardNth22alpha))));
    
    CCTK_REAL_VEC Ats23 = 
      kmadd(Gt123,PDstandardNth1alpha,kmadd(alphaL,R23,ksub(kmadd(PDstandardNth3alpha,kmadd(cdphi2,ToReal(2),Gt323),kmul(PDstandardNth2alpha,kmadd(cdphi3,ToReal(2),Gt223))),PDstandardNth23alpha)));
    
    CCTK_REAL_VEC Ats33 = 
      kmadd(Gt133,PDstandardNth1alpha,kmadd(Gt233,PDstandardNth2alpha,kmadd(alphaL,R33,kmsub(PDstandardNth3alpha,kmadd(cdphi3,ToReal(4),Gt333),PDstandardNth33alpha))));
    
    CCTK_REAL_VEC trAts = 
      kmadd(Ats11,gu11,kmadd(Ats22,gu22,kmadd(Ats33,gu33,kmul(kmadd(Ats12,gu12,kmadd(Ats13,gu13,kmul(Ats23,gu23))),ToReal(2)))));
    
    CCTK_REAL_VEC At11rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(em4phi,kmsub(Ats11,ToReal(3),kmul(g11,trAts)),kmadd(alphaL,kmadd(kmadd(At12L,Atm21,kmul(At13L,Atm31)),ToReal(-6),kmul(At11L,kmadd(Atm11,ToReal(-6),kmul(trKL,ToReal(3))))),kmadd(At11L,kmadd(kadd(PDstandardNth2beta2,PDstandardNth3beta3),ToReal(-2),kmul(PDstandardNth1beta1,ToReal(4))),kmul(kmadd(At12L,PDstandardNth1beta2,kmul(At13L,PDstandardNth1beta3)),ToReal(6))))));
    
    CCTK_REAL_VEC At12rhsL = 
      kmadd(At22L,PDstandardNth1beta2,kmadd(At23L,PDstandardNth1beta3,kmadd(At11L,PDstandardNth2beta1,kmadd(At13L,PDstandardNth2beta3,kmadd(alphaL,kmadd(kmadd(At11L,Atm12,kmul(At13L,Atm32)),ToReal(-2),kmul(At12L,kmadd(Atm22,ToReal(-2),trKL))),kmadd(em4phi,kmadd(g12,kmul(trAts,ToReal(-0.333333333333333333333333333333)),Ats12),kmul(At12L,kmadd(PDstandardNth3beta3,ToReal(-0.666666666666666666666666666667),kmul(kadd(PDstandardNth1beta1,PDstandardNth2beta2),ToReal(0.333333333333333333333333333333))))))))));
    
    CCTK_REAL_VEC At13rhsL = 
      kmadd(At23L,PDstandardNth1beta2,kmadd(At33L,PDstandardNth1beta3,kmadd(At11L,PDstandardNth3beta1,kmadd(At12L,PDstandardNth3beta2,kmadd(alphaL,kmadd(kmadd(At11L,Atm13,kmul(At12L,Atm23)),ToReal(-2),kmul(At13L,kmadd(Atm33,ToReal(-2),trKL))),kmadd(em4phi,kmadd(g13,kmul(trAts,ToReal(-0.333333333333333333333333333333)),Ats13),kmul(At13L,kmadd(PDstandardNth2beta2,ToReal(-0.666666666666666666666666666667),kmul(kadd(PDstandardNth1beta1,PDstandardNth3beta3),ToReal(0.333333333333333333333333333333))))))))));
    
    CCTK_REAL_VEC At22rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(em4phi,kmsub(Ats22,ToReal(3),kmul(g22,trAts)),kmadd(alphaL,kmadd(kmadd(At12L,Atm12,kmul(At23L,Atm32)),ToReal(-6),kmul(At22L,kmadd(Atm22,ToReal(-6),kmul(trKL,ToReal(3))))),kmadd(At22L,kmadd(kadd(PDstandardNth1beta1,PDstandardNth3beta3),ToReal(-2),kmul(PDstandardNth2beta2,ToReal(4))),kmul(kmadd(At12L,PDstandardNth2beta1,kmul(At23L,PDstandardNth2beta3)),ToReal(6))))));
    
    CCTK_REAL_VEC At23rhsL = 
      kmadd(At13L,PDstandardNth2beta1,kmadd(At33L,PDstandardNth2beta3,kmadd(At12L,PDstandardNth3beta1,kmadd(At22L,PDstandardNth3beta2,kmadd(alphaL,kmadd(kmadd(At12L,Atm13,kmul(At22L,Atm23)),ToReal(-2),kmul(At23L,kmadd(Atm33,ToReal(-2),trKL))),kmadd(em4phi,kmadd(g23,kmul(trAts,ToReal(-0.333333333333333333333333333333)),Ats23),kmul(At23L,kmadd(PDstandardNth1beta1,ToReal(-0.666666666666666666666666666667),kmul(kadd(PDstandardNth2beta2,PDstandardNth3beta3),ToReal(0.333333333333333333333333333333))))))))));
    
    CCTK_REAL_VEC At33rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(em4phi,kmsub(Ats33,ToReal(3),kmul(g33,trAts)),kmadd(alphaL,kmadd(kmadd(At13L,Atm13,kmul(At23L,Atm23)),ToReal(-6),kmul(At33L,kmadd(Atm33,ToReal(-6),kmul(trKL,ToReal(3))))),kmadd(At33L,kmadd(kadd(PDstandardNth1beta1,PDstandardNth2beta2),ToReal(-2),kmul(PDstandardNth3beta3,ToReal(4))),kmul(kmadd(At13L,PDstandardNth3beta1,kmul(At23L,PDstandardNth3beta2)),ToReal(6))))));
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,lc_imin,lc_imax);
    vec_store_nta_partial(At11rhs[index],At11rhsL);
    vec_store_nta_partial(At12rhs[index],At12rhsL);
    vec_store_nta_partial(At13rhs[index],At13rhsL);
    vec_store_nta_partial(At22rhs[index],At22rhsL);
    vec_store_nta_partial(At23rhs[index],At23rhsL);
    vec_store_nta_partial(At33rhs[index],At33rhsL);
  }
  LC_ENDLOOP3VEC(ML_BSSN_Host_RHS2);
}

extern "C" void ML_BSSN_Host_RHS2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_Host_RHS2_Body");
  }
  
  if (cctk_iteration % ML_BSSN_Host_RHS2_calc_every != ML_BSSN_Host_RHS2_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "ML_BSSN_Host::ML_curv",
    "ML_BSSN_Host::ML_curvrhs",
    "ML_BSSN_Host::ML_Gamma",
    "ML_BSSN_Host::ML_lapse",
    "ML_BSSN_Host::ML_log_confac",
    "ML_BSSN_Host::ML_metric",
    "ML_BSSN_Host::ML_shift",
    "ML_BSSN_Host::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_RHS2", 8, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Host_RHS2", 4, 4, 4);
  
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_Host_RHS2_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_RHS2_Body");
  }
}
