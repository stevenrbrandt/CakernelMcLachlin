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
#include "vectors.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define QAD(x) (SQR(SQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))

extern "C" void ML_BSSN_O2_constraints1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O2::ML_Ham","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O2::ML_Ham.");
  return;
}

static void ML_BSSN_O2_constraints1_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O2_constraints1_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O2_constraints1_calc_every != ML_BSSN_O2_constraints1_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_BSSN_O2::ML_curv","ML_BSSN_O2::ML_Gamma","ML_BSSN_O2::ML_Ham","ML_BSSN_O2::ML_lapse","ML_BSSN_O2::ML_log_confac","ML_BSSN_O2::ML_metric","ML_BSSN_O2::ML_shift","ML_BSSN_O2::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_O2_constraints1", 8, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_O2_constraints1", 1, 1, 1);
  
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
  CCTK_REAL_VEC const p1o16dx = kmul(INV(dx),ToReal(0.0625));
  CCTK_REAL_VEC const p1o16dy = kmul(INV(dy),ToReal(0.0625));
  CCTK_REAL_VEC const p1o16dz = kmul(INV(dz),ToReal(0.0625));
  CCTK_REAL_VEC const p1o2dx = kmul(INV(dx),ToReal(0.5));
  CCTK_REAL_VEC const p1o2dy = kmul(INV(dy),ToReal(0.5));
  CCTK_REAL_VEC const p1o2dz = kmul(INV(dz),ToReal(0.5));
  CCTK_REAL_VEC const p1o4dx = kmul(INV(dx),ToReal(0.25));
  CCTK_REAL_VEC const p1o4dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(0.25)));
  CCTK_REAL_VEC const p1o4dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(0.25)));
  CCTK_REAL_VEC const p1o4dy = kmul(INV(dy),ToReal(0.25));
  CCTK_REAL_VEC const p1o4dydz = kmul(INV(dy),kmul(INV(dz),ToReal(0.25)));
  CCTK_REAL_VEC const p1o4dz = kmul(INV(dz),ToReal(0.25));
  CCTK_REAL_VEC const p1odx = INV(dx);
  CCTK_REAL_VEC const p1odx2 = INV(SQR(dx));
  CCTK_REAL_VEC const p1ody = INV(dy);
  CCTK_REAL_VEC const p1ody2 = INV(SQR(dy));
  CCTK_REAL_VEC const p1odz = INV(dz);
  CCTK_REAL_VEC const p1odz2 = INV(SQR(dz));
  CCTK_REAL_VEC const pm1o2dx = kmul(INV(dx),ToReal(-0.5));
  CCTK_REAL_VEC const pm1o2dy = kmul(INV(dy),ToReal(-0.5));
  CCTK_REAL_VEC const pm1o2dz = kmul(INV(dz),ToReal(-0.5));
  CCTK_REAL_VEC const pm1o4dx = kmul(INV(dx),ToReal(-0.25));
  CCTK_REAL_VEC const pm1o4dy = kmul(INV(dy),ToReal(-0.25));
  CCTK_REAL_VEC const pm1o4dz = kmul(INV(dz),ToReal(-0.25));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC (ML_BSSN_O2_constraints1,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
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
    
    CCTK_REAL_VEC eTttL, eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL;
    
    if (*stress_energy_state)
    {
      eTttL = vec_load(eTtt[index]);
      eTtxL = vec_load(eTtx[index]);
      eTtyL = vec_load(eTty[index]);
      eTtzL = vec_load(eTtz[index]);
      eTxxL = vec_load(eTxx[index]);
      eTxyL = vec_load(eTxy[index]);
      eTxzL = vec_load(eTxz[index]);
      eTyyL = vec_load(eTyy[index]);
      eTyzL = vec_load(eTyz[index]);
      eTzzL = vec_load(eTzz[index]);
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
    CCTK_REAL_VEC detgt = ToReal(1);
    
    CCTK_REAL_VEC gtu11 = kmul(INV(detgt),kmsub(gt22L,gt33L,SQR(gt23L)));
    
    CCTK_REAL_VEC gtu12 = 
      kmul(INV(detgt),kmsub(gt13L,gt23L,kmul(gt12L,gt33L)));
    
    CCTK_REAL_VEC gtu13 = 
      kmul(INV(detgt),kmsub(gt12L,gt23L,kmul(gt13L,gt22L)));
    
    CCTK_REAL_VEC gtu22 = kmul(INV(detgt),kmsub(gt11L,gt33L,SQR(gt13L)));
    
    CCTK_REAL_VEC gtu23 = 
      kmul(INV(detgt),kmsub(gt12L,gt13L,kmul(gt11L,gt23L)));
    
    CCTK_REAL_VEC gtu33 = kmul(INV(detgt),kmsub(gt11L,gt22L,SQR(gt12L)));
    
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
      IfThen(conformalMethod,kmul(INV(phiL),ToReal(-0.5)),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 = kmul(fac1,PDstandardNth1phi);
    
    CCTK_REAL_VEC cdphi2 = kmul(fac1,PDstandardNth2phi);
    
    CCTK_REAL_VEC cdphi3 = kmul(fac1,PDstandardNth3phi);
    
    CCTK_REAL_VEC fac2 = 
      IfThen(conformalMethod,kmul(INV(SQR(phiL)),ToReal(0.5)),ToReal(0));
    
    CCTK_REAL_VEC cdphi211 = 
      kmsub(fac2,SQR(PDstandardNth1phi),kmul(fac1,kmadd(Gt111,PDstandardNth1phi,kmadd(Gt211,PDstandardNth2phi,kmsub(Gt311,PDstandardNth3phi,PDstandardNth11phi)))));
    
    CCTK_REAL_VEC cdphi212 = 
      kmsub(fac2,kmul(PDstandardNth1phi,PDstandardNth2phi),kmul(fac1,kmadd(Gt112,PDstandardNth1phi,kmadd(Gt212,PDstandardNth2phi,kmsub(Gt312,PDstandardNth3phi,PDstandardNth12phi)))));
    
    CCTK_REAL_VEC cdphi213 = 
      kmsub(fac2,kmul(PDstandardNth1phi,PDstandardNth3phi),kmul(fac1,kmadd(Gt113,PDstandardNth1phi,kmadd(Gt213,PDstandardNth2phi,kmsub(Gt313,PDstandardNth3phi,PDstandardNth13phi)))));
    
    CCTK_REAL_VEC cdphi222 = 
      kmsub(fac2,SQR(PDstandardNth2phi),kmul(fac1,kmadd(Gt122,PDstandardNth1phi,kmadd(Gt222,PDstandardNth2phi,kmsub(Gt322,PDstandardNth3phi,PDstandardNth22phi)))));
    
    CCTK_REAL_VEC cdphi223 = 
      kmsub(fac2,kmul(PDstandardNth2phi,PDstandardNth3phi),kmul(fac1,kmadd(Gt123,PDstandardNth1phi,kmadd(Gt223,PDstandardNth2phi,kmsub(Gt323,PDstandardNth3phi,PDstandardNth23phi)))));
    
    CCTK_REAL_VEC cdphi233 = 
      kmsub(fac2,SQR(PDstandardNth3phi),kmul(fac1,kmadd(Gt133,PDstandardNth1phi,kmadd(Gt233,PDstandardNth2phi,kmsub(Gt333,PDstandardNth3phi,PDstandardNth33phi)))));
    
    CCTK_REAL_VEC Rphi11 = 
      kmul(ToReal(-2),kadd(cdphi211,kmadd(SQR(cdphi1),kmul(kmadd(gt11L,gtu11,ToReal(-1)),ToReal(2)),kmul(gt11L,kmadd(cdphi211,gtu11,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu33,SQR(cdphi3))))),ToReal(2),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmul(kmadd(cdphi1,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13)),kmul(cdphi2,kmul(cdphi3,gtu23))),ToReal(4))))))))));
    
    CCTK_REAL_VEC Rphi12 = 
      kmul(ToReal(-2),kadd(cdphi212,kmadd(gt12L,kmadd(cdphi211,gtu11,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu11,SQR(cdphi1))))),ToReal(2),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmadd(gtu33,kmadd(SQR(cdphi3),ToReal(2),cdphi233),kmul(cdphi2,kmul(cdphi3,kmul(gtu23,ToReal(4)))))))),kmul(cdphi1,kmadd(cdphi3,kmul(gt12L,kmul(gtu13,ToReal(4))),kmul(cdphi2,kmadd(gt12L,kmul(gtu12,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi13 = 
      kmul(ToReal(-2),kadd(cdphi213,kmadd(gt13L,kmadd(cdphi211,gtu11,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu11,SQR(cdphi1))))),ToReal(2),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmadd(gtu33,kmadd(SQR(cdphi3),ToReal(2),cdphi233),kmul(cdphi2,kmul(cdphi3,kmul(gtu23,ToReal(4)))))))),kmul(cdphi1,kmadd(cdphi2,kmul(gt13L,kmul(gtu12,ToReal(4))),kmul(cdphi3,kmadd(gt13L,kmul(gtu13,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi22 = 
      kmul(ToReal(-2),kadd(cdphi222,kmadd(SQR(cdphi2),kmul(kmadd(gt22L,gtu22,ToReal(-1)),ToReal(2)),kmul(gt22L,kmadd(cdphi222,gtu22,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu33,SQR(cdphi3))))),ToReal(2),kmadd(gtu11,kmadd(SQR(cdphi1),ToReal(2),cdphi211),kmul(kmadd(cdphi1,kmul(cdphi3,gtu13),kmul(cdphi2,kmadd(cdphi1,gtu12,kmul(cdphi3,gtu23)))),ToReal(4))))))))));
    
    CCTK_REAL_VEC Rphi23 = 
      kmul(ToReal(-2),kadd(cdphi223,kmadd(gt23L,kmadd(cdphi222,gtu22,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu22,SQR(cdphi2))))),ToReal(2),kmadd(gtu11,kmadd(SQR(cdphi1),ToReal(2),cdphi211),kmadd(gtu33,kmadd(SQR(cdphi3),ToReal(2),cdphi233),kmul(cdphi1,kmul(cdphi3,kmul(gtu13,ToReal(4)))))))),kmul(cdphi2,kmadd(cdphi1,kmul(gt23L,kmul(gtu12,ToReal(4))),kmul(cdphi3,kmadd(gt23L,kmul(gtu23,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi33 = 
      kmul(ToReal(-2),kadd(cdphi233,kmadd(SQR(cdphi3),kmul(kmadd(gt33L,gtu33,ToReal(-1)),ToReal(2)),kmul(gt33L,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi213,gtu13,kmul(cdphi223,gtu23)),ToReal(2),kmadd(gtu11,kmadd(SQR(cdphi1),ToReal(2),cdphi211),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmadd(cdphi3,kmul(kmadd(cdphi1,gtu13,kmul(cdphi2,gtu23)),ToReal(4)),kmul(gtu12,kmadd(cdphi212,ToReal(2),kmul(cdphi1,kmul(cdphi2,ToReal(4))))))))))))));
    
    CCTK_REAL_VEC e4phi = 
      IfThen(conformalMethod,INV(SQR(phiL)),kexp(kmul(phiL,ToReal(4))));
    
    CCTK_REAL_VEC em4phi = INV(e4phi);
    
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
    
    CCTK_REAL_VEC trR = 
      kmadd(gu11,R11,kmadd(gu22,R22,kmadd(gu33,R33,kmul(kmadd(gu12,R12,kmadd(gu13,R13,kmul(gu23,R23))),ToReal(2)))));
    
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
    
    CCTK_REAL_VEC rho = 
      kmul(INV(SQR(alphaL)),kadd(eTttL,kmadd(eTxxL,SQR(beta1L),kmadd(eTyyL,SQR(beta2L),kmadd(eTzzL,SQR(beta3L),kmadd(kmadd(beta2L,eTtyL,kmul(beta3L,eTtzL)),ToReal(-2),kmul(kmadd(beta2L,kmul(beta3L,eTyzL),kmul(beta1L,kmadd(beta2L,eTxyL,kmsub(beta3L,eTxzL,eTtxL)))),ToReal(2))))))));
    
    CCTK_REAL_VEC HL = 
      kadd(trR,kmadd(rho,ToReal(-50.26548245743669181540229413247204614715),kmadd(kmadd(Atm12,Atm21,kmadd(Atm13,Atm31,kmul(Atm23,Atm32))),ToReal(-2.),kmadd(kadd(SQR(Atm11),kadd(SQR(Atm22),SQR(Atm33))),ToReal(-1.),kmul(SQR(trKL),ToReal(0.6666666666666666666666666666666666666667))))));
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(H[index],HL,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(H[index],HL,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(H[index],HL,elt_count);
      break;
    }
    
    /* Copy local copies back to grid functions */
    vec_store_nta(H[index],HL);
  }
  LC_ENDLOOP3VEC (ML_BSSN_O2_constraints1);
}

extern "C" void ML_BSSN_O2_constraints1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_O2_constraints1_Body);
}
