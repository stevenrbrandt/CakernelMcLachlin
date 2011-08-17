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

extern "C" void ML_BSSN_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_trace_curvrhs.");
  return;
}

static void ML_BSSN_RHS1_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_RHS1_Body");
  }
  
  if (cctk_iteration % ML_BSSN_RHS1_calc_every != ML_BSSN_RHS1_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"grid::coordinates","Grid::coordinates","ML_BSSN::ML_curv","ML_BSSN::ML_dtlapse","ML_BSSN::ML_dtlapserhs","ML_BSSN::ML_dtshift","ML_BSSN::ML_dtshiftrhs","ML_BSSN::ML_Gamma","ML_BSSN::ML_Gammarhs","ML_BSSN::ML_lapse","ML_BSSN::ML_lapserhs","ML_BSSN::ML_log_confac","ML_BSSN::ML_log_confacrhs","ML_BSSN::ML_metric","ML_BSSN::ML_metricrhs","ML_BSSN::ML_shift","ML_BSSN::ML_shiftrhs","ML_BSSN::ML_trace_curv","ML_BSSN::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_RHS1", 19, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_RHS1", 2, 2, 2);
  
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
  CCTK_REAL_VEC const p1o12dx = kmul(INV(dx),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o12dy = kmul(INV(dy),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o12dz = kmul(INV(dz),ToReal(0.0833333333333333333333333333333));
  CCTK_REAL_VEC const p1o144dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o144dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o144dydz = kmul(INV(dy),kmul(INV(dz),ToReal(0.00694444444444444444444444444444)));
  CCTK_REAL_VEC const p1o24dx = kmul(INV(dx),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o24dy = kmul(INV(dy),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o24dz = kmul(INV(dz),ToReal(0.0416666666666666666666666666667));
  CCTK_REAL_VEC const p1o64dx = kmul(INV(dx),ToReal(0.015625));
  CCTK_REAL_VEC const p1o64dy = kmul(INV(dy),ToReal(0.015625));
  CCTK_REAL_VEC const p1o64dz = kmul(INV(dz),ToReal(0.015625));
  CCTK_REAL_VEC const p1odx = INV(dx);
  CCTK_REAL_VEC const p1ody = INV(dy);
  CCTK_REAL_VEC const p1odz = INV(dz);
  CCTK_REAL_VEC const pm1o12dx2 = kmul(INV(SQR(dx)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dy2 = kmul(INV(SQR(dy)),ToReal(-0.0833333333333333333333333333333));
  CCTK_REAL_VEC const pm1o12dz2 = kmul(INV(SQR(dz)),ToReal(-0.0833333333333333333333333333333));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC (ML_BSSN_RHS1,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC AL = vec_load(A[index]);
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC At11L = vec_load(At11[index]);
    CCTK_REAL_VEC At12L = vec_load(At12[index]);
    CCTK_REAL_VEC At13L = vec_load(At13[index]);
    CCTK_REAL_VEC At22L = vec_load(At22[index]);
    CCTK_REAL_VEC At23L = vec_load(At23[index]);
    CCTK_REAL_VEC At33L = vec_load(At33[index]);
    CCTK_REAL_VEC B1L = vec_load(B1[index]);
    CCTK_REAL_VEC B2L = vec_load(B2[index]);
    CCTK_REAL_VEC B3L = vec_load(B3[index]);
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
    CCTK_REAL_VEC rL = vec_load(r[index]);
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
    CCTK_REAL_VEC const PDstandardNth11beta1 = PDstandardNth11(&beta1[index]);
    CCTK_REAL_VEC const PDstandardNth22beta1 = PDstandardNth22(&beta1[index]);
    CCTK_REAL_VEC const PDstandardNth33beta1 = PDstandardNth33(&beta1[index]);
    CCTK_REAL_VEC const PDstandardNth12beta1 = PDstandardNth12(&beta1[index]);
    CCTK_REAL_VEC const PDstandardNth13beta1 = PDstandardNth13(&beta1[index]);
    CCTK_REAL_VEC const PDstandardNth23beta1 = PDstandardNth23(&beta1[index]);
    CCTK_REAL_VEC const PDstandardNth1beta2 = PDstandardNth1(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth2beta2 = PDstandardNth2(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth3beta2 = PDstandardNth3(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth11beta2 = PDstandardNth11(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth22beta2 = PDstandardNth22(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth33beta2 = PDstandardNth33(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth12beta2 = PDstandardNth12(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth13beta2 = PDstandardNth13(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth23beta2 = PDstandardNth23(&beta2[index]);
    CCTK_REAL_VEC const PDstandardNth1beta3 = PDstandardNth1(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth2beta3 = PDstandardNth2(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth3beta3 = PDstandardNth3(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth11beta3 = PDstandardNth11(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth22beta3 = PDstandardNth22(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth33beta3 = PDstandardNth33(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth12beta3 = PDstandardNth12(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth13beta3 = PDstandardNth13(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth23beta3 = PDstandardNth23(&beta3[index]);
    CCTK_REAL_VEC const PDstandardNth1gt11 = PDstandardNth1(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth2gt11 = PDstandardNth2(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth3gt11 = PDstandardNth3(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth1gt12 = PDstandardNth1(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth2gt12 = PDstandardNth2(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth3gt12 = PDstandardNth3(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth1gt13 = PDstandardNth1(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth2gt13 = PDstandardNth2(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth3gt13 = PDstandardNth3(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth1gt22 = PDstandardNth1(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth2gt22 = PDstandardNth2(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth3gt22 = PDstandardNth3(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth1gt23 = PDstandardNth1(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth2gt23 = PDstandardNth2(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth3gt23 = PDstandardNth3(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth1gt33 = PDstandardNth1(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth2gt33 = PDstandardNth2(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth3gt33 = PDstandardNth3(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth1phi = PDstandardNth1(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth2phi = PDstandardNth2(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth3phi = PDstandardNth3(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth1trK = PDstandardNth1(&trK[index]);
    CCTK_REAL_VEC const PDstandardNth2trK = PDstandardNth2(&trK[index]);
    CCTK_REAL_VEC const PDstandardNth3trK = PDstandardNth3(&trK[index]);
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
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
    
    CCTK_REAL_VEC fac1 = 
      IfThen(conformalMethod,kmul(INV(phiL),ToReal(-0.5)),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 = kmul(fac1,PDstandardNth1phi);
    
    CCTK_REAL_VEC cdphi2 = kmul(fac1,PDstandardNth2phi);
    
    CCTK_REAL_VEC cdphi3 = kmul(fac1,PDstandardNth3phi);
    
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
    
    CCTK_REAL_VEC Atu11 = 
      kmadd(Atm11,gtu11,kmadd(Atm12,gtu12,kmul(Atm13,gtu13)));
    
    CCTK_REAL_VEC Atu12 = 
      kmadd(Atm11,gtu12,kmadd(Atm12,gtu22,kmul(Atm13,gtu23)));
    
    CCTK_REAL_VEC Atu13 = 
      kmadd(Atm11,gtu13,kmadd(Atm12,gtu23,kmul(Atm13,gtu33)));
    
    CCTK_REAL_VEC Atu22 = 
      kmadd(Atm21,gtu12,kmadd(Atm22,gtu22,kmul(Atm23,gtu23)));
    
    CCTK_REAL_VEC Atu23 = 
      kmadd(Atm21,gtu13,kmadd(Atm22,gtu23,kmul(Atm23,gtu33)));
    
    CCTK_REAL_VEC Atu33 = 
      kmadd(Atm31,gtu13,kmadd(Atm32,gtu23,kmul(Atm33,gtu33)));
    
    CCTK_REAL_VEC e4phi = 
      IfThen(conformalMethod,INV(SQR(phiL)),kexp(kmul(phiL,ToReal(4))));
    
    CCTK_REAL_VEC em4phi = INV(e4phi);
    
    CCTK_REAL_VEC rho = 
      kmul(INV(SQR(alphaL)),kadd(eTttL,kmadd(eTxxL,SQR(beta1L),kmadd(eTyyL,SQR(beta2L),kmadd(eTzzL,SQR(beta3L),kmadd(kmadd(beta2L,eTtyL,kmul(beta3L,eTtzL)),ToReal(-2),kmul(kmadd(beta2L,kmul(beta3L,eTyzL),kmul(beta1L,kmadd(beta2L,eTxyL,kmsub(beta3L,eTxzL,eTtxL)))),ToReal(2))))))));
    
    CCTK_REAL_VEC S1 = 
      kmul(INV(alphaL),kmadd(beta1L,eTxxL,kmadd(beta2L,eTxyL,kmsub(beta3L,eTxzL,eTtxL))));
    
    CCTK_REAL_VEC S2 = 
      kmul(INV(alphaL),kmadd(beta1L,eTxyL,kmadd(beta2L,eTyyL,kmsub(beta3L,eTyzL,eTtyL))));
    
    CCTK_REAL_VEC S3 = 
      kmul(INV(alphaL),kmadd(beta1L,eTxzL,kmadd(beta2L,eTyzL,kmsub(beta3L,eTzzL,eTtzL))));
    
    CCTK_REAL_VEC trS = 
      kmul(em4phi,kmadd(eTxxL,gtu11,kmadd(eTyyL,gtu22,kmadd(eTzzL,gtu33,kmul(kmadd(eTxyL,gtu12,kmadd(eTxzL,gtu13,kmul(eTyzL,gtu23))),ToReal(2))))));
    
    CCTK_REAL_VEC phirhsL = 
      IfThen(conformalMethod,kmul(phiL,kmadd(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(-0.333333333333333333333333333333),kmul(alphaL,kmul(trKL,ToReal(0.333333333333333333333333333333))))),kmadd(alphaL,kmul(trKL,ToReal(-0.166666666666666666666666666667)),kmul(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(0.166666666666666666666666666667))));
    
    CCTK_REAL_VEC gt11rhsL = 
      kmul(ToReal(-0.666666666666666666666666666667),kmadd(kmadd(gt12L,PDstandardNth1beta2,kmul(gt13L,PDstandardNth1beta3)),ToReal(-3),kmadd(gt11L,kadd(PDstandardNth2beta2,kmadd(PDstandardNth1beta1,ToReal(-2),PDstandardNth3beta3)),kmul(alphaL,kmul(At11L,ToReal(3))))));
    
    CCTK_REAL_VEC gt12rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At12L,ToReal(-6)),kmadd(gt12L,kadd(PDstandardNth1beta1,kmadd(PDstandardNth3beta3,ToReal(-2),PDstandardNth2beta2)),kmul(kmadd(gt22L,PDstandardNth1beta2,kmadd(gt23L,PDstandardNth1beta3,kmadd(gt11L,PDstandardNth2beta1,kmul(gt13L,PDstandardNth2beta3)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt13rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At13L,ToReal(-6)),kmadd(gt13L,kadd(PDstandardNth1beta1,kmadd(PDstandardNth2beta2,ToReal(-2),PDstandardNth3beta3)),kmul(kmadd(gt23L,PDstandardNth1beta2,kmadd(gt33L,PDstandardNth1beta3,kmadd(gt11L,PDstandardNth3beta1,kmul(gt12L,PDstandardNth3beta2)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt22rhsL = 
      kmul(ToReal(-0.666666666666666666666666666667),kmadd(kmadd(gt12L,PDstandardNth2beta1,kmul(gt23L,PDstandardNth2beta3)),ToReal(-3),kmadd(gt22L,kadd(PDstandardNth1beta1,kmadd(PDstandardNth2beta2,ToReal(-2),PDstandardNth3beta3)),kmul(alphaL,kmul(At22L,ToReal(3))))));
    
    CCTK_REAL_VEC gt23rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At23L,ToReal(-6)),kmadd(gt23L,kadd(PDstandardNth2beta2,kmadd(PDstandardNth1beta1,ToReal(-2),PDstandardNth3beta3)),kmul(kmadd(gt13L,PDstandardNth2beta1,kmadd(gt33L,PDstandardNth2beta3,kmadd(gt12L,PDstandardNth3beta1,kmul(gt22L,PDstandardNth3beta2)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt33rhsL = 
      kmul(ToReal(-0.666666666666666666666666666667),kmadd(kmadd(gt13L,PDstandardNth3beta1,kmul(gt23L,PDstandardNth3beta2)),ToReal(-3),kmadd(gt33L,kadd(PDstandardNth1beta1,kmadd(PDstandardNth3beta3,ToReal(-2),PDstandardNth2beta2)),kmul(alphaL,kmul(At33L,ToReal(3))))));
    
    CCTK_REAL_VEC dotXt1 = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(gtu12,kadd(PDstandardNth22beta2,kmadd(alphaL,kmul(S2,ToReal(-150.7964473723100754462068823974161384415)),PDstandardNth23beta3)),kmadd(gtu13,kadd(PDstandardNth23beta2,kmadd(alphaL,kmul(S3,ToReal(-150.7964473723100754462068823974161384415)),PDstandardNth33beta3)),kmadd(kmadd(Atu11,PDstandardNth1alpha,kmadd(Atu12,PDstandardNth2alpha,kmul(Atu13,PDstandardNth3alpha))),ToReal(-6),kmadd(kmadd(PDstandardNth2beta1,Xtn2,kmul(PDstandardNth3beta1,Xtn3)),ToReal(-3),kmadd(Xtn1,kmsub(PDstandardNth3beta3,ToReal(2),PDstandardNth1beta1),kmadd(kmadd(gtu22,PDstandardNth22beta1,kmul(gtu33,PDstandardNth33beta1)),ToReal(3),kmadd(gtu11,kadd(PDstandardNth12beta2,kadd(PDstandardNth13beta3,kmadd(alphaL,kmul(S1,ToReal(-150.7964473723100754462068823974161384415)),kmul(PDstandardNth11beta1,ToReal(4))))),kmadd(gtu23,kmul(PDstandardNth23beta1,ToReal(6)),kmadd(kmadd(gtu12,PDstandardNth12beta1,kmul(gtu13,PDstandardNth13beta1)),ToReal(7),kmul(ToReal(2),kmadd(PDstandardNth2beta2,Xtn1,kmul(alphaL,kmadd(kmadd(gtu11,PDstandardNth1trK,kmadd(gtu12,PDstandardNth2trK,kmul(gtu13,PDstandardNth3trK))),ToReal(-2),kmadd(kmadd(Atu11,Gt111,kmadd(Atu22,Gt122,kmul(Atu33,Gt133))),ToReal(3),kmadd(kmadd(Atu12,Gt112,kmadd(Atu13,Gt113,kmul(Atu23,Gt123))),ToReal(6),kmul(kmadd(Atu11,cdphi1,kmadd(Atu12,cdphi2,kmul(Atu13,cdphi3))),ToReal(18))))))))))))))))));
    
    CCTK_REAL_VEC dotXt2 = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(kmadd(Atu12,PDstandardNth1alpha,kmadd(Atu22,PDstandardNth2alpha,kmul(Atu23,PDstandardNth3alpha))),ToReal(-6),kmadd(kmadd(PDstandardNth1beta2,Xtn1,kmul(PDstandardNth3beta2,Xtn3)),ToReal(-3),kmadd(Xtn2,kmsub(PDstandardNth3beta3,ToReal(2),PDstandardNth2beta2),kmadd(kmadd(gtu11,PDstandardNth11beta2,kmul(gtu33,PDstandardNth33beta2)),ToReal(3),kmadd(gtu22,kadd(PDstandardNth12beta1,kadd(PDstandardNth23beta3,kmadd(alphaL,kmul(S2,ToReal(-150.7964473723100754462068823974161384415)),kmul(PDstandardNth22beta2,ToReal(4))))),kmadd(gtu13,kmul(PDstandardNth13beta2,ToReal(6)),kmadd(gtu12,kadd(PDstandardNth11beta1,kadd(PDstandardNth13beta3,kmadd(alphaL,kmul(S1,ToReal(-150.7964473723100754462068823974161384415)),kmul(PDstandardNth12beta2,ToReal(7))))),kmadd(gtu23,kadd(PDstandardNth13beta1,kadd(PDstandardNth33beta3,kmadd(alphaL,kmul(S3,ToReal(-150.7964473723100754462068823974161384415)),kmul(PDstandardNth23beta2,ToReal(7))))),kmul(ToReal(2),kmadd(PDstandardNth1beta1,Xtn2,kmul(alphaL,kmadd(kmadd(gtu12,PDstandardNth1trK,kmadd(gtu22,PDstandardNth2trK,kmul(gtu23,PDstandardNth3trK))),ToReal(-2),kmadd(kmadd(Atu11,Gt211,kmadd(Atu22,Gt222,kmul(Atu33,Gt233))),ToReal(3),kmadd(kmadd(Atu12,Gt212,kmadd(Atu13,Gt213,kmul(Atu23,Gt223))),ToReal(6),kmul(kmadd(Atu12,cdphi1,kmadd(Atu22,cdphi2,kmul(Atu23,cdphi3))),ToReal(18)))))))))))))))));
    
    CCTK_REAL_VEC dotXt3 = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(kmadd(Atu13,PDstandardNth1alpha,kmadd(Atu23,PDstandardNth2alpha,kmul(Atu33,PDstandardNth3alpha))),ToReal(-6),kmadd(kmadd(PDstandardNth1beta3,Xtn1,kmul(PDstandardNth2beta3,Xtn2)),ToReal(-3),kmadd(Xtn3,kmsub(PDstandardNth2beta2,ToReal(2),PDstandardNth3beta3),kmadd(kmadd(gtu11,PDstandardNth11beta3,kmul(gtu22,PDstandardNth22beta3)),ToReal(3),kmadd(gtu33,kadd(PDstandardNth13beta1,kadd(PDstandardNth23beta2,kmadd(alphaL,kmul(S3,ToReal(-150.7964473723100754462068823974161384415)),kmul(PDstandardNth33beta3,ToReal(4))))),kmadd(gtu12,kmul(PDstandardNth12beta3,ToReal(6)),kmadd(gtu13,kadd(PDstandardNth11beta1,kadd(PDstandardNth12beta2,kmadd(alphaL,kmul(S1,ToReal(-150.7964473723100754462068823974161384415)),kmul(PDstandardNth13beta3,ToReal(7))))),kmadd(gtu23,kadd(PDstandardNth12beta1,kadd(PDstandardNth22beta2,kmadd(alphaL,kmul(S2,ToReal(-150.7964473723100754462068823974161384415)),kmul(PDstandardNth23beta3,ToReal(7))))),kmul(ToReal(2),kmadd(PDstandardNth1beta1,Xtn3,kmul(alphaL,kmadd(kmadd(gtu13,PDstandardNth1trK,kmadd(gtu23,PDstandardNth2trK,kmul(gtu33,PDstandardNth3trK))),ToReal(-2),kmadd(kmadd(Atu11,Gt311,kmadd(Atu22,Gt322,kmul(Atu33,Gt333))),ToReal(3),kmadd(kmadd(Atu12,Gt312,kmadd(Atu13,Gt313,kmul(Atu23,Gt323))),ToReal(6),kmul(kmadd(Atu13,cdphi1,kmadd(Atu23,cdphi2,kmul(Atu33,cdphi3))),ToReal(18)))))))))))))))));
    
    CCTK_REAL_VEC Xt1rhsL = dotXt1;
    
    CCTK_REAL_VEC Xt2rhsL = dotXt2;
    
    CCTK_REAL_VEC Xt3rhsL = dotXt3;
    
    CCTK_REAL_VEC dottrK = 
      kmsub(alphaL,kadd(SQR(Atm11),kadd(SQR(Atm22),kadd(SQR(Atm33),kmadd(SQR(trKL),ToReal(0.333333333333333333333333333333),kmadd(kmadd(Atm12,Atm21,kmadd(Atm13,Atm31,kmul(Atm23,Atm32))),ToReal(2),kmul(kadd(rho,trS),ToReal(12.56637061435917295385057353311801153679))))))),kmul(em4phi,kmadd(gtu11,PDstandardNth11alpha,kmadd(gtu22,PDstandardNth22alpha,knmsub(PDstandardNth3alpha,Xtn3,kmadd(kmadd(gtu12,PDstandardNth12alpha,kmadd(gtu13,kmadd(cdphi1,PDstandardNth3alpha,PDstandardNth13alpha),kmul(gtu23,kmadd(cdphi2,PDstandardNth3alpha,PDstandardNth23alpha)))),ToReal(2),kmadd(PDstandardNth1alpha,kmsub(kmadd(cdphi1,gtu11,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13))),ToReal(2),Xtn1),kmadd(PDstandardNth2alpha,kmsub(kmadd(cdphi1,gtu12,kmadd(cdphi2,gtu22,kmul(cdphi3,gtu23))),ToReal(2),Xtn2),kmul(gtu33,kmadd(cdphi3,kmul(PDstandardNth3alpha,ToReal(2)),PDstandardNth33alpha))))))))));
    
    CCTK_REAL_VEC trKrhsL = dottrK;
    
    CCTK_REAL_VEC alpharhsL = 
      kneg(kmul(kpow(alphaL,harmonicN),kmul(ToReal(harmonicF),kmadd(ksub(AL,trKL),ToReal(LapseACoeff),trKL))));
    
    CCTK_REAL_VEC ArhsL = 
      kmul(knmsub(AL,ToReal(AlphaDriver),dottrK),ToReal(LapseACoeff));
    
    CCTK_REAL_VEC eta = 
      kfmin(ToReal(1),kmul(INV(rL),ToReal(SpatialBetaDriverRadius)));
    
    CCTK_REAL_VEC theta = 
      kfmin(ToReal(1),kexp(knmsub(rL,INV(ToReal(SpatialShiftGammaCoeffRadius)),ToReal(1))));
    
    CCTK_REAL_VEC beta1rhsL = 
      kmul(theta,kmul(kadd(Xt1L,kmadd(beta1L,kmul(eta,ToReal(BetaDriver*(-1 + 
      ShiftBCoeff))),kmul(ksub(B1L,Xt1L),ToReal(ShiftBCoeff)))),ToReal(ShiftGammaCoeff)));
    
    CCTK_REAL_VEC beta2rhsL = 
      kmul(theta,kmul(kadd(Xt2L,kmadd(beta2L,kmul(eta,ToReal(BetaDriver*(-1 + 
      ShiftBCoeff))),kmul(ksub(B2L,Xt2L),ToReal(ShiftBCoeff)))),ToReal(ShiftGammaCoeff)));
    
    CCTK_REAL_VEC beta3rhsL = 
      kmul(theta,kmul(kadd(Xt3L,kmadd(beta3L,kmul(eta,ToReal(BetaDriver*(-1 + 
      ShiftBCoeff))),kmul(ksub(B3L,Xt3L),ToReal(ShiftBCoeff)))),ToReal(ShiftGammaCoeff)));
    
    CCTK_REAL_VEC B1rhsL = 
      kmul(knmsub(B1L,kmul(eta,ToReal(BetaDriver)),dotXt1),ToReal(ShiftBCoeff));
    
    CCTK_REAL_VEC B2rhsL = 
      kmul(knmsub(B2L,kmul(eta,ToReal(BetaDriver)),dotXt2),ToReal(ShiftBCoeff));
    
    CCTK_REAL_VEC B3rhsL = 
      kmul(knmsub(B3L,kmul(eta,ToReal(BetaDriver)),dotXt3),ToReal(ShiftBCoeff));
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(alpharhs[index],alpharhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Arhs[index],ArhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B1rhs[index],B1rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B2rhs[index],B2rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(B3rhs[index],B3rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta1rhs[index],beta1rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta2rhs[index],beta2rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(beta3rhs[index],beta3rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt11rhs[index],gt11rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt12rhs[index],gt12rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt13rhs[index],gt13rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt22rhs[index],gt22rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt23rhs[index],gt23rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(gt33rhs[index],gt33rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(phirhs[index],phirhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(trKrhs[index],trKrhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt1rhs[index],Xt1rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt2rhs[index],Xt2rhsL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(Xt3rhs[index],Xt3rhsL,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(alpharhs[index],alpharhsL,elt_count);
      vec_store_nta_partial_hi(Arhs[index],ArhsL,elt_count);
      vec_store_nta_partial_hi(B1rhs[index],B1rhsL,elt_count);
      vec_store_nta_partial_hi(B2rhs[index],B2rhsL,elt_count);
      vec_store_nta_partial_hi(B3rhs[index],B3rhsL,elt_count);
      vec_store_nta_partial_hi(beta1rhs[index],beta1rhsL,elt_count);
      vec_store_nta_partial_hi(beta2rhs[index],beta2rhsL,elt_count);
      vec_store_nta_partial_hi(beta3rhs[index],beta3rhsL,elt_count);
      vec_store_nta_partial_hi(gt11rhs[index],gt11rhsL,elt_count);
      vec_store_nta_partial_hi(gt12rhs[index],gt12rhsL,elt_count);
      vec_store_nta_partial_hi(gt13rhs[index],gt13rhsL,elt_count);
      vec_store_nta_partial_hi(gt22rhs[index],gt22rhsL,elt_count);
      vec_store_nta_partial_hi(gt23rhs[index],gt23rhsL,elt_count);
      vec_store_nta_partial_hi(gt33rhs[index],gt33rhsL,elt_count);
      vec_store_nta_partial_hi(phirhs[index],phirhsL,elt_count);
      vec_store_nta_partial_hi(trKrhs[index],trKrhsL,elt_count);
      vec_store_nta_partial_hi(Xt1rhs[index],Xt1rhsL,elt_count);
      vec_store_nta_partial_hi(Xt2rhs[index],Xt2rhsL,elt_count);
      vec_store_nta_partial_hi(Xt3rhs[index],Xt3rhsL,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(alpharhs[index],alpharhsL,elt_count);
      vec_store_nta_partial_lo(Arhs[index],ArhsL,elt_count);
      vec_store_nta_partial_lo(B1rhs[index],B1rhsL,elt_count);
      vec_store_nta_partial_lo(B2rhs[index],B2rhsL,elt_count);
      vec_store_nta_partial_lo(B3rhs[index],B3rhsL,elt_count);
      vec_store_nta_partial_lo(beta1rhs[index],beta1rhsL,elt_count);
      vec_store_nta_partial_lo(beta2rhs[index],beta2rhsL,elt_count);
      vec_store_nta_partial_lo(beta3rhs[index],beta3rhsL,elt_count);
      vec_store_nta_partial_lo(gt11rhs[index],gt11rhsL,elt_count);
      vec_store_nta_partial_lo(gt12rhs[index],gt12rhsL,elt_count);
      vec_store_nta_partial_lo(gt13rhs[index],gt13rhsL,elt_count);
      vec_store_nta_partial_lo(gt22rhs[index],gt22rhsL,elt_count);
      vec_store_nta_partial_lo(gt23rhs[index],gt23rhsL,elt_count);
      vec_store_nta_partial_lo(gt33rhs[index],gt33rhsL,elt_count);
      vec_store_nta_partial_lo(phirhs[index],phirhsL,elt_count);
      vec_store_nta_partial_lo(trKrhs[index],trKrhsL,elt_count);
      vec_store_nta_partial_lo(Xt1rhs[index],Xt1rhsL,elt_count);
      vec_store_nta_partial_lo(Xt2rhs[index],Xt2rhsL,elt_count);
      vec_store_nta_partial_lo(Xt3rhs[index],Xt3rhsL,elt_count);
      break;
    }
    
    /* Copy local copies back to grid functions */
    vec_store_nta(alpharhs[index],alpharhsL);
    vec_store_nta(Arhs[index],ArhsL);
    vec_store_nta(B1rhs[index],B1rhsL);
    vec_store_nta(B2rhs[index],B2rhsL);
    vec_store_nta(B3rhs[index],B3rhsL);
    vec_store_nta(beta1rhs[index],beta1rhsL);
    vec_store_nta(beta2rhs[index],beta2rhsL);
    vec_store_nta(beta3rhs[index],beta3rhsL);
    vec_store_nta(gt11rhs[index],gt11rhsL);
    vec_store_nta(gt12rhs[index],gt12rhsL);
    vec_store_nta(gt13rhs[index],gt13rhsL);
    vec_store_nta(gt22rhs[index],gt22rhsL);
    vec_store_nta(gt23rhs[index],gt23rhsL);
    vec_store_nta(gt33rhs[index],gt33rhsL);
    vec_store_nta(phirhs[index],phirhsL);
    vec_store_nta(trKrhs[index],trKrhsL);
    vec_store_nta(Xt1rhs[index],Xt1rhsL);
    vec_store_nta(Xt2rhs[index],Xt2rhsL);
    vec_store_nta(Xt3rhs[index],Xt3rhsL);
  }
  LC_ENDLOOP3VEC (ML_BSSN_RHS1);
}

extern "C" void ML_BSSN_RHS1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_RHS1_Body);
}
