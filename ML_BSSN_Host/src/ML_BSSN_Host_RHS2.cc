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
#define QAD(x) (SQR(SQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))

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
  CCTK_REAL_VEC const p1o1024dx = kmul(INV(dx),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1024dy = kmul(INV(dy),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1024dz = kmul(INV(dz),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1680dx = kmul(INV(dx),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o1680dy = kmul(INV(dy),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o1680dz = kmul(INV(dz),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o5040dx2 = kmul(INV(SQR(dx)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dy2 = kmul(INV(SQR(dy)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dz2 = kmul(INV(SQR(dz)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o560dx = kmul(INV(dx),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o560dy = kmul(INV(dy),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o560dz = kmul(INV(dz),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o705600dxdy = kmul(INV(kmul(dx,dy)),ToReal(1.41723356009070294784580498866e-6));
  CCTK_REAL_VEC const p1o705600dxdz = kmul(INV(kmul(dx,dz)),ToReal(1.41723356009070294784580498866e-6));
  CCTK_REAL_VEC const p1o705600dydz = kmul(INV(kmul(dy,dz)),ToReal(1.41723356009070294784580498866e-6));
  CCTK_REAL_VEC const p1o840dx = kmul(INV(dx),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dy = kmul(INV(dy),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dz = kmul(INV(dz),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1odx = INV(dx);
  CCTK_REAL_VEC const p1ody = INV(dy);
  CCTK_REAL_VEC const p1odz = INV(dz);
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC(ML_BSSN_Host_RHS2,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
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
    CCTK_REAL_VEC DPDstandardNthalpha1L = vec_load(DPDstandardNthalpha1[index]);
    CCTK_REAL_VEC DPDstandardNthalpha11L = vec_load(DPDstandardNthalpha11[index]);
    CCTK_REAL_VEC DPDstandardNthalpha12L = vec_load(DPDstandardNthalpha12[index]);
    CCTK_REAL_VEC DPDstandardNthalpha13L = vec_load(DPDstandardNthalpha13[index]);
    CCTK_REAL_VEC DPDstandardNthalpha2L = vec_load(DPDstandardNthalpha2[index]);
    CCTK_REAL_VEC DPDstandardNthalpha22L = vec_load(DPDstandardNthalpha22[index]);
    CCTK_REAL_VEC DPDstandardNthalpha23L = vec_load(DPDstandardNthalpha23[index]);
    CCTK_REAL_VEC DPDstandardNthalpha3L = vec_load(DPDstandardNthalpha3[index]);
    CCTK_REAL_VEC DPDstandardNthalpha33L = vec_load(DPDstandardNthalpha33[index]);
    CCTK_REAL_VEC DPDstandardNthbeta11L = vec_load(DPDstandardNthbeta11[index]);
    CCTK_REAL_VEC DPDstandardNthbeta12L = vec_load(DPDstandardNthbeta12[index]);
    CCTK_REAL_VEC DPDstandardNthbeta13L = vec_load(DPDstandardNthbeta13[index]);
    CCTK_REAL_VEC DPDstandardNthbeta21L = vec_load(DPDstandardNthbeta21[index]);
    CCTK_REAL_VEC DPDstandardNthbeta22L = vec_load(DPDstandardNthbeta22[index]);
    CCTK_REAL_VEC DPDstandardNthbeta23L = vec_load(DPDstandardNthbeta23[index]);
    CCTK_REAL_VEC DPDstandardNthbeta31L = vec_load(DPDstandardNthbeta31[index]);
    CCTK_REAL_VEC DPDstandardNthbeta32L = vec_load(DPDstandardNthbeta32[index]);
    CCTK_REAL_VEC DPDstandardNthbeta33L = vec_load(DPDstandardNthbeta33[index]);
    CCTK_REAL_VEC DPDstandardNthgt111L = vec_load(DPDstandardNthgt111[index]);
    CCTK_REAL_VEC DPDstandardNthgt1111L = vec_load(DPDstandardNthgt1111[index]);
    CCTK_REAL_VEC DPDstandardNthgt1112L = vec_load(DPDstandardNthgt1112[index]);
    CCTK_REAL_VEC DPDstandardNthgt1113L = vec_load(DPDstandardNthgt1113[index]);
    CCTK_REAL_VEC DPDstandardNthgt112L = vec_load(DPDstandardNthgt112[index]);
    CCTK_REAL_VEC DPDstandardNthgt1122L = vec_load(DPDstandardNthgt1122[index]);
    CCTK_REAL_VEC DPDstandardNthgt1123L = vec_load(DPDstandardNthgt1123[index]);
    CCTK_REAL_VEC DPDstandardNthgt113L = vec_load(DPDstandardNthgt113[index]);
    CCTK_REAL_VEC DPDstandardNthgt1133L = vec_load(DPDstandardNthgt1133[index]);
    CCTK_REAL_VEC DPDstandardNthgt121L = vec_load(DPDstandardNthgt121[index]);
    CCTK_REAL_VEC DPDstandardNthgt1211L = vec_load(DPDstandardNthgt1211[index]);
    CCTK_REAL_VEC DPDstandardNthgt1212L = vec_load(DPDstandardNthgt1212[index]);
    CCTK_REAL_VEC DPDstandardNthgt1213L = vec_load(DPDstandardNthgt1213[index]);
    CCTK_REAL_VEC DPDstandardNthgt122L = vec_load(DPDstandardNthgt122[index]);
    CCTK_REAL_VEC DPDstandardNthgt1222L = vec_load(DPDstandardNthgt1222[index]);
    CCTK_REAL_VEC DPDstandardNthgt1223L = vec_load(DPDstandardNthgt1223[index]);
    CCTK_REAL_VEC DPDstandardNthgt123L = vec_load(DPDstandardNthgt123[index]);
    CCTK_REAL_VEC DPDstandardNthgt1233L = vec_load(DPDstandardNthgt1233[index]);
    CCTK_REAL_VEC DPDstandardNthgt131L = vec_load(DPDstandardNthgt131[index]);
    CCTK_REAL_VEC DPDstandardNthgt1311L = vec_load(DPDstandardNthgt1311[index]);
    CCTK_REAL_VEC DPDstandardNthgt1312L = vec_load(DPDstandardNthgt1312[index]);
    CCTK_REAL_VEC DPDstandardNthgt1313L = vec_load(DPDstandardNthgt1313[index]);
    CCTK_REAL_VEC DPDstandardNthgt132L = vec_load(DPDstandardNthgt132[index]);
    CCTK_REAL_VEC DPDstandardNthgt1322L = vec_load(DPDstandardNthgt1322[index]);
    CCTK_REAL_VEC DPDstandardNthgt1323L = vec_load(DPDstandardNthgt1323[index]);
    CCTK_REAL_VEC DPDstandardNthgt133L = vec_load(DPDstandardNthgt133[index]);
    CCTK_REAL_VEC DPDstandardNthgt1333L = vec_load(DPDstandardNthgt1333[index]);
    CCTK_REAL_VEC DPDstandardNthgt221L = vec_load(DPDstandardNthgt221[index]);
    CCTK_REAL_VEC DPDstandardNthgt2211L = vec_load(DPDstandardNthgt2211[index]);
    CCTK_REAL_VEC DPDstandardNthgt2212L = vec_load(DPDstandardNthgt2212[index]);
    CCTK_REAL_VEC DPDstandardNthgt2213L = vec_load(DPDstandardNthgt2213[index]);
    CCTK_REAL_VEC DPDstandardNthgt222L = vec_load(DPDstandardNthgt222[index]);
    CCTK_REAL_VEC DPDstandardNthgt2222L = vec_load(DPDstandardNthgt2222[index]);
    CCTK_REAL_VEC DPDstandardNthgt2223L = vec_load(DPDstandardNthgt2223[index]);
    CCTK_REAL_VEC DPDstandardNthgt223L = vec_load(DPDstandardNthgt223[index]);
    CCTK_REAL_VEC DPDstandardNthgt2233L = vec_load(DPDstandardNthgt2233[index]);
    CCTK_REAL_VEC DPDstandardNthgt231L = vec_load(DPDstandardNthgt231[index]);
    CCTK_REAL_VEC DPDstandardNthgt2311L = vec_load(DPDstandardNthgt2311[index]);
    CCTK_REAL_VEC DPDstandardNthgt2312L = vec_load(DPDstandardNthgt2312[index]);
    CCTK_REAL_VEC DPDstandardNthgt2313L = vec_load(DPDstandardNthgt2313[index]);
    CCTK_REAL_VEC DPDstandardNthgt232L = vec_load(DPDstandardNthgt232[index]);
    CCTK_REAL_VEC DPDstandardNthgt2322L = vec_load(DPDstandardNthgt2322[index]);
    CCTK_REAL_VEC DPDstandardNthgt2323L = vec_load(DPDstandardNthgt2323[index]);
    CCTK_REAL_VEC DPDstandardNthgt233L = vec_load(DPDstandardNthgt233[index]);
    CCTK_REAL_VEC DPDstandardNthgt2333L = vec_load(DPDstandardNthgt2333[index]);
    CCTK_REAL_VEC DPDstandardNthgt331L = vec_load(DPDstandardNthgt331[index]);
    CCTK_REAL_VEC DPDstandardNthgt3311L = vec_load(DPDstandardNthgt3311[index]);
    CCTK_REAL_VEC DPDstandardNthgt3312L = vec_load(DPDstandardNthgt3312[index]);
    CCTK_REAL_VEC DPDstandardNthgt3313L = vec_load(DPDstandardNthgt3313[index]);
    CCTK_REAL_VEC DPDstandardNthgt332L = vec_load(DPDstandardNthgt332[index]);
    CCTK_REAL_VEC DPDstandardNthgt3322L = vec_load(DPDstandardNthgt3322[index]);
    CCTK_REAL_VEC DPDstandardNthgt3323L = vec_load(DPDstandardNthgt3323[index]);
    CCTK_REAL_VEC DPDstandardNthgt333L = vec_load(DPDstandardNthgt333[index]);
    CCTK_REAL_VEC DPDstandardNthgt3333L = vec_load(DPDstandardNthgt3333[index]);
    CCTK_REAL_VEC DPDstandardNthphi1L = vec_load(DPDstandardNthphi1[index]);
    CCTK_REAL_VEC DPDstandardNthphi11L = vec_load(DPDstandardNthphi11[index]);
    CCTK_REAL_VEC DPDstandardNthphi12L = vec_load(DPDstandardNthphi12[index]);
    CCTK_REAL_VEC DPDstandardNthphi13L = vec_load(DPDstandardNthphi13[index]);
    CCTK_REAL_VEC DPDstandardNthphi2L = vec_load(DPDstandardNthphi2[index]);
    CCTK_REAL_VEC DPDstandardNthphi22L = vec_load(DPDstandardNthphi22[index]);
    CCTK_REAL_VEC DPDstandardNthphi23L = vec_load(DPDstandardNthphi23[index]);
    CCTK_REAL_VEC DPDstandardNthphi3L = vec_load(DPDstandardNthphi3[index]);
    CCTK_REAL_VEC DPDstandardNthphi33L = vec_load(DPDstandardNthphi33[index]);
    CCTK_REAL_VEC DPDstandardNthXt11L = vec_load(DPDstandardNthXt11[index]);
    CCTK_REAL_VEC DPDstandardNthXt12L = vec_load(DPDstandardNthXt12[index]);
    CCTK_REAL_VEC DPDstandardNthXt13L = vec_load(DPDstandardNthXt13[index]);
    CCTK_REAL_VEC DPDstandardNthXt21L = vec_load(DPDstandardNthXt21[index]);
    CCTK_REAL_VEC DPDstandardNthXt22L = vec_load(DPDstandardNthXt22[index]);
    CCTK_REAL_VEC DPDstandardNthXt23L = vec_load(DPDstandardNthXt23[index]);
    CCTK_REAL_VEC DPDstandardNthXt31L = vec_load(DPDstandardNthXt31[index]);
    CCTK_REAL_VEC DPDstandardNthXt32L = vec_load(DPDstandardNthXt32[index]);
    CCTK_REAL_VEC DPDstandardNthXt33L = vec_load(DPDstandardNthXt33[index]);
    CCTK_REAL_VEC gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC phiL = vec_load(phi[index]);
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    CCTK_REAL_VEC detgt = ToReal(1);
    
    CCTK_REAL_VEC gtu11 = 
      kmul(INV(detgt),kmsub(gt22L,gt33L,SQR(gt23L)));
    
    CCTK_REAL_VEC gtu12 = 
      kmul(INV(detgt),kmsub(gt13L,gt23L,kmul(gt12L,gt33L)));
    
    CCTK_REAL_VEC gtu13 = 
      kmul(INV(detgt),kmsub(gt12L,gt23L,kmul(gt13L,gt22L)));
    
    CCTK_REAL_VEC gtu22 = 
      kmul(INV(detgt),kmsub(gt11L,gt33L,SQR(gt13L)));
    
    CCTK_REAL_VEC gtu23 = 
      kmul(INV(detgt),kmsub(gt12L,gt13L,kmul(gt11L,gt23L)));
    
    CCTK_REAL_VEC gtu33 = 
      kmul(INV(detgt),kmsub(gt11L,gt22L,SQR(gt12L)));
    
    CCTK_REAL_VEC Gtl111 = kmul(DPDstandardNthgt111L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl112 = kmul(DPDstandardNthgt112L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl113 = kmul(DPDstandardNthgt113L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl122 = 
      kmadd(DPDstandardNthgt221L,ToReal(-0.5),DPDstandardNthgt122L);
    
    CCTK_REAL_VEC Gtl123 = 
      kmul(kadd(DPDstandardNthgt123L,ksub(DPDstandardNthgt132L,DPDstandardNthgt231L)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl133 = 
      kmadd(DPDstandardNthgt331L,ToReal(-0.5),DPDstandardNthgt133L);
    
    CCTK_REAL_VEC Gtl211 = 
      kmadd(DPDstandardNthgt112L,ToReal(-0.5),DPDstandardNthgt121L);
    
    CCTK_REAL_VEC Gtl212 = kmul(DPDstandardNthgt221L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl213 = 
      kmul(kadd(DPDstandardNthgt123L,ksub(DPDstandardNthgt231L,DPDstandardNthgt132L)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl222 = kmul(DPDstandardNthgt222L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl223 = kmul(DPDstandardNthgt223L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl233 = 
      kmadd(DPDstandardNthgt332L,ToReal(-0.5),DPDstandardNthgt233L);
    
    CCTK_REAL_VEC Gtl311 = 
      kmadd(DPDstandardNthgt113L,ToReal(-0.5),DPDstandardNthgt131L);
    
    CCTK_REAL_VEC Gtl312 = 
      kmul(kadd(DPDstandardNthgt132L,ksub(DPDstandardNthgt231L,DPDstandardNthgt123L)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl313 = kmul(DPDstandardNthgt331L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl322 = 
      kmadd(DPDstandardNthgt223L,ToReal(-0.5),DPDstandardNthgt232L);
    
    CCTK_REAL_VEC Gtl323 = kmul(DPDstandardNthgt332L,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl333 = kmul(DPDstandardNthgt333L,ToReal(0.5));
    
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
      kmadd(DPDstandardNthXt11L,gt11L,kmadd(DPDstandardNthXt21L,gt12L,kmadd(DPDstandardNthXt31L,gt13L,kmadd(Gt211,Gtlu211,kmadd(Gt212,Gtlu212,kmadd(Gt213,Gtlu213,kmadd(Gt311,Gtlu311,kmadd(Gt312,Gtlu312,kmadd(Gt313,Gtlu313,kmadd(Gtl111,Xtn1,kmadd(Gtl112,Xtn2,kmadd(Gtl113,Xtn3,knmsub(DPDstandardNthgt1112L,gtu12,knmsub(DPDstandardNthgt1113L,gtu13,knmsub(DPDstandardNthgt1123L,gtu23,kmadd(kmadd(DPDstandardNthgt1111L,gtu11,kmul(DPDstandardNthgt1122L,gtu22)),ToReal(-0.5),kmadd(DPDstandardNthgt1133L,kmul(gtu33,ToReal(-0.5)),kmadd(kmadd(Gt211,Gtlu121,kmadd(Gt212,Gtlu122,kmadd(Gt213,Gtlu123,kmadd(Gt311,Gtlu131,kmadd(Gt312,Gtlu132,kmul(Gt313,Gtlu133)))))),ToReal(2),kmul(kmadd(Gt111,Gtlu111,kmadd(Gt112,Gtlu112,kmul(Gt113,Gtlu113))),ToReal(3))))))))))))))))))));
    
    CCTK_REAL_VEC Rt12 = 
      kmul(ToReal(0.5),kmadd(DPDstandardNthXt12L,gt11L,kmadd(kadd(DPDstandardNthXt11L,DPDstandardNthXt22L),gt12L,kmadd(DPDstandardNthXt32L,gt13L,kmadd(DPDstandardNthXt21L,gt22L,kmadd(DPDstandardNthXt31L,gt23L,kmadd(Gtl112,Xtn1,kmadd(Gtl211,Xtn1,kmadd(Gtl122,Xtn2,kmadd(Gtl212,Xtn2,kmadd(Gtl123,Xtn3,kmadd(Gtl213,Xtn3,kmadd(DPDstandardNthgt1212L,kmul(gtu12,ToReal(-2)),kmadd(DPDstandardNthgt1213L,kmul(gtu13,ToReal(-2)),kmadd(DPDstandardNthgt1223L,kmul(gtu23,ToReal(-2)),knmsub(DPDstandardNthgt1211L,gtu11,knmsub(DPDstandardNthgt1222L,gtu22,knmsub(DPDstandardNthgt1233L,gtu33,kmadd(kmadd(Gt112,Gtlu111,kmadd(Gt122,Gtlu112,kmadd(Gt123,Gtlu113,kmadd(Gt111,Gtlu121,kmadd(Gt212,Gtlu121,kmadd(Gt112,Gtlu122,kmadd(Gt222,Gtlu122,kmadd(Gt113,Gtlu123,kmadd(Gt223,Gtlu123,kmadd(Gt312,Gtlu131,kmadd(Gt322,Gtlu132,kmadd(Gt323,Gtlu133,kmadd(Gt111,Gtlu211,kmadd(Gt112,Gtlu212,kmadd(Gt113,Gtlu213,kmadd(Gt311,Gtlu231,kmadd(Gt312,Gtlu232,kmadd(Gt313,Gtlu233,kmadd(Gt311,Gtlu321,kmadd(Gt312,Gtlu322,kmul(Gt313,Gtlu323))))))))))))))))))))),ToReal(2),kmul(kmadd(Gt211,Gtlu221,kmadd(Gt212,Gtlu222,kmul(Gt213,Gtlu223))),ToReal(4)))))))))))))))))))));
    
    CCTK_REAL_VEC Rt13 = 
      kmul(ToReal(0.5),kmadd(DPDstandardNthXt13L,gt11L,kmadd(DPDstandardNthXt23L,gt12L,kmadd(kadd(DPDstandardNthXt11L,DPDstandardNthXt33L),gt13L,kmadd(DPDstandardNthXt21L,gt23L,kmadd(DPDstandardNthXt31L,gt33L,kmadd(Gtl113,Xtn1,kmadd(Gtl311,Xtn1,kmadd(Gtl123,Xtn2,kmadd(Gtl312,Xtn2,kmadd(Gtl133,Xtn3,kmadd(Gtl313,Xtn3,kmadd(DPDstandardNthgt1312L,kmul(gtu12,ToReal(-2)),kmadd(DPDstandardNthgt1313L,kmul(gtu13,ToReal(-2)),kmadd(DPDstandardNthgt1323L,kmul(gtu23,ToReal(-2)),knmsub(DPDstandardNthgt1311L,gtu11,knmsub(DPDstandardNthgt1322L,gtu22,knmsub(DPDstandardNthgt1333L,gtu33,kmadd(kmadd(Gt113,Gtlu111,kmadd(Gt123,Gtlu112,kmadd(Gt133,Gtlu113,kmadd(Gt213,Gtlu121,kmadd(Gt223,Gtlu122,kmadd(Gt233,Gtlu123,kmadd(Gt111,Gtlu131,kmadd(Gt313,Gtlu131,kmadd(Gt112,Gtlu132,kmadd(Gt323,Gtlu132,kmadd(Gt113,Gtlu133,kmadd(Gt333,Gtlu133,kmadd(Gt211,Gtlu231,kmadd(Gt212,Gtlu232,kmadd(Gt213,Gtlu233,kmadd(Gt111,Gtlu311,kmadd(Gt112,Gtlu312,kmadd(Gt113,Gtlu313,kmadd(Gt211,Gtlu321,kmadd(Gt212,Gtlu322,kmul(Gt213,Gtlu323))))))))))))))))))))),ToReal(2),kmul(kmadd(Gt311,Gtlu331,kmadd(Gt312,Gtlu332,kmul(Gt313,Gtlu333))),ToReal(4)))))))))))))))))))));
    
    CCTK_REAL_VEC Rt22 = 
      kmadd(DPDstandardNthXt12L,gt12L,kmadd(DPDstandardNthXt22L,gt22L,kmadd(DPDstandardNthXt32L,gt23L,kmadd(Gt312,Gtlu321,kmadd(Gt322,Gtlu322,kmadd(Gt323,Gtlu323,kmadd(Gtl212,Xtn1,kmadd(Gtl222,Xtn2,kmadd(Gtl223,Xtn3,knmsub(DPDstandardNthgt2212L,gtu12,knmsub(DPDstandardNthgt2213L,gtu13,knmsub(DPDstandardNthgt2223L,gtu23,kmadd(kmadd(DPDstandardNthgt2211L,gtu11,kmul(DPDstandardNthgt2222L,gtu22)),ToReal(-0.5),kmadd(DPDstandardNthgt2233L,kmul(gtu33,ToReal(-0.5)),kmadd(kmadd(Gt312,Gtlu231,kmadd(Gt322,Gtlu232,kmul(Gt323,Gtlu233))),ToReal(2),kmadd(Gt112,kmadd(Gtlu211,ToReal(2),Gtlu121),kmadd(Gt122,kmadd(Gtlu212,ToReal(2),Gtlu122),kmadd(Gt123,kmadd(Gtlu213,ToReal(2),Gtlu123),kmul(kmadd(Gt212,Gtlu221,kmadd(Gt222,Gtlu222,kmul(Gt223,Gtlu223))),ToReal(3))))))))))))))))))));
    
    CCTK_REAL_VEC Rt23 = 
      kmul(ToReal(0.5),kmadd(DPDstandardNthXt13L,gt12L,kmadd(DPDstandardNthXt12L,gt13L,kmadd(DPDstandardNthXt23L,gt22L,kmadd(kadd(DPDstandardNthXt22L,DPDstandardNthXt33L),gt23L,kmadd(DPDstandardNthXt32L,gt33L,kmadd(Gtl213,Xtn1,kmadd(Gtl312,Xtn1,kmadd(Gtl223,Xtn2,kmadd(Gtl322,Xtn2,kmadd(Gtl233,Xtn3,kmadd(Gtl323,Xtn3,kmadd(DPDstandardNthgt2312L,kmul(gtu12,ToReal(-2)),kmadd(DPDstandardNthgt2313L,kmul(gtu13,ToReal(-2)),kmadd(DPDstandardNthgt2323L,kmul(gtu23,ToReal(-2)),knmsub(DPDstandardNthgt2311L,gtu11,knmsub(DPDstandardNthgt2322L,gtu22,knmsub(DPDstandardNthgt2333L,gtu33,kmadd(kmadd(Gt112,Gtlu131,kmadd(Gt122,Gtlu132,kmadd(Gt123,Gtlu133,kmadd(Gt113,Gtlu211,kmadd(Gt123,Gtlu212,kmadd(Gt133,Gtlu213,kmadd(Gt213,Gtlu221,kmadd(Gt223,Gtlu222,kmadd(Gt233,Gtlu223,kmadd(Gt212,Gtlu231,kmadd(Gt313,Gtlu231,kmadd(Gt222,Gtlu232,kmadd(Gt323,Gtlu232,kmadd(Gt223,Gtlu233,kmadd(Gt333,Gtlu233,kmadd(Gt112,Gtlu311,kmadd(Gt122,Gtlu312,kmadd(Gt123,Gtlu313,kmadd(Gt212,Gtlu321,kmadd(Gt222,Gtlu322,kmul(Gt223,Gtlu323))))))))))))))))))))),ToReal(2),kmul(kmadd(Gt312,Gtlu331,kmadd(Gt322,Gtlu332,kmul(Gt323,Gtlu333))),ToReal(4)))))))))))))))))))));
    
    CCTK_REAL_VEC Rt33 = 
      kmadd(DPDstandardNthXt13L,gt13L,kmadd(DPDstandardNthXt23L,gt23L,kmadd(DPDstandardNthXt33L,gt33L,kmadd(Gtl313,Xtn1,kmadd(Gtl323,Xtn2,kmadd(Gtl333,Xtn3,knmsub(DPDstandardNthgt3312L,gtu12,knmsub(DPDstandardNthgt3313L,gtu13,knmsub(DPDstandardNthgt3323L,gtu23,kmadd(kmadd(DPDstandardNthgt3311L,gtu11,kmul(DPDstandardNthgt3322L,gtu22)),ToReal(-0.5),kmadd(DPDstandardNthgt3333L,kmul(gtu33,ToReal(-0.5)),kmadd(Gt113,kmadd(Gtlu311,ToReal(2),Gtlu131),kmadd(Gt123,kmadd(Gtlu312,ToReal(2),Gtlu132),kmadd(Gt133,kmadd(Gtlu313,ToReal(2),Gtlu133),kmadd(Gt213,kmadd(Gtlu321,ToReal(2),Gtlu231),kmadd(Gt223,kmadd(Gtlu322,ToReal(2),Gtlu232),kmadd(Gt233,kmadd(Gtlu323,ToReal(2),Gtlu233),kmul(kmadd(Gt313,Gtlu331,kmadd(Gt323,Gtlu332,kmul(Gt333,Gtlu333))),ToReal(3)))))))))))))))))));
    
    CCTK_REAL_VEC fac1 = 
      IfThen(conformalMethod,kmul(INV(phiL),ToReal(-0.5)),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 = kmul(DPDstandardNthphi1L,fac1);
    
    CCTK_REAL_VEC cdphi2 = kmul(DPDstandardNthphi2L,fac1);
    
    CCTK_REAL_VEC cdphi3 = kmul(DPDstandardNthphi3L,fac1);
    
    CCTK_REAL_VEC fac2 = 
      IfThen(conformalMethod,kmul(INV(SQR(phiL)),ToReal(0.5)),ToReal(0));
    
    CCTK_REAL_VEC cdphi211 = 
      kmadd(fac2,SQR(DPDstandardNthphi1L),kmul(fac1,ksub(DPDstandardNthphi11L,kmadd(DPDstandardNthphi1L,Gt111,kmadd(DPDstandardNthphi3L,Gt311,kmul(DPDstandardNthphi2L,Gt211))))));
    
    CCTK_REAL_VEC cdphi212 = 
      kmadd(DPDstandardNthphi1L,kmul(DPDstandardNthphi2L,fac2),kmul(fac1,ksub(DPDstandardNthphi12L,kmadd(DPDstandardNthphi1L,Gt112,kmadd(DPDstandardNthphi3L,Gt312,kmul(DPDstandardNthphi2L,Gt212))))));
    
    CCTK_REAL_VEC cdphi213 = 
      kmadd(DPDstandardNthphi1L,kmul(DPDstandardNthphi3L,fac2),kmul(fac1,ksub(DPDstandardNthphi13L,kmadd(DPDstandardNthphi1L,Gt113,kmadd(DPDstandardNthphi3L,Gt313,kmul(DPDstandardNthphi2L,Gt213))))));
    
    CCTK_REAL_VEC cdphi222 = 
      kmadd(fac2,SQR(DPDstandardNthphi2L),kmul(fac1,ksub(DPDstandardNthphi22L,kmadd(DPDstandardNthphi1L,Gt122,kmadd(DPDstandardNthphi3L,Gt322,kmul(DPDstandardNthphi2L,Gt222))))));
    
    CCTK_REAL_VEC cdphi223 = 
      kmadd(DPDstandardNthphi2L,kmul(DPDstandardNthphi3L,fac2),kmul(fac1,ksub(DPDstandardNthphi23L,kmadd(DPDstandardNthphi1L,Gt123,kmadd(DPDstandardNthphi3L,Gt323,kmul(DPDstandardNthphi2L,Gt223))))));
    
    CCTK_REAL_VEC cdphi233 = 
      kmadd(fac2,SQR(DPDstandardNthphi3L),kmul(fac1,ksub(DPDstandardNthphi33L,kmadd(DPDstandardNthphi1L,Gt133,kmadd(DPDstandardNthphi3L,Gt333,kmul(DPDstandardNthphi2L,Gt233))))));
    
    CCTK_REAL_VEC Rphi11 = 
      kmul(ToReal(-2),kadd(cdphi211,kmadd(SQR(cdphi1),kmul(kmadd(gt11L,gtu11,ToReal(-1)),ToReal(2)),kmul(gt11L,kmadd(cdphi211,gtu11,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu33,SQR(cdphi3))))),ToReal(2),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmul(kmadd(cdphi1,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13)),kmul(cdphi2,kmul(cdphi3,gtu23))),ToReal(4))))))))));
    
    CCTK_REAL_VEC Rphi12 = 
      kmul(ToReal(-2),kadd(cdphi212,kmadd(gt12L,kmadd(cdphi211,gtu11,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu11,SQR(cdphi1))))),ToReal(2),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmadd(gtu33,kmadd(SQR(cdphi3),ToReal(2),cdphi233),kmul(cdphi2,kmul(cdphi3,kmul(gtu23,ToReal(4)))))))),kmul(cdphi1,kmadd(gt12L,kmul(cdphi3,kmul(gtu13,ToReal(4))),kmul(cdphi2,kmadd(gt12L,kmul(gtu12,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi13 = 
      kmul(ToReal(-2),kadd(cdphi213,kmadd(gt13L,kmadd(cdphi211,gtu11,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu11,SQR(cdphi1))))),ToReal(2),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmadd(gtu33,kmadd(SQR(cdphi3),ToReal(2),cdphi233),kmul(cdphi2,kmul(cdphi3,kmul(gtu23,ToReal(4)))))))),kmul(cdphi1,kmadd(gt13L,kmul(cdphi2,kmul(gtu12,ToReal(4))),kmul(cdphi3,kmadd(gt13L,kmul(gtu13,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi22 = 
      kmul(ToReal(-2),kadd(cdphi222,kmadd(SQR(cdphi2),kmul(kmadd(gt22L,gtu22,ToReal(-1)),ToReal(2)),kmul(gt22L,kmadd(cdphi222,gtu22,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu33,SQR(cdphi3))))),ToReal(2),kmadd(gtu11,kmadd(SQR(cdphi1),ToReal(2),cdphi211),kmul(kmadd(cdphi1,kmul(cdphi3,gtu13),kmul(cdphi2,kmadd(cdphi1,gtu12,kmul(cdphi3,gtu23)))),ToReal(4))))))))));
    
    CCTK_REAL_VEC Rphi23 = 
      kmul(ToReal(-2),kadd(cdphi223,kmadd(gt23L,kmadd(cdphi222,gtu22,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu22,SQR(cdphi2))))),ToReal(2),kmadd(gtu11,kmadd(SQR(cdphi1),ToReal(2),cdphi211),kmadd(gtu33,kmadd(SQR(cdphi3),ToReal(2),cdphi233),kmul(cdphi1,kmul(cdphi3,kmul(gtu13,ToReal(4)))))))),kmul(cdphi2,kmadd(gt23L,kmul(cdphi1,kmul(gtu12,ToReal(4))),kmul(cdphi3,kmadd(gt23L,kmul(gtu23,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi33 = 
      kmul(ToReal(-2),kadd(cdphi233,kmadd(SQR(cdphi3),kmul(kmadd(gt33L,gtu33,ToReal(-1)),ToReal(2)),kmul(gt33L,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi213,gtu13,kmul(cdphi223,gtu23)),ToReal(2),kmadd(gtu11,kmadd(SQR(cdphi1),ToReal(2),cdphi211),kmadd(gtu22,kmadd(SQR(cdphi2),ToReal(2),cdphi222),kmadd(cdphi3,kmul(kmadd(cdphi1,gtu13,kmul(cdphi2,gtu23)),ToReal(4)),kmul(gtu12,kmadd(cdphi212,ToReal(2),kmul(cdphi1,kmul(cdphi2,ToReal(4))))))))))))));
    
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
      IfThen(conformalMethod,INV(SQR(phiL)),kexp(kmul(phiL,ToReal(4))));
    
    CCTK_REAL_VEC em4phi = INV(e4phi);
    
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
      kmadd(DPDstandardNthalpha2L,Gt211,kmadd(DPDstandardNthalpha3L,Gt311,kmadd(alphaL,R11,kmsub(DPDstandardNthalpha1L,kmadd(cdphi1,ToReal(4),Gt111),DPDstandardNthalpha11L))));
    
    CCTK_REAL_VEC Ats12 = 
      kmadd(DPDstandardNthalpha1L,Gt112,kmadd(DPDstandardNthalpha2L,Gt212,kmadd(DPDstandardNthalpha3L,Gt312,kmadd(alphaL,R12,kmsub(kmadd(DPDstandardNthalpha2L,cdphi1,kmul(DPDstandardNthalpha1L,cdphi2)),ToReal(2),DPDstandardNthalpha12L)))));
    
    CCTK_REAL_VEC Ats13 = 
      kmadd(DPDstandardNthalpha1L,Gt113,kmadd(DPDstandardNthalpha2L,Gt213,kmadd(DPDstandardNthalpha3L,Gt313,kmadd(alphaL,R13,kmsub(kmadd(DPDstandardNthalpha3L,cdphi1,kmul(DPDstandardNthalpha1L,cdphi3)),ToReal(2),DPDstandardNthalpha13L)))));
    
    CCTK_REAL_VEC Ats22 = 
      kmadd(DPDstandardNthalpha1L,Gt122,kmadd(DPDstandardNthalpha3L,Gt322,kmadd(alphaL,R22,kmsub(DPDstandardNthalpha2L,kmadd(cdphi2,ToReal(4),Gt222),DPDstandardNthalpha22L))));
    
    CCTK_REAL_VEC Ats23 = 
      kmadd(DPDstandardNthalpha1L,Gt123,kmadd(DPDstandardNthalpha2L,Gt223,kmadd(DPDstandardNthalpha3L,Gt323,kmadd(alphaL,R23,kmsub(kmadd(DPDstandardNthalpha3L,cdphi2,kmul(DPDstandardNthalpha2L,cdphi3)),ToReal(2),DPDstandardNthalpha23L)))));
    
    CCTK_REAL_VEC Ats33 = 
      kmadd(DPDstandardNthalpha1L,Gt133,kmadd(DPDstandardNthalpha2L,Gt233,kmadd(alphaL,R33,kmsub(DPDstandardNthalpha3L,kmadd(cdphi3,ToReal(4),Gt333),DPDstandardNthalpha33L))));
    
    CCTK_REAL_VEC trAts = 
      kmadd(Ats11,gu11,kmadd(Ats22,gu22,kmadd(Ats33,gu33,kmul(kmadd(Ats12,gu12,kmadd(Ats13,gu13,kmul(Ats23,gu23))),ToReal(2)))));
    
    CCTK_REAL_VEC At11rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(ToReal(-3),kmadd(kmadd(At12L,Atm21,kmul(At13L,Atm31)),ToReal(2),kmul(At11L,kmsub(Atm11,ToReal(2),trKL)))),kmadd(em4phi,kmsub(Ats11,ToReal(3),kmul(g11,trAts)),kmadd(At11L,kmadd(kadd(DPDstandardNthbeta22L,DPDstandardNthbeta33L),ToReal(-2),kmul(DPDstandardNthbeta11L,ToReal(4))),kmul(kmadd(At12L,DPDstandardNthbeta21L,kmul(At13L,DPDstandardNthbeta31L)),ToReal(6))))));
    
    CCTK_REAL_VEC At12rhsL = 
      kmadd(alphaL,kmadd(kmadd(At11L,Atm12,kmul(At13L,Atm32)),ToReal(-2),kmul(At12L,kmadd(Atm22,ToReal(-2),trKL))),kmul(ToReal(0.333333333333333333333333333333),kmadd(At12L,kadd(DPDstandardNthbeta11L,kmadd(DPDstandardNthbeta33L,ToReal(-2),DPDstandardNthbeta22L)),kmsub(kmadd(At11L,DPDstandardNthbeta12L,kmadd(At22L,DPDstandardNthbeta21L,kmadd(At23L,DPDstandardNthbeta31L,kmadd(At13L,DPDstandardNthbeta32L,kmul(Ats12,em4phi))))),ToReal(3),kmul(em4phi,kmul(g12,trAts))))));
    
    CCTK_REAL_VEC At13rhsL = 
      kmadd(alphaL,kmadd(kmadd(At11L,Atm13,kmul(At12L,Atm23)),ToReal(-2),kmul(At13L,kmadd(Atm33,ToReal(-2),trKL))),kmul(ToReal(0.333333333333333333333333333333),kmadd(At13L,kadd(DPDstandardNthbeta11L,kmadd(DPDstandardNthbeta22L,ToReal(-2),DPDstandardNthbeta33L)),kmsub(kmadd(At11L,DPDstandardNthbeta13L,kmadd(At23L,DPDstandardNthbeta21L,kmadd(At12L,DPDstandardNthbeta23L,kmadd(At33L,DPDstandardNthbeta31L,kmul(Ats13,em4phi))))),ToReal(3),kmul(em4phi,kmul(g13,trAts))))));
    
    CCTK_REAL_VEC At22rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(At22L,kmul(ToReal(-2),kadd(DPDstandardNthbeta11L,kmadd(DPDstandardNthbeta22L,ToReal(-2),DPDstandardNthbeta33L))),kmadd(alphaL,kmul(ToReal(-3),kmadd(kmadd(At12L,Atm12,kmul(At23L,Atm32)),ToReal(2),kmul(At22L,kmsub(Atm22,ToReal(2),trKL)))),kmadd(em4phi,kmsub(Ats22,ToReal(3),kmul(g22,trAts)),kmul(kmadd(At12L,DPDstandardNthbeta12L,kmul(At23L,DPDstandardNthbeta32L)),ToReal(6))))));
    
    CCTK_REAL_VEC At23rhsL = 
      kmadd(alphaL,kmadd(kmadd(At12L,Atm13,kmul(At22L,Atm23)),ToReal(-2),kmul(At23L,kmadd(Atm33,ToReal(-2),trKL))),kmul(ToReal(0.333333333333333333333333333333),kmadd(At23L,kadd(DPDstandardNthbeta22L,kmadd(DPDstandardNthbeta11L,ToReal(-2),DPDstandardNthbeta33L)),kmsub(kmadd(At13L,DPDstandardNthbeta12L,kmadd(At12L,DPDstandardNthbeta13L,kmadd(At22L,DPDstandardNthbeta23L,kmadd(At33L,DPDstandardNthbeta32L,kmul(Ats23,em4phi))))),ToReal(3),kmul(em4phi,kmul(g23,trAts))))));
    
    CCTK_REAL_VEC At33rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(At33L,kmul(ToReal(-2),kadd(DPDstandardNthbeta11L,kmadd(DPDstandardNthbeta33L,ToReal(-2),DPDstandardNthbeta22L))),kmadd(alphaL,kmul(ToReal(-3),kmadd(kmadd(At13L,Atm13,kmul(At23L,Atm23)),ToReal(2),kmul(At33L,kmsub(Atm33,ToReal(2),trKL)))),kmadd(em4phi,kmsub(Ats33,ToReal(3),kmul(g33,trAts)),kmul(kmadd(At13L,DPDstandardNthbeta13L,kmul(At23L,DPDstandardNthbeta23L)),ToReal(6))))));
    
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
    "ML_BSSN_Host::DPDstandardNthalpha11_group",
    "ML_BSSN_Host::DPDstandardNthalpha12_group",
    "ML_BSSN_Host::DPDstandardNthalpha13_group",
    "ML_BSSN_Host::DPDstandardNthalpha1_group",
    "ML_BSSN_Host::DPDstandardNthalpha22_group",
    "ML_BSSN_Host::DPDstandardNthalpha23_group",
    "ML_BSSN_Host::DPDstandardNthalpha2_group",
    "ML_BSSN_Host::DPDstandardNthalpha33_group",
    "ML_BSSN_Host::DPDstandardNthalpha3_group",
    "ML_BSSN_Host::DPDstandardNthbeta11_group",
    "ML_BSSN_Host::DPDstandardNthbeta12_group",
    "ML_BSSN_Host::DPDstandardNthbeta13_group",
    "ML_BSSN_Host::DPDstandardNthbeta21_group",
    "ML_BSSN_Host::DPDstandardNthbeta22_group",
    "ML_BSSN_Host::DPDstandardNthbeta23_group",
    "ML_BSSN_Host::DPDstandardNthbeta31_group",
    "ML_BSSN_Host::DPDstandardNthbeta32_group",
    "ML_BSSN_Host::DPDstandardNthbeta33_group",
    "ML_BSSN_Host::DPDstandardNthgt1111_group",
    "ML_BSSN_Host::DPDstandardNthgt1112_group",
    "ML_BSSN_Host::DPDstandardNthgt1113_group",
    "ML_BSSN_Host::DPDstandardNthgt111_group",
    "ML_BSSN_Host::DPDstandardNthgt1122_group",
    "ML_BSSN_Host::DPDstandardNthgt1123_group",
    "ML_BSSN_Host::DPDstandardNthgt112_group",
    "ML_BSSN_Host::DPDstandardNthgt1133_group",
    "ML_BSSN_Host::DPDstandardNthgt113_group",
    "ML_BSSN_Host::DPDstandardNthgt1211_group",
    "ML_BSSN_Host::DPDstandardNthgt1212_group",
    "ML_BSSN_Host::DPDstandardNthgt1213_group",
    "ML_BSSN_Host::DPDstandardNthgt121_group",
    "ML_BSSN_Host::DPDstandardNthgt1222_group",
    "ML_BSSN_Host::DPDstandardNthgt1223_group",
    "ML_BSSN_Host::DPDstandardNthgt122_group",
    "ML_BSSN_Host::DPDstandardNthgt1233_group",
    "ML_BSSN_Host::DPDstandardNthgt123_group",
    "ML_BSSN_Host::DPDstandardNthgt1311_group",
    "ML_BSSN_Host::DPDstandardNthgt1312_group",
    "ML_BSSN_Host::DPDstandardNthgt1313_group",
    "ML_BSSN_Host::DPDstandardNthgt131_group",
    "ML_BSSN_Host::DPDstandardNthgt1322_group",
    "ML_BSSN_Host::DPDstandardNthgt1323_group",
    "ML_BSSN_Host::DPDstandardNthgt132_group",
    "ML_BSSN_Host::DPDstandardNthgt1333_group",
    "ML_BSSN_Host::DPDstandardNthgt133_group",
    "ML_BSSN_Host::DPDstandardNthgt2211_group",
    "ML_BSSN_Host::DPDstandardNthgt2212_group",
    "ML_BSSN_Host::DPDstandardNthgt2213_group",
    "ML_BSSN_Host::DPDstandardNthgt221_group",
    "ML_BSSN_Host::DPDstandardNthgt2222_group",
    "ML_BSSN_Host::DPDstandardNthgt2223_group",
    "ML_BSSN_Host::DPDstandardNthgt222_group",
    "ML_BSSN_Host::DPDstandardNthgt2233_group",
    "ML_BSSN_Host::DPDstandardNthgt223_group",
    "ML_BSSN_Host::DPDstandardNthgt2311_group",
    "ML_BSSN_Host::DPDstandardNthgt2312_group",
    "ML_BSSN_Host::DPDstandardNthgt2313_group",
    "ML_BSSN_Host::DPDstandardNthgt231_group",
    "ML_BSSN_Host::DPDstandardNthgt2322_group",
    "ML_BSSN_Host::DPDstandardNthgt2323_group",
    "ML_BSSN_Host::DPDstandardNthgt232_group",
    "ML_BSSN_Host::DPDstandardNthgt2333_group",
    "ML_BSSN_Host::DPDstandardNthgt233_group",
    "ML_BSSN_Host::DPDstandardNthgt3311_group",
    "ML_BSSN_Host::DPDstandardNthgt3312_group",
    "ML_BSSN_Host::DPDstandardNthgt3313_group",
    "ML_BSSN_Host::DPDstandardNthgt331_group",
    "ML_BSSN_Host::DPDstandardNthgt3322_group",
    "ML_BSSN_Host::DPDstandardNthgt3323_group",
    "ML_BSSN_Host::DPDstandardNthgt332_group",
    "ML_BSSN_Host::DPDstandardNthgt3333_group",
    "ML_BSSN_Host::DPDstandardNthgt333_group",
    "ML_BSSN_Host::DPDstandardNthphi11_group",
    "ML_BSSN_Host::DPDstandardNthphi12_group",
    "ML_BSSN_Host::DPDstandardNthphi13_group",
    "ML_BSSN_Host::DPDstandardNthphi1_group",
    "ML_BSSN_Host::DPDstandardNthphi22_group",
    "ML_BSSN_Host::DPDstandardNthphi23_group",
    "ML_BSSN_Host::DPDstandardNthphi2_group",
    "ML_BSSN_Host::DPDstandardNthphi33_group",
    "ML_BSSN_Host::DPDstandardNthphi3_group",
    "ML_BSSN_Host::DPDstandardNthXt11_group",
    "ML_BSSN_Host::DPDstandardNthXt12_group",
    "ML_BSSN_Host::DPDstandardNthXt13_group",
    "ML_BSSN_Host::DPDstandardNthXt21_group",
    "ML_BSSN_Host::DPDstandardNthXt22_group",
    "ML_BSSN_Host::DPDstandardNthXt23_group",
    "ML_BSSN_Host::DPDstandardNthXt31_group",
    "ML_BSSN_Host::DPDstandardNthXt32_group",
    "ML_BSSN_Host::DPDstandardNthXt33_group",
    "ML_BSSN_Host::ML_curv",
    "ML_BSSN_Host::ML_curvrhs",
    "ML_BSSN_Host::ML_lapse",
    "ML_BSSN_Host::ML_log_confac",
    "ML_BSSN_Host::ML_metric",
    "ML_BSSN_Host::ML_shift",
    "ML_BSSN_Host::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_RHS2", 97, groups);
  
  
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_Host_RHS2_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_RHS2_Body");
  }
}
