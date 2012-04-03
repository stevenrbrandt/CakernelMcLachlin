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

extern "C" void ML_BSSN_Host_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_trace_curvrhs.");
  return;
}

static void ML_BSSN_Host_RHS1_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  LC_LOOP3VEC(ML_BSSN_Host_RHS1,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
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
    CCTK_REAL_VEC DPDstandardNthbeta111L = vec_load(DPDstandardNthbeta111[index]);
    CCTK_REAL_VEC DPDstandardNthbeta112L = vec_load(DPDstandardNthbeta112[index]);
    CCTK_REAL_VEC DPDstandardNthbeta113L = vec_load(DPDstandardNthbeta113[index]);
    CCTK_REAL_VEC DPDstandardNthbeta12L = vec_load(DPDstandardNthbeta12[index]);
    CCTK_REAL_VEC DPDstandardNthbeta122L = vec_load(DPDstandardNthbeta122[index]);
    CCTK_REAL_VEC DPDstandardNthbeta123L = vec_load(DPDstandardNthbeta123[index]);
    CCTK_REAL_VEC DPDstandardNthbeta13L = vec_load(DPDstandardNthbeta13[index]);
    CCTK_REAL_VEC DPDstandardNthbeta133L = vec_load(DPDstandardNthbeta133[index]);
    CCTK_REAL_VEC DPDstandardNthbeta21L = vec_load(DPDstandardNthbeta21[index]);
    CCTK_REAL_VEC DPDstandardNthbeta211L = vec_load(DPDstandardNthbeta211[index]);
    CCTK_REAL_VEC DPDstandardNthbeta212L = vec_load(DPDstandardNthbeta212[index]);
    CCTK_REAL_VEC DPDstandardNthbeta213L = vec_load(DPDstandardNthbeta213[index]);
    CCTK_REAL_VEC DPDstandardNthbeta22L = vec_load(DPDstandardNthbeta22[index]);
    CCTK_REAL_VEC DPDstandardNthbeta222L = vec_load(DPDstandardNthbeta222[index]);
    CCTK_REAL_VEC DPDstandardNthbeta223L = vec_load(DPDstandardNthbeta223[index]);
    CCTK_REAL_VEC DPDstandardNthbeta23L = vec_load(DPDstandardNthbeta23[index]);
    CCTK_REAL_VEC DPDstandardNthbeta233L = vec_load(DPDstandardNthbeta233[index]);
    CCTK_REAL_VEC DPDstandardNthbeta31L = vec_load(DPDstandardNthbeta31[index]);
    CCTK_REAL_VEC DPDstandardNthbeta311L = vec_load(DPDstandardNthbeta311[index]);
    CCTK_REAL_VEC DPDstandardNthbeta312L = vec_load(DPDstandardNthbeta312[index]);
    CCTK_REAL_VEC DPDstandardNthbeta313L = vec_load(DPDstandardNthbeta313[index]);
    CCTK_REAL_VEC DPDstandardNthbeta32L = vec_load(DPDstandardNthbeta32[index]);
    CCTK_REAL_VEC DPDstandardNthbeta322L = vec_load(DPDstandardNthbeta322[index]);
    CCTK_REAL_VEC DPDstandardNthbeta323L = vec_load(DPDstandardNthbeta323[index]);
    CCTK_REAL_VEC DPDstandardNthbeta33L = vec_load(DPDstandardNthbeta33[index]);
    CCTK_REAL_VEC DPDstandardNthbeta333L = vec_load(DPDstandardNthbeta333[index]);
    CCTK_REAL_VEC DPDstandardNthgt111L = vec_load(DPDstandardNthgt111[index]);
    CCTK_REAL_VEC DPDstandardNthgt112L = vec_load(DPDstandardNthgt112[index]);
    CCTK_REAL_VEC DPDstandardNthgt113L = vec_load(DPDstandardNthgt113[index]);
    CCTK_REAL_VEC DPDstandardNthgt121L = vec_load(DPDstandardNthgt121[index]);
    CCTK_REAL_VEC DPDstandardNthgt122L = vec_load(DPDstandardNthgt122[index]);
    CCTK_REAL_VEC DPDstandardNthgt123L = vec_load(DPDstandardNthgt123[index]);
    CCTK_REAL_VEC DPDstandardNthgt131L = vec_load(DPDstandardNthgt131[index]);
    CCTK_REAL_VEC DPDstandardNthgt132L = vec_load(DPDstandardNthgt132[index]);
    CCTK_REAL_VEC DPDstandardNthgt133L = vec_load(DPDstandardNthgt133[index]);
    CCTK_REAL_VEC DPDstandardNthgt221L = vec_load(DPDstandardNthgt221[index]);
    CCTK_REAL_VEC DPDstandardNthgt222L = vec_load(DPDstandardNthgt222[index]);
    CCTK_REAL_VEC DPDstandardNthgt223L = vec_load(DPDstandardNthgt223[index]);
    CCTK_REAL_VEC DPDstandardNthgt231L = vec_load(DPDstandardNthgt231[index]);
    CCTK_REAL_VEC DPDstandardNthgt232L = vec_load(DPDstandardNthgt232[index]);
    CCTK_REAL_VEC DPDstandardNthgt233L = vec_load(DPDstandardNthgt233[index]);
    CCTK_REAL_VEC DPDstandardNthgt331L = vec_load(DPDstandardNthgt331[index]);
    CCTK_REAL_VEC DPDstandardNthgt332L = vec_load(DPDstandardNthgt332[index]);
    CCTK_REAL_VEC DPDstandardNthgt333L = vec_load(DPDstandardNthgt333[index]);
    CCTK_REAL_VEC DPDstandardNthphi1L = vec_load(DPDstandardNthphi1[index]);
    CCTK_REAL_VEC DPDstandardNthphi2L = vec_load(DPDstandardNthphi2[index]);
    CCTK_REAL_VEC DPDstandardNthphi3L = vec_load(DPDstandardNthphi3[index]);
    CCTK_REAL_VEC DPDstandardNthtrK1L = vec_load(DPDstandardNthtrK1[index]);
    CCTK_REAL_VEC DPDstandardNthtrK2L = vec_load(DPDstandardNthtrK2[index]);
    CCTK_REAL_VEC DPDstandardNthtrK3L = vec_load(DPDstandardNthtrK3[index]);
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
    
    CCTK_REAL_VEC cdphi1 = kmul(DPDstandardNthphi1L,fac1);
    
    CCTK_REAL_VEC cdphi2 = kmul(DPDstandardNthphi2L,fac1);
    
    CCTK_REAL_VEC cdphi3 = kmul(DPDstandardNthphi3L,fac1);
    
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
    
    CCTK_REAL_VEC phirhsL = 
      IfThen(conformalMethod,kmul(phiL,kmadd(kadd(DPDstandardNthbeta11L,kadd(DPDstandardNthbeta22L,DPDstandardNthbeta33L)),ToReal(-0.333333333333333333333333333333),kmul(alphaL,kmul(trKL,ToReal(0.333333333333333333333333333333))))),kmadd(alphaL,kmul(trKL,ToReal(-0.166666666666666666666666666667)),kmul(kadd(DPDstandardNthbeta11L,kadd(DPDstandardNthbeta22L,DPDstandardNthbeta33L)),ToReal(0.166666666666666666666666666667))));
    
    CCTK_REAL_VEC gt11rhsL = 
      kmul(ToReal(-0.666666666666666666666666666667),kmadd(kmadd(DPDstandardNthbeta21L,gt12L,kmul(DPDstandardNthbeta31L,gt13L)),ToReal(-3),kmadd(gt11L,kadd(DPDstandardNthbeta22L,kmadd(DPDstandardNthbeta11L,ToReal(-2),DPDstandardNthbeta33L)),kmul(alphaL,kmul(At11L,ToReal(3))))));
    
    CCTK_REAL_VEC gt12rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At12L,ToReal(-6)),kmadd(gt12L,kadd(DPDstandardNthbeta11L,kmadd(DPDstandardNthbeta33L,ToReal(-2),DPDstandardNthbeta22L)),kmul(kmadd(DPDstandardNthbeta12L,gt11L,kmadd(DPDstandardNthbeta32L,gt13L,kmadd(DPDstandardNthbeta21L,gt22L,kmul(DPDstandardNthbeta31L,gt23L)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt13rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At13L,ToReal(-6)),kmadd(gt13L,kadd(DPDstandardNthbeta11L,kmadd(DPDstandardNthbeta22L,ToReal(-2),DPDstandardNthbeta33L)),kmul(kmadd(DPDstandardNthbeta13L,gt11L,kmadd(DPDstandardNthbeta23L,gt12L,kmadd(DPDstandardNthbeta21L,gt23L,kmul(DPDstandardNthbeta31L,gt33L)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt22rhsL = 
      kmul(ToReal(-0.666666666666666666666666666667),kmadd(kmadd(DPDstandardNthbeta12L,gt12L,kmul(DPDstandardNthbeta32L,gt23L)),ToReal(-3),kmadd(gt22L,kadd(DPDstandardNthbeta11L,kmadd(DPDstandardNthbeta22L,ToReal(-2),DPDstandardNthbeta33L)),kmul(alphaL,kmul(At22L,ToReal(3))))));
    
    CCTK_REAL_VEC gt23rhsL = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At23L,ToReal(-6)),kmadd(gt23L,kadd(DPDstandardNthbeta22L,kmadd(DPDstandardNthbeta11L,ToReal(-2),DPDstandardNthbeta33L)),kmul(kmadd(DPDstandardNthbeta13L,gt12L,kmadd(DPDstandardNthbeta12L,gt13L,kmadd(DPDstandardNthbeta23L,gt22L,kmul(DPDstandardNthbeta32L,gt33L)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt33rhsL = 
      kmul(ToReal(-0.666666666666666666666666666667),kmadd(kmadd(DPDstandardNthbeta13L,gt13L,kmul(DPDstandardNthbeta23L,gt23L)),ToReal(-3),kmadd(gt33L,kadd(DPDstandardNthbeta11L,kmadd(DPDstandardNthbeta33L,ToReal(-2),DPDstandardNthbeta22L)),kmul(alphaL,kmul(At33L,ToReal(3))))));
    
    CCTK_REAL_VEC dotXt1 = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(kmadd(DPDstandardNthalpha1L,Atu11,kmadd(DPDstandardNthalpha2L,Atu12,kmul(DPDstandardNthalpha3L,Atu13))),ToReal(-6),kmadd(kmadd(DPDstandardNthbeta12L,Xtn2,kmul(DPDstandardNthbeta13L,Xtn3)),ToReal(-3),kmadd(Xtn1,kmsub(DPDstandardNthbeta33L,ToReal(2),DPDstandardNthbeta11L),kmadd(kmadd(DPDstandardNthbeta122L,gtu22,kmul(DPDstandardNthbeta133L,gtu33)),ToReal(3),kmadd(gtu11,kadd(DPDstandardNthbeta212L,kmadd(DPDstandardNthbeta111L,ToReal(4),DPDstandardNthbeta313L)),kmadd(DPDstandardNthbeta123L,kmul(gtu23,ToReal(6)),kmadd(gtu12,kadd(DPDstandardNthbeta222L,kmadd(DPDstandardNthbeta112L,ToReal(7),DPDstandardNthbeta323L)),kmadd(gtu13,kadd(DPDstandardNthbeta223L,kmadd(DPDstandardNthbeta113L,ToReal(7),DPDstandardNthbeta333L)),kmul(ToReal(2),kmadd(DPDstandardNthbeta22L,Xtn1,kmul(alphaL,kmadd(kmadd(DPDstandardNthtrK1L,gtu11,kmadd(DPDstandardNthtrK2L,gtu12,kmul(DPDstandardNthtrK3L,gtu13))),ToReal(-2),kmadd(kmadd(Atu23,Gt123,kmul(Atu12,kmadd(cdphi2,ToReal(3),Gt112))),ToReal(6),kmadd(ToReal(3),kmadd(Atu22,Gt122,kmadd(Atu33,Gt133,kmul(Atu11,kmadd(cdphi1,ToReal(6),Gt111)))),kmul(Atu13,kmadd(Gt113,ToReal(6),kmul(cdphi3,ToReal(18)))))))))))))))))));
    
    CCTK_REAL_VEC dotXt2 = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(kmadd(DPDstandardNthalpha1L,Atu12,kmadd(DPDstandardNthalpha2L,Atu22,kmul(DPDstandardNthalpha3L,Atu23))),ToReal(-6),kmadd(kmadd(DPDstandardNthbeta21L,Xtn1,kmul(DPDstandardNthbeta23L,Xtn3)),ToReal(-3),kmadd(Xtn2,kmsub(DPDstandardNthbeta33L,ToReal(2),DPDstandardNthbeta22L),kmadd(kmadd(DPDstandardNthbeta211L,gtu11,kmul(DPDstandardNthbeta233L,gtu33)),ToReal(3),kmadd(gtu22,kadd(DPDstandardNthbeta112L,kmadd(DPDstandardNthbeta222L,ToReal(4),DPDstandardNthbeta323L)),kmadd(DPDstandardNthbeta213L,kmul(gtu13,ToReal(6)),kmadd(gtu12,kadd(DPDstandardNthbeta111L,kmadd(DPDstandardNthbeta212L,ToReal(7),DPDstandardNthbeta313L)),kmadd(gtu23,kadd(DPDstandardNthbeta113L,kmadd(DPDstandardNthbeta223L,ToReal(7),DPDstandardNthbeta333L)),kmul(ToReal(2),kmadd(DPDstandardNthbeta11L,Xtn2,kmul(alphaL,kmadd(kmadd(DPDstandardNthtrK1L,gtu12,kmadd(DPDstandardNthtrK2L,gtu22,kmul(DPDstandardNthtrK3L,gtu23))),ToReal(-2),kmadd(kmadd(Atu13,Gt213,kmul(Atu12,kmadd(cdphi1,ToReal(3),Gt212))),ToReal(6),kmadd(ToReal(3),kmadd(Atu11,Gt211,kmadd(Atu33,Gt233,kmul(Atu22,kmadd(cdphi2,ToReal(6),Gt222)))),kmul(Atu23,kmadd(Gt223,ToReal(6),kmul(cdphi3,ToReal(18)))))))))))))))))));
    
    CCTK_REAL_VEC dotXt3 = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(kmadd(DPDstandardNthalpha1L,Atu13,kmadd(DPDstandardNthalpha2L,Atu23,kmul(DPDstandardNthalpha3L,Atu33))),ToReal(-6),kmadd(kmadd(DPDstandardNthbeta31L,Xtn1,kmul(DPDstandardNthbeta32L,Xtn2)),ToReal(-3),kmadd(Xtn3,kmsub(DPDstandardNthbeta22L,ToReal(2),DPDstandardNthbeta33L),kmadd(kmadd(DPDstandardNthbeta311L,gtu11,kmul(DPDstandardNthbeta322L,gtu22)),ToReal(3),kmadd(gtu33,kadd(DPDstandardNthbeta113L,kmadd(DPDstandardNthbeta333L,ToReal(4),DPDstandardNthbeta223L)),kmadd(DPDstandardNthbeta312L,kmul(gtu12,ToReal(6)),kmadd(gtu13,kadd(DPDstandardNthbeta111L,kmadd(DPDstandardNthbeta313L,ToReal(7),DPDstandardNthbeta212L)),kmadd(gtu23,kadd(DPDstandardNthbeta112L,kmadd(DPDstandardNthbeta323L,ToReal(7),DPDstandardNthbeta222L)),kmul(ToReal(2),kmadd(DPDstandardNthbeta11L,Xtn3,kmul(alphaL,kmadd(kmadd(DPDstandardNthtrK1L,gtu13,kmadd(DPDstandardNthtrK2L,gtu23,kmul(DPDstandardNthtrK3L,gtu33))),ToReal(-2),kmadd(kmadd(Atu11,Gt311,kmul(Atu22,Gt322)),ToReal(3),kmadd(kmadd(Atu12,Gt312,kmadd(Atu13,kmadd(cdphi1,ToReal(3),Gt313),kmul(Atu23,kmadd(cdphi2,ToReal(3),Gt323)))),ToReal(6),kmul(Atu33,kmadd(Gt333,ToReal(3),kmul(cdphi3,ToReal(18)))))))))))))))))));
    
    CCTK_REAL_VEC Xt1rhsL = dotXt1;
    
    CCTK_REAL_VEC Xt2rhsL = dotXt2;
    
    CCTK_REAL_VEC Xt3rhsL = dotXt3;
    
    CCTK_REAL_VEC dottrK = 
      kmsub(alphaL,kadd(SQR(Atm11),kadd(SQR(Atm22),kadd(SQR(Atm33),kmadd(SQR(trKL),ToReal(0.333333333333333333333333333333),kmul(kmadd(Atm12,Atm21,kmadd(Atm13,Atm31,kmul(Atm23,Atm32))),ToReal(2)))))),kmul(em4phi,kmadd(DPDstandardNthalpha11L,gtu11,kmadd(DPDstandardNthalpha22L,gtu22,kmadd(DPDstandardNthalpha33L,gtu33,knmsub(DPDstandardNthalpha1L,Xtn1,knmsub(DPDstandardNthalpha2L,Xtn2,kmsub(kmadd(DPDstandardNthalpha12L,gtu12,kmadd(DPDstandardNthalpha1L,kmul(cdphi2,gtu12),kmadd(DPDstandardNthalpha13L,gtu13,kmadd(DPDstandardNthalpha1L,kmul(cdphi3,gtu13),kmadd(cdphi1,kmadd(DPDstandardNthalpha1L,gtu11,kmadd(DPDstandardNthalpha2L,gtu12,kmul(DPDstandardNthalpha3L,gtu13))),kmadd(DPDstandardNthalpha2L,kmul(cdphi2,gtu22),kmadd(DPDstandardNthalpha23L,gtu23,kmadd(DPDstandardNthalpha3L,kmul(cdphi2,gtu23),kmadd(DPDstandardNthalpha2L,kmul(cdphi3,gtu23),kmul(DPDstandardNthalpha3L,kmul(cdphi3,gtu33))))))))))),ToReal(2),kmul(DPDstandardNthalpha3L,Xtn3)))))))));
    
    CCTK_REAL_VEC trKrhsL = dottrK;
    
    CCTK_REAL_VEC alpharhsL = 
      kneg(kmul(kpow(alphaL,harmonicN),kmul(ToReal(harmonicF),kmadd(AL,ToReal(LapseACoeff),kmul(kmadd(kadd(alphaL,ToReal(-1)),ToReal(AlphaDriver),trKL),ksub(ToReal(1),ToReal(LapseACoeff)))))));
    
    CCTK_REAL_VEC ArhsL = 
      kmul(knmsub(AL,ToReal(AlphaDriver),dottrK),ToReal(LapseACoeff));
    
    CCTK_REAL_VEC eta = ToReal(1);
    
    CCTK_REAL_VEC theta = ToReal(1);
    
    CCTK_REAL_VEC beta1rhsL;
    CCTK_REAL_VEC beta2rhsL;
    CCTK_REAL_VEC beta3rhsL;
    
    if (harmonicShift)
    {
      beta1rhsL = 
        kmul(alphaL,kmul(phiL,kmul(ToReal(0.5),kmadd(phiL,kmul(kmadd(DPDstandardNthalpha1L,gtu11,kmadd(DPDstandardNthalpha2L,gtu12,kmul(DPDstandardNthalpha3L,gtu13))),ToReal(-2)),kmul(alphaL,kmadd(kmadd(DPDstandardNthphi2L,gtu12,kmadd(DPDstandardNthphi3L,gtu13,kmadd(gtu11,kmadd(phiL,kmadd(DPDstandardNthgt122L,gtu22,kmadd(DPDstandardNthgt123L,gtu23,kmul(DPDstandardNthgt133L,gtu33))),DPDstandardNthphi1L),kmul(DPDstandardNthgt221L,kmul(phiL,SQR(gtu12)))))),ToReal(2),kmul(phiL,kmadd(DPDstandardNthgt111L,SQR(gtu11),kmadd(DPDstandardNthgt331L,kmul(SQR(gtu13),ToReal(2)),kmadd(gtu12,kmul(gtu33,kmsub(DPDstandardNthgt233L,ToReal(2),DPDstandardNthgt332L)),kmadd(gtu23,kmadd(DPDstandardNthgt231L,kmul(gtu11,ToReal(-2)),kmul(DPDstandardNthgt223L,kmul(gtu12,ToReal(2)))),kmadd(gtu22,kmadd(DPDstandardNthgt222L,gtu12,kmul(DPDstandardNthgt232L,kmul(gtu13,ToReal(2)))),kmadd(gtu11,kmadd(DPDstandardNthgt113L,gtu13,knmsub(DPDstandardNthgt221L,gtu22,knmsub(DPDstandardNthgt331L,gtu33,kmadd(DPDstandardNthgt132L,kmul(gtu23,ToReal(2)),kmul(gtu12,kmadd(DPDstandardNthgt121L,ToReal(2),DPDstandardNthgt112L)))))),kmul(gtu13,kmadd(DPDstandardNthgt333L,gtu33,knmsub(DPDstandardNthgt223L,gtu22,kmadd(kmadd(DPDstandardNthgt131L,gtu11,kmul(DPDstandardNthgt332L,gtu23)),ToReal(2),kmul(DPDstandardNthgt231L,kmul(gtu12,ToReal(4))))))))))))))))))));
      
      beta2rhsL = 
        kmul(alphaL,kmul(phiL,kmul(ToReal(0.5),kmadd(phiL,kmul(kmadd(DPDstandardNthalpha1L,gtu12,kmadd(DPDstandardNthalpha2L,gtu22,kmul(DPDstandardNthalpha3L,gtu23))),ToReal(-2)),kmul(alphaL,kmadd(kmadd(kmadd(phiL,kmadd(DPDstandardNthgt131L,gtu11,kmadd(DPDstandardNthgt331L,gtu13,kmul(DPDstandardNthgt232L,gtu22))),DPDstandardNthphi3L),gtu23,kmul(phiL,kmadd(DPDstandardNthgt112L,SQR(gtu12),kmul(DPDstandardNthgt332L,SQR(gtu23))))),ToReal(2),kmadd(gtu22,kmadd(phiL,kmadd(DPDstandardNthgt221L,gtu12,kmsub(DPDstandardNthgt223L,gtu23,kmul(DPDstandardNthgt112L,gtu11))),kmul(kmadd(phiL,kmadd(DPDstandardNthgt121L,gtu11,kmadd(DPDstandardNthgt123L,gtu13,kmul(DPDstandardNthgt233L,gtu33))),DPDstandardNthphi2L),ToReal(2))),kmadd(phiL,kmadd(DPDstandardNthgt222L,SQR(gtu22),kmadd(gtu23,kmsub(DPDstandardNthgt333L,gtu33,kmul(DPDstandardNthgt113L,gtu11)),kmsub(gtu13,kmadd(DPDstandardNthgt132L,kmul(gtu22,ToReal(-2)),kmul(DPDstandardNthgt231L,kmul(gtu22,ToReal(2)))),kmul(gtu33,kmadd(DPDstandardNthgt332L,gtu22,kmul(DPDstandardNthgt331L,gtu12)))))),kmul(gtu12,kmadd(kmadd(phiL,kmadd(DPDstandardNthgt113L,gtu13,kmul(DPDstandardNthgt133L,gtu33)),DPDstandardNthphi1L),ToReal(2),kmul(phiL,kmadd(DPDstandardNthgt111L,gtu11,kmadd(DPDstandardNthgt122L,kmul(gtu22,ToReal(2)),kmul(DPDstandardNthgt132L,kmul(gtu23,ToReal(4))))))))))))))));
      
      beta3rhsL = 
        kmul(alphaL,kmul(phiL,kmul(ToReal(0.5),kmadd(phiL,kmul(kmadd(DPDstandardNthalpha1L,gtu13,kmadd(DPDstandardNthalpha2L,gtu23,kmul(DPDstandardNthalpha3L,gtu33))),ToReal(-2)),kmul(alphaL,kmadd(gtu33,kmadd(phiL,kmadd(DPDstandardNthgt331L,gtu13,kmadd(DPDstandardNthgt332L,gtu23,kmsub(DPDstandardNthgt123L,kmul(gtu12,ToReal(-2)),kmadd(DPDstandardNthgt223L,gtu22,kmul(DPDstandardNthgt113L,gtu11))))),kmul(kmadd(phiL,kmadd(DPDstandardNthgt131L,gtu11,kmul(DPDstandardNthgt232L,gtu22)),DPDstandardNthphi3L),ToReal(2))),kmadd(phiL,kmadd(DPDstandardNthgt333L,SQR(gtu33),kmul(kmadd(kadd(DPDstandardNthgt132L,DPDstandardNthgt231L),kmul(gtu12,gtu33),kmadd(DPDstandardNthgt113L,SQR(gtu13),kmul(DPDstandardNthgt223L,SQR(gtu23)))),ToReal(2))),kmadd(gtu23,kmadd(kmadd(phiL,kmadd(DPDstandardNthgt121L,gtu11,kmul(DPDstandardNthgt233L,gtu33)),DPDstandardNthphi2L),ToReal(2),kmul(phiL,kmadd(DPDstandardNthgt222L,gtu22,kmsub(DPDstandardNthgt221L,kmul(gtu12,ToReal(2)),kmul(DPDstandardNthgt112L,gtu11))))),kmul(gtu13,kmadd(kmadd(phiL,kmadd(DPDstandardNthgt112L,gtu12,kmul(DPDstandardNthgt133L,gtu33)),DPDstandardNthphi1L),ToReal(2),kmul(phiL,kmadd(DPDstandardNthgt111L,gtu11,kmadd(gtu22,kmsub(DPDstandardNthgt122L,ToReal(2),DPDstandardNthgt221L),kmul(DPDstandardNthgt123L,kmul(gtu23,ToReal(4))))))))))))))));
    }
    else
    {
      beta1rhsL = 
        kmul(theta,kmul(kadd(Xt1L,kmadd(ksub(B1L,Xt1L),ToReal(ShiftBCoeff),kmul(beta1L,kmul(eta,kmul(ToReal(BetaDriver),kadd(ToReal(-1),ToReal(ShiftBCoeff))))))),ToReal(ShiftGammaCoeff)));
      
      beta2rhsL = 
        kmul(theta,kmul(kadd(Xt2L,kmadd(ksub(B2L,Xt2L),ToReal(ShiftBCoeff),kmul(beta2L,kmul(eta,kmul(ToReal(BetaDriver),kadd(ToReal(-1),ToReal(ShiftBCoeff))))))),ToReal(ShiftGammaCoeff)));
      
      beta3rhsL = 
        kmul(theta,kmul(kadd(Xt3L,kmadd(ksub(B3L,Xt3L),ToReal(ShiftBCoeff),kmul(beta3L,kmul(eta,kmul(ToReal(BetaDriver),kadd(ToReal(-1),ToReal(ShiftBCoeff))))))),ToReal(ShiftGammaCoeff)));
    }
    
    CCTK_REAL_VEC B1rhsL = 
      kmul(knmsub(B1L,kmul(eta,ToReal(BetaDriver)),dotXt1),ToReal(ShiftBCoeff));
    
    CCTK_REAL_VEC B2rhsL = 
      kmul(knmsub(B2L,kmul(eta,ToReal(BetaDriver)),dotXt2),ToReal(ShiftBCoeff));
    
    CCTK_REAL_VEC B3rhsL = 
      kmul(knmsub(B3L,kmul(eta,ToReal(BetaDriver)),dotXt3),ToReal(ShiftBCoeff));
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,lc_imin,lc_imax);
    vec_store_nta_partial(alpharhs[index],alpharhsL);
    vec_store_nta_partial(Arhs[index],ArhsL);
    vec_store_nta_partial(B1rhs[index],B1rhsL);
    vec_store_nta_partial(B2rhs[index],B2rhsL);
    vec_store_nta_partial(B3rhs[index],B3rhsL);
    vec_store_nta_partial(beta1rhs[index],beta1rhsL);
    vec_store_nta_partial(beta2rhs[index],beta2rhsL);
    vec_store_nta_partial(beta3rhs[index],beta3rhsL);
    vec_store_nta_partial(gt11rhs[index],gt11rhsL);
    vec_store_nta_partial(gt12rhs[index],gt12rhsL);
    vec_store_nta_partial(gt13rhs[index],gt13rhsL);
    vec_store_nta_partial(gt22rhs[index],gt22rhsL);
    vec_store_nta_partial(gt23rhs[index],gt23rhsL);
    vec_store_nta_partial(gt33rhs[index],gt33rhsL);
    vec_store_nta_partial(phirhs[index],phirhsL);
    vec_store_nta_partial(trKrhs[index],trKrhsL);
    vec_store_nta_partial(Xt1rhs[index],Xt1rhsL);
    vec_store_nta_partial(Xt2rhs[index],Xt2rhsL);
    vec_store_nta_partial(Xt3rhs[index],Xt3rhsL);
  }
  LC_ENDLOOP3VEC(ML_BSSN_Host_RHS1);
}

extern "C" void ML_BSSN_Host_RHS1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_Host_RHS1_Body");
  }
  
  if (cctk_iteration % ML_BSSN_Host_RHS1_calc_every != ML_BSSN_Host_RHS1_calc_offset)
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
    "ML_BSSN_Host::DPDstandardNthbeta111_group",
    "ML_BSSN_Host::DPDstandardNthbeta112_group",
    "ML_BSSN_Host::DPDstandardNthbeta113_group",
    "ML_BSSN_Host::DPDstandardNthbeta11_group",
    "ML_BSSN_Host::DPDstandardNthbeta122_group",
    "ML_BSSN_Host::DPDstandardNthbeta123_group",
    "ML_BSSN_Host::DPDstandardNthbeta12_group",
    "ML_BSSN_Host::DPDstandardNthbeta133_group",
    "ML_BSSN_Host::DPDstandardNthbeta13_group",
    "ML_BSSN_Host::DPDstandardNthbeta211_group",
    "ML_BSSN_Host::DPDstandardNthbeta212_group",
    "ML_BSSN_Host::DPDstandardNthbeta213_group",
    "ML_BSSN_Host::DPDstandardNthbeta21_group",
    "ML_BSSN_Host::DPDstandardNthbeta222_group",
    "ML_BSSN_Host::DPDstandardNthbeta223_group",
    "ML_BSSN_Host::DPDstandardNthbeta22_group",
    "ML_BSSN_Host::DPDstandardNthbeta233_group",
    "ML_BSSN_Host::DPDstandardNthbeta23_group",
    "ML_BSSN_Host::DPDstandardNthbeta311_group",
    "ML_BSSN_Host::DPDstandardNthbeta312_group",
    "ML_BSSN_Host::DPDstandardNthbeta313_group",
    "ML_BSSN_Host::DPDstandardNthbeta31_group",
    "ML_BSSN_Host::DPDstandardNthbeta322_group",
    "ML_BSSN_Host::DPDstandardNthbeta323_group",
    "ML_BSSN_Host::DPDstandardNthbeta32_group",
    "ML_BSSN_Host::DPDstandardNthbeta333_group",
    "ML_BSSN_Host::DPDstandardNthbeta33_group",
    "ML_BSSN_Host::DPDstandardNthgt111_group",
    "ML_BSSN_Host::DPDstandardNthgt112_group",
    "ML_BSSN_Host::DPDstandardNthgt113_group",
    "ML_BSSN_Host::DPDstandardNthgt121_group",
    "ML_BSSN_Host::DPDstandardNthgt122_group",
    "ML_BSSN_Host::DPDstandardNthgt123_group",
    "ML_BSSN_Host::DPDstandardNthgt131_group",
    "ML_BSSN_Host::DPDstandardNthgt132_group",
    "ML_BSSN_Host::DPDstandardNthgt133_group",
    "ML_BSSN_Host::DPDstandardNthgt221_group",
    "ML_BSSN_Host::DPDstandardNthgt222_group",
    "ML_BSSN_Host::DPDstandardNthgt223_group",
    "ML_BSSN_Host::DPDstandardNthgt231_group",
    "ML_BSSN_Host::DPDstandardNthgt232_group",
    "ML_BSSN_Host::DPDstandardNthgt233_group",
    "ML_BSSN_Host::DPDstandardNthgt331_group",
    "ML_BSSN_Host::DPDstandardNthgt332_group",
    "ML_BSSN_Host::DPDstandardNthgt333_group",
    "ML_BSSN_Host::DPDstandardNthphi1_group",
    "ML_BSSN_Host::DPDstandardNthphi2_group",
    "ML_BSSN_Host::DPDstandardNthphi3_group",
    "ML_BSSN_Host::DPDstandardNthtrK1_group",
    "ML_BSSN_Host::DPDstandardNthtrK2_group",
    "ML_BSSN_Host::DPDstandardNthtrK3_group",
    "ML_BSSN_Host::ML_curv",
    "ML_BSSN_Host::ML_dtlapse",
    "ML_BSSN_Host::ML_dtlapserhs",
    "ML_BSSN_Host::ML_dtshift",
    "ML_BSSN_Host::ML_dtshiftrhs",
    "ML_BSSN_Host::ML_Gamma",
    "ML_BSSN_Host::ML_Gammarhs",
    "ML_BSSN_Host::ML_lapse",
    "ML_BSSN_Host::ML_lapserhs",
    "ML_BSSN_Host::ML_log_confac",
    "ML_BSSN_Host::ML_log_confacrhs",
    "ML_BSSN_Host::ML_metric",
    "ML_BSSN_Host::ML_metricrhs",
    "ML_BSSN_Host::ML_shift",
    "ML_BSSN_Host::ML_shiftrhs",
    "ML_BSSN_Host::ML_trace_curv",
    "ML_BSSN_Host::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_RHS1", 77, groups);
  
  
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_Host_RHS1_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_RHS1_Body");
  }
}
