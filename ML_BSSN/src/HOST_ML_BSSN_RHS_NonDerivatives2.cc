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

/* Define macros used in calculations */
#define INITVALUE (42)
#define INV(x) ((CCTK_REAL)1.0 / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * SQR(x))
#define QAD(x) (SQR(SQR(x)))

extern "C" void HOST_ML_BSSN_RHS_NonDerivatives2_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % HOST_ML_BSSN_RHS_NonDerivatives2_calc_every != HOST_ML_BSSN_RHS_NonDerivatives2_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_curvrhs.");
  return;
}

static void HOST_ML_BSSN_RHS_NonDerivatives2_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(0));
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(1));
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(2));
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_TIME);
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = INV(dx);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = INV(dy);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = INV(dz);
  const CCTK_REAL khalf CCTK_ATTRIBUTE_UNUSED = 0.5;
  const CCTK_REAL kthird CCTK_ATTRIBUTE_UNUSED = 
    0.333333333333333333333333333333;
  const CCTK_REAL ktwothird CCTK_ATTRIBUTE_UNUSED = 
    0.666666666666666666666666666667;
  const CCTK_REAL kfourthird CCTK_ATTRIBUTE_UNUSED = 
    1.33333333333333333333333333333;
  const CCTK_REAL hdxi CCTK_ATTRIBUTE_UNUSED = 0.5*dxi;
  const CCTK_REAL hdyi CCTK_ATTRIBUTE_UNUSED = 0.5*dyi;
  const CCTK_REAL hdzi CCTK_ATTRIBUTE_UNUSED = 0.5*dzi;
  
  /* Initialize predefined quantities */
  const CCTK_REAL p1o1024dx CCTK_ATTRIBUTE_UNUSED = 0.0009765625*INV(dx);
  const CCTK_REAL p1o1024dy CCTK_ATTRIBUTE_UNUSED = 0.0009765625*INV(dy);
  const CCTK_REAL p1o1024dz CCTK_ATTRIBUTE_UNUSED = 0.0009765625*INV(dz);
  const CCTK_REAL p1o1680dx CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*INV(dx);
  const CCTK_REAL p1o1680dy CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*INV(dy);
  const CCTK_REAL p1o1680dz CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*INV(dz);
  const CCTK_REAL p1o2dx CCTK_ATTRIBUTE_UNUSED = 0.5*INV(dx);
  const CCTK_REAL p1o2dy CCTK_ATTRIBUTE_UNUSED = 0.5*INV(dy);
  const CCTK_REAL p1o2dz CCTK_ATTRIBUTE_UNUSED = 0.5*INV(dz);
  const CCTK_REAL p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*INV(SQR(dx));
  const CCTK_REAL p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*INV(SQR(dy));
  const CCTK_REAL p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*INV(SQR(dz));
  const CCTK_REAL p1o560dx CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*INV(dx);
  const CCTK_REAL p1o560dy CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*INV(dy);
  const CCTK_REAL p1o560dz CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*INV(dz);
  const CCTK_REAL p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*INV(dx*dy);
  const CCTK_REAL p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*INV(dx*dz);
  const CCTK_REAL p1o705600dydz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*INV(dy*dz);
  const CCTK_REAL p1o840dx CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*INV(dx);
  const CCTK_REAL p1o840dy CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*INV(dy);
  const CCTK_REAL p1o840dz CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*INV(dz);
  const CCTK_REAL p1odx CCTK_ATTRIBUTE_UNUSED = INV(dx);
  const CCTK_REAL p1ody CCTK_ATTRIBUTE_UNUSED = INV(dy);
  const CCTK_REAL p1odz CCTK_ATTRIBUTE_UNUSED = INV(dz);
  const CCTK_REAL pm1o2dx CCTK_ATTRIBUTE_UNUSED = -0.5*INV(dx);
  const CCTK_REAL pm1o2dy CCTK_ATTRIBUTE_UNUSED = -0.5*INV(dy);
  const CCTK_REAL pm1o2dz CCTK_ATTRIBUTE_UNUSED = -0.5*INV(dz);
  const CCTK_REAL pm1o840dx CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*INV(dx);
  const CCTK_REAL pm1o840dy CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*INV(dy);
  const CCTK_REAL pm1o840dz CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*INV(dz);
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  #pragma omp parallel // reduction(+: vec_iter_counter, vec_op_counter, vec_mem_counter)
  CCTK_LOOP3(HOST_ML_BSSN_RHS_NonDerivatives2,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    // ++vec_iter_counter;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpha[index];
    CCTK_REAL At11L CCTK_ATTRIBUTE_UNUSED = At11[index];
    CCTK_REAL At12L CCTK_ATTRIBUTE_UNUSED = At12[index];
    CCTK_REAL At13L CCTK_ATTRIBUTE_UNUSED = At13[index];
    CCTK_REAL At22L CCTK_ATTRIBUTE_UNUSED = At22[index];
    CCTK_REAL At23L CCTK_ATTRIBUTE_UNUSED = At23[index];
    CCTK_REAL At33L CCTK_ATTRIBUTE_UNUSED = At33[index];
    CCTK_REAL beta1L CCTK_ATTRIBUTE_UNUSED = beta1[index];
    CCTK_REAL beta2L CCTK_ATTRIBUTE_UNUSED = beta2[index];
    CCTK_REAL beta3L CCTK_ATTRIBUTE_UNUSED = beta3[index];
    CCTK_REAL DPDstandardNthalpha1L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthalpha1[index];
    CCTK_REAL DPDstandardNthalpha11L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthalpha11[index];
    CCTK_REAL DPDstandardNthalpha12L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthalpha12[index];
    CCTK_REAL DPDstandardNthalpha13L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthalpha13[index];
    CCTK_REAL DPDstandardNthalpha2L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthalpha2[index];
    CCTK_REAL DPDstandardNthalpha22L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthalpha22[index];
    CCTK_REAL DPDstandardNthalpha23L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthalpha23[index];
    CCTK_REAL DPDstandardNthalpha3L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthalpha3[index];
    CCTK_REAL DPDstandardNthalpha33L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthalpha33[index];
    CCTK_REAL DPDstandardNthbeta11L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthbeta11[index];
    CCTK_REAL DPDstandardNthbeta12L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthbeta12[index];
    CCTK_REAL DPDstandardNthbeta13L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthbeta13[index];
    CCTK_REAL DPDstandardNthbeta21L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthbeta21[index];
    CCTK_REAL DPDstandardNthbeta22L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthbeta22[index];
    CCTK_REAL DPDstandardNthbeta23L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthbeta23[index];
    CCTK_REAL DPDstandardNthbeta31L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthbeta31[index];
    CCTK_REAL DPDstandardNthbeta32L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthbeta32[index];
    CCTK_REAL DPDstandardNthbeta33L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthbeta33[index];
    CCTK_REAL DPDstandardNthgt111L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt111[index];
    CCTK_REAL DPDstandardNthgt1111L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1111[index];
    CCTK_REAL DPDstandardNthgt1112L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1112[index];
    CCTK_REAL DPDstandardNthgt1113L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1113[index];
    CCTK_REAL DPDstandardNthgt112L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt112[index];
    CCTK_REAL DPDstandardNthgt1122L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1122[index];
    CCTK_REAL DPDstandardNthgt1123L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1123[index];
    CCTK_REAL DPDstandardNthgt113L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt113[index];
    CCTK_REAL DPDstandardNthgt1133L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1133[index];
    CCTK_REAL DPDstandardNthgt121L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt121[index];
    CCTK_REAL DPDstandardNthgt1211L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1211[index];
    CCTK_REAL DPDstandardNthgt1212L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1212[index];
    CCTK_REAL DPDstandardNthgt1213L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1213[index];
    CCTK_REAL DPDstandardNthgt122L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt122[index];
    CCTK_REAL DPDstandardNthgt1222L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1222[index];
    CCTK_REAL DPDstandardNthgt1223L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1223[index];
    CCTK_REAL DPDstandardNthgt123L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt123[index];
    CCTK_REAL DPDstandardNthgt1233L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1233[index];
    CCTK_REAL DPDstandardNthgt131L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt131[index];
    CCTK_REAL DPDstandardNthgt1311L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1311[index];
    CCTK_REAL DPDstandardNthgt1312L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1312[index];
    CCTK_REAL DPDstandardNthgt1313L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1313[index];
    CCTK_REAL DPDstandardNthgt132L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt132[index];
    CCTK_REAL DPDstandardNthgt1322L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1322[index];
    CCTK_REAL DPDstandardNthgt1323L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1323[index];
    CCTK_REAL DPDstandardNthgt133L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt133[index];
    CCTK_REAL DPDstandardNthgt1333L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt1333[index];
    CCTK_REAL DPDstandardNthgt221L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt221[index];
    CCTK_REAL DPDstandardNthgt2211L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2211[index];
    CCTK_REAL DPDstandardNthgt2212L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2212[index];
    CCTK_REAL DPDstandardNthgt2213L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2213[index];
    CCTK_REAL DPDstandardNthgt222L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt222[index];
    CCTK_REAL DPDstandardNthgt2222L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2222[index];
    CCTK_REAL DPDstandardNthgt2223L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2223[index];
    CCTK_REAL DPDstandardNthgt223L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt223[index];
    CCTK_REAL DPDstandardNthgt2233L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2233[index];
    CCTK_REAL DPDstandardNthgt231L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt231[index];
    CCTK_REAL DPDstandardNthgt2311L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2311[index];
    CCTK_REAL DPDstandardNthgt2312L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2312[index];
    CCTK_REAL DPDstandardNthgt2313L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2313[index];
    CCTK_REAL DPDstandardNthgt232L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt232[index];
    CCTK_REAL DPDstandardNthgt2322L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2322[index];
    CCTK_REAL DPDstandardNthgt2323L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2323[index];
    CCTK_REAL DPDstandardNthgt233L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt233[index];
    CCTK_REAL DPDstandardNthgt2333L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt2333[index];
    CCTK_REAL DPDstandardNthgt331L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt331[index];
    CCTK_REAL DPDstandardNthgt3311L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt3311[index];
    CCTK_REAL DPDstandardNthgt3312L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt3312[index];
    CCTK_REAL DPDstandardNthgt3313L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt3313[index];
    CCTK_REAL DPDstandardNthgt332L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt332[index];
    CCTK_REAL DPDstandardNthgt3322L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt3322[index];
    CCTK_REAL DPDstandardNthgt3323L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt3323[index];
    CCTK_REAL DPDstandardNthgt333L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt333[index];
    CCTK_REAL DPDstandardNthgt3333L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt3333[index];
    CCTK_REAL DPDstandardNthphi1L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi1[index];
    CCTK_REAL DPDstandardNthphi11L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi11[index];
    CCTK_REAL DPDstandardNthphi12L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi12[index];
    CCTK_REAL DPDstandardNthphi13L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi13[index];
    CCTK_REAL DPDstandardNthphi2L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi2[index];
    CCTK_REAL DPDstandardNthphi22L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi22[index];
    CCTK_REAL DPDstandardNthphi23L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi23[index];
    CCTK_REAL DPDstandardNthphi3L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi3[index];
    CCTK_REAL DPDstandardNthphi33L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi33[index];
    CCTK_REAL DPDstandardNthXt11L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt11[index];
    CCTK_REAL DPDstandardNthXt12L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt12[index];
    CCTK_REAL DPDstandardNthXt13L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt13[index];
    CCTK_REAL DPDstandardNthXt21L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt21[index];
    CCTK_REAL DPDstandardNthXt22L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt22[index];
    CCTK_REAL DPDstandardNthXt23L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt23[index];
    CCTK_REAL DPDstandardNthXt31L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt31[index];
    CCTK_REAL DPDstandardNthXt32L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt32[index];
    CCTK_REAL DPDstandardNthXt33L CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt33[index];
    CCTK_REAL gt11L CCTK_ATTRIBUTE_UNUSED = gt11[index];
    CCTK_REAL gt12L CCTK_ATTRIBUTE_UNUSED = gt12[index];
    CCTK_REAL gt13L CCTK_ATTRIBUTE_UNUSED = gt13[index];
    CCTK_REAL gt22L CCTK_ATTRIBUTE_UNUSED = gt22[index];
    CCTK_REAL gt23L CCTK_ATTRIBUTE_UNUSED = gt23[index];
    CCTK_REAL gt33L CCTK_ATTRIBUTE_UNUSED = gt33[index];
    CCTK_REAL phiL CCTK_ATTRIBUTE_UNUSED = phi[index];
    CCTK_REAL trKL CCTK_ATTRIBUTE_UNUSED = trK[index];
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = isgn(beta1L);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = isgn(beta2L);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = isgn(beta3L);
    
    CCTK_REAL detgt CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL gtu11 CCTK_ATTRIBUTE_UNUSED = INV(detgt)*(gt22L*gt33L - 
      SQR(gt23L));
    
    CCTK_REAL gtu12 CCTK_ATTRIBUTE_UNUSED = (gt13L*gt23L - 
      gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu13 CCTK_ATTRIBUTE_UNUSED = (-(gt13L*gt22L) + 
      gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 CCTK_ATTRIBUTE_UNUSED = INV(detgt)*(gt11L*gt33L - 
      SQR(gt13L));
    
    CCTK_REAL gtu23 CCTK_ATTRIBUTE_UNUSED = (gt12L*gt13L - 
      gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 CCTK_ATTRIBUTE_UNUSED = INV(detgt)*(gt11L*gt22L - 
      SQR(gt12L));
    
    CCTK_REAL Gtl111 CCTK_ATTRIBUTE_UNUSED = 0.5*DPDstandardNthgt111L;
    
    CCTK_REAL Gtl112 CCTK_ATTRIBUTE_UNUSED = 0.5*DPDstandardNthgt112L;
    
    CCTK_REAL Gtl113 CCTK_ATTRIBUTE_UNUSED = 0.5*DPDstandardNthgt113L;
    
    CCTK_REAL Gtl122 CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt122L - 
      0.5*DPDstandardNthgt221L;
    
    CCTK_REAL Gtl123 CCTK_ATTRIBUTE_UNUSED = 0.5*(DPDstandardNthgt123L + 
      DPDstandardNthgt132L - DPDstandardNthgt231L);
    
    CCTK_REAL Gtl133 CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt133L - 
      0.5*DPDstandardNthgt331L;
    
    CCTK_REAL Gtl211 CCTK_ATTRIBUTE_UNUSED = -0.5*DPDstandardNthgt112L + 
      DPDstandardNthgt121L;
    
    CCTK_REAL Gtl212 CCTK_ATTRIBUTE_UNUSED = 0.5*DPDstandardNthgt221L;
    
    CCTK_REAL Gtl213 CCTK_ATTRIBUTE_UNUSED = 0.5*(DPDstandardNthgt123L - 
      DPDstandardNthgt132L + DPDstandardNthgt231L);
    
    CCTK_REAL Gtl222 CCTK_ATTRIBUTE_UNUSED = 0.5*DPDstandardNthgt222L;
    
    CCTK_REAL Gtl223 CCTK_ATTRIBUTE_UNUSED = 0.5*DPDstandardNthgt223L;
    
    CCTK_REAL Gtl233 CCTK_ATTRIBUTE_UNUSED = DPDstandardNthgt233L - 
      0.5*DPDstandardNthgt332L;
    
    CCTK_REAL Gtl311 CCTK_ATTRIBUTE_UNUSED = -0.5*DPDstandardNthgt113L + 
      DPDstandardNthgt131L;
    
    CCTK_REAL Gtl312 CCTK_ATTRIBUTE_UNUSED = 0.5*(-DPDstandardNthgt123L + 
      DPDstandardNthgt132L + DPDstandardNthgt231L);
    
    CCTK_REAL Gtl313 CCTK_ATTRIBUTE_UNUSED = 0.5*DPDstandardNthgt331L;
    
    CCTK_REAL Gtl322 CCTK_ATTRIBUTE_UNUSED = -0.5*DPDstandardNthgt223L + 
      DPDstandardNthgt232L;
    
    CCTK_REAL Gtl323 CCTK_ATTRIBUTE_UNUSED = 0.5*DPDstandardNthgt332L;
    
    CCTK_REAL Gtl333 CCTK_ATTRIBUTE_UNUSED = 0.5*DPDstandardNthgt333L;
    
    CCTK_REAL Gtlu111 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu11 + Gtl112*gtu12 
      + Gtl113*gtu13;
    
    CCTK_REAL Gtlu112 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu12 + Gtl112*gtu22 
      + Gtl113*gtu23;
    
    CCTK_REAL Gtlu113 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu13 + Gtl112*gtu23 
      + Gtl113*gtu33;
    
    CCTK_REAL Gtlu121 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu11 + Gtl122*gtu12 
      + Gtl123*gtu13;
    
    CCTK_REAL Gtlu122 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu12 + Gtl122*gtu22 
      + Gtl123*gtu23;
    
    CCTK_REAL Gtlu123 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu13 + Gtl122*gtu23 
      + Gtl123*gtu33;
    
    CCTK_REAL Gtlu131 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu11 + Gtl123*gtu12 
      + Gtl133*gtu13;
    
    CCTK_REAL Gtlu132 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu12 + Gtl123*gtu22 
      + Gtl133*gtu23;
    
    CCTK_REAL Gtlu133 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu13 + Gtl123*gtu23 
      + Gtl133*gtu33;
    
    CCTK_REAL Gtlu211 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu11 + Gtl212*gtu12 
      + Gtl213*gtu13;
    
    CCTK_REAL Gtlu212 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu12 + Gtl212*gtu22 
      + Gtl213*gtu23;
    
    CCTK_REAL Gtlu213 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu13 + Gtl212*gtu23 
      + Gtl213*gtu33;
    
    CCTK_REAL Gtlu221 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu11 + Gtl222*gtu12 
      + Gtl223*gtu13;
    
    CCTK_REAL Gtlu222 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu12 + Gtl222*gtu22 
      + Gtl223*gtu23;
    
    CCTK_REAL Gtlu223 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu13 + Gtl222*gtu23 
      + Gtl223*gtu33;
    
    CCTK_REAL Gtlu231 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu11 + Gtl223*gtu12 
      + Gtl233*gtu13;
    
    CCTK_REAL Gtlu232 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu12 + Gtl223*gtu22 
      + Gtl233*gtu23;
    
    CCTK_REAL Gtlu233 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu13 + Gtl223*gtu23 
      + Gtl233*gtu33;
    
    CCTK_REAL Gtlu311 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu11 + Gtl312*gtu12 
      + Gtl313*gtu13;
    
    CCTK_REAL Gtlu312 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu12 + Gtl312*gtu22 
      + Gtl313*gtu23;
    
    CCTK_REAL Gtlu313 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu13 + Gtl312*gtu23 
      + Gtl313*gtu33;
    
    CCTK_REAL Gtlu321 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu11 + Gtl322*gtu12 
      + Gtl323*gtu13;
    
    CCTK_REAL Gtlu322 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu12 + Gtl322*gtu22 
      + Gtl323*gtu23;
    
    CCTK_REAL Gtlu323 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu13 + Gtl322*gtu23 
      + Gtl323*gtu33;
    
    CCTK_REAL Gtlu331 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu11 + Gtl323*gtu12 
      + Gtl333*gtu13;
    
    CCTK_REAL Gtlu332 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu12 + Gtl323*gtu22 
      + Gtl333*gtu23;
    
    CCTK_REAL Gtlu333 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu13 + Gtl323*gtu23 
      + Gtl333*gtu33;
    
    CCTK_REAL Gt111 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu11 + Gtl211*gtu12 + 
      Gtl311*gtu13;
    
    CCTK_REAL Gt211 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu12 + Gtl211*gtu22 + 
      Gtl311*gtu23;
    
    CCTK_REAL Gt311 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu13 + Gtl211*gtu23 + 
      Gtl311*gtu33;
    
    CCTK_REAL Gt112 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu11 + Gtl212*gtu12 + 
      Gtl312*gtu13;
    
    CCTK_REAL Gt212 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu12 + Gtl212*gtu22 + 
      Gtl312*gtu23;
    
    CCTK_REAL Gt312 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu13 + Gtl212*gtu23 + 
      Gtl312*gtu33;
    
    CCTK_REAL Gt113 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu11 + Gtl213*gtu12 + 
      Gtl313*gtu13;
    
    CCTK_REAL Gt213 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu12 + Gtl213*gtu22 + 
      Gtl313*gtu23;
    
    CCTK_REAL Gt313 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu13 + Gtl213*gtu23 + 
      Gtl313*gtu33;
    
    CCTK_REAL Gt122 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu11 + Gtl222*gtu12 + 
      Gtl322*gtu13;
    
    CCTK_REAL Gt222 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu12 + Gtl222*gtu22 + 
      Gtl322*gtu23;
    
    CCTK_REAL Gt322 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu13 + Gtl222*gtu23 + 
      Gtl322*gtu33;
    
    CCTK_REAL Gt123 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu11 + Gtl223*gtu12 + 
      Gtl323*gtu13;
    
    CCTK_REAL Gt223 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu12 + Gtl223*gtu22 + 
      Gtl323*gtu23;
    
    CCTK_REAL Gt323 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu13 + Gtl223*gtu23 + 
      Gtl323*gtu33;
    
    CCTK_REAL Gt133 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu11 + Gtl233*gtu12 + 
      Gtl333*gtu13;
    
    CCTK_REAL Gt233 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu12 + Gtl233*gtu22 + 
      Gtl333*gtu23;
    
    CCTK_REAL Gt333 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu13 + Gtl233*gtu23 + 
      Gtl333*gtu33;
    
    CCTK_REAL Xtn1 CCTK_ATTRIBUTE_UNUSED = Gt111*gtu11 + Gt122*gtu22 + 
      2*(Gt112*gtu12 + Gt113*gtu13 + Gt123*gtu23) + Gt133*gtu33;
    
    CCTK_REAL Xtn2 CCTK_ATTRIBUTE_UNUSED = Gt211*gtu11 + Gt222*gtu22 + 
      2*(Gt212*gtu12 + Gt213*gtu13 + Gt223*gtu23) + Gt233*gtu33;
    
    CCTK_REAL Xtn3 CCTK_ATTRIBUTE_UNUSED = Gt311*gtu11 + Gt322*gtu22 + 
      2*(Gt312*gtu12 + Gt313*gtu13 + Gt323*gtu23) + Gt333*gtu33;
    
    CCTK_REAL Rt11 CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt11L*gt11L + 
      DPDstandardNthXt21L*gt12L + DPDstandardNthXt31L*gt13L + 
      3*(Gt111*Gtlu111 + Gt112*Gtlu112 + Gt113*Gtlu113) + 2*(Gt211*Gtlu121 + 
      Gt212*Gtlu122 + Gt213*Gtlu123 + Gt311*Gtlu131 + Gt312*Gtlu132 + 
      Gt313*Gtlu133) + Gt211*Gtlu211 + Gt212*Gtlu212 + Gt213*Gtlu213 + 
      Gt311*Gtlu311 + Gt312*Gtlu312 + Gt313*Gtlu313 + 
      0.5*(-(DPDstandardNthgt1111L*gtu11) - 2*DPDstandardNthgt1112L*gtu12 - 
      2*DPDstandardNthgt1113L*gtu13 - DPDstandardNthgt1122L*gtu22 - 
      2*DPDstandardNthgt1123L*gtu23 - DPDstandardNthgt1133L*gtu33) + 
      Gtl111*Xtn1 + Gtl112*Xtn2 + Gtl113*Xtn3;
    
    CCTK_REAL Rt12 CCTK_ATTRIBUTE_UNUSED = 0.5*(DPDstandardNthXt12L*gt11L 
      + (DPDstandardNthXt11L + DPDstandardNthXt22L)*gt12L + 
      DPDstandardNthXt32L*gt13L + DPDstandardNthXt21L*gt22L + 
      DPDstandardNthXt31L*gt23L + 4*(Gt211*Gtlu221 + Gt212*Gtlu222 + 
      Gt213*Gtlu223) + 2*(Gt112*Gtlu111 + Gt122*Gtlu112 + Gt123*Gtlu113 + 
      Gt111*Gtlu121 + Gt212*Gtlu121 + Gt112*Gtlu122 + Gt222*Gtlu122 + 
      Gt113*Gtlu123 + Gt223*Gtlu123 + Gt312*Gtlu131 + Gt322*Gtlu132 + 
      Gt323*Gtlu133 + Gt111*Gtlu211 + Gt112*Gtlu212 + Gt113*Gtlu213 + 
      Gt311*Gtlu231 + Gt312*Gtlu232 + Gt313*Gtlu233 + Gt311*Gtlu321 + 
      Gt312*Gtlu322 + Gt313*Gtlu323) - DPDstandardNthgt1211L*gtu11 - 
      2*DPDstandardNthgt1212L*gtu12 - 2*DPDstandardNthgt1213L*gtu13 - 
      DPDstandardNthgt1222L*gtu22 - 2*DPDstandardNthgt1223L*gtu23 - 
      DPDstandardNthgt1233L*gtu33 + Gtl112*Xtn1 + Gtl211*Xtn1 + Gtl122*Xtn2 + 
      Gtl212*Xtn2 + Gtl123*Xtn3 + Gtl213*Xtn3);
    
    CCTK_REAL Rt13 CCTK_ATTRIBUTE_UNUSED = 0.5*(DPDstandardNthXt13L*gt11L 
      + DPDstandardNthXt23L*gt12L + (DPDstandardNthXt11L + 
      DPDstandardNthXt33L)*gt13L + DPDstandardNthXt21L*gt23L + 
      DPDstandardNthXt31L*gt33L + 2*(Gt113*Gtlu111 + Gt123*Gtlu112 + 
      Gt133*Gtlu113 + Gt213*Gtlu121 + Gt223*Gtlu122 + Gt233*Gtlu123 + 
      Gt111*Gtlu131 + Gt313*Gtlu131 + Gt112*Gtlu132 + Gt323*Gtlu132 + 
      Gt113*Gtlu133 + Gt333*Gtlu133 + Gt211*Gtlu231 + Gt212*Gtlu232 + 
      Gt213*Gtlu233 + Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*Gtlu313 + 
      Gt211*Gtlu321 + Gt212*Gtlu322 + Gt213*Gtlu323) + 4*(Gt311*Gtlu331 + 
      Gt312*Gtlu332 + Gt313*Gtlu333) - DPDstandardNthgt1311L*gtu11 - 
      2*DPDstandardNthgt1312L*gtu12 - 2*DPDstandardNthgt1313L*gtu13 - 
      DPDstandardNthgt1322L*gtu22 - 2*DPDstandardNthgt1323L*gtu23 - 
      DPDstandardNthgt1333L*gtu33 + Gtl113*Xtn1 + Gtl311*Xtn1 + Gtl123*Xtn2 + 
      Gtl312*Xtn2 + Gtl133*Xtn3 + Gtl313*Xtn3);
    
    CCTK_REAL Rt22 CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt12L*gt12L + 
      DPDstandardNthXt22L*gt22L + DPDstandardNthXt32L*gt23L + Gt112*(Gtlu121 
      + 2*Gtlu211) + Gt122*(Gtlu122 + 2*Gtlu212) + Gt123*(Gtlu123 + 
      2*Gtlu213) + 3*(Gt212*Gtlu221 + Gt222*Gtlu222 + Gt223*Gtlu223) + 
      2*(Gt312*Gtlu231 + Gt322*Gtlu232 + Gt323*Gtlu233) + Gt312*Gtlu321 + 
      Gt322*Gtlu322 + Gt323*Gtlu323 + 0.5*(-(DPDstandardNthgt2211L*gtu11) - 
      2*DPDstandardNthgt2212L*gtu12 - 2*DPDstandardNthgt2213L*gtu13 - 
      DPDstandardNthgt2222L*gtu22 - 2*DPDstandardNthgt2223L*gtu23 - 
      DPDstandardNthgt2233L*gtu33) + Gtl212*Xtn1 + Gtl222*Xtn2 + Gtl223*Xtn3;
    
    CCTK_REAL Rt23 CCTK_ATTRIBUTE_UNUSED = 0.5*(DPDstandardNthXt13L*gt12L 
      + DPDstandardNthXt12L*gt13L + DPDstandardNthXt23L*gt22L + 
      (DPDstandardNthXt22L + DPDstandardNthXt33L)*gt23L + 
      DPDstandardNthXt32L*gt33L + 2*(Gt112*Gtlu131 + Gt122*Gtlu132 + 
      Gt123*Gtlu133 + Gt113*Gtlu211 + Gt123*Gtlu212 + Gt133*Gtlu213 + 
      Gt213*Gtlu221 + Gt223*Gtlu222 + Gt233*Gtlu223 + Gt212*Gtlu231 + 
      Gt313*Gtlu231 + Gt222*Gtlu232 + Gt323*Gtlu232 + Gt223*Gtlu233 + 
      Gt333*Gtlu233 + Gt112*Gtlu311 + Gt122*Gtlu312 + Gt123*Gtlu313 + 
      Gt212*Gtlu321 + Gt222*Gtlu322 + Gt223*Gtlu323) + 4*(Gt312*Gtlu331 + 
      Gt322*Gtlu332 + Gt323*Gtlu333) - DPDstandardNthgt2311L*gtu11 - 
      2*DPDstandardNthgt2312L*gtu12 - 2*DPDstandardNthgt2313L*gtu13 - 
      DPDstandardNthgt2322L*gtu22 - 2*DPDstandardNthgt2323L*gtu23 - 
      DPDstandardNthgt2333L*gtu33 + Gtl213*Xtn1 + Gtl312*Xtn1 + Gtl223*Xtn2 + 
      Gtl322*Xtn2 + Gtl233*Xtn3 + Gtl323*Xtn3);
    
    CCTK_REAL Rt33 CCTK_ATTRIBUTE_UNUSED = DPDstandardNthXt13L*gt13L + 
      DPDstandardNthXt23L*gt23L + DPDstandardNthXt33L*gt33L + Gt113*(Gtlu131 
      + 2*Gtlu311) + Gt123*(Gtlu132 + 2*Gtlu312) + Gt133*(Gtlu133 + 
      2*Gtlu313) + Gt213*(Gtlu231 + 2*Gtlu321) + Gt223*(Gtlu232 + 2*Gtlu322) 
      + Gt233*(Gtlu233 + 2*Gtlu323) + 3*(Gt313*Gtlu331 + Gt323*Gtlu332 + 
      Gt333*Gtlu333) + 0.5*(-(DPDstandardNthgt3311L*gtu11) - 
      2*DPDstandardNthgt3312L*gtu12 - 2*DPDstandardNthgt3313L*gtu13 - 
      DPDstandardNthgt3322L*gtu22 - 2*DPDstandardNthgt3323L*gtu23 - 
      DPDstandardNthgt3333L*gtu33) + Gtl313*Xtn1 + Gtl323*Xtn2 + Gtl333*Xtn3;
    
    CCTK_REAL fac1 CCTK_ATTRIBUTE_UNUSED = 
      IfThen(conformalMethod,-0.5*INV(phiL),1);
    
    CCTK_REAL cdphi1 CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi1L*fac1;
    
    CCTK_REAL cdphi2 CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi2L*fac1;
    
    CCTK_REAL cdphi3 CCTK_ATTRIBUTE_UNUSED = DPDstandardNthphi3L*fac1;
    
    CCTK_REAL fac2 CCTK_ATTRIBUTE_UNUSED = 
      IfThen(conformalMethod,0.5*INV(SQR(phiL)),0);
    
    CCTK_REAL cdphi211 CCTK_ATTRIBUTE_UNUSED = fac1*(DPDstandardNthphi11L 
      - DPDstandardNthphi1L*Gt111 - DPDstandardNthphi2L*Gt211 - 
      DPDstandardNthphi3L*Gt311) + fac2*SQR(DPDstandardNthphi1L);
    
    CCTK_REAL cdphi212 CCTK_ATTRIBUTE_UNUSED = 
      DPDstandardNthphi1L*DPDstandardNthphi2L*fac2 + 
      fac1*(DPDstandardNthphi12L - DPDstandardNthphi1L*Gt112 - 
      DPDstandardNthphi2L*Gt212 - DPDstandardNthphi3L*Gt312);
    
    CCTK_REAL cdphi213 CCTK_ATTRIBUTE_UNUSED = 
      DPDstandardNthphi1L*DPDstandardNthphi3L*fac2 + 
      fac1*(DPDstandardNthphi13L - DPDstandardNthphi1L*Gt113 - 
      DPDstandardNthphi2L*Gt213 - DPDstandardNthphi3L*Gt313);
    
    CCTK_REAL cdphi222 CCTK_ATTRIBUTE_UNUSED = fac1*(DPDstandardNthphi22L 
      - DPDstandardNthphi1L*Gt122 - DPDstandardNthphi2L*Gt222 - 
      DPDstandardNthphi3L*Gt322) + fac2*SQR(DPDstandardNthphi2L);
    
    CCTK_REAL cdphi223 CCTK_ATTRIBUTE_UNUSED = 
      DPDstandardNthphi2L*DPDstandardNthphi3L*fac2 + 
      fac1*(DPDstandardNthphi23L - DPDstandardNthphi1L*Gt123 - 
      DPDstandardNthphi2L*Gt223 - DPDstandardNthphi3L*Gt323);
    
    CCTK_REAL cdphi233 CCTK_ATTRIBUTE_UNUSED = fac1*(DPDstandardNthphi33L 
      - DPDstandardNthphi1L*Gt133 - DPDstandardNthphi2L*Gt233 - 
      DPDstandardNthphi3L*Gt333) + fac2*SQR(DPDstandardNthphi3L);
    
    CCTK_REAL Rphi11 CCTK_ATTRIBUTE_UNUSED = -2*(cdphi211 + 2*(-1 + 
      gt11L*gtu11)*SQR(cdphi1) + gt11L*(cdphi211*gtu11 + 
      4*(cdphi1*(cdphi2*gtu12 + cdphi3*gtu13) + cdphi2*cdphi3*gtu23) + 
      cdphi233*gtu33 + gtu22*(cdphi222 + 2*SQR(cdphi2)) + 2*(cdphi212*gtu12 + 
      cdphi213*gtu13 + cdphi223*gtu23 + gtu33*SQR(cdphi3))));
    
    CCTK_REAL Rphi12 CCTK_ATTRIBUTE_UNUSED = -2*(cdphi212 + 
      cdphi1*(cdphi2*(-2 + 4*gt12L*gtu12) + 4*gt12L*cdphi3*gtu13) + 
      gt12L*(cdphi211*gtu11 + 4*cdphi2*cdphi3*gtu23 + 2*(cdphi212*gtu12 + 
      cdphi213*gtu13 + cdphi223*gtu23 + gtu11*SQR(cdphi1)) + gtu22*(cdphi222 
      + 2*SQR(cdphi2)) + gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL Rphi13 CCTK_ATTRIBUTE_UNUSED = -2*(cdphi213 + 
      cdphi1*(4*gt13L*cdphi2*gtu12 + cdphi3*(-2 + 4*gt13L*gtu13)) + 
      gt13L*(cdphi211*gtu11 + 4*cdphi2*cdphi3*gtu23 + 2*(cdphi212*gtu12 + 
      cdphi213*gtu13 + cdphi223*gtu23 + gtu11*SQR(cdphi1)) + gtu22*(cdphi222 
      + 2*SQR(cdphi2)) + gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL Rphi22 CCTK_ATTRIBUTE_UNUSED = -2*(cdphi222 + 2*(-1 + 
      gt22L*gtu22)*SQR(cdphi2) + gt22L*(cdphi222*gtu22 + 
      4*(cdphi1*cdphi3*gtu13 + cdphi2*(cdphi1*gtu12 + cdphi3*gtu23)) + 
      cdphi233*gtu33 + gtu11*(cdphi211 + 2*SQR(cdphi1)) + 2*(cdphi212*gtu12 + 
      cdphi213*gtu13 + cdphi223*gtu23 + gtu33*SQR(cdphi3))));
    
    CCTK_REAL Rphi23 CCTK_ATTRIBUTE_UNUSED = -2*(cdphi223 + 
      cdphi2*(4*gt23L*cdphi1*gtu12 + cdphi3*(-2 + 4*gt23L*gtu23)) + 
      gt23L*(4*cdphi1*cdphi3*gtu13 + cdphi222*gtu22 + gtu11*(cdphi211 + 
      2*SQR(cdphi1)) + 2*(cdphi212*gtu12 + cdphi213*gtu13 + cdphi223*gtu23 + 
      gtu22*SQR(cdphi2)) + gtu33*(cdphi233 + 2*SQR(cdphi3))));
    
    CCTK_REAL Rphi33 CCTK_ATTRIBUTE_UNUSED = -2*(cdphi233 + 
      gt33L*((4*cdphi1*cdphi2 + 2*cdphi212)*gtu12 + 4*cdphi3*(cdphi1*gtu13 + 
      cdphi2*gtu23) + 2*(cdphi213*gtu13 + cdphi223*gtu23) + cdphi233*gtu33 + 
      gtu11*(cdphi211 + 2*SQR(cdphi1)) + gtu22*(cdphi222 + 2*SQR(cdphi2))) + 
      2*(-1 + gt33L*gtu33)*SQR(cdphi3));
    
    CCTK_REAL Atm11 CCTK_ATTRIBUTE_UNUSED = At11L*gtu11 + At12L*gtu12 + 
      At13L*gtu13;
    
    CCTK_REAL Atm21 CCTK_ATTRIBUTE_UNUSED = At11L*gtu12 + At12L*gtu22 + 
      At13L*gtu23;
    
    CCTK_REAL Atm31 CCTK_ATTRIBUTE_UNUSED = At11L*gtu13 + At12L*gtu23 + 
      At13L*gtu33;
    
    CCTK_REAL Atm12 CCTK_ATTRIBUTE_UNUSED = At12L*gtu11 + At22L*gtu12 + 
      At23L*gtu13;
    
    CCTK_REAL Atm22 CCTK_ATTRIBUTE_UNUSED = At12L*gtu12 + At22L*gtu22 + 
      At23L*gtu23;
    
    CCTK_REAL Atm32 CCTK_ATTRIBUTE_UNUSED = At12L*gtu13 + At22L*gtu23 + 
      At23L*gtu33;
    
    CCTK_REAL Atm13 CCTK_ATTRIBUTE_UNUSED = At13L*gtu11 + At23L*gtu12 + 
      At33L*gtu13;
    
    CCTK_REAL Atm23 CCTK_ATTRIBUTE_UNUSED = At13L*gtu12 + At23L*gtu22 + 
      At33L*gtu23;
    
    CCTK_REAL Atm33 CCTK_ATTRIBUTE_UNUSED = At13L*gtu13 + At23L*gtu23 + 
      At33L*gtu33;
    
    CCTK_REAL e4phi CCTK_ATTRIBUTE_UNUSED = 
      IfThen(conformalMethod,INV(SQR(phiL)),exp(4*phiL));
    
    CCTK_REAL em4phi CCTK_ATTRIBUTE_UNUSED = INV(e4phi);
    
    CCTK_REAL g11 CCTK_ATTRIBUTE_UNUSED = gt11L*e4phi;
    
    CCTK_REAL g12 CCTK_ATTRIBUTE_UNUSED = gt12L*e4phi;
    
    CCTK_REAL g13 CCTK_ATTRIBUTE_UNUSED = gt13L*e4phi;
    
    CCTK_REAL g22 CCTK_ATTRIBUTE_UNUSED = gt22L*e4phi;
    
    CCTK_REAL g23 CCTK_ATTRIBUTE_UNUSED = gt23L*e4phi;
    
    CCTK_REAL g33 CCTK_ATTRIBUTE_UNUSED = gt33L*e4phi;
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu11;
    
    CCTK_REAL gu12 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu12;
    
    CCTK_REAL gu13 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu13;
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu22;
    
    CCTK_REAL gu23 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu23;
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu33;
    
    CCTK_REAL R11 CCTK_ATTRIBUTE_UNUSED = Rphi11 + Rt11;
    
    CCTK_REAL R12 CCTK_ATTRIBUTE_UNUSED = Rphi12 + Rt12;
    
    CCTK_REAL R13 CCTK_ATTRIBUTE_UNUSED = Rphi13 + Rt13;
    
    CCTK_REAL R22 CCTK_ATTRIBUTE_UNUSED = Rphi22 + Rt22;
    
    CCTK_REAL R23 CCTK_ATTRIBUTE_UNUSED = Rphi23 + Rt23;
    
    CCTK_REAL R33 CCTK_ATTRIBUTE_UNUSED = Rphi33 + Rt33;
    
    CCTK_REAL Ats11 CCTK_ATTRIBUTE_UNUSED = -DPDstandardNthalpha11L + 
      DPDstandardNthalpha1L*(4*cdphi1 + Gt111) + DPDstandardNthalpha2L*Gt211 
      + DPDstandardNthalpha3L*Gt311 + alphaL*R11;
    
    CCTK_REAL Ats12 CCTK_ATTRIBUTE_UNUSED = -DPDstandardNthalpha12L + 
      2*(DPDstandardNthalpha2L*cdphi1 + DPDstandardNthalpha1L*cdphi2) + 
      DPDstandardNthalpha1L*Gt112 + DPDstandardNthalpha2L*Gt212 + 
      DPDstandardNthalpha3L*Gt312 + alphaL*R12;
    
    CCTK_REAL Ats13 CCTK_ATTRIBUTE_UNUSED = -DPDstandardNthalpha13L + 
      2*(DPDstandardNthalpha3L*cdphi1 + DPDstandardNthalpha1L*cdphi3) + 
      DPDstandardNthalpha1L*Gt113 + DPDstandardNthalpha2L*Gt213 + 
      DPDstandardNthalpha3L*Gt313 + alphaL*R13;
    
    CCTK_REAL Ats22 CCTK_ATTRIBUTE_UNUSED = -DPDstandardNthalpha22L + 
      DPDstandardNthalpha1L*Gt122 + DPDstandardNthalpha2L*(4*cdphi2 + Gt222) 
      + DPDstandardNthalpha3L*Gt322 + alphaL*R22;
    
    CCTK_REAL Ats23 CCTK_ATTRIBUTE_UNUSED = -DPDstandardNthalpha23L + 
      2*(DPDstandardNthalpha3L*cdphi2 + DPDstandardNthalpha2L*cdphi3) + 
      DPDstandardNthalpha1L*Gt123 + DPDstandardNthalpha2L*Gt223 + 
      DPDstandardNthalpha3L*Gt323 + alphaL*R23;
    
    CCTK_REAL Ats33 CCTK_ATTRIBUTE_UNUSED = -DPDstandardNthalpha33L + 
      DPDstandardNthalpha1L*Gt133 + DPDstandardNthalpha2L*Gt233 + 
      DPDstandardNthalpha3L*(4*cdphi3 + Gt333) + alphaL*R33;
    
    CCTK_REAL trAts CCTK_ATTRIBUTE_UNUSED = Ats11*gu11 + Ats22*gu22 + 
      2*(Ats12*gu12 + Ats13*gu13 + Ats23*gu23) + Ats33*gu33;
    
    CCTK_REAL At11rhsL CCTK_ATTRIBUTE_UNUSED = 
      2*(At11L*DPDstandardNthbeta11L + At12L*DPDstandardNthbeta21L + 
      At13L*DPDstandardNthbeta31L) - 
      0.666666666666666666666666666667*At11L*(DPDstandardNthbeta11L + 
      DPDstandardNthbeta22L + DPDstandardNthbeta33L) + alphaL*(At11L*trKL - 
      2*(At11L*Atm11 + At12L*Atm21 + At13L*Atm31)) + em4phi*(Ats11 - 
      0.333333333333333333333333333333*g11*trAts);
    
    CCTK_REAL At12rhsL CCTK_ATTRIBUTE_UNUSED = At11L*DPDstandardNthbeta12L 
      + At22L*DPDstandardNthbeta21L + At23L*DPDstandardNthbeta31L + 
      At13L*DPDstandardNthbeta32L + At12L*(DPDstandardNthbeta11L + 
      DPDstandardNthbeta22L - 
      0.666666666666666666666666666667*(DPDstandardNthbeta11L + 
      DPDstandardNthbeta22L + DPDstandardNthbeta33L)) + alphaL*(At12L*trKL - 
      2*(At11L*Atm12 + At12L*Atm22 + At13L*Atm32)) + em4phi*(Ats12 - 
      0.333333333333333333333333333333*g12*trAts);
    
    CCTK_REAL At13rhsL CCTK_ATTRIBUTE_UNUSED = At11L*DPDstandardNthbeta13L 
      + At23L*DPDstandardNthbeta21L + At12L*DPDstandardNthbeta23L + 
      At33L*DPDstandardNthbeta31L + At13L*(DPDstandardNthbeta11L + 
      DPDstandardNthbeta33L - 
      0.666666666666666666666666666667*(DPDstandardNthbeta11L + 
      DPDstandardNthbeta22L + DPDstandardNthbeta33L)) + alphaL*(At13L*trKL - 
      2*(At11L*Atm13 + At12L*Atm23 + At13L*Atm33)) + em4phi*(Ats13 - 
      0.333333333333333333333333333333*g13*trAts);
    
    CCTK_REAL At22rhsL CCTK_ATTRIBUTE_UNUSED = 
      2*(At12L*DPDstandardNthbeta12L + At22L*DPDstandardNthbeta22L + 
      At23L*DPDstandardNthbeta32L) - 
      0.666666666666666666666666666667*At22L*(DPDstandardNthbeta11L + 
      DPDstandardNthbeta22L + DPDstandardNthbeta33L) + alphaL*(At22L*trKL - 
      2*(At12L*Atm12 + At22L*Atm22 + At23L*Atm32)) + em4phi*(Ats22 - 
      0.333333333333333333333333333333*g22*trAts);
    
    CCTK_REAL At23rhsL CCTK_ATTRIBUTE_UNUSED = At13L*DPDstandardNthbeta12L 
      + At12L*DPDstandardNthbeta13L + At22L*DPDstandardNthbeta23L + 
      At33L*DPDstandardNthbeta32L + At23L*(DPDstandardNthbeta22L + 
      DPDstandardNthbeta33L - 
      0.666666666666666666666666666667*(DPDstandardNthbeta11L + 
      DPDstandardNthbeta22L + DPDstandardNthbeta33L)) + alphaL*(At23L*trKL - 
      2*(At12L*Atm13 + At22L*Atm23 + At23L*Atm33)) + em4phi*(Ats23 - 
      0.333333333333333333333333333333*g23*trAts);
    
    CCTK_REAL At33rhsL CCTK_ATTRIBUTE_UNUSED = 
      -0.666666666666666666666666666667*At33L*(DPDstandardNthbeta11L + 
      DPDstandardNthbeta22L + DPDstandardNthbeta33L) + 
      2*(At13L*DPDstandardNthbeta13L + At23L*DPDstandardNthbeta23L + 
      At33L*DPDstandardNthbeta33L) + alphaL*(At33L*trKL - 2*(At13L*Atm13 + 
      At23L*Atm23 + At33L*Atm33)) + em4phi*(Ats33 - 
      0.333333333333333333333333333333*g33*trAts);
    
    /* Copy local copies back to grid functions */
    At11rhs[index] = At11rhsL;
    At12rhs[index] = At12rhsL;
    At13rhs[index] = At13rhsL;
    At22rhs[index] = At22rhsL;
    At23rhs[index] = At23rhsL;
    At33rhs[index] = At33rhsL;
  }
  CCTK_ENDLOOP3(HOST_ML_BSSN_RHS_NonDerivatives2);
}

extern "C" void HOST_ML_BSSN_RHS_NonDerivatives2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering HOST_ML_BSSN_RHS_NonDerivatives2_Body");
  }
  
  if (cctk_iteration % HOST_ML_BSSN_RHS_NonDerivatives2_calc_every != HOST_ML_BSSN_RHS_NonDerivatives2_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN::DPDstandardNthalpha11_group",
    "ML_BSSN::DPDstandardNthalpha12_group",
    "ML_BSSN::DPDstandardNthalpha13_group",
    "ML_BSSN::DPDstandardNthalpha1_group",
    "ML_BSSN::DPDstandardNthalpha22_group",
    "ML_BSSN::DPDstandardNthalpha23_group",
    "ML_BSSN::DPDstandardNthalpha2_group",
    "ML_BSSN::DPDstandardNthalpha33_group",
    "ML_BSSN::DPDstandardNthalpha3_group",
    "ML_BSSN::DPDstandardNthbeta11_group",
    "ML_BSSN::DPDstandardNthbeta12_group",
    "ML_BSSN::DPDstandardNthbeta13_group",
    "ML_BSSN::DPDstandardNthbeta21_group",
    "ML_BSSN::DPDstandardNthbeta22_group",
    "ML_BSSN::DPDstandardNthbeta23_group",
    "ML_BSSN::DPDstandardNthbeta31_group",
    "ML_BSSN::DPDstandardNthbeta32_group",
    "ML_BSSN::DPDstandardNthbeta33_group",
    "ML_BSSN::DPDstandardNthgt1111_group",
    "ML_BSSN::DPDstandardNthgt1112_group",
    "ML_BSSN::DPDstandardNthgt1113_group",
    "ML_BSSN::DPDstandardNthgt111_group",
    "ML_BSSN::DPDstandardNthgt1122_group",
    "ML_BSSN::DPDstandardNthgt1123_group",
    "ML_BSSN::DPDstandardNthgt112_group",
    "ML_BSSN::DPDstandardNthgt1133_group",
    "ML_BSSN::DPDstandardNthgt113_group",
    "ML_BSSN::DPDstandardNthgt1211_group",
    "ML_BSSN::DPDstandardNthgt1212_group",
    "ML_BSSN::DPDstandardNthgt1213_group",
    "ML_BSSN::DPDstandardNthgt121_group",
    "ML_BSSN::DPDstandardNthgt1222_group",
    "ML_BSSN::DPDstandardNthgt1223_group",
    "ML_BSSN::DPDstandardNthgt122_group",
    "ML_BSSN::DPDstandardNthgt1233_group",
    "ML_BSSN::DPDstandardNthgt123_group",
    "ML_BSSN::DPDstandardNthgt1311_group",
    "ML_BSSN::DPDstandardNthgt1312_group",
    "ML_BSSN::DPDstandardNthgt1313_group",
    "ML_BSSN::DPDstandardNthgt131_group",
    "ML_BSSN::DPDstandardNthgt1322_group",
    "ML_BSSN::DPDstandardNthgt1323_group",
    "ML_BSSN::DPDstandardNthgt132_group",
    "ML_BSSN::DPDstandardNthgt1333_group",
    "ML_BSSN::DPDstandardNthgt133_group",
    "ML_BSSN::DPDstandardNthgt2211_group",
    "ML_BSSN::DPDstandardNthgt2212_group",
    "ML_BSSN::DPDstandardNthgt2213_group",
    "ML_BSSN::DPDstandardNthgt221_group",
    "ML_BSSN::DPDstandardNthgt2222_group",
    "ML_BSSN::DPDstandardNthgt2223_group",
    "ML_BSSN::DPDstandardNthgt222_group",
    "ML_BSSN::DPDstandardNthgt2233_group",
    "ML_BSSN::DPDstandardNthgt223_group",
    "ML_BSSN::DPDstandardNthgt2311_group",
    "ML_BSSN::DPDstandardNthgt2312_group",
    "ML_BSSN::DPDstandardNthgt2313_group",
    "ML_BSSN::DPDstandardNthgt231_group",
    "ML_BSSN::DPDstandardNthgt2322_group",
    "ML_BSSN::DPDstandardNthgt2323_group",
    "ML_BSSN::DPDstandardNthgt232_group",
    "ML_BSSN::DPDstandardNthgt2333_group",
    "ML_BSSN::DPDstandardNthgt233_group",
    "ML_BSSN::DPDstandardNthgt3311_group",
    "ML_BSSN::DPDstandardNthgt3312_group",
    "ML_BSSN::DPDstandardNthgt3313_group",
    "ML_BSSN::DPDstandardNthgt331_group",
    "ML_BSSN::DPDstandardNthgt3322_group",
    "ML_BSSN::DPDstandardNthgt3323_group",
    "ML_BSSN::DPDstandardNthgt332_group",
    "ML_BSSN::DPDstandardNthgt3333_group",
    "ML_BSSN::DPDstandardNthgt333_group",
    "ML_BSSN::DPDstandardNthphi11_group",
    "ML_BSSN::DPDstandardNthphi12_group",
    "ML_BSSN::DPDstandardNthphi13_group",
    "ML_BSSN::DPDstandardNthphi1_group",
    "ML_BSSN::DPDstandardNthphi22_group",
    "ML_BSSN::DPDstandardNthphi23_group",
    "ML_BSSN::DPDstandardNthphi2_group",
    "ML_BSSN::DPDstandardNthphi33_group",
    "ML_BSSN::DPDstandardNthphi3_group",
    "ML_BSSN::DPDstandardNthXt11_group",
    "ML_BSSN::DPDstandardNthXt12_group",
    "ML_BSSN::DPDstandardNthXt13_group",
    "ML_BSSN::DPDstandardNthXt21_group",
    "ML_BSSN::DPDstandardNthXt22_group",
    "ML_BSSN::DPDstandardNthXt23_group",
    "ML_BSSN::DPDstandardNthXt31_group",
    "ML_BSSN::DPDstandardNthXt32_group",
    "ML_BSSN::DPDstandardNthXt33_group",
    "ML_BSSN::ML_curv",
    "ML_BSSN::ML_curvrhs",
    "ML_BSSN::ML_lapse",
    "ML_BSSN::ML_log_confac",
    "ML_BSSN::ML_metric",
    "ML_BSSN::ML_shift",
    "ML_BSSN::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "HOST_ML_BSSN_RHS_NonDerivatives2", 97, groups);
  
  
  GenericFD_LoopOverInterior(cctkGH, HOST_ML_BSSN_RHS_NonDerivatives2_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving HOST_ML_BSSN_RHS_NonDerivatives2_Body");
  }
}
