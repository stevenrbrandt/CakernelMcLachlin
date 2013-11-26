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

extern "C" void ML_BSSN_convertFromADMBaseGamma_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_convertFromADMBaseGamma_calc_every != ML_BSSN_convertFromADMBaseGamma_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_dtlapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_dtshift.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_Gamma.");
  return;
}

static void ML_BSSN_convertFromADMBaseGamma_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(ML_BSSN_convertFromADMBaseGamma,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    // ++vec_iter_counter;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpha[index];
    CCTK_REAL beta1L CCTK_ATTRIBUTE_UNUSED = beta1[index];
    CCTK_REAL beta2L CCTK_ATTRIBUTE_UNUSED = beta2[index];
    CCTK_REAL beta3L CCTK_ATTRIBUTE_UNUSED = beta3[index];
    CCTK_REAL dtalpL CCTK_ATTRIBUTE_UNUSED = dtalp[index];
    CCTK_REAL dtbetaxL CCTK_ATTRIBUTE_UNUSED = dtbetax[index];
    CCTK_REAL dtbetayL CCTK_ATTRIBUTE_UNUSED = dtbetay[index];
    CCTK_REAL dtbetazL CCTK_ATTRIBUTE_UNUSED = dtbetaz[index];
    CCTK_REAL gt11L CCTK_ATTRIBUTE_UNUSED = gt11[index];
    CCTK_REAL gt12L CCTK_ATTRIBUTE_UNUSED = gt12[index];
    CCTK_REAL gt13L CCTK_ATTRIBUTE_UNUSED = gt13[index];
    CCTK_REAL gt22L CCTK_ATTRIBUTE_UNUSED = gt22[index];
    CCTK_REAL gt23L CCTK_ATTRIBUTE_UNUSED = gt23[index];
    CCTK_REAL gt33L CCTK_ATTRIBUTE_UNUSED = gt33[index];
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    const CCTK_REAL PDupwindNthAnti1alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti1(&alpha[index]);
    const CCTK_REAL PDupwindNthSymm1alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm1(&alpha[index]);
    const CCTK_REAL PDupwindNthAnti2alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti2(&alpha[index]);
    const CCTK_REAL PDupwindNthSymm2alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm2(&alpha[index]);
    const CCTK_REAL PDupwindNthAnti3alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti3(&alpha[index]);
    const CCTK_REAL PDupwindNthSymm3alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm3(&alpha[index]);
    const CCTK_REAL PDupwindNthAnti1beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti1(&beta1[index]);
    const CCTK_REAL PDupwindNthSymm1beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm1(&beta1[index]);
    const CCTK_REAL PDupwindNthAnti2beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti2(&beta1[index]);
    const CCTK_REAL PDupwindNthSymm2beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm2(&beta1[index]);
    const CCTK_REAL PDupwindNthAnti3beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti3(&beta1[index]);
    const CCTK_REAL PDupwindNthSymm3beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm3(&beta1[index]);
    const CCTK_REAL PDupwindNthAnti1beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti1(&beta2[index]);
    const CCTK_REAL PDupwindNthSymm1beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm1(&beta2[index]);
    const CCTK_REAL PDupwindNthAnti2beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti2(&beta2[index]);
    const CCTK_REAL PDupwindNthSymm2beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm2(&beta2[index]);
    const CCTK_REAL PDupwindNthAnti3beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti3(&beta2[index]);
    const CCTK_REAL PDupwindNthSymm3beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm3(&beta2[index]);
    const CCTK_REAL PDupwindNthAnti1beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti1(&beta3[index]);
    const CCTK_REAL PDupwindNthSymm1beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm1(&beta3[index]);
    const CCTK_REAL PDupwindNthAnti2beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti2(&beta3[index]);
    const CCTK_REAL PDupwindNthSymm2beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm2(&beta3[index]);
    const CCTK_REAL PDupwindNthAnti3beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti3(&beta3[index]);
    const CCTK_REAL PDupwindNthSymm3beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm3(&beta3[index]);
    const CCTK_REAL PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt11[index]);
    const CCTK_REAL PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt11[index]);
    const CCTK_REAL PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt11[index]);
    const CCTK_REAL PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt12[index]);
    const CCTK_REAL PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt12[index]);
    const CCTK_REAL PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt12[index]);
    const CCTK_REAL PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt13[index]);
    const CCTK_REAL PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt13[index]);
    const CCTK_REAL PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt13[index]);
    const CCTK_REAL PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt22[index]);
    const CCTK_REAL PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt22[index]);
    const CCTK_REAL PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt22[index]);
    const CCTK_REAL PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt23[index]);
    const CCTK_REAL PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt23[index]);
    const CCTK_REAL PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt23[index]);
    const CCTK_REAL PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt33[index]);
    const CCTK_REAL PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt33[index]);
    const CCTK_REAL PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt33[index]);
    
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
    
    CCTK_REAL Gt111 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu11*PDstandardNth1gt11 
      + 2*(gtu12*PDstandardNth1gt12 + gtu13*PDstandardNth1gt13) - 
      gtu12*PDstandardNth2gt11 - gtu13*PDstandardNth3gt11);
    
    CCTK_REAL Gt211 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu12*PDstandardNth1gt11 
      + 2*(gtu22*PDstandardNth1gt12 + gtu23*PDstandardNth1gt13) - 
      gtu22*PDstandardNth2gt11 - gtu23*PDstandardNth3gt11);
    
    CCTK_REAL Gt311 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu13*PDstandardNth1gt11 
      + 2*(gtu23*PDstandardNth1gt12 + gtu33*PDstandardNth1gt13) - 
      gtu23*PDstandardNth2gt11 - gtu33*PDstandardNth3gt11);
    
    CCTK_REAL Gt112 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu12*PDstandardNth1gt22 
      + gtu11*PDstandardNth2gt11 + gtu13*(PDstandardNth1gt23 + 
      PDstandardNth2gt13 - PDstandardNth3gt12));
    
    CCTK_REAL Gt212 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu22*PDstandardNth1gt22 
      + gtu12*PDstandardNth2gt11 + gtu23*(PDstandardNth1gt23 + 
      PDstandardNth2gt13 - PDstandardNth3gt12));
    
    CCTK_REAL Gt312 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu23*PDstandardNth1gt22 
      + gtu13*PDstandardNth2gt11 + gtu33*(PDstandardNth1gt23 + 
      PDstandardNth2gt13 - PDstandardNth3gt12));
    
    CCTK_REAL Gt113 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu13*PDstandardNth1gt33 
      + gtu11*PDstandardNth3gt11 + gtu12*(PDstandardNth1gt23 - 
      PDstandardNth2gt13 + PDstandardNth3gt12));
    
    CCTK_REAL Gt213 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu23*PDstandardNth1gt33 
      + gtu12*PDstandardNth3gt11 + gtu22*(PDstandardNth1gt23 - 
      PDstandardNth2gt13 + PDstandardNth3gt12));
    
    CCTK_REAL Gt313 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu33*PDstandardNth1gt33 
      + gtu13*PDstandardNth3gt11 + gtu23*(PDstandardNth1gt23 - 
      PDstandardNth2gt13 + PDstandardNth3gt12));
    
    CCTK_REAL Gt122 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu11*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + 
      gtu12*PDstandardNth2gt22 + gtu13*(2*PDstandardNth2gt23 - 
      PDstandardNth3gt22));
    
    CCTK_REAL Gt222 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu12*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + 
      gtu22*PDstandardNth2gt22 + gtu23*(2*PDstandardNth2gt23 - 
      PDstandardNth3gt22));
    
    CCTK_REAL Gt322 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu13*(-PDstandardNth1gt22 + 2*PDstandardNth2gt12) + 
      gtu23*PDstandardNth2gt22 + gtu33*(2*PDstandardNth2gt23 - 
      PDstandardNth3gt22));
    
    CCTK_REAL Gt123 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu13*PDstandardNth2gt33 
      + gtu11*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) 
      + gtu12*PDstandardNth3gt22);
    
    CCTK_REAL Gt223 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu23*PDstandardNth2gt33 
      + gtu12*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) 
      + gtu22*PDstandardNth3gt22);
    
    CCTK_REAL Gt323 CCTK_ATTRIBUTE_UNUSED = 0.5*(gtu33*PDstandardNth2gt33 
      + gtu13*(-PDstandardNth1gt23 + PDstandardNth2gt13 + PDstandardNth3gt12) 
      + gtu23*PDstandardNth3gt22);
    
    CCTK_REAL Gt133 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu11*(-PDstandardNth1gt33 + 2*PDstandardNth3gt13) + 
      gtu12*(-PDstandardNth2gt33 + 2*PDstandardNth3gt23) + 
      gtu13*PDstandardNth3gt33);
    
    CCTK_REAL Gt233 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu12*(-PDstandardNth1gt33 + 2*PDstandardNth3gt13) + 
      gtu22*(-PDstandardNth2gt33 + 2*PDstandardNth3gt23) + 
      gtu23*PDstandardNth3gt33);
    
    CCTK_REAL Gt333 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gtu13*(-PDstandardNth1gt33 + 2*PDstandardNth3gt13) + 
      gtu23*(-PDstandardNth2gt33 + 2*PDstandardNth3gt23) + 
      gtu33*PDstandardNth3gt33);
    
    CCTK_REAL Xt1L CCTK_ATTRIBUTE_UNUSED = Gt111*gtu11 + Gt122*gtu22 + 
      2*(Gt112*gtu12 + Gt113*gtu13 + Gt123*gtu23) + Gt133*gtu33;
    
    CCTK_REAL Xt2L CCTK_ATTRIBUTE_UNUSED = Gt211*gtu11 + Gt222*gtu22 + 
      2*(Gt212*gtu12 + Gt213*gtu13 + Gt223*gtu23) + Gt233*gtu33;
    
    CCTK_REAL Xt3L CCTK_ATTRIBUTE_UNUSED = Gt311*gtu11 + Gt322*gtu22 + 
      2*(Gt312*gtu12 + Gt313*gtu13 + Gt323*gtu23) + Gt333*gtu33;
    
    CCTK_REAL AL CCTK_ATTRIBUTE_UNUSED = IfThen(LapseACoeff != 
      0,-(INV(ToReal(harmonicF))*pow(alphaL,-ToReal(harmonicN))*(dtalpL - 
      (beta1L*PDupwindNthAnti1alpha + beta2L*PDupwindNthAnti2alpha + 
      beta3L*PDupwindNthAnti3alpha + PDupwindNthSymm1alpha*fabs(beta1L) + 
      PDupwindNthSymm2alpha*fabs(beta2L) + 
      PDupwindNthSymm3alpha*fabs(beta3L))*ToReal(LapseAdvectionCoeff))),0);
    
    CCTK_REAL theta CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL B1L CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL B2L CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL B3L CCTK_ATTRIBUTE_UNUSED;
    
    if (ShiftBCoeff*ShiftGammaCoeff != 0)
    {
      B1L = INV(theta*ToReal(ShiftGammaCoeff))*(dtbetaxL - 
        (beta1L*PDupwindNthAnti1beta1 + beta2L*PDupwindNthAnti2beta1 + 
        beta3L*PDupwindNthAnti3beta1 + PDupwindNthSymm1beta1*fabs(beta1L) + 
        PDupwindNthSymm2beta1*fabs(beta2L) + 
        PDupwindNthSymm3beta1*fabs(beta3L))*ToReal(ShiftAdvectionCoeff));
      
      B2L = INV(theta*ToReal(ShiftGammaCoeff))*(dtbetayL - 
        (beta1L*PDupwindNthAnti1beta2 + beta2L*PDupwindNthAnti2beta2 + 
        beta3L*PDupwindNthAnti3beta2 + PDupwindNthSymm1beta2*fabs(beta1L) + 
        PDupwindNthSymm2beta2*fabs(beta2L) + 
        PDupwindNthSymm3beta2*fabs(beta3L))*ToReal(ShiftAdvectionCoeff));
      
      B3L = INV(theta*ToReal(ShiftGammaCoeff))*(dtbetazL - 
        (beta1L*PDupwindNthAnti1beta3 + beta2L*PDupwindNthAnti2beta3 + 
        beta3L*PDupwindNthAnti3beta3 + PDupwindNthSymm1beta3*fabs(beta1L) + 
        PDupwindNthSymm2beta3*fabs(beta2L) + 
        PDupwindNthSymm3beta3*fabs(beta3L))*ToReal(ShiftAdvectionCoeff));
    }
    else
    {
      B1L = 0;
      
      B2L = 0;
      
      B3L = 0;
    }
    
    /* Copy local copies back to grid functions */
    A[index] = AL;
    B1[index] = B1L;
    B2[index] = B2L;
    B3[index] = B3L;
    Xt1[index] = Xt1L;
    Xt2[index] = Xt2L;
    Xt3[index] = Xt3L;
  }
  CCTK_ENDLOOP3(ML_BSSN_convertFromADMBaseGamma);
}

extern "C" void ML_BSSN_convertFromADMBaseGamma(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_convertFromADMBaseGamma_Body");
  }
  
  if (cctk_iteration % ML_BSSN_convertFromADMBaseGamma_calc_every != ML_BSSN_convertFromADMBaseGamma_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::dtlapse",
    "ADMBase::dtshift",
    "ML_BSSN::ML_dtlapse",
    "ML_BSSN::ML_dtshift",
    "ML_BSSN::ML_Gamma",
    "ML_BSSN::ML_lapse",
    "ML_BSSN::ML_metric",
    "ML_BSSN::ML_shift"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_convertFromADMBaseGamma", 8, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_convertFromADMBaseGamma", 5, 5, 5);
  
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_convertFromADMBaseGamma_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_convertFromADMBaseGamma_Body");
  }
}
