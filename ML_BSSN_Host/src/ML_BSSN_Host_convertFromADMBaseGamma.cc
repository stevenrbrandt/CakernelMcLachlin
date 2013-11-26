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

extern "C" void ML_BSSN_Host_convertFromADMBaseGamma_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_Host_convertFromADMBaseGamma_calc_every != ML_BSSN_Host_convertFromADMBaseGamma_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_dtlapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_dtshift.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_Gamma.");
  return;
}

static void ML_BSSN_Host_convertFromADMBaseGamma_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL_VEC dx CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(0));
  const CCTK_REAL_VEC dy CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(1));
  const CCTK_REAL_VEC dz CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(2));
  const CCTK_REAL_VEC dt CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_TIME);
  const CCTK_REAL_VEC t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);
  const CCTK_REAL_VEC dxi CCTK_ATTRIBUTE_UNUSED = INV(dx);
  const CCTK_REAL_VEC dyi CCTK_ATTRIBUTE_UNUSED = INV(dy);
  const CCTK_REAL_VEC dzi CCTK_ATTRIBUTE_UNUSED = INV(dz);
  const CCTK_REAL_VEC khalf CCTK_ATTRIBUTE_UNUSED = ToReal(0.5);
  const CCTK_REAL_VEC kthird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(0.333333333333333333333333333333);
  const CCTK_REAL_VEC ktwothird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(0.666666666666666666666666666667);
  const CCTK_REAL_VEC kfourthird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(1.33333333333333333333333333333);
  const CCTK_REAL_VEC hdxi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dxi,ToReal(0.5));
  const CCTK_REAL_VEC hdyi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dyi,ToReal(0.5));
  const CCTK_REAL_VEC hdzi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dzi,ToReal(0.5));
  
  /* Initialize predefined quantities */
  const CCTK_REAL_VEC p1o1024dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dx);
  const CCTK_REAL_VEC p1o1024dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dy);
  const CCTK_REAL_VEC p1o1024dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dz);
  const CCTK_REAL_VEC p1o1680dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dx);
  const CCTK_REAL_VEC p1o1680dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dy);
  const CCTK_REAL_VEC p1o1680dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dz);
  const CCTK_REAL_VEC p1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dx);
  const CCTK_REAL_VEC p1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dy);
  const CCTK_REAL_VEC p1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dz);
  const CCTK_REAL_VEC p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dx,dx));
  const CCTK_REAL_VEC p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dy,dy));
  const CCTK_REAL_VEC p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dz,dz));
  const CCTK_REAL_VEC p1o560dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dx);
  const CCTK_REAL_VEC p1o560dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dy);
  const CCTK_REAL_VEC p1o560dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dz);
  const CCTK_REAL_VEC p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dy,dx));
  const CCTK_REAL_VEC p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dz,dx));
  const CCTK_REAL_VEC p1o705600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dz,dy));
  const CCTK_REAL_VEC p1o840dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dx);
  const CCTK_REAL_VEC p1o840dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dy);
  const CCTK_REAL_VEC p1o840dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dz);
  const CCTK_REAL_VEC p1odx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dx);
  const CCTK_REAL_VEC p1ody CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dy);
  const CCTK_REAL_VEC p1odz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dz);
  const CCTK_REAL_VEC pm1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dx);
  const CCTK_REAL_VEC pm1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dy);
  const CCTK_REAL_VEC pm1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dz);
  
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
  CCTK_LOOP3STR(ML_BSSN_Host_convertFromADMBaseGamma,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    // vec_iter_counter+=CCTK_REAL_VEC_SIZE;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = vec_load(alpha[index]);
    CCTK_REAL_VEC beta1L CCTK_ATTRIBUTE_UNUSED = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L CCTK_ATTRIBUTE_UNUSED = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L CCTK_ATTRIBUTE_UNUSED = vec_load(beta3[index]);
    CCTK_REAL_VEC dtalpL CCTK_ATTRIBUTE_UNUSED = vec_load(dtalp[index]);
    CCTK_REAL_VEC dtbetaxL CCTK_ATTRIBUTE_UNUSED = vec_load(dtbetax[index]);
    CCTK_REAL_VEC dtbetayL CCTK_ATTRIBUTE_UNUSED = vec_load(dtbetay[index]);
    CCTK_REAL_VEC dtbetazL CCTK_ATTRIBUTE_UNUSED = vec_load(dtbetaz[index]);
    CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = vec_load(gt33[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    const CCTK_REAL_VEC PDupwindNthAnti1alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti1(&alpha[index]);
    const CCTK_REAL_VEC PDupwindNthSymm1alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm1(&alpha[index]);
    const CCTK_REAL_VEC PDupwindNthAnti2alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti2(&alpha[index]);
    const CCTK_REAL_VEC PDupwindNthSymm2alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm2(&alpha[index]);
    const CCTK_REAL_VEC PDupwindNthAnti3alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti3(&alpha[index]);
    const CCTK_REAL_VEC PDupwindNthSymm3alpha CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm3(&alpha[index]);
    const CCTK_REAL_VEC PDupwindNthAnti1beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti1(&beta1[index]);
    const CCTK_REAL_VEC PDupwindNthSymm1beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm1(&beta1[index]);
    const CCTK_REAL_VEC PDupwindNthAnti2beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti2(&beta1[index]);
    const CCTK_REAL_VEC PDupwindNthSymm2beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm2(&beta1[index]);
    const CCTK_REAL_VEC PDupwindNthAnti3beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti3(&beta1[index]);
    const CCTK_REAL_VEC PDupwindNthSymm3beta1 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm3(&beta1[index]);
    const CCTK_REAL_VEC PDupwindNthAnti1beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti1(&beta2[index]);
    const CCTK_REAL_VEC PDupwindNthSymm1beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm1(&beta2[index]);
    const CCTK_REAL_VEC PDupwindNthAnti2beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti2(&beta2[index]);
    const CCTK_REAL_VEC PDupwindNthSymm2beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm2(&beta2[index]);
    const CCTK_REAL_VEC PDupwindNthAnti3beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti3(&beta2[index]);
    const CCTK_REAL_VEC PDupwindNthSymm3beta2 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm3(&beta2[index]);
    const CCTK_REAL_VEC PDupwindNthAnti1beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti1(&beta3[index]);
    const CCTK_REAL_VEC PDupwindNthSymm1beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm1(&beta3[index]);
    const CCTK_REAL_VEC PDupwindNthAnti2beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti2(&beta3[index]);
    const CCTK_REAL_VEC PDupwindNthSymm2beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm2(&beta3[index]);
    const CCTK_REAL_VEC PDupwindNthAnti3beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthAnti3(&beta3[index]);
    const CCTK_REAL_VEC PDupwindNthSymm3beta3 CCTK_ATTRIBUTE_UNUSED = PDupwindNthSymm3(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt33[index]);
    const CCTK_REAL_VEC PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt33[index]);
    const CCTK_REAL_VEC PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt33[index]);
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = kisgn(beta1L);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = kisgn(beta2L);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = kisgn(beta3L);
    
    CCTK_REAL_VEC detgt CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC gtu11 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt22L,gt33L,kmul(gt23L,gt23L)),detgt);
    
    CCTK_REAL_VEC gtu12 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt13L,gt23L,kmul(gt12L,gt33L)),detgt);
    
    CCTK_REAL_VEC gtu13 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt12L,gt23L,kmul(gt13L,gt22L)),detgt);
    
    CCTK_REAL_VEC gtu22 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt11L,gt33L,kmul(gt13L,gt13L)),detgt);
    
    CCTK_REAL_VEC gtu23 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt12L,gt13L,kmul(gt11L,gt23L)),detgt);
    
    CCTK_REAL_VEC gtu33 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt11L,gt22L,kmul(gt12L,gt12L)),detgt);
    
    CCTK_REAL_VEC Gt111 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gtu11,PDstandardNth1gt11,knmsub(gtu12,PDstandardNth2gt11,kmsub(kmadd(gtu12,PDstandardNth1gt12,kmul(gtu13,PDstandardNth1gt13)),ToReal(2),kmul(gtu13,PDstandardNth3gt11)))));
    
    CCTK_REAL_VEC Gt211 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gtu12,PDstandardNth1gt11,knmsub(gtu22,PDstandardNth2gt11,kmsub(kmadd(gtu22,PDstandardNth1gt12,kmul(gtu23,PDstandardNth1gt13)),ToReal(2),kmul(gtu23,PDstandardNth3gt11)))));
    
    CCTK_REAL_VEC Gt311 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gtu13,PDstandardNth1gt11,knmsub(gtu23,PDstandardNth2gt11,kmsub(kmadd(gtu23,PDstandardNth1gt12,kmul(gtu33,PDstandardNth1gt13)),ToReal(2),kmul(gtu33,PDstandardNth3gt11)))));
    
    CCTK_REAL_VEC Gt112 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu12,PDstandardNth1gt22,kmadd(gtu11,PDstandardNth2gt11,kmul(gtu13,kadd(PDstandardNth1gt23,ksub(PDstandardNth2gt13,PDstandardNth3gt12))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt212 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu22,PDstandardNth1gt22,kmadd(gtu12,PDstandardNth2gt11,kmul(gtu23,kadd(PDstandardNth1gt23,ksub(PDstandardNth2gt13,PDstandardNth3gt12))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt312 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu23,PDstandardNth1gt22,kmadd(gtu13,PDstandardNth2gt11,kmul(gtu33,kadd(PDstandardNth1gt23,ksub(PDstandardNth2gt13,PDstandardNth3gt12))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt113 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu13,PDstandardNth1gt33,kmadd(gtu11,PDstandardNth3gt11,kmul(gtu12,kadd(PDstandardNth1gt23,ksub(PDstandardNth3gt12,PDstandardNth2gt13))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt213 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu23,PDstandardNth1gt33,kmadd(gtu12,PDstandardNth3gt11,kmul(gtu22,kadd(PDstandardNth1gt23,ksub(PDstandardNth3gt12,PDstandardNth2gt13))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt313 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu33,PDstandardNth1gt33,kmadd(gtu13,PDstandardNth3gt11,kmul(gtu23,kadd(PDstandardNth1gt23,ksub(PDstandardNth3gt12,PDstandardNth2gt13))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt122 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gtu12,PDstandardNth2gt22,kmadd(gtu11,kmsub(PDstandardNth2gt12,ToReal(2),PDstandardNth1gt22),kmul(gtu13,kmsub(PDstandardNth2gt23,ToReal(2),PDstandardNth3gt22)))));
    
    CCTK_REAL_VEC Gt222 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gtu22,PDstandardNth2gt22,kmadd(gtu12,kmsub(PDstandardNth2gt12,ToReal(2),PDstandardNth1gt22),kmul(gtu23,kmsub(PDstandardNth2gt23,ToReal(2),PDstandardNth3gt22)))));
    
    CCTK_REAL_VEC Gt322 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gtu23,PDstandardNth2gt22,kmadd(gtu13,kmsub(PDstandardNth2gt12,ToReal(2),PDstandardNth1gt22),kmul(gtu33,kmsub(PDstandardNth2gt23,ToReal(2),PDstandardNth3gt22)))));
    
    CCTK_REAL_VEC Gt123 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu13,PDstandardNth2gt33,kmadd(gtu12,PDstandardNth3gt22,kmul(gtu11,kadd(PDstandardNth2gt13,ksub(PDstandardNth3gt12,PDstandardNth1gt23))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt223 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu23,PDstandardNth2gt33,kmadd(gtu22,PDstandardNth3gt22,kmul(gtu12,kadd(PDstandardNth2gt13,ksub(PDstandardNth3gt12,PDstandardNth1gt23))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt323 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmadd(gtu33,PDstandardNth2gt33,kmadd(gtu23,PDstandardNth3gt22,kmul(gtu13,kadd(PDstandardNth2gt13,ksub(PDstandardNth3gt12,PDstandardNth1gt23))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt133 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gtu13,PDstandardNth3gt33,kmadd(gtu11,kmsub(PDstandardNth3gt13,ToReal(2),PDstandardNth1gt33),kmul(gtu12,kmsub(PDstandardNth3gt23,ToReal(2),PDstandardNth2gt33)))));
    
    CCTK_REAL_VEC Gt233 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gtu23,PDstandardNth3gt33,kmadd(gtu12,kmsub(PDstandardNth3gt13,ToReal(2),PDstandardNth1gt33),kmul(gtu22,kmsub(PDstandardNth3gt23,ToReal(2),PDstandardNth2gt33)))));
    
    CCTK_REAL_VEC Gt333 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gtu33,PDstandardNth3gt33,kmadd(gtu13,kmsub(PDstandardNth3gt13,ToReal(2),PDstandardNth1gt33),kmul(gtu23,kmsub(PDstandardNth3gt23,ToReal(2),PDstandardNth2gt33)))));
    
    CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmadd(Gt133,gtu33,kmul(kmadd(Gt112,gtu12,kmadd(Gt113,gtu13,kmul(Gt123,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmadd(Gt233,gtu33,kmul(kmadd(Gt212,gtu12,kmadd(Gt213,gtu13,kmul(Gt223,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmadd(Gt333,gtu33,kmul(kmadd(Gt312,gtu12,kmadd(Gt313,gtu13,kmul(Gt323,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC AL CCTK_ATTRIBUTE_UNUSED = IfThen(LapseACoeff != 
      0,kneg(kmul(kmul(kpow(alphaL,-harmonicN),knmsub(kmadd(beta1L,PDupwindNthAnti1alpha,kmadd(beta2L,PDupwindNthAnti2alpha,kmadd(beta3L,PDupwindNthAnti3alpha,kmadd(PDupwindNthSymm1alpha,kfabs(beta1L),kmadd(PDupwindNthSymm2alpha,kfabs(beta2L),kmul(PDupwindNthSymm3alpha,kfabs(beta3L))))))),ToReal(LapseAdvectionCoeff),dtalpL)),ToReal(ScalarINV(harmonicF)))),ToReal(0));
    
    CCTK_REAL_VEC theta CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC B1L CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC B2L CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC B3L CCTK_ATTRIBUTE_UNUSED;
    
    if (ShiftBCoeff*ShiftGammaCoeff != 0)
    {
      B1L = 
        kdiv(knmsub(kmadd(beta1L,PDupwindNthAnti1beta1,kmadd(beta2L,PDupwindNthAnti2beta1,kmadd(beta3L,PDupwindNthAnti3beta1,kmadd(PDupwindNthSymm1beta1,kfabs(beta1L),kmadd(PDupwindNthSymm2beta1,kfabs(beta2L),kmul(PDupwindNthSymm3beta1,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),dtbetaxL),kmul(ToReal(ShiftGammaCoeff),theta));
      
      B2L = 
        kdiv(knmsub(kmadd(beta1L,PDupwindNthAnti1beta2,kmadd(beta2L,PDupwindNthAnti2beta2,kmadd(beta3L,PDupwindNthAnti3beta2,kmadd(PDupwindNthSymm1beta2,kfabs(beta1L),kmadd(PDupwindNthSymm2beta2,kfabs(beta2L),kmul(PDupwindNthSymm3beta2,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),dtbetayL),kmul(ToReal(ShiftGammaCoeff),theta));
      
      B3L = 
        kdiv(knmsub(kmadd(beta1L,PDupwindNthAnti1beta3,kmadd(beta2L,PDupwindNthAnti2beta3,kmadd(beta3L,PDupwindNthAnti3beta3,kmadd(PDupwindNthSymm1beta3,kfabs(beta1L),kmadd(PDupwindNthSymm2beta3,kfabs(beta2L),kmul(PDupwindNthSymm3beta3,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),dtbetazL),kmul(ToReal(ShiftGammaCoeff),theta));
    }
    else
    {
      B1L = ToReal(0);
      
      B2L = ToReal(0);
      
      B3L = ToReal(0);
    }
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(A[index],AL);
    vec_store_nta_partial(B1[index],B1L);
    vec_store_nta_partial(B2[index],B2L);
    vec_store_nta_partial(B3[index],B3L);
    vec_store_nta_partial(Xt1[index],Xt1L);
    vec_store_nta_partial(Xt2[index],Xt2L);
    vec_store_nta_partial(Xt3[index],Xt3L);
  }
  CCTK_ENDLOOP3STR(ML_BSSN_Host_convertFromADMBaseGamma);
}

extern "C" void ML_BSSN_Host_convertFromADMBaseGamma(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_Host_convertFromADMBaseGamma_Body");
  }
  
  if (cctk_iteration % ML_BSSN_Host_convertFromADMBaseGamma_calc_every != ML_BSSN_Host_convertFromADMBaseGamma_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::dtlapse",
    "ADMBase::dtshift",
    "ML_BSSN_Host::ML_dtlapse",
    "ML_BSSN_Host::ML_dtshift",
    "ML_BSSN_Host::ML_Gamma",
    "ML_BSSN_Host::ML_lapse",
    "ML_BSSN_Host::ML_metric",
    "ML_BSSN_Host::ML_shift"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_convertFromADMBaseGamma", 8, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Host_convertFromADMBaseGamma", 5, 5, 5);
  
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_Host_convertFromADMBaseGamma_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_convertFromADMBaseGamma_Body");
  }
}
