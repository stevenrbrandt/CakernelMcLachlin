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
  
  CCTK_INT ierr = 0;
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

static void ML_BSSN_Host_convertFromADMBaseGamma_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  LC_LOOP3VEC(ML_BSSN_Host_convertFromADMBaseGamma,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC dtalpL = vec_load(dtalp[index]);
    CCTK_REAL_VEC dtbetaxL = vec_load(dtbetax[index]);
    CCTK_REAL_VEC dtbetayL = vec_load(dtbetay[index]);
    CCTK_REAL_VEC dtbetazL = vec_load(dtbetaz[index]);
    CCTK_REAL_VEC gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L = vec_load(gt33[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC const PDupwindNthAnti1alpha = PDupwindNthAnti1(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1alpha = PDupwindNthSymm1(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2alpha = PDupwindNthAnti2(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2alpha = PDupwindNthSymm2(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3alpha = PDupwindNthAnti3(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3alpha = PDupwindNthSymm3(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1beta1 = PDupwindNthAnti1(&beta1[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1beta1 = PDupwindNthSymm1(&beta1[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2beta1 = PDupwindNthAnti2(&beta1[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2beta1 = PDupwindNthSymm2(&beta1[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3beta1 = PDupwindNthAnti3(&beta1[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3beta1 = PDupwindNthSymm3(&beta1[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1beta2 = PDupwindNthAnti1(&beta2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1beta2 = PDupwindNthSymm1(&beta2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2beta2 = PDupwindNthAnti2(&beta2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2beta2 = PDupwindNthSymm2(&beta2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3beta2 = PDupwindNthAnti3(&beta2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3beta2 = PDupwindNthSymm3(&beta2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1beta3 = PDupwindNthAnti1(&beta3[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1beta3 = PDupwindNthSymm1(&beta3[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2beta3 = PDupwindNthAnti2(&beta3[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2beta3 = PDupwindNthSymm2(&beta3[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3beta3 = PDupwindNthAnti3(&beta3[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3beta3 = PDupwindNthSymm3(&beta3[index]);
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
    
    CCTK_REAL_VEC Gt111 = 
      kmul(ToReal(0.5),kmadd(gtu11,PDstandardNth1gt11,knmsub(gtu12,PDstandardNth2gt11,kmsub(kmadd(gtu12,PDstandardNth1gt12,kmul(gtu13,PDstandardNth1gt13)),ToReal(2),kmul(gtu13,PDstandardNth3gt11)))));
    
    CCTK_REAL_VEC Gt211 = 
      kmul(ToReal(0.5),kmadd(gtu12,PDstandardNth1gt11,knmsub(gtu22,PDstandardNth2gt11,kmsub(kmadd(gtu22,PDstandardNth1gt12,kmul(gtu23,PDstandardNth1gt13)),ToReal(2),kmul(gtu23,PDstandardNth3gt11)))));
    
    CCTK_REAL_VEC Gt311 = 
      kmul(ToReal(0.5),kmadd(gtu13,PDstandardNth1gt11,knmsub(gtu23,PDstandardNth2gt11,kmsub(kmadd(gtu23,PDstandardNth1gt12,kmul(gtu33,PDstandardNth1gt13)),ToReal(2),kmul(gtu33,PDstandardNth3gt11)))));
    
    CCTK_REAL_VEC Gt112 = 
      kmul(kmadd(gtu12,PDstandardNth1gt22,kmadd(gtu11,PDstandardNth2gt11,kmul(gtu13,kadd(PDstandardNth1gt23,ksub(PDstandardNth2gt13,PDstandardNth3gt12))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt212 = 
      kmul(kmadd(gtu22,PDstandardNth1gt22,kmadd(gtu12,PDstandardNth2gt11,kmul(gtu23,kadd(PDstandardNth1gt23,ksub(PDstandardNth2gt13,PDstandardNth3gt12))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt312 = 
      kmul(kmadd(gtu23,PDstandardNth1gt22,kmadd(gtu13,PDstandardNth2gt11,kmul(gtu33,kadd(PDstandardNth1gt23,ksub(PDstandardNth2gt13,PDstandardNth3gt12))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt113 = 
      kmul(kmadd(gtu13,PDstandardNth1gt33,kmadd(gtu11,PDstandardNth3gt11,kmul(gtu12,kadd(PDstandardNth1gt23,ksub(PDstandardNth3gt12,PDstandardNth2gt13))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt213 = 
      kmul(kmadd(gtu23,PDstandardNth1gt33,kmadd(gtu12,PDstandardNth3gt11,kmul(gtu22,kadd(PDstandardNth1gt23,ksub(PDstandardNth3gt12,PDstandardNth2gt13))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt313 = 
      kmul(kmadd(gtu33,PDstandardNth1gt33,kmadd(gtu13,PDstandardNth3gt11,kmul(gtu23,kadd(PDstandardNth1gt23,ksub(PDstandardNth3gt12,PDstandardNth2gt13))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt122 = 
      kmul(ToReal(0.5),kmadd(gtu12,PDstandardNth2gt22,kmadd(gtu11,kmsub(PDstandardNth2gt12,ToReal(2),PDstandardNth1gt22),kmul(gtu13,kmsub(PDstandardNth2gt23,ToReal(2),PDstandardNth3gt22)))));
    
    CCTK_REAL_VEC Gt222 = 
      kmul(ToReal(0.5),kmadd(gtu22,PDstandardNth2gt22,kmadd(gtu12,kmsub(PDstandardNth2gt12,ToReal(2),PDstandardNth1gt22),kmul(gtu23,kmsub(PDstandardNth2gt23,ToReal(2),PDstandardNth3gt22)))));
    
    CCTK_REAL_VEC Gt322 = 
      kmul(ToReal(0.5),kmadd(gtu23,PDstandardNth2gt22,kmadd(gtu13,kmsub(PDstandardNth2gt12,ToReal(2),PDstandardNth1gt22),kmul(gtu33,kmsub(PDstandardNth2gt23,ToReal(2),PDstandardNth3gt22)))));
    
    CCTK_REAL_VEC Gt123 = 
      kmul(kmadd(gtu13,PDstandardNth2gt33,kmadd(gtu12,PDstandardNth3gt22,kmul(gtu11,kadd(PDstandardNth2gt13,ksub(PDstandardNth3gt12,PDstandardNth1gt23))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt223 = 
      kmul(kmadd(gtu23,PDstandardNth2gt33,kmadd(gtu22,PDstandardNth3gt22,kmul(gtu12,kadd(PDstandardNth2gt13,ksub(PDstandardNth3gt12,PDstandardNth1gt23))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt323 = 
      kmul(kmadd(gtu33,PDstandardNth2gt33,kmadd(gtu23,PDstandardNth3gt22,kmul(gtu13,kadd(PDstandardNth2gt13,ksub(PDstandardNth3gt12,PDstandardNth1gt23))))),ToReal(0.5));
    
    CCTK_REAL_VEC Gt133 = 
      kmul(ToReal(0.5),kmadd(gtu13,PDstandardNth3gt33,kmadd(gtu11,kmsub(PDstandardNth3gt13,ToReal(2),PDstandardNth1gt33),kmul(gtu12,kmsub(PDstandardNth3gt23,ToReal(2),PDstandardNth2gt33)))));
    
    CCTK_REAL_VEC Gt233 = 
      kmul(ToReal(0.5),kmadd(gtu23,PDstandardNth3gt33,kmadd(gtu12,kmsub(PDstandardNth3gt13,ToReal(2),PDstandardNth1gt33),kmul(gtu22,kmsub(PDstandardNth3gt23,ToReal(2),PDstandardNth2gt33)))));
    
    CCTK_REAL_VEC Gt333 = 
      kmul(ToReal(0.5),kmadd(gtu33,PDstandardNth3gt33,kmadd(gtu13,kmsub(PDstandardNth3gt13,ToReal(2),PDstandardNth1gt33),kmul(gtu23,kmsub(PDstandardNth3gt23,ToReal(2),PDstandardNth2gt33)))));
    
    CCTK_REAL_VEC Xt1L = 
      kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmadd(Gt133,gtu33,kmul(kmadd(Gt112,gtu12,kmadd(Gt113,gtu13,kmul(Gt123,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Xt2L = 
      kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmadd(Gt233,gtu33,kmul(kmadd(Gt212,gtu12,kmadd(Gt213,gtu13,kmul(Gt223,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Xt3L = 
      kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmadd(Gt333,gtu33,kmul(kmadd(Gt312,gtu12,kmadd(Gt313,gtu13,kmul(Gt323,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC AL = IfThen(LapseACoeff != 
      0,kneg(kmul(kmul(kpow(alphaL,-harmonicN),knmsub(kmadd(beta1L,PDupwindNthAnti1alpha,kmadd(beta2L,PDupwindNthAnti2alpha,kmadd(beta3L,PDupwindNthAnti3alpha,kmadd(PDupwindNthSymm1alpha,kfabs(beta1L),kmadd(PDupwindNthSymm2alpha,kfabs(beta2L),kmul(PDupwindNthSymm3alpha,kfabs(beta3L))))))),ToReal(LapseAdvectionCoeff),dtalpL)),ToReal(ScalarINV(harmonicF)))),ToReal(0));
    
    CCTK_REAL_VEC theta = ToReal(1);
    
    CCTK_REAL_VEC B1L;
    CCTK_REAL_VEC B2L;
    CCTK_REAL_VEC B3L;
    
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
    vec_store_partial_prepare(i,lc_imin,lc_imax);
    vec_store_nta_partial(A[index],AL);
    vec_store_nta_partial(B1[index],B1L);
    vec_store_nta_partial(B2[index],B2L);
    vec_store_nta_partial(B3[index],B3L);
    vec_store_nta_partial(Xt1[index],Xt1L);
    vec_store_nta_partial(Xt2[index],Xt2L);
    vec_store_nta_partial(Xt3[index],Xt3L);
  }
  LC_ENDLOOP3VEC(ML_BSSN_Host_convertFromADMBaseGamma);
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
  
  const char *const groups[] = {
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
