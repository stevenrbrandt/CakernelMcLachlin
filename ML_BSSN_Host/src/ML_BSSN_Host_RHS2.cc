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
  
  if (cctk_iteration % ML_BSSN_Host_RHS2_calc_every != ML_BSSN_Host_RHS2_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_curvrhs.");
  return;
}

static void ML_BSSN_Host_RHS2_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3STR(ML_BSSN_Host_RHS2,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    // vec_iter_counter+=CCTK_REAL_VEC_SIZE;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = vec_load(alpha[index]);
    CCTK_REAL_VEC At11L CCTK_ATTRIBUTE_UNUSED = vec_load(At11[index]);
    CCTK_REAL_VEC At12L CCTK_ATTRIBUTE_UNUSED = vec_load(At12[index]);
    CCTK_REAL_VEC At13L CCTK_ATTRIBUTE_UNUSED = vec_load(At13[index]);
    CCTK_REAL_VEC At22L CCTK_ATTRIBUTE_UNUSED = vec_load(At22[index]);
    CCTK_REAL_VEC At23L CCTK_ATTRIBUTE_UNUSED = vec_load(At23[index]);
    CCTK_REAL_VEC At33L CCTK_ATTRIBUTE_UNUSED = vec_load(At33[index]);
    CCTK_REAL_VEC beta1L CCTK_ATTRIBUTE_UNUSED = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L CCTK_ATTRIBUTE_UNUSED = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L CCTK_ATTRIBUTE_UNUSED = vec_load(beta3[index]);
    CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = vec_load(gt33[index]);
    CCTK_REAL_VEC phiL CCTK_ATTRIBUTE_UNUSED = vec_load(phi[index]);
    CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = vec_load(trK[index]);
    CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt3[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    const CCTK_REAL_VEC PDstandardNth1alpha CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&alpha[index]);
    const CCTK_REAL_VEC PDstandardNth2alpha CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&alpha[index]);
    const CCTK_REAL_VEC PDstandardNth3alpha CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&alpha[index]);
    const CCTK_REAL_VEC PDstandardNth11alpha CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&alpha[index]);
    const CCTK_REAL_VEC PDstandardNth22alpha CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&alpha[index]);
    const CCTK_REAL_VEC PDstandardNth33alpha CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&alpha[index]);
    const CCTK_REAL_VEC PDstandardNth12alpha CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&alpha[index]);
    const CCTK_REAL_VEC PDstandardNth13alpha CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&alpha[index]);
    const CCTK_REAL_VEC PDstandardNth23alpha CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&alpha[index]);
    const CCTK_REAL_VEC PDstandardNth1beta1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&beta1[index]);
    const CCTK_REAL_VEC PDstandardNth2beta1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&beta1[index]);
    const CCTK_REAL_VEC PDstandardNth3beta1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&beta1[index]);
    const CCTK_REAL_VEC PDstandardNth1beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth2beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth3beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth1beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth2beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth3beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth11gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth22gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth33gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth12gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth13gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth23gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt11[index]);
    const CCTK_REAL_VEC PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth11gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth22gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth33gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth12gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth13gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth23gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt12[index]);
    const CCTK_REAL_VEC PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth11gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth22gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth33gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth12gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth13gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth23gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt13[index]);
    const CCTK_REAL_VEC PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth11gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth22gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth33gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth12gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth13gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth23gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt22[index]);
    const CCTK_REAL_VEC PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth11gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth22gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth33gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth12gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth13gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth23gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt23[index]);
    const CCTK_REAL_VEC PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt33[index]);
    const CCTK_REAL_VEC PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt33[index]);
    const CCTK_REAL_VEC PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt33[index]);
    const CCTK_REAL_VEC PDstandardNth11gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt33[index]);
    const CCTK_REAL_VEC PDstandardNth22gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt33[index]);
    const CCTK_REAL_VEC PDstandardNth33gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt33[index]);
    const CCTK_REAL_VEC PDstandardNth12gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt33[index]);
    const CCTK_REAL_VEC PDstandardNth13gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt33[index]);
    const CCTK_REAL_VEC PDstandardNth23gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt33[index]);
    const CCTK_REAL_VEC PDstandardNth1phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth2phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth3phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth11phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth22phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth33phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth12phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth13phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth23phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth1Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&Xt1[index]);
    const CCTK_REAL_VEC PDstandardNth2Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&Xt1[index]);
    const CCTK_REAL_VEC PDstandardNth3Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&Xt1[index]);
    const CCTK_REAL_VEC PDstandardNth1Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&Xt2[index]);
    const CCTK_REAL_VEC PDstandardNth2Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&Xt2[index]);
    const CCTK_REAL_VEC PDstandardNth3Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&Xt2[index]);
    const CCTK_REAL_VEC PDstandardNth1Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&Xt3[index]);
    const CCTK_REAL_VEC PDstandardNth2Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&Xt3[index]);
    const CCTK_REAL_VEC PDstandardNth3Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&Xt3[index]);
    
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
    
    CCTK_REAL_VEC Gtl111 CCTK_ATTRIBUTE_UNUSED = 
      kmul(PDstandardNth1gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl112 CCTK_ATTRIBUTE_UNUSED = 
      kmul(PDstandardNth2gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl113 CCTK_ATTRIBUTE_UNUSED = 
      kmul(PDstandardNth3gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl122 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(PDstandardNth1gt22,ToReal(-0.5),PDstandardNth2gt12);
    
    CCTK_REAL_VEC Gtl123 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kadd(PDstandardNth2gt13,ksub(PDstandardNth3gt12,PDstandardNth1gt23)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl133 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(PDstandardNth1gt33,ToReal(-0.5),PDstandardNth3gt13);
    
    CCTK_REAL_VEC Gtl211 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(PDstandardNth2gt11,ToReal(-0.5),PDstandardNth1gt12);
    
    CCTK_REAL_VEC Gtl212 CCTK_ATTRIBUTE_UNUSED = 
      kmul(PDstandardNth1gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl213 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kadd(PDstandardNth1gt23,ksub(PDstandardNth3gt12,PDstandardNth2gt13)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl222 CCTK_ATTRIBUTE_UNUSED = 
      kmul(PDstandardNth2gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl223 CCTK_ATTRIBUTE_UNUSED = 
      kmul(PDstandardNth3gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl233 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(PDstandardNth2gt33,ToReal(-0.5),PDstandardNth3gt23);
    
    CCTK_REAL_VEC Gtl311 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(PDstandardNth3gt11,ToReal(-0.5),PDstandardNth1gt13);
    
    CCTK_REAL_VEC Gtl312 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kadd(PDstandardNth1gt23,ksub(PDstandardNth2gt13,PDstandardNth3gt12)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl313 CCTK_ATTRIBUTE_UNUSED = 
      kmul(PDstandardNth1gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl322 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(PDstandardNth3gt22,ToReal(-0.5),PDstandardNth2gt23);
    
    CCTK_REAL_VEC Gtl323 CCTK_ATTRIBUTE_UNUSED = 
      kmul(PDstandardNth2gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl333 CCTK_ATTRIBUTE_UNUSED = 
      kmul(PDstandardNth3gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtlu111 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu11,kmadd(Gtl112,gtu12,kmul(Gtl113,gtu13)));
    
    CCTK_REAL_VEC Gtlu112 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu12,kmadd(Gtl112,gtu22,kmul(Gtl113,gtu23)));
    
    CCTK_REAL_VEC Gtlu113 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu13,kmadd(Gtl112,gtu23,kmul(Gtl113,gtu33)));
    
    CCTK_REAL_VEC Gtlu121 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu11,kmadd(Gtl122,gtu12,kmul(Gtl123,gtu13)));
    
    CCTK_REAL_VEC Gtlu122 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu12,kmadd(Gtl122,gtu22,kmul(Gtl123,gtu23)));
    
    CCTK_REAL_VEC Gtlu123 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu13,kmadd(Gtl122,gtu23,kmul(Gtl123,gtu33)));
    
    CCTK_REAL_VEC Gtlu131 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu11,kmadd(Gtl123,gtu12,kmul(Gtl133,gtu13)));
    
    CCTK_REAL_VEC Gtlu132 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu12,kmadd(Gtl123,gtu22,kmul(Gtl133,gtu23)));
    
    CCTK_REAL_VEC Gtlu133 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu13,kmadd(Gtl123,gtu23,kmul(Gtl133,gtu33)));
    
    CCTK_REAL_VEC Gtlu211 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl211,gtu11,kmadd(Gtl212,gtu12,kmul(Gtl213,gtu13)));
    
    CCTK_REAL_VEC Gtlu212 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl211,gtu12,kmadd(Gtl212,gtu22,kmul(Gtl213,gtu23)));
    
    CCTK_REAL_VEC Gtlu213 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl211,gtu13,kmadd(Gtl212,gtu23,kmul(Gtl213,gtu33)));
    
    CCTK_REAL_VEC Gtlu221 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl212,gtu11,kmadd(Gtl222,gtu12,kmul(Gtl223,gtu13)));
    
    CCTK_REAL_VEC Gtlu222 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl212,gtu12,kmadd(Gtl222,gtu22,kmul(Gtl223,gtu23)));
    
    CCTK_REAL_VEC Gtlu223 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl212,gtu13,kmadd(Gtl222,gtu23,kmul(Gtl223,gtu33)));
    
    CCTK_REAL_VEC Gtlu231 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl213,gtu11,kmadd(Gtl223,gtu12,kmul(Gtl233,gtu13)));
    
    CCTK_REAL_VEC Gtlu232 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl213,gtu12,kmadd(Gtl223,gtu22,kmul(Gtl233,gtu23)));
    
    CCTK_REAL_VEC Gtlu233 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl213,gtu13,kmadd(Gtl223,gtu23,kmul(Gtl233,gtu33)));
    
    CCTK_REAL_VEC Gtlu311 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl311,gtu11,kmadd(Gtl312,gtu12,kmul(Gtl313,gtu13)));
    
    CCTK_REAL_VEC Gtlu312 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl311,gtu12,kmadd(Gtl312,gtu22,kmul(Gtl313,gtu23)));
    
    CCTK_REAL_VEC Gtlu313 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl311,gtu13,kmadd(Gtl312,gtu23,kmul(Gtl313,gtu33)));
    
    CCTK_REAL_VEC Gtlu321 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl312,gtu11,kmadd(Gtl322,gtu12,kmul(Gtl323,gtu13)));
    
    CCTK_REAL_VEC Gtlu322 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl312,gtu12,kmadd(Gtl322,gtu22,kmul(Gtl323,gtu23)));
    
    CCTK_REAL_VEC Gtlu323 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl312,gtu13,kmadd(Gtl322,gtu23,kmul(Gtl323,gtu33)));
    
    CCTK_REAL_VEC Gtlu331 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl313,gtu11,kmadd(Gtl323,gtu12,kmul(Gtl333,gtu13)));
    
    CCTK_REAL_VEC Gtlu332 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl313,gtu12,kmadd(Gtl323,gtu22,kmul(Gtl333,gtu23)));
    
    CCTK_REAL_VEC Gtlu333 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl313,gtu13,kmadd(Gtl323,gtu23,kmul(Gtl333,gtu33)));
    
    CCTK_REAL_VEC Gt111 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu11,kmadd(Gtl211,gtu12,kmul(Gtl311,gtu13)));
    
    CCTK_REAL_VEC Gt211 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu12,kmadd(Gtl211,gtu22,kmul(Gtl311,gtu23)));
    
    CCTK_REAL_VEC Gt311 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu13,kmadd(Gtl211,gtu23,kmul(Gtl311,gtu33)));
    
    CCTK_REAL_VEC Gt112 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu11,kmadd(Gtl212,gtu12,kmul(Gtl312,gtu13)));
    
    CCTK_REAL_VEC Gt212 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu12,kmadd(Gtl212,gtu22,kmul(Gtl312,gtu23)));
    
    CCTK_REAL_VEC Gt312 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu13,kmadd(Gtl212,gtu23,kmul(Gtl312,gtu33)));
    
    CCTK_REAL_VEC Gt113 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu11,kmadd(Gtl213,gtu12,kmul(Gtl313,gtu13)));
    
    CCTK_REAL_VEC Gt213 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu12,kmadd(Gtl213,gtu22,kmul(Gtl313,gtu23)));
    
    CCTK_REAL_VEC Gt313 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu13,kmadd(Gtl213,gtu23,kmul(Gtl313,gtu33)));
    
    CCTK_REAL_VEC Gt122 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl122,gtu11,kmadd(Gtl222,gtu12,kmul(Gtl322,gtu13)));
    
    CCTK_REAL_VEC Gt222 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl122,gtu12,kmadd(Gtl222,gtu22,kmul(Gtl322,gtu23)));
    
    CCTK_REAL_VEC Gt322 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl122,gtu13,kmadd(Gtl222,gtu23,kmul(Gtl322,gtu33)));
    
    CCTK_REAL_VEC Gt123 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl123,gtu11,kmadd(Gtl223,gtu12,kmul(Gtl323,gtu13)));
    
    CCTK_REAL_VEC Gt223 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl123,gtu12,kmadd(Gtl223,gtu22,kmul(Gtl323,gtu23)));
    
    CCTK_REAL_VEC Gt323 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl123,gtu13,kmadd(Gtl223,gtu23,kmul(Gtl323,gtu33)));
    
    CCTK_REAL_VEC Gt133 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl133,gtu11,kmadd(Gtl233,gtu12,kmul(Gtl333,gtu13)));
    
    CCTK_REAL_VEC Gt233 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl133,gtu12,kmadd(Gtl233,gtu22,kmul(Gtl333,gtu23)));
    
    CCTK_REAL_VEC Gt333 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl133,gtu13,kmadd(Gtl233,gtu23,kmul(Gtl333,gtu33)));
    
    CCTK_REAL_VEC Xtn1 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmadd(Gt133,gtu33,kmul(kmadd(Gt112,gtu12,kmadd(Gt113,gtu13,kmul(Gt123,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Xtn2 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmadd(Gt233,gtu33,kmul(kmadd(Gt212,gtu12,kmadd(Gt213,gtu13,kmul(Gt223,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Xtn3 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmadd(Gt333,gtu33,kmul(kmadd(Gt312,gtu12,kmadd(Gt313,gtu13,kmul(Gt323,gtu23))),ToReal(2)))));
    
    CCTK_REAL_VEC Rt11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt211,Gtlu211,kmadd(Gt212,Gtlu212,kmadd(Gt213,Gtlu213,kmadd(Gt311,Gtlu311,kmadd(Gt312,Gtlu312,kmadd(Gt313,Gtlu313,kmadd(gt11L,PDstandardNth1Xt1,kmadd(gt12L,PDstandardNth1Xt2,kmadd(gt13L,PDstandardNth1Xt3,kmadd(Gtl111,Xtn1,kmadd(Gtl112,Xtn2,kmadd(Gtl113,Xtn3,kmadd(kmadd(gtu12,kmul(PDstandardNth12gt11,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt11,ToReal(-2)),kmsub(gtu23,kmul(PDstandardNth23gt11,ToReal(-2)),kmadd(gtu11,PDstandardNth11gt11,kmadd(gtu33,PDstandardNth33gt11,kmul(gtu22,PDstandardNth22gt11)))))),ToReal(0.5),kmadd(kmadd(Gt211,Gtlu121,kmadd(Gt212,Gtlu122,kmadd(Gt213,Gtlu123,kmadd(Gt311,Gtlu131,kmadd(Gt312,Gtlu132,kmul(Gt313,Gtlu133)))))),ToReal(2),kmul(kmadd(Gt111,Gtlu111,kmadd(Gt112,Gtlu112,kmul(Gt113,Gtlu113))),ToReal(3))))))))))))))));
    
    CCTK_REAL_VEC Rt12 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gt12L,PDstandardNth1Xt1,kmadd(gt22L,PDstandardNth1Xt2,kmadd(gt23L,PDstandardNth1Xt3,kmadd(gt11L,PDstandardNth2Xt1,kmadd(gt12L,PDstandardNth2Xt2,kmadd(gt13L,PDstandardNth2Xt3,kmadd(Gtl112,Xtn1,kmadd(Gtl211,Xtn1,kmadd(Gtl122,Xtn2,kmadd(Gtl212,Xtn2,kmadd(Gtl123,Xtn3,kmadd(Gtl213,Xtn3,kmadd(gtu12,kmul(PDstandardNth12gt12,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt12,ToReal(-2)),kmadd(gtu23,kmul(PDstandardNth23gt12,ToReal(-2)),knmsub(gtu11,PDstandardNth11gt12,knmsub(gtu22,PDstandardNth22gt12,knmsub(gtu33,PDstandardNth33gt12,kmadd(kmadd(Gt112,Gtlu111,kmadd(Gt122,Gtlu112,kmadd(Gt123,Gtlu113,kmadd(Gt111,Gtlu121,kmadd(Gt212,Gtlu121,kmadd(Gt112,Gtlu122,kmadd(Gt222,Gtlu122,kmadd(Gt113,Gtlu123,kmadd(Gt223,Gtlu123,kmadd(Gt312,Gtlu131,kmadd(Gt322,Gtlu132,kmadd(Gt323,Gtlu133,kmadd(Gt111,Gtlu211,kmadd(Gt112,Gtlu212,kmadd(Gt113,Gtlu213,kmadd(Gt311,Gtlu231,kmadd(Gt312,Gtlu232,kmadd(Gt313,Gtlu233,kmadd(Gt311,Gtlu321,kmadd(Gt312,Gtlu322,kmul(Gt313,Gtlu323))))))))))))))))))))),ToReal(2),kmul(kmadd(Gt211,Gtlu221,kmadd(Gt212,Gtlu222,kmul(Gt213,Gtlu223))),ToReal(4))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt13 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gt13L,PDstandardNth1Xt1,kmadd(gt23L,PDstandardNth1Xt2,kmadd(gt33L,PDstandardNth1Xt3,kmadd(gt11L,PDstandardNth3Xt1,kmadd(gt12L,PDstandardNth3Xt2,kmadd(gt13L,PDstandardNth3Xt3,kmadd(Gtl113,Xtn1,kmadd(Gtl311,Xtn1,kmadd(Gtl123,Xtn2,kmadd(Gtl312,Xtn2,kmadd(Gtl133,Xtn3,kmadd(Gtl313,Xtn3,kmadd(gtu12,kmul(PDstandardNth12gt13,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt13,ToReal(-2)),kmadd(gtu23,kmul(PDstandardNth23gt13,ToReal(-2)),knmsub(gtu11,PDstandardNth11gt13,knmsub(gtu22,PDstandardNth22gt13,knmsub(gtu33,PDstandardNth33gt13,kmadd(kmadd(Gt113,Gtlu111,kmadd(Gt123,Gtlu112,kmadd(Gt133,Gtlu113,kmadd(Gt213,Gtlu121,kmadd(Gt223,Gtlu122,kmadd(Gt233,Gtlu123,kmadd(Gt111,Gtlu131,kmadd(Gt313,Gtlu131,kmadd(Gt112,Gtlu132,kmadd(Gt323,Gtlu132,kmadd(Gt113,Gtlu133,kmadd(Gt333,Gtlu133,kmadd(Gt211,Gtlu231,kmadd(Gt212,Gtlu232,kmadd(Gt213,Gtlu233,kmadd(Gt111,Gtlu311,kmadd(Gt112,Gtlu312,kmadd(Gt113,Gtlu313,kmadd(Gt211,Gtlu321,kmadd(Gt212,Gtlu322,kmul(Gt213,Gtlu323))))))))))))))))))))),ToReal(2),kmul(kmadd(Gt311,Gtlu331,kmadd(Gt312,Gtlu332,kmul(Gt313,Gtlu333))),ToReal(4))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt312,Gtlu321,kmadd(Gt322,Gtlu322,kmadd(Gt323,Gtlu323,kmadd(gt12L,PDstandardNth2Xt1,kmadd(gt22L,PDstandardNth2Xt2,kmadd(gt23L,PDstandardNth2Xt3,kmadd(Gtl212,Xtn1,kmadd(Gtl222,Xtn2,kmadd(Gtl223,Xtn3,kmadd(kmadd(gtu12,kmul(PDstandardNth12gt22,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt22,ToReal(-2)),kmsub(gtu23,kmul(PDstandardNth23gt22,ToReal(-2)),kmadd(gtu11,PDstandardNth11gt22,kmadd(gtu33,PDstandardNth33gt22,kmul(gtu22,PDstandardNth22gt22)))))),ToReal(0.5),kmadd(kmadd(Gt312,Gtlu231,kmadd(Gt322,Gtlu232,kmul(Gt323,Gtlu233))),ToReal(2),kmadd(Gt112,kmadd(Gtlu211,ToReal(2),Gtlu121),kmadd(Gt122,kmadd(Gtlu212,ToReal(2),Gtlu122),kmadd(Gt123,kmadd(Gtlu213,ToReal(2),Gtlu123),kmul(kmadd(Gt212,Gtlu221,kmadd(Gt222,Gtlu222,kmul(Gt223,Gtlu223))),ToReal(3))))))))))))))));
    
    CCTK_REAL_VEC Rt23 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmadd(gt13L,PDstandardNth2Xt1,kmadd(gt23L,PDstandardNth2Xt2,kmadd(gt33L,PDstandardNth2Xt3,kmadd(gt12L,PDstandardNth3Xt1,kmadd(gt22L,PDstandardNth3Xt2,kmadd(gt23L,PDstandardNth3Xt3,kmadd(Gtl213,Xtn1,kmadd(Gtl312,Xtn1,kmadd(Gtl223,Xtn2,kmadd(Gtl322,Xtn2,kmadd(Gtl233,Xtn3,kmadd(Gtl323,Xtn3,kmadd(gtu12,kmul(PDstandardNth12gt23,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt23,ToReal(-2)),kmadd(gtu23,kmul(PDstandardNth23gt23,ToReal(-2)),knmsub(gtu11,PDstandardNth11gt23,knmsub(gtu22,PDstandardNth22gt23,knmsub(gtu33,PDstandardNth33gt23,kmadd(kmadd(Gt112,Gtlu131,kmadd(Gt122,Gtlu132,kmadd(Gt123,Gtlu133,kmadd(Gt113,Gtlu211,kmadd(Gt123,Gtlu212,kmadd(Gt133,Gtlu213,kmadd(Gt213,Gtlu221,kmadd(Gt223,Gtlu222,kmadd(Gt233,Gtlu223,kmadd(Gt212,Gtlu231,kmadd(Gt313,Gtlu231,kmadd(Gt222,Gtlu232,kmadd(Gt323,Gtlu232,kmadd(Gt223,Gtlu233,kmadd(Gt333,Gtlu233,kmadd(Gt112,Gtlu311,kmadd(Gt122,Gtlu312,kmadd(Gt123,Gtlu313,kmadd(Gt212,Gtlu321,kmadd(Gt222,Gtlu322,kmul(Gt223,Gtlu323))))))))))))))))))))),ToReal(2),kmul(kmadd(Gt312,Gtlu331,kmadd(Gt322,Gtlu332,kmul(Gt323,Gtlu333))),ToReal(4))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt13L,PDstandardNth3Xt1,kmadd(gt23L,PDstandardNth3Xt2,kmadd(gt33L,PDstandardNth3Xt3,kmadd(Gtl313,Xtn1,kmadd(Gtl323,Xtn2,kmadd(Gtl333,Xtn3,kmadd(kmadd(gtu12,kmul(PDstandardNth12gt33,ToReal(-2)),kmadd(gtu13,kmul(PDstandardNth13gt33,ToReal(-2)),kmsub(gtu23,kmul(PDstandardNth23gt33,ToReal(-2)),kmadd(gtu11,PDstandardNth11gt33,kmadd(gtu33,PDstandardNth33gt33,kmul(gtu22,PDstandardNth22gt33)))))),ToReal(0.5),kmadd(Gt113,kmadd(Gtlu311,ToReal(2),Gtlu131),kmadd(Gt123,kmadd(Gtlu312,ToReal(2),Gtlu132),kmadd(Gt133,kmadd(Gtlu313,ToReal(2),Gtlu133),kmadd(Gt213,kmadd(Gtlu321,ToReal(2),Gtlu231),kmadd(Gt223,kmadd(Gtlu322,ToReal(2),Gtlu232),kmadd(Gt233,kmadd(Gtlu323,ToReal(2),Gtlu233),kmul(kmadd(Gt313,Gtlu331,kmadd(Gt323,Gtlu332,kmul(Gt333,Gtlu333))),ToReal(3)))))))))))))));
    
    CCTK_REAL_VEC fac1 CCTK_ATTRIBUTE_UNUSED = 
      IfThen(conformalMethod,kdiv(ToReal(-0.5),phiL),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,PDstandardNth1phi);
    
    CCTK_REAL_VEC cdphi2 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,PDstandardNth2phi);
    
    CCTK_REAL_VEC cdphi3 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,PDstandardNth3phi);
    
    CCTK_REAL_VEC fac2 CCTK_ATTRIBUTE_UNUSED = 
      IfThen(conformalMethod,kdiv(ToReal(0.5),kmul(phiL,phiL)),ToReal(0));
    
    CCTK_REAL_VEC cdphi211 CCTK_ATTRIBUTE_UNUSED = 
      kmsub(fac2,kmul(PDstandardNth1phi,PDstandardNth1phi),kmul(fac1,kmadd(Gt111,PDstandardNth1phi,kmadd(Gt211,PDstandardNth2phi,kmsub(Gt311,PDstandardNth3phi,PDstandardNth11phi)))));
    
    CCTK_REAL_VEC cdphi212 CCTK_ATTRIBUTE_UNUSED = 
      kmsub(fac2,kmul(PDstandardNth1phi,PDstandardNth2phi),kmul(fac1,kmadd(Gt112,PDstandardNth1phi,kmadd(Gt212,PDstandardNth2phi,kmsub(Gt312,PDstandardNth3phi,PDstandardNth12phi)))));
    
    CCTK_REAL_VEC cdphi213 CCTK_ATTRIBUTE_UNUSED = 
      kmsub(fac2,kmul(PDstandardNth1phi,PDstandardNth3phi),kmul(fac1,kmadd(Gt113,PDstandardNth1phi,kmadd(Gt213,PDstandardNth2phi,kmsub(Gt313,PDstandardNth3phi,PDstandardNth13phi)))));
    
    CCTK_REAL_VEC cdphi222 CCTK_ATTRIBUTE_UNUSED = 
      kmsub(fac2,kmul(PDstandardNth2phi,PDstandardNth2phi),kmul(fac1,kmadd(Gt122,PDstandardNth1phi,kmadd(Gt222,PDstandardNth2phi,kmsub(Gt322,PDstandardNth3phi,PDstandardNth22phi)))));
    
    CCTK_REAL_VEC cdphi223 CCTK_ATTRIBUTE_UNUSED = 
      kmsub(fac2,kmul(PDstandardNth2phi,PDstandardNth3phi),kmul(fac1,kmadd(Gt123,PDstandardNth1phi,kmadd(Gt223,PDstandardNth2phi,kmsub(Gt323,PDstandardNth3phi,PDstandardNth23phi)))));
    
    CCTK_REAL_VEC cdphi233 CCTK_ATTRIBUTE_UNUSED = 
      kmsub(fac2,kmul(PDstandardNth3phi,PDstandardNth3phi),kmul(fac1,kmadd(Gt133,PDstandardNth1phi,kmadd(Gt233,PDstandardNth2phi,kmsub(Gt333,PDstandardNth3phi,PDstandardNth33phi)))));
    
    CCTK_REAL_VEC Rphi11 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(-2),kadd(cdphi211,kmadd(kmul(cdphi1,cdphi1),kmul(kmadd(gt11L,gtu11,ToReal(-1)),ToReal(2)),kmul(gt11L,kmadd(cdphi211,gtu11,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu33,kmul(cdphi3,cdphi3))))),ToReal(2),kmadd(gtu22,kmadd(kmul(cdphi2,cdphi2),ToReal(2),cdphi222),kmul(kmadd(cdphi1,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13)),kmul(cdphi2,kmul(cdphi3,gtu23))),ToReal(4))))))))));
    
    CCTK_REAL_VEC Rphi12 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(-2),kadd(cdphi212,kmadd(gt12L,kmadd(cdphi211,gtu11,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu11,kmul(cdphi1,cdphi1))))),ToReal(2),kmadd(gtu22,kmadd(kmul(cdphi2,cdphi2),ToReal(2),cdphi222),kmadd(gtu33,kmadd(kmul(cdphi3,cdphi3),ToReal(2),cdphi233),kmul(cdphi2,kmul(cdphi3,kmul(gtu23,ToReal(4)))))))),kmul(cdphi1,kmadd(gt12L,kmul(cdphi3,kmul(gtu13,ToReal(4))),kmul(cdphi2,kmadd(gt12L,kmul(gtu12,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi13 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(-2),kadd(cdphi213,kmadd(gt13L,kmadd(cdphi211,gtu11,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu11,kmul(cdphi1,cdphi1))))),ToReal(2),kmadd(gtu22,kmadd(kmul(cdphi2,cdphi2),ToReal(2),cdphi222),kmadd(gtu33,kmadd(kmul(cdphi3,cdphi3),ToReal(2),cdphi233),kmul(cdphi2,kmul(cdphi3,kmul(gtu23,ToReal(4)))))))),kmul(cdphi1,kmadd(gt13L,kmul(cdphi2,kmul(gtu12,ToReal(4))),kmul(cdphi3,kmadd(gt13L,kmul(gtu13,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi22 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(-2),kadd(cdphi222,kmadd(kmul(cdphi2,cdphi2),kmul(kmadd(gt22L,gtu22,ToReal(-1)),ToReal(2)),kmul(gt22L,kmadd(cdphi222,gtu22,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu33,kmul(cdphi3,cdphi3))))),ToReal(2),kmadd(gtu11,kmadd(kmul(cdphi1,cdphi1),ToReal(2),cdphi211),kmul(kmadd(cdphi1,kmul(cdphi3,gtu13),kmul(cdphi2,kmadd(cdphi1,gtu12,kmul(cdphi3,gtu23)))),ToReal(4))))))))));
    
    CCTK_REAL_VEC Rphi23 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(-2),kadd(cdphi223,kmadd(gt23L,kmadd(cdphi222,gtu22,kmadd(kmadd(cdphi212,gtu12,kmadd(cdphi213,gtu13,kmadd(cdphi223,gtu23,kmul(gtu22,kmul(cdphi2,cdphi2))))),ToReal(2),kmadd(gtu11,kmadd(kmul(cdphi1,cdphi1),ToReal(2),cdphi211),kmadd(gtu33,kmadd(kmul(cdphi3,cdphi3),ToReal(2),cdphi233),kmul(cdphi1,kmul(cdphi3,kmul(gtu13,ToReal(4)))))))),kmul(cdphi2,kmadd(gt23L,kmul(cdphi1,kmul(gtu12,ToReal(4))),kmul(cdphi3,kmadd(gt23L,kmul(gtu23,ToReal(4)),ToReal(-2))))))));
    
    CCTK_REAL_VEC Rphi33 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(-2),kadd(cdphi233,kmadd(kmul(cdphi3,cdphi3),kmul(kmadd(gt33L,gtu33,ToReal(-1)),ToReal(2)),kmul(gt33L,kmadd(cdphi233,gtu33,kmadd(kmadd(cdphi213,gtu13,kmul(cdphi223,gtu23)),ToReal(2),kmadd(gtu11,kmadd(kmul(cdphi1,cdphi1),ToReal(2),cdphi211),kmadd(gtu22,kmadd(kmul(cdphi2,cdphi2),ToReal(2),cdphi222),kmadd(cdphi3,kmul(kmadd(cdphi1,gtu13,kmul(cdphi2,gtu23)),ToReal(4)),kmul(gtu12,kmadd(cdphi212,ToReal(2),kmul(cdphi1,kmul(cdphi2,ToReal(4))))))))))))));
    
    CCTK_REAL_VEC Atm11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At11L,gtu11,kmadd(At12L,gtu12,kmul(At13L,gtu13)));
    
    CCTK_REAL_VEC Atm21 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At11L,gtu12,kmadd(At12L,gtu22,kmul(At13L,gtu23)));
    
    CCTK_REAL_VEC Atm31 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At11L,gtu13,kmadd(At12L,gtu23,kmul(At13L,gtu33)));
    
    CCTK_REAL_VEC Atm12 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At12L,gtu11,kmadd(At22L,gtu12,kmul(At23L,gtu13)));
    
    CCTK_REAL_VEC Atm22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At12L,gtu12,kmadd(At22L,gtu22,kmul(At23L,gtu23)));
    
    CCTK_REAL_VEC Atm32 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At12L,gtu13,kmadd(At22L,gtu23,kmul(At23L,gtu33)));
    
    CCTK_REAL_VEC Atm13 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At13L,gtu11,kmadd(At23L,gtu12,kmul(At33L,gtu13)));
    
    CCTK_REAL_VEC Atm23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At13L,gtu12,kmadd(At23L,gtu22,kmul(At33L,gtu23)));
    
    CCTK_REAL_VEC Atm33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At13L,gtu13,kmadd(At23L,gtu23,kmul(At33L,gtu33)));
    
    CCTK_REAL_VEC e4phi CCTK_ATTRIBUTE_UNUSED = 
      IfThen(conformalMethod,kdiv(ToReal(1),kmul(phiL,phiL)),kexp(kmul(phiL,ToReal(4))));
    
    CCTK_REAL_VEC em4phi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),e4phi);
    
    CCTK_REAL_VEC g11 CCTK_ATTRIBUTE_UNUSED = kmul(gt11L,e4phi);
    
    CCTK_REAL_VEC g12 CCTK_ATTRIBUTE_UNUSED = kmul(gt12L,e4phi);
    
    CCTK_REAL_VEC g13 CCTK_ATTRIBUTE_UNUSED = kmul(gt13L,e4phi);
    
    CCTK_REAL_VEC g22 CCTK_ATTRIBUTE_UNUSED = kmul(gt22L,e4phi);
    
    CCTK_REAL_VEC g23 CCTK_ATTRIBUTE_UNUSED = kmul(gt23L,e4phi);
    
    CCTK_REAL_VEC g33 CCTK_ATTRIBUTE_UNUSED = kmul(gt33L,e4phi);
    
    CCTK_REAL_VEC gu11 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu11);
    
    CCTK_REAL_VEC gu12 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu12);
    
    CCTK_REAL_VEC gu13 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu13);
    
    CCTK_REAL_VEC gu22 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu22);
    
    CCTK_REAL_VEC gu23 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu23);
    
    CCTK_REAL_VEC gu33 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu33);
    
    CCTK_REAL_VEC R11 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi11,Rt11);
    
    CCTK_REAL_VEC R12 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi12,Rt12);
    
    CCTK_REAL_VEC R13 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi13,Rt13);
    
    CCTK_REAL_VEC R22 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi22,Rt22);
    
    CCTK_REAL_VEC R23 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi23,Rt23);
    
    CCTK_REAL_VEC R33 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi33,Rt33);
    
    CCTK_REAL_VEC Ats11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt211,PDstandardNth2alpha,kmadd(Gt311,PDstandardNth3alpha,kmadd(alphaL,R11,kmsub(PDstandardNth1alpha,kmadd(cdphi1,ToReal(4),Gt111),PDstandardNth11alpha))));
    
    CCTK_REAL_VEC Ats12 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt312,PDstandardNth3alpha,kmadd(alphaL,R12,ksub(kmadd(PDstandardNth2alpha,kmadd(cdphi1,ToReal(2),Gt212),kmul(PDstandardNth1alpha,kmadd(cdphi2,ToReal(2),Gt112))),PDstandardNth12alpha)));
    
    CCTK_REAL_VEC Ats13 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt213,PDstandardNth2alpha,kmadd(alphaL,R13,ksub(kmadd(PDstandardNth3alpha,kmadd(cdphi1,ToReal(2),Gt313),kmul(PDstandardNth1alpha,kmadd(cdphi3,ToReal(2),Gt113))),PDstandardNth13alpha)));
    
    CCTK_REAL_VEC Ats22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt122,PDstandardNth1alpha,kmadd(Gt322,PDstandardNth3alpha,kmadd(alphaL,R22,kmsub(PDstandardNth2alpha,kmadd(cdphi2,ToReal(4),Gt222),PDstandardNth22alpha))));
    
    CCTK_REAL_VEC Ats23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt123,PDstandardNth1alpha,kmadd(alphaL,R23,ksub(kmadd(PDstandardNth3alpha,kmadd(cdphi2,ToReal(2),Gt323),kmul(PDstandardNth2alpha,kmadd(cdphi3,ToReal(2),Gt223))),PDstandardNth23alpha)));
    
    CCTK_REAL_VEC Ats33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt133,PDstandardNth1alpha,kmadd(Gt233,PDstandardNth2alpha,kmadd(alphaL,R33,kmsub(PDstandardNth3alpha,kmadd(cdphi3,ToReal(4),Gt333),PDstandardNth33alpha))));
    
    CCTK_REAL_VEC trAts CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Ats11,gu11,kmadd(Ats22,gu22,kmadd(Ats33,gu33,kmul(kmadd(Ats12,gu12,kmadd(Ats13,gu13,kmul(Ats23,gu23))),ToReal(2)))));
    
    CCTK_REAL_VEC At11rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmadd(At11L,trKL,kmul(kmadd(At11L,Atm11,kmadd(At12L,Atm21,kmul(At13L,Atm31))),ToReal(-2))),kmadd(At11L,kmul(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667)),kmadd(em4phi,kmadd(g11,kmul(trAts,ToReal(-0.333333333333333333333333333333)),Ats11),kmul(kmadd(At11L,PDstandardNth1beta1,kmadd(At12L,PDstandardNth1beta2,kmul(At13L,PDstandardNth1beta3))),ToReal(2)))));
    
    CCTK_REAL_VEC At12rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At22L,PDstandardNth1beta2,kmadd(At23L,PDstandardNth1beta3,kmadd(At11L,PDstandardNth2beta1,kmadd(At13L,PDstandardNth2beta3,kmadd(alphaL,kmadd(At12L,trKL,kmul(kmadd(At11L,Atm12,kmadd(At12L,Atm22,kmul(At13L,Atm32))),ToReal(-2))),kmadd(At12L,kadd(PDstandardNth1beta1,kmadd(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667),PDstandardNth2beta2)),kmul(em4phi,kmadd(g12,kmul(trAts,ToReal(-0.333333333333333333333333333333)),Ats12))))))));
    
    CCTK_REAL_VEC At13rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At23L,PDstandardNth1beta2,kmadd(At33L,PDstandardNth1beta3,kmadd(At11L,PDstandardNth3beta1,kmadd(At12L,PDstandardNth3beta2,kmadd(alphaL,kmadd(At13L,trKL,kmul(kmadd(At11L,Atm13,kmadd(At12L,Atm23,kmul(At13L,Atm33))),ToReal(-2))),kmadd(At13L,kadd(PDstandardNth1beta1,kmadd(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667),PDstandardNth3beta3)),kmul(em4phi,kmadd(g13,kmul(trAts,ToReal(-0.333333333333333333333333333333)),Ats13))))))));
    
    CCTK_REAL_VEC At22rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmadd(At22L,trKL,kmul(kmadd(At12L,Atm12,kmadd(At22L,Atm22,kmul(At23L,Atm32))),ToReal(-2))),kmadd(At22L,kmul(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667)),kmadd(em4phi,kmadd(g22,kmul(trAts,ToReal(-0.333333333333333333333333333333)),Ats22),kmul(kmadd(At12L,PDstandardNth2beta1,kmadd(At22L,PDstandardNth2beta2,kmul(At23L,PDstandardNth2beta3))),ToReal(2)))));
    
    CCTK_REAL_VEC At23rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At13L,PDstandardNth2beta1,kmadd(At33L,PDstandardNth2beta3,kmadd(At12L,PDstandardNth3beta1,kmadd(At22L,PDstandardNth3beta2,kmadd(alphaL,kmadd(At23L,trKL,kmul(kmadd(At12L,Atm13,kmadd(At22L,Atm23,kmul(At23L,Atm33))),ToReal(-2))),kmadd(At23L,kadd(PDstandardNth2beta2,kmadd(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667),PDstandardNth3beta3)),kmul(em4phi,kmadd(g23,kmul(trAts,ToReal(-0.333333333333333333333333333333)),Ats23))))))));
    
    CCTK_REAL_VEC At33rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmadd(At33L,trKL,kmul(kmadd(At13L,Atm13,kmadd(At23L,Atm23,kmul(At33L,Atm33))),ToReal(-2))),kmadd(At33L,kmul(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667)),kmadd(em4phi,kmadd(g33,kmul(trAts,ToReal(-0.333333333333333333333333333333)),Ats33),kmul(kmadd(At13L,PDstandardNth3beta1,kmadd(At23L,PDstandardNth3beta2,kmul(At33L,PDstandardNth3beta3))),ToReal(2)))));
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(At11rhs[index],At11rhsL);
    vec_store_nta_partial(At12rhs[index],At12rhsL);
    vec_store_nta_partial(At13rhs[index],At13rhsL);
    vec_store_nta_partial(At22rhs[index],At22rhsL);
    vec_store_nta_partial(At23rhs[index],At23rhsL);
    vec_store_nta_partial(At33rhs[index],At33rhsL);
  }
  CCTK_ENDLOOP3STR(ML_BSSN_Host_RHS2);
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
  
  const char* const groups[] = {
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
