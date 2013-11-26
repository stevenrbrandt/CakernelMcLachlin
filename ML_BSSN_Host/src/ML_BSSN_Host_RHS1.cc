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

extern "C" void ML_BSSN_Host_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_Host_RHS1_calc_every != ML_BSSN_Host_RHS1_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
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

static void ML_BSSN_Host_RHS1_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3STR(ML_BSSN_Host_RHS1,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    // vec_iter_counter+=CCTK_REAL_VEC_SIZE;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC AL CCTK_ATTRIBUTE_UNUSED = vec_load(A[index]);
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = vec_load(alpha[index]);
    CCTK_REAL_VEC At11L CCTK_ATTRIBUTE_UNUSED = vec_load(At11[index]);
    CCTK_REAL_VEC At12L CCTK_ATTRIBUTE_UNUSED = vec_load(At12[index]);
    CCTK_REAL_VEC At13L CCTK_ATTRIBUTE_UNUSED = vec_load(At13[index]);
    CCTK_REAL_VEC At22L CCTK_ATTRIBUTE_UNUSED = vec_load(At22[index]);
    CCTK_REAL_VEC At23L CCTK_ATTRIBUTE_UNUSED = vec_load(At23[index]);
    CCTK_REAL_VEC At33L CCTK_ATTRIBUTE_UNUSED = vec_load(At33[index]);
    CCTK_REAL_VEC B1L CCTK_ATTRIBUTE_UNUSED = vec_load(B1[index]);
    CCTK_REAL_VEC B2L CCTK_ATTRIBUTE_UNUSED = vec_load(B2[index]);
    CCTK_REAL_VEC B3L CCTK_ATTRIBUTE_UNUSED = vec_load(B3[index]);
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
    const CCTK_REAL_VEC PDstandardNth11beta1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&beta1[index]);
    const CCTK_REAL_VEC PDstandardNth22beta1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&beta1[index]);
    const CCTK_REAL_VEC PDstandardNth33beta1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&beta1[index]);
    const CCTK_REAL_VEC PDstandardNth12beta1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&beta1[index]);
    const CCTK_REAL_VEC PDstandardNth13beta1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&beta1[index]);
    const CCTK_REAL_VEC PDstandardNth23beta1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&beta1[index]);
    const CCTK_REAL_VEC PDstandardNth1beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth2beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth3beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth11beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth22beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth33beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth12beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth13beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth23beta2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&beta2[index]);
    const CCTK_REAL_VEC PDstandardNth1beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth2beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth3beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth11beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth22beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth33beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth12beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth13beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&beta3[index]);
    const CCTK_REAL_VEC PDstandardNth23beta3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&beta3[index]);
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
    const CCTK_REAL_VEC PDstandardNth1phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth2phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth3phi CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&phi[index]);
    const CCTK_REAL_VEC PDstandardNth1trK CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&trK[index]);
    const CCTK_REAL_VEC PDstandardNth2trK CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&trK[index]);
    const CCTK_REAL_VEC PDstandardNth3trK CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&trK[index]);
    
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
    
    CCTK_REAL_VEC fac1 CCTK_ATTRIBUTE_UNUSED = 
      IfThen(conformalMethod,kdiv(ToReal(-0.5),phiL),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,PDstandardNth1phi);
    
    CCTK_REAL_VEC cdphi2 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,PDstandardNth2phi);
    
    CCTK_REAL_VEC cdphi3 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,PDstandardNth3phi);
    
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
    
    CCTK_REAL_VEC Atu11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm11,gtu11,kmadd(Atm12,gtu12,kmul(Atm13,gtu13)));
    
    CCTK_REAL_VEC Atu12 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm11,gtu12,kmadd(Atm12,gtu22,kmul(Atm13,gtu23)));
    
    CCTK_REAL_VEC Atu13 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm11,gtu13,kmadd(Atm12,gtu23,kmul(Atm13,gtu33)));
    
    CCTK_REAL_VEC Atu22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm21,gtu12,kmadd(Atm22,gtu22,kmul(Atm23,gtu23)));
    
    CCTK_REAL_VEC Atu23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm21,gtu13,kmadd(Atm22,gtu23,kmul(Atm23,gtu33)));
    
    CCTK_REAL_VEC Atu33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm31,gtu13,kmadd(Atm32,gtu23,kmul(Atm33,gtu33)));
    
    CCTK_REAL_VEC e4phi CCTK_ATTRIBUTE_UNUSED = 
      IfThen(conformalMethod,kdiv(ToReal(1),kmul(phiL,phiL)),kexp(kmul(phiL,ToReal(4))));
    
    CCTK_REAL_VEC em4phi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),e4phi);
    
    CCTK_REAL_VEC phirhsL CCTK_ATTRIBUTE_UNUSED = 
      IfThen(conformalMethod,kmul(phiL,kmadd(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(-0.333333333333333333333333333333),kmul(alphaL,kmul(trKL,ToReal(0.333333333333333333333333333333))))),kmadd(alphaL,kmul(trKL,ToReal(-0.166666666666666666666666666667)),kmul(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(0.166666666666666666666666666667))));
    
    CCTK_REAL_VEC gt11rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmul(At11L,ToReal(-2)),kmadd(gt11L,kmul(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667)),kmul(kmadd(gt11L,PDstandardNth1beta1,kmadd(gt12L,PDstandardNth1beta2,kmul(gt13L,PDstandardNth1beta3))),ToReal(2))));
    
    CCTK_REAL_VEC gt12rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At12L,ToReal(-6)),kmadd(gt12L,kadd(PDstandardNth1beta1,kmadd(PDstandardNth3beta3,ToReal(-2),PDstandardNth2beta2)),kmul(kmadd(gt22L,PDstandardNth1beta2,kmadd(gt23L,PDstandardNth1beta3,kmadd(gt11L,PDstandardNth2beta1,kmul(gt13L,PDstandardNth2beta3)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt13rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At13L,ToReal(-6)),kmadd(gt13L,kadd(PDstandardNth1beta1,kmadd(PDstandardNth2beta2,ToReal(-2),PDstandardNth3beta3)),kmul(kmadd(gt23L,PDstandardNth1beta2,kmadd(gt33L,PDstandardNth1beta3,kmadd(gt11L,PDstandardNth3beta1,kmul(gt12L,PDstandardNth3beta2)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt22rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmul(At22L,ToReal(-2)),kmadd(gt22L,kmul(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667)),kmul(kmadd(gt12L,PDstandardNth2beta1,kmadd(gt22L,PDstandardNth2beta2,kmul(gt23L,PDstandardNth2beta3))),ToReal(2))));
    
    CCTK_REAL_VEC gt23rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.333333333333333333333333333333),kmadd(alphaL,kmul(At23L,ToReal(-6)),kmadd(gt23L,kadd(PDstandardNth2beta2,kmadd(PDstandardNth1beta1,ToReal(-2),PDstandardNth3beta3)),kmul(kmadd(gt13L,PDstandardNth2beta1,kmadd(gt33L,PDstandardNth2beta3,kmadd(gt12L,PDstandardNth3beta1,kmul(gt22L,PDstandardNth3beta2)))),ToReal(3)))));
    
    CCTK_REAL_VEC gt33rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmul(At33L,ToReal(-2)),kmadd(gt33L,kmul(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667)),kmul(kmadd(gt13L,PDstandardNth3beta1,kmadd(gt23L,PDstandardNth3beta2,kmul(gt33L,PDstandardNth3beta3))),ToReal(2))));
    
    CCTK_REAL_VEC dotXt1 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,PDstandardNth11beta1,kmadd(gtu22,PDstandardNth22beta1,kmadd(gtu33,PDstandardNth33beta1,kmadd(kmadd(Atu11,PDstandardNth1alpha,kmadd(Atu12,PDstandardNth2alpha,kmul(Atu13,PDstandardNth3alpha))),ToReal(-2),knmsub(PDstandardNth2beta1,Xtn2,knmsub(PDstandardNth3beta1,Xtn3,kmadd(kmadd(gtu11,kadd(PDstandardNth11beta1,kadd(PDstandardNth12beta2,PDstandardNth13beta3)),kmadd(gtu12,kadd(PDstandardNth12beta1,kadd(PDstandardNth22beta2,PDstandardNth23beta3)),kmul(gtu13,kadd(PDstandardNth13beta1,kadd(PDstandardNth23beta2,PDstandardNth33beta3))))),ToReal(0.333333333333333333333333333333),kmadd(Xtn1,kmsub(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(0.666666666666666666666666666667),PDstandardNth1beta1),kmul(ToReal(2),kmadd(gtu12,PDstandardNth12beta1,kmadd(gtu13,PDstandardNth13beta1,kmadd(gtu23,PDstandardNth23beta1,kmul(alphaL,kmadd(Atu11,Gt111,kmadd(Atu22,Gt122,kmadd(Atu33,Gt133,kmadd(kmadd(gtu11,PDstandardNth1trK,kmadd(gtu12,PDstandardNth2trK,kmul(gtu13,PDstandardNth3trK))),ToReal(-0.666666666666666666666666666667),kmadd(Atu12,kmul(Gt112,ToReal(2)),kmadd(Atu13,kmul(Gt113,ToReal(2)),kmadd(Atu23,kmul(Gt123,ToReal(2)),kmul(kmadd(Atu11,cdphi1,kmadd(Atu12,cdphi2,kmul(Atu13,cdphi3))),ToReal(6))))))))))))))))))))));
    
    CCTK_REAL_VEC dotXt2 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,PDstandardNth11beta2,kmadd(gtu22,PDstandardNth22beta2,kmadd(gtu33,PDstandardNth33beta2,kmadd(kmadd(Atu12,PDstandardNth1alpha,kmadd(Atu22,PDstandardNth2alpha,kmul(Atu23,PDstandardNth3alpha))),ToReal(-2),knmsub(PDstandardNth1beta2,Xtn1,knmsub(PDstandardNth3beta2,Xtn3,kmadd(kmadd(gtu12,kadd(PDstandardNth11beta1,kadd(PDstandardNth12beta2,PDstandardNth13beta3)),kmadd(gtu22,kadd(PDstandardNth12beta1,kadd(PDstandardNth22beta2,PDstandardNth23beta3)),kmul(gtu23,kadd(PDstandardNth13beta1,kadd(PDstandardNth23beta2,PDstandardNth33beta3))))),ToReal(0.333333333333333333333333333333),kmadd(Xtn2,kmsub(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(0.666666666666666666666666666667),PDstandardNth2beta2),kmul(ToReal(2),kmadd(gtu12,PDstandardNth12beta2,kmadd(gtu13,PDstandardNth13beta2,kmadd(gtu23,PDstandardNth23beta2,kmul(alphaL,kmadd(Atu11,Gt211,kmadd(Atu22,Gt222,kmadd(Atu33,Gt233,kmadd(kmadd(gtu12,PDstandardNth1trK,kmadd(gtu22,PDstandardNth2trK,kmul(gtu23,PDstandardNth3trK))),ToReal(-0.666666666666666666666666666667),kmadd(Atu12,kmul(Gt212,ToReal(2)),kmadd(Atu13,kmul(Gt213,ToReal(2)),kmadd(Atu23,kmul(Gt223,ToReal(2)),kmul(kmadd(Atu12,cdphi1,kmadd(Atu22,cdphi2,kmul(Atu23,cdphi3))),ToReal(6))))))))))))))))))))));
    
    CCTK_REAL_VEC dotXt3 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,PDstandardNth11beta3,kmadd(gtu22,PDstandardNth22beta3,kmadd(gtu33,PDstandardNth33beta3,kmadd(kmadd(Atu13,PDstandardNth1alpha,kmadd(Atu23,PDstandardNth2alpha,kmul(Atu33,PDstandardNth3alpha))),ToReal(-2),knmsub(PDstandardNth1beta3,Xtn1,knmsub(PDstandardNth2beta3,Xtn2,kmadd(kmadd(gtu13,kadd(PDstandardNth11beta1,kadd(PDstandardNth12beta2,PDstandardNth13beta3)),kmadd(gtu23,kadd(PDstandardNth12beta1,kadd(PDstandardNth22beta2,PDstandardNth23beta3)),kmul(gtu33,kadd(PDstandardNth13beta1,kadd(PDstandardNth23beta2,PDstandardNth33beta3))))),ToReal(0.333333333333333333333333333333),kmadd(Xtn3,kmsub(kadd(PDstandardNth1beta1,kadd(PDstandardNth2beta2,PDstandardNth3beta3)),ToReal(0.666666666666666666666666666667),PDstandardNth3beta3),kmul(ToReal(2),kmadd(gtu12,PDstandardNth12beta3,kmadd(gtu13,PDstandardNth13beta3,kmadd(gtu23,PDstandardNth23beta3,kmul(alphaL,kmadd(Atu11,Gt311,kmadd(Atu22,Gt322,kmadd(Atu33,Gt333,kmadd(kmadd(gtu13,PDstandardNth1trK,kmadd(gtu23,PDstandardNth2trK,kmul(gtu33,PDstandardNth3trK))),ToReal(-0.666666666666666666666666666667),kmadd(Atu12,kmul(Gt312,ToReal(2)),kmadd(Atu13,kmul(Gt313,ToReal(2)),kmadd(Atu23,kmul(Gt323,ToReal(2)),kmul(kmadd(Atu13,cdphi1,kmadd(Atu23,cdphi2,kmul(Atu33,cdphi3))),ToReal(6))))))))))))))))))))));
    
    CCTK_REAL_VEC Xt1rhsL CCTK_ATTRIBUTE_UNUSED = dotXt1;
    
    CCTK_REAL_VEC Xt2rhsL CCTK_ATTRIBUTE_UNUSED = dotXt2;
    
    CCTK_REAL_VEC Xt3rhsL CCTK_ATTRIBUTE_UNUSED = dotXt3;
    
    CCTK_REAL_VEC dottrK CCTK_ATTRIBUTE_UNUSED = 
      kmsub(alphaL,kmadd(Atm11,Atm11,kmadd(Atm22,Atm22,kmadd(Atm33,Atm33,kmadd(kmul(trKL,trKL),ToReal(0.333333333333333333333333333333),kmul(kmadd(Atm12,Atm21,kmadd(Atm13,Atm31,kmul(Atm23,Atm32))),ToReal(2)))))),kmul(em4phi,kmadd(gtu11,PDstandardNth11alpha,kmadd(gtu22,PDstandardNth22alpha,knmsub(PDstandardNth3alpha,Xtn3,kmadd(kmadd(gtu12,PDstandardNth12alpha,kmadd(gtu13,kmadd(cdphi1,PDstandardNth3alpha,PDstandardNth13alpha),kmul(gtu23,kmadd(cdphi2,PDstandardNth3alpha,PDstandardNth23alpha)))),ToReal(2),kmadd(PDstandardNth1alpha,kmsub(kmadd(cdphi1,gtu11,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13))),ToReal(2),Xtn1),kmadd(PDstandardNth2alpha,kmsub(kmadd(cdphi1,gtu12,kmadd(cdphi2,gtu22,kmul(cdphi3,gtu23))),ToReal(2),Xtn2),kmul(gtu33,kmadd(cdphi3,kmul(PDstandardNth3alpha,ToReal(2)),PDstandardNth33alpha))))))))));
    
    CCTK_REAL_VEC trKrhsL CCTK_ATTRIBUTE_UNUSED = dottrK;
    
    CCTK_REAL_VEC alpharhsL CCTK_ATTRIBUTE_UNUSED = 
      kneg(kmul(kpow(alphaL,harmonicN),kmul(ToReal(harmonicF),kmadd(AL,ToReal(LapseACoeff),kmul(kmadd(kadd(alphaL,ToReal(-1)),ToReal(AlphaDriver),trKL),ToReal(1 
      - LapseACoeff))))));
    
    CCTK_REAL_VEC ArhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(knmsub(AL,ToReal(AlphaDriver),dottrK),ToReal(LapseACoeff));
    
    CCTK_REAL_VEC eta CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC theta CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC beta1rhsL CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC beta2rhsL CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC beta3rhsL CCTK_ATTRIBUTE_UNUSED;
    
    if (harmonicShift)
    {
      beta1rhsL = 
        kmul(alphaL,kmul(phiL,kmul(ToReal(0.5),kmadd(alphaL,kmadd(kmadd(gtu12,kmadd(phiL,kmul(gtu23,PDstandardNth3gt22),PDstandardNth2phi),kmul(gtu13,PDstandardNth3phi)),ToReal(2),kmadd(gtu11,kmadd(phiL,kmadd(gtu12,PDstandardNth2gt11,kmul(gtu13,PDstandardNth3gt11)),kmul(kmadd(phiL,kmul(gtu23,PDstandardNth2gt13),PDstandardNth1phi),ToReal(2))),kmul(phiL,kmadd(gtu33,kmsub(gtu13,PDstandardNth3gt33,kmul(gtu12,PDstandardNth2gt33)),kmul(kmadd(gtu13,kmul(gtu22,PDstandardNth2gt23),kmadd(gtu11,kmul(gtu33,PDstandardNth3gt13),kmul(PDstandardNth1gt22,kmul(gtu12,gtu12)))),ToReal(2)))))),kmul(phiL,kmadd(kmadd(gtu12,PDstandardNth2alpha,kmul(gtu13,PDstandardNth3alpha)),ToReal(-2),kmadd(gtu11,kmadd(kmadd(alphaL,kmul(gtu23,PDstandardNth1gt23),PDstandardNth1alpha),ToReal(-2),kmul(alphaL,kmsub(kmadd(gtu12,PDstandardNth1gt12,kmadd(gtu13,PDstandardNth1gt13,kmadd(gtu22,PDstandardNth2gt12,kmul(gtu23,PDstandardNth3gt12)))),ToReal(2),kmul(gtu22,PDstandardNth1gt22)))),kmul(alphaL,kmadd(PDstandardNth1gt11,kmul(gtu11,gtu11),kmadd(PDstandardNth1gt33,kmul(kmul(gtu13,gtu13),ToReal(2)),kmadd(gtu13,kmsub(gtu23,kmul(PDstandardNth2gt33,ToReal(2)),kmul(gtu22,PDstandardNth3gt22)),kmadd(gtu33,kmsub(gtu12,kmul(PDstandardNth3gt23,ToReal(2)),kmul(gtu11,PDstandardNth1gt33)),kmul(gtu12,kmadd(gtu22,PDstandardNth2gt22,kmul(gtu13,kmul(PDstandardNth1gt23,ToReal(4)))))))))))))))));
      
      beta2rhsL = 
        kmul(alphaL,kmul(phiL,kmul(ToReal(0.5),kmadd(alphaL,kmadd(kmadd(gtu12,PDstandardNth1phi,kmadd(gtu22,kmadd(phiL,kmadd(gtu23,PDstandardNth2gt23,kmul(gtu13,PDstandardNth3gt12)),PDstandardNth2phi),kmul(gtu23,PDstandardNth3phi))),ToReal(2),kmul(phiL,kmadd(gtu23,kmul(gtu33,PDstandardNth3gt33),kmadd(PDstandardNth2gt22,kmul(gtu22,gtu22),kmadd(gtu22,kmadd(gtu23,PDstandardNth3gt22,kmsub(gtu13,kmul(PDstandardNth2gt13,ToReal(-2)),kmul(gtu11,PDstandardNth2gt11))),kmul(gtu12,kmul(gtu13,kmul(PDstandardNth3gt11,ToReal(2))))))))),kmul(phiL,kmadd(gtu23,kmsub(PDstandardNth3alpha,ToReal(-2),kmul(alphaL,kmul(gtu11,PDstandardNth3gt11))),kmadd(alphaL,kmul(kmadd(gtu11,kmadd(gtu22,PDstandardNth1gt12,kmul(gtu23,PDstandardNth1gt13)),kmadd(gtu13,kmadd(gtu22,PDstandardNth1gt23,kmul(gtu23,PDstandardNth1gt33)),kmadd(gtu12,kmadd(gtu22,PDstandardNth2gt12,kmul(gtu33,PDstandardNth3gt13)),kmadd(PDstandardNth2gt11,kmul(gtu12,gtu12),kmul(PDstandardNth2gt33,kmul(gtu23,gtu23)))))),ToReal(2)),kmadd(gtu22,kmadd(PDstandardNth2alpha,ToReal(-2),kmul(alphaL,kmul(gtu33,kmsub(PDstandardNth3gt23,ToReal(2),PDstandardNth2gt33)))),kmul(gtu12,kmadd(PDstandardNth1alpha,ToReal(-2),kmul(alphaL,kmadd(gtu11,PDstandardNth1gt11,kmadd(gtu22,PDstandardNth1gt22,kmsub(gtu23,kmul(PDstandardNth2gt13,ToReal(4)),kmul(gtu33,PDstandardNth1gt33)))))))))))))));
      
      beta3rhsL = 
        kmul(alphaL,kmul(phiL,kmul(ToReal(0.5),kmadd(alphaL,kmul(kmadd(gtu23,PDstandardNth2phi,kmadd(gtu13,kmadd(phiL,kmadd(gtu12,PDstandardNth2gt11,kmadd(gtu22,PDstandardNth2gt12,kmul(gtu33,PDstandardNth3gt13))),PDstandardNth1phi),kmadd(gtu33,PDstandardNth3phi,kmul(phiL,kmadd(gtu33,kmadd(gtu22,PDstandardNth2gt23,kmul(gtu23,PDstandardNth3gt23)),kmul(PDstandardNth3gt22,kmul(gtu23,gtu23))))))),ToReal(2)),kmul(phiL,kmadd(gtu23,kmadd(PDstandardNth2alpha,ToReal(-2),kmul(alphaL,kmadd(gtu22,PDstandardNth2gt22,kmsub(gtu33,PDstandardNth2gt33,kmul(gtu11,PDstandardNth2gt11))))),kmadd(gtu33,kmsub(kmadd(alphaL,kmul(gtu12,PDstandardNth3gt12),PDstandardNth3alpha),ToReal(-2),kmul(alphaL,kmadd(gtu22,PDstandardNth3gt22,kmul(gtu11,PDstandardNth3gt11)))),kmadd(alphaL,kmadd(PDstandardNth3gt33,kmul(gtu33,gtu33),kmul(kmadd(gtu11,kmadd(gtu23,PDstandardNth1gt12,kmul(gtu33,PDstandardNth1gt13)),kmadd(gtu12,kmadd(gtu23,PDstandardNth1gt22,kmul(gtu33,kadd(PDstandardNth1gt23,PDstandardNth2gt13))),kmul(PDstandardNth3gt11,kmul(gtu13,gtu13)))),ToReal(2))),kmul(gtu13,kmadd(PDstandardNth1alpha,ToReal(-2),kmul(alphaL,kmadd(gtu11,PDstandardNth1gt11,kmadd(gtu33,PDstandardNth1gt33,kmsub(gtu23,kmul(PDstandardNth3gt12,ToReal(4)),kmul(gtu22,PDstandardNth1gt22)))))))))))))));
    }
    else
    {
      beta1rhsL = 
        kmul(theta,kmul(kadd(Xt1L,kmadd(ksub(B1L,Xt1L),ToReal(ShiftBCoeff),kmul(beta1L,kmul(eta,ToReal(BetaDriver*(-1 
        + ShiftBCoeff)))))),ToReal(ShiftGammaCoeff)));
      
      beta2rhsL = 
        kmul(theta,kmul(kadd(Xt2L,kmadd(ksub(B2L,Xt2L),ToReal(ShiftBCoeff),kmul(beta2L,kmul(eta,ToReal(BetaDriver*(-1 
        + ShiftBCoeff)))))),ToReal(ShiftGammaCoeff)));
      
      beta3rhsL = 
        kmul(theta,kmul(kadd(Xt3L,kmadd(ksub(B3L,Xt3L),ToReal(ShiftBCoeff),kmul(beta3L,kmul(eta,ToReal(BetaDriver*(-1 
        + ShiftBCoeff)))))),ToReal(ShiftGammaCoeff)));
    }
    
    CCTK_REAL_VEC B1rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(knmsub(B1L,kmul(eta,ToReal(BetaDriver)),dotXt1),ToReal(ShiftBCoeff));
    
    CCTK_REAL_VEC B2rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(knmsub(B2L,kmul(eta,ToReal(BetaDriver)),dotXt2),ToReal(ShiftBCoeff));
    
    CCTK_REAL_VEC B3rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(knmsub(B3L,kmul(eta,ToReal(BetaDriver)),dotXt3),ToReal(ShiftBCoeff));
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
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
  CCTK_ENDLOOP3STR(ML_BSSN_Host_RHS1);
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
  
  const char* const groups[] = {
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
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_RHS1", 17, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Host_RHS1", 4, 4, 4);
  
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_Host_RHS1_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_RHS1_Body");
  }
}
