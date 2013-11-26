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

extern "C" void ML_BSSN_Host_RHSRadiativeBoundary_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_Host_RHSRadiativeBoundary_calc_every != ML_BSSN_Host_RHSRadiativeBoundary_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_curvrhs.");
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

static void ML_BSSN_Host_RHSRadiativeBoundary_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3STR(ML_BSSN_Host_RHSRadiativeBoundary,
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
    CCTK_REAL_VEC rCopyL CCTK_ATTRIBUTE_UNUSED = vec_load(rCopy[index]);
    CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = vec_load(trK[index]);
    CCTK_REAL_VEC xCopyL CCTK_ATTRIBUTE_UNUSED = vec_load(xCopy[index]);
    CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt3[index]);
    CCTK_REAL_VEC yCopyL CCTK_ATTRIBUTE_UNUSED = vec_load(yCopy[index]);
    CCTK_REAL_VEC zCopyL CCTK_ATTRIBUTE_UNUSED = vec_load(zCopy[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    const CCTK_REAL_VEC PDstandard2nd1A CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&A[index]);
    const CCTK_REAL_VEC PDstandard2nd2A CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&A[index]);
    const CCTK_REAL_VEC PDstandard2nd3A CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&A[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1A CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&A[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1A CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&A[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2A CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&A[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2A CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&A[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3A CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&A[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3A CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&A[index]);
    const CCTK_REAL_VEC PDstandard2nd1alpha CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&alpha[index]);
    const CCTK_REAL_VEC PDstandard2nd2alpha CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&alpha[index]);
    const CCTK_REAL_VEC PDstandard2nd3alpha CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&alpha[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&alpha[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&alpha[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&alpha[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&alpha[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&alpha[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&alpha[index]);
    const CCTK_REAL_VEC PDstandard2nd1At11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At11[index]);
    const CCTK_REAL_VEC PDstandard2nd2At11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At11[index]);
    const CCTK_REAL_VEC PDstandard2nd3At11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At11[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At11[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At11[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At11[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At11[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At11[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At11[index]);
    const CCTK_REAL_VEC PDstandard2nd1At12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At12[index]);
    const CCTK_REAL_VEC PDstandard2nd2At12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At12[index]);
    const CCTK_REAL_VEC PDstandard2nd3At12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At12[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At12[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At12[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At12[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At12[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At12[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At12[index]);
    const CCTK_REAL_VEC PDstandard2nd1At13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At13[index]);
    const CCTK_REAL_VEC PDstandard2nd2At13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At13[index]);
    const CCTK_REAL_VEC PDstandard2nd3At13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At13[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At13[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At13[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At13[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At13[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At13[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At13[index]);
    const CCTK_REAL_VEC PDstandard2nd1At22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At22[index]);
    const CCTK_REAL_VEC PDstandard2nd2At22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At22[index]);
    const CCTK_REAL_VEC PDstandard2nd3At22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At22[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At22[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At22[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At22[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At22[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At22[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At22[index]);
    const CCTK_REAL_VEC PDstandard2nd1At23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At23[index]);
    const CCTK_REAL_VEC PDstandard2nd2At23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At23[index]);
    const CCTK_REAL_VEC PDstandard2nd3At23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At23[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At23[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At23[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At23[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At23[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At23[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At23[index]);
    const CCTK_REAL_VEC PDstandard2nd1At33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At33[index]);
    const CCTK_REAL_VEC PDstandard2nd2At33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At33[index]);
    const CCTK_REAL_VEC PDstandard2nd3At33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At33[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At33[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At33[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At33[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At33[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At33[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At33[index]);
    const CCTK_REAL_VEC PDstandard2nd1B1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&B1[index]);
    const CCTK_REAL_VEC PDstandard2nd2B1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&B1[index]);
    const CCTK_REAL_VEC PDstandard2nd3B1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&B1[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&B1[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&B1[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&B1[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&B1[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&B1[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&B1[index]);
    const CCTK_REAL_VEC PDstandard2nd1B2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&B2[index]);
    const CCTK_REAL_VEC PDstandard2nd2B2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&B2[index]);
    const CCTK_REAL_VEC PDstandard2nd3B2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&B2[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&B2[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&B2[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&B2[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&B2[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&B2[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&B2[index]);
    const CCTK_REAL_VEC PDstandard2nd1B3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&B3[index]);
    const CCTK_REAL_VEC PDstandard2nd2B3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&B3[index]);
    const CCTK_REAL_VEC PDstandard2nd3B3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&B3[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&B3[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&B3[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&B3[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&B3[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&B3[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&B3[index]);
    const CCTK_REAL_VEC PDstandard2nd1beta1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&beta1[index]);
    const CCTK_REAL_VEC PDstandard2nd2beta1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&beta1[index]);
    const CCTK_REAL_VEC PDstandard2nd3beta1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&beta1[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&beta1[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&beta1[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&beta1[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&beta1[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&beta1[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&beta1[index]);
    const CCTK_REAL_VEC PDstandard2nd1beta2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&beta2[index]);
    const CCTK_REAL_VEC PDstandard2nd2beta2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&beta2[index]);
    const CCTK_REAL_VEC PDstandard2nd3beta2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&beta2[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&beta2[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&beta2[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&beta2[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&beta2[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&beta2[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&beta2[index]);
    const CCTK_REAL_VEC PDstandard2nd1beta3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&beta3[index]);
    const CCTK_REAL_VEC PDstandard2nd2beta3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&beta3[index]);
    const CCTK_REAL_VEC PDstandard2nd3beta3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&beta3[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&beta3[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&beta3[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&beta3[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&beta3[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&beta3[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&beta3[index]);
    const CCTK_REAL_VEC PDstandard2nd1gt11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt11[index]);
    const CCTK_REAL_VEC PDstandard2nd2gt11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt11[index]);
    const CCTK_REAL_VEC PDstandard2nd3gt11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt11[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt11[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt11[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt11[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt11[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt11[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt11[index]);
    const CCTK_REAL_VEC PDstandard2nd1gt12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt12[index]);
    const CCTK_REAL_VEC PDstandard2nd2gt12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt12[index]);
    const CCTK_REAL_VEC PDstandard2nd3gt12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt12[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt12[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt12[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt12[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt12[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt12[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt12[index]);
    const CCTK_REAL_VEC PDstandard2nd1gt13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt13[index]);
    const CCTK_REAL_VEC PDstandard2nd2gt13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt13[index]);
    const CCTK_REAL_VEC PDstandard2nd3gt13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt13[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt13[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt13[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt13[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt13[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt13[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt13[index]);
    const CCTK_REAL_VEC PDstandard2nd1gt22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt22[index]);
    const CCTK_REAL_VEC PDstandard2nd2gt22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt22[index]);
    const CCTK_REAL_VEC PDstandard2nd3gt22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt22[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt22[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt22[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt22[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt22[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt22[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt22[index]);
    const CCTK_REAL_VEC PDstandard2nd1gt23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt23[index]);
    const CCTK_REAL_VEC PDstandard2nd2gt23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt23[index]);
    const CCTK_REAL_VEC PDstandard2nd3gt23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt23[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt23[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt23[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt23[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt23[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt23[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt23[index]);
    const CCTK_REAL_VEC PDstandard2nd1gt33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt33[index]);
    const CCTK_REAL_VEC PDstandard2nd2gt33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt33[index]);
    const CCTK_REAL_VEC PDstandard2nd3gt33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt33[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt33[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt33[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt33[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt33[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt33[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt33[index]);
    const CCTK_REAL_VEC PDstandard2nd1phi CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&phi[index]);
    const CCTK_REAL_VEC PDstandard2nd2phi CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&phi[index]);
    const CCTK_REAL_VEC PDstandard2nd3phi CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&phi[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1phi CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&phi[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1phi CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&phi[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2phi CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&phi[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2phi CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&phi[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3phi CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&phi[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3phi CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&phi[index]);
    const CCTK_REAL_VEC PDstandard2nd1trK CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&trK[index]);
    const CCTK_REAL_VEC PDstandard2nd2trK CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&trK[index]);
    const CCTK_REAL_VEC PDstandard2nd3trK CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&trK[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1trK CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&trK[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1trK CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&trK[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2trK CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&trK[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2trK CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&trK[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3trK CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&trK[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3trK CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&trK[index]);
    const CCTK_REAL_VEC PDstandard2nd1Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&Xt1[index]);
    const CCTK_REAL_VEC PDstandard2nd2Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&Xt1[index]);
    const CCTK_REAL_VEC PDstandard2nd3Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&Xt1[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&Xt1[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&Xt1[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&Xt1[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&Xt1[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&Xt1[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&Xt1[index]);
    const CCTK_REAL_VEC PDstandard2nd1Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&Xt2[index]);
    const CCTK_REAL_VEC PDstandard2nd2Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&Xt2[index]);
    const CCTK_REAL_VEC PDstandard2nd3Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&Xt2[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&Xt2[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&Xt2[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&Xt2[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&Xt2[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&Xt2[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&Xt2[index]);
    const CCTK_REAL_VEC PDstandard2nd1Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&Xt3[index]);
    const CCTK_REAL_VEC PDstandard2nd2Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&Xt3[index]);
    const CCTK_REAL_VEC PDstandard2nd3Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&Xt3[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd1Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&Xt3[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd1Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&Xt3[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd2Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&Xt3[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd2Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&Xt3[index]);
    const CCTK_REAL_VEC PDonesidedPlus2nd3Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&Xt3[index]);
    const CCTK_REAL_VEC PDonesidedMinus2nd3Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC rn1 CCTK_ATTRIBUTE_UNUSED = kneg(kdiv(xCopyL,rCopyL));
    
    CCTK_REAL_VEC rn2 CCTK_ATTRIBUTE_UNUSED = kneg(kdiv(yCopyL,rCopyL));
    
    CCTK_REAL_VEC rn3 CCTK_ATTRIBUTE_UNUSED = kneg(kdiv(zCopyL,rCopyL));
    
    CCTK_REAL_VEC phi0 CCTK_ATTRIBUTE_UNUSED = 
      IfThen(conformalMethod,ToReal(1),ToReal(0));
    
    CCTK_REAL_VEC v0 CCTK_ATTRIBUTE_UNUSED = ksqrt(ToReal(harmonicF));
    
    CCTK_REAL_VEC phirhsL CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmul(v0,kadd(phi0,kmsub(rCopyL,kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1phi,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1phi,PDstandard2nd1phi)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2phi,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2phi,PDstandard2nd2phi)),kmul(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3phi,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3phi,PDstandard2nd3phi))))),phiL))),rCopyL);
    
    CCTK_REAL_VEC gt11rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1gt11,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1gt11,PDstandard2nd1gt11)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt11,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt11,PDstandard2nd2gt11)),kmadd(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt11,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt11,PDstandard2nd3gt11)),kdiv(ksub(ToReal(1),gt11L),rCopyL))));
    
    CCTK_REAL_VEC gt12rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1gt12,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1gt12,PDstandard2nd1gt12)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt12,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt12,PDstandard2nd2gt12)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt12,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt12,PDstandard2nd3gt12)),kdiv(gt12L,rCopyL))));
    
    CCTK_REAL_VEC gt13rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1gt13,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1gt13,PDstandard2nd1gt13)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt13,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt13,PDstandard2nd2gt13)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt13,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt13,PDstandard2nd3gt13)),kdiv(gt13L,rCopyL))));
    
    CCTK_REAL_VEC gt22rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1gt22,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1gt22,PDstandard2nd1gt22)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt22,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt22,PDstandard2nd2gt22)),kmadd(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt22,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt22,PDstandard2nd3gt22)),kdiv(ksub(ToReal(1),gt22L),rCopyL))));
    
    CCTK_REAL_VEC gt23rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1gt23,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1gt23,PDstandard2nd1gt23)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt23,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt23,PDstandard2nd2gt23)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt23,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt23,PDstandard2nd3gt23)),kdiv(gt23L,rCopyL))));
    
    CCTK_REAL_VEC gt33rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1gt33,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1gt33,PDstandard2nd1gt33)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt33,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt33,PDstandard2nd2gt33)),kmadd(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt33,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt33,PDstandard2nd3gt33)),kdiv(ksub(ToReal(1),gt33L),rCopyL))));
    
    CCTK_REAL_VEC trKrhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(v0,kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1trK,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1trK,PDstandard2nd1trK)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2trK,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2trK,PDstandard2nd2trK)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3trK,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3trK,PDstandard2nd3trK)),kdiv(trKL,rCopyL)))));
    
    CCTK_REAL_VEC At11rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1At11,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1At11,PDstandard2nd1At11)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At11,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At11,PDstandard2nd2At11)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At11,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At11,PDstandard2nd3At11)),kdiv(At11L,rCopyL))));
    
    CCTK_REAL_VEC At12rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1At12,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1At12,PDstandard2nd1At12)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At12,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At12,PDstandard2nd2At12)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At12,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At12,PDstandard2nd3At12)),kdiv(At12L,rCopyL))));
    
    CCTK_REAL_VEC At13rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1At13,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1At13,PDstandard2nd1At13)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At13,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At13,PDstandard2nd2At13)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At13,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At13,PDstandard2nd3At13)),kdiv(At13L,rCopyL))));
    
    CCTK_REAL_VEC At22rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1At22,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1At22,PDstandard2nd1At22)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At22,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At22,PDstandard2nd2At22)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At22,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At22,PDstandard2nd3At22)),kdiv(At22L,rCopyL))));
    
    CCTK_REAL_VEC At23rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1At23,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1At23,PDstandard2nd1At23)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At23,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At23,PDstandard2nd2At23)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At23,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At23,PDstandard2nd3At23)),kdiv(At23L,rCopyL))));
    
    CCTK_REAL_VEC At33rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1At33,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1At33,PDstandard2nd1At33)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At33,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At33,PDstandard2nd2At33)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At33,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At33,PDstandard2nd3At33)),kdiv(At33L,rCopyL))));
    
    CCTK_REAL_VEC Xt1rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1Xt1,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt1,PDstandard2nd1Xt1)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2Xt1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt1,PDstandard2nd2Xt1)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3Xt1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt1,PDstandard2nd3Xt1)),kdiv(Xt1L,rCopyL))));
    
    CCTK_REAL_VEC Xt2rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1Xt2,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt2,PDstandard2nd1Xt2)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2Xt2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt2,PDstandard2nd2Xt2)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3Xt2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt2,PDstandard2nd3Xt2)),kdiv(Xt2L,rCopyL))));
    
    CCTK_REAL_VEC Xt3rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1Xt3,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt3,PDstandard2nd1Xt3)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2Xt3,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt3,PDstandard2nd2Xt3)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3Xt3,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt3,PDstandard2nd3Xt3)),kdiv(Xt3L,rCopyL))));
    
    CCTK_REAL_VEC alpharhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(v0,kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1alpha,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1alpha,PDstandard2nd1alpha)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2alpha,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2alpha,PDstandard2nd2alpha)),kmadd(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3alpha,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3alpha,PDstandard2nd3alpha)),kdiv(ksub(ToReal(1),alphaL),rCopyL)))));
    
    CCTK_REAL_VEC ArhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(v0,kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1A,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1A,PDstandard2nd1A)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2A,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2A,PDstandard2nd2A)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3A,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3A,PDstandard2nd3A)),kdiv(AL,rCopyL)))));
    
    CCTK_REAL_VEC beta1rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1beta1,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1beta1,PDstandard2nd1beta1)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2beta1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta1,PDstandard2nd2beta1)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3beta1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta1,PDstandard2nd3beta1)),kdiv(beta1L,rCopyL))));
    
    CCTK_REAL_VEC beta2rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1beta2,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1beta2,PDstandard2nd1beta2)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2beta2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta2,PDstandard2nd2beta2)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3beta2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta2,PDstandard2nd3beta2)),kdiv(beta2L,rCopyL))));
    
    CCTK_REAL_VEC beta3rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1beta3,IfThen(normal[0] 
      > 
      0,PDonesidedMinus2nd1beta3,PDstandard2nd1beta3)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2beta3,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta3,PDstandard2nd2beta3)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3beta3,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta3,PDstandard2nd3beta3)),kdiv(beta3L,rCopyL))));
    
    CCTK_REAL_VEC B1rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1B1,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1B1,PDstandard2nd1B1)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2B1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2B1,PDstandard2nd2B1)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3B1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3B1,PDstandard2nd3B1)),kdiv(B1L,rCopyL))));
    
    CCTK_REAL_VEC B2rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1B2,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1B2,PDstandard2nd1B2)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2B2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2B2,PDstandard2nd2B2)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3B2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3B2,PDstandard2nd3B2)),kdiv(B2L,rCopyL))));
    
    CCTK_REAL_VEC B3rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(rn1,IfThen(normal[0] < 0,PDonesidedPlus2nd1B3,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1B3,PDstandard2nd1B3)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2B3,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2B3,PDstandard2nd2B3)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3B3,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3B3,PDstandard2nd3B3)),kdiv(B3L,rCopyL))));
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(alpharhs[index],alpharhsL);
    vec_store_nta_partial(Arhs[index],ArhsL);
    vec_store_nta_partial(At11rhs[index],At11rhsL);
    vec_store_nta_partial(At12rhs[index],At12rhsL);
    vec_store_nta_partial(At13rhs[index],At13rhsL);
    vec_store_nta_partial(At22rhs[index],At22rhsL);
    vec_store_nta_partial(At23rhs[index],At23rhsL);
    vec_store_nta_partial(At33rhs[index],At33rhsL);
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
  CCTK_ENDLOOP3STR(ML_BSSN_Host_RHSRadiativeBoundary);
}

extern "C" void ML_BSSN_Host_RHSRadiativeBoundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_Host_RHSRadiativeBoundary_Body");
  }
  
  if (cctk_iteration % ML_BSSN_Host_RHSRadiativeBoundary_calc_every != ML_BSSN_Host_RHSRadiativeBoundary_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN_Host::coords",
    "ML_BSSN_Host::ML_curv",
    "ML_BSSN_Host::ML_curvrhs",
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
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_RHSRadiativeBoundary", 19, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Host_RHSRadiativeBoundary", 2, 2, 2);
  
  GenericFD_LoopOverBoundary(cctkGH, ML_BSSN_Host_RHSRadiativeBoundary_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_RHSRadiativeBoundary_Body");
  }
}
