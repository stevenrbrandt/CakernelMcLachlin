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

extern "C" void HOST_ML_BSSN_RHSRadiativeBoundary_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % HOST_ML_BSSN_RHSRadiativeBoundary_calc_every != HOST_ML_BSSN_RHSRadiativeBoundary_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_curvrhs.");
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

static void HOST_ML_BSSN_RHSRadiativeBoundary_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(HOST_ML_BSSN_RHSRadiativeBoundary,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    // ++vec_iter_counter;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL AL CCTK_ATTRIBUTE_UNUSED = A[index];
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpha[index];
    CCTK_REAL At11L CCTK_ATTRIBUTE_UNUSED = At11[index];
    CCTK_REAL At12L CCTK_ATTRIBUTE_UNUSED = At12[index];
    CCTK_REAL At13L CCTK_ATTRIBUTE_UNUSED = At13[index];
    CCTK_REAL At22L CCTK_ATTRIBUTE_UNUSED = At22[index];
    CCTK_REAL At23L CCTK_ATTRIBUTE_UNUSED = At23[index];
    CCTK_REAL At33L CCTK_ATTRIBUTE_UNUSED = At33[index];
    CCTK_REAL B1L CCTK_ATTRIBUTE_UNUSED = B1[index];
    CCTK_REAL B2L CCTK_ATTRIBUTE_UNUSED = B2[index];
    CCTK_REAL B3L CCTK_ATTRIBUTE_UNUSED = B3[index];
    CCTK_REAL beta1L CCTK_ATTRIBUTE_UNUSED = beta1[index];
    CCTK_REAL beta2L CCTK_ATTRIBUTE_UNUSED = beta2[index];
    CCTK_REAL beta3L CCTK_ATTRIBUTE_UNUSED = beta3[index];
    CCTK_REAL gt11L CCTK_ATTRIBUTE_UNUSED = gt11[index];
    CCTK_REAL gt12L CCTK_ATTRIBUTE_UNUSED = gt12[index];
    CCTK_REAL gt13L CCTK_ATTRIBUTE_UNUSED = gt13[index];
    CCTK_REAL gt22L CCTK_ATTRIBUTE_UNUSED = gt22[index];
    CCTK_REAL gt23L CCTK_ATTRIBUTE_UNUSED = gt23[index];
    CCTK_REAL gt33L CCTK_ATTRIBUTE_UNUSED = gt33[index];
    CCTK_REAL phiL CCTK_ATTRIBUTE_UNUSED = phi[index];
    CCTK_REAL rCopyL CCTK_ATTRIBUTE_UNUSED = rCopy[index];
    CCTK_REAL trKL CCTK_ATTRIBUTE_UNUSED = trK[index];
    CCTK_REAL xCopyL CCTK_ATTRIBUTE_UNUSED = xCopy[index];
    CCTK_REAL Xt1L CCTK_ATTRIBUTE_UNUSED = Xt1[index];
    CCTK_REAL Xt2L CCTK_ATTRIBUTE_UNUSED = Xt2[index];
    CCTK_REAL Xt3L CCTK_ATTRIBUTE_UNUSED = Xt3[index];
    CCTK_REAL yCopyL CCTK_ATTRIBUTE_UNUSED = yCopy[index];
    CCTK_REAL zCopyL CCTK_ATTRIBUTE_UNUSED = zCopy[index];
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    const CCTK_REAL PDstandard2nd1A CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&A[index]);
    const CCTK_REAL PDstandard2nd2A CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&A[index]);
    const CCTK_REAL PDstandard2nd3A CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&A[index]);
    const CCTK_REAL PDonesidedPlus2nd1A CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&A[index]);
    const CCTK_REAL PDonesidedMinus2nd1A CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&A[index]);
    const CCTK_REAL PDonesidedPlus2nd2A CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&A[index]);
    const CCTK_REAL PDonesidedMinus2nd2A CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&A[index]);
    const CCTK_REAL PDonesidedPlus2nd3A CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&A[index]);
    const CCTK_REAL PDonesidedMinus2nd3A CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&A[index]);
    const CCTK_REAL PDstandard2nd1alpha CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&alpha[index]);
    const CCTK_REAL PDstandard2nd2alpha CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&alpha[index]);
    const CCTK_REAL PDstandard2nd3alpha CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&alpha[index]);
    const CCTK_REAL PDonesidedPlus2nd1alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&alpha[index]);
    const CCTK_REAL PDonesidedMinus2nd1alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&alpha[index]);
    const CCTK_REAL PDonesidedPlus2nd2alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&alpha[index]);
    const CCTK_REAL PDonesidedMinus2nd2alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&alpha[index]);
    const CCTK_REAL PDonesidedPlus2nd3alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&alpha[index]);
    const CCTK_REAL PDonesidedMinus2nd3alpha CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&alpha[index]);
    const CCTK_REAL PDstandard2nd1At11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At11[index]);
    const CCTK_REAL PDstandard2nd2At11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At11[index]);
    const CCTK_REAL PDstandard2nd3At11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At11[index]);
    const CCTK_REAL PDonesidedPlus2nd1At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At11[index]);
    const CCTK_REAL PDonesidedMinus2nd1At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At11[index]);
    const CCTK_REAL PDonesidedPlus2nd2At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At11[index]);
    const CCTK_REAL PDonesidedMinus2nd2At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At11[index]);
    const CCTK_REAL PDonesidedPlus2nd3At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At11[index]);
    const CCTK_REAL PDonesidedMinus2nd3At11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At11[index]);
    const CCTK_REAL PDstandard2nd1At12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At12[index]);
    const CCTK_REAL PDstandard2nd2At12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At12[index]);
    const CCTK_REAL PDstandard2nd3At12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At12[index]);
    const CCTK_REAL PDonesidedPlus2nd1At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At12[index]);
    const CCTK_REAL PDonesidedMinus2nd1At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At12[index]);
    const CCTK_REAL PDonesidedPlus2nd2At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At12[index]);
    const CCTK_REAL PDonesidedMinus2nd2At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At12[index]);
    const CCTK_REAL PDonesidedPlus2nd3At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At12[index]);
    const CCTK_REAL PDonesidedMinus2nd3At12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At12[index]);
    const CCTK_REAL PDstandard2nd1At13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At13[index]);
    const CCTK_REAL PDstandard2nd2At13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At13[index]);
    const CCTK_REAL PDstandard2nd3At13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At13[index]);
    const CCTK_REAL PDonesidedPlus2nd1At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At13[index]);
    const CCTK_REAL PDonesidedMinus2nd1At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At13[index]);
    const CCTK_REAL PDonesidedPlus2nd2At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At13[index]);
    const CCTK_REAL PDonesidedMinus2nd2At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At13[index]);
    const CCTK_REAL PDonesidedPlus2nd3At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At13[index]);
    const CCTK_REAL PDonesidedMinus2nd3At13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At13[index]);
    const CCTK_REAL PDstandard2nd1At22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At22[index]);
    const CCTK_REAL PDstandard2nd2At22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At22[index]);
    const CCTK_REAL PDstandard2nd3At22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At22[index]);
    const CCTK_REAL PDonesidedPlus2nd1At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At22[index]);
    const CCTK_REAL PDonesidedMinus2nd1At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At22[index]);
    const CCTK_REAL PDonesidedPlus2nd2At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At22[index]);
    const CCTK_REAL PDonesidedMinus2nd2At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At22[index]);
    const CCTK_REAL PDonesidedPlus2nd3At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At22[index]);
    const CCTK_REAL PDonesidedMinus2nd3At22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At22[index]);
    const CCTK_REAL PDstandard2nd1At23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At23[index]);
    const CCTK_REAL PDstandard2nd2At23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At23[index]);
    const CCTK_REAL PDstandard2nd3At23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At23[index]);
    const CCTK_REAL PDonesidedPlus2nd1At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At23[index]);
    const CCTK_REAL PDonesidedMinus2nd1At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At23[index]);
    const CCTK_REAL PDonesidedPlus2nd2At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At23[index]);
    const CCTK_REAL PDonesidedMinus2nd2At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At23[index]);
    const CCTK_REAL PDonesidedPlus2nd3At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At23[index]);
    const CCTK_REAL PDonesidedMinus2nd3At23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At23[index]);
    const CCTK_REAL PDstandard2nd1At33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&At33[index]);
    const CCTK_REAL PDstandard2nd2At33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&At33[index]);
    const CCTK_REAL PDstandard2nd3At33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&At33[index]);
    const CCTK_REAL PDonesidedPlus2nd1At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&At33[index]);
    const CCTK_REAL PDonesidedMinus2nd1At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&At33[index]);
    const CCTK_REAL PDonesidedPlus2nd2At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&At33[index]);
    const CCTK_REAL PDonesidedMinus2nd2At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&At33[index]);
    const CCTK_REAL PDonesidedPlus2nd3At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&At33[index]);
    const CCTK_REAL PDonesidedMinus2nd3At33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&At33[index]);
    const CCTK_REAL PDstandard2nd1B1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&B1[index]);
    const CCTK_REAL PDstandard2nd2B1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&B1[index]);
    const CCTK_REAL PDstandard2nd3B1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&B1[index]);
    const CCTK_REAL PDonesidedPlus2nd1B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&B1[index]);
    const CCTK_REAL PDonesidedMinus2nd1B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&B1[index]);
    const CCTK_REAL PDonesidedPlus2nd2B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&B1[index]);
    const CCTK_REAL PDonesidedMinus2nd2B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&B1[index]);
    const CCTK_REAL PDonesidedPlus2nd3B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&B1[index]);
    const CCTK_REAL PDonesidedMinus2nd3B1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&B1[index]);
    const CCTK_REAL PDstandard2nd1B2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&B2[index]);
    const CCTK_REAL PDstandard2nd2B2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&B2[index]);
    const CCTK_REAL PDstandard2nd3B2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&B2[index]);
    const CCTK_REAL PDonesidedPlus2nd1B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&B2[index]);
    const CCTK_REAL PDonesidedMinus2nd1B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&B2[index]);
    const CCTK_REAL PDonesidedPlus2nd2B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&B2[index]);
    const CCTK_REAL PDonesidedMinus2nd2B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&B2[index]);
    const CCTK_REAL PDonesidedPlus2nd3B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&B2[index]);
    const CCTK_REAL PDonesidedMinus2nd3B2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&B2[index]);
    const CCTK_REAL PDstandard2nd1B3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&B3[index]);
    const CCTK_REAL PDstandard2nd2B3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&B3[index]);
    const CCTK_REAL PDstandard2nd3B3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&B3[index]);
    const CCTK_REAL PDonesidedPlus2nd1B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&B3[index]);
    const CCTK_REAL PDonesidedMinus2nd1B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&B3[index]);
    const CCTK_REAL PDonesidedPlus2nd2B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&B3[index]);
    const CCTK_REAL PDonesidedMinus2nd2B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&B3[index]);
    const CCTK_REAL PDonesidedPlus2nd3B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&B3[index]);
    const CCTK_REAL PDonesidedMinus2nd3B3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&B3[index]);
    const CCTK_REAL PDstandard2nd1beta1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&beta1[index]);
    const CCTK_REAL PDstandard2nd2beta1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&beta1[index]);
    const CCTK_REAL PDstandard2nd3beta1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&beta1[index]);
    const CCTK_REAL PDonesidedPlus2nd1beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&beta1[index]);
    const CCTK_REAL PDonesidedMinus2nd1beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&beta1[index]);
    const CCTK_REAL PDonesidedPlus2nd2beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&beta1[index]);
    const CCTK_REAL PDonesidedMinus2nd2beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&beta1[index]);
    const CCTK_REAL PDonesidedPlus2nd3beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&beta1[index]);
    const CCTK_REAL PDonesidedMinus2nd3beta1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&beta1[index]);
    const CCTK_REAL PDstandard2nd1beta2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&beta2[index]);
    const CCTK_REAL PDstandard2nd2beta2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&beta2[index]);
    const CCTK_REAL PDstandard2nd3beta2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&beta2[index]);
    const CCTK_REAL PDonesidedPlus2nd1beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&beta2[index]);
    const CCTK_REAL PDonesidedMinus2nd1beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&beta2[index]);
    const CCTK_REAL PDonesidedPlus2nd2beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&beta2[index]);
    const CCTK_REAL PDonesidedMinus2nd2beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&beta2[index]);
    const CCTK_REAL PDonesidedPlus2nd3beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&beta2[index]);
    const CCTK_REAL PDonesidedMinus2nd3beta2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&beta2[index]);
    const CCTK_REAL PDstandard2nd1beta3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&beta3[index]);
    const CCTK_REAL PDstandard2nd2beta3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&beta3[index]);
    const CCTK_REAL PDstandard2nd3beta3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&beta3[index]);
    const CCTK_REAL PDonesidedPlus2nd1beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&beta3[index]);
    const CCTK_REAL PDonesidedMinus2nd1beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&beta3[index]);
    const CCTK_REAL PDonesidedPlus2nd2beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&beta3[index]);
    const CCTK_REAL PDonesidedMinus2nd2beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&beta3[index]);
    const CCTK_REAL PDonesidedPlus2nd3beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&beta3[index]);
    const CCTK_REAL PDonesidedMinus2nd3beta3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&beta3[index]);
    const CCTK_REAL PDstandard2nd1gt11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt11[index]);
    const CCTK_REAL PDstandard2nd2gt11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt11[index]);
    const CCTK_REAL PDstandard2nd3gt11 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt11[index]);
    const CCTK_REAL PDonesidedPlus2nd1gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt11[index]);
    const CCTK_REAL PDonesidedMinus2nd1gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt11[index]);
    const CCTK_REAL PDonesidedPlus2nd2gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt11[index]);
    const CCTK_REAL PDonesidedMinus2nd2gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt11[index]);
    const CCTK_REAL PDonesidedPlus2nd3gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt11[index]);
    const CCTK_REAL PDonesidedMinus2nd3gt11 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt11[index]);
    const CCTK_REAL PDstandard2nd1gt12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt12[index]);
    const CCTK_REAL PDstandard2nd2gt12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt12[index]);
    const CCTK_REAL PDstandard2nd3gt12 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt12[index]);
    const CCTK_REAL PDonesidedPlus2nd1gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt12[index]);
    const CCTK_REAL PDonesidedMinus2nd1gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt12[index]);
    const CCTK_REAL PDonesidedPlus2nd2gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt12[index]);
    const CCTK_REAL PDonesidedMinus2nd2gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt12[index]);
    const CCTK_REAL PDonesidedPlus2nd3gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt12[index]);
    const CCTK_REAL PDonesidedMinus2nd3gt12 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt12[index]);
    const CCTK_REAL PDstandard2nd1gt13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt13[index]);
    const CCTK_REAL PDstandard2nd2gt13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt13[index]);
    const CCTK_REAL PDstandard2nd3gt13 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt13[index]);
    const CCTK_REAL PDonesidedPlus2nd1gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt13[index]);
    const CCTK_REAL PDonesidedMinus2nd1gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt13[index]);
    const CCTK_REAL PDonesidedPlus2nd2gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt13[index]);
    const CCTK_REAL PDonesidedMinus2nd2gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt13[index]);
    const CCTK_REAL PDonesidedPlus2nd3gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt13[index]);
    const CCTK_REAL PDonesidedMinus2nd3gt13 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt13[index]);
    const CCTK_REAL PDstandard2nd1gt22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt22[index]);
    const CCTK_REAL PDstandard2nd2gt22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt22[index]);
    const CCTK_REAL PDstandard2nd3gt22 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt22[index]);
    const CCTK_REAL PDonesidedPlus2nd1gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt22[index]);
    const CCTK_REAL PDonesidedMinus2nd1gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt22[index]);
    const CCTK_REAL PDonesidedPlus2nd2gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt22[index]);
    const CCTK_REAL PDonesidedMinus2nd2gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt22[index]);
    const CCTK_REAL PDonesidedPlus2nd3gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt22[index]);
    const CCTK_REAL PDonesidedMinus2nd3gt22 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt22[index]);
    const CCTK_REAL PDstandard2nd1gt23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt23[index]);
    const CCTK_REAL PDstandard2nd2gt23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt23[index]);
    const CCTK_REAL PDstandard2nd3gt23 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt23[index]);
    const CCTK_REAL PDonesidedPlus2nd1gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt23[index]);
    const CCTK_REAL PDonesidedMinus2nd1gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt23[index]);
    const CCTK_REAL PDonesidedPlus2nd2gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt23[index]);
    const CCTK_REAL PDonesidedMinus2nd2gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt23[index]);
    const CCTK_REAL PDonesidedPlus2nd3gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt23[index]);
    const CCTK_REAL PDonesidedMinus2nd3gt23 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt23[index]);
    const CCTK_REAL PDstandard2nd1gt33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&gt33[index]);
    const CCTK_REAL PDstandard2nd2gt33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&gt33[index]);
    const CCTK_REAL PDstandard2nd3gt33 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&gt33[index]);
    const CCTK_REAL PDonesidedPlus2nd1gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&gt33[index]);
    const CCTK_REAL PDonesidedMinus2nd1gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&gt33[index]);
    const CCTK_REAL PDonesidedPlus2nd2gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&gt33[index]);
    const CCTK_REAL PDonesidedMinus2nd2gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&gt33[index]);
    const CCTK_REAL PDonesidedPlus2nd3gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&gt33[index]);
    const CCTK_REAL PDonesidedMinus2nd3gt33 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&gt33[index]);
    const CCTK_REAL PDstandard2nd1phi CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&phi[index]);
    const CCTK_REAL PDstandard2nd2phi CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&phi[index]);
    const CCTK_REAL PDstandard2nd3phi CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&phi[index]);
    const CCTK_REAL PDonesidedPlus2nd1phi CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&phi[index]);
    const CCTK_REAL PDonesidedMinus2nd1phi CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&phi[index]);
    const CCTK_REAL PDonesidedPlus2nd2phi CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&phi[index]);
    const CCTK_REAL PDonesidedMinus2nd2phi CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&phi[index]);
    const CCTK_REAL PDonesidedPlus2nd3phi CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&phi[index]);
    const CCTK_REAL PDonesidedMinus2nd3phi CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&phi[index]);
    const CCTK_REAL PDstandard2nd1trK CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&trK[index]);
    const CCTK_REAL PDstandard2nd2trK CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&trK[index]);
    const CCTK_REAL PDstandard2nd3trK CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&trK[index]);
    const CCTK_REAL PDonesidedPlus2nd1trK CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&trK[index]);
    const CCTK_REAL PDonesidedMinus2nd1trK CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&trK[index]);
    const CCTK_REAL PDonesidedPlus2nd2trK CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&trK[index]);
    const CCTK_REAL PDonesidedMinus2nd2trK CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&trK[index]);
    const CCTK_REAL PDonesidedPlus2nd3trK CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&trK[index]);
    const CCTK_REAL PDonesidedMinus2nd3trK CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&trK[index]);
    const CCTK_REAL PDstandard2nd1Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&Xt1[index]);
    const CCTK_REAL PDstandard2nd2Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&Xt1[index]);
    const CCTK_REAL PDstandard2nd3Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&Xt1[index]);
    const CCTK_REAL PDonesidedPlus2nd1Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&Xt1[index]);
    const CCTK_REAL PDonesidedMinus2nd1Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&Xt1[index]);
    const CCTK_REAL PDonesidedPlus2nd2Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&Xt1[index]);
    const CCTK_REAL PDonesidedMinus2nd2Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&Xt1[index]);
    const CCTK_REAL PDonesidedPlus2nd3Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&Xt1[index]);
    const CCTK_REAL PDonesidedMinus2nd3Xt1 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&Xt1[index]);
    const CCTK_REAL PDstandard2nd1Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&Xt2[index]);
    const CCTK_REAL PDstandard2nd2Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&Xt2[index]);
    const CCTK_REAL PDstandard2nd3Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&Xt2[index]);
    const CCTK_REAL PDonesidedPlus2nd1Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&Xt2[index]);
    const CCTK_REAL PDonesidedMinus2nd1Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&Xt2[index]);
    const CCTK_REAL PDonesidedPlus2nd2Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&Xt2[index]);
    const CCTK_REAL PDonesidedMinus2nd2Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&Xt2[index]);
    const CCTK_REAL PDonesidedPlus2nd3Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&Xt2[index]);
    const CCTK_REAL PDonesidedMinus2nd3Xt2 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&Xt2[index]);
    const CCTK_REAL PDstandard2nd1Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&Xt3[index]);
    const CCTK_REAL PDstandard2nd2Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&Xt3[index]);
    const CCTK_REAL PDstandard2nd3Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&Xt3[index]);
    const CCTK_REAL PDonesidedPlus2nd1Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd1(&Xt3[index]);
    const CCTK_REAL PDonesidedMinus2nd1Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd1(&Xt3[index]);
    const CCTK_REAL PDonesidedPlus2nd2Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd2(&Xt3[index]);
    const CCTK_REAL PDonesidedMinus2nd2Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd2(&Xt3[index]);
    const CCTK_REAL PDonesidedPlus2nd3Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedPlus2nd3(&Xt3[index]);
    const CCTK_REAL PDonesidedMinus2nd3Xt3 CCTK_ATTRIBUTE_UNUSED = PDonesidedMinus2nd3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL rn1 CCTK_ATTRIBUTE_UNUSED = -(xCopyL*INV(rCopyL));
    
    CCTK_REAL rn2 CCTK_ATTRIBUTE_UNUSED = -(yCopyL*INV(rCopyL));
    
    CCTK_REAL rn3 CCTK_ATTRIBUTE_UNUSED = -(zCopyL*INV(rCopyL));
    
    CCTK_REAL phi0 CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod,1,0);
    
    CCTK_REAL v0 CCTK_ATTRIBUTE_UNUSED = sqrt(ToReal(harmonicF));
    
    CCTK_REAL phirhsL CCTK_ATTRIBUTE_UNUSED = v0*(-phiL + phi0 + 
      rCopyL*(rn1*IfThen(normal[0] < 0,PDonesidedPlus2nd1phi,IfThen(normal[0] 
      > 0,PDonesidedMinus2nd1phi,PDstandard2nd1phi)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2phi,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2phi,PDstandard2nd2phi)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3phi,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3phi,PDstandard2nd3phi))))*INV(rCopyL);
    
    CCTK_REAL gt11rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt11,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt11,PDstandard2nd1gt11)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt11,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt11,PDstandard2nd2gt11)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt11,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt11,PDstandard2nd3gt11)) + (1 - 
      gt11L)*INV(rCopyL);
    
    CCTK_REAL gt12rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt12,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt12,PDstandard2nd1gt12)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt12,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt12,PDstandard2nd2gt12)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt12,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt12,PDstandard2nd3gt12)) - gt12L*INV(rCopyL);
    
    CCTK_REAL gt13rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt13,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt13,PDstandard2nd1gt13)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt13,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt13,PDstandard2nd2gt13)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt13,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt13,PDstandard2nd3gt13)) - gt13L*INV(rCopyL);
    
    CCTK_REAL gt22rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt22,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt22,PDstandard2nd1gt22)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt22,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt22,PDstandard2nd2gt22)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt22,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt22,PDstandard2nd3gt22)) + (1 - 
      gt22L)*INV(rCopyL);
    
    CCTK_REAL gt23rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt23,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt23,PDstandard2nd1gt23)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt23,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt23,PDstandard2nd2gt23)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt23,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt23,PDstandard2nd3gt23)) - gt23L*INV(rCopyL);
    
    CCTK_REAL gt33rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt33,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt33,PDstandard2nd1gt33)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt33,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt33,PDstandard2nd2gt33)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt33,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt33,PDstandard2nd3gt33)) + (1 - 
      gt33L)*INV(rCopyL);
    
    CCTK_REAL trKrhsL CCTK_ATTRIBUTE_UNUSED = v0*(rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1trK,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1trK,PDstandard2nd1trK)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2trK,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2trK,PDstandard2nd2trK)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3trK,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3trK,PDstandard2nd3trK)) - trKL*INV(rCopyL));
    
    CCTK_REAL At11rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At11,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At11,PDstandard2nd1At11)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At11,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At11,PDstandard2nd2At11)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At11,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At11,PDstandard2nd3At11)) - At11L*INV(rCopyL);
    
    CCTK_REAL At12rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At12,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At12,PDstandard2nd1At12)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At12,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At12,PDstandard2nd2At12)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At12,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At12,PDstandard2nd3At12)) - At12L*INV(rCopyL);
    
    CCTK_REAL At13rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At13,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At13,PDstandard2nd1At13)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At13,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At13,PDstandard2nd2At13)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At13,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At13,PDstandard2nd3At13)) - At13L*INV(rCopyL);
    
    CCTK_REAL At22rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At22,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At22,PDstandard2nd1At22)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At22,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At22,PDstandard2nd2At22)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At22,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At22,PDstandard2nd3At22)) - At22L*INV(rCopyL);
    
    CCTK_REAL At23rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At23,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At23,PDstandard2nd1At23)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At23,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At23,PDstandard2nd2At23)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At23,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At23,PDstandard2nd3At23)) - At23L*INV(rCopyL);
    
    CCTK_REAL At33rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At33,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At33,PDstandard2nd1At33)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At33,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At33,PDstandard2nd2At33)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At33,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At33,PDstandard2nd3At33)) - At33L*INV(rCopyL);
    
    CCTK_REAL Xt1rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1Xt1,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt1,PDstandard2nd1Xt1)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2Xt1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt1,PDstandard2nd2Xt1)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3Xt1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt1,PDstandard2nd3Xt1)) - Xt1L*INV(rCopyL);
    
    CCTK_REAL Xt2rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1Xt2,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt2,PDstandard2nd1Xt2)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2Xt2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt2,PDstandard2nd2Xt2)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3Xt2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt2,PDstandard2nd3Xt2)) - Xt2L*INV(rCopyL);
    
    CCTK_REAL Xt3rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1Xt3,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt3,PDstandard2nd1Xt3)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2Xt3,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt3,PDstandard2nd2Xt3)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3Xt3,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt3,PDstandard2nd3Xt3)) - Xt3L*INV(rCopyL);
    
    CCTK_REAL alpharhsL CCTK_ATTRIBUTE_UNUSED = v0*(rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1alpha,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1alpha,PDstandard2nd1alpha)) + rn2*IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2alpha,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2alpha,PDstandard2nd2alpha)) + rn3*IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3alpha,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3alpha,PDstandard2nd3alpha)) + (1 - 
      alphaL)*INV(rCopyL));
    
    CCTK_REAL ArhsL CCTK_ATTRIBUTE_UNUSED = v0*(rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1A,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1A,PDstandard2nd1A)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2A,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2A,PDstandard2nd2A)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3A,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3A,PDstandard2nd3A)) - AL*INV(rCopyL));
    
    CCTK_REAL beta1rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1beta1,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1beta1,PDstandard2nd1beta1)) + rn2*IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2beta1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta1,PDstandard2nd2beta1)) + rn3*IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3beta1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta1,PDstandard2nd3beta1)) - beta1L*INV(rCopyL);
    
    CCTK_REAL beta2rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1beta2,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1beta2,PDstandard2nd1beta2)) + rn2*IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2beta2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta2,PDstandard2nd2beta2)) + rn3*IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3beta2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta2,PDstandard2nd3beta2)) - beta2L*INV(rCopyL);
    
    CCTK_REAL beta3rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1beta3,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1beta3,PDstandard2nd1beta3)) + rn2*IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2beta3,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta3,PDstandard2nd2beta3)) + rn3*IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3beta3,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta3,PDstandard2nd3beta3)) - beta3L*INV(rCopyL);
    
    CCTK_REAL B1rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1B1,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1B1,PDstandard2nd1B1)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2B1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2B1,PDstandard2nd2B1)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3B1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3B1,PDstandard2nd3B1)) - B1L*INV(rCopyL);
    
    CCTK_REAL B2rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1B2,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1B2,PDstandard2nd1B2)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2B2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2B2,PDstandard2nd2B2)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3B2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3B2,PDstandard2nd3B2)) - B2L*INV(rCopyL);
    
    CCTK_REAL B3rhsL CCTK_ATTRIBUTE_UNUSED = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1B3,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1B3,PDstandard2nd1B3)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2B3,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2B3,PDstandard2nd2B3)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3B3,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3B3,PDstandard2nd3B3)) - B3L*INV(rCopyL);
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    Arhs[index] = ArhsL;
    At11rhs[index] = At11rhsL;
    At12rhs[index] = At12rhsL;
    At13rhs[index] = At13rhsL;
    At22rhs[index] = At22rhsL;
    At23rhs[index] = At23rhsL;
    At33rhs[index] = At33rhsL;
    B1rhs[index] = B1rhsL;
    B2rhs[index] = B2rhsL;
    B3rhs[index] = B3rhsL;
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    gt11rhs[index] = gt11rhsL;
    gt12rhs[index] = gt12rhsL;
    gt13rhs[index] = gt13rhsL;
    gt22rhs[index] = gt22rhsL;
    gt23rhs[index] = gt23rhsL;
    gt33rhs[index] = gt33rhsL;
    phirhs[index] = phirhsL;
    trKrhs[index] = trKrhsL;
    Xt1rhs[index] = Xt1rhsL;
    Xt2rhs[index] = Xt2rhsL;
    Xt3rhs[index] = Xt3rhsL;
  }
  CCTK_ENDLOOP3(HOST_ML_BSSN_RHSRadiativeBoundary);
}

extern "C" void HOST_ML_BSSN_RHSRadiativeBoundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering HOST_ML_BSSN_RHSRadiativeBoundary_Body");
  }
  
  if (cctk_iteration % HOST_ML_BSSN_RHSRadiativeBoundary_calc_every != HOST_ML_BSSN_RHSRadiativeBoundary_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN::coords",
    "ML_BSSN::ML_curv",
    "ML_BSSN::ML_curvrhs",
    "ML_BSSN::ML_dtlapse",
    "ML_BSSN::ML_dtlapserhs",
    "ML_BSSN::ML_dtshift",
    "ML_BSSN::ML_dtshiftrhs",
    "ML_BSSN::ML_Gamma",
    "ML_BSSN::ML_Gammarhs",
    "ML_BSSN::ML_lapse",
    "ML_BSSN::ML_lapserhs",
    "ML_BSSN::ML_log_confac",
    "ML_BSSN::ML_log_confacrhs",
    "ML_BSSN::ML_metric",
    "ML_BSSN::ML_metricrhs",
    "ML_BSSN::ML_shift",
    "ML_BSSN::ML_shiftrhs",
    "ML_BSSN::ML_trace_curv",
    "ML_BSSN::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "HOST_ML_BSSN_RHSRadiativeBoundary", 19, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "HOST_ML_BSSN_RHSRadiativeBoundary", 2, 2, 2);
  
  GenericFD_LoopOverBoundary(cctkGH, HOST_ML_BSSN_RHSRadiativeBoundary_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving HOST_ML_BSSN_RHSRadiativeBoundary_Body");
  }
}
