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
  
  CCTK_INT ierr = 0;
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

static void ML_BSSN_Host_RHSRadiativeBoundary_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  LC_LOOP3VEC(ML_BSSN_Host_RHSRadiativeBoundary,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
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
    CCTK_REAL_VEC gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC phiL = vec_load(phi[index]);
    CCTK_REAL_VEC rCopyL = vec_load(rCopy[index]);
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    CCTK_REAL_VEC xCopyL = vec_load(xCopy[index]);
    CCTK_REAL_VEC Xt1L = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L = vec_load(Xt3[index]);
    CCTK_REAL_VEC yCopyL = vec_load(yCopy[index]);
    CCTK_REAL_VEC zCopyL = vec_load(zCopy[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC const PDstandard2nd1A = PDstandard2nd1(&A[index]);
    CCTK_REAL_VEC const PDstandard2nd2A = PDstandard2nd2(&A[index]);
    CCTK_REAL_VEC const PDstandard2nd3A = PDstandard2nd3(&A[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1A = PDonesidedPlus2nd1(&A[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1A = PDonesidedMinus2nd1(&A[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2A = PDonesidedPlus2nd2(&A[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2A = PDonesidedMinus2nd2(&A[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3A = PDonesidedPlus2nd3(&A[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3A = PDonesidedMinus2nd3(&A[index]);
    CCTK_REAL_VEC const PDstandard2nd1alpha = PDstandard2nd1(&alpha[index]);
    CCTK_REAL_VEC const PDstandard2nd2alpha = PDstandard2nd2(&alpha[index]);
    CCTK_REAL_VEC const PDstandard2nd3alpha = PDstandard2nd3(&alpha[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1alpha = PDonesidedPlus2nd1(&alpha[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1alpha = PDonesidedMinus2nd1(&alpha[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2alpha = PDonesidedPlus2nd2(&alpha[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2alpha = PDonesidedMinus2nd2(&alpha[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3alpha = PDonesidedPlus2nd3(&alpha[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3alpha = PDonesidedMinus2nd3(&alpha[index]);
    CCTK_REAL_VEC const PDstandard2nd1At11 = PDstandard2nd1(&At11[index]);
    CCTK_REAL_VEC const PDstandard2nd2At11 = PDstandard2nd2(&At11[index]);
    CCTK_REAL_VEC const PDstandard2nd3At11 = PDstandard2nd3(&At11[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1At11 = PDonesidedPlus2nd1(&At11[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1At11 = PDonesidedMinus2nd1(&At11[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2At11 = PDonesidedPlus2nd2(&At11[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2At11 = PDonesidedMinus2nd2(&At11[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3At11 = PDonesidedPlus2nd3(&At11[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3At11 = PDonesidedMinus2nd3(&At11[index]);
    CCTK_REAL_VEC const PDstandard2nd1At12 = PDstandard2nd1(&At12[index]);
    CCTK_REAL_VEC const PDstandard2nd2At12 = PDstandard2nd2(&At12[index]);
    CCTK_REAL_VEC const PDstandard2nd3At12 = PDstandard2nd3(&At12[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1At12 = PDonesidedPlus2nd1(&At12[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1At12 = PDonesidedMinus2nd1(&At12[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2At12 = PDonesidedPlus2nd2(&At12[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2At12 = PDonesidedMinus2nd2(&At12[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3At12 = PDonesidedPlus2nd3(&At12[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3At12 = PDonesidedMinus2nd3(&At12[index]);
    CCTK_REAL_VEC const PDstandard2nd1At13 = PDstandard2nd1(&At13[index]);
    CCTK_REAL_VEC const PDstandard2nd2At13 = PDstandard2nd2(&At13[index]);
    CCTK_REAL_VEC const PDstandard2nd3At13 = PDstandard2nd3(&At13[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1At13 = PDonesidedPlus2nd1(&At13[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1At13 = PDonesidedMinus2nd1(&At13[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2At13 = PDonesidedPlus2nd2(&At13[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2At13 = PDonesidedMinus2nd2(&At13[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3At13 = PDonesidedPlus2nd3(&At13[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3At13 = PDonesidedMinus2nd3(&At13[index]);
    CCTK_REAL_VEC const PDstandard2nd1At22 = PDstandard2nd1(&At22[index]);
    CCTK_REAL_VEC const PDstandard2nd2At22 = PDstandard2nd2(&At22[index]);
    CCTK_REAL_VEC const PDstandard2nd3At22 = PDstandard2nd3(&At22[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1At22 = PDonesidedPlus2nd1(&At22[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1At22 = PDonesidedMinus2nd1(&At22[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2At22 = PDonesidedPlus2nd2(&At22[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2At22 = PDonesidedMinus2nd2(&At22[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3At22 = PDonesidedPlus2nd3(&At22[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3At22 = PDonesidedMinus2nd3(&At22[index]);
    CCTK_REAL_VEC const PDstandard2nd1At23 = PDstandard2nd1(&At23[index]);
    CCTK_REAL_VEC const PDstandard2nd2At23 = PDstandard2nd2(&At23[index]);
    CCTK_REAL_VEC const PDstandard2nd3At23 = PDstandard2nd3(&At23[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1At23 = PDonesidedPlus2nd1(&At23[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1At23 = PDonesidedMinus2nd1(&At23[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2At23 = PDonesidedPlus2nd2(&At23[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2At23 = PDonesidedMinus2nd2(&At23[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3At23 = PDonesidedPlus2nd3(&At23[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3At23 = PDonesidedMinus2nd3(&At23[index]);
    CCTK_REAL_VEC const PDstandard2nd1At33 = PDstandard2nd1(&At33[index]);
    CCTK_REAL_VEC const PDstandard2nd2At33 = PDstandard2nd2(&At33[index]);
    CCTK_REAL_VEC const PDstandard2nd3At33 = PDstandard2nd3(&At33[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1At33 = PDonesidedPlus2nd1(&At33[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1At33 = PDonesidedMinus2nd1(&At33[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2At33 = PDonesidedPlus2nd2(&At33[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2At33 = PDonesidedMinus2nd2(&At33[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3At33 = PDonesidedPlus2nd3(&At33[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3At33 = PDonesidedMinus2nd3(&At33[index]);
    CCTK_REAL_VEC const PDstandard2nd1B1 = PDstandard2nd1(&B1[index]);
    CCTK_REAL_VEC const PDstandard2nd2B1 = PDstandard2nd2(&B1[index]);
    CCTK_REAL_VEC const PDstandard2nd3B1 = PDstandard2nd3(&B1[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1B1 = PDonesidedPlus2nd1(&B1[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1B1 = PDonesidedMinus2nd1(&B1[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2B1 = PDonesidedPlus2nd2(&B1[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2B1 = PDonesidedMinus2nd2(&B1[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3B1 = PDonesidedPlus2nd3(&B1[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3B1 = PDonesidedMinus2nd3(&B1[index]);
    CCTK_REAL_VEC const PDstandard2nd1B2 = PDstandard2nd1(&B2[index]);
    CCTK_REAL_VEC const PDstandard2nd2B2 = PDstandard2nd2(&B2[index]);
    CCTK_REAL_VEC const PDstandard2nd3B2 = PDstandard2nd3(&B2[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1B2 = PDonesidedPlus2nd1(&B2[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1B2 = PDonesidedMinus2nd1(&B2[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2B2 = PDonesidedPlus2nd2(&B2[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2B2 = PDonesidedMinus2nd2(&B2[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3B2 = PDonesidedPlus2nd3(&B2[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3B2 = PDonesidedMinus2nd3(&B2[index]);
    CCTK_REAL_VEC const PDstandard2nd1B3 = PDstandard2nd1(&B3[index]);
    CCTK_REAL_VEC const PDstandard2nd2B3 = PDstandard2nd2(&B3[index]);
    CCTK_REAL_VEC const PDstandard2nd3B3 = PDstandard2nd3(&B3[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1B3 = PDonesidedPlus2nd1(&B3[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1B3 = PDonesidedMinus2nd1(&B3[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2B3 = PDonesidedPlus2nd2(&B3[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2B3 = PDonesidedMinus2nd2(&B3[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3B3 = PDonesidedPlus2nd3(&B3[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3B3 = PDonesidedMinus2nd3(&B3[index]);
    CCTK_REAL_VEC const PDstandard2nd1beta1 = PDstandard2nd1(&beta1[index]);
    CCTK_REAL_VEC const PDstandard2nd2beta1 = PDstandard2nd2(&beta1[index]);
    CCTK_REAL_VEC const PDstandard2nd3beta1 = PDstandard2nd3(&beta1[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1beta1 = PDonesidedPlus2nd1(&beta1[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1beta1 = PDonesidedMinus2nd1(&beta1[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2beta1 = PDonesidedPlus2nd2(&beta1[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2beta1 = PDonesidedMinus2nd2(&beta1[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3beta1 = PDonesidedPlus2nd3(&beta1[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3beta1 = PDonesidedMinus2nd3(&beta1[index]);
    CCTK_REAL_VEC const PDstandard2nd1beta2 = PDstandard2nd1(&beta2[index]);
    CCTK_REAL_VEC const PDstandard2nd2beta2 = PDstandard2nd2(&beta2[index]);
    CCTK_REAL_VEC const PDstandard2nd3beta2 = PDstandard2nd3(&beta2[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1beta2 = PDonesidedPlus2nd1(&beta2[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1beta2 = PDonesidedMinus2nd1(&beta2[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2beta2 = PDonesidedPlus2nd2(&beta2[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2beta2 = PDonesidedMinus2nd2(&beta2[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3beta2 = PDonesidedPlus2nd3(&beta2[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3beta2 = PDonesidedMinus2nd3(&beta2[index]);
    CCTK_REAL_VEC const PDstandard2nd1beta3 = PDstandard2nd1(&beta3[index]);
    CCTK_REAL_VEC const PDstandard2nd2beta3 = PDstandard2nd2(&beta3[index]);
    CCTK_REAL_VEC const PDstandard2nd3beta3 = PDstandard2nd3(&beta3[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1beta3 = PDonesidedPlus2nd1(&beta3[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1beta3 = PDonesidedMinus2nd1(&beta3[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2beta3 = PDonesidedPlus2nd2(&beta3[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2beta3 = PDonesidedMinus2nd2(&beta3[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3beta3 = PDonesidedPlus2nd3(&beta3[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3beta3 = PDonesidedMinus2nd3(&beta3[index]);
    CCTK_REAL_VEC const PDstandard2nd1gt11 = PDstandard2nd1(&gt11[index]);
    CCTK_REAL_VEC const PDstandard2nd2gt11 = PDstandard2nd2(&gt11[index]);
    CCTK_REAL_VEC const PDstandard2nd3gt11 = PDstandard2nd3(&gt11[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1gt11 = PDonesidedPlus2nd1(&gt11[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1gt11 = PDonesidedMinus2nd1(&gt11[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2gt11 = PDonesidedPlus2nd2(&gt11[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2gt11 = PDonesidedMinus2nd2(&gt11[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3gt11 = PDonesidedPlus2nd3(&gt11[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3gt11 = PDonesidedMinus2nd3(&gt11[index]);
    CCTK_REAL_VEC const PDstandard2nd1gt12 = PDstandard2nd1(&gt12[index]);
    CCTK_REAL_VEC const PDstandard2nd2gt12 = PDstandard2nd2(&gt12[index]);
    CCTK_REAL_VEC const PDstandard2nd3gt12 = PDstandard2nd3(&gt12[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1gt12 = PDonesidedPlus2nd1(&gt12[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1gt12 = PDonesidedMinus2nd1(&gt12[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2gt12 = PDonesidedPlus2nd2(&gt12[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2gt12 = PDonesidedMinus2nd2(&gt12[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3gt12 = PDonesidedPlus2nd3(&gt12[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3gt12 = PDonesidedMinus2nd3(&gt12[index]);
    CCTK_REAL_VEC const PDstandard2nd1gt13 = PDstandard2nd1(&gt13[index]);
    CCTK_REAL_VEC const PDstandard2nd2gt13 = PDstandard2nd2(&gt13[index]);
    CCTK_REAL_VEC const PDstandard2nd3gt13 = PDstandard2nd3(&gt13[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1gt13 = PDonesidedPlus2nd1(&gt13[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1gt13 = PDonesidedMinus2nd1(&gt13[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2gt13 = PDonesidedPlus2nd2(&gt13[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2gt13 = PDonesidedMinus2nd2(&gt13[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3gt13 = PDonesidedPlus2nd3(&gt13[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3gt13 = PDonesidedMinus2nd3(&gt13[index]);
    CCTK_REAL_VEC const PDstandard2nd1gt22 = PDstandard2nd1(&gt22[index]);
    CCTK_REAL_VEC const PDstandard2nd2gt22 = PDstandard2nd2(&gt22[index]);
    CCTK_REAL_VEC const PDstandard2nd3gt22 = PDstandard2nd3(&gt22[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1gt22 = PDonesidedPlus2nd1(&gt22[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1gt22 = PDonesidedMinus2nd1(&gt22[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2gt22 = PDonesidedPlus2nd2(&gt22[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2gt22 = PDonesidedMinus2nd2(&gt22[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3gt22 = PDonesidedPlus2nd3(&gt22[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3gt22 = PDonesidedMinus2nd3(&gt22[index]);
    CCTK_REAL_VEC const PDstandard2nd1gt23 = PDstandard2nd1(&gt23[index]);
    CCTK_REAL_VEC const PDstandard2nd2gt23 = PDstandard2nd2(&gt23[index]);
    CCTK_REAL_VEC const PDstandard2nd3gt23 = PDstandard2nd3(&gt23[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1gt23 = PDonesidedPlus2nd1(&gt23[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1gt23 = PDonesidedMinus2nd1(&gt23[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2gt23 = PDonesidedPlus2nd2(&gt23[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2gt23 = PDonesidedMinus2nd2(&gt23[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3gt23 = PDonesidedPlus2nd3(&gt23[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3gt23 = PDonesidedMinus2nd3(&gt23[index]);
    CCTK_REAL_VEC const PDstandard2nd1gt33 = PDstandard2nd1(&gt33[index]);
    CCTK_REAL_VEC const PDstandard2nd2gt33 = PDstandard2nd2(&gt33[index]);
    CCTK_REAL_VEC const PDstandard2nd3gt33 = PDstandard2nd3(&gt33[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1gt33 = PDonesidedPlus2nd1(&gt33[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1gt33 = PDonesidedMinus2nd1(&gt33[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2gt33 = PDonesidedPlus2nd2(&gt33[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2gt33 = PDonesidedMinus2nd2(&gt33[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3gt33 = PDonesidedPlus2nd3(&gt33[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3gt33 = PDonesidedMinus2nd3(&gt33[index]);
    CCTK_REAL_VEC const PDstandard2nd1phi = PDstandard2nd1(&phi[index]);
    CCTK_REAL_VEC const PDstandard2nd2phi = PDstandard2nd2(&phi[index]);
    CCTK_REAL_VEC const PDstandard2nd3phi = PDstandard2nd3(&phi[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1phi = PDonesidedPlus2nd1(&phi[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1phi = PDonesidedMinus2nd1(&phi[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2phi = PDonesidedPlus2nd2(&phi[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2phi = PDonesidedMinus2nd2(&phi[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3phi = PDonesidedPlus2nd3(&phi[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3phi = PDonesidedMinus2nd3(&phi[index]);
    CCTK_REAL_VEC const PDstandard2nd1trK = PDstandard2nd1(&trK[index]);
    CCTK_REAL_VEC const PDstandard2nd2trK = PDstandard2nd2(&trK[index]);
    CCTK_REAL_VEC const PDstandard2nd3trK = PDstandard2nd3(&trK[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1trK = PDonesidedPlus2nd1(&trK[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1trK = PDonesidedMinus2nd1(&trK[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2trK = PDonesidedPlus2nd2(&trK[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2trK = PDonesidedMinus2nd2(&trK[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3trK = PDonesidedPlus2nd3(&trK[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3trK = PDonesidedMinus2nd3(&trK[index]);
    CCTK_REAL_VEC const PDstandard2nd1Xt1 = PDstandard2nd1(&Xt1[index]);
    CCTK_REAL_VEC const PDstandard2nd2Xt1 = PDstandard2nd2(&Xt1[index]);
    CCTK_REAL_VEC const PDstandard2nd3Xt1 = PDstandard2nd3(&Xt1[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1Xt1 = PDonesidedPlus2nd1(&Xt1[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1Xt1 = PDonesidedMinus2nd1(&Xt1[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2Xt1 = PDonesidedPlus2nd2(&Xt1[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2Xt1 = PDonesidedMinus2nd2(&Xt1[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3Xt1 = PDonesidedPlus2nd3(&Xt1[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3Xt1 = PDonesidedMinus2nd3(&Xt1[index]);
    CCTK_REAL_VEC const PDstandard2nd1Xt2 = PDstandard2nd1(&Xt2[index]);
    CCTK_REAL_VEC const PDstandard2nd2Xt2 = PDstandard2nd2(&Xt2[index]);
    CCTK_REAL_VEC const PDstandard2nd3Xt2 = PDstandard2nd3(&Xt2[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1Xt2 = PDonesidedPlus2nd1(&Xt2[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1Xt2 = PDonesidedMinus2nd1(&Xt2[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2Xt2 = PDonesidedPlus2nd2(&Xt2[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2Xt2 = PDonesidedMinus2nd2(&Xt2[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3Xt2 = PDonesidedPlus2nd3(&Xt2[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3Xt2 = PDonesidedMinus2nd3(&Xt2[index]);
    CCTK_REAL_VEC const PDstandard2nd1Xt3 = PDstandard2nd1(&Xt3[index]);
    CCTK_REAL_VEC const PDstandard2nd2Xt3 = PDstandard2nd2(&Xt3[index]);
    CCTK_REAL_VEC const PDstandard2nd3Xt3 = PDstandard2nd3(&Xt3[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd1Xt3 = PDonesidedPlus2nd1(&Xt3[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd1Xt3 = PDonesidedMinus2nd1(&Xt3[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd2Xt3 = PDonesidedPlus2nd2(&Xt3[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd2Xt3 = PDonesidedMinus2nd2(&Xt3[index]);
    CCTK_REAL_VEC const PDonesidedPlus2nd3Xt3 = PDonesidedPlus2nd3(&Xt3[index]);
    CCTK_REAL_VEC const PDonesidedMinus2nd3Xt3 = PDonesidedMinus2nd3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC rn1 = kneg(kdiv(xCopyL,rCopyL));
    
    CCTK_REAL_VEC rn2 = kneg(kdiv(yCopyL,rCopyL));
    
    CCTK_REAL_VEC rn3 = kneg(kdiv(zCopyL,rCopyL));
    
    CCTK_REAL_VEC phi0 = IfThen(conformalMethod,ToReal(1),ToReal(0));
    
    CCTK_REAL_VEC v0 = ksqrt(ToReal(harmonicF));
    
    CCTK_REAL_VEC phirhsL = 
      kdiv(kmul(v0,kadd(phi0,kmsub(rCopyL,kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1phi,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1phi,PDstandard2nd1phi)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2phi,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2phi,PDstandard2nd2phi)),kmul(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3phi,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3phi,PDstandard2nd3phi))))),phiL))),rCopyL);
    
    CCTK_REAL_VEC gt11rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt11,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt11,PDstandard2nd1gt11)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt11,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt11,PDstandard2nd2gt11)),kmadd(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt11,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt11,PDstandard2nd3gt11)),kdiv(ksub(ToReal(1),gt11L),rCopyL))));
    
    CCTK_REAL_VEC gt12rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt12,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt12,PDstandard2nd1gt12)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt12,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt12,PDstandard2nd2gt12)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt12,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt12,PDstandard2nd3gt12)),kdiv(gt12L,rCopyL))));
    
    CCTK_REAL_VEC gt13rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt13,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt13,PDstandard2nd1gt13)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt13,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt13,PDstandard2nd2gt13)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt13,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt13,PDstandard2nd3gt13)),kdiv(gt13L,rCopyL))));
    
    CCTK_REAL_VEC gt22rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt22,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt22,PDstandard2nd1gt22)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt22,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt22,PDstandard2nd2gt22)),kmadd(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt22,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt22,PDstandard2nd3gt22)),kdiv(ksub(ToReal(1),gt22L),rCopyL))));
    
    CCTK_REAL_VEC gt23rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt23,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt23,PDstandard2nd1gt23)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt23,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt23,PDstandard2nd2gt23)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt23,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt23,PDstandard2nd3gt23)),kdiv(gt23L,rCopyL))));
    
    CCTK_REAL_VEC gt33rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt33,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt33,PDstandard2nd1gt33)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2gt33,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt33,PDstandard2nd2gt33)),kmadd(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3gt33,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt33,PDstandard2nd3gt33)),kdiv(ksub(ToReal(1),gt33L),rCopyL))));
    
    CCTK_REAL_VEC trKrhsL = kmul(v0,kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1trK,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1trK,PDstandard2nd1trK)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2trK,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2trK,PDstandard2nd2trK)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3trK,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3trK,PDstandard2nd3trK)),kdiv(trKL,rCopyL)))));
    
    CCTK_REAL_VEC At11rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At11,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At11,PDstandard2nd1At11)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At11,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At11,PDstandard2nd2At11)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At11,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At11,PDstandard2nd3At11)),kdiv(At11L,rCopyL))));
    
    CCTK_REAL_VEC At12rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At12,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At12,PDstandard2nd1At12)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At12,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At12,PDstandard2nd2At12)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At12,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At12,PDstandard2nd3At12)),kdiv(At12L,rCopyL))));
    
    CCTK_REAL_VEC At13rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At13,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At13,PDstandard2nd1At13)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At13,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At13,PDstandard2nd2At13)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At13,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At13,PDstandard2nd3At13)),kdiv(At13L,rCopyL))));
    
    CCTK_REAL_VEC At22rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At22,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At22,PDstandard2nd1At22)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At22,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At22,PDstandard2nd2At22)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At22,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At22,PDstandard2nd3At22)),kdiv(At22L,rCopyL))));
    
    CCTK_REAL_VEC At23rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At23,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At23,PDstandard2nd1At23)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At23,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At23,PDstandard2nd2At23)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At23,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At23,PDstandard2nd3At23)),kdiv(At23L,rCopyL))));
    
    CCTK_REAL_VEC At33rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At33,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At33,PDstandard2nd1At33)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2At33,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At33,PDstandard2nd2At33)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3At33,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At33,PDstandard2nd3At33)),kdiv(At33L,rCopyL))));
    
    CCTK_REAL_VEC Xt1rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1Xt1,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt1,PDstandard2nd1Xt1)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2Xt1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt1,PDstandard2nd2Xt1)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3Xt1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt1,PDstandard2nd3Xt1)),kdiv(Xt1L,rCopyL))));
    
    CCTK_REAL_VEC Xt2rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1Xt2,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt2,PDstandard2nd1Xt2)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2Xt2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt2,PDstandard2nd2Xt2)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3Xt2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt2,PDstandard2nd3Xt2)),kdiv(Xt2L,rCopyL))));
    
    CCTK_REAL_VEC Xt3rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1Xt3,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt3,PDstandard2nd1Xt3)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2Xt3,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt3,PDstandard2nd2Xt3)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3Xt3,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt3,PDstandard2nd3Xt3)),kdiv(Xt3L,rCopyL))));
    
    CCTK_REAL_VEC alpharhsL = kmul(v0,kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1alpha,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1alpha,PDstandard2nd1alpha)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2alpha,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2alpha,PDstandard2nd2alpha)),kmadd(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3alpha,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3alpha,PDstandard2nd3alpha)),kdiv(ksub(ToReal(1),alphaL),rCopyL)))));
    
    CCTK_REAL_VEC ArhsL = kmul(v0,kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1A,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1A,PDstandard2nd1A)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2A,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2A,PDstandard2nd2A)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3A,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3A,PDstandard2nd3A)),kdiv(AL,rCopyL)))));
    
    CCTK_REAL_VEC beta1rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1beta1,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1beta1,PDstandard2nd1beta1)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2beta1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta1,PDstandard2nd2beta1)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3beta1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta1,PDstandard2nd3beta1)),kdiv(beta1L,rCopyL))));
    
    CCTK_REAL_VEC beta2rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1beta2,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1beta2,PDstandard2nd1beta2)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2beta2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta2,PDstandard2nd2beta2)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3beta2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta2,PDstandard2nd3beta2)),kdiv(beta2L,rCopyL))));
    
    CCTK_REAL_VEC beta3rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1beta3,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1beta3,PDstandard2nd1beta3)),kmadd(rn2,IfThen(normal[1] 
      < 0,PDonesidedPlus2nd2beta3,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta3,PDstandard2nd2beta3)),kmsub(rn3,IfThen(normal[2] 
      < 0,PDonesidedPlus2nd3beta3,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta3,PDstandard2nd3beta3)),kdiv(beta3L,rCopyL))));
    
    CCTK_REAL_VEC B1rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1B1,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1B1,PDstandard2nd1B1)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2B1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2B1,PDstandard2nd2B1)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3B1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3B1,PDstandard2nd3B1)),kdiv(B1L,rCopyL))));
    
    CCTK_REAL_VEC B2rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1B2,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1B2,PDstandard2nd1B2)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2B2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2B2,PDstandard2nd2B2)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3B2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3B2,PDstandard2nd3B2)),kdiv(B2L,rCopyL))));
    
    CCTK_REAL_VEC B3rhsL = kmadd(rn1,IfThen(normal[0] < 
      0,PDonesidedPlus2nd1B3,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1B3,PDstandard2nd1B3)),kmadd(rn2,IfThen(normal[1] < 
      0,PDonesidedPlus2nd2B3,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2B3,PDstandard2nd2B3)),kmsub(rn3,IfThen(normal[2] < 
      0,PDonesidedPlus2nd3B3,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3B3,PDstandard2nd3B3)),kdiv(B3L,rCopyL))));
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,lc_imin,lc_imax);
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
  LC_ENDLOOP3VEC(ML_BSSN_Host_RHSRadiativeBoundary);
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
  
  const char *const groups[] = {
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
