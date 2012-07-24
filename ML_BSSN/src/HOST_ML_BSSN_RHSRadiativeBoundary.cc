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
  
  CCTK_INT ierr = 0;
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

static void HOST_ML_BSSN_RHSRadiativeBoundary_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_REAL const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL const dz = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL const dt = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL const t = ToReal(cctk_time);
  CCTK_REAL const dxi = INV(dx);
  CCTK_REAL const dyi = INV(dy);
  CCTK_REAL const dzi = INV(dz);
  CCTK_REAL const khalf = 0.5;
  CCTK_REAL const kthird = 1/3.0;
  CCTK_REAL const ktwothird = 2.0/3.0;
  CCTK_REAL const kfourthird = 4.0/3.0;
  CCTK_REAL const keightthird = 8.0/3.0;
  CCTK_REAL const hdxi = 0.5 * dxi;
  CCTK_REAL const hdyi = 0.5 * dyi;
  CCTK_REAL const hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  CCTK_REAL const p1o1024dx = 0.0009765625*INV(dx);
  CCTK_REAL const p1o1024dy = 0.0009765625*INV(dy);
  CCTK_REAL const p1o1024dz = 0.0009765625*INV(dz);
  CCTK_REAL const p1o1680dx = 0.000595238095238095238095238095238*INV(dx);
  CCTK_REAL const p1o1680dy = 0.000595238095238095238095238095238*INV(dy);
  CCTK_REAL const p1o1680dz = 0.000595238095238095238095238095238*INV(dz);
  CCTK_REAL const p1o2dx = 0.5*INV(dx);
  CCTK_REAL const p1o2dy = 0.5*INV(dy);
  CCTK_REAL const p1o2dz = 0.5*INV(dz);
  CCTK_REAL const p1o5040dx2 = 0.000198412698412698412698412698413*INV(SQR(dx));
  CCTK_REAL const p1o5040dy2 = 0.000198412698412698412698412698413*INV(SQR(dy));
  CCTK_REAL const p1o5040dz2 = 0.000198412698412698412698412698413*INV(SQR(dz));
  CCTK_REAL const p1o560dx = 0.00178571428571428571428571428571*INV(dx);
  CCTK_REAL const p1o560dy = 0.00178571428571428571428571428571*INV(dy);
  CCTK_REAL const p1o560dz = 0.00178571428571428571428571428571*INV(dz);
  CCTK_REAL const p1o705600dxdy = 1.41723356009070294784580498866e-6*INV(dx*dy);
  CCTK_REAL const p1o705600dxdz = 1.41723356009070294784580498866e-6*INV(dx*dz);
  CCTK_REAL const p1o705600dydz = 1.41723356009070294784580498866e-6*INV(dy*dz);
  CCTK_REAL const p1o840dx = 0.00119047619047619047619047619048*INV(dx);
  CCTK_REAL const p1o840dy = 0.00119047619047619047619047619048*INV(dy);
  CCTK_REAL const p1o840dz = 0.00119047619047619047619047619048*INV(dz);
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o2dx = -0.5*INV(dx);
  CCTK_REAL const pm1o2dy = -0.5*INV(dy);
  CCTK_REAL const pm1o2dz = -0.5*INV(dz);
  CCTK_REAL const pm1o840dx = -0.00119047619047619047619047619048*INV(dx);
  CCTK_REAL const pm1o840dy = -0.00119047619047619047619047619048*INV(dy);
  CCTK_REAL const pm1o840dz = -0.00119047619047619047619047619048*INV(dz);
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  CCTK_LOOP3(HOST_ML_BSSN_RHSRadiativeBoundary,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL AL = A[index];
    CCTK_REAL alphaL = alpha[index];
    CCTK_REAL At11L = At11[index];
    CCTK_REAL At12L = At12[index];
    CCTK_REAL At13L = At13[index];
    CCTK_REAL At22L = At22[index];
    CCTK_REAL At23L = At23[index];
    CCTK_REAL At33L = At33[index];
    CCTK_REAL B1L = B1[index];
    CCTK_REAL B2L = B2[index];
    CCTK_REAL B3L = B3[index];
    CCTK_REAL beta1L = beta1[index];
    CCTK_REAL beta2L = beta2[index];
    CCTK_REAL beta3L = beta3[index];
    CCTK_REAL gt11L = gt11[index];
    CCTK_REAL gt12L = gt12[index];
    CCTK_REAL gt13L = gt13[index];
    CCTK_REAL gt22L = gt22[index];
    CCTK_REAL gt23L = gt23[index];
    CCTK_REAL gt33L = gt33[index];
    CCTK_REAL phiL = phi[index];
    CCTK_REAL rCopyL = rCopy[index];
    CCTK_REAL trKL = trK[index];
    CCTK_REAL xCopyL = xCopy[index];
    CCTK_REAL Xt1L = Xt1[index];
    CCTK_REAL Xt2L = Xt2[index];
    CCTK_REAL Xt3L = Xt3[index];
    CCTK_REAL yCopyL = yCopy[index];
    CCTK_REAL zCopyL = zCopy[index];
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDstandard2nd1A = PDstandard2nd1(&A[index]);
    CCTK_REAL const PDstandard2nd2A = PDstandard2nd2(&A[index]);
    CCTK_REAL const PDstandard2nd3A = PDstandard2nd3(&A[index]);
    CCTK_REAL const PDonesidedPlus2nd1A = PDonesidedPlus2nd1(&A[index]);
    CCTK_REAL const PDonesidedMinus2nd1A = PDonesidedMinus2nd1(&A[index]);
    CCTK_REAL const PDonesidedPlus2nd2A = PDonesidedPlus2nd2(&A[index]);
    CCTK_REAL const PDonesidedMinus2nd2A = PDonesidedMinus2nd2(&A[index]);
    CCTK_REAL const PDonesidedPlus2nd3A = PDonesidedPlus2nd3(&A[index]);
    CCTK_REAL const PDonesidedMinus2nd3A = PDonesidedMinus2nd3(&A[index]);
    CCTK_REAL const PDstandard2nd1alpha = PDstandard2nd1(&alpha[index]);
    CCTK_REAL const PDstandard2nd2alpha = PDstandard2nd2(&alpha[index]);
    CCTK_REAL const PDstandard2nd3alpha = PDstandard2nd3(&alpha[index]);
    CCTK_REAL const PDonesidedPlus2nd1alpha = PDonesidedPlus2nd1(&alpha[index]);
    CCTK_REAL const PDonesidedMinus2nd1alpha = PDonesidedMinus2nd1(&alpha[index]);
    CCTK_REAL const PDonesidedPlus2nd2alpha = PDonesidedPlus2nd2(&alpha[index]);
    CCTK_REAL const PDonesidedMinus2nd2alpha = PDonesidedMinus2nd2(&alpha[index]);
    CCTK_REAL const PDonesidedPlus2nd3alpha = PDonesidedPlus2nd3(&alpha[index]);
    CCTK_REAL const PDonesidedMinus2nd3alpha = PDonesidedMinus2nd3(&alpha[index]);
    CCTK_REAL const PDstandard2nd1At11 = PDstandard2nd1(&At11[index]);
    CCTK_REAL const PDstandard2nd2At11 = PDstandard2nd2(&At11[index]);
    CCTK_REAL const PDstandard2nd3At11 = PDstandard2nd3(&At11[index]);
    CCTK_REAL const PDonesidedPlus2nd1At11 = PDonesidedPlus2nd1(&At11[index]);
    CCTK_REAL const PDonesidedMinus2nd1At11 = PDonesidedMinus2nd1(&At11[index]);
    CCTK_REAL const PDonesidedPlus2nd2At11 = PDonesidedPlus2nd2(&At11[index]);
    CCTK_REAL const PDonesidedMinus2nd2At11 = PDonesidedMinus2nd2(&At11[index]);
    CCTK_REAL const PDonesidedPlus2nd3At11 = PDonesidedPlus2nd3(&At11[index]);
    CCTK_REAL const PDonesidedMinus2nd3At11 = PDonesidedMinus2nd3(&At11[index]);
    CCTK_REAL const PDstandard2nd1At12 = PDstandard2nd1(&At12[index]);
    CCTK_REAL const PDstandard2nd2At12 = PDstandard2nd2(&At12[index]);
    CCTK_REAL const PDstandard2nd3At12 = PDstandard2nd3(&At12[index]);
    CCTK_REAL const PDonesidedPlus2nd1At12 = PDonesidedPlus2nd1(&At12[index]);
    CCTK_REAL const PDonesidedMinus2nd1At12 = PDonesidedMinus2nd1(&At12[index]);
    CCTK_REAL const PDonesidedPlus2nd2At12 = PDonesidedPlus2nd2(&At12[index]);
    CCTK_REAL const PDonesidedMinus2nd2At12 = PDonesidedMinus2nd2(&At12[index]);
    CCTK_REAL const PDonesidedPlus2nd3At12 = PDonesidedPlus2nd3(&At12[index]);
    CCTK_REAL const PDonesidedMinus2nd3At12 = PDonesidedMinus2nd3(&At12[index]);
    CCTK_REAL const PDstandard2nd1At13 = PDstandard2nd1(&At13[index]);
    CCTK_REAL const PDstandard2nd2At13 = PDstandard2nd2(&At13[index]);
    CCTK_REAL const PDstandard2nd3At13 = PDstandard2nd3(&At13[index]);
    CCTK_REAL const PDonesidedPlus2nd1At13 = PDonesidedPlus2nd1(&At13[index]);
    CCTK_REAL const PDonesidedMinus2nd1At13 = PDonesidedMinus2nd1(&At13[index]);
    CCTK_REAL const PDonesidedPlus2nd2At13 = PDonesidedPlus2nd2(&At13[index]);
    CCTK_REAL const PDonesidedMinus2nd2At13 = PDonesidedMinus2nd2(&At13[index]);
    CCTK_REAL const PDonesidedPlus2nd3At13 = PDonesidedPlus2nd3(&At13[index]);
    CCTK_REAL const PDonesidedMinus2nd3At13 = PDonesidedMinus2nd3(&At13[index]);
    CCTK_REAL const PDstandard2nd1At22 = PDstandard2nd1(&At22[index]);
    CCTK_REAL const PDstandard2nd2At22 = PDstandard2nd2(&At22[index]);
    CCTK_REAL const PDstandard2nd3At22 = PDstandard2nd3(&At22[index]);
    CCTK_REAL const PDonesidedPlus2nd1At22 = PDonesidedPlus2nd1(&At22[index]);
    CCTK_REAL const PDonesidedMinus2nd1At22 = PDonesidedMinus2nd1(&At22[index]);
    CCTK_REAL const PDonesidedPlus2nd2At22 = PDonesidedPlus2nd2(&At22[index]);
    CCTK_REAL const PDonesidedMinus2nd2At22 = PDonesidedMinus2nd2(&At22[index]);
    CCTK_REAL const PDonesidedPlus2nd3At22 = PDonesidedPlus2nd3(&At22[index]);
    CCTK_REAL const PDonesidedMinus2nd3At22 = PDonesidedMinus2nd3(&At22[index]);
    CCTK_REAL const PDstandard2nd1At23 = PDstandard2nd1(&At23[index]);
    CCTK_REAL const PDstandard2nd2At23 = PDstandard2nd2(&At23[index]);
    CCTK_REAL const PDstandard2nd3At23 = PDstandard2nd3(&At23[index]);
    CCTK_REAL const PDonesidedPlus2nd1At23 = PDonesidedPlus2nd1(&At23[index]);
    CCTK_REAL const PDonesidedMinus2nd1At23 = PDonesidedMinus2nd1(&At23[index]);
    CCTK_REAL const PDonesidedPlus2nd2At23 = PDonesidedPlus2nd2(&At23[index]);
    CCTK_REAL const PDonesidedMinus2nd2At23 = PDonesidedMinus2nd2(&At23[index]);
    CCTK_REAL const PDonesidedPlus2nd3At23 = PDonesidedPlus2nd3(&At23[index]);
    CCTK_REAL const PDonesidedMinus2nd3At23 = PDonesidedMinus2nd3(&At23[index]);
    CCTK_REAL const PDstandard2nd1At33 = PDstandard2nd1(&At33[index]);
    CCTK_REAL const PDstandard2nd2At33 = PDstandard2nd2(&At33[index]);
    CCTK_REAL const PDstandard2nd3At33 = PDstandard2nd3(&At33[index]);
    CCTK_REAL const PDonesidedPlus2nd1At33 = PDonesidedPlus2nd1(&At33[index]);
    CCTK_REAL const PDonesidedMinus2nd1At33 = PDonesidedMinus2nd1(&At33[index]);
    CCTK_REAL const PDonesidedPlus2nd2At33 = PDonesidedPlus2nd2(&At33[index]);
    CCTK_REAL const PDonesidedMinus2nd2At33 = PDonesidedMinus2nd2(&At33[index]);
    CCTK_REAL const PDonesidedPlus2nd3At33 = PDonesidedPlus2nd3(&At33[index]);
    CCTK_REAL const PDonesidedMinus2nd3At33 = PDonesidedMinus2nd3(&At33[index]);
    CCTK_REAL const PDstandard2nd1B1 = PDstandard2nd1(&B1[index]);
    CCTK_REAL const PDstandard2nd2B1 = PDstandard2nd2(&B1[index]);
    CCTK_REAL const PDstandard2nd3B1 = PDstandard2nd3(&B1[index]);
    CCTK_REAL const PDonesidedPlus2nd1B1 = PDonesidedPlus2nd1(&B1[index]);
    CCTK_REAL const PDonesidedMinus2nd1B1 = PDonesidedMinus2nd1(&B1[index]);
    CCTK_REAL const PDonesidedPlus2nd2B1 = PDonesidedPlus2nd2(&B1[index]);
    CCTK_REAL const PDonesidedMinus2nd2B1 = PDonesidedMinus2nd2(&B1[index]);
    CCTK_REAL const PDonesidedPlus2nd3B1 = PDonesidedPlus2nd3(&B1[index]);
    CCTK_REAL const PDonesidedMinus2nd3B1 = PDonesidedMinus2nd3(&B1[index]);
    CCTK_REAL const PDstandard2nd1B2 = PDstandard2nd1(&B2[index]);
    CCTK_REAL const PDstandard2nd2B2 = PDstandard2nd2(&B2[index]);
    CCTK_REAL const PDstandard2nd3B2 = PDstandard2nd3(&B2[index]);
    CCTK_REAL const PDonesidedPlus2nd1B2 = PDonesidedPlus2nd1(&B2[index]);
    CCTK_REAL const PDonesidedMinus2nd1B2 = PDonesidedMinus2nd1(&B2[index]);
    CCTK_REAL const PDonesidedPlus2nd2B2 = PDonesidedPlus2nd2(&B2[index]);
    CCTK_REAL const PDonesidedMinus2nd2B2 = PDonesidedMinus2nd2(&B2[index]);
    CCTK_REAL const PDonesidedPlus2nd3B2 = PDonesidedPlus2nd3(&B2[index]);
    CCTK_REAL const PDonesidedMinus2nd3B2 = PDonesidedMinus2nd3(&B2[index]);
    CCTK_REAL const PDstandard2nd1B3 = PDstandard2nd1(&B3[index]);
    CCTK_REAL const PDstandard2nd2B3 = PDstandard2nd2(&B3[index]);
    CCTK_REAL const PDstandard2nd3B3 = PDstandard2nd3(&B3[index]);
    CCTK_REAL const PDonesidedPlus2nd1B3 = PDonesidedPlus2nd1(&B3[index]);
    CCTK_REAL const PDonesidedMinus2nd1B3 = PDonesidedMinus2nd1(&B3[index]);
    CCTK_REAL const PDonesidedPlus2nd2B3 = PDonesidedPlus2nd2(&B3[index]);
    CCTK_REAL const PDonesidedMinus2nd2B3 = PDonesidedMinus2nd2(&B3[index]);
    CCTK_REAL const PDonesidedPlus2nd3B3 = PDonesidedPlus2nd3(&B3[index]);
    CCTK_REAL const PDonesidedMinus2nd3B3 = PDonesidedMinus2nd3(&B3[index]);
    CCTK_REAL const PDstandard2nd1beta1 = PDstandard2nd1(&beta1[index]);
    CCTK_REAL const PDstandard2nd2beta1 = PDstandard2nd2(&beta1[index]);
    CCTK_REAL const PDstandard2nd3beta1 = PDstandard2nd3(&beta1[index]);
    CCTK_REAL const PDonesidedPlus2nd1beta1 = PDonesidedPlus2nd1(&beta1[index]);
    CCTK_REAL const PDonesidedMinus2nd1beta1 = PDonesidedMinus2nd1(&beta1[index]);
    CCTK_REAL const PDonesidedPlus2nd2beta1 = PDonesidedPlus2nd2(&beta1[index]);
    CCTK_REAL const PDonesidedMinus2nd2beta1 = PDonesidedMinus2nd2(&beta1[index]);
    CCTK_REAL const PDonesidedPlus2nd3beta1 = PDonesidedPlus2nd3(&beta1[index]);
    CCTK_REAL const PDonesidedMinus2nd3beta1 = PDonesidedMinus2nd3(&beta1[index]);
    CCTK_REAL const PDstandard2nd1beta2 = PDstandard2nd1(&beta2[index]);
    CCTK_REAL const PDstandard2nd2beta2 = PDstandard2nd2(&beta2[index]);
    CCTK_REAL const PDstandard2nd3beta2 = PDstandard2nd3(&beta2[index]);
    CCTK_REAL const PDonesidedPlus2nd1beta2 = PDonesidedPlus2nd1(&beta2[index]);
    CCTK_REAL const PDonesidedMinus2nd1beta2 = PDonesidedMinus2nd1(&beta2[index]);
    CCTK_REAL const PDonesidedPlus2nd2beta2 = PDonesidedPlus2nd2(&beta2[index]);
    CCTK_REAL const PDonesidedMinus2nd2beta2 = PDonesidedMinus2nd2(&beta2[index]);
    CCTK_REAL const PDonesidedPlus2nd3beta2 = PDonesidedPlus2nd3(&beta2[index]);
    CCTK_REAL const PDonesidedMinus2nd3beta2 = PDonesidedMinus2nd3(&beta2[index]);
    CCTK_REAL const PDstandard2nd1beta3 = PDstandard2nd1(&beta3[index]);
    CCTK_REAL const PDstandard2nd2beta3 = PDstandard2nd2(&beta3[index]);
    CCTK_REAL const PDstandard2nd3beta3 = PDstandard2nd3(&beta3[index]);
    CCTK_REAL const PDonesidedPlus2nd1beta3 = PDonesidedPlus2nd1(&beta3[index]);
    CCTK_REAL const PDonesidedMinus2nd1beta3 = PDonesidedMinus2nd1(&beta3[index]);
    CCTK_REAL const PDonesidedPlus2nd2beta3 = PDonesidedPlus2nd2(&beta3[index]);
    CCTK_REAL const PDonesidedMinus2nd2beta3 = PDonesidedMinus2nd2(&beta3[index]);
    CCTK_REAL const PDonesidedPlus2nd3beta3 = PDonesidedPlus2nd3(&beta3[index]);
    CCTK_REAL const PDonesidedMinus2nd3beta3 = PDonesidedMinus2nd3(&beta3[index]);
    CCTK_REAL const PDstandard2nd1gt11 = PDstandard2nd1(&gt11[index]);
    CCTK_REAL const PDstandard2nd2gt11 = PDstandard2nd2(&gt11[index]);
    CCTK_REAL const PDstandard2nd3gt11 = PDstandard2nd3(&gt11[index]);
    CCTK_REAL const PDonesidedPlus2nd1gt11 = PDonesidedPlus2nd1(&gt11[index]);
    CCTK_REAL const PDonesidedMinus2nd1gt11 = PDonesidedMinus2nd1(&gt11[index]);
    CCTK_REAL const PDonesidedPlus2nd2gt11 = PDonesidedPlus2nd2(&gt11[index]);
    CCTK_REAL const PDonesidedMinus2nd2gt11 = PDonesidedMinus2nd2(&gt11[index]);
    CCTK_REAL const PDonesidedPlus2nd3gt11 = PDonesidedPlus2nd3(&gt11[index]);
    CCTK_REAL const PDonesidedMinus2nd3gt11 = PDonesidedMinus2nd3(&gt11[index]);
    CCTK_REAL const PDstandard2nd1gt12 = PDstandard2nd1(&gt12[index]);
    CCTK_REAL const PDstandard2nd2gt12 = PDstandard2nd2(&gt12[index]);
    CCTK_REAL const PDstandard2nd3gt12 = PDstandard2nd3(&gt12[index]);
    CCTK_REAL const PDonesidedPlus2nd1gt12 = PDonesidedPlus2nd1(&gt12[index]);
    CCTK_REAL const PDonesidedMinus2nd1gt12 = PDonesidedMinus2nd1(&gt12[index]);
    CCTK_REAL const PDonesidedPlus2nd2gt12 = PDonesidedPlus2nd2(&gt12[index]);
    CCTK_REAL const PDonesidedMinus2nd2gt12 = PDonesidedMinus2nd2(&gt12[index]);
    CCTK_REAL const PDonesidedPlus2nd3gt12 = PDonesidedPlus2nd3(&gt12[index]);
    CCTK_REAL const PDonesidedMinus2nd3gt12 = PDonesidedMinus2nd3(&gt12[index]);
    CCTK_REAL const PDstandard2nd1gt13 = PDstandard2nd1(&gt13[index]);
    CCTK_REAL const PDstandard2nd2gt13 = PDstandard2nd2(&gt13[index]);
    CCTK_REAL const PDstandard2nd3gt13 = PDstandard2nd3(&gt13[index]);
    CCTK_REAL const PDonesidedPlus2nd1gt13 = PDonesidedPlus2nd1(&gt13[index]);
    CCTK_REAL const PDonesidedMinus2nd1gt13 = PDonesidedMinus2nd1(&gt13[index]);
    CCTK_REAL const PDonesidedPlus2nd2gt13 = PDonesidedPlus2nd2(&gt13[index]);
    CCTK_REAL const PDonesidedMinus2nd2gt13 = PDonesidedMinus2nd2(&gt13[index]);
    CCTK_REAL const PDonesidedPlus2nd3gt13 = PDonesidedPlus2nd3(&gt13[index]);
    CCTK_REAL const PDonesidedMinus2nd3gt13 = PDonesidedMinus2nd3(&gt13[index]);
    CCTK_REAL const PDstandard2nd1gt22 = PDstandard2nd1(&gt22[index]);
    CCTK_REAL const PDstandard2nd2gt22 = PDstandard2nd2(&gt22[index]);
    CCTK_REAL const PDstandard2nd3gt22 = PDstandard2nd3(&gt22[index]);
    CCTK_REAL const PDonesidedPlus2nd1gt22 = PDonesidedPlus2nd1(&gt22[index]);
    CCTK_REAL const PDonesidedMinus2nd1gt22 = PDonesidedMinus2nd1(&gt22[index]);
    CCTK_REAL const PDonesidedPlus2nd2gt22 = PDonesidedPlus2nd2(&gt22[index]);
    CCTK_REAL const PDonesidedMinus2nd2gt22 = PDonesidedMinus2nd2(&gt22[index]);
    CCTK_REAL const PDonesidedPlus2nd3gt22 = PDonesidedPlus2nd3(&gt22[index]);
    CCTK_REAL const PDonesidedMinus2nd3gt22 = PDonesidedMinus2nd3(&gt22[index]);
    CCTK_REAL const PDstandard2nd1gt23 = PDstandard2nd1(&gt23[index]);
    CCTK_REAL const PDstandard2nd2gt23 = PDstandard2nd2(&gt23[index]);
    CCTK_REAL const PDstandard2nd3gt23 = PDstandard2nd3(&gt23[index]);
    CCTK_REAL const PDonesidedPlus2nd1gt23 = PDonesidedPlus2nd1(&gt23[index]);
    CCTK_REAL const PDonesidedMinus2nd1gt23 = PDonesidedMinus2nd1(&gt23[index]);
    CCTK_REAL const PDonesidedPlus2nd2gt23 = PDonesidedPlus2nd2(&gt23[index]);
    CCTK_REAL const PDonesidedMinus2nd2gt23 = PDonesidedMinus2nd2(&gt23[index]);
    CCTK_REAL const PDonesidedPlus2nd3gt23 = PDonesidedPlus2nd3(&gt23[index]);
    CCTK_REAL const PDonesidedMinus2nd3gt23 = PDonesidedMinus2nd3(&gt23[index]);
    CCTK_REAL const PDstandard2nd1gt33 = PDstandard2nd1(&gt33[index]);
    CCTK_REAL const PDstandard2nd2gt33 = PDstandard2nd2(&gt33[index]);
    CCTK_REAL const PDstandard2nd3gt33 = PDstandard2nd3(&gt33[index]);
    CCTK_REAL const PDonesidedPlus2nd1gt33 = PDonesidedPlus2nd1(&gt33[index]);
    CCTK_REAL const PDonesidedMinus2nd1gt33 = PDonesidedMinus2nd1(&gt33[index]);
    CCTK_REAL const PDonesidedPlus2nd2gt33 = PDonesidedPlus2nd2(&gt33[index]);
    CCTK_REAL const PDonesidedMinus2nd2gt33 = PDonesidedMinus2nd2(&gt33[index]);
    CCTK_REAL const PDonesidedPlus2nd3gt33 = PDonesidedPlus2nd3(&gt33[index]);
    CCTK_REAL const PDonesidedMinus2nd3gt33 = PDonesidedMinus2nd3(&gt33[index]);
    CCTK_REAL const PDstandard2nd1phi = PDstandard2nd1(&phi[index]);
    CCTK_REAL const PDstandard2nd2phi = PDstandard2nd2(&phi[index]);
    CCTK_REAL const PDstandard2nd3phi = PDstandard2nd3(&phi[index]);
    CCTK_REAL const PDonesidedPlus2nd1phi = PDonesidedPlus2nd1(&phi[index]);
    CCTK_REAL const PDonesidedMinus2nd1phi = PDonesidedMinus2nd1(&phi[index]);
    CCTK_REAL const PDonesidedPlus2nd2phi = PDonesidedPlus2nd2(&phi[index]);
    CCTK_REAL const PDonesidedMinus2nd2phi = PDonesidedMinus2nd2(&phi[index]);
    CCTK_REAL const PDonesidedPlus2nd3phi = PDonesidedPlus2nd3(&phi[index]);
    CCTK_REAL const PDonesidedMinus2nd3phi = PDonesidedMinus2nd3(&phi[index]);
    CCTK_REAL const PDstandard2nd1trK = PDstandard2nd1(&trK[index]);
    CCTK_REAL const PDstandard2nd2trK = PDstandard2nd2(&trK[index]);
    CCTK_REAL const PDstandard2nd3trK = PDstandard2nd3(&trK[index]);
    CCTK_REAL const PDonesidedPlus2nd1trK = PDonesidedPlus2nd1(&trK[index]);
    CCTK_REAL const PDonesidedMinus2nd1trK = PDonesidedMinus2nd1(&trK[index]);
    CCTK_REAL const PDonesidedPlus2nd2trK = PDonesidedPlus2nd2(&trK[index]);
    CCTK_REAL const PDonesidedMinus2nd2trK = PDonesidedMinus2nd2(&trK[index]);
    CCTK_REAL const PDonesidedPlus2nd3trK = PDonesidedPlus2nd3(&trK[index]);
    CCTK_REAL const PDonesidedMinus2nd3trK = PDonesidedMinus2nd3(&trK[index]);
    CCTK_REAL const PDstandard2nd1Xt1 = PDstandard2nd1(&Xt1[index]);
    CCTK_REAL const PDstandard2nd2Xt1 = PDstandard2nd2(&Xt1[index]);
    CCTK_REAL const PDstandard2nd3Xt1 = PDstandard2nd3(&Xt1[index]);
    CCTK_REAL const PDonesidedPlus2nd1Xt1 = PDonesidedPlus2nd1(&Xt1[index]);
    CCTK_REAL const PDonesidedMinus2nd1Xt1 = PDonesidedMinus2nd1(&Xt1[index]);
    CCTK_REAL const PDonesidedPlus2nd2Xt1 = PDonesidedPlus2nd2(&Xt1[index]);
    CCTK_REAL const PDonesidedMinus2nd2Xt1 = PDonesidedMinus2nd2(&Xt1[index]);
    CCTK_REAL const PDonesidedPlus2nd3Xt1 = PDonesidedPlus2nd3(&Xt1[index]);
    CCTK_REAL const PDonesidedMinus2nd3Xt1 = PDonesidedMinus2nd3(&Xt1[index]);
    CCTK_REAL const PDstandard2nd1Xt2 = PDstandard2nd1(&Xt2[index]);
    CCTK_REAL const PDstandard2nd2Xt2 = PDstandard2nd2(&Xt2[index]);
    CCTK_REAL const PDstandard2nd3Xt2 = PDstandard2nd3(&Xt2[index]);
    CCTK_REAL const PDonesidedPlus2nd1Xt2 = PDonesidedPlus2nd1(&Xt2[index]);
    CCTK_REAL const PDonesidedMinus2nd1Xt2 = PDonesidedMinus2nd1(&Xt2[index]);
    CCTK_REAL const PDonesidedPlus2nd2Xt2 = PDonesidedPlus2nd2(&Xt2[index]);
    CCTK_REAL const PDonesidedMinus2nd2Xt2 = PDonesidedMinus2nd2(&Xt2[index]);
    CCTK_REAL const PDonesidedPlus2nd3Xt2 = PDonesidedPlus2nd3(&Xt2[index]);
    CCTK_REAL const PDonesidedMinus2nd3Xt2 = PDonesidedMinus2nd3(&Xt2[index]);
    CCTK_REAL const PDstandard2nd1Xt3 = PDstandard2nd1(&Xt3[index]);
    CCTK_REAL const PDstandard2nd2Xt3 = PDstandard2nd2(&Xt3[index]);
    CCTK_REAL const PDstandard2nd3Xt3 = PDstandard2nd3(&Xt3[index]);
    CCTK_REAL const PDonesidedPlus2nd1Xt3 = PDonesidedPlus2nd1(&Xt3[index]);
    CCTK_REAL const PDonesidedMinus2nd1Xt3 = PDonesidedMinus2nd1(&Xt3[index]);
    CCTK_REAL const PDonesidedPlus2nd2Xt3 = PDonesidedPlus2nd2(&Xt3[index]);
    CCTK_REAL const PDonesidedMinus2nd2Xt3 = PDonesidedMinus2nd2(&Xt3[index]);
    CCTK_REAL const PDonesidedPlus2nd3Xt3 = PDonesidedPlus2nd3(&Xt3[index]);
    CCTK_REAL const PDonesidedMinus2nd3Xt3 = PDonesidedMinus2nd3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL rn1 = -(xCopyL*INV(rCopyL));
    
    CCTK_REAL rn2 = -(yCopyL*INV(rCopyL));
    
    CCTK_REAL rn3 = -(zCopyL*INV(rCopyL));
    
    CCTK_REAL phi0 = IfThen(conformalMethod,1,0);
    
    CCTK_REAL v0 = sqrt(ToReal(harmonicF));
    
    CCTK_REAL phirhsL = v0*(-phiL + phi0 + rCopyL*(rn1*IfThen(normal[0] 
      < 0,PDonesidedPlus2nd1phi,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1phi,PDstandard2nd1phi)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2phi,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2phi,PDstandard2nd2phi)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3phi,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3phi,PDstandard2nd3phi))))*INV(rCopyL);
    
    CCTK_REAL gt11rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt11,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt11,PDstandard2nd1gt11)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt11,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt11,PDstandard2nd2gt11)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt11,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt11,PDstandard2nd3gt11)) + (1 - 
      gt11L)*INV(rCopyL);
    
    CCTK_REAL gt12rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt12,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt12,PDstandard2nd1gt12)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt12,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt12,PDstandard2nd2gt12)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt12,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt12,PDstandard2nd3gt12)) - gt12L*INV(rCopyL);
    
    CCTK_REAL gt13rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt13,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt13,PDstandard2nd1gt13)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt13,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt13,PDstandard2nd2gt13)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt13,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt13,PDstandard2nd3gt13)) - gt13L*INV(rCopyL);
    
    CCTK_REAL gt22rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt22,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt22,PDstandard2nd1gt22)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt22,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt22,PDstandard2nd2gt22)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt22,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt22,PDstandard2nd3gt22)) + (1 - 
      gt22L)*INV(rCopyL);
    
    CCTK_REAL gt23rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt23,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt23,PDstandard2nd1gt23)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt23,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt23,PDstandard2nd2gt23)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt23,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt23,PDstandard2nd3gt23)) - gt23L*INV(rCopyL);
    
    CCTK_REAL gt33rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1gt33,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1gt33,PDstandard2nd1gt33)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2gt33,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2gt33,PDstandard2nd2gt33)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3gt33,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3gt33,PDstandard2nd3gt33)) + (1 - 
      gt33L)*INV(rCopyL);
    
    CCTK_REAL trKrhsL = v0*(rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1trK,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1trK,PDstandard2nd1trK)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2trK,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2trK,PDstandard2nd2trK)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3trK,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3trK,PDstandard2nd3trK)) - trKL*INV(rCopyL));
    
    CCTK_REAL At11rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At11,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At11,PDstandard2nd1At11)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At11,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At11,PDstandard2nd2At11)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At11,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At11,PDstandard2nd3At11)) - At11L*INV(rCopyL);
    
    CCTK_REAL At12rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At12,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At12,PDstandard2nd1At12)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At12,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At12,PDstandard2nd2At12)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At12,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At12,PDstandard2nd3At12)) - At12L*INV(rCopyL);
    
    CCTK_REAL At13rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At13,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At13,PDstandard2nd1At13)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At13,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At13,PDstandard2nd2At13)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At13,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At13,PDstandard2nd3At13)) - At13L*INV(rCopyL);
    
    CCTK_REAL At22rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At22,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At22,PDstandard2nd1At22)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At22,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At22,PDstandard2nd2At22)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At22,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At22,PDstandard2nd3At22)) - At22L*INV(rCopyL);
    
    CCTK_REAL At23rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At23,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At23,PDstandard2nd1At23)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At23,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At23,PDstandard2nd2At23)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At23,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At23,PDstandard2nd3At23)) - At23L*INV(rCopyL);
    
    CCTK_REAL At33rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1At33,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1At33,PDstandard2nd1At33)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2At33,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2At33,PDstandard2nd2At33)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3At33,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3At33,PDstandard2nd3At33)) - At33L*INV(rCopyL);
    
    CCTK_REAL Xt1rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1Xt1,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt1,PDstandard2nd1Xt1)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2Xt1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt1,PDstandard2nd2Xt1)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3Xt1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt1,PDstandard2nd3Xt1)) - Xt1L*INV(rCopyL);
    
    CCTK_REAL Xt2rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1Xt2,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt2,PDstandard2nd1Xt2)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2Xt2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt2,PDstandard2nd2Xt2)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3Xt2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt2,PDstandard2nd3Xt2)) - Xt2L*INV(rCopyL);
    
    CCTK_REAL Xt3rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1Xt3,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1Xt3,PDstandard2nd1Xt3)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2Xt3,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2Xt3,PDstandard2nd2Xt3)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3Xt3,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3Xt3,PDstandard2nd3Xt3)) - Xt3L*INV(rCopyL);
    
    CCTK_REAL alpharhsL = v0*(rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1alpha,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1alpha,PDstandard2nd1alpha)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2alpha,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2alpha,PDstandard2nd2alpha)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3alpha,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3alpha,PDstandard2nd3alpha)) + (1 - 
      alphaL)*INV(rCopyL));
    
    CCTK_REAL ArhsL = v0*(rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1A,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1A,PDstandard2nd1A)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2A,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2A,PDstandard2nd2A)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3A,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3A,PDstandard2nd3A)) - AL*INV(rCopyL));
    
    CCTK_REAL beta1rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1beta1,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1beta1,PDstandard2nd1beta1)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2beta1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta1,PDstandard2nd2beta1)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3beta1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta1,PDstandard2nd3beta1)) - 
      beta1L*INV(rCopyL);
    
    CCTK_REAL beta2rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1beta2,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1beta2,PDstandard2nd1beta2)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2beta2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta2,PDstandard2nd2beta2)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3beta2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta2,PDstandard2nd3beta2)) - 
      beta2L*INV(rCopyL);
    
    CCTK_REAL beta3rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1beta3,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1beta3,PDstandard2nd1beta3)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2beta3,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2beta3,PDstandard2nd2beta3)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3beta3,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3beta3,PDstandard2nd3beta3)) - 
      beta3L*INV(rCopyL);
    
    CCTK_REAL B1rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1B1,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1B1,PDstandard2nd1B1)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2B1,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2B1,PDstandard2nd2B1)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3B1,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3B1,PDstandard2nd3B1)) - B1L*INV(rCopyL);
    
    CCTK_REAL B2rhsL = rn1*IfThen(normal[0] < 
      0,PDonesidedPlus2nd1B2,IfThen(normal[0] > 
      0,PDonesidedMinus2nd1B2,PDstandard2nd1B2)) + rn2*IfThen(normal[1] < 
      0,PDonesidedPlus2nd2B2,IfThen(normal[1] > 
      0,PDonesidedMinus2nd2B2,PDstandard2nd2B2)) + rn3*IfThen(normal[2] < 
      0,PDonesidedPlus2nd3B2,IfThen(normal[2] > 
      0,PDonesidedMinus2nd3B2,PDstandard2nd3B2)) - B2L*INV(rCopyL);
    
    CCTK_REAL B3rhsL = rn1*IfThen(normal[0] < 
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
  
  const char *const groups[] = {
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
