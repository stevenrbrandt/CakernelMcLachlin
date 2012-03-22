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
#define QAD(x) (SQR(SQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))

extern "C" void ML_BSSN_Host_Advect_SelectBCs(CCTK_ARGUMENTS)
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

static void ML_BSSN_Host_Advect_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_REAL_VEC const p1o1024dx = kmul(INV(dx),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1024dy = kmul(INV(dy),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1024dz = kmul(INV(dz),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1680dx = kmul(INV(dx),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o1680dy = kmul(INV(dy),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o1680dz = kmul(INV(dz),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o5040dx2 = kmul(INV(SQR(dx)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dy2 = kmul(INV(SQR(dy)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dz2 = kmul(INV(SQR(dz)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o560dx = kmul(INV(dx),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o560dy = kmul(INV(dy),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o560dz = kmul(INV(dz),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o705600dxdy = kmul(INV(kmul(dx,dy)),ToReal(1.41723356009070294784580498866e-6));
  CCTK_REAL_VEC const p1o705600dxdz = kmul(INV(kmul(dx,dz)),ToReal(1.41723356009070294784580498866e-6));
  CCTK_REAL_VEC const p1o705600dydz = kmul(INV(kmul(dy,dz)),ToReal(1.41723356009070294784580498866e-6));
  CCTK_REAL_VEC const p1o840dx = kmul(INV(dx),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dy = kmul(INV(dy),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dz = kmul(INV(dz),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1odx = INV(dx);
  CCTK_REAL_VEC const p1ody = INV(dy);
  CCTK_REAL_VEC const p1odz = INV(dz);
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC(ML_BSSN_Host_Advect,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC AL = vec_load(A[index]);
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC alpharhsL = vec_load(alpharhs[index]);
    CCTK_REAL_VEC ArhsL = vec_load(Arhs[index]);
    CCTK_REAL_VEC At11L = vec_load(At11[index]);
    CCTK_REAL_VEC At11rhsL = vec_load(At11rhs[index]);
    CCTK_REAL_VEC At12L = vec_load(At12[index]);
    CCTK_REAL_VEC At12rhsL = vec_load(At12rhs[index]);
    CCTK_REAL_VEC At13L = vec_load(At13[index]);
    CCTK_REAL_VEC At13rhsL = vec_load(At13rhs[index]);
    CCTK_REAL_VEC At22L = vec_load(At22[index]);
    CCTK_REAL_VEC At22rhsL = vec_load(At22rhs[index]);
    CCTK_REAL_VEC At23L = vec_load(At23[index]);
    CCTK_REAL_VEC At23rhsL = vec_load(At23rhs[index]);
    CCTK_REAL_VEC At33L = vec_load(At33[index]);
    CCTK_REAL_VEC At33rhsL = vec_load(At33rhs[index]);
    CCTK_REAL_VEC B1L = vec_load(B1[index]);
    CCTK_REAL_VEC B1rhsL = vec_load(B1rhs[index]);
    CCTK_REAL_VEC B2L = vec_load(B2[index]);
    CCTK_REAL_VEC B2rhsL = vec_load(B2rhs[index]);
    CCTK_REAL_VEC B3L = vec_load(B3[index]);
    CCTK_REAL_VEC B3rhsL = vec_load(B3rhs[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta1rhsL = vec_load(beta1rhs[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta2rhsL = vec_load(beta2rhs[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC beta3rhsL = vec_load(beta3rhs[index]);
    CCTK_REAL_VEC gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC gt11rhsL = vec_load(gt11rhs[index]);
    CCTK_REAL_VEC gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC gt12rhsL = vec_load(gt12rhs[index]);
    CCTK_REAL_VEC gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC gt13rhsL = vec_load(gt13rhs[index]);
    CCTK_REAL_VEC gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC gt22rhsL = vec_load(gt22rhs[index]);
    CCTK_REAL_VEC gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC gt23rhsL = vec_load(gt23rhs[index]);
    CCTK_REAL_VEC gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC gt33rhsL = vec_load(gt33rhs[index]);
    CCTK_REAL_VEC phiL = vec_load(phi[index]);
    CCTK_REAL_VEC phirhsL = vec_load(phirhs[index]);
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    CCTK_REAL_VEC trKrhsL = vec_load(trKrhs[index]);
    CCTK_REAL_VEC Xt1L = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt1rhsL = vec_load(Xt1rhs[index]);
    CCTK_REAL_VEC Xt2L = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt2rhsL = vec_load(Xt2rhs[index]);
    CCTK_REAL_VEC Xt3L = vec_load(Xt3[index]);
    CCTK_REAL_VEC Xt3rhsL = vec_load(Xt3rhs[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC const PDupwindNthAnti1A = PDupwindNthAnti1(&A[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1A = PDupwindNthSymm1(&A[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2A = PDupwindNthAnti2(&A[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2A = PDupwindNthSymm2(&A[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3A = PDupwindNthAnti3(&A[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3A = PDupwindNthSymm3(&A[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1alpha = PDupwindNthAnti1(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1alpha = PDupwindNthSymm1(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2alpha = PDupwindNthAnti2(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2alpha = PDupwindNthSymm2(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3alpha = PDupwindNthAnti3(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3alpha = PDupwindNthSymm3(&alpha[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1At11 = PDupwindNthAnti1(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At11 = PDupwindNthSymm1(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At11 = PDupwindNthAnti2(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At11 = PDupwindNthSymm2(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At11 = PDupwindNthAnti3(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At11 = PDupwindNthSymm3(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1At12 = PDupwindNthAnti1(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At12 = PDupwindNthSymm1(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At12 = PDupwindNthAnti2(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At12 = PDupwindNthSymm2(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At12 = PDupwindNthAnti3(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At12 = PDupwindNthSymm3(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1At13 = PDupwindNthAnti1(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At13 = PDupwindNthSymm1(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At13 = PDupwindNthAnti2(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At13 = PDupwindNthSymm2(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At13 = PDupwindNthAnti3(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At13 = PDupwindNthSymm3(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1At22 = PDupwindNthAnti1(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At22 = PDupwindNthSymm1(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At22 = PDupwindNthAnti2(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At22 = PDupwindNthSymm2(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At22 = PDupwindNthAnti3(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At22 = PDupwindNthSymm3(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1At23 = PDupwindNthAnti1(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At23 = PDupwindNthSymm1(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At23 = PDupwindNthAnti2(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At23 = PDupwindNthSymm2(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At23 = PDupwindNthAnti3(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At23 = PDupwindNthSymm3(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1At33 = PDupwindNthAnti1(&At33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At33 = PDupwindNthSymm1(&At33[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At33 = PDupwindNthAnti2(&At33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At33 = PDupwindNthSymm2(&At33[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At33 = PDupwindNthAnti3(&At33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At33 = PDupwindNthSymm3(&At33[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1B1 = PDupwindNthAnti1(&B1[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1B1 = PDupwindNthSymm1(&B1[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2B1 = PDupwindNthAnti2(&B1[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2B1 = PDupwindNthSymm2(&B1[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3B1 = PDupwindNthAnti3(&B1[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3B1 = PDupwindNthSymm3(&B1[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1B2 = PDupwindNthAnti1(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1B2 = PDupwindNthSymm1(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2B2 = PDupwindNthAnti2(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2B2 = PDupwindNthSymm2(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3B2 = PDupwindNthAnti3(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3B2 = PDupwindNthSymm3(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1B3 = PDupwindNthAnti1(&B3[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1B3 = PDupwindNthSymm1(&B3[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2B3 = PDupwindNthAnti2(&B3[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2B3 = PDupwindNthSymm2(&B3[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3B3 = PDupwindNthAnti3(&B3[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3B3 = PDupwindNthSymm3(&B3[index]);
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
    CCTK_REAL_VEC const PDupwindNthAnti1gt11 = PDupwindNthAnti1(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt11 = PDupwindNthSymm1(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt11 = PDupwindNthAnti2(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt11 = PDupwindNthSymm2(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt11 = PDupwindNthAnti3(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt11 = PDupwindNthSymm3(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1gt12 = PDupwindNthAnti1(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt12 = PDupwindNthSymm1(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt12 = PDupwindNthAnti2(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt12 = PDupwindNthSymm2(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt12 = PDupwindNthAnti3(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt12 = PDupwindNthSymm3(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1gt13 = PDupwindNthAnti1(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt13 = PDupwindNthSymm1(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt13 = PDupwindNthAnti2(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt13 = PDupwindNthSymm2(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt13 = PDupwindNthAnti3(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt13 = PDupwindNthSymm3(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1gt22 = PDupwindNthAnti1(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt22 = PDupwindNthSymm1(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt22 = PDupwindNthAnti2(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt22 = PDupwindNthSymm2(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt22 = PDupwindNthAnti3(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt22 = PDupwindNthSymm3(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1gt23 = PDupwindNthAnti1(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt23 = PDupwindNthSymm1(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt23 = PDupwindNthAnti2(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt23 = PDupwindNthSymm2(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt23 = PDupwindNthAnti3(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt23 = PDupwindNthSymm3(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1gt33 = PDupwindNthAnti1(&gt33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt33 = PDupwindNthSymm1(&gt33[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt33 = PDupwindNthAnti2(&gt33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt33 = PDupwindNthSymm2(&gt33[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt33 = PDupwindNthAnti3(&gt33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt33 = PDupwindNthSymm3(&gt33[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1phi = PDupwindNthAnti1(&phi[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1phi = PDupwindNthSymm1(&phi[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2phi = PDupwindNthAnti2(&phi[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2phi = PDupwindNthSymm2(&phi[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3phi = PDupwindNthAnti3(&phi[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3phi = PDupwindNthSymm3(&phi[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1trK = PDupwindNthAnti1(&trK[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1trK = PDupwindNthSymm1(&trK[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2trK = PDupwindNthAnti2(&trK[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2trK = PDupwindNthSymm2(&trK[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3trK = PDupwindNthAnti3(&trK[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3trK = PDupwindNthSymm3(&trK[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1Xt1 = PDupwindNthAnti1(&Xt1[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1Xt1 = PDupwindNthSymm1(&Xt1[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2Xt1 = PDupwindNthAnti2(&Xt1[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2Xt1 = PDupwindNthSymm2(&Xt1[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3Xt1 = PDupwindNthAnti3(&Xt1[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3Xt1 = PDupwindNthSymm3(&Xt1[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1Xt2 = PDupwindNthAnti1(&Xt2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1Xt2 = PDupwindNthSymm1(&Xt2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2Xt2 = PDupwindNthAnti2(&Xt2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2Xt2 = PDupwindNthSymm2(&Xt2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3Xt2 = PDupwindNthAnti3(&Xt2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3Xt2 = PDupwindNthSymm3(&Xt2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1Xt3 = PDupwindNthAnti1(&Xt3[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1Xt3 = PDupwindNthSymm1(&Xt3[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2Xt3 = PDupwindNthAnti2(&Xt3[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2Xt3 = PDupwindNthSymm2(&Xt3[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3Xt3 = PDupwindNthAnti3(&Xt3[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3Xt3 = PDupwindNthSymm3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    phirhsL = 
      kadd(phirhsL,kmadd(beta1L,PDupwindNthAnti1phi,kmadd(beta2L,PDupwindNthAnti2phi,kmadd(beta3L,PDupwindNthAnti3phi,kmadd(PDupwindNthSymm1phi,kfabs(beta1L),kmadd(PDupwindNthSymm2phi,kfabs(beta2L),kmul(PDupwindNthSymm3phi,kfabs(beta3L))))))));
    
    gt11rhsL = 
      kadd(gt11rhsL,kmadd(beta1L,PDupwindNthAnti1gt11,kmadd(beta2L,PDupwindNthAnti2gt11,kmadd(beta3L,PDupwindNthAnti3gt11,kmadd(PDupwindNthSymm1gt11,kfabs(beta1L),kmadd(PDupwindNthSymm2gt11,kfabs(beta2L),kmul(PDupwindNthSymm3gt11,kfabs(beta3L))))))));
    
    gt12rhsL = 
      kadd(gt12rhsL,kmadd(beta1L,PDupwindNthAnti1gt12,kmadd(beta2L,PDupwindNthAnti2gt12,kmadd(beta3L,PDupwindNthAnti3gt12,kmadd(PDupwindNthSymm1gt12,kfabs(beta1L),kmadd(PDupwindNthSymm2gt12,kfabs(beta2L),kmul(PDupwindNthSymm3gt12,kfabs(beta3L))))))));
    
    gt13rhsL = 
      kadd(gt13rhsL,kmadd(beta1L,PDupwindNthAnti1gt13,kmadd(beta2L,PDupwindNthAnti2gt13,kmadd(beta3L,PDupwindNthAnti3gt13,kmadd(PDupwindNthSymm1gt13,kfabs(beta1L),kmadd(PDupwindNthSymm2gt13,kfabs(beta2L),kmul(PDupwindNthSymm3gt13,kfabs(beta3L))))))));
    
    gt22rhsL = 
      kadd(gt22rhsL,kmadd(beta1L,PDupwindNthAnti1gt22,kmadd(beta2L,PDupwindNthAnti2gt22,kmadd(beta3L,PDupwindNthAnti3gt22,kmadd(PDupwindNthSymm1gt22,kfabs(beta1L),kmadd(PDupwindNthSymm2gt22,kfabs(beta2L),kmul(PDupwindNthSymm3gt22,kfabs(beta3L))))))));
    
    gt23rhsL = 
      kadd(gt23rhsL,kmadd(beta1L,PDupwindNthAnti1gt23,kmadd(beta2L,PDupwindNthAnti2gt23,kmadd(beta3L,PDupwindNthAnti3gt23,kmadd(PDupwindNthSymm1gt23,kfabs(beta1L),kmadd(PDupwindNthSymm2gt23,kfabs(beta2L),kmul(PDupwindNthSymm3gt23,kfabs(beta3L))))))));
    
    gt33rhsL = 
      kadd(gt33rhsL,kmadd(beta1L,PDupwindNthAnti1gt33,kmadd(beta2L,PDupwindNthAnti2gt33,kmadd(beta3L,PDupwindNthAnti3gt33,kmadd(PDupwindNthSymm1gt33,kfabs(beta1L),kmadd(PDupwindNthSymm2gt33,kfabs(beta2L),kmul(PDupwindNthSymm3gt33,kfabs(beta3L))))))));
    
    Xt1rhsL = 
      kadd(Xt1rhsL,kmadd(beta1L,PDupwindNthAnti1Xt1,kmadd(beta2L,PDupwindNthAnti2Xt1,kmadd(beta3L,PDupwindNthAnti3Xt1,kmadd(PDupwindNthSymm1Xt1,kfabs(beta1L),kmadd(PDupwindNthSymm2Xt1,kfabs(beta2L),kmul(PDupwindNthSymm3Xt1,kfabs(beta3L))))))));
    
    Xt2rhsL = 
      kadd(Xt2rhsL,kmadd(beta1L,PDupwindNthAnti1Xt2,kmadd(beta2L,PDupwindNthAnti2Xt2,kmadd(beta3L,PDupwindNthAnti3Xt2,kmadd(PDupwindNthSymm1Xt2,kfabs(beta1L),kmadd(PDupwindNthSymm2Xt2,kfabs(beta2L),kmul(PDupwindNthSymm3Xt2,kfabs(beta3L))))))));
    
    Xt3rhsL = 
      kadd(Xt3rhsL,kmadd(beta1L,PDupwindNthAnti1Xt3,kmadd(beta2L,PDupwindNthAnti2Xt3,kmadd(beta3L,PDupwindNthAnti3Xt3,kmadd(PDupwindNthSymm1Xt3,kfabs(beta1L),kmadd(PDupwindNthSymm2Xt3,kfabs(beta2L),kmul(PDupwindNthSymm3Xt3,kfabs(beta3L))))))));
    
    trKrhsL = 
      kadd(trKrhsL,kmadd(beta1L,PDupwindNthAnti1trK,kmadd(beta2L,PDupwindNthAnti2trK,kmadd(beta3L,PDupwindNthAnti3trK,kmadd(PDupwindNthSymm1trK,kfabs(beta1L),kmadd(PDupwindNthSymm2trK,kfabs(beta2L),kmul(PDupwindNthSymm3trK,kfabs(beta3L))))))));
    
    At11rhsL = 
      kadd(At11rhsL,kmadd(beta1L,PDupwindNthAnti1At11,kmadd(beta2L,PDupwindNthAnti2At11,kmadd(beta3L,PDupwindNthAnti3At11,kmadd(PDupwindNthSymm1At11,kfabs(beta1L),kmadd(PDupwindNthSymm2At11,kfabs(beta2L),kmul(PDupwindNthSymm3At11,kfabs(beta3L))))))));
    
    At12rhsL = 
      kadd(At12rhsL,kmadd(beta1L,PDupwindNthAnti1At12,kmadd(beta2L,PDupwindNthAnti2At12,kmadd(beta3L,PDupwindNthAnti3At12,kmadd(PDupwindNthSymm1At12,kfabs(beta1L),kmadd(PDupwindNthSymm2At12,kfabs(beta2L),kmul(PDupwindNthSymm3At12,kfabs(beta3L))))))));
    
    At13rhsL = 
      kadd(At13rhsL,kmadd(beta1L,PDupwindNthAnti1At13,kmadd(beta2L,PDupwindNthAnti2At13,kmadd(beta3L,PDupwindNthAnti3At13,kmadd(PDupwindNthSymm1At13,kfabs(beta1L),kmadd(PDupwindNthSymm2At13,kfabs(beta2L),kmul(PDupwindNthSymm3At13,kfabs(beta3L))))))));
    
    At22rhsL = 
      kadd(At22rhsL,kmadd(beta1L,PDupwindNthAnti1At22,kmadd(beta2L,PDupwindNthAnti2At22,kmadd(beta3L,PDupwindNthAnti3At22,kmadd(PDupwindNthSymm1At22,kfabs(beta1L),kmadd(PDupwindNthSymm2At22,kfabs(beta2L),kmul(PDupwindNthSymm3At22,kfabs(beta3L))))))));
    
    At23rhsL = 
      kadd(At23rhsL,kmadd(beta1L,PDupwindNthAnti1At23,kmadd(beta2L,PDupwindNthAnti2At23,kmadd(beta3L,PDupwindNthAnti3At23,kmadd(PDupwindNthSymm1At23,kfabs(beta1L),kmadd(PDupwindNthSymm2At23,kfabs(beta2L),kmul(PDupwindNthSymm3At23,kfabs(beta3L))))))));
    
    At33rhsL = 
      kadd(At33rhsL,kmadd(beta1L,PDupwindNthAnti1At33,kmadd(beta2L,PDupwindNthAnti2At33,kmadd(beta3L,PDupwindNthAnti3At33,kmadd(PDupwindNthSymm1At33,kfabs(beta1L),kmadd(PDupwindNthSymm2At33,kfabs(beta2L),kmul(PDupwindNthSymm3At33,kfabs(beta3L))))))));
    
    alpharhsL = 
      kmadd(kmadd(beta1L,PDupwindNthAnti1alpha,kmadd(beta2L,PDupwindNthAnti2alpha,kmadd(beta3L,PDupwindNthAnti3alpha,kmadd(PDupwindNthSymm1alpha,kfabs(beta1L),kmadd(PDupwindNthSymm2alpha,kfabs(beta2L),kmul(PDupwindNthSymm3alpha,kfabs(beta3L))))))),ToReal(LapseAdvectionCoeff),alpharhsL);
    
    ArhsL = 
      kmadd(ToReal(LapseACoeff),kmsub(kmadd(beta1L,PDupwindNthAnti1A,kmadd(beta2L,PDupwindNthAnti2A,kmadd(beta3L,PDupwindNthAnti3A,kmadd(PDupwindNthSymm1A,kfabs(beta1L),kmadd(PDupwindNthSymm2A,kfabs(beta2L),kmul(PDupwindNthSymm3A,kfabs(beta3L))))))),ToReal(LapseAdvectionCoeff),kmul(kmadd(beta1L,PDupwindNthAnti1trK,kmadd(beta2L,PDupwindNthAnti2trK,kmadd(beta3L,PDupwindNthAnti3trK,kmadd(PDupwindNthSymm1trK,kfabs(beta1L),kmadd(PDupwindNthSymm2trK,kfabs(beta2L),kmul(PDupwindNthSymm3trK,kfabs(beta3L))))))),kadd(ToReal(-1),ToReal(LapseAdvectionCoeff)))),ArhsL);
    
    beta1rhsL = 
      kmadd(kmadd(beta1L,PDupwindNthAnti1beta1,kmadd(beta2L,PDupwindNthAnti2beta1,kmadd(beta3L,PDupwindNthAnti3beta1,kmadd(PDupwindNthSymm1beta1,kfabs(beta1L),kmadd(PDupwindNthSymm2beta1,kfabs(beta2L),kmul(PDupwindNthSymm3beta1,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),beta1rhsL);
    
    beta2rhsL = 
      kmadd(kmadd(beta1L,PDupwindNthAnti1beta2,kmadd(beta2L,PDupwindNthAnti2beta2,kmadd(beta3L,PDupwindNthAnti3beta2,kmadd(PDupwindNthSymm1beta2,kfabs(beta1L),kmadd(PDupwindNthSymm2beta2,kfabs(beta2L),kmul(PDupwindNthSymm3beta2,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),beta2rhsL);
    
    beta3rhsL = 
      kmadd(kmadd(beta1L,PDupwindNthAnti1beta3,kmadd(beta2L,PDupwindNthAnti2beta3,kmadd(beta3L,PDupwindNthAnti3beta3,kmadd(PDupwindNthSymm1beta3,kfabs(beta1L),kmadd(PDupwindNthSymm2beta3,kfabs(beta2L),kmul(PDupwindNthSymm3beta3,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),beta3rhsL);
    
    B1rhsL = 
      kmadd(kmsub(kmadd(beta1L,PDupwindNthAnti1B1,kmadd(beta2L,PDupwindNthAnti2B1,kmadd(beta3L,PDupwindNthAnti3B1,kmadd(PDupwindNthSymm1B1,kfabs(beta1L),kmadd(PDupwindNthSymm2B1,kfabs(beta2L),kmul(PDupwindNthSymm3B1,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),kmul(kmadd(beta1L,PDupwindNthAnti1Xt1,kmadd(beta2L,PDupwindNthAnti2Xt1,kmadd(beta3L,PDupwindNthAnti3Xt1,kmadd(PDupwindNthSymm1Xt1,kfabs(beta1L),kmadd(PDupwindNthSymm2Xt1,kfabs(beta2L),kmul(PDupwindNthSymm3Xt1,kfabs(beta3L))))))),kadd(ToReal(-1),ToReal(ShiftAdvectionCoeff)))),ToReal(ShiftBCoeff),B1rhsL);
    
    B2rhsL = 
      kmadd(kmsub(kmadd(beta1L,PDupwindNthAnti1B2,kmadd(beta2L,PDupwindNthAnti2B2,kmadd(beta3L,PDupwindNthAnti3B2,kmadd(PDupwindNthSymm1B2,kfabs(beta1L),kmadd(PDupwindNthSymm2B2,kfabs(beta2L),kmul(PDupwindNthSymm3B2,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),kmul(kmadd(beta1L,PDupwindNthAnti1Xt2,kmadd(beta2L,PDupwindNthAnti2Xt2,kmadd(beta3L,PDupwindNthAnti3Xt2,kmadd(PDupwindNthSymm1Xt2,kfabs(beta1L),kmadd(PDupwindNthSymm2Xt2,kfabs(beta2L),kmul(PDupwindNthSymm3Xt2,kfabs(beta3L))))))),kadd(ToReal(-1),ToReal(ShiftAdvectionCoeff)))),ToReal(ShiftBCoeff),B2rhsL);
    
    B3rhsL = 
      kmadd(kmsub(kmadd(beta1L,PDupwindNthAnti1B3,kmadd(beta2L,PDupwindNthAnti2B3,kmadd(beta3L,PDupwindNthAnti3B3,kmadd(PDupwindNthSymm1B3,kfabs(beta1L),kmadd(PDupwindNthSymm2B3,kfabs(beta2L),kmul(PDupwindNthSymm3B3,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),kmul(kmadd(beta1L,PDupwindNthAnti1Xt3,kmadd(beta2L,PDupwindNthAnti2Xt3,kmadd(beta3L,PDupwindNthAnti3Xt3,kmadd(PDupwindNthSymm1Xt3,kfabs(beta1L),kmadd(PDupwindNthSymm2Xt3,kfabs(beta2L),kmul(PDupwindNthSymm3Xt3,kfabs(beta3L))))))),kadd(ToReal(-1),ToReal(ShiftAdvectionCoeff)))),ToReal(ShiftBCoeff),B3rhsL);
    
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
  LC_ENDLOOP3VEC(ML_BSSN_Host_Advect);
}

extern "C" void ML_BSSN_Host_Advect(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_Host_Advect_Body");
  }
  
  if (cctk_iteration % ML_BSSN_Host_Advect_calc_every != ML_BSSN_Host_Advect_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
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
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_Advect", 18, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Host_Advect", 5, 5, 5);
  
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_Host_Advect_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_Advect_Body");
  }
}
