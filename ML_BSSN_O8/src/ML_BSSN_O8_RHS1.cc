/*  File produced by Kranc */

#define KRANC_C

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"
#include "Differencing.h"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define QAD(x) (SQR(SQR(x)))
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))

extern "C" void ML_BSSN_O8_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_O8::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_O8::ML_trace_curvrhs.");
  return;
}

static void ML_BSSN_O8_RHS1_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_O8_RHS1_Body");
  }
  
  if (cctk_iteration % ML_BSSN_O8_RHS1_calc_every != ML_BSSN_O8_RHS1_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"grid::coordinates","Grid::coordinates","ML_BSSN_O8::ML_curv","ML_BSSN_O8::ML_dtlapse","ML_BSSN_O8::ML_dtlapserhs","ML_BSSN_O8::ML_dtshift","ML_BSSN_O8::ML_dtshiftrhs","ML_BSSN_O8::ML_Gamma","ML_BSSN_O8::ML_Gammarhs","ML_BSSN_O8::ML_lapse","ML_BSSN_O8::ML_lapserhs","ML_BSSN_O8::ML_log_confac","ML_BSSN_O8::ML_log_confacrhs","ML_BSSN_O8::ML_metric","ML_BSSN_O8::ML_metricrhs","ML_BSSN_O8::ML_shift","ML_BSSN_O8::ML_shiftrhs","ML_BSSN_O8::ML_trace_curv","ML_BSSN_O8::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_O8_RHS1", 19, groups);
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  ptrdiff_t const di = 1;
  ptrdiff_t const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  CCTK_REAL const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL const dz = ToReal(CCTK_DELTA_SPACE(2));
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
  CCTK_REAL const p1o5040dx2 = 0.000198412698412698412698412698413*INV(SQR(dx));
  CCTK_REAL const p1o5040dy2 = 0.000198412698412698412698412698413*INV(SQR(dy));
  CCTK_REAL const p1o5040dz2 = 0.000198412698412698412698412698413*INV(SQR(dz));
  CCTK_REAL const p1o560dx = 0.00178571428571428571428571428571*INV(dx);
  CCTK_REAL const p1o560dy = 0.00178571428571428571428571428571*INV(dy);
  CCTK_REAL const p1o560dz = 0.00178571428571428571428571428571*INV(dz);
  CCTK_REAL const p1o705600dxdy = 1.41723356009070294784580498866e-6*INV(dx)*INV(dy);
  CCTK_REAL const p1o705600dxdz = 1.41723356009070294784580498866e-6*INV(dx)*INV(dz);
  CCTK_REAL const p1o705600dydz = 1.41723356009070294784580498866e-6*INV(dy)*INV(dz);
  CCTK_REAL const p1o840dx = 0.00119047619047619047619047619048*INV(dx);
  CCTK_REAL const p1o840dy = 0.00119047619047619047619047619048*INV(dy);
  CCTK_REAL const p1o840dz = 0.00119047619047619047619047619048*INV(dz);
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o840dx = -0.00119047619047619047619047619048*INV(dx);
  CCTK_REAL const pm1o840dy = -0.00119047619047619047619047619048*INV(dy);
  CCTK_REAL const pm1o840dz = -0.00119047619047619047619047619048*INV(dz);
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_O8_RHS1,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
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
    CCTK_REAL eTttL = (*stress_energy_state) ? eTtt[index] : ToReal(0.0);
    CCTK_REAL eTtxL = (*stress_energy_state) ? eTtx[index] : ToReal(0.0);
    CCTK_REAL eTtyL = (*stress_energy_state) ? eTty[index] : ToReal(0.0);
    CCTK_REAL eTtzL = (*stress_energy_state) ? eTtz[index] : ToReal(0.0);
    CCTK_REAL eTxxL = (*stress_energy_state) ? eTxx[index] : ToReal(0.0);
    CCTK_REAL eTxyL = (*stress_energy_state) ? eTxy[index] : ToReal(0.0);
    CCTK_REAL eTxzL = (*stress_energy_state) ? eTxz[index] : ToReal(0.0);
    CCTK_REAL eTyyL = (*stress_energy_state) ? eTyy[index] : ToReal(0.0);
    CCTK_REAL eTyzL = (*stress_energy_state) ? eTyz[index] : ToReal(0.0);
    CCTK_REAL eTzzL = (*stress_energy_state) ? eTzz[index] : ToReal(0.0);
    CCTK_REAL gt11L = gt11[index];
    CCTK_REAL gt12L = gt12[index];
    CCTK_REAL gt13L = gt13[index];
    CCTK_REAL gt22L = gt22[index];
    CCTK_REAL gt23L = gt23[index];
    CCTK_REAL gt33L = gt33[index];
    CCTK_REAL phiL = phi[index];
    CCTK_REAL rL = r[index];
    CCTK_REAL trKL = trK[index];
    CCTK_REAL Xt1L = Xt1[index];
    CCTK_REAL Xt2L = Xt2[index];
    CCTK_REAL Xt3L = Xt3[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDdissipationNth1A = PDdissipationNth1(&A[index]);
    CCTK_REAL const PDdissipationNth2A = PDdissipationNth2(&A[index]);
    CCTK_REAL const PDdissipationNth3A = PDdissipationNth3(&A[index]);
    CCTK_REAL const PDupwindNthAnti1A = PDupwindNthAnti1(&A[index]);
    CCTK_REAL const PDupwindNthSymm1A = PDupwindNthSymm1(&A[index]);
    CCTK_REAL const PDupwindNthAnti2A = PDupwindNthAnti2(&A[index]);
    CCTK_REAL const PDupwindNthSymm2A = PDupwindNthSymm2(&A[index]);
    CCTK_REAL const PDupwindNthAnti3A = PDupwindNthAnti3(&A[index]);
    CCTK_REAL const PDupwindNthSymm3A = PDupwindNthSymm3(&A[index]);
    CCTK_REAL const PDstandardNth1alpha = PDstandardNth1(&alpha[index]);
    CCTK_REAL const PDstandardNth2alpha = PDstandardNth2(&alpha[index]);
    CCTK_REAL const PDstandardNth3alpha = PDstandardNth3(&alpha[index]);
    CCTK_REAL const PDstandardNth11alpha = PDstandardNth11(&alpha[index]);
    CCTK_REAL const PDstandardNth22alpha = PDstandardNth22(&alpha[index]);
    CCTK_REAL const PDstandardNth33alpha = PDstandardNth33(&alpha[index]);
    CCTK_REAL const PDstandardNth12alpha = PDstandardNth12(&alpha[index]);
    CCTK_REAL const PDstandardNth13alpha = PDstandardNth13(&alpha[index]);
    CCTK_REAL const PDstandardNth23alpha = PDstandardNth23(&alpha[index]);
    CCTK_REAL const PDdissipationNth1alpha = PDdissipationNth1(&alpha[index]);
    CCTK_REAL const PDdissipationNth2alpha = PDdissipationNth2(&alpha[index]);
    CCTK_REAL const PDdissipationNth3alpha = PDdissipationNth3(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti1alpha = PDupwindNthAnti1(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm1alpha = PDupwindNthSymm1(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti2alpha = PDupwindNthAnti2(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm2alpha = PDupwindNthSymm2(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti3alpha = PDupwindNthAnti3(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm3alpha = PDupwindNthSymm3(&alpha[index]);
    CCTK_REAL const PDdissipationNth1B1 = PDdissipationNth1(&B1[index]);
    CCTK_REAL const PDdissipationNth2B1 = PDdissipationNth2(&B1[index]);
    CCTK_REAL const PDdissipationNth3B1 = PDdissipationNth3(&B1[index]);
    CCTK_REAL const PDupwindNthAnti1B1 = PDupwindNthAnti1(&B1[index]);
    CCTK_REAL const PDupwindNthSymm1B1 = PDupwindNthSymm1(&B1[index]);
    CCTK_REAL const PDupwindNthAnti2B1 = PDupwindNthAnti2(&B1[index]);
    CCTK_REAL const PDupwindNthSymm2B1 = PDupwindNthSymm2(&B1[index]);
    CCTK_REAL const PDupwindNthAnti3B1 = PDupwindNthAnti3(&B1[index]);
    CCTK_REAL const PDupwindNthSymm3B1 = PDupwindNthSymm3(&B1[index]);
    CCTK_REAL const PDdissipationNth1B2 = PDdissipationNth1(&B2[index]);
    CCTK_REAL const PDdissipationNth2B2 = PDdissipationNth2(&B2[index]);
    CCTK_REAL const PDdissipationNth3B2 = PDdissipationNth3(&B2[index]);
    CCTK_REAL const PDupwindNthAnti1B2 = PDupwindNthAnti1(&B2[index]);
    CCTK_REAL const PDupwindNthSymm1B2 = PDupwindNthSymm1(&B2[index]);
    CCTK_REAL const PDupwindNthAnti2B2 = PDupwindNthAnti2(&B2[index]);
    CCTK_REAL const PDupwindNthSymm2B2 = PDupwindNthSymm2(&B2[index]);
    CCTK_REAL const PDupwindNthAnti3B2 = PDupwindNthAnti3(&B2[index]);
    CCTK_REAL const PDupwindNthSymm3B2 = PDupwindNthSymm3(&B2[index]);
    CCTK_REAL const PDdissipationNth1B3 = PDdissipationNth1(&B3[index]);
    CCTK_REAL const PDdissipationNth2B3 = PDdissipationNth2(&B3[index]);
    CCTK_REAL const PDdissipationNth3B3 = PDdissipationNth3(&B3[index]);
    CCTK_REAL const PDupwindNthAnti1B3 = PDupwindNthAnti1(&B3[index]);
    CCTK_REAL const PDupwindNthSymm1B3 = PDupwindNthSymm1(&B3[index]);
    CCTK_REAL const PDupwindNthAnti2B3 = PDupwindNthAnti2(&B3[index]);
    CCTK_REAL const PDupwindNthSymm2B3 = PDupwindNthSymm2(&B3[index]);
    CCTK_REAL const PDupwindNthAnti3B3 = PDupwindNthAnti3(&B3[index]);
    CCTK_REAL const PDupwindNthSymm3B3 = PDupwindNthSymm3(&B3[index]);
    CCTK_REAL const PDstandardNth1beta1 = PDstandardNth1(&beta1[index]);
    CCTK_REAL const PDstandardNth2beta1 = PDstandardNth2(&beta1[index]);
    CCTK_REAL const PDstandardNth3beta1 = PDstandardNth3(&beta1[index]);
    CCTK_REAL const PDstandardNth11beta1 = PDstandardNth11(&beta1[index]);
    CCTK_REAL const PDstandardNth22beta1 = PDstandardNth22(&beta1[index]);
    CCTK_REAL const PDstandardNth33beta1 = PDstandardNth33(&beta1[index]);
    CCTK_REAL const PDstandardNth12beta1 = PDstandardNth12(&beta1[index]);
    CCTK_REAL const PDstandardNth13beta1 = PDstandardNth13(&beta1[index]);
    CCTK_REAL const PDstandardNth23beta1 = PDstandardNth23(&beta1[index]);
    CCTK_REAL const PDdissipationNth1beta1 = PDdissipationNth1(&beta1[index]);
    CCTK_REAL const PDdissipationNth2beta1 = PDdissipationNth2(&beta1[index]);
    CCTK_REAL const PDdissipationNth3beta1 = PDdissipationNth3(&beta1[index]);
    CCTK_REAL const PDupwindNthAnti1beta1 = PDupwindNthAnti1(&beta1[index]);
    CCTK_REAL const PDupwindNthSymm1beta1 = PDupwindNthSymm1(&beta1[index]);
    CCTK_REAL const PDupwindNthAnti2beta1 = PDupwindNthAnti2(&beta1[index]);
    CCTK_REAL const PDupwindNthSymm2beta1 = PDupwindNthSymm2(&beta1[index]);
    CCTK_REAL const PDupwindNthAnti3beta1 = PDupwindNthAnti3(&beta1[index]);
    CCTK_REAL const PDupwindNthSymm3beta1 = PDupwindNthSymm3(&beta1[index]);
    CCTK_REAL const PDstandardNth1beta2 = PDstandardNth1(&beta2[index]);
    CCTK_REAL const PDstandardNth2beta2 = PDstandardNth2(&beta2[index]);
    CCTK_REAL const PDstandardNth3beta2 = PDstandardNth3(&beta2[index]);
    CCTK_REAL const PDstandardNth11beta2 = PDstandardNth11(&beta2[index]);
    CCTK_REAL const PDstandardNth22beta2 = PDstandardNth22(&beta2[index]);
    CCTK_REAL const PDstandardNth33beta2 = PDstandardNth33(&beta2[index]);
    CCTK_REAL const PDstandardNth12beta2 = PDstandardNth12(&beta2[index]);
    CCTK_REAL const PDstandardNth13beta2 = PDstandardNth13(&beta2[index]);
    CCTK_REAL const PDstandardNth23beta2 = PDstandardNth23(&beta2[index]);
    CCTK_REAL const PDdissipationNth1beta2 = PDdissipationNth1(&beta2[index]);
    CCTK_REAL const PDdissipationNth2beta2 = PDdissipationNth2(&beta2[index]);
    CCTK_REAL const PDdissipationNth3beta2 = PDdissipationNth3(&beta2[index]);
    CCTK_REAL const PDupwindNthAnti1beta2 = PDupwindNthAnti1(&beta2[index]);
    CCTK_REAL const PDupwindNthSymm1beta2 = PDupwindNthSymm1(&beta2[index]);
    CCTK_REAL const PDupwindNthAnti2beta2 = PDupwindNthAnti2(&beta2[index]);
    CCTK_REAL const PDupwindNthSymm2beta2 = PDupwindNthSymm2(&beta2[index]);
    CCTK_REAL const PDupwindNthAnti3beta2 = PDupwindNthAnti3(&beta2[index]);
    CCTK_REAL const PDupwindNthSymm3beta2 = PDupwindNthSymm3(&beta2[index]);
    CCTK_REAL const PDstandardNth1beta3 = PDstandardNth1(&beta3[index]);
    CCTK_REAL const PDstandardNth2beta3 = PDstandardNth2(&beta3[index]);
    CCTK_REAL const PDstandardNth3beta3 = PDstandardNth3(&beta3[index]);
    CCTK_REAL const PDstandardNth11beta3 = PDstandardNth11(&beta3[index]);
    CCTK_REAL const PDstandardNth22beta3 = PDstandardNth22(&beta3[index]);
    CCTK_REAL const PDstandardNth33beta3 = PDstandardNth33(&beta3[index]);
    CCTK_REAL const PDstandardNth12beta3 = PDstandardNth12(&beta3[index]);
    CCTK_REAL const PDstandardNth13beta3 = PDstandardNth13(&beta3[index]);
    CCTK_REAL const PDstandardNth23beta3 = PDstandardNth23(&beta3[index]);
    CCTK_REAL const PDdissipationNth1beta3 = PDdissipationNth1(&beta3[index]);
    CCTK_REAL const PDdissipationNth2beta3 = PDdissipationNth2(&beta3[index]);
    CCTK_REAL const PDdissipationNth3beta3 = PDdissipationNth3(&beta3[index]);
    CCTK_REAL const PDupwindNthAnti1beta3 = PDupwindNthAnti1(&beta3[index]);
    CCTK_REAL const PDupwindNthSymm1beta3 = PDupwindNthSymm1(&beta3[index]);
    CCTK_REAL const PDupwindNthAnti2beta3 = PDupwindNthAnti2(&beta3[index]);
    CCTK_REAL const PDupwindNthSymm2beta3 = PDupwindNthSymm2(&beta3[index]);
    CCTK_REAL const PDupwindNthAnti3beta3 = PDupwindNthAnti3(&beta3[index]);
    CCTK_REAL const PDupwindNthSymm3beta3 = PDupwindNthSymm3(&beta3[index]);
    CCTK_REAL const PDstandardNth1gt11 = PDstandardNth1(&gt11[index]);
    CCTK_REAL const PDstandardNth2gt11 = PDstandardNth2(&gt11[index]);
    CCTK_REAL const PDstandardNth3gt11 = PDstandardNth3(&gt11[index]);
    CCTK_REAL const PDdissipationNth1gt11 = PDdissipationNth1(&gt11[index]);
    CCTK_REAL const PDdissipationNth2gt11 = PDdissipationNth2(&gt11[index]);
    CCTK_REAL const PDdissipationNth3gt11 = PDdissipationNth3(&gt11[index]);
    CCTK_REAL const PDupwindNthAnti1gt11 = PDupwindNthAnti1(&gt11[index]);
    CCTK_REAL const PDupwindNthSymm1gt11 = PDupwindNthSymm1(&gt11[index]);
    CCTK_REAL const PDupwindNthAnti2gt11 = PDupwindNthAnti2(&gt11[index]);
    CCTK_REAL const PDupwindNthSymm2gt11 = PDupwindNthSymm2(&gt11[index]);
    CCTK_REAL const PDupwindNthAnti3gt11 = PDupwindNthAnti3(&gt11[index]);
    CCTK_REAL const PDupwindNthSymm3gt11 = PDupwindNthSymm3(&gt11[index]);
    CCTK_REAL const PDstandardNth1gt12 = PDstandardNth1(&gt12[index]);
    CCTK_REAL const PDstandardNth2gt12 = PDstandardNth2(&gt12[index]);
    CCTK_REAL const PDstandardNth3gt12 = PDstandardNth3(&gt12[index]);
    CCTK_REAL const PDdissipationNth1gt12 = PDdissipationNth1(&gt12[index]);
    CCTK_REAL const PDdissipationNth2gt12 = PDdissipationNth2(&gt12[index]);
    CCTK_REAL const PDdissipationNth3gt12 = PDdissipationNth3(&gt12[index]);
    CCTK_REAL const PDupwindNthAnti1gt12 = PDupwindNthAnti1(&gt12[index]);
    CCTK_REAL const PDupwindNthSymm1gt12 = PDupwindNthSymm1(&gt12[index]);
    CCTK_REAL const PDupwindNthAnti2gt12 = PDupwindNthAnti2(&gt12[index]);
    CCTK_REAL const PDupwindNthSymm2gt12 = PDupwindNthSymm2(&gt12[index]);
    CCTK_REAL const PDupwindNthAnti3gt12 = PDupwindNthAnti3(&gt12[index]);
    CCTK_REAL const PDupwindNthSymm3gt12 = PDupwindNthSymm3(&gt12[index]);
    CCTK_REAL const PDstandardNth1gt13 = PDstandardNth1(&gt13[index]);
    CCTK_REAL const PDstandardNth2gt13 = PDstandardNth2(&gt13[index]);
    CCTK_REAL const PDstandardNth3gt13 = PDstandardNth3(&gt13[index]);
    CCTK_REAL const PDdissipationNth1gt13 = PDdissipationNth1(&gt13[index]);
    CCTK_REAL const PDdissipationNth2gt13 = PDdissipationNth2(&gt13[index]);
    CCTK_REAL const PDdissipationNth3gt13 = PDdissipationNth3(&gt13[index]);
    CCTK_REAL const PDupwindNthAnti1gt13 = PDupwindNthAnti1(&gt13[index]);
    CCTK_REAL const PDupwindNthSymm1gt13 = PDupwindNthSymm1(&gt13[index]);
    CCTK_REAL const PDupwindNthAnti2gt13 = PDupwindNthAnti2(&gt13[index]);
    CCTK_REAL const PDupwindNthSymm2gt13 = PDupwindNthSymm2(&gt13[index]);
    CCTK_REAL const PDupwindNthAnti3gt13 = PDupwindNthAnti3(&gt13[index]);
    CCTK_REAL const PDupwindNthSymm3gt13 = PDupwindNthSymm3(&gt13[index]);
    CCTK_REAL const PDstandardNth1gt22 = PDstandardNth1(&gt22[index]);
    CCTK_REAL const PDstandardNth2gt22 = PDstandardNth2(&gt22[index]);
    CCTK_REAL const PDstandardNth3gt22 = PDstandardNth3(&gt22[index]);
    CCTK_REAL const PDdissipationNth1gt22 = PDdissipationNth1(&gt22[index]);
    CCTK_REAL const PDdissipationNth2gt22 = PDdissipationNth2(&gt22[index]);
    CCTK_REAL const PDdissipationNth3gt22 = PDdissipationNth3(&gt22[index]);
    CCTK_REAL const PDupwindNthAnti1gt22 = PDupwindNthAnti1(&gt22[index]);
    CCTK_REAL const PDupwindNthSymm1gt22 = PDupwindNthSymm1(&gt22[index]);
    CCTK_REAL const PDupwindNthAnti2gt22 = PDupwindNthAnti2(&gt22[index]);
    CCTK_REAL const PDupwindNthSymm2gt22 = PDupwindNthSymm2(&gt22[index]);
    CCTK_REAL const PDupwindNthAnti3gt22 = PDupwindNthAnti3(&gt22[index]);
    CCTK_REAL const PDupwindNthSymm3gt22 = PDupwindNthSymm3(&gt22[index]);
    CCTK_REAL const PDstandardNth1gt23 = PDstandardNth1(&gt23[index]);
    CCTK_REAL const PDstandardNth2gt23 = PDstandardNth2(&gt23[index]);
    CCTK_REAL const PDstandardNth3gt23 = PDstandardNth3(&gt23[index]);
    CCTK_REAL const PDdissipationNth1gt23 = PDdissipationNth1(&gt23[index]);
    CCTK_REAL const PDdissipationNth2gt23 = PDdissipationNth2(&gt23[index]);
    CCTK_REAL const PDdissipationNth3gt23 = PDdissipationNth3(&gt23[index]);
    CCTK_REAL const PDupwindNthAnti1gt23 = PDupwindNthAnti1(&gt23[index]);
    CCTK_REAL const PDupwindNthSymm1gt23 = PDupwindNthSymm1(&gt23[index]);
    CCTK_REAL const PDupwindNthAnti2gt23 = PDupwindNthAnti2(&gt23[index]);
    CCTK_REAL const PDupwindNthSymm2gt23 = PDupwindNthSymm2(&gt23[index]);
    CCTK_REAL const PDupwindNthAnti3gt23 = PDupwindNthAnti3(&gt23[index]);
    CCTK_REAL const PDupwindNthSymm3gt23 = PDupwindNthSymm3(&gt23[index]);
    CCTK_REAL const PDstandardNth1gt33 = PDstandardNth1(&gt33[index]);
    CCTK_REAL const PDstandardNth2gt33 = PDstandardNth2(&gt33[index]);
    CCTK_REAL const PDstandardNth3gt33 = PDstandardNth3(&gt33[index]);
    CCTK_REAL const PDdissipationNth1gt33 = PDdissipationNth1(&gt33[index]);
    CCTK_REAL const PDdissipationNth2gt33 = PDdissipationNth2(&gt33[index]);
    CCTK_REAL const PDdissipationNth3gt33 = PDdissipationNth3(&gt33[index]);
    CCTK_REAL const PDupwindNthAnti1gt33 = PDupwindNthAnti1(&gt33[index]);
    CCTK_REAL const PDupwindNthSymm1gt33 = PDupwindNthSymm1(&gt33[index]);
    CCTK_REAL const PDupwindNthAnti2gt33 = PDupwindNthAnti2(&gt33[index]);
    CCTK_REAL const PDupwindNthSymm2gt33 = PDupwindNthSymm2(&gt33[index]);
    CCTK_REAL const PDupwindNthAnti3gt33 = PDupwindNthAnti3(&gt33[index]);
    CCTK_REAL const PDupwindNthSymm3gt33 = PDupwindNthSymm3(&gt33[index]);
    CCTK_REAL const PDstandardNth1phi = PDstandardNth1(&phi[index]);
    CCTK_REAL const PDstandardNth2phi = PDstandardNth2(&phi[index]);
    CCTK_REAL const PDstandardNth3phi = PDstandardNth3(&phi[index]);
    CCTK_REAL const PDdissipationNth1phi = PDdissipationNth1(&phi[index]);
    CCTK_REAL const PDdissipationNth2phi = PDdissipationNth2(&phi[index]);
    CCTK_REAL const PDdissipationNth3phi = PDdissipationNth3(&phi[index]);
    CCTK_REAL const PDupwindNthAnti1phi = PDupwindNthAnti1(&phi[index]);
    CCTK_REAL const PDupwindNthSymm1phi = PDupwindNthSymm1(&phi[index]);
    CCTK_REAL const PDupwindNthAnti2phi = PDupwindNthAnti2(&phi[index]);
    CCTK_REAL const PDupwindNthSymm2phi = PDupwindNthSymm2(&phi[index]);
    CCTK_REAL const PDupwindNthAnti3phi = PDupwindNthAnti3(&phi[index]);
    CCTK_REAL const PDupwindNthSymm3phi = PDupwindNthSymm3(&phi[index]);
    CCTK_REAL const PDstandardNth1trK = PDstandardNth1(&trK[index]);
    CCTK_REAL const PDstandardNth2trK = PDstandardNth2(&trK[index]);
    CCTK_REAL const PDstandardNth3trK = PDstandardNth3(&trK[index]);
    CCTK_REAL const PDdissipationNth1trK = PDdissipationNth1(&trK[index]);
    CCTK_REAL const PDdissipationNth2trK = PDdissipationNth2(&trK[index]);
    CCTK_REAL const PDdissipationNth3trK = PDdissipationNth3(&trK[index]);
    CCTK_REAL const PDupwindNthAnti1trK = PDupwindNthAnti1(&trK[index]);
    CCTK_REAL const PDupwindNthSymm1trK = PDupwindNthSymm1(&trK[index]);
    CCTK_REAL const PDupwindNthAnti2trK = PDupwindNthAnti2(&trK[index]);
    CCTK_REAL const PDupwindNthSymm2trK = PDupwindNthSymm2(&trK[index]);
    CCTK_REAL const PDupwindNthAnti3trK = PDupwindNthAnti3(&trK[index]);
    CCTK_REAL const PDupwindNthSymm3trK = PDupwindNthSymm3(&trK[index]);
    CCTK_REAL const PDdissipationNth1Xt1 = PDdissipationNth1(&Xt1[index]);
    CCTK_REAL const PDdissipationNth2Xt1 = PDdissipationNth2(&Xt1[index]);
    CCTK_REAL const PDdissipationNth3Xt1 = PDdissipationNth3(&Xt1[index]);
    CCTK_REAL const PDupwindNthAnti1Xt1 = PDupwindNthAnti1(&Xt1[index]);
    CCTK_REAL const PDupwindNthSymm1Xt1 = PDupwindNthSymm1(&Xt1[index]);
    CCTK_REAL const PDupwindNthAnti2Xt1 = PDupwindNthAnti2(&Xt1[index]);
    CCTK_REAL const PDupwindNthSymm2Xt1 = PDupwindNthSymm2(&Xt1[index]);
    CCTK_REAL const PDupwindNthAnti3Xt1 = PDupwindNthAnti3(&Xt1[index]);
    CCTK_REAL const PDupwindNthSymm3Xt1 = PDupwindNthSymm3(&Xt1[index]);
    CCTK_REAL const PDdissipationNth1Xt2 = PDdissipationNth1(&Xt2[index]);
    CCTK_REAL const PDdissipationNth2Xt2 = PDdissipationNth2(&Xt2[index]);
    CCTK_REAL const PDdissipationNth3Xt2 = PDdissipationNth3(&Xt2[index]);
    CCTK_REAL const PDupwindNthAnti1Xt2 = PDupwindNthAnti1(&Xt2[index]);
    CCTK_REAL const PDupwindNthSymm1Xt2 = PDupwindNthSymm1(&Xt2[index]);
    CCTK_REAL const PDupwindNthAnti2Xt2 = PDupwindNthAnti2(&Xt2[index]);
    CCTK_REAL const PDupwindNthSymm2Xt2 = PDupwindNthSymm2(&Xt2[index]);
    CCTK_REAL const PDupwindNthAnti3Xt2 = PDupwindNthAnti3(&Xt2[index]);
    CCTK_REAL const PDupwindNthSymm3Xt2 = PDupwindNthSymm3(&Xt2[index]);
    CCTK_REAL const PDdissipationNth1Xt3 = PDdissipationNth1(&Xt3[index]);
    CCTK_REAL const PDdissipationNth2Xt3 = PDdissipationNth2(&Xt3[index]);
    CCTK_REAL const PDdissipationNth3Xt3 = PDdissipationNth3(&Xt3[index]);
    CCTK_REAL const PDupwindNthAnti1Xt3 = PDupwindNthAnti1(&Xt3[index]);
    CCTK_REAL const PDupwindNthSymm1Xt3 = PDupwindNthSymm1(&Xt3[index]);
    CCTK_REAL const PDupwindNthAnti2Xt3 = PDupwindNthAnti2(&Xt3[index]);
    CCTK_REAL const PDupwindNthSymm2Xt3 = PDupwindNthSymm2(&Xt3[index]);
    CCTK_REAL const PDupwindNthAnti3Xt3 = PDupwindNthAnti3(&Xt3[index]);
    CCTK_REAL const PDupwindNthSymm3Xt3 = PDupwindNthSymm3(&Xt3[index]);
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    CCTK_REAL epsdiss1 = ToReal(EpsDiss);
    
    CCTK_REAL epsdiss2 = ToReal(EpsDiss);
    
    CCTK_REAL epsdiss3 = ToReal(EpsDiss);
    
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu12 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu13 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu23 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL Gtl111 = 0.5*PDstandardNth1gt11;
    
    CCTK_REAL Gtl112 = 0.5*PDstandardNth2gt11;
    
    CCTK_REAL Gtl113 = 0.5*PDstandardNth3gt11;
    
    CCTK_REAL Gtl122 = -0.5*PDstandardNth1gt22 + PDstandardNth2gt12;
    
    CCTK_REAL Gtl123 = 0.5*(-PDstandardNth1gt23 + PDstandardNth2gt13 + 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl133 = -0.5*PDstandardNth1gt33 + PDstandardNth3gt13;
    
    CCTK_REAL Gtl211 = PDstandardNth1gt12 - 0.5*PDstandardNth2gt11;
    
    CCTK_REAL Gtl212 = 0.5*PDstandardNth1gt22;
    
    CCTK_REAL Gtl213 = 0.5*(PDstandardNth1gt23 - PDstandardNth2gt13 + 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl222 = 0.5*PDstandardNth2gt22;
    
    CCTK_REAL Gtl223 = 0.5*PDstandardNth3gt22;
    
    CCTK_REAL Gtl233 = -0.5*PDstandardNth2gt33 + PDstandardNth3gt23;
    
    CCTK_REAL Gtl311 = PDstandardNth1gt13 - 0.5*PDstandardNth3gt11;
    
    CCTK_REAL Gtl312 = 0.5*(PDstandardNth1gt23 + PDstandardNth2gt13 - 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl313 = 0.5*PDstandardNth1gt33;
    
    CCTK_REAL Gtl322 = PDstandardNth2gt23 - 0.5*PDstandardNth3gt22;
    
    CCTK_REAL Gtl323 = 0.5*PDstandardNth2gt33;
    
    CCTK_REAL Gtl333 = 0.5*PDstandardNth3gt33;
    
    CCTK_REAL Gt111 = Gtl111*gtu11 + Gtl211*gtu12 + Gtl311*gtu13;
    
    CCTK_REAL Gt211 = Gtl111*gtu12 + Gtl211*gtu22 + Gtl311*gtu23;
    
    CCTK_REAL Gt311 = Gtl111*gtu13 + Gtl211*gtu23 + Gtl311*gtu33;
    
    CCTK_REAL Gt112 = Gtl112*gtu11 + Gtl212*gtu12 + Gtl312*gtu13;
    
    CCTK_REAL Gt212 = Gtl112*gtu12 + Gtl212*gtu22 + Gtl312*gtu23;
    
    CCTK_REAL Gt312 = Gtl112*gtu13 + Gtl212*gtu23 + Gtl312*gtu33;
    
    CCTK_REAL Gt113 = Gtl113*gtu11 + Gtl213*gtu12 + Gtl313*gtu13;
    
    CCTK_REAL Gt213 = Gtl113*gtu12 + Gtl213*gtu22 + Gtl313*gtu23;
    
    CCTK_REAL Gt313 = Gtl113*gtu13 + Gtl213*gtu23 + Gtl313*gtu33;
    
    CCTK_REAL Gt122 = Gtl122*gtu11 + Gtl222*gtu12 + Gtl322*gtu13;
    
    CCTK_REAL Gt222 = Gtl122*gtu12 + Gtl222*gtu22 + Gtl322*gtu23;
    
    CCTK_REAL Gt322 = Gtl122*gtu13 + Gtl222*gtu23 + Gtl322*gtu33;
    
    CCTK_REAL Gt123 = Gtl123*gtu11 + Gtl223*gtu12 + Gtl323*gtu13;
    
    CCTK_REAL Gt223 = Gtl123*gtu12 + Gtl223*gtu22 + Gtl323*gtu23;
    
    CCTK_REAL Gt323 = Gtl123*gtu13 + Gtl223*gtu23 + Gtl323*gtu33;
    
    CCTK_REAL Gt133 = Gtl133*gtu11 + Gtl233*gtu12 + Gtl333*gtu13;
    
    CCTK_REAL Gt233 = Gtl133*gtu12 + Gtl233*gtu22 + Gtl333*gtu23;
    
    CCTK_REAL Gt333 = Gtl133*gtu13 + Gtl233*gtu23 + Gtl333*gtu33;
    
    CCTK_REAL Xtn1 = Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu12 + 
      Gt113*gtu13 + Gt123*gtu23) + Gt133*gtu33;
    
    CCTK_REAL Xtn2 = Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu12 + 
      Gt213*gtu13 + Gt223*gtu23) + Gt233*gtu33;
    
    CCTK_REAL Xtn3 = Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu12 + 
      Gt313*gtu13 + Gt323*gtu23) + Gt333*gtu33;
    
    CCTK_REAL fac1 = IfThen(ToReal(conformalMethod),-0.5*INV(phiL),1);
    
    CCTK_REAL cdphi1 = fac1*PDstandardNth1phi;
    
    CCTK_REAL cdphi2 = fac1*PDstandardNth2phi;
    
    CCTK_REAL cdphi3 = fac1*PDstandardNth3phi;
    
    CCTK_REAL Atm11 = At11L*gtu11 + At12L*gtu12 + At13L*gtu13;
    
    CCTK_REAL Atm21 = At11L*gtu12 + At12L*gtu22 + At13L*gtu23;
    
    CCTK_REAL Atm31 = At11L*gtu13 + At12L*gtu23 + At13L*gtu33;
    
    CCTK_REAL Atm12 = At12L*gtu11 + At22L*gtu12 + At23L*gtu13;
    
    CCTK_REAL Atm22 = At12L*gtu12 + At22L*gtu22 + At23L*gtu23;
    
    CCTK_REAL Atm32 = At12L*gtu13 + At22L*gtu23 + At23L*gtu33;
    
    CCTK_REAL Atm13 = At13L*gtu11 + At23L*gtu12 + At33L*gtu13;
    
    CCTK_REAL Atm23 = At13L*gtu12 + At23L*gtu22 + At33L*gtu23;
    
    CCTK_REAL Atm33 = At13L*gtu13 + At23L*gtu23 + At33L*gtu33;
    
    CCTK_REAL Atu11 = Atm11*gtu11 + Atm12*gtu12 + Atm13*gtu13;
    
    CCTK_REAL Atu12 = Atm11*gtu12 + Atm12*gtu22 + Atm13*gtu23;
    
    CCTK_REAL Atu13 = Atm11*gtu13 + Atm12*gtu23 + Atm13*gtu33;
    
    CCTK_REAL Atu22 = Atm21*gtu12 + Atm22*gtu22 + Atm23*gtu23;
    
    CCTK_REAL Atu23 = Atm21*gtu13 + Atm22*gtu23 + Atm23*gtu33;
    
    CCTK_REAL Atu33 = Atm31*gtu13 + Atm32*gtu23 + Atm33*gtu33;
    
    CCTK_REAL e4phi = 
      IfThen(ToReal(conformalMethod),INV(SQR(phiL)),exp(4*phiL));
    
    CCTK_REAL em4phi = INV(e4phi);
    
    CCTK_REAL rho = INV(SQR(alphaL))*(eTttL - 2*(beta2L*eTtyL + 
      beta3L*eTtzL) + 2*(beta1L*(-eTtxL + beta2L*eTxyL + beta3L*eTxzL) + 
      beta2L*beta3L*eTyzL) + eTxxL*SQR(beta1L) + eTyyL*SQR(beta2L) + 
      eTzzL*SQR(beta3L));
    
    CCTK_REAL S1 = (-eTtxL + beta1L*eTxxL + beta2L*eTxyL + 
      beta3L*eTxzL)*INV(alphaL);
    
    CCTK_REAL S2 = (-eTtyL + beta1L*eTxyL + beta2L*eTyyL + 
      beta3L*eTyzL)*INV(alphaL);
    
    CCTK_REAL S3 = (-eTtzL + beta1L*eTxzL + beta2L*eTyzL + 
      beta3L*eTzzL)*INV(alphaL);
    
    CCTK_REAL trS = em4phi*(eTxxL*gtu11 + eTyyL*gtu22 + 2*(eTxyL*gtu12 + 
      eTxzL*gtu13 + eTyzL*gtu23) + eTzzL*gtu33);
    
    CCTK_REAL phirhsL = epsdiss1*PDdissipationNth1phi + 
      epsdiss2*PDdissipationNth2phi + epsdiss3*PDdissipationNth3phi + 
      beta1L*PDupwindNthAnti1phi + beta2L*PDupwindNthAnti2phi + 
      beta3L*PDupwindNthAnti3phi + PDupwindNthSymm1phi*Abs(beta1L) + 
      PDupwindNthSymm2phi*Abs(beta2L) + PDupwindNthSymm3phi*Abs(beta3L) + 
      IfThen(ToReal(conformalMethod),phiL*(-0.333333333333333333333333333333*(PDstandardNth1beta1 
      + PDstandardNth2beta2 + PDstandardNth3beta3) + 
      0.333333333333333333333333333333*alphaL*trKL),0.166666666666666666666666666667*(PDstandardNth1beta1 
      + PDstandardNth2beta2 + PDstandardNth3beta3) - 
      0.166666666666666666666666666667*alphaL*trKL);
    
    CCTK_REAL gt11rhsL = -2*alphaL*At11L + epsdiss1*PDdissipationNth1gt11 
      + epsdiss2*PDdissipationNth2gt11 + epsdiss3*PDdissipationNth3gt11 + 
      2*(gt12L*PDstandardNth1beta2 + gt13L*PDstandardNth1beta3) + 
      gt11L*(1.33333333333333333333333333333*PDstandardNth1beta1 - 
      0.666666666666666666666666666667*(PDstandardNth2beta2 + 
      PDstandardNth3beta3)) + beta1L*PDupwindNthAnti1gt11 + 
      beta2L*PDupwindNthAnti2gt11 + beta3L*PDupwindNthAnti3gt11 + 
      PDupwindNthSymm1gt11*Abs(beta1L) + PDupwindNthSymm2gt11*Abs(beta2L) + 
      PDupwindNthSymm3gt11*Abs(beta3L);
    
    CCTK_REAL gt12rhsL = -2*alphaL*At12L + epsdiss1*PDdissipationNth1gt12 
      + epsdiss2*PDdissipationNth2gt12 + epsdiss3*PDdissipationNth3gt12 + 
      gt22L*PDstandardNth1beta2 + gt23L*PDstandardNth1beta3 + 
      gt11L*PDstandardNth2beta1 + gt13L*PDstandardNth2beta3 + 
      gt12L*(0.333333333333333333333333333333*(PDstandardNth1beta1 + 
      PDstandardNth2beta2) - 
      0.666666666666666666666666666667*PDstandardNth3beta3) + 
      beta1L*PDupwindNthAnti1gt12 + beta2L*PDupwindNthAnti2gt12 + 
      beta3L*PDupwindNthAnti3gt12 + PDupwindNthSymm1gt12*Abs(beta1L) + 
      PDupwindNthSymm2gt12*Abs(beta2L) + PDupwindNthSymm3gt12*Abs(beta3L);
    
    CCTK_REAL gt13rhsL = -2*alphaL*At13L + epsdiss1*PDdissipationNth1gt13 
      + epsdiss2*PDdissipationNth2gt13 + epsdiss3*PDdissipationNth3gt13 + 
      gt23L*PDstandardNth1beta2 + gt33L*PDstandardNth1beta3 + 
      gt11L*PDstandardNth3beta1 + gt12L*PDstandardNth3beta2 + 
      gt13L*(-0.666666666666666666666666666667*PDstandardNth2beta2 + 
      0.333333333333333333333333333333*(PDstandardNth1beta1 + 
      PDstandardNth3beta3)) + beta1L*PDupwindNthAnti1gt13 + 
      beta2L*PDupwindNthAnti2gt13 + beta3L*PDupwindNthAnti3gt13 + 
      PDupwindNthSymm1gt13*Abs(beta1L) + PDupwindNthSymm2gt13*Abs(beta2L) + 
      PDupwindNthSymm3gt13*Abs(beta3L);
    
    CCTK_REAL gt22rhsL = -2*alphaL*At22L + epsdiss1*PDdissipationNth1gt22 
      + epsdiss2*PDdissipationNth2gt22 + epsdiss3*PDdissipationNth3gt22 + 
      2*(gt12L*PDstandardNth2beta1 + gt23L*PDstandardNth2beta3) + 
      gt22L*(1.33333333333333333333333333333*PDstandardNth2beta2 - 
      0.666666666666666666666666666667*(PDstandardNth1beta1 + 
      PDstandardNth3beta3)) + beta1L*PDupwindNthAnti1gt22 + 
      beta2L*PDupwindNthAnti2gt22 + beta3L*PDupwindNthAnti3gt22 + 
      PDupwindNthSymm1gt22*Abs(beta1L) + PDupwindNthSymm2gt22*Abs(beta2L) + 
      PDupwindNthSymm3gt22*Abs(beta3L);
    
    CCTK_REAL gt23rhsL = -2*alphaL*At23L + epsdiss1*PDdissipationNth1gt23 
      + epsdiss2*PDdissipationNth2gt23 + epsdiss3*PDdissipationNth3gt23 + 
      gt13L*PDstandardNth2beta1 + gt33L*PDstandardNth2beta3 + 
      gt12L*PDstandardNth3beta1 + gt22L*PDstandardNth3beta2 + 
      gt23L*(-0.666666666666666666666666666667*PDstandardNth1beta1 + 
      0.333333333333333333333333333333*(PDstandardNth2beta2 + 
      PDstandardNth3beta3)) + beta1L*PDupwindNthAnti1gt23 + 
      beta2L*PDupwindNthAnti2gt23 + beta3L*PDupwindNthAnti3gt23 + 
      PDupwindNthSymm1gt23*Abs(beta1L) + PDupwindNthSymm2gt23*Abs(beta2L) + 
      PDupwindNthSymm3gt23*Abs(beta3L);
    
    CCTK_REAL gt33rhsL = -2*alphaL*At33L + epsdiss1*PDdissipationNth1gt33 
      + epsdiss2*PDdissipationNth2gt33 + epsdiss3*PDdissipationNth3gt33 + 
      2*(gt13L*PDstandardNth3beta1 + gt23L*PDstandardNth3beta2) + 
      gt33L*(-0.666666666666666666666666666667*(PDstandardNth1beta1 + 
      PDstandardNth2beta2) + 
      1.33333333333333333333333333333*PDstandardNth3beta3) + 
      beta1L*PDupwindNthAnti1gt33 + beta2L*PDupwindNthAnti2gt33 + 
      beta3L*PDupwindNthAnti3gt33 + PDupwindNthSymm1gt33*Abs(beta1L) + 
      PDupwindNthSymm2gt33*Abs(beta2L) + PDupwindNthSymm3gt33*Abs(beta3L);
    
    CCTK_REAL dotXt1 = 
      0.333333333333333333333333333333*(7*(gtu12*PDstandardNth12beta1 + 
      gtu13*PDstandardNth13beta1) + gtu11*(4*PDstandardNth11beta1 + 
      PDstandardNth12beta2 + PDstandardNth13beta3) + 
      gtu12*(PDstandardNth22beta2 + PDstandardNth23beta3) + 
      gtu13*(PDstandardNth23beta2 + PDstandardNth33beta3) - 
      6*(Atu11*PDstandardNth1alpha + Atu12*PDstandardNth2alpha + 
      Atu13*PDstandardNth3alpha) + 6*(gtu23*PDstandardNth23beta1 + 
      alphaL*(6*(Atu11*cdphi1 + Atu12*cdphi2 + Atu13*cdphi3) + Atu11*Gt111 + 
      Atu22*Gt122 + 2*(Atu12*Gt112 + Atu13*Gt113 + Atu23*Gt123) + Atu33*Gt133 
      - 0.666666666666666666666666666667*(gtu11*PDstandardNth1trK + 
      gtu12*PDstandardNth2trK + gtu13*PDstandardNth3trK))) - 
      150.7964473723100754462068823974161384415*alphaL*(gtu11*S1 + gtu12*S2 + 
      gtu13*S3) + (-3*PDstandardNth1beta1 + 2*(PDstandardNth1beta1 + 
      PDstandardNth2beta2 + PDstandardNth3beta3))*Xtn1 - 
      3*(PDstandardNth2beta1*Xtn2 + PDstandardNth3beta1*Xtn3) + 
      3*(epsdiss1*PDdissipationNth1Xt1 + epsdiss2*PDdissipationNth2Xt1 + 
      epsdiss3*PDdissipationNth3Xt1 + gtu22*PDstandardNth22beta1 + 
      gtu33*PDstandardNth33beta1 + beta1L*PDupwindNthAnti1Xt1 + 
      beta2L*PDupwindNthAnti2Xt1 + beta3L*PDupwindNthAnti3Xt1 + 
      PDupwindNthSymm1Xt1*Abs(beta1L) + PDupwindNthSymm2Xt1*Abs(beta2L) + 
      PDupwindNthSymm3Xt1*Abs(beta3L)));
    
    CCTK_REAL dotXt2 = 
      0.333333333333333333333333333333*(gtu12*(PDstandardNth11beta1 + 
      7*PDstandardNth12beta2 + PDstandardNth13beta3) + 
      gtu22*(PDstandardNth12beta1 + 4*PDstandardNth22beta2 + 
      PDstandardNth23beta3) + gtu23*(PDstandardNth13beta1 + 
      7*PDstandardNth23beta2 + PDstandardNth33beta3) - 
      6*(Atu12*PDstandardNth1alpha + Atu22*PDstandardNth2alpha + 
      Atu23*PDstandardNth3alpha) + 6*(gtu13*PDstandardNth13beta2 + 
      alphaL*(6*(Atu12*cdphi1 + Atu22*cdphi2 + Atu23*cdphi3) + Atu11*Gt211 + 
      Atu22*Gt222 + 2*(Atu12*Gt212 + Atu13*Gt213 + Atu23*Gt223) + Atu33*Gt233 
      - 0.666666666666666666666666666667*(gtu12*PDstandardNth1trK + 
      gtu22*PDstandardNth2trK + gtu23*PDstandardNth3trK))) - 
      150.7964473723100754462068823974161384415*alphaL*(gtu12*S1 + gtu22*S2 + 
      gtu23*S3) + 2*(PDstandardNth1beta1 + PDstandardNth2beta2 + 
      PDstandardNth3beta3)*Xtn2 - 3*(PDstandardNth1beta2*Xtn1 + 
      PDstandardNth2beta2*Xtn2 + PDstandardNth3beta2*Xtn3) + 
      3*(epsdiss1*PDdissipationNth1Xt2 + epsdiss2*PDdissipationNth2Xt2 + 
      epsdiss3*PDdissipationNth3Xt2 + gtu11*PDstandardNth11beta2 + 
      gtu33*PDstandardNth33beta2 + beta1L*PDupwindNthAnti1Xt2 + 
      beta2L*PDupwindNthAnti2Xt2 + beta3L*PDupwindNthAnti3Xt2 + 
      PDupwindNthSymm1Xt2*Abs(beta1L) + PDupwindNthSymm2Xt2*Abs(beta2L) + 
      PDupwindNthSymm3Xt2*Abs(beta3L)));
    
    CCTK_REAL dotXt3 = 
      0.333333333333333333333333333333*(gtu13*(PDstandardNth11beta1 + 
      PDstandardNth12beta2 + 7*PDstandardNth13beta3) + 
      gtu23*(PDstandardNth12beta1 + PDstandardNth22beta2 + 
      7*PDstandardNth23beta3) + gtu33*(PDstandardNth13beta1 + 
      PDstandardNth23beta2 + 4*PDstandardNth33beta3) - 
      6*(Atu13*PDstandardNth1alpha + Atu23*PDstandardNth2alpha + 
      Atu33*PDstandardNth3alpha) + 6*(gtu12*PDstandardNth12beta3 + 
      alphaL*(6*(Atu13*cdphi1 + Atu23*cdphi2 + Atu33*cdphi3) + Atu11*Gt311 + 
      Atu22*Gt322 + 2*(Atu12*Gt312 + Atu13*Gt313 + Atu23*Gt323) + Atu33*Gt333 
      - 0.666666666666666666666666666667*(gtu13*PDstandardNth1trK + 
      gtu23*PDstandardNth2trK + gtu33*PDstandardNth3trK))) - 
      150.7964473723100754462068823974161384415*alphaL*(gtu13*S1 + gtu23*S2 + 
      gtu33*S3) + 2*(PDstandardNth1beta1 + PDstandardNth2beta2 + 
      PDstandardNth3beta3)*Xtn3 - 3*(PDstandardNth1beta3*Xtn1 + 
      PDstandardNth2beta3*Xtn2 + PDstandardNth3beta3*Xtn3) + 
      3*(epsdiss1*PDdissipationNth1Xt3 + epsdiss2*PDdissipationNth2Xt3 + 
      epsdiss3*PDdissipationNth3Xt3 + gtu11*PDstandardNth11beta3 + 
      gtu22*PDstandardNth22beta3 + beta1L*PDupwindNthAnti1Xt3 + 
      beta2L*PDupwindNthAnti2Xt3 + beta3L*PDupwindNthAnti3Xt3 + 
      PDupwindNthSymm1Xt3*Abs(beta1L) + PDupwindNthSymm2Xt3*Abs(beta2L) + 
      PDupwindNthSymm3Xt3*Abs(beta3L)));
    
    CCTK_REAL Xt1rhsL = dotXt1;
    
    CCTK_REAL Xt2rhsL = dotXt2;
    
    CCTK_REAL Xt3rhsL = dotXt3;
    
    CCTK_REAL dottrK = epsdiss1*PDdissipationNth1trK + 
      epsdiss2*PDdissipationNth2trK + epsdiss3*PDdissipationNth3trK + 
      beta1L*PDupwindNthAnti1trK + beta2L*PDupwindNthAnti2trK + 
      beta3L*PDupwindNthAnti3trK - em4phi*(gtu11*PDstandardNth11alpha + 
      gtu22*PDstandardNth22alpha + gtu33*(PDstandardNth33alpha + 
      2*cdphi3*PDstandardNth3alpha) + 2*(gtu12*PDstandardNth12alpha + 
      gtu13*(PDstandardNth13alpha + cdphi1*PDstandardNth3alpha) + 
      gtu23*(PDstandardNth23alpha + cdphi2*PDstandardNth3alpha)) + 
      PDstandardNth1alpha*(2*(cdphi1*gtu11 + cdphi2*gtu12 + cdphi3*gtu13) - 
      Xtn1) + PDstandardNth2alpha*(2*(cdphi1*gtu12 + cdphi2*gtu22 + 
      cdphi3*gtu23) - Xtn2) - PDstandardNth3alpha*Xtn3) + 
      PDupwindNthSymm1trK*Abs(beta1L) + PDupwindNthSymm2trK*Abs(beta2L) + 
      PDupwindNthSymm3trK*Abs(beta3L) + alphaL*(2*(Atm12*Atm21 + Atm13*Atm31 
      + Atm23*Atm32) + 12.56637061435917295385057353311801153679*(rho + trS) 
      + SQR(Atm11) + SQR(Atm22) + SQR(Atm33) + 
      0.333333333333333333333333333333*SQR(trKL));
    
    CCTK_REAL trKrhsL = dottrK;
    
    CCTK_REAL alpharhsL = epsdiss1*PDdissipationNth1alpha + 
      epsdiss2*PDdissipationNth2alpha + epsdiss3*PDdissipationNth3alpha - 
      pow(alphaL,ToReal(harmonicN))*ToReal(harmonicF)*(trKL + (AL - 
      trKL)*ToReal(LapseACoeff)) + (beta1L*PDupwindNthAnti1alpha + 
      beta2L*PDupwindNthAnti2alpha + beta3L*PDupwindNthAnti3alpha + 
      PDupwindNthSymm1alpha*Abs(beta1L) + PDupwindNthSymm2alpha*Abs(beta2L) + 
      PDupwindNthSymm3alpha*Abs(beta3L))*ToReal(LapseAdvectionCoeff);
    
    CCTK_REAL ArhsL = epsdiss1*PDdissipationNth1A + 
      epsdiss2*PDdissipationNth2A + epsdiss3*PDdissipationNth3A + (dottrK - 
      AL*ToReal(AlphaDriver))*ToReal(LapseACoeff) + (beta1L*PDupwindNthAnti1A 
      + beta2L*PDupwindNthAnti2A + beta3L*PDupwindNthAnti3A + 
      PDupwindNthSymm1A*Abs(beta1L) + PDupwindNthSymm2A*Abs(beta2L) + 
      PDupwindNthSymm3A*Abs(beta3L))*ToReal(LapseAdvectionCoeff);
    
    CCTK_REAL eta = fmin(1,INV(rL)*ToReal(SpatialBetaDriverRadius));
    
    CCTK_REAL theta = fmin(1,exp(1 - 
      rL*INV(ToReal(SpatialShiftGammaCoeffRadius))));
    
    CCTK_REAL beta1rhsL = epsdiss1*PDdissipationNth1beta1 + 
      epsdiss2*PDdissipationNth2beta1 + epsdiss3*PDdissipationNth3beta1 + 
      (beta1L*PDupwindNthAnti1beta1 + beta2L*PDupwindNthAnti2beta1 + 
      beta3L*PDupwindNthAnti3beta1 + PDupwindNthSymm1beta1*Abs(beta1L) + 
      PDupwindNthSymm2beta1*Abs(beta2L) + 
      PDupwindNthSymm3beta1*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      theta*(Xt1L + beta1L*eta*ToReal(BetaDriver)*(-1 + ToReal(ShiftBCoeff)) 
      + (B1L - Xt1L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    CCTK_REAL beta2rhsL = epsdiss1*PDdissipationNth1beta2 + 
      epsdiss2*PDdissipationNth2beta2 + epsdiss3*PDdissipationNth3beta2 + 
      (beta1L*PDupwindNthAnti1beta2 + beta2L*PDupwindNthAnti2beta2 + 
      beta3L*PDupwindNthAnti3beta2 + PDupwindNthSymm1beta2*Abs(beta1L) + 
      PDupwindNthSymm2beta2*Abs(beta2L) + 
      PDupwindNthSymm3beta2*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      theta*(Xt2L + beta2L*eta*ToReal(BetaDriver)*(-1 + ToReal(ShiftBCoeff)) 
      + (B2L - Xt2L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    CCTK_REAL beta3rhsL = epsdiss1*PDdissipationNth1beta3 + 
      epsdiss2*PDdissipationNth2beta3 + epsdiss3*PDdissipationNth3beta3 + 
      (beta1L*PDupwindNthAnti1beta3 + beta2L*PDupwindNthAnti2beta3 + 
      beta3L*PDupwindNthAnti3beta3 + PDupwindNthSymm1beta3*Abs(beta1L) + 
      PDupwindNthSymm2beta3*Abs(beta2L) + 
      PDupwindNthSymm3beta3*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      theta*(Xt3L + beta3L*eta*ToReal(BetaDriver)*(-1 + ToReal(ShiftBCoeff)) 
      + (B3L - Xt3L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    CCTK_REAL B1rhsL = epsdiss1*PDdissipationNth1B1 + 
      epsdiss2*PDdissipationNth2B1 + epsdiss3*PDdissipationNth3B1 + 
      (beta1L*(PDupwindNthAnti1B1 - PDupwindNthAnti1Xt1) + 
      beta2L*(PDupwindNthAnti2B1 - PDupwindNthAnti2Xt1) + 
      beta3L*(PDupwindNthAnti3B1 - PDupwindNthAnti3Xt1) + (PDupwindNthSymm1B1 
      - PDupwindNthSymm1Xt1)*Abs(beta1L) + (PDupwindNthSymm2B1 - 
      PDupwindNthSymm2Xt1)*Abs(beta2L) + (PDupwindNthSymm3B1 - 
      PDupwindNthSymm3Xt1)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + (dotXt1 
      - B1L*eta*ToReal(BetaDriver))*ToReal(ShiftBCoeff);
    
    CCTK_REAL B2rhsL = epsdiss1*PDdissipationNth1B2 + 
      epsdiss2*PDdissipationNth2B2 + epsdiss3*PDdissipationNth3B2 + 
      (beta1L*(PDupwindNthAnti1B2 - PDupwindNthAnti1Xt2) + 
      beta2L*(PDupwindNthAnti2B2 - PDupwindNthAnti2Xt2) + 
      beta3L*(PDupwindNthAnti3B2 - PDupwindNthAnti3Xt2) + (PDupwindNthSymm1B2 
      - PDupwindNthSymm1Xt2)*Abs(beta1L) + (PDupwindNthSymm2B2 - 
      PDupwindNthSymm2Xt2)*Abs(beta2L) + (PDupwindNthSymm3B2 - 
      PDupwindNthSymm3Xt2)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + (dotXt2 
      - B2L*eta*ToReal(BetaDriver))*ToReal(ShiftBCoeff);
    
    CCTK_REAL B3rhsL = epsdiss1*PDdissipationNth1B3 + 
      epsdiss2*PDdissipationNth2B3 + epsdiss3*PDdissipationNth3B3 + 
      (beta1L*(PDupwindNthAnti1B3 - PDupwindNthAnti1Xt3) + 
      beta2L*(PDupwindNthAnti2B3 - PDupwindNthAnti2Xt3) + 
      beta3L*(PDupwindNthAnti3B3 - PDupwindNthAnti3Xt3) + (PDupwindNthSymm1B3 
      - PDupwindNthSymm1Xt3)*Abs(beta1L) + (PDupwindNthSymm2B3 - 
      PDupwindNthSymm2Xt3)*Abs(beta2L) + (PDupwindNthSymm3B3 - 
      PDupwindNthSymm3Xt3)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + (dotXt3 
      - B3L*eta*ToReal(BetaDriver))*ToReal(ShiftBCoeff);
    
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    Arhs[index] = ArhsL;
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
  LC_ENDLOOP3 (ML_BSSN_O8_RHS1);
}

extern "C" void ML_BSSN_O8_RHS1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_O8_RHS1_Body);
}
