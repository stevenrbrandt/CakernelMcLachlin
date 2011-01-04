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
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

void ML_BSSN_UPW_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_UPW::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_UPW::ML_trace_curvrhs.");
  return;
}

void ML_BSSN_UPW_RHS1_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_UPW_RHS1_Body");
  }
  
  if (cctk_iteration % ML_BSSN_UPW_RHS1_calc_every != ML_BSSN_UPW_RHS1_calc_offset)
  {
    return;
  }
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  CCTK_REAL const dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL const dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL const dz = CCTK_DELTA_SPACE(2);
  int const di = 1;
  int const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  int const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  CCTK_REAL const dxi = 1.0 / dx;
  CCTK_REAL const dyi = 1.0 / dy;
  CCTK_REAL const dzi = 1.0 / dz;
  CCTK_REAL const khalf = 0.5;
  CCTK_REAL const kthird = 1/3.0;
  CCTK_REAL const ktwothird = 2.0/3.0;
  CCTK_REAL const kfourthird = 4.0/3.0;
  CCTK_REAL const keightthird = 8.0/3.0;
  CCTK_REAL const hdxi = 0.5 * dxi;
  CCTK_REAL const hdyi = 0.5 * dyi;
  CCTK_REAL const hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  CCTK_REAL const p1o12dx = INV(dx)/12.;
  CCTK_REAL const p1o12dy = INV(dy)/12.;
  CCTK_REAL const p1o12dz = INV(dz)/12.;
  CCTK_REAL const p1o144dxdy = (INV(dx)*INV(dy))/144.;
  CCTK_REAL const p1o144dxdz = (INV(dx)*INV(dz))/144.;
  CCTK_REAL const p1o144dydz = (INV(dy)*INV(dz))/144.;
  CCTK_REAL const p1o24dx = INV(dx)/24.;
  CCTK_REAL const p1o24dy = INV(dy)/24.;
  CCTK_REAL const p1o24dz = INV(dz)/24.;
  CCTK_REAL const p1o64dx = INV(dx)/64.;
  CCTK_REAL const p1o64dy = INV(dy)/64.;
  CCTK_REAL const p1o64dz = INV(dz)/64.;
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o12dx2 = -pow(dx,-2)/12.;
  CCTK_REAL const pm1o12dy2 = -pow(dy,-2)/12.;
  CCTK_REAL const pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_UPW_RHS1,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL PDdissipationNth1A = INITVALUE;
    // CCTK_REAL PDdissipationNth2A = INITVALUE;
    // CCTK_REAL PDdissipationNth3A = INITVALUE;
    // CCTK_REAL PDstandardNth1alpha = INITVALUE;
    // CCTK_REAL PDstandardNth2alpha = INITVALUE;
    // CCTK_REAL PDstandardNth3alpha = INITVALUE;
    // CCTK_REAL PDstandardNth11alpha = INITVALUE;
    // CCTK_REAL PDstandardNth22alpha = INITVALUE;
    // CCTK_REAL PDstandardNth33alpha = INITVALUE;
    // CCTK_REAL PDstandardNth12alpha = INITVALUE;
    // CCTK_REAL PDstandardNth13alpha = INITVALUE;
    // CCTK_REAL PDstandardNth23alpha = INITVALUE;
    // CCTK_REAL PDdissipationNth1alpha = INITVALUE;
    // CCTK_REAL PDdissipationNth2alpha = INITVALUE;
    // CCTK_REAL PDdissipationNth3alpha = INITVALUE;
    // CCTK_REAL PDdissipationNth1B1 = INITVALUE;
    // CCTK_REAL PDdissipationNth2B1 = INITVALUE;
    // CCTK_REAL PDdissipationNth3B1 = INITVALUE;
    // CCTK_REAL PDdissipationNth1B2 = INITVALUE;
    // CCTK_REAL PDdissipationNth2B2 = INITVALUE;
    // CCTK_REAL PDdissipationNth3B2 = INITVALUE;
    // CCTK_REAL PDdissipationNth1B3 = INITVALUE;
    // CCTK_REAL PDdissipationNth2B3 = INITVALUE;
    // CCTK_REAL PDdissipationNth3B3 = INITVALUE;
    // CCTK_REAL PDstandardNth1beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth2beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth3beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth11beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth22beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth33beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth12beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth13beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth23beta1 = INITVALUE;
    // CCTK_REAL PDdissipationNth1beta1 = INITVALUE;
    // CCTK_REAL PDdissipationNth2beta1 = INITVALUE;
    // CCTK_REAL PDdissipationNth3beta1 = INITVALUE;
    // CCTK_REAL PDstandardNth1beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth2beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth3beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth11beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth22beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth33beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth12beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth13beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth23beta2 = INITVALUE;
    // CCTK_REAL PDdissipationNth1beta2 = INITVALUE;
    // CCTK_REAL PDdissipationNth2beta2 = INITVALUE;
    // CCTK_REAL PDdissipationNth3beta2 = INITVALUE;
    // CCTK_REAL PDstandardNth1beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth2beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth3beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth11beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth22beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth33beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth12beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth13beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth23beta3 = INITVALUE;
    // CCTK_REAL PDdissipationNth1beta3 = INITVALUE;
    // CCTK_REAL PDdissipationNth2beta3 = INITVALUE;
    // CCTK_REAL PDdissipationNth3beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt11 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt11 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt11 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt12 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt12 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt12 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt13 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt13 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt13 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt22 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt22 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt22 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt23 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt23 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt23 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt33 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt33 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt33 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth1phi = INITVALUE;
    // CCTK_REAL PDstandardNth2phi = INITVALUE;
    // CCTK_REAL PDstandardNth3phi = INITVALUE;
    // CCTK_REAL PDdissipationNth1phi = INITVALUE;
    // CCTK_REAL PDdissipationNth2phi = INITVALUE;
    // CCTK_REAL PDdissipationNth3phi = INITVALUE;
    // CCTK_REAL PDstandardNth1trK = INITVALUE;
    // CCTK_REAL PDstandardNth2trK = INITVALUE;
    // CCTK_REAL PDstandardNth3trK = INITVALUE;
    // CCTK_REAL PDdissipationNth1trK = INITVALUE;
    // CCTK_REAL PDdissipationNth2trK = INITVALUE;
    // CCTK_REAL PDdissipationNth3trK = INITVALUE;
    // CCTK_REAL PDdissipationNth1Xt1 = INITVALUE;
    // CCTK_REAL PDdissipationNth2Xt1 = INITVALUE;
    // CCTK_REAL PDdissipationNth3Xt1 = INITVALUE;
    // CCTK_REAL PDdissipationNth1Xt2 = INITVALUE;
    // CCTK_REAL PDdissipationNth2Xt2 = INITVALUE;
    // CCTK_REAL PDdissipationNth3Xt2 = INITVALUE;
    // CCTK_REAL PDdissipationNth1Xt3 = INITVALUE;
    // CCTK_REAL PDdissipationNth2Xt3 = INITVALUE;
    // CCTK_REAL PDdissipationNth3Xt3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL  AL = A[index];
    CCTK_REAL  alphaL = alpha[index];
    CCTK_REAL  At11L = At11[index];
    CCTK_REAL  At12L = At12[index];
    CCTK_REAL  At13L = At13[index];
    CCTK_REAL  At22L = At22[index];
    CCTK_REAL  At23L = At23[index];
    CCTK_REAL  At33L = At33[index];
    CCTK_REAL  B1L = B1[index];
    CCTK_REAL  B2L = B2[index];
    CCTK_REAL  B3L = B3[index];
    CCTK_REAL  beta1L = beta1[index];
    CCTK_REAL  beta2L = beta2[index];
    CCTK_REAL  beta3L = beta3[index];
    CCTK_REAL  eTttL = (*stress_energy_state) ? (eTtt[index]) : 0.0;
    CCTK_REAL  eTtxL = (*stress_energy_state) ? (eTtx[index]) : 0.0;
    CCTK_REAL  eTtyL = (*stress_energy_state) ? (eTty[index]) : 0.0;
    CCTK_REAL  eTtzL = (*stress_energy_state) ? (eTtz[index]) : 0.0;
    CCTK_REAL  eTxxL = (*stress_energy_state) ? (eTxx[index]) : 0.0;
    CCTK_REAL  eTxyL = (*stress_energy_state) ? (eTxy[index]) : 0.0;
    CCTK_REAL  eTxzL = (*stress_energy_state) ? (eTxz[index]) : 0.0;
    CCTK_REAL  eTyyL = (*stress_energy_state) ? (eTyy[index]) : 0.0;
    CCTK_REAL  eTyzL = (*stress_energy_state) ? (eTyz[index]) : 0.0;
    CCTK_REAL  eTzzL = (*stress_energy_state) ? (eTzz[index]) : 0.0;
    CCTK_REAL  gt11L = gt11[index];
    CCTK_REAL  gt12L = gt12[index];
    CCTK_REAL  gt13L = gt13[index];
    CCTK_REAL  gt22L = gt22[index];
    CCTK_REAL  gt23L = gt23[index];
    CCTK_REAL  gt33L = gt33[index];
    CCTK_REAL  phiL = phi[index];
    CCTK_REAL  rL = r[index];
    CCTK_REAL  trKL = trK[index];
    CCTK_REAL  Xt1L = Xt1[index];
    CCTK_REAL  Xt2L = Xt2[index];
    CCTK_REAL  Xt3L = Xt3[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDdissipationNth1A = PDdissipationNth1(A, i, j, k);
    CCTK_REAL const PDdissipationNth2A = PDdissipationNth2(A, i, j, k);
    CCTK_REAL const PDdissipationNth3A = PDdissipationNth3(A, i, j, k);
    CCTK_REAL const PDstandardNth1alpha = PDstandardNth1(alpha, i, j, k);
    CCTK_REAL const PDstandardNth2alpha = PDstandardNth2(alpha, i, j, k);
    CCTK_REAL const PDstandardNth3alpha = PDstandardNth3(alpha, i, j, k);
    CCTK_REAL const PDstandardNth11alpha = PDstandardNth11(alpha, i, j, k);
    CCTK_REAL const PDstandardNth22alpha = PDstandardNth22(alpha, i, j, k);
    CCTK_REAL const PDstandardNth33alpha = PDstandardNth33(alpha, i, j, k);
    CCTK_REAL const PDstandardNth12alpha = PDstandardNth12(alpha, i, j, k);
    CCTK_REAL const PDstandardNth13alpha = PDstandardNth13(alpha, i, j, k);
    CCTK_REAL const PDstandardNth23alpha = PDstandardNth23(alpha, i, j, k);
    CCTK_REAL const PDdissipationNth1alpha = PDdissipationNth1(alpha, i, j, k);
    CCTK_REAL const PDdissipationNth2alpha = PDdissipationNth2(alpha, i, j, k);
    CCTK_REAL const PDdissipationNth3alpha = PDdissipationNth3(alpha, i, j, k);
    CCTK_REAL const PDdissipationNth1B1 = PDdissipationNth1(B1, i, j, k);
    CCTK_REAL const PDdissipationNth2B1 = PDdissipationNth2(B1, i, j, k);
    CCTK_REAL const PDdissipationNth3B1 = PDdissipationNth3(B1, i, j, k);
    CCTK_REAL const PDdissipationNth1B2 = PDdissipationNth1(B2, i, j, k);
    CCTK_REAL const PDdissipationNth2B2 = PDdissipationNth2(B2, i, j, k);
    CCTK_REAL const PDdissipationNth3B2 = PDdissipationNth3(B2, i, j, k);
    CCTK_REAL const PDdissipationNth1B3 = PDdissipationNth1(B3, i, j, k);
    CCTK_REAL const PDdissipationNth2B3 = PDdissipationNth2(B3, i, j, k);
    CCTK_REAL const PDdissipationNth3B3 = PDdissipationNth3(B3, i, j, k);
    CCTK_REAL const PDstandardNth1beta1 = PDstandardNth1(beta1, i, j, k);
    CCTK_REAL const PDstandardNth2beta1 = PDstandardNth2(beta1, i, j, k);
    CCTK_REAL const PDstandardNth3beta1 = PDstandardNth3(beta1, i, j, k);
    CCTK_REAL const PDstandardNth11beta1 = PDstandardNth11(beta1, i, j, k);
    CCTK_REAL const PDstandardNth22beta1 = PDstandardNth22(beta1, i, j, k);
    CCTK_REAL const PDstandardNth33beta1 = PDstandardNth33(beta1, i, j, k);
    CCTK_REAL const PDstandardNth12beta1 = PDstandardNth12(beta1, i, j, k);
    CCTK_REAL const PDstandardNth13beta1 = PDstandardNth13(beta1, i, j, k);
    CCTK_REAL const PDstandardNth23beta1 = PDstandardNth23(beta1, i, j, k);
    CCTK_REAL const PDdissipationNth1beta1 = PDdissipationNth1(beta1, i, j, k);
    CCTK_REAL const PDdissipationNth2beta1 = PDdissipationNth2(beta1, i, j, k);
    CCTK_REAL const PDdissipationNth3beta1 = PDdissipationNth3(beta1, i, j, k);
    CCTK_REAL const PDstandardNth1beta2 = PDstandardNth1(beta2, i, j, k);
    CCTK_REAL const PDstandardNth2beta2 = PDstandardNth2(beta2, i, j, k);
    CCTK_REAL const PDstandardNth3beta2 = PDstandardNth3(beta2, i, j, k);
    CCTK_REAL const PDstandardNth11beta2 = PDstandardNth11(beta2, i, j, k);
    CCTK_REAL const PDstandardNth22beta2 = PDstandardNth22(beta2, i, j, k);
    CCTK_REAL const PDstandardNth33beta2 = PDstandardNth33(beta2, i, j, k);
    CCTK_REAL const PDstandardNth12beta2 = PDstandardNth12(beta2, i, j, k);
    CCTK_REAL const PDstandardNth13beta2 = PDstandardNth13(beta2, i, j, k);
    CCTK_REAL const PDstandardNth23beta2 = PDstandardNth23(beta2, i, j, k);
    CCTK_REAL const PDdissipationNth1beta2 = PDdissipationNth1(beta2, i, j, k);
    CCTK_REAL const PDdissipationNth2beta2 = PDdissipationNth2(beta2, i, j, k);
    CCTK_REAL const PDdissipationNth3beta2 = PDdissipationNth3(beta2, i, j, k);
    CCTK_REAL const PDstandardNth1beta3 = PDstandardNth1(beta3, i, j, k);
    CCTK_REAL const PDstandardNth2beta3 = PDstandardNth2(beta3, i, j, k);
    CCTK_REAL const PDstandardNth3beta3 = PDstandardNth3(beta3, i, j, k);
    CCTK_REAL const PDstandardNth11beta3 = PDstandardNth11(beta3, i, j, k);
    CCTK_REAL const PDstandardNth22beta3 = PDstandardNth22(beta3, i, j, k);
    CCTK_REAL const PDstandardNth33beta3 = PDstandardNth33(beta3, i, j, k);
    CCTK_REAL const PDstandardNth12beta3 = PDstandardNth12(beta3, i, j, k);
    CCTK_REAL const PDstandardNth13beta3 = PDstandardNth13(beta3, i, j, k);
    CCTK_REAL const PDstandardNth23beta3 = PDstandardNth23(beta3, i, j, k);
    CCTK_REAL const PDdissipationNth1beta3 = PDdissipationNth1(beta3, i, j, k);
    CCTK_REAL const PDdissipationNth2beta3 = PDdissipationNth2(beta3, i, j, k);
    CCTK_REAL const PDdissipationNth3beta3 = PDdissipationNth3(beta3, i, j, k);
    CCTK_REAL const PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    CCTK_REAL const PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    CCTK_REAL const PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    CCTK_REAL const PDdissipationNth1gt11 = PDdissipationNth1(gt11, i, j, k);
    CCTK_REAL const PDdissipationNth2gt11 = PDdissipationNth2(gt11, i, j, k);
    CCTK_REAL const PDdissipationNth3gt11 = PDdissipationNth3(gt11, i, j, k);
    CCTK_REAL const PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    CCTK_REAL const PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    CCTK_REAL const PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    CCTK_REAL const PDdissipationNth1gt12 = PDdissipationNth1(gt12, i, j, k);
    CCTK_REAL const PDdissipationNth2gt12 = PDdissipationNth2(gt12, i, j, k);
    CCTK_REAL const PDdissipationNth3gt12 = PDdissipationNth3(gt12, i, j, k);
    CCTK_REAL const PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    CCTK_REAL const PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    CCTK_REAL const PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    CCTK_REAL const PDdissipationNth1gt13 = PDdissipationNth1(gt13, i, j, k);
    CCTK_REAL const PDdissipationNth2gt13 = PDdissipationNth2(gt13, i, j, k);
    CCTK_REAL const PDdissipationNth3gt13 = PDdissipationNth3(gt13, i, j, k);
    CCTK_REAL const PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    CCTK_REAL const PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    CCTK_REAL const PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    CCTK_REAL const PDdissipationNth1gt22 = PDdissipationNth1(gt22, i, j, k);
    CCTK_REAL const PDdissipationNth2gt22 = PDdissipationNth2(gt22, i, j, k);
    CCTK_REAL const PDdissipationNth3gt22 = PDdissipationNth3(gt22, i, j, k);
    CCTK_REAL const PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    CCTK_REAL const PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    CCTK_REAL const PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    CCTK_REAL const PDdissipationNth1gt23 = PDdissipationNth1(gt23, i, j, k);
    CCTK_REAL const PDdissipationNth2gt23 = PDdissipationNth2(gt23, i, j, k);
    CCTK_REAL const PDdissipationNth3gt23 = PDdissipationNth3(gt23, i, j, k);
    CCTK_REAL const PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    CCTK_REAL const PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    CCTK_REAL const PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    CCTK_REAL const PDdissipationNth1gt33 = PDdissipationNth1(gt33, i, j, k);
    CCTK_REAL const PDdissipationNth2gt33 = PDdissipationNth2(gt33, i, j, k);
    CCTK_REAL const PDdissipationNth3gt33 = PDdissipationNth3(gt33, i, j, k);
    CCTK_REAL const PDstandardNth1phi = PDstandardNth1(phi, i, j, k);
    CCTK_REAL const PDstandardNth2phi = PDstandardNth2(phi, i, j, k);
    CCTK_REAL const PDstandardNth3phi = PDstandardNth3(phi, i, j, k);
    CCTK_REAL const PDdissipationNth1phi = PDdissipationNth1(phi, i, j, k);
    CCTK_REAL const PDdissipationNth2phi = PDdissipationNth2(phi, i, j, k);
    CCTK_REAL const PDdissipationNth3phi = PDdissipationNth3(phi, i, j, k);
    CCTK_REAL const PDstandardNth1trK = PDstandardNth1(trK, i, j, k);
    CCTK_REAL const PDstandardNth2trK = PDstandardNth2(trK, i, j, k);
    CCTK_REAL const PDstandardNth3trK = PDstandardNth3(trK, i, j, k);
    CCTK_REAL const PDdissipationNth1trK = PDdissipationNth1(trK, i, j, k);
    CCTK_REAL const PDdissipationNth2trK = PDdissipationNth2(trK, i, j, k);
    CCTK_REAL const PDdissipationNth3trK = PDdissipationNth3(trK, i, j, k);
    CCTK_REAL const PDdissipationNth1Xt1 = PDdissipationNth1(Xt1, i, j, k);
    CCTK_REAL const PDdissipationNth2Xt1 = PDdissipationNth2(Xt1, i, j, k);
    CCTK_REAL const PDdissipationNth3Xt1 = PDdissipationNth3(Xt1, i, j, k);
    CCTK_REAL const PDdissipationNth1Xt2 = PDdissipationNth1(Xt2, i, j, k);
    CCTK_REAL const PDdissipationNth2Xt2 = PDdissipationNth2(Xt2, i, j, k);
    CCTK_REAL const PDdissipationNth3Xt2 = PDdissipationNth3(Xt2, i, j, k);
    CCTK_REAL const PDdissipationNth1Xt3 = PDdissipationNth1(Xt3, i, j, k);
    CCTK_REAL const PDdissipationNth2Xt3 = PDdissipationNth2(Xt3, i, j, k);
    CCTK_REAL const PDdissipationNth3Xt3 = PDdissipationNth3(Xt3, i, j, k);
    
    /* Calculate temporaries and grid functions */
    int dir1 = Sign(beta1L);
    
    int dir2 = Sign(beta2L);
    
    int dir3 = Sign(beta3L);
    
    CCTK_REAL epsdiss1 = EpsDiss;
    
    CCTK_REAL epsdiss2 = EpsDiss;
    
    CCTK_REAL epsdiss3 = EpsDiss;
    
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu12 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu13 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu23 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL Gtl111 = khalf*PDstandardNth1gt11;
    
    CCTK_REAL Gtl112 = khalf*PDstandardNth2gt11;
    
    CCTK_REAL Gtl113 = khalf*PDstandardNth3gt11;
    
    CCTK_REAL Gtl122 = -(khalf*PDstandardNth1gt22) + PDstandardNth2gt12;
    
    CCTK_REAL Gtl123 = khalf*(-PDstandardNth1gt23 + PDstandardNth2gt13 + 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl133 = -(khalf*PDstandardNth1gt33) + PDstandardNth3gt13;
    
    CCTK_REAL Gtl211 = PDstandardNth1gt12 - khalf*PDstandardNth2gt11;
    
    CCTK_REAL Gtl212 = khalf*PDstandardNth1gt22;
    
    CCTK_REAL Gtl213 = khalf*(PDstandardNth1gt23 - PDstandardNth2gt13 + 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl222 = khalf*PDstandardNth2gt22;
    
    CCTK_REAL Gtl223 = khalf*PDstandardNth3gt22;
    
    CCTK_REAL Gtl233 = -(khalf*PDstandardNth2gt33) + PDstandardNth3gt23;
    
    CCTK_REAL Gtl311 = PDstandardNth1gt13 - khalf*PDstandardNth3gt11;
    
    CCTK_REAL Gtl312 = khalf*(PDstandardNth1gt23 + PDstandardNth2gt13 - 
      PDstandardNth3gt12);
    
    CCTK_REAL Gtl313 = khalf*PDstandardNth1gt33;
    
    CCTK_REAL Gtl322 = PDstandardNth2gt23 - khalf*PDstandardNth3gt22;
    
    CCTK_REAL Gtl323 = khalf*PDstandardNth2gt33;
    
    CCTK_REAL Gtl333 = khalf*PDstandardNth3gt33;
    
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
    
    CCTK_REAL fac1 = IfThen(conformalMethod,-(khalf*INV(phiL)),1);
    
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
    
    CCTK_REAL e4phi = IfThen(conformalMethod,pow(phiL,-2),exp(4*phiL));
    
    CCTK_REAL em4phi = INV(e4phi);
    
    CCTK_REAL rho = pow(alphaL,-2)*(eTttL - 2*(beta2L*eTtyL + 
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
    
    CCTK_REAL phirhsL = PDupwindNth1(phi, i, j, k)*beta1L + 
      PDupwindNth2(phi, i, j, k)*beta2L + PDupwindNth3(phi, i, j, 
      k)*beta3L + epsdiss1*PDdissipationNth1phi + 
      epsdiss2*PDdissipationNth2phi + epsdiss3*PDdissipationNth3phi + 
      (PDstandardNth1beta1 + PDstandardNth2beta2 + 
      PDstandardNth3beta3)*IfThen(conformalMethod,-(kthird*phiL),0.16666666666666666) 
      + alphaL*trKL*IfThen(conformalMethod,kthird*phiL,-0.16666666666666666);
    
    CCTK_REAL gt11rhsL = -2*alphaL*At11L + PDupwindNth1(gt11, i, j, 
      k)*beta1L + PDupwindNth2(gt11, i, j, k)*beta2L + PDupwindNth3(gt11, 
      i, j, k)*beta3L + epsdiss1*PDdissipationNth1gt11 + 
      epsdiss2*PDdissipationNth2gt11 + epsdiss3*PDdissipationNth3gt11 + 
      2*(gt12L*PDstandardNth1beta2 + gt13L*PDstandardNth1beta3) + 
      gt11L*(kfourthird*PDstandardNth1beta1 - ktwothird*PDstandardNth2beta2 - 
      ktwothird*PDstandardNth3beta3);
    
    CCTK_REAL gt12rhsL = -2*alphaL*At12L + PDupwindNth1(gt12, i, j, 
      k)*beta1L + PDupwindNth2(gt12, i, j, k)*beta2L + PDupwindNth3(gt12, 
      i, j, k)*beta3L + epsdiss1*PDdissipationNth1gt12 + 
      epsdiss2*PDdissipationNth2gt12 + epsdiss3*PDdissipationNth3gt12 + 
      gt22L*PDstandardNth1beta2 + gt23L*PDstandardNth1beta3 + 
      gt11L*PDstandardNth2beta1 + gt13L*PDstandardNth2beta3 + 
      gt12L*(kthird*(PDstandardNth1beta1 + PDstandardNth2beta2) - 
      ktwothird*PDstandardNth3beta3);
    
    CCTK_REAL gt13rhsL = -2*alphaL*At13L + PDupwindNth1(gt13, i, j, 
      k)*beta1L + PDupwindNth2(gt13, i, j, k)*beta2L + PDupwindNth3(gt13, 
      i, j, k)*beta3L + epsdiss1*PDdissipationNth1gt13 + 
      epsdiss2*PDdissipationNth2gt13 + epsdiss3*PDdissipationNth3gt13 + 
      gt23L*PDstandardNth1beta2 + gt33L*PDstandardNth1beta3 + 
      gt11L*PDstandardNth3beta1 + gt12L*PDstandardNth3beta2 + 
      gt13L*(-(ktwothird*PDstandardNth2beta2) + kthird*(PDstandardNth1beta1 + 
      PDstandardNth3beta3));
    
    CCTK_REAL gt22rhsL = -2*alphaL*At22L + PDupwindNth1(gt22, i, j, 
      k)*beta1L + PDupwindNth2(gt22, i, j, k)*beta2L + PDupwindNth3(gt22, 
      i, j, k)*beta3L + epsdiss1*PDdissipationNth1gt22 + 
      epsdiss2*PDdissipationNth2gt22 + epsdiss3*PDdissipationNth3gt22 + 
      2*(gt12L*PDstandardNth2beta1 + gt23L*PDstandardNth2beta3) + 
      gt22L*(-(ktwothird*PDstandardNth1beta1) + 
      kfourthird*PDstandardNth2beta2 - ktwothird*PDstandardNth3beta3);
    
    CCTK_REAL gt23rhsL = -2*alphaL*At23L + PDupwindNth1(gt23, i, j, 
      k)*beta1L + PDupwindNth2(gt23, i, j, k)*beta2L + PDupwindNth3(gt23, 
      i, j, k)*beta3L + epsdiss1*PDdissipationNth1gt23 + 
      epsdiss2*PDdissipationNth2gt23 + epsdiss3*PDdissipationNth3gt23 + 
      gt13L*PDstandardNth2beta1 + gt33L*PDstandardNth2beta3 + 
      gt12L*PDstandardNth3beta1 + gt22L*PDstandardNth3beta2 + 
      gt23L*(-(ktwothird*PDstandardNth1beta1) + kthird*(PDstandardNth2beta2 + 
      PDstandardNth3beta3));
    
    CCTK_REAL gt33rhsL = -2*alphaL*At33L + PDupwindNth1(gt33, i, j, 
      k)*beta1L + PDupwindNth2(gt33, i, j, k)*beta2L + PDupwindNth3(gt33, 
      i, j, k)*beta3L + epsdiss1*PDdissipationNth1gt33 + 
      epsdiss2*PDdissipationNth2gt33 + epsdiss3*PDdissipationNth3gt33 - 
      gt33L*ktwothird*PDstandardNth1beta1 - 
      gt33L*ktwothird*PDstandardNth2beta2 + 2*gt13L*PDstandardNth3beta1 + 
      2*gt23L*PDstandardNth3beta2 + gt33L*kfourthird*PDstandardNth3beta3;
    
    CCTK_REAL dotXt1 = kthird*(7*(gtu12*PDstandardNth12beta1 + 
      gtu13*PDstandardNth13beta1) + gtu11*(4*PDstandardNth11beta1 + 
      PDstandardNth12beta2 + PDstandardNth13beta3) + 
      gtu12*(PDstandardNth22beta2 + PDstandardNth23beta3) + 
      3*(PDupwindNth1(Xt1, i, j, k)*beta1L + PDupwindNth2(Xt1, i, j, 
      k)*beta2L + PDupwindNth3(Xt1, i, j, k)*beta3L + 
      epsdiss1*PDdissipationNth1Xt1 + epsdiss2*PDdissipationNth2Xt1 + 
      epsdiss3*PDdissipationNth3Xt1 + gtu22*PDstandardNth22beta1 + 
      gtu33*PDstandardNth33beta1) + gtu13*(PDstandardNth23beta2 + 
      PDstandardNth33beta3) - 6*(Atu11*PDstandardNth1alpha + 
      Atu12*PDstandardNth2alpha + Atu13*PDstandardNth3alpha) + 
      6*(gtu23*PDstandardNth23beta1 + alphaL*(6*(Atu11*cdphi1 + Atu12*cdphi2 
      + Atu13*cdphi3) + Atu11*Gt111 + Atu22*Gt122 + 2*(Atu12*Gt112 + 
      Atu13*Gt113 + Atu23*Gt123) + Atu33*Gt133 - 
      ktwothird*(gtu11*PDstandardNth1trK + gtu12*PDstandardNth2trK + 
      gtu13*PDstandardNth3trK))) - 
      150.7964473723100754462068823974161384415*alphaL*(gtu11*S1 + gtu12*S2 + 
      gtu13*S3) + (-3*PDstandardNth1beta1 + 2*(PDstandardNth1beta1 + 
      PDstandardNth2beta2 + PDstandardNth3beta3))*Xtn1 - 
      3*(PDstandardNth2beta1*Xtn2 + PDstandardNth3beta1*Xtn3));
    
    CCTK_REAL dotXt2 = kthird*(gtu12*(PDstandardNth11beta1 + 
      7*PDstandardNth12beta2 + PDstandardNth13beta3) + 
      gtu22*(PDstandardNth12beta1 + 4*PDstandardNth22beta2 + 
      PDstandardNth23beta3) + 3*(PDupwindNth1(Xt2, i, j, k)*beta1L + 
      PDupwindNth2(Xt2, i, j, k)*beta2L + PDupwindNth3(Xt2, i, j, 
      k)*beta3L + epsdiss1*PDdissipationNth1Xt2 + 
      epsdiss2*PDdissipationNth2Xt2 + epsdiss3*PDdissipationNth3Xt2 + 
      gtu11*PDstandardNth11beta2 + gtu33*PDstandardNth33beta2) + 
      gtu23*(PDstandardNth13beta1 + 7*PDstandardNth23beta2 + 
      PDstandardNth33beta3) - 6*(Atu12*PDstandardNth1alpha + 
      Atu22*PDstandardNth2alpha + Atu23*PDstandardNth3alpha) + 
      6*(gtu13*PDstandardNth13beta2 + alphaL*(6*(Atu12*cdphi1 + Atu22*cdphi2 
      + Atu23*cdphi3) + Atu11*Gt211 + Atu22*Gt222 + 2*(Atu12*Gt212 + 
      Atu13*Gt213 + Atu23*Gt223) + Atu33*Gt233 - 
      ktwothird*(gtu12*PDstandardNth1trK + gtu22*PDstandardNth2trK + 
      gtu23*PDstandardNth3trK))) - 
      150.7964473723100754462068823974161384415*alphaL*(gtu12*S1 + gtu22*S2 + 
      gtu23*S3) + 2*(PDstandardNth1beta1 + PDstandardNth2beta2 + 
      PDstandardNth3beta3)*Xtn2 - 3*(PDstandardNth1beta2*Xtn1 + 
      PDstandardNth2beta2*Xtn2 + PDstandardNth3beta2*Xtn3));
    
    CCTK_REAL dotXt3 = kthird*(gtu13*(PDstandardNth11beta1 + 
      PDstandardNth12beta2 + 7*PDstandardNth13beta3) + 3*(PDupwindNth1(Xt3, 
      i, j, k)*beta1L + PDupwindNth2(Xt3, i, j, k)*beta2L + 
      PDupwindNth3(Xt3, i, j, k)*beta3L + epsdiss1*PDdissipationNth1Xt3 + 
      epsdiss2*PDdissipationNth2Xt3 + epsdiss3*PDdissipationNth3Xt3 + 
      gtu11*PDstandardNth11beta3 + gtu22*PDstandardNth22beta3) + 
      gtu23*(PDstandardNth12beta1 + PDstandardNth22beta2 + 
      7*PDstandardNth23beta3) + gtu33*(PDstandardNth13beta1 + 
      PDstandardNth23beta2 + 4*PDstandardNth33beta3) - 
      6*(Atu13*PDstandardNth1alpha + Atu23*PDstandardNth2alpha + 
      Atu33*PDstandardNth3alpha) + 6*(gtu12*PDstandardNth12beta3 + 
      alphaL*(6*(Atu13*cdphi1 + Atu23*cdphi2 + Atu33*cdphi3) + Atu11*Gt311 + 
      Atu22*Gt322 + 2*(Atu12*Gt312 + Atu13*Gt313 + Atu23*Gt323) + Atu33*Gt333 
      - ktwothird*(gtu13*PDstandardNth1trK + gtu23*PDstandardNth2trK + 
      gtu33*PDstandardNth3trK))) - 
      150.7964473723100754462068823974161384415*alphaL*(gtu13*S1 + gtu23*S2 + 
      gtu33*S3) + 2*(PDstandardNth1beta1 + PDstandardNth2beta2 + 
      PDstandardNth3beta3)*Xtn3 - 3*(PDstandardNth1beta3*Xtn1 + 
      PDstandardNth2beta3*Xtn2 + PDstandardNth3beta3*Xtn3));
    
    CCTK_REAL Xt1rhsL = dotXt1;
    
    CCTK_REAL Xt2rhsL = dotXt2;
    
    CCTK_REAL Xt3rhsL = dotXt3;
    
    CCTK_REAL dottrK = PDupwindNth1(trK, i, j, k)*beta1L + 
      PDupwindNth2(trK, i, j, k)*beta2L + PDupwindNth3(trK, i, j, 
      k)*beta3L + epsdiss1*PDdissipationNth1trK + 
      epsdiss2*PDdissipationNth2trK + epsdiss3*PDdissipationNth3trK - 
      em4phi*(gtu11*PDstandardNth11alpha + gtu22*PDstandardNth22alpha + 
      gtu33*(PDstandardNth33alpha + 2*cdphi3*PDstandardNth3alpha) + 
      2*(gtu12*PDstandardNth12alpha + gtu13*(PDstandardNth13alpha + 
      cdphi1*PDstandardNth3alpha) + gtu23*(PDstandardNth23alpha + 
      cdphi2*PDstandardNth3alpha)) + PDstandardNth1alpha*(2*(cdphi1*gtu11 + 
      cdphi2*gtu12 + cdphi3*gtu13) - Xtn1) + 
      PDstandardNth2alpha*(2*(cdphi1*gtu12 + cdphi2*gtu22 + cdphi3*gtu23) - 
      Xtn2) - PDstandardNth3alpha*Xtn3) + alphaL*(2*(Atm12*Atm21 + 
      Atm13*Atm31 + Atm23*Atm32) + 
      12.56637061435917295385057353311801153679*(rho + trS) + SQR(Atm11) + 
      SQR(Atm22) + SQR(Atm33) + kthird*SQR(trKL));
    
    CCTK_REAL trKrhsL = dottrK;
    
    CCTK_REAL alpharhsL = (PDupwindNth1(alpha, i, j, k)*beta1L + 
      PDupwindNth2(alpha, i, j, k)*beta2L + PDupwindNth3(alpha, i, j, 
      k)*beta3L)*LapseAdvectionCoeff + epsdiss1*PDdissipationNth1alpha + 
      epsdiss2*PDdissipationNth2alpha + epsdiss3*PDdissipationNth3alpha - 
      harmonicF*(LapseACoeff*(AL - trKL) + trKL)*pow(alphaL,harmonicN);
    
    CCTK_REAL ArhsL = (-(AL*AlphaDriver) + dottrK)*LapseACoeff + 
      (PDupwindNth1(A, i, j, k)*beta1L + PDupwindNth2(A, i, j, k)*beta2L 
      + PDupwindNth3(A, i, j, k)*beta3L)*LapseAdvectionCoeff + 
      epsdiss1*PDdissipationNth1A + epsdiss2*PDdissipationNth2A + 
      epsdiss3*PDdissipationNth3A;
    
    CCTK_REAL eta = fmin(1,SpatialBetaDriverRadius*INV(rL));
    
    CCTK_REAL theta = fmin(1,exp(1 - 
      rL*INV(SpatialShiftGammaCoeffRadius)));
    
    CCTK_REAL beta1rhsL = epsdiss1*PDdissipationNth1beta1 + 
      epsdiss2*PDdissipationNth2beta1 + epsdiss3*PDdissipationNth3beta1 + 
      (PDupwindNth1(beta1, i, j, k)*beta1L + PDupwindNth2(beta1, i, j, 
      k)*beta2L + PDupwindNth3(beta1, i, j, k)*beta3L)*ShiftAdvectionCoeff 
      + ShiftGammaCoeff*theta*(beta1L*BetaDriver*eta*(-1 + ShiftBCoeff) + 
      ShiftBCoeff*(B1L - Xt1L) + Xt1L);
    
    CCTK_REAL beta2rhsL = epsdiss1*PDdissipationNth1beta2 + 
      epsdiss2*PDdissipationNth2beta2 + epsdiss3*PDdissipationNth3beta2 + 
      (PDupwindNth1(beta2, i, j, k)*beta1L + PDupwindNth2(beta2, i, j, 
      k)*beta2L + PDupwindNth3(beta2, i, j, k)*beta3L)*ShiftAdvectionCoeff 
      + ShiftGammaCoeff*theta*(beta2L*BetaDriver*eta*(-1 + ShiftBCoeff) + 
      ShiftBCoeff*(B2L - Xt2L) + Xt2L);
    
    CCTK_REAL beta3rhsL = epsdiss1*PDdissipationNth1beta3 + 
      epsdiss2*PDdissipationNth2beta3 + epsdiss3*PDdissipationNth3beta3 + 
      (PDupwindNth1(beta3, i, j, k)*beta1L + PDupwindNth2(beta3, i, j, 
      k)*beta2L + PDupwindNth3(beta3, i, j, k)*beta3L)*ShiftAdvectionCoeff 
      + ShiftGammaCoeff*theta*(beta3L*BetaDriver*eta*(-1 + ShiftBCoeff) + 
      ShiftBCoeff*(B3L - Xt3L) + Xt3L);
    
    CCTK_REAL B1rhsL = epsdiss1*PDdissipationNth1B1 + 
      epsdiss2*PDdissipationNth2B1 + epsdiss3*PDdissipationNth3B1 + 
      ((PDupwindNth1(B1, i, j, k) - PDupwindNth1(Xt1, i, j, k))*beta1L + 
      (PDupwindNth2(B1, i, j, k) - PDupwindNth2(Xt1, i, j, k))*beta2L + 
      (PDupwindNth3(B1, i, j, k) - PDupwindNth3(Xt1, i, j, 
      k))*beta3L)*ShiftAdvectionCoeff + (dotXt1 - 
      B1L*BetaDriver*eta)*ShiftBCoeff;
    
    CCTK_REAL B2rhsL = epsdiss1*PDdissipationNth1B2 + 
      epsdiss2*PDdissipationNth2B2 + epsdiss3*PDdissipationNth3B2 + 
      ((PDupwindNth1(B2, i, j, k) - PDupwindNth1(Xt2, i, j, k))*beta1L + 
      (PDupwindNth2(B2, i, j, k) - PDupwindNth2(Xt2, i, j, k))*beta2L + 
      (PDupwindNth3(B2, i, j, k) - PDupwindNth3(Xt2, i, j, 
      k))*beta3L)*ShiftAdvectionCoeff + (dotXt2 - 
      B2L*BetaDriver*eta)*ShiftBCoeff;
    
    CCTK_REAL B3rhsL = epsdiss1*PDdissipationNth1B3 + 
      epsdiss2*PDdissipationNth2B3 + epsdiss3*PDdissipationNth3B3 + 
      ((PDupwindNth1(B3, i, j, k) - PDupwindNth1(Xt3, i, j, k))*beta1L + 
      (PDupwindNth2(B3, i, j, k) - PDupwindNth2(Xt3, i, j, k))*beta2L + 
      (PDupwindNth3(B3, i, j, k) - PDupwindNth3(Xt3, i, j, 
      k))*beta3L)*ShiftAdvectionCoeff + (dotXt3 - 
      B3L*BetaDriver*eta)*ShiftBCoeff;
    
    
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
  LC_ENDLOOP3 (ML_BSSN_UPW_RHS1);
}

void ML_BSSN_UPW_RHS1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_UPW_RHS1_Body);
}
