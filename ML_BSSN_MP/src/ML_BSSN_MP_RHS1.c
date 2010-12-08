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

void ML_BSSN_MP_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_trace_curvrhs.");
  return;
}

void ML_BSSN_MP_RHS1_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_RHS1_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_RHS1_calc_every != ML_BSSN_MP_RHS1_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"Coordinates::jacobian","Coordinates::jacobian2","grid::coordinates","Grid::coordinates","ML_BSSN_MP::ML_curv","ML_BSSN_MP::ML_dtlapse","ML_BSSN_MP::ML_dtlapserhs","ML_BSSN_MP::ML_dtshift","ML_BSSN_MP::ML_dtshiftrhs","ML_BSSN_MP::ML_Gamma","ML_BSSN_MP::ML_Gammarhs","ML_BSSN_MP::ML_lapse","ML_BSSN_MP::ML_lapserhs","ML_BSSN_MP::ML_log_confac","ML_BSSN_MP::ML_log_confacrhs","ML_BSSN_MP::ML_metric","ML_BSSN_MP::ML_metricrhs","ML_BSSN_MP::ML_shift","ML_BSSN_MP::ML_shiftrhs","ML_BSSN_MP::ML_trace_curv","ML_BSSN_MP::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_RHS1", 21, groups);
  
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
  LC_LOOP3 (ML_BSSN_MP_RHS1,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL PDdissipationNth1A = INITVALUE;
    // CCTK_REAL PDdissipationNth2A = INITVALUE;
    // CCTK_REAL PDdissipationNth3A = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1A = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1A = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2A = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2A = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3A = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3A = INITVALUE;
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
    // CCTK_REAL PDupwindNthAnti1alpha = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1alpha = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2alpha = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2alpha = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3alpha = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3alpha = INITVALUE;
    // CCTK_REAL PDdissipationNth1B1 = INITVALUE;
    // CCTK_REAL PDdissipationNth2B1 = INITVALUE;
    // CCTK_REAL PDdissipationNth3B1 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1B1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1B1 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2B1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2B1 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3B1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3B1 = INITVALUE;
    // CCTK_REAL PDdissipationNth1B2 = INITVALUE;
    // CCTK_REAL PDdissipationNth2B2 = INITVALUE;
    // CCTK_REAL PDdissipationNth3B2 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1B2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1B2 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2B2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2B2 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3B2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3B2 = INITVALUE;
    // CCTK_REAL PDdissipationNth1B3 = INITVALUE;
    // CCTK_REAL PDdissipationNth2B3 = INITVALUE;
    // CCTK_REAL PDdissipationNth3B3 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1B3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1B3 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2B3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2B3 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3B3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3B3 = INITVALUE;
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
    // CCTK_REAL PDupwindNthAnti1beta1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1beta1 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2beta1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2beta1 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3beta1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3beta1 = INITVALUE;
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
    // CCTK_REAL PDupwindNthAnti1beta2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1beta2 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2beta2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2beta2 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3beta2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3beta2 = INITVALUE;
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
    // CCTK_REAL PDupwindNthAnti1beta3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1beta3 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2beta3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2beta3 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3beta3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt11 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt11 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt11 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt11 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1gt11 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1gt11 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2gt11 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2gt11 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3gt11 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt12 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt12 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt12 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt12 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1gt12 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1gt12 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2gt12 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2gt12 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3gt12 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt13 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt13 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt13 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt13 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1gt13 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1gt13 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2gt13 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2gt13 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3gt13 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt22 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt22 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt22 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt22 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1gt22 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1gt22 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2gt22 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2gt22 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3gt22 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt23 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt23 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt23 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt23 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1gt23 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1gt23 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2gt23 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2gt23 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3gt23 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt33 = INITVALUE;
    // CCTK_REAL PDdissipationNth1gt33 = INITVALUE;
    // CCTK_REAL PDdissipationNth2gt33 = INITVALUE;
    // CCTK_REAL PDdissipationNth3gt33 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1gt33 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1gt33 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2gt33 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2gt33 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3gt33 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth1phi = INITVALUE;
    // CCTK_REAL PDstandardNth2phi = INITVALUE;
    // CCTK_REAL PDstandardNth3phi = INITVALUE;
    // CCTK_REAL PDdissipationNth1phi = INITVALUE;
    // CCTK_REAL PDdissipationNth2phi = INITVALUE;
    // CCTK_REAL PDdissipationNth3phi = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1phi = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1phi = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2phi = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2phi = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3phi = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3phi = INITVALUE;
    // CCTK_REAL PDstandardNth1trK = INITVALUE;
    // CCTK_REAL PDstandardNth2trK = INITVALUE;
    // CCTK_REAL PDstandardNth3trK = INITVALUE;
    // CCTK_REAL PDdissipationNth1trK = INITVALUE;
    // CCTK_REAL PDdissipationNth2trK = INITVALUE;
    // CCTK_REAL PDdissipationNth3trK = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1trK = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1trK = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2trK = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2trK = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3trK = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3trK = INITVALUE;
    // CCTK_REAL PDdissipationNth1Xt1 = INITVALUE;
    // CCTK_REAL PDdissipationNth2Xt1 = INITVALUE;
    // CCTK_REAL PDdissipationNth3Xt1 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1Xt1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1Xt1 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2Xt1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2Xt1 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3Xt1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3Xt1 = INITVALUE;
    // CCTK_REAL PDdissipationNth1Xt2 = INITVALUE;
    // CCTK_REAL PDdissipationNth2Xt2 = INITVALUE;
    // CCTK_REAL PDdissipationNth3Xt2 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1Xt2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1Xt2 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2Xt2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2Xt2 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3Xt2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3Xt2 = INITVALUE;
    // CCTK_REAL PDdissipationNth1Xt3 = INITVALUE;
    // CCTK_REAL PDdissipationNth2Xt3 = INITVALUE;
    // CCTK_REAL PDdissipationNth3Xt3 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1Xt3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1Xt3 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2Xt3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2Xt3 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3Xt3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3Xt3 = INITVALUE;
    
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
    CCTK_REAL  dJ111L = dJ111[index];
    CCTK_REAL  dJ112L = dJ112[index];
    CCTK_REAL  dJ113L = dJ113[index];
    CCTK_REAL  dJ122L = dJ122[index];
    CCTK_REAL  dJ123L = dJ123[index];
    CCTK_REAL  dJ133L = dJ133[index];
    CCTK_REAL  dJ211L = dJ211[index];
    CCTK_REAL  dJ212L = dJ212[index];
    CCTK_REAL  dJ213L = dJ213[index];
    CCTK_REAL  dJ222L = dJ222[index];
    CCTK_REAL  dJ223L = dJ223[index];
    CCTK_REAL  dJ233L = dJ233[index];
    CCTK_REAL  dJ311L = dJ311[index];
    CCTK_REAL  dJ312L = dJ312[index];
    CCTK_REAL  dJ313L = dJ313[index];
    CCTK_REAL  dJ322L = dJ322[index];
    CCTK_REAL  dJ323L = dJ323[index];
    CCTK_REAL  dJ333L = dJ333[index];
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
    CCTK_REAL  J11L = J11[index];
    CCTK_REAL  J12L = J12[index];
    CCTK_REAL  J13L = J13[index];
    CCTK_REAL  J21L = J21[index];
    CCTK_REAL  J22L = J22[index];
    CCTK_REAL  J23L = J23[index];
    CCTK_REAL  J31L = J31[index];
    CCTK_REAL  J32L = J32[index];
    CCTK_REAL  J33L = J33[index];
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
    CCTK_REAL const PDupwindNthAnti1A = PDupwindNthAnti1(A, i, j, k);
    CCTK_REAL const PDupwindNthSymm1A = PDupwindNthSymm1(A, i, j, k);
    CCTK_REAL const PDupwindNthAnti2A = PDupwindNthAnti2(A, i, j, k);
    CCTK_REAL const PDupwindNthSymm2A = PDupwindNthSymm2(A, i, j, k);
    CCTK_REAL const PDupwindNthAnti3A = PDupwindNthAnti3(A, i, j, k);
    CCTK_REAL const PDupwindNthSymm3A = PDupwindNthSymm3(A, i, j, k);
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
    CCTK_REAL const PDupwindNthAnti1alpha = PDupwindNthAnti1(alpha, i, j, k);
    CCTK_REAL const PDupwindNthSymm1alpha = PDupwindNthSymm1(alpha, i, j, k);
    CCTK_REAL const PDupwindNthAnti2alpha = PDupwindNthAnti2(alpha, i, j, k);
    CCTK_REAL const PDupwindNthSymm2alpha = PDupwindNthSymm2(alpha, i, j, k);
    CCTK_REAL const PDupwindNthAnti3alpha = PDupwindNthAnti3(alpha, i, j, k);
    CCTK_REAL const PDupwindNthSymm3alpha = PDupwindNthSymm3(alpha, i, j, k);
    CCTK_REAL const PDdissipationNth1B1 = PDdissipationNth1(B1, i, j, k);
    CCTK_REAL const PDdissipationNth2B1 = PDdissipationNth2(B1, i, j, k);
    CCTK_REAL const PDdissipationNth3B1 = PDdissipationNth3(B1, i, j, k);
    CCTK_REAL const PDupwindNthAnti1B1 = PDupwindNthAnti1(B1, i, j, k);
    CCTK_REAL const PDupwindNthSymm1B1 = PDupwindNthSymm1(B1, i, j, k);
    CCTK_REAL const PDupwindNthAnti2B1 = PDupwindNthAnti2(B1, i, j, k);
    CCTK_REAL const PDupwindNthSymm2B1 = PDupwindNthSymm2(B1, i, j, k);
    CCTK_REAL const PDupwindNthAnti3B1 = PDupwindNthAnti3(B1, i, j, k);
    CCTK_REAL const PDupwindNthSymm3B1 = PDupwindNthSymm3(B1, i, j, k);
    CCTK_REAL const PDdissipationNth1B2 = PDdissipationNth1(B2, i, j, k);
    CCTK_REAL const PDdissipationNth2B2 = PDdissipationNth2(B2, i, j, k);
    CCTK_REAL const PDdissipationNth3B2 = PDdissipationNth3(B2, i, j, k);
    CCTK_REAL const PDupwindNthAnti1B2 = PDupwindNthAnti1(B2, i, j, k);
    CCTK_REAL const PDupwindNthSymm1B2 = PDupwindNthSymm1(B2, i, j, k);
    CCTK_REAL const PDupwindNthAnti2B2 = PDupwindNthAnti2(B2, i, j, k);
    CCTK_REAL const PDupwindNthSymm2B2 = PDupwindNthSymm2(B2, i, j, k);
    CCTK_REAL const PDupwindNthAnti3B2 = PDupwindNthAnti3(B2, i, j, k);
    CCTK_REAL const PDupwindNthSymm3B2 = PDupwindNthSymm3(B2, i, j, k);
    CCTK_REAL const PDdissipationNth1B3 = PDdissipationNth1(B3, i, j, k);
    CCTK_REAL const PDdissipationNth2B3 = PDdissipationNth2(B3, i, j, k);
    CCTK_REAL const PDdissipationNth3B3 = PDdissipationNth3(B3, i, j, k);
    CCTK_REAL const PDupwindNthAnti1B3 = PDupwindNthAnti1(B3, i, j, k);
    CCTK_REAL const PDupwindNthSymm1B3 = PDupwindNthSymm1(B3, i, j, k);
    CCTK_REAL const PDupwindNthAnti2B3 = PDupwindNthAnti2(B3, i, j, k);
    CCTK_REAL const PDupwindNthSymm2B3 = PDupwindNthSymm2(B3, i, j, k);
    CCTK_REAL const PDupwindNthAnti3B3 = PDupwindNthAnti3(B3, i, j, k);
    CCTK_REAL const PDupwindNthSymm3B3 = PDupwindNthSymm3(B3, i, j, k);
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
    CCTK_REAL const PDupwindNthAnti1beta1 = PDupwindNthAnti1(beta1, i, j, k);
    CCTK_REAL const PDupwindNthSymm1beta1 = PDupwindNthSymm1(beta1, i, j, k);
    CCTK_REAL const PDupwindNthAnti2beta1 = PDupwindNthAnti2(beta1, i, j, k);
    CCTK_REAL const PDupwindNthSymm2beta1 = PDupwindNthSymm2(beta1, i, j, k);
    CCTK_REAL const PDupwindNthAnti3beta1 = PDupwindNthAnti3(beta1, i, j, k);
    CCTK_REAL const PDupwindNthSymm3beta1 = PDupwindNthSymm3(beta1, i, j, k);
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
    CCTK_REAL const PDupwindNthAnti1beta2 = PDupwindNthAnti1(beta2, i, j, k);
    CCTK_REAL const PDupwindNthSymm1beta2 = PDupwindNthSymm1(beta2, i, j, k);
    CCTK_REAL const PDupwindNthAnti2beta2 = PDupwindNthAnti2(beta2, i, j, k);
    CCTK_REAL const PDupwindNthSymm2beta2 = PDupwindNthSymm2(beta2, i, j, k);
    CCTK_REAL const PDupwindNthAnti3beta2 = PDupwindNthAnti3(beta2, i, j, k);
    CCTK_REAL const PDupwindNthSymm3beta2 = PDupwindNthSymm3(beta2, i, j, k);
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
    CCTK_REAL const PDupwindNthAnti1beta3 = PDupwindNthAnti1(beta3, i, j, k);
    CCTK_REAL const PDupwindNthSymm1beta3 = PDupwindNthSymm1(beta3, i, j, k);
    CCTK_REAL const PDupwindNthAnti2beta3 = PDupwindNthAnti2(beta3, i, j, k);
    CCTK_REAL const PDupwindNthSymm2beta3 = PDupwindNthSymm2(beta3, i, j, k);
    CCTK_REAL const PDupwindNthAnti3beta3 = PDupwindNthAnti3(beta3, i, j, k);
    CCTK_REAL const PDupwindNthSymm3beta3 = PDupwindNthSymm3(beta3, i, j, k);
    CCTK_REAL const PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    CCTK_REAL const PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    CCTK_REAL const PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    CCTK_REAL const PDdissipationNth1gt11 = PDdissipationNth1(gt11, i, j, k);
    CCTK_REAL const PDdissipationNth2gt11 = PDdissipationNth2(gt11, i, j, k);
    CCTK_REAL const PDdissipationNth3gt11 = PDdissipationNth3(gt11, i, j, k);
    CCTK_REAL const PDupwindNthAnti1gt11 = PDupwindNthAnti1(gt11, i, j, k);
    CCTK_REAL const PDupwindNthSymm1gt11 = PDupwindNthSymm1(gt11, i, j, k);
    CCTK_REAL const PDupwindNthAnti2gt11 = PDupwindNthAnti2(gt11, i, j, k);
    CCTK_REAL const PDupwindNthSymm2gt11 = PDupwindNthSymm2(gt11, i, j, k);
    CCTK_REAL const PDupwindNthAnti3gt11 = PDupwindNthAnti3(gt11, i, j, k);
    CCTK_REAL const PDupwindNthSymm3gt11 = PDupwindNthSymm3(gt11, i, j, k);
    CCTK_REAL const PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    CCTK_REAL const PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    CCTK_REAL const PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    CCTK_REAL const PDdissipationNth1gt12 = PDdissipationNth1(gt12, i, j, k);
    CCTK_REAL const PDdissipationNth2gt12 = PDdissipationNth2(gt12, i, j, k);
    CCTK_REAL const PDdissipationNth3gt12 = PDdissipationNth3(gt12, i, j, k);
    CCTK_REAL const PDupwindNthAnti1gt12 = PDupwindNthAnti1(gt12, i, j, k);
    CCTK_REAL const PDupwindNthSymm1gt12 = PDupwindNthSymm1(gt12, i, j, k);
    CCTK_REAL const PDupwindNthAnti2gt12 = PDupwindNthAnti2(gt12, i, j, k);
    CCTK_REAL const PDupwindNthSymm2gt12 = PDupwindNthSymm2(gt12, i, j, k);
    CCTK_REAL const PDupwindNthAnti3gt12 = PDupwindNthAnti3(gt12, i, j, k);
    CCTK_REAL const PDupwindNthSymm3gt12 = PDupwindNthSymm3(gt12, i, j, k);
    CCTK_REAL const PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    CCTK_REAL const PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    CCTK_REAL const PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    CCTK_REAL const PDdissipationNth1gt13 = PDdissipationNth1(gt13, i, j, k);
    CCTK_REAL const PDdissipationNth2gt13 = PDdissipationNth2(gt13, i, j, k);
    CCTK_REAL const PDdissipationNth3gt13 = PDdissipationNth3(gt13, i, j, k);
    CCTK_REAL const PDupwindNthAnti1gt13 = PDupwindNthAnti1(gt13, i, j, k);
    CCTK_REAL const PDupwindNthSymm1gt13 = PDupwindNthSymm1(gt13, i, j, k);
    CCTK_REAL const PDupwindNthAnti2gt13 = PDupwindNthAnti2(gt13, i, j, k);
    CCTK_REAL const PDupwindNthSymm2gt13 = PDupwindNthSymm2(gt13, i, j, k);
    CCTK_REAL const PDupwindNthAnti3gt13 = PDupwindNthAnti3(gt13, i, j, k);
    CCTK_REAL const PDupwindNthSymm3gt13 = PDupwindNthSymm3(gt13, i, j, k);
    CCTK_REAL const PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    CCTK_REAL const PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    CCTK_REAL const PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    CCTK_REAL const PDdissipationNth1gt22 = PDdissipationNth1(gt22, i, j, k);
    CCTK_REAL const PDdissipationNth2gt22 = PDdissipationNth2(gt22, i, j, k);
    CCTK_REAL const PDdissipationNth3gt22 = PDdissipationNth3(gt22, i, j, k);
    CCTK_REAL const PDupwindNthAnti1gt22 = PDupwindNthAnti1(gt22, i, j, k);
    CCTK_REAL const PDupwindNthSymm1gt22 = PDupwindNthSymm1(gt22, i, j, k);
    CCTK_REAL const PDupwindNthAnti2gt22 = PDupwindNthAnti2(gt22, i, j, k);
    CCTK_REAL const PDupwindNthSymm2gt22 = PDupwindNthSymm2(gt22, i, j, k);
    CCTK_REAL const PDupwindNthAnti3gt22 = PDupwindNthAnti3(gt22, i, j, k);
    CCTK_REAL const PDupwindNthSymm3gt22 = PDupwindNthSymm3(gt22, i, j, k);
    CCTK_REAL const PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    CCTK_REAL const PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    CCTK_REAL const PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    CCTK_REAL const PDdissipationNth1gt23 = PDdissipationNth1(gt23, i, j, k);
    CCTK_REAL const PDdissipationNth2gt23 = PDdissipationNth2(gt23, i, j, k);
    CCTK_REAL const PDdissipationNth3gt23 = PDdissipationNth3(gt23, i, j, k);
    CCTK_REAL const PDupwindNthAnti1gt23 = PDupwindNthAnti1(gt23, i, j, k);
    CCTK_REAL const PDupwindNthSymm1gt23 = PDupwindNthSymm1(gt23, i, j, k);
    CCTK_REAL const PDupwindNthAnti2gt23 = PDupwindNthAnti2(gt23, i, j, k);
    CCTK_REAL const PDupwindNthSymm2gt23 = PDupwindNthSymm2(gt23, i, j, k);
    CCTK_REAL const PDupwindNthAnti3gt23 = PDupwindNthAnti3(gt23, i, j, k);
    CCTK_REAL const PDupwindNthSymm3gt23 = PDupwindNthSymm3(gt23, i, j, k);
    CCTK_REAL const PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    CCTK_REAL const PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    CCTK_REAL const PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    CCTK_REAL const PDdissipationNth1gt33 = PDdissipationNth1(gt33, i, j, k);
    CCTK_REAL const PDdissipationNth2gt33 = PDdissipationNth2(gt33, i, j, k);
    CCTK_REAL const PDdissipationNth3gt33 = PDdissipationNth3(gt33, i, j, k);
    CCTK_REAL const PDupwindNthAnti1gt33 = PDupwindNthAnti1(gt33, i, j, k);
    CCTK_REAL const PDupwindNthSymm1gt33 = PDupwindNthSymm1(gt33, i, j, k);
    CCTK_REAL const PDupwindNthAnti2gt33 = PDupwindNthAnti2(gt33, i, j, k);
    CCTK_REAL const PDupwindNthSymm2gt33 = PDupwindNthSymm2(gt33, i, j, k);
    CCTK_REAL const PDupwindNthAnti3gt33 = PDupwindNthAnti3(gt33, i, j, k);
    CCTK_REAL const PDupwindNthSymm3gt33 = PDupwindNthSymm3(gt33, i, j, k);
    CCTK_REAL const PDstandardNth1phi = PDstandardNth1(phi, i, j, k);
    CCTK_REAL const PDstandardNth2phi = PDstandardNth2(phi, i, j, k);
    CCTK_REAL const PDstandardNth3phi = PDstandardNth3(phi, i, j, k);
    CCTK_REAL const PDdissipationNth1phi = PDdissipationNth1(phi, i, j, k);
    CCTK_REAL const PDdissipationNth2phi = PDdissipationNth2(phi, i, j, k);
    CCTK_REAL const PDdissipationNth3phi = PDdissipationNth3(phi, i, j, k);
    CCTK_REAL const PDupwindNthAnti1phi = PDupwindNthAnti1(phi, i, j, k);
    CCTK_REAL const PDupwindNthSymm1phi = PDupwindNthSymm1(phi, i, j, k);
    CCTK_REAL const PDupwindNthAnti2phi = PDupwindNthAnti2(phi, i, j, k);
    CCTK_REAL const PDupwindNthSymm2phi = PDupwindNthSymm2(phi, i, j, k);
    CCTK_REAL const PDupwindNthAnti3phi = PDupwindNthAnti3(phi, i, j, k);
    CCTK_REAL const PDupwindNthSymm3phi = PDupwindNthSymm3(phi, i, j, k);
    CCTK_REAL const PDstandardNth1trK = PDstandardNth1(trK, i, j, k);
    CCTK_REAL const PDstandardNth2trK = PDstandardNth2(trK, i, j, k);
    CCTK_REAL const PDstandardNth3trK = PDstandardNth3(trK, i, j, k);
    CCTK_REAL const PDdissipationNth1trK = PDdissipationNth1(trK, i, j, k);
    CCTK_REAL const PDdissipationNth2trK = PDdissipationNth2(trK, i, j, k);
    CCTK_REAL const PDdissipationNth3trK = PDdissipationNth3(trK, i, j, k);
    CCTK_REAL const PDupwindNthAnti1trK = PDupwindNthAnti1(trK, i, j, k);
    CCTK_REAL const PDupwindNthSymm1trK = PDupwindNthSymm1(trK, i, j, k);
    CCTK_REAL const PDupwindNthAnti2trK = PDupwindNthAnti2(trK, i, j, k);
    CCTK_REAL const PDupwindNthSymm2trK = PDupwindNthSymm2(trK, i, j, k);
    CCTK_REAL const PDupwindNthAnti3trK = PDupwindNthAnti3(trK, i, j, k);
    CCTK_REAL const PDupwindNthSymm3trK = PDupwindNthSymm3(trK, i, j, k);
    CCTK_REAL const PDdissipationNth1Xt1 = PDdissipationNth1(Xt1, i, j, k);
    CCTK_REAL const PDdissipationNth2Xt1 = PDdissipationNth2(Xt1, i, j, k);
    CCTK_REAL const PDdissipationNth3Xt1 = PDdissipationNth3(Xt1, i, j, k);
    CCTK_REAL const PDupwindNthAnti1Xt1 = PDupwindNthAnti1(Xt1, i, j, k);
    CCTK_REAL const PDupwindNthSymm1Xt1 = PDupwindNthSymm1(Xt1, i, j, k);
    CCTK_REAL const PDupwindNthAnti2Xt1 = PDupwindNthAnti2(Xt1, i, j, k);
    CCTK_REAL const PDupwindNthSymm2Xt1 = PDupwindNthSymm2(Xt1, i, j, k);
    CCTK_REAL const PDupwindNthAnti3Xt1 = PDupwindNthAnti3(Xt1, i, j, k);
    CCTK_REAL const PDupwindNthSymm3Xt1 = PDupwindNthSymm3(Xt1, i, j, k);
    CCTK_REAL const PDdissipationNth1Xt2 = PDdissipationNth1(Xt2, i, j, k);
    CCTK_REAL const PDdissipationNth2Xt2 = PDdissipationNth2(Xt2, i, j, k);
    CCTK_REAL const PDdissipationNth3Xt2 = PDdissipationNth3(Xt2, i, j, k);
    CCTK_REAL const PDupwindNthAnti1Xt2 = PDupwindNthAnti1(Xt2, i, j, k);
    CCTK_REAL const PDupwindNthSymm1Xt2 = PDupwindNthSymm1(Xt2, i, j, k);
    CCTK_REAL const PDupwindNthAnti2Xt2 = PDupwindNthAnti2(Xt2, i, j, k);
    CCTK_REAL const PDupwindNthSymm2Xt2 = PDupwindNthSymm2(Xt2, i, j, k);
    CCTK_REAL const PDupwindNthAnti3Xt2 = PDupwindNthAnti3(Xt2, i, j, k);
    CCTK_REAL const PDupwindNthSymm3Xt2 = PDupwindNthSymm3(Xt2, i, j, k);
    CCTK_REAL const PDdissipationNth1Xt3 = PDdissipationNth1(Xt3, i, j, k);
    CCTK_REAL const PDdissipationNth2Xt3 = PDdissipationNth2(Xt3, i, j, k);
    CCTK_REAL const PDdissipationNth3Xt3 = PDdissipationNth3(Xt3, i, j, k);
    CCTK_REAL const PDupwindNthAnti1Xt3 = PDupwindNthAnti1(Xt3, i, j, k);
    CCTK_REAL const PDupwindNthSymm1Xt3 = PDupwindNthSymm1(Xt3, i, j, k);
    CCTK_REAL const PDupwindNthAnti2Xt3 = PDupwindNthAnti2(Xt3, i, j, k);
    CCTK_REAL const PDupwindNthSymm2Xt3 = PDupwindNthSymm2(Xt3, i, j, k);
    CCTK_REAL const PDupwindNthAnti3Xt3 = PDupwindNthAnti3(Xt3, i, j, k);
    CCTK_REAL const PDupwindNthSymm3Xt3 = PDupwindNthSymm3(Xt3, i, j, k);
    
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
    
    CCTK_REAL Gtl111 = khalf*(J11L*PDstandardNth1gt11 + 
      J21L*PDstandardNth2gt11 + J31L*PDstandardNth3gt11);
    
    CCTK_REAL Gtl112 = khalf*(J12L*PDstandardNth1gt11 + 
      J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11);
    
    CCTK_REAL Gtl113 = khalf*(J13L*PDstandardNth1gt11 + 
      J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11);
    
    CCTK_REAL Gtl122 = J12L*PDstandardNth1gt12 - 
      J11L*khalf*PDstandardNth1gt22 + J22L*PDstandardNth2gt12 - 
      J21L*khalf*PDstandardNth2gt22 + J32L*PDstandardNth3gt12 - 
      J31L*khalf*PDstandardNth3gt22;
    
    CCTK_REAL Gtl123 = khalf*(J13L*PDstandardNth1gt12 + 
      J12L*PDstandardNth1gt13 - J11L*PDstandardNth1gt23 + 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 - 
      J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 - J31L*PDstandardNth3gt23);
    
    CCTK_REAL Gtl133 = J13L*PDstandardNth1gt13 - 
      J11L*khalf*PDstandardNth1gt33 + J23L*PDstandardNth2gt13 - 
      J21L*khalf*PDstandardNth2gt33 + J33L*PDstandardNth3gt13 - 
      J31L*khalf*PDstandardNth3gt33;
    
    CCTK_REAL Gtl211 = khalf*(-(J12L*PDstandardNth1gt11) + 
      2*J11L*PDstandardNth1gt12 - J22L*PDstandardNth2gt11 + 
      2*J21L*PDstandardNth2gt12 - J32L*PDstandardNth3gt11 + 
      2*J31L*PDstandardNth3gt12);
    
    CCTK_REAL Gtl212 = khalf*(J11L*PDstandardNth1gt22 + 
      J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22);
    
    CCTK_REAL Gtl213 = khalf*(J13L*PDstandardNth1gt12 - 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 + 
      J23L*PDstandardNth2gt12 - J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 + J33L*PDstandardNth3gt12 - 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23);
    
    CCTK_REAL Gtl222 = khalf*(J12L*PDstandardNth1gt22 + 
      J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22);
    
    CCTK_REAL Gtl223 = khalf*(J13L*PDstandardNth1gt22 + 
      J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22);
    
    CCTK_REAL Gtl233 = J13L*PDstandardNth1gt23 - 
      J12L*khalf*PDstandardNth1gt33 + J23L*PDstandardNth2gt23 - 
      J22L*khalf*PDstandardNth2gt33 + J33L*PDstandardNth3gt23 - 
      J32L*khalf*PDstandardNth3gt33;
    
    CCTK_REAL Gtl311 = khalf*(-(J13L*PDstandardNth1gt11) + 
      2*J11L*PDstandardNth1gt13 - J23L*PDstandardNth2gt11 + 
      2*J21L*PDstandardNth2gt13 - J33L*PDstandardNth3gt11 + 
      2*J31L*PDstandardNth3gt13);
    
    CCTK_REAL Gtl312 = khalf*(-(J13L*PDstandardNth1gt12) + 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23);
    
    CCTK_REAL Gtl313 = khalf*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33);
    
    CCTK_REAL Gtl322 = khalf*(-(J13L*PDstandardNth1gt22) + 
      2*J12L*PDstandardNth1gt23 - J23L*PDstandardNth2gt22 + 
      2*J22L*PDstandardNth2gt23 - J33L*PDstandardNth3gt22 + 
      2*J32L*PDstandardNth3gt23);
    
    CCTK_REAL Gtl323 = khalf*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33);
    
    CCTK_REAL Gtl333 = khalf*(J13L*PDstandardNth1gt33 + 
      J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33);
    
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
    
    CCTK_REAL cdphi1 = fac1*(J11L*PDstandardNth1phi + 
      J21L*PDstandardNth2phi + J31L*PDstandardNth3phi);
    
    CCTK_REAL cdphi2 = fac1*(J12L*PDstandardNth1phi + 
      J22L*PDstandardNth2phi + J32L*PDstandardNth3phi);
    
    CCTK_REAL cdphi3 = fac1*(J13L*PDstandardNth1phi + 
      J23L*PDstandardNth2phi + J33L*PDstandardNth3phi);
    
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
    
    CCTK_REAL phirhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1phi + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2phi + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3phi + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1phi + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2phi + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3phi + (J11L*PDupwindNthSymm1phi + 
      J21L*PDupwindNthSymm2phi + J31L*PDupwindNthSymm3phi)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1phi + J22L*PDupwindNthSymm2phi + 
      J32L*PDupwindNthSymm3phi)*Abs(beta2L) + (J13L*PDupwindNthSymm1phi + 
      J23L*PDupwindNthSymm2phi + J33L*PDupwindNthSymm3phi)*Abs(beta3L) + 
      (J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3)*IfThen(conformalMethod,-(kthird*phiL),0.16666666666666666) 
      + alphaL*trKL*IfThen(conformalMethod,kthird*phiL,-0.16666666666666666);
    
    CCTK_REAL gt11rhsL = -2*alphaL*At11L + (epsdiss1*J11L + epsdiss2*J12L 
      + epsdiss3*J13L)*PDdissipationNth1gt11 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2gt11 + (epsdiss1*J31L + epsdiss2*J32L 
      + epsdiss3*J33L)*PDdissipationNth3gt11 - 
      gt11L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + 2*(J11L*(gt11L*PDstandardNth1beta1 + 
      gt12L*PDstandardNth1beta2 + gt13L*PDstandardNth1beta3) + 
      J21L*(gt11L*PDstandardNth2beta1 + gt12L*PDstandardNth2beta2 + 
      gt13L*PDstandardNth2beta3) + J31L*(gt11L*PDstandardNth3beta1 + 
      gt12L*PDstandardNth3beta2 + gt13L*PDstandardNth3beta3)) + (beta1L*J11L 
      + beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1gt11 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2gt11 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3gt11 + 
      (J11L*PDupwindNthSymm1gt11 + J21L*PDupwindNthSymm2gt11 + 
      J31L*PDupwindNthSymm3gt11)*Abs(beta1L) + (J12L*PDupwindNthSymm1gt11 + 
      J22L*PDupwindNthSymm2gt11 + J32L*PDupwindNthSymm3gt11)*Abs(beta2L) + 
      (J13L*PDupwindNthSymm1gt11 + J23L*PDupwindNthSymm2gt11 + 
      J33L*PDupwindNthSymm3gt11)*Abs(beta3L);
    
    CCTK_REAL gt12rhsL = -2*alphaL*At12L + (epsdiss1*J11L + epsdiss2*J12L 
      + epsdiss3*J13L)*PDdissipationNth1gt12 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2gt12 + (epsdiss1*J31L + epsdiss2*J32L 
      + epsdiss3*J33L)*PDdissipationNth3gt12 + (gt12L*J11L + 
      gt11L*J12L)*PDstandardNth1beta1 + (gt22L*J11L + 
      gt12L*J12L)*PDstandardNth1beta2 + (gt23L*J11L + 
      gt13L*J12L)*PDstandardNth1beta3 + (gt12L*J21L + 
      gt11L*J22L)*PDstandardNth2beta1 + (gt22L*J21L + 
      gt12L*J22L)*PDstandardNth2beta2 + (gt23L*J21L + 
      gt13L*J22L)*PDstandardNth2beta3 + (gt12L*J31L + 
      gt11L*J32L)*PDstandardNth3beta1 + (gt22L*J31L + 
      gt12L*J32L)*PDstandardNth3beta2 + (gt23L*J31L + 
      gt13L*J32L)*PDstandardNth3beta3 - 
      gt12L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1gt12 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2gt12 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3gt12 + (J11L*PDupwindNthSymm1gt12 + 
      J21L*PDupwindNthSymm2gt12 + J31L*PDupwindNthSymm3gt12)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1gt12 + J22L*PDupwindNthSymm2gt12 + 
      J32L*PDupwindNthSymm3gt12)*Abs(beta2L) + (J13L*PDupwindNthSymm1gt12 + 
      J23L*PDupwindNthSymm2gt12 + J33L*PDupwindNthSymm3gt12)*Abs(beta3L);
    
    CCTK_REAL gt13rhsL = -2*alphaL*At13L + (epsdiss1*J11L + epsdiss2*J12L 
      + epsdiss3*J13L)*PDdissipationNth1gt13 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2gt13 + (epsdiss1*J31L + epsdiss2*J32L 
      + epsdiss3*J33L)*PDdissipationNth3gt13 + (gt13L*J11L + 
      gt11L*J13L)*PDstandardNth1beta1 + (gt23L*J11L + 
      gt12L*J13L)*PDstandardNth1beta2 + (gt33L*J11L + 
      gt13L*J13L)*PDstandardNth1beta3 + (gt13L*J21L + 
      gt11L*J23L)*PDstandardNth2beta1 + (gt23L*J21L + 
      gt12L*J23L)*PDstandardNth2beta2 + (gt33L*J21L + 
      gt13L*J23L)*PDstandardNth2beta3 + (gt13L*J31L + 
      gt11L*J33L)*PDstandardNth3beta1 + (gt23L*J31L + 
      gt12L*J33L)*PDstandardNth3beta2 + (gt33L*J31L + 
      gt13L*J33L)*PDstandardNth3beta3 - 
      gt13L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1gt13 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2gt13 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3gt13 + (J11L*PDupwindNthSymm1gt13 + 
      J21L*PDupwindNthSymm2gt13 + J31L*PDupwindNthSymm3gt13)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1gt13 + J22L*PDupwindNthSymm2gt13 + 
      J32L*PDupwindNthSymm3gt13)*Abs(beta2L) + (J13L*PDupwindNthSymm1gt13 + 
      J23L*PDupwindNthSymm2gt13 + J33L*PDupwindNthSymm3gt13)*Abs(beta3L);
    
    CCTK_REAL gt22rhsL = -2*alphaL*At22L + (epsdiss1*J11L + epsdiss2*J12L 
      + epsdiss3*J13L)*PDdissipationNth1gt22 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2gt22 + (epsdiss1*J31L + epsdiss2*J32L 
      + epsdiss3*J33L)*PDdissipationNth3gt22 - 
      gt22L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + 2*(J12L*(gt12L*PDstandardNth1beta1 + 
      gt22L*PDstandardNth1beta2 + gt23L*PDstandardNth1beta3) + 
      J22L*(gt12L*PDstandardNth2beta1 + gt22L*PDstandardNth2beta2 + 
      gt23L*PDstandardNth2beta3) + J32L*(gt12L*PDstandardNth3beta1 + 
      gt22L*PDstandardNth3beta2 + gt23L*PDstandardNth3beta3)) + (beta1L*J11L 
      + beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1gt22 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2gt22 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3gt22 + 
      (J11L*PDupwindNthSymm1gt22 + J21L*PDupwindNthSymm2gt22 + 
      J31L*PDupwindNthSymm3gt22)*Abs(beta1L) + (J12L*PDupwindNthSymm1gt22 + 
      J22L*PDupwindNthSymm2gt22 + J32L*PDupwindNthSymm3gt22)*Abs(beta2L) + 
      (J13L*PDupwindNthSymm1gt22 + J23L*PDupwindNthSymm2gt22 + 
      J33L*PDupwindNthSymm3gt22)*Abs(beta3L);
    
    CCTK_REAL gt23rhsL = -2*alphaL*At23L + (epsdiss1*J11L + epsdiss2*J12L 
      + epsdiss3*J13L)*PDdissipationNth1gt23 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2gt23 + (epsdiss1*J31L + epsdiss2*J32L 
      + epsdiss3*J33L)*PDdissipationNth3gt23 + (gt13L*J12L + 
      gt12L*J13L)*PDstandardNth1beta1 + (gt23L*J12L + 
      gt22L*J13L)*PDstandardNth1beta2 + (gt33L*J12L + 
      gt23L*J13L)*PDstandardNth1beta3 + (gt13L*J22L + 
      gt12L*J23L)*PDstandardNth2beta1 + (gt23L*J22L + 
      gt22L*J23L)*PDstandardNth2beta2 + (gt33L*J22L + 
      gt23L*J23L)*PDstandardNth2beta3 + (gt13L*J32L + 
      gt12L*J33L)*PDstandardNth3beta1 + (gt23L*J32L + 
      gt22L*J33L)*PDstandardNth3beta2 + (gt33L*J32L + 
      gt23L*J33L)*PDstandardNth3beta3 - 
      gt23L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1gt23 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2gt23 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3gt23 + (J11L*PDupwindNthSymm1gt23 + 
      J21L*PDupwindNthSymm2gt23 + J31L*PDupwindNthSymm3gt23)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1gt23 + J22L*PDupwindNthSymm2gt23 + 
      J32L*PDupwindNthSymm3gt23)*Abs(beta2L) + (J13L*PDupwindNthSymm1gt23 + 
      J23L*PDupwindNthSymm2gt23 + J33L*PDupwindNthSymm3gt23)*Abs(beta3L);
    
    CCTK_REAL gt33rhsL = -2*alphaL*At33L + (epsdiss1*J11L + epsdiss2*J12L 
      + epsdiss3*J13L)*PDdissipationNth1gt33 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2gt33 + (epsdiss1*J31L + epsdiss2*J32L 
      + epsdiss3*J33L)*PDdissipationNth3gt33 - 
      gt33L*ktwothird*(J11L*PDstandardNth1beta1 + J12L*PDstandardNth1beta2 + 
      J13L*PDstandardNth1beta3 + J21L*PDstandardNth2beta1 + 
      J22L*PDstandardNth2beta2 + J23L*PDstandardNth2beta3 + 
      J31L*PDstandardNth3beta1 + J32L*PDstandardNth3beta2 + 
      J33L*PDstandardNth3beta3) + 2*(J13L*(gt13L*PDstandardNth1beta1 + 
      gt23L*PDstandardNth1beta2 + gt33L*PDstandardNth1beta3) + 
      J23L*(gt13L*PDstandardNth2beta1 + gt23L*PDstandardNth2beta2 + 
      gt33L*PDstandardNth2beta3) + J33L*(gt13L*PDstandardNth3beta1 + 
      gt23L*PDstandardNth3beta2 + gt33L*PDstandardNth3beta3)) + (beta1L*J11L 
      + beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1gt33 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2gt33 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3gt33 + 
      (J11L*PDupwindNthSymm1gt33 + J21L*PDupwindNthSymm2gt33 + 
      J31L*PDupwindNthSymm3gt33)*Abs(beta1L) + (J12L*PDupwindNthSymm1gt33 + 
      J22L*PDupwindNthSymm2gt33 + J32L*PDupwindNthSymm3gt33)*Abs(beta2L) + 
      (J13L*PDupwindNthSymm1gt33 + J23L*PDupwindNthSymm2gt33 + 
      J33L*PDupwindNthSymm3gt33)*Abs(beta3L);
    
    CCTK_REAL dotXt1 = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1Xt1 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2Xt1 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3Xt1 - 2*((Atu11*J11L + Atu12*J12L + 
      Atu13*J13L)*PDstandardNth1alpha + (Atu11*J21L + Atu12*J22L + 
      Atu13*J23L)*PDstandardNth2alpha + (Atu11*J31L + Atu12*J32L + 
      Atu13*J33L)*PDstandardNth3alpha) + 
      2*(gtu12*(J11L*J12L*PDstandardNth11beta1 + 
      J12L*J21L*PDstandardNth12beta1 + J11L*J22L*PDstandardNth12beta1 + 
      J12L*J31L*PDstandardNth13beta1 + J11L*J32L*PDstandardNth13beta1 + 
      dJ112L*PDstandardNth1beta1 + J21L*J22L*PDstandardNth22beta1 + 
      J22L*J31L*PDstandardNth23beta1 + J21L*J32L*PDstandardNth23beta1 + 
      dJ212L*PDstandardNth2beta1 + J31L*J32L*PDstandardNth33beta1 + 
      dJ312L*PDstandardNth3beta1) + gtu13*(J11L*J13L*PDstandardNth11beta1 + 
      J13L*J21L*PDstandardNth12beta1 + J11L*J23L*PDstandardNth12beta1 + 
      J13L*J31L*PDstandardNth13beta1 + J11L*J33L*PDstandardNth13beta1 + 
      dJ113L*PDstandardNth1beta1 + J21L*J23L*PDstandardNth22beta1 + 
      J23L*J31L*PDstandardNth23beta1 + J21L*J33L*PDstandardNth23beta1 + 
      dJ213L*PDstandardNth2beta1 + J31L*J33L*PDstandardNth33beta1 + 
      dJ313L*PDstandardNth3beta1) + gtu23*(J12L*J13L*PDstandardNth11beta1 + 
      J13L*J22L*PDstandardNth12beta1 + J12L*J23L*PDstandardNth12beta1 + 
      J13L*J32L*PDstandardNth13beta1 + J12L*J33L*PDstandardNth13beta1 + 
      dJ123L*PDstandardNth1beta1 + J22L*J23L*PDstandardNth22beta1 + 
      J23L*J32L*PDstandardNth23beta1 + J22L*J33L*PDstandardNth23beta1 + 
      dJ223L*PDstandardNth2beta1 + J32L*J33L*PDstandardNth33beta1 + 
      dJ323L*PDstandardNth3beta1) + alphaL*(6*(Atu11*cdphi1 + Atu12*cdphi2 + 
      Atu13*cdphi3) + Atu11*Gt111 + 2*Atu12*Gt112 + 2*Atu13*Gt113 + 
      Atu22*Gt122 + 2*Atu23*Gt123 + Atu33*Gt133 - ktwothird*((gtu11*J11L + 
      gtu12*J12L + gtu13*J13L)*PDstandardNth1trK + (gtu11*J21L + gtu12*J22L + 
      gtu13*J23L)*PDstandardNth2trK + (gtu11*J31L + gtu12*J32L + 
      gtu13*J33L)*PDstandardNth3trK))) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1Xt1 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2Xt1 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3Xt1 - 
      50.26548245743669181540229413247204614715*alphaL*(gtu11*S1 + gtu12*S2 + 
      gtu13*S3) + ktwothird*(J11L*PDstandardNth1beta1 + 
      J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
      J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + 
      J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
      J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3)*Xtn1 - 
      PDstandardNth1beta1*(J11L*Xtn1 + J12L*Xtn2 + J13L*Xtn3) - 
      PDstandardNth2beta1*(J21L*Xtn1 + J22L*Xtn2 + J23L*Xtn3) - 
      PDstandardNth3beta1*(J31L*Xtn1 + J32L*Xtn2 + J33L*Xtn3) + 
      (J11L*PDupwindNthSymm1Xt1 + J21L*PDupwindNthSymm2Xt1 + 
      J31L*PDupwindNthSymm3Xt1)*Abs(beta1L) + (J12L*PDupwindNthSymm1Xt1 + 
      J22L*PDupwindNthSymm2Xt1 + J32L*PDupwindNthSymm3Xt1)*Abs(beta2L) + 
      (J13L*PDupwindNthSymm1Xt1 + J23L*PDupwindNthSymm2Xt1 + 
      J33L*PDupwindNthSymm3Xt1)*Abs(beta3L) + 
      gtu11*(2*J11L*J21L*PDstandardNth12beta1 + 
      2*J11L*J31L*PDstandardNth13beta1 + dJ111L*PDstandardNth1beta1 + 
      2*J21L*J31L*PDstandardNth23beta1 + dJ211L*PDstandardNth2beta1 + 
      dJ311L*PDstandardNth3beta1 + PDstandardNth11beta1*SQR(J11L) + 
      PDstandardNth22beta1*SQR(J21L) + PDstandardNth33beta1*SQR(J31L)) + 
      gtu22*(2*J12L*J22L*PDstandardNth12beta1 + 
      2*J12L*J32L*PDstandardNth13beta1 + dJ122L*PDstandardNth1beta1 + 
      2*J22L*J32L*PDstandardNth23beta1 + dJ222L*PDstandardNth2beta1 + 
      dJ322L*PDstandardNth3beta1 + PDstandardNth11beta1*SQR(J12L) + 
      PDstandardNth22beta1*SQR(J22L) + PDstandardNth33beta1*SQR(J32L)) + 
      gtu33*(2*J13L*J23L*PDstandardNth12beta1 + 
      2*J13L*J33L*PDstandardNth13beta1 + dJ133L*PDstandardNth1beta1 + 
      2*J23L*J33L*PDstandardNth23beta1 + dJ233L*PDstandardNth2beta1 + 
      dJ333L*PDstandardNth3beta1 + PDstandardNth11beta1*SQR(J13L) + 
      PDstandardNth22beta1*SQR(J23L) + PDstandardNth33beta1*SQR(J33L)) + 
      kthird*(gtu11*(J11L*J12L*PDstandardNth11beta2 + 
      J11L*J13L*PDstandardNth11beta3 + 2*J11L*J21L*PDstandardNth12beta1 + 
      J12L*J21L*PDstandardNth12beta2 + J11L*J22L*PDstandardNth12beta2 + 
      J13L*J21L*PDstandardNth12beta3 + J11L*J23L*PDstandardNth12beta3 + 
      2*J11L*J31L*PDstandardNth13beta1 + J12L*J31L*PDstandardNth13beta2 + 
      J11L*J32L*PDstandardNth13beta2 + J13L*J31L*PDstandardNth13beta3 + 
      J11L*J33L*PDstandardNth13beta3 + dJ111L*PDstandardNth1beta1 + 
      dJ112L*PDstandardNth1beta2 + dJ113L*PDstandardNth1beta3 + 
      J21L*J22L*PDstandardNth22beta2 + J21L*J23L*PDstandardNth22beta3 + 
      2*J21L*J31L*PDstandardNth23beta1 + J22L*J31L*PDstandardNth23beta2 + 
      J21L*J32L*PDstandardNth23beta2 + J23L*J31L*PDstandardNth23beta3 + 
      J21L*J33L*PDstandardNth23beta3 + dJ211L*PDstandardNth2beta1 + 
      dJ212L*PDstandardNth2beta2 + dJ213L*PDstandardNth2beta3 + 
      J31L*J32L*PDstandardNth33beta2 + J31L*J33L*PDstandardNth33beta3 + 
      dJ311L*PDstandardNth3beta1 + dJ312L*PDstandardNth3beta2 + 
      dJ313L*PDstandardNth3beta3 + PDstandardNth11beta1*SQR(J11L) + 
      PDstandardNth22beta1*SQR(J21L) + PDstandardNth33beta1*SQR(J31L)) + 
      gtu12*(J11L*J12L*PDstandardNth11beta1 + J12L*J13L*PDstandardNth11beta3 
      + J12L*J21L*PDstandardNth12beta1 + J11L*J22L*PDstandardNth12beta1 + 
      2*J12L*J22L*PDstandardNth12beta2 + J13L*J22L*PDstandardNth12beta3 + 
      J12L*J23L*PDstandardNth12beta3 + J12L*J31L*PDstandardNth13beta1 + 
      J11L*J32L*PDstandardNth13beta1 + 2*J12L*J32L*PDstandardNth13beta2 + 
      J13L*J32L*PDstandardNth13beta3 + J12L*J33L*PDstandardNth13beta3 + 
      dJ112L*PDstandardNth1beta1 + dJ122L*PDstandardNth1beta2 + 
      dJ123L*PDstandardNth1beta3 + J21L*J22L*PDstandardNth22beta1 + 
      J22L*J23L*PDstandardNth22beta3 + J22L*J31L*PDstandardNth23beta1 + 
      J21L*J32L*PDstandardNth23beta1 + 2*J22L*J32L*PDstandardNth23beta2 + 
      J23L*J32L*PDstandardNth23beta3 + J22L*J33L*PDstandardNth23beta3 + 
      dJ212L*PDstandardNth2beta1 + dJ222L*PDstandardNth2beta2 + 
      dJ223L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta1 + 
      J32L*J33L*PDstandardNth33beta3 + dJ312L*PDstandardNth3beta1 + 
      dJ322L*PDstandardNth3beta2 + dJ323L*PDstandardNth3beta3 + 
      PDstandardNth11beta2*SQR(J12L) + PDstandardNth22beta2*SQR(J22L) + 
      PDstandardNth33beta2*SQR(J32L)) + gtu13*(J11L*J13L*PDstandardNth11beta1 
      + J12L*J13L*PDstandardNth11beta2 + J13L*J21L*PDstandardNth12beta1 + 
      J11L*J23L*PDstandardNth12beta1 + J13L*J22L*PDstandardNth12beta2 + 
      J12L*J23L*PDstandardNth12beta2 + 2*J13L*J23L*PDstandardNth12beta3 + 
      J13L*J31L*PDstandardNth13beta1 + J11L*J33L*PDstandardNth13beta1 + 
      J13L*J32L*PDstandardNth13beta2 + J12L*J33L*PDstandardNth13beta2 + 
      2*J13L*J33L*PDstandardNth13beta3 + dJ113L*PDstandardNth1beta1 + 
      dJ123L*PDstandardNth1beta2 + dJ133L*PDstandardNth1beta3 + 
      J21L*J23L*PDstandardNth22beta1 + J22L*J23L*PDstandardNth22beta2 + 
      J23L*J31L*PDstandardNth23beta1 + J21L*J33L*PDstandardNth23beta1 + 
      J23L*J32L*PDstandardNth23beta2 + J22L*J33L*PDstandardNth23beta2 + 
      2*J23L*J33L*PDstandardNth23beta3 + dJ213L*PDstandardNth2beta1 + 
      dJ223L*PDstandardNth2beta2 + dJ233L*PDstandardNth2beta3 + 
      J31L*J33L*PDstandardNth33beta1 + J32L*J33L*PDstandardNth33beta2 + 
      dJ313L*PDstandardNth3beta1 + dJ323L*PDstandardNth3beta2 + 
      dJ333L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J13L) + 
      PDstandardNth22beta3*SQR(J23L) + PDstandardNth33beta3*SQR(J33L)));
    
    CCTK_REAL dotXt2 = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1Xt2 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2Xt2 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3Xt2 - 2*((Atu12*J11L + Atu22*J12L + 
      Atu23*J13L)*PDstandardNth1alpha + (Atu12*J21L + Atu22*J22L + 
      Atu23*J23L)*PDstandardNth2alpha + (Atu12*J31L + Atu22*J32L + 
      Atu23*J33L)*PDstandardNth3alpha) + 
      2*(gtu12*(J11L*J12L*PDstandardNth11beta2 + 
      J12L*J21L*PDstandardNth12beta2 + J11L*J22L*PDstandardNth12beta2 + 
      J12L*J31L*PDstandardNth13beta2 + J11L*J32L*PDstandardNth13beta2 + 
      dJ112L*PDstandardNth1beta2 + J21L*J22L*PDstandardNth22beta2 + 
      J22L*J31L*PDstandardNth23beta2 + J21L*J32L*PDstandardNth23beta2 + 
      dJ212L*PDstandardNth2beta2 + J31L*J32L*PDstandardNth33beta2 + 
      dJ312L*PDstandardNth3beta2) + gtu13*(J11L*J13L*PDstandardNth11beta2 + 
      J13L*J21L*PDstandardNth12beta2 + J11L*J23L*PDstandardNth12beta2 + 
      J13L*J31L*PDstandardNth13beta2 + J11L*J33L*PDstandardNth13beta2 + 
      dJ113L*PDstandardNth1beta2 + J21L*J23L*PDstandardNth22beta2 + 
      J23L*J31L*PDstandardNth23beta2 + J21L*J33L*PDstandardNth23beta2 + 
      dJ213L*PDstandardNth2beta2 + J31L*J33L*PDstandardNth33beta2 + 
      dJ313L*PDstandardNth3beta2) + gtu23*(J12L*J13L*PDstandardNth11beta2 + 
      J13L*J22L*PDstandardNth12beta2 + J12L*J23L*PDstandardNth12beta2 + 
      J13L*J32L*PDstandardNth13beta2 + J12L*J33L*PDstandardNth13beta2 + 
      dJ123L*PDstandardNth1beta2 + J22L*J23L*PDstandardNth22beta2 + 
      J23L*J32L*PDstandardNth23beta2 + J22L*J33L*PDstandardNth23beta2 + 
      dJ223L*PDstandardNth2beta2 + J32L*J33L*PDstandardNth33beta2 + 
      dJ323L*PDstandardNth3beta2) + alphaL*(6*(Atu12*cdphi1 + Atu22*cdphi2 + 
      Atu23*cdphi3) + Atu11*Gt211 + 2*Atu12*Gt212 + 2*Atu13*Gt213 + 
      Atu22*Gt222 + 2*Atu23*Gt223 + Atu33*Gt233 - ktwothird*((gtu12*J11L + 
      gtu22*J12L + gtu23*J13L)*PDstandardNth1trK + (gtu12*J21L + gtu22*J22L + 
      gtu23*J23L)*PDstandardNth2trK + (gtu12*J31L + gtu22*J32L + 
      gtu23*J33L)*PDstandardNth3trK))) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1Xt2 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2Xt2 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3Xt2 - 
      50.26548245743669181540229413247204614715*alphaL*(gtu12*S1 + gtu22*S2 + 
      gtu23*S3) + ktwothird*(J11L*PDstandardNth1beta1 + 
      J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
      J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + 
      J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
      J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3)*Xtn2 - 
      PDstandardNth1beta2*(J11L*Xtn1 + J12L*Xtn2 + J13L*Xtn3) - 
      PDstandardNth2beta2*(J21L*Xtn1 + J22L*Xtn2 + J23L*Xtn3) - 
      PDstandardNth3beta2*(J31L*Xtn1 + J32L*Xtn2 + J33L*Xtn3) + 
      (J11L*PDupwindNthSymm1Xt2 + J21L*PDupwindNthSymm2Xt2 + 
      J31L*PDupwindNthSymm3Xt2)*Abs(beta1L) + (J12L*PDupwindNthSymm1Xt2 + 
      J22L*PDupwindNthSymm2Xt2 + J32L*PDupwindNthSymm3Xt2)*Abs(beta2L) + 
      (J13L*PDupwindNthSymm1Xt2 + J23L*PDupwindNthSymm2Xt2 + 
      J33L*PDupwindNthSymm3Xt2)*Abs(beta3L) + 
      gtu11*(2*J11L*J21L*PDstandardNth12beta2 + 
      2*J11L*J31L*PDstandardNth13beta2 + dJ111L*PDstandardNth1beta2 + 
      2*J21L*J31L*PDstandardNth23beta2 + dJ211L*PDstandardNth2beta2 + 
      dJ311L*PDstandardNth3beta2 + PDstandardNth11beta2*SQR(J11L) + 
      PDstandardNth22beta2*SQR(J21L) + PDstandardNth33beta2*SQR(J31L)) + 
      gtu22*(2*J12L*J22L*PDstandardNth12beta2 + 
      2*J12L*J32L*PDstandardNth13beta2 + dJ122L*PDstandardNth1beta2 + 
      2*J22L*J32L*PDstandardNth23beta2 + dJ222L*PDstandardNth2beta2 + 
      dJ322L*PDstandardNth3beta2 + PDstandardNth11beta2*SQR(J12L) + 
      PDstandardNth22beta2*SQR(J22L) + PDstandardNth33beta2*SQR(J32L)) + 
      gtu33*(2*J13L*J23L*PDstandardNth12beta2 + 
      2*J13L*J33L*PDstandardNth13beta2 + dJ133L*PDstandardNth1beta2 + 
      2*J23L*J33L*PDstandardNth23beta2 + dJ233L*PDstandardNth2beta2 + 
      dJ333L*PDstandardNth3beta2 + PDstandardNth11beta2*SQR(J13L) + 
      PDstandardNth22beta2*SQR(J23L) + PDstandardNth33beta2*SQR(J33L)) + 
      kthird*(gtu12*(J11L*J12L*PDstandardNth11beta2 + 
      J11L*J13L*PDstandardNth11beta3 + 2*J11L*J21L*PDstandardNth12beta1 + 
      J12L*J21L*PDstandardNth12beta2 + J11L*J22L*PDstandardNth12beta2 + 
      J13L*J21L*PDstandardNth12beta3 + J11L*J23L*PDstandardNth12beta3 + 
      2*J11L*J31L*PDstandardNth13beta1 + J12L*J31L*PDstandardNth13beta2 + 
      J11L*J32L*PDstandardNth13beta2 + J13L*J31L*PDstandardNth13beta3 + 
      J11L*J33L*PDstandardNth13beta3 + dJ111L*PDstandardNth1beta1 + 
      dJ112L*PDstandardNth1beta2 + dJ113L*PDstandardNth1beta3 + 
      J21L*J22L*PDstandardNth22beta2 + J21L*J23L*PDstandardNth22beta3 + 
      2*J21L*J31L*PDstandardNth23beta1 + J22L*J31L*PDstandardNth23beta2 + 
      J21L*J32L*PDstandardNth23beta2 + J23L*J31L*PDstandardNth23beta3 + 
      J21L*J33L*PDstandardNth23beta3 + dJ211L*PDstandardNth2beta1 + 
      dJ212L*PDstandardNth2beta2 + dJ213L*PDstandardNth2beta3 + 
      J31L*J32L*PDstandardNth33beta2 + J31L*J33L*PDstandardNth33beta3 + 
      dJ311L*PDstandardNth3beta1 + dJ312L*PDstandardNth3beta2 + 
      dJ313L*PDstandardNth3beta3 + PDstandardNth11beta1*SQR(J11L) + 
      PDstandardNth22beta1*SQR(J21L) + PDstandardNth33beta1*SQR(J31L)) + 
      gtu22*(J11L*J12L*PDstandardNth11beta1 + J12L*J13L*PDstandardNth11beta3 
      + J12L*J21L*PDstandardNth12beta1 + J11L*J22L*PDstandardNth12beta1 + 
      2*J12L*J22L*PDstandardNth12beta2 + J13L*J22L*PDstandardNth12beta3 + 
      J12L*J23L*PDstandardNth12beta3 + J12L*J31L*PDstandardNth13beta1 + 
      J11L*J32L*PDstandardNth13beta1 + 2*J12L*J32L*PDstandardNth13beta2 + 
      J13L*J32L*PDstandardNth13beta3 + J12L*J33L*PDstandardNth13beta3 + 
      dJ112L*PDstandardNth1beta1 + dJ122L*PDstandardNth1beta2 + 
      dJ123L*PDstandardNth1beta3 + J21L*J22L*PDstandardNth22beta1 + 
      J22L*J23L*PDstandardNth22beta3 + J22L*J31L*PDstandardNth23beta1 + 
      J21L*J32L*PDstandardNth23beta1 + 2*J22L*J32L*PDstandardNth23beta2 + 
      J23L*J32L*PDstandardNth23beta3 + J22L*J33L*PDstandardNth23beta3 + 
      dJ212L*PDstandardNth2beta1 + dJ222L*PDstandardNth2beta2 + 
      dJ223L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta1 + 
      J32L*J33L*PDstandardNth33beta3 + dJ312L*PDstandardNth3beta1 + 
      dJ322L*PDstandardNth3beta2 + dJ323L*PDstandardNth3beta3 + 
      PDstandardNth11beta2*SQR(J12L) + PDstandardNth22beta2*SQR(J22L) + 
      PDstandardNth33beta2*SQR(J32L)) + gtu23*(J11L*J13L*PDstandardNth11beta1 
      + J12L*J13L*PDstandardNth11beta2 + J13L*J21L*PDstandardNth12beta1 + 
      J11L*J23L*PDstandardNth12beta1 + J13L*J22L*PDstandardNth12beta2 + 
      J12L*J23L*PDstandardNth12beta2 + 2*J13L*J23L*PDstandardNth12beta3 + 
      J13L*J31L*PDstandardNth13beta1 + J11L*J33L*PDstandardNth13beta1 + 
      J13L*J32L*PDstandardNth13beta2 + J12L*J33L*PDstandardNth13beta2 + 
      2*J13L*J33L*PDstandardNth13beta3 + dJ113L*PDstandardNth1beta1 + 
      dJ123L*PDstandardNth1beta2 + dJ133L*PDstandardNth1beta3 + 
      J21L*J23L*PDstandardNth22beta1 + J22L*J23L*PDstandardNth22beta2 + 
      J23L*J31L*PDstandardNth23beta1 + J21L*J33L*PDstandardNth23beta1 + 
      J23L*J32L*PDstandardNth23beta2 + J22L*J33L*PDstandardNth23beta2 + 
      2*J23L*J33L*PDstandardNth23beta3 + dJ213L*PDstandardNth2beta1 + 
      dJ223L*PDstandardNth2beta2 + dJ233L*PDstandardNth2beta3 + 
      J31L*J33L*PDstandardNth33beta1 + J32L*J33L*PDstandardNth33beta2 + 
      dJ313L*PDstandardNth3beta1 + dJ323L*PDstandardNth3beta2 + 
      dJ333L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J13L) + 
      PDstandardNth22beta3*SQR(J23L) + PDstandardNth33beta3*SQR(J33L)));
    
    CCTK_REAL dotXt3 = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1Xt3 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2Xt3 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3Xt3 - 2*((Atu13*J11L + Atu23*J12L + 
      Atu33*J13L)*PDstandardNth1alpha + (Atu13*J21L + Atu23*J22L + 
      Atu33*J23L)*PDstandardNth2alpha + (Atu13*J31L + Atu23*J32L + 
      Atu33*J33L)*PDstandardNth3alpha) + 
      2*(gtu12*(J11L*J12L*PDstandardNth11beta3 + 
      J12L*J21L*PDstandardNth12beta3 + J11L*J22L*PDstandardNth12beta3 + 
      J12L*J31L*PDstandardNth13beta3 + J11L*J32L*PDstandardNth13beta3 + 
      dJ112L*PDstandardNth1beta3 + J21L*J22L*PDstandardNth22beta3 + 
      J22L*J31L*PDstandardNth23beta3 + J21L*J32L*PDstandardNth23beta3 + 
      dJ212L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta3 + 
      dJ312L*PDstandardNth3beta3) + gtu13*(J11L*J13L*PDstandardNth11beta3 + 
      J13L*J21L*PDstandardNth12beta3 + J11L*J23L*PDstandardNth12beta3 + 
      J13L*J31L*PDstandardNth13beta3 + J11L*J33L*PDstandardNth13beta3 + 
      dJ113L*PDstandardNth1beta3 + J21L*J23L*PDstandardNth22beta3 + 
      J23L*J31L*PDstandardNth23beta3 + J21L*J33L*PDstandardNth23beta3 + 
      dJ213L*PDstandardNth2beta3 + J31L*J33L*PDstandardNth33beta3 + 
      dJ313L*PDstandardNth3beta3) + gtu23*(J12L*J13L*PDstandardNth11beta3 + 
      J13L*J22L*PDstandardNth12beta3 + J12L*J23L*PDstandardNth12beta3 + 
      J13L*J32L*PDstandardNth13beta3 + J12L*J33L*PDstandardNth13beta3 + 
      dJ123L*PDstandardNth1beta3 + J22L*J23L*PDstandardNth22beta3 + 
      J23L*J32L*PDstandardNth23beta3 + J22L*J33L*PDstandardNth23beta3 + 
      dJ223L*PDstandardNth2beta3 + J32L*J33L*PDstandardNth33beta3 + 
      dJ323L*PDstandardNth3beta3) + alphaL*(6*(Atu13*cdphi1 + Atu23*cdphi2 + 
      Atu33*cdphi3) + Atu11*Gt311 + 2*Atu12*Gt312 + 2*Atu13*Gt313 + 
      Atu22*Gt322 + 2*Atu23*Gt323 + Atu33*Gt333 - ktwothird*((gtu13*J11L + 
      gtu23*J12L + gtu33*J13L)*PDstandardNth1trK + (gtu13*J21L + gtu23*J22L + 
      gtu33*J23L)*PDstandardNth2trK + (gtu13*J31L + gtu23*J32L + 
      gtu33*J33L)*PDstandardNth3trK))) + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1Xt3 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2Xt3 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3Xt3 - 
      50.26548245743669181540229413247204614715*alphaL*(gtu13*S1 + gtu23*S2 + 
      gtu33*S3) + ktwothird*(J11L*PDstandardNth1beta1 + 
      J12L*PDstandardNth1beta2 + J13L*PDstandardNth1beta3 + 
      J21L*PDstandardNth2beta1 + J22L*PDstandardNth2beta2 + 
      J23L*PDstandardNth2beta3 + J31L*PDstandardNth3beta1 + 
      J32L*PDstandardNth3beta2 + J33L*PDstandardNth3beta3)*Xtn3 - 
      PDstandardNth1beta3*(J11L*Xtn1 + J12L*Xtn2 + J13L*Xtn3) - 
      PDstandardNth2beta3*(J21L*Xtn1 + J22L*Xtn2 + J23L*Xtn3) - 
      PDstandardNth3beta3*(J31L*Xtn1 + J32L*Xtn2 + J33L*Xtn3) + 
      (J11L*PDupwindNthSymm1Xt3 + J21L*PDupwindNthSymm2Xt3 + 
      J31L*PDupwindNthSymm3Xt3)*Abs(beta1L) + (J12L*PDupwindNthSymm1Xt3 + 
      J22L*PDupwindNthSymm2Xt3 + J32L*PDupwindNthSymm3Xt3)*Abs(beta2L) + 
      (J13L*PDupwindNthSymm1Xt3 + J23L*PDupwindNthSymm2Xt3 + 
      J33L*PDupwindNthSymm3Xt3)*Abs(beta3L) + 
      gtu11*(2*J11L*J21L*PDstandardNth12beta3 + 
      2*J11L*J31L*PDstandardNth13beta3 + dJ111L*PDstandardNth1beta3 + 
      2*J21L*J31L*PDstandardNth23beta3 + dJ211L*PDstandardNth2beta3 + 
      dJ311L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J11L) + 
      PDstandardNth22beta3*SQR(J21L) + PDstandardNth33beta3*SQR(J31L)) + 
      gtu22*(2*J12L*J22L*PDstandardNth12beta3 + 
      2*J12L*J32L*PDstandardNth13beta3 + dJ122L*PDstandardNth1beta3 + 
      2*J22L*J32L*PDstandardNth23beta3 + dJ222L*PDstandardNth2beta3 + 
      dJ322L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J12L) + 
      PDstandardNth22beta3*SQR(J22L) + PDstandardNth33beta3*SQR(J32L)) + 
      gtu33*(2*J13L*J23L*PDstandardNth12beta3 + 
      2*J13L*J33L*PDstandardNth13beta3 + dJ133L*PDstandardNth1beta3 + 
      2*J23L*J33L*PDstandardNth23beta3 + dJ233L*PDstandardNth2beta3 + 
      dJ333L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J13L) + 
      PDstandardNth22beta3*SQR(J23L) + PDstandardNth33beta3*SQR(J33L)) + 
      kthird*(gtu13*(J11L*J12L*PDstandardNth11beta2 + 
      J11L*J13L*PDstandardNth11beta3 + 2*J11L*J21L*PDstandardNth12beta1 + 
      J12L*J21L*PDstandardNth12beta2 + J11L*J22L*PDstandardNth12beta2 + 
      J13L*J21L*PDstandardNth12beta3 + J11L*J23L*PDstandardNth12beta3 + 
      2*J11L*J31L*PDstandardNth13beta1 + J12L*J31L*PDstandardNth13beta2 + 
      J11L*J32L*PDstandardNth13beta2 + J13L*J31L*PDstandardNth13beta3 + 
      J11L*J33L*PDstandardNth13beta3 + dJ111L*PDstandardNth1beta1 + 
      dJ112L*PDstandardNth1beta2 + dJ113L*PDstandardNth1beta3 + 
      J21L*J22L*PDstandardNth22beta2 + J21L*J23L*PDstandardNth22beta3 + 
      2*J21L*J31L*PDstandardNth23beta1 + J22L*J31L*PDstandardNth23beta2 + 
      J21L*J32L*PDstandardNth23beta2 + J23L*J31L*PDstandardNth23beta3 + 
      J21L*J33L*PDstandardNth23beta3 + dJ211L*PDstandardNth2beta1 + 
      dJ212L*PDstandardNth2beta2 + dJ213L*PDstandardNth2beta3 + 
      J31L*J32L*PDstandardNth33beta2 + J31L*J33L*PDstandardNth33beta3 + 
      dJ311L*PDstandardNth3beta1 + dJ312L*PDstandardNth3beta2 + 
      dJ313L*PDstandardNth3beta3 + PDstandardNth11beta1*SQR(J11L) + 
      PDstandardNth22beta1*SQR(J21L) + PDstandardNth33beta1*SQR(J31L)) + 
      gtu23*(J11L*J12L*PDstandardNth11beta1 + J12L*J13L*PDstandardNth11beta3 
      + J12L*J21L*PDstandardNth12beta1 + J11L*J22L*PDstandardNth12beta1 + 
      2*J12L*J22L*PDstandardNth12beta2 + J13L*J22L*PDstandardNth12beta3 + 
      J12L*J23L*PDstandardNth12beta3 + J12L*J31L*PDstandardNth13beta1 + 
      J11L*J32L*PDstandardNth13beta1 + 2*J12L*J32L*PDstandardNth13beta2 + 
      J13L*J32L*PDstandardNth13beta3 + J12L*J33L*PDstandardNth13beta3 + 
      dJ112L*PDstandardNth1beta1 + dJ122L*PDstandardNth1beta2 + 
      dJ123L*PDstandardNth1beta3 + J21L*J22L*PDstandardNth22beta1 + 
      J22L*J23L*PDstandardNth22beta3 + J22L*J31L*PDstandardNth23beta1 + 
      J21L*J32L*PDstandardNth23beta1 + 2*J22L*J32L*PDstandardNth23beta2 + 
      J23L*J32L*PDstandardNth23beta3 + J22L*J33L*PDstandardNth23beta3 + 
      dJ212L*PDstandardNth2beta1 + dJ222L*PDstandardNth2beta2 + 
      dJ223L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta1 + 
      J32L*J33L*PDstandardNth33beta3 + dJ312L*PDstandardNth3beta1 + 
      dJ322L*PDstandardNth3beta2 + dJ323L*PDstandardNth3beta3 + 
      PDstandardNth11beta2*SQR(J12L) + PDstandardNth22beta2*SQR(J22L) + 
      PDstandardNth33beta2*SQR(J32L)) + gtu33*(J11L*J13L*PDstandardNth11beta1 
      + J12L*J13L*PDstandardNth11beta2 + J13L*J21L*PDstandardNth12beta1 + 
      J11L*J23L*PDstandardNth12beta1 + J13L*J22L*PDstandardNth12beta2 + 
      J12L*J23L*PDstandardNth12beta2 + 2*J13L*J23L*PDstandardNth12beta3 + 
      J13L*J31L*PDstandardNth13beta1 + J11L*J33L*PDstandardNth13beta1 + 
      J13L*J32L*PDstandardNth13beta2 + J12L*J33L*PDstandardNth13beta2 + 
      2*J13L*J33L*PDstandardNth13beta3 + dJ113L*PDstandardNth1beta1 + 
      dJ123L*PDstandardNth1beta2 + dJ133L*PDstandardNth1beta3 + 
      J21L*J23L*PDstandardNth22beta1 + J22L*J23L*PDstandardNth22beta2 + 
      J23L*J31L*PDstandardNth23beta1 + J21L*J33L*PDstandardNth23beta1 + 
      J23L*J32L*PDstandardNth23beta2 + J22L*J33L*PDstandardNth23beta2 + 
      2*J23L*J33L*PDstandardNth23beta3 + dJ213L*PDstandardNth2beta1 + 
      dJ223L*PDstandardNth2beta2 + dJ233L*PDstandardNth2beta3 + 
      J31L*J33L*PDstandardNth33beta1 + J32L*J33L*PDstandardNth33beta2 + 
      dJ313L*PDstandardNth3beta1 + dJ323L*PDstandardNth3beta2 + 
      dJ333L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J13L) + 
      PDstandardNth22beta3*SQR(J23L) + PDstandardNth33beta3*SQR(J33L)));
    
    CCTK_REAL Xt1rhsL = dotXt1;
    
    CCTK_REAL Xt2rhsL = dotXt2;
    
    CCTK_REAL Xt3rhsL = dotXt3;
    
    CCTK_REAL dottrK = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1trK + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2trK + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3trK + (beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1trK + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2trK + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3trK + (J11L*PDupwindNthSymm1trK + 
      J21L*PDupwindNthSymm2trK + J31L*PDupwindNthSymm3trK)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1trK + J22L*PDupwindNthSymm2trK + 
      J32L*PDupwindNthSymm3trK)*Abs(beta2L) + (J13L*PDupwindNthSymm1trK + 
      J23L*PDupwindNthSymm2trK + J33L*PDupwindNthSymm3trK)*Abs(beta3L) - 
      em4phi*(2*((gtu12*J11L*J22L + gtu13*(J13L*J21L + 
      J11L*J23L))*PDstandardNth12alpha + (gtu12*J11L*J32L + gtu13*(J13L*J31L 
      + J11L*J33L))*PDstandardNth13alpha + J11L*((gtu12*J12L + 
      gtu13*J13L)*PDstandardNth11alpha + gtu11*(J21L*PDstandardNth12alpha + 
      J31L*PDstandardNth13alpha)) + J12L*(gtu23*J13L*PDstandardNth11alpha + 
      gtu12*(J21L*PDstandardNth12alpha + J31L*PDstandardNth13alpha)) + 
      (gtu11*J21L*J31L + (gtu22*J22L + gtu23*J23L)*J32L + (gtu23*J22L + 
      gtu33*J23L)*J33L)*PDstandardNth23alpha + J22L*((gtu22*J12L + 
      gtu23*J13L)*PDstandardNth12alpha + (gtu12*J21L + 
      gtu23*J23L)*PDstandardNth22alpha + gtu12*J31L*PDstandardNth23alpha) + 
      J23L*((gtu23*J12L + gtu33*J13L)*PDstandardNth12alpha + 
      gtu13*(J21L*PDstandardNth22alpha + J31L*PDstandardNth23alpha)) + 
      (gtu12*(dJ312L + cdphi2*J31L) + cdphi3*(gtu13*J31L + 
      gtu23*J32L))*PDstandardNth3alpha + J32L*((gtu22*J12L + 
      gtu23*J13L)*PDstandardNth13alpha + gtu23*J33L*PDstandardNth33alpha + 
      gtu12*(J21L*PDstandardNth23alpha + J31L*PDstandardNth33alpha) + 
      cdphi2*gtu22*PDstandardNth3alpha) + J33L*((gtu23*J12L + 
      gtu33*J13L)*PDstandardNth13alpha + gtu13*(J21L*PDstandardNth23alpha + 
      J31L*PDstandardNth33alpha) + (cdphi2*gtu23 + 
      cdphi3*gtu33)*PDstandardNth3alpha)) + 
      PDstandardNth1alpha*(gtu11*(dJ111L + 2*cdphi1*J11L) + gtu22*(dJ122L + 
      2*cdphi2*J12L) + gtu33*(dJ133L + 2*cdphi3*J13L) + 2*(dJ112L*gtu12 + 
      dJ113L*gtu13 + dJ123L*gtu23 + cdphi2*gtu12*J11L + cdphi3*gtu13*J11L + 
      cdphi1*gtu12*J12L + cdphi3*gtu23*J12L + cdphi1*gtu13*J13L + 
      cdphi2*gtu23*J13L) - J11L*Xtn1 - J12L*Xtn2 - J13L*Xtn3) + 
      PDstandardNth2alpha*(gtu11*(dJ211L + 2*cdphi1*J21L) + gtu22*(dJ222L + 
      2*cdphi2*J22L) + gtu33*(dJ233L + 2*cdphi3*J23L) + 2*(dJ212L*gtu12 + 
      dJ213L*gtu13 + dJ223L*gtu23 + cdphi2*gtu12*J21L + cdphi3*gtu13*J21L + 
      cdphi1*gtu12*J22L + cdphi3*gtu23*J22L + cdphi1*gtu13*J23L + 
      cdphi2*gtu23*J23L) - J21L*Xtn1 - J22L*Xtn2 - J23L*Xtn3) + 
      PDstandardNth3alpha*(dJ322L*gtu22 + dJ333L*gtu33 + gtu11*(dJ311L + 
      2*cdphi1*J31L) + 2*(dJ313L*gtu13 + dJ323L*gtu23 + cdphi1*gtu12*J32L + 
      cdphi1*gtu13*J33L) - J31L*Xtn1 - J32L*Xtn2 - J33L*Xtn3) + 
      PDstandardNth11alpha*(gtu11*SQR(J11L) + gtu22*SQR(J12L) + 
      gtu33*SQR(J13L)) + PDstandardNth22alpha*(gtu11*SQR(J21L) + 
      gtu22*SQR(J22L) + gtu33*SQR(J23L)) + 
      PDstandardNth33alpha*(gtu11*SQR(J31L) + gtu22*SQR(J32L) + 
      gtu33*SQR(J33L))) + alphaL*(2*(Atm12*Atm21 + Atm13*Atm31 + Atm23*Atm32) 
      + 12.56637061435917295385057353311801153679*(rho + trS) + SQR(Atm11) + 
      SQR(Atm22) + SQR(Atm33) + kthird*SQR(trKL));
    
    CCTK_REAL trKrhsL = dottrK;
    
    CCTK_REAL alpharhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1alpha + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2alpha + (epsdiss1*J31L + 
      epsdiss2*J32L + epsdiss3*J33L)*PDdissipationNth3alpha + 
      LapseAdvectionCoeff*((beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1alpha + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2alpha + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3alpha + (J11L*PDupwindNthSymm1alpha + 
      J21L*PDupwindNthSymm2alpha + J31L*PDupwindNthSymm3alpha)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1alpha + J22L*PDupwindNthSymm2alpha + 
      J32L*PDupwindNthSymm3alpha)*Abs(beta2L) + (J13L*PDupwindNthSymm1alpha + 
      J23L*PDupwindNthSymm2alpha + J33L*PDupwindNthSymm3alpha)*Abs(beta3L)) - 
      harmonicF*(LapseACoeff*(AL - trKL) + trKL)*pow(alphaL,harmonicN);
    
    CCTK_REAL ArhsL = (-(AL*AlphaDriver) + dottrK)*LapseACoeff + 
      (epsdiss1*J11L + epsdiss2*J12L + epsdiss3*J13L)*PDdissipationNth1A + 
      (epsdiss1*J21L + epsdiss2*J22L + epsdiss3*J23L)*PDdissipationNth2A + 
      (epsdiss1*J31L + epsdiss2*J32L + epsdiss3*J33L)*PDdissipationNth3A + 
      LapseAdvectionCoeff*((beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1A + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2A + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3A + (J11L*PDupwindNthSymm1A + 
      J21L*PDupwindNthSymm2A + J31L*PDupwindNthSymm3A)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1A + J22L*PDupwindNthSymm2A + 
      J32L*PDupwindNthSymm3A)*Abs(beta2L) + (J13L*PDupwindNthSymm1A + 
      J23L*PDupwindNthSymm2A + J33L*PDupwindNthSymm3A)*Abs(beta3L));
    
    CCTK_REAL eta = fmin(1,SpatialBetaDriverRadius*INV(rL));
    
    CCTK_REAL theta = fmin(1,exp(1 - 
      rL*INV(SpatialShiftGammaCoeffRadius)));
    
    CCTK_REAL beta1rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1beta1 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2beta1 + (epsdiss1*J31L + 
      epsdiss2*J32L + epsdiss3*J33L)*PDdissipationNth3beta1 + 
      ShiftGammaCoeff*theta*(beta1L*BetaDriver*eta*(-1 + ShiftBCoeff) + 
      ShiftBCoeff*(B1L - Xt1L) + Xt1L) + ShiftAdvectionCoeff*((beta1L*J11L + 
      beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1beta1 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2beta1 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3beta1 + 
      (J11L*PDupwindNthSymm1beta1 + J21L*PDupwindNthSymm2beta1 + 
      J31L*PDupwindNthSymm3beta1)*Abs(beta1L) + (J12L*PDupwindNthSymm1beta1 + 
      J22L*PDupwindNthSymm2beta1 + J32L*PDupwindNthSymm3beta1)*Abs(beta2L) + 
      (J13L*PDupwindNthSymm1beta1 + J23L*PDupwindNthSymm2beta1 + 
      J33L*PDupwindNthSymm3beta1)*Abs(beta3L));
    
    CCTK_REAL beta2rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1beta2 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2beta2 + (epsdiss1*J31L + 
      epsdiss2*J32L + epsdiss3*J33L)*PDdissipationNth3beta2 + 
      ShiftGammaCoeff*theta*(beta2L*BetaDriver*eta*(-1 + ShiftBCoeff) + 
      ShiftBCoeff*(B2L - Xt2L) + Xt2L) + ShiftAdvectionCoeff*((beta1L*J11L + 
      beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1beta2 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2beta2 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3beta2 + 
      (J11L*PDupwindNthSymm1beta2 + J21L*PDupwindNthSymm2beta2 + 
      J31L*PDupwindNthSymm3beta2)*Abs(beta1L) + (J12L*PDupwindNthSymm1beta2 + 
      J22L*PDupwindNthSymm2beta2 + J32L*PDupwindNthSymm3beta2)*Abs(beta2L) + 
      (J13L*PDupwindNthSymm1beta2 + J23L*PDupwindNthSymm2beta2 + 
      J33L*PDupwindNthSymm3beta2)*Abs(beta3L));
    
    CCTK_REAL beta3rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1beta3 + (epsdiss1*J21L + epsdiss2*J22L 
      + epsdiss3*J23L)*PDdissipationNth2beta3 + (epsdiss1*J31L + 
      epsdiss2*J32L + epsdiss3*J33L)*PDdissipationNth3beta3 + 
      ShiftGammaCoeff*theta*(beta3L*BetaDriver*eta*(-1 + ShiftBCoeff) + 
      ShiftBCoeff*(B3L - Xt3L) + Xt3L) + ShiftAdvectionCoeff*((beta1L*J11L + 
      beta2L*J12L + beta3L*J13L)*PDupwindNthAnti1beta3 + (beta1L*J21L + 
      beta2L*J22L + beta3L*J23L)*PDupwindNthAnti2beta3 + (beta1L*J31L + 
      beta2L*J32L + beta3L*J33L)*PDupwindNthAnti3beta3 + 
      (J11L*PDupwindNthSymm1beta3 + J21L*PDupwindNthSymm2beta3 + 
      J31L*PDupwindNthSymm3beta3)*Abs(beta1L) + (J12L*PDupwindNthSymm1beta3 + 
      J22L*PDupwindNthSymm2beta3 + J32L*PDupwindNthSymm3beta3)*Abs(beta2L) + 
      (J13L*PDupwindNthSymm1beta3 + J23L*PDupwindNthSymm2beta3 + 
      J33L*PDupwindNthSymm3beta3)*Abs(beta3L));
    
    CCTK_REAL B1rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1B1 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2B1 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3B1 + (dotXt1 - 
      B1L*BetaDriver*eta)*ShiftBCoeff + ShiftAdvectionCoeff*((beta1L*J11L + 
      beta2L*J12L + beta3L*J13L)*(PDupwindNthAnti1B1 - PDupwindNthAnti1Xt1) + 
      (beta1L*J21L + beta2L*J22L + beta3L*J23L)*(PDupwindNthAnti2B1 - 
      PDupwindNthAnti2Xt1) + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*(PDupwindNthAnti3B1 - PDupwindNthAnti3Xt1) + 
      (J11L*(PDupwindNthSymm1B1 - PDupwindNthSymm1Xt1) + 
      J21L*(PDupwindNthSymm2B1 - PDupwindNthSymm2Xt1) + 
      J31L*(PDupwindNthSymm3B1 - PDupwindNthSymm3Xt1))*Abs(beta1L) + 
      (J12L*(PDupwindNthSymm1B1 - PDupwindNthSymm1Xt1) + 
      J22L*(PDupwindNthSymm2B1 - PDupwindNthSymm2Xt1) + 
      J32L*(PDupwindNthSymm3B1 - PDupwindNthSymm3Xt1))*Abs(beta2L) + 
      (J13L*(PDupwindNthSymm1B1 - PDupwindNthSymm1Xt1) + 
      J23L*(PDupwindNthSymm2B1 - PDupwindNthSymm2Xt1) + 
      J33L*(PDupwindNthSymm3B1 - PDupwindNthSymm3Xt1))*Abs(beta3L));
    
    CCTK_REAL B2rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1B2 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2B2 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3B2 + (dotXt2 - 
      B2L*BetaDriver*eta)*ShiftBCoeff + ShiftAdvectionCoeff*((beta1L*J11L + 
      beta2L*J12L + beta3L*J13L)*(PDupwindNthAnti1B2 - PDupwindNthAnti1Xt2) + 
      (beta1L*J21L + beta2L*J22L + beta3L*J23L)*(PDupwindNthAnti2B2 - 
      PDupwindNthAnti2Xt2) + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*(PDupwindNthAnti3B2 - PDupwindNthAnti3Xt2) + 
      (J11L*(PDupwindNthSymm1B2 - PDupwindNthSymm1Xt2) + 
      J21L*(PDupwindNthSymm2B2 - PDupwindNthSymm2Xt2) + 
      J31L*(PDupwindNthSymm3B2 - PDupwindNthSymm3Xt2))*Abs(beta1L) + 
      (J12L*(PDupwindNthSymm1B2 - PDupwindNthSymm1Xt2) + 
      J22L*(PDupwindNthSymm2B2 - PDupwindNthSymm2Xt2) + 
      J32L*(PDupwindNthSymm3B2 - PDupwindNthSymm3Xt2))*Abs(beta2L) + 
      (J13L*(PDupwindNthSymm1B2 - PDupwindNthSymm1Xt2) + 
      J23L*(PDupwindNthSymm2B2 - PDupwindNthSymm2Xt2) + 
      J33L*(PDupwindNthSymm3B2 - PDupwindNthSymm3Xt2))*Abs(beta3L));
    
    CCTK_REAL B3rhsL = (epsdiss1*J11L + epsdiss2*J12L + 
      epsdiss3*J13L)*PDdissipationNth1B3 + (epsdiss1*J21L + epsdiss2*J22L + 
      epsdiss3*J23L)*PDdissipationNth2B3 + (epsdiss1*J31L + epsdiss2*J32L + 
      epsdiss3*J33L)*PDdissipationNth3B3 + (dotXt3 - 
      B3L*BetaDriver*eta)*ShiftBCoeff + ShiftAdvectionCoeff*((beta1L*J11L + 
      beta2L*J12L + beta3L*J13L)*(PDupwindNthAnti1B3 - PDupwindNthAnti1Xt3) + 
      (beta1L*J21L + beta2L*J22L + beta3L*J23L)*(PDupwindNthAnti2B3 - 
      PDupwindNthAnti2Xt3) + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*(PDupwindNthAnti3B3 - PDupwindNthAnti3Xt3) + 
      (J11L*(PDupwindNthSymm1B3 - PDupwindNthSymm1Xt3) + 
      J21L*(PDupwindNthSymm2B3 - PDupwindNthSymm2Xt3) + 
      J31L*(PDupwindNthSymm3B3 - PDupwindNthSymm3Xt3))*Abs(beta1L) + 
      (J12L*(PDupwindNthSymm1B3 - PDupwindNthSymm1Xt3) + 
      J22L*(PDupwindNthSymm2B3 - PDupwindNthSymm2Xt3) + 
      J32L*(PDupwindNthSymm3B3 - PDupwindNthSymm3Xt3))*Abs(beta2L) + 
      (J13L*(PDupwindNthSymm1B3 - PDupwindNthSymm1Xt3) + 
      J23L*(PDupwindNthSymm2B3 - PDupwindNthSymm2Xt3) + 
      J33L*(PDupwindNthSymm3B3 - PDupwindNthSymm3Xt3))*Abs(beta3L));
    
    
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
  LC_ENDLOOP3 (ML_BSSN_MP_RHS1);
}

void ML_BSSN_MP_RHS1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_RHS1_Body);
}
