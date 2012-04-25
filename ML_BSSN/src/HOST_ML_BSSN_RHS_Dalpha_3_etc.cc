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
#define QAD(x) (SQR(SQR(x)))
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))

extern "C" void HOST_ML_BSSN_RHS_Dalpha_3_etc_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthalpha33_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthalpha33_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthalpha3_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthalpha3_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta133_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta133_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta13_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta13_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta233_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta233_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta23_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta23_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta333_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta333_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta33_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta33_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt1133_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt1133_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt113_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt113_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt1233_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt1233_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt123_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt123_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt1333_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt1333_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt133_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt133_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt2233_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt2233_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt223_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt223_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt2333_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt2333_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt233_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt233_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt3333_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt3333_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt333_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt333_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthphi33_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthphi33_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthphi3_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthphi3_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthtrK3_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthtrK3_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthXt13_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthXt13_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthXt23_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthXt23_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthXt33_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthXt33_group.");
  return;
}

static void HOST_ML_BSSN_RHS_Dalpha_3_etc_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_REAL const pm1o840dx = -0.00119047619047619047619047619048*INV(dx);
  CCTK_REAL const pm1o840dy = -0.00119047619047619047619047619048*INV(dy);
  CCTK_REAL const pm1o840dz = -0.00119047619047619047619047619048*INV(dz);
  
  /* Loop over the grid points */
  #pragma omp parallel
  CCTK_LOOP3(HOST_ML_BSSN_RHS_Dalpha_3_etc,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Calculate temporaries and grid functions */
    DPDstandardNthalpha3[index] = PDstandardNth3(&alpha[index]);
    DPDstandardNthalpha33[index] = PDstandardNth33(&alpha[index]);
    DPDstandardNthbeta13[index] = PDstandardNth3(&beta1[index]);
    DPDstandardNthbeta133[index] = PDstandardNth33(&beta1[index]);
    DPDstandardNthbeta23[index] = PDstandardNth3(&beta2[index]);
    DPDstandardNthbeta233[index] = PDstandardNth33(&beta2[index]);
    DPDstandardNthbeta33[index] = PDstandardNth3(&beta3[index]);
    DPDstandardNthbeta333[index] = PDstandardNth33(&beta3[index]);
    DPDstandardNthgt113[index] = PDstandardNth3(&gt11[index]);
    DPDstandardNthgt1133[index] = PDstandardNth33(&gt11[index]);
    DPDstandardNthgt123[index] = PDstandardNth3(&gt12[index]);
    DPDstandardNthgt1233[index] = PDstandardNth33(&gt12[index]);
    DPDstandardNthgt133[index] = PDstandardNth3(&gt13[index]);
    DPDstandardNthgt1333[index] = PDstandardNth33(&gt13[index]);
    DPDstandardNthgt223[index] = PDstandardNth3(&gt22[index]);
    DPDstandardNthgt2233[index] = PDstandardNth33(&gt22[index]);
    DPDstandardNthgt233[index] = PDstandardNth3(&gt23[index]);
    DPDstandardNthgt2333[index] = PDstandardNth33(&gt23[index]);
    DPDstandardNthgt333[index] = PDstandardNth3(&gt33[index]);
    DPDstandardNthgt3333[index] = PDstandardNth33(&gt33[index]);
    DPDstandardNthphi3[index] = PDstandardNth3(&phi[index]);
    DPDstandardNthphi33[index] = PDstandardNth33(&phi[index]);
    DPDstandardNthtrK3[index] = PDstandardNth3(&trK[index]);
    DPDstandardNthXt13[index] = PDstandardNth3(&Xt1[index]);
    DPDstandardNthXt23[index] = PDstandardNth3(&Xt2[index]);
    DPDstandardNthXt33[index] = PDstandardNth3(&Xt3[index]);
  }
  CCTK_ENDLOOP3(HOST_ML_BSSN_RHS_Dalpha_3_etc);
}

extern "C" void HOST_ML_BSSN_RHS_Dalpha_3_etc(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering HOST_ML_BSSN_RHS_Dalpha_3_etc_Body");
  }
  
  if (cctk_iteration % HOST_ML_BSSN_RHS_Dalpha_3_etc_calc_every != HOST_ML_BSSN_RHS_Dalpha_3_etc_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "ML_BSSN::DPDstandardNthalpha33_group",
    "ML_BSSN::DPDstandardNthalpha3_group",
    "ML_BSSN::DPDstandardNthbeta133_group",
    "ML_BSSN::DPDstandardNthbeta13_group",
    "ML_BSSN::DPDstandardNthbeta233_group",
    "ML_BSSN::DPDstandardNthbeta23_group",
    "ML_BSSN::DPDstandardNthbeta333_group",
    "ML_BSSN::DPDstandardNthbeta33_group",
    "ML_BSSN::DPDstandardNthgt1133_group",
    "ML_BSSN::DPDstandardNthgt113_group",
    "ML_BSSN::DPDstandardNthgt1233_group",
    "ML_BSSN::DPDstandardNthgt123_group",
    "ML_BSSN::DPDstandardNthgt1333_group",
    "ML_BSSN::DPDstandardNthgt133_group",
    "ML_BSSN::DPDstandardNthgt2233_group",
    "ML_BSSN::DPDstandardNthgt223_group",
    "ML_BSSN::DPDstandardNthgt2333_group",
    "ML_BSSN::DPDstandardNthgt233_group",
    "ML_BSSN::DPDstandardNthgt3333_group",
    "ML_BSSN::DPDstandardNthgt333_group",
    "ML_BSSN::DPDstandardNthphi33_group",
    "ML_BSSN::DPDstandardNthphi3_group",
    "ML_BSSN::DPDstandardNthtrK3_group",
    "ML_BSSN::DPDstandardNthXt13_group",
    "ML_BSSN::DPDstandardNthXt23_group",
    "ML_BSSN::DPDstandardNthXt33_group",
    "ML_BSSN::ML_Gamma",
    "ML_BSSN::ML_lapse",
    "ML_BSSN::ML_log_confac",
    "ML_BSSN::ML_metric",
    "ML_BSSN::ML_shift",
    "ML_BSSN::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "HOST_ML_BSSN_RHS_Dalpha_3_etc", 32, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "HOST_ML_BSSN_RHS_Dalpha_3_etc", 4, 4, 4);
  
  GenericFD_LoopOverInterior(cctkGH, HOST_ML_BSSN_RHS_Dalpha_3_etc_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving HOST_ML_BSSN_RHS_Dalpha_3_etc_Body");
  }
}
