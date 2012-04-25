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

extern "C" void HOST_ML_BSSN_RHS_Dalpha_2_etc_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthalpha22_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthalpha22_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthalpha2_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthalpha2_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta122_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta122_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta12_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta12_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta222_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta222_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta22_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta22_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta322_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta322_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthbeta32_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthbeta32_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt1122_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt1122_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt112_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt112_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt1222_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt1222_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt122_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt122_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt1322_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt1322_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt132_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt132_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt2222_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt2222_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt222_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt222_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt2322_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt2322_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt232_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt232_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt3322_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt3322_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthgt332_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthgt332_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthphi22_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthphi22_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthphi2_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthphi2_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthtrK2_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthtrK2_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthXt12_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthXt12_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthXt22_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthXt22_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::DPDstandardNthXt32_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::DPDstandardNthXt32_group.");
  return;
}

static void HOST_ML_BSSN_RHS_Dalpha_2_etc_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_LOOP3(HOST_ML_BSSN_RHS_Dalpha_2_etc,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Calculate temporaries and grid functions */
    DPDstandardNthalpha2[index] = PDstandardNth2(&alpha[index]);
    DPDstandardNthalpha22[index] = PDstandardNth22(&alpha[index]);
    DPDstandardNthbeta12[index] = PDstandardNth2(&beta1[index]);
    DPDstandardNthbeta122[index] = PDstandardNth22(&beta1[index]);
    DPDstandardNthbeta22[index] = PDstandardNth2(&beta2[index]);
    DPDstandardNthbeta222[index] = PDstandardNth22(&beta2[index]);
    DPDstandardNthbeta32[index] = PDstandardNth2(&beta3[index]);
    DPDstandardNthbeta322[index] = PDstandardNth22(&beta3[index]);
    DPDstandardNthgt112[index] = PDstandardNth2(&gt11[index]);
    DPDstandardNthgt1122[index] = PDstandardNth22(&gt11[index]);
    DPDstandardNthgt122[index] = PDstandardNth2(&gt12[index]);
    DPDstandardNthgt1222[index] = PDstandardNth22(&gt12[index]);
    DPDstandardNthgt132[index] = PDstandardNth2(&gt13[index]);
    DPDstandardNthgt1322[index] = PDstandardNth22(&gt13[index]);
    DPDstandardNthgt222[index] = PDstandardNth2(&gt22[index]);
    DPDstandardNthgt2222[index] = PDstandardNth22(&gt22[index]);
    DPDstandardNthgt232[index] = PDstandardNth2(&gt23[index]);
    DPDstandardNthgt2322[index] = PDstandardNth22(&gt23[index]);
    DPDstandardNthgt332[index] = PDstandardNth2(&gt33[index]);
    DPDstandardNthgt3322[index] = PDstandardNth22(&gt33[index]);
    DPDstandardNthphi2[index] = PDstandardNth2(&phi[index]);
    DPDstandardNthphi22[index] = PDstandardNth22(&phi[index]);
    DPDstandardNthtrK2[index] = PDstandardNth2(&trK[index]);
    DPDstandardNthXt12[index] = PDstandardNth2(&Xt1[index]);
    DPDstandardNthXt22[index] = PDstandardNth2(&Xt2[index]);
    DPDstandardNthXt32[index] = PDstandardNth2(&Xt3[index]);
  }
  CCTK_ENDLOOP3(HOST_ML_BSSN_RHS_Dalpha_2_etc);
}

extern "C" void HOST_ML_BSSN_RHS_Dalpha_2_etc(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering HOST_ML_BSSN_RHS_Dalpha_2_etc_Body");
  }
  
  if (cctk_iteration % HOST_ML_BSSN_RHS_Dalpha_2_etc_calc_every != HOST_ML_BSSN_RHS_Dalpha_2_etc_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "ML_BSSN::DPDstandardNthalpha22_group",
    "ML_BSSN::DPDstandardNthalpha2_group",
    "ML_BSSN::DPDstandardNthbeta122_group",
    "ML_BSSN::DPDstandardNthbeta12_group",
    "ML_BSSN::DPDstandardNthbeta222_group",
    "ML_BSSN::DPDstandardNthbeta22_group",
    "ML_BSSN::DPDstandardNthbeta322_group",
    "ML_BSSN::DPDstandardNthbeta32_group",
    "ML_BSSN::DPDstandardNthgt1122_group",
    "ML_BSSN::DPDstandardNthgt112_group",
    "ML_BSSN::DPDstandardNthgt1222_group",
    "ML_BSSN::DPDstandardNthgt122_group",
    "ML_BSSN::DPDstandardNthgt1322_group",
    "ML_BSSN::DPDstandardNthgt132_group",
    "ML_BSSN::DPDstandardNthgt2222_group",
    "ML_BSSN::DPDstandardNthgt222_group",
    "ML_BSSN::DPDstandardNthgt2322_group",
    "ML_BSSN::DPDstandardNthgt232_group",
    "ML_BSSN::DPDstandardNthgt3322_group",
    "ML_BSSN::DPDstandardNthgt332_group",
    "ML_BSSN::DPDstandardNthphi22_group",
    "ML_BSSN::DPDstandardNthphi2_group",
    "ML_BSSN::DPDstandardNthtrK2_group",
    "ML_BSSN::DPDstandardNthXt12_group",
    "ML_BSSN::DPDstandardNthXt22_group",
    "ML_BSSN::DPDstandardNthXt32_group",
    "ML_BSSN::ML_Gamma",
    "ML_BSSN::ML_lapse",
    "ML_BSSN::ML_log_confac",
    "ML_BSSN::ML_metric",
    "ML_BSSN::ML_shift",
    "ML_BSSN::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "HOST_ML_BSSN_RHS_Dalpha_2_etc", 32, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "HOST_ML_BSSN_RHS_Dalpha_2_etc", 4, 4, 4);
  
  GenericFD_LoopOverInterior(cctkGH, HOST_ML_BSSN_RHS_Dalpha_2_etc_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving HOST_ML_BSSN_RHS_Dalpha_2_etc_Body");
  }
}
