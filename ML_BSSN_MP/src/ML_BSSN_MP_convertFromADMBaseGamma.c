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

void ML_BSSN_MP_convertFromADMBaseGamma_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_dtlapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_dtshift.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP::ML_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP::ML_Gamma.");
  return;
}

void ML_BSSN_MP_convertFromADMBaseGamma_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_convertFromADMBaseGamma_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_convertFromADMBaseGamma_calc_every != ML_BSSN_MP_convertFromADMBaseGamma_calc_offset)
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
  LC_LOOP3 (ML_BSSN_MP_convertFromADMBaseGamma,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    // CCTK_REAL PDupwindNthAnti1alpha = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1alpha = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2alpha = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2alpha = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3alpha = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3alpha = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1beta1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1beta1 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2beta1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2beta1 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3beta1 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3beta1 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1beta2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1beta2 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2beta2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2beta2 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3beta2 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3beta2 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti1beta3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm1beta3 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti2beta3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm2beta3 = INITVALUE;
    // CCTK_REAL PDupwindNthAnti3beta3 = INITVALUE;
    // CCTK_REAL PDupwindNthSymm3beta3 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt11 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt12 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt13 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt22 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt23 = INITVALUE;
    // CCTK_REAL PDstandardNth1gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth2gt33 = INITVALUE;
    // CCTK_REAL PDstandardNth3gt33 = INITVALUE;
    
    /* Assign local copies of grid functions */
    CCTK_REAL  alphaL = alpha[index];
    CCTK_REAL  beta1L = beta1[index];
    CCTK_REAL  beta2L = beta2[index];
    CCTK_REAL  beta3L = beta3[index];
    CCTK_REAL  dtalpL = dtalp[index];
    CCTK_REAL  dtbetaxL = dtbetax[index];
    CCTK_REAL  dtbetayL = dtbetay[index];
    CCTK_REAL  dtbetazL = dtbetaz[index];
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
    CCTK_REAL  rL = r[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDupwindNthAnti1alpha = PDupwindNthAnti1(alpha, i, j, k);
    CCTK_REAL const PDupwindNthSymm1alpha = PDupwindNthSymm1(alpha, i, j, k);
    CCTK_REAL const PDupwindNthAnti2alpha = PDupwindNthAnti2(alpha, i, j, k);
    CCTK_REAL const PDupwindNthSymm2alpha = PDupwindNthSymm2(alpha, i, j, k);
    CCTK_REAL const PDupwindNthAnti3alpha = PDupwindNthAnti3(alpha, i, j, k);
    CCTK_REAL const PDupwindNthSymm3alpha = PDupwindNthSymm3(alpha, i, j, k);
    CCTK_REAL const PDupwindNthAnti1beta1 = PDupwindNthAnti1(beta1, i, j, k);
    CCTK_REAL const PDupwindNthSymm1beta1 = PDupwindNthSymm1(beta1, i, j, k);
    CCTK_REAL const PDupwindNthAnti2beta1 = PDupwindNthAnti2(beta1, i, j, k);
    CCTK_REAL const PDupwindNthSymm2beta1 = PDupwindNthSymm2(beta1, i, j, k);
    CCTK_REAL const PDupwindNthAnti3beta1 = PDupwindNthAnti3(beta1, i, j, k);
    CCTK_REAL const PDupwindNthSymm3beta1 = PDupwindNthSymm3(beta1, i, j, k);
    CCTK_REAL const PDupwindNthAnti1beta2 = PDupwindNthAnti1(beta2, i, j, k);
    CCTK_REAL const PDupwindNthSymm1beta2 = PDupwindNthSymm1(beta2, i, j, k);
    CCTK_REAL const PDupwindNthAnti2beta2 = PDupwindNthAnti2(beta2, i, j, k);
    CCTK_REAL const PDupwindNthSymm2beta2 = PDupwindNthSymm2(beta2, i, j, k);
    CCTK_REAL const PDupwindNthAnti3beta2 = PDupwindNthAnti3(beta2, i, j, k);
    CCTK_REAL const PDupwindNthSymm3beta2 = PDupwindNthSymm3(beta2, i, j, k);
    CCTK_REAL const PDupwindNthAnti1beta3 = PDupwindNthAnti1(beta3, i, j, k);
    CCTK_REAL const PDupwindNthSymm1beta3 = PDupwindNthSymm1(beta3, i, j, k);
    CCTK_REAL const PDupwindNthAnti2beta3 = PDupwindNthAnti2(beta3, i, j, k);
    CCTK_REAL const PDupwindNthSymm2beta3 = PDupwindNthSymm2(beta3, i, j, k);
    CCTK_REAL const PDupwindNthAnti3beta3 = PDupwindNthAnti3(beta3, i, j, k);
    CCTK_REAL const PDupwindNthSymm3beta3 = PDupwindNthSymm3(beta3, i, j, k);
    CCTK_REAL const PDstandardNth1gt11 = PDstandardNth1(gt11, i, j, k);
    CCTK_REAL const PDstandardNth2gt11 = PDstandardNth2(gt11, i, j, k);
    CCTK_REAL const PDstandardNth3gt11 = PDstandardNth3(gt11, i, j, k);
    CCTK_REAL const PDstandardNth1gt12 = PDstandardNth1(gt12, i, j, k);
    CCTK_REAL const PDstandardNth2gt12 = PDstandardNth2(gt12, i, j, k);
    CCTK_REAL const PDstandardNth3gt12 = PDstandardNth3(gt12, i, j, k);
    CCTK_REAL const PDstandardNth1gt13 = PDstandardNth1(gt13, i, j, k);
    CCTK_REAL const PDstandardNth2gt13 = PDstandardNth2(gt13, i, j, k);
    CCTK_REAL const PDstandardNth3gt13 = PDstandardNth3(gt13, i, j, k);
    CCTK_REAL const PDstandardNth1gt22 = PDstandardNth1(gt22, i, j, k);
    CCTK_REAL const PDstandardNth2gt22 = PDstandardNth2(gt22, i, j, k);
    CCTK_REAL const PDstandardNth3gt22 = PDstandardNth3(gt22, i, j, k);
    CCTK_REAL const PDstandardNth1gt23 = PDstandardNth1(gt23, i, j, k);
    CCTK_REAL const PDstandardNth2gt23 = PDstandardNth2(gt23, i, j, k);
    CCTK_REAL const PDstandardNth3gt23 = PDstandardNth3(gt23, i, j, k);
    CCTK_REAL const PDstandardNth1gt33 = PDstandardNth1(gt33, i, j, k);
    CCTK_REAL const PDstandardNth2gt33 = PDstandardNth2(gt33, i, j, k);
    CCTK_REAL const PDstandardNth3gt33 = PDstandardNth3(gt33, i, j, k);
    
    /* Calculate temporaries and grid functions */
    int dir1 = Sign(beta1L);
    
    int dir2 = Sign(beta2L);
    
    int dir3 = Sign(beta3L);
    
    CCTK_REAL detgt = 1;
    
    CCTK_REAL gtu11 = INV(detgt)*(gt22L*gt33L - SQR(gt23L));
    
    CCTK_REAL gtu12 = (gt13L*gt23L - gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu13 = (-(gt13L*gt22L) + gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 = INV(detgt)*(gt11L*gt33L - SQR(gt13L));
    
    CCTK_REAL gtu23 = (gt12L*gt13L - gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 = INV(detgt)*(gt11L*gt22L - SQR(gt12L));
    
    CCTK_REAL Gt111 = khalf*((gtu11*J11L - gtu12*J12L - 
      gtu13*J13L)*PDstandardNth1gt11 + (gtu11*J21L - gtu12*J22L - 
      gtu13*J23L)*PDstandardNth2gt11 + (gtu11*J31L - gtu12*J32L - 
      gtu13*J33L)*PDstandardNth3gt11 + 2*(J11L*(gtu12*PDstandardNth1gt12 + 
      gtu13*PDstandardNth1gt13) + J21L*(gtu12*PDstandardNth2gt12 + 
      gtu13*PDstandardNth2gt13) + J31L*(gtu12*PDstandardNth3gt12 + 
      gtu13*PDstandardNth3gt13)));
    
    CCTK_REAL Gt211 = khalf*((gtu12*J11L - gtu22*J12L - 
      gtu23*J13L)*PDstandardNth1gt11 + (gtu12*J21L - gtu22*J22L - 
      gtu23*J23L)*PDstandardNth2gt11 + (gtu12*J31L - gtu22*J32L - 
      gtu23*J33L)*PDstandardNth3gt11 + 2*(J11L*(gtu22*PDstandardNth1gt12 + 
      gtu23*PDstandardNth1gt13) + J21L*(gtu22*PDstandardNth2gt12 + 
      gtu23*PDstandardNth2gt13) + J31L*(gtu22*PDstandardNth3gt12 + 
      gtu23*PDstandardNth3gt13)));
    
    CCTK_REAL Gt311 = khalf*((gtu13*J11L - gtu23*J12L - 
      gtu33*J13L)*PDstandardNth1gt11 + (gtu13*J21L - gtu23*J22L - 
      gtu33*J23L)*PDstandardNth2gt11 + (gtu13*J31L - gtu23*J32L - 
      gtu33*J33L)*PDstandardNth3gt11 + 2*(J11L*(gtu23*PDstandardNth1gt12 + 
      gtu33*PDstandardNth1gt13) + J21L*(gtu23*PDstandardNth2gt12 + 
      gtu33*PDstandardNth2gt13) + J31L*(gtu23*PDstandardNth3gt12 + 
      gtu33*PDstandardNth3gt13)));
    
    CCTK_REAL Gt112 = khalf*(gtu11*(J12L*PDstandardNth1gt11 + 
      J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
      gtu12*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
      J31L*PDstandardNth3gt22) + gtu13*(-(J13L*PDstandardNth1gt12) + 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23));
    
    CCTK_REAL Gt212 = khalf*(gtu12*(J12L*PDstandardNth1gt11 + 
      J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
      gtu22*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
      J31L*PDstandardNth3gt22) + gtu23*(-(J13L*PDstandardNth1gt12) + 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23));
    
    CCTK_REAL Gt312 = khalf*(gtu13*(J12L*PDstandardNth1gt11 + 
      J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11) + 
      gtu23*(J11L*PDstandardNth1gt22 + J21L*PDstandardNth2gt22 + 
      J31L*PDstandardNth3gt22) + gtu33*(-(J13L*PDstandardNth1gt12) + 
      J12L*PDstandardNth1gt13 + J11L*PDstandardNth1gt23 - 
      J23L*PDstandardNth2gt12 + J22L*PDstandardNth2gt13 + 
      J21L*PDstandardNth2gt23 - J33L*PDstandardNth3gt12 + 
      J32L*PDstandardNth3gt13 + J31L*PDstandardNth3gt23));
    
    CCTK_REAL Gt113 = khalf*(gtu11*(J13L*PDstandardNth1gt11 + 
      J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
      gtu12*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
      J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
      J31L*PDstandardNth3gt23) + gtu13*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL Gt213 = khalf*(gtu12*(J13L*PDstandardNth1gt11 + 
      J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
      gtu22*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
      J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
      J31L*PDstandardNth3gt23) + gtu23*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL Gt313 = khalf*(gtu13*(J13L*PDstandardNth1gt11 + 
      J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11) + 
      gtu23*(J13L*PDstandardNth1gt12 - J12L*PDstandardNth1gt13 + 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 - 
      J22L*PDstandardNth2gt13 + J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 - J32L*PDstandardNth3gt13 + 
      J31L*PDstandardNth3gt23) + gtu33*(J11L*PDstandardNth1gt33 + 
      J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33));
    
    CCTK_REAL Gt122 = khalf*(gtu11*(-(J11L*PDstandardNth1gt22) + 
      2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
      J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
      J31L*PDstandardNth3gt22) + gtu12*(J12L*PDstandardNth1gt22 + 
      J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
      gtu13*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
      J33L*PDstandardNth3gt22 - 2*(J12L*PDstandardNth1gt23 + 
      J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL Gt222 = khalf*(gtu12*(-(J11L*PDstandardNth1gt22) + 
      2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
      J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
      J31L*PDstandardNth3gt22) + gtu22*(J12L*PDstandardNth1gt22 + 
      J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
      gtu23*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
      J33L*PDstandardNth3gt22 - 2*(J12L*PDstandardNth1gt23 + 
      J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL Gt322 = khalf*(gtu13*(-(J11L*PDstandardNth1gt22) + 
      2*(J12L*PDstandardNth1gt12 + J22L*PDstandardNth2gt12) - 
      J21L*PDstandardNth2gt22 + 2*J32L*PDstandardNth3gt12 - 
      J31L*PDstandardNth3gt22) + gtu23*(J12L*PDstandardNth1gt22 + 
      J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22) - 
      gtu33*(J13L*PDstandardNth1gt22 + J23L*PDstandardNth2gt22 + 
      J33L*PDstandardNth3gt22 - 2*(J12L*PDstandardNth1gt23 + 
      J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23)));
    
    CCTK_REAL Gt123 = khalf*(gtu12*(J13L*PDstandardNth1gt22 + 
      J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
      gtu11*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
      J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt23) + gtu13*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL Gt223 = khalf*(gtu22*(J13L*PDstandardNth1gt22 + 
      J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
      gtu12*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
      J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt23) + gtu23*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL Gt323 = khalf*(gtu23*(J13L*PDstandardNth1gt22 + 
      J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22) + 
      gtu13*(J13L*PDstandardNth1gt12 + J12L*PDstandardNth1gt13 - 
      J11L*PDstandardNth1gt23 + J23L*PDstandardNth2gt12 + 
      J22L*PDstandardNth2gt13 - J21L*PDstandardNth2gt23 + 
      J33L*PDstandardNth3gt12 + J32L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt23) + gtu33*(J12L*PDstandardNth1gt33 + 
      J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33));
    
    CCTK_REAL Gt133 = khalf*(gtu11*(-(J11L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
      J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt33) + gtu12*(-(J12L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
      J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
      J32L*PDstandardNth3gt33) + gtu13*(J13L*PDstandardNth1gt33 + 
      J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL Gt233 = khalf*(gtu12*(-(J11L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
      J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt33) + gtu22*(-(J12L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
      J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
      J32L*PDstandardNth3gt33) + gtu23*(J13L*PDstandardNth1gt33 + 
      J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL Gt333 = khalf*(gtu13*(-(J11L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt13 + J23L*PDstandardNth2gt13) - 
      J21L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt13 - 
      J31L*PDstandardNth3gt33) + gtu23*(-(J12L*PDstandardNth1gt33) + 
      2*(J13L*PDstandardNth1gt23 + J23L*PDstandardNth2gt23) - 
      J22L*PDstandardNth2gt33 + 2*J33L*PDstandardNth3gt23 - 
      J32L*PDstandardNth3gt33) + gtu33*(J13L*PDstandardNth1gt33 + 
      J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33));
    
    CCTK_REAL Xt1L = Gt111*gtu11 + Gt122*gtu22 + 2*(Gt112*gtu12 + 
      Gt113*gtu13 + Gt123*gtu23) + Gt133*gtu33;
    
    CCTK_REAL Xt2L = Gt211*gtu11 + Gt222*gtu22 + 2*(Gt212*gtu12 + 
      Gt213*gtu13 + Gt223*gtu23) + Gt233*gtu33;
    
    CCTK_REAL Xt3L = Gt311*gtu11 + Gt322*gtu22 + 2*(Gt312*gtu12 + 
      Gt313*gtu13 + Gt323*gtu23) + Gt333*gtu33;
    
    CCTK_REAL AL = IfThen(LapseACoeff != 0,(-dtalpL + 
      J11L*LapseAdvectionCoeff*(beta1L*PDupwindNthAnti1alpha + 
      PDupwindNthSymm1alpha*Abs(beta1L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta1L*J21L*PDupwindNthAnti2alpha + 
      J11L*PDupwindNthSymm1alpha*Abs(beta1L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta1L*J31L*PDupwindNthAnti3alpha + 
      J11L*PDupwindNthSymm1alpha*Abs(beta1L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      J21L*LapseAdvectionCoeff*(beta1L*PDupwindNthAnti2alpha + 
      PDupwindNthSymm2alpha*Abs(beta1L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta1L*J11L*PDupwindNthAnti1alpha + 
      J21L*PDupwindNthSymm2alpha*Abs(beta1L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta1L*J31L*PDupwindNthAnti3alpha + 
      J21L*PDupwindNthSymm2alpha*Abs(beta1L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      J31L*LapseAdvectionCoeff*(beta1L*PDupwindNthAnti3alpha + 
      PDupwindNthSymm3alpha*Abs(beta1L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta1L*J11L*PDupwindNthAnti1alpha + 
      J31L*PDupwindNthSymm3alpha*Abs(beta1L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta1L*J21L*PDupwindNthAnti2alpha + 
      J31L*PDupwindNthSymm3alpha*Abs(beta1L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      J12L*LapseAdvectionCoeff*(beta2L*PDupwindNthAnti1alpha + 
      PDupwindNthSymm1alpha*Abs(beta2L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta2L*J22L*PDupwindNthAnti2alpha + 
      J12L*PDupwindNthSymm1alpha*Abs(beta2L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta2L*J32L*PDupwindNthAnti3alpha + 
      J12L*PDupwindNthSymm1alpha*Abs(beta2L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      J22L*LapseAdvectionCoeff*(beta2L*PDupwindNthAnti2alpha + 
      PDupwindNthSymm2alpha*Abs(beta2L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta2L*J12L*PDupwindNthAnti1alpha + 
      J22L*PDupwindNthSymm2alpha*Abs(beta2L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta2L*J32L*PDupwindNthAnti3alpha + 
      J22L*PDupwindNthSymm2alpha*Abs(beta2L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      J32L*LapseAdvectionCoeff*(beta2L*PDupwindNthAnti3alpha + 
      PDupwindNthSymm3alpha*Abs(beta2L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta2L*J12L*PDupwindNthAnti1alpha + 
      J32L*PDupwindNthSymm3alpha*Abs(beta2L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta2L*J22L*PDupwindNthAnti2alpha + 
      J32L*PDupwindNthSymm3alpha*Abs(beta2L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      J13L*LapseAdvectionCoeff*(beta3L*PDupwindNthAnti1alpha + 
      PDupwindNthSymm1alpha*Abs(beta3L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta3L*J23L*PDupwindNthAnti2alpha + 
      J13L*PDupwindNthSymm1alpha*Abs(beta3L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta3L*J33L*PDupwindNthAnti3alpha + 
      J13L*PDupwindNthSymm1alpha*Abs(beta3L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      J23L*LapseAdvectionCoeff*(beta3L*PDupwindNthAnti2alpha + 
      PDupwindNthSymm2alpha*Abs(beta3L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta3L*J13L*PDupwindNthAnti1alpha + 
      J23L*PDupwindNthSymm2alpha*Abs(beta3L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta3L*J33L*PDupwindNthAnti3alpha + 
      J23L*PDupwindNthSymm2alpha*Abs(beta3L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      J33L*LapseAdvectionCoeff*(beta3L*PDupwindNthAnti3alpha + 
      PDupwindNthSymm3alpha*Abs(beta3L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta3L*J13L*PDupwindNthAnti1alpha + 
      J33L*PDupwindNthSymm3alpha*Abs(beta3L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0) 
      + IfThen(LapseACoeff != 0,(-dtalpL + 
      LapseAdvectionCoeff*(beta3L*J23L*PDupwindNthAnti2alpha + 
      J33L*PDupwindNthSymm3alpha*Abs(beta3L)))*INV(harmonicF)*pow(alphaL,-harmonicN),0);
    
    CCTK_REAL theta = fmin(1,exp(1 - 
      rL*INV(SpatialShiftGammaCoeffRadius)));
    
    CCTK_REAL B1L = IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      J11L*ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti1beta1 + 
      PDupwindNthSymm1beta1*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta1L*J21L*PDupwindNthAnti2beta1 + 
      J11L*PDupwindNthSymm1beta1*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta1L*J31L*PDupwindNthAnti3beta1 + 
      J11L*PDupwindNthSymm1beta1*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      J21L*ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti2beta1 + 
      PDupwindNthSymm2beta1*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta1L*J11L*PDupwindNthAnti1beta1 + 
      J21L*PDupwindNthSymm2beta1*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta1L*J31L*PDupwindNthAnti3beta1 + 
      J21L*PDupwindNthSymm2beta1*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      J31L*ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti3beta1 + 
      PDupwindNthSymm3beta1*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta1L*J11L*PDupwindNthAnti1beta1 + 
      J31L*PDupwindNthSymm3beta1*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta1L*J21L*PDupwindNthAnti2beta1 + 
      J31L*PDupwindNthSymm3beta1*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      J12L*ShiftAdvectionCoeff*(beta2L*PDupwindNthAnti1beta1 + 
      PDupwindNthSymm1beta1*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta2L*J22L*PDupwindNthAnti2beta1 + 
      J12L*PDupwindNthSymm1beta1*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta2L*J32L*PDupwindNthAnti3beta1 + 
      J12L*PDupwindNthSymm1beta1*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      J22L*ShiftAdvectionCoeff*(beta2L*PDupwindNthAnti2beta1 + 
      PDupwindNthSymm2beta1*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta2L*J12L*PDupwindNthAnti1beta1 + 
      J22L*PDupwindNthSymm2beta1*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta2L*J32L*PDupwindNthAnti3beta1 + 
      J22L*PDupwindNthSymm2beta1*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      J32L*ShiftAdvectionCoeff*(beta2L*PDupwindNthAnti3beta1 + 
      PDupwindNthSymm3beta1*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta2L*J12L*PDupwindNthAnti1beta1 + 
      J32L*PDupwindNthSymm3beta1*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta2L*J22L*PDupwindNthAnti2beta1 + 
      J32L*PDupwindNthSymm3beta1*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      J13L*ShiftAdvectionCoeff*(beta3L*PDupwindNthAnti1beta1 + 
      PDupwindNthSymm1beta1*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta3L*J23L*PDupwindNthAnti2beta1 + 
      J13L*PDupwindNthSymm1beta1*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta3L*J33L*PDupwindNthAnti3beta1 + 
      J13L*PDupwindNthSymm1beta1*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      J23L*ShiftAdvectionCoeff*(beta3L*PDupwindNthAnti2beta1 + 
      PDupwindNthSymm2beta1*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta3L*J13L*PDupwindNthAnti1beta1 + 
      J23L*PDupwindNthSymm2beta1*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta3L*J33L*PDupwindNthAnti3beta1 + 
      J23L*PDupwindNthSymm2beta1*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      J33L*ShiftAdvectionCoeff*(beta3L*PDupwindNthAnti3beta1 + 
      PDupwindNthSymm3beta1*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta3L*J13L*PDupwindNthAnti1beta1 + 
      J33L*PDupwindNthSymm3beta1*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetaxL - 
      ShiftAdvectionCoeff*(beta3L*J23L*PDupwindNthAnti2beta1 + 
      J33L*PDupwindNthSymm3beta1*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0);
    
    CCTK_REAL B2L = IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      J11L*ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti1beta2 + 
      PDupwindNthSymm1beta2*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta1L*J21L*PDupwindNthAnti2beta2 + 
      J11L*PDupwindNthSymm1beta2*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta1L*J31L*PDupwindNthAnti3beta2 + 
      J11L*PDupwindNthSymm1beta2*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      J21L*ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti2beta2 + 
      PDupwindNthSymm2beta2*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta1L*J11L*PDupwindNthAnti1beta2 + 
      J21L*PDupwindNthSymm2beta2*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta1L*J31L*PDupwindNthAnti3beta2 + 
      J21L*PDupwindNthSymm2beta2*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      J31L*ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti3beta2 + 
      PDupwindNthSymm3beta2*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta1L*J11L*PDupwindNthAnti1beta2 + 
      J31L*PDupwindNthSymm3beta2*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta1L*J21L*PDupwindNthAnti2beta2 + 
      J31L*PDupwindNthSymm3beta2*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      J12L*ShiftAdvectionCoeff*(beta2L*PDupwindNthAnti1beta2 + 
      PDupwindNthSymm1beta2*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta2L*J22L*PDupwindNthAnti2beta2 + 
      J12L*PDupwindNthSymm1beta2*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta2L*J32L*PDupwindNthAnti3beta2 + 
      J12L*PDupwindNthSymm1beta2*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      J22L*ShiftAdvectionCoeff*(beta2L*PDupwindNthAnti2beta2 + 
      PDupwindNthSymm2beta2*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta2L*J12L*PDupwindNthAnti1beta2 + 
      J22L*PDupwindNthSymm2beta2*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta2L*J32L*PDupwindNthAnti3beta2 + 
      J22L*PDupwindNthSymm2beta2*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      J32L*ShiftAdvectionCoeff*(beta2L*PDupwindNthAnti3beta2 + 
      PDupwindNthSymm3beta2*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta2L*J12L*PDupwindNthAnti1beta2 + 
      J32L*PDupwindNthSymm3beta2*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta2L*J22L*PDupwindNthAnti2beta2 + 
      J32L*PDupwindNthSymm3beta2*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      J13L*ShiftAdvectionCoeff*(beta3L*PDupwindNthAnti1beta2 + 
      PDupwindNthSymm1beta2*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta3L*J23L*PDupwindNthAnti2beta2 + 
      J13L*PDupwindNthSymm1beta2*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta3L*J33L*PDupwindNthAnti3beta2 + 
      J13L*PDupwindNthSymm1beta2*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      J23L*ShiftAdvectionCoeff*(beta3L*PDupwindNthAnti2beta2 + 
      PDupwindNthSymm2beta2*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta3L*J13L*PDupwindNthAnti1beta2 + 
      J23L*PDupwindNthSymm2beta2*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta3L*J33L*PDupwindNthAnti3beta2 + 
      J23L*PDupwindNthSymm2beta2*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      J33L*ShiftAdvectionCoeff*(beta3L*PDupwindNthAnti3beta2 + 
      PDupwindNthSymm3beta2*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta3L*J13L*PDupwindNthAnti1beta2 + 
      J33L*PDupwindNthSymm3beta2*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetayL - 
      ShiftAdvectionCoeff*(beta3L*J23L*PDupwindNthAnti2beta2 + 
      J33L*PDupwindNthSymm3beta2*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0);
    
    CCTK_REAL B3L = IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      J11L*ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti1beta3 + 
      PDupwindNthSymm1beta3*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta1L*J21L*PDupwindNthAnti2beta3 + 
      J11L*PDupwindNthSymm1beta3*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta1L*J31L*PDupwindNthAnti3beta3 + 
      J11L*PDupwindNthSymm1beta3*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      J21L*ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti2beta3 + 
      PDupwindNthSymm2beta3*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta1L*J11L*PDupwindNthAnti1beta3 + 
      J21L*PDupwindNthSymm2beta3*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta1L*J31L*PDupwindNthAnti3beta3 + 
      J21L*PDupwindNthSymm2beta3*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      J31L*ShiftAdvectionCoeff*(beta1L*PDupwindNthAnti3beta3 + 
      PDupwindNthSymm3beta3*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta1L*J11L*PDupwindNthAnti1beta3 + 
      J31L*PDupwindNthSymm3beta3*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta1L*J21L*PDupwindNthAnti2beta3 + 
      J31L*PDupwindNthSymm3beta3*Abs(beta1L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      J12L*ShiftAdvectionCoeff*(beta2L*PDupwindNthAnti1beta3 + 
      PDupwindNthSymm1beta3*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta2L*J22L*PDupwindNthAnti2beta3 + 
      J12L*PDupwindNthSymm1beta3*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta2L*J32L*PDupwindNthAnti3beta3 + 
      J12L*PDupwindNthSymm1beta3*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      J22L*ShiftAdvectionCoeff*(beta2L*PDupwindNthAnti2beta3 + 
      PDupwindNthSymm2beta3*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta2L*J12L*PDupwindNthAnti1beta3 + 
      J22L*PDupwindNthSymm2beta3*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta2L*J32L*PDupwindNthAnti3beta3 + 
      J22L*PDupwindNthSymm2beta3*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      J32L*ShiftAdvectionCoeff*(beta2L*PDupwindNthAnti3beta3 + 
      PDupwindNthSymm3beta3*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta2L*J12L*PDupwindNthAnti1beta3 + 
      J32L*PDupwindNthSymm3beta3*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta2L*J22L*PDupwindNthAnti2beta3 + 
      J32L*PDupwindNthSymm3beta3*Abs(beta2L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      J13L*ShiftAdvectionCoeff*(beta3L*PDupwindNthAnti1beta3 + 
      PDupwindNthSymm1beta3*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta3L*J23L*PDupwindNthAnti2beta3 + 
      J13L*PDupwindNthSymm1beta3*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta3L*J33L*PDupwindNthAnti3beta3 + 
      J13L*PDupwindNthSymm1beta3*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      J23L*ShiftAdvectionCoeff*(beta3L*PDupwindNthAnti2beta3 + 
      PDupwindNthSymm2beta3*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta3L*J13L*PDupwindNthAnti1beta3 + 
      J23L*PDupwindNthSymm2beta3*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta3L*J33L*PDupwindNthAnti3beta3 + 
      J23L*PDupwindNthSymm2beta3*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      J33L*ShiftAdvectionCoeff*(beta3L*PDupwindNthAnti3beta3 + 
      PDupwindNthSymm3beta3*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta3L*J13L*PDupwindNthAnti1beta3 + 
      J33L*PDupwindNthSymm3beta3*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0) 
      + IfThen(ShiftBCoeff*ShiftGammaCoeff != 0,(dtbetazL - 
      ShiftAdvectionCoeff*(beta3L*J23L*PDupwindNthAnti2beta3 + 
      J33L*PDupwindNthSymm3beta3*Abs(beta3L)))*INV(ShiftGammaCoeff)*INV(theta),0);
    
    
    /* Copy local copies back to grid functions */
    A[index] = AL;
    B1[index] = B1L;
    B2[index] = B2L;
    B3[index] = B3L;
    Xt1[index] = Xt1L;
    Xt2[index] = Xt2L;
    Xt3[index] = Xt3L;
  }
  LC_ENDLOOP3 (ML_BSSN_MP_convertFromADMBaseGamma);
}

void ML_BSSN_MP_convertFromADMBaseGamma(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_convertFromADMBaseGamma_Body);
}
