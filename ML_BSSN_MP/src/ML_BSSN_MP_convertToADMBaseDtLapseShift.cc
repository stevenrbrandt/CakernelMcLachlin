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

extern "C" void ML_BSSN_MP_convertToADMBaseDtLapseShift_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ADMBase::dtlapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ADMBase::dtshift.");
  return;
}

static void ML_BSSN_MP_convertToADMBaseDtLapseShift_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_convertToADMBaseDtLapseShift_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_convertToADMBaseDtLapseShift_calc_every != ML_BSSN_MP_convertToADMBaseDtLapseShift_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ADMBase::dtlapse","ADMBase::dtshift","Coordinates::jacobian","grid::coordinates","Grid::coordinates","ML_BSSN_MP::ML_dtlapse","ML_BSSN_MP::ML_dtshift","ML_BSSN_MP::ML_Gamma","ML_BSSN_MP::ML_lapse","ML_BSSN_MP::ML_shift","ML_BSSN_MP::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_convertToADMBaseDtLapseShift", 11, groups);
  
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
  CCTK_REAL const p1o12dx = 0.0833333333333333333333333333333*INV(dx);
  CCTK_REAL const p1o12dy = 0.0833333333333333333333333333333*INV(dy);
  CCTK_REAL const p1o12dz = 0.0833333333333333333333333333333*INV(dz);
  CCTK_REAL const p1o144dxdy = 0.00694444444444444444444444444444*INV(dx)*INV(dy);
  CCTK_REAL const p1o144dxdz = 0.00694444444444444444444444444444*INV(dx)*INV(dz);
  CCTK_REAL const p1o144dydz = 0.00694444444444444444444444444444*INV(dy)*INV(dz);
  CCTK_REAL const p1o24dx = 0.0416666666666666666666666666667*INV(dx);
  CCTK_REAL const p1o24dy = 0.0416666666666666666666666666667*INV(dy);
  CCTK_REAL const p1o24dz = 0.0416666666666666666666666666667*INV(dz);
  CCTK_REAL const p1o64dx = 0.015625*INV(dx);
  CCTK_REAL const p1o64dy = 0.015625*INV(dy);
  CCTK_REAL const p1o64dz = 0.015625*INV(dz);
  CCTK_REAL const p1odx = INV(dx);
  CCTK_REAL const p1ody = INV(dy);
  CCTK_REAL const p1odz = INV(dz);
  CCTK_REAL const pm1o12dx2 = -0.0833333333333333333333333333333*INV(SQR(dx));
  CCTK_REAL const pm1o12dy2 = -0.0833333333333333333333333333333*INV(SQR(dy));
  CCTK_REAL const pm1o12dz2 = -0.0833333333333333333333333333333*INV(SQR(dz));
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSN_MP_convertToADMBaseDtLapseShift,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    CCTK_REAL AL = A[index];
    CCTK_REAL alphaL = alpha[index];
    CCTK_REAL B1L = B1[index];
    CCTK_REAL B2L = B2[index];
    CCTK_REAL B3L = B3[index];
    CCTK_REAL beta1L = beta1[index];
    CCTK_REAL beta2L = beta2[index];
    CCTK_REAL beta3L = beta3[index];
    CCTK_REAL J11L = J11[index];
    CCTK_REAL J12L = J12[index];
    CCTK_REAL J13L = J13[index];
    CCTK_REAL J21L = J21[index];
    CCTK_REAL J22L = J22[index];
    CCTK_REAL J23L = J23[index];
    CCTK_REAL J31L = J31[index];
    CCTK_REAL J32L = J32[index];
    CCTK_REAL J33L = J33[index];
    CCTK_REAL rL = r[index];
    CCTK_REAL trKL = trK[index];
    CCTK_REAL Xt1L = Xt1[index];
    CCTK_REAL Xt2L = Xt2[index];
    CCTK_REAL Xt3L = Xt3[index];
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL const PDupwindNthAnti1alpha = PDupwindNthAnti1(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm1alpha = PDupwindNthSymm1(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti2alpha = PDupwindNthAnti2(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm2alpha = PDupwindNthSymm2(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti3alpha = PDupwindNthAnti3(&alpha[index]);
    CCTK_REAL const PDupwindNthSymm3alpha = PDupwindNthSymm3(&alpha[index]);
    CCTK_REAL const PDupwindNthAnti1beta1 = PDupwindNthAnti1(&beta1[index]);
    CCTK_REAL const PDupwindNthSymm1beta1 = PDupwindNthSymm1(&beta1[index]);
    CCTK_REAL const PDupwindNthAnti2beta1 = PDupwindNthAnti2(&beta1[index]);
    CCTK_REAL const PDupwindNthSymm2beta1 = PDupwindNthSymm2(&beta1[index]);
    CCTK_REAL const PDupwindNthAnti3beta1 = PDupwindNthAnti3(&beta1[index]);
    CCTK_REAL const PDupwindNthSymm3beta1 = PDupwindNthSymm3(&beta1[index]);
    CCTK_REAL const PDupwindNthAnti1beta2 = PDupwindNthAnti1(&beta2[index]);
    CCTK_REAL const PDupwindNthSymm1beta2 = PDupwindNthSymm1(&beta2[index]);
    CCTK_REAL const PDupwindNthAnti2beta2 = PDupwindNthAnti2(&beta2[index]);
    CCTK_REAL const PDupwindNthSymm2beta2 = PDupwindNthSymm2(&beta2[index]);
    CCTK_REAL const PDupwindNthAnti3beta2 = PDupwindNthAnti3(&beta2[index]);
    CCTK_REAL const PDupwindNthSymm3beta2 = PDupwindNthSymm3(&beta2[index]);
    CCTK_REAL const PDupwindNthAnti1beta3 = PDupwindNthAnti1(&beta3[index]);
    CCTK_REAL const PDupwindNthSymm1beta3 = PDupwindNthSymm1(&beta3[index]);
    CCTK_REAL const PDupwindNthAnti2beta3 = PDupwindNthAnti2(&beta3[index]);
    CCTK_REAL const PDupwindNthSymm2beta3 = PDupwindNthSymm2(&beta3[index]);
    CCTK_REAL const PDupwindNthAnti3beta3 = PDupwindNthAnti3(&beta3[index]);
    CCTK_REAL const PDupwindNthSymm3beta3 = PDupwindNthSymm3(&beta3[index]);
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    CCTK_REAL eta = fmin(1,INV(rL)*ToReal(SpatialBetaDriverRadius));
    
    CCTK_REAL theta = fmin(1,exp(1 - 
      rL*INV(ToReal(SpatialShiftGammaCoeffRadius))));
    
    CCTK_REAL dtalpL = 
      -(pow(alphaL,ToReal(harmonicN))*ToReal(harmonicF)*(trKL + (AL - 
      trKL)*ToReal(LapseACoeff))) + ((beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1alpha + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2alpha + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3alpha + (J11L*PDupwindNthSymm1alpha + 
      J21L*PDupwindNthSymm2alpha + J31L*PDupwindNthSymm3alpha)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1alpha + J22L*PDupwindNthSymm2alpha + 
      J32L*PDupwindNthSymm3alpha)*Abs(beta2L) + (J13L*PDupwindNthSymm1alpha + 
      J23L*PDupwindNthSymm2alpha + 
      J33L*PDupwindNthSymm3alpha)*Abs(beta3L))*ToReal(LapseAdvectionCoeff);
    
    CCTK_REAL dtbetaxL = ((beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1beta1 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2beta1 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3beta1 + (J11L*PDupwindNthSymm1beta1 + 
      J21L*PDupwindNthSymm2beta1 + J31L*PDupwindNthSymm3beta1)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1beta1 + J22L*PDupwindNthSymm2beta1 + 
      J32L*PDupwindNthSymm3beta1)*Abs(beta2L) + (J13L*PDupwindNthSymm1beta1 + 
      J23L*PDupwindNthSymm2beta1 + 
      J33L*PDupwindNthSymm3beta1)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      theta*(Xt1L + beta1L*eta*ToReal(BetaDriver)*(-1 + ToReal(ShiftBCoeff)) 
      + (B1L - Xt1L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    CCTK_REAL dtbetayL = ((beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1beta2 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2beta2 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3beta2 + (J11L*PDupwindNthSymm1beta2 + 
      J21L*PDupwindNthSymm2beta2 + J31L*PDupwindNthSymm3beta2)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1beta2 + J22L*PDupwindNthSymm2beta2 + 
      J32L*PDupwindNthSymm3beta2)*Abs(beta2L) + (J13L*PDupwindNthSymm1beta2 + 
      J23L*PDupwindNthSymm2beta2 + 
      J33L*PDupwindNthSymm3beta2)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      theta*(Xt2L + beta2L*eta*ToReal(BetaDriver)*(-1 + ToReal(ShiftBCoeff)) 
      + (B2L - Xt2L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    CCTK_REAL dtbetazL = ((beta1L*J11L + beta2L*J12L + 
      beta3L*J13L)*PDupwindNthAnti1beta3 + (beta1L*J21L + beta2L*J22L + 
      beta3L*J23L)*PDupwindNthAnti2beta3 + (beta1L*J31L + beta2L*J32L + 
      beta3L*J33L)*PDupwindNthAnti3beta3 + (J11L*PDupwindNthSymm1beta3 + 
      J21L*PDupwindNthSymm2beta3 + J31L*PDupwindNthSymm3beta3)*Abs(beta1L) + 
      (J12L*PDupwindNthSymm1beta3 + J22L*PDupwindNthSymm2beta3 + 
      J32L*PDupwindNthSymm3beta3)*Abs(beta2L) + (J13L*PDupwindNthSymm1beta3 + 
      J23L*PDupwindNthSymm2beta3 + 
      J33L*PDupwindNthSymm3beta3)*Abs(beta3L))*ToReal(ShiftAdvectionCoeff) + 
      theta*(Xt3L + beta3L*eta*ToReal(BetaDriver)*(-1 + ToReal(ShiftBCoeff)) 
      + (B3L - Xt3L)*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    
    
    /* Copy local copies back to grid functions */
    dtalp[index] = dtalpL;
    dtbetax[index] = dtbetaxL;
    dtbetay[index] = dtbetayL;
    dtbetaz[index] = dtbetazL;
  }
  LC_ENDLOOP3 (ML_BSSN_MP_convertToADMBaseDtLapseShift);
}

extern "C" void ML_BSSN_MP_convertToADMBaseDtLapseShift(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_convertToADMBaseDtLapseShift_Body);
}
