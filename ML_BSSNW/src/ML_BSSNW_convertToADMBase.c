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

void ML_BSSNW_convertToADMBase_Body(cGH const * const cctkGH, CCTK_INT const dir, CCTK_INT const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], CCTK_INT const min[3], CCTK_INT const max[3], CCTK_INT const n_subblock_gfs, CCTK_REAL * const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  CCTK_REAL dx = INITVALUE, dy = INITVALUE, dz = INITVALUE;
  CCTK_REAL dxi = INITVALUE, dyi = INITVALUE, dzi = INITVALUE;
  CCTK_REAL khalf = INITVALUE, kthird = INITVALUE, ktwothird = INITVALUE, kfourthird = INITVALUE, keightthird = INITVALUE;
  CCTK_REAL hdxi = INITVALUE, hdyi = INITVALUE, hdzi = INITVALUE;
  
  
  /* Declare predefined quantities */
  CCTK_REAL p1o12dx = INITVALUE;
  CCTK_REAL p1o12dy = INITVALUE;
  CCTK_REAL p1o12dz = INITVALUE;
  CCTK_REAL p1o144dxdy = INITVALUE;
  CCTK_REAL p1o144dxdz = INITVALUE;
  CCTK_REAL p1o144dydz = INITVALUE;
  CCTK_REAL pm1o12dx2 = INITVALUE;
  CCTK_REAL pm1o12dy2 = INITVALUE;
  CCTK_REAL pm1o12dz2 = INITVALUE;
  CCTK_REAL Differencing`Private`liName$23747 = INITVALUE;
  CCTK_REAL Differencing`Private`liName$23763 = INITVALUE;
  CCTK_REAL Differencing`Private`liName$23779 = INITVALUE;
  CCTK_REAL Differencing`Private`liName$23795 = INITVALUE;
  CCTK_REAL Differencing`Private`liName$23811 = INITVALUE;
  CCTK_REAL Differencing`Private`liName$23827 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSNW_convertToADMBase_Body");
  }
  
  if (cctk_iteration % ML_BSSNW_convertToADMBase_calc_every != ML_BSSNW_convertToADMBase_calc_offset)
  {
    return;
  }
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  dx = CCTK_DELTA_SPACE(0);
  dy = CCTK_DELTA_SPACE(1);
  dz = CCTK_DELTA_SPACE(2);
  dxi = 1.0 / dx;
  dyi = 1.0 / dy;
  dzi = 1.0 / dz;
  khalf = 0.5;
  kthird = 1/3.0;
  ktwothird = 2.0/3.0;
  kfourthird = 4.0/3.0;
  keightthird = 8.0/3.0;
  hdxi = 0.5 * dxi;
  hdyi = 0.5 * dyi;
  hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  p1o12dx = INV(dx)/12.;
  p1o12dy = INV(dy)/12.;
  p1o12dz = INV(dz)/12.;
  p1o144dxdy = (INV(dx)*INV(dy))/144.;
  p1o144dxdz = (INV(dx)*INV(dz))/144.;
  p1o144dydz = (INV(dy)*INV(dz))/144.;
  pm1o12dx2 = -pow(dx,-2)/12.;
  pm1o12dy2 = -pow(dy,-2)/12.;
  pm1o12dz2 = -pow(dz,-2)/12.;
  Differencing`Private`liName$23747 = Differencing_Private_num$23747*Differencing_Private_ss$23747*INV(Differencing_Private_den$23747);
  Differencing`Private`liName$23763 = Differencing_Private_num$23763*Differencing_Private_ss$23763*INV(Differencing_Private_den$23763);
  Differencing`Private`liName$23779 = Differencing_Private_num$23779*Differencing_Private_ss$23779*INV(Differencing_Private_den$23779);
  Differencing`Private`liName$23795 = Differencing_Private_num$23795*Differencing_Private_ss$23795*INV(Differencing_Private_den$23795);
  Differencing`Private`liName$23811 = Differencing_Private_num$23811*Differencing_Private_ss$23811*INV(Differencing_Private_den$23811);
  Differencing`Private`liName$23827 = Differencing_Private_num$23827*Differencing_Private_ss$23827*INV(Differencing_Private_den$23827);
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3 (ML_BSSNW_convertToADMBase,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lssh[CCTK_LSSH_IDX(0,0)],cctk_lssh[CCTK_LSSH_IDX(0,1)],cctk_lssh[CCTK_LSSH_IDX(0,2)])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL betam1 = INITVALUE, betam2 = INITVALUE, betam3 = INITVALUE;
    CCTK_REAL betap1 = INITVALUE, betap2 = INITVALUE, betap3 = INITVALUE;
    CCTK_REAL g11 = INITVALUE, g12 = INITVALUE, g13 = INITVALUE, g22 = INITVALUE, g23 = INITVALUE, g33 = INITVALUE;
    CCTK_REAL invW2 = INITVALUE;
    CCTK_REAL K11 = INITVALUE, K12 = INITVALUE, K13 = INITVALUE, K22 = INITVALUE, K23 = INITVALUE, K33 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL AL = INITVALUE;
    CCTK_REAL alpL = INITVALUE;
    CCTK_REAL alphaL = INITVALUE;
    CCTK_REAL At11L = INITVALUE, At12L = INITVALUE, At13L = INITVALUE, At22L = INITVALUE, At23L = INITVALUE, At33L = INITVALUE;
    CCTK_REAL B1L = INITVALUE, B2L = INITVALUE, B3L = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta2L = INITVALUE, beta3L = INITVALUE;
    CCTK_REAL betaxL = INITVALUE;
    CCTK_REAL betayL = INITVALUE;
    CCTK_REAL betazL = INITVALUE;
    CCTK_REAL dtalpL = INITVALUE;
    CCTK_REAL dtbetaxL = INITVALUE;
    CCTK_REAL dtbetayL = INITVALUE;
    CCTK_REAL dtbetazL = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt12L = INITVALUE, gt13L = INITVALUE, gt22L = INITVALUE, gt23L = INITVALUE, gt33L = INITVALUE;
    CCTK_REAL gxxL = INITVALUE;
    CCTK_REAL gxyL = INITVALUE;
    CCTK_REAL gxzL = INITVALUE;
    CCTK_REAL gyyL = INITVALUE;
    CCTK_REAL gyzL = INITVALUE;
    CCTK_REAL gzzL = INITVALUE;
    CCTK_REAL kxxL = INITVALUE;
    CCTK_REAL kxyL = INITVALUE;
    CCTK_REAL kxzL = INITVALUE;
    CCTK_REAL kyyL = INITVALUE;
    CCTK_REAL kyzL = INITVALUE;
    CCTK_REAL kzzL = INITVALUE;
    CCTK_REAL trKL = INITVALUE;
    CCTK_REAL WL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDupwindpNth1alpha = INITVALUE;
    CCTK_REAL PDupwindpNth2alpha = INITVALUE;
    CCTK_REAL PDupwindpNth3alpha = INITVALUE;
    CCTK_REAL PDupwindmNth1alpha = INITVALUE;
    CCTK_REAL PDupwindmNth2alpha = INITVALUE;
    CCTK_REAL PDupwindmNth3alpha = INITVALUE;
    CCTK_REAL PDupwindpNth1beta1 = INITVALUE;
    CCTK_REAL PDupwindpNth2beta1 = INITVALUE;
    CCTK_REAL PDupwindpNth3beta1 = INITVALUE;
    CCTK_REAL PDupwindmNth1beta1 = INITVALUE;
    CCTK_REAL PDupwindmNth2beta1 = INITVALUE;
    CCTK_REAL PDupwindmNth3beta1 = INITVALUE;
    
    /* Assign local copies of grid functions */
    AL = A[index];
    alphaL = alpha[index];
    At11L = At11[index];
    At12L = At12[index];
    At13L = At13[index];
    At22L = At22[index];
    At23L = At23[index];
    At33L = At33[index];
    B1L = B1[index];
    B2L = B2[index];
    B3L = B3[index];
    beta1L = beta1[index];
    beta2L = beta2[index];
    beta3L = beta3[index];
    gt11L = gt11[index];
    gt12L = gt12[index];
    gt13L = gt13[index];
    gt22L = gt22[index];
    gt23L = gt23[index];
    gt33L = gt33[index];
    trKL = trK[index];
    WL = W[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDupwindpNth1alpha = PDupwindpNth1(alpha, i, j, k);
    PDupwindpNth2alpha = PDupwindpNth2(alpha, i, j, k);
    PDupwindpNth3alpha = PDupwindpNth3(alpha, i, j, k);
    PDupwindmNth1alpha = PDupwindmNth1(alpha, i, j, k);
    PDupwindmNth2alpha = PDupwindmNth2(alpha, i, j, k);
    PDupwindmNth3alpha = PDupwindmNth3(alpha, i, j, k);
    PDupwindpNth1beta1 = PDupwindpNth1(beta1, i, j, k);
    PDupwindpNth2beta1 = PDupwindpNth2(beta1, i, j, k);
    PDupwindpNth3beta1 = PDupwindpNth3(beta1, i, j, k);
    PDupwindmNth1beta1 = PDupwindmNth1(beta1, i, j, k);
    PDupwindmNth2beta1 = PDupwindmNth2(beta1, i, j, k);
    PDupwindmNth3beta1 = PDupwindmNth3(beta1, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    betam1  =  khalf*(beta1L - Abs(beta1L));
    
    betam2  =  khalf*(beta2L - Abs(beta2L));
    
    betam3  =  khalf*(beta3L - Abs(beta3L));
    
    betap1  =  khalf*(beta1L + Abs(beta1L));
    
    betap2  =  khalf*(beta2L + Abs(beta2L));
    
    betap3  =  khalf*(beta3L + Abs(beta3L));
    
    invW2  =  pow(WL,-2);
    
    g11  =  gt11L*invW2;
    
    g12  =  gt12L*invW2;
    
    g13  =  gt13L*invW2;
    
    g22  =  gt22L*invW2;
    
    g23  =  gt23L*invW2;
    
    g33  =  gt33L*invW2;
    
    gxxL  =  g11;
    
    gxyL  =  g12;
    
    gxzL  =  g13;
    
    gyyL  =  g22;
    
    gyzL  =  g23;
    
    gzzL  =  g33;
    
    K11  =  At11L*invW2 + g11*kthird*trKL;
    
    K12  =  At12L*invW2 + g12*kthird*trKL;
    
    K13  =  At13L*invW2 + g13*kthird*trKL;
    
    K22  =  At22L*invW2 + g22*kthird*trKL;
    
    K23  =  At23L*invW2 + g23*kthird*trKL;
    
    K33  =  At33L*invW2 + g33*kthird*trKL;
    
    kxxL  =  K11;
    
    kxyL  =  K12;
    
    kxzL  =  K13;
    
    kyyL  =  K22;
    
    kyzL  =  K23;
    
    kzzL  =  K33;
    
    alpL  =  alphaL;
    
    betaxL  =  beta1L;
    
    betayL  =  beta2L;
    
    betazL  =  beta3L;
    
    dtalpL  =  LapseAdvectionCoeff*(betam1*PDupwindmNth1alpha + betam2*PDupwindmNth2alpha + betam3*PDupwindmNth3alpha + 
           betap1*PDupwindpNth1alpha + betap2*PDupwindpNth2alpha + betap3*PDupwindpNth3alpha) + 
        harmonicF*(AL*(-1 + LapseAdvectionCoeff) - LapseAdvectionCoeff*trKL)*pow(alphaL,harmonicN);
    
    dtbetaxL  =  (betam1*PDupwindmNth1beta1 + betam2*PDupwindmNth2beta1 + betam3*PDupwindmNth3beta1 + 
           betap1*PDupwindpNth1beta1 + betap2*PDupwindpNth2beta1 + betap3*PDupwindpNth3beta1)*ShiftAdvectionCoeff + 
        B1L*ShiftGammaCoeff;
    
    dtbetayL  =  (betam1*PDupwindmNth1beta1 + betam2*PDupwindmNth2beta1 + betam3*PDupwindmNth3beta1 + 
           betap1*PDupwindpNth1beta1 + betap2*PDupwindpNth2beta1 + betap3*PDupwindpNth3beta1)*ShiftAdvectionCoeff + 
        B2L*ShiftGammaCoeff;
    
    dtbetazL  =  (betam1*PDupwindmNth1beta1 + betam2*PDupwindmNth2beta1 + betam3*PDupwindmNth3beta1 + 
           betap1*PDupwindpNth1beta1 + betap2*PDupwindpNth2beta1 + betap3*PDupwindpNth3beta1)*ShiftAdvectionCoeff + 
        B3L*ShiftGammaCoeff;
    
    
    /* Copy local copies back to grid functions */
    alp[index] = alpL;
    betax[index] = betaxL;
    betay[index] = betayL;
    betaz[index] = betazL;
    dtalp[index] = dtalpL;
    dtbetax[index] = dtbetaxL;
    dtbetay[index] = dtbetayL;
    dtbetaz[index] = dtbetazL;
    gxx[index] = gxxL;
    gxy[index] = gxyL;
    gxz[index] = gxzL;
    gyy[index] = gyyL;
    gyz[index] = gyzL;
    gzz[index] = gzzL;
    kxx[index] = kxxL;
    kxy[index] = kxyL;
    kxz[index] = kxzL;
    kyy[index] = kyyL;
    kyz[index] = kyzL;
    kzz[index] = kzzL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (ML_BSSNW_convertToADMBase);
}

void ML_BSSNW_convertToADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSNW_convertToADMBase_Body);
}
