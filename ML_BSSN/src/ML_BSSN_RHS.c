/*  File produced by user eschnett */
/*  Produced with Mathematica Version 6.0 for Mac OS X x86 (32-bit) (April 20, 2007) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#define KRANC_C

#include <math.h>
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

void ML_BSSN_RHS_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  
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
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_RHS_Body");
  }
  
  if (cctk_iteration % ML_BSSN_RHS_calc_every != ML_BSSN_RHS_calc_offset)
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
  
  /* Loop over the grid points */
  LC_LOOP3 (somename,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL Atm11 = INITVALUE, Atm12 = INITVALUE, Atm13 = INITVALUE, Atm21 = INITVALUE, Atm22 = INITVALUE, Atm23 = INITVALUE;
    CCTK_REAL Atm31 = INITVALUE, Atm32 = INITVALUE, Atm33 = INITVALUE;
    CCTK_REAL Atu11 = INITVALUE, Atu12 = INITVALUE, Atu13 = INITVALUE, Atu22 = INITVALUE, Atu23 = INITVALUE, Atu33 = INITVALUE;
    CCTK_REAL ddetg1 = INITVALUE, ddetg2 = INITVALUE, ddetg3 = INITVALUE;
    CCTK_REAL ddetgt1 = INITVALUE, ddetgt2 = INITVALUE, ddetgt3 = INITVALUE;
    CCTK_REAL ddgtu1111 = INITVALUE, ddgtu1121 = INITVALUE, ddgtu1131 = INITVALUE, ddgtu1211 = INITVALUE, ddgtu1221 = INITVALUE, ddgtu1222 = INITVALUE;
    CCTK_REAL ddgtu1231 = INITVALUE, ddgtu1232 = INITVALUE, ddgtu1311 = INITVALUE, ddgtu1321 = INITVALUE, ddgtu1331 = INITVALUE, ddgtu1332 = INITVALUE;
    CCTK_REAL ddgtu1333 = INITVALUE, ddgtu2221 = INITVALUE, ddgtu2222 = INITVALUE, ddgtu2232 = INITVALUE, ddgtu2321 = INITVALUE, ddgtu2322 = INITVALUE;
    CCTK_REAL ddgtu2331 = INITVALUE, ddgtu2332 = INITVALUE, ddgtu2333 = INITVALUE, ddgtu3331 = INITVALUE, ddgtu3332 = INITVALUE, ddgtu3333 = INITVALUE;
    CCTK_REAL detg = INITVALUE;
    CCTK_REAL detgt = INITVALUE;
    CCTK_REAL dgtu111 = INITVALUE, dgtu112 = INITVALUE, dgtu113 = INITVALUE, dgtu121 = INITVALUE, dgtu122 = INITVALUE, dgtu123 = INITVALUE;
    CCTK_REAL dgtu131 = INITVALUE, dgtu132 = INITVALUE, dgtu133 = INITVALUE, dgtu221 = INITVALUE, dgtu222 = INITVALUE, dgtu223 = INITVALUE;
    CCTK_REAL dgtu231 = INITVALUE, dgtu232 = INITVALUE, dgtu233 = INITVALUE, dgtu331 = INITVALUE, dgtu332 = INITVALUE, dgtu333 = INITVALUE;
    CCTK_REAL e4phi = INITVALUE;
    CCTK_REAL em4phi = INITVALUE;
    CCTK_REAL g11 = INITVALUE;
    CCTK_REAL G111 = INITVALUE, G121 = INITVALUE, G122 = INITVALUE, G131 = INITVALUE, G132 = INITVALUE, G133 = INITVALUE;
    CCTK_REAL g21 = INITVALUE;
    CCTK_REAL G211 = INITVALUE;
    CCTK_REAL g22 = INITVALUE;
    CCTK_REAL G221 = INITVALUE, G222 = INITVALUE, G231 = INITVALUE, G232 = INITVALUE, G233 = INITVALUE;
    CCTK_REAL g31 = INITVALUE;
    CCTK_REAL G311 = INITVALUE;
    CCTK_REAL g32 = INITVALUE;
    CCTK_REAL G321 = INITVALUE, G322 = INITVALUE;
    CCTK_REAL g33 = INITVALUE;
    CCTK_REAL G331 = INITVALUE, G332 = INITVALUE, G333 = INITVALUE;
    CCTK_REAL Gt111 = INITVALUE, Gt121 = INITVALUE, Gt122 = INITVALUE, Gt131 = INITVALUE, Gt132 = INITVALUE, Gt133 = INITVALUE;
    CCTK_REAL Gt211 = INITVALUE, Gt221 = INITVALUE, Gt222 = INITVALUE, Gt231 = INITVALUE, Gt232 = INITVALUE, Gt233 = INITVALUE;
    CCTK_REAL Gt311 = INITVALUE, Gt321 = INITVALUE, Gt322 = INITVALUE, Gt331 = INITVALUE, Gt332 = INITVALUE, Gt333 = INITVALUE;
    CCTK_REAL gtu11 = INITVALUE, gtu12 = INITVALUE, gtu13 = INITVALUE, gtu22 = INITVALUE, gtu23 = INITVALUE, gtu33 = INITVALUE;
    CCTK_REAL gu11 = INITVALUE, gu12 = INITVALUE, gu13 = INITVALUE, gu22 = INITVALUE, gu23 = INITVALUE, gu33 = INITVALUE;
    CCTK_REAL R11 = INITVALUE, R21 = INITVALUE, R22 = INITVALUE, R31 = INITVALUE, R32 = INITVALUE, R33 = INITVALUE;
    CCTK_REAL Rphi11 = INITVALUE, Rphi21 = INITVALUE, Rphi22 = INITVALUE, Rphi31 = INITVALUE, Rphi32 = INITVALUE, Rphi33 = INITVALUE;
    CCTK_REAL Rt11 = INITVALUE, Rt21 = INITVALUE, Rt22 = INITVALUE, Rt31 = INITVALUE, Rt32 = INITVALUE, Rt33 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL alphaL = INITVALUE, alpharhsL = INITVALUE;
    CCTK_REAL At11L = INITVALUE, At11rhsL = INITVALUE, At21L = INITVALUE, At21rhsL = INITVALUE, At22L = INITVALUE, At22rhsL = INITVALUE;
    CCTK_REAL At31L = INITVALUE, At31rhsL = INITVALUE, At32L = INITVALUE, At32rhsL = INITVALUE, At33L = INITVALUE, At33rhsL = INITVALUE;
    CCTK_REAL beta1L = INITVALUE, beta1rhsL = INITVALUE, beta2L = INITVALUE, beta2rhsL = INITVALUE, beta3L = INITVALUE, beta3rhsL = INITVALUE;
    CCTK_REAL dtalphaL = INITVALUE, dtalpharhsL = INITVALUE;
    CCTK_REAL dtbeta1L = INITVALUE, dtbeta1rhsL = INITVALUE, dtbeta2L = INITVALUE, dtbeta2rhsL = INITVALUE, dtbeta3L = INITVALUE, dtbeta3rhsL = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt11rhsL = INITVALUE, gt21L = INITVALUE, gt21rhsL = INITVALUE, gt22L = INITVALUE, gt22rhsL = INITVALUE;
    CCTK_REAL gt31L = INITVALUE, gt31rhsL = INITVALUE, gt32L = INITVALUE, gt32rhsL = INITVALUE, gt33L = INITVALUE, gt33rhsL = INITVALUE;
    CCTK_REAL phiL = INITVALUE, phirhsL = INITVALUE;
    CCTK_REAL trKL = INITVALUE, trKrhsL = INITVALUE;
    CCTK_REAL Xt1L = INITVALUE, Xt1rhsL = INITVALUE, Xt2L = INITVALUE, Xt2rhsL = INITVALUE, Xt3L = INITVALUE, Xt3rhsL = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandard4th1alpha = INITVALUE;
    CCTK_REAL PDstandard4th2alpha = INITVALUE;
    CCTK_REAL PDstandard4th3alpha = INITVALUE;
    CCTK_REAL PDstandard4th11alpha = INITVALUE;
    CCTK_REAL PDstandard4th22alpha = INITVALUE;
    CCTK_REAL PDstandard4th33alpha = INITVALUE;
    CCTK_REAL PDstandard4th12alpha = INITVALUE;
    CCTK_REAL PDstandard4th13alpha = INITVALUE;
    CCTK_REAL PDstandard4th21alpha = INITVALUE;
    CCTK_REAL PDstandard4th23alpha = INITVALUE;
    CCTK_REAL PDstandard4th31alpha = INITVALUE;
    CCTK_REAL PDstandard4th32alpha = INITVALUE;
    CCTK_REAL PDstandard4th1At11 = INITVALUE;
    CCTK_REAL PDstandard4th2At11 = INITVALUE;
    CCTK_REAL PDstandard4th3At11 = INITVALUE;
    CCTK_REAL PDstandard4th1At21 = INITVALUE;
    CCTK_REAL PDstandard4th2At21 = INITVALUE;
    CCTK_REAL PDstandard4th3At21 = INITVALUE;
    CCTK_REAL PDstandard4th1At22 = INITVALUE;
    CCTK_REAL PDstandard4th2At22 = INITVALUE;
    CCTK_REAL PDstandard4th3At22 = INITVALUE;
    CCTK_REAL PDstandard4th1At31 = INITVALUE;
    CCTK_REAL PDstandard4th2At31 = INITVALUE;
    CCTK_REAL PDstandard4th3At31 = INITVALUE;
    CCTK_REAL PDstandard4th1At32 = INITVALUE;
    CCTK_REAL PDstandard4th2At32 = INITVALUE;
    CCTK_REAL PDstandard4th3At32 = INITVALUE;
    CCTK_REAL PDstandard4th1At33 = INITVALUE;
    CCTK_REAL PDstandard4th2At33 = INITVALUE;
    CCTK_REAL PDstandard4th3At33 = INITVALUE;
    CCTK_REAL PDstandard4th1beta1 = INITVALUE;
    CCTK_REAL PDstandard4th2beta1 = INITVALUE;
    CCTK_REAL PDstandard4th3beta1 = INITVALUE;
    CCTK_REAL PDstandard4th11beta1 = INITVALUE;
    CCTK_REAL PDstandard4th22beta1 = INITVALUE;
    CCTK_REAL PDstandard4th33beta1 = INITVALUE;
    CCTK_REAL PDstandard4th12beta1 = INITVALUE;
    CCTK_REAL PDstandard4th13beta1 = INITVALUE;
    CCTK_REAL PDstandard4th21beta1 = INITVALUE;
    CCTK_REAL PDstandard4th23beta1 = INITVALUE;
    CCTK_REAL PDstandard4th31beta1 = INITVALUE;
    CCTK_REAL PDstandard4th32beta1 = INITVALUE;
    CCTK_REAL PDstandard4th1beta2 = INITVALUE;
    CCTK_REAL PDstandard4th2beta2 = INITVALUE;
    CCTK_REAL PDstandard4th3beta2 = INITVALUE;
    CCTK_REAL PDstandard4th11beta2 = INITVALUE;
    CCTK_REAL PDstandard4th22beta2 = INITVALUE;
    CCTK_REAL PDstandard4th33beta2 = INITVALUE;
    CCTK_REAL PDstandard4th12beta2 = INITVALUE;
    CCTK_REAL PDstandard4th13beta2 = INITVALUE;
    CCTK_REAL PDstandard4th21beta2 = INITVALUE;
    CCTK_REAL PDstandard4th23beta2 = INITVALUE;
    CCTK_REAL PDstandard4th31beta2 = INITVALUE;
    CCTK_REAL PDstandard4th32beta2 = INITVALUE;
    CCTK_REAL PDstandard4th1beta3 = INITVALUE;
    CCTK_REAL PDstandard4th2beta3 = INITVALUE;
    CCTK_REAL PDstandard4th3beta3 = INITVALUE;
    CCTK_REAL PDstandard4th11beta3 = INITVALUE;
    CCTK_REAL PDstandard4th22beta3 = INITVALUE;
    CCTK_REAL PDstandard4th33beta3 = INITVALUE;
    CCTK_REAL PDstandard4th12beta3 = INITVALUE;
    CCTK_REAL PDstandard4th13beta3 = INITVALUE;
    CCTK_REAL PDstandard4th21beta3 = INITVALUE;
    CCTK_REAL PDstandard4th23beta3 = INITVALUE;
    CCTK_REAL PDstandard4th31beta3 = INITVALUE;
    CCTK_REAL PDstandard4th32beta3 = INITVALUE;
    CCTK_REAL PDstandard4th1gt11 = INITVALUE;
    CCTK_REAL PDstandard4th2gt11 = INITVALUE;
    CCTK_REAL PDstandard4th3gt11 = INITVALUE;
    CCTK_REAL PDstandard4th11gt11 = INITVALUE;
    CCTK_REAL PDstandard4th22gt11 = INITVALUE;
    CCTK_REAL PDstandard4th33gt11 = INITVALUE;
    CCTK_REAL PDstandard4th12gt11 = INITVALUE;
    CCTK_REAL PDstandard4th13gt11 = INITVALUE;
    CCTK_REAL PDstandard4th21gt11 = INITVALUE;
    CCTK_REAL PDstandard4th23gt11 = INITVALUE;
    CCTK_REAL PDstandard4th31gt11 = INITVALUE;
    CCTK_REAL PDstandard4th32gt11 = INITVALUE;
    CCTK_REAL PDstandard4th1gt21 = INITVALUE;
    CCTK_REAL PDstandard4th2gt21 = INITVALUE;
    CCTK_REAL PDstandard4th3gt21 = INITVALUE;
    CCTK_REAL PDstandard4th11gt21 = INITVALUE;
    CCTK_REAL PDstandard4th22gt21 = INITVALUE;
    CCTK_REAL PDstandard4th33gt21 = INITVALUE;
    CCTK_REAL PDstandard4th12gt21 = INITVALUE;
    CCTK_REAL PDstandard4th13gt21 = INITVALUE;
    CCTK_REAL PDstandard4th21gt21 = INITVALUE;
    CCTK_REAL PDstandard4th23gt21 = INITVALUE;
    CCTK_REAL PDstandard4th31gt21 = INITVALUE;
    CCTK_REAL PDstandard4th32gt21 = INITVALUE;
    CCTK_REAL PDstandard4th1gt22 = INITVALUE;
    CCTK_REAL PDstandard4th2gt22 = INITVALUE;
    CCTK_REAL PDstandard4th3gt22 = INITVALUE;
    CCTK_REAL PDstandard4th11gt22 = INITVALUE;
    CCTK_REAL PDstandard4th22gt22 = INITVALUE;
    CCTK_REAL PDstandard4th33gt22 = INITVALUE;
    CCTK_REAL PDstandard4th12gt22 = INITVALUE;
    CCTK_REAL PDstandard4th13gt22 = INITVALUE;
    CCTK_REAL PDstandard4th21gt22 = INITVALUE;
    CCTK_REAL PDstandard4th23gt22 = INITVALUE;
    CCTK_REAL PDstandard4th31gt22 = INITVALUE;
    CCTK_REAL PDstandard4th32gt22 = INITVALUE;
    CCTK_REAL PDstandard4th1gt31 = INITVALUE;
    CCTK_REAL PDstandard4th2gt31 = INITVALUE;
    CCTK_REAL PDstandard4th3gt31 = INITVALUE;
    CCTK_REAL PDstandard4th11gt31 = INITVALUE;
    CCTK_REAL PDstandard4th22gt31 = INITVALUE;
    CCTK_REAL PDstandard4th33gt31 = INITVALUE;
    CCTK_REAL PDstandard4th12gt31 = INITVALUE;
    CCTK_REAL PDstandard4th13gt31 = INITVALUE;
    CCTK_REAL PDstandard4th21gt31 = INITVALUE;
    CCTK_REAL PDstandard4th23gt31 = INITVALUE;
    CCTK_REAL PDstandard4th31gt31 = INITVALUE;
    CCTK_REAL PDstandard4th32gt31 = INITVALUE;
    CCTK_REAL PDstandard4th1gt32 = INITVALUE;
    CCTK_REAL PDstandard4th2gt32 = INITVALUE;
    CCTK_REAL PDstandard4th3gt32 = INITVALUE;
    CCTK_REAL PDstandard4th11gt32 = INITVALUE;
    CCTK_REAL PDstandard4th22gt32 = INITVALUE;
    CCTK_REAL PDstandard4th33gt32 = INITVALUE;
    CCTK_REAL PDstandard4th12gt32 = INITVALUE;
    CCTK_REAL PDstandard4th13gt32 = INITVALUE;
    CCTK_REAL PDstandard4th21gt32 = INITVALUE;
    CCTK_REAL PDstandard4th23gt32 = INITVALUE;
    CCTK_REAL PDstandard4th31gt32 = INITVALUE;
    CCTK_REAL PDstandard4th32gt32 = INITVALUE;
    CCTK_REAL PDstandard4th1gt33 = INITVALUE;
    CCTK_REAL PDstandard4th2gt33 = INITVALUE;
    CCTK_REAL PDstandard4th3gt33 = INITVALUE;
    CCTK_REAL PDstandard4th11gt33 = INITVALUE;
    CCTK_REAL PDstandard4th22gt33 = INITVALUE;
    CCTK_REAL PDstandard4th33gt33 = INITVALUE;
    CCTK_REAL PDstandard4th12gt33 = INITVALUE;
    CCTK_REAL PDstandard4th13gt33 = INITVALUE;
    CCTK_REAL PDstandard4th21gt33 = INITVALUE;
    CCTK_REAL PDstandard4th23gt33 = INITVALUE;
    CCTK_REAL PDstandard4th31gt33 = INITVALUE;
    CCTK_REAL PDstandard4th32gt33 = INITVALUE;
    CCTK_REAL PDstandard4th1phi = INITVALUE;
    CCTK_REAL PDstandard4th2phi = INITVALUE;
    CCTK_REAL PDstandard4th3phi = INITVALUE;
    CCTK_REAL PDstandard4th11phi = INITVALUE;
    CCTK_REAL PDstandard4th22phi = INITVALUE;
    CCTK_REAL PDstandard4th33phi = INITVALUE;
    CCTK_REAL PDstandard4th12phi = INITVALUE;
    CCTK_REAL PDstandard4th13phi = INITVALUE;
    CCTK_REAL PDstandard4th21phi = INITVALUE;
    CCTK_REAL PDstandard4th23phi = INITVALUE;
    CCTK_REAL PDstandard4th31phi = INITVALUE;
    CCTK_REAL PDstandard4th32phi = INITVALUE;
    CCTK_REAL PDstandard4th1trK = INITVALUE;
    CCTK_REAL PDstandard4th2trK = INITVALUE;
    CCTK_REAL PDstandard4th3trK = INITVALUE;
    CCTK_REAL PDstandard4th1Xt1 = INITVALUE;
    CCTK_REAL PDstandard4th2Xt1 = INITVALUE;
    CCTK_REAL PDstandard4th3Xt1 = INITVALUE;
    CCTK_REAL PDstandard4th1Xt2 = INITVALUE;
    CCTK_REAL PDstandard4th2Xt2 = INITVALUE;
    CCTK_REAL PDstandard4th3Xt2 = INITVALUE;
    CCTK_REAL PDstandard4th1Xt3 = INITVALUE;
    CCTK_REAL PDstandard4th2Xt3 = INITVALUE;
    CCTK_REAL PDstandard4th3Xt3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    alphaL = alpha[index];
    At11L = At11[index];
    At21L = At21[index];
    At22L = At22[index];
    At31L = At31[index];
    At32L = At32[index];
    At33L = At33[index];
    beta1L = beta1[index];
    beta2L = beta2[index];
    beta3L = beta3[index];
    dtalphaL = dtalpha[index];
    dtbeta1L = dtbeta1[index];
    dtbeta2L = dtbeta2[index];
    dtbeta3L = dtbeta3[index];
    gt11L = gt11[index];
    gt21L = gt21[index];
    gt22L = gt22[index];
    gt31L = gt31[index];
    gt32L = gt32[index];
    gt33L = gt33[index];
    phiL = phi[index];
    trKL = trK[index];
    trKrhsL = trKrhs[index];
    Xt1L = Xt1[index];
    Xt1rhsL = Xt1rhs[index];
    Xt2L = Xt2[index];
    Xt2rhsL = Xt2rhs[index];
    Xt3L = Xt3[index];
    Xt3rhsL = Xt3rhs[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandard4th1alpha = PDstandard4th1(alpha, i, j, k);
    PDstandard4th2alpha = PDstandard4th2(alpha, i, j, k);
    PDstandard4th3alpha = PDstandard4th3(alpha, i, j, k);
    PDstandard4th11alpha = PDstandard4th11(alpha, i, j, k);
    PDstandard4th22alpha = PDstandard4th22(alpha, i, j, k);
    PDstandard4th33alpha = PDstandard4th33(alpha, i, j, k);
    PDstandard4th12alpha = PDstandard4th12(alpha, i, j, k);
    PDstandard4th13alpha = PDstandard4th13(alpha, i, j, k);
    PDstandard4th23alpha = PDstandard4th23(alpha, i, j, k);
    PDstandard4th1At11 = PDstandard4th1(At11, i, j, k);
    PDstandard4th2At11 = PDstandard4th2(At11, i, j, k);
    PDstandard4th3At11 = PDstandard4th3(At11, i, j, k);
    PDstandard4th1At21 = PDstandard4th1(At21, i, j, k);
    PDstandard4th2At21 = PDstandard4th2(At21, i, j, k);
    PDstandard4th3At21 = PDstandard4th3(At21, i, j, k);
    PDstandard4th1At22 = PDstandard4th1(At22, i, j, k);
    PDstandard4th2At22 = PDstandard4th2(At22, i, j, k);
    PDstandard4th3At22 = PDstandard4th3(At22, i, j, k);
    PDstandard4th1At31 = PDstandard4th1(At31, i, j, k);
    PDstandard4th2At31 = PDstandard4th2(At31, i, j, k);
    PDstandard4th3At31 = PDstandard4th3(At31, i, j, k);
    PDstandard4th1At32 = PDstandard4th1(At32, i, j, k);
    PDstandard4th2At32 = PDstandard4th2(At32, i, j, k);
    PDstandard4th3At32 = PDstandard4th3(At32, i, j, k);
    PDstandard4th1At33 = PDstandard4th1(At33, i, j, k);
    PDstandard4th2At33 = PDstandard4th2(At33, i, j, k);
    PDstandard4th3At33 = PDstandard4th3(At33, i, j, k);
    PDstandard4th1beta1 = PDstandard4th1(beta1, i, j, k);
    PDstandard4th2beta1 = PDstandard4th2(beta1, i, j, k);
    PDstandard4th3beta1 = PDstandard4th3(beta1, i, j, k);
    PDstandard4th11beta1 = PDstandard4th11(beta1, i, j, k);
    PDstandard4th22beta1 = PDstandard4th22(beta1, i, j, k);
    PDstandard4th33beta1 = PDstandard4th33(beta1, i, j, k);
    PDstandard4th12beta1 = PDstandard4th12(beta1, i, j, k);
    PDstandard4th13beta1 = PDstandard4th13(beta1, i, j, k);
    PDstandard4th23beta1 = PDstandard4th23(beta1, i, j, k);
    PDstandard4th1beta2 = PDstandard4th1(beta2, i, j, k);
    PDstandard4th2beta2 = PDstandard4th2(beta2, i, j, k);
    PDstandard4th3beta2 = PDstandard4th3(beta2, i, j, k);
    PDstandard4th11beta2 = PDstandard4th11(beta2, i, j, k);
    PDstandard4th22beta2 = PDstandard4th22(beta2, i, j, k);
    PDstandard4th33beta2 = PDstandard4th33(beta2, i, j, k);
    PDstandard4th12beta2 = PDstandard4th12(beta2, i, j, k);
    PDstandard4th13beta2 = PDstandard4th13(beta2, i, j, k);
    PDstandard4th23beta2 = PDstandard4th23(beta2, i, j, k);
    PDstandard4th1beta3 = PDstandard4th1(beta3, i, j, k);
    PDstandard4th2beta3 = PDstandard4th2(beta3, i, j, k);
    PDstandard4th3beta3 = PDstandard4th3(beta3, i, j, k);
    PDstandard4th11beta3 = PDstandard4th11(beta3, i, j, k);
    PDstandard4th22beta3 = PDstandard4th22(beta3, i, j, k);
    PDstandard4th33beta3 = PDstandard4th33(beta3, i, j, k);
    PDstandard4th12beta3 = PDstandard4th12(beta3, i, j, k);
    PDstandard4th13beta3 = PDstandard4th13(beta3, i, j, k);
    PDstandard4th23beta3 = PDstandard4th23(beta3, i, j, k);
    PDstandard4th1gt11 = PDstandard4th1(gt11, i, j, k);
    PDstandard4th2gt11 = PDstandard4th2(gt11, i, j, k);
    PDstandard4th3gt11 = PDstandard4th3(gt11, i, j, k);
    PDstandard4th11gt11 = PDstandard4th11(gt11, i, j, k);
    PDstandard4th22gt11 = PDstandard4th22(gt11, i, j, k);
    PDstandard4th33gt11 = PDstandard4th33(gt11, i, j, k);
    PDstandard4th12gt11 = PDstandard4th12(gt11, i, j, k);
    PDstandard4th13gt11 = PDstandard4th13(gt11, i, j, k);
    PDstandard4th23gt11 = PDstandard4th23(gt11, i, j, k);
    PDstandard4th1gt21 = PDstandard4th1(gt21, i, j, k);
    PDstandard4th2gt21 = PDstandard4th2(gt21, i, j, k);
    PDstandard4th3gt21 = PDstandard4th3(gt21, i, j, k);
    PDstandard4th11gt21 = PDstandard4th11(gt21, i, j, k);
    PDstandard4th22gt21 = PDstandard4th22(gt21, i, j, k);
    PDstandard4th33gt21 = PDstandard4th33(gt21, i, j, k);
    PDstandard4th12gt21 = PDstandard4th12(gt21, i, j, k);
    PDstandard4th13gt21 = PDstandard4th13(gt21, i, j, k);
    PDstandard4th23gt21 = PDstandard4th23(gt21, i, j, k);
    PDstandard4th1gt22 = PDstandard4th1(gt22, i, j, k);
    PDstandard4th2gt22 = PDstandard4th2(gt22, i, j, k);
    PDstandard4th3gt22 = PDstandard4th3(gt22, i, j, k);
    PDstandard4th11gt22 = PDstandard4th11(gt22, i, j, k);
    PDstandard4th22gt22 = PDstandard4th22(gt22, i, j, k);
    PDstandard4th33gt22 = PDstandard4th33(gt22, i, j, k);
    PDstandard4th12gt22 = PDstandard4th12(gt22, i, j, k);
    PDstandard4th13gt22 = PDstandard4th13(gt22, i, j, k);
    PDstandard4th23gt22 = PDstandard4th23(gt22, i, j, k);
    PDstandard4th1gt31 = PDstandard4th1(gt31, i, j, k);
    PDstandard4th2gt31 = PDstandard4th2(gt31, i, j, k);
    PDstandard4th3gt31 = PDstandard4th3(gt31, i, j, k);
    PDstandard4th11gt31 = PDstandard4th11(gt31, i, j, k);
    PDstandard4th22gt31 = PDstandard4th22(gt31, i, j, k);
    PDstandard4th33gt31 = PDstandard4th33(gt31, i, j, k);
    PDstandard4th12gt31 = PDstandard4th12(gt31, i, j, k);
    PDstandard4th13gt31 = PDstandard4th13(gt31, i, j, k);
    PDstandard4th23gt31 = PDstandard4th23(gt31, i, j, k);
    PDstandard4th1gt32 = PDstandard4th1(gt32, i, j, k);
    PDstandard4th2gt32 = PDstandard4th2(gt32, i, j, k);
    PDstandard4th3gt32 = PDstandard4th3(gt32, i, j, k);
    PDstandard4th11gt32 = PDstandard4th11(gt32, i, j, k);
    PDstandard4th22gt32 = PDstandard4th22(gt32, i, j, k);
    PDstandard4th33gt32 = PDstandard4th33(gt32, i, j, k);
    PDstandard4th12gt32 = PDstandard4th12(gt32, i, j, k);
    PDstandard4th13gt32 = PDstandard4th13(gt32, i, j, k);
    PDstandard4th23gt32 = PDstandard4th23(gt32, i, j, k);
    PDstandard4th1gt33 = PDstandard4th1(gt33, i, j, k);
    PDstandard4th2gt33 = PDstandard4th2(gt33, i, j, k);
    PDstandard4th3gt33 = PDstandard4th3(gt33, i, j, k);
    PDstandard4th11gt33 = PDstandard4th11(gt33, i, j, k);
    PDstandard4th22gt33 = PDstandard4th22(gt33, i, j, k);
    PDstandard4th33gt33 = PDstandard4th33(gt33, i, j, k);
    PDstandard4th12gt33 = PDstandard4th12(gt33, i, j, k);
    PDstandard4th13gt33 = PDstandard4th13(gt33, i, j, k);
    PDstandard4th23gt33 = PDstandard4th23(gt33, i, j, k);
    PDstandard4th1phi = PDstandard4th1(phi, i, j, k);
    PDstandard4th2phi = PDstandard4th2(phi, i, j, k);
    PDstandard4th3phi = PDstandard4th3(phi, i, j, k);
    PDstandard4th11phi = PDstandard4th11(phi, i, j, k);
    PDstandard4th22phi = PDstandard4th22(phi, i, j, k);
    PDstandard4th33phi = PDstandard4th33(phi, i, j, k);
    PDstandard4th12phi = PDstandard4th12(phi, i, j, k);
    PDstandard4th13phi = PDstandard4th13(phi, i, j, k);
    PDstandard4th23phi = PDstandard4th23(phi, i, j, k);
    PDstandard4th1trK = PDstandard4th1(trK, i, j, k);
    PDstandard4th2trK = PDstandard4th2(trK, i, j, k);
    PDstandard4th3trK = PDstandard4th3(trK, i, j, k);
    PDstandard4th1Xt1 = PDstandard4th1(Xt1, i, j, k);
    PDstandard4th2Xt1 = PDstandard4th2(Xt1, i, j, k);
    PDstandard4th3Xt1 = PDstandard4th3(Xt1, i, j, k);
    PDstandard4th1Xt2 = PDstandard4th1(Xt2, i, j, k);
    PDstandard4th2Xt2 = PDstandard4th2(Xt2, i, j, k);
    PDstandard4th3Xt2 = PDstandard4th3(Xt2, i, j, k);
    PDstandard4th1Xt3 = PDstandard4th1(Xt3, i, j, k);
    PDstandard4th2Xt3 = PDstandard4th2(Xt3, i, j, k);
    PDstandard4th3Xt3 = PDstandard4th3(Xt3, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    detgt  =  1;
    
    ddetgt1  =  0;
    
    ddetgt2  =  0;
    
    ddetgt3  =  0;
    
    gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt32L));
    
    gtu12  =  (gt31L*gt32L - gt21L*gt33L)*INV(detgt);
    
    gtu13  =  (-(gt22L*gt31L) + gt21L*gt32L)*INV(detgt);
    
    gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt31L));
    
    gtu23  =  (gt21L*gt31L - gt11L*gt32L)*INV(detgt);
    
    gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt21L));
    
    dgtu111  =  -2*(gtu11*gtu12*PDstandard4th1gt21 + gtu11*gtu13*PDstandard4th1gt31 + gtu12*gtu13*PDstandard4th1gt32) - 
        PDstandard4th1gt11*SQR(gtu11) - PDstandard4th1gt22*SQR(gtu12) - PDstandard4th1gt33*SQR(gtu13);
    
    dgtu121  =  -(gtu11*gtu12*PDstandard4th1gt11) - gtu11*gtu22*PDstandard4th1gt21 - gtu12*gtu22*PDstandard4th1gt22 - 
        gtu12*gtu13*PDstandard4th1gt31 - gtu11*gtu23*PDstandard4th1gt31 - gtu13*gtu22*PDstandard4th1gt32 - 
        gtu12*gtu23*PDstandard4th1gt32 - gtu13*gtu23*PDstandard4th1gt33 - PDstandard4th1gt21*SQR(gtu12);
    
    dgtu131  =  -(gtu11*gtu13*PDstandard4th1gt11) - gtu12*gtu13*PDstandard4th1gt21 - gtu11*gtu23*PDstandard4th1gt21 - 
        gtu12*gtu23*PDstandard4th1gt22 - gtu11*gtu33*PDstandard4th1gt31 - gtu13*gtu23*PDstandard4th1gt32 - 
        gtu12*gtu33*PDstandard4th1gt32 - gtu13*gtu33*PDstandard4th1gt33 - PDstandard4th1gt31*SQR(gtu13);
    
    dgtu221  =  -2*(gtu12*gtu22*PDstandard4th1gt21 + gtu12*gtu23*PDstandard4th1gt31 + gtu22*gtu23*PDstandard4th1gt32) - 
        PDstandard4th1gt11*SQR(gtu12) - PDstandard4th1gt22*SQR(gtu22) - PDstandard4th1gt33*SQR(gtu23);
    
    dgtu231  =  -(gtu12*gtu13*PDstandard4th1gt11) - gtu13*gtu22*PDstandard4th1gt21 - gtu12*gtu23*PDstandard4th1gt21 - 
        gtu22*gtu23*PDstandard4th1gt22 - gtu13*gtu23*PDstandard4th1gt31 - gtu12*gtu33*PDstandard4th1gt31 - 
        gtu22*gtu33*PDstandard4th1gt32 - gtu23*gtu33*PDstandard4th1gt33 - PDstandard4th1gt32*SQR(gtu23);
    
    dgtu331  =  -2*(gtu13*gtu23*PDstandard4th1gt21 + gtu13*gtu33*PDstandard4th1gt31 + gtu23*gtu33*PDstandard4th1gt32) - 
        PDstandard4th1gt11*SQR(gtu13) - PDstandard4th1gt22*SQR(gtu23) - PDstandard4th1gt33*SQR(gtu33);
    
    dgtu112  =  -2*(gtu11*gtu12*PDstandard4th2gt21 + gtu11*gtu13*PDstandard4th2gt31 + gtu12*gtu13*PDstandard4th2gt32) - 
        PDstandard4th2gt11*SQR(gtu11) - PDstandard4th2gt22*SQR(gtu12) - PDstandard4th2gt33*SQR(gtu13);
    
    dgtu122  =  -(gtu11*gtu12*PDstandard4th2gt11) - gtu11*gtu22*PDstandard4th2gt21 - gtu12*gtu22*PDstandard4th2gt22 - 
        gtu12*gtu13*PDstandard4th2gt31 - gtu11*gtu23*PDstandard4th2gt31 - gtu13*gtu22*PDstandard4th2gt32 - 
        gtu12*gtu23*PDstandard4th2gt32 - gtu13*gtu23*PDstandard4th2gt33 - PDstandard4th2gt21*SQR(gtu12);
    
    dgtu132  =  -(gtu11*gtu13*PDstandard4th2gt11) - gtu12*gtu13*PDstandard4th2gt21 - gtu11*gtu23*PDstandard4th2gt21 - 
        gtu12*gtu23*PDstandard4th2gt22 - gtu11*gtu33*PDstandard4th2gt31 - gtu13*gtu23*PDstandard4th2gt32 - 
        gtu12*gtu33*PDstandard4th2gt32 - gtu13*gtu33*PDstandard4th2gt33 - PDstandard4th2gt31*SQR(gtu13);
    
    dgtu222  =  -2*(gtu12*gtu22*PDstandard4th2gt21 + gtu12*gtu23*PDstandard4th2gt31 + gtu22*gtu23*PDstandard4th2gt32) - 
        PDstandard4th2gt11*SQR(gtu12) - PDstandard4th2gt22*SQR(gtu22) - PDstandard4th2gt33*SQR(gtu23);
    
    dgtu232  =  -(gtu12*gtu13*PDstandard4th2gt11) - gtu13*gtu22*PDstandard4th2gt21 - gtu12*gtu23*PDstandard4th2gt21 - 
        gtu22*gtu23*PDstandard4th2gt22 - gtu13*gtu23*PDstandard4th2gt31 - gtu12*gtu33*PDstandard4th2gt31 - 
        gtu22*gtu33*PDstandard4th2gt32 - gtu23*gtu33*PDstandard4th2gt33 - PDstandard4th2gt32*SQR(gtu23);
    
    dgtu332  =  -2*(gtu13*gtu23*PDstandard4th2gt21 + gtu13*gtu33*PDstandard4th2gt31 + gtu23*gtu33*PDstandard4th2gt32) - 
        PDstandard4th2gt11*SQR(gtu13) - PDstandard4th2gt22*SQR(gtu23) - PDstandard4th2gt33*SQR(gtu33);
    
    dgtu113  =  -2*(gtu11*gtu12*PDstandard4th3gt21 + gtu11*gtu13*PDstandard4th3gt31 + gtu12*gtu13*PDstandard4th3gt32) - 
        PDstandard4th3gt11*SQR(gtu11) - PDstandard4th3gt22*SQR(gtu12) - PDstandard4th3gt33*SQR(gtu13);
    
    dgtu123  =  -(gtu11*gtu12*PDstandard4th3gt11) - gtu11*gtu22*PDstandard4th3gt21 - gtu12*gtu22*PDstandard4th3gt22 - 
        gtu12*gtu13*PDstandard4th3gt31 - gtu11*gtu23*PDstandard4th3gt31 - gtu13*gtu22*PDstandard4th3gt32 - 
        gtu12*gtu23*PDstandard4th3gt32 - gtu13*gtu23*PDstandard4th3gt33 - PDstandard4th3gt21*SQR(gtu12);
    
    dgtu133  =  -(gtu11*gtu13*PDstandard4th3gt11) - gtu12*gtu13*PDstandard4th3gt21 - gtu11*gtu23*PDstandard4th3gt21 - 
        gtu12*gtu23*PDstandard4th3gt22 - gtu11*gtu33*PDstandard4th3gt31 - gtu13*gtu23*PDstandard4th3gt32 - 
        gtu12*gtu33*PDstandard4th3gt32 - gtu13*gtu33*PDstandard4th3gt33 - PDstandard4th3gt31*SQR(gtu13);
    
    dgtu223  =  -2*(gtu12*gtu22*PDstandard4th3gt21 + gtu12*gtu23*PDstandard4th3gt31 + gtu22*gtu23*PDstandard4th3gt32) - 
        PDstandard4th3gt11*SQR(gtu12) - PDstandard4th3gt22*SQR(gtu22) - PDstandard4th3gt33*SQR(gtu23);
    
    dgtu233  =  -(gtu12*gtu13*PDstandard4th3gt11) - gtu13*gtu22*PDstandard4th3gt21 - gtu12*gtu23*PDstandard4th3gt21 - 
        gtu22*gtu23*PDstandard4th3gt22 - gtu13*gtu23*PDstandard4th3gt31 - gtu12*gtu33*PDstandard4th3gt31 - 
        gtu22*gtu33*PDstandard4th3gt32 - gtu23*gtu33*PDstandard4th3gt33 - PDstandard4th3gt32*SQR(gtu23);
    
    dgtu333  =  -2*(gtu13*gtu23*PDstandard4th3gt21 + gtu13*gtu33*PDstandard4th3gt31 + gtu23*gtu33*PDstandard4th3gt32) - 
        PDstandard4th3gt11*SQR(gtu13) - PDstandard4th3gt22*SQR(gtu23) - PDstandard4th3gt33*SQR(gtu33);
    
    ddgtu1111  =  -2*(gtu11*gtu12*PDstandard4th11gt21 + gtu11*gtu13*PDstandard4th11gt31 + gtu12*gtu13*PDstandard4th11gt32 + 
           dgtu111*gtu11*PDstandard4th1gt11 + dgtu121*gtu11*PDstandard4th1gt21 + dgtu111*gtu12*PDstandard4th1gt21 + 
           dgtu121*gtu12*PDstandard4th1gt22 + dgtu131*gtu11*PDstandard4th1gt31 + dgtu111*gtu13*PDstandard4th1gt31 + 
           dgtu131*gtu12*PDstandard4th1gt32 + dgtu121*gtu13*PDstandard4th1gt32 + dgtu131*gtu13*PDstandard4th1gt33) - 
        PDstandard4th11gt11*SQR(gtu11) - PDstandard4th11gt22*SQR(gtu12) - PDstandard4th11gt33*SQR(gtu13);
    
    ddgtu1211  =  -(gtu11*gtu12*PDstandard4th11gt11) - gtu11*gtu22*PDstandard4th11gt21 - gtu12*gtu22*PDstandard4th11gt22 - 
        gtu12*gtu13*PDstandard4th11gt31 - gtu11*gtu23*PDstandard4th11gt31 - gtu13*gtu22*PDstandard4th11gt32 - 
        gtu12*gtu23*PDstandard4th11gt32 - gtu13*gtu23*PDstandard4th11gt33 - dgtu121*gtu11*PDstandard4th1gt11 - 
        dgtu111*gtu12*PDstandard4th1gt11 - dgtu221*gtu11*PDstandard4th1gt21 - 2*dgtu121*gtu12*PDstandard4th1gt21 - 
        dgtu111*gtu22*PDstandard4th1gt21 - dgtu221*gtu12*PDstandard4th1gt22 - dgtu121*gtu22*PDstandard4th1gt22 - 
        dgtu231*gtu11*PDstandard4th1gt31 - dgtu131*gtu12*PDstandard4th1gt31 - dgtu121*gtu13*PDstandard4th1gt31 - 
        dgtu111*gtu23*PDstandard4th1gt31 - dgtu231*gtu12*PDstandard4th1gt32 - dgtu221*gtu13*PDstandard4th1gt32 - 
        dgtu131*gtu22*PDstandard4th1gt32 - dgtu121*gtu23*PDstandard4th1gt32 - dgtu231*gtu13*PDstandard4th1gt33 - 
        dgtu131*gtu23*PDstandard4th1gt33 - PDstandard4th11gt21*SQR(gtu12);
    
    ddgtu1311  =  -(gtu11*gtu13*PDstandard4th11gt11) - gtu12*gtu13*PDstandard4th11gt21 - gtu11*gtu23*PDstandard4th11gt21 - 
        gtu12*gtu23*PDstandard4th11gt22 - gtu11*gtu33*PDstandard4th11gt31 - gtu13*gtu23*PDstandard4th11gt32 - 
        gtu12*gtu33*PDstandard4th11gt32 - gtu13*gtu33*PDstandard4th11gt33 - dgtu131*gtu11*PDstandard4th1gt11 - 
        dgtu111*gtu13*PDstandard4th1gt11 - dgtu231*gtu11*PDstandard4th1gt21 - dgtu131*gtu12*PDstandard4th1gt21 - 
        dgtu121*gtu13*PDstandard4th1gt21 - dgtu111*gtu23*PDstandard4th1gt21 - dgtu231*gtu12*PDstandard4th1gt22 - 
        dgtu121*gtu23*PDstandard4th1gt22 - dgtu331*gtu11*PDstandard4th1gt31 - 2*dgtu131*gtu13*PDstandard4th1gt31 - 
        dgtu111*gtu33*PDstandard4th1gt31 - dgtu331*gtu12*PDstandard4th1gt32 - dgtu231*gtu13*PDstandard4th1gt32 - 
        dgtu131*gtu23*PDstandard4th1gt32 - dgtu121*gtu33*PDstandard4th1gt32 - dgtu331*gtu13*PDstandard4th1gt33 - 
        dgtu131*gtu33*PDstandard4th1gt33 - PDstandard4th11gt31*SQR(gtu13);
    
    ddgtu1121  =  -2*(gtu11*gtu12*PDstandard4th12gt21 + gtu11*gtu13*PDstandard4th12gt31 + gtu12*gtu13*PDstandard4th12gt32 + 
           dgtu112*gtu11*PDstandard4th1gt11 + dgtu122*gtu11*PDstandard4th1gt21 + dgtu112*gtu12*PDstandard4th1gt21 + 
           dgtu122*gtu12*PDstandard4th1gt22 + dgtu132*gtu11*PDstandard4th1gt31 + dgtu112*gtu13*PDstandard4th1gt31 + 
           dgtu132*gtu12*PDstandard4th1gt32 + dgtu122*gtu13*PDstandard4th1gt32 + dgtu132*gtu13*PDstandard4th1gt33) - 
        PDstandard4th12gt11*SQR(gtu11) - PDstandard4th12gt22*SQR(gtu12) - PDstandard4th12gt33*SQR(gtu13);
    
    ddgtu1221  =  -(gtu11*gtu12*PDstandard4th12gt11) - gtu11*gtu22*PDstandard4th12gt21 - gtu12*gtu22*PDstandard4th12gt22 - 
        gtu12*gtu13*PDstandard4th12gt31 - gtu11*gtu23*PDstandard4th12gt31 - gtu13*gtu22*PDstandard4th12gt32 - 
        gtu12*gtu23*PDstandard4th12gt32 - gtu13*gtu23*PDstandard4th12gt33 - dgtu122*gtu11*PDstandard4th1gt11 - 
        dgtu112*gtu12*PDstandard4th1gt11 - dgtu222*gtu11*PDstandard4th1gt21 - 2*dgtu122*gtu12*PDstandard4th1gt21 - 
        dgtu112*gtu22*PDstandard4th1gt21 - dgtu222*gtu12*PDstandard4th1gt22 - dgtu122*gtu22*PDstandard4th1gt22 - 
        dgtu232*gtu11*PDstandard4th1gt31 - dgtu132*gtu12*PDstandard4th1gt31 - dgtu122*gtu13*PDstandard4th1gt31 - 
        dgtu112*gtu23*PDstandard4th1gt31 - dgtu232*gtu12*PDstandard4th1gt32 - dgtu222*gtu13*PDstandard4th1gt32 - 
        dgtu132*gtu22*PDstandard4th1gt32 - dgtu122*gtu23*PDstandard4th1gt32 - dgtu232*gtu13*PDstandard4th1gt33 - 
        dgtu132*gtu23*PDstandard4th1gt33 - PDstandard4th12gt21*SQR(gtu12);
    
    ddgtu1321  =  -(gtu11*gtu13*PDstandard4th12gt11) - gtu12*gtu13*PDstandard4th12gt21 - gtu11*gtu23*PDstandard4th12gt21 - 
        gtu12*gtu23*PDstandard4th12gt22 - gtu11*gtu33*PDstandard4th12gt31 - gtu13*gtu23*PDstandard4th12gt32 - 
        gtu12*gtu33*PDstandard4th12gt32 - gtu13*gtu33*PDstandard4th12gt33 - dgtu132*gtu11*PDstandard4th1gt11 - 
        dgtu112*gtu13*PDstandard4th1gt11 - dgtu232*gtu11*PDstandard4th1gt21 - dgtu132*gtu12*PDstandard4th1gt21 - 
        dgtu122*gtu13*PDstandard4th1gt21 - dgtu112*gtu23*PDstandard4th1gt21 - dgtu232*gtu12*PDstandard4th1gt22 - 
        dgtu122*gtu23*PDstandard4th1gt22 - dgtu332*gtu11*PDstandard4th1gt31 - 2*dgtu132*gtu13*PDstandard4th1gt31 - 
        dgtu112*gtu33*PDstandard4th1gt31 - dgtu332*gtu12*PDstandard4th1gt32 - dgtu232*gtu13*PDstandard4th1gt32 - 
        dgtu132*gtu23*PDstandard4th1gt32 - dgtu122*gtu33*PDstandard4th1gt32 - dgtu332*gtu13*PDstandard4th1gt33 - 
        dgtu132*gtu33*PDstandard4th1gt33 - PDstandard4th12gt31*SQR(gtu13);
    
    ddgtu2221  =  -2*(gtu12*gtu22*PDstandard4th12gt21 + gtu12*gtu23*PDstandard4th12gt31 + gtu22*gtu23*PDstandard4th12gt32 + 
           dgtu122*gtu12*PDstandard4th1gt11 + dgtu222*gtu12*PDstandard4th1gt21 + dgtu122*gtu22*PDstandard4th1gt21 + 
           dgtu222*gtu22*PDstandard4th1gt22 + dgtu232*gtu12*PDstandard4th1gt31 + dgtu122*gtu23*PDstandard4th1gt31 + 
           dgtu232*gtu22*PDstandard4th1gt32 + dgtu222*gtu23*PDstandard4th1gt32 + dgtu232*gtu23*PDstandard4th1gt33) - 
        PDstandard4th12gt11*SQR(gtu12) - PDstandard4th12gt22*SQR(gtu22) - PDstandard4th12gt33*SQR(gtu23);
    
    ddgtu2321  =  -(gtu12*gtu13*PDstandard4th12gt11) - gtu13*gtu22*PDstandard4th12gt21 - gtu12*gtu23*PDstandard4th12gt21 - 
        gtu22*gtu23*PDstandard4th12gt22 - gtu13*gtu23*PDstandard4th12gt31 - gtu12*gtu33*PDstandard4th12gt31 - 
        gtu22*gtu33*PDstandard4th12gt32 - gtu23*gtu33*PDstandard4th12gt33 - dgtu132*gtu12*PDstandard4th1gt11 - 
        dgtu122*gtu13*PDstandard4th1gt11 - dgtu232*gtu12*PDstandard4th1gt21 - dgtu222*gtu13*PDstandard4th1gt21 - 
        dgtu132*gtu22*PDstandard4th1gt21 - dgtu122*gtu23*PDstandard4th1gt21 - dgtu232*gtu22*PDstandard4th1gt22 - 
        dgtu222*gtu23*PDstandard4th1gt22 - dgtu332*gtu12*PDstandard4th1gt31 - dgtu232*gtu13*PDstandard4th1gt31 - 
        dgtu132*gtu23*PDstandard4th1gt31 - dgtu122*gtu33*PDstandard4th1gt31 - dgtu332*gtu22*PDstandard4th1gt32 - 
        2*dgtu232*gtu23*PDstandard4th1gt32 - dgtu222*gtu33*PDstandard4th1gt32 - dgtu332*gtu23*PDstandard4th1gt33 - 
        dgtu232*gtu33*PDstandard4th1gt33 - PDstandard4th12gt32*SQR(gtu23);
    
    ddgtu1131  =  -2*(gtu11*gtu12*PDstandard4th13gt21 + gtu11*gtu13*PDstandard4th13gt31 + gtu12*gtu13*PDstandard4th13gt32 + 
           dgtu113*gtu11*PDstandard4th1gt11 + dgtu123*gtu11*PDstandard4th1gt21 + dgtu113*gtu12*PDstandard4th1gt21 + 
           dgtu123*gtu12*PDstandard4th1gt22 + dgtu133*gtu11*PDstandard4th1gt31 + dgtu113*gtu13*PDstandard4th1gt31 + 
           dgtu133*gtu12*PDstandard4th1gt32 + dgtu123*gtu13*PDstandard4th1gt32 + dgtu133*gtu13*PDstandard4th1gt33) - 
        PDstandard4th13gt11*SQR(gtu11) - PDstandard4th13gt22*SQR(gtu12) - PDstandard4th13gt33*SQR(gtu13);
    
    ddgtu1231  =  -(gtu11*gtu12*PDstandard4th13gt11) - gtu11*gtu22*PDstandard4th13gt21 - gtu12*gtu22*PDstandard4th13gt22 - 
        gtu12*gtu13*PDstandard4th13gt31 - gtu11*gtu23*PDstandard4th13gt31 - gtu13*gtu22*PDstandard4th13gt32 - 
        gtu12*gtu23*PDstandard4th13gt32 - gtu13*gtu23*PDstandard4th13gt33 - dgtu123*gtu11*PDstandard4th1gt11 - 
        dgtu113*gtu12*PDstandard4th1gt11 - dgtu223*gtu11*PDstandard4th1gt21 - 2*dgtu123*gtu12*PDstandard4th1gt21 - 
        dgtu113*gtu22*PDstandard4th1gt21 - dgtu223*gtu12*PDstandard4th1gt22 - dgtu123*gtu22*PDstandard4th1gt22 - 
        dgtu233*gtu11*PDstandard4th1gt31 - dgtu133*gtu12*PDstandard4th1gt31 - dgtu123*gtu13*PDstandard4th1gt31 - 
        dgtu113*gtu23*PDstandard4th1gt31 - dgtu233*gtu12*PDstandard4th1gt32 - dgtu223*gtu13*PDstandard4th1gt32 - 
        dgtu133*gtu22*PDstandard4th1gt32 - dgtu123*gtu23*PDstandard4th1gt32 - dgtu233*gtu13*PDstandard4th1gt33 - 
        dgtu133*gtu23*PDstandard4th1gt33 - PDstandard4th13gt21*SQR(gtu12);
    
    ddgtu1331  =  -(gtu11*gtu13*PDstandard4th13gt11) - gtu12*gtu13*PDstandard4th13gt21 - gtu11*gtu23*PDstandard4th13gt21 - 
        gtu12*gtu23*PDstandard4th13gt22 - gtu11*gtu33*PDstandard4th13gt31 - gtu13*gtu23*PDstandard4th13gt32 - 
        gtu12*gtu33*PDstandard4th13gt32 - gtu13*gtu33*PDstandard4th13gt33 - dgtu133*gtu11*PDstandard4th1gt11 - 
        dgtu113*gtu13*PDstandard4th1gt11 - dgtu233*gtu11*PDstandard4th1gt21 - dgtu133*gtu12*PDstandard4th1gt21 - 
        dgtu123*gtu13*PDstandard4th1gt21 - dgtu113*gtu23*PDstandard4th1gt21 - dgtu233*gtu12*PDstandard4th1gt22 - 
        dgtu123*gtu23*PDstandard4th1gt22 - dgtu333*gtu11*PDstandard4th1gt31 - 2*dgtu133*gtu13*PDstandard4th1gt31 - 
        dgtu113*gtu33*PDstandard4th1gt31 - dgtu333*gtu12*PDstandard4th1gt32 - dgtu233*gtu13*PDstandard4th1gt32 - 
        dgtu133*gtu23*PDstandard4th1gt32 - dgtu123*gtu33*PDstandard4th1gt32 - dgtu333*gtu13*PDstandard4th1gt33 - 
        dgtu133*gtu33*PDstandard4th1gt33 - PDstandard4th13gt31*SQR(gtu13);
    
    ddgtu2331  =  -(gtu12*gtu13*PDstandard4th13gt11) - gtu13*gtu22*PDstandard4th13gt21 - gtu12*gtu23*PDstandard4th13gt21 - 
        gtu22*gtu23*PDstandard4th13gt22 - gtu13*gtu23*PDstandard4th13gt31 - gtu12*gtu33*PDstandard4th13gt31 - 
        gtu22*gtu33*PDstandard4th13gt32 - gtu23*gtu33*PDstandard4th13gt33 - dgtu133*gtu12*PDstandard4th1gt11 - 
        dgtu123*gtu13*PDstandard4th1gt11 - dgtu233*gtu12*PDstandard4th1gt21 - dgtu223*gtu13*PDstandard4th1gt21 - 
        dgtu133*gtu22*PDstandard4th1gt21 - dgtu123*gtu23*PDstandard4th1gt21 - dgtu233*gtu22*PDstandard4th1gt22 - 
        dgtu223*gtu23*PDstandard4th1gt22 - dgtu333*gtu12*PDstandard4th1gt31 - dgtu233*gtu13*PDstandard4th1gt31 - 
        dgtu133*gtu23*PDstandard4th1gt31 - dgtu123*gtu33*PDstandard4th1gt31 - dgtu333*gtu22*PDstandard4th1gt32 - 
        2*dgtu233*gtu23*PDstandard4th1gt32 - dgtu223*gtu33*PDstandard4th1gt32 - dgtu333*gtu23*PDstandard4th1gt33 - 
        dgtu233*gtu33*PDstandard4th1gt33 - PDstandard4th13gt32*SQR(gtu23);
    
    ddgtu3331  =  -2*(gtu13*gtu23*PDstandard4th13gt21 + gtu13*gtu33*PDstandard4th13gt31 + gtu23*gtu33*PDstandard4th13gt32 + 
           dgtu133*gtu13*PDstandard4th1gt11 + dgtu233*gtu13*PDstandard4th1gt21 + dgtu133*gtu23*PDstandard4th1gt21 + 
           dgtu233*gtu23*PDstandard4th1gt22 + dgtu333*gtu13*PDstandard4th1gt31 + dgtu133*gtu33*PDstandard4th1gt31 + 
           dgtu333*gtu23*PDstandard4th1gt32 + dgtu233*gtu33*PDstandard4th1gt32 + dgtu333*gtu33*PDstandard4th1gt33) - 
        PDstandard4th13gt11*SQR(gtu13) - PDstandard4th13gt22*SQR(gtu23) - PDstandard4th13gt33*SQR(gtu33);
    
    ddgtu1222  =  -(gtu11*gtu12*PDstandard4th22gt11) - gtu11*gtu22*PDstandard4th22gt21 - gtu12*gtu22*PDstandard4th22gt22 - 
        gtu12*gtu13*PDstandard4th22gt31 - gtu11*gtu23*PDstandard4th22gt31 - gtu13*gtu22*PDstandard4th22gt32 - 
        gtu12*gtu23*PDstandard4th22gt32 - gtu13*gtu23*PDstandard4th22gt33 - dgtu122*gtu11*PDstandard4th2gt11 - 
        dgtu112*gtu12*PDstandard4th2gt11 - dgtu222*gtu11*PDstandard4th2gt21 - 2*dgtu122*gtu12*PDstandard4th2gt21 - 
        dgtu112*gtu22*PDstandard4th2gt21 - dgtu222*gtu12*PDstandard4th2gt22 - dgtu122*gtu22*PDstandard4th2gt22 - 
        dgtu232*gtu11*PDstandard4th2gt31 - dgtu132*gtu12*PDstandard4th2gt31 - dgtu122*gtu13*PDstandard4th2gt31 - 
        dgtu112*gtu23*PDstandard4th2gt31 - dgtu232*gtu12*PDstandard4th2gt32 - dgtu222*gtu13*PDstandard4th2gt32 - 
        dgtu132*gtu22*PDstandard4th2gt32 - dgtu122*gtu23*PDstandard4th2gt32 - dgtu232*gtu13*PDstandard4th2gt33 - 
        dgtu132*gtu23*PDstandard4th2gt33 - PDstandard4th22gt21*SQR(gtu12);
    
    ddgtu2222  =  -2*(gtu12*gtu22*PDstandard4th22gt21 + gtu12*gtu23*PDstandard4th22gt31 + gtu22*gtu23*PDstandard4th22gt32 + 
           dgtu122*gtu12*PDstandard4th2gt11 + dgtu222*gtu12*PDstandard4th2gt21 + dgtu122*gtu22*PDstandard4th2gt21 + 
           dgtu222*gtu22*PDstandard4th2gt22 + dgtu232*gtu12*PDstandard4th2gt31 + dgtu122*gtu23*PDstandard4th2gt31 + 
           dgtu232*gtu22*PDstandard4th2gt32 + dgtu222*gtu23*PDstandard4th2gt32 + dgtu232*gtu23*PDstandard4th2gt33) - 
        PDstandard4th22gt11*SQR(gtu12) - PDstandard4th22gt22*SQR(gtu22) - PDstandard4th22gt33*SQR(gtu23);
    
    ddgtu2322  =  -(gtu12*gtu13*PDstandard4th22gt11) - gtu13*gtu22*PDstandard4th22gt21 - gtu12*gtu23*PDstandard4th22gt21 - 
        gtu22*gtu23*PDstandard4th22gt22 - gtu13*gtu23*PDstandard4th22gt31 - gtu12*gtu33*PDstandard4th22gt31 - 
        gtu22*gtu33*PDstandard4th22gt32 - gtu23*gtu33*PDstandard4th22gt33 - dgtu132*gtu12*PDstandard4th2gt11 - 
        dgtu122*gtu13*PDstandard4th2gt11 - dgtu232*gtu12*PDstandard4th2gt21 - dgtu222*gtu13*PDstandard4th2gt21 - 
        dgtu132*gtu22*PDstandard4th2gt21 - dgtu122*gtu23*PDstandard4th2gt21 - dgtu232*gtu22*PDstandard4th2gt22 - 
        dgtu222*gtu23*PDstandard4th2gt22 - dgtu332*gtu12*PDstandard4th2gt31 - dgtu232*gtu13*PDstandard4th2gt31 - 
        dgtu132*gtu23*PDstandard4th2gt31 - dgtu122*gtu33*PDstandard4th2gt31 - dgtu332*gtu22*PDstandard4th2gt32 - 
        2*dgtu232*gtu23*PDstandard4th2gt32 - dgtu222*gtu33*PDstandard4th2gt32 - dgtu332*gtu23*PDstandard4th2gt33 - 
        dgtu232*gtu33*PDstandard4th2gt33 - PDstandard4th22gt32*SQR(gtu23);
    
    ddgtu1232  =  -(gtu11*gtu12*PDstandard4th23gt11) - gtu11*gtu22*PDstandard4th23gt21 - gtu12*gtu22*PDstandard4th23gt22 - 
        gtu12*gtu13*PDstandard4th23gt31 - gtu11*gtu23*PDstandard4th23gt31 - gtu13*gtu22*PDstandard4th23gt32 - 
        gtu12*gtu23*PDstandard4th23gt32 - gtu13*gtu23*PDstandard4th23gt33 - dgtu123*gtu11*PDstandard4th2gt11 - 
        dgtu113*gtu12*PDstandard4th2gt11 - dgtu223*gtu11*PDstandard4th2gt21 - 2*dgtu123*gtu12*PDstandard4th2gt21 - 
        dgtu113*gtu22*PDstandard4th2gt21 - dgtu223*gtu12*PDstandard4th2gt22 - dgtu123*gtu22*PDstandard4th2gt22 - 
        dgtu233*gtu11*PDstandard4th2gt31 - dgtu133*gtu12*PDstandard4th2gt31 - dgtu123*gtu13*PDstandard4th2gt31 - 
        dgtu113*gtu23*PDstandard4th2gt31 - dgtu233*gtu12*PDstandard4th2gt32 - dgtu223*gtu13*PDstandard4th2gt32 - 
        dgtu133*gtu22*PDstandard4th2gt32 - dgtu123*gtu23*PDstandard4th2gt32 - dgtu233*gtu13*PDstandard4th2gt33 - 
        dgtu133*gtu23*PDstandard4th2gt33 - PDstandard4th23gt21*SQR(gtu12);
    
    ddgtu1332  =  -(gtu11*gtu13*PDstandard4th23gt11) - gtu12*gtu13*PDstandard4th23gt21 - gtu11*gtu23*PDstandard4th23gt21 - 
        gtu12*gtu23*PDstandard4th23gt22 - gtu11*gtu33*PDstandard4th23gt31 - gtu13*gtu23*PDstandard4th23gt32 - 
        gtu12*gtu33*PDstandard4th23gt32 - gtu13*gtu33*PDstandard4th23gt33 - dgtu133*gtu11*PDstandard4th2gt11 - 
        dgtu113*gtu13*PDstandard4th2gt11 - dgtu233*gtu11*PDstandard4th2gt21 - dgtu133*gtu12*PDstandard4th2gt21 - 
        dgtu123*gtu13*PDstandard4th2gt21 - dgtu113*gtu23*PDstandard4th2gt21 - dgtu233*gtu12*PDstandard4th2gt22 - 
        dgtu123*gtu23*PDstandard4th2gt22 - dgtu333*gtu11*PDstandard4th2gt31 - 2*dgtu133*gtu13*PDstandard4th2gt31 - 
        dgtu113*gtu33*PDstandard4th2gt31 - dgtu333*gtu12*PDstandard4th2gt32 - dgtu233*gtu13*PDstandard4th2gt32 - 
        dgtu133*gtu23*PDstandard4th2gt32 - dgtu123*gtu33*PDstandard4th2gt32 - dgtu333*gtu13*PDstandard4th2gt33 - 
        dgtu133*gtu33*PDstandard4th2gt33 - PDstandard4th23gt31*SQR(gtu13);
    
    ddgtu2232  =  -2*(gtu12*gtu22*PDstandard4th23gt21 + gtu12*gtu23*PDstandard4th23gt31 + gtu22*gtu23*PDstandard4th23gt32 + 
           dgtu123*gtu12*PDstandard4th2gt11 + dgtu223*gtu12*PDstandard4th2gt21 + dgtu123*gtu22*PDstandard4th2gt21 + 
           dgtu223*gtu22*PDstandard4th2gt22 + dgtu233*gtu12*PDstandard4th2gt31 + dgtu123*gtu23*PDstandard4th2gt31 + 
           dgtu233*gtu22*PDstandard4th2gt32 + dgtu223*gtu23*PDstandard4th2gt32 + dgtu233*gtu23*PDstandard4th2gt33) - 
        PDstandard4th23gt11*SQR(gtu12) - PDstandard4th23gt22*SQR(gtu22) - PDstandard4th23gt33*SQR(gtu23);
    
    ddgtu2332  =  -(gtu12*gtu13*PDstandard4th23gt11) - gtu13*gtu22*PDstandard4th23gt21 - gtu12*gtu23*PDstandard4th23gt21 - 
        gtu22*gtu23*PDstandard4th23gt22 - gtu13*gtu23*PDstandard4th23gt31 - gtu12*gtu33*PDstandard4th23gt31 - 
        gtu22*gtu33*PDstandard4th23gt32 - gtu23*gtu33*PDstandard4th23gt33 - dgtu133*gtu12*PDstandard4th2gt11 - 
        dgtu123*gtu13*PDstandard4th2gt11 - dgtu233*gtu12*PDstandard4th2gt21 - dgtu223*gtu13*PDstandard4th2gt21 - 
        dgtu133*gtu22*PDstandard4th2gt21 - dgtu123*gtu23*PDstandard4th2gt21 - dgtu233*gtu22*PDstandard4th2gt22 - 
        dgtu223*gtu23*PDstandard4th2gt22 - dgtu333*gtu12*PDstandard4th2gt31 - dgtu233*gtu13*PDstandard4th2gt31 - 
        dgtu133*gtu23*PDstandard4th2gt31 - dgtu123*gtu33*PDstandard4th2gt31 - dgtu333*gtu22*PDstandard4th2gt32 - 
        2*dgtu233*gtu23*PDstandard4th2gt32 - dgtu223*gtu33*PDstandard4th2gt32 - dgtu333*gtu23*PDstandard4th2gt33 - 
        dgtu233*gtu33*PDstandard4th2gt33 - PDstandard4th23gt32*SQR(gtu23);
    
    ddgtu3332  =  -2*(gtu13*gtu23*PDstandard4th23gt21 + gtu13*gtu33*PDstandard4th23gt31 + gtu23*gtu33*PDstandard4th23gt32 + 
           dgtu133*gtu13*PDstandard4th2gt11 + dgtu233*gtu13*PDstandard4th2gt21 + dgtu133*gtu23*PDstandard4th2gt21 + 
           dgtu233*gtu23*PDstandard4th2gt22 + dgtu333*gtu13*PDstandard4th2gt31 + dgtu133*gtu33*PDstandard4th2gt31 + 
           dgtu333*gtu23*PDstandard4th2gt32 + dgtu233*gtu33*PDstandard4th2gt32 + dgtu333*gtu33*PDstandard4th2gt33) - 
        PDstandard4th23gt11*SQR(gtu13) - PDstandard4th23gt22*SQR(gtu23) - PDstandard4th23gt33*SQR(gtu33);
    
    ddgtu1333  =  -(gtu11*gtu13*PDstandard4th33gt11) - gtu12*gtu13*PDstandard4th33gt21 - gtu11*gtu23*PDstandard4th33gt21 - 
        gtu12*gtu23*PDstandard4th33gt22 - gtu11*gtu33*PDstandard4th33gt31 - gtu13*gtu23*PDstandard4th33gt32 - 
        gtu12*gtu33*PDstandard4th33gt32 - gtu13*gtu33*PDstandard4th33gt33 - dgtu133*gtu11*PDstandard4th3gt11 - 
        dgtu113*gtu13*PDstandard4th3gt11 - dgtu233*gtu11*PDstandard4th3gt21 - dgtu133*gtu12*PDstandard4th3gt21 - 
        dgtu123*gtu13*PDstandard4th3gt21 - dgtu113*gtu23*PDstandard4th3gt21 - dgtu233*gtu12*PDstandard4th3gt22 - 
        dgtu123*gtu23*PDstandard4th3gt22 - dgtu333*gtu11*PDstandard4th3gt31 - 2*dgtu133*gtu13*PDstandard4th3gt31 - 
        dgtu113*gtu33*PDstandard4th3gt31 - dgtu333*gtu12*PDstandard4th3gt32 - dgtu233*gtu13*PDstandard4th3gt32 - 
        dgtu133*gtu23*PDstandard4th3gt32 - dgtu123*gtu33*PDstandard4th3gt32 - dgtu333*gtu13*PDstandard4th3gt33 - 
        dgtu133*gtu33*PDstandard4th3gt33 - PDstandard4th33gt31*SQR(gtu13);
    
    ddgtu2333  =  -(gtu12*gtu13*PDstandard4th33gt11) - gtu13*gtu22*PDstandard4th33gt21 - gtu12*gtu23*PDstandard4th33gt21 - 
        gtu22*gtu23*PDstandard4th33gt22 - gtu13*gtu23*PDstandard4th33gt31 - gtu12*gtu33*PDstandard4th33gt31 - 
        gtu22*gtu33*PDstandard4th33gt32 - gtu23*gtu33*PDstandard4th33gt33 - dgtu133*gtu12*PDstandard4th3gt11 - 
        dgtu123*gtu13*PDstandard4th3gt11 - dgtu233*gtu12*PDstandard4th3gt21 - dgtu223*gtu13*PDstandard4th3gt21 - 
        dgtu133*gtu22*PDstandard4th3gt21 - dgtu123*gtu23*PDstandard4th3gt21 - dgtu233*gtu22*PDstandard4th3gt22 - 
        dgtu223*gtu23*PDstandard4th3gt22 - dgtu333*gtu12*PDstandard4th3gt31 - dgtu233*gtu13*PDstandard4th3gt31 - 
        dgtu133*gtu23*PDstandard4th3gt31 - dgtu123*gtu33*PDstandard4th3gt31 - dgtu333*gtu22*PDstandard4th3gt32 - 
        2*dgtu233*gtu23*PDstandard4th3gt32 - dgtu223*gtu33*PDstandard4th3gt32 - dgtu333*gtu23*PDstandard4th3gt33 - 
        dgtu233*gtu33*PDstandard4th3gt33 - PDstandard4th33gt32*SQR(gtu23);
    
    ddgtu3333  =  -2*(gtu13*gtu23*PDstandard4th33gt21 + gtu13*gtu33*PDstandard4th33gt31 + gtu23*gtu33*PDstandard4th33gt32 + 
           dgtu133*gtu13*PDstandard4th3gt11 + dgtu233*gtu13*PDstandard4th3gt21 + dgtu133*gtu23*PDstandard4th3gt21 + 
           dgtu233*gtu23*PDstandard4th3gt22 + dgtu333*gtu13*PDstandard4th3gt31 + dgtu133*gtu33*PDstandard4th3gt31 + 
           dgtu333*gtu23*PDstandard4th3gt32 + dgtu233*gtu33*PDstandard4th3gt32 + dgtu333*gtu33*PDstandard4th3gt33) - 
        PDstandard4th33gt11*SQR(gtu13) - PDstandard4th33gt22*SQR(gtu23) - PDstandard4th33gt33*SQR(gtu33);
    
    Gt111  =  khalf*(gtu11*PDstandard4th1gt11 + 2*(gtu12*PDstandard4th1gt21 + gtu13*PDstandard4th1gt31) - 
          gtu12*PDstandard4th2gt11 - gtu13*PDstandard4th3gt11);
    
    Gt211  =  khalf*(gtu12*PDstandard4th1gt11 + 2*(gtu22*PDstandard4th1gt21 + gtu23*PDstandard4th1gt31) - 
          gtu22*PDstandard4th2gt11 - gtu23*PDstandard4th3gt11);
    
    Gt311  =  khalf*(gtu13*PDstandard4th1gt11 + 2*(gtu23*PDstandard4th1gt21 + gtu33*PDstandard4th1gt31) - 
          gtu23*PDstandard4th2gt11 - gtu33*PDstandard4th3gt11);
    
    Gt121  =  khalf*(gtu12*PDstandard4th1gt22 + gtu11*PDstandard4th2gt11 + 
          gtu13*(PDstandard4th1gt32 + PDstandard4th2gt31 - PDstandard4th3gt21));
    
    Gt221  =  khalf*(gtu22*PDstandard4th1gt22 + gtu12*PDstandard4th2gt11 + 
          gtu23*(PDstandard4th1gt32 + PDstandard4th2gt31 - PDstandard4th3gt21));
    
    Gt321  =  khalf*(gtu23*PDstandard4th1gt22 + gtu13*PDstandard4th2gt11 + 
          gtu33*(PDstandard4th1gt32 + PDstandard4th2gt31 - PDstandard4th3gt21));
    
    Gt131  =  khalf*(gtu13*PDstandard4th1gt33 + gtu11*PDstandard4th3gt11 + 
          gtu12*(PDstandard4th1gt32 - PDstandard4th2gt31 + PDstandard4th3gt21));
    
    Gt231  =  khalf*(gtu23*PDstandard4th1gt33 + gtu12*PDstandard4th3gt11 + 
          gtu22*(PDstandard4th1gt32 - PDstandard4th2gt31 + PDstandard4th3gt21));
    
    Gt331  =  khalf*(gtu33*PDstandard4th1gt33 + gtu13*PDstandard4th3gt11 + 
          gtu23*(PDstandard4th1gt32 - PDstandard4th2gt31 + PDstandard4th3gt21));
    
    Gt122  =  khalf*(gtu11*(-PDstandard4th1gt22 + 2*PDstandard4th2gt21) + gtu12*PDstandard4th2gt22 + 
          gtu13*(2*PDstandard4th2gt32 - PDstandard4th3gt22));
    
    Gt222  =  khalf*(gtu12*(-PDstandard4th1gt22 + 2*PDstandard4th2gt21) + gtu22*PDstandard4th2gt22 + 
          gtu23*(2*PDstandard4th2gt32 - PDstandard4th3gt22));
    
    Gt322  =  khalf*(gtu13*(-PDstandard4th1gt22 + 2*PDstandard4th2gt21) + gtu23*PDstandard4th2gt22 + 
          gtu33*(2*PDstandard4th2gt32 - PDstandard4th3gt22));
    
    Gt132  =  khalf*(gtu13*PDstandard4th2gt33 + gtu11*(-PDstandard4th1gt32 + PDstandard4th2gt31 + PDstandard4th3gt21) + 
          gtu12*PDstandard4th3gt22);
    
    Gt232  =  khalf*(gtu23*PDstandard4th2gt33 + gtu12*(-PDstandard4th1gt32 + PDstandard4th2gt31 + PDstandard4th3gt21) + 
          gtu22*PDstandard4th3gt22);
    
    Gt332  =  khalf*(gtu33*PDstandard4th2gt33 + gtu13*(-PDstandard4th1gt32 + PDstandard4th2gt31 + PDstandard4th3gt21) + 
          gtu23*PDstandard4th3gt22);
    
    Gt133  =  khalf*(-(gtu11*PDstandard4th1gt33) - gtu12*PDstandard4th2gt33 + 2*gtu11*PDstandard4th3gt31 + 
          2*gtu12*PDstandard4th3gt32 + gtu13*PDstandard4th3gt33);
    
    Gt233  =  khalf*(-(gtu12*PDstandard4th1gt33) - gtu22*PDstandard4th2gt33 + 2*gtu12*PDstandard4th3gt31 + 
          2*gtu22*PDstandard4th3gt32 + gtu23*PDstandard4th3gt33);
    
    Gt333  =  khalf*(-(gtu13*PDstandard4th1gt33) - gtu23*PDstandard4th2gt33 + 2*gtu13*PDstandard4th3gt31 + 
          2*gtu23*PDstandard4th3gt32 + gtu33*PDstandard4th3gt33);
    
    Rt11  =  -(gtu11*khalf*PDstandard4th11gt11) + gtu12*
         (2*Gt211*Gt221*gt22L + 4*gt21L*Gt232*Gt311 + 6*Gt121*Gt311*gt31L + 4*gt11L*Gt131*Gt321 + 4*Gt221*gt31L*Gt321 + 
           4*Gt211*gt31L*Gt322 + 2*Gt221*Gt311*gt32L + 2*Gt211*Gt321*gt32L + 4*gt31L*Gt321*Gt331 + 4*Gt311*gt31L*Gt332 + 
           2*Gt311*Gt321*gt33L - PDstandard4th12gt11) + 
        gtu13*(2*Gt231*Gt311*gt32L + 4*gt11L*Gt131*Gt331 + 2*Gt211*gt32L*Gt331 + 4*Gt211*gt31L*Gt332 + 
           4*Gt311*gt31L*Gt333 + 2*Gt311*Gt331*gt33L - PDstandard4th13gt11) + 2*gt11L*PDstandard4th1Xt1 + 
        gt21L*(6*Gt111*Gt221*gtu12 + 4*Gt211*Gt222*gtu12 + 6*Gt131*Gt211*gtu13 + 6*Gt121*Gt221*gtu22 + 
           6*Gt131*Gt221*gtu23 + 6*Gt131*Gt231*gtu33 + 2*PDstandard4th1Xt2) + 2*gt31L*PDstandard4th1Xt3 - 
        gtu22*khalf*PDstandard4th22gt11 - gtu23*PDstandard4th23gt11 - gtu33*khalf*PDstandard4th33gt11 + 
        Gt111*(10*gt11L*Gt131*gtu13 + 6*gt21L*Gt231*gtu13 + 2*gt11L*Xt1L) + 
        Gt211*(4*gt11L*Gt121*gtu11 + 6*Gt111*gt21L*gtu11 + 4*gt31L*Gt321*gtu11 + 2*Gt311*gt32L*gtu11 + 
           4*gt11L*Gt122*gtu12 + 4*gt11L*Gt132*gtu13 + 2*gt22L*Gt231*gtu13 + 4*gt21L*Gt232*gtu13 + 2*gt21L*Xt1L) + 
        Gt311*(4*gt21L*Gt231*gtu11 + 6*Gt111*gt31L*gtu11 + 4*gt11L*Gt132*gtu12 + 4*gt11L*Gt133*gtu13 + 
           4*gt21L*Gt233*gtu13 + 6*Gt131*gt31L*gtu13 + 2*gt31L*Xt1L) + 
        Gt121*(10*Gt111*gt11L*gtu12 + 6*Gt211*gt21L*gtu12 + 4*gt11L*Gt231*gtu13 + 6*gt31L*Gt321*gtu22 + 
           10*gt11L*Gt131*gtu23 + 2*gt11L*Xt2L) + Gt221*
         (4*gt31L*Gt322*gtu22 + 2*Gt321*gt32L*gtu22 + 2*gt22L*Gt231*gtu23 + 2*gt32L*Gt331*gtu23 + 4*gt31L*Gt332*gtu23 + 
           2*gt21L*Xt2L) + Gt321*(4*gt21L*Gt231*gtu12 + 6*Gt111*gt31L*gtu12 + 4*Gt231*gt31L*gtu13 + 4*gt11L*Gt132*gtu22 + 
           4*gt31L*Gt332*gtu22 + 4*gt11L*Gt133*gtu23 + 4*gt21L*Gt233*gtu23 + 6*Gt131*gt31L*gtu23 + 2*gt31L*Xt2L) + 
        2*gt11L*Gt131*Xt3L + Gt231*(4*gt11L*Gt122*gtu23 + 6*Gt121*gt21L*gtu23 + 4*gt31L*Gt322*gtu23 + 2*Gt321*gt32L*gtu23 + 
           4*gt11L*Gt132*gtu33 + 2*gt32L*Gt331*gtu33 + 4*gt31L*Gt332*gtu33 + 2*gt21L*Xt3L) + 
        Gt331*(4*gt21L*Gt231*gtu13 + 6*Gt111*gt31L*gtu13 + 4*gt11L*Gt132*gtu23 + 4*gt21L*Gt232*gtu23 + 
           6*Gt121*gt31L*gtu23 + 2*Gt321*gt33L*gtu23 + 4*gt11L*Gt133*gtu33 + 4*gt21L*Gt233*gtu33 + 6*Gt131*gt31L*gtu33 + 
           4*gt31L*Gt333*gtu33 + 2*gt31L*Xt3L) + 5*gt11L*gtu11*SQR(Gt111) + 5*gt11L*gtu22*SQR(Gt121) + 
        5*gt11L*gtu33*SQR(Gt131) + gt22L*gtu11*SQR(Gt211) + gt22L*gtu22*SQR(Gt221) + 
        4*(Gt211*gt21L*Gt221*gtu11 + gt11L*Gt131*Gt311*gtu11 + Gt311*gt31L*Gt331*gtu11 + gt11L*Gt121*Gt221*gtu12 + 
           gt21L*Gt221*Gt231*gtu13 + gt11L*Gt122*Gt221*gtu22 + gt21L*Gt221*Gt222*gtu22 + gt21L*Gt232*Gt321*gtu22 + 
           gt11L*Gt132*Gt221*gtu23 + gt21L*Gt222*Gt231*gtu23 + gt21L*Gt221*Gt232*gtu23 + gt31L*Gt331*Gt332*gtu23 + 
           gt31L*Gt321*Gt333*gtu23 + gt21L*Gt231*Gt232*gtu33 + gt21L*gtu12*SQR(Gt221)) + gt22L*gtu33*SQR(Gt231) + 
        gt33L*gtu11*SQR(Gt311) + gt33L*gtu22*SQR(Gt321) + 4*gt31L*gtu13*SQR(Gt331) + gt33L*gtu33*SQR(Gt331);
    
    Rt21  =  gt22L*PDstandard4th1Xt2 + gt32L*PDstandard4th1Xt3 + gt11L*PDstandard4th2Xt1 + 
        gt21L*(PDstandard4th1Xt1 + PDstandard4th2Xt2) + gt31L*PDstandard4th2Xt3 + 
        khalf*(-(gtu11*PDstandard4th11gt21) - 2*gtu12*PDstandard4th12gt21 - 2*gtu13*PDstandard4th13gt21 - 
           gtu22*PDstandard4th22gt21 - 2*gtu23*PDstandard4th23gt21 - gtu33*PDstandard4th33gt21) + 
        (gt11L*Gt121 + gt21L*(Gt111 + Gt221) + Gt211*gt22L + gt31L*Gt321 + Gt311*gt32L)*Xt1L + 
        (gt11L*Gt122 + gt21L*(Gt121 + Gt222) + Gt221*gt22L + gt31L*Gt322 + Gt321*gt32L)*Xt2L + 
        (gt11L*Gt132 + gt22L*Gt231 + gt21L*(Gt131 + Gt232) + gt32L*Gt331 + gt31L*Gt332)*Xt3L + 
        gtu11*(gt11L*(Gt121*(3*Gt111 + 2*Gt221) + 2*Gt131*Gt321) + 
           Gt111*(gt21L*Gt221 + gt31L*Gt321 + 2*(Gt211*gt22L + Gt311*gt32L)) + 
           Gt211*(5*Gt121*gt21L + 3*(Gt221*gt22L + Gt321*gt32L)) + gt31L*(3*Gt121*Gt311 + 2*Gt321*(Gt221 + Gt331)) + 
           Gt311*(Gt221*gt32L + Gt321*gt33L) + 2*(Gt311*(Gt131*gt21L + gt22L*Gt231 + gt32L*Gt331) + 
              gt21L*(Gt231*Gt321 + SQR(Gt111) + SQR(Gt221)))) + 
        gtu22*(gt11L*(Gt122*(3*Gt121 + 2*Gt222) + 2*Gt132*Gt322) + 
           Gt121*(gt21L*Gt222 + gt31L*Gt322 + 2*(Gt221*gt22L + Gt321*gt32L)) + 
           Gt221*(5*Gt122*gt21L + 3*(Gt222*gt22L + Gt322*gt32L)) + gt31L*(3*Gt122*Gt321 + 2*Gt322*(Gt222 + Gt332)) + 
           Gt321*(Gt222*gt32L + Gt322*gt33L) + 2*(Gt321*(Gt132*gt21L + gt22L*Gt232 + gt32L*Gt332) + 
              gt21L*(Gt232*Gt322 + SQR(Gt121) + SQR(Gt222)))) + 
        gtu33*(gt11L*(Gt132*(3*Gt131 + 2*Gt232) + 2*Gt133*Gt332) + 
           Gt131*(gt21L*Gt232 + 2*(gt22L*Gt231 + gt32L*Gt331) + gt31L*Gt332) + 
           Gt231*(5*Gt132*gt21L + 3*(gt22L*Gt232 + gt32L*Gt332)) + 
           Gt331*(3*Gt132*gt31L + gt32L*(Gt232 + 2*Gt333) + Gt332*gt33L) + 
           2*((Gt133*gt21L + gt22L*Gt233)*Gt331 + Gt332*(gt21L*Gt233 + gt31L*(Gt232 + Gt333)) + 
              gt21L*(SQR(Gt131) + SQR(Gt232)))) + gtu12*
         (Gt122*(Gt111*gt11L + 3*Gt211*gt21L + Gt311*gt31L) + 
           gt21L*(Gt111*Gt222 + Gt221*(4*Gt121 + 2*Gt222) + 2*(Gt132*Gt311 + Gt232*Gt321)) + 
           gt31L*((Gt111 + 2*Gt221)*Gt322 + Gt321*(4*Gt121 + 2*Gt332)) + 
           gt32L*(2*Gt221*Gt321 + 3*Gt211*Gt322 + Gt311*(Gt222 + 2*Gt332)) + 
           gt11L*(2*(Gt122*Gt221 + Gt132*Gt321) + 3*SQR(Gt121)) + gt22L*(3*Gt211*Gt222 + 2*Gt232*Gt311 + SQR(Gt221)) + 
           2*(gt11L*(Gt111*Gt122 + Gt121*Gt222 + Gt131*Gt322) + 
              gt21L*(Gt122*Gt211 + Gt121*(Gt111 + Gt221) + Gt221*Gt222 + Gt131*Gt321 + Gt231*Gt322) + 
              Gt111*(Gt121*gt21L + Gt221*gt22L + Gt321*gt32L) + gt31L*(Gt122*Gt311 + Gt222*Gt321 + Gt322*Gt331) + 
              gt32L*(Gt121*Gt311 + Gt321*(Gt221 + Gt331)) + gt22L*(Gt121*Gt211 + Gt231*Gt321 + SQR(Gt221))) + 
           gt33L*(Gt311*Gt322 + SQR(Gt321))) + gtu13*
         (gt22L*(3*Gt211*Gt232 + 2*Gt233*Gt311) + Gt132*(Gt111*gt11L + 3*Gt211*gt21L + Gt311*gt31L) + 
           gt11L*(3*Gt121*Gt131 + 2*(Gt132*Gt221 + Gt133*Gt321)) + 
           gt21L*(Gt131*Gt221 + 3*Gt121*Gt231 + Gt111*Gt232 + 2*(Gt133*Gt311 + Gt233*Gt321)) + 
           Gt221*(gt22L*Gt231 + 2*gt21L*Gt232 + gt32L*Gt331) + gt31L*(3*Gt121*Gt331 + Gt111*Gt332) + 
           gt32L*(3*Gt211*Gt332 + Gt311*(Gt232 + 2*Gt333)) + Gt311*Gt332*gt33L + 
           Gt321*(Gt231*gt32L + gt31L*(Gt131 + 2*Gt333) + Gt331*gt33L) + 
           2*(Gt111*(Gt131*gt21L + gt22L*Gt231 + gt32L*Gt331) + gt22L*(Gt131*Gt211 + Gt231*(Gt221 + Gt331)) + 
              gt11L*(Gt111*Gt132 + Gt121*Gt232 + Gt131*Gt332) + gt31L*(Gt132*Gt311 + Gt232*Gt321 + (Gt221 + Gt331)*Gt332) + 
              gt21L*(Gt132*Gt211 + Gt221*Gt232 + Gt131*(Gt111 + Gt331) + Gt231*(Gt121 + Gt332)) + 
              gt32L*(Gt131*Gt311 + Gt231*Gt321 + SQR(Gt331)))) + 
        gtu23*(gt22L*(3*Gt221*Gt232 + 2*Gt233*Gt321) + Gt132*(gt11L*Gt121 + 3*gt21L*Gt221 + gt31L*Gt321) + 
           gt11L*(3*Gt122*Gt131 + 2*(Gt132*Gt222 + Gt133*Gt322)) + 
           gt21L*(Gt131*Gt222 + 3*Gt122*Gt231 + Gt121*Gt232 + 2*(Gt133*Gt321 + Gt233*Gt322)) + 
           Gt222*(gt22L*Gt231 + 2*gt21L*Gt232 + gt32L*Gt331) + gt31L*(3*Gt122*Gt331 + Gt121*Gt332) + 
           gt32L*(3*Gt221*Gt332 + Gt321*(Gt232 + 2*Gt333)) + Gt321*Gt332*gt33L + 
           Gt322*(Gt231*gt32L + gt31L*(Gt131 + 2*Gt333) + Gt331*gt33L) + 
           2*(Gt131*(Gt121*gt21L + Gt221*gt22L + Gt321*gt32L) + Gt231*(Gt122*gt21L + Gt222*gt22L + Gt322*gt32L) + 
              Gt121*(Gt131*gt21L + gt22L*Gt231 + gt32L*Gt331) + Gt331*(Gt132*gt21L + gt22L*Gt232 + gt32L*Gt332) + 
              gt11L*(Gt122*Gt232 + Gt132*(Gt121 + Gt332)) + gt21L*(Gt132*Gt221 + Gt232*(Gt222 + Gt332)) + 
              gt31L*(Gt132*Gt321 + Gt232*Gt322 + Gt222*Gt332 + SQR(Gt332))));
    
    Rt31  =  gt32L*PDstandard4th1Xt2 + gt33L*PDstandard4th1Xt3 + 
        khalf*(-(gtu11*PDstandard4th11gt31) - 2*gtu12*PDstandard4th12gt31 - 2*gtu13*PDstandard4th13gt31 - 
           gtu22*PDstandard4th22gt31 - 2*gtu23*PDstandard4th23gt31 - gtu33*PDstandard4th33gt31) + gt11L*PDstandard4th3Xt1 + 
        gt21L*PDstandard4th3Xt2 + gt31L*(PDstandard4th1Xt1 + PDstandard4th3Xt3) + 
        (gt11L*Gt131 + gt21L*Gt231 + Gt211*gt32L + gt31L*(Gt111 + Gt331) + Gt311*gt33L)*Xt1L + 
        (gt11L*Gt132 + gt21L*Gt232 + Gt221*gt32L + gt31L*(Gt121 + Gt332) + Gt321*gt33L)*Xt2L + 
        (gt11L*Gt133 + gt21L*Gt233 + Gt231*gt32L + gt31L*(Gt131 + Gt333) + Gt331*gt33L)*Xt3L + 
        gtu12*(gt22L*(Gt221*Gt231 + Gt211*Gt232) + Gt132*(Gt111*gt11L + Gt211*gt21L + 3*Gt311*gt31L) + 
           gt21L*((Gt121 + 2*Gt222)*Gt231 + Gt232*(Gt111 + 2*Gt331)) + 
           gt11L*(3*Gt121*Gt131 + 2*(Gt122*Gt231 + Gt132*Gt331)) + gt31L*(Gt121*Gt331 + Gt111*Gt332) + 
           gt32L*(Gt231*Gt321 + Gt221*Gt331 + Gt211*(2*Gt222 + Gt332)) + (2*Gt211*Gt322 + Gt321*Gt331)*gt33L + 
           3*(Gt131*(gt21L*Gt221 + gt31L*Gt321) + Gt311*(Gt232*gt32L + Gt332*gt33L)) + 
           2*(gt11L*(Gt111*Gt132 + Gt121*Gt232 + Gt131*Gt332) + gt21L*(Gt132*Gt211 + Gt221*Gt232 + Gt231*Gt332) + 
              gt31L*(Gt122*Gt211 + Gt121*(2*Gt111 + Gt221) + Gt132*Gt311 + (Gt131 + Gt232)*Gt321 + Gt231*Gt322 + 
                 2*Gt331*Gt332) + (Gt121*Gt311 + Gt321*(Gt221 + Gt331))*gt33L + Gt111*(Gt221*gt32L + Gt321*gt33L) + 
              gt32L*(Gt121*Gt211 + Gt231*Gt321 + SQR(Gt221)))) + 
        gtu23*(Gt233*(Gt221*gt22L + 2*gt31L*Gt322 + Gt321*gt32L) + Gt132*(gt11L*Gt131 + gt21L*Gt231 + 3*gt31L*Gt331) + 
           Gt232*(Gt131*gt21L + gt22L*Gt231 + 3*gt32L*Gt331) + (Gt131*gt31L + Gt231*gt32L)*Gt332 + 
           gt21L*((Gt121 + 2*Gt222)*Gt233 + 2*Gt232*Gt333) + gt11L*(3*Gt121*Gt133 + 2*(Gt122*Gt233 + Gt132*Gt333)) + 
           Gt333*(Gt121*gt31L + Gt221*gt32L + Gt321*gt33L) + 3*(Gt133*(gt21L*Gt221 + gt31L*Gt321) + Gt331*Gt332*gt33L) + 
           2*((Gt222*Gt231 + Gt221*(Gt131 + Gt232) + Gt233*Gt321)*gt32L + gt11L*(Gt132*(Gt131 + Gt232) + Gt133*Gt332) + 
              gt31L*(Gt121*Gt131 + Gt122*Gt231 + Gt133*Gt321 + Gt132*(Gt221 + Gt331) + Gt332*(Gt232 + 2*Gt333)) + 
              (Gt231*Gt322 + Gt221*Gt332 + Gt321*(Gt131 + Gt333))*gt33L + Gt121*(Gt131*gt31L + Gt231*gt32L + Gt331*gt33L) + 
              gt21L*(Gt132*Gt231 + Gt233*Gt332 + SQR(Gt232)))) + 
        gtu11*(Gt231*(Gt211*gt22L + 2*(gt21L*Gt221 + gt31L*Gt321) + 3*Gt311*gt32L) + 
           gt21L*(3*Gt131*Gt211 + Gt231*(Gt111 + 2*Gt331)) + gt11L*(2*Gt121*Gt231 + Gt131*(3*Gt111 + 2*Gt331)) + 
           3*Gt311*Gt331*gt33L + Gt211*(gt32L*Gt331 + 2*Gt321*gt33L) + 
           gt31L*(5*Gt131*Gt311 + Gt111*Gt331 + 2*(Gt121*Gt211 + SQR(Gt111))) + 
           2*(Gt211*(Gt111 + Gt221)*gt32L + Gt111*Gt311*gt33L + gt31L*SQR(Gt331))) + 
        gtu13*(Gt133*(Gt111*gt11L + Gt211*gt21L + 3*Gt311*gt31L) + gt32L*(3*Gt233*Gt311 + 2*Gt231*(Gt111 + Gt331)) + 
           gt21L*(Gt111*Gt233 + 2*(Gt231*Gt232 + Gt233*Gt331)) + Gt111*(gt31L*Gt333 + 2*Gt331*gt33L) + 
           Gt211*(gt22L*Gt233 + gt32L*Gt333 + 2*(Gt132*gt31L + Gt232*gt32L + Gt332*gt33L)) + 
           gt11L*(2*(Gt132*Gt231 + Gt133*Gt331) + 3*SQR(Gt131)) + gt22L*SQR(Gt231) + gt33L*(3*Gt311*Gt333 + SQR(Gt331)) + 
           2*(gt32L*(Gt131*Gt211 + Gt231*(Gt221 + Gt331)) + Gt131*(Gt111*gt31L + 2*(gt21L*Gt231 + gt31L*Gt331)) + 
              gt11L*(Gt111*Gt133 + Gt121*Gt233 + Gt131*Gt333) + gt21L*(Gt133*Gt211 + Gt221*Gt233 + Gt231*Gt333) + 
              gt31L*(Gt133*Gt311 + Gt233*Gt321 + Gt131*(Gt111 + Gt331) + Gt231*(Gt121 + Gt332) + 2*Gt331*Gt333) + 
              gt33L*(Gt131*Gt311 + Gt231*Gt321 + SQR(Gt331)))) + 
        gtu22*(Gt232*(Gt221*gt22L + 2*(gt21L*Gt222 + gt31L*Gt322) + 3*Gt321*gt32L) + 
           gt21L*(3*Gt132*Gt221 + Gt232*(Gt121 + 2*Gt332)) + gt11L*(2*Gt122*Gt232 + Gt132*(3*Gt121 + 2*Gt332)) + 
           3*Gt321*Gt332*gt33L + Gt221*(gt32L*Gt332 + 2*Gt322*gt33L) + 
           gt31L*(5*Gt132*Gt321 + Gt121*Gt332 + 2*(Gt122*Gt221 + SQR(Gt121))) + 
           2*(Gt221*(Gt121 + Gt222)*gt32L + Gt121*Gt321*gt33L + gt31L*SQR(Gt332))) + 
        gtu33*(Gt233*(gt22L*Gt231 + 2*gt21L*Gt232 + 3*gt32L*Gt331) + gt21L*(3*Gt133*Gt231 + Gt233*(Gt131 + 2*Gt333)) + 
           gt11L*(2*Gt132*Gt233 + Gt133*(3*Gt131 + 2*Gt333)) + 3*Gt331*Gt333*gt33L + Gt231*(gt32L*Gt333 + 2*Gt332*gt33L) + 
           gt31L*(5*Gt133*Gt331 + Gt131*Gt333 + 2*SQR(Gt131)) + 
           2*(Gt231*(Gt132*gt31L + (Gt131 + Gt232)*gt32L) + Gt131*Gt331*gt33L + gt31L*(Gt233*Gt332 + SQR(Gt333))));
    
    Rt22  =  6*(Gt122*gt21L*Gt221*gtu12 + Gt121*gt21L*Gt222*gtu12 + Gt222*Gt321*gt32L*gtu12 + Gt221*Gt322*gt32L*gtu12 + 
           Gt132*gt21L*Gt221*gtu13 + Gt122*gt21L*Gt222*gtu22 + Gt132*gt21L*Gt222*gtu23 + Gt232*gt32L*Gt332*gtu33) - 
        gtu11*khalf*PDstandard4th11gt22 + gtu12*(2*Gt122*gt31L*Gt321 + 4*Gt131*gt21L*Gt322 + 4*Gt121*Gt321*gt32L + 
           4*Gt322*gt32L*Gt331 + 4*Gt321*gt32L*Gt332 + 2*Gt321*Gt322*gt33L - PDstandard4th12gt22) - 
        gtu13*PDstandard4th13gt22 - gtu22*khalf*PDstandard4th22gt22 + 
        gtu23*(2*Gt132*gt31L*Gt322 + 4*Gt122*gt32L*Gt331 + 4*gt22L*Gt232*Gt332 + 4*Gt322*gt32L*Gt333 + 
           2*Gt322*Gt332*gt33L - PDstandard4th23gt22) + 2*gt21L*PDstandard4th2Xt1 + 
        gt22L*(10*Gt221*Gt222*gtu12 + 4*Gt232*Gt321*gtu12 + 4*Gt132*Gt211*gtu13 + 4*Gt122*Gt221*gtu22 + 
           4*Gt132*Gt221*gtu23 + 10*Gt222*Gt232*gtu23 + 4*Gt233*Gt322*gtu23 + 4*Gt132*Gt231*gtu33 + 2*PDstandard4th2Xt2) + 
        gt32L*(4*Gt121*Gt311*gtu11 + 6*Gt221*Gt321*gtu11 + 4*Gt122*Gt311*gtu12 + 4*Gt132*Gt311*gtu13 + 
           6*Gt232*Gt321*gtu13 + 4*Gt121*Gt331*gtu13 + 4*Gt122*Gt321*gtu22 + 6*Gt232*Gt322*gtu23 + 4*Gt132*Gt331*gtu33 + 
           2*PDstandard4th2Xt3) - gtu33*khalf*PDstandard4th33gt22 + 2*Gt221*gt22L*Xt1L + 
        Gt121*(4*Gt111*gt21L*gtu11 + 6*gt21L*Gt221*gtu11 + 2*gt11L*Gt122*gtu12 + 2*gt11L*Gt132*gtu13 + 
           4*Gt122*gt21L*gtu22 + 4*Gt132*gt21L*gtu23 + 2*gt21L*Xt1L) + 
        Gt321*(4*gt22L*Gt231*gtu11 + 2*Gt121*gt31L*gtu11 + 4*gt32L*Gt331*gtu11 + 4*Gt132*gt21L*gtu12 + 
           4*Gt133*gt21L*gtu13 + 4*gt22L*Gt233*gtu13 + 2*Gt132*gt31L*gtu13 + 2*Gt332*gt33L*gtu13 + 4*Gt132*gt32L*gtu23 + 
           2*gt32L*Xt1L) + 2*Gt222*gt22L*Xt2L + Gt122*(2*gt11L*Gt132*gtu23 + 4*gt22L*Gt231*gtu23 + 2*gt21L*Xt2L) + 
        Gt322*(4*gt22L*Gt231*gtu12 + 2*Gt121*gt31L*gtu12 + 4*Gt132*gt21L*gtu22 + 2*Gt122*gt31L*gtu22 + 
           6*Gt222*gt32L*gtu22 + 4*Gt133*gt21L*gtu23 + 2*gt32L*Xt2L) + 2*Gt132*gt21L*Xt3L + 
        Gt232*(6*Gt121*gt21L*gtu13 + 10*Gt221*gt22L*gtu13 + 4*gt22L*Gt322*gtu22 + 6*Gt122*gt21L*gtu23 + 
           6*Gt132*gt21L*gtu33 + 2*gt22L*Xt3L) + Gt332*
         (2*Gt121*gt31L*gtu13 + 6*Gt221*gt32L*gtu13 + 4*Gt322*gt32L*gtu22 + 4*Gt132*gt21L*gtu23 + 2*Gt122*gt31L*gtu23 + 
           6*Gt222*gt32L*gtu23 + 4*Gt133*gt21L*gtu33 + 4*gt22L*Gt233*gtu33 + 2*Gt132*gt31L*gtu33 + 2*gt32L*Xt3L) + 
        gt11L*gtu11*SQR(Gt121) + 4*(Gt121*Gt211*gt22L*gtu11 + Gt131*gt21L*Gt321*gtu11 + Gt111*Gt122*gt21L*gtu12 + 
           Gt122*Gt211*gt22L*gtu12 + Gt121*Gt221*gt22L*gtu12 + Gt121*Gt131*gt21L*gtu13 + Gt111*Gt132*gt21L*gtu13 + 
           Gt121*gt22L*Gt231*gtu13 + Gt131*gt21L*Gt332*gtu13 + gt22L*Gt231*Gt332*gtu13 + gt32L*Gt331*Gt332*gtu13 + 
           Gt321*gt32L*Gt333*gtu13 + Gt122*Gt131*gt21L*gtu23 + Gt131*Gt132*gt21L*gtu33 + gt32L*Gt332*Gt333*gtu33 + 
           gt21L*gtu12*SQR(Gt121)) + gt11L*gtu22*SQR(Gt122) + gt11L*gtu33*SQR(Gt132) + 5*gt22L*gtu11*SQR(Gt221) + 
        5*gt22L*gtu22*SQR(Gt222) + 5*gt22L*gtu33*SQR(Gt232) + gt33L*gtu11*SQR(Gt321) + gt33L*gtu22*SQR(Gt322) + 
        4*gt32L*gtu23*SQR(Gt332) + gt33L*gtu33*SQR(Gt332);
    
    Rt32  =  gt31L*PDstandard4th2Xt1 + gt33L*PDstandard4th2Xt3 + 
        khalf*(-(gtu11*PDstandard4th11gt32) - 2*gtu12*PDstandard4th12gt32 - 2*gtu13*PDstandard4th13gt32 - 
           gtu22*PDstandard4th22gt32 - 2*gtu23*PDstandard4th23gt32 - gtu33*PDstandard4th33gt32) + gt21L*PDstandard4th3Xt1 + 
        gt22L*PDstandard4th3Xt2 + gt32L*(PDstandard4th2Xt2 + PDstandard4th3Xt3) + 
        (Gt131*gt21L + gt22L*Gt231 + Gt121*gt31L + gt32L*(Gt221 + Gt331) + Gt321*gt33L)*Xt1L + 
        (Gt132*gt21L + gt22L*Gt232 + Gt122*gt31L + gt32L*(Gt222 + Gt332) + Gt322*gt33L)*Xt2L + 
        (Gt133*gt21L + gt22L*Gt233 + Gt132*gt31L + gt32L*(Gt232 + Gt333) + Gt332*gt33L)*Xt3L + 
        gtu12*(gt11L*(Gt122*Gt131 + Gt121*Gt132) + gt32L*
            (2*(Gt132*Gt311 + Gt131*Gt321) + 5*(Gt232*Gt321 + Gt231*Gt322) + Gt222*(4*Gt221 + Gt331) + 
              (Gt221 + 4*Gt331)*Gt332) + gt21L*(3*Gt122*Gt231 + Gt132*(Gt221 + 2*Gt331) + Gt131*(Gt222 + 2*Gt332)) + 
           3*(Gt121*gt21L*Gt232 + gt22L*(Gt222*Gt231 + Gt221*Gt232) + gt31L*(Gt132*Gt321 + Gt131*Gt322) + 
              (Gt322*Gt331 + Gt321*Gt332)*gt33L) + 2*
            ((Gt122*Gt221 + Gt121*Gt222)*gt31L + Gt111*(Gt132*gt21L + Gt122*gt31L) + 
              gt22L*(Gt132*Gt211 + Gt131*Gt221 + Gt232*Gt331 + Gt231*Gt332) + (Gt222*Gt321 + Gt221*Gt322)*gt33L + 
              Gt122*(Gt211*gt32L + Gt311*gt33L) + Gt121*(Gt131*gt21L + Gt221*gt32L + Gt321*gt33L)) + 
           gt31L*(Gt122*Gt331 + Gt121*Gt332 + 2*SQR(Gt121))) + 
        gtu13*(gt11L*(Gt131*Gt132 + Gt121*Gt133) + Gt133*(gt21L*Gt221 + 3*gt31L*Gt321 + 2*(Gt211*gt22L + Gt311*gt32L)) + 
           gt32L*(5*Gt233*Gt321 + Gt232*(4*Gt221 + Gt331) + Gt231*(2*Gt121 + 5*Gt332) + Gt221*Gt333) + 
           Gt331*(2*gt22L*Gt233 + Gt132*gt31L + 4*gt32L*Gt333) + 
           Gt131*(2*(gt22L*Gt231 + Gt121*gt31L + gt32L*Gt331) + gt31L*Gt332 + gt21L*(Gt232 + 2*Gt333)) + 
           Gt121*(gt31L*Gt333 + 2*Gt331*gt33L) + 3*(Gt231*(Gt132*gt21L + gt22L*Gt232) + (Gt121*gt21L + Gt221*gt22L)*Gt233 + 
              (Gt331*Gt332 + Gt321*Gt333)*gt33L) + 2*
            (Gt111*(Gt133*gt21L + Gt132*gt31L) + gt31L*(Gt132*Gt221 + Gt121*Gt232 + Gt131*Gt332) + gt22L*Gt231*Gt333 + 
              (Gt232*Gt321 + Gt221*Gt332)*gt33L + Gt132*(Gt211*gt32L + Gt311*gt33L) + gt21L*(Gt133*Gt331 + SQR(Gt131)))) + 
        gtu11*(Gt131*(gt11L*Gt121 + gt21L*Gt221 + 3*gt31L*Gt321 + 2*(Gt211*gt22L + Gt311*gt32L)) + 
           (Gt121*gt31L + Gt221*gt32L)*Gt331 + Gt231*(5*Gt321*gt32L + 2*gt22L*Gt331) + 
           3*((Gt121*gt21L + Gt221*gt22L)*Gt231 + Gt321*Gt331*gt33L) + 
           2*(Gt111*(Gt131*gt21L + Gt121*gt31L) + Gt131*gt21L*Gt331 + Gt221*Gt321*gt33L + 
              Gt121*(Gt221*gt31L + Gt211*gt32L + Gt311*gt33L) + gt32L*(SQR(Gt221) + SQR(Gt331)))) + 
        gtu23*(Gt233*(Gt122*gt21L + Gt222*gt22L + 3*Gt322*gt32L) + 
           Gt133*(gt11L*Gt122 + 3*gt31L*Gt322 + gt21L*(Gt222 + 2*Gt332)) + 
           Gt132*(4*gt21L*Gt232 + 2*(gt22L*Gt231 + gt32L*Gt331 + gt31L*(Gt222 + Gt332))) + 
           gt32L*(4*Gt232*Gt332 + Gt222*Gt333) + Gt122*(gt31L*Gt333 + 2*Gt331*gt33L) + gt11L*SQR(Gt132) + 
           3*(Gt322*Gt333*gt33L + gt22L*SQR(Gt232)) + gt33L*SQR(Gt332) + 
           2*(Gt131*(Gt132*gt21L + Gt122*gt31L) + Gt133*(Gt121*gt21L + Gt221*gt22L + Gt321*gt32L) + 
              gt31L*(Gt122*Gt232 + Gt132*(Gt121 + Gt332)) + Gt233*(Gt122*gt21L + Gt322*gt32L + gt22L*(Gt222 + Gt332)) + 
              (Gt132*gt21L + gt22L*Gt232)*Gt333 + gt32L*
               (Gt132*Gt221 + Gt122*Gt231 + Gt232*(2*Gt222 + Gt332) + 2*Gt332*Gt333) + 
              gt33L*(Gt132*Gt321 + Gt232*Gt322 + Gt222*Gt332 + SQR(Gt332)))) + 
        gtu22*(Gt132*(gt11L*Gt122 + gt21L*Gt222 + 3*gt31L*Gt322 + 2*(Gt221*gt22L + Gt321*gt32L)) + 
           (Gt122*gt31L + Gt222*gt32L)*Gt332 + Gt232*(5*Gt322*gt32L + 2*gt22L*Gt332) + 
           3*((Gt122*gt21L + Gt222*gt22L)*Gt232 + Gt322*Gt332*gt33L) + 
           2*(Gt121*(Gt132*gt21L + Gt122*gt31L) + Gt132*gt21L*Gt332 + Gt222*Gt322*gt33L + 
              Gt122*(Gt222*gt31L + Gt221*gt32L + Gt321*gt33L) + gt32L*(SQR(Gt222) + SQR(Gt332)))) + 
        gtu33*(Gt133*(gt11L*Gt132 + gt21L*Gt232 + 2*(gt22L*Gt231 + gt32L*Gt331)) + (Gt132*gt31L + Gt232*gt32L)*Gt333 + 
           Gt233*(5*gt32L*Gt332 + 2*gt22L*Gt333) + 3*
            ((Gt132*gt21L + gt22L*Gt232)*Gt233 + Gt332*(Gt133*gt31L + Gt333*gt33L)) + 
           2*(Gt131*(Gt133*gt21L + Gt132*gt31L) + Gt133*gt21L*Gt333 + Gt232*Gt332*gt33L + 
              Gt132*(Gt232*gt31L + Gt231*gt32L + Gt331*gt33L) + gt32L*(SQR(Gt232) + SQR(Gt333))));
    
    Rt33  =  6*(Gt133*gt31L*Gt331*gtu13 + Gt233*gt32L*Gt331*gtu13 + Gt131*gt31L*Gt333*gtu13 + Gt231*gt32L*Gt333*gtu13 + 
           Gt132*gt31L*Gt332*gtu22 + Gt133*gt31L*Gt332*gtu23 + Gt132*gt31L*Gt333*gtu23 + Gt232*gt32L*Gt333*gtu23 + 
           Gt133*gt31L*Gt333*gtu33) + gtu12*(2*gt22L*Gt231*Gt232 + 4*Gt111*Gt132*gt31L + 4*Gt132*Gt211*gt32L + 
           4*Gt221*Gt232*gt32L + 6*Gt132*gt31L*Gt331 + 6*Gt232*gt32L*Gt331 + 6*Gt131*gt31L*Gt332 + 6*Gt231*gt32L*Gt332 + 
           4*Gt132*Gt311*gt33L - PDstandard4th12gt33) - gtu13*PDstandard4th13gt33 - gtu22*khalf*PDstandard4th22gt33 - 
        gtu23*PDstandard4th23gt33 - gtu33*khalf*PDstandard4th33gt33 + 
        2*(Gt132*gt21L*Gt231*gtu12 + Gt131*gt21L*Gt232*gtu12 + Gt133*gt21L*Gt231*gtu13 + Gt131*gt21L*Gt233*gtu13 + 
           gt22L*Gt231*Gt233*gtu13 + Gt132*gt21L*Gt232*gtu22 + gt11L*Gt132*Gt133*gtu23 + Gt132*gt21L*Gt233*gtu23 + 
           gt22L*Gt232*Gt233*gtu23 + Gt133*gt21L*Gt233*gtu33 + gt31L*PDstandard4th3Xt1) + 2*gt32L*PDstandard4th3Xt2 + 
        gt33L*(4*Gt231*Gt322*gtu12 + 10*Gt331*Gt332*gtu12 + 4*Gt133*Gt311*gtu13 + 4*Gt231*Gt332*gtu13 + 
           10*Gt331*Gt333*gtu13 + 4*Gt132*Gt321*gtu22 + 4*Gt133*Gt321*gtu23 + 4*Gt133*Gt331*gtu33 + 2*PDstandard4th3Xt3) + 
        2*Gt231*gt32L*Xt1L + 2*Gt331*gt33L*Xt1L + Gt131*
         (2*gt21L*Gt231*gtu11 + 4*Gt111*gt31L*gtu11 + 2*gt11L*Gt132*gtu12 + 2*gt11L*Gt133*gtu13 + 4*Gt132*gt31L*gtu23 + 
           4*Gt133*gt31L*gtu33 + 2*gt31L*Xt1L) + 2*Gt132*gt31L*Xt2L + 
        Gt232*(4*Gt222*gt32L*gtu22 + 6*gt32L*Gt332*gtu22 + 2*Gt133*gt21L*gtu23 + 4*Gt233*gt32L*gtu33 + 2*gt32L*Xt2L) + 
        Gt332*(4*Gt232*gt33L*gtu23 + 10*Gt333*gt33L*gtu23 + 4*Gt233*gt33L*gtu33 + 2*gt33L*Xt2L) + 2*Gt133*gt31L*Xt3L + 
        2*Gt333*gt33L*Xt3L + Gt233*(4*Gt222*gt32L*gtu23 + 6*gt32L*Gt332*gtu23 + 4*Gt132*gt31L*gtu33 + 6*gt32L*Gt333*gtu33 + 
           2*gt32L*Xt3L) + gtu11*(4*Gt221*Gt231*gt32L + 6*Gt131*gt31L*Gt331 + 6*Gt231*gt32L*Gt331 + 4*Gt131*Gt311*gt33L + 
           4*Gt231*Gt321*gt33L - khalf*PDstandard4th11gt33 + gt11L*SQR(Gt131)) + 
        4*(Gt121*Gt231*gt31L*gtu11 + Gt131*Gt211*gt32L*gtu11 + Gt121*Gt131*gt31L*gtu12 + Gt122*Gt231*gt31L*gtu12 + 
           Gt121*Gt232*gt31L*gtu12 + Gt131*Gt221*gt32L*gtu12 + Gt222*Gt231*gt32L*gtu12 + Gt131*Gt321*gt33L*gtu12 + 
           Gt232*Gt321*gt33L*gtu12 + Gt111*Gt133*gt31L*gtu13 + Gt132*Gt231*gt31L*gtu13 + Gt121*Gt233*gt31L*gtu13 + 
           Gt133*Gt211*gt32L*gtu13 + Gt131*Gt231*gt32L*gtu13 + Gt231*Gt232*gt32L*gtu13 + Gt221*Gt233*gt32L*gtu13 + 
           Gt233*Gt321*gt33L*gtu13 + Gt131*Gt331*gt33L*gtu13 + Gt121*Gt132*gt31L*gtu22 + Gt122*Gt232*gt31L*gtu22 + 
           Gt132*Gt221*gt32L*gtu22 + Gt232*Gt322*gt33L*gtu22 + Gt121*Gt133*gt31L*gtu23 + Gt132*Gt232*gt31L*gtu23 + 
           Gt122*Gt233*gt31L*gtu23 + Gt133*Gt221*gt32L*gtu23 + Gt132*Gt231*gt32L*gtu23 + Gt233*Gt322*gt33L*gtu23 + 
           Gt132*Gt331*gt33L*gtu23 + Gt133*Gt231*gt32L*gtu33 + gt31L*gtu13*SQR(Gt131)) + gt11L*gtu22*SQR(Gt132) + 
        gt11L*gtu33*SQR(Gt133) + gt22L*gtu11*SQR(Gt231) + gt22L*gtu22*SQR(Gt232) + 4*gt32L*gtu23*SQR(Gt232) + 
        gt22L*gtu33*SQR(Gt233) + 5*gt33L*gtu11*SQR(Gt331) + 5*gt33L*gtu22*SQR(Gt332) + 5*gt33L*gtu33*SQR(Gt333);
    
    Rphi11  =  2*(-PDstandard4th11phi - gt11L*gtu11*PDstandard4th11phi - 2*gt11L*gtu12*PDstandard4th12phi - 
          2*gt11L*gtu13*PDstandard4th13phi - gt11L*gtu22*PDstandard4th22phi - 2*gt11L*gtu23*PDstandard4th23phi - 
          gt11L*gtu33*PDstandard4th33phi + Gt311*PDstandard4th3phi + gt11L*Gt311*gtu11*PDstandard4th3phi + 
          2*gt11L*Gt321*gtu12*PDstandard4th3phi + 2*gt11L*Gt331*gtu13*PDstandard4th3phi + 
          gt11L*Gt322*gtu22*PDstandard4th3phi + 2*gt11L*Gt332*gtu23*PDstandard4th3phi + 
          gt11L*Gt333*gtu33*PDstandard4th3phi + PDstandard4th1phi*
           (Gt111 + Gt111*gt11L*gtu11 + 2*gt11L*Gt121*gtu12 + 2*gt11L*Gt131*gtu13 + gt11L*Gt122*gtu22 + 
             2*gt11L*Gt132*gtu23 + gt11L*Gt133*gtu33 - 4*gt11L*gtu12*PDstandard4th2phi - 4*gt11L*gtu13*PDstandard4th3phi) + 
          PDstandard4th2phi*(Gt211 + gt11L*Gt211*gtu11 + 
             gt11L*(2*Gt221*gtu12 + 2*Gt231*gtu13 + Gt222*gtu22 + 2*Gt232*gtu23 + Gt233*gtu33) - 
             4*gt11L*gtu23*PDstandard4th3phi) + (2 - 2*gt11L*gtu11)*SQR(PDstandard4th1phi) - 
          2*gt11L*gtu22*SQR(PDstandard4th2phi) - 2*gt11L*gtu33*SQR(PDstandard4th3phi));
    
    Rphi21  =  2*(-(gt21L*gtu11*PDstandard4th11phi) - PDstandard4th12phi - 2*gt21L*gtu12*PDstandard4th12phi - 
          2*gt21L*gtu13*PDstandard4th13phi - gt21L*gtu22*PDstandard4th22phi - 2*gt21L*gtu23*PDstandard4th23phi - 
          gt21L*gtu33*PDstandard4th33phi + Gt321*PDstandard4th3phi + gt21L*Gt311*gtu11*PDstandard4th3phi + 
          2*gt21L*Gt321*gtu12*PDstandard4th3phi + 2*gt21L*Gt331*gtu13*PDstandard4th3phi + 
          gt21L*Gt322*gtu22*PDstandard4th3phi + 2*gt21L*Gt332*gtu23*PDstandard4th3phi + 
          gt21L*Gt333*gtu33*PDstandard4th3phi + PDstandard4th1phi*
           (Gt121 + Gt111*gt21L*gtu11 + 2*Gt121*gt21L*gtu12 + 2*Gt131*gt21L*gtu13 + Gt122*gt21L*gtu22 + 
             2*Gt132*gt21L*gtu23 + Gt133*gt21L*gtu33 + (2 - 4*gt21L*gtu12)*PDstandard4th2phi - 
             4*gt21L*gtu13*PDstandard4th3phi) + PDstandard4th2phi*
           (Gt221 + 2*gt21L*Gt221*gtu12 + gt21L*(Gt211*gtu11 + 2*Gt231*gtu13 + Gt222*gtu22 + 2*Gt232*gtu23 + Gt233*gtu33) - 
             4*gt21L*gtu23*PDstandard4th3phi) - 2*gt21L*gtu11*SQR(PDstandard4th1phi) - 
          2*gt21L*gtu22*SQR(PDstandard4th2phi) - 2*gt21L*gtu33*SQR(PDstandard4th3phi));
    
    Rphi31  =  2*(-PDstandard4th13phi + gt31L*(-(gtu11*PDstandard4th11phi) - 2*gtu12*PDstandard4th12phi - 
             2*gtu13*PDstandard4th13phi) - gt31L*gtu22*PDstandard4th22phi - 2*gt31L*gtu23*PDstandard4th23phi - 
          gt31L*gtu33*PDstandard4th33phi + Gt331*PDstandard4th3phi + Gt311*gt31L*gtu11*PDstandard4th3phi + 
          2*gt31L*Gt321*gtu12*PDstandard4th3phi + 2*gt31L*Gt331*gtu13*PDstandard4th3phi + 
          gt31L*Gt322*gtu22*PDstandard4th3phi + 2*gt31L*Gt332*gtu23*PDstandard4th3phi + 
          gt31L*Gt333*gtu33*PDstandard4th3phi + PDstandard4th1phi*
           (Gt131 + Gt111*gt31L*gtu11 + 2*Gt121*gt31L*gtu12 + 2*Gt131*gt31L*gtu13 + Gt122*gt31L*gtu22 + 
             2*Gt132*gt31L*gtu23 + Gt133*gt31L*gtu33 - 4*gt31L*gtu12*PDstandard4th2phi + 
             (2 - 4*gt31L*gtu13)*PDstandard4th3phi) + 
          PDstandard4th2phi*(Gt231 + 2*Gt231*gt31L*gtu13 + 
             gt31L*(Gt211*gtu11 + 2*Gt221*gtu12 + Gt222*gtu22 + 2*Gt232*gtu23 + Gt233*gtu33) - 
             4*gt31L*gtu23*PDstandard4th3phi) - 2*gt31L*gtu11*SQR(PDstandard4th1phi) - 
          2*gt31L*gtu22*SQR(PDstandard4th2phi) - 2*gt31L*gtu33*SQR(PDstandard4th3phi));
    
    Rphi22  =  2*(-PDstandard4th22phi + gt22L*(-(gtu11*PDstandard4th11phi) - 2*gtu12*PDstandard4th12phi - 
             2*gtu13*PDstandard4th13phi - gtu22*PDstandard4th22phi) - 2*gt22L*gtu23*PDstandard4th23phi - 
          gt22L*gtu33*PDstandard4th33phi + Gt322*PDstandard4th3phi + gt22L*Gt311*gtu11*PDstandard4th3phi + 
          2*gt22L*Gt321*gtu12*PDstandard4th3phi + 2*gt22L*Gt331*gtu13*PDstandard4th3phi + 
          gt22L*Gt322*gtu22*PDstandard4th3phi + 2*gt22L*Gt332*gtu23*PDstandard4th3phi + 
          gt22L*Gt333*gtu33*PDstandard4th3phi + PDstandard4th1phi*
           (Gt122 + Gt111*gt22L*gtu11 + 2*Gt121*gt22L*gtu12 + 2*Gt131*gt22L*gtu13 + Gt122*gt22L*gtu22 + 
             2*Gt132*gt22L*gtu23 + Gt133*gt22L*gtu33 - 4*gt22L*gtu12*PDstandard4th2phi - 4*gt22L*gtu13*PDstandard4th3phi) + 
          PDstandard4th2phi*(Gt222 + Gt222*gt22L*gtu22 + 
             gt22L*(Gt211*gtu11 + 2*Gt221*gtu12 + 2*Gt231*gtu13 + 2*Gt232*gtu23 + Gt233*gtu33) - 
             4*gt22L*gtu23*PDstandard4th3phi) - 2*gt22L*gtu11*SQR(PDstandard4th1phi) + 
          (2 - 2*gt22L*gtu22)*SQR(PDstandard4th2phi) - 2*gt22L*gtu33*SQR(PDstandard4th3phi));
    
    Rphi32  =  2*(-PDstandard4th23phi + gt32L*(-(gtu11*PDstandard4th11phi) - 2*gtu12*PDstandard4th12phi - 
             2*gtu13*PDstandard4th13phi - gtu22*PDstandard4th22phi - 2*gtu23*PDstandard4th23phi) - 
          gt32L*gtu33*PDstandard4th33phi + Gt332*PDstandard4th3phi + Gt311*gt32L*gtu11*PDstandard4th3phi + 
          2*Gt321*gt32L*gtu12*PDstandard4th3phi + 2*gt32L*Gt331*gtu13*PDstandard4th3phi + 
          Gt322*gt32L*gtu22*PDstandard4th3phi + 2*gt32L*Gt332*gtu23*PDstandard4th3phi + 
          gt32L*Gt333*gtu33*PDstandard4th3phi + PDstandard4th1phi*
           (Gt132 + Gt111*gt32L*gtu11 + 2*Gt121*gt32L*gtu12 + 2*Gt131*gt32L*gtu13 + Gt122*gt32L*gtu22 + 
             2*Gt132*gt32L*gtu23 + Gt133*gt32L*gtu33 - 4*gt32L*gtu12*PDstandard4th2phi - 4*gt32L*gtu13*PDstandard4th3phi) + 
          PDstandard4th2phi*(Gt232 + 2*Gt232*gt32L*gtu23 + 
             gt32L*(Gt211*gtu11 + 2*Gt221*gtu12 + 2*Gt231*gtu13 + Gt222*gtu22 + Gt233*gtu33) + 
             (2 - 4*gt32L*gtu23)*PDstandard4th3phi) - 2*gt32L*gtu11*SQR(PDstandard4th1phi) - 
          2*gt32L*gtu22*SQR(PDstandard4th2phi) - 2*gt32L*gtu33*SQR(PDstandard4th3phi));
    
    Rphi33  =  2*(-PDstandard4th33phi + (Gt333 + gt33L*
              (Gt322*gtu22 + 2*(Gt321*gtu12 + Gt331*gtu13 + Gt332*gtu23) + Gt333*gtu33))*PDstandard4th3phi + 
          PDstandard4th2phi*(Gt233 + gt33L*(Gt211*gtu11 + Gt222*gtu22 + 2*(Gt221*gtu12 + Gt231*gtu13 + Gt232*gtu23) + 
                Gt233*gtu33 - 4*gtu23*PDstandard4th3phi)) + 
          PDstandard4th1phi*(Gt133 + gt33L*(Gt111*gtu11 + Gt122*gtu22 + 2*(Gt121*gtu12 + Gt131*gtu13 + Gt132*gtu23) + 
                Gt133*gtu33 - 4*(gtu12*PDstandard4th2phi + gtu13*PDstandard4th3phi))) + 2*SQR(PDstandard4th3phi) + 
          gt33L*(-(gtu11*PDstandard4th11phi) - 2*gtu12*PDstandard4th12phi - 2*gtu13*PDstandard4th13phi - 
             gtu22*PDstandard4th22phi - 2*gtu23*PDstandard4th23phi - gtu33*PDstandard4th33phi + 
             Gt311*gtu11*PDstandard4th3phi - 2*gtu11*SQR(PDstandard4th1phi) - 2*gtu22*SQR(PDstandard4th2phi) - 
             2*gtu33*SQR(PDstandard4th3phi)));
    
    R11  =  Rphi11 + Rt11;
    
    R21  =  Rphi21 + Rt21;
    
    R31  =  Rphi31 + Rt31;
    
    R22  =  Rphi22 + Rt22;
    
    R32  =  Rphi32 + Rt32;
    
    R33  =  Rphi33 + Rt33;
    
    Atm11  =  At11L*gtu11 + At21L*gtu12 + At31L*gtu13;
    
    Atm21  =  At11L*gtu12 + At21L*gtu22 + At31L*gtu23;
    
    Atm31  =  At11L*gtu13 + At21L*gtu23 + At31L*gtu33;
    
    Atm12  =  At21L*gtu11 + At22L*gtu12 + At32L*gtu13;
    
    Atm22  =  At21L*gtu12 + At22L*gtu22 + At32L*gtu23;
    
    Atm32  =  At21L*gtu13 + At22L*gtu23 + At32L*gtu33;
    
    Atm13  =  At31L*gtu11 + At32L*gtu12 + At33L*gtu13;
    
    Atm23  =  At31L*gtu12 + At32L*gtu22 + At33L*gtu23;
    
    Atm33  =  At31L*gtu13 + At32L*gtu23 + At33L*gtu33;
    
    Atu11  =  Atm11*gtu11 + Atm12*gtu12 + Atm13*gtu13;
    
    Atu12  =  Atm11*gtu12 + Atm12*gtu22 + Atm13*gtu23;
    
    Atu13  =  Atm11*gtu13 + Atm12*gtu23 + Atm13*gtu33;
    
    Atu22  =  Atm21*gtu12 + Atm22*gtu22 + Atm23*gtu23;
    
    Atu23  =  Atm21*gtu13 + Atm22*gtu23 + Atm23*gtu33;
    
    Atu33  =  Atm31*gtu13 + Atm32*gtu23 + Atm33*gtu33;
    
    e4phi  =  exp(4*phiL);
    
    em4phi  =  INV(e4phi);
    
    g11  =  e4phi*gt11L;
    
    g21  =  e4phi*gt21L;
    
    g31  =  e4phi*gt31L;
    
    g22  =  e4phi*gt22L;
    
    g32  =  e4phi*gt32L;
    
    g33  =  e4phi*gt33L;
    
    detg  =  2*g21*g31*g32 + g33*(g11*g22 - SQR(g21)) - g22*SQR(g31) - g11*SQR(g32);
    
    gu11  =  em4phi*gtu11;
    
    gu12  =  em4phi*gtu12;
    
    gu13  =  em4phi*gtu13;
    
    gu22  =  em4phi*gtu22;
    
    gu23  =  em4phi*gtu23;
    
    gu33  =  em4phi*gtu33;
    
    ddetg1  =  e4phi*(ddetgt1 + 4*detgt*PDstandard4th1phi);
    
    ddetg2  =  e4phi*(ddetgt2 + 4*detgt*PDstandard4th2phi);
    
    ddetg3  =  e4phi*(ddetgt3 + 4*detgt*PDstandard4th3phi);
    
    G111  =  -((-6*detg*Gt111 + ddetg1*(-6 + g11*gu11) + g11*(ddetg2*gu12 + ddetg3*gu13))*INV(detg))/6.;
    
    G211  =  ((6*detg*Gt211 - g11*(ddetg1*gu12 + ddetg2*gu22 + ddetg3*gu23))*INV(detg))/6.;
    
    G311  =  ((6*detg*Gt311 - g11*(ddetg1*gu13 + ddetg2*gu23 + ddetg3*gu33))*INV(detg))/6.;
    
    G121  =  -((-6*detg*Gt121 + ddetg2*(-3 + g21*gu12) + g21*(ddetg1*gu11 + ddetg3*gu13))*INV(detg))/6.;
    
    G221  =  -((-6*detg*Gt221 + ddetg1*(-3 + g21*gu12) + g21*(ddetg2*gu22 + ddetg3*gu23))*INV(detg))/6.;
    
    G321  =  ((6*detg*Gt321 - g21*(ddetg1*gu13 + ddetg2*gu23 + ddetg3*gu33))*INV(detg))/6.;
    
    G131  =  -((-6*detg*Gt131 + g31*(ddetg1*gu11 + ddetg2*gu12) + ddetg3*(-3 + g31*gu13))*INV(detg))/6.;
    
    G231  =  ((6*detg*Gt231 - g31*(ddetg1*gu12 + ddetg2*gu22 + ddetg3*gu23))*INV(detg))/6.;
    
    G331  =  -((-6*detg*Gt331 + ddetg1*(-3 + g31*gu13) + g31*(ddetg2*gu23 + ddetg3*gu33))*INV(detg))/6.;
    
    G122  =  ((6*detg*Gt122 - g22*(ddetg1*gu11 + ddetg2*gu12 + ddetg3*gu13))*INV(detg))/6.;
    
    G222  =  -((-6*detg*Gt222 + ddetg2*(-6 + g22*gu22) + g22*(ddetg1*gu12 + ddetg3*gu23))*INV(detg))/6.;
    
    G322  =  ((6*detg*Gt322 - g22*(ddetg1*gu13 + ddetg2*gu23 + ddetg3*gu33))*INV(detg))/6.;
    
    G132  =  ((6*detg*Gt132 - g32*(ddetg1*gu11 + ddetg2*gu12 + ddetg3*gu13))*INV(detg))/6.;
    
    G232  =  -((-6*detg*Gt232 + g32*(ddetg1*gu12 + ddetg2*gu22) + ddetg3*(-3 + g32*gu23))*INV(detg))/6.;
    
    G332  =  -((-6*detg*Gt332 + ddetg2*(-3 + g32*gu23) + g32*(ddetg1*gu13 + ddetg3*gu33))*INV(detg))/6.;
    
    G133  =  ((6*detg*Gt133 - g33*(ddetg1*gu11 + ddetg2*gu12 + ddetg3*gu13))*INV(detg))/6.;
    
    G233  =  ((6*detg*Gt233 - g33*(ddetg1*gu12 + ddetg2*gu22 + ddetg3*gu23))*INV(detg))/6.;
    
    G333  =  -((-6*detg*Gt333 + g33*(ddetg1*gu13 + ddetg2*gu23) + ddetg3*(-6 + g33*gu33))*INV(detg))/6.;
    
    phirhsL  =  beta1L*PDstandard4th1phi + beta2L*PDstandard4th2phi + beta3L*PDstandard4th3phi + 
        ((PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3)*phiL)/6. - (alphaL*trKL)/6.;
    
    gt11rhsL  =  -2*alphaL*At11L + 2*(gt11L*PDstandard4th1beta1 + gt21L*PDstandard4th1beta2 + gt31L*PDstandard4th1beta3) + 
        beta1L*PDstandard4th1gt11 + beta2L*PDstandard4th2gt11 - 
        gt11L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) + beta3L*PDstandard4th3gt11;
    
    gt21rhsL  =  -2*alphaL*At21L + gt22L*PDstandard4th1beta2 + gt32L*PDstandard4th1beta3 + beta1L*PDstandard4th1gt21 + 
        gt11L*PDstandard4th2beta1 + gt31L*PDstandard4th2beta3 + beta2L*PDstandard4th2gt21 + 
        gt21L*(kthird*(PDstandard4th1beta1 + PDstandard4th2beta2) - ktwothird*PDstandard4th3beta3) + 
        beta3L*PDstandard4th3gt21;
    
    gt31rhsL  =  -2*alphaL*At31L + gt32L*PDstandard4th1beta2 + gt33L*PDstandard4th1beta3 + beta1L*PDstandard4th1gt31 + 
        beta2L*PDstandard4th2gt31 + gt11L*PDstandard4th3beta1 + gt21L*PDstandard4th3beta2 + 
        gt31L*(-(ktwothird*PDstandard4th2beta2) + kthird*(PDstandard4th1beta1 + PDstandard4th3beta3)) + 
        beta3L*PDstandard4th3gt31;
    
    gt22rhsL  =  -2*alphaL*At22L + beta1L*PDstandard4th1gt22 + 
        2*(gt21L*PDstandard4th2beta1 + gt22L*PDstandard4th2beta2 + gt32L*PDstandard4th2beta3) + beta2L*PDstandard4th2gt22 - 
        gt22L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) + beta3L*PDstandard4th3gt22;
    
    gt32rhsL  =  -2*alphaL*At32L + beta1L*PDstandard4th1gt32 + gt31L*PDstandard4th2beta1 + gt33L*PDstandard4th2beta3 + 
        beta2L*PDstandard4th2gt32 + gt21L*PDstandard4th3beta1 + gt22L*PDstandard4th3beta2 + 
        gt32L*(-(ktwothird*PDstandard4th1beta1) + kthird*(PDstandard4th2beta2 + PDstandard4th3beta3)) + 
        beta3L*PDstandard4th3gt32;
    
    gt33rhsL  =  -2*alphaL*At33L + beta1L*PDstandard4th1gt33 + beta2L*PDstandard4th2gt33 - 
        gt33L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) + 
        2*(gt31L*PDstandard4th3beta1 + gt32L*PDstandard4th3beta2 + gt33L*PDstandard4th3beta3) + beta3L*PDstandard4th3gt33;
    
    Xt1rhsL  =  kthird*(-3*(beta1L*(ddgtu1111 + ddgtu1221 + ddgtu1331) + beta2L*(ddgtu1121 + ddgtu1222 + ddgtu1332) + 
             beta3L*(ddgtu1131 + ddgtu1232 + ddgtu1333)) + 16*(gtu12*PDstandard4th12beta1 + gtu13*PDstandard4th13beta1) + 
          (7*dgtu111 + 4*(dgtu122 + dgtu133))*PDstandard4th1beta1 + 
          gtu11*(10*PDstandard4th11beta1 + 4*(PDstandard4th12beta2 + PDstandard4th13beta3) - 4*alphaL*PDstandard4th1trK) + 
          12*gtu23*PDstandard4th23beta1 + (9*dgtu121 + 6*dgtu233)*PDstandard4th2beta1 + 
          (dgtu122 - 2*(dgtu111 + dgtu133))*PDstandard4th2beta2 + 
          4*(gtu12*(PDstandard4th22beta2 + PDstandard4th23beta3) + gtu13*(PDstandard4th23beta2 + PDstandard4th33beta3)) - 
          6*(Atu11*PDstandard4th1alpha + Atu12*PDstandard4th2alpha + Atu13*PDstandard4th3alpha) + 
          (9*dgtu131 + 6*dgtu333)*PDstandard4th3beta1 + 
          6*(gtu22*PDstandard4th22beta1 + dgtu222*PDstandard4th2beta1 + gtu33*PDstandard4th33beta1 + 
             dgtu232*PDstandard4th3beta1) + 3*(dgtu112*PDstandard4th1beta2 + dgtu113*PDstandard4th1beta3 + 
             dgtu123*PDstandard4th2beta3 + dgtu132*PDstandard4th3beta2) + 
          (-2*(dgtu111 + dgtu122) + dgtu133)*PDstandard4th3beta3 + 
          alphaL*(12*(Atu12*Gt121 + Atu13*Gt131 + Atu23*Gt132) + 6*(Atu11*Gt111 + Atu22*Gt122 + Atu33*Gt133) + 
             36*(Atu11*PDstandard4th1phi + Atu12*PDstandard4th2phi + Atu13*PDstandard4th3phi) - 
             4*(gtu12*PDstandard4th2trK + gtu13*PDstandard4th3trK)));
    
    Xt2rhsL  =  kthird*(-3*(beta1L*(ddgtu1211 + ddgtu2221 + ddgtu2331) + beta2L*(ddgtu1221 + ddgtu2222 + ddgtu2332) + 
             beta3L*(ddgtu1231 + ddgtu2232 + ddgtu2333)) + 12*gtu13*PDstandard4th13beta2 + 
          (dgtu121 - 2*(dgtu222 + dgtu233))*PDstandard4th1beta1 + (9*dgtu122 + 6*dgtu133)*PDstandard4th1beta2 + 
          gtu12*(16*PDstandard4th12beta2 + 4*PDstandard4th13beta3 - 4*alphaL*PDstandard4th1trK) + 
          (7*dgtu222 + 4*dgtu233)*PDstandard4th2beta2 + 
          4*(gtu12*PDstandard4th11beta1 + gtu23*PDstandard4th13beta1 + 
             gtu22*(PDstandard4th12beta1 + PDstandard4th23beta3) + dgtu121*PDstandard4th2beta2) + 
          gtu22*(10*PDstandard4th22beta2 - 4*alphaL*PDstandard4th2trK) - 
          6*(Atu12*PDstandard4th1alpha + Atu22*PDstandard4th2alpha + Atu23*PDstandard4th3alpha) + 
          3*(dgtu123*PDstandard4th1beta3 + dgtu221*PDstandard4th2beta1 + dgtu223*PDstandard4th2beta3 + 
             dgtu231*PDstandard4th3beta1) + (9*dgtu232 + 6*dgtu333)*PDstandard4th3beta2 + 
          6*(gtu11*PDstandard4th11beta2 + dgtu111*PDstandard4th1beta2 + gtu33*PDstandard4th33beta2 + 
             dgtu131*PDstandard4th3beta2) + (-2*(dgtu121 + dgtu222) + dgtu233)*PDstandard4th3beta3 + 
          alphaL*(12*(Atu12*Gt221 + Atu13*Gt231 + Atu23*Gt232) + 6*(Atu11*Gt211 + Atu22*Gt222 + Atu33*Gt233) + 
             36*(Atu12*PDstandard4th1phi + Atu22*PDstandard4th2phi + Atu23*PDstandard4th3phi)) + 
          gtu23*(16*PDstandard4th23beta2 + 4*PDstandard4th33beta3 - 4*alphaL*PDstandard4th3trK));
    
    Xt3rhsL  =  kthird*(-3*(beta1L*(ddgtu1311 + ddgtu2321 + ddgtu3331) + beta2L*(ddgtu1321 + ddgtu2322 + ddgtu3332) + 
             beta3L*(ddgtu1331 + ddgtu2332 + ddgtu3333)) + 12*gtu12*PDstandard4th12beta3 + 
          (dgtu131 - 2*(dgtu232 + dgtu333))*PDstandard4th1beta1 + (6*dgtu122 + 9*dgtu133)*PDstandard4th1beta3 + 
          gtu13*(16*PDstandard4th13beta3 - 4*alphaL*PDstandard4th1trK) + 
          (dgtu232 - 2*(dgtu131 + dgtu333))*PDstandard4th2beta2 + (6*dgtu222 + 9*dgtu233)*PDstandard4th2beta3 + 
          6*(gtu11*PDstandard4th11beta3 + dgtu111*PDstandard4th1beta3 + gtu22*PDstandard4th22beta3 + 
             dgtu121*PDstandard4th2beta3) + gtu23*(16*PDstandard4th23beta3 - 4*alphaL*PDstandard4th2trK) - 
          6*(Atu13*PDstandard4th1alpha + Atu23*PDstandard4th2alpha + Atu33*PDstandard4th3alpha) + 
          3*(dgtu132*PDstandard4th1beta2 + dgtu231*PDstandard4th2beta1 + dgtu331*PDstandard4th3beta1 + 
             dgtu332*PDstandard4th3beta2) + (4*dgtu232 + 7*dgtu333)*PDstandard4th3beta3 + 
          4*(gtu13*(PDstandard4th11beta1 + PDstandard4th12beta2) + gtu23*(PDstandard4th12beta1 + PDstandard4th22beta2) + 
             gtu33*(PDstandard4th13beta1 + PDstandard4th23beta2) + dgtu131*PDstandard4th3beta3) + 
          alphaL*(12*(Atu12*Gt321 + Atu13*Gt331 + Atu23*Gt332) + 6*(Atu11*Gt311 + Atu22*Gt322 + Atu33*Gt333) + 
             36*(Atu13*PDstandard4th1phi + Atu23*PDstandard4th2phi + Atu33*PDstandard4th3phi)) + 
          gtu33*(10*PDstandard4th33beta3 - 4*alphaL*PDstandard4th3trK));
    
    trKrhsL  =  -(gu11*PDstandard4th11alpha) - 2*gu12*PDstandard4th12alpha - 2*gu13*PDstandard4th13alpha + 
        (G111*gu11 + 2*G121*gu12 + 2*G131*gu13 + G122*gu22 + 2*G132*gu23 + G133*gu33)*PDstandard4th1alpha + 
        beta1L*PDstandard4th1trK - gu22*PDstandard4th22alpha - 2*gu23*PDstandard4th23alpha + 
        (G211*gu11 + 2*G221*gu12 + 2*G231*gu13 + G222*gu22 + 2*G232*gu23 + G233*gu33)*PDstandard4th2alpha + 
        beta2L*PDstandard4th2trK - gu33*PDstandard4th33alpha + G311*gu11*PDstandard4th3alpha + 
        2*G331*gu13*PDstandard4th3alpha + G322*gu22*PDstandard4th3alpha + 2*G332*gu23*PDstandard4th3alpha + 
        G333*gu33*PDstandard4th3alpha + 2*(alphaL*Atm12*Atm21 + alphaL*Atm13*Atm31 + alphaL*Atm23*Atm32 + 
           G321*gu12*PDstandard4th3alpha) + beta3L*PDstandard4th3trK + alphaL*SQR(Atm11) + alphaL*SQR(Atm22) + 
        alphaL*SQR(Atm33) + alphaL*kthird*SQR(trKL);
    
    At11rhsL  =  beta1L*PDstandard4th1At11 + 2*(At11L*PDstandard4th1beta1 + At21L*PDstandard4th1beta2 + 
           At31L*PDstandard4th1beta3) + beta2L*PDstandard4th2At11 + beta3L*PDstandard4th3At11 - 
        At11L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) - 
        em4phi*kthird*(3*PDstandard4th11alpha - g11*gu11*PDstandard4th11alpha - 2*g11*gu12*PDstandard4th12alpha - 
           2*g11*gu13*PDstandard4th13alpha + (G111*(-3 + g11*gu11) + 
              g11*(2*G121*gu12 + 2*G131*gu13 + G122*gu22 + 2*G132*gu23 + G133*gu33))*PDstandard4th1alpha - 
           g11*gu22*PDstandard4th22alpha - 2*g11*gu23*PDstandard4th23alpha + 
           (G211*(-3 + g11*gu11) + g11*(2*G221*gu12 + 2*G231*gu13 + G222*gu22 + 2*G232*gu23 + G233*gu33))*
            PDstandard4th2alpha - g11*gu33*PDstandard4th33alpha - 3*G311*PDstandard4th3alpha + 
           g11*G311*gu11*PDstandard4th3alpha + 2*g11*G321*gu12*PDstandard4th3alpha + 2*g11*G331*gu13*PDstandard4th3alpha + 
           g11*G322*gu22*PDstandard4th3alpha + 2*g11*G332*gu23*PDstandard4th3alpha + g11*G333*gu33*PDstandard4th3alpha - 
           3*alphaL*R11 + alphaL*g11*gu11*R11 + 2*alphaL*g11*gu12*R21 + alphaL*g11*gu22*R22 + 2*alphaL*g11*gu13*R31 + 
           2*alphaL*g11*gu23*R32 + alphaL*g11*gu33*R33) + alphaL*(-2*(At11L*Atm11 + At21L*Atm21 + At31L*Atm31) + At11L*trKL);
    
    At21rhsL  =  beta1L*PDstandard4th1At21 + At22L*PDstandard4th1beta2 + At32L*PDstandard4th1beta3 + 
        beta2L*PDstandard4th2At21 + At11L*PDstandard4th2beta1 + At31L*PDstandard4th2beta3 + beta3L*PDstandard4th3At21 + 
        At21L*(PDstandard4th1beta1 + PDstandard4th2beta2 - 
           ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3)) - 
        em4phi*kthird*(3*PDstandard4th12alpha + (G121*(-3 + 2*g21*gu12) + 
              g21*(G111*gu11 + G122*gu22 + 2*(G131*gu13 + G132*gu23) + G133*gu33))*PDstandard4th1alpha + 
           (G221*(-3 + 2*g21*gu12) + g21*(G211*gu11 + G222*gu22 + 2*(G231*gu13 + G232*gu23) + G233*gu33))*
            PDstandard4th2alpha + G321*(-3 + 2*g21*gu12)*PDstandard4th3alpha + 
           alphaL*(-3*R21 + g21*(gu22*R22 + 2*gu23*R32)) + 
           g21*(-(gu11*PDstandard4th11alpha) - 2*gu12*PDstandard4th12alpha - 2*gu13*PDstandard4th13alpha - 
              gu22*PDstandard4th22alpha - 2*gu23*PDstandard4th23alpha - gu33*PDstandard4th33alpha + 
              G311*gu11*PDstandard4th3alpha + 2*G331*gu13*PDstandard4th3alpha + G322*gu22*PDstandard4th3alpha + 
              2*G332*gu23*PDstandard4th3alpha + G333*gu33*PDstandard4th3alpha + alphaL*gu11*R11 + 2*alphaL*gu12*R21 + 
              2*alphaL*gu13*R31 + alphaL*gu33*R33)) + alphaL*(-2*(At11L*Atm12 + At21L*Atm22 + At31L*Atm32) + At21L*trKL);
    
    At31rhsL  =  beta1L*PDstandard4th1At31 + At32L*PDstandard4th1beta2 + At33L*PDstandard4th1beta3 + 
        beta2L*PDstandard4th2At31 + beta3L*PDstandard4th3At31 + At11L*PDstandard4th3beta1 + At21L*PDstandard4th3beta2 + 
        At31L*(PDstandard4th1beta1 + PDstandard4th3beta3 - 
           ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3)) - 
        em4phi*kthird*(3*PDstandard4th13alpha + (G131*(-3 + 2*g31*gu13) + 
              g31*(G111*gu11 + G122*gu22 + 2*(G121*gu12 + G132*gu23) + G133*gu33))*PDstandard4th1alpha + 
           (G231*(-3 + 2*g31*gu13) + g31*(G211*gu11 + G222*gu22 + 2*(G221*gu12 + G232*gu23) + G233*gu33))*
            PDstandard4th2alpha + (G331*(-3 + 2*g31*gu13) + g31*(G322*gu22 + 2*(G321*gu12 + G332*gu23) + G333*gu33))*
            PDstandard4th3alpha + alphaL*(-3*R31 + 2*g31*gu23*R32) + 
           g31*(-(gu11*PDstandard4th11alpha) - 2*gu12*PDstandard4th12alpha - 2*gu13*PDstandard4th13alpha - 
              gu22*PDstandard4th22alpha - 2*gu23*PDstandard4th23alpha - gu33*PDstandard4th33alpha + 
              G311*gu11*PDstandard4th3alpha + alphaL*gu11*R11 + 2*alphaL*gu12*R21 + alphaL*gu22*R22 + 2*alphaL*gu13*R31 + 
              alphaL*gu33*R33)) + alphaL*(-2*(At11L*Atm13 + At21L*Atm23 + At31L*Atm33) + At31L*trKL);
    
    At22rhsL  =  beta1L*PDstandard4th1At22 + beta2L*PDstandard4th2At22 + 
        2*(At21L*PDstandard4th2beta1 + At22L*PDstandard4th2beta2 + At32L*PDstandard4th2beta3) + beta3L*PDstandard4th3At22 - 
        At22L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) - 
        em4phi*kthird*(-(g22*gu11*PDstandard4th11alpha) - 2*g22*gu12*PDstandard4th12alpha - 
           2*g22*gu13*PDstandard4th13alpha + (G122*(-3 + g22*gu22) + 
              g22*(G111*gu11 + 2*G121*gu12 + 2*G131*gu13 + 2*G132*gu23 + G133*gu33))*PDstandard4th1alpha + 
           3*PDstandard4th22alpha - g22*gu22*PDstandard4th22alpha - 2*g22*gu23*PDstandard4th23alpha + 
           (G222*(-3 + g22*gu22) + g22*(G211*gu11 + 2*G221*gu12 + 2*G231*gu13 + 2*G232*gu23 + G233*gu33))*
            PDstandard4th2alpha - g22*gu33*PDstandard4th33alpha - 3*G322*PDstandard4th3alpha + 
           g22*G311*gu11*PDstandard4th3alpha + 2*g22*G321*gu12*PDstandard4th3alpha + 2*g22*G331*gu13*PDstandard4th3alpha + 
           g22*G322*gu22*PDstandard4th3alpha + 2*g22*G332*gu23*PDstandard4th3alpha + g22*G333*gu33*PDstandard4th3alpha + 
           alphaL*g22*gu11*R11 + 2*alphaL*g22*gu12*R21 - 3*alphaL*R22 + alphaL*g22*gu22*R22 + 2*alphaL*g22*gu13*R31 + 
           2*alphaL*g22*gu23*R32 + alphaL*g22*gu33*R33) + alphaL*(-2*(At21L*Atm12 + At22L*Atm22 + At32L*Atm32) + At22L*trKL);
    
    At32rhsL  =  beta1L*PDstandard4th1At32 + beta2L*PDstandard4th2At32 + At31L*PDstandard4th2beta1 + 
        At33L*PDstandard4th2beta3 + beta3L*PDstandard4th3At32 + At21L*PDstandard4th3beta1 + At22L*PDstandard4th3beta2 + 
        At32L*(PDstandard4th2beta2 + PDstandard4th3beta3 - 
           ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3)) - 
        em4phi*kthird*((G132*(-3 + 2*g32*gu23) + g32*(G111*gu11 + 2*(G121*gu12 + G131*gu13) + G122*gu22 + G133*gu33))*
            PDstandard4th1alpha + 3*PDstandard4th23alpha + 
           (G232*(-3 + 2*g32*gu23) + g32*(G211*gu11 + 2*(G221*gu12 + G231*gu13) + G222*gu22 + G233*gu33))*
            PDstandard4th2alpha + (G332*(-3 + 2*g32*gu23) + g32*(2*(G321*gu12 + G331*gu13) + G322*gu22 + G333*gu33))*
            PDstandard4th3alpha + g32*(-(gu11*PDstandard4th11alpha) - 2*gu12*PDstandard4th12alpha - 
              2*gu13*PDstandard4th13alpha - gu22*PDstandard4th22alpha - 2*gu23*PDstandard4th23alpha - 
              gu33*PDstandard4th33alpha + G311*gu11*PDstandard4th3alpha + alphaL*gu11*R11 + 2*alphaL*gu12*R21 + 
              alphaL*gu22*R22 + 2*alphaL*gu13*R31 + 2*alphaL*gu23*R32) + alphaL*(-3*R32 + g32*gu33*R33)) + 
        alphaL*(-2*(At21L*Atm13 + At22L*Atm23 + At32L*Atm33) + At32L*trKL);
    
    At33rhsL  =  beta1L*PDstandard4th1At33 + beta2L*PDstandard4th2At33 + beta3L*PDstandard4th3At33 - 
        At33L*ktwothird*(PDstandard4th1beta1 + PDstandard4th2beta2 + PDstandard4th3beta3) + 
        2*(At31L*PDstandard4th3beta1 + At32L*PDstandard4th3beta2 + At33L*PDstandard4th3beta3) - 
        em4phi*kthird*(-(g33*gu11*PDstandard4th11alpha) - 2*g33*gu12*PDstandard4th12alpha - 
           2*g33*gu13*PDstandard4th13alpha + (g33*(G111*gu11 + 2*G121*gu12 + 2*G131*gu13 + G122*gu22 + 2*G132*gu23) + 
              G133*(-3 + g33*gu33))*PDstandard4th1alpha - g33*gu22*PDstandard4th22alpha - 2*g33*gu23*PDstandard4th23alpha + 
           (g33*(G211*gu11 + 2*G221*gu12 + 2*G231*gu13 + G222*gu22 + 2*G232*gu23) + G233*(-3 + g33*gu33))*
            PDstandard4th2alpha + 3*PDstandard4th33alpha - g33*gu33*PDstandard4th33alpha - 3*G333*PDstandard4th3alpha + 
           G311*g33*gu11*PDstandard4th3alpha + 2*G321*g33*gu12*PDstandard4th3alpha + 2*g33*G331*gu13*PDstandard4th3alpha + 
           G322*g33*gu22*PDstandard4th3alpha + 2*g33*G332*gu23*PDstandard4th3alpha + g33*G333*gu33*PDstandard4th3alpha + 
           alphaL*g33*gu11*R11 + 2*alphaL*g33*gu12*R21 + alphaL*g33*gu22*R22 + 2*alphaL*g33*gu13*R31 + 
           2*alphaL*g33*gu23*R32 - 3*alphaL*R33 + alphaL*g33*gu33*R33) + 
        alphaL*(-2*(At31L*Atm13 + At32L*Atm23 + At33L*Atm33) + At33L*trKL);
    
    alpharhsL  =  beta1L*PDstandard4th1alpha + beta2L*PDstandard4th2alpha + beta3L*PDstandard4th3alpha - 
        dtalphaL*harmonicF*pow(alphaL,harmonicN);
    
    dtalpharhsL  =  -(AlphaDriver*dtalphaL) + trKrhsL;
    
    beta1rhsL  =  dtbeta1L*ShiftGammaCoeff*pow(alphaL,ShiftAlphaPower);
    
    beta2rhsL  =  dtbeta2L*ShiftGammaCoeff*pow(alphaL,ShiftAlphaPower);
    
    beta3rhsL  =  dtbeta3L*ShiftGammaCoeff*pow(alphaL,ShiftAlphaPower);
    
    dtbeta1rhsL  =  -(BetaDriver*dtbeta1L) + Xt1rhsL;
    
    dtbeta2rhsL  =  -(BetaDriver*dtbeta2L) + Xt2rhsL;
    
    dtbeta3rhsL  =  -(BetaDriver*dtbeta3L) + Xt3rhsL;
    
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    At11rhs[index] = At11rhsL;
    At21rhs[index] = At21rhsL;
    At22rhs[index] = At22rhsL;
    At31rhs[index] = At31rhsL;
    At32rhs[index] = At32rhsL;
    At33rhs[index] = At33rhsL;
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    dtalpharhs[index] = dtalpharhsL;
    dtbeta1rhs[index] = dtbeta1rhsL;
    dtbeta2rhs[index] = dtbeta2rhsL;
    dtbeta3rhs[index] = dtbeta3rhsL;
    gt11rhs[index] = gt11rhsL;
    gt21rhs[index] = gt21rhsL;
    gt22rhs[index] = gt22rhsL;
    gt31rhs[index] = gt31rhsL;
    gt32rhs[index] = gt32rhsL;
    gt33rhs[index] = gt33rhsL;
    phirhs[index] = phirhsL;
    trKrhs[index] = trKrhsL;
    Xt1rhs[index] = Xt1rhsL;
    Xt2rhs[index] = Xt2rhsL;
    Xt3rhs[index] = Xt3rhsL;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (somename);
}

void ML_BSSN_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_RHS_Body);
}
