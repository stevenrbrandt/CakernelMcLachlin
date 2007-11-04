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

void BSSN_constraints_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering BSSN_constraints_Body");
  }
  
  if (cctk_iteration % BSSN_constraints_calc_every != BSSN_constraints_calc_offset)
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
            max[0]-min[0],max[1]-min[1],max[2]-min[2])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL ddetg1 = INITVALUE, ddetg2 = INITVALUE, ddetg3 = INITVALUE;
    CCTK_REAL ddetgt1 = INITVALUE, ddetgt2 = INITVALUE, ddetgt3 = INITVALUE;
    CCTK_REAL detg = INITVALUE;
    CCTK_REAL detgt = INITVALUE;
    CCTK_REAL e4phi = INITVALUE;
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
    CCTK_REAL gK112 = INITVALUE, gK113 = INITVALUE, gK211 = INITVALUE, gK212 = INITVALUE, gK213 = INITVALUE, gK221 = INITVALUE;
    CCTK_REAL gK223 = INITVALUE, gK311 = INITVALUE, gK312 = INITVALUE, gK313 = INITVALUE, gK321 = INITVALUE, gK322 = INITVALUE;
    CCTK_REAL gK323 = INITVALUE, gK331 = INITVALUE, gK332 = INITVALUE;
    CCTK_REAL Gt111 = INITVALUE, Gt121 = INITVALUE, Gt122 = INITVALUE, Gt131 = INITVALUE, Gt132 = INITVALUE, Gt133 = INITVALUE;
    CCTK_REAL Gt211 = INITVALUE, Gt221 = INITVALUE, Gt222 = INITVALUE, Gt231 = INITVALUE, Gt232 = INITVALUE, Gt233 = INITVALUE;
    CCTK_REAL Gt311 = INITVALUE, Gt321 = INITVALUE, Gt322 = INITVALUE, Gt331 = INITVALUE, Gt332 = INITVALUE, Gt333 = INITVALUE;
    CCTK_REAL gtu11 = INITVALUE, gtu12 = INITVALUE, gtu13 = INITVALUE, gtu22 = INITVALUE, gtu23 = INITVALUE, gtu33 = INITVALUE;
    CCTK_REAL gu11 = INITVALUE, gu12 = INITVALUE, gu13 = INITVALUE, gu22 = INITVALUE, gu23 = INITVALUE, gu33 = INITVALUE;
    CCTK_REAL K11 = INITVALUE, K21 = INITVALUE, K22 = INITVALUE, K31 = INITVALUE, K32 = INITVALUE, K33 = INITVALUE;
    CCTK_REAL Km11 = INITVALUE, Km12 = INITVALUE, Km13 = INITVALUE, Km21 = INITVALUE, Km22 = INITVALUE, Km23 = INITVALUE;
    CCTK_REAL Km31 = INITVALUE, Km32 = INITVALUE, Km33 = INITVALUE;
    CCTK_REAL R11 = INITVALUE, R21 = INITVALUE, R22 = INITVALUE, R31 = INITVALUE, R32 = INITVALUE, R33 = INITVALUE;
    CCTK_REAL Rphi11 = INITVALUE, Rphi21 = INITVALUE, Rphi22 = INITVALUE, Rphi31 = INITVALUE, Rphi32 = INITVALUE, Rphi33 = INITVALUE;
    CCTK_REAL Rt11 = INITVALUE, Rt21 = INITVALUE, Rt22 = INITVALUE, Rt31 = INITVALUE, Rt32 = INITVALUE, Rt33 = INITVALUE;
    CCTK_REAL trR = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL At11L = INITVALUE, At21L = INITVALUE, At22L = INITVALUE, At31L = INITVALUE, At32L = INITVALUE, At33L = INITVALUE;
    CCTK_REAL cAL = INITVALUE;
    CCTK_REAL cSL = INITVALUE;
    CCTK_REAL cXt1L = INITVALUE, cXt2L = INITVALUE, cXt3L = INITVALUE;
    CCTK_REAL gt11L = INITVALUE, gt21L = INITVALUE, gt22L = INITVALUE, gt31L = INITVALUE, gt32L = INITVALUE, gt33L = INITVALUE;
    CCTK_REAL HL = INITVALUE;
    CCTK_REAL M1L = INITVALUE, M2L = INITVALUE, M3L = INITVALUE;
    CCTK_REAL phiL = INITVALUE;
    CCTK_REAL trKL = INITVALUE;
    CCTK_REAL Xt1L = INITVALUE, Xt2L = INITVALUE, Xt3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandard4th2At11 = INITVALUE;
    CCTK_REAL PDstandard4th3At11 = INITVALUE;
    CCTK_REAL PDstandard4th1At21 = INITVALUE;
    CCTK_REAL PDstandard4th2At21 = INITVALUE;
    CCTK_REAL PDstandard4th3At21 = INITVALUE;
    CCTK_REAL PDstandard4th1At22 = INITVALUE;
    CCTK_REAL PDstandard4th3At22 = INITVALUE;
    CCTK_REAL PDstandard4th1At31 = INITVALUE;
    CCTK_REAL PDstandard4th2At31 = INITVALUE;
    CCTK_REAL PDstandard4th3At31 = INITVALUE;
    CCTK_REAL PDstandard4th1At32 = INITVALUE;
    CCTK_REAL PDstandard4th2At32 = INITVALUE;
    CCTK_REAL PDstandard4th3At32 = INITVALUE;
    CCTK_REAL PDstandard4th1At33 = INITVALUE;
    CCTK_REAL PDstandard4th2At33 = INITVALUE;
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
    At11L = At11[index];
    At21L = At21[index];
    At22L = At22[index];
    At31L = At31[index];
    At32L = At32[index];
    At33L = At33[index];
    gt11L = gt11[index];
    gt21L = gt21[index];
    gt22L = gt22[index];
    gt31L = gt31[index];
    gt32L = gt32[index];
    gt33L = gt33[index];
    phiL = phi[index];
    trKL = trK[index];
    Xt1L = Xt1[index];
    Xt2L = Xt2[index];
    Xt3L = Xt3[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandard4th2At11 = PDstandard4th2(At11, i, j, k);
    PDstandard4th3At11 = PDstandard4th3(At11, i, j, k);
    PDstandard4th1At21 = PDstandard4th1(At21, i, j, k);
    PDstandard4th2At21 = PDstandard4th2(At21, i, j, k);
    PDstandard4th3At21 = PDstandard4th3(At21, i, j, k);
    PDstandard4th1At22 = PDstandard4th1(At22, i, j, k);
    PDstandard4th3At22 = PDstandard4th3(At22, i, j, k);
    PDstandard4th1At31 = PDstandard4th1(At31, i, j, k);
    PDstandard4th2At31 = PDstandard4th2(At31, i, j, k);
    PDstandard4th3At31 = PDstandard4th3(At31, i, j, k);
    PDstandard4th1At32 = PDstandard4th1(At32, i, j, k);
    PDstandard4th2At32 = PDstandard4th2(At32, i, j, k);
    PDstandard4th3At32 = PDstandard4th3(At32, i, j, k);
    PDstandard4th1At33 = PDstandard4th1(At33, i, j, k);
    PDstandard4th2At33 = PDstandard4th2(At33, i, j, k);
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
    detgt  =  2*gt21L*gt31L*gt32L + gt33L*(gt11L*gt22L - SQR(gt21L)) - gt22L*SQR(gt31L) - gt11L*SQR(gt32L);
    
    ddetgt1  =  (-2*gt22L*gt31L + 2*gt21L*gt32L)*PDstandard4th1gt31 + 
        2*((gt31L*gt32L - gt21L*gt33L)*PDstandard4th1gt21 + (gt21L*gt31L - gt11L*gt32L)*PDstandard4th1gt32) + 
        PDstandard4th1gt33*(gt11L*gt22L - SQR(gt21L)) + PDstandard4th1gt22*(gt11L*gt33L - SQR(gt31L)) + 
        PDstandard4th1gt11*(gt22L*gt33L - SQR(gt32L));
    
    ddetgt2  =  (-2*gt22L*gt31L + 2*gt21L*gt32L)*PDstandard4th2gt31 + 
        2*((gt31L*gt32L - gt21L*gt33L)*PDstandard4th2gt21 + (gt21L*gt31L - gt11L*gt32L)*PDstandard4th2gt32) + 
        PDstandard4th2gt33*(gt11L*gt22L - SQR(gt21L)) + PDstandard4th2gt22*(gt11L*gt33L - SQR(gt31L)) + 
        PDstandard4th2gt11*(gt22L*gt33L - SQR(gt32L));
    
    ddetgt3  =  (-2*gt22L*gt31L + 2*gt21L*gt32L)*PDstandard4th3gt31 + 
        2*((gt31L*gt32L - gt21L*gt33L)*PDstandard4th3gt21 + (gt21L*gt31L - gt11L*gt32L)*PDstandard4th3gt32) + 
        PDstandard4th3gt33*(gt11L*gt22L - SQR(gt21L)) + PDstandard4th3gt22*(gt11L*gt33L - SQR(gt31L)) + 
        PDstandard4th3gt11*(gt22L*gt33L - SQR(gt32L));
    
    gtu11  =  INV(detgt)*(gt22L*gt33L - SQR(gt32L));
    
    gtu12  =  (gt31L*gt32L - gt21L*gt33L)*INV(detgt);
    
    gtu13  =  (-(gt22L*gt31L) + gt21L*gt32L)*INV(detgt);
    
    gtu22  =  INV(detgt)*(gt11L*gt33L - SQR(gt31L));
    
    gtu23  =  (gt21L*gt31L - gt11L*gt32L)*INV(detgt);
    
    gtu33  =  INV(detgt)*(gt11L*gt22L - SQR(gt21L));
    
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
    
    e4phi  =  exp(4*phiL);
    
    g11  =  e4phi*gt11L;
    
    g21  =  e4phi*gt21L;
    
    g31  =  e4phi*gt31L;
    
    g22  =  e4phi*gt22L;
    
    g32  =  e4phi*gt32L;
    
    g33  =  e4phi*gt33L;
    
    detg  =  2*g21*g31*g32 + g33*(g11*g22 - SQR(g21)) - g22*SQR(g31) - g11*SQR(g32);
    
    gu11  =  INV(detg)*(g22*g33 - SQR(g32));
    
    gu12  =  (g31*g32 - g21*g33)*INV(detg);
    
    gu13  =  (-(g22*g31) + g21*g32)*INV(detg);
    
    gu22  =  INV(detg)*(g11*g33 - SQR(g31));
    
    gu23  =  (g21*g31 - g11*g32)*INV(detg);
    
    gu33  =  INV(detg)*(g11*g22 - SQR(g21));
    
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
    
    trR  =  gu11*R11 + gu22*R22 + 2*(gu12*R21 + gu13*R31 + gu23*R32) + gu33*R33;
    
    K11  =  At11L*e4phi + g11*kthird*trKL;
    
    K21  =  At21L*e4phi + g21*kthird*trKL;
    
    K31  =  At31L*e4phi + g31*kthird*trKL;
    
    K22  =  At22L*e4phi + g22*kthird*trKL;
    
    K32  =  At32L*e4phi + g32*kthird*trKL;
    
    K33  =  At33L*e4phi + g33*kthird*trKL;
    
    Km11  =  gu11*K11 + gu12*K21 + gu13*K31;
    
    Km21  =  gu12*K11 + gu22*K21 + gu23*K31;
    
    Km31  =  gu13*K11 + gu23*K21 + gu33*K31;
    
    Km12  =  gu11*K21 + gu12*K22 + gu13*K32;
    
    Km22  =  gu12*K21 + gu22*K22 + gu23*K32;
    
    Km32  =  gu13*K21 + gu23*K22 + gu33*K32;
    
    Km13  =  gu11*K31 + gu12*K32 + gu13*K33;
    
    Km23  =  gu12*K31 + gu22*K32 + gu23*K33;
    
    Km33  =  gu13*K31 + gu23*K32 + gu33*K33;
    
    HL  =  -2*(Km12*Km21 + Km13*Km31 + Km23*Km32) + trR - SQR(Km11) - SQR(Km22) - SQR(Km33) + SQR(trKL);
    
    gK112  =  e4phi*(-2*(At11L*G121 + At21L*G221 + At31L*G321) + PDstandard4th2At11 + 4*At11L*PDstandard4th2phi) + 
        g11*kthird*PDstandard4th2trK;
    
    gK113  =  e4phi*(-2*(At11L*G131 + At21L*G231 + At31L*G331) + PDstandard4th3At11 + 4*At11L*PDstandard4th3phi) + 
        g11*kthird*PDstandard4th3trK;
    
    gK211  =  e4phi*(-(At11L*G121) - At22L*G211 - At21L*(G111 + G221) - At32L*G311 - At31L*G321 + PDstandard4th1At21 + 
           4*At21L*PDstandard4th1phi) + g21*kthird*PDstandard4th1trK;
    
    gK212  =  e4phi*(-(At11L*G122) - At22L*G221 - At21L*(G121 + G222) - At32L*G321 - At31L*G322 + PDstandard4th2At21 + 
           4*At21L*PDstandard4th2phi) + g21*kthird*PDstandard4th2trK;
    
    gK213  =  e4phi*(-(At11L*G132) - At22L*G231 - At21L*(G131 + G232) - At32L*G331 - At31L*G332 + PDstandard4th3At21 + 
           4*At21L*PDstandard4th3phi) + g21*kthird*PDstandard4th3trK;
    
    gK311  =  e4phi*(-(At11L*G131) - At32L*G211 - At21L*G231 - At33L*G311 - At31L*(G111 + G331) + PDstandard4th1At31 + 
           4*At31L*PDstandard4th1phi) + g31*kthird*PDstandard4th1trK;
    
    gK312  =  e4phi*(-(At11L*G132) - At32L*G221 - At21L*G232 - At33L*G321 - At31L*(G121 + G332) + PDstandard4th2At31 + 
           4*At31L*PDstandard4th2phi) + g31*kthird*PDstandard4th2trK;
    
    gK313  =  e4phi*(-(At11L*G133) - At32L*G231 - At21L*G233 - At33L*G331 - At31L*(G131 + G333) + PDstandard4th3At31 + 
           4*At31L*PDstandard4th3phi) + g31*kthird*PDstandard4th3trK;
    
    gK221  =  e4phi*(-2*(At21L*G121 + At22L*G221 + At32L*G321) + PDstandard4th1At22 + 4*At22L*PDstandard4th1phi) + 
        g22*kthird*PDstandard4th1trK;
    
    gK223  =  e4phi*(-2*(At21L*G132 + At22L*G232 + At32L*G332) + PDstandard4th3At22 + 4*At22L*PDstandard4th3phi) + 
        g22*kthird*PDstandard4th3trK;
    
    gK321  =  e4phi*(-(At31L*G121) - At21L*G131 - At32L*G221 - At22L*G231 - At33L*G321 - At32L*G331 + PDstandard4th1At32 + 
           4*At32L*PDstandard4th1phi) + g32*kthird*PDstandard4th1trK;
    
    gK322  =  e4phi*(-(At31L*G122) - At21L*G132 - At32L*G222 - At22L*G232 - At33L*G322 - At32L*G332 + PDstandard4th2At32 + 
           4*At32L*PDstandard4th2phi) + g32*kthird*PDstandard4th2trK;
    
    gK323  =  e4phi*(-(At31L*G132) - At21L*G133 - At32L*G232 - At22L*G233 - At33L*G332 - At32L*G333 + PDstandard4th3At32 + 
           4*At32L*PDstandard4th3phi) + g32*kthird*PDstandard4th3trK;
    
    gK331  =  e4phi*(-2*(At31L*G131 + At32L*G231 + At33L*G331) + PDstandard4th1At33 + 4*At33L*PDstandard4th1phi) + 
        g33*kthird*PDstandard4th1trK;
    
    gK332  =  e4phi*(-2*(At31L*G132 + At32L*G232 + At33L*G332) + PDstandard4th2At33 + 4*At33L*PDstandard4th2phi) + 
        g33*kthird*PDstandard4th2trK;
    
    M1L  =  (gK112 - gK211)*gu12 + (gK113 - gK311)*gu13 + (gK212 - gK221)*gu22 + (gK213 + gK312 - 2*gK321)*gu23 + 
        (gK313 - gK331)*gu33;
    
    M2L  =  (-gK112 + gK211)*gu11 + (-gK212 + gK221)*gu12 + (gK213 - 2*gK312 + gK321)*gu13 + (gK223 - gK322)*gu23 + 
        (gK323 - gK332)*gu33;
    
    M3L  =  (-gK113 + gK311)*gu11 + (-2*gK213 + gK312 + gK321)*gu12 + (-gK313 + gK331)*gu13 + (-gK223 + gK322)*gu22 + 
        (-gK323 + gK332)*gu23;
    
    cSL  =  Log(detgt);
    
    cXt1L  =  Gt111*gtu11 + Gt122*gtu22 + 2*(Gt121*gtu12 + Gt131*gtu13 + Gt132*gtu23) + Gt133*gtu33 - Xt1L;
    
    cXt2L  =  Gt211*gtu11 + Gt222*gtu22 + 2*(Gt221*gtu12 + Gt231*gtu13 + Gt232*gtu23) + Gt233*gtu33 - Xt2L;
    
    cXt3L  =  Gt311*gtu11 + Gt322*gtu22 + 2*(Gt321*gtu12 + Gt331*gtu13 + Gt332*gtu23) + Gt333*gtu33 - Xt3L;
    
    cAL  =  At11L*gtu11 + At22L*gtu22 + 2*(At21L*gtu12 + At31L*gtu13 + At32L*gtu23) + At33L*gtu33;
    
    
    /* Copy local copies back to grid functions */
    cA[index] = cAL;
    cS[index] = cSL;
    cXt1[index] = cXt1L;
    cXt2[index] = cXt2L;
    cXt3[index] = cXt3L;
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (somename);
}

void BSSN_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverEverything(cctkGH, &BSSN_constraints_Body);
}
