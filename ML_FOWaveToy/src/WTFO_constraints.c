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

void WTFO_constraints_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
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
  CCTK_REAL p1o144dx2dy = INITVALUE;
  CCTK_REAL p1o144dx2dz = INITVALUE;
  CCTK_REAL p1o144dxdy = INITVALUE;
  CCTK_REAL p1o144dxdy2 = INITVALUE;
  CCTK_REAL p1o144dxdz = INITVALUE;
  CCTK_REAL p1o144dxdz2 = INITVALUE;
  CCTK_REAL p1o144dy2dz = INITVALUE;
  CCTK_REAL p1o144dydz = INITVALUE;
  CCTK_REAL p1o144dydz2 = INITVALUE;
  CCTK_REAL p1o1728dxdydz = INITVALUE;
  CCTK_REAL p1o2dx3 = INITVALUE;
  CCTK_REAL p1o2dy3 = INITVALUE;
  CCTK_REAL p1o2dz3 = INITVALUE;
  CCTK_REAL pm1o12dx2 = INITVALUE;
  CCTK_REAL pm1o12dy2 = INITVALUE;
  CCTK_REAL pm1o12dz2 = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WTFO_constraints_Body");
  }
  
  if (cctk_iteration % WTFO_constraints_calc_every != WTFO_constraints_calc_offset)
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
  p1o144dx2dy = (INV(dy)*pow(dx,-2))/144.;
  p1o144dx2dz = (INV(dz)*pow(dx,-2))/144.;
  p1o144dxdy = (INV(dx)*INV(dy))/144.;
  p1o144dxdy2 = (INV(dx)*pow(dy,-2))/144.;
  p1o144dxdz = (INV(dx)*INV(dz))/144.;
  p1o144dxdz2 = (INV(dx)*pow(dz,-2))/144.;
  p1o144dy2dz = (INV(dz)*pow(dy,-2))/144.;
  p1o144dydz = (INV(dy)*INV(dz))/144.;
  p1o144dydz2 = (INV(dy)*pow(dz,-2))/144.;
  p1o1728dxdydz = (INV(dx)*INV(dy)*INV(dz))/1728.;
  p1o2dx3 = khalf*pow(dx,-3);
  p1o2dy3 = khalf*pow(dy,-3);
  p1o2dz3 = khalf*pow(dz,-3);
  pm1o12dx2 = -pow(dx,-2)/12.;
  pm1o12dy2 = -pow(dy,-2)/12.;
  pm1o12dz2 = -pow(dz,-2)/12.;
  
  /* Loop over the grid points */
  _Pragma ("omp parallel")
  LC_LOOP3 (WTFO_constraints,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    int index = INITVALUE;
    int subblock_index = INITVALUE;
    index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    subblock_index = i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2]));
    
    /* Declare shorthands */
    CCTK_REAL Jinv11 = INITVALUE, Jinv12 = INITVALUE, Jinv13 = INITVALUE, Jinv21 = INITVALUE, Jinv22 = INITVALUE, Jinv23 = INITVALUE;
    CCTK_REAL Jinv31 = INITVALUE, Jinv32 = INITVALUE, Jinv33 = INITVALUE;
    
    /* Declare local copies of grid functions */
    CCTK_REAL dadxL = INITVALUE;
    CCTK_REAL dadyL = INITVALUE;
    CCTK_REAL dadzL = INITVALUE;
    CCTK_REAL dbdxL = INITVALUE;
    CCTK_REAL dbdyL = INITVALUE;
    CCTK_REAL dbdzL = INITVALUE;
    CCTK_REAL dcdxL = INITVALUE;
    CCTK_REAL dcdyL = INITVALUE;
    CCTK_REAL dcdzL = INITVALUE;
    CCTK_REAL v1L = INITVALUE, v2L = INITVALUE, v3L = INITVALUE;
    CCTK_REAL w1L = INITVALUE, w2L = INITVALUE, w3L = INITVALUE;
    /* Declare precomputed derivatives*/
    
    /* Declare derivatives */
    CCTK_REAL PDstandardNth1v1 = INITVALUE;
    CCTK_REAL PDstandardNth2v1 = INITVALUE;
    CCTK_REAL PDstandardNth3v1 = INITVALUE;
    CCTK_REAL PDstandardNth1v2 = INITVALUE;
    CCTK_REAL PDstandardNth2v2 = INITVALUE;
    CCTK_REAL PDstandardNth3v2 = INITVALUE;
    CCTK_REAL PDstandardNth1v3 = INITVALUE;
    CCTK_REAL PDstandardNth2v3 = INITVALUE;
    CCTK_REAL PDstandardNth3v3 = INITVALUE;
    
    /* Assign local copies of grid functions */
    dadxL = dadx[index];
    dadyL = dady[index];
    dadzL = dadz[index];
    dbdxL = dbdx[index];
    dbdyL = dbdy[index];
    dbdzL = dbdz[index];
    dcdxL = dcdx[index];
    dcdyL = dcdy[index];
    dcdzL = dcdz[index];
    v1L = v1[index];
    v2L = v2[index];
    v3L = v3[index];
    
    /* Assign local copies of subblock grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives (new style) */
    PDstandardNth1v1 = PDstandardNth1(v1, i, j, k);
    PDstandardNth2v1 = PDstandardNth2(v1, i, j, k);
    PDstandardNth3v1 = PDstandardNth3(v1, i, j, k);
    PDstandardNth1v2 = PDstandardNth1(v2, i, j, k);
    PDstandardNth2v2 = PDstandardNth2(v2, i, j, k);
    PDstandardNth3v2 = PDstandardNth3(v2, i, j, k);
    PDstandardNth1v3 = PDstandardNth1(v3, i, j, k);
    PDstandardNth2v3 = PDstandardNth2(v3, i, j, k);
    PDstandardNth3v3 = PDstandardNth3(v3, i, j, k);
    
    /* Precompute derivatives (old style) */
    
    /* Calculate temporaries and grid functions */
    Jinv11  =  dadxL;
    
    Jinv12  =  dadyL;
    
    Jinv13  =  dadzL;
    
    Jinv21  =  dbdxL;
    
    Jinv22  =  dbdyL;
    
    Jinv23  =  dbdzL;
    
    Jinv31  =  dcdxL;
    
    Jinv32  =  dcdyL;
    
    Jinv33  =  dcdzL;
    
    w1L  =  Jinv13*PDstandardNth1v2 - Jinv12*PDstandardNth1v3 + Jinv23*PDstandardNth2v2 - Jinv22*PDstandardNth2v3 + 
        Jinv33*PDstandardNth3v2 - Jinv32*PDstandardNth3v3;
    
    w2L  =  -(Jinv13*PDstandardNth1v1) + Jinv11*PDstandardNth1v3 - Jinv23*PDstandardNth2v1 + Jinv21*PDstandardNth2v3 - 
        Jinv33*PDstandardNth3v1 + Jinv31*PDstandardNth3v3;
    
    w3L  =  Jinv12*PDstandardNth1v1 - Jinv11*PDstandardNth1v2 + Jinv22*PDstandardNth2v1 - Jinv21*PDstandardNth2v2 + 
        Jinv32*PDstandardNth3v1 - Jinv31*PDstandardNth3v2;
    
    
    /* Copy local copies back to grid functions */
    w1[index] = w1L;
    w2[index] = w2L;
    w3[index] = w3L;
    
    /* Copy local copies back to subblock grid functions */
  }
  LC_ENDLOOP3 (WTFO_constraints);
}

void WTFO_constraints(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &WTFO_constraints_Body);
}
