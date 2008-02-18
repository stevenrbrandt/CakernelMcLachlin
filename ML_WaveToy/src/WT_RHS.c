/*  File produced by user diener */
/*  Produced with Mathematica Version 6.0 for Linux x86 (32-bit) (April 20, 2007) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

#define KRANC_C

#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"
#include "Differencing.h"

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

void WT_RHS_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangentA[3], CCTK_REAL tangentB[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  
  /* Declare the variables used for looping over grid points */
  CCTK_INT i = INITVALUE, j = INITVALUE, k = INITVALUE;
  CCTK_INT index = INITVALUE;
  CCTK_INT subblock_index = INITVALUE;
  
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
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WT_RHS_Body");
  }
  
  if (cctk_iteration % WT_RHS_calc_every != WT_RHS_calc_offset)
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
  for (k = min[2]; k < max[2]; k++)
  {
    for (j = min[1]; j < max[1]; j++)
    {
      for (i = min[0]; i < max[0]; i++)
      {
         index  =  CCTK_GFINDEX3D(cctkGH,i,j,k) ;
         subblock_index  =  i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2])) ;
        
        /* Declare shorthands */
        CCTK_REAL dJinv111 = INITVALUE, dJinv122 = INITVALUE, dJinv133 = INITVALUE, dJinv211 = INITVALUE, dJinv222 = INITVALUE, dJinv233 = INITVALUE;
        CCTK_REAL dJinv311 = INITVALUE, dJinv322 = INITVALUE, dJinv333 = INITVALUE;
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
        CCTK_REAL ddadxdxL = INITVALUE;
        CCTK_REAL ddadydyL = INITVALUE;
        CCTK_REAL ddadzdzL = INITVALUE;
        CCTK_REAL rhoL = INITVALUE, rhorhsL = INITVALUE;
        CCTK_REAL uL = INITVALUE, urhsL = INITVALUE;
        /* Declare precomputed derivatives*/
        
        /* Declare derivatives */
        CCTK_REAL PDstandardNth1u = INITVALUE;
        CCTK_REAL PDstandardNth2u = INITVALUE;
        CCTK_REAL PDstandardNth3u = INITVALUE;
        CCTK_REAL PDstandardNth11u = INITVALUE;
        CCTK_REAL PDstandardNth22u = INITVALUE;
        CCTK_REAL PDstandardNth33u = INITVALUE;
        CCTK_REAL PDstandardNth12u = INITVALUE;
        CCTK_REAL PDstandardNth13u = INITVALUE;
        CCTK_REAL PDstandardNth21u = INITVALUE;
        CCTK_REAL PDstandardNth23u = INITVALUE;
        CCTK_REAL PDstandardNth31u = INITVALUE;
        CCTK_REAL PDstandardNth32u = INITVALUE;
        
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
        ddadxdxL = ddadxdx[index];
        ddadydyL = ddadydy[index];
        ddadzdzL = ddadzdz[index];
        rhoL = rho[index];
        uL = u[index];
        
        /* Assign local copies of subblock grid functions */
        
        /* Include user supplied include files */
        
        /* Precompute derivatives (new style) */
        PDstandardNth1u = PDstandardNth1(u, i, j, k);
        PDstandardNth2u = PDstandardNth2(u, i, j, k);
        PDstandardNth3u = PDstandardNth3(u, i, j, k);
        PDstandardNth11u = PDstandardNth11(u, i, j, k);
        PDstandardNth22u = PDstandardNth22(u, i, j, k);
        PDstandardNth33u = PDstandardNth33(u, i, j, k);
        PDstandardNth12u = PDstandardNth12(u, i, j, k);
        PDstandardNth13u = PDstandardNth13(u, i, j, k);
        PDstandardNth23u = PDstandardNth23(u, i, j, k);
        
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
        
        dJinv111  =  ddadxdxL;
        
        dJinv122  =  ddadydyL;
        
        dJinv133  =  ddadzdzL;
        
        dJinv211  =  ddadxdxL;
        
        dJinv222  =  ddadydyL;
        
        dJinv233  =  ddadzdzL;
        
        dJinv311  =  ddadxdxL;
        
        dJinv322  =  ddadydyL;
        
        dJinv333  =  ddadzdzL;
        
        urhsL  =  rhoL;
        
        rhorhsL  =  (dJinv111 + dJinv122 + dJinv133)*PDstandardNth1u + 
            2*((Jinv11*Jinv21 + Jinv12*Jinv22 + Jinv13*Jinv23)*PDstandardNth12u + 
               (Jinv11*Jinv31 + Jinv12*Jinv32 + Jinv13*Jinv33)*PDstandardNth13u + 
               (Jinv21*Jinv31 + Jinv22*Jinv32 + Jinv23*Jinv33)*PDstandardNth23u) + 
            (dJinv211 + dJinv222 + dJinv233)*PDstandardNth2u + (dJinv311 + dJinv322 + dJinv333)*PDstandardNth3u + 
            PDstandardNth11u*(SQR(Jinv11) + SQR(Jinv12) + SQR(Jinv13)) + 
            PDstandardNth22u*(SQR(Jinv21) + SQR(Jinv22) + SQR(Jinv23)) + 
            PDstandardNth33u*(SQR(Jinv31) + SQR(Jinv32) + SQR(Jinv33));
        
        
        /* Copy local copies back to grid functions */
        rhorhs[index] = rhorhsL;
        urhs[index] = urhsL;
        
        /* Copy local copies back to subblock grid functions */
      }
    }
  }
}

void WT_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &WT_RHS_Body);
}
