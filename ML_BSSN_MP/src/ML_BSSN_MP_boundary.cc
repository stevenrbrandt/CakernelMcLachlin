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
#include "Vectors.hh"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

static void ML_BSSN_MP_boundary_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_boundary_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_boundary_calc_every != ML_BSSN_MP_boundary_calc_offset)
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
  LC_LOOP3 (ML_BSSN_MP_boundary,
            i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
            cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    // int index = INITVALUE;
    int const index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    /* Declare derivatives */
    
    /* Assign local copies of grid functions */
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC phiL = IfThen(conformalMethod,1,0);
    
    CCTK_REAL_VEC gt11L = 1;
    
    CCTK_REAL_VEC gt12L = 0;
    
    CCTK_REAL_VEC gt13L = 0;
    
    CCTK_REAL_VEC gt22L = 1;
    
    CCTK_REAL_VEC gt23L = 0;
    
    CCTK_REAL_VEC gt33L = 1;
    
    CCTK_REAL_VEC trKL = 0;
    
    CCTK_REAL_VEC At11L = 0;
    
    CCTK_REAL_VEC At12L = 0;
    
    CCTK_REAL_VEC At13L = 0;
    
    CCTK_REAL_VEC At22L = 0;
    
    CCTK_REAL_VEC At23L = 0;
    
    CCTK_REAL_VEC At33L = 0;
    
    CCTK_REAL_VEC Xt1L = 0;
    
    CCTK_REAL_VEC Xt2L = 0;
    
    CCTK_REAL_VEC Xt3L = 0;
    
    CCTK_REAL_VEC alphaL = 1;
    
    CCTK_REAL_VEC AL = 0;
    
    CCTK_REAL_VEC beta1L = 0;
    
    CCTK_REAL_VEC beta2L = 0;
    
    CCTK_REAL_VEC beta3L = 0;
    
    CCTK_REAL_VEC B1L = 0;
    
    CCTK_REAL_VEC B2L = 0;
    
    CCTK_REAL_VEC B3L = 0;
    
    
    /* Copy local copies back to grid functions */
    vec_store_nta(A[index],AL);
    vec_store_nta(alpha[index],alphaL);
    vec_store_nta(At11[index],At11L);
    vec_store_nta(At12[index],At12L);
    vec_store_nta(At13[index],At13L);
    vec_store_nta(At22[index],At22L);
    vec_store_nta(At23[index],At23L);
    vec_store_nta(At33[index],At33L);
    vec_store_nta(B1[index],B1L);
    vec_store_nta(B2[index],B2L);
    vec_store_nta(B3[index],B3L);
    vec_store_nta(beta1[index],beta1L);
    vec_store_nta(beta2[index],beta2L);
    vec_store_nta(beta3[index],beta3L);
    vec_store_nta(gt11[index],gt11L);
    vec_store_nta(gt12[index],gt12L);
    vec_store_nta(gt13[index],gt13L);
    vec_store_nta(gt22[index],gt22L);
    vec_store_nta(gt23[index],gt23L);
    vec_store_nta(gt33[index],gt33L);
    vec_store_nta(phi[index],phiL);
    vec_store_nta(trK[index],trKL);
    vec_store_nta(Xt1[index],Xt1L);
    vec_store_nta(Xt2[index],Xt2L);
    vec_store_nta(Xt3[index],Xt3L);
  i += CCTK_REAL_VEC_SIZE-1;
  }
  LC_ENDLOOP3 (ML_BSSN_MP_boundary);
}

extern "C" void ML_BSSN_MP_boundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverBoundaryWithGhosts(cctkGH, &ML_BSSN_MP_boundary_Body);
}
