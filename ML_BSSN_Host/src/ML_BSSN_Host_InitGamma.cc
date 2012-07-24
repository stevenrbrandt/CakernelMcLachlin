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
#include "vectors.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define ScalarINV(x) ((CCTK_REAL)1.0 / (x))
#define ScalarSQR(x) ((x) * (x))
#define ScalarCUB(x) ((x) * ScalarSQR(x))
#define ScalarQAD(x) (ScalarSQR(ScalarSQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))
#define QAD(x) (SQR(SQR(x)))

static void ML_BSSN_Host_InitGamma_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_REAL_VEC const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL_VEC const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL_VEC const dz = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL_VEC const dt = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL_VEC const t = ToReal(cctk_time);
  CCTK_REAL_VEC const dxi = INV(dx);
  CCTK_REAL_VEC const dyi = INV(dy);
  CCTK_REAL_VEC const dzi = INV(dz);
  CCTK_REAL_VEC const khalf = ToReal(0.5);
  CCTK_REAL_VEC const kthird = ToReal(1.0/3.0);
  CCTK_REAL_VEC const ktwothird = ToReal(2.0/3.0);
  CCTK_REAL_VEC const kfourthird = ToReal(4.0/3.0);
  CCTK_REAL_VEC const keightthird = ToReal(8.0/3.0);
  CCTK_REAL_VEC const hdxi = kmul(ToReal(0.5), dxi);
  CCTK_REAL_VEC const hdyi = kmul(ToReal(0.5), dyi);
  CCTK_REAL_VEC const hdzi = kmul(ToReal(0.5), dzi);
  
  /* Initialize predefined quantities */
  CCTK_REAL_VEC const p1o1024dx = kdiv(ToReal(0.0009765625),dx);
  CCTK_REAL_VEC const p1o1024dy = kdiv(ToReal(0.0009765625),dy);
  CCTK_REAL_VEC const p1o1024dz = kdiv(ToReal(0.0009765625),dz);
  CCTK_REAL_VEC const p1o1680dx = kdiv(ToReal(0.000595238095238095238095238095238),dx);
  CCTK_REAL_VEC const p1o1680dy = kdiv(ToReal(0.000595238095238095238095238095238),dy);
  CCTK_REAL_VEC const p1o1680dz = kdiv(ToReal(0.000595238095238095238095238095238),dz);
  CCTK_REAL_VEC const p1o2dx = kdiv(ToReal(0.5),dx);
  CCTK_REAL_VEC const p1o2dy = kdiv(ToReal(0.5),dy);
  CCTK_REAL_VEC const p1o2dz = kdiv(ToReal(0.5),dz);
  CCTK_REAL_VEC const p1o5040dx2 = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dx,dx));
  CCTK_REAL_VEC const p1o5040dy2 = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dy,dy));
  CCTK_REAL_VEC const p1o5040dz2 = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dz,dz));
  CCTK_REAL_VEC const p1o560dx = kdiv(ToReal(0.00178571428571428571428571428571),dx);
  CCTK_REAL_VEC const p1o560dy = kdiv(ToReal(0.00178571428571428571428571428571),dy);
  CCTK_REAL_VEC const p1o560dz = kdiv(ToReal(0.00178571428571428571428571428571),dz);
  CCTK_REAL_VEC const p1o705600dxdy = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dy,dx));
  CCTK_REAL_VEC const p1o705600dxdz = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dz,dx));
  CCTK_REAL_VEC const p1o705600dydz = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dz,dy));
  CCTK_REAL_VEC const p1o840dx = kdiv(ToReal(0.00119047619047619047619047619048),dx);
  CCTK_REAL_VEC const p1o840dy = kdiv(ToReal(0.00119047619047619047619047619048),dy);
  CCTK_REAL_VEC const p1o840dz = kdiv(ToReal(0.00119047619047619047619047619048),dz);
  CCTK_REAL_VEC const p1odx = kdiv(ToReal(1),dx);
  CCTK_REAL_VEC const p1ody = kdiv(ToReal(1),dy);
  CCTK_REAL_VEC const p1odz = kdiv(ToReal(1),dz);
  CCTK_REAL_VEC const pm1o2dx = kdiv(ToReal(-0.5),dx);
  CCTK_REAL_VEC const pm1o2dy = kdiv(ToReal(-0.5),dy);
  CCTK_REAL_VEC const pm1o2dz = kdiv(ToReal(-0.5),dz);
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC(ML_BSSN_Host_InitGamma,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC Xt1L = ToReal(0);
    
    CCTK_REAL_VEC Xt2L = ToReal(0);
    
    CCTK_REAL_VEC Xt3L = ToReal(0);
    
    CCTK_REAL_VEC AL = ToReal(0);
    
    CCTK_REAL_VEC B1L = ToReal(0);
    
    CCTK_REAL_VEC B2L = ToReal(0);
    
    CCTK_REAL_VEC B3L = ToReal(0);
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,lc_imin,lc_imax);
    vec_store_nta_partial(A[index],AL);
    vec_store_nta_partial(B1[index],B1L);
    vec_store_nta_partial(B2[index],B2L);
    vec_store_nta_partial(B3[index],B3L);
    vec_store_nta_partial(Xt1[index],Xt1L);
    vec_store_nta_partial(Xt2[index],Xt2L);
    vec_store_nta_partial(Xt3[index],Xt3L);
  }
  LC_ENDLOOP3VEC(ML_BSSN_Host_InitGamma);
}

extern "C" void ML_BSSN_Host_InitGamma(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_Host_InitGamma_Body");
  }
  
  if (cctk_iteration % ML_BSSN_Host_InitGamma_calc_every != ML_BSSN_Host_InitGamma_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "ML_BSSN_Host::ML_dtlapse",
    "ML_BSSN_Host::ML_dtshift",
    "ML_BSSN_Host::ML_Gamma"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_InitGamma", 3, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, ML_BSSN_Host_InitGamma_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_InitGamma_Body");
  }
}
