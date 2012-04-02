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
#define QAD(x) (SQR(SQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))

extern "C" void ML_BSSN_Host_Advect_gtlalb_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_metricrhs.");
  return;
}

static void ML_BSSN_Host_Advect_gtlalb_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_REAL_VEC const p1o1024dx = kmul(INV(dx),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1024dy = kmul(INV(dy),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1024dz = kmul(INV(dz),ToReal(0.0009765625));
  CCTK_REAL_VEC const p1o1680dx = kmul(INV(dx),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o1680dy = kmul(INV(dy),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o1680dz = kmul(INV(dz),ToReal(0.000595238095238095238095238095238));
  CCTK_REAL_VEC const p1o5040dx2 = kmul(INV(SQR(dx)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dy2 = kmul(INV(SQR(dy)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o5040dz2 = kmul(INV(SQR(dz)),ToReal(0.000198412698412698412698412698413));
  CCTK_REAL_VEC const p1o560dx = kmul(INV(dx),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o560dy = kmul(INV(dy),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o560dz = kmul(INV(dz),ToReal(0.00178571428571428571428571428571));
  CCTK_REAL_VEC const p1o705600dxdy = kmul(INV(kmul(dx,dy)),ToReal(1.41723356009070294784580498866e-6));
  CCTK_REAL_VEC const p1o705600dxdz = kmul(INV(kmul(dx,dz)),ToReal(1.41723356009070294784580498866e-6));
  CCTK_REAL_VEC const p1o705600dydz = kmul(INV(kmul(dy,dz)),ToReal(1.41723356009070294784580498866e-6));
  CCTK_REAL_VEC const p1o840dx = kmul(INV(dx),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dy = kmul(INV(dy),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dz = kmul(INV(dz),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1odx = INV(dx);
  CCTK_REAL_VEC const p1ody = INV(dy);
  CCTK_REAL_VEC const p1odz = INV(dz);
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC(ML_BSSN_Host_Advect_gtlalb,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC gt11rhsL = vec_load(gt11rhs[index]);
    CCTK_REAL_VEC gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC gt12rhsL = vec_load(gt12rhs[index]);
    CCTK_REAL_VEC gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC gt13rhsL = vec_load(gt13rhs[index]);
    CCTK_REAL_VEC gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC gt22rhsL = vec_load(gt22rhs[index]);
    CCTK_REAL_VEC gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC gt23rhsL = vec_load(gt23rhs[index]);
    CCTK_REAL_VEC gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC gt33rhsL = vec_load(gt33rhs[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC const PDupwindNthAnti1gt11 = PDupwindNthAnti1(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt11 = PDupwindNthSymm1(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt11 = PDupwindNthAnti2(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt11 = PDupwindNthSymm2(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt11 = PDupwindNthAnti3(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt11 = PDupwindNthSymm3(&gt11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1gt12 = PDupwindNthAnti1(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt12 = PDupwindNthSymm1(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt12 = PDupwindNthAnti2(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt12 = PDupwindNthSymm2(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt12 = PDupwindNthAnti3(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt12 = PDupwindNthSymm3(&gt12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1gt13 = PDupwindNthAnti1(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt13 = PDupwindNthSymm1(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt13 = PDupwindNthAnti2(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt13 = PDupwindNthSymm2(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt13 = PDupwindNthAnti3(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt13 = PDupwindNthSymm3(&gt13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1gt22 = PDupwindNthAnti1(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt22 = PDupwindNthSymm1(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt22 = PDupwindNthAnti2(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt22 = PDupwindNthSymm2(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt22 = PDupwindNthAnti3(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt22 = PDupwindNthSymm3(&gt22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1gt23 = PDupwindNthAnti1(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt23 = PDupwindNthSymm1(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt23 = PDupwindNthAnti2(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt23 = PDupwindNthSymm2(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt23 = PDupwindNthAnti3(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt23 = PDupwindNthSymm3(&gt23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1gt33 = PDupwindNthAnti1(&gt33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1gt33 = PDupwindNthSymm1(&gt33[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2gt33 = PDupwindNthAnti2(&gt33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2gt33 = PDupwindNthSymm2(&gt33[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3gt33 = PDupwindNthAnti3(&gt33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3gt33 = PDupwindNthSymm3(&gt33[index]);
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    gt11rhsL = 
      kadd(gt11rhsL,kmadd(beta1L,PDupwindNthAnti1gt11,kmadd(beta2L,PDupwindNthAnti2gt11,kmadd(beta3L,PDupwindNthAnti3gt11,kmadd(PDupwindNthSymm1gt11,kfabs(beta1L),kmadd(PDupwindNthSymm2gt11,kfabs(beta2L),kmul(PDupwindNthSymm3gt11,kfabs(beta3L))))))));
    
    gt12rhsL = 
      kadd(gt12rhsL,kmadd(beta1L,PDupwindNthAnti1gt12,kmadd(beta2L,PDupwindNthAnti2gt12,kmadd(beta3L,PDupwindNthAnti3gt12,kmadd(PDupwindNthSymm1gt12,kfabs(beta1L),kmadd(PDupwindNthSymm2gt12,kfabs(beta2L),kmul(PDupwindNthSymm3gt12,kfabs(beta3L))))))));
    
    gt13rhsL = 
      kadd(gt13rhsL,kmadd(beta1L,PDupwindNthAnti1gt13,kmadd(beta2L,PDupwindNthAnti2gt13,kmadd(beta3L,PDupwindNthAnti3gt13,kmadd(PDupwindNthSymm1gt13,kfabs(beta1L),kmadd(PDupwindNthSymm2gt13,kfabs(beta2L),kmul(PDupwindNthSymm3gt13,kfabs(beta3L))))))));
    
    gt22rhsL = 
      kadd(gt22rhsL,kmadd(beta1L,PDupwindNthAnti1gt22,kmadd(beta2L,PDupwindNthAnti2gt22,kmadd(beta3L,PDupwindNthAnti3gt22,kmadd(PDupwindNthSymm1gt22,kfabs(beta1L),kmadd(PDupwindNthSymm2gt22,kfabs(beta2L),kmul(PDupwindNthSymm3gt22,kfabs(beta3L))))))));
    
    gt23rhsL = 
      kadd(gt23rhsL,kmadd(beta1L,PDupwindNthAnti1gt23,kmadd(beta2L,PDupwindNthAnti2gt23,kmadd(beta3L,PDupwindNthAnti3gt23,kmadd(PDupwindNthSymm1gt23,kfabs(beta1L),kmadd(PDupwindNthSymm2gt23,kfabs(beta2L),kmul(PDupwindNthSymm3gt23,kfabs(beta3L))))))));
    
    gt33rhsL = 
      kadd(gt33rhsL,kmadd(beta1L,PDupwindNthAnti1gt33,kmadd(beta2L,PDupwindNthAnti2gt33,kmadd(beta3L,PDupwindNthAnti3gt33,kmadd(PDupwindNthSymm1gt33,kfabs(beta1L),kmadd(PDupwindNthSymm2gt33,kfabs(beta2L),kmul(PDupwindNthSymm3gt33,kfabs(beta3L))))))));
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,lc_imin,lc_imax);
    vec_store_nta_partial(gt11rhs[index],gt11rhsL);
    vec_store_nta_partial(gt12rhs[index],gt12rhsL);
    vec_store_nta_partial(gt13rhs[index],gt13rhsL);
    vec_store_nta_partial(gt22rhs[index],gt22rhsL);
    vec_store_nta_partial(gt23rhs[index],gt23rhsL);
    vec_store_nta_partial(gt33rhs[index],gt33rhsL);
  }
  LC_ENDLOOP3VEC(ML_BSSN_Host_Advect_gtlalb);
}

extern "C" void ML_BSSN_Host_Advect_gtlalb(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_Host_Advect_gtlalb_Body");
  }
  
  if (cctk_iteration % ML_BSSN_Host_Advect_gtlalb_calc_every != ML_BSSN_Host_Advect_gtlalb_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "ML_BSSN_Host::ML_metric",
    "ML_BSSN_Host::ML_metricrhs",
    "ML_BSSN_Host::ML_shift"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_Advect_gtlalb", 3, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Host_Advect_gtlalb", 5, 5, 5);
  
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_Host_Advect_gtlalb_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_Advect_gtlalb_Body");
  }
}
