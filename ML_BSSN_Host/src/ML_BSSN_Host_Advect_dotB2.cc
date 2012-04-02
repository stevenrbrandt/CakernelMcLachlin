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

extern "C" void ML_BSSN_Host_Advect_dotB2_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_dtshiftrhs.");
  return;
}

static void ML_BSSN_Host_Advect_dotB2_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  LC_LOOP3VEC(ML_BSSN_Host_Advect_dotB2,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC B2L = vec_load(B2[index]);
    CCTK_REAL_VEC B2rhsL = vec_load(B2rhs[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC Xt2L = vec_load(Xt2[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC const PDupwindNthAnti1B2 = PDupwindNthAnti1(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1B2 = PDupwindNthSymm1(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2B2 = PDupwindNthAnti2(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2B2 = PDupwindNthSymm2(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3B2 = PDupwindNthAnti3(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3B2 = PDupwindNthSymm3(&B2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1Xt2 = PDupwindNthAnti1(&Xt2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1Xt2 = PDupwindNthSymm1(&Xt2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2Xt2 = PDupwindNthAnti2(&Xt2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2Xt2 = PDupwindNthSymm2(&Xt2[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3Xt2 = PDupwindNthAnti3(&Xt2[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3Xt2 = PDupwindNthSymm3(&Xt2[index]);
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    B2rhsL = 
      kmadd(kmsub(kmadd(beta1L,PDupwindNthAnti1B2,kmadd(beta2L,PDupwindNthAnti2B2,kmadd(beta3L,PDupwindNthAnti3B2,kmadd(PDupwindNthSymm1B2,kfabs(beta1L),kmadd(PDupwindNthSymm2B2,kfabs(beta2L),kmul(PDupwindNthSymm3B2,kfabs(beta3L))))))),ToReal(ShiftAdvectionCoeff),kmul(kmadd(beta1L,PDupwindNthAnti1Xt2,kmadd(beta2L,PDupwindNthAnti2Xt2,kmadd(beta3L,PDupwindNthAnti3Xt2,kmadd(PDupwindNthSymm1Xt2,kfabs(beta1L),kmadd(PDupwindNthSymm2Xt2,kfabs(beta2L),kmul(PDupwindNthSymm3Xt2,kfabs(beta3L))))))),kadd(ToReal(-1),ToReal(ShiftAdvectionCoeff)))),ToReal(ShiftBCoeff),B2rhsL);
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,lc_imin,lc_imax);
    vec_store_nta_partial(B2rhs[index],B2rhsL);
  }
  LC_ENDLOOP3VEC(ML_BSSN_Host_Advect_dotB2);
}

extern "C" void ML_BSSN_Host_Advect_dotB2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_Host_Advect_dotB2_Body");
  }
  
  if (cctk_iteration % ML_BSSN_Host_Advect_dotB2_calc_every != ML_BSSN_Host_Advect_dotB2_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "ML_BSSN_Host::ML_dtshift",
    "ML_BSSN_Host::ML_dtshiftrhs",
    "ML_BSSN_Host::ML_Gamma",
    "ML_BSSN_Host::ML_shift"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_Advect_dotB2", 4, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Host_Advect_dotB2", 5, 5, 5);
  
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_Host_Advect_dotB2_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_Advect_dotB2_Body");
  }
}
