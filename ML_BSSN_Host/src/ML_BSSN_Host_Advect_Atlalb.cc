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

extern "C" void ML_BSSN_Host_Advect_Atlalb_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_Host::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_Host::ML_curvrhs.");
  return;
}

static void ML_BSSN_Host_Advect_Atlalb_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  LC_LOOP3VEC(ML_BSSN_Host_Advect_Atlalb,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC At11L = vec_load(At11[index]);
    CCTK_REAL_VEC At11rhsL = vec_load(At11rhs[index]);
    CCTK_REAL_VEC At12L = vec_load(At12[index]);
    CCTK_REAL_VEC At12rhsL = vec_load(At12rhs[index]);
    CCTK_REAL_VEC At13L = vec_load(At13[index]);
    CCTK_REAL_VEC At13rhsL = vec_load(At13rhs[index]);
    CCTK_REAL_VEC At22L = vec_load(At22[index]);
    CCTK_REAL_VEC At22rhsL = vec_load(At22rhs[index]);
    CCTK_REAL_VEC At23L = vec_load(At23[index]);
    CCTK_REAL_VEC At23rhsL = vec_load(At23rhs[index]);
    CCTK_REAL_VEC At33L = vec_load(At33[index]);
    CCTK_REAL_VEC At33rhsL = vec_load(At33rhs[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC const PDupwindNthAnti1At11 = PDupwindNthAnti1(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At11 = PDupwindNthSymm1(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At11 = PDupwindNthAnti2(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At11 = PDupwindNthSymm2(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At11 = PDupwindNthAnti3(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At11 = PDupwindNthSymm3(&At11[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1At12 = PDupwindNthAnti1(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At12 = PDupwindNthSymm1(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At12 = PDupwindNthAnti2(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At12 = PDupwindNthSymm2(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At12 = PDupwindNthAnti3(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At12 = PDupwindNthSymm3(&At12[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1At13 = PDupwindNthAnti1(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At13 = PDupwindNthSymm1(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At13 = PDupwindNthAnti2(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At13 = PDupwindNthSymm2(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At13 = PDupwindNthAnti3(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At13 = PDupwindNthSymm3(&At13[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1At22 = PDupwindNthAnti1(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At22 = PDupwindNthSymm1(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At22 = PDupwindNthAnti2(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At22 = PDupwindNthSymm2(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At22 = PDupwindNthAnti3(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At22 = PDupwindNthSymm3(&At22[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1At23 = PDupwindNthAnti1(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At23 = PDupwindNthSymm1(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At23 = PDupwindNthAnti2(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At23 = PDupwindNthSymm2(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At23 = PDupwindNthAnti3(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At23 = PDupwindNthSymm3(&At23[index]);
    CCTK_REAL_VEC const PDupwindNthAnti1At33 = PDupwindNthAnti1(&At33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm1At33 = PDupwindNthSymm1(&At33[index]);
    CCTK_REAL_VEC const PDupwindNthAnti2At33 = PDupwindNthAnti2(&At33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm2At33 = PDupwindNthSymm2(&At33[index]);
    CCTK_REAL_VEC const PDupwindNthAnti3At33 = PDupwindNthAnti3(&At33[index]);
    CCTK_REAL_VEC const PDupwindNthSymm3At33 = PDupwindNthSymm3(&At33[index]);
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 = Sign(beta1L);
    
    ptrdiff_t dir2 = Sign(beta2L);
    
    ptrdiff_t dir3 = Sign(beta3L);
    
    At11rhsL = 
      kadd(At11rhsL,kmadd(beta1L,PDupwindNthAnti1At11,kmadd(beta2L,PDupwindNthAnti2At11,kmadd(beta3L,PDupwindNthAnti3At11,kmadd(PDupwindNthSymm1At11,kfabs(beta1L),kmadd(PDupwindNthSymm2At11,kfabs(beta2L),kmul(PDupwindNthSymm3At11,kfabs(beta3L))))))));
    
    At12rhsL = 
      kadd(At12rhsL,kmadd(beta1L,PDupwindNthAnti1At12,kmadd(beta2L,PDupwindNthAnti2At12,kmadd(beta3L,PDupwindNthAnti3At12,kmadd(PDupwindNthSymm1At12,kfabs(beta1L),kmadd(PDupwindNthSymm2At12,kfabs(beta2L),kmul(PDupwindNthSymm3At12,kfabs(beta3L))))))));
    
    At13rhsL = 
      kadd(At13rhsL,kmadd(beta1L,PDupwindNthAnti1At13,kmadd(beta2L,PDupwindNthAnti2At13,kmadd(beta3L,PDupwindNthAnti3At13,kmadd(PDupwindNthSymm1At13,kfabs(beta1L),kmadd(PDupwindNthSymm2At13,kfabs(beta2L),kmul(PDupwindNthSymm3At13,kfabs(beta3L))))))));
    
    At22rhsL = 
      kadd(At22rhsL,kmadd(beta1L,PDupwindNthAnti1At22,kmadd(beta2L,PDupwindNthAnti2At22,kmadd(beta3L,PDupwindNthAnti3At22,kmadd(PDupwindNthSymm1At22,kfabs(beta1L),kmadd(PDupwindNthSymm2At22,kfabs(beta2L),kmul(PDupwindNthSymm3At22,kfabs(beta3L))))))));
    
    At23rhsL = 
      kadd(At23rhsL,kmadd(beta1L,PDupwindNthAnti1At23,kmadd(beta2L,PDupwindNthAnti2At23,kmadd(beta3L,PDupwindNthAnti3At23,kmadd(PDupwindNthSymm1At23,kfabs(beta1L),kmadd(PDupwindNthSymm2At23,kfabs(beta2L),kmul(PDupwindNthSymm3At23,kfabs(beta3L))))))));
    
    At33rhsL = 
      kadd(At33rhsL,kmadd(beta1L,PDupwindNthAnti1At33,kmadd(beta2L,PDupwindNthAnti2At33,kmadd(beta3L,PDupwindNthAnti3At33,kmadd(PDupwindNthSymm1At33,kfabs(beta1L),kmadd(PDupwindNthSymm2At33,kfabs(beta2L),kmul(PDupwindNthSymm3At33,kfabs(beta3L))))))));
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,lc_imin,lc_imax);
    vec_store_nta_partial(At11rhs[index],At11rhsL);
    vec_store_nta_partial(At12rhs[index],At12rhsL);
    vec_store_nta_partial(At13rhs[index],At13rhsL);
    vec_store_nta_partial(At22rhs[index],At22rhsL);
    vec_store_nta_partial(At23rhs[index],At23rhsL);
    vec_store_nta_partial(At33rhs[index],At33rhsL);
  }
  LC_ENDLOOP3VEC(ML_BSSN_Host_Advect_Atlalb);
}

extern "C" void ML_BSSN_Host_Advect_Atlalb(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_Host_Advect_Atlalb_Body");
  }
  
  if (cctk_iteration % ML_BSSN_Host_Advect_Atlalb_calc_every != ML_BSSN_Host_Advect_Atlalb_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "ML_BSSN_Host::ML_curv",
    "ML_BSSN_Host::ML_curvrhs",
    "ML_BSSN_Host::ML_shift"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_Advect_Atlalb", 3, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_Host_Advect_Atlalb", 5, 5, 5);
  
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_Host_Advect_Atlalb_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_Advect_Atlalb_Body");
  }
}
