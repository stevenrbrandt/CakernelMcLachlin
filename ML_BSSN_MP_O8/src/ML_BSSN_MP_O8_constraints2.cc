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
#include "loopcontrol.h"
#include "vectors.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define QAD(x) (SQR(SQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))

extern "C" void ML_BSSN_MP_O8_constraints2_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_cons_detg","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_cons_detg.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_cons_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_cons_Gamma.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_cons_traceA","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_cons_traceA.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_MP_O8::ML_mom","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_MP_O8::ML_mom.");
  return;
}

static void ML_BSSN_MP_O8_constraints2_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_MP_O8_constraints2_Body");
  }
  
  if (cctk_iteration % ML_BSSN_MP_O8_constraints2_calc_every != ML_BSSN_MP_O8_constraints2_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"ML_BSSN_MP_O8::ML_cons_detg","ML_BSSN_MP_O8::ML_cons_Gamma","ML_BSSN_MP_O8::ML_cons_traceA","ML_BSSN_MP_O8::ML_curv","ML_BSSN_MP_O8::ML_Gamma","ML_BSSN_MP_O8::ML_lapse","ML_BSSN_MP_O8::ML_log_confac","ML_BSSN_MP_O8::ML_metric","ML_BSSN_MP_O8::ML_mom","ML_BSSN_MP_O8::ML_shift","ML_BSSN_MP_O8::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_MP_O8_constraints2", 11, groups);
  
  GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_MP_O8_constraints2", 4, 4, 4);
  
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
  CCTK_REAL_VEC const p1o705600dxdy = kmul(INV(dx),kmul(INV(dy),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o705600dxdz = kmul(INV(dx),kmul(INV(dz),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o705600dydz = kmul(INV(dy),kmul(INV(dz),ToReal(1.41723356009070294784580498866e-6)));
  CCTK_REAL_VEC const p1o840dx = kmul(INV(dx),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dy = kmul(INV(dy),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1o840dz = kmul(INV(dz),ToReal(0.00119047619047619047619047619048));
  CCTK_REAL_VEC const p1odx = INV(dx);
  CCTK_REAL_VEC const p1ody = INV(dy);
  CCTK_REAL_VEC const p1odz = INV(dz);
  
  /* Jacobian variable pointers */
  bool const use_jacobian = (!CCTK_IsFunctionAliased("MultiPatch_GetMap") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)
                       && strlen(jacobian_group) > 0;
  if (use_jacobian && strlen(jacobian_derivative_group) == 0)
  {
    CCTK_WARN (1, "GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names");
  }
  
  CCTK_REAL const *restrict jacobian_ptrs[9];
  if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_group,
                                                9, jacobian_ptrs);
  
  CCTK_REAL const *restrict const J11 = use_jacobian ? jacobian_ptrs[0] : 0;
  CCTK_REAL const *restrict const J12 = use_jacobian ? jacobian_ptrs[1] : 0;
  CCTK_REAL const *restrict const J13 = use_jacobian ? jacobian_ptrs[2] : 0;
  CCTK_REAL const *restrict const J21 = use_jacobian ? jacobian_ptrs[3] : 0;
  CCTK_REAL const *restrict const J22 = use_jacobian ? jacobian_ptrs[4] : 0;
  CCTK_REAL const *restrict const J23 = use_jacobian ? jacobian_ptrs[5] : 0;
  CCTK_REAL const *restrict const J31 = use_jacobian ? jacobian_ptrs[6] : 0;
  CCTK_REAL const *restrict const J32 = use_jacobian ? jacobian_ptrs[7] : 0;
  CCTK_REAL const *restrict const J33 = use_jacobian ? jacobian_ptrs[8] : 0;
  
  CCTK_REAL const *restrict jacobian_derivative_ptrs[18];
  if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_derivative_group,
                                                18, jacobian_derivative_ptrs);
  
  CCTK_REAL const *restrict const dJ111 = use_jacobian ? jacobian_derivative_ptrs[0] : 0;
  CCTK_REAL const *restrict const dJ112 = use_jacobian ? jacobian_derivative_ptrs[1] : 0;
  CCTK_REAL const *restrict const dJ113 = use_jacobian ? jacobian_derivative_ptrs[2] : 0;
  CCTK_REAL const *restrict const dJ122 = use_jacobian ? jacobian_derivative_ptrs[3] : 0;
  CCTK_REAL const *restrict const dJ123 = use_jacobian ? jacobian_derivative_ptrs[4] : 0;
  CCTK_REAL const *restrict const dJ133 = use_jacobian ? jacobian_derivative_ptrs[5] : 0;
  CCTK_REAL const *restrict const dJ211 = use_jacobian ? jacobian_derivative_ptrs[6] : 0;
  CCTK_REAL const *restrict const dJ212 = use_jacobian ? jacobian_derivative_ptrs[7] : 0;
  CCTK_REAL const *restrict const dJ213 = use_jacobian ? jacobian_derivative_ptrs[8] : 0;
  CCTK_REAL const *restrict const dJ222 = use_jacobian ? jacobian_derivative_ptrs[9] : 0;
  CCTK_REAL const *restrict const dJ223 = use_jacobian ? jacobian_derivative_ptrs[10] : 0;
  CCTK_REAL const *restrict const dJ233 = use_jacobian ? jacobian_derivative_ptrs[11] : 0;
  CCTK_REAL const *restrict const dJ311 = use_jacobian ? jacobian_derivative_ptrs[12] : 0;
  CCTK_REAL const *restrict const dJ312 = use_jacobian ? jacobian_derivative_ptrs[13] : 0;
  CCTK_REAL const *restrict const dJ313 = use_jacobian ? jacobian_derivative_ptrs[14] : 0;
  CCTK_REAL const *restrict const dJ322 = use_jacobian ? jacobian_derivative_ptrs[15] : 0;
  CCTK_REAL const *restrict const dJ323 = use_jacobian ? jacobian_derivative_ptrs[16] : 0;
  CCTK_REAL const *restrict const dJ333 = use_jacobian ? jacobian_derivative_ptrs[17] : 0;
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC (ML_BSSN_MP_O8_constraints2,
    i,j,k, min[0],min[1],min[2], max[0],max[1],max[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL = vec_load(alpha[index]);
    CCTK_REAL_VEC At11L = vec_load(At11[index]);
    CCTK_REAL_VEC At12L = vec_load(At12[index]);
    CCTK_REAL_VEC At13L = vec_load(At13[index]);
    CCTK_REAL_VEC At22L = vec_load(At22[index]);
    CCTK_REAL_VEC At23L = vec_load(At23[index]);
    CCTK_REAL_VEC At33L = vec_load(At33[index]);
    CCTK_REAL_VEC beta1L = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L = vec_load(beta3[index]);
    CCTK_REAL_VEC gt11L = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L = vec_load(gt33[index]);
    CCTK_REAL_VEC phiL = vec_load(phi[index]);
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    CCTK_REAL_VEC Xt1L = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L = vec_load(Xt3[index]);
    
    CCTK_REAL_VEC eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL;
    
    if (*stress_energy_state)
    {
      eTtxL = vec_load(eTtx[index]);
      eTtyL = vec_load(eTty[index]);
      eTtzL = vec_load(eTtz[index]);
      eTxxL = vec_load(eTxx[index]);
      eTxyL = vec_load(eTxy[index]);
      eTxzL = vec_load(eTxz[index]);
      eTyyL = vec_load(eTyy[index]);
      eTyzL = vec_load(eTyz[index]);
      eTzzL = vec_load(eTzz[index]);
    }
    else
    {
      eTtxL = ToReal(0.0);
      eTtyL = ToReal(0.0);
      eTtzL = ToReal(0.0);
      eTxxL = ToReal(0.0);
      eTxyL = ToReal(0.0);
      eTxzL = ToReal(0.0);
      eTyyL = ToReal(0.0);
      eTyzL = ToReal(0.0);
      eTzzL = ToReal(0.0);
    }
    
    CCTK_REAL_VEC J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L;
    
    if (use_jacobian)
    {
      J11L = vec_load(J11[index]);
      J12L = vec_load(J12[index]);
      J13L = vec_load(J13[index]);
      J21L = vec_load(J21[index]);
      J22L = vec_load(J22[index]);
      J23L = vec_load(J23[index]);
      J31L = vec_load(J31[index]);
      J32L = vec_load(J32[index]);
      J33L = vec_load(J33[index]);
    }
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC const PDstandardNth1At11 = PDstandardNth1(&At11[index]);
    CCTK_REAL_VEC const PDstandardNth2At11 = PDstandardNth2(&At11[index]);
    CCTK_REAL_VEC const PDstandardNth3At11 = PDstandardNth3(&At11[index]);
    CCTK_REAL_VEC const PDstandardNth1At12 = PDstandardNth1(&At12[index]);
    CCTK_REAL_VEC const PDstandardNth2At12 = PDstandardNth2(&At12[index]);
    CCTK_REAL_VEC const PDstandardNth3At12 = PDstandardNth3(&At12[index]);
    CCTK_REAL_VEC const PDstandardNth1At13 = PDstandardNth1(&At13[index]);
    CCTK_REAL_VEC const PDstandardNth2At13 = PDstandardNth2(&At13[index]);
    CCTK_REAL_VEC const PDstandardNth3At13 = PDstandardNth3(&At13[index]);
    CCTK_REAL_VEC const PDstandardNth1At22 = PDstandardNth1(&At22[index]);
    CCTK_REAL_VEC const PDstandardNth2At22 = PDstandardNth2(&At22[index]);
    CCTK_REAL_VEC const PDstandardNth3At22 = PDstandardNth3(&At22[index]);
    CCTK_REAL_VEC const PDstandardNth1At23 = PDstandardNth1(&At23[index]);
    CCTK_REAL_VEC const PDstandardNth2At23 = PDstandardNth2(&At23[index]);
    CCTK_REAL_VEC const PDstandardNth3At23 = PDstandardNth3(&At23[index]);
    CCTK_REAL_VEC const PDstandardNth1At33 = PDstandardNth1(&At33[index]);
    CCTK_REAL_VEC const PDstandardNth2At33 = PDstandardNth2(&At33[index]);
    CCTK_REAL_VEC const PDstandardNth3At33 = PDstandardNth3(&At33[index]);
    CCTK_REAL_VEC const PDstandardNth1gt11 = PDstandardNth1(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth2gt11 = PDstandardNth2(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth3gt11 = PDstandardNth3(&gt11[index]);
    CCTK_REAL_VEC const PDstandardNth1gt12 = PDstandardNth1(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth2gt12 = PDstandardNth2(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth3gt12 = PDstandardNth3(&gt12[index]);
    CCTK_REAL_VEC const PDstandardNth1gt13 = PDstandardNth1(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth2gt13 = PDstandardNth2(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth3gt13 = PDstandardNth3(&gt13[index]);
    CCTK_REAL_VEC const PDstandardNth1gt22 = PDstandardNth1(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth2gt22 = PDstandardNth2(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth3gt22 = PDstandardNth3(&gt22[index]);
    CCTK_REAL_VEC const PDstandardNth1gt23 = PDstandardNth1(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth2gt23 = PDstandardNth2(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth3gt23 = PDstandardNth3(&gt23[index]);
    CCTK_REAL_VEC const PDstandardNth1gt33 = PDstandardNth1(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth2gt33 = PDstandardNth2(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth3gt33 = PDstandardNth3(&gt33[index]);
    CCTK_REAL_VEC const PDstandardNth1phi = PDstandardNth1(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth2phi = PDstandardNth2(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth3phi = PDstandardNth3(&phi[index]);
    CCTK_REAL_VEC const PDstandardNth1trK = PDstandardNth1(&trK[index]);
    CCTK_REAL_VEC const PDstandardNth2trK = PDstandardNth2(&trK[index]);
    CCTK_REAL_VEC const PDstandardNth3trK = PDstandardNth3(&trK[index]);
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC JacPDstandardNth1At11;
    CCTK_REAL_VEC JacPDstandardNth1At12;
    CCTK_REAL_VEC JacPDstandardNth1At13;
    CCTK_REAL_VEC JacPDstandardNth1At22;
    CCTK_REAL_VEC JacPDstandardNth1At23;
    CCTK_REAL_VEC JacPDstandardNth1At33;
    CCTK_REAL_VEC JacPDstandardNth1gt11;
    CCTK_REAL_VEC JacPDstandardNth1gt12;
    CCTK_REAL_VEC JacPDstandardNth1gt13;
    CCTK_REAL_VEC JacPDstandardNth1gt22;
    CCTK_REAL_VEC JacPDstandardNth1gt23;
    CCTK_REAL_VEC JacPDstandardNth1gt33;
    CCTK_REAL_VEC JacPDstandardNth1phi;
    CCTK_REAL_VEC JacPDstandardNth1trK;
    CCTK_REAL_VEC JacPDstandardNth2At11;
    CCTK_REAL_VEC JacPDstandardNth2At12;
    CCTK_REAL_VEC JacPDstandardNth2At13;
    CCTK_REAL_VEC JacPDstandardNth2At22;
    CCTK_REAL_VEC JacPDstandardNth2At23;
    CCTK_REAL_VEC JacPDstandardNth2At33;
    CCTK_REAL_VEC JacPDstandardNth2gt11;
    CCTK_REAL_VEC JacPDstandardNth2gt12;
    CCTK_REAL_VEC JacPDstandardNth2gt13;
    CCTK_REAL_VEC JacPDstandardNth2gt22;
    CCTK_REAL_VEC JacPDstandardNth2gt23;
    CCTK_REAL_VEC JacPDstandardNth2gt33;
    CCTK_REAL_VEC JacPDstandardNth2phi;
    CCTK_REAL_VEC JacPDstandardNth2trK;
    CCTK_REAL_VEC JacPDstandardNth3At11;
    CCTK_REAL_VEC JacPDstandardNth3At12;
    CCTK_REAL_VEC JacPDstandardNth3At13;
    CCTK_REAL_VEC JacPDstandardNth3At22;
    CCTK_REAL_VEC JacPDstandardNth3At23;
    CCTK_REAL_VEC JacPDstandardNth3At33;
    CCTK_REAL_VEC JacPDstandardNth3gt11;
    CCTK_REAL_VEC JacPDstandardNth3gt12;
    CCTK_REAL_VEC JacPDstandardNth3gt13;
    CCTK_REAL_VEC JacPDstandardNth3gt22;
    CCTK_REAL_VEC JacPDstandardNth3gt23;
    CCTK_REAL_VEC JacPDstandardNth3gt33;
    CCTK_REAL_VEC JacPDstandardNth3phi;
    CCTK_REAL_VEC JacPDstandardNth3trK;
    
    if (use_jacobian)
    {
      JacPDstandardNth1At11 = 
        kmadd(J11L,PDstandardNth1At11,kmadd(J21L,PDstandardNth2At11,kmul(J31L,PDstandardNth3At11)));
      
      JacPDstandardNth1At12 = 
        kmadd(J11L,PDstandardNth1At12,kmadd(J21L,PDstandardNth2At12,kmul(J31L,PDstandardNth3At12)));
      
      JacPDstandardNth1At13 = 
        kmadd(J11L,PDstandardNth1At13,kmadd(J21L,PDstandardNth2At13,kmul(J31L,PDstandardNth3At13)));
      
      JacPDstandardNth1At22 = 
        kmadd(J11L,PDstandardNth1At22,kmadd(J21L,PDstandardNth2At22,kmul(J31L,PDstandardNth3At22)));
      
      JacPDstandardNth1At23 = 
        kmadd(J11L,PDstandardNth1At23,kmadd(J21L,PDstandardNth2At23,kmul(J31L,PDstandardNth3At23)));
      
      JacPDstandardNth1At33 = 
        kmadd(J11L,PDstandardNth1At33,kmadd(J21L,PDstandardNth2At33,kmul(J31L,PDstandardNth3At33)));
      
      JacPDstandardNth1gt11 = 
        kmadd(J11L,PDstandardNth1gt11,kmadd(J21L,PDstandardNth2gt11,kmul(J31L,PDstandardNth3gt11)));
      
      JacPDstandardNth1gt12 = 
        kmadd(J11L,PDstandardNth1gt12,kmadd(J21L,PDstandardNth2gt12,kmul(J31L,PDstandardNth3gt12)));
      
      JacPDstandardNth1gt13 = 
        kmadd(J11L,PDstandardNth1gt13,kmadd(J21L,PDstandardNth2gt13,kmul(J31L,PDstandardNth3gt13)));
      
      JacPDstandardNth1gt22 = 
        kmadd(J11L,PDstandardNth1gt22,kmadd(J21L,PDstandardNth2gt22,kmul(J31L,PDstandardNth3gt22)));
      
      JacPDstandardNth1gt23 = 
        kmadd(J11L,PDstandardNth1gt23,kmadd(J21L,PDstandardNth2gt23,kmul(J31L,PDstandardNth3gt23)));
      
      JacPDstandardNth1gt33 = 
        kmadd(J11L,PDstandardNth1gt33,kmadd(J21L,PDstandardNth2gt33,kmul(J31L,PDstandardNth3gt33)));
      
      JacPDstandardNth1phi = 
        kmadd(J11L,PDstandardNth1phi,kmadd(J21L,PDstandardNth2phi,kmul(J31L,PDstandardNth3phi)));
      
      JacPDstandardNth1trK = 
        kmadd(J11L,PDstandardNth1trK,kmadd(J21L,PDstandardNth2trK,kmul(J31L,PDstandardNth3trK)));
      
      JacPDstandardNth2At11 = 
        kmadd(J12L,PDstandardNth1At11,kmadd(J22L,PDstandardNth2At11,kmul(J32L,PDstandardNth3At11)));
      
      JacPDstandardNth2At12 = 
        kmadd(J12L,PDstandardNth1At12,kmadd(J22L,PDstandardNth2At12,kmul(J32L,PDstandardNth3At12)));
      
      JacPDstandardNth2At13 = 
        kmadd(J12L,PDstandardNth1At13,kmadd(J22L,PDstandardNth2At13,kmul(J32L,PDstandardNth3At13)));
      
      JacPDstandardNth2At22 = 
        kmadd(J12L,PDstandardNth1At22,kmadd(J22L,PDstandardNth2At22,kmul(J32L,PDstandardNth3At22)));
      
      JacPDstandardNth2At23 = 
        kmadd(J12L,PDstandardNth1At23,kmadd(J22L,PDstandardNth2At23,kmul(J32L,PDstandardNth3At23)));
      
      JacPDstandardNth2At33 = 
        kmadd(J12L,PDstandardNth1At33,kmadd(J22L,PDstandardNth2At33,kmul(J32L,PDstandardNth3At33)));
      
      JacPDstandardNth2gt11 = 
        kmadd(J12L,PDstandardNth1gt11,kmadd(J22L,PDstandardNth2gt11,kmul(J32L,PDstandardNth3gt11)));
      
      JacPDstandardNth2gt12 = 
        kmadd(J12L,PDstandardNth1gt12,kmadd(J22L,PDstandardNth2gt12,kmul(J32L,PDstandardNth3gt12)));
      
      JacPDstandardNth2gt13 = 
        kmadd(J12L,PDstandardNth1gt13,kmadd(J22L,PDstandardNth2gt13,kmul(J32L,PDstandardNth3gt13)));
      
      JacPDstandardNth2gt22 = 
        kmadd(J12L,PDstandardNth1gt22,kmadd(J22L,PDstandardNth2gt22,kmul(J32L,PDstandardNth3gt22)));
      
      JacPDstandardNth2gt23 = 
        kmadd(J12L,PDstandardNth1gt23,kmadd(J22L,PDstandardNth2gt23,kmul(J32L,PDstandardNth3gt23)));
      
      JacPDstandardNth2gt33 = 
        kmadd(J12L,PDstandardNth1gt33,kmadd(J22L,PDstandardNth2gt33,kmul(J32L,PDstandardNth3gt33)));
      
      JacPDstandardNth2phi = 
        kmadd(J12L,PDstandardNth1phi,kmadd(J22L,PDstandardNth2phi,kmul(J32L,PDstandardNth3phi)));
      
      JacPDstandardNth2trK = 
        kmadd(J12L,PDstandardNth1trK,kmadd(J22L,PDstandardNth2trK,kmul(J32L,PDstandardNth3trK)));
      
      JacPDstandardNth3At11 = 
        kmadd(J13L,PDstandardNth1At11,kmadd(J23L,PDstandardNth2At11,kmul(J33L,PDstandardNth3At11)));
      
      JacPDstandardNth3At12 = 
        kmadd(J13L,PDstandardNth1At12,kmadd(J23L,PDstandardNth2At12,kmul(J33L,PDstandardNth3At12)));
      
      JacPDstandardNth3At13 = 
        kmadd(J13L,PDstandardNth1At13,kmadd(J23L,PDstandardNth2At13,kmul(J33L,PDstandardNth3At13)));
      
      JacPDstandardNth3At22 = 
        kmadd(J13L,PDstandardNth1At22,kmadd(J23L,PDstandardNth2At22,kmul(J33L,PDstandardNth3At22)));
      
      JacPDstandardNth3At23 = 
        kmadd(J13L,PDstandardNth1At23,kmadd(J23L,PDstandardNth2At23,kmul(J33L,PDstandardNth3At23)));
      
      JacPDstandardNth3At33 = 
        kmadd(J13L,PDstandardNth1At33,kmadd(J23L,PDstandardNth2At33,kmul(J33L,PDstandardNth3At33)));
      
      JacPDstandardNth3gt11 = 
        kmadd(J13L,PDstandardNth1gt11,kmadd(J23L,PDstandardNth2gt11,kmul(J33L,PDstandardNth3gt11)));
      
      JacPDstandardNth3gt12 = 
        kmadd(J13L,PDstandardNth1gt12,kmadd(J23L,PDstandardNth2gt12,kmul(J33L,PDstandardNth3gt12)));
      
      JacPDstandardNth3gt13 = 
        kmadd(J13L,PDstandardNth1gt13,kmadd(J23L,PDstandardNth2gt13,kmul(J33L,PDstandardNth3gt13)));
      
      JacPDstandardNth3gt22 = 
        kmadd(J13L,PDstandardNth1gt22,kmadd(J23L,PDstandardNth2gt22,kmul(J33L,PDstandardNth3gt22)));
      
      JacPDstandardNth3gt23 = 
        kmadd(J13L,PDstandardNth1gt23,kmadd(J23L,PDstandardNth2gt23,kmul(J33L,PDstandardNth3gt23)));
      
      JacPDstandardNth3gt33 = 
        kmadd(J13L,PDstandardNth1gt33,kmadd(J23L,PDstandardNth2gt33,kmul(J33L,PDstandardNth3gt33)));
      
      JacPDstandardNth3phi = 
        kmadd(J13L,PDstandardNth1phi,kmadd(J23L,PDstandardNth2phi,kmul(J33L,PDstandardNth3phi)));
      
      JacPDstandardNth3trK = 
        kmadd(J13L,PDstandardNth1trK,kmadd(J23L,PDstandardNth2trK,kmul(J33L,PDstandardNth3trK)));
    }
    else
    {
      JacPDstandardNth1At11 = PDstandardNth1At11;
      
      JacPDstandardNth1At12 = PDstandardNth1At12;
      
      JacPDstandardNth1At13 = PDstandardNth1At13;
      
      JacPDstandardNth1At22 = PDstandardNth1At22;
      
      JacPDstandardNth1At23 = PDstandardNth1At23;
      
      JacPDstandardNth1At33 = PDstandardNth1At33;
      
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth1phi = PDstandardNth1phi;
      
      JacPDstandardNth1trK = PDstandardNth1trK;
      
      JacPDstandardNth2At11 = PDstandardNth2At11;
      
      JacPDstandardNth2At12 = PDstandardNth2At12;
      
      JacPDstandardNth2At13 = PDstandardNth2At13;
      
      JacPDstandardNth2At22 = PDstandardNth2At22;
      
      JacPDstandardNth2At23 = PDstandardNth2At23;
      
      JacPDstandardNth2At33 = PDstandardNth2At33;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth2phi = PDstandardNth2phi;
      
      JacPDstandardNth2trK = PDstandardNth2trK;
      
      JacPDstandardNth3At11 = PDstandardNth3At11;
      
      JacPDstandardNth3At12 = PDstandardNth3At12;
      
      JacPDstandardNth3At13 = PDstandardNth3At13;
      
      JacPDstandardNth3At22 = PDstandardNth3At22;
      
      JacPDstandardNth3At23 = PDstandardNth3At23;
      
      JacPDstandardNth3At33 = PDstandardNth3At33;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDstandardNth3phi = PDstandardNth3phi;
      
      JacPDstandardNth3trK = PDstandardNth3trK;
    }
    
    CCTK_REAL_VEC detgt = ToReal(1);
    
    CCTK_REAL_VEC gtu11 = kmul(INV(detgt),kmsub(gt22L,gt33L,SQR(gt23L)));
    
    CCTK_REAL_VEC gtu12 = 
      kmul(INV(detgt),kmsub(gt13L,gt23L,kmul(gt12L,gt33L)));
    
    CCTK_REAL_VEC gtu13 = 
      kmul(INV(detgt),kmsub(gt12L,gt23L,kmul(gt13L,gt22L)));
    
    CCTK_REAL_VEC gtu22 = kmul(INV(detgt),kmsub(gt11L,gt33L,SQR(gt13L)));
    
    CCTK_REAL_VEC gtu23 = 
      kmul(INV(detgt),kmsub(gt12L,gt13L,kmul(gt11L,gt23L)));
    
    CCTK_REAL_VEC gtu33 = kmul(INV(detgt),kmsub(gt11L,gt22L,SQR(gt12L)));
    
    CCTK_REAL_VEC Gtl111 = kmul(JacPDstandardNth1gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl112 = kmul(JacPDstandardNth2gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl113 = kmul(JacPDstandardNth3gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl122 = 
      kmadd(JacPDstandardNth1gt22,ToReal(-0.5),JacPDstandardNth2gt12);
    
    CCTK_REAL_VEC Gtl123 = 
      kmul(kadd(JacPDstandardNth2gt13,ksub(JacPDstandardNth3gt12,JacPDstandardNth1gt23)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl133 = 
      kmadd(JacPDstandardNth1gt33,ToReal(-0.5),JacPDstandardNth3gt13);
    
    CCTK_REAL_VEC Gtl211 = 
      kmadd(JacPDstandardNth2gt11,ToReal(-0.5),JacPDstandardNth1gt12);
    
    CCTK_REAL_VEC Gtl212 = kmul(JacPDstandardNth1gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl213 = 
      kmul(kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth3gt12,JacPDstandardNth2gt13)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl222 = kmul(JacPDstandardNth2gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl223 = kmul(JacPDstandardNth3gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl233 = 
      kmadd(JacPDstandardNth2gt33,ToReal(-0.5),JacPDstandardNth3gt23);
    
    CCTK_REAL_VEC Gtl311 = 
      kmadd(JacPDstandardNth3gt11,ToReal(-0.5),JacPDstandardNth1gt13);
    
    CCTK_REAL_VEC Gtl312 = 
      kmul(kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth2gt13,JacPDstandardNth3gt12)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl313 = kmul(JacPDstandardNth1gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl322 = 
      kmadd(JacPDstandardNth3gt22,ToReal(-0.5),JacPDstandardNth2gt23);
    
    CCTK_REAL_VEC Gtl323 = kmul(JacPDstandardNth2gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl333 = kmul(JacPDstandardNth3gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gt111 = 
      kmadd(Gtl111,gtu11,kmadd(Gtl211,gtu12,kmul(Gtl311,gtu13)));
    
    CCTK_REAL_VEC Gt211 = 
      kmadd(Gtl111,gtu12,kmadd(Gtl211,gtu22,kmul(Gtl311,gtu23)));
    
    CCTK_REAL_VEC Gt311 = 
      kmadd(Gtl111,gtu13,kmadd(Gtl211,gtu23,kmul(Gtl311,gtu33)));
    
    CCTK_REAL_VEC Gt112 = 
      kmadd(Gtl112,gtu11,kmadd(Gtl212,gtu12,kmul(Gtl312,gtu13)));
    
    CCTK_REAL_VEC Gt212 = 
      kmadd(Gtl112,gtu12,kmadd(Gtl212,gtu22,kmul(Gtl312,gtu23)));
    
    CCTK_REAL_VEC Gt312 = 
      kmadd(Gtl112,gtu13,kmadd(Gtl212,gtu23,kmul(Gtl312,gtu33)));
    
    CCTK_REAL_VEC Gt113 = 
      kmadd(Gtl113,gtu11,kmadd(Gtl213,gtu12,kmul(Gtl313,gtu13)));
    
    CCTK_REAL_VEC Gt213 = 
      kmadd(Gtl113,gtu12,kmadd(Gtl213,gtu22,kmul(Gtl313,gtu23)));
    
    CCTK_REAL_VEC Gt313 = 
      kmadd(Gtl113,gtu13,kmadd(Gtl213,gtu23,kmul(Gtl313,gtu33)));
    
    CCTK_REAL_VEC Gt122 = 
      kmadd(Gtl122,gtu11,kmadd(Gtl222,gtu12,kmul(Gtl322,gtu13)));
    
    CCTK_REAL_VEC Gt222 = 
      kmadd(Gtl122,gtu12,kmadd(Gtl222,gtu22,kmul(Gtl322,gtu23)));
    
    CCTK_REAL_VEC Gt322 = 
      kmadd(Gtl122,gtu13,kmadd(Gtl222,gtu23,kmul(Gtl322,gtu33)));
    
    CCTK_REAL_VEC Gt123 = 
      kmadd(Gtl123,gtu11,kmadd(Gtl223,gtu12,kmul(Gtl323,gtu13)));
    
    CCTK_REAL_VEC Gt223 = 
      kmadd(Gtl123,gtu12,kmadd(Gtl223,gtu22,kmul(Gtl323,gtu23)));
    
    CCTK_REAL_VEC Gt323 = 
      kmadd(Gtl123,gtu13,kmadd(Gtl223,gtu23,kmul(Gtl323,gtu33)));
    
    CCTK_REAL_VEC Gt133 = 
      kmadd(Gtl133,gtu11,kmadd(Gtl233,gtu12,kmul(Gtl333,gtu13)));
    
    CCTK_REAL_VEC Gt233 = 
      kmadd(Gtl133,gtu12,kmadd(Gtl233,gtu22,kmul(Gtl333,gtu23)));
    
    CCTK_REAL_VEC Gt333 = 
      kmadd(Gtl133,gtu13,kmadd(Gtl233,gtu23,kmul(Gtl333,gtu33)));
    
    CCTK_REAL_VEC fac1 = 
      IfThen(conformalMethod,kmul(INV(phiL),ToReal(-0.5)),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 = kmul(fac1,JacPDstandardNth1phi);
    
    CCTK_REAL_VEC cdphi2 = kmul(fac1,JacPDstandardNth2phi);
    
    CCTK_REAL_VEC cdphi3 = kmul(fac1,JacPDstandardNth3phi);
    
    CCTK_REAL_VEC S1 = 
      kmul(INV(alphaL),kmadd(beta1L,eTxxL,kmadd(beta2L,eTxyL,kmsub(beta3L,eTxzL,eTtxL))));
    
    CCTK_REAL_VEC S2 = 
      kmul(INV(alphaL),kmadd(beta1L,eTxyL,kmadd(beta2L,eTyyL,kmsub(beta3L,eTyzL,eTtyL))));
    
    CCTK_REAL_VEC S3 = 
      kmul(INV(alphaL),kmadd(beta1L,eTxzL,kmadd(beta2L,eTyzL,kmsub(beta3L,eTzzL,eTtzL))));
    
    CCTK_REAL_VEC M1L = 
      kmadd(gtu11,JacPDstandardNth1At11,kmadd(gtu22,JacPDstandardNth2At12,kmadd(gtu23,kadd(JacPDstandardNth2At13,JacPDstandardNth3At12),kmadd(gtu33,JacPDstandardNth3At13,kmadd(S1,ToReal(-25.13274122871834590770114706623602307358),kmadd(kmadd(kmadd(At22L,Gt212,kmul(At23L,Gt312)),gtu22,kmadd(kmadd(At13L,Gt112,kmadd(At22L,Gt213,kmadd(At33L,Gt312,kmul(At23L,kadd(Gt212,Gt313))))),gtu23,kmul(kmadd(At13L,Gt113,kmadd(At23L,Gt213,kmul(At33L,Gt313))),gtu33))),ToReal(-1.),kmadd(gtu12,kadd(JacPDstandardNth1At12,kadd(JacPDstandardNth2At11,kmadd(At13L,kmul(Gt312,ToReal(-3.)),kmul(At22L,kmul(Gt211,ToReal(-1.)))))),kmadd(gtu13,kadd(JacPDstandardNth1At13,kadd(JacPDstandardNth3At11,kmadd(At13L,kmul(Gt313,ToReal(-3.)),kmul(At23L,kmul(Gt211,ToReal(-1.)))))),kmadd(Gt311,kmadd(At13L,kmul(gtu11,ToReal(-2.)),kmul(kmadd(At23L,gtu12,kmul(At33L,gtu13)),ToReal(-1.))),kmadd(JacPDstandardNth1trK,ToReal(-0.6666666666666666666666666666666666666667),kmadd(At13L,kmadd(kmadd(Gt322,gtu22,kmul(Gt333,gtu33)),ToReal(-1.),kmadd(cdphi3,kmul(gtu33,ToReal(6.)),kmadd(gtu13,kmadd(Gt111,ToReal(-1.),kmul(cdphi1,ToReal(6.))),kmul(gtu23,kmadd(Gt323,ToReal(-2.),kmul(cdphi2,ToReal(6.))))))),kmadd(At11L,kmadd(Gt123,kmul(gtu23,ToReal(-2.)),kmadd(kmadd(Gt122,gtu22,kmul(Gt133,gtu33)),ToReal(-1.),kmadd(gtu11,kmadd(Gt111,ToReal(-2.),kmul(cdphi1,ToReal(6.))),kmadd(gtu12,kmadd(Gt112,ToReal(-3.),kmul(cdphi2,ToReal(6.))),kmul(gtu13,kmadd(Gt113,ToReal(-3.),kmul(cdphi3,ToReal(6.)))))))),kmul(At12L,kmadd(Gt213,kmul(gtu13,ToReal(-3.)),kmadd(kmadd(Gt211,gtu11,kmul(Gt223,gtu23)),ToReal(-2.),kmadd(Gt233,kmul(gtu33,ToReal(-1.)),kmadd(gtu12,kmadd(Gt212,ToReal(-3.),kmadd(Gt111,ToReal(-1.),kmul(cdphi1,ToReal(6.)))),kmadd(gtu22,kmadd(kadd(Gt112,Gt222),ToReal(-1.),kmul(cdphi2,ToReal(6.))),kmul(gtu23,kmadd(Gt113,ToReal(-1.),kmul(cdphi3,ToReal(6.))))))))))))))))))))));
    
    CCTK_REAL_VEC M2L = 
      kmadd(gtu11,JacPDstandardNth1At12,kmadd(gtu12,kadd(JacPDstandardNth1At22,JacPDstandardNth2At12),kmadd(gtu22,JacPDstandardNth2At22,kmadd(gtu33,JacPDstandardNth3At23,kmadd(S2,ToReal(-25.13274122871834590770114706623602307358),kmadd(kmadd(kmadd(At22L,Gt211,kmadd(At23L,Gt311,kmul(At13L,Gt312))),gtu11,kmadd(kmadd(At23L,Gt212,kmul(At33L,Gt312)),gtu13,kmadd(At11L,kmadd(Gt112,gtu11,kmadd(Gt122,gtu12,kmul(Gt123,gtu13))),kmadd(kmadd(At23L,Gt223,kmul(At33L,Gt323)),gtu33,kmul(At13L,kmadd(Gt322,gtu12,kmadd(Gt112,gtu13,kmadd(Gt122,gtu23,kmul(Gt123,gtu33))))))))),ToReal(-1.),kmadd(gtu23,kadd(JacPDstandardNth2At23,kadd(JacPDstandardNth3At22,kmadd(kmadd(At22L,Gt223,kmul(At23L,Gt323)),ToReal(-3.),kmul(kmadd(At23L,Gt222,kmul(At33L,Gt322)),ToReal(-1.))))),kmadd(gtu13,kadd(JacPDstandardNth1At23,kadd(JacPDstandardNth3At12,kmadd(At23L,kmul(Gt313,ToReal(-2.)),kmul(At13L,kmul(Gt323,ToReal(-1.)))))),kmadd(JacPDstandardNth2trK,ToReal(-0.6666666666666666666666666666666666666667),kmadd(At23L,kmadd(Gt312,kmul(gtu12,ToReal(-3.)),kmadd(Gt322,kmul(gtu22,ToReal(-2.)),kmadd(Gt333,kmul(gtu33,ToReal(-1.)),kmul(kmadd(cdphi1,gtu13,kmadd(cdphi2,gtu23,kmul(cdphi3,gtu33))),ToReal(6.))))),kmadd(At22L,kmadd(kmadd(Gt213,gtu13,kmul(Gt222,gtu22)),ToReal(-2.),kmadd(Gt233,kmul(gtu33,ToReal(-1.)),kmadd(kmadd(cdphi2,gtu22,kmul(cdphi3,gtu23)),ToReal(6.),kmul(gtu12,kmadd(Gt212,ToReal(-3.),kmul(cdphi1,ToReal(6.))))))),kmul(At12L,kmadd(Gt123,kmul(gtu23,ToReal(-3.)),kmadd(Gt122,kmul(gtu22,ToReal(-2.)),kmadd(Gt133,kmul(gtu33,ToReal(-1.)),kmadd(gtu11,kmadd(kadd(Gt111,Gt212),ToReal(-1.),kmul(cdphi1,ToReal(6.))),kmadd(gtu12,kmadd(Gt112,ToReal(-3.),kmadd(Gt222,ToReal(-1.),kmul(cdphi2,ToReal(6.)))),kmul(gtu13,kmadd(Gt113,ToReal(-2.),kmadd(Gt223,ToReal(-1.),kmul(cdphi3,ToReal(6.))))))))))))))))))))));
    
    CCTK_REAL_VEC M3L = 
      kmadd(gtu11,JacPDstandardNth1At13,kmadd(gtu22,JacPDstandardNth2At23,kmadd(gtu13,kadd(JacPDstandardNth1At33,JacPDstandardNth3At13),kmadd(gtu33,JacPDstandardNth3At33,kmadd(S3,ToReal(-25.13274122871834590770114706623602307358),kmadd(kmadd(kmadd(At23L,Gt211,kmadd(At12L,Gt213,kmul(At33L,Gt311))),gtu11,kmadd(kmadd(At22L,Gt213,kmul(At12L,kadd(Gt113,Gt223))),gtu12,kmadd(At11L,kmadd(Gt113,gtu11,kmadd(Gt123,gtu12,kmul(Gt133,gtu13))),kmadd(kmadd(At23L,Gt222,kmul(At22L,Gt223)),gtu22,kmul(At12L,kmadd(Gt233,gtu13,kmadd(Gt123,gtu22,kmul(Gt133,gtu23)))))))),ToReal(-1.),kmadd(gtu12,kadd(JacPDstandardNth1At23,kadd(JacPDstandardNth2At13,kmadd(At33L,kmul(Gt312,ToReal(-2.)),kmul(At23L,kmul(Gt313,ToReal(-1.)))))),kmadd(gtu23,kadd(JacPDstandardNth2At33,kadd(JacPDstandardNth3At23,kmadd(kmadd(At23L,Gt223,kmul(At33L,Gt323)),ToReal(-3.),kmul(kmadd(At22L,Gt233,kmul(At23L,Gt333)),ToReal(-1.))))),kmadd(JacPDstandardNth3trK,ToReal(-0.6666666666666666666666666666666666666667),kmadd(At33L,kmadd(Gt333,kmul(gtu33,ToReal(-2.)),kmadd(Gt322,kmul(gtu22,ToReal(-1.)),kmadd(kmadd(cdphi2,gtu23,kmul(cdphi3,gtu33)),ToReal(6.),kmul(gtu13,kmadd(Gt313,ToReal(-3.),kmul(cdphi1,ToReal(6.))))))),kmadd(At23L,kmadd(Gt213,kmul(gtu13,ToReal(-3.)),kmadd(Gt233,kmul(gtu33,ToReal(-2.)),kmadd(cdphi3,kmul(gtu23,ToReal(6.)),kmadd(gtu12,kmadd(Gt212,ToReal(-2.),kmul(cdphi1,ToReal(6.))),kmul(gtu22,kmadd(Gt323,ToReal(-1.),kmul(cdphi2,ToReal(6.)))))))),kmul(At13L,kmadd(Gt123,kmul(gtu23,ToReal(-3.)),kmadd(Gt133,kmul(gtu33,ToReal(-2.)),kmadd(Gt122,kmul(gtu22,ToReal(-1.)),kmadd(gtu11,kmadd(kadd(Gt111,Gt313),ToReal(-1.),kmul(cdphi1,ToReal(6.))),kmadd(gtu12,kmadd(Gt112,ToReal(-2.),kmadd(Gt323,ToReal(-1.),kmul(cdphi2,ToReal(6.)))),kmul(gtu13,kmadd(Gt113,ToReal(-3.),kmadd(Gt333,ToReal(-1.),kmul(cdphi3,ToReal(6.))))))))))))))))))))));
    
    CCTK_REAL_VEC cSL = klog(detgt);
    
    CCTK_REAL_VEC cXt1L = 
      kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmadd(Gt133,gtu33,kmsub(kmadd(Gt112,gtu12,kmadd(Gt113,gtu13,kmul(Gt123,gtu23))),ToReal(2),Xt1L))));
    
    CCTK_REAL_VEC cXt2L = 
      kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmadd(Gt233,gtu33,kmsub(kmadd(Gt212,gtu12,kmadd(Gt213,gtu13,kmul(Gt223,gtu23))),ToReal(2),Xt2L))));
    
    CCTK_REAL_VEC cXt3L = 
      kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmadd(Gt333,gtu33,kmsub(kmadd(Gt312,gtu12,kmadd(Gt313,gtu13,kmul(Gt323,gtu23))),ToReal(2),Xt3L))));
    
    CCTK_REAL_VEC cAL = 
      kmadd(At11L,gtu11,kmadd(At22L,gtu22,kmadd(At33L,gtu33,kmul(kmadd(At12L,gtu12,kmadd(At13L,gtu13,kmul(At23L,gtu23))),ToReal(2)))));
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 2 && CCTK_BUILTIN_EXPECT(i < lc_imin && i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count_lo = lc_imin-i;
      ptrdiff_t const elt_count_hi = lc_imax-i;
      vec_store_nta_partial_mid(cA[index],cAL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(cS[index],cSL,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(cXt1[index],cXt1L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(cXt2[index],cXt2L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(cXt3[index],cXt3L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(M1[index],M1L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(M2[index],M2L,elt_count_lo,elt_count_hi);
      vec_store_nta_partial_mid(M3[index],M3L,elt_count_lo,elt_count_hi);
      break;
    }
    
    /* If necessary, store only partial vectors after the first iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i < lc_imin, 0))
    {
      ptrdiff_t const elt_count = lc_imin-i;
      vec_store_nta_partial_hi(cA[index],cAL,elt_count);
      vec_store_nta_partial_hi(cS[index],cSL,elt_count);
      vec_store_nta_partial_hi(cXt1[index],cXt1L,elt_count);
      vec_store_nta_partial_hi(cXt2[index],cXt2L,elt_count);
      vec_store_nta_partial_hi(cXt3[index],cXt3L,elt_count);
      vec_store_nta_partial_hi(M1[index],M1L,elt_count);
      vec_store_nta_partial_hi(M2[index],M2L,elt_count);
      vec_store_nta_partial_hi(M3[index],M3L,elt_count);
      continue;
    }
    
    /* If necessary, store only partial vectors after the last iteration */
    
    if (CCTK_REAL_VEC_SIZE > 1 && CCTK_BUILTIN_EXPECT(i+CCTK_REAL_VEC_SIZE > lc_imax, 0))
    {
      ptrdiff_t const elt_count = lc_imax-i;
      vec_store_nta_partial_lo(cA[index],cAL,elt_count);
      vec_store_nta_partial_lo(cS[index],cSL,elt_count);
      vec_store_nta_partial_lo(cXt1[index],cXt1L,elt_count);
      vec_store_nta_partial_lo(cXt2[index],cXt2L,elt_count);
      vec_store_nta_partial_lo(cXt3[index],cXt3L,elt_count);
      vec_store_nta_partial_lo(M1[index],M1L,elt_count);
      vec_store_nta_partial_lo(M2[index],M2L,elt_count);
      vec_store_nta_partial_lo(M3[index],M3L,elt_count);
      break;
    }
    
    /* Copy local copies back to grid functions */
    vec_store_nta(cA[index],cAL);
    vec_store_nta(cS[index],cSL);
    vec_store_nta(cXt1[index],cXt1L);
    vec_store_nta(cXt2[index],cXt2L);
    vec_store_nta(cXt3[index],cXt3L);
    vec_store_nta(M1[index],M1L);
    vec_store_nta(M2[index],M2L);
    vec_store_nta(M3[index],M3L);
  }
  LC_ENDLOOP3VEC (ML_BSSN_MP_O8_constraints2);
}

extern "C" void ML_BSSN_MP_O8_constraints2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverInterior(cctkGH, &ML_BSSN_MP_O8_constraints2_Body);
}
