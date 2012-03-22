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

static void ML_BSSN_Host_convertFromADMBase_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  LC_LOOP3VEC(ML_BSSN_Host_convertFromADMBase,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alpL = vec_load(alp[index]);
    CCTK_REAL_VEC betaxL = vec_load(betax[index]);
    CCTK_REAL_VEC betayL = vec_load(betay[index]);
    CCTK_REAL_VEC betazL = vec_load(betaz[index]);
    CCTK_REAL_VEC gxxL = vec_load(gxx[index]);
    CCTK_REAL_VEC gxyL = vec_load(gxy[index]);
    CCTK_REAL_VEC gxzL = vec_load(gxz[index]);
    CCTK_REAL_VEC gyyL = vec_load(gyy[index]);
    CCTK_REAL_VEC gyzL = vec_load(gyz[index]);
    CCTK_REAL_VEC gzzL = vec_load(gzz[index]);
    CCTK_REAL_VEC kxxL = vec_load(kxx[index]);
    CCTK_REAL_VEC kxyL = vec_load(kxy[index]);
    CCTK_REAL_VEC kxzL = vec_load(kxz[index]);
    CCTK_REAL_VEC kyyL = vec_load(kyy[index]);
    CCTK_REAL_VEC kyzL = vec_load(kyz[index]);
    CCTK_REAL_VEC kzzL = vec_load(kzz[index]);
    CCTK_REAL_VEC phiL = vec_load(phi[index]);
    CCTK_REAL_VEC rL = vec_load(r[index]);
    CCTK_REAL_VEC trKL = vec_load(trK[index]);
    CCTK_REAL_VEC xL = vec_load(x[index]);
    CCTK_REAL_VEC yL = vec_load(y[index]);
    CCTK_REAL_VEC zL = vec_load(z[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC g11 = gxxL;
    
    CCTK_REAL_VEC g12 = gxyL;
    
    CCTK_REAL_VEC g13 = gxzL;
    
    CCTK_REAL_VEC g22 = gyyL;
    
    CCTK_REAL_VEC g23 = gyzL;
    
    CCTK_REAL_VEC g33 = gzzL;
    
    CCTK_REAL_VEC detg = 
      knmsub(g22,SQR(g13),knmsub(g11,SQR(g23),kmadd(g33,kmsub(g11,g22,SQR(g12)),kmul(g12,kmul(g13,kmul(g23,ToReal(2)))))));
    
    CCTK_REAL_VEC gu11 = kmul(INV(detg),kmsub(g22,g33,SQR(g23)));
    
    CCTK_REAL_VEC gu12 = kmul(INV(detg),kmsub(g13,g23,kmul(g12,g33)));
    
    CCTK_REAL_VEC gu13 = kmul(INV(detg),kmsub(g12,g23,kmul(g13,g22)));
    
    CCTK_REAL_VEC gu22 = kmul(INV(detg),kmsub(g11,g33,SQR(g13)));
    
    CCTK_REAL_VEC gu23 = kmul(INV(detg),kmsub(g12,g13,kmul(g11,g23)));
    
    CCTK_REAL_VEC gu33 = kmul(INV(detg),kmsub(g11,g22,SQR(g12)));
    
    CCTK_REAL_VEC em4phi;
    
    if (conformalMethod)
    {
      phiL = kpow(detg,-0.166666666666666666666666666667);
      
      em4phi = SQR(phiL);
    }
    else
    {
      phiL = kmul(klog(detg),ToReal(0.0833333333333333333333333333333));
      
      em4phi = kexp(kmul(phiL,ToReal(-4)));
    }
    
    CCTK_REAL_VEC gt11L = kmul(em4phi,g11);
    
    CCTK_REAL_VEC gt12L = kmul(em4phi,g12);
    
    CCTK_REAL_VEC gt13L = kmul(em4phi,g13);
    
    CCTK_REAL_VEC gt22L = kmul(em4phi,g22);
    
    CCTK_REAL_VEC gt23L = kmul(em4phi,g23);
    
    CCTK_REAL_VEC gt33L = kmul(em4phi,g33);
    
    trKL = 
      kmadd(kxxL,gu11,kmadd(kyyL,gu22,kmadd(kzzL,gu33,kmul(kmadd(kxyL,gu12,kmadd(kxzL,gu13,kmul(kyzL,gu23))),ToReal(2)))));
    
    CCTK_REAL_VEC At11L = 
      kmul(em4phi,kmadd(trKL,kmul(g11,ToReal(-0.333333333333333333333333333333)),kxxL));
    
    CCTK_REAL_VEC At12L = 
      kmul(em4phi,kmadd(trKL,kmul(g12,ToReal(-0.333333333333333333333333333333)),kxyL));
    
    CCTK_REAL_VEC At13L = 
      kmul(em4phi,kmadd(trKL,kmul(g13,ToReal(-0.333333333333333333333333333333)),kxzL));
    
    CCTK_REAL_VEC At22L = 
      kmul(em4phi,kmadd(trKL,kmul(g22,ToReal(-0.333333333333333333333333333333)),kyyL));
    
    CCTK_REAL_VEC At23L = 
      kmul(em4phi,kmadd(trKL,kmul(g23,ToReal(-0.333333333333333333333333333333)),kyzL));
    
    CCTK_REAL_VEC At33L = 
      kmul(em4phi,kmadd(trKL,kmul(g33,ToReal(-0.333333333333333333333333333333)),kzzL));
    
    CCTK_REAL_VEC alphaL = alpL;
    
    CCTK_REAL_VEC beta1L = betaxL;
    
    CCTK_REAL_VEC beta2L = betayL;
    
    CCTK_REAL_VEC beta3L = betazL;
    
    CCTK_REAL_VEC xCopyL = xL;
    
    CCTK_REAL_VEC yCopyL = yL;
    
    CCTK_REAL_VEC zCopyL = zL;
    
    CCTK_REAL_VEC rCopyL = rL;
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,lc_imin,lc_imax);
    vec_store_nta_partial(alpha[index],alphaL);
    vec_store_nta_partial(At11[index],At11L);
    vec_store_nta_partial(At12[index],At12L);
    vec_store_nta_partial(At13[index],At13L);
    vec_store_nta_partial(At22[index],At22L);
    vec_store_nta_partial(At23[index],At23L);
    vec_store_nta_partial(At33[index],At33L);
    vec_store_nta_partial(beta1[index],beta1L);
    vec_store_nta_partial(beta2[index],beta2L);
    vec_store_nta_partial(beta3[index],beta3L);
    vec_store_nta_partial(gt11[index],gt11L);
    vec_store_nta_partial(gt12[index],gt12L);
    vec_store_nta_partial(gt13[index],gt13L);
    vec_store_nta_partial(gt22[index],gt22L);
    vec_store_nta_partial(gt23[index],gt23L);
    vec_store_nta_partial(gt33[index],gt33L);
    vec_store_nta_partial(phi[index],phiL);
    vec_store_nta_partial(rCopy[index],rCopyL);
    vec_store_nta_partial(trK[index],trKL);
    vec_store_nta_partial(xCopy[index],xCopyL);
    vec_store_nta_partial(yCopy[index],yCopyL);
    vec_store_nta_partial(zCopy[index],zCopyL);
  }
  LC_ENDLOOP3VEC(ML_BSSN_Host_convertFromADMBase);
}

extern "C" void ML_BSSN_Host_convertFromADMBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_Host_convertFromADMBase_Body");
  }
  
  if (cctk_iteration % ML_BSSN_Host_convertFromADMBase_calc_every != ML_BSSN_Host_convertFromADMBase_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "ADMBase::curv",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "ML_BSSN_Host::coords",
    "grid::coordinates",
    "Grid::coordinates",
    "ML_BSSN_Host::ML_curv",
    "ML_BSSN_Host::ML_lapse",
    "ML_BSSN_Host::ML_log_confac",
    "ML_BSSN_Host::ML_metric",
    "ML_BSSN_Host::ML_shift",
    "ML_BSSN_Host::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_convertFromADMBase", 13, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, ML_BSSN_Host_convertFromADMBase_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_convertFromADMBase_Body");
  }
}
