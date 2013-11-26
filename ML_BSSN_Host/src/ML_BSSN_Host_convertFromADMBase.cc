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

static void ML_BSSN_Host_convertFromADMBase_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const CCTK_REAL_VEC dx CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(0));
  const CCTK_REAL_VEC dy CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(1));
  const CCTK_REAL_VEC dz CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(2));
  const CCTK_REAL_VEC dt CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_TIME);
  const CCTK_REAL_VEC t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);
  const CCTK_REAL_VEC dxi CCTK_ATTRIBUTE_UNUSED = INV(dx);
  const CCTK_REAL_VEC dyi CCTK_ATTRIBUTE_UNUSED = INV(dy);
  const CCTK_REAL_VEC dzi CCTK_ATTRIBUTE_UNUSED = INV(dz);
  const CCTK_REAL_VEC khalf CCTK_ATTRIBUTE_UNUSED = ToReal(0.5);
  const CCTK_REAL_VEC kthird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(0.333333333333333333333333333333);
  const CCTK_REAL_VEC ktwothird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(0.666666666666666666666666666667);
  const CCTK_REAL_VEC kfourthird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(1.33333333333333333333333333333);
  const CCTK_REAL_VEC hdxi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dxi,ToReal(0.5));
  const CCTK_REAL_VEC hdyi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dyi,ToReal(0.5));
  const CCTK_REAL_VEC hdzi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dzi,ToReal(0.5));
  
  /* Initialize predefined quantities */
  const CCTK_REAL_VEC p1o1024dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dx);
  const CCTK_REAL_VEC p1o1024dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dy);
  const CCTK_REAL_VEC p1o1024dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0009765625),dz);
  const CCTK_REAL_VEC p1o1680dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dx);
  const CCTK_REAL_VEC p1o1680dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dy);
  const CCTK_REAL_VEC p1o1680dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000595238095238095238095238095238),dz);
  const CCTK_REAL_VEC p1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dx);
  const CCTK_REAL_VEC p1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dy);
  const CCTK_REAL_VEC p1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dz);
  const CCTK_REAL_VEC p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dx,dx));
  const CCTK_REAL_VEC p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dy,dy));
  const CCTK_REAL_VEC p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.000198412698412698412698412698413),kmul(dz,dz));
  const CCTK_REAL_VEC p1o560dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dx);
  const CCTK_REAL_VEC p1o560dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dy);
  const CCTK_REAL_VEC p1o560dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00178571428571428571428571428571),dz);
  const CCTK_REAL_VEC p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dy,dx));
  const CCTK_REAL_VEC p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dz,dx));
  const CCTK_REAL_VEC p1o705600dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1.41723356009070294784580498866e-6),kmul(dz,dy));
  const CCTK_REAL_VEC p1o840dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dx);
  const CCTK_REAL_VEC p1o840dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dy);
  const CCTK_REAL_VEC p1o840dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00119047619047619047619047619048),dz);
  const CCTK_REAL_VEC p1odx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dx);
  const CCTK_REAL_VEC p1ody CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dy);
  const CCTK_REAL_VEC p1odz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dz);
  const CCTK_REAL_VEC pm1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dx);
  const CCTK_REAL_VEC pm1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dy);
  const CCTK_REAL_VEC pm1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dz);
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  #pragma omp parallel // reduction(+: vec_iter_counter, vec_op_counter, vec_mem_counter)
  CCTK_LOOP3STR(ML_BSSN_Host_convertFromADMBase,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    // vec_iter_counter+=CCTK_REAL_VEC_SIZE;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alpL CCTK_ATTRIBUTE_UNUSED = vec_load(alp[index]);
    CCTK_REAL_VEC betaxL CCTK_ATTRIBUTE_UNUSED = vec_load(betax[index]);
    CCTK_REAL_VEC betayL CCTK_ATTRIBUTE_UNUSED = vec_load(betay[index]);
    CCTK_REAL_VEC betazL CCTK_ATTRIBUTE_UNUSED = vec_load(betaz[index]);
    CCTK_REAL_VEC gxxL CCTK_ATTRIBUTE_UNUSED = vec_load(gxx[index]);
    CCTK_REAL_VEC gxyL CCTK_ATTRIBUTE_UNUSED = vec_load(gxy[index]);
    CCTK_REAL_VEC gxzL CCTK_ATTRIBUTE_UNUSED = vec_load(gxz[index]);
    CCTK_REAL_VEC gyyL CCTK_ATTRIBUTE_UNUSED = vec_load(gyy[index]);
    CCTK_REAL_VEC gyzL CCTK_ATTRIBUTE_UNUSED = vec_load(gyz[index]);
    CCTK_REAL_VEC gzzL CCTK_ATTRIBUTE_UNUSED = vec_load(gzz[index]);
    CCTK_REAL_VEC kxxL CCTK_ATTRIBUTE_UNUSED = vec_load(kxx[index]);
    CCTK_REAL_VEC kxyL CCTK_ATTRIBUTE_UNUSED = vec_load(kxy[index]);
    CCTK_REAL_VEC kxzL CCTK_ATTRIBUTE_UNUSED = vec_load(kxz[index]);
    CCTK_REAL_VEC kyyL CCTK_ATTRIBUTE_UNUSED = vec_load(kyy[index]);
    CCTK_REAL_VEC kyzL CCTK_ATTRIBUTE_UNUSED = vec_load(kyz[index]);
    CCTK_REAL_VEC kzzL CCTK_ATTRIBUTE_UNUSED = vec_load(kzz[index]);
    CCTK_REAL_VEC phiL CCTK_ATTRIBUTE_UNUSED = vec_load(phi[index]);
    CCTK_REAL_VEC rL CCTK_ATTRIBUTE_UNUSED = vec_load(r[index]);
    CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = vec_load(trK[index]);
    CCTK_REAL_VEC xL CCTK_ATTRIBUTE_UNUSED = vec_load(x[index]);
    CCTK_REAL_VEC yL CCTK_ATTRIBUTE_UNUSED = vec_load(y[index]);
    CCTK_REAL_VEC zL CCTK_ATTRIBUTE_UNUSED = vec_load(z[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC g11 CCTK_ATTRIBUTE_UNUSED = gxxL;
    
    CCTK_REAL_VEC g12 CCTK_ATTRIBUTE_UNUSED = gxyL;
    
    CCTK_REAL_VEC g13 CCTK_ATTRIBUTE_UNUSED = gxzL;
    
    CCTK_REAL_VEC g22 CCTK_ATTRIBUTE_UNUSED = gyyL;
    
    CCTK_REAL_VEC g23 CCTK_ATTRIBUTE_UNUSED = gyzL;
    
    CCTK_REAL_VEC g33 CCTK_ATTRIBUTE_UNUSED = gzzL;
    
    CCTK_REAL_VEC detg CCTK_ATTRIBUTE_UNUSED = 
      knmsub(g22,kmul(g13,g13),knmsub(g11,kmul(g23,g23),kmadd(g33,kmsub(g11,g22,kmul(g12,g12)),kmul(g12,kmul(g13,kmul(g23,ToReal(2)))))));
    
    CCTK_REAL_VEC gu11 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g22,g33,kmul(g23,g23)),detg);
    
    CCTK_REAL_VEC gu12 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g13,g23,kmul(g12,g33)),detg);
    
    CCTK_REAL_VEC gu13 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g12,g23,kmul(g13,g22)),detg);
    
    CCTK_REAL_VEC gu22 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g11,g33,kmul(g13,g13)),detg);
    
    CCTK_REAL_VEC gu23 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g12,g13,kmul(g11,g23)),detg);
    
    CCTK_REAL_VEC gu33 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(g11,g22,kmul(g12,g12)),detg);
    
    CCTK_REAL_VEC em4phi CCTK_ATTRIBUTE_UNUSED;
    
    if (conformalMethod)
    {
      phiL = kpow(detg,-0.166666666666666666666666666667);
      
      em4phi = kmul(phiL,phiL);
    }
    else
    {
      phiL = kmul(klog(detg),ToReal(0.0833333333333333333333333333333));
      
      em4phi = kexp(kmul(phiL,ToReal(-4)));
    }
    
    CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,g11);
    
    CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,g12);
    
    CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,g13);
    
    CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,g22);
    
    CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,g23);
    
    CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,g33);
    
    trKL = 
      kmadd(kxxL,gu11,kmadd(kyyL,gu22,kmadd(kzzL,gu33,kmul(kmadd(kxyL,gu12,kmadd(kxzL,gu13,kmul(kyzL,gu23))),ToReal(2)))));
    
    CCTK_REAL_VEC At11L CCTK_ATTRIBUTE_UNUSED = 
      kmul(em4phi,kmadd(trKL,kmul(g11,ToReal(-0.333333333333333333333333333333)),kxxL));
    
    CCTK_REAL_VEC At12L CCTK_ATTRIBUTE_UNUSED = 
      kmul(em4phi,kmadd(trKL,kmul(g12,ToReal(-0.333333333333333333333333333333)),kxyL));
    
    CCTK_REAL_VEC At13L CCTK_ATTRIBUTE_UNUSED = 
      kmul(em4phi,kmadd(trKL,kmul(g13,ToReal(-0.333333333333333333333333333333)),kxzL));
    
    CCTK_REAL_VEC At22L CCTK_ATTRIBUTE_UNUSED = 
      kmul(em4phi,kmadd(trKL,kmul(g22,ToReal(-0.333333333333333333333333333333)),kyyL));
    
    CCTK_REAL_VEC At23L CCTK_ATTRIBUTE_UNUSED = 
      kmul(em4phi,kmadd(trKL,kmul(g23,ToReal(-0.333333333333333333333333333333)),kyzL));
    
    CCTK_REAL_VEC At33L CCTK_ATTRIBUTE_UNUSED = 
      kmul(em4phi,kmadd(trKL,kmul(g33,ToReal(-0.333333333333333333333333333333)),kzzL));
    
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = alpL;
    
    CCTK_REAL_VEC beta1L CCTK_ATTRIBUTE_UNUSED = betaxL;
    
    CCTK_REAL_VEC beta2L CCTK_ATTRIBUTE_UNUSED = betayL;
    
    CCTK_REAL_VEC beta3L CCTK_ATTRIBUTE_UNUSED = betazL;
    
    CCTK_REAL_VEC xCopyL CCTK_ATTRIBUTE_UNUSED = xL;
    
    CCTK_REAL_VEC yCopyL CCTK_ATTRIBUTE_UNUSED = yL;
    
    CCTK_REAL_VEC zCopyL CCTK_ATTRIBUTE_UNUSED = zL;
    
    CCTK_REAL_VEC rCopyL CCTK_ATTRIBUTE_UNUSED = rL;
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
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
  CCTK_ENDLOOP3STR(ML_BSSN_Host_convertFromADMBase);
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
  
  const char* const groups[] = {
    "ADMBase::curv",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "ML_BSSN_Host::coords",
    "grid::coordinates",
    "ML_BSSN_Host::ML_curv",
    "ML_BSSN_Host::ML_lapse",
    "ML_BSSN_Host::ML_log_confac",
    "ML_BSSN_Host::ML_metric",
    "ML_BSSN_Host::ML_shift",
    "ML_BSSN_Host::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_Host_convertFromADMBase", 12, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, ML_BSSN_Host_convertFromADMBase_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_Host_convertFromADMBase_Body");
  }
}
