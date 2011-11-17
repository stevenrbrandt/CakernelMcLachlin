#include <assert.h>
#include "vectors.h"

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder21(u) (kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0))))
#else
#  define PDstandardNthfdOrder21(u) (PDstandardNthfdOrder21_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder22(u) (kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0))))
#else
#  define PDstandardNthfdOrder22(u) (PDstandardNthfdOrder22_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder23(u) (kmul(p1o2dz,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,-1))))
#else
#  define PDstandardNthfdOrder23(u) (PDstandardNthfdOrder23_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o2dz,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,-1)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder41(u) (kmul(p1o12dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-8),kmsub(KRANC_GFOFFSET3D(u,1,0,0),ToReal(8),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDstandardNthfdOrder41(u) (PDstandardNthfdOrder41_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder41_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder41_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o12dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-8),kmsub(KRANC_GFOFFSET3D(u,1,0,0),ToReal(8),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder42(u) (kmul(p1o12dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-8),kmsub(KRANC_GFOFFSET3D(u,0,1,0),ToReal(8),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDstandardNthfdOrder42(u) (PDstandardNthfdOrder42_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder42_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder42_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o12dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-8),kmsub(KRANC_GFOFFSET3D(u,0,1,0),ToReal(8),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder43(u) (kmul(p1o12dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-8),kmsub(KRANC_GFOFFSET3D(u,0,0,1),ToReal(8),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDstandardNthfdOrder43(u) (PDstandardNthfdOrder43_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder43_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder43_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o12dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-8),kmsub(KRANC_GFOFFSET3D(u,0,0,1),ToReal(8),KRANC_GFOFFSET3D(u,0,0,2)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder61(u) (kmul(p1o60dx,kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-45),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-9),ksub(kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(9),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(45))),KRANC_GFOFFSET3D(u,-3,0,0)))))))
#else
#  define PDstandardNthfdOrder61(u) (PDstandardNthfdOrder61_impl(u,p1o60dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder61_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o60dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder61_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o60dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o60dx,kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-45),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-9),ksub(kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(9),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(45))),KRANC_GFOFFSET3D(u,-3,0,0))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder62(u) (kmul(p1o60dy,kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-45),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-9),ksub(kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(9),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(45))),KRANC_GFOFFSET3D(u,0,-3,0)))))))
#else
#  define PDstandardNthfdOrder62(u) (PDstandardNthfdOrder62_impl(u,p1o60dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder62_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o60dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder62_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o60dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o60dy,kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-45),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-9),ksub(kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(9),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(45))),KRANC_GFOFFSET3D(u,0,-3,0))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder63(u) (kmul(p1o60dz,kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-45),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-9),ksub(kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(9),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(45))),KRANC_GFOFFSET3D(u,0,0,-3)))))))
#else
#  define PDstandardNthfdOrder63(u) (PDstandardNthfdOrder63_impl(u,p1o60dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder63_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o60dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder63_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o60dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o60dz,kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-45),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-9),ksub(kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(9),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(45))),KRANC_GFOFFSET3D(u,0,0,-3))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder81(u) (kmul(p1o840dx,kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,-3,0,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,4,0,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,-4,0,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(168),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(672)))))))))))
#else
#  define PDstandardNthfdOrder81(u) (PDstandardNthfdOrder81_impl(u,p1o840dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder81_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder81_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o840dx,kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,-3,0,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,4,0,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,-4,0,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(168),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(672))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder82(u) (kmul(p1o840dy,kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,0,-3,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,4,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,-4,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(168),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(672)))))))))))
#else
#  define PDstandardNthfdOrder82(u) (PDstandardNthfdOrder82_impl(u,p1o840dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder82_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder82_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o840dy,kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,0,-3,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,4,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,-4,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(168),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(672))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder83(u) (kmul(p1o840dz,kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,0,0,-3),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,0,4),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,0,-4),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(168),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(672)))))))))))
#else
#  define PDstandardNthfdOrder83(u) (PDstandardNthfdOrder83_impl(u,p1o840dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder83_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder83_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o840dz,kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,0,0,-3),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,0,4),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,0,-4),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(168),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(672))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder211(u) (kmul(p1odx2,kadd(KRANC_GFOFFSET3D(u,-1,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-2),KRANC_GFOFFSET3D(u,1,0,0)))))
#else
#  define PDstandardNthfdOrder211(u) (PDstandardNthfdOrder211_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder211_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odx2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder211_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odx2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1odx2,kadd(KRANC_GFOFFSET3D(u,-1,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-2),KRANC_GFOFFSET3D(u,1,0,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder222(u) (kmul(p1ody2,kadd(KRANC_GFOFFSET3D(u,0,-1,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-2),KRANC_GFOFFSET3D(u,0,1,0)))))
#else
#  define PDstandardNthfdOrder222(u) (PDstandardNthfdOrder222_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder222_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1ody2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder222_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1ody2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1ody2,kadd(KRANC_GFOFFSET3D(u,0,-1,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-2),KRANC_GFOFFSET3D(u,0,1,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder233(u) (kmul(p1odz2,kadd(KRANC_GFOFFSET3D(u,0,0,-1),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-2),KRANC_GFOFFSET3D(u,0,0,1)))))
#else
#  define PDstandardNthfdOrder233(u) (PDstandardNthfdOrder233_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder233_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odz2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder233_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odz2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1odz2,kadd(KRANC_GFOFFSET3D(u,0,0,-1),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-2),KRANC_GFOFFSET3D(u,0,0,1))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder411(u) (kmul(pm1o12dx2,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-16),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(30)))))))
#else
#  define PDstandardNthfdOrder411(u) (PDstandardNthfdOrder411_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder411_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dx2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder411_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dx2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o12dx2,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-16),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(30))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder422(u) (kmul(pm1o12dy2,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-16),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(30)))))))
#else
#  define PDstandardNthfdOrder422(u) (PDstandardNthfdOrder422_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder422_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dy2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder422_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dy2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o12dy2,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-16),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(30))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder433(u) (kmul(pm1o12dz2,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-16),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(30)))))))
#else
#  define PDstandardNthfdOrder433(u) (PDstandardNthfdOrder433_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder433_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dz2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder433_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dz2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o12dz2,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-16),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(30))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder611(u) (kmul(p1o180dx2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-490),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-27),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(2),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(270)))))))
#else
#  define PDstandardNthfdOrder611(u) (PDstandardNthfdOrder611_impl(u,p1o180dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder611_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o180dx2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder611_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o180dx2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o180dx2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-490),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-27),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(2),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(270))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder622(u) (kmul(p1o180dy2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-490),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-27),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(2),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(270)))))))
#else
#  define PDstandardNthfdOrder622(u) (PDstandardNthfdOrder622_impl(u,p1o180dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder622_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o180dy2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder622_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o180dy2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o180dy2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-490),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-27),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(2),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(270))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder633(u) (kmul(p1o180dz2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-490),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-27),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(2),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(270)))))))
#else
#  define PDstandardNthfdOrder633(u) (PDstandardNthfdOrder633_impl(u,p1o180dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder633_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o180dz2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder633_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o180dz2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o180dz2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-490),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-27),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(2),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(270))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder811(u) (kmul(p1o5040dx2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(8064))))))))
#else
#  define PDstandardNthfdOrder811(u) (PDstandardNthfdOrder811_impl(u,p1o5040dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder811_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dx2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder811_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dx2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o5040dx2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(8064)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder822(u) (kmul(p1o5040dy2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(8064))))))))
#else
#  define PDstandardNthfdOrder822(u) (PDstandardNthfdOrder822_impl(u,p1o5040dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder822_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dy2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder822_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dy2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o5040dy2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(8064)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder833(u) (kmul(p1o5040dz2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(8064))))))))
#else
#  define PDstandardNthfdOrder833(u) (PDstandardNthfdOrder833_impl(u,p1o5040dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder833_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dz2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder833_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dz2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o5040dz2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(8064)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder212(u) (kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(KRANC_GFOFFSET3D(u,1,1,0),kadd(KRANC_GFOFFSET3D(u,1,-1,0),KRANC_GFOFFSET3D(u,-1,1,0))))))
#else
#  define PDstandardNthfdOrder212(u) (PDstandardNthfdOrder212_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder212_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder212_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(KRANC_GFOFFSET3D(u,1,1,0),kadd(KRANC_GFOFFSET3D(u,1,-1,0),KRANC_GFOFFSET3D(u,-1,1,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder213(u) (kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(KRANC_GFOFFSET3D(u,1,0,1),kadd(KRANC_GFOFFSET3D(u,1,0,-1),KRANC_GFOFFSET3D(u,-1,0,1))))))
#else
#  define PDstandardNthfdOrder213(u) (PDstandardNthfdOrder213_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder213_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder213_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(KRANC_GFOFFSET3D(u,1,0,1),kadd(KRANC_GFOFFSET3D(u,1,0,-1),KRANC_GFOFFSET3D(u,-1,0,1)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder221(u) (kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(KRANC_GFOFFSET3D(u,1,1,0),kadd(KRANC_GFOFFSET3D(u,1,-1,0),KRANC_GFOFFSET3D(u,-1,1,0))))))
#else
#  define PDstandardNthfdOrder221(u) (PDstandardNthfdOrder221_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder221_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder221_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(KRANC_GFOFFSET3D(u,1,1,0),kadd(KRANC_GFOFFSET3D(u,1,-1,0),KRANC_GFOFFSET3D(u,-1,1,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder223(u) (kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(KRANC_GFOFFSET3D(u,0,1,1),kadd(KRANC_GFOFFSET3D(u,0,1,-1),KRANC_GFOFFSET3D(u,0,-1,1))))))
#else
#  define PDstandardNthfdOrder223(u) (PDstandardNthfdOrder223_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder223_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder223_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(KRANC_GFOFFSET3D(u,0,1,1),kadd(KRANC_GFOFFSET3D(u,0,1,-1),KRANC_GFOFFSET3D(u,0,-1,1)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder231(u) (kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(KRANC_GFOFFSET3D(u,1,0,1),kadd(KRANC_GFOFFSET3D(u,1,0,-1),KRANC_GFOFFSET3D(u,-1,0,1))))))
#else
#  define PDstandardNthfdOrder231(u) (PDstandardNthfdOrder231_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder231_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder231_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(KRANC_GFOFFSET3D(u,1,0,1),kadd(KRANC_GFOFFSET3D(u,1,0,-1),KRANC_GFOFFSET3D(u,-1,0,1)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder232(u) (kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(KRANC_GFOFFSET3D(u,0,1,1),kadd(KRANC_GFOFFSET3D(u,0,1,-1),KRANC_GFOFFSET3D(u,0,-1,1))))))
#else
#  define PDstandardNthfdOrder232(u) (PDstandardNthfdOrder232_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder232_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder232_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(KRANC_GFOFFSET3D(u,0,1,1),kadd(KRANC_GFOFFSET3D(u,0,1,-1),KRANC_GFOFFSET3D(u,0,-1,1)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder412(u) (kmul(p1o144dxdy,kadd(KRANC_GFOFFSET3D(u,-2,-2,0),kadd(KRANC_GFOFFSET3D(u,2,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(64))),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0))))))))
#else
#  define PDstandardNthfdOrder412(u) (PDstandardNthfdOrder412_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder412_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder412_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dxdy,kadd(KRANC_GFOFFSET3D(u,-2,-2,0),kadd(KRANC_GFOFFSET3D(u,2,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(64))),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder413(u) (kmul(p1o144dxdz,kadd(KRANC_GFOFFSET3D(u,-2,0,-2),kadd(KRANC_GFOFFSET3D(u,2,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(64))),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2))))))))
#else
#  define PDstandardNthfdOrder413(u) (PDstandardNthfdOrder413_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder413_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder413_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dxdz,kadd(KRANC_GFOFFSET3D(u,-2,0,-2),kadd(KRANC_GFOFFSET3D(u,2,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(64))),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder421(u) (kmul(p1o144dxdy,kadd(KRANC_GFOFFSET3D(u,-2,-2,0),kadd(KRANC_GFOFFSET3D(u,2,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(64))),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0))))))))
#else
#  define PDstandardNthfdOrder421(u) (PDstandardNthfdOrder421_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder421_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder421_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dxdy,kadd(KRANC_GFOFFSET3D(u,-2,-2,0),kadd(KRANC_GFOFFSET3D(u,2,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(64))),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder423(u) (kmul(p1o144dydz,kadd(KRANC_GFOFFSET3D(u,0,-2,-2),kadd(KRANC_GFOFFSET3D(u,0,2,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(64))),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2))))))))
#else
#  define PDstandardNthfdOrder423(u) (PDstandardNthfdOrder423_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder423_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder423_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dydz,kadd(KRANC_GFOFFSET3D(u,0,-2,-2),kadd(KRANC_GFOFFSET3D(u,0,2,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(64))),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder431(u) (kmul(p1o144dxdz,kadd(KRANC_GFOFFSET3D(u,-2,0,-2),kadd(KRANC_GFOFFSET3D(u,2,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(64))),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2))))))))
#else
#  define PDstandardNthfdOrder431(u) (PDstandardNthfdOrder431_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder431_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder431_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dxdz,kadd(KRANC_GFOFFSET3D(u,-2,0,-2),kadd(KRANC_GFOFFSET3D(u,2,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(64))),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder432(u) (kmul(p1o144dydz,kadd(KRANC_GFOFFSET3D(u,0,-2,-2),kadd(KRANC_GFOFFSET3D(u,0,2,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(64))),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2))))))))
#else
#  define PDstandardNthfdOrder432(u) (PDstandardNthfdOrder432_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder432_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder432_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dydz,kadd(KRANC_GFOFFSET3D(u,0,-2,-2),kadd(KRANC_GFOFFSET3D(u,0,2,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(64))),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder612(u) (kmul(p1o3600dxdy,kadd(KRANC_GFOFFSET3D(u,-3,-3,0),kadd(KRANC_GFOFFSET3D(u,3,3,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,3,-3,0)),KRANC_GFOFFSET3D(u,-3,3,0)))))))))))
#else
#  define PDstandardNthfdOrder612(u) (PDstandardNthfdOrder612_impl(u,p1o3600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder612_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder612_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o3600dxdy,kadd(KRANC_GFOFFSET3D(u,-3,-3,0),kadd(KRANC_GFOFFSET3D(u,3,3,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,3,-3,0)),KRANC_GFOFFSET3D(u,-3,3,0))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder613(u) (kmul(p1o3600dxdz,kadd(KRANC_GFOFFSET3D(u,-3,0,-3),kadd(KRANC_GFOFFSET3D(u,3,0,3),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,3,0,-3)),KRANC_GFOFFSET3D(u,-3,0,3)))))))))))
#else
#  define PDstandardNthfdOrder613(u) (PDstandardNthfdOrder613_impl(u,p1o3600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder613_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder613_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o3600dxdz,kadd(KRANC_GFOFFSET3D(u,-3,0,-3),kadd(KRANC_GFOFFSET3D(u,3,0,3),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,3,0,-3)),KRANC_GFOFFSET3D(u,-3,0,3))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder621(u) (kmul(p1o3600dxdy,kadd(KRANC_GFOFFSET3D(u,-3,-3,0),kadd(KRANC_GFOFFSET3D(u,3,3,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,3,-3,0)),KRANC_GFOFFSET3D(u,-3,3,0)))))))))))
#else
#  define PDstandardNthfdOrder621(u) (PDstandardNthfdOrder621_impl(u,p1o3600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder621_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder621_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o3600dxdy,kadd(KRANC_GFOFFSET3D(u,-3,-3,0),kadd(KRANC_GFOFFSET3D(u,3,3,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,3,-3,0)),KRANC_GFOFFSET3D(u,-3,3,0))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder623(u) (kmul(p1o3600dydz,kadd(KRANC_GFOFFSET3D(u,0,-3,-3),kadd(KRANC_GFOFFSET3D(u,0,3,3),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,0,3,-3)),KRANC_GFOFFSET3D(u,0,-3,3)))))))))))
#else
#  define PDstandardNthfdOrder623(u) (PDstandardNthfdOrder623_impl(u,p1o3600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder623_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder623_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o3600dydz,kadd(KRANC_GFOFFSET3D(u,0,-3,-3),kadd(KRANC_GFOFFSET3D(u,0,3,3),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,0,3,-3)),KRANC_GFOFFSET3D(u,0,-3,3))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder631(u) (kmul(p1o3600dxdz,kadd(KRANC_GFOFFSET3D(u,-3,0,-3),kadd(KRANC_GFOFFSET3D(u,3,0,3),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,3,0,-3)),KRANC_GFOFFSET3D(u,-3,0,3)))))))))))
#else
#  define PDstandardNthfdOrder631(u) (PDstandardNthfdOrder631_impl(u,p1o3600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder631_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder631_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o3600dxdz,kadd(KRANC_GFOFFSET3D(u,-3,0,-3),kadd(KRANC_GFOFFSET3D(u,3,0,3),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,3,0,-3)),KRANC_GFOFFSET3D(u,-3,0,3))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder632(u) (kmul(p1o3600dydz,kadd(KRANC_GFOFFSET3D(u,0,-3,-3),kadd(KRANC_GFOFFSET3D(u,0,3,3),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,0,3,-3)),KRANC_GFOFFSET3D(u,0,-3,3)))))))))))
#else
#  define PDstandardNthfdOrder632(u) (PDstandardNthfdOrder632_impl(u,p1o3600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder632_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder632_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o3600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o3600dydz,kadd(KRANC_GFOFFSET3D(u,0,-3,-3),kadd(KRANC_GFOFFSET3D(u,0,3,3),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-2025),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-405),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-81),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-45),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-9),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(45),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(81),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(405),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(2025)))))),KRANC_GFOFFSET3D(u,0,3,-3)),KRANC_GFOFFSET3D(u,0,-3,3))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder812(u) (kmul(p1o705600dxdy,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder812(u) (PDstandardNthfdOrder812_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder812_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder812_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dxdy,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(451584))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder813(u) (kmul(p1o705600dxdz,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder813(u) (PDstandardNthfdOrder813_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder813_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder813_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dxdz,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(451584))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder821(u) (kmul(p1o705600dxdy,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder821(u) (PDstandardNthfdOrder821_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder821_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder821_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dxdy,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(451584))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder823(u) (kmul(p1o705600dydz,kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder823(u) (PDstandardNthfdOrder823_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder823_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder823_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dydz,kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(451584))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder831(u) (kmul(p1o705600dxdz,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder831(u) (PDstandardNthfdOrder831_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder831_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder831_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dxdz,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(451584))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder832(u) (kmul(p1o705600dydz,kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNthfdOrder832(u) (PDstandardNthfdOrder832_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder832_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder832_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dydz,kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(451584))))))))))))))))))))));
}
#endif

