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

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder21(u) (kmul(p1o16dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDdissipationNthfdOrder21(u) (PDdissipationNthfdOrder21_impl(u,p1o16dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o16dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o16dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o16dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder22(u) (kmul(p1o16dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDdissipationNthfdOrder22(u) (PDdissipationNthfdOrder22_impl(u,p1o16dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o16dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o16dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o16dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder23(u) (kmul(p1o16dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDdissipationNthfdOrder23(u) (PDdissipationNthfdOrder23_impl(u,p1o16dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o16dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o16dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o16dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder41(u) (kmul(p1o64dx,kadd(KRANC_GFOFFSET3D(u,-3,0,0),kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(15))))))))
#else
#  define PDdissipationNthfdOrder41(u) (PDdissipationNthfdOrder41_impl(u,p1o64dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder41_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder41_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o64dx,kadd(KRANC_GFOFFSET3D(u,-3,0,0),kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder42(u) (kmul(p1o64dy,kadd(KRANC_GFOFFSET3D(u,0,-3,0),kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(15))))))))
#else
#  define PDdissipationNthfdOrder42(u) (PDdissipationNthfdOrder42_impl(u,p1o64dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder42_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder42_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o64dy,kadd(KRANC_GFOFFSET3D(u,0,-3,0),kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder43(u) (kmul(p1o64dz,kadd(KRANC_GFOFFSET3D(u,0,0,-3),kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(15))))))))
#else
#  define PDdissipationNthfdOrder43(u) (PDdissipationNthfdOrder43_impl(u,p1o64dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder43_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder43_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o64dz,kadd(KRANC_GFOFFSET3D(u,0,0,-3),kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder61(u) (kmul(p1o256dx,kadd(KRANC_GFOFFSET3D(u,-4,0,0),kadd(KRANC_GFOFFSET3D(u,4,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70)))))))))
#else
#  define PDdissipationNthfdOrder61(u) (PDdissipationNthfdOrder61_impl(u,p1o256dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder61_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o256dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder61_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o256dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o256dx,kadd(KRANC_GFOFFSET3D(u,-4,0,0),kadd(KRANC_GFOFFSET3D(u,4,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder62(u) (kmul(p1o256dy,kadd(KRANC_GFOFFSET3D(u,0,-4,0),kadd(KRANC_GFOFFSET3D(u,0,4,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70)))))))))
#else
#  define PDdissipationNthfdOrder62(u) (PDdissipationNthfdOrder62_impl(u,p1o256dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder62_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o256dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder62_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o256dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o256dy,kadd(KRANC_GFOFFSET3D(u,0,-4,0),kadd(KRANC_GFOFFSET3D(u,0,4,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder63(u) (kmul(p1o256dz,kadd(KRANC_GFOFFSET3D(u,0,0,-4),kadd(KRANC_GFOFFSET3D(u,0,0,4),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70)))))))))
#else
#  define PDdissipationNthfdOrder63(u) (PDdissipationNthfdOrder63_impl(u,p1o256dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder63_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o256dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder63_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o256dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o256dz,kadd(KRANC_GFOFFSET3D(u,0,0,-4),kadd(KRANC_GFOFFSET3D(u,0,0,4),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder81(u) (kmul(p1o1024dx,kadd(KRANC_GFOFFSET3D(u,-5,0,0),kadd(KRANC_GFOFFSET3D(u,5,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(210))))))))))
#else
#  define PDdissipationNthfdOrder81(u) (PDdissipationNthfdOrder81_impl(u,p1o1024dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder81_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder81_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1024dx,kadd(KRANC_GFOFFSET3D(u,-5,0,0),kadd(KRANC_GFOFFSET3D(u,5,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(210)))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder82(u) (kmul(p1o1024dy,kadd(KRANC_GFOFFSET3D(u,0,-5,0),kadd(KRANC_GFOFFSET3D(u,0,5,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(210))))))))))
#else
#  define PDdissipationNthfdOrder82(u) (PDdissipationNthfdOrder82_impl(u,p1o1024dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder82_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder82_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1024dy,kadd(KRANC_GFOFFSET3D(u,0,-5,0),kadd(KRANC_GFOFFSET3D(u,0,5,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(210)))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder83(u) (kmul(p1o1024dz,kadd(KRANC_GFOFFSET3D(u,0,0,-5),kadd(KRANC_GFOFFSET3D(u,0,0,5),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(210))))))))))
#else
#  define PDdissipationNthfdOrder83(u) (PDdissipationNthfdOrder83_impl(u,p1o1024dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder83_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder83_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1024dz,kadd(KRANC_GFOFFSET3D(u,0,0,-5),kadd(KRANC_GFOFFSET3D(u,0,0,5),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(210)))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder21(u) (kmul(pm1o2dx,kmul(dir1,kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(KRANC_GFOFFSET3D(u,1,0,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))))))
#else
#  define PDupwindNthfdOrder21(u) (PDupwindNthfdOrder21_impl(u,pm1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder41(u) (kmul(p1o12dx,kmul(dir1,kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-10),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-6),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-3),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(18)))))))))
#else
#  define PDupwindNthfdOrder41(u) (PDupwindNthfdOrder41_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder41_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder41_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder61(u) (kmul(pm1o60dx,kmul(dir1,kadd(KRANC_GFOFFSET3D(u,4,0,0),kmadd(KRANC_GFOFFSET3D(u,1,0,0),ToReal(-80),kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(-8),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(-2),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(24),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(30),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(35)))))))))))
#else
#  define PDupwindNthfdOrder61(u) (PDupwindNthfdOrder61_impl(u,pm1o60dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder61_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o60dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder61_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o60dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder81(u) (kmul(p1o840dx,kmul(dir1,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-378),kmadd(KRANC_GFOFFSET3D(u,5,0,0),ToReal(3),kmul(ToReal(-5),kadd(KRANC_GFOFFSET3D(u,-3,0,0),kmadd(KRANC_GFOFFSET3D(u,1,0,0),ToReal(-210),kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(-28),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(-12),kmadd(KRANC_GFOFFSET3D(u,4,0,0),ToReal(6),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(84)))))))))))))
#else
#  define PDupwindNthfdOrder81(u) (PDupwindNthfdOrder81_impl(u,p1o840dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder81_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder81_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder21(u) (kmul(p1o4dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-4),kmsub(KRANC_GFOFFSET3D(u,1,0,0),ToReal(4),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDupwindNthAntifdOrder21(u) (PDupwindNthAntifdOrder21_impl(u,p1o4dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o4dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-4),kmsub(KRANC_GFOFFSET3D(u,1,0,0),ToReal(4),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder41(u) (kmul(p1o24dx,kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-21),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-6),ksub(kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(6),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(21))),KRANC_GFOFFSET3D(u,-3,0,0)))))))
#else
#  define PDupwindNthAntifdOrder41(u) (PDupwindNthAntifdOrder41_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder41_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder41_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dx,kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-21),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-6),ksub(kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(6),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(21))),KRANC_GFOFFSET3D(u,-3,0,0))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder61(u) (kmul(p1o120dx,kadd(KRANC_GFOFFSET3D(u,-4,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-104),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,-3,0,0),ToReal(-8),ksub(kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(8),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(32),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(104)))),KRANC_GFOFFSET3D(u,4,0,0))))))))
#else
#  define PDupwindNthAntifdOrder61(u) (PDupwindNthAntifdOrder61_impl(u,p1o120dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder61_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o120dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder61_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o120dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o120dx,kadd(KRANC_GFOFFSET3D(u,-4,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-104),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,-3,0,0),ToReal(-8),ksub(kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(8),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(32),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(104)))),KRANC_GFOFFSET3D(u,4,0,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder81(u) (kmul(p1o1680dx,kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,-3,0,0),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,4,0,0),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,-5,0,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,5,0,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,-4,0,0),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(480),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(1470)))))))))))))
#else
#  define PDupwindNthAntifdOrder81(u) (PDupwindNthAntifdOrder81_impl(u,p1o1680dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder81_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder81_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1680dx,kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,-3,0,0),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,4,0,0),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,-5,0,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,5,0,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,-4,0,0),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(480),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(1470))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder21(u) (kmul(pm1o4dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDupwindNthSymmfdOrder21(u) (PDupwindNthSymmfdOrder21_impl(u,pm1o4dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o4dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o4dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o4dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder41(u) (kmul(p1o24dx,kadd(KRANC_GFOFFSET3D(u,-3,0,0),kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(15))))))))
#else
#  define PDupwindNthSymmfdOrder41(u) (PDupwindNthSymmfdOrder41_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder41_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder41_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dx,kadd(KRANC_GFOFFSET3D(u,-3,0,0),kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder61(u) (kmul(pm1o120dx,kadd(KRANC_GFOFFSET3D(u,-4,0,0),kadd(KRANC_GFOFFSET3D(u,4,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70)))))))))
#else
#  define PDupwindNthSymmfdOrder61(u) (PDupwindNthSymmfdOrder61_impl(u,pm1o120dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder61_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o120dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder61_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o120dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o120dx,kadd(KRANC_GFOFFSET3D(u,-4,0,0),kadd(KRANC_GFOFFSET3D(u,4,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder81(u) (kmul(p1o560dx,kadd(KRANC_GFOFFSET3D(u,-5,0,0),kadd(KRANC_GFOFFSET3D(u,5,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(210))))))))))
#else
#  define PDupwindNthSymmfdOrder81(u) (PDupwindNthSymmfdOrder81_impl(u,p1o560dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder81_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder81_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o560dx,kadd(KRANC_GFOFFSET3D(u,-5,0,0),kadd(KRANC_GFOFFSET3D(u,5,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(210)))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided1(u) (kmul(p1odx,kmul(dir1,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,0,0,0)))))
#else
#  define PDonesided1(u) (PDonesided1_impl(u,p1odx,cdj,cdk))
static CCTK_REAL_VEC PDonesided1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesided1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder22(u) (kmul(pm1o2dy,kmul(dir2,kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(KRANC_GFOFFSET3D(u,0,1,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))))))
#else
#  define PDupwindNthfdOrder22(u) (PDupwindNthfdOrder22_impl(u,pm1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder42(u) (kmul(p1o12dy,kmul(dir2,kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-10),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-6),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-3),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(18)))))))))
#else
#  define PDupwindNthfdOrder42(u) (PDupwindNthfdOrder42_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder42_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder42_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder62(u) (kmul(pm1o60dy,kmul(dir2,kadd(KRANC_GFOFFSET3D(u,0,4,0),kmadd(KRANC_GFOFFSET3D(u,0,1,0),ToReal(-80),kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(-8),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(-2),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(24),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(30),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(35)))))))))))
#else
#  define PDupwindNthfdOrder62(u) (PDupwindNthfdOrder62_impl(u,pm1o60dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder62_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o60dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder62_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o60dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder82(u) (kmul(p1o840dy,kmul(dir2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-378),kmadd(KRANC_GFOFFSET3D(u,0,5,0),ToReal(3),kmul(ToReal(-5),kadd(KRANC_GFOFFSET3D(u,0,-3,0),kmadd(KRANC_GFOFFSET3D(u,0,1,0),ToReal(-210),kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(-28),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(-12),kmadd(KRANC_GFOFFSET3D(u,0,4,0),ToReal(6),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(84)))))))))))))
#else
#  define PDupwindNthfdOrder82(u) (PDupwindNthfdOrder82_impl(u,p1o840dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder82_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder82_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder22(u) (kmul(p1o4dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-4),kmsub(KRANC_GFOFFSET3D(u,0,1,0),ToReal(4),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDupwindNthAntifdOrder22(u) (PDupwindNthAntifdOrder22_impl(u,p1o4dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o4dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-4),kmsub(KRANC_GFOFFSET3D(u,0,1,0),ToReal(4),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder42(u) (kmul(p1o24dy,kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-21),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-6),ksub(kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(6),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(21))),KRANC_GFOFFSET3D(u,0,-3,0)))))))
#else
#  define PDupwindNthAntifdOrder42(u) (PDupwindNthAntifdOrder42_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder42_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder42_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dy,kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-21),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-6),ksub(kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(6),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(21))),KRANC_GFOFFSET3D(u,0,-3,0))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder62(u) (kmul(p1o120dy,kadd(KRANC_GFOFFSET3D(u,0,-4,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-104),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,-3,0),ToReal(-8),ksub(kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(8),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(32),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(104)))),KRANC_GFOFFSET3D(u,0,4,0))))))))
#else
#  define PDupwindNthAntifdOrder62(u) (PDupwindNthAntifdOrder62_impl(u,p1o120dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder62_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o120dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder62_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o120dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o120dy,kadd(KRANC_GFOFFSET3D(u,0,-4,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-104),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,-3,0),ToReal(-8),ksub(kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(8),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(32),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(104)))),KRANC_GFOFFSET3D(u,0,4,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder82(u) (kmul(p1o1680dy,kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,0,-3,0),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,0,4,0),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,0,-5,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,5,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,-4,0),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(480),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(1470)))))))))))))
#else
#  define PDupwindNthAntifdOrder82(u) (PDupwindNthAntifdOrder82_impl(u,p1o1680dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder82_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder82_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1680dy,kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,0,-3,0),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,0,4,0),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,0,-5,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,5,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,-4,0),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(480),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(1470))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder22(u) (kmul(pm1o4dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDupwindNthSymmfdOrder22(u) (PDupwindNthSymmfdOrder22_impl(u,pm1o4dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o4dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o4dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o4dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder42(u) (kmul(p1o24dy,kadd(KRANC_GFOFFSET3D(u,0,-3,0),kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(15))))))))
#else
#  define PDupwindNthSymmfdOrder42(u) (PDupwindNthSymmfdOrder42_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder42_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder42_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dy,kadd(KRANC_GFOFFSET3D(u,0,-3,0),kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder62(u) (kmul(pm1o120dy,kadd(KRANC_GFOFFSET3D(u,0,-4,0),kadd(KRANC_GFOFFSET3D(u,0,4,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70)))))))))
#else
#  define PDupwindNthSymmfdOrder62(u) (PDupwindNthSymmfdOrder62_impl(u,pm1o120dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder62_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o120dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder62_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o120dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o120dy,kadd(KRANC_GFOFFSET3D(u,0,-4,0),kadd(KRANC_GFOFFSET3D(u,0,4,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder82(u) (kmul(p1o560dy,kadd(KRANC_GFOFFSET3D(u,0,-5,0),kadd(KRANC_GFOFFSET3D(u,0,5,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(210))))))))))
#else
#  define PDupwindNthSymmfdOrder82(u) (PDupwindNthSymmfdOrder82_impl(u,p1o560dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder82_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder82_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o560dy,kadd(KRANC_GFOFFSET3D(u,0,-5,0),kadd(KRANC_GFOFFSET3D(u,0,5,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(210)))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided2(u) (kmul(p1ody,kmul(dir2,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,0,0)))))
#else
#  define PDonesided2(u) (PDonesided2_impl(u,p1ody,cdj,cdk))
static CCTK_REAL_VEC PDonesided2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesided2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder23(u) (kmul(pm1o2dz,kmul(dir3,kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(KRANC_GFOFFSET3D(u,0,0,1),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))))))
#else
#  define PDupwindNthfdOrder23(u) (PDupwindNthfdOrder23_impl(u,pm1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder43(u) (kmul(p1o12dz,kmul(dir3,kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-10),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-6),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-3),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(18)))))))))
#else
#  define PDupwindNthfdOrder43(u) (PDupwindNthfdOrder43_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder43_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder43_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder63(u) (kmul(pm1o60dz,kmul(dir3,kadd(KRANC_GFOFFSET3D(u,0,0,4),kmadd(KRANC_GFOFFSET3D(u,0,0,1),ToReal(-80),kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(-8),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(-2),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(24),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(30),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(35)))))))))))
#else
#  define PDupwindNthfdOrder63(u) (PDupwindNthfdOrder63_impl(u,pm1o60dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder63_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o60dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder63_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o60dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder83(u) (kmul(p1o840dz,kmul(dir3,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-378),kmadd(KRANC_GFOFFSET3D(u,0,0,5),ToReal(3),kmul(ToReal(-5),kadd(KRANC_GFOFFSET3D(u,0,0,-3),kmadd(KRANC_GFOFFSET3D(u,0,0,1),ToReal(-210),kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(-28),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(-12),kmadd(KRANC_GFOFFSET3D(u,0,0,4),ToReal(6),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(84)))))))))))))
#else
#  define PDupwindNthfdOrder83(u) (PDupwindNthfdOrder83_impl(u,p1o840dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder83_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder83_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder23(u) (kmul(p1o4dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-4),kmsub(KRANC_GFOFFSET3D(u,0,0,1),ToReal(4),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDupwindNthAntifdOrder23(u) (PDupwindNthAntifdOrder23_impl(u,p1o4dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o4dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o4dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-4),kmsub(KRANC_GFOFFSET3D(u,0,0,1),ToReal(4),KRANC_GFOFFSET3D(u,0,0,2)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder43(u) (kmul(p1o24dz,kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-21),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-6),ksub(kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(6),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(21))),KRANC_GFOFFSET3D(u,0,0,-3)))))))
#else
#  define PDupwindNthAntifdOrder43(u) (PDupwindNthAntifdOrder43_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder43_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder43_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dz,kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-21),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-6),ksub(kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(6),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(21))),KRANC_GFOFFSET3D(u,0,0,-3))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder63(u) (kmul(p1o120dz,kadd(KRANC_GFOFFSET3D(u,0,0,-4),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-104),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,0,-3),ToReal(-8),ksub(kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(8),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(32),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(104)))),KRANC_GFOFFSET3D(u,0,0,4))))))))
#else
#  define PDupwindNthAntifdOrder63(u) (PDupwindNthAntifdOrder63_impl(u,p1o120dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder63_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o120dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder63_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o120dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o120dz,kadd(KRANC_GFOFFSET3D(u,0,0,-4),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-104),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,0,-3),ToReal(-8),ksub(kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(8),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(32),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(104)))),KRANC_GFOFFSET3D(u,0,0,4)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder83(u) (kmul(p1o1680dz,kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,0,0,-3),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,0,0,4),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,0,0,-5),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,0,5),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,0,-4),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(480),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(1470)))))))))))))
#else
#  define PDupwindNthAntifdOrder83(u) (PDupwindNthAntifdOrder83_impl(u,p1o1680dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder83_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder83_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1680dz,kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,0,0,-3),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,0,0,4),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,0,0,-5),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,0,5),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,0,-4),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(480),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(1470))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder23(u) (kmul(pm1o4dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDupwindNthSymmfdOrder23(u) (PDupwindNthSymmfdOrder23_impl(u,pm1o4dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o4dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o4dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o4dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder43(u) (kmul(p1o24dz,kadd(KRANC_GFOFFSET3D(u,0,0,-3),kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(15))))))))
#else
#  define PDupwindNthSymmfdOrder43(u) (PDupwindNthSymmfdOrder43_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder43_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder43_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dz,kadd(KRANC_GFOFFSET3D(u,0,0,-3),kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder63(u) (kmul(pm1o120dz,kadd(KRANC_GFOFFSET3D(u,0,0,-4),kadd(KRANC_GFOFFSET3D(u,0,0,4),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70)))))))))
#else
#  define PDupwindNthSymmfdOrder63(u) (PDupwindNthSymmfdOrder63_impl(u,pm1o120dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder63_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o120dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder63_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o120dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o120dz,kadd(KRANC_GFOFFSET3D(u,0,0,-4),kadd(KRANC_GFOFFSET3D(u,0,0,4),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-56),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(-8),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(28),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(70))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder83(u) (kmul(p1o560dz,kadd(KRANC_GFOFFSET3D(u,0,0,-5),kadd(KRANC_GFOFFSET3D(u,0,0,5),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(210))))))))))
#else
#  define PDupwindNthSymmfdOrder83(u) (PDupwindNthSymmfdOrder83_impl(u,p1o560dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder83_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder83_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o560dz,kadd(KRANC_GFOFFSET3D(u,0,0,-5),kadd(KRANC_GFOFFSET3D(u,0,0,5),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(210)))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided3(u) (kmul(p1odz,kmul(dir3,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,0)))))
#else
#  define PDonesided3(u) (PDonesided3_impl(u,p1odz,cdj,cdk))
static CCTK_REAL_VEC PDonesided3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesided3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

