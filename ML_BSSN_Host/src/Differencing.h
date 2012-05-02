#include <assert.h>
#include "vectors.h"

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth1(u) (kmul(p1o840dx,kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,-3,0,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,4,0,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,-4,0,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(168),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(672)))))))))))
#else
#  define PDstandardNth1(u) (PDstandardNth1_impl(u,p1o840dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o840dx,kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,-3,0,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,4,0,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,-4,0,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(168),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(672))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth2(u) (kmul(p1o840dy,kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,0,-3,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,4,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,-4,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(168),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(672)))))))))))
#else
#  define PDstandardNth2(u) (PDstandardNth2_impl(u,p1o840dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o840dy,kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,0,-3,0),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,4,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,-4,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(168),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(672))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth3(u) (kmul(p1o840dz,kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,0,0,-3),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,0,4),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,0,-4),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(168),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(672)))))))))))
#else
#  define PDstandardNth3(u) (PDstandardNth3_impl(u,p1o840dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o840dz,kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-672),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-168),kmadd(KRANC_GFOFFSET3D(u,0,0,-3),ToReal(-32),kmadd(KRANC_GFOFFSET3D(u,0,0,4),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,0,-4),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(32),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(168),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(672))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth11(u) (kmul(p1o5040dx2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(8064))))))))
#else
#  define PDstandardNth11(u) (PDstandardNth11_impl(u,p1o5040dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth11_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dx2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth11_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dx2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o5040dx2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(8064)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth22(u) (kmul(p1o5040dy2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(8064))))))))
#else
#  define PDstandardNth22(u) (PDstandardNth22_impl(u,p1o5040dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dy2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dy2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o5040dy2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(8064)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth33(u) (kmul(p1o5040dz2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(8064))))))))
#else
#  define PDstandardNth33(u) (PDstandardNth33_impl(u,p1o5040dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth33_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dz2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth33_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o5040dz2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o5040dz2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-14350),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-1008),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(128),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(8064)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth12(u) (kmul(p1o705600dxdy,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNth12(u) (PDstandardNth12_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth12_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth12_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dxdy,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(451584))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth13(u) (kmul(p1o705600dxdz,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNth13(u) (PDstandardNth13_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth13_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth13_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dxdz,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(451584))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth21(u) (kmul(p1o705600dxdy,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNth21(u) (PDstandardNth21_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dxdy,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,3,0),kadd(KRANC_GFOFFSET3D(u,1,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,1,0),KRANC_GFOFFSET3D(u,3,-1,0)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-3,0),kadd(KRANC_GFOFFSET3D(u,2,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-2,0),KRANC_GFOFFSET3D(u,3,2,0)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-4,0),kadd(KRANC_GFOFFSET3D(u,1,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-1,0),KRANC_GFOFFSET3D(u,4,1,0)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,3,0),KRANC_GFOFFSET3D(u,3,-3,0)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,4,0),kadd(KRANC_GFOFFSET3D(u,2,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,2,0),KRANC_GFOFFSET3D(u,4,-2,0)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-4,0),kadd(KRANC_GFOFFSET3D(u,3,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-3,0),KRANC_GFOFFSET3D(u,4,3,0)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,4,0),KRANC_GFOFFSET3D(u,4,-4,0)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,-4,0),KRANC_GFOFFSET3D(u,4,4,0)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,4,0),kadd(KRANC_GFOFFSET3D(u,3,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,3,0),KRANC_GFOFFSET3D(u,4,-3,0)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-4,0),kadd(KRANC_GFOFFSET3D(u,2,4,0),kadd(KRANC_GFOFFSET3D(u,-4,-2,0),KRANC_GFOFFSET3D(u,4,2,0)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,-3,0),KRANC_GFOFFSET3D(u,3,3,0)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,4,0),kadd(KRANC_GFOFFSET3D(u,1,-4,0),kadd(KRANC_GFOFFSET3D(u,-4,1,0),KRANC_GFOFFSET3D(u,4,-1,0)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,3,0),kadd(KRANC_GFOFFSET3D(u,2,-3,0),kadd(KRANC_GFOFFSET3D(u,-3,2,0),KRANC_GFOFFSET3D(u,3,-2,0)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-3,0),kadd(KRANC_GFOFFSET3D(u,1,3,0),kadd(KRANC_GFOFFSET3D(u,-3,-1,0),KRANC_GFOFFSET3D(u,3,1,0)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,-2,0),KRANC_GFOFFSET3D(u,2,2,0)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(451584))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth23(u) (kmul(p1o705600dydz,kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNth23(u) (PDstandardNth23_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dydz,kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(451584))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth31(u) (kmul(p1o705600dxdz,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNth31(u) (PDstandardNth31_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth31_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth31_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dxdz,kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,3),kadd(KRANC_GFOFFSET3D(u,1,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,1),KRANC_GFOFFSET3D(u,3,0,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-3),kadd(KRANC_GFOFFSET3D(u,2,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-2),KRANC_GFOFFSET3D(u,3,0,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-4),kadd(KRANC_GFOFFSET3D(u,1,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-1),KRANC_GFOFFSET3D(u,4,0,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,3),KRANC_GFOFFSET3D(u,3,0,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,4),kadd(KRANC_GFOFFSET3D(u,2,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,2),KRANC_GFOFFSET3D(u,4,0,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-4),kadd(KRANC_GFOFFSET3D(u,3,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-3),KRANC_GFOFFSET3D(u,4,0,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,4),KRANC_GFOFFSET3D(u,4,0,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,-4),KRANC_GFOFFSET3D(u,4,0,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,4),kadd(KRANC_GFOFFSET3D(u,3,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,3),KRANC_GFOFFSET3D(u,4,0,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-4),kadd(KRANC_GFOFFSET3D(u,2,0,4),kadd(KRANC_GFOFFSET3D(u,-4,0,-2),KRANC_GFOFFSET3D(u,4,0,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,-3),KRANC_GFOFFSET3D(u,3,0,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,4),kadd(KRANC_GFOFFSET3D(u,1,0,-4),kadd(KRANC_GFOFFSET3D(u,-4,0,1),KRANC_GFOFFSET3D(u,4,0,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,3),kadd(KRANC_GFOFFSET3D(u,2,0,-3),kadd(KRANC_GFOFFSET3D(u,-3,0,2),KRANC_GFOFFSET3D(u,3,0,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-3),kadd(KRANC_GFOFFSET3D(u,1,0,3),kadd(KRANC_GFOFFSET3D(u,-3,0,-1),KRANC_GFOFFSET3D(u,3,0,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,-2),KRANC_GFOFFSET3D(u,2,0,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(451584))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth32(u) (kmul(p1o705600dydz,kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(451584)))))))))))))))))))))))
#else
#  define PDstandardNth32(u) (PDstandardNth32_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth32_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth32_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o705600dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o705600dydz,kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-451584),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-112896),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,2),KRANC_GFOFFSET3D(u,0,2,-2)),ToReal(-28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,3),kadd(KRANC_GFOFFSET3D(u,0,1,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,1),KRANC_GFOFFSET3D(u,0,3,-1)))),ToReal(-21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-3),kadd(KRANC_GFOFFSET3D(u,0,2,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-2),KRANC_GFOFFSET3D(u,0,3,2)))),ToReal(-5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-4),kadd(KRANC_GFOFFSET3D(u,0,1,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-1),KRANC_GFOFFSET3D(u,0,4,1)))),ToReal(-2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,3),KRANC_GFOFFSET3D(u,0,3,-3)),ToReal(-1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,4),kadd(KRANC_GFOFFSET3D(u,0,2,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,2),KRANC_GFOFFSET3D(u,0,4,-2)))),ToReal(-504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-4),kadd(KRANC_GFOFFSET3D(u,0,3,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-3),KRANC_GFOFFSET3D(u,0,4,3)))),ToReal(-96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,4),KRANC_GFOFFSET3D(u,0,4,-4)),ToReal(-9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,-4),KRANC_GFOFFSET3D(u,0,4,4)),ToReal(9),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,4),kadd(KRANC_GFOFFSET3D(u,0,3,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,3),KRANC_GFOFFSET3D(u,0,4,-3)))),ToReal(96),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-4),kadd(KRANC_GFOFFSET3D(u,0,2,4),kadd(KRANC_GFOFFSET3D(u,0,-4,-2),KRANC_GFOFFSET3D(u,0,4,2)))),ToReal(504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,-3),KRANC_GFOFFSET3D(u,0,3,3)),ToReal(1024),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,4),kadd(KRANC_GFOFFSET3D(u,0,1,-4),kadd(KRANC_GFOFFSET3D(u,0,-4,1),KRANC_GFOFFSET3D(u,0,4,-1)))),ToReal(2016),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,3),kadd(KRANC_GFOFFSET3D(u,0,2,-3),kadd(KRANC_GFOFFSET3D(u,0,-3,2),KRANC_GFOFFSET3D(u,0,3,-2)))),ToReal(5376),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-3),kadd(KRANC_GFOFFSET3D(u,0,1,3),kadd(KRANC_GFOFFSET3D(u,0,-3,-1),KRANC_GFOFFSET3D(u,0,3,1)))),ToReal(21504),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,-2),KRANC_GFOFFSET3D(u,0,2,2)),ToReal(28224),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(112896),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(451584))))))))))))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd1(u) (kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0))))
#else
#  define PDstandard2nd1(u) (PDstandard2nd1_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd2(u) (kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0))))
#else
#  define PDstandard2nd2(u) (PDstandard2nd2_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd3(u) (kmul(p1o2dz,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,-1))))
#else
#  define PDstandard2nd3(u) (PDstandard2nd3_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDstandard2nd3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandard2nd3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o2dz,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,-1)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth1(u) (kmul(p1o1024dx,kadd(KRANC_GFOFFSET3D(u,-5,0,0),kadd(KRANC_GFOFFSET3D(u,5,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(210))))))))))
#else
#  define PDdissipationNth1(u) (PDdissipationNth1_impl(u,p1o1024dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1024dx,kadd(KRANC_GFOFFSET3D(u,-5,0,0),kadd(KRANC_GFOFFSET3D(u,5,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(210)))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth2(u) (kmul(p1o1024dy,kadd(KRANC_GFOFFSET3D(u,0,-5,0),kadd(KRANC_GFOFFSET3D(u,0,5,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(210))))))))))
#else
#  define PDdissipationNth2(u) (PDdissipationNth2_impl(u,p1o1024dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1024dy,kadd(KRANC_GFOFFSET3D(u,0,-5,0),kadd(KRANC_GFOFFSET3D(u,0,5,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(210)))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth3(u) (kmul(p1o1024dz,kadd(KRANC_GFOFFSET3D(u,0,0,-5),kadd(KRANC_GFOFFSET3D(u,0,0,5),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(210))))))))))
#else
#  define PDdissipationNth3(u) (PDdissipationNth3_impl(u,p1o1024dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1024dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1024dz,kadd(KRANC_GFOFFSET3D(u,0,0,-5),kadd(KRANC_GFOFFSET3D(u,0,0,5),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(210)))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth1(u) (kmul(p1o840dx,kmul(dir1,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-378),kmadd(KRANC_GFOFFSET3D(u,5,0,0),ToReal(3),kmul(ToReal(-5),kadd(KRANC_GFOFFSET3D(u,-3,0,0),kmadd(KRANC_GFOFFSET3D(u,1,0,0),ToReal(-210),kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(-28),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(-12),kmadd(KRANC_GFOFFSET3D(u,4,0,0),ToReal(6),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(84)))))))))))))
#else
#  define PDupwindNth1(u) (PDupwindNth1_impl(u,p1o840dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti1(u) (kmul(p1o1680dx,kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,-3,0,0),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,4,0,0),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,-5,0,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,5,0,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,-4,0,0),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(480),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(1470)))))))))))))
#else
#  define PDupwindNthAnti1(u) (PDupwindNthAnti1_impl(u,p1o1680dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAnti1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAnti1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1680dx,kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,-3,0,0),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,4,0,0),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,-5,0,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,5,0,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,-4,0,0),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,3,0,0),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(480),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(1470))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm1(u) (kmul(p1o560dx,kadd(KRANC_GFOFFSET3D(u,-5,0,0),kadd(KRANC_GFOFFSET3D(u,5,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,-4,0,0),KRANC_GFOFFSET3D(u,4,0,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,-3,0,0),KRANC_GFOFFSET3D(u,3,0,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(210))))))))))
#else
#  define PDupwindNthSymm1(u) (PDupwindNthSymm1_impl(u,p1o560dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymm1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymm1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
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
#  define PDonesidedPlus2nd1(u) (kmul(pm1o2dx,kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(KRANC_GFOFFSET3D(u,1,0,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3))))))
#else
#  define PDonesidedPlus2nd1(u) (PDonesidedPlus2nd1_impl(u,pm1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDonesidedPlus2nd1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedPlus2nd1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o2dx,kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(KRANC_GFOFFSET3D(u,1,0,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedMinus2nd1(u) (kmul(p1o2dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3))))))
#else
#  define PDonesidedMinus2nd1(u) (PDonesidedMinus2nd1_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDonesidedMinus2nd1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedMinus2nd1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o2dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth2(u) (kmul(p1o840dy,kmul(dir2,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-378),kmadd(KRANC_GFOFFSET3D(u,0,5,0),ToReal(3),kmul(ToReal(-5),kadd(KRANC_GFOFFSET3D(u,0,-3,0),kmadd(KRANC_GFOFFSET3D(u,0,1,0),ToReal(-210),kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(-28),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(-12),kmadd(KRANC_GFOFFSET3D(u,0,4,0),ToReal(6),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(84)))))))))))))
#else
#  define PDupwindNth2(u) (PDupwindNth2_impl(u,p1o840dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti2(u) (kmul(p1o1680dy,kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,0,-3,0),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,0,4,0),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,0,-5,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,5,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,-4,0),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(480),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(1470)))))))))))))
#else
#  define PDupwindNthAnti2(u) (PDupwindNthAnti2_impl(u,p1o1680dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAnti2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAnti2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1680dy,kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,0,-3,0),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,0,4,0),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,0,-5,0),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,5,0),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,-4,0),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,0,3,0),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(480),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(1470))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm2(u) (kmul(p1o560dy,kadd(KRANC_GFOFFSET3D(u,0,-5,0),kadd(KRANC_GFOFFSET3D(u,0,5,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-4,0),KRANC_GFOFFSET3D(u,0,4,0)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-3,0),KRANC_GFOFFSET3D(u,0,3,0)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(210))))))))))
#else
#  define PDupwindNthSymm2(u) (PDupwindNthSymm2_impl(u,p1o560dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymm2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymm2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
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
#  define PDonesidedPlus2nd2(u) (kmul(pm1o2dy,kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(KRANC_GFOFFSET3D(u,0,1,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3))))))
#else
#  define PDonesidedPlus2nd2(u) (PDonesidedPlus2nd2_impl(u,pm1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDonesidedPlus2nd2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedPlus2nd2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o2dy,kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(KRANC_GFOFFSET3D(u,0,1,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedMinus2nd2(u) (kmul(p1o2dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3))))))
#else
#  define PDonesidedMinus2nd2(u) (PDonesidedMinus2nd2_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDonesidedMinus2nd2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedMinus2nd2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o2dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth3(u) (kmul(p1o840dz,kmul(dir3,kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-378),kmadd(KRANC_GFOFFSET3D(u,0,0,5),ToReal(3),kmul(ToReal(-5),kadd(KRANC_GFOFFSET3D(u,0,0,-3),kmadd(KRANC_GFOFFSET3D(u,0,0,1),ToReal(-210),kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(-28),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(-12),kmadd(KRANC_GFOFFSET3D(u,0,0,4),ToReal(6),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(84)))))))))))))
#else
#  define PDupwindNth3(u) (PDupwindNth3_impl(u,p1o840dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o840dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti3(u) (kmul(p1o1680dz,kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,0,0,-3),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,0,0,4),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,0,0,-5),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,0,5),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,0,-4),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(480),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(1470)))))))))))))
#else
#  define PDupwindNthAnti3(u) (PDupwindNthAnti3_impl(u,p1o1680dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAnti3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAnti3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o1680dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o1680dz,kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-1470),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-480),kmadd(KRANC_GFOFFSET3D(u,0,0,-3),ToReal(-145),kmadd(KRANC_GFOFFSET3D(u,0,0,4),ToReal(-30),kmadd(KRANC_GFOFFSET3D(u,0,0,-5),ToReal(-3),kmadd(KRANC_GFOFFSET3D(u,0,0,5),ToReal(3),kmadd(KRANC_GFOFFSET3D(u,0,0,-4),ToReal(30),kmadd(KRANC_GFOFFSET3D(u,0,0,3),ToReal(145),kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(480),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(1470))))))))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm3(u) (kmul(p1o560dz,kadd(KRANC_GFOFFSET3D(u,0,0,-5),kadd(KRANC_GFOFFSET3D(u,0,0,5),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-252),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-120),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-4),KRANC_GFOFFSET3D(u,0,0,4)),ToReal(-10),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-3),KRANC_GFOFFSET3D(u,0,0,3)),ToReal(45),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(210))))))))))
#else
#  define PDupwindNthSymm3(u) (PDupwindNthSymm3_impl(u,p1o560dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymm3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymm3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o560dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
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

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedPlus2nd3(u) (kmul(pm1o2dz,kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(KRANC_GFOFFSET3D(u,0,0,1),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3))))))
#else
#  define PDonesidedPlus2nd3(u) (PDonesidedPlus2nd3_impl(u,pm1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDonesidedPlus2nd3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedPlus2nd3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o2dz,kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(KRANC_GFOFFSET3D(u,0,0,1),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedMinus2nd3(u) (kmul(p1o2dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3))))))
#else
#  define PDonesidedMinus2nd3(u) (PDonesidedMinus2nd3_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDonesidedMinus2nd3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesidedMinus2nd3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o2dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDPlus1(u) (kmul(p1odx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,0,0,0))))
#else
#  define PDPlus1(u) (PDPlus1_impl(u,p1odx,cdj,cdk))
static CCTK_REAL_VEC PDPlus1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDPlus1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1odx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,0,0,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDPlus2(u) (kmul(p1ody,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,0,0))))
#else
#  define PDPlus2(u) (PDPlus2_impl(u,p1ody,cdj,cdk))
static CCTK_REAL_VEC PDPlus2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDPlus2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1ody,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,0,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDPlus3(u) (kmul(p1odz,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,0))))
#else
#  define PDPlus3(u) (PDPlus3_impl(u,p1odz,cdj,cdk))
static CCTK_REAL_VEC PDPlus3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDPlus3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1odz,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDMinus1(u) (kmul(p1odx,ksub(KRANC_GFOFFSET3D(u,0,0,0),KRANC_GFOFFSET3D(u,-1,0,0))))
#else
#  define PDMinus1(u) (PDMinus1_impl(u,p1odx,cdj,cdk))
static CCTK_REAL_VEC PDMinus1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDMinus1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1odx,ksub(KRANC_GFOFFSET3D(u,0,0,0),KRANC_GFOFFSET3D(u,-1,0,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDMinus2(u) (kmul(p1ody,ksub(KRANC_GFOFFSET3D(u,0,0,0),KRANC_GFOFFSET3D(u,0,-1,0))))
#else
#  define PDMinus2(u) (PDMinus2_impl(u,p1ody,cdj,cdk))
static CCTK_REAL_VEC PDMinus2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDMinus2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1ody,ksub(KRANC_GFOFFSET3D(u,0,0,0),KRANC_GFOFFSET3D(u,0,-1,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDMinus3(u) (kmul(p1odz,ksub(KRANC_GFOFFSET3D(u,0,0,0),KRANC_GFOFFSET3D(u,0,0,-1))))
#else
#  define PDMinus3(u) (PDMinus3_impl(u,p1odz,cdj,cdk))
static CCTK_REAL_VEC PDMinus3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDMinus3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1odz,ksub(KRANC_GFOFFSET3D(u,0,0,0),KRANC_GFOFFSET3D(u,0,0,-1)));
}
#endif

