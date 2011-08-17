#include <assert.h>
#include "vectors.h"

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth1(u) (kmul(p1o12dx,kadd(vec_loadu_maybe3(-2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]),kmadd(vec_loadu_maybe3(-1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]),ToReal(-8),kmsub(vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]),ToReal(8),vec_loadu_maybe3(2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]))))))
#else
#  define PDstandardNth1(u) (PDstandardNth1_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o12dx,kadd(vec_loadu_maybe3(-2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]),kmadd(vec_loadu_maybe3(-1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]),ToReal(-8),kmsub(vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]),ToReal(8),vec_loadu_maybe3(2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)])))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth2(u) (kmul(p1o12dy,kadd(vec_loadu_maybe3(0,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]),ToReal(-8),kmsub(vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]),ToReal(8),vec_loadu_maybe3(0,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]))))))
#else
#  define PDstandardNth2(u) (PDstandardNth2_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o12dy,kadd(vec_loadu_maybe3(0,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]),ToReal(-8),kmsub(vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]),ToReal(8),vec_loadu_maybe3(0,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)])))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth3(u) (kmul(p1o12dz,kadd(vec_loadu_maybe3(0,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]),kmadd(vec_loadu_maybe3(0,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]),ToReal(-8),kmsub(vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]),ToReal(8),vec_loadu_maybe3(0,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]))))))
#else
#  define PDstandardNth3(u) (PDstandardNth3_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o12dz,kadd(vec_loadu_maybe3(0,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]),kmadd(vec_loadu_maybe3(0,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]),ToReal(-8),kmsub(vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]),ToReal(8),vec_loadu_maybe3(0,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)])))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth11(u) (kmul(pm1o12dx2,kadd(vec_loadu_maybe3(-2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]),kadd(vec_loadu_maybe3(2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]),kmadd(kadd(vec_loadu_maybe3(-1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]),vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])),ToReal(-16),kmul(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(30)))))))
#else
#  define PDstandardNth11(u) (PDstandardNth11_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth11_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dx2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth11_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dx2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o12dx2,kadd(vec_loadu_maybe3(-2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]),kadd(vec_loadu_maybe3(2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]),kmadd(kadd(vec_loadu_maybe3(-1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]),vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])),ToReal(-16),kmul(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(30))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth22(u) (kmul(pm1o12dy2,kadd(vec_loadu_maybe3(0,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(0,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]),kmadd(kadd(vec_loadu_maybe3(0,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])),ToReal(-16),kmul(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(30)))))))
#else
#  define PDstandardNth22(u) (PDstandardNth22_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dy2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth22_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dy2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o12dy2,kadd(vec_loadu_maybe3(0,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(0,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]),kmadd(kadd(vec_loadu_maybe3(0,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])),ToReal(-16),kmul(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(30))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth33(u) (kmul(pm1o12dz2,kadd(vec_loadu_maybe3(0,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]),kmadd(kadd(vec_loadu_maybe3(0,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])),ToReal(-16),kmul(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(30)))))))
#else
#  define PDstandardNth33(u) (PDstandardNth33_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth33_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dz2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth33_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const pm1o12dz2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(pm1o12dz2,kadd(vec_loadu_maybe3(0,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]),kmadd(kadd(vec_loadu_maybe3(0,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])),ToReal(-16),kmul(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(30))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth12(u) (kmul(p1o144dxdy,kadd(vec_loadu_maybe3(-2,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(2,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(2)+cdk*(0)]),kmadd(kadd(vec_loadu_maybe3(-1,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]),vec_loadu_maybe3(1,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(-1,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(1,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(2)+cdk*(0)]),kadd(vec_loadu_maybe3(-2,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(2,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(1)+cdk*(0)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(-1,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(2)+cdk*(0)]),kadd(vec_loadu_maybe3(1,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(-2,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(1)+cdk*(0)]),vec_loadu_maybe3(2,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-1)+cdk*(0)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(-1,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(1,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)])),ToReal(64))),vec_loadu_maybe3(2,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-2)+cdk*(0)])),vec_loadu_maybe3(-2,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(2)+cdk*(0)]))))))))
#else
#  define PDstandardNth12(u) (PDstandardNth12_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth12_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth12_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dxdy,kadd(vec_loadu_maybe3(-2,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(2,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(2)+cdk*(0)]),kmadd(kadd(vec_loadu_maybe3(-1,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]),vec_loadu_maybe3(1,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(-1,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(1,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(2)+cdk*(0)]),kadd(vec_loadu_maybe3(-2,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(2,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(1)+cdk*(0)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(-1,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(2)+cdk*(0)]),kadd(vec_loadu_maybe3(1,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(-2,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(1)+cdk*(0)]),vec_loadu_maybe3(2,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-1)+cdk*(0)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(-1,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(1,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)])),ToReal(64))),vec_loadu_maybe3(2,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-2)+cdk*(0)])),vec_loadu_maybe3(-2,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(2)+cdk*(0)])))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth13(u) (kmul(p1o144dxdz,kadd(vec_loadu_maybe3(-2,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(2,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(2)]),kmadd(kadd(vec_loadu_maybe3(-1,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]),vec_loadu_maybe3(1,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(-1,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(1,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(2)]),kadd(vec_loadu_maybe3(-2,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(2,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(1)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(-1,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(2)]),kadd(vec_loadu_maybe3(1,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(-2,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(1)]),vec_loadu_maybe3(2,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-1)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(-1,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(1,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)])),ToReal(64))),vec_loadu_maybe3(2,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-2)])),vec_loadu_maybe3(-2,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(2)]))))))))
#else
#  define PDstandardNth13(u) (PDstandardNth13_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth13_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth13_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dxdz,kadd(vec_loadu_maybe3(-2,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(2,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(2)]),kmadd(kadd(vec_loadu_maybe3(-1,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]),vec_loadu_maybe3(1,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(-1,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(1,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(2)]),kadd(vec_loadu_maybe3(-2,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(2,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(1)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(-1,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(2)]),kadd(vec_loadu_maybe3(1,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(-2,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(1)]),vec_loadu_maybe3(2,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-1)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(-1,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(1,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)])),ToReal(64))),vec_loadu_maybe3(2,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-2)])),vec_loadu_maybe3(-2,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(2)])))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth21(u) (kmul(p1o144dxdy,kadd(vec_loadu_maybe3(-2,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(2,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(2)+cdk*(0)]),kmadd(kadd(vec_loadu_maybe3(-1,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]),vec_loadu_maybe3(1,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(-1,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(1,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(2)+cdk*(0)]),kadd(vec_loadu_maybe3(-2,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(2,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(1)+cdk*(0)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(-1,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(2)+cdk*(0)]),kadd(vec_loadu_maybe3(1,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(-2,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(1)+cdk*(0)]),vec_loadu_maybe3(2,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-1)+cdk*(0)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(-1,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(1,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)])),ToReal(64))),vec_loadu_maybe3(2,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-2)+cdk*(0)])),vec_loadu_maybe3(-2,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(2)+cdk*(0)]))))))))
#else
#  define PDstandardNth21(u) (PDstandardNth21_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth21_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dxdy,kadd(vec_loadu_maybe3(-2,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(2,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(2)+cdk*(0)]),kmadd(kadd(vec_loadu_maybe3(-1,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]),vec_loadu_maybe3(1,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(-1,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(1,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(2)+cdk*(0)]),kadd(vec_loadu_maybe3(-2,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(2,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(1)+cdk*(0)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(-1,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(2)+cdk*(0)]),kadd(vec_loadu_maybe3(1,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-2)+cdk*(0)]),kadd(vec_loadu_maybe3(-2,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(1)+cdk*(0)]),vec_loadu_maybe3(2,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-1)+cdk*(0)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(-1,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(1,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)])),ToReal(64))),vec_loadu_maybe3(2,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-2)+cdk*(0)])),vec_loadu_maybe3(-2,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(2)+cdk*(0)])))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth23(u) (kmul(p1o144dydz,kadd(vec_loadu_maybe3(0,-2,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,2,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(2)]),kmadd(kadd(vec_loadu_maybe3(0,-1,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]),vec_loadu_maybe3(0,1,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(0,-1,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,1,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(2)]),kadd(vec_loadu_maybe3(0,-2,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-1)]),vec_loadu_maybe3(0,2,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(1)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(0,-1,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(2)]),kadd(vec_loadu_maybe3(0,1,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,-2,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(1)]),vec_loadu_maybe3(0,2,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-1)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(0,-1,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]),vec_loadu_maybe3(0,1,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)])),ToReal(64))),vec_loadu_maybe3(0,2,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-2)])),vec_loadu_maybe3(0,-2,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(2)]))))))))
#else
#  define PDstandardNth23(u) (PDstandardNth23_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth23_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dydz,kadd(vec_loadu_maybe3(0,-2,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,2,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(2)]),kmadd(kadd(vec_loadu_maybe3(0,-1,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]),vec_loadu_maybe3(0,1,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(0,-1,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,1,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(2)]),kadd(vec_loadu_maybe3(0,-2,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-1)]),vec_loadu_maybe3(0,2,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(1)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(0,-1,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(2)]),kadd(vec_loadu_maybe3(0,1,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,-2,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(1)]),vec_loadu_maybe3(0,2,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-1)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(0,-1,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]),vec_loadu_maybe3(0,1,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)])),ToReal(64))),vec_loadu_maybe3(0,2,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-2)])),vec_loadu_maybe3(0,-2,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(2)])))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth31(u) (kmul(p1o144dxdz,kadd(vec_loadu_maybe3(-2,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(2,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(2)]),kmadd(kadd(vec_loadu_maybe3(-1,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]),vec_loadu_maybe3(1,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(-1,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(1,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(2)]),kadd(vec_loadu_maybe3(-2,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(2,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(1)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(-1,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(2)]),kadd(vec_loadu_maybe3(1,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(-2,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(1)]),vec_loadu_maybe3(2,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-1)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(-1,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(1,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)])),ToReal(64))),vec_loadu_maybe3(2,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-2)])),vec_loadu_maybe3(-2,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(2)]))))))))
#else
#  define PDstandardNth31(u) (PDstandardNth31_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth31_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth31_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dxdz,kadd(vec_loadu_maybe3(-2,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(2,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(2)]),kmadd(kadd(vec_loadu_maybe3(-1,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]),vec_loadu_maybe3(1,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(-1,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(1,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(2)]),kadd(vec_loadu_maybe3(-2,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(2,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(1)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(-1,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(2)]),kadd(vec_loadu_maybe3(1,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-2)]),kadd(vec_loadu_maybe3(-2,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(1)]),vec_loadu_maybe3(2,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-1)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(-1,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(1,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)])),ToReal(64))),vec_loadu_maybe3(2,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-2)])),vec_loadu_maybe3(-2,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(2)])))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth32(u) (kmul(p1o144dydz,kadd(vec_loadu_maybe3(0,-2,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,2,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(2)]),kmadd(kadd(vec_loadu_maybe3(0,-1,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]),vec_loadu_maybe3(0,1,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(0,-1,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,1,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(2)]),kadd(vec_loadu_maybe3(0,-2,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-1)]),vec_loadu_maybe3(0,2,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(1)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(0,-1,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(2)]),kadd(vec_loadu_maybe3(0,1,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,-2,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(1)]),vec_loadu_maybe3(0,2,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-1)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(0,-1,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]),vec_loadu_maybe3(0,1,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)])),ToReal(64))),vec_loadu_maybe3(0,2,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-2)])),vec_loadu_maybe3(0,-2,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(2)]))))))))
#else
#  define PDstandardNth32(u) (PDstandardNth32_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNth32_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNth32_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o144dydz,kadd(vec_loadu_maybe3(0,-2,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,2,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(2)]),kmadd(kadd(vec_loadu_maybe3(0,-1,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]),vec_loadu_maybe3(0,1,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)])),ToReal(-64),kmadd(kadd(vec_loadu_maybe3(0,-1,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,1,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(2)]),kadd(vec_loadu_maybe3(0,-2,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-1)]),vec_loadu_maybe3(0,2,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(1)])))),ToReal(-8),ksub(ksub(kmadd(kadd(vec_loadu_maybe3(0,-1,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(2)]),kadd(vec_loadu_maybe3(0,1,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-2)]),kadd(vec_loadu_maybe3(0,-2,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(1)]),vec_loadu_maybe3(0,2,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-1)])))),ToReal(8),kmul(kadd(vec_loadu_maybe3(0,-1,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]),vec_loadu_maybe3(0,1,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)])),ToReal(64))),vec_loadu_maybe3(0,2,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-2)])),vec_loadu_maybe3(0,-2,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(2)])))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth1(u) (kmul(p1o64dx,kadd(vec_loadu_maybe3(-3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-3)+cdj*(0)+cdk*(0)]),kadd(vec_loadu_maybe3(3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(3)+cdj*(0)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(-2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]),vec_loadu_maybe3(2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(-1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]),vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])),ToReal(15))))))))
#else
#  define PDdissipationNth1(u) (PDdissipationNth1_impl(u,p1o64dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o64dx,kadd(vec_loadu_maybe3(-3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-3)+cdj*(0)+cdk*(0)]),kadd(vec_loadu_maybe3(3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(3)+cdj*(0)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(-2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]),vec_loadu_maybe3(2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(-1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]),vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth2(u) (kmul(p1o64dy,kadd(vec_loadu_maybe3(0,-3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-3)+cdk*(0)]),kadd(vec_loadu_maybe3(0,3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(3)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(0,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]),vec_loadu_maybe3(0,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(0,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])),ToReal(15))))))))
#else
#  define PDdissipationNth2(u) (PDdissipationNth2_impl(u,p1o64dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o64dy,kadd(vec_loadu_maybe3(0,-3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-3)+cdk*(0)]),kadd(vec_loadu_maybe3(0,3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(3)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(0,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]),vec_loadu_maybe3(0,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(0,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth3(u) (kmul(p1o64dz,kadd(vec_loadu_maybe3(0,0,-3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-3)]),kadd(vec_loadu_maybe3(0,0,3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(3)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(0,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]),vec_loadu_maybe3(0,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(0,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])),ToReal(15))))))))
#else
#  define PDdissipationNth3(u) (PDdissipationNth3_impl(u,p1o64dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o64dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o64dz,kadd(vec_loadu_maybe3(0,0,-3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-3)]),kadd(vec_loadu_maybe3(0,0,3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(3)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(0,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]),vec_loadu_maybe3(0,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(0,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth1(u) (kmul(p1o12dx,kmul(dir1,kadd(vec_loadu_maybe3(3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(3*dir1)+cdj*(0)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-10),kmadd(vec_loadu_maybe3(2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2*dir1)+cdj*(0)+cdk*(0)]),ToReal(-6),kmadd(vec_loadu_maybe3(-1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-dir1)+cdj*(0)+cdk*(0)]),ToReal(-3),kmul(vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(dir1)+cdj*(0)+cdk*(0)]),ToReal(18)))))))))
#else
#  define PDupwindNth1(u) (PDupwindNth1_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth2(u) (kmul(p1o12dy,kmul(dir2,kadd(vec_loadu_maybe3(0,3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(3*dir2)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-10),kmadd(vec_loadu_maybe3(0,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2*dir2)+cdk*(0)]),ToReal(-6),kmadd(vec_loadu_maybe3(0,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-dir2)+cdk*(0)]),ToReal(-3),kmul(vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(dir2)+cdk*(0)]),ToReal(18)))))))))
#else
#  define PDupwindNth2(u) (PDupwindNth2_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth3(u) (kmul(p1o12dz,kmul(dir3,kadd(vec_loadu_maybe3(0,0,3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(3*dir3)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-10),kmadd(vec_loadu_maybe3(0,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2*dir3)]),ToReal(-6),kmadd(vec_loadu_maybe3(0,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-dir3)]),ToReal(-3),kmul(vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(dir3)]),ToReal(18)))))))))
#else
#  define PDupwindNth3(u) (PDupwindNth3_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided1(u) (kmul(p1odx,kmul(dir1,ksub(vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(dir1)+cdj*(0)+cdk*(0)]),vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)])))))
#else
#  define PDonesided1(u) (PDonesided1_impl(u,p1odx,cdj,cdk))
static CCTK_REAL_VEC PDonesided1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesided1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided2(u) (kmul(p1ody,kmul(dir2,ksub(vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(dir2)+cdk*(0)]),vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)])))))
#else
#  define PDonesided2(u) (PDonesided2_impl(u,p1ody,cdj,cdk))
static CCTK_REAL_VEC PDonesided2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesided2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided3(u) (kmul(p1odz,kmul(dir3,ksub(vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(dir3)]),vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)])))))
#else
#  define PDonesided3(u) (PDonesided3_impl(u,p1odz,cdj,cdk))
static CCTK_REAL_VEC PDonesided3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesided3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti1(u) (kmul(p1o24dx,kadd(vec_loadu_maybe3(3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(3)+cdj*(0)+cdk*(0)]),kmadd(vec_loadu_maybe3(-1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]),ToReal(-21),kmadd(vec_loadu_maybe3(2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]),ToReal(-6),ksub(kmadd(vec_loadu_maybe3(-2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]),ToReal(6),kmul(vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]),ToReal(21))),vec_loadu_maybe3(-3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-3)+cdj*(0)+cdk*(0)])))))))
#else
#  define PDupwindNthAnti1(u) (PDupwindNthAnti1_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAnti1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAnti1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dx,kadd(vec_loadu_maybe3(3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(3)+cdj*(0)+cdk*(0)]),kmadd(vec_loadu_maybe3(-1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]),ToReal(-21),kmadd(vec_loadu_maybe3(2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]),ToReal(-6),ksub(kmadd(vec_loadu_maybe3(-2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]),ToReal(6),kmul(vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]),ToReal(21))),vec_loadu_maybe3(-3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-3)+cdj*(0)+cdk*(0)]))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm1(u) (kmul(p1o24dx,kadd(vec_loadu_maybe3(-3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-3)+cdj*(0)+cdk*(0)]),kadd(vec_loadu_maybe3(3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(3)+cdj*(0)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(-2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]),vec_loadu_maybe3(2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(-1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]),vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])),ToReal(15))))))))
#else
#  define PDupwindNthSymm1(u) (PDupwindNthSymm1_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymm1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymm1_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dx,kadd(vec_loadu_maybe3(-3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-3)+cdj*(0)+cdk*(0)]),kadd(vec_loadu_maybe3(3,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(3)+cdj*(0)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(-2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]),vec_loadu_maybe3(2,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(-1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]),vec_loadu_maybe3(1,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti2(u) (kmul(p1o24dy,kadd(vec_loadu_maybe3(0,3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(3)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]),ToReal(-21),kmadd(vec_loadu_maybe3(0,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]),ToReal(-6),ksub(kmadd(vec_loadu_maybe3(0,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]),ToReal(6),kmul(vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]),ToReal(21))),vec_loadu_maybe3(0,-3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-3)+cdk*(0)])))))))
#else
#  define PDupwindNthAnti2(u) (PDupwindNthAnti2_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAnti2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAnti2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dy,kadd(vec_loadu_maybe3(0,3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(3)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]),ToReal(-21),kmadd(vec_loadu_maybe3(0,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]),ToReal(-6),ksub(kmadd(vec_loadu_maybe3(0,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]),ToReal(6),kmul(vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]),ToReal(21))),vec_loadu_maybe3(0,-3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-3)+cdk*(0)]))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm2(u) (kmul(p1o24dy,kadd(vec_loadu_maybe3(0,-3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-3)+cdk*(0)]),kadd(vec_loadu_maybe3(0,3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(3)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(0,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]),vec_loadu_maybe3(0,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(0,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])),ToReal(15))))))))
#else
#  define PDupwindNthSymm2(u) (PDupwindNthSymm2_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymm2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymm2_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dy,kadd(vec_loadu_maybe3(0,-3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-3)+cdk*(0)]),kadd(vec_loadu_maybe3(0,3,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(3)+cdk*(0)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(0,-2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]),vec_loadu_maybe3(0,2,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(0,-1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]),vec_loadu_maybe3(0,1,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti3(u) (kmul(p1o24dz,kadd(vec_loadu_maybe3(0,0,3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(3)]),kmadd(vec_loadu_maybe3(0,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]),ToReal(-21),kmadd(vec_loadu_maybe3(0,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]),ToReal(-6),ksub(kmadd(vec_loadu_maybe3(0,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]),ToReal(6),kmul(vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]),ToReal(21))),vec_loadu_maybe3(0,0,-3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-3)])))))))
#else
#  define PDupwindNthAnti3(u) (PDupwindNthAnti3_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAnti3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAnti3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dz,kadd(vec_loadu_maybe3(0,0,3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(3)]),kmadd(vec_loadu_maybe3(0,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]),ToReal(-21),kmadd(vec_loadu_maybe3(0,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]),ToReal(-6),ksub(kmadd(vec_loadu_maybe3(0,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]),ToReal(6),kmul(vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]),ToReal(21))),vec_loadu_maybe3(0,0,-3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-3)]))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm3(u) (kmul(p1o24dz,kadd(vec_loadu_maybe3(0,0,-3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-3)]),kadd(vec_loadu_maybe3(0,0,3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(3)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(0,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]),vec_loadu_maybe3(0,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(0,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])),ToReal(15))))))))
#else
#  define PDupwindNthSymm3(u) (PDupwindNthSymm3_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymm3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymm3_impl(CCTK_REAL const* restrict const u, CCTK_REAL_VEC const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return kmul(p1o24dz,kadd(vec_loadu_maybe3(0,0,-3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-3)]),kadd(vec_loadu_maybe3(0,0,3,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(3)]),kmadd(vec_loadu_maybe3(0,0,0,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]),ToReal(-20),kmadd(kadd(vec_loadu_maybe3(0,0,-2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]),vec_loadu_maybe3(0,0,2,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)])),ToReal(-6),kmul(kadd(vec_loadu_maybe3(0,0,-1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]),vec_loadu_maybe3(0,0,1,*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])),ToReal(15)))))));
}
#endif

