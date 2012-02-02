#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth1(u) ((-8*KRANC_GFOFFSET3D(u,-1,0,0) + 8*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o12dx)
#else
#  define PDstandardNth1(u) (PDstandardNth1_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL PDstandardNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-8*KRANC_GFOFFSET3D(u,-1,0,0) + 8*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o12dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth2(u) ((-8*KRANC_GFOFFSET3D(u,0,-1,0) + 8*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o12dy)
#else
#  define PDstandardNth2(u) (PDstandardNth2_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL PDstandardNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-8*KRANC_GFOFFSET3D(u,0,-1,0) + 8*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o12dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth3(u) ((-8*KRANC_GFOFFSET3D(u,0,0,-1) + 8*KRANC_GFOFFSET3D(u,0,0,1) + KRANC_GFOFFSET3D(u,0,0,-2) - KRANC_GFOFFSET3D(u,0,0,2))*p1o12dz)
#else
#  define PDstandardNth3(u) (PDstandardNth3_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL PDstandardNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-8*KRANC_GFOFFSET3D(u,0,0,-1) + 8*KRANC_GFOFFSET3D(u,0,0,1) + KRANC_GFOFFSET3D(u,0,0,-2) - KRANC_GFOFFSET3D(u,0,0,2))*p1o12dz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth11(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o12dx2)
#else
#  define PDstandardNth11(u) (PDstandardNth11_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL PDstandardNth11_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dx2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth11_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dx2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o12dx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth22(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o12dy2)
#else
#  define PDstandardNth22(u) (PDstandardNth22_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL PDstandardNth22_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dy2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth22_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dy2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o12dy2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth33(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) + KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2))*pm1o12dz2)
#else
#  define PDstandardNth33(u) (PDstandardNth33_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL PDstandardNth33_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dz2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth33_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dz2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) + KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2))*pm1o12dz2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth12(u) ((-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy)
#else
#  define PDstandardNth12(u) (PDstandardNth12_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandardNth12_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth12_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth13(u) ((-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz)
#else
#  define PDstandardNth13(u) (PDstandardNth13_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandardNth13_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth13_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth21(u) ((-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy)
#else
#  define PDstandardNth21(u) (PDstandardNth21_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandardNth21_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth21_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth23(u) ((-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz)
#else
#  define PDstandardNth23(u) (PDstandardNth23_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandardNth23_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth23_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth31(u) ((-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz)
#else
#  define PDstandardNth31(u) (PDstandardNth31_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandardNth31_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth31_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth32(u) ((-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz)
#else
#  define PDstandardNth32(u) (PDstandardNth32_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandardNth32_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth32_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth1(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 6*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o64dx)
#else
#  define PDdissipationNth1(u) (PDdissipationNth1_impl(u,p1o64dx,cdj,cdk))
static CCTK_REAL PDdissipationNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o64dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o64dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 6*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o64dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth2(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 6*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o64dy)
#else
#  define PDdissipationNth2(u) (PDdissipationNth2_impl(u,p1o64dy,cdj,cdk))
static CCTK_REAL PDdissipationNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o64dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o64dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 6*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o64dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth3(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 6*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3))*p1o64dz)
#else
#  define PDdissipationNth3(u) (PDdissipationNth3_impl(u,p1o64dz,cdj,cdk))
static CCTK_REAL PDdissipationNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o64dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o64dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 6*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3))*p1o64dz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth1(u) ((-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,2*dir1,0,0) + KRANC_GFOFFSET3D(u,3*dir1,0,0) - 3*KRANC_GFOFFSET3D(u,-dir1,0,0) + 18*KRANC_GFOFFSET3D(u,dir1,0,0))*p1o12dx*dir1)
#else
#  define PDupwindNth1(u) (PDupwindNth1_impl(u,p1o12dx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,2*dir1,0,0) + KRANC_GFOFFSET3D(u,3*dir1,0,0) - 3*KRANC_GFOFFSET3D(u,-dir1,0,0) + 18*KRANC_GFOFFSET3D(u,dir1,0,0))*p1o12dx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti1(u) ((-21*KRANC_GFOFFSET3D(u,-1,0,0) + 21*KRANC_GFOFFSET3D(u,1,0,0) + 6*KRANC_GFOFFSET3D(u,-2,0,0) - 6*KRANC_GFOFFSET3D(u,2,0,0) - KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o24dx)
#else
#  define PDupwindNthAnti1(u) (PDupwindNthAnti1_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL PDupwindNthAnti1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAnti1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-21*KRANC_GFOFFSET3D(u,-1,0,0) + 21*KRANC_GFOFFSET3D(u,1,0,0) + 6*KRANC_GFOFFSET3D(u,-2,0,0) - 6*KRANC_GFOFFSET3D(u,2,0,0) - KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o24dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm1(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 6*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o24dx)
#else
#  define PDupwindNthSymm1(u) (PDupwindNthSymm1_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL PDupwindNthSymm1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymm1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 6*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o24dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided1(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,dir1,0,0))*p1odx*dir1)
#else
#  define PDonesided1(u) (PDonesided1_impl(u,p1odx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,dir1,0,0))*p1odx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth2(u) ((-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,0,2*dir2,0) + KRANC_GFOFFSET3D(u,0,3*dir2,0) - 3*KRANC_GFOFFSET3D(u,0,-dir2,0) + 18*KRANC_GFOFFSET3D(u,0,dir2,0))*p1o12dy*dir2)
#else
#  define PDupwindNth2(u) (PDupwindNth2_impl(u,p1o12dy,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,0,2*dir2,0) + KRANC_GFOFFSET3D(u,0,3*dir2,0) - 3*KRANC_GFOFFSET3D(u,0,-dir2,0) + 18*KRANC_GFOFFSET3D(u,0,dir2,0))*p1o12dy*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti2(u) ((-21*KRANC_GFOFFSET3D(u,0,-1,0) + 21*KRANC_GFOFFSET3D(u,0,1,0) + 6*KRANC_GFOFFSET3D(u,0,-2,0) - 6*KRANC_GFOFFSET3D(u,0,2,0) - KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o24dy)
#else
#  define PDupwindNthAnti2(u) (PDupwindNthAnti2_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL PDupwindNthAnti2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAnti2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-21*KRANC_GFOFFSET3D(u,0,-1,0) + 21*KRANC_GFOFFSET3D(u,0,1,0) + 6*KRANC_GFOFFSET3D(u,0,-2,0) - 6*KRANC_GFOFFSET3D(u,0,2,0) - KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o24dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm2(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 6*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o24dy)
#else
#  define PDupwindNthSymm2(u) (PDupwindNthSymm2_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL PDupwindNthSymm2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymm2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 6*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o24dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided2(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,dir2,0))*p1ody*dir2)
#else
#  define PDonesided2(u) (PDonesided2_impl(u,p1ody,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,dir2,0))*p1ody*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth3(u) ((-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,0,0,2*dir3) + KRANC_GFOFFSET3D(u,0,0,3*dir3) - 3*KRANC_GFOFFSET3D(u,0,0,-dir3) + 18*KRANC_GFOFFSET3D(u,0,0,dir3))*p1o12dz*dir3)
#else
#  define PDupwindNth3(u) (PDupwindNth3_impl(u,p1o12dz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,0,0,2*dir3) + KRANC_GFOFFSET3D(u,0,0,3*dir3) - 3*KRANC_GFOFFSET3D(u,0,0,-dir3) + 18*KRANC_GFOFFSET3D(u,0,0,dir3))*p1o12dz*dir3;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti3(u) ((-21*KRANC_GFOFFSET3D(u,0,0,-1) + 21*KRANC_GFOFFSET3D(u,0,0,1) + 6*KRANC_GFOFFSET3D(u,0,0,-2) - 6*KRANC_GFOFFSET3D(u,0,0,2) - KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3))*p1o24dz)
#else
#  define PDupwindNthAnti3(u) (PDupwindNthAnti3_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL PDupwindNthAnti3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAnti3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-21*KRANC_GFOFFSET3D(u,0,0,-1) + 21*KRANC_GFOFFSET3D(u,0,0,1) + 6*KRANC_GFOFFSET3D(u,0,0,-2) - 6*KRANC_GFOFFSET3D(u,0,0,2) - KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3))*p1o24dz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm3(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 6*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3))*p1o24dz)
#else
#  define PDupwindNthSymm3(u) (PDupwindNthSymm3_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL PDupwindNthSymm3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymm3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o24dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 6*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3))*p1o24dz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided3(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,dir3))*p1odz*dir3)
#else
#  define PDonesided3(u) (PDonesided3_impl(u,p1odz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,dir3))*p1odz*dir3;
}
#endif

