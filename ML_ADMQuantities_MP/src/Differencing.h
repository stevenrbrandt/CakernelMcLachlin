#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth1(u) ((-8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + 8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]))*p1o12dx)
#else
#  define PDstandardNth1(u) (PDstandardNth1_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL PDstandardNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + 8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]))*p1o12dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth2(u) ((-8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + 8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]))*p1o12dy)
#else
#  define PDstandardNth2(u) (PDstandardNth2_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL PDstandardNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + 8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]))*p1o12dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth3(u) ((-8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + 8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]))*p1o12dz)
#else
#  define PDstandardNth3(u) (PDstandardNth3_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL PDstandardNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + 8*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]))*p1o12dz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth11(u) ((30*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 16*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]))*pm1o12dx2)
#else
#  define PDstandardNth11(u) (PDstandardNth11_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL PDstandardNth11_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dx2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth11_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dx2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (30*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 16*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]))*pm1o12dx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth22(u) ((30*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 16*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]))*pm1o12dy2)
#else
#  define PDstandardNth22(u) (PDstandardNth22_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL PDstandardNth22_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dy2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth22_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dy2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (30*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 16*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]))*pm1o12dy2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth33(u) ((30*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 16*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]))*pm1o12dz2)
#else
#  define PDstandardNth33(u) (PDstandardNth33_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL PDstandardNth33_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dz2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth33_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o12dz2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (30*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 16*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]))*pm1o12dz2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth12(u) ((-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-1)+cdk*(0)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(1)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(2)+cdk*(0)]))*p1o144dxdy)
#else
#  define PDstandardNth12(u) (PDstandardNth12_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandardNth12_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth12_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-1)+cdk*(0)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(1)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(2)+cdk*(0)]))*p1o144dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth13(u) ((-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-1)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(2)]))*p1o144dxdz)
#else
#  define PDstandardNth13(u) (PDstandardNth13_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandardNth13_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth13_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-1)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(2)]))*p1o144dxdz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth21(u) ((-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-1)+cdk*(0)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(1)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(2)+cdk*(0)]))*p1o144dxdy)
#else
#  define PDstandardNth21(u) (PDstandardNth21_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandardNth21_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth21_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-1)+cdk*(0)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(1)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(-2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(2)+cdk*(0)]))*p1o144dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth23(u) ((-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-1)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(2)]))*p1o144dydz)
#else
#  define PDstandardNth23(u) (PDstandardNth23_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandardNth23_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth23_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-1)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(2)]))*p1o144dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth31(u) ((-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-1)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(2)]))*p1o144dxdz)
#else
#  define PDstandardNth31(u) (PDstandardNth31_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandardNth31_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth31_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-1)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(2)]))*p1o144dxdz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth32(u) ((-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-1)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(2)]))*p1o144dydz)
#else
#  define PDstandardNth32(u) (PDstandardNth32_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandardNth32_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth32_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o144dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)])) + 64*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)])) + 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-1)])) - 8*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(2)]))*p1o144dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth1(u) ((-10*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(2*dir1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(3*dir1)+cdj*(0)+cdk*(0)]) - 3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(-dir1)+cdj*(0)+cdk*(0)]) + 18*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(dir1)+cdj*(0)+cdk*(0)]))*p1o12dx*dir1)
#else
#  define PDupwindNth1(u) (PDupwindNth1_impl(u,p1o12dx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dx, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-10*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(2*dir1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(3*dir1)+cdj*(0)+cdk*(0)]) - 3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(-dir1)+cdj*(0)+cdk*(0)]) + 18*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(dir1)+cdj*(0)+cdk*(0)]))*p1o12dx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth2(u) ((-10*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2*dir2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(3*dir2)+cdk*(0)]) - 3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-dir2)+cdk*(0)]) + 18*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(dir2)+cdk*(0)]))*p1o12dy*dir2)
#else
#  define PDupwindNth2(u) (PDupwindNth2_impl(u,p1o12dy,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dy, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-10*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2*dir2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(3*dir2)+cdk*(0)]) - 3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-dir2)+cdk*(0)]) + 18*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(dir2)+cdk*(0)]))*p1o12dy*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth3(u) ((-10*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2*dir3)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(3*dir3)]) - 3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-dir3)]) + 18*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(dir3)]))*p1o12dz*dir3)
#else
#  define PDupwindNth3(u) (PDupwindNth3_impl(u,p1o12dz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o12dz, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-10*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2*dir3)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(3*dir3)]) - 3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-dir3)]) + 18*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(dir3)]))*p1o12dz*dir3;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided1(u) ((-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(dir1)+cdj*(0)+cdk*(0)]))*p1odx*dir1)
#else
#  define PDonesided1(u) (PDonesided1_impl(u,p1odx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odx, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(dir1)+cdj*(0)+cdk*(0)]))*p1odx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided2(u) ((-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(dir2)+cdk*(0)]))*p1ody*dir2)
#else
#  define PDonesided2(u) (PDonesided2_impl(u,p1ody,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1ody, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(dir2)+cdk*(0)]))*p1ody*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided3(u) ((-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(dir3)]))*p1odz*dir3)
#else
#  define PDonesided3(u) (PDonesided3_impl(u,p1odz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odz, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(dir3)]))*p1odz*dir3;
}
#endif

