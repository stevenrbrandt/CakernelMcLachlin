#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth1(u) ((-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]))*p1o2dx)
#else
#  define PDstandardNth1(u) (PDstandardNth1_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL PDstandardNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]))*p1o2dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth2(u) ((-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]))*p1o2dy)
#else
#  define PDstandardNth2(u) (PDstandardNth2_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL PDstandardNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]))*p1o2dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth3(u) ((-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]))*p1o2dz)
#else
#  define PDstandardNth3(u) (PDstandardNth3_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL PDstandardNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]))*p1o2dz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth11(u) ((-2*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]))*p1odx2)
#else
#  define PDstandardNth11(u) (PDstandardNth11_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL PDstandardNth11_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odx2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth11_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odx2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-2*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]))*p1odx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth22(u) ((-2*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]))*p1ody2)
#else
#  define PDstandardNth22(u) (PDstandardNth22_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL PDstandardNth22_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1ody2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth22_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1ody2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-2*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]))*p1ody2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth33(u) ((-2*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]))*p1odz2)
#else
#  define PDstandardNth33(u) (PDstandardNth33_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL PDstandardNth33_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odz2, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth33_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1odz2, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-2*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]))*p1odz2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth12(u) (((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)]))*p1o4dxdy)
#else
#  define PDstandardNth12(u) (PDstandardNth12_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandardNth12_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth12_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return ((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)]))*p1o4dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth13(u) (((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)]))*p1o4dxdz)
#else
#  define PDstandardNth13(u) (PDstandardNth13_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandardNth13_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth13_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return ((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)]))*p1o4dxdz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth21(u) (((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)]))*p1o4dxdy)
#else
#  define PDstandardNth21(u) (PDstandardNth21_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandardNth21_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth21_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dxdy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return ((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(-1)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(1)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(1)+cdk*(0)]))*p1o4dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth23(u) (((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)]))*p1o4dydz)
#else
#  define PDstandardNth23(u) (PDstandardNth23_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandardNth23_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth23_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return ((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)]))*p1o4dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth31(u) (((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)]))*p1o4dxdz)
#else
#  define PDstandardNth31(u) (PDstandardNth31_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandardNth31_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth31_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dxdz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return ((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(-1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(1)]))*p1o4dxdz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth32(u) (((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)]))*p1o4dydz)
#else
#  define PDstandardNth32(u) (PDstandardNth32_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandardNth32_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dydz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth32_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dydz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return ((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(-1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(1)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(1)]))*p1o4dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth1(u) ((6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]))*p1o16dx)
#else
#  define PDdissipationNth1(u) (PDdissipationNth1_impl(u,p1o16dx,cdj,cdk))
static CCTK_REAL PDdissipationNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o16dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o16dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]))*p1o16dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth2(u) ((6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]))*p1o16dy)
#else
#  define PDdissipationNth2(u) (PDdissipationNth2_impl(u,p1o16dy,cdj,cdk))
static CCTK_REAL PDdissipationNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o16dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o16dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]))*p1o16dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth3(u) ((6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]))*p1o16dz)
#else
#  define PDdissipationNth3(u) (PDdissipationNth3_impl(u,p1o16dz,cdj,cdk))
static CCTK_REAL PDdissipationNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o16dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o16dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]))*p1o16dz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth1(u) ((3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2*dir1)+cdj*(0)+cdk*(0)]) - 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(dir1)+cdj*(0)+cdk*(0)]))*pm1o2dx*dir1)
#else
#  define PDupwindNth1(u) (PDupwindNth1_impl(u,pm1o2dx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o2dx, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2*dir1)+cdj*(0)+cdk*(0)]) - 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(dir1)+cdj*(0)+cdk*(0)]))*pm1o2dx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth2(u) ((3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2*dir2)+cdk*(0)]) - 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(dir2)+cdk*(0)]))*pm1o2dy*dir2)
#else
#  define PDupwindNth2(u) (PDupwindNth2_impl(u,pm1o2dy,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o2dy, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2*dir2)+cdk*(0)]) - 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(dir2)+cdk*(0)]))*pm1o2dy*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth3(u) ((3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2*dir3)]) - 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(dir3)]))*pm1o2dz*dir3)
#else
#  define PDupwindNth3(u) (PDupwindNth3_impl(u,pm1o2dz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o2dz, ptrdiff_t const cdj, ptrdiff_t const cdk, ptrdiff_t const dir1, ptrdiff_t const dir2, ptrdiff_t const dir3)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (3*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2*dir3)]) - 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(dir3)]))*pm1o2dz*dir3;
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

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti1(u) ((-4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]))*p1o4dx)
#else
#  define PDupwindNthAnti1(u) (PDupwindNthAnti1_impl(u,p1o4dx,cdj,cdk))
static CCTK_REAL PDupwindNthAnti1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAnti1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]))*p1o4dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm1(u) ((6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]))*pm1o4dx)
#else
#  define PDupwindNthSymm1(u) (PDupwindNthSymm1_impl(u,pm1o4dx,cdj,cdk))
static CCTK_REAL PDupwindNthSymm1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o4dx, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymm1_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o4dx, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(-1)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(1)+cdj*(0)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(-2)+cdj*(0)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(2)+cdj*(0)+cdk*(0)]))*pm1o4dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti2(u) ((-4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]))*p1o4dy)
#else
#  define PDupwindNthAnti2(u) (PDupwindNthAnti2_impl(u,p1o4dy,cdj,cdk))
static CCTK_REAL PDupwindNthAnti2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAnti2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]))*p1o4dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm2(u) ((6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]))*pm1o4dy)
#else
#  define PDupwindNthSymm2(u) (PDupwindNthSymm2_impl(u,pm1o4dy,cdj,cdk))
static CCTK_REAL PDupwindNthSymm2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o4dy, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymm2_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o4dy, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-1)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(1)+cdk*(0)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(-2)+cdk*(0)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(2)+cdk*(0)]))*pm1o4dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti3(u) ((-4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]))*p1o4dz)
#else
#  define PDupwindNthAnti3(u) (PDupwindNthAnti3_impl(u,p1o4dz,cdj,cdk))
static CCTK_REAL PDupwindNthAnti3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAnti3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const p1o4dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (-4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + 4*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]) - (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]))*p1o4dz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm3(u) ((6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]))*pm1o4dz)
#else
#  define PDupwindNthSymm3(u) (PDupwindNthSymm3_impl(u,pm1o4dz,cdj,cdk))
static CCTK_REAL PDupwindNthSymm3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o4dz, ptrdiff_t const cdj, ptrdiff_t const cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymm3_impl(CCTK_REAL const* restrict const u, CCTK_REAL const pm1o4dz, ptrdiff_t const cdj, ptrdiff_t const cdk)
{
  ptrdiff_t const cdi=sizeof(CCTK_REAL);
  return (6*(*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(0)]) - 4*((*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-1)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(1)])) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(-2)]) + (*(CCTK_REAL const*)&((char const*)(u))[cdi*(0)+cdj*(0)+cdk*(2)]))*pm1o4dz;
}
#endif

