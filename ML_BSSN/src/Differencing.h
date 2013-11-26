#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth1(u) ((-672*KRANC_GFOFFSET3D(u,-1,0,0) + 672*KRANC_GFOFFSET3D(u,1,0,0) + 168*KRANC_GFOFFSET3D(u,-2,0,0) - 168*KRANC_GFOFFSET3D(u,2,0,0) - 32*KRANC_GFOFFSET3D(u,-3,0,0) + 32*KRANC_GFOFFSET3D(u,3,0,0) + 3*KRANC_GFOFFSET3D(u,-4,0,0) - 3*KRANC_GFOFFSET3D(u,4,0,0))*p1o840dx)
#else
#  define PDstandardNth1(u) (PDstandardNth1_impl(u,p1o840dx,cdj,cdk))
static CCTK_REAL PDstandardNth1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-672*KRANC_GFOFFSET3D(u,-1,0,0) + 672*KRANC_GFOFFSET3D(u,1,0,0) + 168*KRANC_GFOFFSET3D(u,-2,0,0) - 168*KRANC_GFOFFSET3D(u,2,0,0) - 32*KRANC_GFOFFSET3D(u,-3,0,0) + 32*KRANC_GFOFFSET3D(u,3,0,0) + 3*KRANC_GFOFFSET3D(u,-4,0,0) - 3*KRANC_GFOFFSET3D(u,4,0,0))*p1o840dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth2(u) ((-672*KRANC_GFOFFSET3D(u,0,-1,0) + 672*KRANC_GFOFFSET3D(u,0,1,0) + 168*KRANC_GFOFFSET3D(u,0,-2,0) - 168*KRANC_GFOFFSET3D(u,0,2,0) - 32*KRANC_GFOFFSET3D(u,0,-3,0) + 32*KRANC_GFOFFSET3D(u,0,3,0) + 3*KRANC_GFOFFSET3D(u,0,-4,0) - 3*KRANC_GFOFFSET3D(u,0,4,0))*p1o840dy)
#else
#  define PDstandardNth2(u) (PDstandardNth2_impl(u,p1o840dy,cdj,cdk))
static CCTK_REAL PDstandardNth2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-672*KRANC_GFOFFSET3D(u,0,-1,0) + 672*KRANC_GFOFFSET3D(u,0,1,0) + 168*KRANC_GFOFFSET3D(u,0,-2,0) - 168*KRANC_GFOFFSET3D(u,0,2,0) - 32*KRANC_GFOFFSET3D(u,0,-3,0) + 32*KRANC_GFOFFSET3D(u,0,3,0) + 3*KRANC_GFOFFSET3D(u,0,-4,0) - 3*KRANC_GFOFFSET3D(u,0,4,0))*p1o840dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth3(u) ((-672*KRANC_GFOFFSET3D(u,0,0,-1) + 672*KRANC_GFOFFSET3D(u,0,0,1) + 168*KRANC_GFOFFSET3D(u,0,0,-2) - 168*KRANC_GFOFFSET3D(u,0,0,2) - 32*KRANC_GFOFFSET3D(u,0,0,-3) + 32*KRANC_GFOFFSET3D(u,0,0,3) + 3*KRANC_GFOFFSET3D(u,0,0,-4) - 3*KRANC_GFOFFSET3D(u,0,0,4))*p1o840dz)
#else
#  define PDstandardNth3(u) (PDstandardNth3_impl(u,p1o840dz,cdj,cdk))
static CCTK_REAL PDstandardNth3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth2_impl(u, p1o840dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth11(u) ((-14350*KRANC_GFOFFSET3D(u,0,0,0) + 8064*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 1008*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 128*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)) - 9*(KRANC_GFOFFSET3D(u,-4,0,0) + KRANC_GFOFFSET3D(u,4,0,0)))*p1o5040dx2)
#else
#  define PDstandardNth11(u) (PDstandardNth11_impl(u,p1o5040dx2,cdj,cdk))
static CCTK_REAL PDstandardNth11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-14350*KRANC_GFOFFSET3D(u,0,0,0) + 8064*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 1008*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 128*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)) - 9*(KRANC_GFOFFSET3D(u,-4,0,0) + KRANC_GFOFFSET3D(u,4,0,0)))*p1o5040dx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth22(u) ((-14350*KRANC_GFOFFSET3D(u,0,0,0) + 8064*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 1008*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 128*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)) - 9*(KRANC_GFOFFSET3D(u,0,-4,0) + KRANC_GFOFFSET3D(u,0,4,0)))*p1o5040dy2)
#else
#  define PDstandardNth22(u) (PDstandardNth22_impl(u,p1o5040dy2,cdj,cdk))
static CCTK_REAL PDstandardNth22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-14350*KRANC_GFOFFSET3D(u,0,0,0) + 8064*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 1008*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 128*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)) - 9*(KRANC_GFOFFSET3D(u,0,-4,0) + KRANC_GFOFFSET3D(u,0,4,0)))*p1o5040dy2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth33(u) ((-14350*KRANC_GFOFFSET3D(u,0,0,0) + 8064*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 1008*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + 128*(KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3)) - 9*(KRANC_GFOFFSET3D(u,0,0,-4) + KRANC_GFOFFSET3D(u,0,0,4)))*p1o5040dz2)
#else
#  define PDstandardNth33(u) (PDstandardNth33_impl(u,p1o5040dz2,cdj,cdk))
static CCTK_REAL PDstandardNth33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o5040dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth22_impl(u, p1o5040dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth12(u) ((-451584*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 451584*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 112896*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 112896*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 28224*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 28224*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 21504*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 21504*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 5376*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 5376*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) - 1024*(KRANC_GFOFFSET3D(u,-3,3,0) + KRANC_GFOFFSET3D(u,3,-3,0)) + 1024*(KRANC_GFOFFSET3D(u,-3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0)) + 2016*(KRANC_GFOFFSET3D(u,-1,4,0) + KRANC_GFOFFSET3D(u,1,-4,0) + KRANC_GFOFFSET3D(u,-4,1,0) + KRANC_GFOFFSET3D(u,4,-1,0)) - 2016*(KRANC_GFOFFSET3D(u,-1,-4,0) + KRANC_GFOFFSET3D(u,1,4,0) + KRANC_GFOFFSET3D(u,-4,-1,0) + KRANC_GFOFFSET3D(u,4,1,0)) - 504*(KRANC_GFOFFSET3D(u,-2,4,0) + KRANC_GFOFFSET3D(u,2,-4,0) + KRANC_GFOFFSET3D(u,-4,2,0) + KRANC_GFOFFSET3D(u,4,-2,0)) + 504*(KRANC_GFOFFSET3D(u,-2,-4,0) + KRANC_GFOFFSET3D(u,2,4,0) + KRANC_GFOFFSET3D(u,-4,-2,0) + KRANC_GFOFFSET3D(u,4,2,0)) + 96*(KRANC_GFOFFSET3D(u,-3,4,0) + KRANC_GFOFFSET3D(u,3,-4,0) + KRANC_GFOFFSET3D(u,-4,3,0) + KRANC_GFOFFSET3D(u,4,-3,0)) - 96*(KRANC_GFOFFSET3D(u,-3,-4,0) + KRANC_GFOFFSET3D(u,3,4,0) + KRANC_GFOFFSET3D(u,-4,-3,0) + KRANC_GFOFFSET3D(u,4,3,0)) - 9*(KRANC_GFOFFSET3D(u,-4,4,0) + KRANC_GFOFFSET3D(u,4,-4,0)) + 9*(KRANC_GFOFFSET3D(u,-4,-4,0) + KRANC_GFOFFSET3D(u,4,4,0)))*p1o705600dxdy)
#else
#  define PDstandardNth12(u) (PDstandardNth12_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL PDstandardNth12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-451584*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 451584*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 112896*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 112896*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 28224*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 28224*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 21504*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 21504*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 5376*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 5376*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) - 1024*(KRANC_GFOFFSET3D(u,-3,3,0) + KRANC_GFOFFSET3D(u,3,-3,0)) + 1024*(KRANC_GFOFFSET3D(u,-3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0)) + 2016*(KRANC_GFOFFSET3D(u,-1,4,0) + KRANC_GFOFFSET3D(u,1,-4,0) + KRANC_GFOFFSET3D(u,-4,1,0) + KRANC_GFOFFSET3D(u,4,-1,0)) - 2016*(KRANC_GFOFFSET3D(u,-1,-4,0) + KRANC_GFOFFSET3D(u,1,4,0) + KRANC_GFOFFSET3D(u,-4,-1,0) + KRANC_GFOFFSET3D(u,4,1,0)) - 504*(KRANC_GFOFFSET3D(u,-2,4,0) + KRANC_GFOFFSET3D(u,2,-4,0) + KRANC_GFOFFSET3D(u,-4,2,0) + KRANC_GFOFFSET3D(u,4,-2,0)) + 504*(KRANC_GFOFFSET3D(u,-2,-4,0) + KRANC_GFOFFSET3D(u,2,4,0) + KRANC_GFOFFSET3D(u,-4,-2,0) + KRANC_GFOFFSET3D(u,4,2,0)) + 96*(KRANC_GFOFFSET3D(u,-3,4,0) + KRANC_GFOFFSET3D(u,3,-4,0) + KRANC_GFOFFSET3D(u,-4,3,0) + KRANC_GFOFFSET3D(u,4,-3,0)) - 96*(KRANC_GFOFFSET3D(u,-3,-4,0) + KRANC_GFOFFSET3D(u,3,4,0) + KRANC_GFOFFSET3D(u,-4,-3,0) + KRANC_GFOFFSET3D(u,4,3,0)) - 9*(KRANC_GFOFFSET3D(u,-4,4,0) + KRANC_GFOFFSET3D(u,4,-4,0)) + 9*(KRANC_GFOFFSET3D(u,-4,-4,0) + KRANC_GFOFFSET3D(u,4,4,0)))*p1o705600dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth13(u) ((-451584*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 451584*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 112896*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 112896*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) - 28224*(KRANC_GFOFFSET3D(u,-2,0,2) + KRANC_GFOFFSET3D(u,2,0,-2)) + 28224*(KRANC_GFOFFSET3D(u,-2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2)) - 21504*(KRANC_GFOFFSET3D(u,-1,0,3) + KRANC_GFOFFSET3D(u,1,0,-3) + KRANC_GFOFFSET3D(u,-3,0,1) + KRANC_GFOFFSET3D(u,3,0,-1)) + 21504*(KRANC_GFOFFSET3D(u,-1,0,-3) + KRANC_GFOFFSET3D(u,1,0,3) + KRANC_GFOFFSET3D(u,-3,0,-1) + KRANC_GFOFFSET3D(u,3,0,1)) + 5376*(KRANC_GFOFFSET3D(u,-2,0,3) + KRANC_GFOFFSET3D(u,2,0,-3) + KRANC_GFOFFSET3D(u,-3,0,2) + KRANC_GFOFFSET3D(u,3,0,-2)) - 5376*(KRANC_GFOFFSET3D(u,-2,0,-3) + KRANC_GFOFFSET3D(u,2,0,3) + KRANC_GFOFFSET3D(u,-3,0,-2) + KRANC_GFOFFSET3D(u,3,0,2)) - 1024*(KRANC_GFOFFSET3D(u,-3,0,3) + KRANC_GFOFFSET3D(u,3,0,-3)) + 1024*(KRANC_GFOFFSET3D(u,-3,0,-3) + KRANC_GFOFFSET3D(u,3,0,3)) + 2016*(KRANC_GFOFFSET3D(u,-1,0,4) + KRANC_GFOFFSET3D(u,1,0,-4) + KRANC_GFOFFSET3D(u,-4,0,1) + KRANC_GFOFFSET3D(u,4,0,-1)) - 2016*(KRANC_GFOFFSET3D(u,-1,0,-4) + KRANC_GFOFFSET3D(u,1,0,4) + KRANC_GFOFFSET3D(u,-4,0,-1) + KRANC_GFOFFSET3D(u,4,0,1)) - 504*(KRANC_GFOFFSET3D(u,-2,0,4) + KRANC_GFOFFSET3D(u,2,0,-4) + KRANC_GFOFFSET3D(u,-4,0,2) + KRANC_GFOFFSET3D(u,4,0,-2)) + 504*(KRANC_GFOFFSET3D(u,-2,0,-4) + KRANC_GFOFFSET3D(u,2,0,4) + KRANC_GFOFFSET3D(u,-4,0,-2) + KRANC_GFOFFSET3D(u,4,0,2)) + 96*(KRANC_GFOFFSET3D(u,-3,0,4) + KRANC_GFOFFSET3D(u,3,0,-4) + KRANC_GFOFFSET3D(u,-4,0,3) + KRANC_GFOFFSET3D(u,4,0,-3)) - 96*(KRANC_GFOFFSET3D(u,-3,0,-4) + KRANC_GFOFFSET3D(u,3,0,4) + KRANC_GFOFFSET3D(u,-4,0,-3) + KRANC_GFOFFSET3D(u,4,0,3)) - 9*(KRANC_GFOFFSET3D(u,-4,0,4) + KRANC_GFOFFSET3D(u,4,0,-4)) + 9*(KRANC_GFOFFSET3D(u,-4,0,-4) + KRANC_GFOFFSET3D(u,4,0,4)))*p1o705600dxdz)
#else
#  define PDstandardNth13(u) (PDstandardNth13_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL PDstandardNth13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth12_impl(u, p1o705600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth21(u) ((-451584*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 451584*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 112896*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 112896*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 28224*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 28224*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 21504*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 21504*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 5376*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 5376*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) - 1024*(KRANC_GFOFFSET3D(u,-3,3,0) + KRANC_GFOFFSET3D(u,3,-3,0)) + 1024*(KRANC_GFOFFSET3D(u,-3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0)) + 2016*(KRANC_GFOFFSET3D(u,-1,4,0) + KRANC_GFOFFSET3D(u,1,-4,0) + KRANC_GFOFFSET3D(u,-4,1,0) + KRANC_GFOFFSET3D(u,4,-1,0)) - 2016*(KRANC_GFOFFSET3D(u,-1,-4,0) + KRANC_GFOFFSET3D(u,1,4,0) + KRANC_GFOFFSET3D(u,-4,-1,0) + KRANC_GFOFFSET3D(u,4,1,0)) - 504*(KRANC_GFOFFSET3D(u,-2,4,0) + KRANC_GFOFFSET3D(u,2,-4,0) + KRANC_GFOFFSET3D(u,-4,2,0) + KRANC_GFOFFSET3D(u,4,-2,0)) + 504*(KRANC_GFOFFSET3D(u,-2,-4,0) + KRANC_GFOFFSET3D(u,2,4,0) + KRANC_GFOFFSET3D(u,-4,-2,0) + KRANC_GFOFFSET3D(u,4,2,0)) + 96*(KRANC_GFOFFSET3D(u,-3,4,0) + KRANC_GFOFFSET3D(u,3,-4,0) + KRANC_GFOFFSET3D(u,-4,3,0) + KRANC_GFOFFSET3D(u,4,-3,0)) - 96*(KRANC_GFOFFSET3D(u,-3,-4,0) + KRANC_GFOFFSET3D(u,3,4,0) + KRANC_GFOFFSET3D(u,-4,-3,0) + KRANC_GFOFFSET3D(u,4,3,0)) - 9*(KRANC_GFOFFSET3D(u,-4,4,0) + KRANC_GFOFFSET3D(u,4,-4,0)) + 9*(KRANC_GFOFFSET3D(u,-4,-4,0) + KRANC_GFOFFSET3D(u,4,4,0)))*p1o705600dxdy)
#else
#  define PDstandardNth21(u) (PDstandardNth21_impl(u,p1o705600dxdy,cdj,cdk))
static CCTK_REAL PDstandardNth21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth12_impl(u, p1o705600dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth23(u) ((-451584*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 451584*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 112896*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 112896*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 28224*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 28224*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 21504*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 21504*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 5376*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 5376*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) - 1024*(KRANC_GFOFFSET3D(u,0,-3,3) + KRANC_GFOFFSET3D(u,0,3,-3)) + 1024*(KRANC_GFOFFSET3D(u,0,-3,-3) + KRANC_GFOFFSET3D(u,0,3,3)) + 2016*(KRANC_GFOFFSET3D(u,0,-1,4) + KRANC_GFOFFSET3D(u,0,1,-4) + KRANC_GFOFFSET3D(u,0,-4,1) + KRANC_GFOFFSET3D(u,0,4,-1)) - 2016*(KRANC_GFOFFSET3D(u,0,-1,-4) + KRANC_GFOFFSET3D(u,0,1,4) + KRANC_GFOFFSET3D(u,0,-4,-1) + KRANC_GFOFFSET3D(u,0,4,1)) - 504*(KRANC_GFOFFSET3D(u,0,-2,4) + KRANC_GFOFFSET3D(u,0,2,-4) + KRANC_GFOFFSET3D(u,0,-4,2) + KRANC_GFOFFSET3D(u,0,4,-2)) + 504*(KRANC_GFOFFSET3D(u,0,-2,-4) + KRANC_GFOFFSET3D(u,0,2,4) + KRANC_GFOFFSET3D(u,0,-4,-2) + KRANC_GFOFFSET3D(u,0,4,2)) + 96*(KRANC_GFOFFSET3D(u,0,-3,4) + KRANC_GFOFFSET3D(u,0,3,-4) + KRANC_GFOFFSET3D(u,0,-4,3) + KRANC_GFOFFSET3D(u,0,4,-3)) - 96*(KRANC_GFOFFSET3D(u,0,-3,-4) + KRANC_GFOFFSET3D(u,0,3,4) + KRANC_GFOFFSET3D(u,0,-4,-3) + KRANC_GFOFFSET3D(u,0,4,3)) - 9*(KRANC_GFOFFSET3D(u,0,-4,4) + KRANC_GFOFFSET3D(u,0,4,-4)) + 9*(KRANC_GFOFFSET3D(u,0,-4,-4) + KRANC_GFOFFSET3D(u,0,4,4)))*p1o705600dydz)
#else
#  define PDstandardNth23(u) (PDstandardNth23_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL PDstandardNth23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-451584*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 451584*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 112896*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 112896*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 28224*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 28224*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 21504*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 21504*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 5376*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 5376*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) - 1024*(KRANC_GFOFFSET3D(u,0,-3,3) + KRANC_GFOFFSET3D(u,0,3,-3)) + 1024*(KRANC_GFOFFSET3D(u,0,-3,-3) + KRANC_GFOFFSET3D(u,0,3,3)) + 2016*(KRANC_GFOFFSET3D(u,0,-1,4) + KRANC_GFOFFSET3D(u,0,1,-4) + KRANC_GFOFFSET3D(u,0,-4,1) + KRANC_GFOFFSET3D(u,0,4,-1)) - 2016*(KRANC_GFOFFSET3D(u,0,-1,-4) + KRANC_GFOFFSET3D(u,0,1,4) + KRANC_GFOFFSET3D(u,0,-4,-1) + KRANC_GFOFFSET3D(u,0,4,1)) - 504*(KRANC_GFOFFSET3D(u,0,-2,4) + KRANC_GFOFFSET3D(u,0,2,-4) + KRANC_GFOFFSET3D(u,0,-4,2) + KRANC_GFOFFSET3D(u,0,4,-2)) + 504*(KRANC_GFOFFSET3D(u,0,-2,-4) + KRANC_GFOFFSET3D(u,0,2,4) + KRANC_GFOFFSET3D(u,0,-4,-2) + KRANC_GFOFFSET3D(u,0,4,2)) + 96*(KRANC_GFOFFSET3D(u,0,-3,4) + KRANC_GFOFFSET3D(u,0,3,-4) + KRANC_GFOFFSET3D(u,0,-4,3) + KRANC_GFOFFSET3D(u,0,4,-3)) - 96*(KRANC_GFOFFSET3D(u,0,-3,-4) + KRANC_GFOFFSET3D(u,0,3,4) + KRANC_GFOFFSET3D(u,0,-4,-3) + KRANC_GFOFFSET3D(u,0,4,3)) - 9*(KRANC_GFOFFSET3D(u,0,-4,4) + KRANC_GFOFFSET3D(u,0,4,-4)) + 9*(KRANC_GFOFFSET3D(u,0,-4,-4) + KRANC_GFOFFSET3D(u,0,4,4)))*p1o705600dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth31(u) ((-451584*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 451584*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 112896*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 112896*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) - 28224*(KRANC_GFOFFSET3D(u,-2,0,2) + KRANC_GFOFFSET3D(u,2,0,-2)) + 28224*(KRANC_GFOFFSET3D(u,-2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2)) - 21504*(KRANC_GFOFFSET3D(u,-1,0,3) + KRANC_GFOFFSET3D(u,1,0,-3) + KRANC_GFOFFSET3D(u,-3,0,1) + KRANC_GFOFFSET3D(u,3,0,-1)) + 21504*(KRANC_GFOFFSET3D(u,-1,0,-3) + KRANC_GFOFFSET3D(u,1,0,3) + KRANC_GFOFFSET3D(u,-3,0,-1) + KRANC_GFOFFSET3D(u,3,0,1)) + 5376*(KRANC_GFOFFSET3D(u,-2,0,3) + KRANC_GFOFFSET3D(u,2,0,-3) + KRANC_GFOFFSET3D(u,-3,0,2) + KRANC_GFOFFSET3D(u,3,0,-2)) - 5376*(KRANC_GFOFFSET3D(u,-2,0,-3) + KRANC_GFOFFSET3D(u,2,0,3) + KRANC_GFOFFSET3D(u,-3,0,-2) + KRANC_GFOFFSET3D(u,3,0,2)) - 1024*(KRANC_GFOFFSET3D(u,-3,0,3) + KRANC_GFOFFSET3D(u,3,0,-3)) + 1024*(KRANC_GFOFFSET3D(u,-3,0,-3) + KRANC_GFOFFSET3D(u,3,0,3)) + 2016*(KRANC_GFOFFSET3D(u,-1,0,4) + KRANC_GFOFFSET3D(u,1,0,-4) + KRANC_GFOFFSET3D(u,-4,0,1) + KRANC_GFOFFSET3D(u,4,0,-1)) - 2016*(KRANC_GFOFFSET3D(u,-1,0,-4) + KRANC_GFOFFSET3D(u,1,0,4) + KRANC_GFOFFSET3D(u,-4,0,-1) + KRANC_GFOFFSET3D(u,4,0,1)) - 504*(KRANC_GFOFFSET3D(u,-2,0,4) + KRANC_GFOFFSET3D(u,2,0,-4) + KRANC_GFOFFSET3D(u,-4,0,2) + KRANC_GFOFFSET3D(u,4,0,-2)) + 504*(KRANC_GFOFFSET3D(u,-2,0,-4) + KRANC_GFOFFSET3D(u,2,0,4) + KRANC_GFOFFSET3D(u,-4,0,-2) + KRANC_GFOFFSET3D(u,4,0,2)) + 96*(KRANC_GFOFFSET3D(u,-3,0,4) + KRANC_GFOFFSET3D(u,3,0,-4) + KRANC_GFOFFSET3D(u,-4,0,3) + KRANC_GFOFFSET3D(u,4,0,-3)) - 96*(KRANC_GFOFFSET3D(u,-3,0,-4) + KRANC_GFOFFSET3D(u,3,0,4) + KRANC_GFOFFSET3D(u,-4,0,-3) + KRANC_GFOFFSET3D(u,4,0,3)) - 9*(KRANC_GFOFFSET3D(u,-4,0,4) + KRANC_GFOFFSET3D(u,4,0,-4)) + 9*(KRANC_GFOFFSET3D(u,-4,0,-4) + KRANC_GFOFFSET3D(u,4,0,4)))*p1o705600dxdz)
#else
#  define PDstandardNth31(u) (PDstandardNth31_impl(u,p1o705600dxdz,cdj,cdk))
static CCTK_REAL PDstandardNth31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth12_impl(u, p1o705600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNth32(u) ((-451584*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 451584*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 112896*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 112896*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 28224*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 28224*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 21504*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 21504*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 5376*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 5376*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) - 1024*(KRANC_GFOFFSET3D(u,0,-3,3) + KRANC_GFOFFSET3D(u,0,3,-3)) + 1024*(KRANC_GFOFFSET3D(u,0,-3,-3) + KRANC_GFOFFSET3D(u,0,3,3)) + 2016*(KRANC_GFOFFSET3D(u,0,-1,4) + KRANC_GFOFFSET3D(u,0,1,-4) + KRANC_GFOFFSET3D(u,0,-4,1) + KRANC_GFOFFSET3D(u,0,4,-1)) - 2016*(KRANC_GFOFFSET3D(u,0,-1,-4) + KRANC_GFOFFSET3D(u,0,1,4) + KRANC_GFOFFSET3D(u,0,-4,-1) + KRANC_GFOFFSET3D(u,0,4,1)) - 504*(KRANC_GFOFFSET3D(u,0,-2,4) + KRANC_GFOFFSET3D(u,0,2,-4) + KRANC_GFOFFSET3D(u,0,-4,2) + KRANC_GFOFFSET3D(u,0,4,-2)) + 504*(KRANC_GFOFFSET3D(u,0,-2,-4) + KRANC_GFOFFSET3D(u,0,2,4) + KRANC_GFOFFSET3D(u,0,-4,-2) + KRANC_GFOFFSET3D(u,0,4,2)) + 96*(KRANC_GFOFFSET3D(u,0,-3,4) + KRANC_GFOFFSET3D(u,0,3,-4) + KRANC_GFOFFSET3D(u,0,-4,3) + KRANC_GFOFFSET3D(u,0,4,-3)) - 96*(KRANC_GFOFFSET3D(u,0,-3,-4) + KRANC_GFOFFSET3D(u,0,3,4) + KRANC_GFOFFSET3D(u,0,-4,-3) + KRANC_GFOFFSET3D(u,0,4,3)) - 9*(KRANC_GFOFFSET3D(u,0,-4,4) + KRANC_GFOFFSET3D(u,0,4,-4)) + 9*(KRANC_GFOFFSET3D(u,0,-4,-4) + KRANC_GFOFFSET3D(u,0,4,4)))*p1o705600dydz)
#else
#  define PDstandardNth32(u) (PDstandardNth32_impl(u,p1o705600dydz,cdj,cdk))
static CCTK_REAL PDstandardNth32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNth32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o705600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNth23_impl(u, p1o705600dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd1(u) ((-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx)
#else
#  define PDstandard2nd1(u) (PDstandard2nd1_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL PDstandard2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd2(u) ((-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy)
#else
#  define PDstandard2nd2(u) (PDstandard2nd2_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL PDstandard2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd3(u) ((-KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1))*p1o2dz)
#else
#  define PDstandard2nd3(u) (PDstandard2nd3_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL PDstandard2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd2_impl(u, p1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth1(u) ((-252*KRANC_GFOFFSET3D(u,0,0,0) + 210*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 120*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 45*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)) - 10*(KRANC_GFOFFSET3D(u,-4,0,0) + KRANC_GFOFFSET3D(u,4,0,0)) + KRANC_GFOFFSET3D(u,-5,0,0) + KRANC_GFOFFSET3D(u,5,0,0))*p1o1024dx)
#else
#  define PDdissipationNth1(u) (PDdissipationNth1_impl(u,p1o1024dx,cdj,cdk))
static CCTK_REAL PDdissipationNth1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1024dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNth1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1024dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-252*KRANC_GFOFFSET3D(u,0,0,0) + 210*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 120*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 45*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)) - 10*(KRANC_GFOFFSET3D(u,-4,0,0) + KRANC_GFOFFSET3D(u,4,0,0)) + KRANC_GFOFFSET3D(u,-5,0,0) + KRANC_GFOFFSET3D(u,5,0,0))*p1o1024dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth2(u) ((-252*KRANC_GFOFFSET3D(u,0,0,0) + 210*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 120*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 45*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)) - 10*(KRANC_GFOFFSET3D(u,0,-4,0) + KRANC_GFOFFSET3D(u,0,4,0)) + KRANC_GFOFFSET3D(u,0,-5,0) + KRANC_GFOFFSET3D(u,0,5,0))*p1o1024dy)
#else
#  define PDdissipationNth2(u) (PDdissipationNth2_impl(u,p1o1024dy,cdj,cdk))
static CCTK_REAL PDdissipationNth2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1024dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNth2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1024dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-252*KRANC_GFOFFSET3D(u,0,0,0) + 210*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 120*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 45*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)) - 10*(KRANC_GFOFFSET3D(u,0,-4,0) + KRANC_GFOFFSET3D(u,0,4,0)) + KRANC_GFOFFSET3D(u,0,-5,0) + KRANC_GFOFFSET3D(u,0,5,0))*p1o1024dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNth3(u) ((-252*KRANC_GFOFFSET3D(u,0,0,0) + 210*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 120*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + 45*(KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3)) - 10*(KRANC_GFOFFSET3D(u,0,0,-4) + KRANC_GFOFFSET3D(u,0,0,4)) + KRANC_GFOFFSET3D(u,0,0,-5) + KRANC_GFOFFSET3D(u,0,0,5))*p1o1024dz)
#else
#  define PDdissipationNth3(u) (PDdissipationNth3_impl(u,p1o1024dz,cdj,cdk))
static CCTK_REAL PDdissipationNth3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1024dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNth3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1024dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDdissipationNth2_impl(u, p1o1024dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth1(u) ((378*KRANC_GFOFFSET3D(u,0,0,0) - 60*KRANC_GFOFFSET3D(u,-2*dir1,0,0) + 5*KRANC_GFOFFSET3D(u,-3*dir1,0,0) - 140*KRANC_GFOFFSET3D(u,3*dir1,0,0) + 30*KRANC_GFOFFSET3D(u,4*dir1,0,0) - 3*KRANC_GFOFFSET3D(u,5*dir1,0,0) + 420*(KRANC_GFOFFSET3D(u,2*dir1,0,0) + KRANC_GFOFFSET3D(u,-dir1,0,0)) - 1050*KRANC_GFOFFSET3D(u,dir1,0,0))*pm1o840dx*dir1)
#else
#  define PDupwindNth1(u) (PDupwindNth1_impl(u,pm1o840dx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o840dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (378*KRANC_GFOFFSET3D(u,0,0,0) - 60*KRANC_GFOFFSET3D(u,-2*dir1,0,0) + 5*KRANC_GFOFFSET3D(u,-3*dir1,0,0) - 140*KRANC_GFOFFSET3D(u,3*dir1,0,0) + 30*KRANC_GFOFFSET3D(u,4*dir1,0,0) - 3*KRANC_GFOFFSET3D(u,5*dir1,0,0) + 420*(KRANC_GFOFFSET3D(u,2*dir1,0,0) + KRANC_GFOFFSET3D(u,-dir1,0,0)) - 1050*KRANC_GFOFFSET3D(u,dir1,0,0))*pm1o840dx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti1(u) ((-1470*KRANC_GFOFFSET3D(u,-1,0,0) + 1470*KRANC_GFOFFSET3D(u,1,0,0) + 480*KRANC_GFOFFSET3D(u,-2,0,0) - 480*KRANC_GFOFFSET3D(u,2,0,0) - 145*KRANC_GFOFFSET3D(u,-3,0,0) + 145*KRANC_GFOFFSET3D(u,3,0,0) + 30*KRANC_GFOFFSET3D(u,-4,0,0) - 30*KRANC_GFOFFSET3D(u,4,0,0) - 3*KRANC_GFOFFSET3D(u,-5,0,0) + 3*KRANC_GFOFFSET3D(u,5,0,0))*p1o1680dx)
#else
#  define PDupwindNthAnti1(u) (PDupwindNthAnti1_impl(u,p1o1680dx,cdj,cdk))
static CCTK_REAL PDupwindNthAnti1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1680dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAnti1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1680dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-1470*KRANC_GFOFFSET3D(u,-1,0,0) + 1470*KRANC_GFOFFSET3D(u,1,0,0) + 480*KRANC_GFOFFSET3D(u,-2,0,0) - 480*KRANC_GFOFFSET3D(u,2,0,0) - 145*KRANC_GFOFFSET3D(u,-3,0,0) + 145*KRANC_GFOFFSET3D(u,3,0,0) + 30*KRANC_GFOFFSET3D(u,-4,0,0) - 30*KRANC_GFOFFSET3D(u,4,0,0) - 3*KRANC_GFOFFSET3D(u,-5,0,0) + 3*KRANC_GFOFFSET3D(u,5,0,0))*p1o1680dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm1(u) ((-252*KRANC_GFOFFSET3D(u,0,0,0) + 210*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 120*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 45*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)) - 10*(KRANC_GFOFFSET3D(u,-4,0,0) + KRANC_GFOFFSET3D(u,4,0,0)) + KRANC_GFOFFSET3D(u,-5,0,0) + KRANC_GFOFFSET3D(u,5,0,0))*p1o560dx)
#else
#  define PDupwindNthSymm1(u) (PDupwindNthSymm1_impl(u,p1o560dx,cdj,cdk))
static CCTK_REAL PDupwindNthSymm1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o560dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymm1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o560dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-252*KRANC_GFOFFSET3D(u,0,0,0) + 210*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 120*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 45*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)) - 10*(KRANC_GFOFFSET3D(u,-4,0,0) + KRANC_GFOFFSET3D(u,4,0,0)) + KRANC_GFOFFSET3D(u,-5,0,0) + KRANC_GFOFFSET3D(u,5,0,0))*p1o560dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided1(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,dir1,0,0))*p1odx*dir1)
#else
#  define PDonesided1(u) (PDonesided1_impl(u,p1odx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,dir1,0,0))*p1odx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedPlus2nd1(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) - 4*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o2dx)
#else
#  define PDonesidedPlus2nd1(u) (PDonesidedPlus2nd1_impl(u,pm1o2dx,cdj,cdk))
static CCTK_REAL PDonesidedPlus2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesidedPlus2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (3*KRANC_GFOFFSET3D(u,0,0,0) - 4*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o2dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedMinus2nd1(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) - 4*KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0))*p1o2dx)
#else
#  define PDonesidedMinus2nd1(u) (PDonesidedMinus2nd1_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL PDonesidedMinus2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesidedMinus2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (3*KRANC_GFOFFSET3D(u,0,0,0) - 4*KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0))*p1o2dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth2(u) ((378*KRANC_GFOFFSET3D(u,0,0,0) - 60*KRANC_GFOFFSET3D(u,0,-2*dir2,0) + 5*KRANC_GFOFFSET3D(u,0,-3*dir2,0) - 140*KRANC_GFOFFSET3D(u,0,3*dir2,0) + 30*KRANC_GFOFFSET3D(u,0,4*dir2,0) - 3*KRANC_GFOFFSET3D(u,0,5*dir2,0) + 420*(KRANC_GFOFFSET3D(u,0,2*dir2,0) + KRANC_GFOFFSET3D(u,0,-dir2,0)) - 1050*KRANC_GFOFFSET3D(u,0,dir2,0))*pm1o840dy*dir2)
#else
#  define PDupwindNth2(u) (PDupwindNth2_impl(u,pm1o840dy,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o840dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (378*KRANC_GFOFFSET3D(u,0,0,0) - 60*KRANC_GFOFFSET3D(u,0,-2*dir2,0) + 5*KRANC_GFOFFSET3D(u,0,-3*dir2,0) - 140*KRANC_GFOFFSET3D(u,0,3*dir2,0) + 30*KRANC_GFOFFSET3D(u,0,4*dir2,0) - 3*KRANC_GFOFFSET3D(u,0,5*dir2,0) + 420*(KRANC_GFOFFSET3D(u,0,2*dir2,0) + KRANC_GFOFFSET3D(u,0,-dir2,0)) - 1050*KRANC_GFOFFSET3D(u,0,dir2,0))*pm1o840dy*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti2(u) ((-1470*KRANC_GFOFFSET3D(u,0,-1,0) + 1470*KRANC_GFOFFSET3D(u,0,1,0) + 480*KRANC_GFOFFSET3D(u,0,-2,0) - 480*KRANC_GFOFFSET3D(u,0,2,0) - 145*KRANC_GFOFFSET3D(u,0,-3,0) + 145*KRANC_GFOFFSET3D(u,0,3,0) + 30*KRANC_GFOFFSET3D(u,0,-4,0) - 30*KRANC_GFOFFSET3D(u,0,4,0) - 3*KRANC_GFOFFSET3D(u,0,-5,0) + 3*KRANC_GFOFFSET3D(u,0,5,0))*p1o1680dy)
#else
#  define PDupwindNthAnti2(u) (PDupwindNthAnti2_impl(u,p1o1680dy,cdj,cdk))
static CCTK_REAL PDupwindNthAnti2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1680dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAnti2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1680dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-1470*KRANC_GFOFFSET3D(u,0,-1,0) + 1470*KRANC_GFOFFSET3D(u,0,1,0) + 480*KRANC_GFOFFSET3D(u,0,-2,0) - 480*KRANC_GFOFFSET3D(u,0,2,0) - 145*KRANC_GFOFFSET3D(u,0,-3,0) + 145*KRANC_GFOFFSET3D(u,0,3,0) + 30*KRANC_GFOFFSET3D(u,0,-4,0) - 30*KRANC_GFOFFSET3D(u,0,4,0) - 3*KRANC_GFOFFSET3D(u,0,-5,0) + 3*KRANC_GFOFFSET3D(u,0,5,0))*p1o1680dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm2(u) ((-252*KRANC_GFOFFSET3D(u,0,0,0) + 210*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 120*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 45*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)) - 10*(KRANC_GFOFFSET3D(u,0,-4,0) + KRANC_GFOFFSET3D(u,0,4,0)) + KRANC_GFOFFSET3D(u,0,-5,0) + KRANC_GFOFFSET3D(u,0,5,0))*p1o560dy)
#else
#  define PDupwindNthSymm2(u) (PDupwindNthSymm2_impl(u,p1o560dy,cdj,cdk))
static CCTK_REAL PDupwindNthSymm2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o560dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymm2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o560dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-252*KRANC_GFOFFSET3D(u,0,0,0) + 210*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 120*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 45*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)) - 10*(KRANC_GFOFFSET3D(u,0,-4,0) + KRANC_GFOFFSET3D(u,0,4,0)) + KRANC_GFOFFSET3D(u,0,-5,0) + KRANC_GFOFFSET3D(u,0,5,0))*p1o560dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided2(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,dir2,0))*p1ody*dir2)
#else
#  define PDonesided2(u) (PDonesided2_impl(u,p1ody,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,dir2,0))*p1ody*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedPlus2nd2(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) - 4*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o2dy)
#else
#  define PDonesidedPlus2nd2(u) (PDonesidedPlus2nd2_impl(u,pm1o2dy,cdj,cdk))
static CCTK_REAL PDonesidedPlus2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesidedPlus2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (3*KRANC_GFOFFSET3D(u,0,0,0) - 4*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o2dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedMinus2nd2(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) - 4*KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,-2,0))*p1o2dy)
#else
#  define PDonesidedMinus2nd2(u) (PDonesidedMinus2nd2_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL PDonesidedMinus2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesidedMinus2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (3*KRANC_GFOFFSET3D(u,0,0,0) - 4*KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,-2,0))*p1o2dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNth3(u) ((378*KRANC_GFOFFSET3D(u,0,0,0) - 60*KRANC_GFOFFSET3D(u,0,0,-2*dir3) + 5*KRANC_GFOFFSET3D(u,0,0,-3*dir3) - 140*KRANC_GFOFFSET3D(u,0,0,3*dir3) + 30*KRANC_GFOFFSET3D(u,0,0,4*dir3) - 3*KRANC_GFOFFSET3D(u,0,0,5*dir3) + 420*(KRANC_GFOFFSET3D(u,0,0,2*dir3) + KRANC_GFOFFSET3D(u,0,0,-dir3)) - 1050*KRANC_GFOFFSET3D(u,0,0,dir3))*pm1o840dz*dir3)
#else
#  define PDupwindNth3(u) (PDupwindNth3_impl(u,pm1o840dz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNth3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNth3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o840dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNth2_impl(u, pm1o840dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAnti3(u) ((-1470*KRANC_GFOFFSET3D(u,0,0,-1) + 1470*KRANC_GFOFFSET3D(u,0,0,1) + 480*KRANC_GFOFFSET3D(u,0,0,-2) - 480*KRANC_GFOFFSET3D(u,0,0,2) - 145*KRANC_GFOFFSET3D(u,0,0,-3) + 145*KRANC_GFOFFSET3D(u,0,0,3) + 30*KRANC_GFOFFSET3D(u,0,0,-4) - 30*KRANC_GFOFFSET3D(u,0,0,4) - 3*KRANC_GFOFFSET3D(u,0,0,-5) + 3*KRANC_GFOFFSET3D(u,0,0,5))*p1o1680dz)
#else
#  define PDupwindNthAnti3(u) (PDupwindNthAnti3_impl(u,p1o1680dz,cdj,cdk))
static CCTK_REAL PDupwindNthAnti3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1680dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAnti3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o1680dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthAnti2_impl(u, p1o1680dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymm3(u) ((-252*KRANC_GFOFFSET3D(u,0,0,0) + 210*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 120*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + 45*(KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3)) - 10*(KRANC_GFOFFSET3D(u,0,0,-4) + KRANC_GFOFFSET3D(u,0,0,4)) + KRANC_GFOFFSET3D(u,0,0,-5) + KRANC_GFOFFSET3D(u,0,0,5))*p1o560dz)
#else
#  define PDupwindNthSymm3(u) (PDupwindNthSymm3_impl(u,p1o560dz,cdj,cdk))
static CCTK_REAL PDupwindNthSymm3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o560dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymm3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o560dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthSymm2_impl(u, p1o560dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided3(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,dir3))*p1odz*dir3)
#else
#  define PDonesided3(u) (PDonesided3_impl(u,p1odz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDonesided2_impl(u, p1odz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedPlus2nd3(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) - 4*KRANC_GFOFFSET3D(u,0,0,1) + KRANC_GFOFFSET3D(u,0,0,2))*pm1o2dz)
#else
#  define PDonesidedPlus2nd3(u) (PDonesidedPlus2nd3_impl(u,pm1o2dz,cdj,cdk))
static CCTK_REAL PDonesidedPlus2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesidedPlus2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDonesidedPlus2nd2_impl(u, pm1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesidedMinus2nd3(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) - 4*KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,-2))*p1o2dz)
#else
#  define PDonesidedMinus2nd3(u) (PDonesidedMinus2nd3_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL PDonesidedMinus2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesidedMinus2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDonesidedMinus2nd2_impl(u, p1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDPlus1(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx)
#else
#  define PDPlus1(u) (PDPlus1_impl(u,p1odx,cdj,cdk))
static CCTK_REAL PDPlus1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDPlus1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDPlus2(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody)
#else
#  define PDPlus2(u) (PDPlus2_impl(u,p1ody,cdj,cdk))
static CCTK_REAL PDPlus2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDPlus2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDPlus3(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,1))*p1odz)
#else
#  define PDPlus3(u) (PDPlus3_impl(u,p1odz,cdj,cdk))
static CCTK_REAL PDPlus3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDPlus3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDPlus2_impl(u, p1odz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDMinus1(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,-1,0,0))*p1odx)
#else
#  define PDMinus1(u) (PDMinus1_impl(u,p1odx,cdj,cdk))
static CCTK_REAL PDMinus1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDMinus1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,-1,0,0))*p1odx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDMinus2(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,-1,0))*p1ody)
#else
#  define PDMinus2(u) (PDMinus2_impl(u,p1ody,cdj,cdk))
static CCTK_REAL PDMinus2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDMinus2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,-1,0))*p1ody;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDMinus3(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,0,-1))*p1odz)
#else
#  define PDMinus3(u) (PDMinus3_impl(u,p1odz,cdj,cdk))
static CCTK_REAL PDMinus3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDMinus3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDMinus2_impl(u, p1odz, cdk, cdj);
}
#endif

