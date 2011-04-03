#define PDstandardNth1(u) (p1o12dx*(-8*(u)[di*(-1)+dj*(0)+dk*(0)] + 8*(u)[di*(1)+dj*(0)+dk*(0)] + (u)[di*(-2)+dj*(0)+dk*(0)] - (u)[di*(2)+dj*(0)+dk*(0)]))
#define PDstandardNth2(u) (p1o12dy*(-8*(u)[di*(0)+dj*(-1)+dk*(0)] + 8*(u)[di*(0)+dj*(1)+dk*(0)] + (u)[di*(0)+dj*(-2)+dk*(0)] - (u)[di*(0)+dj*(2)+dk*(0)]))
#define PDstandardNth3(u) (p1o12dz*(-8*(u)[di*(0)+dj*(0)+dk*(-1)] + 8*(u)[di*(0)+dj*(0)+dk*(1)] + (u)[di*(0)+dj*(0)+dk*(-2)] - (u)[di*(0)+dj*(0)+dk*(2)]))
#define PDstandardNth11(u) (pm1o12dx2*(30*(u)[di*(0)+dj*(0)+dk*(0)] - 16*((u)[di*(-1)+dj*(0)+dk*(0)] + (u)[di*(1)+dj*(0)+dk*(0)]) + (u)[di*(-2)+dj*(0)+dk*(0)] + (u)[di*(2)+dj*(0)+dk*(0)]))
#define PDstandardNth22(u) (pm1o12dy2*(30*(u)[di*(0)+dj*(0)+dk*(0)] - 16*((u)[di*(0)+dj*(-1)+dk*(0)] + (u)[di*(0)+dj*(1)+dk*(0)]) + (u)[di*(0)+dj*(-2)+dk*(0)] + (u)[di*(0)+dj*(2)+dk*(0)]))
#define PDstandardNth33(u) (pm1o12dz2*(30*(u)[di*(0)+dj*(0)+dk*(0)] - 16*((u)[di*(0)+dj*(0)+dk*(-1)] + (u)[di*(0)+dj*(0)+dk*(1)]) + (u)[di*(0)+dj*(0)+dk*(-2)] + (u)[di*(0)+dj*(0)+dk*(2)]))
#define PDstandardNth12(u) (p1o144dxdy*(-64*((u)[di*(-1)+dj*(1)+dk*(0)] + (u)[di*(1)+dj*(-1)+dk*(0)]) + 64*((u)[di*(-1)+dj*(-1)+dk*(0)] + (u)[di*(1)+dj*(1)+dk*(0)]) + 8*((u)[di*(-1)+dj*(2)+dk*(0)] + (u)[di*(1)+dj*(-2)+dk*(0)] + (u)[di*(-2)+dj*(1)+dk*(0)] + (u)[di*(2)+dj*(-1)+dk*(0)]) - 8*((u)[di*(-1)+dj*(-2)+dk*(0)] + (u)[di*(1)+dj*(2)+dk*(0)] + (u)[di*(-2)+dj*(-1)+dk*(0)] + (u)[di*(2)+dj*(1)+dk*(0)]) + (u)[di*(-2)+dj*(-2)+dk*(0)] - (u)[di*(-2)+dj*(2)+dk*(0)] - (u)[di*(2)+dj*(-2)+dk*(0)] + (u)[di*(2)+dj*(2)+dk*(0)]))
#define PDstandardNth13(u) (p1o144dxdz*(-64*((u)[di*(-1)+dj*(0)+dk*(1)] + (u)[di*(1)+dj*(0)+dk*(-1)]) + 64*((u)[di*(-1)+dj*(0)+dk*(-1)] + (u)[di*(1)+dj*(0)+dk*(1)]) + 8*((u)[di*(-1)+dj*(0)+dk*(2)] + (u)[di*(1)+dj*(0)+dk*(-2)] + (u)[di*(-2)+dj*(0)+dk*(1)] + (u)[di*(2)+dj*(0)+dk*(-1)]) - 8*((u)[di*(-1)+dj*(0)+dk*(-2)] + (u)[di*(1)+dj*(0)+dk*(2)] + (u)[di*(-2)+dj*(0)+dk*(-1)] + (u)[di*(2)+dj*(0)+dk*(1)]) + (u)[di*(-2)+dj*(0)+dk*(-2)] - (u)[di*(-2)+dj*(0)+dk*(2)] - (u)[di*(2)+dj*(0)+dk*(-2)] + (u)[di*(2)+dj*(0)+dk*(2)]))
#define PDstandardNth21(u) (p1o144dxdy*(-64*((u)[di*(-1)+dj*(1)+dk*(0)] + (u)[di*(1)+dj*(-1)+dk*(0)]) + 64*((u)[di*(-1)+dj*(-1)+dk*(0)] + (u)[di*(1)+dj*(1)+dk*(0)]) + 8*((u)[di*(-1)+dj*(2)+dk*(0)] + (u)[di*(1)+dj*(-2)+dk*(0)] + (u)[di*(-2)+dj*(1)+dk*(0)] + (u)[di*(2)+dj*(-1)+dk*(0)]) - 8*((u)[di*(-1)+dj*(-2)+dk*(0)] + (u)[di*(1)+dj*(2)+dk*(0)] + (u)[di*(-2)+dj*(-1)+dk*(0)] + (u)[di*(2)+dj*(1)+dk*(0)]) + (u)[di*(-2)+dj*(-2)+dk*(0)] - (u)[di*(-2)+dj*(2)+dk*(0)] - (u)[di*(2)+dj*(-2)+dk*(0)] + (u)[di*(2)+dj*(2)+dk*(0)]))
#define PDstandardNth23(u) (p1o144dydz*(-64*((u)[di*(0)+dj*(-1)+dk*(1)] + (u)[di*(0)+dj*(1)+dk*(-1)]) + 64*((u)[di*(0)+dj*(-1)+dk*(-1)] + (u)[di*(0)+dj*(1)+dk*(1)]) + 8*((u)[di*(0)+dj*(-1)+dk*(2)] + (u)[di*(0)+dj*(1)+dk*(-2)] + (u)[di*(0)+dj*(-2)+dk*(1)] + (u)[di*(0)+dj*(2)+dk*(-1)]) - 8*((u)[di*(0)+dj*(-1)+dk*(-2)] + (u)[di*(0)+dj*(1)+dk*(2)] + (u)[di*(0)+dj*(-2)+dk*(-1)] + (u)[di*(0)+dj*(2)+dk*(1)]) + (u)[di*(0)+dj*(-2)+dk*(-2)] - (u)[di*(0)+dj*(-2)+dk*(2)] - (u)[di*(0)+dj*(2)+dk*(-2)] + (u)[di*(0)+dj*(2)+dk*(2)]))
#define PDstandardNth31(u) (p1o144dxdz*(-64*((u)[di*(-1)+dj*(0)+dk*(1)] + (u)[di*(1)+dj*(0)+dk*(-1)]) + 64*((u)[di*(-1)+dj*(0)+dk*(-1)] + (u)[di*(1)+dj*(0)+dk*(1)]) + 8*((u)[di*(-1)+dj*(0)+dk*(2)] + (u)[di*(1)+dj*(0)+dk*(-2)] + (u)[di*(-2)+dj*(0)+dk*(1)] + (u)[di*(2)+dj*(0)+dk*(-1)]) - 8*((u)[di*(-1)+dj*(0)+dk*(-2)] + (u)[di*(1)+dj*(0)+dk*(2)] + (u)[di*(-2)+dj*(0)+dk*(-1)] + (u)[di*(2)+dj*(0)+dk*(1)]) + (u)[di*(-2)+dj*(0)+dk*(-2)] - (u)[di*(-2)+dj*(0)+dk*(2)] - (u)[di*(2)+dj*(0)+dk*(-2)] + (u)[di*(2)+dj*(0)+dk*(2)]))
#define PDstandardNth32(u) (p1o144dydz*(-64*((u)[di*(0)+dj*(-1)+dk*(1)] + (u)[di*(0)+dj*(1)+dk*(-1)]) + 64*((u)[di*(0)+dj*(-1)+dk*(-1)] + (u)[di*(0)+dj*(1)+dk*(1)]) + 8*((u)[di*(0)+dj*(-1)+dk*(2)] + (u)[di*(0)+dj*(1)+dk*(-2)] + (u)[di*(0)+dj*(-2)+dk*(1)] + (u)[di*(0)+dj*(2)+dk*(-1)]) - 8*((u)[di*(0)+dj*(-1)+dk*(-2)] + (u)[di*(0)+dj*(1)+dk*(2)] + (u)[di*(0)+dj*(-2)+dk*(-1)] + (u)[di*(0)+dj*(2)+dk*(1)]) + (u)[di*(0)+dj*(-2)+dk*(-2)] - (u)[di*(0)+dj*(-2)+dk*(2)] - (u)[di*(0)+dj*(2)+dk*(-2)] + (u)[di*(0)+dj*(2)+dk*(2)]))
#define PDdissipationNth1(u) (p1o64dx*(-20*(u)[di*(0)+dj*(0)+dk*(0)] + 15*((u)[di*(-1)+dj*(0)+dk*(0)] + (u)[di*(1)+dj*(0)+dk*(0)]) - 6*((u)[di*(-2)+dj*(0)+dk*(0)] + (u)[di*(2)+dj*(0)+dk*(0)]) + (u)[di*(-3)+dj*(0)+dk*(0)] + (u)[di*(3)+dj*(0)+dk*(0)]))
#define PDdissipationNth2(u) (p1o64dy*(-20*(u)[di*(0)+dj*(0)+dk*(0)] + 15*((u)[di*(0)+dj*(-1)+dk*(0)] + (u)[di*(0)+dj*(1)+dk*(0)]) - 6*((u)[di*(0)+dj*(-2)+dk*(0)] + (u)[di*(0)+dj*(2)+dk*(0)]) + (u)[di*(0)+dj*(-3)+dk*(0)] + (u)[di*(0)+dj*(3)+dk*(0)]))
#define PDdissipationNth3(u) (p1o64dz*(-20*(u)[di*(0)+dj*(0)+dk*(0)] + 15*((u)[di*(0)+dj*(0)+dk*(-1)] + (u)[di*(0)+dj*(0)+dk*(1)]) - 6*((u)[di*(0)+dj*(0)+dk*(-2)] + (u)[di*(0)+dj*(0)+dk*(2)]) + (u)[di*(0)+dj*(0)+dk*(-3)] + (u)[di*(0)+dj*(0)+dk*(3)]))
#define PDupwindNth1(u) (p1o12dx*(-10*(u)[di*(0)+dj*(0)+dk*(0)] - 6*(u)[di*(2*dir1)+dj*(0)+dk*(0)] + (u)[di*(3*dir1)+dj*(0)+dk*(0)] - 3*(u)[di*(-dir1)+dj*(0)+dk*(0)] + 18*(u)[di*(dir1)+dj*(0)+dk*(0)])*dir1)
#define PDupwindNth2(u) (p1o12dy*(-10*(u)[di*(0)+dj*(0)+dk*(0)] - 6*(u)[di*(0)+dj*(2*dir2)+dk*(0)] + (u)[di*(0)+dj*(3*dir2)+dk*(0)] - 3*(u)[di*(0)+dj*(-dir2)+dk*(0)] + 18*(u)[di*(0)+dj*(dir2)+dk*(0)])*dir2)
#define PDupwindNth3(u) (p1o12dz*(-10*(u)[di*(0)+dj*(0)+dk*(0)] - 6*(u)[di*(0)+dj*(0)+dk*(2*dir3)] + (u)[di*(0)+dj*(0)+dk*(3*dir3)] - 3*(u)[di*(0)+dj*(0)+dk*(-dir3)] + 18*(u)[di*(0)+dj*(0)+dk*(dir3)])*dir3)
#define PDonesided1(u) (p1odx*(-(u)[di*(0)+dj*(0)+dk*(0)] + (u)[di*(dir1)+dj*(0)+dk*(0)])*dir1)
#define PDonesided2(u) (p1ody*(-(u)[di*(0)+dj*(0)+dk*(0)] + (u)[di*(0)+dj*(dir2)+dk*(0)])*dir2)
#define PDonesided3(u) (p1odz*(-(u)[di*(0)+dj*(0)+dk*(0)] + (u)[di*(0)+dj*(0)+dk*(dir3)])*dir3)
#define PDupwindNthAnti1(u) (p1o24dx*(-21*(u)[di*(-1)+dj*(0)+dk*(0)] + 21*(u)[di*(1)+dj*(0)+dk*(0)] + 6*(u)[di*(-2)+dj*(0)+dk*(0)] - 6*(u)[di*(2)+dj*(0)+dk*(0)] - (u)[di*(-3)+dj*(0)+dk*(0)] + (u)[di*(3)+dj*(0)+dk*(0)]))
#define PDupwindNthSymm1(u) (p1o24dx*(-20*(u)[di*(0)+dj*(0)+dk*(0)] + 15*((u)[di*(-1)+dj*(0)+dk*(0)] + (u)[di*(1)+dj*(0)+dk*(0)]) - 6*((u)[di*(-2)+dj*(0)+dk*(0)] + (u)[di*(2)+dj*(0)+dk*(0)]) + (u)[di*(-3)+dj*(0)+dk*(0)] + (u)[di*(3)+dj*(0)+dk*(0)]))
#define PDupwindNthAnti2(u) (p1o24dy*(-21*(u)[di*(0)+dj*(-1)+dk*(0)] + 21*(u)[di*(0)+dj*(1)+dk*(0)] + 6*(u)[di*(0)+dj*(-2)+dk*(0)] - 6*(u)[di*(0)+dj*(2)+dk*(0)] - (u)[di*(0)+dj*(-3)+dk*(0)] + (u)[di*(0)+dj*(3)+dk*(0)]))
#define PDupwindNthSymm2(u) (p1o24dy*(-20*(u)[di*(0)+dj*(0)+dk*(0)] + 15*((u)[di*(0)+dj*(-1)+dk*(0)] + (u)[di*(0)+dj*(1)+dk*(0)]) - 6*((u)[di*(0)+dj*(-2)+dk*(0)] + (u)[di*(0)+dj*(2)+dk*(0)]) + (u)[di*(0)+dj*(-3)+dk*(0)] + (u)[di*(0)+dj*(3)+dk*(0)]))
#define PDupwindNthAnti3(u) (p1o24dz*(-21*(u)[di*(0)+dj*(0)+dk*(-1)] + 21*(u)[di*(0)+dj*(0)+dk*(1)] + 6*(u)[di*(0)+dj*(0)+dk*(-2)] - 6*(u)[di*(0)+dj*(0)+dk*(2)] - (u)[di*(0)+dj*(0)+dk*(-3)] + (u)[di*(0)+dj*(0)+dk*(3)]))
#define PDupwindNthSymm3(u) (p1o24dz*(-20*(u)[di*(0)+dj*(0)+dk*(0)] + 15*((u)[di*(0)+dj*(0)+dk*(-1)] + (u)[di*(0)+dj*(0)+dk*(1)]) - 6*((u)[di*(0)+dj*(0)+dk*(-2)] + (u)[di*(0)+dj*(0)+dk*(2)]) + (u)[di*(0)+dj*(0)+dk*(-3)] + (u)[di*(0)+dj*(0)+dk*(3)]))
