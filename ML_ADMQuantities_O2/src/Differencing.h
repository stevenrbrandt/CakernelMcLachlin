#define PDstandardNth1(u,i,j,k) (p1o2dx*(-(u)[index+di*(-1)+dj*(0)+dk*(0)] + (u)[index+di*(1)+dj*(0)+dk*(0)]))
#define PDstandardNth2(u,i,j,k) (p1o2dy*(-(u)[index+di*(0)+dj*(-1)+dk*(0)] + (u)[index+di*(0)+dj*(1)+dk*(0)]))
#define PDstandardNth3(u,i,j,k) (p1o2dz*(-(u)[index+di*(0)+dj*(0)+dk*(-1)] + (u)[index+di*(0)+dj*(0)+dk*(1)]))
#define PDstandardNth11(u,i,j,k) (p1odx2*(-2*(u)[index+di*(0)+dj*(0)+dk*(0)] + (u)[index+di*(-1)+dj*(0)+dk*(0)] + (u)[index+di*(1)+dj*(0)+dk*(0)]))
#define PDstandardNth22(u,i,j,k) (p1ody2*(-2*(u)[index+di*(0)+dj*(0)+dk*(0)] + (u)[index+di*(0)+dj*(-1)+dk*(0)] + (u)[index+di*(0)+dj*(1)+dk*(0)]))
#define PDstandardNth33(u,i,j,k) (p1odz2*(-2*(u)[index+di*(0)+dj*(0)+dk*(0)] + (u)[index+di*(0)+dj*(0)+dk*(-1)] + (u)[index+di*(0)+dj*(0)+dk*(1)]))
#define PDstandardNth12(u,i,j,k) (p1o4dxdy*((u)[index+di*(-1)+dj*(-1)+dk*(0)] - (u)[index+di*(-1)+dj*(1)+dk*(0)] - (u)[index+di*(1)+dj*(-1)+dk*(0)] + (u)[index+di*(1)+dj*(1)+dk*(0)]))
#define PDstandardNth13(u,i,j,k) (p1o4dxdz*((u)[index+di*(-1)+dj*(0)+dk*(-1)] - (u)[index+di*(-1)+dj*(0)+dk*(1)] - (u)[index+di*(1)+dj*(0)+dk*(-1)] + (u)[index+di*(1)+dj*(0)+dk*(1)]))
#define PDstandardNth21(u,i,j,k) (p1o4dxdy*((u)[index+di*(-1)+dj*(-1)+dk*(0)] - (u)[index+di*(-1)+dj*(1)+dk*(0)] - (u)[index+di*(1)+dj*(-1)+dk*(0)] + (u)[index+di*(1)+dj*(1)+dk*(0)]))
#define PDstandardNth23(u,i,j,k) (p1o4dydz*((u)[index+di*(0)+dj*(-1)+dk*(-1)] - (u)[index+di*(0)+dj*(-1)+dk*(1)] - (u)[index+di*(0)+dj*(1)+dk*(-1)] + (u)[index+di*(0)+dj*(1)+dk*(1)]))
#define PDstandardNth31(u,i,j,k) (p1o4dxdz*((u)[index+di*(-1)+dj*(0)+dk*(-1)] - (u)[index+di*(-1)+dj*(0)+dk*(1)] - (u)[index+di*(1)+dj*(0)+dk*(-1)] + (u)[index+di*(1)+dj*(0)+dk*(1)]))
#define PDstandardNth32(u,i,j,k) (p1o4dydz*((u)[index+di*(0)+dj*(-1)+dk*(-1)] - (u)[index+di*(0)+dj*(-1)+dk*(1)] - (u)[index+di*(0)+dj*(1)+dk*(-1)] + (u)[index+di*(0)+dj*(1)+dk*(1)]))
#define PDupwindNth1(u,i,j,k) (pm1o2dx*(3*(u)[index+di*(0)+dj*(0)+dk*(0)] + (u)[index+di*(2*dir1)+dj*(0)+dk*(0)] - 4*(u)[index+di*(dir1)+dj*(0)+dk*(0)])*dir1)
#define PDupwindNth2(u,i,j,k) (pm1o2dy*(3*(u)[index+di*(0)+dj*(0)+dk*(0)] + (u)[index+di*(0)+dj*(2*dir2)+dk*(0)] - 4*(u)[index+di*(0)+dj*(dir2)+dk*(0)])*dir2)
#define PDupwindNth3(u,i,j,k) (pm1o2dz*(3*(u)[index+di*(0)+dj*(0)+dk*(0)] + (u)[index+di*(0)+dj*(0)+dk*(2*dir3)] - 4*(u)[index+di*(0)+dj*(0)+dk*(dir3)])*dir3)
#define PDonesided1(u,i,j,k) (p1odx*(-(u)[index+di*(0)+dj*(0)+dk*(0)] + (u)[index+di*(dir1)+dj*(0)+dk*(0)])*dir1)
#define PDonesided2(u,i,j,k) (p1ody*(-(u)[index+di*(0)+dj*(0)+dk*(0)] + (u)[index+di*(0)+dj*(dir2)+dk*(0)])*dir2)
#define PDonesided3(u,i,j,k) (p1odz*(-(u)[index+di*(0)+dj*(0)+dk*(0)] + (u)[index+di*(0)+dj*(0)+dk*(dir3)])*dir3)
