#define PDstandardNth1(u,i,j,k) (p1o12dx*(-8*vec_loadu((u)[index+di*(-1)+dj*(0)+dk*(0)]) + 8*vec_loadu((u)[index+di*(1)+dj*(0)+dk*(0)]) + vec_loadu((u)[index+di*(-2)+dj*(0)+dk*(0)]) - vec_loadu((u)[index+di*(2)+dj*(0)+dk*(0)])))
#define PDstandardNth2(u,i,j,k) (p1o12dy*(-8*vec_loadu((u)[index+di*(0)+dj*(-1)+dk*(0)]) + 8*vec_loadu((u)[index+di*(0)+dj*(1)+dk*(0)]) + vec_loadu((u)[index+di*(0)+dj*(-2)+dk*(0)]) - vec_loadu((u)[index+di*(0)+dj*(2)+dk*(0)])))
#define PDstandardNth3(u,i,j,k) (p1o12dz*(-8*vec_loadu((u)[index+di*(0)+dj*(0)+dk*(-1)]) + 8*vec_loadu((u)[index+di*(0)+dj*(0)+dk*(1)]) + vec_loadu((u)[index+di*(0)+dj*(0)+dk*(-2)]) - vec_loadu((u)[index+di*(0)+dj*(0)+dk*(2)])))
#define PDstandardNth11(u,i,j,k) (pm1o12dx2*(30*vec_loadu((u)[index+di*(0)+dj*(0)+dk*(0)]) - 16*(vec_loadu((u)[index+di*(-1)+dj*(0)+dk*(0)]) + vec_loadu((u)[index+di*(1)+dj*(0)+dk*(0)])) + vec_loadu((u)[index+di*(-2)+dj*(0)+dk*(0)]) + vec_loadu((u)[index+di*(2)+dj*(0)+dk*(0)])))
#define PDstandardNth22(u,i,j,k) (pm1o12dy2*(30*vec_loadu((u)[index+di*(0)+dj*(0)+dk*(0)]) - 16*(vec_loadu((u)[index+di*(0)+dj*(-1)+dk*(0)]) + vec_loadu((u)[index+di*(0)+dj*(1)+dk*(0)])) + vec_loadu((u)[index+di*(0)+dj*(-2)+dk*(0)]) + vec_loadu((u)[index+di*(0)+dj*(2)+dk*(0)])))
#define PDstandardNth33(u,i,j,k) (pm1o12dz2*(30*vec_loadu((u)[index+di*(0)+dj*(0)+dk*(0)]) - 16*(vec_loadu((u)[index+di*(0)+dj*(0)+dk*(-1)]) + vec_loadu((u)[index+di*(0)+dj*(0)+dk*(1)])) + vec_loadu((u)[index+di*(0)+dj*(0)+dk*(-2)]) + vec_loadu((u)[index+di*(0)+dj*(0)+dk*(2)])))
#define PDstandardNth12(u,i,j,k) (p1o144dxdy*(-64*(vec_loadu((u)[index+di*(-1)+dj*(1)+dk*(0)]) + vec_loadu((u)[index+di*(1)+dj*(-1)+dk*(0)])) + 64*(vec_loadu((u)[index+di*(-1)+dj*(-1)+dk*(0)]) + vec_loadu((u)[index+di*(1)+dj*(1)+dk*(0)])) + 8*(vec_loadu((u)[index+di*(-1)+dj*(2)+dk*(0)]) + vec_loadu((u)[index+di*(1)+dj*(-2)+dk*(0)]) + vec_loadu((u)[index+di*(-2)+dj*(1)+dk*(0)]) + vec_loadu((u)[index+di*(2)+dj*(-1)+dk*(0)])) - 8*(vec_loadu((u)[index+di*(-1)+dj*(-2)+dk*(0)]) + vec_loadu((u)[index+di*(1)+dj*(2)+dk*(0)]) + vec_loadu((u)[index+di*(-2)+dj*(-1)+dk*(0)]) + vec_loadu((u)[index+di*(2)+dj*(1)+dk*(0)])) + vec_loadu((u)[index+di*(-2)+dj*(-2)+dk*(0)]) - vec_loadu((u)[index+di*(-2)+dj*(2)+dk*(0)]) - vec_loadu((u)[index+di*(2)+dj*(-2)+dk*(0)]) + vec_loadu((u)[index+di*(2)+dj*(2)+dk*(0)])))
#define PDstandardNth13(u,i,j,k) (p1o144dxdz*(-64*(vec_loadu((u)[index+di*(-1)+dj*(0)+dk*(1)]) + vec_loadu((u)[index+di*(1)+dj*(0)+dk*(-1)])) + 64*(vec_loadu((u)[index+di*(-1)+dj*(0)+dk*(-1)]) + vec_loadu((u)[index+di*(1)+dj*(0)+dk*(1)])) + 8*(vec_loadu((u)[index+di*(-1)+dj*(0)+dk*(2)]) + vec_loadu((u)[index+di*(1)+dj*(0)+dk*(-2)]) + vec_loadu((u)[index+di*(-2)+dj*(0)+dk*(1)]) + vec_loadu((u)[index+di*(2)+dj*(0)+dk*(-1)])) - 8*(vec_loadu((u)[index+di*(-1)+dj*(0)+dk*(-2)]) + vec_loadu((u)[index+di*(1)+dj*(0)+dk*(2)]) + vec_loadu((u)[index+di*(-2)+dj*(0)+dk*(-1)]) + vec_loadu((u)[index+di*(2)+dj*(0)+dk*(1)])) + vec_loadu((u)[index+di*(-2)+dj*(0)+dk*(-2)]) - vec_loadu((u)[index+di*(-2)+dj*(0)+dk*(2)]) - vec_loadu((u)[index+di*(2)+dj*(0)+dk*(-2)]) + vec_loadu((u)[index+di*(2)+dj*(0)+dk*(2)])))
#define PDstandardNth21(u,i,j,k) (p1o144dxdy*(-64*(vec_loadu((u)[index+di*(-1)+dj*(1)+dk*(0)]) + vec_loadu((u)[index+di*(1)+dj*(-1)+dk*(0)])) + 64*(vec_loadu((u)[index+di*(-1)+dj*(-1)+dk*(0)]) + vec_loadu((u)[index+di*(1)+dj*(1)+dk*(0)])) + 8*(vec_loadu((u)[index+di*(-1)+dj*(2)+dk*(0)]) + vec_loadu((u)[index+di*(1)+dj*(-2)+dk*(0)]) + vec_loadu((u)[index+di*(-2)+dj*(1)+dk*(0)]) + vec_loadu((u)[index+di*(2)+dj*(-1)+dk*(0)])) - 8*(vec_loadu((u)[index+di*(-1)+dj*(-2)+dk*(0)]) + vec_loadu((u)[index+di*(1)+dj*(2)+dk*(0)]) + vec_loadu((u)[index+di*(-2)+dj*(-1)+dk*(0)]) + vec_loadu((u)[index+di*(2)+dj*(1)+dk*(0)])) + vec_loadu((u)[index+di*(-2)+dj*(-2)+dk*(0)]) - vec_loadu((u)[index+di*(-2)+dj*(2)+dk*(0)]) - vec_loadu((u)[index+di*(2)+dj*(-2)+dk*(0)]) + vec_loadu((u)[index+di*(2)+dj*(2)+dk*(0)])))
#define PDstandardNth23(u,i,j,k) (p1o144dydz*(-64*(vec_loadu((u)[index+di*(0)+dj*(-1)+dk*(1)]) + vec_loadu((u)[index+di*(0)+dj*(1)+dk*(-1)])) + 64*(vec_loadu((u)[index+di*(0)+dj*(-1)+dk*(-1)]) + vec_loadu((u)[index+di*(0)+dj*(1)+dk*(1)])) + 8*(vec_loadu((u)[index+di*(0)+dj*(-1)+dk*(2)]) + vec_loadu((u)[index+di*(0)+dj*(1)+dk*(-2)]) + vec_loadu((u)[index+di*(0)+dj*(-2)+dk*(1)]) + vec_loadu((u)[index+di*(0)+dj*(2)+dk*(-1)])) - 8*(vec_loadu((u)[index+di*(0)+dj*(-1)+dk*(-2)]) + vec_loadu((u)[index+di*(0)+dj*(1)+dk*(2)]) + vec_loadu((u)[index+di*(0)+dj*(-2)+dk*(-1)]) + vec_loadu((u)[index+di*(0)+dj*(2)+dk*(1)])) + vec_loadu((u)[index+di*(0)+dj*(-2)+dk*(-2)]) - vec_loadu((u)[index+di*(0)+dj*(-2)+dk*(2)]) - vec_loadu((u)[index+di*(0)+dj*(2)+dk*(-2)]) + vec_loadu((u)[index+di*(0)+dj*(2)+dk*(2)])))
#define PDstandardNth31(u,i,j,k) (p1o144dxdz*(-64*(vec_loadu((u)[index+di*(-1)+dj*(0)+dk*(1)]) + vec_loadu((u)[index+di*(1)+dj*(0)+dk*(-1)])) + 64*(vec_loadu((u)[index+di*(-1)+dj*(0)+dk*(-1)]) + vec_loadu((u)[index+di*(1)+dj*(0)+dk*(1)])) + 8*(vec_loadu((u)[index+di*(-1)+dj*(0)+dk*(2)]) + vec_loadu((u)[index+di*(1)+dj*(0)+dk*(-2)]) + vec_loadu((u)[index+di*(-2)+dj*(0)+dk*(1)]) + vec_loadu((u)[index+di*(2)+dj*(0)+dk*(-1)])) - 8*(vec_loadu((u)[index+di*(-1)+dj*(0)+dk*(-2)]) + vec_loadu((u)[index+di*(1)+dj*(0)+dk*(2)]) + vec_loadu((u)[index+di*(-2)+dj*(0)+dk*(-1)]) + vec_loadu((u)[index+di*(2)+dj*(0)+dk*(1)])) + vec_loadu((u)[index+di*(-2)+dj*(0)+dk*(-2)]) - vec_loadu((u)[index+di*(-2)+dj*(0)+dk*(2)]) - vec_loadu((u)[index+di*(2)+dj*(0)+dk*(-2)]) + vec_loadu((u)[index+di*(2)+dj*(0)+dk*(2)])))
#define PDstandardNth32(u,i,j,k) (p1o144dydz*(-64*(vec_loadu((u)[index+di*(0)+dj*(-1)+dk*(1)]) + vec_loadu((u)[index+di*(0)+dj*(1)+dk*(-1)])) + 64*(vec_loadu((u)[index+di*(0)+dj*(-1)+dk*(-1)]) + vec_loadu((u)[index+di*(0)+dj*(1)+dk*(1)])) + 8*(vec_loadu((u)[index+di*(0)+dj*(-1)+dk*(2)]) + vec_loadu((u)[index+di*(0)+dj*(1)+dk*(-2)]) + vec_loadu((u)[index+di*(0)+dj*(-2)+dk*(1)]) + vec_loadu((u)[index+di*(0)+dj*(2)+dk*(-1)])) - 8*(vec_loadu((u)[index+di*(0)+dj*(-1)+dk*(-2)]) + vec_loadu((u)[index+di*(0)+dj*(1)+dk*(2)]) + vec_loadu((u)[index+di*(0)+dj*(-2)+dk*(-1)]) + vec_loadu((u)[index+di*(0)+dj*(2)+dk*(1)])) + vec_loadu((u)[index+di*(0)+dj*(-2)+dk*(-2)]) - vec_loadu((u)[index+di*(0)+dj*(-2)+dk*(2)]) - vec_loadu((u)[index+di*(0)+dj*(2)+dk*(-2)]) + vec_loadu((u)[index+di*(0)+dj*(2)+dk*(2)])))
