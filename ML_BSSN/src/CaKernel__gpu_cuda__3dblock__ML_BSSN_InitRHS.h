#include"CaKernel__gpu_cuda__shared.h"
/*@@
 * @file    CaKernel__gpu_cuda__3dblock__ML_BSSN_InitRHS.h
 * @date    Thu Feb 02 02:59:33 CST 2012
 * @author  Marek Blazewicz
 * @desc
 * The prototype of the CaCUDA computational schema. It contains macros
 * which enable to declare, define and launch kernels as well as to copy
 * the data required for proper computations. The macros presented
 * in this file in the future will be automatically generated by the
 * Cactus parser depending on the input in interfaces.ccl file.
 * @enddesc
 * @version  $Header$
 *
 @@*/


#ifndef CAKERNEL__GPU_CUDA__3DBLOCK__ML_BSSN_INITRHS_H
#define CAKERNEL__GPU_CUDA__3DBLOCK__ML_BSSN_INITRHS_H
#include <algorithm>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

                                       


                                       




#ifdef __CUDACC__

/// !!!!!!!!!!!!! BEGIN of global definitions (not auto generated) !!!!!!!!!!!!!!!!

#define CAKERNEL_Threadsx 8
#define CAKERNEL_Threadsy 8
#define CAKERNEL_Threadsz 1

/* JT: 16x16x16 failed to compile on spider. 8x8x8 is ok to compile.
 * Let's make 8x8x8 as the default setting temporarily.
 * We will need to estimate the best configuration based on
 * the number of variables and the memory available.
 * */
#define CAKERNEL_Tilex 8
#define CAKERNEL_Tiley 8
#define CAKERNEL_Tilez 8

/// !!!!!!!!!!!!! END of global definitions (not auto generated) !!!!!!!!!!!!!!!!

/// !!!!!!!!!!!!!!!!!!!!!!!!! BEGIN ML_BSSN_InitRHS Kernel macors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef stncl_xn
#error "You can't include two header file in one execution file"
#endif

#define stncl_xn 0
#define stncl_xp 0
#define stncl_yn 0
#define stncl_yp 0
#define stncl_zn 0
#define stncl_zp 0

  
#define CAKERNEL_GFINDEX3D_ML_BSSN_InitRHS(ptr, i, j, k)                       \
 ptr[(i + gi) + params.cagh_ni * ((j + gj) + params.cagh_nj * (k + gk))]

#define CAKERNEL_GFINDEX3DV_ML_BSSN_InitRHS(ptr, i, j, k, vind)                \
 CAKERNEL_GFINDEX3D_ML_BSSN_InitRHS((&ptr[vind * vstride]), i, j, k)

#define I3DV CAKERNEL_GFINDEX3DV_ML_BSSN_InitRHS
#define I3D CAKERNEL_GFINDEX3D_ML_BSSN_InitRHS

   
/** Element  */
#define CAKERNEL_GFINDEX3D_ML_BSSN_InitRHS_l(ptr, i, j, k)                     \
  (ptr##_sh[k + lk][j + lj][i + li])

#define CAKERNEL_GFINDEX3DV_ML_BSSN_InitRHS_l(ptr, i, j, k, vind)              \
  (ptr##_sh[vind][k + lk][j + lj][i + li])

#define I3D_l CAKERNEL_GFINDEX3D_ML_BSSN_InitRHS_l
#define I3DV_l CAKERNEL_GFINDEX3DV_ML_BSSN_InitRHS_l
 


#define CAKERNEL_ML_BSSN_InitRHS_Declare_Begin_s                               \
template<class params_type>                                                    \
__global__ void CAKERNEL_ML_BSSN_InitRHS(                                      \
                                                                               \
                                                                               \
 CCTK_REAL * alpharhs, CCTK_REAL * Arhs, CCTK_REAL * At11rhs, CCTK_REAL * At12rhs, CCTK_REAL * At13rhs, CCTK_REAL * At22rhs, CCTK_REAL * At23rhs, CCTK_REAL * At33rhs, CCTK_REAL * B1rhs, CCTK_REAL * B2rhs, CCTK_REAL * B3rhs, CCTK_REAL * beta1rhs, CCTK_REAL * beta2rhs, CCTK_REAL * beta3rhs, CCTK_REAL * gt11rhs, CCTK_REAL * gt12rhs, CCTK_REAL * gt13rhs, CCTK_REAL * gt22rhs, CCTK_REAL * gt23rhs, CCTK_REAL * gt33rhs, CCTK_REAL * phirhs, CCTK_REAL * trKrhs, CCTK_REAL * Xt1rhs, CCTK_REAL * Xt2rhs, CCTK_REAL * Xt3rhs,\
                                                                               \
                                                                               \
                                                                               \
                                                                               \
                                                                               \
                                                                               \
/** Statically added variables to each kernel: */                              \
const params_type params)                                                      \
{                                                                           
# define CAKERNEL_ML_BSSN_InitRHS_Declare_Cached_Variables_s                   \
/** Kernel specific variables declaration. */                                  \
                                                                               \
                                                                               \
                                                                               \
                                                                               \

/** Common variables declaration; values are kernel specific. */            
# define CAKERNEL_ML_BSSN_InitRHS_Declare_Flow_Variables_s                     \
  int li = threadIdx.x;                                                        \
  int lj = threadIdx.y;                                                        \
  int lk = threadIdx.z + stncl_zn;                                             \
  int gi = blockIdx.x * (CAKERNEL_Tilex - stncl_xn - stncl_xp) + li;           \
  int gj = (blockIdx.y % params.cagh_blocky) *                                 \
          (CAKERNEL_Tiley - stncl_yn - stncl_yp) + lj;                         \
  int gk2= (blockIdx.y / params.cagh_blocky) *                                 \
          (CAKERNEL_Tilez - stncl_zn - stncl_zp) + lk;                         \
  int gk = gk2;                                                                \
  int vstride = params.cagh_ni * params.cagh_nj * params.cagh_nk;              \
  bool fetch_data = gi < params.cagh_ni && gj < params.cagh_nj;                \
  bool compute = gi < params.cagh_ni - stncl_xp && gj < params.cagh_nj -       \
    stncl_yp && li >= stncl_xn && lj >= stncl_yn &&                            \
    li < CAKERNEL_Tilex - stncl_xp &&                                          \
    lj < CAKERNEL_Tiley - stncl_yp;                                            \
  int tilez_to = min(CAKERNEL_Tilez - stncl_zp - stncl_zn,                     \
                        params.cagh_nk - gk - stncl_zp);                       \
    /** Dynamically set fetching from global memory */                        
# define CAKERNEL_ML_BSSN_InitRHS_Limit_Threads_To_LSH_Begin_s                 \
  if(fetch_data)                                                               \
  {                                                                           
#   define CAKERNEL_ML_BSSN_InitRHS_Fetch_Data_To_Cache_s                      \
                                                                               \

#   define CAKERNEL_ML_BSSN_InitRHS_Computations_Begin_s                       \
    for ( short tmpj = 0; tmpj < tilez_to; tmpj++)                             \
    {                                                                          \
      SYNC_BLOCK();                                                          
#     define CAKERNEL_ML_BSSN_InitRHS_Iterate_Local_Tile_s                     \
                                                                               \
      gk = gk2 + tmpj;                                                      
//    The loop is suppose to iterate local variables, as the tiles 'walks'     \
//      through the z dimension. For cached variables only!                    \
                                                                             
#     define CAKERNEL_ML_BSSN_InitRHS_Fetch_Front_Tile_To_Cache_s              \
                                                                               \
                                                                               \
      SYNC_BLOCK();                                                      

#     define CAKERNEL_ML_BSSN_InitRHS_Limit_Threads_To_Compute_Begin_s         \
      if(compute)                                                              \
      {                                                                        \
      /*if(threadIdx.x == 1 && threadIdx.y == 1)                               \
          printf("3cmpt [%02d, %02d, %02d]\n", gi, gj, gk);*/                  \
         /** TODO Add your computations here */                                \
         /** TODO Store the results to global array ({...}_out)  */          
                                                                             
#     define CAKERNEL_ML_BSSN_InitRHS_Limit_Threads_To_Compute_End_s           \
      }                                                                      
#   define CAKERNEL_ML_BSSN_InitRHS_Computations_End_s                         \
    }                                                                        
# define CAKERNEL_ML_BSSN_InitRHS_Limit_Threads_To_LSH_End_s                   \
  }                                                                          
#define CAKERNEL_ML_BSSN_InitRHS_Declare_End_s                                 \
}

#define CAKERNEL_ML_BSSN_InitRHS_Begin                                         \
CAKERNEL_ML_BSSN_InitRHS_Declare_Begin_s                                       \
  CAKERNEL_ML_BSSN_InitRHS_Declare_Cached_Variables_s                          \
  CAKERNEL_ML_BSSN_InitRHS_Declare_Flow_Variables_s                            \
  CAKERNEL_ML_BSSN_InitRHS_Limit_Threads_To_LSH_Begin_s                        \
    CAKERNEL_ML_BSSN_InitRHS_Fetch_Data_To_Cache_s                                   
                                                                              
#   define CAKERNEL_ML_BSSN_InitRHS_Computations_Begin                         \
    CAKERNEL_ML_BSSN_InitRHS_Computations_Begin_s                              \
      CAKERNEL_ML_BSSN_InitRHS_Iterate_Local_Tile_s                            \
      CAKERNEL_ML_BSSN_InitRHS_Fetch_Front_Tile_To_Cache_s                     \
      CAKERNEL_ML_BSSN_InitRHS_Limit_Threads_To_Compute_Begin_s                      
                                                                              
#   define CAKERNEL_ML_BSSN_InitRHS_Computations_End                           \
      CAKERNEL_ML_BSSN_InitRHS_Limit_Threads_To_Compute_End_s                  \
    CAKERNEL_ML_BSSN_InitRHS_Computations_End_s                                \
  CAKERNEL_ML_BSSN_InitRHS_Limit_Threads_To_LSH_End_s                    
                                                                            
#define CAKERNEL_ML_BSSN_InitRHS_End                                           \
CAKERNEL_ML_BSSN_InitRHS_Declare_End_s


/// !!!!!!!!!!!!!!!!!!!!!!!!! END ML_BSSN_InitRHS Kernel macors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#endif

#endif /* CAKERNEL__GPU_CUDA__3DBLOCK__ML_BSSN_INITRHS_H */
