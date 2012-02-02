#include"CaKernel__gpu_cuda__3dblock__ML_BSSN_InitRHS.h"
#include"CaKernel__ML_BSSN_InitRHS.code"
#include"CaKernel__gpu_cuda__mem.h"



                                       
  

                                       
  


#include<assert.h>

void CAKERNEL_Launch_ML_BSSN_InitRHS(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

#define TS_ERROR(d)                                                            \
    if ( CAKERNEL_Tile##d <= stncl_##d##n + stncl_##d##p )                     \
      CCTK_VWarn                                                               \
       (CCTK_WARN_ABORT, __LINE__,__FILE__, CCTK_THORNSTRING,                  \
        "Tile size for ML_BSSN_InitRHS along %s axis (%d) too small for stencil %d -- %d\n",\
        #d, CAKERNEL_Tile##d, stncl_##d##n, stncl_##d##p);

    TS_ERROR(x); TS_ERROR(y); TS_ERROR(z);
#   undef TS_ERROR

    int vi = 0;
    
    
      assert((vi = CCTK_VarIndex("ML_BSSN::alpharhs"))>=0);
      void * d_alpharhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_alpharhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::Arhs"))>=0);
      void * d_Arhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_Arhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At11rhs"))>=0);
      void * d_At11rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At11rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At12rhs"))>=0);
      void * d_At12rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At12rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At13rhs"))>=0);
      void * d_At13rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At13rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At22rhs"))>=0);
      void * d_At22rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At22rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At23rhs"))>=0);
      void * d_At23rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At23rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At33rhs"))>=0);
      void * d_At33rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At33rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::B1rhs"))>=0);
      void * d_B1rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_B1rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::B2rhs"))>=0);
      void * d_B2rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_B2rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::B3rhs"))>=0);
      void * d_B3rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_B3rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::beta1rhs"))>=0);
      void * d_beta1rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_beta1rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::beta2rhs"))>=0);
      void * d_beta2rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_beta2rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::beta3rhs"))>=0);
      void * d_beta3rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_beta3rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt11rhs"))>=0);
      void * d_gt11rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt11rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt12rhs"))>=0);
      void * d_gt12rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt12rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt13rhs"))>=0);
      void * d_gt13rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt13rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt22rhs"))>=0);
      void * d_gt22rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt22rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt23rhs"))>=0);
      void * d_gt23rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt23rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt33rhs"))>=0);
      void * d_gt33rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt33rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::phirhs"))>=0);
      void * d_phirhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_phirhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::trKrhs"))>=0);
      void * d_trKrhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_trKrhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::Xt1rhs"))>=0);
      void * d_Xt1rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_Xt1rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::Xt2rhs"))>=0);
      void * d_Xt2rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_Xt2rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::Xt3rhs"))>=0);
      void * d_Xt3rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_Xt3rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
    
    size_t datasize = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

    
    
    const int blocky = iDivUp(cctk_lsh[1] - stncl_yn - stncl_yp, CAKERNEL_Tiley - stncl_yn - stncl_yp);

    CaCUDA_Kernel_Launch_Parameters prms(cctk_iteration,
        cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
        cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2], blocky,
        cctk_delta_space[0], cctk_delta_space[1], cctk_delta_space[2],
        cctk_delta_time,
        cctk_origin_space[0], cctk_origin_space[1], cctk_origin_space[2],
        cctk_time);


    CAKERNEL_ML_BSSN_InitRHS<<<                                                    
 dim3(iDivUp(prms.cagh_ni - stncl_xn - stncl_xp, CAKERNEL_Tilex - stncl_xn - stncl_xp), 
      iDivUp(prms.cagh_nk - stncl_zn - stncl_zp, CAKERNEL_Tilez - stncl_zn - stncl_zp) 
          * blocky),
 dim3(CAKERNEL_Tilex, CAKERNEL_Threadsy, CAKERNEL_Threadsz)>>>(
      

(typeof(alpharhs)) d_alpharhs,(typeof(Arhs)) d_Arhs,(typeof(At11rhs)) d_At11rhs,(typeof(At12rhs)) d_At12rhs,(typeof(At13rhs)) d_At13rhs,(typeof(At22rhs)) d_At22rhs,(typeof(At23rhs)) d_At23rhs,(typeof(At33rhs)) d_At33rhs,(typeof(B1rhs)) d_B1rhs,(typeof(B2rhs)) d_B2rhs,(typeof(B3rhs)) d_B3rhs,(typeof(beta1rhs)) d_beta1rhs,(typeof(beta2rhs)) d_beta2rhs,(typeof(beta3rhs)) d_beta3rhs,(typeof(gt11rhs)) d_gt11rhs,(typeof(gt12rhs)) d_gt12rhs,(typeof(gt13rhs)) d_gt13rhs,(typeof(gt22rhs)) d_gt22rhs,(typeof(gt23rhs)) d_gt23rhs,(typeof(gt33rhs)) d_gt33rhs,(typeof(phirhs)) d_phirhs,(typeof(trKrhs)) d_trKrhs,(typeof(Xt1rhs)) d_Xt1rhs,(typeof(Xt2rhs)) d_Xt2rhs,(typeof(Xt3rhs)) d_Xt3rhs,





 prms);
//    cutilCheckMsg("failed while updating the velocity");
    CUDA_SAFE_CALL(cudaThreadSynchronize());
    
    
}
