#include"CaKernel__gpu_cuda__3dblock__ML_BSSN_Dissipation.h"
#include"CaKernel__ML_BSSN_Dissipation.code"
#include"CaKernel__gpu_cuda__mem.h"



                                       
  

                                       
  #undef EpsDiss 



#include<assert.h>

void CAKERNEL_Launch_ML_BSSN_Dissipation(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

#define TS_ERROR(d)                                                            \
    if ( CAKERNEL_Tile##d <= stncl_##d##n + stncl_##d##p )                     \
      CCTK_VWarn                                                               \
       (CCTK_WARN_ABORT, __LINE__,__FILE__, CCTK_THORNSTRING,                  \
        "Tile size for ML_BSSN_Dissipation along %s axis (%d) too small for stencil %d -- %d\n",\
        #d, CAKERNEL_Tile##d, stncl_##d##n, stncl_##d##p);

    TS_ERROR(x); TS_ERROR(y); TS_ERROR(z);
#   undef TS_ERROR

    int vi = 0;
    
    
      assert((vi = CCTK_VarIndex("ML_BSSN::A"))>=0);
      void * d_A = Device_GetVarI(cctkGH, vi, 0); 
      void * d_A_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::alpha"))>=0);
      void * d_alpha = Device_GetVarI(cctkGH, vi, 0); 
      void * d_alpha_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::alpharhs"))>=0);
      void * d_alpharhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_alpharhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::Arhs"))>=0);
      void * d_Arhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_Arhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At11"))>=0);
      void * d_At11 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At11_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At11rhs"))>=0);
      void * d_At11rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At11rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At12"))>=0);
      void * d_At12 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At12_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At12rhs"))>=0);
      void * d_At12rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At12rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At13"))>=0);
      void * d_At13 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At13_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At13rhs"))>=0);
      void * d_At13rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At13rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At22"))>=0);
      void * d_At22 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At22_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At22rhs"))>=0);
      void * d_At22rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At22rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At23"))>=0);
      void * d_At23 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At23_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At23rhs"))>=0);
      void * d_At23rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At23rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At33"))>=0);
      void * d_At33 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At33_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::At33rhs"))>=0);
      void * d_At33rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_At33rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::B1"))>=0);
      void * d_B1 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_B1_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::B1rhs"))>=0);
      void * d_B1rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_B1rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::B2"))>=0);
      void * d_B2 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_B2_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::B2rhs"))>=0);
      void * d_B2rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_B2rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::B3"))>=0);
      void * d_B3 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_B3_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::B3rhs"))>=0);
      void * d_B3rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_B3rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::beta1"))>=0);
      void * d_beta1 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_beta1_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::beta1rhs"))>=0);
      void * d_beta1rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_beta1rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::beta2"))>=0);
      void * d_beta2 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_beta2_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::beta2rhs"))>=0);
      void * d_beta2rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_beta2rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::beta3"))>=0);
      void * d_beta3 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_beta3_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::beta3rhs"))>=0);
      void * d_beta3rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_beta3rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt11"))>=0);
      void * d_gt11 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt11_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt11rhs"))>=0);
      void * d_gt11rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt11rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt12"))>=0);
      void * d_gt12 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt12_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt12rhs"))>=0);
      void * d_gt12rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt12rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt13"))>=0);
      void * d_gt13 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt13_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt13rhs"))>=0);
      void * d_gt13rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt13rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt22"))>=0);
      void * d_gt22 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt22_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt22rhs"))>=0);
      void * d_gt22rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt22rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt23"))>=0);
      void * d_gt23 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt23_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt23rhs"))>=0);
      void * d_gt23rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt23rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt33"))>=0);
      void * d_gt33 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt33_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::gt33rhs"))>=0);
      void * d_gt33rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_gt33rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::phi"))>=0);
      void * d_phi = Device_GetVarI(cctkGH, vi, 0); 
      void * d_phi_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::phirhs"))>=0);
      void * d_phirhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_phirhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::trK"))>=0);
      void * d_trK = Device_GetVarI(cctkGH, vi, 0); 
      void * d_trK_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::trKrhs"))>=0);
      void * d_trKrhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_trKrhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::Xt1"))>=0);
      void * d_Xt1 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_Xt1_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::Xt1rhs"))>=0);
      void * d_Xt1rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_Xt1rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::Xt2"))>=0);
      void * d_Xt2 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_Xt2_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::Xt2rhs"))>=0);
      void * d_Xt2rhs = Device_GetVarI(cctkGH, vi, 0); 
      void * d_Xt2rhs_out = Device_GetVarI(cctkGH, vi, 1);
    
      assert((vi = CCTK_VarIndex("ML_BSSN::Xt3"))>=0);
      void * d_Xt3 = Device_GetVarI(cctkGH, vi, 0); 
      void * d_Xt3_out = Device_GetVarI(cctkGH, vi, 1);
    
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


    CAKERNEL_ML_BSSN_Dissipation<<<                                                    
 dim3(iDivUp(prms.cagh_ni - stncl_xn - stncl_xp, CAKERNEL_Tilex - stncl_xn - stncl_xp), 
      iDivUp(prms.cagh_nk - stncl_zn - stncl_zp, CAKERNEL_Tilez - stncl_zn - stncl_zp) 
          * blocky),
 dim3(CAKERNEL_Tilex, CAKERNEL_Threadsy, CAKERNEL_Threadsz)>>>(
      

(typeof(A)) d_A,(typeof(alpha)) d_alpha,(typeof(alpharhs)) d_alpharhs,(typeof(Arhs)) d_Arhs,(typeof(At11)) d_At11,(typeof(At11rhs)) d_At11rhs,(typeof(At12)) d_At12,(typeof(At12rhs)) d_At12rhs,(typeof(At13)) d_At13,(typeof(At13rhs)) d_At13rhs,(typeof(At22)) d_At22,(typeof(At22rhs)) d_At22rhs,(typeof(At23)) d_At23,(typeof(At23rhs)) d_At23rhs,(typeof(At33)) d_At33,(typeof(At33rhs)) d_At33rhs,(typeof(B1)) d_B1,(typeof(B1rhs)) d_B1rhs,(typeof(B2)) d_B2,(typeof(B2rhs)) d_B2rhs,(typeof(B3)) d_B3,(typeof(B3rhs)) d_B3rhs,(typeof(beta1)) d_beta1,(typeof(beta1rhs)) d_beta1rhs,(typeof(beta2)) d_beta2,(typeof(beta2rhs)) d_beta2rhs,(typeof(beta3)) d_beta3,(typeof(beta3rhs)) d_beta3rhs,(typeof(gt11)) d_gt11,(typeof(gt11rhs)) d_gt11rhs,(typeof(gt12)) d_gt12,(typeof(gt12rhs)) d_gt12rhs,(typeof(gt13)) d_gt13,(typeof(gt13rhs)) d_gt13rhs,(typeof(gt22)) d_gt22,(typeof(gt22rhs)) d_gt22rhs,(typeof(gt23)) d_gt23,(typeof(gt23rhs)) d_gt23rhs,(typeof(gt33)) d_gt33,(typeof(gt33rhs)) d_gt33rhs,(typeof(phi)) d_phi,(typeof(phirhs)) d_phirhs,(typeof(trK)) d_trK,(typeof(trKrhs)) d_trKrhs,(typeof(Xt1)) d_Xt1,(typeof(Xt1rhs)) d_Xt1rhs,(typeof(Xt2)) d_Xt2,(typeof(Xt2rhs)) d_Xt2rhs,(typeof(Xt3)) d_Xt3,(typeof(Xt3rhs)) d_Xt3rhs,



 d_pars.ptr, 
d_pars.EpsDiss_offset, 
 prms);
//    cutilCheckMsg("failed while updating the velocity");
    CUDA_SAFE_CALL(cudaThreadSynchronize());
    
    
}
