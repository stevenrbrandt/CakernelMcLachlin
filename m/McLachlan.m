$Path = Join[$Path, {"~/Calpha/Kranc-devel/Tools/CodeGen",
                     "~/Calpha/Kranc-devel/Tools/MathematicaMisc"}];

Get["KrancThorn`"];

SetEnhancedTimes[False];
SetSourceLanguage["C"];

(******************************************************************************)
(* Derivatives *)
(******************************************************************************)

derivOrder = 4;

derivatives =
{
  (*
  PDstandard2nd[i_]     -> StandardCenteredDifferenceOperator[1,1,i],
  PDstandard2nd[i_, i_] -> StandardCenteredDifferenceOperator[2,1,i],
  PDstandard2nd[i_, j_] -> StandardCenteredDifferenceOperator[1,1,i]
                           StandardCenteredDifferenceOperator[1,1,j]
  *)

  (*
  PDstandard4th[i_]     -> StandardCenteredDifferenceOperator[1,2,i],
  PDstandard4th[i_, i_] -> StandardCenteredDifferenceOperator[2,2,i],
  PDstandard4th[i_, j_] -> StandardCenteredDifferenceOperator[1,2,i]
                           StandardCenteredDifferenceOperator[1,2,j]
  *)

  PDstandardNth[i_]     -> StandardCenteredDifferenceOperator[1,derivOrder/2,i],
  PDstandardNth[i_, i_] -> StandardCenteredDifferenceOperator[2,derivOrder/2,i],
  PDstandardNth[i_, j_] -> StandardCenteredDifferenceOperator[1,derivOrder/2,i]
                           StandardCenteredDifferenceOperator[1,derivOrder/2,j]
};

PD = PDstandardNth;

KD = KroneckerDelta;

(******************************************************************************)
(* Tensors *)
(******************************************************************************)

(* Register the tensor quantities with the TensorTools package *)
Map [DefineTensor,
     {g, K, alpha, beta, H, M, detg, gu, G, R, trR, Km, trK,
      phi, gt, At, Xt, A, B, Atm, Atu, trA, cXt, cS, cA,
      e4phi, em4phi, ddetg, detgt, gtu, ddetgt, dgtu, ddgtu, Gt, Rt, Rphi, gK}];

(* NOTE: It seems as if Lie[.,.] did not take these tensor weights
   into account.  Presumably, CD[.,.] and CDt[.,.] don't do this either.  *)
SetTensorAttribute[phi, TensorWeight, +1/6];
SetTensorAttribute[gt,  TensorWeight, -2/3];
SetTensorAttribute[Xt,  TensorWeight, +2/3];
SetTensorAttribute[At,  TensorWeight, -2/3];
SetTensorAttribute[cXt, TensorWeight, +2/3];
SetTensorAttribute[cS,  TensorWeight, +2  ];

Map [AssertSymmetricIncreasing,
     {g[la,lb], K[la,lb], R[la,lb],
      gt[la,lb], At[la,lb], Rt[la,lb], Rphi[la,lb]}];
AssertSymmetricIncreasing [G[ua,lb,lc], lb, lc];
AssertSymmetricIncreasing [Gt[ua,lb,lc], lb, lc];
AssertSymmetricIncreasing [gK[la,lb,lc], la, lb];
Map [AssertSymmetricDecreasing, {gu[ua,ub], gtu[ua,ub], Atu[ua,ub]}];
AssertSymmetricDecreasing [dgtu[ua,ub,lc], ua, ub];
AssertSymmetricDecreasing [ddgtu[ua,ub,lc,ld], ua, ub];
AssertSymmetricIncreasing [ddgtu[ua,ub,lc,ld], lc, ld];

DefineConnection [CD, PD, G];
DefineConnection [CDt, PD, Gt];

Map [DefineTensor,
     {gxx, gxy, gxz, gyy, gyz, gzz,
      kxx, kxy, kxz, kyy, kyz, kzz,
      alp,
      dtalp,
      betax, betay, betaz,
      dtbetax, dtbetay, dtbetaz}];

(******************************************************************************)
(* Expressions *)
(******************************************************************************)

detgExpr  = Det [MatrixOfComponents [g [la,lb]]];
ddetgExpr[la_] =
  Sum [D[Det[MatrixOfComponents[g[la, lb]]], X] PD[X, la],
       {X, Union[Flatten[MatrixOfComponents[g[la, lb]]]]}];

detgtExpr = Det [MatrixOfComponents [gt[la,lb]]];
ddetgtExpr[la_] =
  Sum [D[Det[MatrixOfComponents[gt[la, lb]]], X] PD[X, la],
       {X, Union[Flatten[MatrixOfComponents[gt[la, lb]]]]}];

(******************************************************************************)
(* Groups *)
(******************************************************************************)

evolvedGroups =
  {SetGroupName [CreateGroupFromTensor [g[la,lb]], "ml_metric"],
   SetGroupName [CreateGroupFromTensor [K[la,lb]], "ml_curv"  ],
   SetGroupName [CreateGroupFromTensor [alpha   ], "ml_lapse" ],
   SetGroupName [CreateGroupFromTensor [beta[ua]], "ml_shift" ]};
evaluatedGroups =
  {SetGroupName [CreateGroupFromTensor [H    ], "Ham"],
   SetGroupName [CreateGroupFromTensor [M[la]], "mom"]};

declaredGroups = Join [evolvedGroups, evaluatedGroups];
declaredGroupNames = Map [First, declaredGroups];

evolvedGroupsBSSN =
  {SetGroupName [CreateGroupFromTensor [phi      ], "ML_log_confac"],
   SetGroupName [CreateGroupFromTensor [gt[la,lb]], "ML_metric"    ],
   SetGroupName [CreateGroupFromTensor [Xt[ua]   ], "ML_Gamma"     ],
   SetGroupName [CreateGroupFromTensor [trK      ], "ML_trace_curv"],
   SetGroupName [CreateGroupFromTensor [At[la,lb]], "ML_curv"      ],
   SetGroupName [CreateGroupFromTensor [alpha    ], "ML_lapse"     ],
   SetGroupName [CreateGroupFromTensor [A        ], "ML_dtlapse"   ],
   SetGroupName [CreateGroupFromTensor [beta[ua] ], "ML_shift"     ],
   SetGroupName [CreateGroupFromTensor [B[ua]    ], "ML_dtshift"   ]};
evaluatedGroupsBSSN =
  {SetGroupName [CreateGroupFromTensor [H      ], "Ham"],
   SetGroupName [CreateGroupFromTensor [M[la]  ], "mom"],
   SetGroupName [CreateGroupFromTensor [cS     ], "cons_detg"],
   SetGroupName [CreateGroupFromTensor [cXt[ua]], "cons_Gamma"],
   SetGroupName [CreateGroupFromTensor [cA     ], "cons_traceA"]};

declaredGroupsBSSN = Join [evolvedGroupsBSSN, evaluatedGroupsBSSN];
declaredGroupNamesBSSN = Map [First, declaredGroupsBSSN];



extraGroups =
  {{"ADMBase::metric",   {gxx, gxy, gxz, gyy, gyz, gzz}},
   {"ADMBase::curv",     {kxx, kxy, kxz, kyy, kyz, kzz}},
   {"ADMBase::lapse",    {alp}},
   {"ADMBase::dtlapse",  {dtalp}},
   {"ADMBase::shift",    {betax, betay, betaz}},
   {"ADMBase::dtshift",  {dtbetax, dtbetay, dtbetaz}}};



groups = Join [declaredGroups, extraGroups];
groupsBSSN = Join [declaredGroupsBSSN, extraGroups];

(******************************************************************************)
(* Initial data *)
(******************************************************************************)

initialCalc =
{
  Name -> "ML_ADM_Minkowski",
  Schedule -> {"IN ADMBase_InitialData"},
  ConditionalOnKeyword -> {"my_initial_data", "Minkowski"},
  (* Where -> Boundary, *)
  (* Where -> Interior, *)
  Equations -> 
  {
    g[la,lb] -> KD[la,lb],
    K[la,lb] -> 0,
    alpha    -> 1,
    beta[ua] -> 0
  }
}

initialCalcBSSN =
{
  Name -> "ML_BSSN_Minkowski",
  Schedule -> {"IN ADMBase_InitialData"},
  ConditionalOnKeyword -> {"my_initial_data", "Minkowski"},
  (* Where -> Boundary, *)
  (* Where -> Interior, *)
  Equations -> 
  {
    phi       -> 0,
    gt[la,lb] -> KD[la,lb],
    trK       -> 0,
    At[la,lb] -> 0,
    Xt[ua]    -> 0,
    alpha     -> 1,
    A         -> 0,
    beta[ua]  -> 0,
    B[ua]     -> 0
  }
}

(******************************************************************************)
(* Convert from ADMBase *)
(******************************************************************************)

convertFromADMBaseCalc =
{
  Name -> "ML_ADM_convertFromADMBase",
  Schedule -> {"AT initial AFTER ADMBase_PostInitial"},
  ConditionalOnKeyword -> {"my_initial_data", "ADMBase"},
  Equations -> 
  {
    g11   -> gxx,
    g12   -> gxy,
    g13   -> gxz,
    g22   -> gyy,
    g23   -> gyz,
    g33   -> gzz,
    K11   -> kxx,
    K12   -> kxy,
    K13   -> kxz,
    K22   -> kyy,
    K23   -> kyz,
    K33   -> kzz,
    (* TODO: this is incomplete; it ignores dtalp and dtbeta^i *)
    alpha -> alp,
    beta1 -> betax,
    beta2 -> betay,
    beta3 -> betaz
  }
}

convertFromADMBaseCalcBSSN =
{
  Name -> "ML_BSSN_convertFromADMBase",
  Schedule -> {"AT initial AFTER ADMBase_PostInitial"},
  ConditionalOnKeyword -> {"my_initial_data", "ADMBase"},
  Shorthands -> {g[la,lb], detg, gu[ua,ub], em4phi, K[la,lb]},
  Equations -> 
  {
    g11 -> gxx,
    g12 -> gxy,
    g13 -> gxz,
    g22 -> gyy,
    g23 -> gyz,
    g33 -> gzz,
    
    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    
    phi       -> Log [detg] / 12,
    em4phi    -> Exp [-4 phi],
    gt[la,lb] -> em4phi g[la,lb],
    
    K11 -> kxx,
    K12 -> kxy,
    K13 -> kxz,
    K22 -> kyy,
    K23 -> kyz,
    K33 -> kzz,
    
    trK       -> gu[ua,ub] K[la,lb],
    At[la,lb] -> em4phi (K[la,lb] - (1/3) g[la,lb] trK),
    
    alpha -> alp,
    (* TODO: this is wrong *)
    A     -> dtalp,
    
    beta1 -> betax,
    beta2 -> betay,
    beta3 -> betaz,
    (* TODO: this is wrong *)
    B1    -> dtbetax,
    B2    -> dtbetay,
    B3    -> dtbetaz
  }
}

convertFromADMBaseCalcBSSNGamma =
{
  Name -> "ML_BSSN_convertFromADMBaseGamma",
  Schedule -> {"AT initial AFTER ML_BSSN_convertFromADMBase"},
  ConditionalOnKeyword -> {"my_initial_data", "ADMBase"},
  Where -> Interior,
  Shorthands -> {detgt, gtu[ua,ub], Gt[ua,lb,lc]},
  Equations -> 
  {
    detgt        -> 1 (* detgtExpr *),
    gtu[ua,ub]   -> 1/detgt detgtExpr MatrixInverse [gt[ua,ub]],
    Gt[ua,lb,lc] -> 1/2 gtu[ua,ud]
                    (PD[gt[lb,ld],lc] + PD[gt[lc,ld],lb] - PD[gt[lb,lc],ld]),
    Xt[ua] -> gtu[ub,uc] Gt[ua,lb,lc]
  }
}

(******************************************************************************)
(* Convert to ADMBase *)
(******************************************************************************)

convertToADMBaseCalc =
{
  Name -> "ML_ADM_convertToADMBase",
  Schedule -> {"IN MoL_PostStep AFTER ML_ADM_ApplyBCs"},
  Equations -> 
  {
    gxx     -> g11,
    gxy     -> g12,
    gxz     -> g13,
    gyy     -> g22,
    gyz     -> g23,
    gzz     -> g33,
    kxx     -> K11,
    kxy     -> K12,
    kxz     -> K13,
    kyy     -> K22,
    kyz     -> K23,
    kzz     -> K33,
    (* TODO: this is wrong; it sets dtalp and dtbeta^i incorrectly *)
    alp     -> alpha,
    dtalp   -> 0,
    betax   -> beta1,
    betay   -> beta2,
    betaz   -> beta3,
    dtbetax -> 0,
    dtbetay -> 0,
    dtbetaz -> 0
  }
}

convertToADMBaseCalcBSSN =
{
  Name -> "ML_BSSN_convertToADMBase",
  Schedule -> {"IN MoL_PostStep AFTER ML_BSSN_ApplyBCs AFTER ML_BSSN_enforce"},
  Shorthands -> {e4phi, g[la,lb], K[la,lb]},
  Equations -> 
  {
    e4phi    -> Exp [4 phi],
    g[la,lb] -> e4phi gt[la,lb],
    gxx      -> g11,
    gxy      -> g12,
    gxz      -> g13,
    gyy      -> g22,
    gyz      -> g23,
    gzz      -> g33,
    K[la,lb] -> e4phi At[la,lb] + (1/3) g[la,lb] trK,
    kxx      -> K11,
    kxy      -> K12,
    kxz      -> K13,
    kyy      -> K22,
    kyz      -> K23,
    kzz      -> K33,
    alp      -> alpha,
    (* TODO: this is wrong *)
    (* TODO: rename dtalp->A, dtbeta->B *)
    dtalp    -> A,
    betax    -> beta1,
    betay    -> beta2,
    betaz    -> beta3,
    (* TODO: this is wrong *)
    dtbetax  -> B1,
    dtbetay  -> B2,
    dtbetaz  -> B3
  }
}

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

evolCalc =
{
  Name -> "ML_ADM_RHS",
  Schedule -> {"IN MoL_CalcRHS", "AT analysis"},
  Where -> Interior,
  Shorthands -> {detg, gu[ua,ub], G[ua,lb,lc], R[la,lb], Km[ua,lb], trK},
  Equations -> 
  {
    detg -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    G[ua,lb,lc] -> 1/2 gu[ua,ud]
                   (PD[g[lb,ld],lc] + PD[g[lc,ld],lb] - PD[g[lb,lc],ld]),
    R[la,lb] -> G[u1,l2,la] G[u2,l1,lb] - G[u1,la,lb] G[u2,l1,l2]
                + 1/2 gu[u1,u2] (- PD[g[l1,l2],la,lb] + PD[g[l1,la],l2,lb]
                                 - PD[g[la,lb],l1,l2] + PD[g[l2,lb],l1,la]),
    Km[ua,lb] -> gu[ua,uc] K[lc,lb],
    trK -> Km[ua,la],

    dot[g[la,lb]] -> -2 alpha K[la,lb]
                     + Lie[g[la,lb], beta],
    dot[K[la,lb]] -> - CD[alpha,la,lb]
                     + alpha (+ R[la,lb] + K[la,lb] trK - 2 K[la,lc] Km[uc,lb])
                     + Lie[K[la,lb], beta],
    dot[alpha]    -> 0,
    dot[beta[ua]] -> 0
  }
}

evolCalcBSSN =
{
  Name -> "ML_BSSN_RHS",
  Schedule -> {"IN MoL_CalcRHS", "AT analysis"},
  Where -> Interior,
  Shorthands -> {detgt, ddetgt[la], gtu[ua,ub],
                 dgtu[ua,ub,lc], ddgtu[ua,ub,lc,ld], Gt[ua,lb,lc],
                 Rt[la,lb], Rphi[la,lb], R[la,lb],
                 Atm[ua,lb], Atu[ua,ub],
                 e4phi, em4phi, g[la,lb], detg,
                 ddetg[la], gu[ua,ub], G[ua,lb,lc]},
  Equations -> 
  {
    detgt        -> 1 (* detgtExpr *),
    ddetgt[la]   -> 0 (* ddetgtExpr[la] *),
    
    (* This leads to simpler code... *)
    gtu[ua,ub]   -> 1/detgt detgtExpr MatrixInverse [gt[ua,ub]],
    dgtu[ua,ub,lc] -> - gtu[ua,ud] gtu[ub,ue] PD[gt[ld,le],lc],
    ddgtu[ua,ub,lc,ld] -> - dgtu[ua,ue,ld] gtu[ub,uf] PD[gt[le,lf],lc]
                          - gtu[ua,ue] dgtu[ub,uf,ld] PD[gt[le,lf],lc]
                          - gtu[ua,ue] gtu[ub,uf] PD[gt[le,lf],lc,ld],
    Gt[ua,lb,lc] -> 1/2 gtu[ua,ud]
                    (PD[gt[lb,ld],lc] + PD[gt[lc,ld],lb] - PD[gt[lb,lc],ld]),
    
    (* PRD 62, 044034 (2000), eqn. (18) *)
    Rt[li,lj] -> - (1/2) gtu[ul,um] PD[gt[li,lj],ll,lm]
                 + (1/2) gt[lk,li] PD[Xt[uk],lj]
                 + (1/2) gt[lk,lj] PD[Xt[uk],li]
                 + (1/2) Xt[uk] gt[li,ln] Gt[un,lj,lk]
                 + (1/2) Xt[uk] gt[lj,ln] Gt[un,li,lk]
                 + gtu[ul,um] (+ Gt[uk,ll,li] gt[lj,ln] Gt[un,lk,lm]
                               + Gt[uk,ll,lj] gt[li,ln] Gt[un,lk,lm]
                               + Gt[uk,li,lm] gt[lk,ln] Gt[un,ll,lj]),
(*    Rt[li,lj] -> (1/2) ( - gtu[ul,um] PD[gt[li,lj],ll,lm]
                           + gt[lk,li] PD[Xt[uk],lj] +
                           + gt[lk,lj] PD[Xt[uk],li] 
                           + Xt[uk] gt[li,ln] Gt[un,lj,lk]
                           + Xt[uk] gt[lj,ln] Gt[un,li,lk] )
                 + gtu[ul,um] (+ Gt[uk,ll,li] gt[lj,ln] Gt[un,lk,lm]
                               + Gt[uk,ll,lj] gt[li,ln] Gt[un,lk,lm]
                               + Gt[uk,li,lm] gt[lk,ln] Gt[un,ll,lj]), *)
    (* PRD 62, 044034 (2000), eqn. (15) *)
    (* TODO: Check that CDt takes the tensor weight of phi into account *)
    Rphi[li,lj] -> - 2 CDt[phi,lj,li]
                   - 2 gt[li,lj] gtu[ul,un] CDt[phi,ll,ln]
                   + 4 CDt[phi,li] CDt[phi,lj]
                   - 4 gt[li,lj] gtu[ul,un] CDt[phi,ln] CDt[phi,ll],
    
    Atm[ua,lb] -> gtu[ua,uc] At[lc,lb],
    Atu[ua,ub] -> Atm[ua,lc] gtu[ub,uc],
    
    e4phi       -> Exp [4 phi],
    em4phi      -> 1 / e4phi,
    g[la,lb]    -> e4phi gt[la,lb],
    detg        -> detgExpr,
    (* gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]], *)
    gu[ua,ub]   -> em4phi gtu[ua,ub],
    ddetg[la]   -> 12 detg PD[phi,la],
    G[ua,lb,lc] -> Gt[ua,lb,lc]
                   + 1/(6 detg) ( KD[ua,lb] ddetg[lc] + KD[ua,lc] ddetg[lb]
                                  - gtu[ua,ud] gt[lb,lc] ddetg[ld] ),
    
    R[la,lb] -> Rt[la,lb] + Rphi[la,lb],
    
    (* PRD 62, 044034 (2000), eqn. (10) *)
    (* PRD 67 084023 (2003), eqn. (16) and (23) *)  
    dot[phi]       -> - (1/6) alpha trK 
                      + Lie[phi, beta] + (1/6) PD[beta[ua],la],
    (* PRD 62, 044034 (2000), eqn. (9) *)
    dot[gt[la,lb]] -> - 2 alpha At[la,lb]
                      + Lie[gt[la,lb], beta] - (2/3) gt[la,lb] PD[beta[uc],lc],
    (* PRD 62, 044034 (2000), eqn. (20) *)
    dot[Xt[ui]]    -> - 2 Atu[ui,uj] PD[alpha,lj]
                      + 2 alpha (+ Gt[ui,lj,lk] Atu[uk,uj]
                                 - (2/3) gtu[ui,uj] PD[trK,lj]
                                 + 6 Atu[ui,uj] PD[phi,lj])
                      - (+ (+ PD[beta[ul],lj] dgtu[ui,uj,ll]
                            + beta[ul] ddgtu[ui,uj,ll,lj])
                         - 2 (+ dgtu[um,uj,lj] PD[beta[ui],lm]
                              + dgtu[um,ui,lj] PD[beta[uj],lm]
                              + gtu[um,uj] PD[beta[ui],lm,lj]
                              + gtu[um,ui] PD[beta[uj],lm,lj])
                         + (2/3) (+ dgtu[ui,uj,lj] PD[beta[ul],ll]
                                  + gtu[ui,uj] PD[beta[ul],ll,lj])),
    
(*    dot[Xt[ui]]    -> - 2 Atu[ui,uj] PD[alpha,lj]
                      + 2 alpha (+ Gt[ui,lj,lk] Atu[uk,uj]
                                 - (2/3) gtu[ui,uj] PD[trK,lj]
                                 + 6 Atu[ui,uj] PD[phi,lj])
                      + gtu[uj,ul] PD[beta[ui],lj,ll]
                      + (1/3) gtu[ui,uj] PD[beta[ul],lj,ll]
                      + beta[uj] PD[Xt[ui],lj]
                      + PD[gtu[ul,uj],ll] PD[beta[ui],lj]
                      - (2/3) PD[gtu[ui,uj],lj] PD[beta[ul],ll], *)
    
    (* PRD 62, 044034 (2000), eqn. (11) *)
    dot[trK]       -> - gu[ua,ub] CD[alpha,la,lb]
                      + alpha (Atm[ua,lb] Atm[ub,la] + (1/3) trK^2)
                      + Lie[trK, beta],

    (* PRD 62, 044034 (2000), eqn. (12) *)
    (* TODO: use Hamiltonian constraint to make tracefree *)
    dot[At[la,lb]] -> + em4phi (+ (- CD[alpha,la,lb] + alpha R[la,lb])
                                - (1/3) g[la,lb] gu[uc,ud]
                                        (- CD[alpha,lc,ld] + alpha R[lc,ld]))
                      + alpha (trK At[la,lb] - 2 At[la,lc] Atm[uc,lb])
                      + Lie[At[la,lb], beta] - (2/3) At[la,lb] PD[beta[uc],lc],
    
    (* dot[alpha] -> - harmonicF alpha^harmonicN trK, *)
    dot[alpha] -> - harmonicF alpha^harmonicN A
                       + Lie[alpha, beta],
    dot[A]    -> dot[trK] - AlphaDriver A,
    (* dot[beta[ua]] -> eta Xt[ua], *)
    dot[beta[ua]] -> ShiftGammaCoeff alpha^ShiftAlphaPower B[ua],
    dot[B[ua]]    -> dot[Xt[ua]] - BetaDriver B[ua]
  }
}

enforceCalcBSSN =
{
  Name -> "ML_BSSN_enforce",
  Schedule -> {"IN MoL_PostStep BEFORE ML_BSSN_BoundConds"},
  Shorthands -> {detgt, gtu[ua,ub], trA},
  Equations -> 
  {
    detgt -> 1 (* detgtExpr *),
    gtu[ua,ub] -> 1/detgt detgtExpr MatrixInverse [gt[ua,ub]],
    
    trA -> gtu[ua,ub] At[la,lb],
    
    At[la,lb] -> At[la,lb] - (1/3) gt[la,lb] trA
  }
}

(******************************************************************************)
(* Constraint equations *)
(******************************************************************************)

constraintsCalc =
{
  Name -> "ML_ADM_constraints",
  Schedule -> {"AT analysis"},
  Where -> Interior,
  Shorthands -> {detg, gu[ua,ub], G[ua,lb,lc], R[la,lb], trR, Km[ua,lb], trK},
  Equations -> 
  {
    detg -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse[g[ua,ub]],
    G[ua,lb,lc] -> 1/2 gu[ua,ud]
                   (PD[g[lb,ld],lc] + PD[g[lc,ld],lb] - PD[g[lb,lc],ld]),
    R[la,lb] -> G[u1,l2,la] G[l1,lb,u2] - G[u1,la,lb] G[l1,l2,u2]
                + 1/2 gu[u1,u2] (- PD[g[l1,l2],la,lb] + PD[g[l1,la],l2,lb]
                                 - PD[g[la,lb],l1,l2] + PD[g[l2,lb],l1,la]),
    trR -> R[la,lb] gu[ua,ub],
    Km[ua,lb] -> gu[ua,uc] K[lc,lb],
    trK -> Km[ua,la],

    H -> trR - Km[ua,lb] Km[ub,la] + trK^2,
    M[la] -> gu[ub,uc] (CD[K[lc,la], lb] - CD[K[lc,lb], la])
  }
}

constraintsCalcBSSN =
{
  Name -> "ML_BSSN_constraints",
  Schedule -> {"AT analysis"},
  Where -> Interior,
  Shorthands -> {detgt, ddetgt[la], gtu[ua,ub], Gt[ua,lb,lc], e4phi, em4phi,
                 g[la,lb], detg, gu[ua,ub], ddetg[la], G[ua,lb,lc],
                 Rt[la,lb], Rphi[la,lb], R[la,lb], trR, Atm[la,lb],
                 gK[la,lb,lc]},
  Equations -> 
  {
    detgt        -> 1 (* detgtExpr *),
    ddetgt[la]   -> 0 (* ddetgtExpr[la] *),
    gtu[ua,ub]   -> 1/detgt detgtExpr MatrixInverse [gt[ua,ub]],
    Gt[ua,lb,lc] -> 1/2 gtu[ua,ud]
                    (PD[gt[lb,ld],lc] + PD[gt[lc,ld],lb] - PD[gt[lb,lc],ld]),
    
    (* PRD 62, 044034 (2000), eqn. (18) *)
    (* Note: This differs from the Goddard formulation,
             which is e.g. described in PRD 70 (2004) 124025, eqn. (6).
             Goddard has a Gamma^k replaced by its definition in terms
             of metric derivatives.  *)
    Rt[li,lj] -> - (1/2) gtu[ul,um] PD[gt[li,lj],ll,lm]
                 + (1/2) gt[lk,li] PD[Xt[uk],lj]
                 + (1/2) gt[lk,lj] PD[Xt[uk],li]
                 + (1/2) Xt[uk] gt[li,ln] Gt[un,lj,lk]
                 + (1/2) Xt[uk] gt[lj,ln] Gt[un,li,lk]
                 + gtu[ul,um] (+ Gt[uk,ll,li] gt[lj,ln] Gt[un,lk,lm]
                               + Gt[uk,ll,lj] gt[li,ln] Gt[un,lk,lm]
                               + Gt[uk,li,lm] gt[lk,ln] Gt[un,ll,lj]),

    (* From the long turducken paper.
       This expression seems to give the same result as the one from 044034.  *)
    (* TODO: symmetrise correctly: (ij) = (1/2) [i+j] *)
(*
    Rt[li,lj] -> - (1/2) gtu[uk,ul] PD[gt[li,lj],lk,ll]
                 + gt[lk,li] PD[Xt[uk],lj] + gt[lk,lj] PD[Xt[uk],li]
                 + gt[li,ln] Gt[un,lj,lk] gtu[um,ua] gtu[uk,ub] PD[gt[la,lb],lm]
                 + gt[lj,ln] Gt[un,li,lk] gtu[um,ua] gtu[uk,ub] PD[gt[la,lb],lm]
                 + gtu[ul,us] (+ 2 Gt[uk,ll,li] gt[lj,ln] Gt[un,lk,ls]
                               + 2 Gt[uk,ll,lj] gt[li,ln] Gt[un,lk,ls]
                               + Gt[uk,li,ls] gt[lk,ln] Gt[un,ll,lj]),
*)

    (* Below would be a straightforward calculation,
       without taking any Gamma^i into account.
       This expression gives a different answer!  *)
(*
    Rt[la,lb] -> + Gt[u1,l2,la] Gt[l1,lb,u2] - Gt[u1,la,lb] Gt[l1,l2,u2]
                 + 1/2 gtu[u1,u2] (- PD[gt[l1,l2],la,lb] + PD[gt[l1,la],l2,lb]
                                   - PD[gt[la,lb],l1,l2] + PD[gt[l2,lb],l1,la]),
*)
    (* PRD 62, 044034 (2000), eqn. (15) *)
    Rphi[li,lj] -> - 2 CDt[phi,lj,li]
                   - 2 gt[li,lj] gtu[ul,un] CDt[phi,ll,ln]
                   + 4 CDt[phi,li] CDt[phi,lj]
                   - 4 gt[li,lj] gtu[ul,un] CDt[phi,ln] CDt[phi,ll],
    
    e4phi       -> Exp [4 phi],
    em4phi      -> 1 / e4phi,
    g[la,lb]    -> e4phi gt[la,lb],
    (* detg      -> detgExpr, *)
    (* gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]], *)
    detg        -> e4phi^3,
    gu[ua,ub]   -> em4phi gtu[ua,ub],
    (* ddetg[la] -> PD[e4phi detg,la], *)
    ddetg[la]   -> e4phi ddetgt[la] + 4 detgt e4phi PD[phi,la],
    (* TODO: check this equation, maybe simplify it by omitting ddetg *)
    G[ua,lb,lc] -> Gt[ua,lb,lc]
                   + 1/(2 detg) (+ KD[ua,lb] ddetg[lc] + KD[ua,lc] ddetg[lb]
                                 - (1/3) g[lb,lc] gu[ua,ud] ddetg[ld]),
    
    R[la,lb] -> Rt[la,lb] + Rphi[la,lb],
    trR      -> gu[ua,ub] R[la,lb],
    
    (* K[la,lb] -> e4phi At[la,lb] + (1/3) g[la,lb] trK, *)
    (* Km[ua,lb] -> gu[ua,uc] K[lc,lb], *)
    Atm[ua,lb] -> gtu[ua,uc] At[lc,lb],
    
    (* H -> trR - Km[ua,lb] Km[ub,la] + trK^2, *)
    (* PRD 67, 084023 (2003), eqn. (19) *)
    H -> trR - Atm[ua,lb] Atm[ub,la] + (2/3) trK^2,
    
(*    (* gK[la,lb,lc] -> CD[K[la,lb],lc], *)
    gK[la,lb,lc] -> + 4 e4phi PD[phi,lc] At[la,lb] + e4phi CD[At[la,lb],lc]
                    + (1/3) g[la,lb] PD[trK,lc],

    M[la] -> gu[ub,uc] (gK[lc,la,lb] - gK[lc,lb,la]), *)

    M[li] -> + gtu[uj,uk] ( CDt[At[li,lj],lk] + 6 At[li,lj] PD[phi,lk] )
             - (2/3) PD[trK,li],
    (* TODO: use PRD 67, 084023 (2003), eqn. (20) *)
    
    (* det gamma-tilde *)
    cS -> Log [detgt],
    
    (* Gamma constraint *)
    cXt[ua] -> gtu[ub,uc] Gt[ua,lb,lc] - Xt[ua],
    
    (* trace A-tilde *)
    cA -> gtu[ua,ub] At[la,lb]
  }
}

(******************************************************************************)
(* Implementations *)
(******************************************************************************)

inheritedImplementations = {"ADMBase"};

(******************************************************************************)
(* Parameters *)
(******************************************************************************)

inheritedKeywordParameters = {"ADMBase::initial_data"};

keywordParameters =
{
  {
    Name -> "my_initial_data",
    (* Visibility -> "restricted", *)
    (* Description -> "ddd", *)
    AllowedValues -> {"ADMBase", "Minkowski"},
    Default -> "ADMBase"
  }
};

intParameters =
{
  {
    Name -> harmonicN,
    Description -> "d/dt alpha = - f alpha^n K  (harmonic=1, 1+log=0)",
    Default -> 1
  },
  {
    Name -> ShiftAlphaPower,
    Default -> 0
  }
};

realParameters =
{
  {
    Name -> harmonicF,
    Description -> "d/dt alpha = - f alpha^n K   (harmonic=1, 1+log=2)",
    Default -> 1
  },
  {
    Name -> AlphaDriver,
    Default -> 0
  },
  {
    Name -> ShiftGammaCoeff,
    Default -> 0
  },
  {
    Name -> BetaDriver,
    Default -> 0
  }
};

(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

calculations = 
{
  initialCalc,
  convertFromADMBaseCalc,
  evolCalc,
  convertToADMBaseCalc,
  constraintsCalc
};

CreateKrancThornTT [groups, ".", "ML_ADM",
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  InheritedImplementations -> inheritedImplementations,
  InheritedKeywordParameters -> inheritedKeywordParameters,
  KeywordParameters -> keywordParameters
];

calculationsBSSN = 
{
  initialCalcBSSN,
  convertFromADMBaseCalcBSSN,
  convertFromADMBaseCalcBSSNGamma,
  evolCalcBSSN,
  enforceCalcBSSN,
  convertToADMBaseCalcBSSN,
  constraintsCalcBSSN
};

CreateKrancThornTT [groupsBSSN, ".", "ML_BSSN",
  Calculations -> calculationsBSSN,
  DeclaredGroups -> declaredGroupNamesBSSN,
  PartialDerivatives -> derivatives,
  InheritedImplementations -> inheritedImplementations,
  InheritedKeywordParameters -> inheritedKeywordParameters,
  KeywordParameters -> keywordParameters,
  IntParameters -> intParameters,
  RealParameters -> realParameters
];
