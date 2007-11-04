(* < McLachlan.m /Applications/Mathematica.app/Contents/MacOS/MathKernel | tee McLachlan.out *)

$Path = Join[$Path, {"~/Calpha/Kranc-devel/Tools/CodeGen",
                     "~/Calpha/Kranc-devel/Tools/MathematicaMisc"}];

Get["KrancThorn`"];

SetEnhancedTimes[False];
SetSourceLanguage["C"];

(******************************************************************************)
(* Derivatives *)
(******************************************************************************)

derivatives =
{
  (*
  PDstandard2nd[i_]     -> StandardCenteredDifferenceOperator[1,1,i],
  PDstandard2nd[i_, i_] -> StandardCenteredDifferenceOperator[2,1,i],
  PDstandard2nd[i_, j_] -> StandardCenteredDifferenceOperator[1,1,i]
                           StandardCenteredDifferenceOperator[1,1,j],
  *)

  PDstandard4th[i_]     -> StandardCenteredDifferenceOperator[1,2,i],
  PDstandard4th[i_, i_] -> StandardCenteredDifferenceOperator[2,2,i],
  PDstandard4th[i_, j_] -> StandardCenteredDifferenceOperator[1,2,i]
                           StandardCenteredDifferenceOperator[1,2,j]

  (*
  PDstandard6th[i_]     -> StandardCenteredDifferenceOperator[1,3,i],
  PDstandard6th[i_, i_] -> StandardCenteredDifferenceOperator[2,3,i],
  PDstandard6th[i_, j_] -> StandardCenteredDifferenceOperator[1,3,i]
                           StandardCenteredDifferenceOperator[1,3,j]
  *)
};

PD = PDstandard4th;

KD = KroneckerDelta;

(******************************************************************************)
(* Tensors *)
(******************************************************************************)

(* Register the tensor quantities with the TensorTools package *)
Map [DefineTensor,
     {g, K, alpha, beta, H, M, detg, gu, G, R, trR, Km, trK,
      phi, gt, At, Xt, dtalpha, dtbeta, trA, cXt, cS, cA,
      e4phi, em4phi, gtu, ddetg, ddetgt, Gt, Rt, Rphi, gK}];
(* SetTensorAttribute[g, TensorWeight, 0]; *)

Map [AssertSymmetricDecreasing,
     {g[la,lb], K[la,lb], R[la,lb],
      gt[la,lb], At[la,lb], Rt[la,lb], Rphi[la,lb]}];
AssertSymmetricDecreasing [G[ua,lb,lc], lb, lc];
AssertSymmetricDecreasing [Gt[ua,lb,lc], lb, lc];
AssertSymmetricDecreasing [gK[la,lb,lc], la, lb];
Map [AssertSymmetricIncreasing, {gu[ua,ub]}];
Map [AssertSymmetricIncreasing, {gtu[ua,ub]}];

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
  {SetGroupName [CreateGroupFromTensor [g[la,lb]], "metric"],
   SetGroupName [CreateGroupFromTensor [K[la,lb]], "curv"  ],
   SetGroupName [CreateGroupFromTensor [alpha   ], "lapse" ],
   SetGroupName [CreateGroupFromTensor [beta[ua]], "shift" ]};
evaluatedGroups =
  {SetGroupName [CreateGroupFromTensor [H    ], "Ham"],
   SetGroupName [CreateGroupFromTensor [M[la]], "mom"]};

declaredGroups = Join [evolvedGroups, evaluatedGroups];
declaredGroupNames = Map [First, declaredGroups];

evolvedGroupsBSSN =
  {SetGroupName [CreateGroupFromTensor [phi       ], "log_confac"],
   SetGroupName [CreateGroupFromTensor [gt[la,lb] ], "metric"    ],
   SetGroupName [CreateGroupFromTensor [Xt[ua]    ], "Gamma"     ],
   SetGroupName [CreateGroupFromTensor [trK       ], "trace_curv"],
   SetGroupName [CreateGroupFromTensor [At[la,lb] ], "curv"      ],
   SetGroupName [CreateGroupFromTensor [alpha     ], "lapse"     ],
   SetGroupName [CreateGroupFromTensor [dtalpha   ], "dtlapse"   ],
   SetGroupName [CreateGroupFromTensor [beta[ua]  ], "shift"     ],
   SetGroupName [CreateGroupFromTensor [dtbeta[ua]], "dtshift"   ]};
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
    phi        -> 0,
    gt[la,lb]  -> KD[la,lb],
    trK        -> 0,
    At[la,lb]  -> 0,
    Xt[ua]     -> 0,
    alpha      -> 1,
    dtalpha    -> 0,
    beta[ua]   -> 0,
    dtbeta[ua] -> 0
  }
}

(******************************************************************************)
(* Convert from ADMBase *)
(******************************************************************************)

convertFromADMBaseCalc =
{
  Name -> "ML_ADM_convertFromADMBase",
  Schedule -> {"AT initial", "AFTER ADMBase_PostInitial"},
  ConditionalOnKeyword -> {"my_initial_data", "ADMBase"},
  Equations -> 
  {
    g11   -> gxx,
    g21   -> gxy,
    g31   -> gxz,
    g22   -> gyy,
    g32   -> gyz,
    g33   -> gzz,
    K11   -> kxx,
    K21   -> kxy,
    K31   -> kxz,
    K22   -> kyy,
    K32   -> kyz,
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
  Schedule -> {"AT initial", "AFTER ADMBase_PostInitial"},
  ConditionalOnKeyword -> {"my_initial_data", "ADMBase"},
  Shorthands -> {g[la,lb], detg, gu[ua,ub], em4phi, K[la,lb], Km[ua,lb]},
  Equations -> 
  {
    g11 -> gxx,
    g21 -> gxy,
    g31 -> gxz,
    g22 -> gyy,
    g32 -> gyz,
    g33 -> gzz,
    
    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    
    (* em4phi -> Exp [-4 phi], *)
    em4phi    -> 1 / detg^3,
    phi        -> Log [detg] / 12,
    gt[la,lb]  -> em4phi g[la,lb],
    
    K11 -> kxx,
    K21 -> kxy,
    K31 -> kxz,
    K22 -> kyy,
    K32 -> kyz,
    K33 -> kzz,
    
    Km[ua,lb] -> gu[ua,uc] K[lc,lb],
    trK       -> Km[ua,la],
    At[la,lb] -> em4phi (K[la,lb] - (1/3) g[la,lb] trK),
    
    alpha   -> alp,
    dtalpha -> dtalp,
    
    beta1   -> betax,
    beta2   -> betay,
    beta3   -> betaz,
    (* TODO: this is wrong *)
    dtbeta1 -> dtbetax,
    dtbeta2 -> dtbetay,
    dtbeta3 -> dtbetaz
  }
}

convertFromADMBaseCalcBSSNGamma =
{
  Name -> "ML_BSSN_convertFromADMBaseGamma",
  Schedule -> {"AT initial", "AFTER ML_BSSN_convertFromADMBase"},
  ConditionalOnKeyword -> {"my_initial_data", "ADMBase"},
  Where -> Interior,
  Shorthands -> {detgt, gtu[ua,ub], Gt[ua,lb,lc]},
  Equations -> 
  {
    detgt        -> detgtExpr,
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
  Schedule -> {"IN MoL_PostStep", "AFTER ADM_ApplyBoundConds"},
  Equations -> 
  {
    gxx     -> g11,
    gxy     -> g21,
    gxz     -> g31,
    gyy     -> g22,
    gyz     -> g32,
    gzz     -> g33,
    kxx     -> K11,
    kxy     -> K21,
    kxz     -> K31,
    kyy     -> K22,
    kyz     -> K32,
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
  Schedule -> {"IN MoL_PostStep", "AFTER ML_ADM_ApplyBoundConds"},
  Shorthands -> {e4phi, g[la,lb], K[la,lb]},
  Equations -> 
  {
    e4phi    -> Exp [4 phi],
    g[la,lb] -> e4phi gt[la,lb],
    gxx      -> g11,
    gxy      -> g21,
    gxz      -> g31,
    gyy      -> g22,
    gyz      -> g32,
    gzz      -> g33,
    K[la,lb] -> e4phi At[la,lb] + (1/3) g[la,lb] trK,
    kxx      -> K11,
    kxy      -> K21,
    kxz      -> K31,
    kyy      -> K22,
    kyz      -> K32,
    kzz      -> K33,
    alp      -> alpha,
    dtalp    -> dtalpha,
    betax    -> beta1,
    betay    -> beta2,
    betaz    -> beta3,
    dtbetax  -> dtbeta1,
    dtbetay  -> dtbeta2,
    dtbetaz  -> dtbeta3
  }
}

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

evolCalc =
{
  Name -> "ML_ADM_RHS",
  Schedule -> {"IN MoL_CalcRHS"},
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

enforceCalcBSSN =
{
  Name -> "ML_BSSN_enforce",
  Schedule -> {"IN MoL_PostStep"},
  Shorthands -> {detgt, gtu[ua,ub], trA},
  Equations -> 
  {
    detgt -> detgtExpr,
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
  Shorthands -> {detgt, ddetgt[la], gtu[ua,ub], Gt[ua,lb,lc], e4phi,
                 g[la,lb], detg, gu[ua,ub], ddetg[la], G[ua,lb,lc],
                 Rt[la,lb], Rphi[la,lb], R[la,lb], trR,
                 K[la,lb], Km[la,lb], gK[la,lb,lc]},
  Equations -> 
  {
    detgt -> detgtExpr,
    ddetgt[la] -> ddetgtExpr[la],
    gtu[ua,ub]   -> 1/detgt detgtExpr MatrixInverse [gt[ua,ub]],
    Gt[ua,lb,lc] -> 1/2 gtu[ua,ud]
                    (PD[gt[lb,ld],lc] + PD[gt[lc,ld],lb] - PD[gt[lb,lc],ld]),
    
    e4phi     -> Exp [4 phi],
    g[la,lb]  -> e4phi gt[la,lb],
    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    (* ddetg[la] -> PD[e4phi detg,la], *)
    ddetg[la] -> e4phi ddetgt[la] + 4 detgt e4phi PD[phi,la],
    G[ua,lb,lc] -> Gt[ua,lb,lc]
                   + 1/(2 detg) (+ KD[ua,lb] ddetg[lc] + KD[ua,lc] ddetg[lb]
                                 - (1/3) g[lb,lc] gu[ua,ud] ddetg[ld]),

    (* PRD 62, 044034 (2000), eqn. (18) *)
    Rt[li,lj] -> - (1/2) gtu[ul,um] PD[gt[li,lj],ll,lm]
                 + gt[lk,li] PD[Xt[uk],lj] + gt[lk,lj] PD[Xt[uk],li]
                 + Xt[uk] gt[li,ln] Gt[un,lj,lk] + Xt[uk] gt[lj,ln] Gt[un,li,lk]
                 + gtu[ul,um] (+ 2 Gt[uk,ll,li] gt[lj,ln] Gt[un,lk,lm]
                               + 2 Gt[uk,ll,lj] gt[li,ln] Gt[un,lk,lm]
                               + Gt[uk,li,lm] gt[lk,ln] Gt[un,ll,lj]),
    (* PRD 62, 044034 (2000), eqn. (15) *)
    (* TODO: Check that CDt takes the tensor weight of phi into account *)
    Rphi[li,lj] -> - 2 CDt[phi,lj,li]
                   - 2 gt[li,lj] gtu[ul,un] CDt[phi,ll,ln]
                   + 4 CDt[phi,li] CDt[phi,lj]
                   - 4 gt[li,lj] gtu[ul,un] CDt[phi,ln] CDt[phi,ll],
    
    R[la,lb] -> Rt[la,lb] + Rphi[la,lb],
    trR -> gu[ua,ub] R[la,lb],
    
    K[la,lb] -> e4phi At[la,lb] + (1/3) g[la,lb] trK,
    Km[ua,lb] -> gu[ua,uc] K[lc,lb],
    
    H -> trR - Km[ua,lb] Km[ub,la] + trK^2,
    
    (* gK[la,lb,lc] -> CD[K[la,lb],lc], *)
    gK[la,lb,lc] -> + 4 e4phi PD[phi,lc] At[la,lb] + e4phi CD[At[la,lb],lc]
                    + (1/3) g[la,lb] PD[trK,lc],
    M[la] -> gu[ub,uc] (gK[lc,la,lb] - gK[lc,lb,la]),
    
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
{{
  Name -> "my_initial_data",
  (* Visibility -> "restricted", *)
  (* Description -> "ddd", *)
  AllowedValues -> {"ADMBase", "Minkowski"},
  Default -> "ADMBase"
}};

realParameters = {};

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
  KeywordParameters -> keywordParameters,
  RealParameters -> realParameters
];

calculationsBSSN = 
{
  initialCalcBSSN,
  convertFromADMBaseCalcBSSN,
  convertFromADMBaseCalcBSSNGamma,
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
  RealParameters -> realParameters
];
