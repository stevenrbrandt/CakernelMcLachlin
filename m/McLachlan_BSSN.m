$Path = Join[$Path, {"~/Calpha/kranc/Tools/CodeGen",
                     "~/Calpha/kranc/Tools/MathematicaMisc"}];

Get["KrancThorn`"];

SetEnhancedTimes[False];
SetSourceLanguage["C"];

(******************************************************************************)
(* Options *)
(******************************************************************************)

createCode[derivOrder_, useGlobalDerivs_, evolutionTimelevels_, addMatter_] :=
Module[{},

prefix = "ML_";
suffix =
  If [useGlobalDerivs, "_MP", ""] <>
  If [derivOrder!=4, "_O" <> ToString[derivOrder], ""] <>
  If [evolutionTimelevels!=3, "_TL" <> ToString[evolutionTimelevels], ""] <>
  If [addMatter==1, "_M", ""];

BSSN = prefix <> "BSSN" <> suffix;

(******************************************************************************)
(* Derivatives *)
(******************************************************************************)

KD = KroneckerDelta;

derivatives =
{
  PDstandardNth[i_]    -> StandardCenteredDifferenceOperator[1,derivOrder/2,i],
  PDstandardNth[i_,i_] -> StandardCenteredDifferenceOperator[2,derivOrder/2,i],
  PDstandardNth[i_,j_] -> StandardCenteredDifferenceOperator[1,derivOrder/2,i] *
                          StandardCenteredDifferenceOperator[1,derivOrder/2,j],
  
(*
  PDPlus    [i_] -> (+1) (-1 + shift[i]^(+1)) / spacing[i],
  PDMinus   [i_] -> (-1) (-1 + shift[i]^(-1)) / spacing[i]
*)
  PDPlus    [i_] -> StandardUpwindDifferenceOperator[1,0,derivOrder/2,i],
  PDMinus   [i_] -> StandardUpwindDifferenceOperator[1,derivOrder/2,0,i]
};

FD = PDstandardNth;

ResetJacobians;
If [useGlobalDerivs,
    DefineJacobian[PD, FD, J, dJ],
    DefineJacobian[PD, FD, KD, Zero3]];



(******************************************************************************)
(* Tensors *)
(******************************************************************************)

(* Register the tensor quantities with the TensorTools package *)
Map [DefineTensor,
     {normal, tangentA, tangentB, dir,
      nn, nu, nlen, nlen2, su, vg,
      xx, rr, th, ph,
      J, dJ,
      g, K, alpha, beta, H, M, detg, gu, G, R, trR, Km, trK,
      phi, gt, At, Xt, Xtn, A, B, Atm, Atu, trA, Ats, trAts, cXt, cS, cA,
      e4phi, em4phi, ddetg, detgt, gtu, ddetgt, dgtu, ddgtu, Gt, Rt, Rphi, gK,
      T00, T0, T, rho, S,
      Psi0re, Psi0im, Psi1re, Psi1im, Psi2re, Psi2im, Psi3re, Psi3im,
      Psi4re, Psi4im,
      er, eth, eph, mm1A, mm1L, mm1, mm2A, mm2B, mm2L, mm2,
      ssA, ssB, ssC, ssL, ss, ss0, tt, ss0, kk, nn, kk0, nn0, mmre, mmim,
      EE, BB}];

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
      gt[la,lb], At[la,lb], Ats[la,lb], Rt[la,lb], Rphi[la,lb], T[la,lb]}];
AssertSymmetricIncreasing [dJ[ua,lb,lc], lb, lc];
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
      dtbetax, dtbetay, dtbetaz,
      eTtt,
      eTtx, eTty, eTtz,
      eTxx, eTxy, eTxz, eTyy, eTyz, eTzz}];

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

pi = N[Pi,40]; 

(******************************************************************************)
(* Groups *)
(******************************************************************************)

SetGroupTimelevels[g_,tl_] = Join[g, {Timelevels -> tl}];

evolvedGroups =
  {SetGroupName [CreateGroupFromTensor [phi      ], prefix <> "log_confac"],
   SetGroupName [CreateGroupFromTensor [gt[la,lb]], prefix <> "metric"    ],
   SetGroupName [CreateGroupFromTensor [Xt[ua]   ], prefix <> "Gamma"     ],
   SetGroupName [CreateGroupFromTensor [trK      ], prefix <> "trace_curv"],
   SetGroupName [CreateGroupFromTensor [At[la,lb]], prefix <> "curv"      ],
   SetGroupName [CreateGroupFromTensor [alpha    ], prefix <> "lapse"     ],
   SetGroupName [CreateGroupFromTensor [A        ], prefix <> "dtlapse"   ],
   SetGroupName [CreateGroupFromTensor [beta[ua] ], prefix <> "shift"     ],
   SetGroupName [CreateGroupFromTensor [B[ua]    ], prefix <> "dtshift"   ]};
evaluatedGroups =
  {SetGroupName [CreateGroupFromTensor [H      ], prefix <> "Ham"],
   SetGroupName [CreateGroupFromTensor [M[la]  ], prefix <> "mom"],
   SetGroupName [CreateGroupFromTensor [cS     ], prefix <> "cons_detg"],
   SetGroupName [CreateGroupFromTensor [cXt[ua]], prefix <> "cons_Gamma"],
   SetGroupName [CreateGroupFromTensor [cA     ], prefix <> "cons_traceA"]};

declaredGroups = Join [evolvedGroups, evaluatedGroups];
declaredGroupNames = Map [First, declaredGroups];



extraGroups =
  {{"ADMBase::metric",   {gxx, gxy, gxz, gyy, gyz, gzz}},
   {"ADMBase::curv",     {kxx, kxy, kxz, kyy, kyz, kzz}},
   {"ADMBase::lapse",    {alp}},
   {"ADMBase::dtlapse",  {dtalp}},
   {"ADMBase::shift",    {betax, betay, betaz}},
   {"ADMBase::dtshift",  {dtbetax, dtbetay, dtbetaz}},
   {"TmunuBase::stress_energy_scalar", {eTtt}},
   {"TmunuBase::stress_energy_vector", {eTtx, eTty, eTtz}},
   {"TmunuBase::stress_energy_tensor", {eTxx, eTxy, eTxz, eTyy, eTyz, eTzz}},
   {"Coordinates::jacobian", {J11, J12, J13, J21, J22, J23, J31, J32, J33}},
   {"Coordinates::jacobian2", {dJ111, dJ112, dJ113, dJ122, dJ123, dJ133,
                               dJ211, dJ212, dJ213, dJ222, dJ223, dJ233,
                               dJ311, dJ312, dJ313, dJ322, dJ323, dJ333}}
};



groups = Join [declaredGroups, extraGroups];

(******************************************************************************)
(* Initial data *)
(******************************************************************************)

initialCalc =
{
  Name -> BSSN <> "_Minkowski",
  Schedule -> {"IN ADMBase_InitialData"},
  ConditionalOnKeyword -> {"my_initial_data", "Minkowski"},
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
};

(*
updateCalc[calc_, name_, value_] :=
  ReplacePart[calc, Position[calc, name][[1]][[1]] -> (name -> value)]

rhs[list_] := Map[#[[1]]&, list]

duplicateID[name_] := ToExpression[ToString[name] <> "copy"]

vars = rhs[Equations /. initialCalc]
newVars = Map[duplicateID, vars]
*)

(* initialCalc = updateCalc[initialCalc, Equations, {phi->1}] *)

(******************************************************************************)
(* Convert from ADMBase *)
(******************************************************************************)

convertFromADMBaseCalc =
{
  Name -> BSSN <> "_convertFromADMBase",
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
    
    phi       -> Log[detg] / 12,
    em4phi    -> Exp[-4 phi],
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
    
    beta1 -> betax,
    beta2 -> betay,
    beta3 -> betaz
  }
};

convertFromADMBaseGammaCalc =
{
  Name -> BSSN <> "_convertFromADMBaseGamma",
  Schedule -> {"AT initial AFTER " <> BSSN <> "_convertFromADMBase"},
  ConditionalOnKeyword -> {"my_initial_data", "ADMBase"},
  Where -> Interior,
  (* should not sync Gamma, since boundary conditions and
     synchronisation are applied later anyway *)
  Shorthands -> {detgt, gtu[ua,ub], Gt[ua,lb,lc]},
  Equations -> 
  {
    detgt        -> 1 (* detgtExpr *),
    gtu[ua,ub]   -> 1/detgt detgtExpr MatrixInverse [gt[ua,ub]],
    Gt[ua,lb,lc] -> 1/2 gtu[ua,ud]
                    (PD[gt[lb,ld],lc] + PD[gt[lc,ld],lb] - PD[gt[lb,lc],ld]),
    Xt[ua] -> gtu[ub,uc] Gt[ua,lb,lc],
    
    A -> - dtalp / (harmonicF alpha^harmonicN) (LapseAdvectionCoeff - 1),
    
    B1 -> 1/ShiftGammaCoeff
          (dtbetax - ShiftAdvectionCoeff beta[ua] PD[beta1,la]),
    B2 -> 1/ShiftGammaCoeff
          (dtbetay - ShiftAdvectionCoeff beta[ua] PD[beta2,la]),
    B3 -> 1/ShiftGammaCoeff
          (dtbetaz - ShiftAdvectionCoeff beta[ua] PD[beta3,la])
  }
};

(******************************************************************************)
(* Convert to ADMBase *)
(******************************************************************************)

convertToADMBaseCalc =
{
  Name -> BSSN <> "_convertToADMBase",
  Schedule -> {"IN " <> BSSN <>"_convertToADMBaseGroup"},
  Where -> Interior,
  Shorthands -> {e4phi, g[la,lb], K[la,lb]},
  Equations -> 
  {
    e4phi    -> Exp[4 phi],
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
    betax    -> beta1,
    betay    -> beta2,
    betaz    -> beta3,
    (* see RHS *)
    dtalp    -> - harmonicF alpha^harmonicN
                  ((1 - LapseAdvectionCoeff) A + LapseAdvectionCoeff trK)
                + LapseAdvectionCoeff beta[ua] PD[alpha,la],
    dtbetax  -> + ShiftGammaCoeff B1
                + ShiftAdvectionCoeff beta[ub] PD[beta[ua],lb],
    dtbetay  -> + ShiftGammaCoeff B2
                + ShiftAdvectionCoeff beta[ub] PD[beta[ua],lb],
    dtbetaz  -> + ShiftGammaCoeff B3
                + ShiftAdvectionCoeff beta[ub] PD[beta[ua],lb]
  }
};

(* TODO: Need to apply a boundary condition to the ADM variables in
   all cases!  *)
(*
boundaryCalcADMBase =
{
  Name -> BSSN <> "_ADMBaseBoundary",
  Schedule -> {"IN " <> BSSN <> "_convertToADMBaseGroup AFTER " <> BSSN <> "_convertToADMBase"},
  ConditionalOnKeyword -> {"my_boundary_condition", "Minkowski"},
  Where -> BoundaryWithGhosts,
  Equations -> 
  {
    gxx     -> 1,
    gxy     -> 0,
    gxz     -> 0,
    gyy     -> 1,
    gyz     -> 0,
    gzz     -> 1,
    kxx     -> 0,
    kxy     -> 0,
    kxz     -> 0,
    kyy     -> 0,
    kyz     -> 0,
    kzz     -> 0,
    alp     -> 1,
    dtalp   -> 0,
    betax   -> 0,
    betay   -> 0,
    betaz   -> 0,
    dtbetax -> 0,
    dtbetay -> 0,
    dtbetaz -> 0
  }
};
*)

boundaryCalcADMBase =
{
  Name -> BSSN <> "_ADMBaseBoundary",
  Schedule -> {"IN " <> BSSN <> "_convertToADMBaseGroup AFTER " <> BSSN <> "_convertToADMBase"},
  Where -> BoundaryWithGhosts,
  Shorthands -> {e4phi, g[la,lb], K[la,lb]},
  Equations -> 
  {
    e4phi    -> Exp[4 phi],
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
    betax    -> beta1,
    betay    -> beta2,
    betaz    -> beta3,
    (* see RHS, but omit derivatives near the boundary *)
    dtalp    -> - harmonicF alpha^harmonicN
                  ((1 - LapseAdvectionCoeff) A + LapseAdvectionCoeff trK),
    dtbetax  -> + ShiftGammaCoeff B1,
    dtbetay  -> + ShiftGammaCoeff B2,
    dtbetaz  -> + ShiftGammaCoeff B3
  }
};

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

evolCalc =
{
  Name -> BSSN <> "_RHS",
  Schedule -> {"IN " <> BSSN <> "_evolCalcGroup"},
  Where -> Interior,
  Shorthands -> {detgt, ddetgt[la], gtu[ua,ub],
                 dgtu[ua,ub,lc], ddgtu[ua,ub,lc,ld], Gt[ua,lb,lc],
                 Xtn[ua], Rt[la,lb], Rphi[la,lb], R[la,lb],
                 Atm[ua,lb], Atu[ua,ub],
                 e4phi, em4phi, g[la,lb], detg,
                 ddetg[la], gu[ua,ub], G[ua,lb,lc], Ats[la,lb], trAts,
                 T00, T0[la], T[la,lb], rho, S[la], trS},
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
    
    (* The conformal connection functions calculated from the conformal metric,
       used instead of Xt where no derivatives of Xt are taken *)
    Xtn[ui]   -> gtu[uj,uk] Gt[ui,lj,lk],

    (* PRD 62, 044034 (2000), eqn. (18) *)
    Rt[li,lj] -> - (1/2) gtu[ul,um] PD[gt[li,lj],ll,lm]
                 + (1/2) gt[lk,li] PD[Xt[uk],lj]
                 + (1/2) gt[lk,lj] PD[Xt[uk],li]
                 + (1/2) Xtn[uk] gt[li,ln] Gt[un,lj,lk]
                 + (1/2) Xtn[uk] gt[lj,ln] Gt[un,li,lk]
                 + gtu[ul,um] (+ Gt[uk,ll,li] gt[lj,ln] Gt[un,lk,lm]
                               + Gt[uk,ll,lj] gt[li,ln] Gt[un,lk,lm]
                               + Gt[uk,li,lm] gt[lk,ln] Gt[un,ll,lj]),
    (* PRD 62, 044034 (2000), eqn. (15) *)
    Rphi[li,lj] -> - 2 CDt[phi,lj,li]
                   - 2 gt[li,lj] gtu[ul,un] CDt[phi,ll,ln]
                   + 4 CDt[phi,li] CDt[phi,lj]
                   - 4 gt[li,lj] gtu[ul,un] CDt[phi,ln] CDt[phi,ll],
    
    Atm[ua,lb] -> gtu[ua,uc] At[lc,lb],
    Atu[ua,ub] -> Atm[ua,lc] gtu[ub,uc],
    
    e4phi       -> Exp[4 phi],
    em4phi      -> 1 / e4phi,
    g[la,lb]    -> e4phi gt[la,lb],
    detg        -> detgExpr,
    (* gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]], *)
    gu[ua,ub]   -> em4phi gtu[ua,ub],
    G[ua,lb,lc] -> Gt[ua,lb,lc]
                   + 2 (KD[ua,lb] PD[phi,lc] + KD[ua,lc] PD[phi,lb] 
                        - gtu[ua,ud] gt[lb,lc] PD[phi,ld]),
    
    R[la,lb] -> Rt[la,lb] + Rphi[la,lb],
    
    (* Matter terms *)
    
    T00 -> addMatter eTtt,
    T01 -> addMatter eTtx,
    T02 -> addMatter eTty,
    T03 -> addMatter eTtz,
    T11 -> addMatter eTxx,
    T12 -> addMatter eTxy,
    T13 -> addMatter eTxz,
    T22 -> addMatter eTyy,
    T23 -> addMatter eTyz,
    T33 -> addMatter eTzz,
    
    (* rho = n^a n^b T_ab *)
    rho -> addMatter
           (1/alpha^2 (T00 - 2 beta[ui] T0[li] + beta[ui] beta[uj] T[li,lj])),
    
    (* S_i = -p^a_i n^b T_ab, where p^a_i = delta^a_i + n^a n_i *)
    S[li] -> addMatter (-1/alpha (T0[li] - beta[uj] T[li,lj])),
    
    (* trS = gamma^ij T_ij  *)
    trS -> addMatter (gu[ui,uj] T[li,lj]),
    
    (* RHS terms *)
    
    (* PRD 62, 044034 (2000), eqn. (10) *)
    (* PRD 67 084023 (2003), eqn. (16) and (23) *)  
    dot[phi]       -> - (1/6) alpha trK 
                      + Lie[phi, beta] + (1/6) PD[beta[ua],la],
    (* PRD 62, 044034 (2000), eqn. (9) *)
    dot[gt[la,lb]] -> - 2 alpha At[la,lb]
                      + Lie[gt[la,lb], beta] - (2/3) gt[la,lb] PD[beta[uc],lc],
    (* PRD 62, 044034 (2000), eqn. (20) *)
    (* PRD 67 084023 (2003), eqn (26) *)    
    dot[Xt[ui]]    -> - 2 Atu[ui,uj] PD[alpha,lj]
                      + 2 alpha (+ Gt[ui,lj,lk] Atu[uk,uj]
                                 - (2/3) gtu[ui,uj] PD[trK,lj]
                                 + 6 Atu[ui,uj] PD[phi,lj])
                      + gtu[uj,ul] PD[beta[ui],lj,ll]
                      + (1/3) gtu[ui,uj] PD[beta[ul],lj,ll]
                      + beta[uj] PD[Xt[ui],lj]
                      - Xtn[uj] PD[beta[ui],lj] 
                      + (2/3) Xtn[ui] PD[beta[uj],lj]
    (* Equation (4.28) in Baumgarte & Shapiro (Phys. Rept. 376 (2003) 41-131) *)
                      + addMatter (- 16 pi alpha gtu[ui,uj] S[lj]),
    
    (* PRD 62, 044034 (2000), eqn. (11) *)
    dot[trK]       -> - gu[ua,ub] CD[alpha,la,lb]
                      + alpha (Atm[ua,lb] Atm[ub,la] + (1/3) trK^2)
                      + Lie[trK, beta]
    (* Equation (4.21) in Baumgarte & Shapiro (Phys. Rept. 376 (2003) 41-131) *)
                      + addMatter (4 pi alpha (rho + trS)),

    (* PRD 62, 044034 (2000), eqn. (12) *)
    (* TODO: Should we use the Hamiltonian constraint to make Rij tracefree? *)
    Ats[la,lb]     -> - CD[alpha,la,lb] + alpha R[la,lb],
    trAts          -> gu[ua,ub] Ats[la,lb],
    dot[At[la,lb]] -> + em4phi (+ Ats[la,lb] - (1/3) g[la,lb] trAts )
                      + alpha (trK At[la,lb] - 2 At[la,lc] Atm[uc,lb])
                      + Lie[At[la,lb], beta] - (2/3) At[la,lb] PD[beta[uc],lc]
    (* Equation (4.23) in Baumgarte & Shapiro (Phys. Rept. 376 (2003) 41-131) *)
                      + addMatter (- em4phi alpha 8 pi
                                     (T[la,lb] - (1/3) g[la,lb] trS)),
    
    (* dot[alpha] -> - harmonicF alpha^harmonicN trK, *)
    (* dot[alpha] -> - harmonicF alpha^harmonicN A + Lie[alpha, beta], *)
    dot[alpha] -> - harmonicF alpha^harmonicN (
                    (1 - LapseAdvectionCoeff) A + LapseAdvectionCoeff trK)
                  + LapseAdvectionCoeff beta[ua] PD[alpha,la],
    dot[A]     -> (1 - LapseAdvectionCoeff) (dot[trK] - AlphaDriver A),
    (* dot[beta[ua]] -> eta Xt[ua], *)
    (* dot[beta[ua]] -> ShiftGammaCoeff alpha^ShiftAlphaPower B[ua], *)

    dot[beta[ua]] -> + ShiftGammaCoeff B[ua]
                     + ShiftAdvectionCoeff beta[ub] PD[beta[ua],lb],

    dot[B[ua]]    -> + dot[Xt[ua]] - BetaDriver B[ua]
                     + ShiftAdvectionCoeff beta[ub] (+ PD[B[ua],lb]
                                                     - PD[Xt[ua],lb])
  }
};

RHSBoundaryCalc =
{
  Name -> BSSN <> "_RHSBoundary",
  Schedule -> {"IN MoL_RHSBoundaries"},
  ConditionalOnKeyword -> {"my_rhs_boundary_condition", "radiative"},
  Where -> Boundary,
  Shorthands -> {detgt, gtu[ua,ub], em4phi, gu[ua,ub],
                 nn[la], nu[ua], nlen, nlen2, su[ua],
                 vg},
  Equations -> 
  {
    detgt      -> 1 (* detgtExpr *),
    gtu[ua,ub] -> 1/detgt detgtExpr MatrixInverse [gt[ua,ub]],
    em4phi     -> Exp[-4 phi],
    gu[ua,ub]  -> em4phi gtu[ua,ub],
    
    nn[la] -> normal[la],
    nu[ua] -> gu[ua,ub] nn[lb],
    nlen2  -> nu[ua] nn[la],
    nlen   -> Sqrt[nlen2],
    su[ua] -> nu[ua] / nlen,
    
    vg -> Sqrt[harmonicF],
    
    (* TODO: Use PDPlus at lower boundary *)
    (* TODO: Use standard PD operators in directions which are not a
       boundary *)
    dot[phi]       -> - vg su[uc] PDMinus[phi      ,lc],
    dot[gt[la,lb]] -> -    su[uc] PDMinus[gt[la,lb],lc],
    dot[trK]       -> - vg su[uc] PDMinus[trK      ,lc],
    dot[At[la,lb]] -> -    su[uc] PDMinus[At[la,lb],lc],
    dot[Xt[ua]]    -> -    su[uc] PDMinus[Xt[ua]   ,lc],
    dot[alpha]     -> - vg su[uc] PDMinus[alpha    ,lc],
    dot[A]         -> - vg su[uc] PDMinus[A        ,lc],
    dot[beta[ua]]  -> -    su[uc] PDMinus[beta[ua] ,lc],
    dot[B[ua]]     -> -    su[uc] PDMinus[B[ua]    ,lc]
  }
};

enforceCalc =
{
  Name -> BSSN <> "_enforce",
  Schedule -> {"IN MoL_PostStep BEFORE " <> BSSN <> "_BoundConds"},
  Shorthands -> {detgt, gtu[ua,ub], trAt},
  Equations -> 
  {
    detgt -> 1 (* detgtExpr *),
    gtu[ua,ub] -> 1/detgt detgtExpr MatrixInverse [gt[ua,ub]],
    
    trAt -> gtu[ua,ub] At[la,lb],
    
    At[la,lb] -> At[la,lb] - (1/3) gt[la,lb] trAt
    
    (* alpha -> Max[alpha, 10^(-10)] *)
  }
};

(******************************************************************************)
(* Boundary conditions *)
(******************************************************************************)

boundaryCalc =
{
  Name -> BSSN <> "_boundary",
  Schedule -> {"IN MoL_PostStep"},
  ConditionalOnKeyword -> {"my_boundary_condition", "Minkowski"},
  Where -> BoundaryWithGhosts,
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
};

(******************************************************************************)
(* Constraint equations *)
(******************************************************************************)

constraintsCalc =
{
  Name -> BSSN <> "_constraints",
  Schedule -> {"IN " <> BSSN <> "_constraintsCalcGroup"},
  Where -> Interior,
  Shorthands -> {detgt, ddetgt[la], gtu[ua,ub], Gt[ua,lb,lc], e4phi, em4phi,
                 g[la,lb], detg, gu[ua,ub], ddetg[la], G[ua,lb,lc],
                 Rt[la,lb], Rphi[la,lb], R[la,lb], trR, Atm[la,lb],
                 gK[la,lb,lc],
                 T00, T0[la], T[la,lb], rho, S[la]},
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
    
    e4phi       -> Exp[4 phi],
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
    
    (* Matter terms *)
    
    T00 -> eTtt,
    T01 -> eTtx,
    T02 -> eTty,
    T03 -> eTtz,
    T11 -> eTxx,
    T12 -> eTxy,
    T13 -> eTxz,
    T22 -> eTyy,
    T23 -> eTyz,
    T33 -> eTzz,
    
    (* rho = n^a n^b T_ab *)
    rho -> 1/alpha^2 (T00 - 2 beta[ui] T0[li] + beta[ui] beta[uj] T[li,lj]),
    
    (* S_i = -p^a_i n^b T_ab, where p^a_i = delta^a_i + n^a n_i *)
    S[li] -> -1/alpha (T0[li] - beta[uj] T[li,lj]),
    
    (* Constraints *)
    
    (* H -> trR - Km[ua,lb] Km[ub,la] + trK^2, *)
    (* PRD 67, 084023 (2003), eqn. (19) *)
    H -> trR - Atm[ua,lb] Atm[ub,la] + (2/3) trK^2 - addMatter 16 pi rho,
    
    (* gK[la,lb,lc] -> CD[K[la,lb],lc], *)
(*    gK[la,lb,lc] -> + 4 e4phi PD[phi,lc] At[la,lb] + e4phi CD[At[la,lb],lc]
                    + (1/3) g[la,lb] PD[trK,lc],

    M[la] -> gu[ub,uc] (gK[lc,la,lb] - gK[lc,lb,la]), *)

    M[li] -> + gtu[uj,uk] (CDt[At[li,lj],lk] + 6 At[li,lj] PD[phi,lk])
             - (2/3) PD[trK,li]
             - addMatter 8 pi S[li],
    (* TODO: use PRD 67, 084023 (2003), eqn. (20) *)
    
    (* det gamma-tilde *)
    cS -> Log[detgt],
    
    (* Gamma constraint *)
    cXt[ua] -> gtu[ub,uc] Gt[ua,lb,lc] - Xt[ua],
    
    (* trace A-tilde *)
    cA -> gtu[ua,ub] At[la,lb]
  }
};

constraintsBoundaryCalc =
{
  Name -> BSSN <> "_constraints_boundary",
  Schedule -> {"IN " <> BSSN <> "_constraintsCalcGroup AFTER " <> BSSN <> "_constraints"},
  Where -> BoundaryWithGhosts,
  Equations -> 
  {
    H     -> 0,
    M[la] -> 0
  }
};

(******************************************************************************)
(* Implementations *)
(******************************************************************************)

inheritedImplementations =
  Join[{"ADMBase"},
       If [addMatter!=0, {"TmunuBase"}, {}],
       If [useGlobalDerivs, {"Coordinates"}, {}]];

(******************************************************************************)
(* Parameters *)
(******************************************************************************)

inheritedKeywordParameters = {};

extendedKeywordParameters =
{
  {
    Name -> "ADMBase::evolution_method",
    AllowedValues -> {BSSN}
  },
  {
    Name -> "ADMBase::lapse_evolution_method",
    AllowedValues -> {BSSN}
  },
  {
    Name -> "ADMBase::shift_evolution_method",
    AllowedValues -> {BSSN}
  }
};

keywordParameters =
{
  {
    Name -> "my_initial_data",
    (* Visibility -> "restricted", *)
    (* Description -> "ddd", *)
    AllowedValues -> {"ADMBase", "Minkowski"},
    Default -> "ADMBase"
  },
  {
    Name -> "my_rhs_boundary_condition",
    (* Visibility -> "restricted", *)
    (* Description -> "ddd", *)
    AllowedValues -> {"none", "radiative"},
    Default -> "none"
  },
  {
    Name -> "my_boundary_condition",
    (* Visibility -> "restricted", *)
    (* Description -> "ddd", *)
    AllowedValues -> {"none", "Minkowski"},
    Default -> "none"
  },
  {
    Name -> "calculate_ADMBase_variables_at",
    Visibility -> "restricted",
    (* Description -> "ddd", *)
    AllowedValues -> {"MoL_PostStep", "CCTK_EVOL", "CCTK_ANALYSIS"},
    Default -> "MoL_PostStep"
  }
};

intParameters =
{
  {
    Name -> harmonicN,
    Description -> "d/dt alpha = - f alpha^n K  (harmonic=2, 1+log=1)",
    Default -> 2
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
  },
  {
    Name -> LapseAdvectionCoeff,
    Description -> "Factor in front of the shift advection terms in 1+log",
    Default -> 1
  },
  {
    Name -> ShiftAdvectionCoeff,
    Description -> "Factor in front of the shift advection terms in gamma driver",
    Default -> 1
  }
};

(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

calculations =
{
  initialCalc,
  convertFromADMBaseCalc,
  convertFromADMBaseGammaCalc,
  evolCalc,
  RHSBoundaryCalc,
  enforceCalc,
  boundaryCalc,
  convertToADMBaseCalc,
  boundaryCalcADMBase,
  constraintsCalc,
  constraintsBoundaryCalc
};

CreateKrancThornTT [groups, ".", BSSN,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  EvolutionTimelevels -> evolutionTimelevels,
  UseLoopControl -> True,
  InheritedImplementations -> inheritedImplementations,
  InheritedKeywordParameters -> inheritedKeywordParameters,
  ExtendedKeywordParameters -> extendedKeywordParameters,
  KeywordParameters -> keywordParameters,
  IntParameters -> intParameters,
  RealParameters -> realParameters
];

];



(******************************************************************************)
(* Options *)
(******************************************************************************)

(* derivative order: 2, 4, 6, 8, ... *)
(* useGlobalDerivs: False or True *)
(* timelevels: 2 or 3
   (keep this at 3; this is better chosen with a run-time parameter) *)
(* matter: 0 or 1
   (matter seems cheap; it should be always enabled) *)

createCode[4, False, 3, 0];
createCode[4, False, 3, 1];
createCode[4, True,  3, 0];
