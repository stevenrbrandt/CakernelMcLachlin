$Path = Join[$Path, {"../../../repos/Kranc/Tools/CodeGen",
                     "../../../repos/Kranc/Tools/MathematicaMisc"}];
 
Get["KrancThorn`"];


SetEnhancedTimes[False];
SetSourceLanguage["C"];

(******************************************************************************)
(* Options *)
(******************************************************************************)

(* derivative order: 2, 4, 6, 8, ... *)
derivOrder = 4;

(* useGlobalDerivs: True or False *)
useGlobalDerivs = False;

(* timelevels: 2 or 3
   (keep this at 3; this is better chosen with a run-time parameter) *)
evolutionTimelevels = 3;

(* matter: 0 or 1 *)
addMatter = 0;



prefix = "ML_";
suffix =
  If [useGlobalDerivs, "_MP", ""] <>
  If [derivOrder!=4, "_O" <> ToString[derivOrder], ""] <>
  If [evolutionTimelevels!=3, "_TL" <> ToString[evolutionTimelevels], ""] <>
  If [addMatter!=0, "_M", ""];

ADM = prefix <> "ADM" <> suffix;

(******************************************************************************)
(* Derivatives *)
(******************************************************************************)

KD = KroneckerDelta;

derivatives =
{
  PDstandardNth[i_]     -> StandardCenteredDifferenceOperator[1,derivOrder/2,i],
  PDstandardNth[i_, i_] -> StandardCenteredDifferenceOperator[2,derivOrder/2,i],
  PDstandardNth[i_, j_] -> StandardCenteredDifferenceOperator[1,derivOrder/2,i]
                           StandardCenteredDifferenceOperator[1,derivOrder/2,j]
};

FD = PDstandardNth;

If [useGlobalDerivs,
    DefineJacobian[PD, FD, J, dJ],
    DefineJacobian[PD, FD, KD, Zero3]];



(******************************************************************************)
(* Tensors *)
(******************************************************************************)

(* Register the tensor quantities with the TensorTools package *)
Map [DefineTensor,
     {J, dJ,
      g, K, alpha, beta, H, M, detg, gu, G, R, trR, Km, trK,
      T00, T0, T, rho, S}];

Map [AssertSymmetricIncreasing,
     {g[la,lb], K[la,lb], R[la,lb],
      T[la,lb]}];
AssertSymmetricIncreasing [dJ[ua,lb,lc], lb, lc];
AssertSymmetricIncreasing [G[ua,lb,lc], lb, lc];

DefineConnection [CD, PD, G];

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

pi = N[Pi,40]; 

(******************************************************************************)
(* Groups *)
(******************************************************************************)

SetGroupTimelevels[g_,tl_] = Join[g, {Timelevels -> tl}];

evolvedGroups =
  {SetGroupName [CreateGroupFromTensor [g[la,lb]], prefix <> "metric"],
   SetGroupName [CreateGroupFromTensor [K[la,lb]], prefix <> "curv"  ],
   SetGroupName [CreateGroupFromTensor [alpha   ], prefix <> "lapse" ],
   SetGroupName [CreateGroupFromTensor [beta[ua]], prefix <> "shift" ]};
evaluatedGroups =
  {SetGroupName [CreateGroupFromTensor [H    ], prefix <> "Ham"],
   SetGroupName [CreateGroupFromTensor [M[la]], prefix <> "mom"]};

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
  Name -> ADM <> "_Minkowski",
  Schedule -> {"IN ADMBase_InitialData"},
  ConditionalOnKeyword -> {"my_initial_data", "Minkowski"},
  Equations -> 
  {
    g[la,lb] -> KD[la,lb],
    K[la,lb] -> 0,
    alpha    -> 1,
    beta[ua] -> 0
  }
};

(******************************************************************************)
(* Convert from ADMBase *)
(******************************************************************************)

convertFromADMBaseCalc =
{
  Name -> ADM <> "_convertFromADMBase",
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
};

(******************************************************************************)
(* Convert to ADMBase *)
(******************************************************************************)

convertToADMBaseCalc =
{
  Name -> ADM <> "_convertToADMBase",
  Schedule -> {"IN MoL_PostStep AFTER " <> ADM <> "_ApplyBCs"},
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
};

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

evolCalc =
{
  Name -> ADM <> "_RHS",
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
};

(******************************************************************************)
(* Boundary conditions *)
(******************************************************************************)

boundaryCalc =
{
  Name -> ADM <> "_boundary",
  Schedule -> {"IN MoL_PostStep"},
  ConditionalOnKeyword -> {"my_boundary_condition", "Minkowski"},
  Where -> BoundaryWithGhosts,
  Equations -> 
  {
    g[la,lb] -> KD[la,lb],
    K[la,lb] -> 0,
    alpha    -> 1,
    beta[ua] -> 0
  }
};

(******************************************************************************)
(* Constraint equations *)
(******************************************************************************)

constraintsCalc =
{
  Name -> ADM <> "_constraints",
  Schedule -> {"AT analysis"},
  TriggerGroups -> {prefix <> "Ham", prefix <> "mom"},
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
};

constraintsBoundaryCalc =
{
  Name -> ADM <> "_constraints_boundary",
  Schedule -> {"AT analysis AFTER " <> ADM <> "_constraints"},
  (* TriggerGroups -> {prefix <> "Ham", prefix <> "mom"}, *)
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
    Name -> "my_boundary_condition",
    (* Visibility -> "restricted", *)
    (* Description -> "ddd", *)
    AllowedValues -> {"none", "Minkowski"},
    Default -> "none"
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
  boundaryCalc,
  convertToADMBaseCalc,
  constraintsCalc,
  constraintsBoundaryCalc
};

CreateKrancThornTT [groups, ".", ADM,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  EvolutionTimelevels -> evolutionTimelevels,
  UseLoopControl -> True,
  InheritedImplementations -> inheritedImplementations,
  KeywordParameters -> keywordParameters
];
