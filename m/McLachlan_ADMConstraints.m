$Path = Join[$Path, {"../../../repos/Kranc/Tools/CodeGen",
                     "../../../repos/Kranc/Tools/MathematicaMisc"}];

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
  ""
  <> If [useGlobalDerivs, "_MP", ""]
  <> If [derivOrder!=4, "_O" <> ToString[derivOrder], ""]
  (* <> If [evolutionTimelevels!=3, "_TL" <> ToString[evolutionTimelevels], ""] *)
  (* <> If [addMatter==1, "_M", ""] *)
  ;

ADMConstraints = prefix <> "ADMConstraints" <> suffix;

(******************************************************************************)
(* Derivatives *)
(******************************************************************************)

KD = KroneckerDelta;

derivatives =
{
  PDstandardNth[i_]    -> StandardCenteredDifferenceOperator[1,derivOrder/2,i],
  PDstandardNth[i_,i_] -> StandardCenteredDifferenceOperator[2,derivOrder/2,i],
  PDstandardNth[i_,j_] -> StandardCenteredDifferenceOperator[1,derivOrder/2,i] *
                          StandardCenteredDifferenceOperator[1,derivOrder/2,j]
};

FD  = PDstandardNth;

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
      J, dJ,
      g, K, alpha, beta, detg, gu, G, R, trR, Km, trK,
      H, M,
      T00, T0, T, rho, S,
      x, y, z, r}];

Map [AssertSymmetricIncreasing,
     {g[la,lb], K[la,lb], R[la,lb], T[la,lb]}];
AssertSymmetricIncreasing [dJ[ua,lb,lc], lb, lc];
AssertSymmetricIncreasing [G[ua,lb,lc], lb, lc];
Map [AssertSymmetricDecreasing,
     {gu[ua,ub]}];

DefineConnection [CD, PD, G];

(* Use the CartGrid3D variable names *)
x1=x; x2=y; x3=z;

(* Use the ADMBase variable names *)
g11=gxx; g12=gxy; g22=gyy; g13=gxz; g23=gyz; g33=gzz;
K11=kxx; K12=kxy; K22=kyy; K13=kxz; K23=kyz; K33=kzz;
alpha=alp;
beta1=betax; beta2=betay; beta3=betaz;

(* Use the TmunuBase variable names *)
T00=eTtt;
T01=eTtx; T02=eTty; T03=eTtz;
T11=eTxx; T12=eTxy; T22=eTyy; T13=eTxz; T23=eTyz; T33=eTzz;

(******************************************************************************)
(* Expressions *)
(******************************************************************************)

detgExpr = Det [MatrixOfComponents [g [la,lb]]];
ddetgExpr[la_] =
  Sum [D[Det[MatrixOfComponents[g[la, lb]]], X] PD[X, la],
       {X, Union[Flatten[MatrixOfComponents[g[la, lb]]]]}];

pi = N[Pi,40]; 

(******************************************************************************)
(* Groups *)
(******************************************************************************)

SetGroupTimelevels[g_,tl_] = Join[g, {Timelevels -> tl}];

evolvedGroups = {};
evaluatedGroups =
  {SetGroupName [CreateGroupFromTensor [H    ], prefix <> "Ham"],
   SetGroupName [CreateGroupFromTensor [M[la]], prefix <> "mom"]};
(*
evolvedGroups =
  {SetGroupName [CreateGroupFromTensor [H    ], prefix <> "Ham"],
   SetGroupName [CreateGroupFromTensor [M[la]], prefix <> "mom"]};
evaluatedGroups = {};
*)

declaredGroups = Join [evolvedGroups, evaluatedGroups];
declaredGroups = Map [AddGroupExtra [#, Timelevels -> evolutionTimelevels] &, declaredGroups];

declaredGroupNames = Map [First, declaredGroups];



extraGroups =
  {{"Grid::coordinates", {x, y, z, r}},
   {"ADMBase::metric",  {gxx, gxy, gxz, gyy, gyz, gzz}},
   {"ADMBase::curv",    {kxx, kxy, kxz, kyy, kyz, kzz}},
   {"ADMBase::lapse",   {alp}},
   {"ADMBase::dtlapse", {dtalp}},
   {"ADMBase::shift",   {betax, betay, betaz}},
   {"ADMBase::dtshift", {dtbetax, dtbetay, dtbetaz}},
   {"TmunuBase::stress_energy_scalar", {eTtt}},
   {"TmunuBase::stress_energy_vector", {eTtx, eTty, eTtz}},
   {"TmunuBase::stress_energy_tensor", {eTxx, eTxy, eTxz, eTyy, eTyz, eTzz}},
   {"Coordinates::jacobian",  {J11, J12, J13, J21, J22, J23, J31, J32, J33}},
   {"Coordinates::jacobian2", {dJ111, dJ112, dJ113, dJ122, dJ123, dJ133,
                               dJ211, dJ212, dJ213, dJ222, dJ223, dJ233,
                               dJ311, dJ312, dJ313, dJ322, dJ323, dJ333}}
};



groups = Join [declaredGroups, extraGroups];

(******************************************************************************)
(* Constraint equations *)
(******************************************************************************)

ADMConstraintsCalc =
{
  Name -> ADMConstraints,
  Schedule -> Automatic,
  Where -> Interior,
  Shorthands -> {detg, gu[ua,ub], G[ua,lb,lc],
                 R[la,lb], trR, Km[la,lb], trK,
                 rho, S[la]},
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

    (* Matter terms *)

    (* rho = n^a n^b T_ab *)
    rho -> 1/alpha^2 (T00 - 2 beta[ui] T0[li] + beta[ui] beta[uj] T[li,lj]),

    (* S_i = -p^a_i n^b T_ab, where p^a_i = delta^a_i + n^a n_i *)
    S[li] -> -1/alpha (T0[li] - beta[uj] T[li,lj]),

    (* ADM constraints *)

    H     -> + trR - Km[ua,lb] Km[ub,la] + trK^2
             - addMatter 16 pi rho,
    M[la] -> + gu[ub,uc] (CD[K[lc,la], lb] - CD[K[lc,lb], la])
             - addMatter 8 pi S[la]
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

(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

calculations =
{
  ADMConstraintsCalc
};

CreateKrancThornTT [groups, ".", ADMConstraints,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  EvolutionTimelevels -> evolutionTimelevels,
  UseLoopControl -> True,
  InheritedImplementations -> inheritedImplementations
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

createCode[2, False, 3, 1];
createCode[4, False, 3, 1];
createCode[4, True,  3, 1];
