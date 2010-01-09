$Path = Join[$Path, {"~/Calpha/kranc/Tools/CodeGen",
                     "~/Calpha/kranc/Tools/MathematicaMisc"}];

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

WeylScalars = prefix <> "WeylScalars" <> suffix;

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
     {xx, rr, th, ph,
      J, dJ,
      g, K, alpha, beta, H, M, detg, gu, G, R, trR, Km, trK,
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

SetTensorAttribute[Psi0re, TensorManualCartesianParities, {+1,+1,+1}];
SetTensorAttribute[Psi2re, TensorManualCartesianParities, {+1,+1,+1}];
SetTensorAttribute[Psi4re, TensorManualCartesianParities, {+1,+1,+1}];
SetTensorAttribute[Psi0im, TensorManualCartesianParities, {-1,-1,-1}];
SetTensorAttribute[Psi2im, TensorManualCartesianParities, {-1,-1,-1}];
SetTensorAttribute[Psi4im, TensorManualCartesianParities, {-1,-1,-1}];
SetTensorAttribute[Psi1re, TensorManualCartesianParities, {+1,+1,-1}];
SetTensorAttribute[Psi3re, TensorManualCartesianParities, {-1,-1,+1}];
SetTensorAttribute[Psi1im, TensorManualCartesianParities, {-1,-1,+1}];
SetTensorAttribute[Psi3im, TensorManualCartesianParities, {-1,-1,+1}];

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
      dtbetax, dtbetay, dtbetaz}];

(******************************************************************************)
(* Expressions *)
(******************************************************************************)

pi = N[Pi,40]; 

(******************************************************************************)
(* Groups *)
(******************************************************************************)

SetGroupTimelevels[g_,tl_] = Join[g, {Timelevels -> tl}];

WeylScalarsVars = {Psi0re, Psi0im, Psi1re, Psi1im, Psi2re, Psi2im,
                   Psi3re, Psi3im, Psi4re, Psi4im};
WeylScalarsGroups = Map[AddGroupExtra[#, Timelevels -> evolutionTimelevels] &,
                        Map[CreateGroupFromTensor,
                            WeylScalarsVars]];

evolvedGroupsWeylScalars = {};
evaluatedGroupsWeylScalars = WeylScalarsGroups;

declaredGroupsWeylScalars =
  Join [evolvedGroupsWeylScalars, evaluatedGroupsWeylScalars];
declaredGroupNamesWeylScalars = Map [First, declaredGroupsWeylScalars];



extraGroups =
  {{"ADMBase::metric",   {gxx, gxy, gxz, gyy, gyz, gzz}},
   {"ADMBase::curv",     {kxx, kxy, kxz, kyy, kyz, kzz}},
   {"ADMBase::lapse",    {alp}},
   {"ADMBase::dtlapse",  {dtalp}},
   {"ADMBase::shift",    {betax, betay, betaz}},
   {"ADMBase::dtshift",  {dtbetax, dtbetay, dtbetaz}},
   {"Coordinates::jacobian", {J11, J12, J13, J21, J22, J23, J31, J32, J33}},
   {"Coordinates::jacobian2", {dJ111, dJ112, dJ113, dJ122, dJ123, dJ133,
                               dJ211, dJ212, dJ213, dJ222, dJ223, dJ233,
                               dJ311, dJ312, dJ313, dJ322, dJ323, dJ333}}
};



groupsWeylScalars = Join [declaredGroupsWeylScalars, extraGroups];

(******************************************************************************)
(* Weyl scalars *)
(******************************************************************************)

WeylScalarsCalc =
{
  Name -> WeylScalars,
  Schedule -> {"IN " <> WeylScalars <> "_Group"},
  Where -> Interior,
  Shorthands -> {xx[ua], rho, rr, th, ph,
                 g[la,lb], detg, gu[ua,ub], G[ua,lb,lc], R[la,lb],
                 K[la,lb], Km[ua,lb], trK,
                 alpha, beta[ua],
                 er[ua], eth[ua], eph[ua],
                 mm1A[ua], mm1L, mm1[ua],
                 mm2A[ua], mm2B[ua], mm2L, mm2[ua],
                 ssA[ua], ssB[ua], ssC[ua], ssL, ss[ua], ss0,
                 tt[ua], tt0,
                 kk[ua], nn[ua], kk0, nn0, mmre[ua], mmim[ua],
                 EE[la,lb], BB[la,lb]},
  Equations ->
  {
    (* current position *)
    xx1 -> x,                      (* r sin theta cos phi *)
    xx2 -> y,                      (* r sin theta sin phi *)
    xx3 -> z,                      (* r cos theta         *)
    rho -> Sqrt [xx1^2 + xx2^2],   (* r sin theta         *)
    rr  -> r,
    th  -> ArcCos [xx3 / rr],
    ph  -> ArcTan [xx1, xx2],
    
    (* metric *)
    g11 -> gxx,
    g12 -> gxy,
    g13 -> gxz,
    g22 -> gyy,
    g23 -> gyz,
    g33 -> gzz,
    
    detg -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    G[ua,lb,lc] -> 1/2 gu[ua,ud]
                   (PD[g[lb,ld],lc] + PD[g[lc,ld],lb] - PD[g[lb,lc],ld]),
    R[la,lb] -> G[u1,l2,la] G[u2,l1,lb] - G[u1,la,lb] G[u2,l1,l2]
                + 1/2 gu[u1,u2] (- PD[g[l1,l2],la,lb] + PD[g[l1,la],l2,lb]
                                 - PD[g[la,lb],l1,l2] + PD[g[l2,lb],l1,la]),
    
    K11 -> Kxx,
    K12 -> Kxy,
    K13 -> Kxz,
    K22 -> Kyy,
    K23 -> Kyz,
    K33 -> Kzz,
    Km[ua,lb] -> gu[ua,uc] K[lc,lb],
    trK -> Km[ua,la],
    
    alpha -> alp,
    beta1 -> betax,
    beta2 -> betay,
    beta3 -> betaz,
    
    (* build a tetrad *)
    er1 -> xx1 / r,
    er2 -> xx2 / r,
    er3 -> xx3 / r,
    
    eth1 ->   xx3 xx1 / (rho r),
    eth2 ->   xx3 xx2 / (rho r),
    eth3 -> - rho / r,
    
    eph1 -> - xx2 / rho,
    eph2 ->   xx1 / rho,
    eph3 ->   0,
    
    mm1A[ua] -> eph[ua],
    mm1L     -> mm1A[ua] mm1A[ub] g[la,lb],
    mm1[ua]  -> mm1A[ua] / Sqrt[mm1L],
    
    mm2A[ua] -> eth[ua],
    mm2B[ua] -> mm2A[ua] - mm1[ua] mm1[ub] mm2A[uc] g[lb,lc],
    mm2L     -> mm2B[ua] mm2B[ub] g[la,lb],
    mm2[ua]  -> mm2B[ua] / Sqrt[mm2L],
    
    ssA[ua] -> er[ua],
    ssB[ua] -> ssA[ua] - mm1[ua] mm1[ub] ssA[uc] g[lb,lc],
    ssC[ua] -> ssB[ua] - mm2[ua] mm2[ub] ssB[uc] g[lb,lc],
    ssL     -> ssC[ua] ssC[ub] g[la,lb],
    ss[ua]  -> ssB[ua] / Sqrt[ssL],
    ss0     -> 0,
    
    tt[ua] -> - beta[ua] / alpha,
    tt0    -> 1 / alpha,
    
    kk[ua]  -> (tt[ua] + ss[ua]) / Sqrt[2],
    nn[ua]  -> (tt[ua] - ss[ua]) / Sqrt[2],
    kk0     -> (tt0 + ss0) / Sqrt[2],
    nn0     -> (tt0 - ss0) / Sqrt[2],
    mmre[ua] -> mm1[ua] / Sqrt[2],
    mmim[ua] -> mm2[ua] / Sqrt[2],
    
    (* Weyl tensor *)
    (* PRD 72 024013 (22) *)
    EE[la,lb] -> R[la,lb] + trK K[la,lb] - K[la,lc] Km[uc,lb] - 1/2 T[la,lb],
    BB[la,lb] -> - Eps[la,uc,ud] CD[K[ld,lb],lc],
    
    (* C[la,lb] -> E[la,lb] + I B[la,lb] *)
    
    (* Weyl scalars *)
    Psi0re -> 0,
    Psi0im -> 0,
    Psi1re -> 0,
    Psi1im -> 0,
    Psi2re -> 0,
    Psi2im -> 0,
    Psi3re -> 0,
    Psi3im -> 0,
    Psi4re -> 0,
    Psi4im -> 0
  }
}

WeylScalarsBoundaryCalc =
{
  Name -> WeylScalars <> "_boundary",
  Schedule -> {"IN " <> WeylScalars <> "_Group AFTER " <> WeylScalars},
  Where -> BoundaryWithGhosts,
  Equations -> 
  {
    Psi0re -> 0,
    Psi0im -> 0,
    Psi1re -> 0,
    Psi1im -> 0,
    Psi2re -> 0,
    Psi2im -> 0,
    Psi3re -> 0,
    Psi3im -> 0,
    Psi4re -> 0,
    Psi4im -> 0
  }
}

(******************************************************************************)
(* Implementations *)
(******************************************************************************)

inheritedImplementations =
  Join[{"ADMBase"},
       If [useGlobalDerivs, {"Coordinates"}, {}]];

(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

calculationsWeylScalars =
{
  WeylScalarsCalc,
  WeylScalarsBoundaryCalc
};

CreateKrancThornTT [groupsWeylScalars, ".", WeylScalars,
  Calculations -> calculationsWeylScalars,
  DeclaredGroups -> declaredGroupNamesWeylScalars,
  PartialDerivatives -> derivatives,
  EvolutionTimelevels -> evolutionTimelevels,
  UseLoopControl -> True
];
