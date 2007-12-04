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
  PDstandardNth[i_]     -> StandardCenteredDifferenceOperator[1,derivOrder/2,i],
  PDstandardNth[i_, i_] -> StandardCenteredDifferenceOperator[2,derivOrder/2,i],
  PDstandardNth[i_, j_] -> StandardCenteredDifferenceOperator[1,derivOrder/2,i]
                           StandardCenteredDifferenceOperator[1,derivOrder/2,j]
};

(* local derivatives *)
PDloc = PDstandardNth;

(* global derivatives *)
PDglob[var_,lx_] := Jinv[u1,lx] PDloc[var,l1];
PDglob[var_,lx_,ly_] :=
  dJinv[u1,lx,ly] PDloc[var,l1] + Jinv[u1,lx] Jinv[u2,ly] PDloc[var,l1,l2];

UseGlobalDerivs = True;
PD := If [UseGlobalDerivs, PDglob, PDloc];

KD = KroneckerDelta;

(******************************************************************************)
(* Tensors *)
(******************************************************************************)

(* Register the tensor quantities with the TensorTools package *)
Map [DefineTensor, {u, rho, v, w, J, Jinv, dJ, dJinv}];

AssertSymmetricIncreasing [dJ[uA,lb,lc], lb,lc];
AssertSymmetricIncreasing [dJinv[ua,lx,ly], lx,ly];

(******************************************************************************)
(* Groups *)
(******************************************************************************)

evolvedGroups =
  {SetGroupName [CreateGroupFromTensor [u  ], "WT_u"  ],
   SetGroupName [CreateGroupFromTensor [rho], "WT_rho"]};
evaluatedGroups = {};

declaredGroups = Join [evolvedGroups, evaluatedGroups];
declaredGroupNames = Map [First, declaredGroups];

evolvedGroupsFO =
  {SetGroupName [CreateGroupFromTensor [u    ], "WT_u"  ],
   SetGroupName [CreateGroupFromTensor [v[la]], "WT_v"  ],
   SetGroupName [CreateGroupFromTensor [rho  ], "WT_rho"]};
evaluatedGroupsFO =
  {SetGroupName [CreateGroupFromTensor [w[ua]], "WT_w"]};

declaredGroupsFO = Join [evolvedGroupsFO, evaluatedGroupsFO];
declaredGroupNamesFO = Map [First, declaredGroupsFO];



extraGroups =
  {{"MultiPatch::transformation",
    {dxda,dxdb,dxdc, dyda,dydb,dydc, dzda,dzdb,dzdc}},
   {"MultiPatch::transformation_inv",
    {dadx,dady,dadz, dbdx,dbdy,dbdz, dcdx,dcdy,dcdz}},
   {"MultiPatch::transformation_derivs",
    {ddxdada,ddxdadb,ddxdadc,ddxdbdb,ddxdbdc,ddxdcdc,
     ddydada,ddydadb,ddydadc,ddydbdb,ddydbdc,ddydcdc,
     ddzdada,ddzdadb,ddzdadc,ddzdbdb,ddzdbdc,ddzdcdc}},
   {"MultiPatch::transformation_inv_derivs",
    {ddadxdx,ddadxdy,ddadxdz,ddadydy,ddadzdz,ddadydz,
     ddbdxdx,ddbdxdy,ddbdxdz,ddbdydy,ddbdzdz,ddbdydz,
     ddcdxdx,ddcdxdy,ddcdxdz,ddcdydy,ddcdzdz,ddcdydz}}};



groups = Join [declaredGroups, extraGroups];
groupsFO = Join [declaredGroupsFO, extraGroups];

(******************************************************************************)
(* Initial data *)
(******************************************************************************)

initialCalc =
{
  Name -> "WT_Gaussian",
  Schedule -> {"AT initial"},
  (* Where -> Boundary, *)
  (* Where -> Interior, *)
  Equations -> 
  {
    u -> 0,
    rho -> 0
  }
}

initialCalcFO =
{
  Name -> "WTFO_Gaussian",
  Schedule -> {"AT initial"},
  (* Where -> Boundary, *)
  (* Where -> Interior, *)
  Equations -> 
  {
    u -> 0,
    v[la] -> 0,
    rho -> 0
  }
}

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

evolCalc =
{
  Name -> "WT_RHS",
  Schedule -> {"IN MoL_CalcRHS", "AT analysis"},
  Where -> Interior,
  Shorthands -> {Jinv[ua,lx], dJinv[ua,lx,ly]},
  Equations -> 
  {
    Jinv11 -> dadx,
    Jinv12 -> dady,
    Jinv13 -> dadz,
    Jinv21 -> dbdx,
    Jinv22 -> dbdy,
    Jinv23 -> dbdz,
    Jinv31 -> dcdx,
    Jinv32 -> dcdy,
    Jinv33 -> dcdz,
    dJinv111 -> ddadxdx,
    dJinv112 -> ddadxdy,
    dJinv113 -> ddadxdz,
    dJinv122 -> ddadydy,
    dJinv123 -> ddadydz,
    dJinv133 -> ddadzdz,
    dJinv211 -> ddadxdx,
    dJinv212 -> ddadxdy,
    dJinv213 -> ddadxdz,
    dJinv222 -> ddadydy,
    dJinv223 -> ddadydz,
    dJinv233 -> ddadzdz,
    dJinv311 -> ddadxdx,
    dJinv312 -> ddadxdy,
    dJinv313 -> ddadxdz,
    dJinv322 -> ddadydy,
    dJinv323 -> ddadydz,
    dJinv333 -> ddadzdz,
    dot[u] -> rho,
    dot[rho] -> KD[ua,ub] PD[u,la,lb]
  }
}

evolCalcFO =
{
  Name -> "WTFO_RHS",
  Schedule -> {"IN MoL_CalcRHS", "AT analysis"},
  Where -> Interior,
  Shorthands -> {Jinv[ua,lx]},
  Equations -> 
  {
    Jinv11 -> dadx,
    Jinv12 -> dady,
    Jinv13 -> dadz,
    Jinv21 -> dbdx,
    Jinv22 -> dbdy,
    Jinv23 -> dbdz,
    Jinv31 -> dcdx,
    Jinv32 -> dcdy,
    Jinv33 -> dcdz,
    dot[u] -> rho,
    dot[rho] -> KD[ua,ub] PD[v[la],lb],
    dot[v[la]] -> PD[rho,la]
  }
}

(******************************************************************************)
(* Constraint equations *)
(******************************************************************************)

constraintsCalcFO =
{
  Name -> "WTFO_constraints",
  Schedule -> {"AT analysis"},
  Where -> Interior,
  Shorthands -> {Jinv[ua,lx]},
  Equations -> 
  {
    Jinv11 -> dadx,
    Jinv12 -> dady,
    Jinv13 -> dadz,
    Jinv21 -> dbdx,
    Jinv22 -> dbdy,
    Jinv23 -> dbdz,
    Jinv31 -> dcdx,
    Jinv32 -> dcdy,
    Jinv33 -> dcdz,
    w[ua] -> Eps[ua,ub,uc] PD[v[lb],lc]
  }
}

(******************************************************************************)
(* Implementations *)
(******************************************************************************)

inheritedImplementations = {"MultiPatch"};

(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

calculations = 
{
  initialCalc,
  evolCalc
};

CreateKrancThornTT [groups, ".", "ML_WaveToy",
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  InheritedImplementations -> inheritedImplementations
];



calculationsFO = 
{
  initialCalcFO,
  evolCalcFO,
  constraintsCalcFO
};

CreateKrancThornTT [groupsFO, ".", "ML_FOWaveToy",
  Calculations -> calculationsFO,
  DeclaredGroups -> declaredGroupNamesFO,
  PartialDerivatives -> derivatives,
  InheritedImplementations -> inheritedImplementations
];
