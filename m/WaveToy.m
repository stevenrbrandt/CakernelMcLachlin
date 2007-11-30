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

PD = PDstandardNth;

KD = KroneckerDelta;

(******************************************************************************)
(* Tensors *)
(******************************************************************************)

(* Register the tensor quantities with the TensorTools package *)
Map [DefineTensor, {u, rho}];

(******************************************************************************)
(* Groups *)
(******************************************************************************)

evolvedGroups =
  {SetGroupName [CreateGroupFromTensor [u  ], "WT_u"  ],
   SetGroupName [CreateGroupFromTensor [rho], "WT_rho"]};
evaluatedGroups = {};

declaredGroups = Join [evolvedGroups, evaluatedGroups];
declaredGroupNames = Map [First, declaredGroups];



groups = Join [declaredGroups];

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

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

evolCalc =
{
  Name -> "WT_RHS",
  Schedule -> {"IN MoL_CalcRHS", "AT analysis"},
  Where -> Interior,
  Equations -> 
  {
    dot[u] -> rho,
    dot[rho] -> KD[ua,ub] PD[u,la,lb]
  }
}

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
  PartialDerivatives -> derivatives
];
