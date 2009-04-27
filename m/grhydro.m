Get["KrancThorn`"];

(*SetDebugLevel[InfoFull];*)

SetEnhancedTimes[False];
SetSourceLanguage["C"];

(****************************************************************************
 Derivatives
****************************************************************************)

derivatives =
{
 PDstandard2nd[i_] -> StandardCenteredDifferenceOperator[1,1,i],
 PDstandard2nd[i_, i_] -> StandardCenteredDifferenceOperator[2,1,i],
 PDstandard2nd[i_, j_] -> StandardCenteredDifferenceOperator[1,1,i] *
   StandardCenteredDifferenceOperator[1,1,j],

 PDstandard4th[i_] -> StandardCenteredDifferenceOperator[1,2,i],
 PDstandard4th[i_, i_] -> StandardCenteredDifferenceOperator[2,2,i],
 PDstandard4th[i_, j_] -> StandardCenteredDifferenceOperator[1,2,i] *
   StandardCenteredDifferenceOperator[1,2,j],

 PDplus[i_] -> DPlus[i],
 PDminus[i_] -> DMinus[i],
 PDplus[i_,j_] -> DPlus[i] DPlus[j],

 PDonesided2nd[1] -> dir[1] (-shift[1]^(2 dir[1]) + 4 shift[1]^dir[1] - 3 )/
   (2 spacing[1]),
 PDonesided2nd[2] -> dir[2] (-shift[2]^(2 dir[2]) + 4 shift[2]^dir[2] - 3 )/
   (2 spacing[2]),
 PDonesided2nd[3] -> dir[3] (-shift[3]^(2 dir[3]) + 4 shift[3]^dir[3] - 3 )/
   (2 spacing[3])
};

(*PD = PDstandard2nd;*)

(****************************************************************************
 Tensors 
****************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  rho, vel, eps, press,
  dens, mom, tau
}];

(* Register the TensorTools symmetries (this is very simplistic) *)
Map[AssertSymmetricDecreasing, 
{
}];

(* Determinants of the metrics in terms of their components
  (Mathematica symbolic expressions) *)
gDet = Det[MatrixOfComponents[g[la,lb]]];

(****************************************************************************
 Groups
****************************************************************************)

SetGroupTimelevels[g_,tl_] = Join[g, {Timelevels -> tl}];

evolvedGroups =
  {CreateGroupFromTensor [dens   ],             (* mass density *)
   CreateGroupFromTensor [mom[la]],             (* momentum density *)
   CreateGroupFromTensor [tau    ]};            (* energy density *)
evaluatedGroups =
  {CreateGroupFromTensor [rho    ],             (* mass density *)
   CreateGroupFromTensor [vel[ua]],             (* velocity *)
   CreateGroupFromTensor [eps    ],             (* specific internal energy *)
   CreateGroupFromTensor [press  ]};            (* pressure *)

declaredGroups = Join [evolvedGroups, evaluatedGroups];
declaredGroupNames = Map [First, declaredGroups];

groups = Join[declaredGroups];

(******************************************************************************)
(* Initial data *)
(******************************************************************************)

initialCalc =
{
  Name -> BSSN <> "_vacuum",
  Schedule -> {"IN ADMBase_InitialData"},
  ConditionalOnKeyword -> {"my_initial_data", "vacuum"},
  Equations -> 
  {
    rho     -> 0,
    vel[ua] -> 0,
    eps     -> 0
  }
};

(******************************************************************************)
(* Convert from primitive to conserved variables *)
(******************************************************************************)

prim2conCalc =
{
  Name -> BSSN <> "_con2prim",
  Schedule -> {"AT initial AFTER ADMBase_PostInitial"},
  Equations -> 
  {
    dens    -> rho,
    mom[la] -> rho vel[ua],
    tau     -> (1/2) rho vel[ua] vel[la] + rho eps
  }
};

(******************************************************************************)
(* Convert from conserved to primitive variables *)
(******************************************************************************)

con2primCalc =
{
  Name -> BSSN <> "_con2prim",
  Schedule -> {"IN " <> BSSN <>"_con2primGroup"},
  Equations -> 
  {
    rho     -> dens,
    vel[ua] -> mom[la] / dens,
    eps     -> tau / dens - (1/2) vel[ua] vel[la],
    
    press   -> Gamma rho eps +
               alpha PD[vel[ua],la]
  }
};

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

evolCalc =
{
  Name -> BSSN <> "_RHS",
  Schedule -> {"IN " <> BSSN <>"_evolCalcGroup"},
  Shorthands -> {rhov[ua], momv[la,ub], tauv[ua]},
  Equations -> 
  {
    (* dt rho + div rho v = 0 *)
    rhov[ua]  -> mom[la],
    dot[dens] -> - PD[rhov[ua],la],
    
    (* dt pi + div (pi v + P) = 0 *)
    momv[la,ub]  -> mom[la] vel[ub] + KD[la,ub] press,
    dot[mom[la]] -> - PD[momv[la,ub],lb],
    
    (* dt tau + div (tauv + P) = 0 *)
    tauv[ua] -> tau vel[ua] + press vel[ua],
    dot[tau] -> - PD[tauv[ua],la]
  }
};
