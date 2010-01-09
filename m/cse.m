expr = {a -> Sin[Sqrt[x]] + Cos[Sqrt[x]],
        b -> Cos[Sqrt[x]]}

pexpr = expr //. {a_, b__} -> CSequence[a, {b}] //. {a_} -> 
    a //. (a_ -> b_) -> CAssign[a, b]

oexp = pexpr // Experimental`OptimizeExpression

{locals, code} = 
 ReleaseHold[(Hold @@ oexp) /. 
   Verbatim[Block][vars_, seq_] :> {vars, Hold[seq]}]

code1 = 
 code /. Hold[CompoundExpression[seq__]] :> Hold[{seq}]

code2 = 
 First[code1 //. 
   Hold[{a___Hold, b_, c___}] /; Head[Unevaluated[b]] =!= Hold :> 
    Hold[{a, Hold[b], c}]]

statements = 
 StringReplace[ToString[CForm[#]], 
    "Hold(" ~~ ShortestMatch[a___] ~~ ")" :> a] & /@ code2

mycsequence = StringJoin @@ Riffle[statements, ";\n"]

replacevar = 
 Rule @@@ Transpose[{ToString[CForm[#]] & /@ locals, 
    StringReplace[StringReplace[ToString[#], {__ ~~ "`" ~~ a_ :> a}], 
       "$" -> "_"] & /@ locals}]

mycsequence1 = StringReplace[mycsequence, replacevar]

final = 
 "{\ndouble " <> StringJoin @@ Riffle[Last /@ replacevar, ","] <> 
  ";\n\n" <> mycsequence1 <> ";\n}\n"
