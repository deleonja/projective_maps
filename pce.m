(* ::Package:: *)

(* ::Section:: *)
(*Component erasing maps package*)


(* ::Input::Initialization:: *)
(*Author: Jos\[EAcute] Alfredo de Le\[OAcute]n*)
BeginPackage["pce`"]

Pauli::usage=
"Pauli[Indices_List] gives the tensor product of Pauli matrices (\!\(\*SubscriptBox[\(\[Sigma]\), \(0\)]\)=\[DoubleStruckOne]) with indices in Indices_List."
Reshuffle::usage=
"Reshuffle[m] applies the reshuffle transformation to the matrix m with dimension \!\(\*SuperscriptBox[\(d\), \(2\)]\)\[Times]\!\(\*SuperscriptBox[\(d\), \(2\)]\)."
PCESuperoperator::usage=
"PCESuperoperator[\[Tau]] gives the superoperator, in computational basis, of the PCE operation characterized by the elements \!\(\*FormBox[SubscriptBox[\(\[Tau]\), \(\*SubscriptBox[\(j\), \(1\)], \[Ellipsis], \*SubscriptBox[\(j\), \(n\)]\)],
TraditionalForm]\) in the 1-D array \[Tau]."
PositivityTest::usage=
"PositivityTest[\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)] evaluates if the matrix \!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\) is positive semidefinite."
Dirac::usage=
"Dirac[vector] returns vector in Dirac notation in computational basis."
Erase::usage=
"Erase[EigInfo,PCEfrom,invariantComponents] erases in all valid forms."
BlochSphTransformation::usage="Returns Bloch Ball transformation of a 1-qubit quantum chanenl given a
list with center point and the factor of x,y and z. 
Example: BlochSphTransformation[{{0,0,0.3},{1,1/2,1/2}}] returns 
a taco-like figure with center at (0,0,0.3).
BlochSphTransformation[coord_List]
"
PCEgenerators::usage=
"PCEgenerators[n] returns the diagonal superoperators of the generators
of PCE quantum channels of n qubits."
TensorPower::usage=
"TensorPower[A,n] gives the n-times Kronecker product of A."
PCEFigures::usage=
"PCEFigures[\[Tau]] generates a geometrical plot to visualize the configuration of the \!\(\*SubscriptBox[\(\[Tau]\), \(\*SubscriptBox[\(j\), \(1\)], \[Ellipsis], \*SubscriptBox[\(j\), \(n\)]\)]\) in the 1-D array \[Tau] of a PCE operation of 1, 2 or 3 qubits. "

Begin["`Private`"]
(*Pauli*)
Pauli[0]=Pauli[{0}]={{1,0},{0,1}};
Pauli[1]=Pauli[{1}]={{0,1},{1,0}};
Pauli[2]=Pauli[{2}]={{0,-I},{I,0}};
Pauli[3]=Pauli[{3}]={{1,0},{0,-1}};
Pauli[Indices_List]:=KroneckerProduct@@(Pauli/@Indices);

(*Reshuffle*)
Reshuffle[m_]:=ArrayFlatten[ArrayFlatten/@Partition[Partition[ArrayReshape[#,{Sqrt[Dimensions[m][[1]]],Sqrt[Dimensions[m][[1]]]}]&/@m,Sqrt[Dimensions[m][[1]]]],Sqrt[Dimensions[m][[1]]]],1];

(*PCESuperoperator*)
PCESuperoperator[pauliDiagonal_List]:=Module[{indices,n,pauliToComp},
n=Log[4,Length[pauliDiagonal]];
indices=Tuples[Range[0,3],n];
pauliToComp=Transpose[Map[Flatten[Pauli[#]]&,indices]];
pauliToComp.DiagonalMatrix[pauliDiagonal].Inverse[pauliToComp]
]

(*PCEFigures*)
PCEFigures[\[Tau]_]:=Which[Length[\[Tau]]==4,ArrayPlot[{#}&/@\[Tau],ImageSize->50],Length[\[Tau]]==16,ArrayPlot[ArrayReshape[\[Tau],{4,4}],ImageSize->80],Length[\[Tau]]==64,Module[{cubeIndices},
cubeIndices=Position[ArrayReshape[\[Tau],{4,4,4}],1]-1;
Graphics3D[{If[Count[#,0]==3,{Black,Cube[#]},
If[Count[#,0]==2,{RGBColor["#CC0000"],Cube[#]},
If[Count[#,0]==1,{RGBColor["#004C99"],Cube[#]},
If[Count[#,0]==0,{RGBColor["#99FF33"],Cube[#]}]]]]&/@cubeIndices,{Thickness[0.012],Line[{{{-0.5,-0.5,-0.5},{-0.5,-0.5,3.5}},{{-0.5,-0.5,-0.5},{-0.5,3.5,-0.5}},{{-0.5,-0.5,-0.5},{3.5,-0.5,-0.5}},{{3.5,-0.5,-0.5},{3.5,-0.5,3.5}},{{-0.5,-0.5,3.5},{3.5,-0.5,3.5}},{{-0.5,3.5,-0.5},{3.5,3.5,-0.5}},{{3.5,3.5,-0.5},{3.5,3.5,3.5}},{{3.5,3.5,3.5},{-0.5,3.5,3.5}},{{-0.5,3.5,3.5},{-0.5,3.5,-0.5}},{{-0.5,3.5,3.5},{-0.5,-0.5,3.5}},{{3.5,3.5,3.5},{3.5,-0.5,3.5}},{{3.5,3.5,-0.5},{3.5,-0.5,-0.5}}}]}},Axes->False,AxesLabel->{"x","y","z"},LabelStyle->Directive[Bold,Medium,Black],PlotRange->{{-0.5,3.5},{-0.5,3.5},{-0.5,3.5}},AxesOrigin->{0.5,0.5,0.5},AxesStyle->Thickness[0.005],ImageSize->Small,ImagePadding->45]
]]

(*PositivityTest*)
PositivityTest[A_]:=Min[Eigenvalues[A]]>=0

(*Dirac*)
Dirac[vector_List]:=(vector[[#]]Ket[IntegerString[(#-1),2,Log[2,Length[vector]]]])&/@Delete[Range[Length[vector]],Position[vector,0]]//Total

(*Erase*)
Erase[EigInfo_,PCEfrom_,invariantComponents_]:=Module[{dimPCE},
dimPCE=Length[Dimensions[PCEfrom][[2;;]]];
If[Count[#//Flatten,1]==invariantComponents,#,Nothing]&/@
DeleteDuplicates[
Flatten[
Table[
DeleteDuplicates[ReplacePart[PCEfrom[[i]]+#-ConstantArray[1,ConstantArray[4,dimPCE]],#->0&/@Position[PCEfrom[[i]]+#-ConstantArray[1,ConstantArray[4,dimPCE]],-1]]&/@EigInfo]
,{i,Length[PCEfrom]}]
,1]
]
]

(*BlochSphTransformation*)
BlochSphTransformation[coord_]:=Module[{x0,y0,z0,a,b,c},
{x0,y0,z0}=coord[[1]];
{a,b,c}=If[#==0,0.02,#]&/@coord[[2]];
Style[Show[ContourPlot3D[x^2+y^2+z^2==1,{x,-1,1},{y,-1,1},{z,-1,1},
ContourStyle->{Yellow,Opacity[0.25]},Mesh->None],
ContourPlot3D[(x-x0)^2/a^2+(y-y0)^2/b^2+(z-z0)^2/c^2==1,{x,-1,1},{y,-1,1},{z,-1,1},
ContourStyle->{Dashed,Pink,Opacity[0.65]},Mesh->None],
Graphics3D[{
        Black, Arrow[Tube[{{0,0,0},1.3*Normalize[{1,0,0}]}],0.05],
        Black,Arrow[Tube[{{0,0,0},1.3*Normalize[{0,1,0}]}],0.05],
        Black,Arrow[Tube[{{0,0,0},1.3*Normalize[{0,0,1}]}],0.05],
Text["x",{1.3,0,0}],Text["y",{0,1.3,0}],Text["z",{0,0,1.3}] },
Boxed->False],Boxed->False,Axes->False,PlotRange->1.3],
RenderingOptions->{"3DRenderingMethod"->"HardwareDepthPeeling"}]
]

(*TensorPower*)
TensorPower[A_, n_] := Nest[KroneckerProduct[A, #] &, A, n - 1]

(*PCEgenerators*)
PCEgenerators[n_]:=Module[{a},
a={{1,1,1,1},{1,1,-1,-1},{1,-1,1,-1},{1,-1,-1,1}};
DiagonalMatrix[#]&/@ReplacePart[TensorPower[a,n],Position[TensorPower[a,n],-1]->0][[2;;]]
]
End[];
EndPackage[]


(* ::InheritFromParent:: *)
(*"PauliToComp[qbitsNum] calculates the change of basis matrix for a qubit-system map from computational to tensor product of Pauli matrices basis."*)
