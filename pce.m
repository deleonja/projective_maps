(* ::Package:: *)

(* ::Section:: *)
(*Pauli component erasing*)


(* ::Input::Initialization:: *)
(*Author: Jos\[EAcute] Alfredo de Le\[OAcute]n*)
BeginPackage["pce`"]

Pauli::usage=
"Pauli[Indices_List] gives the tensor product of Pauli matrices (\!\(\*SubscriptBox[\(\[Sigma]\), \(0\)]\)=\[DoubleStruckOne]) with indices in Indices_List."
PCEsTauConfigurations::usage=
"PCEsTauConfigurations[\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)] gives all configurations of 1's and 0's of the elements \!\(\*FormBox[SubscriptBox[\(\[Tau]\), \(\*SubscriptBox[\(j\), \(1\)], \[Ellipsis], \*SubscriptBox[\(j\), \(n\)]\)],
TraditionalForm]\) of all PCE operations of \!\(\*FormBox[\(n\),
TraditionalForm]\) qubits. "
PCESuperoperator::usage=
"PCESuperoperator[\[Tau]] gives the superoperator, in computational basis, of the PCE operation characterized by the elements \!\(\*FormBox[SubscriptBox[\(\[Tau]\), \(\*SubscriptBox[\(j\), \(1\)], \[Ellipsis], \*SubscriptBox[\(j\), \(n\)]\)],
TraditionalForm]\) in the 1-D array \[Tau]."
Reshuffle::usage=
"Reshuffle[m] applies the reshuffle transformation to the matrix m with dimension \!\(\*SuperscriptBox[\(d\), \(2\)]\)\[Times]\!\(\*SuperscriptBox[\(d\), \(2\)]\)."
PositivityTest::usage=
"PositivityTest[\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)] evaluates if the matrix \!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\) is positive semidefinite."
PCEFigures::usage=
"PCEFigures[\[Tau]] generates a geometrical plot to visualize the configuration of the \!\(\*SubscriptBox[\(\[Tau]\), \(\*SubscriptBox[\(j\), \(1\)], \[Ellipsis], \*SubscriptBox[\(j\), \(n\)]\)]\) in the 1-D array \[Tau] of a PCE operation of 1, 2 or 3 qubits. "
TausOfPCEQuantumChannels::usage=
"TausOfPCEQuantumChannels[] returns all arrays of 1's and 0's of a PCE
	quantum channel."
NumberOfPCEOperations::usage=
"NumberOfPCEOperations[] returns the number of PCE operations."

Begin["`Private`"]
(*Pauli*)
Pauli[0]=Pauli[{0}]={{1,0},{0,1}};
Pauli[1]=Pauli[{1}]={{0,1},{1,0}};
Pauli[2]=Pauli[{2}]={{0,-I},{I,0}};
Pauli[3]=Pauli[{3}]={{1,0},{0,-1}};
Pauli[Indices_List]:=KroneckerProduct@@(Pauli/@Indices);

(*PCEsTauConfigurations*)
PCEsTauConfigurations[n_]:=MemoryConstrained[Normal[SparseArray[{1}~Join~#->ConstantArray[1,Length[{1}~Join~#]],{4^n}]]&/@Subsets[Range[2,4^n],4^n],2000000000,"Maximum usage of memory allowed exceeded"]

PCEsTauConfigurations[n_,k_]:=Normal[SparseArray[{1}~Join~#->ConstantArray[1,Length[{1}~Join~#]],{4^n}]]&/@Subsets[Range[2,4^n],{k-1}]

(*PCESuperoperator*)
PCESuperoperator[pauliDiagonal_List]:=Module[{indices,n,pauliToComp},
n=Log[4,Length[pauliDiagonal]];
indices=Tuples[Range[0,3],n];
pauliToComp=Transpose[Map[Flatten[Pauli[#]]&,indices]];
pauliToComp.DiagonalMatrix[pauliDiagonal].Inverse[pauliToComp]
]

(*Reshuffle*)
Reshuffle[m_]:=ArrayFlatten[ArrayFlatten/@Partition[Partition[ArrayReshape[#,{Sqrt[Dimensions[m][[1]]],Sqrt[Dimensions[m][[1]]]}]&/@m,Sqrt[Dimensions[m][[1]]]],Sqrt[Dimensions[m][[1]]]],1];

(*PositivityTest*)
PositivityTest[A_]:=Min[Eigenvalues[A]]>=0

(*PCEFigures*)
PCEFigures[\[Tau]_]:=Which[Length[\[Tau]]==4,ArrayPlot[{#}&/@\[Tau],ImageSize->50],Length[\[Tau]]==16,ArrayPlot[ArrayReshape[\[Tau],{4,4}],ImageSize->80],Length[\[Tau]]==64,Module[{cubeIndices},
cubeIndices=Position[ArrayReshape[\[Tau],{4,4,4}],1]-1;
Graphics3D[{If[Count[#,0]==3,{Black,Cube[#]},
If[Count[#,0]==2,{RGBColor["#CC0000"],Cube[#]},
If[Count[#,0]==1,{RGBColor["#004C99"],Cube[#]},
If[Count[#,0]==0,{RGBColor["#99FF33"],Cube[#]}]]]]&/@cubeIndices,{Thickness[0.012],Line[{{{-0.5,-0.5,-0.5},{-0.5,-0.5,3.5}},{{-0.5,-0.5,-0.5},{-0.5,3.5,-0.5}},{{-0.5,-0.5,-0.5},{3.5,-0.5,-0.5}},{{3.5,-0.5,-0.5},{3.5,-0.5,3.5}},{{-0.5,-0.5,3.5},{3.5,-0.5,3.5}},{{-0.5,3.5,-0.5},{3.5,3.5,-0.5}},{{3.5,3.5,-0.5},{3.5,3.5,3.5}},{{3.5,3.5,3.5},{-0.5,3.5,3.5}},{{-0.5,3.5,3.5},{-0.5,3.5,-0.5}},{{-0.5,3.5,3.5},{-0.5,-0.5,3.5}},{{3.5,3.5,3.5},{3.5,-0.5,3.5}},{{3.5,3.5,-0.5},{3.5,-0.5,-0.5}}}]}},Axes->False,AxesLabel->{"x","y","z"},LabelStyle->Directive[Bold,Medium,Black],PlotRange->{{-0.5,3.5},{-0.5,3.5},{-0.5,3.5}},AxesOrigin->{0.5,0.5,0.5},AxesStyle->Thickness[0.005],ImageSize->Small,ImagePadding->45]
]]

(*TausOfPCEQuantumChannels*)
TausOfPCEQuantumChannels[n_]:=MemoryConstrained[
Flatten[Table[If[PositivityTest[Reshuffle[PCESuperoperator[#]]]==True,#,Nothing]&/@PCEsTauConfigurations[n,i],{i,1,4^n}],1],2000000000,"Maximum usage of memory allowed exceeded. Try using PCEQuantumChannelsTaus[\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\)] or PCEQuantumChannelsTaus[\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)kmin\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"kmax\",\nFontSlant->\"Italic\"]\)]."]

TausOfPCEQuantumChannels[n_,k_]:=MemoryConstrained[
If[PositivityTest[Reshuffle[PCESuperoperator[#]]]==True,#,Nothing]&/@PCEsTauConfigurations[n,k],2000000000,"Maximum usage of memory allowed exceeded."]

TausOfPCEQuantumChannels[n_,kmin_,kmax_]:=MemoryConstrained[
Flatten[Table[If[PositivityTest[Reshuffle[PCESuperoperator[#]]]==True,#,Nothing]&/@PCEsTauConfigurations[n,i],{i,kmin,kmax}],1],2000000000,"Maximum usage of memory allowed exceeded. Try using PCEQuantumChannelsTaus[\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\)]."]

(*NumberOfPCEOperations*)
NumberOfPCEOperations[n_]:=2^(4^n-1)

NumberOfPCEOperations[n_,k_]:=Binomial[4^n-1,k-1]

End[];
EndPackage[]


(* ::InheritFromParent:: *)
(*"PauliToComp[qbitsNum] calculates the change of basis matrix for a qubit-system map from computational to tensor product of Pauli matrices basis."*)
