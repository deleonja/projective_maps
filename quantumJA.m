(* ::Package:: *)

(*Author: Jos\[EAcute] Alfredo de Le\[OAcute]n*)
(*Date: August 07, 2020*)
BeginPackage["quantumJA`"]

Reshuffle::usage=
"Reshuffle[SqMatrix_List] reshuffles SqMatrix, a matrix of dimension \!\(\*SuperscriptBox[\(4\), \(n\)]\)."
Pauli::usage=
"Pauli[Indices_List] gives the tensor product of Pauli Matrices with indices in Indices_List."
PauliToComp::usage=
"PauliToComp[qbitsNum] constructs the change of basis matrix for a PCE operation from computational to Pauli basis."
PCE::usage=
"PCE[diagElements_List] calculates the matrix representation of a PCE operation in computational basis given its diagonal elements in Pauli basis."
PTest::usage=
"PTest[A_List] evaluates the positive-semidefiniteness of A with its eigenvalues."

Begin["`Private`"]
Reshuffle[SqMatrix_List]:=ArrayFlatten[ArrayFlatten/@Partition[Partition[
ArrayReshape[#,{Sqrt[Dimensions[SqMatrix][[1]]],Sqrt[
Dimensions[SqMatrix][[1]]]}]&/@SqMatrix,Sqrt[
Dimensions[SqMatrix][[1]]]],Sqrt[Dimensions[SqMatrix][[1]]]],1];

Pauli[0]=Pauli[{0}]={{1,0},{0,1}};
Pauli[1]=Pauli[{1}]={{0,1},{1,0}};
Pauli[2]=Pauli[{2}]={{0,-I},{I,0}};
Pauli[3]=Pauli[{3}]={{1,0},{0,-1}};
Pauli[Indices_List]:=KroneckerProduct@@(Pauli/@Indices);

PauliToComp[qbitsNum_Integer]:=
Flatten/@(Pauli[#]&/@Tuples[Range[0,3],qbitsNum])//Transpose

PCE[diagElements_List]:=
PauliToComp[Log[4,Length[diagElements]]].
DiagonalMatrix[diagElements].
Inverse[PauliToComp[Log[4,Length[diagElements]]]]

PTest[A_List]:=(A//Eigenvalues//Min)>=0

End[];
EndPackage[]



