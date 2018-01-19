(* ::Package:: *)

Needs["NDSolve`FEM`"]
BeginPackage["FEMtools`"];


(* ::Section:: *)
(*Compute local functions \[Psi] and the local mass and stiffness matrix*)


(*Initialize \[Psi] on the reference element*)
\[Psi]R = {1-x (3-2 x-4 y)-(3-2 y) y,-(1-2 x) x,-(1-2 y) y,x (4-4 x-4 y),4 x y,-4 x y+(4-4 y) y};

(*Compute transformation to general element*)
Bk[{x1_,y1_},{x2_,y2_},{x3_,y3_}]=({
 {x2-x1, x3-x1},
 {y2-y1, y3-y1}
});
Fk[{x1_,y1_},{x2_,y2_},{x3_,y3_}]=(Bk[{x1,y1},{x2,y2},{x3,y3}].#+{x1,y1}&);
FIk[{x1_,y1_},{x2_,y2_},{x3_,y3_}]=(Inverse[Bk[{x1,y1},{x2,y2},{x3,y3}]].(#-{x1,y1})&);

(*Compute \[Psi] on general element*)
\[Psi]s[{x1_,y1_},{x2_,y2_},{x3_,y3_}]=Simplify[\[Psi]R/.{x->First@FIk[{x1,y1},{x2,y2},{x3,y3}]@{x,y},y->Last@FIk[{x1,y1},{x2,y2},{x3,y3}]@{x,y}}];
\[Psi]s[mesh_,e_]:= \[Psi]s@@MeshElement[mesh,e][[;;3]];

(*Store local mass and stiffness matrix on general element*)
IntMatrix = Compile[{x1,y1,x2,y2,x3,y3},{{1/60 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],-(1/360) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],-(1/360) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],0,-(1/90) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],0},{-(1/360) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],1/60 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],-(1/360) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],0,0,-(1/90) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3]},{-(1/360) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],-(1/360) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],1/60 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],-(1/90) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],0,0},{0,0,-(1/90) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],4/45 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],2/45 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],2/45 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3]},{-(1/90) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],0,0,2/45 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],4/45 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],2/45 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3]},{0,-(1/90) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],0,2/45 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],2/45 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3],4/45 Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3]}}];

IntGradMatrix =  Compile[{x1,y1,x2,y2,x3,y3},{{(((x2-x3)^2+(y2-y3)^2) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(2 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),(((x1-x3) (x2-x3)+(y1-y3) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(6 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),-((((x1-x2) (x2-x3)+(y1-y2) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(6 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),-((2 ((x1-x3) (x2-x3)+(y1-y3) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),0,(2 ((x1-x2) (x2-x3)+(y1-y2) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)},{(((x1-x3) (x2-x3)+(y1-y3) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(6 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),(((x1-x3)^2+(y1-y3)^2) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(2 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),((x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(6 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),-((2 ((x1-x3) (x2-x3)+(y1-y3) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),-((2 (x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),0},{-((((x1-x2) (x2-x3)+(y1-y2) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(6 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),((x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(6 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),(((x1-x2)^2+(y1-y2)^2) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(2 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),0,-((2 (x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),(2 ((x1-x2) (x2-x3)+(y1-y2) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)},{-((2 ((x1-x3) (x2-x3)+(y1-y3) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),-((2 ((x1-x3) (x2-x3)+(y1-y3) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),0,(4 (x1^2+x2^2-x2 x3+x3^2-x1 (x2+x3)+y1^2-y1 y2+y2^2-y1 y3-y2 y3+y3^2) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),(4 ((x1-x2) (x2-x3)+(y1-y2) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),-((4 (x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2))},{0,-((2 (x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),-((2 (x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),(4 ((x1-x2) (x2-x3)+(y1-y2) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),(4 (x1^2+x2^2-x2 x3+x3^2-x1 (x2+x3)+y1^2-y1 y2+y2^2-y1 y3-y2 y3+y3^2) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),-((4 ((x1-x3) (x2-x3)+(y1-y3) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2))},{(2 ((x1-x2) (x2-x3)+(y1-y2) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),0,(2 ((x1-x2) (x2-x3)+(y1-y2) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2),-((4 (x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),-((4 ((x1-x3) (x2-x3)+(y1-y3) (y2-y3)) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)),(4 (x1^2+x2^2-x2 x3+x3^2-x1 (x2+x3)+y1^2-y1 y2+y2^2-y1 y3-y2 y3+y3^2) Abs[-x2 y1+x3 y1+x1 y2-x3 y2-x1 y3+x2 y3])/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))^2)}}];


(* ::Section:: *)
(*Put together global function \[Psi]*)


\[Psi][mesh_,nn_,e_]:=(
\[Psi][mesh,nn,e]= \[Psi]s[mesh,e][[Position[ElementNodes[mesh,e],nn][[1,1]]]];
\[Psi][mesh,nn,e]
)

\[Psi][mesh_,nn_]:=Block[{psi},
psi={\[Psi][mesh,nn,#],{x,y}\[Element]MeshRegions[mesh][[#]]}&/@ NodeElements[mesh,nn];
\[Psi][mesh,nn]=psi;
\[Psi][mesh,nn]
]


(* ::Section:: *)
(*Define functions related to the mesh*)


MeshLength[mesh_]:=(MeshLength[mesh]=Length@mesh["Coordinates"])
MeshLength::usage = "MeshLength[mesh] = Number of nodes in mesh"

MeshArea[mesh_]:=Sum[Total[Integrate[\[Psi][mesh,nn,#],{x,y}\[Element]MeshRegions[mesh][[#]]]&
/@Neighborhood[mesh,nn]],{nn,MeshLength[mesh]}]

MeshRegions[mesh_]:=(MeshRegions[mesh]=Polygon[mesh["Coordinates"][[Riffle@@Partition[#,3]]]]&
/@mesh["MeshElements"][[1,1,;;]])
MeshRegions::usage = "MeshRegions[mesh] = List of mesh elements as polygonal regions."

MeshElement[mesh_,nn_]:=(MeshElement[mesh,nn]=mesh["Coordinates"][[mesh["MeshElements"][[1,1,nn]]]])
MeshElement::usage="MeshElement[mesh,e] = Coordinates of nodes of a particular element."

MeshElements[mesh_]:=mesh["Coordinates"][[#]]&/@mesh["MeshElements"][[1,1]]
MeshElements::usage="MeshElements[mesh] = List of all MeshElement[mesh,e]."

ElementNodes[mesh_,e_]:=(ElementNodes[mesh,e]=mesh["MeshElements"][[1,1,e]])
ElementNodes::usage="ElementNodes[mesh,e] = All node numbers in a particular element."

Unprotect[Line];
Line[a_]/;Length[a]<=1:= EmptyRegion[2];
Protect[Line];

NodeElements[mesh_,nn_]:=(NodeElements[mesh,nn]=
Flatten@Position[Normal[mesh["VertexElementConnectivity"][[nn]]],1])
NodeElements::usage="NodeElements[mesh,node] = All elements including a particular node."

NodeRegion[mesh_,nn_]/;1<=nn<=Length[mesh["MeshElements"][[1,1]]]&&nn\[Element]Integers:=
MeshRegions[mesh][[NodeElements[mesh,nn]]]
NodeRegion::usage="NodeRegion[mesh,node] = All polygonal regions including a particular node."


Neighbors[mesh_,nn_]:=(Neighbors[mesh,nn]=DeleteDuplicates@Flatten[mesh["MeshElements"][[1,1]][[NodeElements[mesh,nn]]]])
Neighbors::usage = "Neighbors[mesh,node] = All nodes in the neighborhood of a node."

Neighborhood[mesh_,nn_]:=(Neighborhood[mesh,nn]=NodeElements[mesh,nn])
Neighborhood::usage = "Neighborhood[mesh,nn] = All elements in the neighborhood of a node."


(* ::Section:: *)
(*Stiffness matrix*)


N\[Psi]N\[Psi]::usage = "N\[Psi]N\[Psi][mesh] = Stiffness matrix for a mesh"
N\[Psi]N\[Psi][mesh_]:=Block[{storage},
storage=SparseArray@Flatten@Table[Table[
{i,j} -> Total[((IntGradMatrix@@Flatten@MeshElement[mesh,#][[;;3]])[[First@Flatten@Position[ElementNodes[mesh,#],i],First@Flatten@Position[ElementNodes[mesh,#],j]]])&
/@(Neighborhood[mesh,i]\[Intersection]Neighborhood[mesh,j])],
{j,DeleteCases[Neighbors[mesh,i],_?(#>i&)]}],{i,Length@mesh["Coordinates"]}];
N\[Psi]N\[Psi][mesh]=storage+storage\[Transpose]-DiagonalMatrix@Diagonal[storage]
];


(* ::Section:: *)
(*Mass matrix*)


\[Psi]\[Psi]::usage = "\[Psi]\[Psi][mesh] = Mass matrix for a mesh"
\[Psi]\[Psi][mesh_]:=Block[{storage},
storage=SparseArray@Flatten@Table[Table[
{i,j} -> Total[((IntMatrix@@Flatten@MeshElement[mesh,#][[;;3]])[[First@Flatten@Position[ElementNodes[mesh,#],i],First@Flatten@Position[ElementNodes[mesh,#],j]]])&
/@(Neighborhood[mesh,i]\[Intersection]Neighborhood[mesh,j])],
{j,DeleteCases[Neighbors[mesh,i],_?(#>i&)]}],{i,Length@mesh["Coordinates"]}];
\[Psi]\[Psi][mesh]=storage+storage\[Transpose]-DiagonalMatrix@Diagonal[storage]
];

\[Psi]E::usage = "\[Psi]E[mesh,phi] = Mass matrix acting on e^(2 phi)"
\[Psi]E[mesh_,phi_]:=\[Psi]\[Psi][mesh].(SparseArray@Table[i->Exp[2 phi[[i]]],{i,MeshLength[mesh]}])

\[Psi]\[Psi]E::usage = "\[Psi]\[Psi]E[mesh,phi] = Mass matrix times e^(2 phi) times the identity matrix"
\[Psi]\[Psi]E[mesh_,phi_]:=\[Psi]\[Psi][mesh].(SparseArray@Table[{i,i}->Exp[2 phi[[i]]],{i,MeshLength[mesh]}])


(* ::Section:: *)
(*K matrix*)


KMatrix[1][{x1_,y1_},{x2_,y2_},{x3_,y3_}]=Norm[{x1,y1}-{x2,y2}]{{((x1-x2) (x2-x3)+(y1-y2) (y2-y3))/(2 \[Sqrt]((x1-x2)^2+(y1-y2)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3))/(6 \[Sqrt]((x1-x2)^2+(y1-y2)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(\[Sqrt]((x1-x2)^2+(y1-y2)^2))/(6 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(2 (x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)))/(3 \[Sqrt]((x1-x2)^2+(y1-y2)^2) (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),0,(2 \[Sqrt]((x1-x2)^2+(y1-y2)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3)))},{((x1-x2) (x2-x3)+(y1-y2) (y2-y3))/(6 \[Sqrt]((x1-x2)^2+(y1-y2)^2) (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3))/(2 \[Sqrt]((x1-x2)^2+(y1-y2)^2) (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(\[Sqrt]((x1-x2)^2+(y1-y2)^2))/(6 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(2 ((x1-x2) (x2-x3)+(y1-y2) (y2-y3)))/(3 \[Sqrt]((x1-x2)^2+(y1-y2)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(2 \[Sqrt]((x1-x2)^2+(y1-y2)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),0},{0,0,0,0,0,0},{(2 ((x1-x2) (x2-x3)+(y1-y2) (y2-y3)))/(3 \[Sqrt]((x1-x2)^2+(y1-y2)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(2 (x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)))/(3 \[Sqrt]((x1-x2)^2+(y1-y2)^2) (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(2 \[Sqrt]((x1-x2)^2+(y1-y2)^2))/(3 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(4 \[Sqrt]((x1-x2)^2+(y1-y2)^2))/(3 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(4 \[Sqrt]((x1-x2)^2+(y1-y2)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(4 \[Sqrt]((x1-x2)^2+(y1-y2)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3)))},{0,0,0,0,0,0},{0,0,0,0,0,0}};

KMatrix[2][{x1_,y1_},{x2_,y2_},{x3_,y3_}]=Norm[{x3,y3}-{x2,y2}]{{0,0,0,0,0,0},{(\[Sqrt]((x2-x3)^2+(y2-y3)^2))/(6 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),((x1-x3) (-x2+x3)+(y1-y3) (-y2+y3))/(2 \[Sqrt]((x2-x3)^2+(y2-y3)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),((x1-x2) (x2-x3)+(y1-y2) (y2-y3))/(6 \[Sqrt]((x2-x3)^2+(y2-y3)^2) (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(2 \[Sqrt]((x2-x3)^2+(y2-y3)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(2 ((x1-x2) (x2-x3)+(y1-y2) (y2-y3)))/(3 \[Sqrt]((x2-x3)^2+(y2-y3)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),0},{(\[Sqrt]((x2-x3)^2+(y2-y3)^2))/(6 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),((x1-x3) (x2-x3)+(y1-y3) (y2-y3))/(6 \[Sqrt]((x2-x3)^2+(y2-y3)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),((x1-x2) (x2-x3)+(y1-y2) (y2-y3))/(2 \[Sqrt]((x2-x3)^2+(y2-y3)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),0,(2 ((x1-x3) (x2-x3)+(y1-y3) (y2-y3)))/(3 \[Sqrt]((x2-x3)^2+(y2-y3)^2) (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(2 \[Sqrt]((x2-x3)^2+(y2-y3)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3)))},{0,0,0,0,0,0},{(2 \[Sqrt]((x2-x3)^2+(y2-y3)^2))/(3 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(2 ((x1-x3) (x2-x3)+(y1-y3) (y2-y3)))/(3 \[Sqrt]((x2-x3)^2+(y2-y3)^2) (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(2 ((x1-x2) (x2-x3)+(y1-y2) (y2-y3)))/(3 \[Sqrt]((x2-x3)^2+(y2-y3)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(4 \[Sqrt]((x2-x3)^2+(y2-y3)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(4 \[Sqrt]((x2-x3)^2+(y2-y3)^2))/(3 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(4 \[Sqrt]((x2-x3)^2+(y2-y3)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3)))},{0,0,0,0,0,0}};
KMatrix[3][{x1_,y1_},{x2_,y2_},{x3_,y3_}]=Norm[{x1,y1}-{x3,y3}]{{((x1-x3) (x2-x3)+(y1-y3) (y2-y3))/(2 \[Sqrt]((x1-x3)^2+(y1-y3)^2) (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(\[Sqrt]((x1-x3)^2+(y1-y3)^2))/(6 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3))/(6 \[Sqrt]((x1-x3)^2+(y1-y3)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(2 \[Sqrt]((x1-x3)^2+(y1-y3)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),0,-((2 (x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)))/(3 \[Sqrt]((x1-x3)^2+(y1-y3)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))))},{0,0,0,0,0,0},{((x1-x3) (x2-x3)+(y1-y3) (y2-y3))/(6 \[Sqrt]((x1-x3)^2+(y1-y3)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(\[Sqrt]((x1-x3)^2+(y1-y3)^2))/(6 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),(-(x1-x2) (x1-x3)-(y1-y2) (y1-y3))/(2 \[Sqrt]((x1-x3)^2+(y1-y3)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),0,(2 \[Sqrt]((x1-x3)^2+(y1-y3)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(2 ((x1-x3) (x2-x3)+(y1-y3) (y2-y3)))/(3 \[Sqrt]((x1-x3)^2+(y1-y3)^2) (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3)))},{0,0,0,0,0,0},{0,0,0,0,0,0},{-((2 ((x1-x3) (x2-x3)+(y1-y3) (y2-y3)))/(3 \[Sqrt]((x1-x3)^2+(y1-y3)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3)))),(2 \[Sqrt]((x1-x3)^2+(y1-y3)^2))/(3 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3))),-((2 (x1^2+x2 x3-x1 (x2+x3)+(y1-y2) (y1-y3)))/(3 \[Sqrt]((x1-x3)^2+(y1-y3)^2) (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3)))),(4 \[Sqrt]((x1-x3)^2+(y1-y3)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(4 \[Sqrt]((x1-x3)^2+(y1-y3)^2))/(3 (x3 (-y1+y2)+x2 (y1-y3)+x1 (-y2+y3))),(4 \[Sqrt]((x1-x3)^2+(y1-y3)^2))/(3 (x3 (y1-y2)+x1 (y2-y3)+x2 (-y1+y3)))}};

K\[Psi]::usage="K\[Psi][mesh,c,r] = Construct global K matrix associated with a boundary circle."
K\[Psi][mesh_,c_,r_]:=Block[{Kglobal,vers,Klocal},
Kglobal=SparseArray[{MeshLength@mesh,MeshLength@mesh}->0];
Do[
verts=ElementNodes[mesh,elem];

Klocal=Which[
MemberQ[CircleNodes[mesh,c,r],verts[[1]]]&&MemberQ[CircleNodes[mesh,c,r],verts[[2]]],
KMatrix[1]@@mesh["Coordinates"][[verts[[;;3]]]],

MemberQ[CircleNodes[mesh,c,r],verts[[2]]]&&MemberQ[CircleNodes[mesh,c,r],verts[[3]]],
KMatrix[2]@@mesh["Coordinates"][[verts[[;;3]]]],

MemberQ[CircleNodes[mesh,c,r],verts[[1]]]&&MemberQ[CircleNodes[mesh,c,r],verts[[3]]],
KMatrix[3]@@mesh["Coordinates"][[verts[[;;3]]]]
];

Do[Kglobal[[verts[[i]],verts[[j]]]]=Kglobal[[verts[[i]],verts[[j]]]]+Klocal[[i,j]],{i,6},{j,6}];
,{elem,CircleElements[mesh,c,r]}];

K\[Psi][mesh,c,r]=Kglobal

]


(* ::Section:: *)
(*Functions for understanding the boundary*)


CircleNodes::usage="CircleNodes[mesh,c,r] = All nodes on a circle of radius r with +/- infty corresponding to x and y axis."

CircleNodes[mesh_,c_,r_]:=(CircleNodes[mesh,c,r]=
Cases[DeleteDuplicates@Flatten[mesh["BoundaryElements"][[1,1]]],_?(Abs[Total[(mesh["Coordinates"][[#]]-c)^2]-r^2]<10^(-5)&),1])

CircleNodes[mesh_,c_,\[Infinity]]:=Cases[DeleteDuplicates@Flatten[mesh["BoundaryElements"][[1,1]]],_?(mesh["Coordinates"][[#,2]]<10^(-8)&&c[[1]]<=mesh["Coordinates"][[#,1]]<=c[[2]]&),1]
CircleNodes[mesh_,c_,-\[Infinity]]:=Cases[DeleteDuplicates@Flatten[mesh["BoundaryElements"][[1,1]]],_?(mesh["Coordinates"][[#,1]]<10^(-8)&&c[[1]]<=mesh["Coordinates"][[#,2]]<=c[[2]]&),1]


CircleNeighborhood::usage="CircleNeighborhood[mesh,c,r] = All nodes in the neighborhood of a circle."
CircleNeighborhood[mesh_,c_,r_]:=(CircleNeighborhood[mesh,c,r]=
Union@@(ElementNodes[mesh,#]&/@Select[Length[ElementNodes[mesh,#]\[Intersection]CircleNodes[mesh,c,r]]>1&]@(Union@@(Neighborhood[mesh,#]&/@CircleNodes[mesh,c,r]))))

CircleElements[mesh_,c_,r_]:=(CircleElements[mesh,c,r]=Select[Length[ElementNodes[mesh,#]\[Intersection]CircleNodes[mesh,c,r]]>1&]@(Union@@(Neighborhood[mesh,#]&/@CircleNodes[mesh,c,r])))


(*Functions for the boundary length of a general element*)

CircLength[{x1_,y1_},{x2_,y2_},{x3_,y3_}]=(
Which[
Abs[(x1-#1[[1]])^2+(y1-#1[[2]])^2-#2^2]<#3&&Abs[(x2-#1[[1]])^2+(y2-#1[[2]])^2-#2^2]<#3,
Sqrt[1+((y2-y1)/(x2-x1))^2]{1/6 Abs[-x1+x2],1/6 Abs[-x1+x2],0,2/3 Abs[x1-x2],0,0},
Abs[(x1-#1[[1]])^2+(y1-#1[[2]])^2-#2^2]<#3&&Abs[(x3-#1[[1]])^2+(y3-#1[[2]])^2-#2^2]<#3,
Sqrt[1+((y1-y3)/(x1-x3))^2]{1/6 Abs[-x1+x3],0,1/6 Abs[-x1+x3],0,0,2/3 Abs[x1-x3]},
Abs[(x3-#1[[1]])^2+(y3-#1[[2]])^2-#2^2]<#3&&Abs[(x2-#1[[1]])^2+(y2-#1[[2]])^2-#2^2]<#3,
Sqrt[1+((y2-y3)/(x2-x3))^2]{0,1/6 Abs[x2-x3],1/6 Abs[x2-x3],0,2/3 Abs[x2-x3],0},
True,
{0,0,0,0,0,0}
]
)&;

XLength[{x1_,y1_},{x2_,y2_},{x3_,y3_}]:=(
Which[
y1<10^-5 && y2<10^-5,
Sqrt[1+((y2-y1)/(x2-x1))^2]{1/6 Abs[-x1+x2],1/6 Abs[-x1+x2],0,2/3 Abs[x1-x2],0,0},
y1<10^-5 && y3<10^-5,
Sqrt[1+((y1-y3)/(x1-x3))^2]{1/6 Abs[-x1+x3],0,1/6 Abs[-x1+x3],0,0,2/3 Abs[x1-x3]},
y3<10^-5 && y2<10^-5,
Sqrt[1+((y2-y3)/(x2-x3))^2]{0,1/6 Abs[x2-x3],1/6 Abs[x2-x3],0,2/3 Abs[x2-x3],0},
True,
{0,0,0,0,0,0}
]
);

YLength[{x1_,y1_},{x2_,y2_},{x3_,y3_}]:=(
Which[
x1<10^-5 && x2<10^-5,
Sqrt[1+((x2-x1)/(y2-y1))^2]{1/6 Abs[-y1+y2],1/6 Abs[-y1+y2],0,2/3 Abs[y1-y2],0,0},
x1<10^-5 && x3<10^-5,
Sqrt[1+((x1-x3)/(y1-y3))^2]{1/6 Abs[-y1+y3],0,1/6 Abs[-y1+y3],0,0,2/3 Abs[y1-y3]},
x3<10^-5 && x2<10^-5,
Sqrt[1+((x2-x3)/(y2-y3))^2]{0,1/6 Abs[y2-y3],1/6 Abs[y2-y3],0,2/3 Abs[y2-y3],0},
True,
{0,0,0,0,0,0}
]
);


(* ::Section:: *)
(*Integrations*)


CircInt::usage = "CircInt[mesh,node,c,r,tol] = Compute integration of a function around a circle using d\[Theta]0 measure."
CircInt[mesh_,node_,c_,r_,tol_]:=(CircInt[mesh,node,c,r,tol]=
1/r Total[(CircLength@@MeshElement[mesh,#][[;;3]])[c,r,tol][[Position[ElementNodes[mesh,#],node][[1,1]]]]&/@
Select[Neighborhood[mesh,node],Length[ElementNodes[mesh,#]\[Intersection]CircleNodes[mesh,c,r]]>1&]])
CircInt[mesh_,c_,r_,tol_:10^-10]:=SparseArray@(Table[nn->CircInt[mesh,nn,c,r,tol],{nn,CircleNeighborhood[mesh,c,r]}]~Join~{MeshLength[mesh]->0})

CircInt[mesh_,node_,c_,r_,tol_]/;r==\[Infinity]:=(CircInt[mesh,node,c,r,tol]=
 Total[(XLength@@MeshElement[mesh,#][[;;3]])[[Position[ElementNodes[mesh,#],node][[1,1]]]]&/@
Select[Neighborhood[mesh,node],Length[ElementNodes[mesh,#]\[Intersection]CircleNodes[mesh,c,r]]>1&]])

CircInt[mesh_,node_,c_,r_,tol_]/;r==-\[Infinity]:=(CircInt[mesh,node,c,r,tol]=
 Total[(YLength@@MeshElement[mesh,#][[;;3]])[[Position[ElementNodes[mesh,#],node][[1,1]]]]&/@
Select[Neighborhood[mesh,node],Length[ElementNodes[mesh,#]\[Intersection]CircleNodes[mesh,c,r]]>1&]])


IntPhiA::usage = "IntPhiA[mesh,phi] = Compute the area of the domain using phi."
IntPhiA[mesh_,phi_]:=Total[\[Psi]\[Psi][mesh].(E^(2#)&/@phi)]

CutLength::usage = "CutLength[mesh,phi,c,r,tol] = Compute the length of a circle using metric."
CutLength[mesh_,phi_,c_,r_,tol_:10^-10]:= r CircInt[mesh,c,r,tol].Exp[phi]

BoundLength::usage = "BoundLength[mesh,phi,{x0,x1}] = Compute the length of an interval on the x-axis using metric."
BoundLength[mesh_,phi_,c_]:= CircInt[mesh,c,\[Infinity]].Exp[phi]

BoundLengthY::usage = "BoundLength[mesh,phi,{y0,y1}] = Compute the length of an interval on the y-axis using metric."
BoundLengthY[mesh_,phi_,c_]:= CircInt[mesh,c,-\[Infinity]].Exp[phi]

IntPhiC::usage = "IntPhiC[mesh,phi,c,r,tol] = Compute the integral of phi around C with measure d\[Theta]0"
IntPhiC[mesh_,phi_,c_,r_,tol_:10^-10]:= CircInt[mesh,c,r,tol].phi

IntPhiE::usage = "IntPhiE[mesh,phi] = Compute bulk integral term in action, i.e. integral of e^(2 phi)"
IntPhiE[mesh_,phi_]:=Total[\[Psi]\[Psi][mesh].(# E^(2#)&/@phi)]


EndPackage[];
