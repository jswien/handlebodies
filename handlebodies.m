(* ::Package:: *)

(* ::Input::Initialization:: *)
SetDirectory[NotebookDirectory[]];
<<"FEMfine.m"


(* ::Title:: *)
(*Function Definitions*)


(* ::Section:: *)
(*Solution functions*)


(* ::Input::Initialization:: *)
PhiSol[]:=PhiSol[1]
PhiSol[n_]/;n<50:=If[Norm[phi[n]-phi[n-1],\[Infinity]]>10^-10,(*Print[{n,Norm[phi[n]-phi[n-1],\[Infinity]],DateString["Time"]}];*)PhiSol[n+1],phi[n]];
PhiSol[a_]/;a>=50:="Too Many Iterations";

PhiSolve[seed_]:=(
ClearAll[phi];
phi[0]=seed;
phi[n_]:=(phi[n]=phi[n-1]+LinearSolve[diffop[phi[n-1]],source[phi[n-1]]]);
PhiSol[])


(* ::Section:: *)
(*Geometric functions*)


(* ::Input::Initialization:: *)
Cen[\[Alpha]_,r_]:=Sin[\[Alpha]]/(Cos[r]-Cos[\[Alpha]]); 
Cen::usage = "Cen[\[Phi],\[Rho]] = Center of circle on equator of Riemann sphere at position \[Phi] and angular radius \[Rho].";

Rad[\[Alpha]_,r_]:=Abs[Sin[r]/(Cos[r]-Cos[\[Alpha]])];
Rad::usage = "Rad[\[Phi],\[Rho]] = Radius of circle on equator of Riemann sphere at position \[Phi] and angular radius \[Rho].";


(* ::Input::Initialization:: *)
L[\[Alpha]1_,\[Alpha]2_,r_]:=-{{-Csc[r] (Sin[(\[Alpha]1-\[Alpha]2)/2]+Cos[r] Sin[(\[Alpha]1+\[Alpha]2)/2]),(Cos[(\[Alpha]1-\[Alpha]2)/2]+Cos[r] Cos[(\[Alpha]1+\[Alpha]2)/2]) Csc[r]},{Cos[(\[Alpha]1+\[Alpha]2)/2] Cot[r]-Cos[(\[Alpha]1-\[Alpha]2)/2] Csc[r],-Csc[r] (Sin[(\[Alpha]1-\[Alpha]2)/2]-Cos[r] Sin[(\[Alpha]1+\[Alpha]2)/2])}};
L::usage="L[\[Alpha]1,\[Alpha]2,\[Rho]] = Loxodromic transformation between circles on the equator of Riemann sphere with angular positions \[Alpha]1 and \[Alpha]2 and radius \[Rho].";

Gen[w_,{\[Alpha]_,\[Beta]_,r_}]:=(L[\[Alpha],\[Beta],r][[1,1]] w+L[\[Alpha],\[Beta],r][[1,2]])/(L[\[Alpha],\[Beta],r][[2,1]] w+L[\[Alpha],\[Beta],r][[2,2]]);
Gen::usage="Gen[w,{\[Alpha]1,\[Alpha]2,\[Rho]}] = Mobius transformation on w between circles on the equator of Riemann sphere with angular positions \[Alpha]1 and \[Alpha]2 and radius \[Rho].";

HLen[g_]:=ArcCosh[(Tr[g]/2)^2+Abs[(Tr[g]/2)^2-1]]
HLen::usage="Hlen[L] = Length of closed bulk geodesic corresponding to transformation L.";


(* ::Title:: *)
(*Handlebody Class*)


(* ::Section:: *)
(*Init Handlebody*)


(* ::Input::Initialization:: *)
InitializeHandlebody[]:=(
circles = {{{0,0},1,"U"}};
symmetries = {"U"};
mcm=0.005;ag=4;
genus=0;
)
InitializeHandlebody::usage = "InitializeHandlebody[] = Set Handlebody parameters to defaults.";


(* ::Input::Initialization:: *)
AddCircle[circ_]:=If[Length@circ==3&&Length@circ[[1]]==2,
If[Not@MemberQ[{"Inv","R","InvR","RU","InvRU"},circ[[3]]],Return["Circle of unknown type"]];(circles=AppendTo[circles,circ];),Return["Circles must be of form {{cx,cy}, rad, type}"]];
AddCircle::usage ="AddCircle[circ] = Add circ = {cen, rad, type} to list of circles in domain.";

AddSymmetry[sym_]:=If[sym=="x" || sym=="y",(symmetries=AppendTo[symmetries,sym];),"Symmetries must be x or y as strings."]
AddSymmetry::usage ="AddSymmetry[sym] = Add sym to list of symmetries in domain, and sym can be either x or y.";


(* ::Input::Initialization:: *)
SolveHandlebody[name_:Handlebody]:=Block[{},
If[Length@circles==1,Return["Please Initialize Handlebody Circles"]];

syms=Length@symmetries;
genus=ComputeGenus[];
If[genus<2, Return["Genus must be bigger at least 2."]];
name["genus"]=genus;
name["circles"]=circles;
name["symmetries"]=symmetries;
name["mcm"]=mcm;
name["ag"]=ag;

If[StringQ@#,Return[#]]&@MakeMesh[];
name["mesh"]=mesh;
name["CError"]=CircleErrors[];
SolvePhi[];
name["sol"]=sol;
If[StringQ@sol,Return[sol]];
name["AError"]=AreaError[];
If[AreaError[]>3 10^-3,Print["High area error: "<>ToString@AreaError[]]];
name["BoundaryLengths"]=BoundaryLengths[];
name["Action"]=ComputeAction[];
]
SolveHandlebody::usage = "SolveHandlebody[name_:Handlebody] = Solve for handlebody solution. Circles and symmetries must be already defined. Creates parameters name[Attribute] with attributes being {Action, AError, BoundaryLengths, CError, Sol, genus, circles, symmetries, mcm, ag, mesh}";


(* ::Section:: *)
(*Make Mesh*)


(* ::Input::Initialization:: *)
MakeMesh[]:=(
(*Check that circles are not intersecting*)
If[1-RegionMeasure[RegionUnion@@Table[Disk@@circ[[;;2]],{circ,circles[[2;;]]}]]/Total[RegionMeasure/@Table[Disk@@circ[[;;2]],{circ,circles[[2;;]]}]]>10^-6,Return["Circles must be non-intersecting"]];

(*Construct geometric regions*)
regions = Table[Disk@@circ[[;;2]],{circ,circles[[2;;]]}]~Join~If[MemberQ[symmetries,"y"],{Rectangle[{-1,-1},{0,1}]},{}]~Join~If[MemberQ[symmetries,"x"],{Rectangle[{-1,-1},{1,0}]},{}];
mregion=RegionDifference[Disk[{0,0},1],RegionUnion@@regions];

(*Make mesh*)
mesh=ToElementMesh[mregion,MaxCellMeasure->mcm,AccuracyGoal->ag];

(*Straighten curved nodes*)
Do[
Do[If[Length[ElementNodes[mesh,e]\[Intersection]CircleNodes[mesh,circ[[1]],circ[[2]]]]==3,coors=ReplacePart[coors,#[[3]]->Mean[mesh["Coordinates"][[#[[1;;2]]]]]&@(ElementNodes[mesh,e]\[Intersection]CircleNodes[mesh,circ[[1]],circ[[2]]])]],{e,Length@MeshElements[mesh]}];
,{circ,circles}];

mesh=ToElementMesh[RegionDifference[Disk[{0,0},1],RegionUnion@@regions],MaxCellMeasure->mcm,AccuracyGoal->ag];
)
MakeMesh::usage = "MakeMesh[] = Make mesh for domain with specified circles and symmetries.";


(* ::Section:: *)
(*Estimate Mesh Errors*)


(* ::Input::Initialization:: *)
\[Theta]InvRU[cx_,r_]=Block[{cx,r},Abs@ArcTan[-((-1+cx^2+r^2)/(cx r)),1/(cx r) (\[Sqrt](-(-1+cx-r) (1+cx-r) (-1+cx+r) (1+cx+r)))]];
\[Theta]InvRUu[cx_,r_]=Block[{cx,r},Abs@ArcTan[(1+cx^2-r^2)/cx,1/cx (\[Sqrt](-(-1+cx-r) (1+cx-r) (-1+cx+r) (1+cx+r)))]];


(* ::Input::Initialization:: *)
CircleErrors[]:=(
(*Analytically compute correct Euclidean circle lengths*)
cints = Table[2\[Pi],Length@circles];

(*Correct circle lengths of type RU*)
Do[
If[circles[[ci,3]]=="RU",cints[[ci]]=2ArcTan[1/circles[[ci,2]]];cints[[1]]=cints[[1]]-4ArcTan[circles[[ci,2]]];];
,{ci,Length@circles}
];

(*Correct circle lengths of type InvRU*)
Do[
If[circles[[ci,3]]=="InvRU",
If[circles[[ci,1,1]]<0,
cints[[ci]]=2\[Theta]InvRU[circles[[ci,1,1]],circles[[ci,2]]];cints[[1]]=4\[Theta]InvRUu[circles[[ci,1,1]],circles[[ci,2]]]-2\[Pi];,
cints[[ci]]=2\[Pi]-2\[Theta]InvRU[circles[[ci,1,1]],circles[[ci,2]]];];
];
,{ci,Length@circles}
];

(*Reduce by symmetries*)
If[MemberQ[symmetries,"y"],
Do[
If[circles[[ci]][[1,1]]==0,cints[[ci]]=cints[[ci]]/2];
,{ci,Length@circles}
];
];
If[MemberQ[symmetries,"x"],
Do[
If[circles[[ci]][[1,2]]==0,cints[[ci]]=cints[[ci]]/2];
,{ci,Length@circles}
];
];

(*Adjust optimal tolerances*)
tols=Table[10^-10,Length@circles];
Do[While[Abs[1-Total@CircInt[mesh,circles[[ci,1]],circles[[ci,2]],tols[[ci]]]/cints[[ci]]]>5 10^-3&&tols[[ci]]<10^-5,
tols[[ci]]=10 tols[[ci]];
];,{ci,Length@circles}];

(*Report Errors*)
Abs@Table[1-Total@CircInt[mesh,circles[[ci,1]],circles[[ci,2]],tols[[ci]]]/cints[[ci]],{ci,Length@circles}]
)

CircleErrors::usage = "CircleErrors[] = Estimate numerical error for mesh.";


(* ::Section:: *)
(*Compute genus*)


(* ::Input::Initialization:: *)
ComputeGenus[]:=(
genus = 0;
(*Contribution for Inv or R type circles*)
Do[
If[circles[[ci,3]]=="Inv"||circles[[ci,3]]=="R",
genus+=Which[
circles[[ci,1]]=={0,0}, 
1,

circles[[ci,1,2]]==0,
If[MemberQ[symmetries,"y"],2,1],

circles[[ci,1,1]]==0,
If[MemberQ[symmetries,"x"],2,1],

True, 2^(syms-1)
];
];
,{ci,Length@circles}];

(*Contribution from RU type circles*)
Do[
If[circles[[ci,3]]=="RU",genus +=1;];
,{ci,Length@circles}];

(*Contribution from InvR type circles*)
Do[
If[circles[[ci,3]]=="InvR",genus +=1;];
,{ci,Length@circles}];

(*Contribution from InvRU type circles*)
Do[
If[circles[[ci,3]]=="InvRU",genus +=1/2;];
,{ci,Length@circles}];

Return[genus];
)
ComputeGenus::usage = "ComputeGenus[] = Compute the handlebody genus.";


(* ::Section:: *)
(*Solve for Phi*)


(* ::Input::Initialization:: *)
SolvePhi[]:=(
N\[Psi]N\[Psi][mesh];
source[phi_]:=N\[Psi]N\[Psi][mesh].phi + \[Psi]E[mesh,phi]-Sum[CircInt[mesh,circles[[ci,1]],circles[[ci,2]],tols[[ci]]],{ci,2,Length@circles}]+CircInt[mesh,{0,0},1,tols[[1]]];
diffop[phi_]:=-(N\[Psi]N\[Psi][mesh]+2\[Psi]\[Psi]E[mesh,phi]);
sol=PhiSolve[SparseArray[MeshLength@mesh->0]];
)
SolvePhi::usage = "SolvePhi[] = Solve for conformal factor phi."

AreaError[]:=Abs[1-IntPhiA[mesh,sol]/(4\[Pi] (genus-1)/2^syms)]
AreaError::usage = "AreaError[] = Estimate numerical error for solution by the Gauss-Bonnet theorem.";


(* ::Section:: *)
(*Lengths of Boundary Segments*)


(* ::Input::Initialization:: *)
BoundaryLengths[]:=(
{Table[CutLength[mesh,sol,circles[[ci,1]],circles[[ci,2]],tols[[ci]]],{ci,Length@circles}],

If[MemberQ[symmetries,"x"],
xintervals=Sort@Flatten[{#[[1,1]]-#[[2]],#[[1,1]]+#[[2]]}&/@Select[circles,#[[1,2]]==0&]];
xintervals=Select[xintervals,#<=1&];
If[MemberQ[symmetries,"y"],
xintervals=Select[xintervals,#>0&];
If[Length@Select[circles,#[[1]]=={0,0}&]==1,xintervals={0}~Join~xintervals~Join~{1}];
];
xintervals=Partition[xintervals,2];
Table[BoundLength[mesh,sol,xinterval],{xinterval,xintervals}],{}
],

If[MemberQ[symmetries,"y"],
yintervals=Sort@Flatten[{#[[1,2]]-#[[2]],#[[1,2]]+#[[2]]}&/@Select[circles,#[[1,1]]==0&]];
yintervals=Select[yintervals,#<=1&];
If[MemberQ[symmetries,"x"],
yintervals=Select[yintervals,#>0&];
If[Length@Select[circles,#[[1]]=={0,0}&]==1,yintervals={0}~Join~yintervals~Join~{1}];
];
yintervals=Partition[yintervals,2];
Table[BoundLengthY[mesh,sol,yinterval],{yinterval,yintervals}],{}
]
}
)
BoundaryLengths::usage = "BoundaryLengths[] = Compute the length of all boundary segments.";


(* ::Section:: *)
(*GR Action*)


(* ::Input::Initialization:: *)
logw[mesh_]:=Log@Norm[#]&/@mesh["Coordinates"]
IntPhiJ[mesh_,func_,jac_,c_,r_,tol_]:=CircInt[mesh,c,r,tol].SparseArray[{MeshLength@mesh->0}~Join~(#->func[[#]]jac@(mesh["Coordinates"][[#]]/._?(Abs[#]<10^-7&)->10^-7)&/@CircleNeighborhood[mesh,c,r])];

J\[Theta]\[Theta]0[cs_,R_]:=Block[{\[Theta]0},(R (R+c Cos[\[Theta]0]))/(c^2+R^2+2 c R Cos[\[Theta]0])/.c->Total[cs]/.\[Theta]0->ArcTan[(#[[1]]-cs[[1]])/R,(#[[2]]-cs[[2]])/R]&];
J\[Theta]I\[Theta]0[cs_,R_]:=Block[{\[Theta]0,d},(R (R-d Cos[\[Theta]0]))/(d^2+R^2-2 d R Cos[\[Theta]0])/.{\[Theta]0->ArcTan[(#[[1]]-cs[[1]])/R,(#[[2]]-cs[[2]])/R],d->R^2/c ((1-c^2+R^2)/(1+c^2-R^2))}/.c->Total[cs]&];
J\[Theta]Inv\[Theta]0[cs_,R_]:=Block[{\[Theta]0},(-c^2+R^2)/(c^2+R^2+2 c R Cos[\[Theta]0])/.c->Total[cs]/.{\[Theta]0->ArcTan[(#[[1]]-cs[[1]])/R,(#[[2]]-cs[[2]])/R]}&];
logarg[cs_,R_]:=Block[{d},(1-R^2/c^2)/(R^2-d^2)/.d->R^2/c ((1-c^2+R^2)/(1+c^2-R^2))/.c->Total[cs]];
logargF[cs_,R_]:=Block[{d},R^2-d^2/.d->R^2/c ((1-c^2+R^2)/(1+c^2-R^2))/.c->Total[cs]]; 

J3x[{c_,R_},{cp_,rp_}]=Block[{c,R,cp,rp,\[Theta],\[Theta]I,d},(cp^2-rp^2)/(cp^2+rp^2+2 cp rp Cos[\[Theta]]) (R (R-d Cos[\[Theta]I]))/(d^2+R^2-2 d R Cos[\[Theta]I])/.d-> R^2/c ((1-c^2+R^2)/(1+c^2-R^2))/.\[Theta]I->I Log[(cp+E^(I \[Theta]) rp)/(cp E^(I \[Theta])+rp)]/.{\[Theta]->ArcTan[(#[[1]]-cp)/rp,#[[2]]/rp]}&];
J3y[{c_,R_},{cp_,rp_}]=Block[{c,R,cp,rp,\[Theta],\[Theta]I,d},(cp^2-rp^2)/(cp^2+rp^2+2 cp rp Cos[\[Theta]]) (R (R-d Cos[\[Theta]I]))/(d^2+R^2-2 d R Cos[\[Theta]I])/.d-> R^2/c ((1-c^2+R^2)/(1+c^2-R^2))/.\[Theta]I->I Log[(cp+E^(I \[Theta]) rp)/(cp E^(I \[Theta])+rp)]/.{\[Theta]->ArcTan[#[[1]]/rp,(#[[2]]-cp)/rp]}&]; 
RadInv[c_,r_]=1/2 Abs[1/(c+r)-1/(c-r)];


(* ::Input::Initialization:: *)
ComputeAction[]:=(
(*Base action contribution from bulk and Unit circle*)
action=-4\[Pi](genus-1)(1-Log[4])-2^syms IntPhiE[mesh,sol]+2^syms IntPhiC[mesh,sol,{0,0},1,tols[[1]]];

(*Contribution from Inv type circles*)
Do[
If[circles[[ci,3]]=="Inv",
action +=- 2^syms IntPhiC[mesh,sol,circles[[ci,1]],circles[[ci,2]],tols[[ci]]];
action+=Which[
circles[[ci,1]]=={0,0},
-8\[Pi] Log[circles[[ci,2]]],

circles[[ci,1,2]]==0,
If[MemberQ[symmetries,"y"],-16\[Pi] Log[circles[[ci,2]]],-8\[Pi] Log[circles[[ci,2]]]],

circles[[ci,1,1]]==0,
If[MemberQ[symmetries,"x"],-16\[Pi] Log[circles[[ci,2]]],-8\[Pi] Log[circles[[ci,2]]]],

True,
-2^(syms-1)8\[Pi] Log[circles[[ci,2]]]
];
];
,{ci,Length@circles}];

(*Contribution from RU type circles*)
Do[
If[circles[[ci,3]]=="RU",
action +=- 2^syms IntPhiC[mesh,sol,circles[[ci,1]],circles[[ci,2]],tols[[ci]]] -8\[Pi] Log[circles[[ci,2]]]+8NIntegrate[x/Sin[x],{x,0,2 ArcTan[circles[[ci,2]]]}];
];
,{ci,Length@circles}];

(*Contribution from InvR type circles*)
Do[
If[circles[[ci,3]]=="InvR",
action +=(
2^syms IntPhiC[mesh,sol,circles[[ci,1]],circles[[ci,2]],tols[[ci]]]
-2^(syms+1) IntPhiJ[mesh,sol ,J\[Theta]\[Theta]0[circles[[ci,1]],circles[[ci,2]]],circles[[ci,1]],circles[[ci,2]],tols[[ci]]]
-2^(syms+1) IntPhiJ[mesh,sol,J\[Theta]I\[Theta]0[circles[[ci,1]],circles[[ci,2]]],circles[[ci,1]],circles[[ci,2]],tols[[ci]]]
+4\[Pi] Log[logarg[circles[[ci,1]],circles[[ci,2]]]]
);
];
,{ci,Length@circles}];

(*Contribution from InvRU type circles*)
If[MemberQ[circles[[;;,3]],"InvRU"],
circsIRUx=(#[[1,;;2]]~Join~#[[{2}]]&)/@(SortBy[Select[{circles,tols}\[Transpose],#[[1,3]]=="InvRU"&&#[[1,1,2]]==0&],#[[1,1,1]]&]);
circsIRUy=(#[[1,;;2]]~Join~#[[{2}]]&)/@(SortBy[Select[{circles,tols}\[Transpose],#[[1,3]]=="InvRU"&&#[[1,1,1]]==0&],#[[1,1,2]]&]);

If[Length@circsIRUx==2,
action+=(-2^(syms+1)IntPhiJ[mesh,sol ,J\[Theta]I\[Theta]0[circsIRUx[[1,1]],circsIRUx[[1,2]]],circsIRUx[[1,1]],circsIRUx[[1,2]],circsIRUx[[1,3]]]
-2^(syms+1) IntPhiJ[mesh,sol ,J\[Theta]\[Theta]0[circsIRUx[[1,1]],circsIRUx[[1,2]]],circsIRUx[[1,1]],circsIRUx[[1,2]],circsIRUx[[1,3]]]
+2^syms IntPhiC[mesh,sol,circsIRUx[[1,1]],circsIRUx[[1,2]],circsIRUx[[1,3]]]

-2^(syms+1) IntPhiJ[mesh,sol, J3x[{circsIRUx[[1,1,1]],circsIRUx[[1,2]]},{circsIRUx[[2,1,1]],circsIRUx[[2,2]]}],circsIRUx[[2,1]],circsIRUx[[2,2]],circsIRUx[[2,3]]]
-2^(syms+1) IntPhiJ[mesh,sol ,J\[Theta]\[Theta]0[circsIRUx[[2,1]],circsIRUx[[2,2]]],circsIRUx[[2,1]],circsIRUx[[2,2]],circsIRUx[[2,3]]]
+2^syms IntPhiC[mesh,sol ,circsIRUx[[2,1]],circsIRUx[[2,2]],circsIRUx[[2,3]]]

-2^(syms+1) IntPhiJ[mesh,logw[mesh],J\[Theta]\[Theta]0[circsIRUx[[1,1]],circsIRUx[[1,2]]],circsIRUx[[1,1]],circsIRUx[[1,2]],circsIRUx[[1,3]]]
-2^(syms+1) IntPhiJ[mesh,logw[mesh],J\[Theta]\[Theta]0[circsIRUx[[2,1]],circsIRUx[[2,2]]],circsIRUx[[2,1]],circsIRUx[[2,2]],circsIRUx[[2,3]]]
-2^(syms+2) IntPhiJ[mesh,logw[mesh] , J3x[{circsIRUx[[1,1,1]],circsIRUx[[1,2]]},{circsIRUx[[2,1,1]],circsIRUx[[2,2]]}],circsIRUx[[2,1]],circsIRUx[[2,2]],circsIRUx[[2,3]]]
-4\[Pi] Log@logargF[circsIRUx[[1,1]],circsIRUx[[1,2]]]);
];

If[Length@circsIRUy==2,
action+=(-2^(syms+1)IntPhiJ[mesh,sol ,J\[Theta]I\[Theta]0[circsIRUy[[1,1]],circsIRUy[[1,2]]],circsIRUy[[1,1]],circsIRUy[[1,2]],circsIRUy[[1,3]]]
-2^(syms+1) IntPhiJ[mesh,sol ,J\[Theta]\[Theta]0[circsIRUy[[1,1]],circsIRUy[[1,2]]],circsIRUy[[1,1]],circsIRUy[[1,2]],circsIRUy[[1,3]]]
+2^syms IntPhiC[mesh,sol,circsIRUy[[1,1]],circsIRUy[[1,2]],circsIRUy[[1,3]]]

-2^(syms+1) IntPhiJ[mesh,sol, J3y[{circsIRUy[[1,1,2]],circsIRUy[[1,2]]},{circsIRUy[[2,1,2]],circsIRUy[[2,2]]}],circsIRUy[[2,1]],circsIRUy[[2,2]],circsIRUy[[2,3]]]
-2^(syms+1) IntPhiJ[mesh,sol ,J\[Theta]\[Theta]0[circsIRUy[[2,1]],circsIRUy[[2,2]]],circsIRUy[[2,1]],circsIRUy[[2,2]],circsIRUy[[2,3]]]
+2^syms IntPhiC[mesh,sol ,circsIRUy[[2,1]],circsIRUy[[2,2]],circsIRUy[[2,3]]]

-2^(syms+1) IntPhiJ[mesh,logw[mesh],J\[Theta]\[Theta]0[circsIRUy[[1,1]],circsIRUy[[1,2]]],circsIRUy[[1,1]],circsIRUy[[1,2]],circsIRUy[[1,3]]]
-2^(syms+1) IntPhiJ[mesh,logw[mesh],J\[Theta]\[Theta]0[circsIRUy[[2,1]],circsIRUy[[2,2]]],circsIRUy[[2,1]],circsIRUy[[2,2]],circsIRUy[[2,3]]]
-2^(syms+2) IntPhiJ[mesh,logw[mesh] , J3y[{circsIRUy[[1,1,2]],circsIRUy[[1,2]]},{circsIRUy[[2,1,2]],circsIRUy[[2,2]]}],circsIRUy[[2,1]],circsIRUy[[2,2]],circsIRUy[[2,3]]]
-4\[Pi] Log@logargF[circsIRUy[[1,1]],circsIRUy[[1,2]]]);
];
];

(*Contribution from R type circles*)
Do[
If[circles[[ci,3]]=="R",
action +=(
-2^(syms+1) IntPhiJ[mesh,sol,J\[Theta]\[Theta]0[circles[[ci,1]],circles[[ci,2]]],circles[[ci,1]],circles[[ci,2]],tols[[ci]]]
+2^syms IntPhiJ[mesh,sol+2logw[mesh],J\[Theta]Inv\[Theta]0[circles[[ci,1]],circles[[ci,2]]],circles[[ci,1]],circles[[ci,2]],tols[[ci]]]
);

action+=If[
circles[[ci,1,2]]==0||circles[[ci,1,1]]==0,
4\[Pi] Log[(-(1/Norm[circles[[ci,1]]]^2)+1/circles[[ci,2]]^2)]+4\[Pi] Log[(1-circles[[ci,2]]^2/Norm[circles[[ci,1]]]^2)]-4\[Pi] Log[RadInv[Norm[circles[[ci,1]]],circles[[ci,2]]]^2],
2(4\[Pi] Log[(-(1/Norm[circles[[ci,1]]]^2)+1/circles[[ci,2]]^2)]+4\[Pi] Log[(1-circles[[ci,2]]^2/Norm[circles[[ci,1]]]^2)]-4\[Pi] Log[RadInv[Norm[circles[[ci,1]]],circles[[ci,2]]]^2])
];

];
,{ci,Length@circles}];

-(1/(24\[Pi]))Chop@action

)
ComputeAction::usage = "ComputeAction[] = Compute the handlebody action.";


(* ::Title:: *)
(*Minimization Functions*)


goldensearch[f_,initial_,accuracy_]:=Block[{list,mid,xpoint,iter,interval},
iter=0;
interval=initial;
While[Abs[f[interval[[2]]]]>accuracy && iter <100,
If[interval[[2]]-interval[[1]]>=interval[[3]]-interval[[2]],xpoint=interval[[1]]+(interval[[2]]-interval[[1]])/GoldenRatio;list=Insert[interval,xpoint,2];,xpoint=interval[[3]]-(interval[[3]]-interval[[2]])/GoldenRatio;list=Insert[interval,xpoint,3];];
mid=1+First@First@Position[f/@list[[2;;3]],Min[f/@list[[2;;3]]]];
++iter;
interval=list[[mid-1;;mid+1]]];
Return[interval[[2]]]
]
goldensearch::usage = "goldensearch[f, init, accuracy] = Do a golden ratio search to minimize function f starting with interval init.";

NM[f_,x_List,\[Delta]_List]:=x-Inverse@Table[1/\[Delta][[col]] ((f@@(x+ReplacePart[Table[0,Length[x]],col->\[Delta][[col]]]))[[row]]-(f@@x)[[row]]),{row,Length[f@@x]},{col,Length[x]}].(f@@x);

NM[f_,x_List,\[Delta]_List,\[Gamma]_]:=x-1/\[Gamma] Inverse@Table[1/\[Delta][[col]] ((f@@(x+ReplacePart[Table[0,Length[x]],col->\[Delta][[col]]]))[[row]]-(f@@x)[[row]]),{row,Length[f@@x]},{col,Length[x]}].(f@@x);
NM::usage = "NM[f, x0, \[Delta], \[Gamma]] = Solve the vector equation f[x] = 0 with Newton's method with initial guess x0. \[Gamma] is an optional parameter that decreases the step taken with each iteration."


(* ::Input::Initialization:: *)
GradSearch[TestF_,pstart_,delta_,gamma_,tol_,file_:False]:=Block[{POINT,TEST,POINTp,TESTP,cache,\[CapitalGamma]},
POINT=pstart;
\[CapitalGamma]=gamma;
cache={};
While[(TEST=TestF@POINT)>tol&&\[CapitalGamma]<10^6,
POINTp=POINT+1/\[CapitalGamma] Table[(TestF@POINT-TestF@(ReplacePart[POINT,pi->POINT[[pi]]+delta[[pi]]]))/(delta[[pi]]+If[delta[[pi]]==0,1,0]),{pi,Length@POINT}];
If[(TESTP=TestF@POINTp)<1.5TEST,
POINT=POINTp;AppendTo[cache,{POINT,TESTP}];If[file=!=False,Export[file,cache]];,
\[CapitalGamma]=\[CapitalGamma] 2.5;Continue[]];
];
Return[POINT]
]
GradSearch::usage="GradSearch[TestF,pstart,\[Delta],\[CapitalGamma],tol,file] = Perform a gradient search to minimize TestF with learning rate \[CapitalGamma]. TestF must be positive everywhere. Write results of each step to file (optional).";
