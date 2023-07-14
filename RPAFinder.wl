(* ::Package:: *)

BeginPackage["RPAFinder`"];


(* Public functions *)

makeReactionSystem::usage="makeReactionSystem[reactions] creates a reaction system from a list of reactions.";

computeInfluenceIndex::usage="computeInfluenceIndex[subnetwork, system] computes the influence index of a subnetwork";

computeFluxInfluenceIndex::usage="computeFluxInfluenceIndex[subnetwork, system] computes the flux influence index of a subnetwork";

enumerateBufferingStructures::usage="enumerateBufferingStructures[system] returns the list of all the buffering structures in a reaction system. The subnetworks are specified by indices.";

enumerateBufferingStructures::usage="enumerateBufferingStructures[system] returns the list of all the buffering structures in a reaction system. The species are specified by names.";

enumerateLabeledBufferingStructuresByIndex::usage="enumerateLabeledBufferingStructuresByIndex[system] returns the list of all the labeled buffering structures in a reaction system. The subnetworks are specified by indices.";

enumerateLabeledBufferingStructures::usage="enumerateLabeledBufferingStructures[system] returns the list of all the labeled buffering structures in a reaction system. The species are specified by names.";

makeMatrixANum::usage="makeMatrixANum[system] creates the matrix A for the system whose nonzero \!\(\*SubscriptBox[\(r\), \(A, i\)]\) components are picked from random numbers"

lBSToBS::usage="lBStoBS[lbs] converts a labeled buffering structure to an ordinary buffering structure";

selectSubnetwork::usage="selectSubnetwork[p1,p2,lbs,system] returns species and reactions that are sensitive to all the elements of p1 and and not sensitive to all the elements of p2";

getEmergentConservedQuantities::usage="getEmergentConservedQuantities[gamma, system] returns emergent conserved quantities in subnetwork gamma";

getEmergentCycles::usage="getEmergentCycles[gamma, system] returns emergent conserved cycles in subnetwork gamma";

getLostConservedQuantities::usage="getLostConservedQuantities[gamma, system] returns lost conserved quantities";




(* Error messages *)

makeReactionSystem::noSpecies="No species found.";
makeReactionSystem::reactionsWrongFormat="The format of reactions is wrong."

computeInfluenceIndex::notReactionSystem="Reaction system object required.";
computeFluxInfluenceIndex::notReactionSystem="Reaction system object required.";
enumerateLabeledBufferingStructuresByIndex::notReactionSystem="Reaction system object required.";

computeInfluenceIndex::argsNotCorrect = "Arguments are not in the correct format.";
computeFluxInfluenceIndex::argsNotCorrect = "Arguments are not in the correct format.";

lBSToBS::arg="The argment is not a labeled buffering structure.";



Begin["`Private`"];

computeResponseMat::ninv = "The matrix A is not invertible";

{iChemicals,i\[Nu],iReactions}={1,2,3};

(*
makeReactionSystem[reactions]
Arguments: 
  reactions: list of reactions

Returns: a reaction system object, which consists of 
 - a list of species,
 - a stoichiometric matrix
 - a list of reactions
*)
makeReactionSystem[reactions_List]:=Module[
	{chemicals=Variables[reactions],\[Nu]},
	If[!AllTrue[reactions,Length@#==2&], Message[makeReactionSystem::reactionsWrongFormat];Return@$Failed ];
	If[Length@chemicals==0, Message[makeReactionSystem::noSpecies];Return@$Failed];
	\[Nu]= makeStoichiometricMatrix[reactions,chemicals];
	Return[{chemicals, \[Nu], reactions}]
];

makeStoichiometricMatrix[reactions_, vs_]:=Module[
	{from,to,\[Nu]},
	Table[
		{from, to} = Coefficient[#,vs]&/@ reactions[[i]];
		to - from,
		{i,1,Length@reactions}
	]//Transpose
];

makeReactantMatrix[reactions_,vs_]:= Module[{},
	Table[Coefficient[reactions[[i,1]],vs],{i,1,Length@reactions}]//Transpose
];

makeReactantMatrix[reactions_]:= makeReactantMatrix[reactions, Variables[reactions]];

(* create projection matrix to sub chemicals *)
projMat[xsub_, xall_]:=DiagonalMatrix[(If[MemberQ[xsub,#],1,0]&)/@xall]; 

projMatComp[xsub_, xall_]:=Module[{},IdentityMatrix[Length@xall] - projMat[xsub, xall]];

dimker[mat_]:=Length@NullSpace[mat];

dimcoker[mat_]:=Length@NullSpace[mat\[Transpose]];

(*
computeInfluenceIndex[xsargs, reactions, system]

Arguments: 
  xsargs: list of species of a subnetwork. It can be either by indices or names
  reactions: list of reactions of a subnetwork
  system: reaction system object

Returns: the influence index of the given subnetwork 
*)
computeInfluenceIndex[xsargs_, rs_, system_]:=Module[
	{xs, xsAllName,xall,\[Nu],rall,pr,pmbar},
	If[Length@system!=3, Message[computeInfluenceIndex::notReactionSystem];Return@$Failed];
	\[Nu] = system[[i\[Nu]]];
	xsAllName = system[[iChemicals]];	
	xs = If[ (Length@xsargs == 0) || NumberQ[xsargs[[1]]], xsargs, 
		Position[xsAllName,#][[1,1]]&/@ xsargs
	];
	
	xall = Range[1, Length@system[[iChemicals]]];
	rall = Range[1, Length@system[[iReactions]]];
	pr = projMat[rs,rall]; 
	pmbar = projMatComp[xs, xall];
	Length@rall + dimcoker@\[Nu] - dimker[\[Nu].pr]- dimcoker[pmbar.\[Nu]]
];

(*
computeInfluenceIndex[gamma, system]

Arguments: 
  gamma: subnetwork, specified by two lists {species, reactions}. Species can be specified either by indices or names. 
  system: reaction system object

Returns: the influence index of the given subnetwork 
*)
computeInfluenceIndex[gamma_, system_]:=Module[{},
	If[Length@gamma==2, Return@computeInfluenceIndex[gamma[[1]], gamma[[2]], system]];
	If[Length@gamma==4, Return@computeInfluenceIndex[lBSToBS[gamma], system]];
	Message[computeInfluenceIndex::argsNotCorrect];
	Return@$Failed
];

(*
computeFluxInfluenceIndex[xsargs, reactions, system]

Arguments: 
  xsargs: list of species of a subnetwork. It can be either by indices or names
  reactions: list of reactions of a subnetwork
  system: reaction system object

Returns: the flux influence index of the given subnetwork 
*)
computeFluxInfluenceIndex[xsargs_, rs_, system_]:=Module[
	{xsAllName,xs,xall,\[Nu],rall,pr,pmbar},
	If[Length@system!=3, Message[computeFluxInfluenceIndex::notReactionSystem]; Return@$Failed];
	\[Nu] = system[[i\[Nu]]];
	xsAllName = system[[iChemicals]];	
	xs = If[ (Length@xsargs == 0) || NumberQ[xsargs[[1]]], 
		xsargs, 
		Position[xsAllName,#][[1,1]]&/@ xsargs
	];

	xall = Range[1, Length@system[[iChemicals]]];
	rall = Range[1, Length@system[[iReactions]]];
	pr = projMat[rs,rall]; 
	pmbar = projMatComp[xs, xall];
	Length@rs + dimcoker@\[Nu] - dimcoker[pmbar.\[Nu]]
];

(*
computeFluxInfluenceIndex[gamma, system]

Arguments: 
  gamma: subnetwork, specified by two lists {species, reactions}. Species can be specified either by indices or names. 
  system: reaction system object

Returns: the flux influence index of the given subnetwork 
*)

computeFluxInfluenceIndex[gamma_, system_]:=Module[{},
	If[Length@gamma==2, Return@computeFluxInfluenceIndex[gamma[[1]], gamma[[2]], system]];
	If[Length@gamma==4, Return@computeFluxInfluenceIndex[lBSToBS[gamma], system]];
	Message[computeFluxInfluenceIndex::argsNotCorrect];
	Return@$Failed
];

(* xs and rs are sets of indices *)
isOutputComplete[xs_,rs_,system_]:=Module[
	{reMat,list,affectedRs},
	reMat=makeReactantMatrix[system[[iReactions]]];
	list=Total@reMat[[xs]];
	affectedRs=Select[Range[1, Length@list], list[[#]]>0&];
	ContainsAll[rs,affectedRs]
];

(* gamma = (xs,rs) is specified by indices *)
makeOutputComplete[gamma_,system_]:=Module[
	{xs,rs,reMat,list,affectedRs},
	{xs,rs}=gamma;
	If[isOutputComplete[xs,rs,system],Return@gamma];
	reMat=makeReactantMatrix[system[[iReactions]]];
	list=Total@reMat[[xs]];
	affectedRs=Select[Range[1, Length@list], list[[#]]>0&];
	Return[{xs, Sort@DeleteDuplicates@Union[rs,affectedRs]}]
];

findReactionsToMakeItOC[gamma_,system_]:=Module[
	{xs,rs,reMat,list,affectedRs},
	{xs,rs}=gamma;
	If[isOutputComplete[xs,rs,system],Return@{}];
	reMat=makeReactantMatrix[system[[iReactions]]];
	list=Total@reMat[[xs]];
	affectedRs=Select[Range[1, Length@list], list[[#]]>0&];
	Complement[affectedRs,rs]
];

(* 
makeMatrixANum[system]

Arguments: 
  system: a reaction system object 
 
Returns: A-matrix for this system where nonzero components of r_{A,i} are given by random numbers
*)
makeMatrixANum[system_, verbose_:False]:=Module[
	{chemicals, \[Nu],reactions, range, drdx, cv,dv,nx, nr, nc, nd,dim},
	{chemicals,\[Nu],reactions} = system[[{iChemicals,i\[Nu],iReactions}]];

	range={1,10}; (* The range of values from which random numbers are picked *)
	drdx=Map[If[#===0,0,RandomReal[range]]&, makeReactantMatrix[system[[iReactions]]]//Transpose,{2}];

	{cv,dv} = NullSpace/@{\[Nu], Transpose@\[Nu]};
	{nx,nr,nc,nd} = Length/@{chemicals,reactions,cv,dv};

	dim=nx+nc;
	Assert[nx+nc == nr+nd ];(*Fredholm's theorem*)
	If[verbose, 
		Print["chemicals: ",chemicals];
		Print["nx,nr,nc,nd: ",{nx,nr,nc,nd}];
		Print["drdx: ",MatrixForm@drdx];Print["cv: ",cv];Print["dv: ",dv];
	];

	Return@Table[
		Which[
			i<= nr && j <= nx, drdx[[i,j]],
			i<= nr && j >= nx, cv[[j-nx,i]],
			i>= nr && j <= nx, dv[[i-nr,j]],
			i>= nr && j >= nx, 0 
		]
	,
	{i,1,dim},{j,1,dim}]
];

makeMatrixANumWithData[system_]:=Module[
	{chemicals,\[Nu],reactions,amat,cv,dv,nx,nr,nc,nd},
	{chemicals,\[Nu],reactions} = system[[{iChemicals,i\[Nu],iReactions}]];
	{cv,dv} = NullSpace/@{\[Nu], Transpose@\[Nu]};
	{nx,nr,nc,nd} = Length/@{chemicals,reactions,cv,dv};
	amat = makeMatrixANum[system];
	{
	amat,{Join[Table[Subscript[r, i],{i,1,nr}],Table[Subscript[d, i],{i,1,nd}]],Join[chemicals,Table[Subscript[c, i],{i,1,nc}]]}
	}
];

computeResponseMat[system_]:=Module[
	{chemicals,\[Nu],reactions,cv,dv,nx,nr,nc,nd,amat,rs,cs,ainv,sResp,rResp},
	{chemicals,\[Nu],reactions} = system[[{iChemicals,i\[Nu],iReactions}]];
	{cv,dv} = NullSpace/@{\[Nu], Transpose@\[Nu]};
	{nx,nr,nc,nd} = Length/@{chemicals,reactions,cv,dv};
	
	amat = makeMatrixANum[system];	
	If[Det[amat]==0, Message[computeResponseMat::ninv]; Throw[amat, computeResponseMat::ninv]];
	
	ainv = Inverse[amat];
	
	sResp = ainv[[Range[nx],All]];
	rResp = cv\[Transpose].ainv[[Range[nx+1,nx+nc],Range[nr+nd]]];
	
	Return[{sResp, rResp}]
];

(* 
enumerateLabeledBufferingStructures[system]

Arguments: 
  system: a reaction system object 

Returns: the list of all the labeled buffering structures of the reaction network (species are represented by indices) 
*)
enumerateLabeledBufferingStructuresByIndex[system_]:=Module[
	{cs,cv,dv,nc,nd,\[Nu],reactions,nx,nr,sResp,rResp,\[Epsilon]=0.000001,sList,rList,bs,lbs,ubs,talliedLbs},
	If[Length@system!=3, Message[enumerateLabeledBufferingStructuresByIndex::notReactionSystem];Return@$Failed];	
	{cs,\[Nu],reactions} = system[[{iChemicals,i\[Nu],iReactions}]];
	{cv,dv} = NullSpace/@{\[Nu], Transpose@\[Nu]};
	{nx,nr,nc,nd} = Length/@{cs,reactions,cv,dv};

	Catch[
		{sResp, rResp}=computeResponseMat[system], computeResponseMat::ninv, Return@$Failed&
	];
	
	bs=Table[
		sList=Flatten@Position[sResp[[All,i]],_?(Abs[#]>\[Epsilon]&)];
		rList=Flatten@Position[rResp[[All,i]],_?(Abs[#]>\[Epsilon]&)];
		{i,sList,rList},
		{i,1,nr+nd}
	]; (* would-be buffering structures *)

	(* Would-be buffering structures obtained here include non-OC ones. They can be made OC by adding reactions *)
	lbs =  {#[[1]],#[[2]],#[[3]],findReactionsToMakeItOC[#[[2;;3]],system]}& /@bs;

	ubs= DeleteDuplicates[lbs[[All, 2;;]]];

	talliedLbs=
	Table[{lbs[[All,1]][[Flatten@Position[lbs[[All, 2;;]],_?(#==ubs[[i]]&)] ]]}~Join~ubs[[i]],{i,1,Length@ubs}];
		
	Assert[computeInfluenceIndex[#,system]==0]& /@ talliedLbs; (* check if the indices are actually zero *)

	talliedLbs
];

(* 
enumerateLabeledBufferingStructures[system]

Arguments: 
  system: a reaction system object 

Returns: the list of all the labeled buffering structures of the reaction network (species are represented by names) 
*)

enumerateLabeledBufferingStructures[system_]:=Module[
	{lbs,cs,\[Nu],reactions,nx,nr,rlabels},
	{cs,\[Nu],reactions} = system[[{iChemicals,i\[Nu],iReactions}]];
	{nx,nr} = Length/@{cs,reactions};

	rlabels=Table[i,{i,1,nr}];
	lbs = enumerateLabeledBufferingStructuresByIndex[system];

	Table[
	{lbs[[i,1]],cs[[#]]&/@lbs[[i,2]],rlabels[[#]]&/@lbs[[i,3]],rlabels[[#]]&/@lbs[[i,4]]},
	{i,1,Length@lbs}]
];

(* 
enumerateBufferingStructuresByIndex[system]

Arguments: 
  system: a reaction system object 

Returns: the list of all the buffering structures of the reaction network (species are represented by indices) 
*)

enumerateBufferingStructuresByIndex[system_]:=Module[{},
	lBSToBS[#]&/@ enumerateLabeledBufferingStructuresByIndex[system]
];

(* 
enumerateBufferingStructures[system]

Arguments: 
  system: a reaction system object 

Returns: the list of all the buffering structures of the reaction network (species are represented by names) 
*)

enumerateBufferingStructures[system_]:=Module[{},
	lBSToBS[#]&/@ enumerateLabeledBufferingStructures[system]
];

(* 
selectSubnetwork[p1,p2,lbs,system]

Arguments: 
  p1: list of reactions 
  p2: list of reactions 
  lbs: list of labeled buffering structures (species are specified by names)
  system: reaction system object
Returns: Species and reactions that do depend on 
         all the element of p1 but not depend on all the elements of p2
*)
selectSubnetwork[p1_,p2_,lbs_,system_]:=Module[
	{xall,rall,xsList,rsList,xs,rs,xnsList,rnsList,xns,rns},
	If[Length@lbs==0,Return@{}];
	xall=system[[1]];
	rall=Range[Length@system[[3]]];
	
	(* sensitive to p1 *)
	If[Length@p1!=0,
		{xsList,rsList} = (#[[{2,3}]]&/@ Select[lbs, ContainsAny[#[[1]],p1]&])//Transpose;
		xs = Intersection@@xsList;
		rs = Intersection@@rsList;
	,
		{xs,rs}={xall,rall};
	];
	
	(* insensitive to p2 *)	
	If[Length@p2!=0,
		{xnsList,rnsList} = 
		({Complement[xall, #[[2]]],Complement[rall, #[[3]]]} &/@ Select[lbs, ContainsAny[#[[1]],p2]&])//Transpose;	
		xns = Intersection@@xnsList;
		rns = Intersection@@rnsList;
	,
		{xns,rns}={xall,rall};
	];
		
	Return[{Intersection[xs,xns], Intersection[rs,rns]}]
];

(*
lBSToBS[lbs]

Argument:
  lbs: a labeled buffering structure

Returns: a buffering structure obtained by forgetting about the labels of the given labeled buffering structure. 
*)
lBSToBS[lbs_]:=Module[{},
	If[Length@lbs != 4, Message[lBSToBS::arg]; Return@$Failed];
	{lbs[[2]], Sort@Union[lbs[[3]],lbs[[4]]]}
];

isZeroVector[v_]:=Apply[And,#] & [ (#==0)&/@v];

(*
getIntersection[l1, l2] 

Argument:
	l1, l2: lists of vectors
Returns: a basis of span{l1} \cap span{l2}
*)
getIntersection[l1_, l2_] :=Module[{n1, kerL1L2, c1},
	If[Length@l1==0 || Length@l2==0, Return@{}];
	kerL1L2 = NullSpace[Transpose[Join[l1,l2]]];
	n1 = Length[l1];
	c1 = #[[1;;n1]] & /@ kerL1L2;
	If[Length@c1==0,Return@{}];
	Select[ #.l1&/@ RowReduce[c1], !isZeroVector[#]&]
];

getOrthogonalComplement[l1_]:=NullSpace[l1];


(* 
getEmergentCycles[gamma, system]

Arguments: 
  gamma: a subnetwork
  system: a reaction system object 

Returns: the list of emergent cycles in gamma
*)
getEmergentCycles[gamma_, system_]:=Module[
	{xs,rs,s11,s21,xsAll,smat,xhash,xsInd,xsCompInd,rsAll,rsComp,reactions,kerS21},
	{xs,rs}=gamma;
	{xsAll, smat,reactions} = system[[{iChemicals,i\[Nu],iReactions}]];
	rsAll=Range[1,Length@reactions];
	
	If[Length@rs==0, Return@{}];
	
	If[Length@xs==Length@xsAll && Length@rs==Length@rsAll, Return@{}]; (* gamma is the whole network *)
	
	xhash=Association@Table[xsAll[[i]]-> i,{i,1,Length@xsAll}];
	
	xsInd=xhash[#]&/@xs;
	xsCompInd=Complement[Range[1,Length@xsAll], xsInd];
	rsComp=Complement[rsAll, rs];
			
	s11=smat[[xsInd, rs]];
	s21=smat[[xsCompInd, rs]];
	
	If[Length@xs==0, 
		Return@getIntersection[
			IdentityMatrix[Length@rs],
			If[
				Length[kerS21=NullSpace@s21]==0, IdentityMatrix[Length@rs], 
				getOrthogonalComplement@kerS21
			]
		]
	];
	
	If[Length@xs==Length@xsAll, Return@{}]; (*ker s21 is the whole space and (ker s21)^perp is empty *)

	getEmergentCyclesFromMat[s11, s21]
];

getEmergentCyclesFromMat[s11_, s21_]:=Module[
	{kerS11=NullSpace@s11,kerS21=NullSpace@s21},
	If[Length@kerS11==0,Return@{}];
	If[Length@kerS21==0,Return@kerS11];
	getIntersection[kerS11, getOrthogonalComplement@kerS21]
];

(* 
getEmergentConservedQuantities[gamma, system]

Arguments: 
  gamma: a subnetwork
  system: a reaction system object 

Returns: the list of emergent conserved quantities in gamma
*)
getEmergentConservedQuantities[gamma_, system_]:=Module[
	{xs,rs,xsAll,smat,xhash,xsInd,xsCompInd,rsAll,rsComp,reactions},
	{xs,rs}=gamma;
	{xsAll, smat, reactions} = system[[{iChemicals,i\[Nu],iReactions}]];
	rsAll=Range[1,Length@reactions];
	
	If[Length@xs==0, Return@{}];
	
	xhash=Association@Table[xsAll[[i]]->i,{i,1,Length@xsAll}];
	
	xsInd=xhash[#]&/@xs;
	xsCompInd=Complement[Range[1,Length@xsAll], xsInd];
	rsComp=Complement[rsAll, rs];

	getEmergentConservedQuantitiesFromMat[
		smat[[ Join[xsInd,xsCompInd], Join[rs,rsComp] ]], 
		Length@xsInd, Length@rs
	]
];

getEmergentConservedQuantitiesFromMat[s_,nrow_,ncolumn_]:=Module[
{nrowAll,ncolumnAll,s11,cokerS11,wholeMat,cokerWholeMat,sol,D11},
	{nrowAll,ncolumnAll}=Dimensions@s;

	s11=s[[ 1;;nrow,1;;ncolumn ]];
	cokerS11=If[ ncolumn==0, IdentityMatrix[nrow], NullSpace@Transpose@s11 ];
	
	If[Length@cokerS11==0,Return@{}];

	(* 
	Make the following matrix : 
	( \[NoBreak]S11	S12	S11 )
    ( S21	S22	0   )
	*)
	wholeMat=ConstantArray[0,{nrowAll,ncolumnAll+ncolumn}];
	wholeMat[[1;;nrowAll, 1;;ncolumnAll]]=s;
	wholeMat[[1;;nrow, ncolumnAll+1;;ncolumnAll+ncolumn]]=s11;

	cokerWholeMat=NullSpace@Transpose@wholeMat;
	
	If[Length@cokerWholeMat==0, Return@cokerS11];
	
	sol=#[[1;;nrow]]&/@cokerWholeMat;

	D11=Select[RowReduce[sol],!isZeroVector[#]&];
	If[Length@D11==0, Return@cokerS11];

	Return@getIntersection[cokerS11, getOrthogonalComplement[D11]]
];

(* 
getLostConservedQuantities[gamma, system]

Arguments: 
  gamma: a subnetwork
  system: a reaction system object 

Returns: the list of lost conserved quantities
*)
getLostConservedQuantities[gamma_, system_]:=Module[
	{xs,rs,xsAll,smat,xhash,xsInd,xsCompInd,rsAll,rsComp,reactions,cokerS},
	{xs,rs}=gamma;
	{xsAll, smat, reactions} = system[[{iChemicals,i\[Nu],iReactions}]];
	rsAll=Range[1,Length@reactions];

	cokerS = NullSpace@Transpose@smat;	
	If[Length@xs==0, Return@{}];
	
	xhash=Association@Table[xsAll[[i]]->i,{i,1,Length@xsAll}];
	
	xsInd=xhash[#]&/@xs;
	xsCompInd=Complement[Range[1,Length@xsAll], xsInd];
	rsComp=Complement[rsAll, rs];

	getLostConservedQuantitiesFromMat[
		smat[[ Join[xsInd,xsCompInd], Join[rs,rsComp] ]], 
		Length@xsInd, Length@rs
	]
];

getLostConservedQuantitiesFromMat[s_,nrow_,ncolumn_]:=Module[
{nrowAll,ncolumnAll,cokerS,s11,cokerS11,wholeMat,spaceX},
	{nrowAll,ncolumnAll}=Dimensions@s;

	s11=s[[ 1;;nrow,1;;ncolumn ]];	
	cokerS11=If[ ncolumn==0, IdentityMatrix[nrow], NullSpace@Transpose@s11 ];
	cokerS=NullSpace@Transpose@s;

	(* 
	Make the following matrix : 
	( \[NoBreak]S11	S12	S11 )
    ( S21	S22	0   )
	*)
	wholeMat=ConstantArray[0,{nrowAll,ncolumnAll+ncolumn}];
	wholeMat[[1;;nrowAll, 1;;ncolumnAll]]=s;
	wholeMat[[1;;nrow, ncolumnAll+1;;ncolumnAll+ncolumn]]=s11;

	spaceX=NullSpace@Transpose@wholeMat;
	
	If[Length@spaceX==0, Return@cokerS];
	Return@getIntersection[cokerS, getOrthogonalComplement[spaceX]]
];




End[];


EndPackage[];
