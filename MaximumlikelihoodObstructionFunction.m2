
newPackage(
    "MaximumlikelihoodObstructionFunction",
    Version => "1.0", 
    Date => "01 Oct 2018",
    Authors => {
   {Name => "Jose Israel Rodriguez",
       Email => "Joisro@UChicago.edu",
       HomePage => "http://home.uchicago.edu/joisro/"}
    },
    Headline => "Produces equations and computes Euler obstructions. ",
    DebuggingMode => true,
    AuxiliaryFiles => false,
    PackageImports => {"SimpleDoc","Bertini"},
  Configuration => { "RandomCoefficients"=>(1,30102),
      "Continuation"=>Bertini },
  CacheExampleOutput => false
)
randomValue=()->random ((options MaximumlikelihoodObstructionFunction).Configuration#"RandomCoefficients")

--Copyright 2017
--  You may redistribute this file under the terms of the GNU General
--  Public License as published by the Free Software Foundation,
--  either version 2 of the License, or any later version.

export { 
    "RemovalMLDegree",
    "MLDegreeWitnessSet",
    "MLDegreeWitnessCollection",
--Standard symbolic formulation
  "MLDegreeVariety",
  "DefiningEquations",
       "Hyperplanes",
       "Data0",
       "Data1",
    "newRemovalMLDegree",
    "newMLDegreeVariety",
--Numerical continuation method
    "newMLDegreeWitnessSet",
    "newMLDegreeWitnessCollection",
----
    "randomVector",  --
    "removalMLDegree",
    "mlObstructionFunction",
    "solveRemovalMLDegree",
    "reclassifyWitnessPoints",
    "homotopyRemovalMLDegree",
    "getWitnessCollection",
    "saveWitnessCollectionConfiguration",
--Keys
    "MLDegrees",
    "TheVariety",
    "ThePoint",
    "TheDimension",
----Keys for numerical computation
    "WitnessPoints",
    "WitnessSets",
    "SortPoints",
    "WitnessRing",
    "VariableGroups",
    "Directory"
        }


exportMutable{
    "randomValue"---Gets a random complex number as the default
    }

--###################################
-- TYPE DEFINITIONS
--###################################
--SubspaceArrangement=new Type of MutableHashTable

RemovalMLDegree=new Type of MutableHashTable
MLDegreeVariety=new Type of MutableHashTable


MLDegreeWitnessSet=new Type of MutableHashTable
MLDegreeWitnessCollection=new Type of RemovalMLDegree

--###################################
-- METHODS
--###################################
randomVector=method(Options=>{		})
randomVector(ZZ):= o->n ->apply(n,i->randomValue())--list of length n of randomValue

-- ###
-- From an ideal, we consider an affine variety. 
-- ###
dim(RemovalMLDegree):=(M)->(
    if not member(TheDimension,keys M#TheVariety) then	error"TheDimension not specified";        
    M#TheVariety#TheDimension
    )
--M.TheVariety=(M#TheVariety#DefiningEquations//dim)

newMLDegreeVariety=method(Options=>{		})

newMLDegreeVariety(Ideal):= o->(I)->(
    numZ:=#gens ring I;
    new MLDegreeVariety from {
      DefiningEquations=> I,
      Hyperplanes=>for i to numZ+1 list randomVector(numZ),
      Data0=>randomVector(numZ),
      Data1=>randomVector(1+numZ)})

--######
newRemovalMLDegree=method(Options=>{   })
--######
--For symbolic compuation
newRemovalMLDegree(MLDegreeVariety):= o->(L)->newRemovalMLDegree(L,apply(#gens ring (L#DefiningEquations),i->1))--newRemovalMLDegree(L,{1,1,...,1})
newRemovalMLDegree(MLDegreeVariety,List):= o->(L,P)->new RemovalMLDegree from {
    ThePoint=>P,
    TheVariety=>L,
    MLDegrees=>{}}    
--For numerical computation
newRemovalMLDegree(MLDegreeWitnessCollection,List):= o->(WC,p)->return new RemovalMLDegree from {TheVariety=>WC,MLDegrees=>{},ThePoint=>p,WitnessSets=>{}}

solveRemovalMLDegree=method(Options=>{ })

--For symbolic computation 
solveRemovalMLDegree(RemovalMLDegree):= o->(M)->solveRemovalMLDegree(M,codim (M#TheVariety#DefiningEquations))---solveRemovalMLDegree(M,codimI) where codimI is the codimension of I. 
solveRemovalMLDegree(RemovalMLDegree,ZZ):= o->(M,codimI)->(
    (M#TheVariety).TheDimension=#gens ring (M#TheVariety#DefiningEquations)-codimI;
    apply(dimI+2,rk->solveRemovalMLDegree(rk,M,codimI)))
    
solveRemovalMLDegree(ZZ,RemovalMLDegree):= o->(rk,M)->(
    theI:=M#TheVariety#DefiningEquations;
--    apply(M#Saturate,i->theI=saturate(theI,i));    
    codimI:=codim theI;
    dimI:=#gens ring (M#TheVariety#DefiningEquations)-codimI;
    (M#TheVariety).TheDimension=dimI;        
    solveRemovalMLDegree(rk,M,codimI))
solveRemovalMLDegree(ZZ,RemovalMLDegree,ZZ):= o->(rk,M,codimI)->(
--    print o;
    P:=M#ThePoint;
    bigF:=(M#TheVariety#DefiningEquations)//gens//entries//flatten;
    zList:=gens ring first bigF;
    mlKK:=coefficientRing ring first bigF;
    mlR:=mlKK[zList,value("y")];
    bigF=apply(bigF,i->sub(i,mlR));
    zList=drop(gens mlR,-1); 
    theY:=last gens mlR;
    zPoint:= apply(zList,P,(z,p)->z-sub(p,mlR));
    print zPoint;
    usedH:=apply(rk-1,i->sum apply((M#TheVariety#Hyperplanes)_(i+1),zPoint,(c,z)->sub(c,mlR)*z  ));
    if rk==0 then littleF:={} else littleF={-theY+sum apply(first(M#TheVariety#Hyperplanes),zPoint,(c,z)->sub(c,mlR)*z)};
    I:=bigF|littleF|usedH;
--    print (toString I);
--    print (toString littleF);
--    print (toString usedH);    
    if rk=!=0 then (
	useVars:=zList|{theY};
	topJac:=matrix{M#TheVariety#Data1}) else if rk==0 
    then (
	useVars=zList;
    	topJac=matrix{M#TheVariety#Data0}); 
 --   print(1,topJac);
    Jac:=(matrix makeJac(I,useVars));   
--    print(2,Jac);
    augM:=sub(topJac,mlR)||(Jac*diagonalMatrix(useVars));
--    print(4,toString augM);
    if rk==0 then minorSize:=codimI else minorSize=codimI+rk;
    LC:=ideal(I)+minors(1+minorSize,augM);
    print ("minors"=>minorSize,Jac,augM); 
    sl:= minors(minorSize,Jac);
    apply(useVars,i->LC=saturate(LC,i));
    LC=saturate(LC,sl);
--    if member(Saturate,keys M) then (print "winKey";      apply(M#Saturate,i->LC=saturate(LC,sub(i,mlR)))	);
    theDegree:=degree LC;
--    print(codim LC,dim LC); 
    M#MLDegrees=append(M#MLDegrees,rk=>theDegree);
--    print theDegree;
    theDegree)
 


--Extracting information
removalMLDegree=method(Options=>{  })
removalMLDegree(RemovalMLDegree):= o->M->(
    if length (M#MLDegrees)<1 then return {};
    maxRK:=M#MLDegrees/toList/first//max;
    outMLDegree:=apply(maxRK+1,i->null);
    apply(M#MLDegrees/toList,i->outMLDegree=replace(i_0,i_1,outMLDegree));
    return outMLDegree
    )
--This method takes a RemovalMLDegree as an input and returns a list of kth removal ML degrees. If the same kth removal ML degree is computed multiple times then only of these will be considered. 

mlObstructionFunction=method(Options=>{  })
{*
mlObstructionFunction(RemovalMLDegree):= o->(M)->(
    if 0=!=length removalMLDegree(M) then error"Use mlObstructionFunction(RemovalMLDegree,ZZ). ";
    solveRemovalMLDegree(M);
    if member(TheDimension,keys(M#TheVariety))
    then theDim:=M#TheVariety#TheDimension
    else   theDim=M#TheVariety#DefiningEquations//dim;
    mlObstructionFunction(M,theDim)
    )
*}
mlObstructionFunction(RemovalMLDegree):= o->(M)->(
    if member(TheDimension,keys(M#TheVariety))
    then theDim:=M#TheVariety#TheDimension    else error"TheDimension not specified. Use mlObstructionFunction(M,d). ";
    mlObstructionFunction(M,theDim)
    )

mlObstructionFunction(RemovalMLDegree,ZZ):= o->(M,d)->(
    theDegrees:=removalMLDegree(M);
    if length theDegrees<1 then return 0;    
    (-1)^d*sum for i to #theDegrees-1 list if theDegrees_i=!= null then (-1)^i*theDegrees_i else 0
    )
--This functon returns an alternating sum of ML degrees \sum_k (-1)^k r_k(1,X) times (-1)^d.
----------------------------------------------------------------------------------------------------------------
checkZero=(aSol,eps)->if aSol/abs//min<eps then false else true
----------------------------------------------------------------------------------------------------------------
sortPointFunction=(aSol)->(if not (apply(aSol,i->{realPart i,imaginaryPart i}/abs//max)//min<1e-8) then true else false	    );
----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
newMLDegreeWitnessCollection=method(Options=>{}) 
----------------------------------------------------------------------------------------------------------------
newMLDegreeWitnessCollection(Ideal,ZZ):= o->(I,d)->(
	s:= temporaryFileName() | "/";                            
	mkdir s;
	newMLDegreeWitnessCollection(I,d,s))
----------------------------------------------------------------------------------------------------------------
newMLDegreeWitnessCollection(Ideal,ZZ,String):= o->(I,d,s)->(
    newKeys:={TheDimension=>d,WitnessSets=>{},Directory=>s,MLDegrees=>{}}; 
    WC:=new MLDegreeWitnessCollection from newKeys|pairs newMLDegreeVariety(I);
    WC.Hyperplanes=apply(WC#Hyperplanes,i->append(i,randomValue()));
    pVars:= gens ring (WC#DefiningEquations);
    lVars:=apply(#pVars+2,i->value("lagV"|i));
    bVars:=apply(d+1,i->value("bV"|i));
    kk:=coefficientRing ring (WC#DefiningEquations);
--    print pVars;
    WC.WitnessRing=(kk[pVars|{"y"}]**kk[lVars]**kk[bVars]);
    WC.DefiningEquations=sub(WC#DefiningEquations,WC#WitnessRing);
--    print keys WC;
    WC.SortPoints=(aSol)->sortPointFunction(aSol);
--    saveWitnessCollectionConfiguration(WC,s);
    return WC
    )  
----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
reclassifyWitnessPoints=method(Options=>{		})
--reclassifyWitnessPoints(ZZ,RemovalMLDegree,FunctionClosure):= o->(rk,M,sf)->()

reclassifyWitnessPoints(MutableHashTable,RR):= o->(WCM,eps)->(
    spf:=(aSol)->(if not(apply(aSol,i->{realPart i,imaginaryPart i}/abs//max)//min<eps) then true else false	    );
    reclassifyWitnessPoints(WCM,spf)    	);
reclassifyWitnessPoints(MutableHashTable,FunctionClosure):= o->(WCM,spf)->(
    apply(#WCM#WitnessSets,p->reclassifyWitnessPoints(WCM,spf,p))
	);    
reclassifyWitnessPoints(MutableHashTable,RR,ZZ):= o->(WCM,eps,p)->(
    spf:=(aSol)->not(if (apply(aSol,i->{realPart i,imaginaryPart i}/abs//max)//min<eps) then true else false	    );
    reclassifyWitnessPoints(WCM,spf,p)
    )
reclassifyWitnessPoints(MutableHashTable,FunctionClosure,ZZ):= o->(WCM,spf,p)->(
    WCM.SortPoints=spf;
    rk:=(WCM#WitnessSets)_p//first;
    WS:=(WCM#WitnessSets)_p//last;
    theSols:=first    (WS#WitnessPoints);
    checkSols:= apply(theSols,i->WCM#SortPoints(i));
--    print checkSols;
    tallySols:=(tally checkSols);
    theDegree:=if member(true,keys tallySols) then  tallySols#true else 0;
--    print 1;
    WS.Degree=theDegree   ; 
--    print theDegree;
    WS.WitnessPoints=(theSols,checkSols);
--    print checkSols;
    WCM.WitnessSets=replace(p,rk=>WS,WCM#WitnessSets);
--    print (p,rk,WS#Degree);
--    print (p,WCM#MLDegrees);
    WCM.MLDegrees=replace(p,rk=>theDegree,WCM#MLDegrees);
--    print (p,WCM#MLDegrees);
    theDegree
    )

newMLDegreeWitnessSet=method(Options=>{	})        
newMLDegreeWitnessSet(MLDegreeWitnessCollection):= o->(WC)->(
    apply(WC#TheDimension+2,rk->newMLDegreeWitnessSet(rk,WC)))    
newMLDegreeWitnessSet(ZZ,MLDegreeWitnessCollection,Boolean):= o->(rk,WC,doRunBertini)->(
    gV:=flatten entries basis({0,0,1},WC#WitnessRing);
    pV:=flatten entries basis({1,0,0},WC#WitnessRing);
    lV:=flatten entries basis({0,1,0},WC#WitnessRing);
    bigF:=(WC#DefiningEquations)//gens//entries//flatten;
    s:=addSlash(WC#Directory)|rk; 
    if not fileExists(s) then mkdir s;
    mlR:=WC#WitnessRing; 
    usedH:=apply(rk,i->-gV_i+sum apply((WC#Hyperplanes)_i,pV,(c,z)->sub(c,mlR)*z  ));    
    I:=bigF|usedH;
    if rk=!=0 then (
      useVars:=pV;
      topJac:=matrix{WC#Data1}) else if rk==0 
    then (
      useVars=drop(pV,-1);
      topJac=matrix{WC#Data0}); 
    Jac:=(matrix makeJac(I,useVars));   
    augM:=sub(topJac,mlR)||(Jac*diagonalMatrix(useVars));
--    print(4,toString augM);
--    print useVars;
    useLagVars:=for i to numrows augM-1 list lV_i;
    useGraphVars:=for i to #usedH-1 list gV_i;    
    WS:=new MLDegreeWitnessSet from {VariableGroups=>{useLagVars,useVars,useGraphVars}};
    allEqs:={bigF,usedH,reverse flatten entries ((matrix{useLagVars})*augM)};
    WS#DefiningEquations=allEqs;
    if doRunBertini 
    then(
      startTime:=currentTime();
      makeB'InputFile(s,
    	B'Polynomials=>flatten allEqs,
    	HomVariableGroup=>useLagVars,
    	AffVariableGroup=>useVars,
	ParameterGroup=>gV,
    	B'Configs=>{"ParameterHomotopy"=>1,"UseRegeneration"=>1,"MPTYPE"=>2,"PrintPathProgress"=>1000}
    	); 
      try runBertini(s,PreparePH2=>true)
      then       (endTime:=currentTime();
      print("Bertini run time ", endTime-startTime))
      else print("Bertini did not run as expected.")      );
    try theSols:=importSolutionsFile(s,NameSolutionsFile=>"start")
      then "Success"       else theSols={};
--    theSols:=importSolutionsFile(s,NameSolutionsFile=>"start");
    checkSols:= apply(theSols,i->WC#SortPoints(i));
--    print checkSols;
    tallySols:=(tally checkSols);
    theDegree:=if member(true,keys tallySols) then  tallySols#true else 0;
    WS#Degree=theDegree   ; 
    WS#WitnessPoints=(theSols,checkSols);
    WC#WitnessSets=append(WC#WitnessSets,rk=>WS);
    WC.MLDegrees=append(WC#MLDegrees,rk=>theDegree);
    theDegree
    )
newMLDegreeWitnessSet(ZZ,MLDegreeWitnessCollection):= o->(rk,WC)->newMLDegreeWitnessSet(rk,WC,true)
----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
homotopyRemovalMLDegree=method(Options=>{		}) 
----------------------------------------------------------------------------------------------------------------
homotopyRemovalMLDegree(RemovalMLDegree):= o->(M)->(
    apply(M#TheVariety#TheDimension+2,rk->homotopyRemovalMLDegree(rk,M)))                  
homotopyRemovalMLDegree(ZZ,RemovalMLDegree):= o->(rk,M)->(    
  WC:=M#TheVariety;
  s:=addSlash(WC#Directory)|rk;
  WS:=null;
  if member(rk,WC#WitnessSets/toList/first) then( 
      indexWS:=position(WC#WitnessSets/toList/first,i->i==rk);
      WS=((WC#WitnessSets)_indexWS)//toList//last
      ) else return {};---Returns an empty list if that witness set is not present.  
  if rk=!=0 
  then(
    theP:=apply(WS#VariableGroups#1,append(M#ThePoint,0),(i,j)->i=>j);
    theSP:=importParameterFile(addSlash(WC#Directory)|rk,NameParameterFile=>"start_parameters");
    bV:=WS#VariableGroups#2;
    theB:=apply(#bV,i->bV_i=>0);
    theFP:=apply(    WS#DefiningEquations#1,i->sub(i,theP|theB));
--    print theFP;
    theFP=theFP|drop(theSP,#theFP);
    writeParameterFile(s,theFP,NameParameterFile=>"final_parameters");
    startTime:=currentTime();
    try runBertini(s)  
    then (
      endTime:=currentTime();
      print("Bertini run time ", endTime-startTime);
      try theSols:=importSolutionsFile(s,NameSolutionsFile=>"nonsingular_solutions")
      then "Success" 
      else theSols={}
      )
    else print("Bertini did not run as expected.")          
  )else if rk===0 then(
    if fileExists(s|"start") 
    then theSols=importSolutionsFile(s,NameSolutionsFile=>"start")
    else theSols={}); 
  checkSols:= apply(theSols,i->WC#SortPoints(i));
--  print checkSols;
  tallySols:=(tally checkSols);
  theDegree:=if member(true,keys tallySols) then  tallySols#true else 0;
  M#MLDegrees=append(M#MLDegrees,rk=>theDegree)   ; 
  WS.WitnessPoints=(theSols,checkSols);
  M#WitnessSets=append(M#WitnessSets,rk=>WS);
  return theDegree
  )
----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
saveKeys={Hyperplanes,TheDimension,Data0,Data1,MLDegrees}
----------------------------------------------------------------------------------------------------------------
saveWitnessCollectionConfiguration=method(Options=>{	})        
----------------------------------------------------------------------------------------------------------------
saveWitnessCollectionConfiguration(MLDegreeWitnessCollection,String):= o->(WC,s)->(
  apply(saveKeys,i->(
    newFile := openOut (addSlash(s)|"save_WC_"|toString i); 
    newFile << toString (WC#i) << endl;
    close newFile    
    ));)                 
----------------------------------------------------------------------------------------------------------------
saveWitnessCollectionConfiguration(MLDegreeWitnessCollection):= o->(WC)->saveWitnessCollectionConfiguration(WC,WC#Directory)
----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
getWitnessCollection=method(Options=>{	})        
----------------------------------------------------------------------------------------------------------------
getWitnessCollection(Ideal,ZZ,String):= o->(I,d,s)->(
  WC:=newMLDegreeWitnessCollection(I,d,s);
  apply(saveKeys,i->(
    savedFile := (addSlash(s)|"save_WC_"|toString i); 
    WC#i=value get savedFile));
  apply(WC#TheDimension+2,rk->newMLDegreeWitnessSet(rk,WC,false));
  return WC
  )                        




{*
restart
loadPackage("MaximumlikelihoodObstructionFunction",Configuration=>{"RandomCoefficients"=>CC},Reload=>true)
R=CC[p1,p2]
I=ideal((p1-1)^3+(p2-1)^4)
WC=newMLDegreeWitnessCollection(I,1,theDir)
newMLDegreeWitnessSet(0,WC)
peek WC
reclassifyWitnessPoints(WC,1e-20,0)
reclassifyWitnessPoints(WC,1e-10,0)
peek 
WC#WitnessSets
removalMLDegree(WC)
WC1=newMLDegreeWitnessCollection(I,1,theDir)
newMLDegreeWitnessSet(WC1)
reclassifyWitnessPoints(WC1,1e-20)
reclassifyWitnessPoints(WC1,1e-9)
M1=newRemovalMLDegree(WC1,{1,1})
peek M1
homotopyRemovalMLDegree(M1)
removalMLDegree(WC1)
removalMLDegree(M1)
mlObstructionFunction(M1)
*}

{*
restart
loadPackage("MaximumlikelihoodObstructionFunction",Configuration=>{"RandomCoefficients"=>CC_300},Reload=>true)
printingPrecision=300
R=CC[p12,p13,p22,p23,p33]
f=det matrix{
    {1,p12,p13},
    {p12,p22,p23},
    {p13,p23,p33}}
I=ideal(f)
WC=newMLDegreeWitnessCollection(I,4,theDir)
newMLDegreeWitnessSet(0,WC)
newMLDegreeWitnessSet(1,WC)
newMLDegreeWitnessSet(2,WC)
newMLDegreeWitnessSet(3,WC)
newMLDegreeWitnessSet(4,WC)
newMLDegreeWitnessSet(5,WC)

M=newRemovalMLDegree(WC,{1,1,1,1,1})
peek WC
keys WC
apply({0,1,2,3,4,5},i-> sum apply(
	(WC#i#WitnessPoints)_0,(WC#i#WitnessPoints)_1,(j,k
	    )->if k then 1 else 0))
{0, -16, 47, -49, 21, -3}//sum
apply({0,1,2,3,4,5},i->homotopyRemovalMLDegree(i,M,WC));
removalMLDegree M
mlObstructionFunction(M,5)
{0, 16, 47, 49, 19, 1}
methods RemovalMLDegree
{0, 16, 47, 49, 19, 1}=>{1,1,1,1,1}
{0, 16, 47, 49, 21, 2}=>{1,1,1,1,2}
*}



{* 
restart
loadPackage("MaximumlikelihoodObstructionFunction",Configuration=>{"RandomCoefficients"=>CC_300},Reload=>true)
printingPrecision=300
R=CC[p1,p2,p3,p4,        p5,	p6]
f=det matrix{
    {p1,p2,p3,1},
    {p2,p3,1,p4},
    {p3,1,p4,p5},
    {1,p4,p5,p6}}
I=ideal(f)
theDim=#gens R-1
WC=newMLDegreeWitnessCollection(I,theDim,theDir)
apply(theDim+1,i->newMLDegreeWitnessSet(i,WC))

newMLDegreeWitnessSet(1,WC)

newMLDegreeWitnessSet(2,WC)
newMLDegreeWitnessSet(3,WC)
newMLDegreeWitnessSet(4,WC)
newMLDegreeWitnessSet(5,WC)
newMLDegreeWitnessSet(6,WC)
newMLDegreeWitnessSet(7,WC)
newMLDegreeWitnessSet(8,WC)

*}

--##########################################################################--
-- INTERNAL METHODS
--##########################################################################--
----------------------------------------
parString=(aString)->("("|toString(aString)|")");
addSlash=(aString)->(
    if aString_-1===" " then error (aString|" cannot end with whitespace.");
    if aString_-1=!="/" then aString=aString|"/";
    return aString    )
newHyperplanes=A->for i to (numColumns A)+1 list randomVector(numRows A)
makeJac=(system,unknowns)->(--it is a list of lists of partial derivatives of a polynomial
         for i in system list for j in unknowns list  diff(j,i))


beginDocumentation()

load "./MLOF_Supplement/DOC_MLOF.m2";

TEST///
--load concatenate(MultiprojectiveWitnessSets#"source directory","./AEO/TST/Example1.tst.m2")
///


end

