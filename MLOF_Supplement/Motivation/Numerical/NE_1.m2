--These files have two parts. 
--The first part gets a witness collection from the Motivation directory. 
--The second part can be used to compute the witness collection.
---Solving. 
loadPackage("MaximumlikelihoodObstructionFunction",Reload=>true,Configuration=>{"RandomCoefficients"=>CC}    )
s=temporaryFileName()|"/"
mkdir s
s="/Users/jo/Documents/GitStuff/ReciprocalLinearSpacesAtInfinity/MLOF_Supplement/Motivation/Numerical/NE_1"
printingPrecision=300
R=CC[p1,p2,p3]
I=ideal {(p1-1)^2-(p2-1)^2*(p3-1)}
d=2
WC=getWitnessCollection(I,d,s)                          
peek WC
thePoints={
  {3,2,1},--S0
  {3,3,2},--S1
  {1,1,2},--S2
  {1,1,1}}--S3
apply(thePoints,P->(
    startTime=currentTime();
    M=newRemovalMLDegree(WC,P);
    win1'=homotopyRemovalMLDegree M; 
    win1=removalMLDegree M;
    L=newMLDegreeVariety(I);
    win2:=mlObstructionFunction M;
    endTime=currentTime();
--	print 1;
     (P=>(endTime-startTime,win1,win2))
     )    )
netList oo

end
--Preprocessing
restart
loadPackage("MaximumlikelihoodObstructionFunction",Reload=>true,Configuration=>{"RandomCoefficients"=>CC}    )
s=temporaryFileName()|"/"
mkdir s
s="/Users/jo/Documents/GitStuff/ReciprocalLinearSpacesAtInfinity/MLOF_Supplement/Motivation/Numerical/NE_1"
printingPrecision=300
R=CC[p1,p2,p3]
I=ideal {(p1-1)^2-(p2-1)^2*(p3-1)}
d=2
WC=newMLDegreeWitnessCollection(I,d,s)                          
WC#Data1={.318344+.730461*ii, .432042+.378417*ii, .782374+.952397*ii,      .873999+.511173*ii}
WC#Data0={.976133+.564431*ii, .470859+.33303*ii, .496929+.577981*ii}
WC#Hyperplanes={{.229183+.879231*ii, .415547+.0727349*ii, .645009+.370209*ii, .751082+.842956*ii}, {.567109+.772207*ii, .612607+.722364*ii, .685995+.629651*ii, .0670266+.320013*ii}, {.206756+.461355*ii, .749389+.536174*ii, .125451+.812638*ii, .254213+.642019*ii}, {.68259+.633135*ii, .148048+.897712*ii, .97991+.54262*ii, .442083+.0530267*ii}, {.159642+.449563*ii, .325694+.764952*ii, .938067+.244003*ii, .819177+.354919*ii}}
--toString WC#Hyperplanes
--toString WC#Data1
--toString WC#Data0
--newMLDegreeWitnessSet(WC)                              
WC.Directory=s
newMLDegreeWitnessSet(WC)
saveWitnessCollectionConfiguration(WC,s)
