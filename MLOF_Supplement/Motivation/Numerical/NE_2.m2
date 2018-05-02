--These files have two parts. 
--The first part gets a witness collection from the Motivation directory. 
--The second part can be used to compute the witness collection.
--{0, 16, 47, 49, 21, 3}
---Solving. 
theExample="NE_2"

loadPackage("MaximumlikelihoodObstructionFunction",Reload=>true,Configuration=>{"RandomCoefficients"=>CC}    )
s=temporaryFileName()|"/"
mkdir s
s="/Users/jo/Documents/GitStuff/ReciprocalLinearSpacesAtInfinity/MLOF_Supplement/Motivation/Numerical/"|theExample

printingPrecision=300
R=CC[p1,p2,p3,p4,p5]
I=ideal {det matrix{{p1,p2,p3},{p2,p4,p5},{p3,p5,1}}}
d=#gens R-numgens I
WC=getWitnessCollection(I,d,s)                          
peek WC
thePoints={
  {1,2,3,5,7},--Rank 3
  {2,1,1,1,1},--Rank 2
  {1,1,1,1,1}}--Rank 1

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
s="/Users/jo/Documents/GitStuff/ReciprocalLinearSpacesAtInfinity/MLOF_Supplement/Motivation/Numerical/"|theExample

printingPrecision=300
R=CC[p1,p2,p3,p4,p5]
I=ideal {det matrix{{p1,p2,p3},{p2,p4,p5},{p3,p5,1}}}
d=#gens R-numgens I
WC=newMLDegreeWitnessCollection(I,d,s)                          
WC#Data1={.19113580965665+.935103300329118*ii, .482388990178692+.843487507389781*ii, .490690501911577+.245366076512194*ii, .251415040009651+.0498365746831353*ii, .166223555598076+.497213293915745*ii, .758876282349089+.883677494789585*ii};
WC#Data0={.675686886661785+.804558923530637*ii, .1342280558682+.144168130915478*ii, .358634198495929+.193358016804975*ii, .714304084371768+.654323025214591*ii, .200749167629783+.0562389581539768*ii};
WC#Hyperplanes={
    {.810501765785929+.699011464747371*ii, .233716556838891+.313873549074706*ii, .357126110039482+.0508569512708977*ii, .970271701838245+.26489734813169*ii, .288003618519728+.548819039319049*ii, .290864686132322+.881971629224713*ii},
{.402297068224863+.632569703259367*ii, .0112559325320295+.523925658316823*ii, .866659426462128+.903392581287414*ii, .331788640417778+.197146722562498*ii, .594599056775829+.934726271695532*ii, .166064374797875+.857101546115526*ii},
{.614255588001682+.365990201954648*ii, .765990878595206+.0598852452689215*ii, .499593140294508+.562675931669175*ii, .345185346906543+.198576252025757*ii, .705016369719636+.922029107313792*ii, .225007931708862+.815969323626468*ii},
{.14888111174523+.2282770725131*ii, .395300490899302+.448793082959414*ii, .0445329393506518+.46339681381017*ii, .659619403605046+.690024750201218*ii, .189428084921205+.377396382893493*ii, .111382611430274+.58303347110279*ii},
{.391802717342185+.259897667292941*ii, .4686439535312+.99080335518113*ii, .159990150707672+.879827994536015*ii, .254764098218317+.550870345950083*ii, .720558449791416+.865795249946525*ii, .601156979292572+.745836482258895*ii},
{.154832405270898+.682615018855798*ii, .942915294245362+.539813040720056*ii, .298176396429141+.334382328157683*ii, .417202601000835+.842327209458092*ii, .520353081450659+.261629766436895*ii, .240885893077631+.834047478006461*ii},
{.938961207848625+.831015181123528*ii, .595241315161505+.142486298001414*ii, .554694615670221+.58678779023901*ii, .0280572709633171+.699132179542595*ii, .709668737999846+.0774629824512825*ii, .0769838358857137+.367033937478421*ii}};
--toString WC#Hyperplanes
--toString WC#Data1
--toString WC#Data0
--newMLDegreeWitnessSet(WC)                              
WC.Directory=s
newMLDegreeWitnessSet(WC)
saveWitnessCollectionConfiguration(WC,s)


