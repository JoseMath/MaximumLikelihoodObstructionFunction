--These files have two parts. 
--The first part gets a witness collection from the Motivation directory. 
--The second part can be used to compute the witness collection.
---Solving. 
theExample="NE_3"

loadPackage("MaximumlikelihoodObstructionFunction",Reload=>true,Configuration=>{"RandomCoefficients"=>CC}    )
s=temporaryFileName()|"/"
mkdir s
s="/Users/jo/Documents/GitStuff/ReciprocalLinearSpacesAtInfinity/MLOF_Supplement/Motivation/Numerical/"|theExample
printingPrecision=300
R=CC[p1,p2,p3,p4]
I=ideal {det matrix{
	{p1,p2,p3},
	{p2,p3,p4},
	{p3,p4,1}}}
d=#gens R-numgens I
WC=getWitnessCollection(I,d,s)                          
thePoints={
  {1,2,3,5},--rank 3
  {2,1,1,1},--Rank 2
  {1,1,1,1}}--Rank 1

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
R=CC[p1,p2,p3,p4]
I=ideal {det matrix{
	{p1,p2,p3},
	{p2,p3,p4},
	{p3,p4,1}}}
d=#gens R-numgens I
WC=newMLDegreeWitnessCollection(I,d,s)                          
--for i in WC#Hyperplanes do print toString i
--print toString WC#Data1
--print toString WC#Data0
WC.Hyperplanes={{.836577757857502+.65138853220081*ii, .602986610469868+.0224031987686349*ii, .261398627747969+.390066126837452*ii, .207177720679944+.170867809996595*ii, .0401624243163546+.326190567058788*ii},
{.139130208774658+.268637478994035*ii, .429502936191858+.359874969021436*ii, .343842650061254+.477763812970088*ii, .454911660663154+.819802744765628*ii, .814728342303367+.285936884667439*ii},
{.508898658021246+.0560685003090858*ii, .118808938293562+.0404839394819848*ii, .823645873473805+.00548421808478361*ii, .517879988209287+.78413364683507*ii, .409582162041713+.966966805630896*ii},
{.0149696309914841+.556001199597431*ii, .0220325778179703+.768979844295788*ii, .207216898881612+.772710514786686*ii, .0780500071158976+.997942662630897*ii, .796622250372538+.538021915294525*ii},
{.00691064885550829+.275887569626361*ii, .263953394714634+.612847042545502*ii, .721763467067246+.642430670795737*ii, .0219409430908032+.196557843995854*ii, .00464058852828719+.357997723832319*ii},
{.484549306075739+.730131458210868*ii, .241882114139105+.818567561879329*ii, .539355721291226+.793181638475794*ii, .753522920114049+.0219779615078567*ii, .368640513239587+.138066989088306*ii}}
WC.Data1={.564018440715025+.27768327769232*ii, .592680231232282+.728689440243153*ii, .187576421016555+.617425043778195*ii, .631231762760909+.338123909353498*ii, .821996612367309+.289545493484167*ii}
WC.DAta0={.45292687726111+.808624388662052*ii, .71339094047815+.662636877852282*ii, .0588198179154032+.853534202028084*ii, .0150499307528988+.462398622139109*ii}
WC.Directory=s
newMLDegreeWitnessSet(WC)
saveWitnessCollectionConfiguration(WC,s)





