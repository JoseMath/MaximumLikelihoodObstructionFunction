--These files have two parts. 
--The first part gets a witness collection from the Motivation directory. 
--The second part can be used to compute the witness collection.
---Solving. 
theExample="NE_4"

loadPackage("MaximumlikelihoodObstructionFunction",Reload=>true,Configuration=>{"RandomCoefficients"=>CC}    )
s=temporaryFileName()|"/"
mkdir s
s="/Users/jo/Documents/GitStuff/ReciprocalLinearSpacesAtInfinity/MLOF_Supplement/Motivation/Numerical/"|theExample
printingPrecision=300
R=CC[p1,p2,p3,p4,p5]
 I=ideal {det matrix{
     {p1,p2,1,p5},
     {p2,p1,p5,1},
     {1,p5,p3,p3},
     {p5,1,p3,p4}}}
d=#gens R-numgens I
WC=getWitnessCollection(I,d,s)                          
thePoints=apply({{p1=>2,p2=>-2,p3=>3,p4=>27/5,p5=>-1},
{p1=>2,p2=>2,p3=>1/2,p4=>3,p5=>1},
{p1=>3,p2=>-2,p3=>2,p4=>2,p5=>1},
{p1=>2,p2=>2,p3=>2,p4=>2,p5=>1},
{p1=>1,p2=>1,p3=>1,p4=>1,p5=>1}},i->apply(i/last,j->sub(j,CC)))

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

(Bertini run time , 113)
Warning: The HomVariableGroup is written first and then the AffVariableGroup is written second.
(Bertini run time , 1030)
Warning: The HomVariableGroup is written first and then the AffVariableGroup is written second.
(Bertini run time , 1237)
Warning: The HomVariableGroup is written first and then the AffVariableGroup is written second.
(Bertini run time , 710)
Warning: The HomVariableGroup is written first and then the AffVariableGroup is written second.
(Bertini run time , 249)
Warning: The HomVariableGroup is written first and then the AffVariableGroup is written second.
(Bertini run time , 17)

      +--------------------------------------------------------+
o12 = |{2, -2, 3, 5.4, -1} => (10, {0, 76, 156, 108, 30, 2}, 0)|
      +--------------------------------------------------------+
      |{2, 2, .5, 3, 1} => (6, {0, 76, 156, 108, 30, 2}, 0)    |
      +--------------------------------------------------------+
      |{3, -2, 2, 2, 1} => (8, {0, 76, 156, 108, 30, 2}, 0)    |
      +--------------------------------------------------------+
      |{2, 2, 2, 2, 1} => (12, {0, 76, 156, 108, 30, 2}, 0)    |
      +--------------------------------------------------------+
      |{1, 1, 1, 1, 1} => (8, {0, 76, 156, 104, 26, 1}, 1)     |
      +--------------------------------------------------------+


 {0, 76, 156, 108, 32, 4}
 {0, -76, 156, -108, 32, -4}
end
--Preprocessing
restart
loadPackage("MaximumlikelihoodObstructionFunction",Reload=>true,Configuration=>{"RandomCoefficients"=>CC}    )
s=temporaryFileName()|"/"
mkdir s
s="/Users/jo/Documents/GitStuff/ReciprocalLinearSpacesAtInfinity/MLOF_Supplement/Motivation/Numerical/"|theExample

printingPrecision=300
R=CC[p1,p2,p3,p4,p5]
I=ideal {det matrix{
	{p1,p2,1,p5},
	{p2,p1,p5,1},
	{1,p5,p3,p3},
	{p5,1,p3,p4}}}
d=#gens R-numgens I
WC=newMLDegreeWitnessCollection(I,d,s)                          
--for i in WC#Hyperplanes do print toString i
--print toString WC#Data1
--print toString WC#Data0
WC.Hyperplanes={
{.404683598354642+.651089307210356*ii, .0270026772626002+.619356716744283*ii, .973091810949759+.87992498076463*ii, .885346201250907+.594002248812794*ii, .308647444252525+.60260058585148*ii, .857908567894742+.516757398639123*ii},
{.296963805294701+.578811580974739*ii, .772972394878659+.721394356570421*ii, .199945036311587+.412238889300797*ii, .0668048165493365+.424116733335587*ii, .967156856963974+.241066160380166*ii, .669895324348012+.900032715290284*ii},
{.519497690568545+.469516471935929*ii, .816412693078858+.769666676542302*ii, .940652190871592+.363676446156691*ii, .153523269641693+.538662771516429*ii, .825360389781554+.518672698437896*ii, .71189515643541+.631384318341153*ii},
{.93954954225564+.689975843230882*ii, .308329992505337+.73547808956304*ii, .0548351958901906+.0858060688518425*ii, .185934784261474+.0946520219771586*ii, .140638485896219+.40687036350765*ii, .241577866587824+.68171037896817*ii},
{.875742914273261+.508406858971348*ii, .94037131707586+.726814800037344*ii, .94473365444556+.814744354748193*ii, .038131200024379+.71223673644495*ii, .329081954656209+.240167844113876*ii, .771693418000762+.648421971087825*ii},
{.569274765294051+.464927971289213*ii, .739724661197907+.248182675566057*ii, .552263731900957+.859806163893949*ii, .363941936808009+.774547395247285*ii, .00116006232393417+.805024434134287*ii, .386033512632351+.455993637286812*ii},
{.700112264199811+.694793143563107*ii, .945580596597809+.702804231410869*ii, .415959839141936+.210040577542976*ii, .51849962165238+.246117527205634*ii, .607798315664427+.992123961883083*ii, .849357265141165+.0814470275027265*ii}}
WC.Data1={.936889268054669+.768092882167579*ii, .220291904750972+.325681596249602*ii, .524906858211409+.980076835830625*ii, .823200306487354+.577754715101639*ii, .14049925332264+.214714128087272*ii, .494901540087479+.594394090432666*ii}
WC.Data0={.600115212414696+.155522867620408*ii, .480768957893592+.524642695340341*ii, .106092633070783+.820858712550693*ii, .550096062509066+.561080279940632*ii, .938108635454372+.114667045465803*ii}
WC.Directory=s
newMLDegreeWitnessSet(WC)
saveWitnessCollectionConfiguration(WC,s)

---We want to understand the Whitney stratification.
R=QQ[p1,p2,p3,p4,p5]
 I=ideal {det matrix{
     {p1,p2,1,p5},
     {p2,p1,p5,1},
     {1,p5,p3,p3},
     {p5,1,p3,p4}}}
decompose I
S2=primaryDecomposition ideal singularLocus I
S2rad=S2/radical
decompose I
S2rad/codim
last S2rad==S3
for i in S2rad list for j in S2rad list i==j
netList S2rad


netList oo
S3=ideal singularLocus S2rad_0

+S2rad_1
decompose S3

S3=ideal (p5 - 1, p1 - p2, - p3 + p4, p2*p3 - 1)
{1,1,1,1,1}


      +---------------------------------------------------------+
      |                                 2                       |
o22 = |ideal (- p5 - 1, - p1 - p2, p2*p3  - p2*p3*p4 - 3p3 - p4)|
      +---------------------------------------------------------+
      |ideal (- p5 + 1, - p1 + p2, - p2*p3 + 1)                 |
      +---------------------------------------------------------+
      |ideal (- p5 + 1, - p3 + p4, - p1*p4 - p2*p4 + 2)         |
      +---------------------------------------------------------+
      |ideal (- p5 + 1, p3 - p4, p1 - p2)                       |


for i in S2rad do print toString i
pSub={p1=>2,p2=>-2,p3=>3,p4=>27/5,p5=>-1}
sub(ideal(-p5-1,-p1-p2,p2*p3^2-p2*p3*p4-3*p3-p4),pSub)

pSub={p1=>2,p2=>2,p3=>1/2,p4=>3,p5=>1}
sub(ideal(-p5+1,-p1+p2,-p2*p3+1),pSub)

pSub={p1=>3,p2=>-2,p3=>2,p4=>2,p5=>1}
sub(ideal(-p5+1,-p3+p4,-p1*p4-p2*p4+2),pSub)

pSub={p1=>2,p2=>2,p3=>2,p4=>2,p5=>1}
sub(ideal(-p5+1,p3-p4,p1-p2),pSub)



