--These files have two parts. 
--The first part gets a witness collection from the Motivation directory. 
--The second part can be used to compute the witness collection.
---Solving. 
theExample="NE_5"

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

end
--Preprocessing
restart
theExample="NE_5"

loadPackage("MaximumlikelihoodObstructionFunction",Reload=>true,Configuration=>{"RandomCoefficients"=>CC}    )
s=temporaryFileName()|"/"
mkdir s
--s="/Users/jo/Documents/GitStuff/ReciprocalLinearSpacesAtInfinity/MLOF_Supplement/Motivation/Numerical/"|theExample

printingPrecision=300
R=CC[p1,p2,p3,p4,p5,p6]
I=ideal {det matrix{
	{p1,p2,p3,p4},
	{p2,p3,p4,p5},
	{p3,p4,p5,p6},
	{p4,p5,p6,1}}}

d=#gens R-numgens I
WC=newMLDegreeWitnessCollection(I,d,s)                          

--for i in WC#Hyperplanes do print toString i
--print toString WC#Data1
--print toString WC#Data0
WC.Hyperplanes={
{.863604635007954+.504469618460703*ii, .558279522564033+.035146543092516*ii, .118308815282904+.294708517275343*ii, .893310479021742+.212978921743777*ii, .491188100158899+.294774838019557*ii, .565170422101381+.650441390170898*ii, .532351711830692+.343491591332896*ii},
{.972781747638505+.802303257318813*ii, .0847418413239819+.470575567102895*ii, .172958780693913+.0142470739108336*ii, .513170316269748+.0300960970883815*ii, .172978233708103+.527536076272898*ii, .537488599688496+.042828839590881*ii, .748966300274089+.942438742232203*ii},
{.871977496949491+.0311056080216155*ii, .510786787103613+.302516638706013*ii, .0634049173998081+.172667074158338*ii, .0140264197921858+.714814808280168*ii, .270400705489284+.16688904079982*ii, .127534773891495+.959129185922611*ii, .339244103619025+.0194106554727191*ii},
{.723795380843429+.402859069046906*ii, .491416916943226+.469757567524431*ii, .379367578018462+.385918012534295*ii, .400891071954386+.887939126712501*ii, .0597487969089885+.388364456069015*ii, .140576672364323+.491146905763802*ii, .154537901036053+.607959700841086*ii},
{.403121739588835+.19075134608003*ii, .108160506591593+.118675756913575*ii, .946382733586289+.55251880954286*ii, .979865128991216+.138104597885224*ii, .601675355596689+.536745857694205*ii, .448182183287909+.762720868070346*ii, .19862960280149+.869323444308189*ii},
{.770181253202099+.778215312629569*ii, .497278381336345+.260115977831259*ii, .917615281545755+.724373617845925*ii, .794868856501123+.0884765414466813*ii, .805009071113024+.454835978006802*ii, .926668861131034+.509616732416782*ii, .440314564377563+.925255629918619*ii},
{.936459606786381+.278310217571028*ii, .392669195530609+.372537824849004*ii, .775861779827556+.230618978107689*ii, .311165937518994+.775958890544416*ii, .240673711757541+.72329176439663*ii, .75524773546996+.336729423956142*ii, .975085571425044+.501734769549083*ii},
{.491444002532311+.130626894698217*ii, .383824326856809+.671532534392353*ii, .766864394732969+.932746201953048*ii, .461622388064771+.00688037503634331*ii, .134822791137508+.459883529109502*ii, .587576963720767+.997835959458902*ii, .206750423750497+.298769238340229*ii}}
WC.Data1={.374107787897095+.482221677769261*ii, .139084080846696+.547582564916353*ii, .154875784374677+.850614485033597*ii, .95208582132038+.539098285127349*ii, .0651143065711171+.815594721549913*ii, .782088465349817+.68174269237269*ii, .707648093003382+.304773823794988*ii}
WC.Data0={.956938827072436+.761815438544562*ii, .526882395011883+.302996979706311*ii, .412548456150462+.270216730299861*ii, .658200463107726+.990655070012695*ii, .664461397109036+.822054293981365*ii, .369548022471534+.588127137459249*ii}

WC.Directory=s
newMLDegreeWitnessSet(WC)
saveWitnessCollectionConfiguration(WC,s)
