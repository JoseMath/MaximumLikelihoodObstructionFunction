loadPackage("MaximumlikelihoodObstructionFunction",Reload=>true)
--Whitney paraguas
R=QQ[p1,p2,p3]
F={(p1-1)^2-(p2-1)^2*(p3-1)}
thePoints={
  {3,2,1},--S0
  {3,3,1},--S1
  {1,1,2},--S2
  {1,1,1}}--S3

apply(thePoints,P->(
--	print 1;
    startTime=currentTime();
--	print 1;
    L=newMLDegreeVariety(ideal F);
--	print 1;
    M=newRemovalMLDegree(L,P);
--	print 1;
    win1:=solveRemovalMLDegree M;
--	print 1;
    win2:=mlObstructionFunction M;
--	print 1;
    endTime=currentTime();
--	print 1;
     (P=>(endTime-startTime,win1,win2))
     )    )
