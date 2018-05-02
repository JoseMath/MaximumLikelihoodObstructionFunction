loadPackage("MaximumlikelihoodObstructionFunction",Reload=>true)
--Symmetric
---Restart
R=QQ[p2,p3,p4,p5,p6]
F={det matrix{{1,p2,p3},{p2,p4,p5},{p3,p5,p6}}}
thePoints={
  {7,5,3,2,1},--Rank 3
  {1,1,1,1,2},--Rank 2
  {1,1,1,1,1}}--Rank 1

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
