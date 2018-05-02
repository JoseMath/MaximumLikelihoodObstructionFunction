
doc /// --MaximumlikelihoodObstructionFunction
    Key
        MaximumlikelihoodObstructionFunction 
    Headline
        Determine Euler obstructions using maximum likelihood degrees.
    Description
      Text
        This package provides several routines for determining Euler obstructions and maximum likelihood degrees.	
		 
      Example
	R=QQ[p1,p2,p3]
	f=sub(p1^2-p2^2*p3,{p1=>p1-1,    p2=>p2-1,    p3=>p3-1})
	L=newMLDegreeVariety(ideal {f})	
	P={2,3,4}
	M1=newRemovalMLDegree(L,P)
	solveRemovalMLDegree M1
	mlObstructionFunction M1
    	peek M
      Example
	R=QQ[p1,p2,p3]
	f=sub(p1^2-p2^2*p3,{p1=>p1-1,    p2=>p2-1,    p3=>p3-1})
	L=newMLDegreeVariety(ideal {f})	
	M=newRemovalMLDegree(L)
    	M#ThePoint--Default point is {1,1,...,1}.
    	solveRemovalMLDegree(0,M)
    	solveRemovalMLDegree(1,M)
    	solveRemovalMLDegree(2,M)
    	solveRemovalMLDegree(3,M)
	M#MLDegrees
	removalMLDegree(M)
	peek M#TheVariety
	mlObstructionFunction(M)
      Example
        loadPackage("MaximumlikelihoodObstructionFunction",Configuration=>{"RandomCoefficients"=>CC},Reload=>true)
	s = temporaryFileName() | "/"                            
	mkdir s
	printingPrecision=300
	R=CC[p1,p2,p3]
	f=sub(p1^2-p2^2*p3,{p1=>p1-1,    p2=>p2-1,    p3=>p3-1})	
	I=ideal(f)
    	d=2;
	WC=newMLDegreeWitnessCollection(I,d,s)
	newMLDegreeWitnessSet(WC)
    	P={1,1,1}
	M=newRemovalMLDegree(WC,P)
	homotopyRemovalMLDegree M
	removalMLDegree M
	mlObstructionFunction M
	reclassifyWitnessPoints(M,1e-300)                                                                                 
	reclassifyWitnessPoints(M,1e10)                                                       
	reclassifyWitnessPoints(M,1e-6)                                                       
	mlObstructionFunction M
      Example
        loadPackage("MaximumlikelihoodObstructionFunction",Configuration=>{"RandomCoefficients"=>CC},Reload=>true)
	s = temporaryFileName() | "/"                            
	mkdir s
	R=CC[p1,p2,p3]
	printingPrecision=300
	f=sub(p1^2-p2^2*p3,{p1=>p1-1,    p2=>p2-1,    p3=>p3-1})	
	I=ideal(f)
    	d=2;
	WC=newMLDegreeWitnessCollection(I,d,s)
    	saveWitnessCollectionConfiguration(WC,s)
	newMLDegreeWitnessSet(0,WC)
	newMLDegreeWitnessSet(1,WC)
	newMLDegreeWitnessSet(2,WC)
	newMLDegreeWitnessSet(3,WC)
	M=newRemovalMLDegree(WC,{1,1,1})
	homotopyRemovalMLDegree(0,M)
	homotopyRemovalMLDegree(1,M)
	homotopyRemovalMLDegree(2,M)
	homotopyRemovalMLDegree(3,M)
	removalMLDegree M
	mlObstructionFunction M
///


--###################################
-- keys and types
--###################################


doc /// --MutableHashTable Keys
    Key
	MLDegrees
	TheVariety
	ThePoint
       	DefiningEquations
       	Hyperplanes
        Data0
        Data1
    	TheDimension
----Keys for numerical computation
        WitnessPoints
    	WitnessSets
    	SortPoints
    	WitnessRing
    	VariableGroups
    	Directory
    Headline
        keys for the new types of MutableHashTable
    Description
      Text
        These are keys for the new types of MutableHashTable in this package. 

    SeeAlso 
        RemovalMLDegree
	MLDegreeVariety
	MLDegreeWitnessSet
	MLDegreeWitnessCollection
///

doc /// --new type of MutableHashTable
    Key
        RemovalMLDegree
	MLDegreeVariety
	MLDegreeWitnessSet
	MLDegreeWitnessCollection
    Headline
        new types of MutableHashTable
    Description
      Text
        new types of MutableHashTable in this package.
	
	RemovalMLDegree is used in both numerical and symbolic computation.
	
	MLDegreeVariety is used in symbolic computation.
	
	MLDegreeWitnessCollection and MLDegreeWitnessSet are used in numerical computation.

    SeeAlso 
        newRemovalMLDegree
	newMLDegreeVariety
	newMLDegreeWitnessSet
	newMLDegreeWitnessCollection    
///



--###################################
-- Random values and vectors--
--###################################
doc /// --RandomVector
    Key
	randomVector
	(randomVector,ZZ)
        randomValue
    Headline
        Random number genreator.
    Description
      Text
        The method randomVector(n) produces a list of length n with elements given by randomValue().
	The method randomVector is used when choosing the normal vectors of the hyperplanes, Data0, and Data1. 
      Example
    	randomValue()
        randomVector(3)
      Text
	The output of randomValue() depends on the configuration of RandomCoefficients of the package. 
      Example
        ((MaximumlikelihoodObstructionFunction//options)#Configuration)
	randomValue()
        loadPackage("MaximumlikelihoodObstructionFunction",Configuration=>{"RandomCoefficients"=>ZZ/30103},Reload=>true)
	randomValue()
        loadPackage("MaximumlikelihoodObstructionFunction",Configuration=>{"RandomCoefficients"=>CC},Reload=>true)
	randomValue()
///

--###################################
-- RemovalMLDegree--
--###################################

doc /// --RemovalMLDegree
    Key
	newRemovalMLDegree
	(newRemovalMLDegree,MLDegreeVariety)
	(newRemovalMLDegree,MLDegreeVariety,List)
	(newRemovalMLDegree,MLDegreeWitnessCollection,List)
    Headline
        newRemovalMLDegree is a method that creates a RemovalMLDegree
    Usage
	newRemovalMLDegree(L)
	newRemovalMLDegree(L,P)
	newRemovalMLDegree(WC,P)
    Inputs
        L:MLDegreeVariety
          used in symbolic computation
	WC:MLDegreeWitnessCollection
	  used in numerical computation
        P:RemovalMLDegree
          specifies the point where we want to compute the Euler obstruction function.	  
    Outputs
        M:RemovalMLDegree
	  Returns a RemovalMLDegree which is used to store the results of our ML degree computations. 
    Description
      Text
        The method newRemovalMLDegree creates a new RemovalMLDegree which is ued to organize the removal ML degrees of a variety to compute the Euler obstruction at a point P. 
	
	The function newRemovalMLDegree(L) takes P to be {1,1,...,1}. To specify the coordinates of P use newRemovalMLDegree(L,P).	

      Text
	The keys of the output M  {TheVariety, ThePoint, MLDegrees, WitnessSets}. 
	
	TheVariety=>Ideal or MLDegreeWitnessCollection that defines the variety.
	
	ThePoint=>List are the coordinates of the point we compute the value of the Euler obstruction at. 
	
    	MLDegrees=>List Each computation of a ML degree appends k=>d to this list. There may be multiple elements with the same value of k. 

    	WitnessSets=>List Each computation of an ML degree appends k=>WS to this list. There may be multiple elements with the same value of k. The ith element of MLDegrees and WitnessSets correspond to one another. 

///

--###################################
-- solveRemovalMLDegree--
--###################################

doc /// --solveRemovalMLDegree
    Key
	solveRemovalMLDegree
--	(solveRemovalMLDegree,RemovalMLDegree)
	(solveRemovalMLDegree,RemovalMLDegree,ZZ)
	(solveRemovalMLDegree,ZZ,RemovalMLDegree)
	(solveRemovalMLDegree,ZZ,RemovalMLDegree,ZZ)
    Usage
  --      solveRemovalMLDegree(M),
        solveRemovalMLDegree(M,c),
        solveRemovalMLDegree(rk,M),
        solveRemovalMLDegree(rk,M,c) 
    Inputs
        rk:ZZ
          Specifies the type of removal ML degree that is computed.
        M:RemovalMLDegree
          A MutableHashTable to store the outputs of computations. 
	c:ZZ
	  Codimension of the variety of interest.
    Outputs
        d:ZZ
	  Returns the computed removal ML degree and the list at the key @ TO MLDegrees @ of M is appended with rk=>d
    Headline
        This method determines ML degrees using symbolic computation.
    Description
      Text
        solveRemovalMLDegree is a method that determines a ML degree. 

	See @ TO MaximumlikelihoodObstructionFunction @ for examples.  	
	
///


doc /// --homotopyRemovalMLDegree
    Key
	homotopyRemovalMLDegree
    	(homotopyRemovalMLDegree,RemovalMLDegree)
    	(homotopyRemovalMLDegree,ZZ,RemovalMLDegree)
    Usage
        homotopyRemovalMLDegree(M),
        homotopyRemovalMLDegree(rk,M) 
    Inputs
        rk:ZZ
          Specifies the type of removal ML degree that is computed.
        M:RemovalMLDegree
          A MutableHashTable to store the outputs of computations. 
    Outputs
        d:ZZ
	  Returns the computed removal ML degree and the list under the key @ TO MLDegrees @ of M is appended with rk=>d
    Headline
        This method determines ML degrees using homotopy continuation.
    Description
      Text
        homotopyRemovalMLDegree is a method that determines a ML degree. 

	See @ TO MaximumlikelihoodObstructionFunction @ for examples.  	
	
///


doc /// --MLDegreeWitnessCollection
    Key
	newMLDegreeWitnessCollection
    	(newMLDegreeWitnessCollection,Ideal,ZZ,String)
    	(newMLDegreeWitnessCollection,Ideal,ZZ)
    Usage
        newMLDegreeWitnessCollection(I,d,s)
        newMLDegreeWitnessCollection(I,d)
    Inputs
        I:Ideal
          Witness system (complete intersection) that contains the variety of interest as an irreducible component.
        d:ZZ
          The dimension of the the variety.
        s:String
	  A string that specifies a directory where auxillary files are saved. 
	  If s is not specified then the director temporaryFileName and mkdir will be used to create s and the directory.
    Outputs
        WC:MLDegreeWitnessCollection
	  stores the witness set information for the euler obstruction at a general point. 
    Headline
        creates a MLDegreeWitnessCollection 	
    Description
      Text
        creates a MLDegreeWitnessCollection that stores the MLDegreeWitnessSet which are used to determine the ML degree using a parameter homotopy. 

	See @ TO MaximumlikelihoodObstructionFunction @ for examples.  	
	
///



doc /// --MLDegreeWitnessCollection
    Key
	newMLDegreeVariety
    	(newMLDegreeVariety,Ideal)
    Usage
        newMLDegreeWitnessCollection(I)
    Inputs
        I:Ideal
          defines the variety of interest.
    Outputs
        L:MLDegreeVariety
	  stores information about our variety. 
    Headline
        creates a MLDegreeVariety 	
    Description
      Text
        creates a MLDegreeVariety  that stores information about the variety of interest. 
	keys include {Hyperplanes, Data0, Data1, DefiningEquations}
	
	Hyperplanes=>List of List that define normal vectors of hyperplanes that go through the point P.

    	Data0=>List of random numbers that correspond to the data

    	Data1=>List of random numbers that correspond to the data
	
	DefiningEquations=>Ideal defining the variety. 

	See @ TO MaximumlikelihoodObstructionFunction @ for examples.  	

	
///



doc /// --MLDegreeWitnessCollection
    Key
	newMLDegreeWitnessSet
	(newMLDegreeWitnessSet,MLDegreeWitnessCollection)
	(newMLDegreeWitnessSet,ZZ,MLDegreeWitnessCollection)    	
    Usage
	newMLDegreeWitnessSet(WC)
	newMLDegreeWitnessSet(rk,WC)    	
    Inputs
        rk:ZZ
          Specifies the removal ML degree
        WC:MLDegreeWitnessCollection
	  the witness set information will be stored here. 
    Outputs
        m:ZZ
	  The rk-th removal ML degree. 
    Headline
        creates a MLDegreeWitnessSet 
    Description
      Text
        creates a MLDegreeWitnessSet WS which is used as the start points of a homotopy.
	This witness set is appended to the List in the WitnessSets key of WC. 
	The ML degree m that correspond to this witness set is appended to the List in the MLDegrees keys of WC. 

	See @ TO MaximumlikelihoodObstructionFunction @ for examples.  	
	
///

doc /// --MLDegreeWitnessCollection
    Key
	mlObstructionFunction
	(mlObstructionFunction,RemovalMLDegree)
 	(mlObstructionFunction,RemovalMLDegree,ZZ)
    Usage
	mlObstructionFunction(M)
	newMLDegreeWitnessSet(M,d)    	
    Inputs
        M:RemovalMLDegree
          A MutableHashTable that is used to store the outputs of the ML degree computations. 
        d:ZZ
	  the dimension of the variety.
    Outputs
        mlo:ZZ
	  The Euler obstruction at the point P=M#ThePoint.
    Headline
        compute the Euler obstruction function at a point.
    Description
      Text
        Takes the alternating sum of computed ML degrees. 

      Text
        The method mlObstructionFunction(M,d) returns the alternating sum of ML degrees \sum_k (-1)^k r_k(1,X) times (-1)^d. It is assumed the ML degrees are correct and complete otherwise incorrect answers will be returned. 
	
	See @ TO MaximumlikelihoodObstructionFunction @ for examples.  	
	
///

doc /// --MLDegreeWitnessCollection
    Key
	removalMLDegree
	(removalMLDegree,RemovalMLDegree)
    Usage
	removalMLDegree(M)
    Inputs
        M:RemovalMLDegree
          A MutableHashTable that is used to store the outputs of the ML degree computations. 
    Outputs
        mld:List
	  List of computed ML degrees.
    Headline
        returns the computed ML degrees.
    Description
      Text
        Returns the values of the ML degrees that have been computed. 
	
	See @ TO MaximumlikelihoodObstructionFunction @ for examples.  	
	
///

doc /// --MLDegreeWitnessCollection
    Key
	reclassifyWitnessPoints
	(reclassifyWitnessPoints,MutableHashTable,FunctionClosure)
	(reclassifyWitnessPoints,MutableHashTable,FunctionClosure,ZZ)
	(reclassifyWitnessPoints,MutableHashTable,RR)
	(reclassifyWitnessPoints,MutableHashTable,RR,ZZ)
    Usage
	reclassifyWitnessPoints(WCM,spf)
	reclassifyWitnessPoints(WCM,spf,p)
	reclassifyWitnessPoints(WCM,epsilon)
	reclassifyWitnessPoints(WCM,epsilon,p)
    Inputs
        WCM:MutableHashTable
          This can be a MLDegreeWitnessCollection or a RemovalMLDegree 
        spf:FunctionClosure
          A function used to determine if the point is in the coordinate hyperplanes. 
        epsilon:RR
          The tolerance used to determine if a point is in a coordinate hyerplane. 
        p:ZZ
          position of the witness set whose points are going to be reclassified. 
    Outputs
        mld:List
	  Computed ML degree.
    Headline
        returns the computed ML degree 	
    Description
      Text
        Returns the the computed ML degree after reclassifying the witness points according to the tolerance epsilon or function spf.
      Text
        The ML degree counts critical points outside the coordinate hyerplanes and singular locus. A default tolerance is set to determine if point is in the coordinate hyperplanes and singular locus. This function reclassifies the computed witness points according to a new tolerance.        	
      Example
        loadPackage("MaximumlikelihoodObstructionFunction",Configuration=>{"RandomCoefficients"=>CC},Reload=>true)
	s = temporaryFileName() | "/"                            
	mkdir s
	R=CC[p1,p2,p3]
	f=sub(p1^2-p2^2*p3,{p1=>p1-1,    p2=>p2-1,    p3=>p3-1})	
	I=ideal(f);    	d=2;
	WC=newMLDegreeWitnessCollection(I,d,s)
	newMLDegreeWitnessSet(WC)
	reclassifyWitnessPoints(WC,10.0)
	reclassifyWitnessPoints(WC,1e-10)	
	M=newRemovalMLDegree(WC,{1,1,1})
	homotopyRemovalMLDegree M
	removalMLDegree M
        reclassifyWitnessPoints(M,1e-8)
        reclassifyWitnessPoints(M,1e-200)
	reclassifyWitnessPoints(WC,10.0)
        reclassifyWitnessPoints(M,1e-8)
	
///


end








