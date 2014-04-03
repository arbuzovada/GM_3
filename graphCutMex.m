% Matlab interface to Vladimir Kolmogorov implementation of min-cut algorithm downloadable from:
% http://www.cs.ucl.ac.uk/staff/V.Kolmogorov/software.html
% 
% The algorithm is described in:
% Yuri Boykov and Vladimir Kolmogorov, 'An Experimental Comparison of Min-Cut/Max-Flow Algorithms for
% Energy Minimization in Vision', IEEE Transactions on Pattern Analysis and Machine Intelligence, vol.
% 26, no. 9, pp. 1124-1137, Sept. 2004.
% 
% 
% Usage:
% [cut] = graphCutMex(termWeights, edgeWeights);
% [cut, labels] = graphCutMex(termWeights, edgeWeights);
% 
% 
% Inputs:
% termWeights	-	the edges connecting the source and the sink with regular nodes (array of type double, size : [numNodes, 2])
% 				termWeights(i, 1) is the weight of the edge connecting the source with node #i
% 				termWeights(i, 2) is the weight of the edge connecting node #i with the sink
% 				numNodes is determined from the size of termWeights.
% edgeWeights	-	the edges connecting regular nodes with each other (array of type double, array size [numEdges, 4])
% 				edgeWeights(i, 3) connects node #edgeWeights(i, 1) to node #edgeWeights(i, 2)
% 				edgeWeights(i, 4) connects node #edgeWeights(i, 2) to node #edgeWeights(i, 1)
% Outputs:
% cut           -	the minimum cut value (type double)
% labels		-	a vector of length numNodes, where labels(i) is 0 or 1 if node #i belongs to S (source) or T (sink) respectively.
% 
% 
% To build the code in matlab choose reasonable compiler and run buildGraphCutMex.m
% 
% by Anton Osokin, firstname.lastname@gmail.com, Spring 2011