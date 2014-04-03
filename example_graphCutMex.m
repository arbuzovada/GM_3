nNodes=4;
%source,sink
terminalWeights=[
    16,0;
    13,0;
    0,20;
    0,4
];

%From,To,Capacity,Rev_Capacity
edgeWeights=[
    1,2,10,4;
    1,3,12,0;
    2,3,0,9;
    2,4,14,0;
    3,4,0,7
    ];

[cut, labels] = graphCutMex(terminalWeights,edgeWeights);
