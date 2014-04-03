function [labels, energy, time] = alphaBetaSwapGridPotts(unary, vertC, ...
    horC, metric, options)
% This function finds errors in message
% INPUT:
%    unary: m-by-n-by-k array of double, unary potentials
%    vertC: (n - 1)-by-m array of double, coefficients for vertical edges
%    horC: n-by-(m - 1) array of double, coefficients for horizontal edges
%    metric: k-by-k array of double, distance between labels
%    options (optional): structure, may contain fields:
%        'maxIter': integer, max number of iterations in algorithm
%            (default: 500)
%        'display': logical, display statistics on each iteration
%            (default: false)
%        'numStart': integer, number of launches for different start points
%            (default: 1)
%        'randOrder': logical, use random order of alpha-beta pairs
%            (default: false)
%
% OUTPUT:
%    labels: n-by-m array of double ?, result labels
%    energy: niter-by-1 array of double, energy on each iteration
%    time: niter-by-1 array of double, time in sec from start to each
%        iteration

    MAX_ITER = 500;
    DISPLAY = false;
    NUM_START = 1;
    RAND_ORDER = false;
    if nargin > 4
        if isfield(options, 'maxIter')
            MAX_ITER = options.maxIter;
        end
        if isfield(options, 'display')
            DISPLAY = options.display;
        end
        if isfield(options, 'numStart')
            NUM_START = options.numStart;
        end
        if isfield(options, 'randOrder')
            RAND_ORDER = options.randOrder;
        end
    end
    
end