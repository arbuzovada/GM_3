function [labels, energy, time] = alphaBetaSwapGridPotts(unary, vertC, ...
    horC, metric, options)
% This function finds errors in message
% INPUT:
%    unary: N-by-M-by-K array of double, unary potentials
%    vertC: (N - 1)-by-M array of double, coefficients for vertical edges
%    horC: N-by-(M - 1) array of double, coefficients for horizontal edges
%    metric: K-by-K array of double, distance between labels
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
%    labels: N-by-M array of double ?, result labels
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
    
    [N, M, K] = size(unary);
    
    % initialization
    labels = randi(K, [N, M]);
    time = zeros(1, MAX_ITER);
    
    for t = 1 : MAX_ITER
        if DISPLAY
%             energy = get_energy(labels, unary, pair_pot);
            energy = get_energy(labels, unary, vertC, horC, metric);
            fprintf('Iteration %d: E = %f\n', t, energy);
        end
        tic;
        for alpha = 1 : K
            for beta = (alpha + 1) : K
                % find alpha and beta indices, (i, j) -> i + (j - 1) * N
                inds = find((labels == alpha) | (labels == beta)); % col
                % build graph
                % source, sink
                terminal_weights = [unary(inds + N * M * (alpha - 1)), ...
                    unary(inds + N * M * (beta - 1))];
                
                % from, to, capacity, reverse_capacity
                mask = zeros(N + 1, M + 1);
                mask_ver = zeros(N + 1, M + 1);
                mask_hor = zeros(N + 1, M + 1);
                mask(1 : (end - 1), 1 : (end - 1)) = ...
                    (labels == alpha) | (labels == beta);
                mask_ver(2 : end, 1 : (end - 1)) = mask;
                mask_hor(1 : (end - 1), 2 : end) = mask;
                inds_ver = find(mask == mask_ver);
                inds_hor = find(mask == mask_hor);
                edge_weights = zeros(length([inds_ver; inds_hor]), 4);
                for i = 1 : length(inds_ver)
                    edge_weights(i, :) = [find(inds == inds_ver(i), 1), ...
                        find(inds == inds_ver(i) - 1, 1), ...
                        pair_pot(alpha, beta) - pair_pot(beta, beta), ...
                        pair_pot(beta, alpha) - pair_pot(alpha, alpha)];
                end
                for i = 1 : length(inds_hor)
                    edge_weights(i + length(inds_ver), :) = ...
                        [find(inds == inds_hor(i), 1), ...
                        find(inds == inds_hor(i) - N, 1), ...
                        pair_pot(alpha, beta) - pair_pot(beta, beta), ...
                        pair_pot(beta, alpha) - pair_pot(alpha, alpha)];
                end
                
                % find labels of min cut
                [~, bin_labels] = graphCutMex(terminal_weights, edge_weights);
                
                % respectively change labels to alpha/beta
                bin_labels(bin_labels == 1) = alpha;
                bin_labels(bin_labels == 0) = beta;
                labels(inds) = bin_labels;
            end
        end
        time(t) = toc;
    end
    time = cumsum(time);
end