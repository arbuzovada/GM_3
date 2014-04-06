function [labels, energy, time] = alphaBetaSwapGridPottsC(unary_pot, ...
    vertC, horC, metric, options)
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
    if nargin > 2
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
    
    [N, M, K] = size(unary_pot);
    
    % initialization
    labels = unary_pot(:, :, 2) + 1;%randi(K, [N, M]);
    time = zeros(1, MAX_ITER);
    mask_ver = zeros(N, M);
    mask_hor = zeros(N, M);
                    
    for t = 1 : MAX_ITER
        if DISPLAY
            energy = get_energyC(labels, unary_pot, vertC, horC, metric);
            fprintf('Iteration %d: E = %f\nNumber of 2nd pic: %d\n', t, ...
                energy, sum(sum(labels == 2)));
        end
        tic;
        for alpha = 1 : K
            for beta = (alpha + 1) : K
                % find alpha and beta indices, (i, j) -> i + (j - 1) * N
                inds = find((labels == alpha) | (labels == beta)); % col
                % build graph
                % source, sink
                terminal_weights = [unary_pot(inds + N * M * (alpha - 1)), ...
                    unary_pot(inds + N * M * (beta - 1))];
                
                % from, to, capacity, reverse_capacity
                mask = (labels == alpha) | (labels == beta);
                mask_ver(2 : end, :) = mask(1 : (end - 1), :);
                mask_hor(:, 2 : end) = mask(:, 1 : (end - 1));
                inds_ver = find(mask == mask_ver);
%                 max(inds_ver)
                inds_hor = find(mask == mask_hor);
%                 max(inds_hor)
                edge_weights = zeros(length([inds_ver; inds_hor]), 4);
                for i = 1 : length(inds_ver)
%                     if isempty(find(inds == inds_ver(i), 1)) || isempty(find(inds == inds_ver(i) - 1, 1))
%                         inds_ver(i)
%                         inds_ver(i) - 1
%                     end
                    edge_weights(i, :) = ...
                        get_pair_pot(find(inds == inds_ver(i), 1), ...
                        find(inds == inds_ver(i) - 1, 1), ...
                        alpha, beta, vertC, horC, metric, true);
                    if any(edge_weights(i, :) < 0)
                        edge_weights(i, :)
                        return
                    end
%                     edge_weights(i, :) = [find(inds == inds_ver(i), 1), ...
%                         find(inds == inds_ver(i) - 1, 1), ...
%                         pair_pot(alpha, beta) - pair_pot(beta, beta), ...
%                         pair_pot(beta, alpha) - pair_pot(alpha, alpha)];
                end
                fprintf('hor...\n');
                for i = 1 : length(inds_hor)
                    edge_weights(i + length(inds_ver), :) = ...
                        get_pair_pot(find(inds == inds_hor(i), 1), ...
                        find(inds == inds_hor(i) - N, 1), ...
                        alpha, beta, vertC, horC, metric, false);
                end
                
                % find labels of min cut
                fprintf('  Start graph cut...')
                [~, bin_labels] = graphCutMex(terminal_weights, edge_weights);
                fprintf(' Done.\n')
                % respectively change labels to alpha / beta
                bin_labels(bin_labels == 1) = alpha;
                bin_labels(bin_labels == 0) = beta;
                labels(inds) = bin_labels;
            end
        end
        time(t) = toc;
    end
    time = cumsum(time);
end

function res = get_pair_pot(ind1, ind2, alpha, beta, vertC, horC, metric, ver)
    [N, ~, ~] = size(horC);
    i2 = mod(ind2, N);
    j2 = fix(ind2 / N) + 1;
    if i2 == 0
        i2 = N;
        j2 = j2 - 1;
    end
    if ver
        pp = vertC(i2, j2) * [metric(alpha, beta), metric(beta, beta), ...
            metric(beta, alpha), metric(alpha, alpha)];
    else
        pp = horC(i2, j2) * [metric(alpha, beta), metric(beta, beta), ...
            metric(beta, alpha), metric(alpha, alpha)];
    end
    res = [ind1, ind2, pp(1) - pp(2), pp(3) - pp(4)];
end