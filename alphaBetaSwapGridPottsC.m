function [best_labels, best_energy, best_time] = alphaBetaSwapGridPottsC(unary_pot, ...
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

    EPS = 1e-1;
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
    
    [N, M, K] = size(unary_pot);
    best_energy = Inf;
    
    for q = 1 : NUM_START
        % initialization
        labels = randi(K, [N, M]);
        time = zeros(1, MAX_ITER);
        energy = zeros(1, MAX_ITER);
        cur_energy = Inf;
        [alphas, betas] = find(triu(1 - eye(K)));
        list = [alphas, betas];

        for t = 1 : MAX_ITER
            tic;
            if mod(t - 1, size(list, 1)) == 0
                if RAND_ORDER
                    list = list(randperm(length(list)), :);
                end
            end
            s = mod(t - 1, size(list, 1)) + 1;
            alpha = list(s, 1);
            beta = list(s, 2);
            
            old_energy = cur_energy;
            cur_energy = get_energyC(labels, unary_pot, vertC, horC, metric);
            energy(t) = cur_energy;
            if DISPLAY
                fprintf('Iteration %d: alpha = %d, beta = %d, E = %f\n', ...
                    t, alpha, beta, cur_energy);
            end
            if (abs(cur_energy - old_energy) < EPS)
                break;
            end
            
            % find alpha and beta indices, (i, j) -> i + (j - 1) * N
            inds = find((labels == alpha) | (labels == beta));
            rev_inds = zeros(1, N * M);
            rev_inds(inds) = 1 : length(inds);

            mask = (labels == alpha) | (labels == beta);
            inds_ver = find([false(1, M); (mask(1 : (end - 1), :) == ...
                mask(2 : end, :))] & (mask == true));
            inds_hor = find([false(N, 1), (mask(:, 1 : (end - 1)) == ...
                mask(:, 2 : end))] & (mask == true));

            inds_edge = find((mask == true) & ...
                ([false(1, M); (~mask(1 : (end - 1), :) == mask(2 : end, :))] | ...
                [(~mask(:, 1 : (end - 1)) == mask(:, 2 : end)), false(N, 1)] | ...
                [(~mask(2 : end, :) == mask(1 : (end - 1), :)); false(1, M)] | ...
                [false(N, 1), (~mask(:, 2 : end) == mask(:, 1 : (end - 1)))]));

            % build graph
            % source, sink
            terminal_weights = [unary_pot(inds + N * M * (alpha - 1)), ...
                unary_pot(inds + N * M * (beta - 1))];
            terminal_weights(rev_inds(inds_edge), :) = ...
                terminal_weights(rev_inds(inds_edge), :) + ...
                    edge_additions(inds_edge, alpha, beta, labels, ...
                    vertC, horC, metric);

            % from, to, capacity, reverse_capacity
            edge_weights = zeros(length([inds_ver; inds_hor]), 4);
            edge_weights(1 : length(inds_ver), :) = ...
                get_edges(rev_inds, inds_ver, ...
                inds_ver - 1, ...
                alpha, beta, vertC, horC, metric, true);
            edge_weights((length(inds_ver) + 1) : end, :) = ...
                get_edges(rev_inds, inds_hor, ...
                inds_hor - N, ...
                alpha, beta, vertC, horC, metric, false);

            % find labels of min cut
            [~, bin_labels] = graphCutMex(terminal_weights, edge_weights);
            % respectively change labels to alpha / beta
            bin_labels(bin_labels == 1) = alpha;
            bin_labels(bin_labels == 0) = beta;
            labels(inds) = bin_labels;
            
            time(t) = toc;
        end
        time = cumsum(time);
        if (cur_energy < best_energy(end))
            best_energy = energy;
            best_labels = labels;
            best_time = time;
        end
    end
end

function res = get_edges(rev_inds, ind1, ind2, alpha, beta, vertC, horC, metric, ver)
    rev_inds = rev_inds(:);
    ind1 = ind1(:);
    ind2 = ind2(:);
    if ver
        % (i, j) -> i + (j - 1) * size_1
        [N, ~, ~] = size(horC);
        i2 = mod(ind2, N);
        j2 = fix(ind2 / N) + 1;
        j2(i2 == 0) = j2(i2 == 0) - 1;
        i2(i2 == 0) = N;
        lin_inds = i2 + (j2 - 1) * size(vertC, 1);
        pp = repmat(vertC(lin_inds), 1, 2) .* ...
            repmat([metric(alpha, beta) - metric(alpha, alpha), ...
            metric(beta, alpha) - metric(beta, beta)], length(lin_inds), 1);
    else
        pp = repmat(horC(ind2), 1, 2) .* ...
            repmat([metric(alpha, beta) - metric(alpha, alpha), ...
            metric(beta, alpha) - metric(beta, beta)], length(ind2), 1);
    end
    res = [rev_inds(ind1), rev_inds(ind2), pp(:, 1), pp(:, 2)];
end

function [res, ngbr] = edge_additions(inds, alpha, beta, labels, vertC, horC, metric)
    [N, ~] = size(horC);
    [~, M] = size(vertC);
    [K, ~] = size(metric);
    inds = inds(:);
    add = zeros(length(inds), 4, 2);
    ngbr = zeros(length(inds), 4);
    ngbr(:, 1) = inds + 1;
    ngbr(:, 2) = inds - 1;
    ngbr(:, 3) = inds + N;
    ngbr(:, 4) = inds - N;
    mask = (ngbr > 0) & (ngbr <= N * M);
    mask(:, 1) = mask(:, 1) & (mod(inds, N) ~= 0);
    mask(:, 2) = mask(:, 2) & (mod(inds - 1, N) ~= 0);
    mask(mask) = (labels(ngbr(mask)) ~= alpha) & (labels(ngbr(mask)) ~= beta);
    % n1
    i = mod(inds, N);
    j = fix(inds / N) + 1;
    j(i == 0) = j(i == 0) - 1;
    i(i == 0) = N;
    lin_inds = i + (j - 1) * size(vertC, 1);
    add(mask(:, 1), 1, 1) = vertC(lin_inds(mask(:, 1))) .* ...
        metric(labels(ngbr(mask(:, 1), 1)) + K * (alpha - 1));
    add(mask(:, 1), 1, 2) = vertC(lin_inds(mask(:, 1))) .* ...
        metric(labels(ngbr(mask(:, 1), 1)) + K * (beta - 1));
    % n2
    i = mod(ngbr(:, 2), N);
    j = fix(ngbr(:, 2) / N) + 1;
    j(i == 0) = j(i == 0) - 1;
    i(i == 0) = N;
    lin_inds = i + (j - 1) * size(vertC, 1);
    add(mask(:, 2), 2, 1) = vertC(lin_inds(mask(:, 2))) .* ...
        metric(alpha + K * (labels(ngbr(mask(:, 2), 2)) - 1));
    add(mask(:, 2), 2, 2) = vertC(lin_inds(mask(:, 2))) .* ...
        metric(beta + K * (labels(ngbr(mask(:, 2), 2)) - 1));
    % in horC the same absolute index as in big matrix
    % n3
    add(mask(:, 3), 3, 1) = horC(inds(mask(:, 3))) .* ...
        metric(labels(ngbr(mask(:, 3), 3)) + K * (alpha - 1));
    add(mask(:, 3), 3, 2) = horC(inds(mask(:, 3))) .* ...
        metric(labels(ngbr(mask(:, 3), 3)) + K * (beta - 1));
    % n4
    add(mask(:, 4), 4, 1) = horC(ngbr(mask(:, 4), 4)) .* ...
        metric(alpha + K * (labels(ngbr(mask(:, 4), 4)) - 1));
    add(mask(:, 4), 4, 2) = horC(ngbr(mask(:, 4), 4)) .* ...
        metric(beta + K * (labels(ngbr(mask(:, 4), 4)) - 1));
    res = squeeze(sum(add, 2));
    ngbr = ngbr(mask);
end