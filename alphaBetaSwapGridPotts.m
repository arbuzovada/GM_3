function [best_labels, best_energy, best_time] = alphaBetaSwapGridPotts(unary_pot, ...
    intensity, options)
% This function finds errors in message
% INPUT:
%    unary: N-by-M-by-K array of double, unary potentials
%    intensity: N-by-M-by-K array of double, intensity of images
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

    EPS = 1;
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
            cur_energy = get_energy(labels, unary_pot, intensity);
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
                    edge_additions(inds_edge, alpha, beta, labels, intensity);

            % from, to, capacity, reverse_capacity
            edge_weights = zeros(length([inds_ver; inds_hor]), 4);
            edge_weights(1 : length(inds_ver), :) = ...
                get_edges(rev_inds, inds_ver, ...
                inds_ver - 1, ...
                alpha, beta, intensity);
            edge_weights((length(inds_ver) + 1) : end, :) = ...
                get_edges(rev_inds, inds_hor, ...
                inds_hor - N, ...
                alpha, beta, intensity);

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

function res = get_edges(rev_inds, inds1, inds2, alpha, beta, intensity)
    rev_inds = rev_inds(:);    
    inds1 = inds1(:);
    inds2 = inds2(:);
    pp = zeros(length(inds1), 2);
    pp(:, 1) = get_pp(inds1, inds2, alpha, beta, intensity);
    pp(:, 2) = get_pp(inds1, inds2, beta, alpha, intensity);
    res = [rev_inds(inds1), rev_inds(inds2), pp(:, 1), pp(:, 2)];
end

function res = get_pp(inds1, inds2, labels1, labels2, intensity)
    [N, M, ~] = size(intensity);
    res = ...
        abs(intensity(inds1 + N * M * (labels1 - 1)) - ...
            intensity(inds1 + N * M * (labels2 - 1))) + ...
        abs(intensity(inds2 + N * M * (labels1 - 1)) - ...
            intensity(inds2 + N * M * (labels2 - 1)));
end

function res = edge_additions(inds, alpha, beta, labels, intensity)
% res: len(inds)-by-2
    [N, M, ~] = size(intensity);
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
    % i, j - inds
    % i + 1, j - ngbr
    add(mask(:, 1), 1, 1) = get_pp(inds(mask(:, 1)), ngbr(mask(:, 1), 1), ...
        alpha, labels(ngbr(mask(:, 1), 1)), intensity);
    add(mask(:, 1), 1, 2) = get_pp(inds(mask(:, 1)), ngbr(mask(:, 1), 1), ...
        beta, labels(ngbr(mask(:, 1), 1)), intensity);
    % n2
    % i, j - ngbr
    % i + 1, j - inds
    add(mask(:, 2), 2, 1) = get_pp(inds(mask(:, 2)), ngbr(mask(:, 2), 2), ...
        alpha, labels(ngbr(mask(:, 2), 2)), intensity);
    add(mask(:, 2), 2, 2) = get_pp(inds(mask(:, 2)), ngbr(mask(:, 2), 2), ...
        beta, labels(ngbr(mask(:, 2), 2)), intensity);
    % in horC the same absolute index as in big matrix
    % n3
    add(mask(:, 3), 3, 1) = get_pp(inds(mask(:, 3)), ngbr(mask(:, 3), 3), ...
        alpha, labels(ngbr(mask(:, 3), 3)), intensity);
    add(mask(:, 3), 3, 2) = get_pp(inds(mask(:, 3)), ngbr(mask(:, 3), 3), ...
        beta, labels(ngbr(mask(:, 3), 3)), intensity);
    % n4
    add(mask(:, 4), 4, 1) = get_pp(inds(mask(:, 4)), ngbr(mask(:, 4), 4), ...
        alpha, labels(ngbr(mask(:, 4), 4)), intensity);
    add(mask(:, 4), 4, 2) = get_pp(inds(mask(:, 4)), ngbr(mask(:, 4), 4), ...
        beta, labels(ngbr(mask(:, 4), 4)), intensity);
    res = squeeze(sum(add, 2)); %[sum(add(:, :, 1), 2), sum(add(:, :, 2), 2)]; 
end