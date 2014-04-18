function [labels, energy, time] = alphaExpansionGridPotts(unary, vertC, horC, metric, options)
% ����
% unary � ������� ����������, ������ ���� double ������� N x M x K, ��� N � ������ �������, M � ������ �������, K � ���������� �����.
% vertC � ������������  c_{pq}, ��������������� ������������ ������, ������ ���� double ������� (N � 1) x M;
% horC � ������������  c_{pq}, ��������������� �������������� ������, ������ ���� double ������� N x (M � 1);
% metric � ���������� ����� ������� �������� ����������, ������ ���� double ������� K x K;
% options � (�������������� ��������) ����� �������������� ����������, ��������� � ������
%   'maxIter' � ����������� ���������� ����� �������� ��������� (�� ��������� = 500);
%   'display' � �������� ���� logical: ���� true, �� ��� ������ ������� ��������� ������� ����� ����� �������� �� ����� ����� ��������, ������ �������������� �����, ������� �������� �������;
%   'numStart' � ���������� �������� �� ������ ��������� �����������;
%   'randOrder' � �������� ���� logical: ���� true, �� ��� ������ ������� ������������ ��������� ������� ����� ? � ?;
% �����
% labels � ��������, ���������� ���������� ��������, ������ ���� double ������� N x M;
% energy � �������� ������� �� ������ ��������, ������ ���� double �����, ������ ���������� �������� ���������;
% time � �����, ���������� � ������ ������ ��������� �� ������ ��������, ������ ���� double �����, ������ ���������� �������� ���������;

    max_iter = 10;
    display = false;
    num_start = 1;
    rand_order = false;
    
    if nargin > 4
        if isfield(options, 'maxIter')
            max_iter = options.maxIter;
        end
        if isfield(options, 'display')
            display = options.display;
        end
        if isfield(options, 'numStart')
            num_start = options.numStart;
        end
        if isfield(options, 'randOrder')
            rand_order = options.randOrder;
        end
    end
    
    
    [N, M, K] = size(unary);
    labels = randi(K, N, M);
    for start = 1:num_start
        en_num = 0;
        for iter = 1:max_iter
            iter
            if rand_order
                order = randperm(K);
            else
                order = 1:K;
            end
            for alpha = order
                if display
                    fprintf('Start � %d, iteration %d, expanding label = %d, processing...\n', start, en_num, alpha);
                end
                if display
                    fprintf('Energy = %f\n', get_energyC(labels, unary, vertC, horC, metric));
                end
                tic;
                en_num = en_num + 1;
                mask = (repmat(labels, [1, 1, K]) == repmat(reshape(1:K, 1, 1, K), N, M));
                tmp = sum(unary .* mask, 3);
                terminal_weights = [reshape(unary(:, :, alpha), N * M, 1), reshape(tmp, N * M, 1)];
                clear edge_weights tmp

                rep = [labels(:, 1), reshape(repmat(labels(:, 2:(end - 1)), 2, 1), N, ...
                    2 * (M - 2)), labels(:, end)];
                hor = (reshape(rep', 2, N * (M - 1)))';
                nums = hor(:, 1) + (hor(:, 2) - 1) * K;
                line = reshape(metric, K * K, 1);
                res = (reshape(line(nums), M - 1, N))';
                tmp1 = labels(:, 1:(end - 1));
                tmp2 = labels(:, 2:end);
                hor_msg = horC .* (reshape(metric(tmp1, alpha), N, M - 1) + reshape(metric(alpha, tmp2), N, M - 1) - res - metric(alpha, alpha));
                inds = 1:(N * (M - 1));
                terminal_weights(inds, 2) = terminal_weights(inds, 2) + ...
                    reshape(horC .* res, [], 1);
                terminal_weights(inds, 1) = terminal_weights(inds, 1) + ...
                    reshape(horC, [], 1) .* metric(alpha, alpha) - ...
                    reshape(horC, [], 1) .* metric(tmp1, alpha) + ...
                    reshape(horC .* res, [], 1);
                terminal_weights(inds + N, 1) = terminal_weights(inds + N, 1) + ...
                    reshape(horC, [], 1) .* metric(tmp1, alpha) - reshape(horC .* res, [], 1);


                rep = [labels(1, :); (reshape((repmat(labels(2:(end - 1), :), 1, 2))', M, 2 * (N - 2)))'; labels(end, :)];
                vert = (reshape(rep, 2, []))';
                nums = vert(:, 1) + (vert(:, 2) - 1) * K;
                res = reshape(line(nums), N - 1, M);

                tmp1 = labels(1:(end - 1), :);
                tmp2 = labels(2:end, :);
                vert_msg = vertC .* (reshape(metric(tmp1, alpha), N - 1, M) + reshape(metric(alpha, tmp2), N - 1, M) - res - metric(alpha, alpha));
                inds = (1:(N * M))';
                inds = inds(rem(inds, N) ~= 0);
                terminal_weights(inds, 2) = terminal_weights(inds, 2) + reshape(vertC .* res, [], 1);
                terminal_weights(inds, 1) = terminal_weights(inds, 1) + ...
                    reshape(vertC, [], 1) .* metric(alpha, alpha) - ...
                    reshape(vertC, [], 1) .* metric(tmp1, alpha) + ...
                    reshape(vertC .* res, [], 1);
                terminal_weights(inds + 1, 1) = terminal_weights(inds + 1, 1) + ...
                    reshape(vertC, [], 1) .* metric(tmp1, alpha) - ...
                    reshape(vertC .* res, [], 1);

                from = (1:(N * M))';
                from = from(rem(from, N) ~= 0);
                to = from + 1;
                edge_weights = [from, to, zeros(size(from)), reshape(vert_msg, [], 1)];
                from = (1:(N * (M - 1)))';
                to = from + N;
                edge_weights = [edge_weights; from, to, zeros(size(from)), reshape(hor_msg, [], 1)];
%                edge_weights = [edge_weights; from, to, zeros(size(from)), zeros(size(from))];

                delta = min(terminal_weights, [], 2);
                terminal_weights = terminal_weights - repmat(delta, 1, 2);
                delta_const = sum(delta);
                [cut, bin_labels] = graphCutMex(terminal_weights,edge_weights);
                labels(logical(reshape(bin_labels, N, M))) = alpha;
                energy_cur(en_num) = cut + delta_const;
                if display
                    fprintf('Energy = %f\n', energy_cur(end));
                end
                time_cur(en_num) = toc;
            end
            if (all(energy_cur((end - K + 1) : end) - energy_cur(end) < 1e-0))
                break;
            end
        end
        if start == 1
            energy = energy_cur;
            time = time_cur;
        elseif energy_cur(end) < energy(end)
            energy = energy_cur;
            time = time_cur;
        end
    end
    
    time = cumsum(time);
end