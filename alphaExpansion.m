function [labels, energy, time] = alphaExpansion(unary, L, options)
% ����
% unary � ������� ����������, ������ ���� double ������� N x M x K, ��� N � ������ �������, M � ������ �������, K � ���������� �����.
% L - ������� �������� (NxMxK)
% options � (�������������� ��������) ����� �������������� ����������, ��������� � ������
%   'maxIter' � ����������� ���������� ����� �������� ��������� (�� ��������� = 500);
%   'display' � �������� ���� logical: ���� true, �� ��� ������ ������� ��������� ������� ����� ����� �������� �� ����� ����� ��������, ������ �������������� �����, ������� �������� �������;
%   'numStart' � ���������� �������� �� ������ ��������� �����������;
%   'randOrder' � �������� ���� logical: ���� true, �� ��� ������ ������� ������������ ��������� ������� ����� ? � ?;
% �����
% labels � ��������, ���������� ���������� ��������, ������ ���� double ������� N x M;
% energy � �������� ������� �� ������ ��������, ������ ���� double �����, ������ ���������� �������� ���������;
% time � �����, ���������� � ������ ������ ��������� �� ������ ��������, ������ ���� double �����, ������ ���������� �������� ���������;

    max_iter = 500;
    display = false;
    num_start = 1;
    rand_order = false;
    
    if nargin > 2
        if isfield(options, 'maxIter')
            max_iter = options.maxIter;
        end
        if isfield(options, 'display')
            fprintf('tu');
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
                %display
                %options.display
                if display
                    fprintf('Start � %d, iteration %d, expanding label = %d, processing...\n', start, en_num, alpha);
                end
                tic;
                en_num = en_num + 1;
                mask = (repmat(labels, [1, 1, K]) == repmat(reshape(1:K, 1, 1, K), N, M));
                tmp = sum(unary .* mask, 3);
                terminal_weights = [reshape(unary(:, :, alpha), N * M, 1), reshape(tmp, N * M, 1)];
                clear edge_weights tmp
                
                mask = (repmat(labels(1:(end - 1), :), [1, 1, K]) == repmat(reshape(1:K, 1, 1, K), N - 1, M));
                tmp1 = sum(L(1:(end - 1), :, :) .* mask, 3);
                mask = (repmat(labels(2:end, :), [1, 1, K]) == repmat(reshape(1:K, 1, 1, K), N - 1, M));
                tmp2 = sum(L(1:(end - 1), :, :) .* mask, 3);
                mask = (repmat(labels(1:(end - 1), :), [1, 1, K]) == repmat(reshape(1:K, 1, 1, K), N - 1, M));
                tmp3 = sum(L(2:end, :, :) .* mask, 3);
                mask = (repmat(labels(2:end, :), [1, 1, K]) == repmat(reshape(1:K, 1, 1, K), N - 1, M));
                tmp4 = sum(L(2:end, :, :) .* mask, 3);
                c = abs(tmp1 - L(1:(end - 1), :, alpha)) + abs(tmp3 - L(2:end, :, alpha));
                d = abs(L(1:(end - 1), :, alpha) - tmp2) + abs(L(2:end, :, alpha) - tmp4);
                a = abs(tmp1 - tmp2) + abs(tmp3 - tmp4);
                vert = c + d - a;
                from = (1:(N * M))';
                from = from(rem(from, N) ~= 0);
                to = from + 1;               
                terminal_weights(from, 2) = terminal_weights(from, 2) + reshape(a, [], 1);
                terminal_weights(from, 1) = terminal_weights(from, 1) + reshape(-c + a, [] ,1);
                terminal_weights(to, 1) = terminal_weights(to, 1) + reshape(c - a, [], 1);
                edge_weights = [from, to, zeros(size(from)), reshape(vert, [], 1)];
                
                mask = (repmat(labels(:, 1:(end - 1)), [1, 1, K]) == repmat(reshape(1:K, 1, 1, K), N, M - 1));
                tmp1 = sum(L(:, 1:(end - 1), :) .* mask, 3);
                mask = (repmat(labels(:, 2:end), [1, 1, K]) == repmat(reshape(1:K, 1, 1, K), N, M - 1));
                tmp2 = sum(L(:, 1:(end - 1), :) .* mask, 3);
                mask = (repmat(labels(:, 1:(end - 1)), [1, 1, K]) == repmat(reshape(1:K, 1, 1, K), N, M - 1));
                tmp3 = sum(L(:, 2:end, :) .* mask, 3);
                mask = (repmat(labels(:, 2:end), [1, 1, K]) == repmat(reshape(1:K, 1, 1, K), N, M - 1));
                tmp4 = sum(L(:, 2:end, :) .* mask, 3);
                
                c = abs(tmp1 - L(:, 1:(end - 1), alpha)) + abs(tmp3 - L(:, 2:end, alpha));
                d = abs(L(:, 1:(end - 1), alpha) - tmp2) + abs(L(:, 2:end, alpha) - tmp4);
                a = abs(tmp1 - tmp2) + abs(tmp3 - tmp4);
                hor = c + d - a;

                from = (1 : (N * (M - 1)))';
                to = from + N;
                terminal_weights(from, 2) = terminal_weights(from, 2) + reshape(a, [], 1);
                terminal_weights(from, 1) = terminal_weights(from, 1) + reshape(-c + a, [] ,1);
                terminal_weights(to, 1) = terminal_weights(to, 1) + reshape(c - a, [], 1);
                edge_weights = [edge_weights; from, to, zeros(size(from)), reshape(hor, [], 1)];

                delta = min(terminal_weights, [], 2);
                terminal_weights = terminal_weights - repmat(delta, 1, 2);
                delta_const = sum(delta);
                [cut, bin_labels] = graphCutMex(terminal_weights,edge_weights);
                labels(logical(reshape(bin_labels, N, M))) = alpha;
                energy_cur(en_num) = cut + delta_const;
                if display
                    fprintf('Energy = %d\n', energy_cur(end));
                end
                time_cur(en_num) = toc;
            end
            if (all(energy_cur((end - K + 1) : end) - energy_cur(end) == 0))
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