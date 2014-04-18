function [labels, energy, time] = alphaExpansionGridPotts(unary, vertC, horC, metric, options)
% ¬’ќƒ
% unary Ч унарные потенциалы, массив типа double размера N x M x K, где N Ч высота решетки, M Ч ширина решетки, K Ч количество меток.
% vertC Ч коэффициенты  c_{pq}, соответствующие вертикальным ребрам, массив типа double размера (N Ч 1) x M;
% horC Ч коэффициенты  c_{pq}, соответствующие горизонтальным ребрам, массив типа double размера N x (M Ч 1);
% metric Ч рассто€ние между метками соседних переменных, массив типа double размера K x K;
% options Ч (необ€зательный аргумент) набор дополнительных параметров, структура с пол€ми
%   'maxIter' Ч максимально допустимое число итераций алгоритма (по умолчанию = 500);
%   'display' Ч параметр типа logical: если true, то при каждом запуске алгоритма разреза графа нужно выводить на экран номер итерации, номера обрабатываемых меток, текущее значение энергии;
%   'numStart' Ч количество запусков из разных начальных приближений;
%   'randOrder' Ч параметр типа logical: если true, то при каждом запуске использовать случайный пор€док меток ? и ?;
% ¬џ’ќƒ
% labels Ч разметка, обладающа€ наименьшей энергией, массив типа double размера N x M;
% energy Ч значени€ энергии на каждой итерации, вектор типа double длины, равной количеству итераций алгоритма;
% time Ч врем€, пройденное с начала работы алгоритма до каждой итерации, вектор типа double длины, равной количеству итераций алгоритма;

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
                    fprintf('Start є %d, iteration %d, expanding label = %d, processing...\n', start, en_num, alpha);
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