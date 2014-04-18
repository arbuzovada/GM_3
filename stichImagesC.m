function [resultImage, resultMask] = stichImagesC(images, seeds)
% This function stiches several images in one
% INPUT:
%    images: K-by-1 cell array of images
%    seeds: K-by-1 cell array of logical masks
%
% OUPUT:
%    resultImage: result image
%    resultMask: height-by-width array of integer, resultMask(i, j) = k, if
%        (i, j)-th pixel was taken from k-th original image

    K = size(images, 1);
    [N, M, ~] = size(images{1});
    
    % initialize unary potentials
    unary_pot = zeros(N, M, K);
    for i = 1 : K
        unary_pot(:, :, i) = seeds{i};
        [rows, cols] = find(seeds{i});
%         length(rows)
%         j = repmat([1 : N]', [1, M, length(rows)]);
%         h = repmat([1 : M], [N, 1, length(rows)]);
%         rows = repmat(rows, [N, M, 1]);
%         cols = repmat(cols, [N, M, 1]); 
        for j = 1 : N
            for h = 1 : M
%                 dist = (rows - j) .^ 2 + (cols - h) .^ 2;
%                 unary_pot(j, h, i) = 1 / (min(dist) + 1);
% %                 unary_pot(:, :, i) = (N ^ 2 + M ^ 2 - min(dist)) / (N ^ 2 + M ^ 2);
            end
        end
% %         unary_pot(:, :, i) = max(max(unary_pot(:, :, i))) - unary_pot(:, :, i) / (max(max(unary_pot(:, :, i))) - min(min(unary_pot(:, :, i))));
    end
    unary_pot = 1e4 * (1 - unary_pot);
    
    % initialize pair potentials
%     pair_pot = zeros(N, M, K);
%     for i = 1 : K
%         pair_pot(:, :, i) = rgb2gray(images{i});
%     end
%     vertC = ones(N - 1, M); 
%     horC = ones(N, M - 1);

    tmp = zeros(N, M, K);
    for i = 1 : K
        tmp(:, :, i) = double(rgb2gray(images{i}));
    end
    vertC = 1 - exp(-sum(abs(tmp(1 : (end - 1), :, :) - tmp(2 : end, :, :)), 3));
    horC = 1 - exp(-sum(abs(tmp(:, 1 : (end - 1), :) - tmp(:, 2 : end, :)), 3));
    metric = 1 - eye(K);
        
    options.display = true;
    options.randOrder = true;
    [resultMask, energy, time] = alphaBetaSwapGridPottsC(unary_pot, ...
        vertC, horC, metric, options);
%     [Elja_resultMask, Elja_energy, Elja_time] = alphaExpansionGridPotts(unary_pot, ...
%         vertC, horC, metric, options);
    save results.mat resultMask energy time
%     save Elja_results.mat Elja_resultMask Elja_energy Elja_time
    figure();
    energy = energy(energy > 0)
    plot(2 : length(energy), energy(2 : end), 'b')
    xlabel('iteration')
    ylabel('energy')
    title('Energy')
    legend('energy')
    figure();
    plot(1 : length(time), time, 'r')
    xlabel('iteration')
    ylabel('time')
    legend('time')
    
    resultImage = zeros(size(images{1}));
    for i = 1 : N
        for j = 1 : M
            resultImage(i, j, :) = images{resultMask(i, j)}(i, j, :);
        end
    end
end