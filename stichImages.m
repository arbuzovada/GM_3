function [resultImage, resultMask] = stichImages(images, seeds)
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
    end
    
    % initialize pair potentials
%     pair_pot = zeros(N, M, K);
%     for i = 1 : K
%         pair_pot(:, :, i) = rgb2gray(images{i});
%     end
    vertC = ones(N - 1, M);
    horC = ones(N, M - 1);
    metric = eye(K);
    
    options.display = true;
    [resultMask, energy, time] = alphaBetaSwapGridPottsC(unary_pot, ...%pair_pot, ...
        vertC, horC, metric, options);
    figure();
    plot(1 : length(energy), energy, 'b')
    xlabel('iteration')
    ylabel('energy')
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