K = 2;
images = cell(K, 1);
images{1} = imread('img1.jpg');
[N, M, ~] = size(images{1})
images{2} = imread('img2.jpg');
seeds = cell(K, 1);
seeds{1} = false([N, M]);
seeds{1}(185 : 264, 225 : 280) = true;
seeds{2} = false([N, M]);
seeds{2}(120 : 165, 215 : 250) = true;

[resultImage, resultMask] = stichImages(images, seeds);
image(resultImage)