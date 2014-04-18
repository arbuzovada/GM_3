K = 3;
images = cell(K, 1);
seeds = cell(K, 1);
for i = 1 : K
    images{i} = imread(['m', num2str(i), '.jpg']);
%     size(images{i})
    cur = imread(['ms', num2str(i), '.jpg']);
    seeds{i} = logical(1 - cur(:, :, 1) / 255);
%     size(seeds{i})
end
save model.mat seeds images
% images{1} = imread('dd1_s.jpg');
% images{1} = imread('mimg1.jpg');
% [N, M, ~] = size(images{1})
% images{2} = imread('dd2_s.jpg');
% images{2} = imread('mimg2.jpg');
% seeds = cell(K, 1);
% cur = imread('sseeds1.jpg');
% seeds{1} = logical(1 - cur(:, :, 1) / 255);
% cur = imread('sseeds2.jpg');
% seeds{2} = logical(1 - cur(:, :, 1) / 255);

% d_s
% seeds{1}(159 : 430, 302 : 430) = true;
% seeds{2}(120 : 165, 215 : 250) = true;

% seeds{1} = false(N, M);
% seeds{2} = false(N, M);

% pillows
% seeds{1}(185 : 264, 225 : 280) = true;
% seeds{2}(120 : 165, 215 : 250) = true;

% seeds{1}(120 : 165, 215 : 250) = true;
% seeds{2}(185 : 264, 225 : 280) = true;
% seeds{2}(140 : 270, 260 : 340) = true;

% seeds{1}(200 : 290, 46 : 82) = true;
% seeds{1}(232 : 300, 44 : 210) = true;
% seeds{1}(210 : 405, 390 : 510) = true;
% seeds{2}(300 : 435, 70 : 123) = true;
% seeds{2}(384 : 450, 66 : 315) = true;

[resultImage, resultMask] = stichImagesC(images, seeds);
figure();
image(resultImage / 255)
imwrite(resultImage / 255, 'res_mC.jpg')
resultMask(resultMask == 2) = 10;
resultMask(resultMask == 3) = 20;
resultMask(resultMask == 4) = 30;
resultMask(resultMask == 5) = 40;
figure()

image(resultMask)
axis equal
print('mask_mC', '-depsc2', '-r300');