clc;
clear;
close all;

X = imread('cameraman.tif');
I = imresize(X, 0.5);
I = im2double(I);

signalPower = mean(abs(I(:)).^2);
SNR1 = 10;
noisePower1 = signalPower / (10^(SNR1/10));
I = imnoise(I,'gaussian',0,noisePower1);
figure(1); imshow(I); title("Original image with 128 \times 128 pixles");  colormap(gray);



blockSize1 = 8;
ratio = 4;
[Y, dctCoeffsCompressed, num_of_non_zero_pixel] = dctCompression(I, ratio, blockSize1);
psnr1 = psnr(Y/max(Y(:)), im2double(I));
compressionRatio1 = (128 * 128) / num_of_non_zero_pixel;
disp("The comprssion ratio is measured as: " + compressionRatio1);

% block size is 16 by 16 - noiseless
blockSize2 = 16;
ratio = 4;
[Y2, dctCoeffsCompressed2, num_of_non_zero_pixel] = dctCompression(I, ratio, blockSize2);
psnr2 = psnr(Y2/max(Y2(:)), im2double(I));
compressionRatio2 = (128 * 128) / num_of_non_zero_pixel;
disp("The comprssion ratio is measured as: " + compressionRatio2);

% block size is 32 by 32 - noiseless
blockSize3 = 32;
ratio = 4;
[Y3, dctCoeffsCompressed3, num_of_non_zero_pixel] = dctCompression(I, ratio, blockSize3);
psnr3 = psnr(Y3/max(Y3(:)), im2double(I));
compressionRatio3 = (128 * 128) / num_of_non_zero_pixel;
disp("The comprssion ratio is measured as: " + compressionRatio1);

% block size is 4 by 4 - noiseless
blockSize4 = 4;
ratio = 4;
[Y4, dctCoeffsCompressed4, num_of_non_zero_pixel] = dctCompression(I, ratio, blockSize4);
psnr4 = psnr(Y4/max(Y4(:)), im2double(I));
compressionRatio4 = (128 * 128) / num_of_non_zero_pixel;
disp("The comprssion ratio is measured as: " + compressionRatio4);


figure(6);
tiledlayout(1,4,'TileSpacing','compact');
% nexttile, imshow(I); title("Original image with 128 \times 128 pixles");
nexttile, imshow(Y4/max(Y4(:))); title(["Reconstructed Image " + ...
    "(" + blockSize4 + " \times " + blockSize4 + " block) " ,...
    "(PSNR = " + sprintf("%.2f", psnr4) + " dB)"]);    colormap(gray);
nexttile, imshow(Y/max(Y(:))); title(["Reconstructed Image " + ...
    "(" + blockSize1 + " \times " + blockSize1 + " block) " ,...
    "(PSNR = " + sprintf("%.2f", psnr1) + " dB)"]);    colormap(gray);
nexttile, imshow(Y2/max(Y2(:))); title(["Reconstructed Image " + ...
    "(" + blockSize2 + " \times " + blockSize2 + " block) " ,...
    "(PSNR = " + sprintf("%.2f", psnr2) + " dB)"]);    colormap(gray);
nexttile, imshow(Y3/max(Y3(:))); title(["Reconstructed Image " + ...
    "(" + blockSize3 + " \times " + blockSize3 + " block) " ,...
    "(PSNR = " + sprintf("%.2f", psnr3) + " dB)"]);    colormap(gray);

function [imgZ, dctCoeffs, num_of_non_zero_pixel] = dctCompression(I, ratio, N)
img = I;
imgX = double(img);
[imgX_rows, imgX_cols] = size(imgX);
% DCT coefficient matrix
% N = 8;
index_row = 1:N:imgX_rows;
index_col = 1:N:imgX_cols;
if index_row(end) ~= imgX_rows
    index_row = [index_row, imgX_rows+1];
end
if index_col(end) ~= imgX_cols
    index_col = [index_col, imgX_cols+1];
end

dctCoeffs = zeros(imgX_rows, imgX_cols);
imgZ = dctCoeffs;

for i = 1:1:length(index_row)-1
    for j = 1:1:length(index_col)-1
        dctMatrix = dct2(imgX(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1));
        dctSorted = sort(abs(dctMatrix(:)), 'descend');
        thr = dctSorted(round(N * N / ratio) + 1);
        if thr <= 0
            disp("smaller than 0")
        end
        dctMatrix = wthresh(dctMatrix,'h',abs(thr));
        dctCoeffs(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1) = dctMatrix;
    end
end

nz = find(abs(dctCoeffs)>0);
num_of_non_zero_pixel = length(nz);

% inverse 2D DCT
for i = 1:1:length(index_row)-1
    for j = 1:1:length(index_col)-1
        imgZ(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1) = idct2(dctCoeffs(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1));
    end
end

end

% global threshold
function [imgZ, dctCoeffs, num_of_non_zero_pixel] = dctCompression2(I, ratio, N)
img = I;
imgX = double(img);
[imgX_rows, imgX_cols] = size(imgX);
% DCT coefficient matrix
% N = 8;
index_row = 1:N:imgX_rows;
index_col = 1:N:imgX_cols;
if index_row(end) ~= imgX_rows
    index_row = [index_row, imgX_rows+1];
end
if index_col(end) ~= imgX_cols
    index_col = [index_col, imgX_cols+1];
end

dctCoeffs = zeros(imgX_rows, imgX_cols);
imgZ = dctCoeffs;

for i = 1:1:length(index_row)-1
    for j = 1:1:length(index_col)-1
        dctMatrix = dct2(imgX(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1));
        dctCoeffs(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1) = dctMatrix;
    end
end

dctSorted = sort(abs(dctCoeffs(:)), 'descend');
thr = dctSorted(round(imgX_rows * imgX_cols / ratio) + 1);
dctCoeffs = wthresh(dctCoeffs,'h',thr);

nz = find(abs(dctCoeffs)>0);
num_of_non_zero_pixel = length(nz);

% inverse 2D DCT
for i = 1:1:length(index_row)-1
    for j = 1:1:length(index_col)-1
        imgZ(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1) = idct2(dctCoeffs(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1));
    end
end

end