clc;
clear;
close all;

X = imread('cameraman.tif');
J = imresize(X, 0.5);
I = double(J);
signalPower = mean(abs(I(:)).^2);

% In the unit of dB
SNR1 = 20;
SNR2 = 10;
SNR3 = 5;
noisePower1 = signalPower / (10^(SNR1/10));
noisePower2 = signalPower / (10^(SNR2/10));
noisePower3 = signalPower / (10^(SNR3/10));

seed = 0;       rng(seed);
Im1 = I + sqrt(noisePower1) * randn(size(I));
% Im1 = Im1 / max(Im1(:));
Im2 = I + sqrt(noisePower2) * randn(size(I));
% Im2 = Im2 / max(Im2(:));
Im3 = I + sqrt(noisePower3) * randn(size(I));
% Im3 = Im3 / max(Im3(:));

figure(1);
tiledlayout(2,2,'TileSpacing','compact');
nexttile, imshow(J); title("Original image with 128 \times 128 pixles");
nexttile, imshow(Im1/max(Im1(:))); title("Image contaminated by " + SNR1 + " dB noise");
nexttile, imshow(Im2/max(Im2(:))); title("Image contaminated by " + SNR2 + " dB noise");
nexttile, imshow(Im3/max(Im3(:))); title("Image contaminated by " + SNR3 + " dB noise");

%% loacl threshold
% block size is 8 by 8 - noisy - local
blockSize1 = 8;
ratio = 4;
[Y1, dctCoeffsCompressed, num_of_non_zero_pixel] = dctCompression(Im1, ratio, blockSize1);
psnr1 = psnr(Y1/max(Y1(:)), im2double(J));
compressionRatio1 = (128 * 128) / num_of_non_zero_pixel;
disp("The comprssion ratio is measured as: " + compressionRatio1);
figure(2); imagesc(dctCoeffsCompressed); title(["Quantized DCT coefficient matrix"+ ...
    "(The block size is " + blockSize1 + " \times " + blockSize1 + ")", ...
    "(The local threshold is used to achieve the compression ratio of " + compressionRatio1 + ")"]); colormap(gray);
figure(3); imshow(Y1/max(Y1(:))); title(["Reconstructed Image (\sigma^2 = " + noisePower1 + ") " + ...
    "(" + blockSize1 + " \times " + blockSize1 + " block) " + ...
    "(PSNR = " + sprintf("%.2f", psnr1) + " dB)", ...
    "(The local threshold is used to achieve the compression ratio of " + compressionRatio1 + ")"]);    colormap(gray);

% block size is 16 by 16 - noisy - local
blockSize2 = 16;
ratio = 4;
[Y2, dctCoeffsCompressed, num_of_non_zero_pixel] = dctCompression(Im1, ratio, blockSize2);
psnr2 = psnr(Y2/max(Y2(:)), im2double(J));
compressionRatio2 = (128 * 128) / num_of_non_zero_pixel;
disp("The comprssion ratio is measured as: " + compressionRatio2);
figure(4); imagesc(dctCoeffsCompressed); title(["Quantized DCT coefficient matrix"+ ...
    "(The block size is " + blockSize2 + " \times " + blockSize2 + ")", ...
    "(The local threshold is used to achieve the compression ratio of " + compressionRatio2 + ")"]); colormap(gray);
figure(5); imshow(Y2/max(Y2(:))); title(["Reconstructed Image (\sigma^2 = " + noisePower1 + ") " + ...
    "(" + blockSize2 + " \times " + blockSize2 + " block) " + ...
    "(PSNR = " + sprintf("%.2f", psnr2) + " dB)", ...
    "(The local threshold is used to achieve the compression ratio of " + compressionRatio2 + ")"]);    colormap(gray);


%% block size is 8 by 8 - noisy - global threshold
ratio = 4;
[Y3, dctCoeffsCompressed, num_of_non_zero_pixel] = dctCompression2(Im1, ratio, blockSize1);
psnr3 = psnr(Y3/max(Y3(:)), im2double(J));
compressionRatio3 = (128 * 128) / num_of_non_zero_pixel;
disp("The comprssion ratio is measured as: " + compressionRatio3);
figure(6); imagesc(dctCoeffsCompressed); title(["Quantized DCT coefficient matrix"+ ...
    "(The block size is " + blockSize1 + " \times " + blockSize1 + ")", ...
    "(The global threshold is used to achieve the compression ratio of " + compressionRatio3 + ")"]); colormap(gray);
figure(7); imshow(Y3/max(Y3(:))); title(["Reconstructed Image (\sigma^2 = " + noisePower1 + ") " + ...
    "(" + blockSize1 + " \times " + blockSize1 + " block) " + ...
    "(PSNR = " + sprintf("%.2f", psnr3) + " dB)", ...
    "(The global threshold is used to achieve the compression ratio of " + compressionRatio3 + ")"]);    colormap(gray);

% block size is 16 by 16 - noisy - global
ratio = 4;
[Y4, dctCoeffsCompressed, num_of_non_zero_pixel] = dctCompression2(Im1, ratio, blockSize2);
psnr4 = psnr(Y4/max(Y4(:)), im2double(J));
compressionRatio4 = (128 * 128) / num_of_non_zero_pixel;
disp("The comprssion ratio is measured as: " + compressionRatio4);
figure(8); imagesc(dctCoeffsCompressed); title(["Quantized DCT coefficient matrix"+ ...
    "(The block size is " + blockSize2 + " \times " + blockSize2 + ")", ...
    "(The global threshold is used to achieve the compression ratio of " + compressionRatio4 + ")"]); colormap(gray);
figure(9); imshow(Y4/max(Y4(:))); title(["Reconstructed Image (\sigma^2 = " + noisePower1 + ") " + ...
    "(" + blockSize2 + " \times " + blockSize2 + " block) " + ...
    "(PSNR = " + sprintf("%.2f", psnr4) + " dB)", ...
    "(The global threshold is used to achieve the compression ratio of " + compressionRatio4 + ")"]);    colormap(gray);

%% Comparison
figure(10);
tiledlayout(2,2,'TileSpacing','compact');
nexttile, imshow(Y1); imshow(Y1/max(Y1(:))); title(["Reconstructed Image (\sigma^2 = " + noisePower1 + ") " + ...
    "(" + blockSize1 + " \times " + blockSize1 + " block) " + ...
    "(PSNR = " + sprintf("%.2f", psnr1) + " dB)", ...
    "(The local threshold is used to achieve the compression ratio of " + compressionRatio1 + ")"]);    colormap(gray);
nexttile, imshow(Y2/max(Y2(:))); title(["Reconstructed Image (\sigma^2 = " + noisePower1 + ") " + ...
    "(" + blockSize2 + " \times " + blockSize2 + " block) " + ...
    "(PSNR = " + sprintf("%.2f", psnr2) + " dB)", ...
    "(The local threshold is used to achieve the compression ratio of " + compressionRatio2 + ")"]);    colormap(gray);
nexttile, imshow(Y3/max(Y3(:))); title(["Reconstructed Image (\sigma^2 = " + noisePower1 + ") " + ...
    "(" + blockSize1 + " \times " + blockSize1 + " block) " + ...
    "(PSNR = " + sprintf("%.2f", psnr3) + " dB)", ...
    "(The global threshold is used to achieve the compression ratio of " + compressionRatio3 + ")"]);    colormap(gray);
nexttile, imshow(Y4/max(Y4(:))); title(["Reconstructed Image (\sigma^2 = " + noisePower1 + ") " + ...
    "(" + blockSize2 + " \times " + blockSize2 + " block) " + ...
    "(PSNR = " + sprintf("%.2f", psnr4) + " dB)", ...
    "(The global threshold is used to achieve the compression ratio of " + compressionRatio4 + ")"]);    colormap(gray);

%% Original Noisy and two local thresholding results
figure(11);
tiledlayout(1,3,'TileSpacing','compact');
nexttile, imshow(Im1/max(Im1(:))); title("Image contaminated by " + SNR1 + " dB noise");
nexttile, imshow(Y1); imshow(Y1/max(Y1(:))); title(["Reconstructed Image" + ...
    "(" + blockSize1 + " \times " + blockSize1 + " block) " + ...
    "(PSNR = " + sprintf("%.2f", psnr1) + " dB)", ...
    "(The local threshold is used to achieve the compression ratio of " + compressionRatio1 + ")"]);    colormap(gray);
nexttile, imshow(Y2/max(Y2(:))); title(["Reconstructed Image" + ...
    "(" + blockSize2 + " \times " + blockSize2 + " block) " + ...
    "(PSNR = " + sprintf("%.2f", psnr2) + " dB)", ...
    "(The local threshold is used to achieve the compression ratio of " + compressionRatio2 + ")"]);    colormap(gray);
%% functions
% local threshold
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