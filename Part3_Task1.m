clc;
clear;
close all;

X = imread('autumn.tif');
I = rgb2gray(X);
% DCT magnitude threshold = 1
th1 = 1;
J = dct2(I);
nz = find(abs(J)<th1);
num_of_zero_pixel1 = length(nz);
J(nz) = zeros(size(nz));
K1 = idct2(J)/255;
psnr1 = psnr(K1, double(I)/255);


% DCT magnitude threshold = 10
th2 = 10;
J = dct2(I);
nz = find(abs(J)<th2);
num_of_zero_pixel2 = length(nz);
J(nz) = zeros(size(nz));
K2 = idct2(J)/255;
psnr2 = psnr(K2, double(I)/255);

% DCT magnitude threshold = 100
th3 = 100;
J = dct2(I);
nz = find(abs(J)<th3);
num_of_zero_pixel3 = length(nz);
J(nz) = zeros(size(nz));
K3 = idct2(J)/255;
psnr3 = psnr(K3, double(I)/255);

J = dct2(I);
disp("The maximum DCT magnitude is: " + max(max(abs(J))));
disp("The mean DCT magnitude is: " + mean(abs(J(:))));
disp("The median DCT magnitude is: " + median(abs(J(:))));

figure(1);
tiledlayout(2,2,'TileSpacing','compact');
nexttile; imshow(double(I)/255), axis off;    title("Original Image");
nexttile; imshow(K1), axis off;    title(["Compressed Image (Threshold of DCT magnitude is " + th1 + ", PSNR = " + sprintf("%.2f", psnr1) + " dB)",...
    "(" + (num_of_zero_pixel1 / (prod(size(I),2))*100) + "% DCT coefficients are set as 0)"]);
nexttile; imshow(K2), axis off;    title(["Compressed Image (Threshold of DCT magnitude is " + th2 + ", PSNR = " + sprintf("%.2f", psnr2) + " dB)",...
    "(" + (num_of_zero_pixel2 / (prod(size(I),2))*100) + "% DCT coefficients are set as 0)"]);
nexttile; imshow(K3), axis off;    title(["Compressed Image (Threshold of DCT magnitude is " + th3 + ", PSNR = " + sprintf("%.2f", psnr3) + " dB)",...
    "(" + (num_of_zero_pixel3 / (prod(size(I),2))*100) + "% DCT coefficients are set as 0)"]);

%%

I = rgb2gray(X);
[imgZ, ~] = dctCompression(I, 0, 64);
th1 = 20;
blockSize1 = 8;
blockSize2 = 16;
blockSize3 = 64;
blockSize4 = 256;
[imgZ1, num_of_zero_pixel1] = dctCompression(I, th1, blockSize1);
PSNR1 = psnr(imgZ1/255, double(I)/255);
[imgZ2, num_of_zero_pixel2] = dctCompression(I, th1, blockSize2);
PSNR2 = psnr(imgZ2/255, double(I)/255);
[imgZ3, num_of_zero_pixel3] = dctCompression(I, th1, blockSize3);
PSNR3 = psnr(imgZ3/255, double(I)/255);
[imgZ4, num_of_zero_pixel4] = dctCompression(I, th1, blockSize4);
PSNR4 = psnr(imgZ4/255, double(I)/255);

figure(2);
tiledlayout(1,4,'TileSpacing','compact');
% nexttile; imshow(double(I)/255), axis off;    title("Original Image");
nexttile; imshow(imgZ1/255), axis off;    title(["(PSNR = " + sprintf("%.2f", PSNR1) + " dB)", "Compressed Image (Thr = " + th1 + ", " + ...
    blockSize1 + " \times " + blockSize1 + " block)"]);

nexttile; imshow(imgZ2/255), axis off;    title(["(PSNR = " + sprintf("%.2f", PSNR2) + " dB)", "Compressed Image (Thr = " + th1 + ", " + ...
    blockSize2 + " \times " + blockSize2 + " block)"]);

nexttile; imshow(imgZ3/255), axis off;    title(["(PSNR = " + sprintf("%.2f", PSNR3) + " dB)", "Compressed Image (Thr = " + th1 + ", " + ...
    blockSize3 + " \times " + blockSize3 + " block)"]);

% title(["Compressed Image (Threshold of = " + th1 + ", PSNR = " + sprintf("%.2f", PSNR3) + " dB)",...
%     "(DCT is applied to " + blockSize3 + " \times " + blockSize3 + " blocks of image individually)",...
%     "(" + (num_of_zero_pixel3 / (prod(size(I),2))*100) + "% DCT coefficients are set as 0)"]);
nexttile; imshow(imgZ4/255), axis off;    title(["(PSNR = " + sprintf("%.2f", PSNR4) + " dB)", "Compressed Image (Thr = " + th1 + ", " + ...
    blockSize4 + " \times " + blockSize4 + " block)"]);

% title(["Compressed Image (Threshold = " + th1 + ", PSNR = " + sprintf("%.2f", PSNR4) + " dB)",...
%     "(DCT is applied to " + blockSize4 + " \times " + blockSize4 + " blocks of image individually)",...
%     "(" + (num_of_zero_pixel4 / (prod(size(I),2))*100) + "% DCT coefficients are set as 0)"]);

function [imgZ, num_of_zero_pixel] = dctCompression(I, thr, N)
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
        dctCoeffs(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1) = dct2(imgX(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1));
    end
end

nz = find(abs(dctCoeffs)<thr);
num_of_zero_pixel = length(nz);
% dctCoeffs(nz) = zeros(size(nz));
dctCoeffs = wthresh(dctCoeffs,'h',thr);

if thr == 0
    disp("The maximum DCT magnitude is: " + max(max(abs(dctCoeffs))));
    disp("The mean DCT magnitude is: " + mean(abs(dctCoeffs(:))));
    disp("The median DCT magnitude is: " + median(abs(dctCoeffs(:))));
end

% inverse 2D DCT
for i = 1:1:length(index_row)-1
    for j = 1:1:length(index_col)-1
        imgZ(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1) = idct2(dctCoeffs(index_row(i):index_row(i+1)-1, index_col(j):index_col(j+1)-1));
    end
end

end