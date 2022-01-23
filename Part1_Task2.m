clc;
clear;
close all;

X = imread("autumn.tif");
I = rgb2gray(X);
figure(1);
imshow(I); title("Original Autumn image");

J = dct2(I);
figure(2);
colormap(jet(64)), imagesc(log(abs(J))), colorbar
title("The logarithmic magnitude of the 2D DCT of autumn");
daspect([1,1,1]);

img_F_mag = abs(fftshift(fft2(I)));
figure(3);
imagesc(log(img_F_mag));
title("The logarithmic magnitude of the 2D FFT of autumn (fftshift)");
daspect([1,1,1]);
colormap(jet(64));
colorbar;
