clc;
clear;
close all;

X = imread("autumn.tif");
I = rgb2gray(X);
J = dct2(I);
figure(1);
colormap(jet(64)), imagesc(log(abs(J))), colorbar
title("The logarithmic magnitude of the 2D DCT of autumn")

img_F_mag = abs((fft2(I)));
figure(2);
imagesc(log(img_F_mag));
title("The logarithmic magnitude of the 2D FFT of autumn (No fftshift)");
colormap(jet(64));
colorbar;
