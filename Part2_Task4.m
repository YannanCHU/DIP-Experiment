clc;
clear;
close all;

X = imread('autumn.tif');
I = rgb2gray(X);
J = imnoise(I,'salt & pepper');
K = medfilt2(J,[3 3]);
K2 = medfilt2(J,[5 5]);
K3 = medfilt2(J,[7 7]);
tiledlayout(2,2,'TileSpacing','compact') 
nexttile, imshow(J);  title("Original Noisy Image (" + size(I,1) + " \times " + size(I,2) + " pixels)");  
nexttile, imshow(K);  title("Median filter with 3 \times 3 neighbourhood size");
nexttile, imshow(K2);  title("Median filter with 5 \times 5 neighbourhood size");
nexttile, imshow(K3);  title("Median filter with 7 \times 7 neighbourhood size");
