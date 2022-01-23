clc;
clear;
close all;

X = imread('autumn.tif');
I = rgb2gray(X);
d = 0.05;
J = imnoise(I,'salt & pepper',d);
K = medfilt2(J,[3 3]);
K2 = medfilt2(J,[5 5]);
K3 = medfilt2(J,[7 7]);
figure();
tiledlayout(1,4,'TileSpacing','compact');
nexttile, imshow(J);  title(["Original Noisy Image (" + size(I,1) + " \times " + size(I,2) + " pixels)", ...
    "(Noise Density = " + d + ")"]);  
nexttile, imshow(K);  title("Median filter with 3 \times 3 neighbourhood size");
nexttile, imshow(K2);  title("Median filter with 5 \times 5 neighbourhood size");
nexttile, imshow(K3);  title("Median filter with 7 \times 7 neighbourhood size");
%%
d = 0.2;
J = imnoise(I,'salt & pepper',d);
K = medfilt2(J,[3 3]);
K2 = medfilt2(J,[5 5]);
K3 = medfilt2(J,[7 7]);
figure();
tiledlayout(1,4,'TileSpacing','compact');
nexttile, imshow(J);  title(["Original Noisy Image (" + size(I,1) + " \times " + size(I,2) + " pixels)", ...
    "(Noise Density = " + d + ")"]);  
nexttile, imshow(K);  title("Median filter with 3 \times 3 neighbourhood size");
nexttile, imshow(K2);  title("Median filter with 5 \times 5 neighbourhood size");
nexttile, imshow(K3);  title("Median filter with 7 \times 7 neighbourhood size");
%%
% SNR in dB
snr = 20;
snr2 = 10 ^ (snr/10);
% compute the signal power
signalPower = mean(abs(double(I(:))).^2);
noiseVar = signalPower / snr2;
J2 = double(I) + sqrt(noiseVar) * randn(size(I));
K = medfilt2(J2,[3 3]);
K2 = medfilt2(J2,[5 5]);
K3 = medfilt2(J2,[7 7]);
figure();
tiledlayout(1,4,'TileSpacing','compact');
nexttile, imshow(J2/255);  title(["Image contaminated by Gaussian Noise", "(SNR = " + snr + " dB)"]);  
nexttile, imshow(K/255);  title("Median filter with 3 \times 3 neighbourhood size");
nexttile, imshow(K2/255);  title("Median filter with 5 \times 5 neighbourhood size");
nexttile, imshow(K3/255);  title("Median filter with 7 \times 7 neighbourhood size");