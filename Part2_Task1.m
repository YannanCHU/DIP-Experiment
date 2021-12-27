clc;
clear;
close all;

I = imread('cameraman.tif');
% A roughly equal number of pixels is mapped to each of the n levels in J, 
% so that the histogram of J is approximately flat. 
% The histogram of J is flatter when n is much smaller than the number of discrete levels in I.
J = histeq(I, 32);
K = histeq(I,4);
figure(1);
subplot(2,2,1), [Jcounts,JbinLocations] = imhist(J,32);
stem(JbinLocations, Jcounts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title(size(I,1) + " \times " + size(I,2) + " pixels are mapped to 32 intensty levels");
subplot(2,2,2), imshow(J);
subplot(2,2,3), [Kcounts,KbinLocations] = imhist(K,32);
stem(KbinLocations, Kcounts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title(size(I,1) + " \times " + size(I,2) + " pixels are mapped to 4 intensty levels");
subplot(2,2,4), imshow(K);

%% Task 2       carry a histogram modification of the intensity of a given image.
% a) a data file
hist_data_1 = 32:-1:1;
hist_data_2 = 1:1:4;

J2 = histeq(I, hist_data_1);
K2 = histeq(I, hist_data_2);
figure(2);
subplot(2,2,1); [J2counts,J2binLocations] = imhist(J2,32);
stem(J2binLocations, J2counts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title({size(I,1) + " \times " + size(I,2) + " pixels are mapped to 32 intensty levels", "(Histogram is given as a data file)"});
subplot(2,2,2); imshow(J2);
subplot(2,2,3); [K2counts,K2binLocations] = imhist(K2,32);
stem(K2binLocations, K2counts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title({size(I,1) + " \times " + size(I,2) + " pixels are mapped to 4 intensty levels", "(Histogram is given as a data file)"});
subplot(2,2,4), imshow(K2);


%%
% b) a function
syms x

f1 = x.^2;
f2 = x;

hist_fun_1 = double(subs(f1, x, {(1:1:32)}));
hist_fun_2 = double(subs(f2, x, {[1,2,3,4]}));

J3 = histeq(I, hist_fun_1);
K3 = histeq(I, hist_fun_2);
figure(3);
subplot(2,2,1); [J3counts,J3binLocations] = imhist(J3,32);
stem(J3binLocations, J3counts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title({size(I,1) + " \times " + size(I,2) + " pixels are mapped to 32 intensty levels", "(Histogram is given as a function)"});
subplot(2,2,2); imshow(J3)
subplot(2,2,3); [K3counts,K3binLocations] = imhist(K3,32)
stem(K3binLocations, K3counts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title({size(I,1) + " \times " + size(I,2) + " pixels are mapped to 4 intensty levels", "(Histogram is given as a function)"});
subplot(2,2,4), imshow(K3);

