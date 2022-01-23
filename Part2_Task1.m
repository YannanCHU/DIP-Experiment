clc;
clear;
close all;

I = imread('cameraman.tif');
% A roughly equal number of pixels is mapped to each of the n levels in J, 
% so that the histogram of J is approximately flat. 
% The histogram of J is flatter when n is much smaller than the number of discrete levels in I.
J = histeq(I, 32);
K = histeq(I,16);
Q = histeq(I, 8);
R = histeq(I,4);

figure();
subplot(2,2,1), [Icounts,IbinLocations] = imhist(I,256);
stem(IbinLocations, Icounts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title("Histogram of the original image");
xlabel("Pixel intensity");  ylabel("Number of pixels");
subplot(2,2,2), imshow(I); title("Original Image")

figure();
subplot(2,2,1), [Jcounts,JbinLocations] = imhist(J,32);
stem(JbinLocations, Jcounts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title(size(I,1) + " \times " + size(I,2) + " pixels are mapped to 32 intensty levels");
xlabel("Pixel intensity");  ylabel("Number of pixels");
subplot(2,2,2), imshow(J); title("Image with 32 intensty levels")

subplot(2,2,3), [Kcounts,KbinLocations] = imhist(K,32);
stem(KbinLocations, Kcounts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title(size(I,1) + " \times " + size(I,2) + " pixels are mapped to 16 intensty levels");
xlabel("Pixel intensity");  ylabel("Number of pixels");
subplot(2,2,4), imshow(K);  title("Image with 16 intensty levels")

figure();
subplot(2,2,1), [Qcounts,QbinLocations] = imhist(Q,32);
stem(QbinLocations, Qcounts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title(size(I,1) + " \times " + size(I,2) + " pixels are mapped to 8 intensty levels");
xlabel("Pixel intensity");  ylabel("Number of pixels");
subplot(2,2,2), imshow(Q);  title("Image with 8 intensty levels")

subplot(2,2,3), [Rcounts,RbinLocations] = imhist(R,32);
stem(RbinLocations, Rcounts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title(size(I,1) + " \times " + size(I,2) + " pixels are mapped to 4 intensty levels");
xlabel("Pixel intensity");  ylabel("Number of pixels");
subplot(2,2,4), imshow(R);  title("Image with 4 intensty levels")

%% Task 2       carry a histogram modification of the intensity of a given image.
% a) a data file
hist_data_1 = 32:-1:1;
hist_data_2 = 1:1:4;
hist_data_1 = hist_data_1 / sum(hist_data_1);
hist_data_2 = hist_data_2 / sum(hist_data_2);

J2 = histeq(I, hist_data_1);
K2 = histeq(I, hist_data_2);
figure();
subplot(2,3,1); stem(hist_data_1, 'LineWidth', 2, 'Marker', 'none');
title("Given data set which specifies the resulting histogram");
xlabel("Index");  ylabel("Data Value");
subplot(2,3,2); [J2counts,J2binLocations] = imhist(J2,32);
stem(J2binLocations, J2counts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title({size(I,1) + " \times " + size(I,2) + " pixels are mapped to 32 intensty levels", "(Histogram is given as a data file)"});
xlabel("Pixel intensity");  ylabel("Number of pixels");
subplot(2,3,3); imshow(J2); title("Image with 32 intensty levels")

subplot(2,3,4); stem(hist_data_2, 'LineWidth', 2, 'Marker', 'none');
title("Given data set which specifies the resulting histogram");
xlabel("Index");  ylabel("Data Value");
subplot(2,3,5); [K2counts,K2binLocations] = imhist(K2,32);
stem(K2binLocations, K2counts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title({size(I,1) + " \times " + size(I,2) + " pixels are mapped to 4 intensty levels", "(Histogram is given as a data file)"});
xlabel("Pixel intensity");  ylabel("Number of pixels");
subplot(2,3,6), imshow(K2); title("Image with 4 intensty levels")


%%
% b) a function
syms x

f1 = x.^2;
f2 = x;

hist_fun_1 = double(subs(f1, x, {linspace(1,32,32)}));
hist_fun_2 = double(subs(f2, x, {linspace(1,32,32)}));
hist_fun_1 = hist_fun_1 / sum(hist_fun_1);
hist_fun_2 = hist_fun_2 / sum(hist_fun_2);

J3 = histeq(I, hist_fun_1);
K3 = histeq(I, hist_fun_2);
figure();
subplot(2,3,1); plot(linspace(1,32,32), hist_fun_1);   
title(["Function which specifies the resulting histogram"]);
xlabel("x");  ylabel("y");
subplot(2,3,2); [J3counts,J3binLocations] = imhist(J3,32);
stem(J3binLocations, J3counts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title({size(I,1) + " \times " + size(I,2) + " pixels are mapped to 32 intensty levels", "(Histogram is given as a function)"});
xlabel("Pixel intensity");  ylabel("Number of pixels");
subplot(2,3,3); imshow(J3); title("Image with 32 intensty levels")

subplot(2,3,4); plot(linspace(1,32,32), hist_fun_2);   
title(["Function which specifies the resulting histogram"]);
xlabel("x");  ylabel("y");
subplot(2,3,5); [K3counts,K3binLocations] = imhist(K3,32);
stem(K3binLocations, K3counts, 'LineWidth', 2, 'Marker', 'none');   xlim([0 256]);
title({size(I,1) + " \times " + size(I,2) + " pixels are mapped to 32 intensty levels", "(Histogram is given as a function)"});
xlabel("Pixel intensity");  ylabel("Number of pixels");
subplot(2,3,6), imshow(K3); title("Image with 32 intensty levels");

