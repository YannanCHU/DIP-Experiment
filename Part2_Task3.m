clc;
clear;
close all;

I = imread('trees.tif');
figure(11);
imshow(I);  title("Original Image");

figure(1);
tiledlayout(2,2,'TileSpacing','compact');
nexttile, imshow(edge(I,'sobel'));    title("sobel operator");
nexttile, imshow(edge(I,'roberts'));  title("roberts operator");
nexttile, imshow(edge(I,'prewitt'));  title("prewitt operator");
nexttile, imshow(edge(I,'Log'));      title("Log operator");
%%
% I = imread('test.png');
I=im2double(I);
[~, thresholdVal]=edge(I, 'sobel');

kernel_neg_45 = [0 1 2; -1 0 1; -2 -1 0] / 8;
kernel_pos_45 = [-2 -1 0; -1 0 1; 0 1 2] / 8;

grad_neg_45 = imfilter(I, kernel_neg_45, 'replicate');
grad_pos_45 = imfilter(I, kernel_pos_45, 'replicate');
grad_diag = abs(grad_neg_45) + abs(grad_pos_45);

figure(2);
subplot(1,3,1);     imshow(bwskel(imbinarize(grad_neg_45)));    title("Edges in an image inclined at - 45 degree");
subplot(1,3,2);     imshow(bwskel(imbinarize(grad_pos_45)));    title("Edges in an image inclined at + 45 degree")
subplot(1,3,3);     
% imshow(grad_diag >= thresholdVal);      
imshow(bwskel(imbinarize(grad_diag)));    
title("Edges in an image inclined at -/+ 45 degree")

figure();
imshow(bwskel(imbinarize(grad_diag)));    
 

% imshow(edge_neg_45);
