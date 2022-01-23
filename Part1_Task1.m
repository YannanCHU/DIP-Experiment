clc;
clear;
close all;

%% a
img = imread("Cameraman.tif");
figure();
imshow(img);    title("Original Image");
% img = im2double(img);

img_F = fftshift(fft2(img));
img_F_mag = abs(img_F);

figure();
imagesc(log(img_F_mag));
title("The logarithmic magnitude of the 2D FFT of Cameraman");
colormap jet
colorbar;

figure();
surf(log(img_F_mag), 'LineStyle', 'none');
title("The logarithmic magnitude of the 2D FFT of Cameraman");
colormap jet
colorbar;

min(log(img_F_mag(:)))

%% b
clc;
clear;
% close  all;

img_grid_impulse = zeros(256, 256);
img_grid_impulse(1:8:256, 1:8:256) = 255;

figure();
imagesc(img_grid_impulse); title("Artificial image of a grid of impulses"); colormap gray;

img_grid_impulse_F = fftshift(fft2(img_grid_impulse));
img_grid_impulse_F_mag = abs(img_grid_impulse_F);

figure();
imagesc(log(img_grid_impulse_F_mag));
title("The logarithmic magnitude of the 2D FFT of grid of impulses");
colormap jet
colorbar;

figure();
surf((img_grid_impulse_F_mag), 'LineStyle', 'none');
title("The magnitude of the 2D FFT of grid of impulses");
colormap jet

% grid lines
img_grid_lines = zeros(256, 256);
img_grid_lines(1:8:256, :) = 255;
img_grid_lines(:, 1:8:256) = 255;

figure();
imagesc(img_grid_lines); title("Artificial image of a grid of lines"); colormap gray;

img_grid_lines_F = fftshift(fft2(img_grid_lines));
img_grid_lines_F_mag = abs(img_grid_lines_F);

figure();
imagesc(log(img_grid_lines_F_mag));
title("The logarithmic magnitude of the 2D FFT of grid of lines");
colormap jet
colorbar;

%% non-periodic
% random noise-like image
img_randi = randi([0, 255], 256, 256);
figure();
imagesc(img_randi); title("Artificial image of pseudo-random integers"); colormap gray;

img_randi_F = fftshift(fft2(img_randi));
img_randi_F_mag = abs(img_randi_F);

figure();
surf(log(img_randi_F_mag), 'LineStyle', 'none');
title("Logarithmic magnitude of the 2D FFT of noise-like image");
colormap jet
colorbar;




index = 2.^(0:1:8);
img_grid_lines2 = zeros(256, 256);
img_grid_lines2(index, :) = 255;
% img_grid_lines2(:, index) = 255;

figure();
imagesc(img_grid_lines2); title("Artificial image of lines with varying distance"); colormap gray;

img_grid_lines2_F = fftshift(fft2(img_grid_lines2));
img_grid_lines2_F_mag = abs(img_grid_lines2_F);

figure();
imagesc(log(img_grid_lines2_F_mag));
title("The logarithmic magnitude of the 2D FFT of lines with varying distance");
colormap jet
colorbar;

% 
%% c
clc;
clear;

imgA = imread("Cameraman.tif");
imgB = imread("rice.png");
figure();
subplot(1,2,1);     imshow(imgA);   title("Image A");
subplot(1,2,2);     imshow(imgB);   title("Image B");

% imgA = im2double(imgA);
% imgB = im2double(imgB);

imgA_F = fft2(imgA);
imgB_F = fft2(imgB);

imgA_F_mag = abs(imgA_F);
imgB_F_mag = abs(imgB_F);

imgA_F_phase = angle(imgA_F);
imgB_F_phase = angle(imgB_F);

img_complex = imgB_F_mag .* exp(1j * imgA_F_phase);
img_ifft2 = real(ifft2(img_complex));

figure();   imshow(img_ifft2/255);   title("Reconstructed Image with A's FFT phase and B's FFT amplitude");  colormap gray;
