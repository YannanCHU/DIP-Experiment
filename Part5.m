clc;
clear;
close all;

I = imread('cameraman.tif');
I = double(I);
signalPower = mean(abs(I(:)).^2);

% 20 dB noise addition
seed = 0;       rng(seed);
SNR1 = 20;
noisePower1 = signalPower / (10^(SNR1/10));
noiseSigma = sqrt(noisePower1);
% noiseSigma = 10;
noise = noiseSigma * randn(size(I));
noisePow_check = 10 * log10(mean(noise(:).^2));
In = I + noise;
% In = awgn(I,20,'measured');

% 2D Gaussian point spread function for image degradation
N1 = 5;
N2 = 7;
sigma = 1;
h1 = fspecial('gaussian', [N1 N1], sigma);
h2 = fspecial('gaussian', [N2 N2], sigma);

Ig1 = imfilter(I,h1,'circular','conv');
Ig2 = imfilter(I,h2,'circular','conv');
% Ig2 = imgaussfilt(I, sigma, 'FilterSize', N);

% add the noise to the blurred image
% noise = 0;
Ign1 = Ig1 + noise;
BSNR1 = bsnr(Ig1, noisePower1);
Ign2 = Ig2 + noise;
BSNR2 = bsnr(Ig2, noisePower1);

figure(1);
tiledlayout(2,2,'TileSpacing','compact');
nexttile, imshow(I/255); title("Original image with 256 \times 256 pixles");
nexttile, imshow(In/255); title("Noisy image with 20 dB noise");
nexttile, imshow(Ig1/255); title("Degraded image by 2-D Gaussian filtering");
nexttile, imshow(Ign1/255); title("Degraded and noisy image");

Ign1_fre = fft2(Ign1);
h1_fre = fft2(h1,size(Ign1,1),size(Ign1,2));
figure();
surf(fftshift(abs(h1_fre)), 'LineStyle', 'None');  title("2D DFT amplitude of degradation model - h");
% Ign_deblurred_f1 = abs(Ign1_fre ./ h1_fre);
% Ign_deblurred1 = real(ifft2(Ign1_fre ./ h1_fre));

Ign_deblurred_inv_1 = inverseFiltering(Ign1, h1);           
        ISNR_inv_1 = insr(I, Ign1, Ign_deblurred_inv_1);
Ign_deblurred_inv_2 = inverseFiltering(Ign2, h2);           
        ISNR_inv_2 = insr(I, Ign2, Ign_deblurred_inv_2);

epsilon = 0.1;
Ign_deblurred_pseudo_1 = pseudoinverseFiltering(Ign1, h1, epsilon);
        ISNR_pseudo_1 = insr(I, Ign1, Ign_deblurred_pseudo_1);
Ign_deblurred_pseudo_2 = pseudoinverseFiltering(Ign2, h2, epsilon);
        ISNR_pseudo_2 = insr(I, Ign2, Ign_deblurred_pseudo_2);

figure();   
tiledlayout(2,3,'TileSpacing','compact');
nexttile, imshow(Ign1 / 255); title(["Origianl noisy image deblurred by a " + N1 + "\times" + N1 + " Gaussian filter (\sigma = " + sigma + ")", ...
    sprintf("BSNR is %.2f dB", BSNR1)]);
nexttile, imshow(Ign_deblurred_inv_1 / 255); title(["Inverse filtering " + ...
    sprintf("(ISNR is %.2f dB)", ISNR_inv_1), ...
    "Origianl image is deblurred by a " + N1 + "\times" + N1 + " Gaussian filter"]);
nexttile, imshow(Ign_deblurred_pseudo_1 / 255); title(["Pseudoinverse filtering (Threshold = " + epsilon + ") " + ...
    sprintf("(ISNR is %.2f dB)", ISNR_pseudo_1), ...
    "Origianl image is deblurred by a " + N1 + "\times" + N1 + " Gaussian filter"]);

nexttile, imshow(Ign2 / 255); title(["Origianl noisy image deblurred by a " + N2 + "\times" + N2 + " Gaussian filter (\sigma = " + sigma + ")", ...
    sprintf("BSNR is %.2f dB", BSNR2)]);
nexttile, imshow(Ign_deblurred_inv_2 / 255); title(["Inverse filtering "+ ...
    sprintf("(ISNR is %.2f dB)", ISNR_inv_2), ...
    "Origianl image is deblurred by a " + N2 + "\times" + N2 + " Gaussian filter"]);
nexttile, imshow(Ign_deblurred_pseudo_2 / 255); title(["Pseudoinverse filtering (Threshold = " + epsilon + ") "+ ...
    sprintf("(ISNR is %.2f dB)", ISNR_pseudo_2), ...
    "Origianl image is deblurred by a " + N2 + "\times" + N2 + " Gaussian filter"]);

%% Wiener Filter
img_deblurred_wnr_1 = wienerFiltering(Ign1, h1, noiseSigma);
        ISNR_wnr_1 = insr(I, Ign1, img_deblurred_wnr_1);
img_deblurred_wnr_2 = wienerFiltering(Ign2, h2, noiseSigma);
        ISNR_wnr_2 = insr(I, Ign2, img_deblurred_wnr_1);
        
figure();
tiledlayout(1,2,'TileSpacing','compact');
nexttile, imshow(img_deblurred_wnr_1/255);  title(["Inverse filtering "+ ...
    sprintf("(ISNR is %.2f dB)", ISNR_wnr_1), ...
    "Origianl image is deblurred by a " + N1 + "\times" + N1 + " Gaussian filter"]);
nexttile, imshow(img_deblurred_wnr_2/255);  title(["Inverse filtering "+ ...
    sprintf("(ISNR is %.2f dB)", ISNR_wnr_2), ...
    "Origianl image is deblurred by a " + N2 + "\times" + N2 + " Gaussian filter"]);


J = deconvwnr(Ign1,h1,numel(Ign1) * noiseSigma^2 ./ (abs(fft2(Ign1)).^2));
% J = deconvwnr(Ign1,h1, noiseSigma^2 / var(Ign1(:)));
figure();   imshow(J/255);
ISNR_wnr_1 = insr(I, Ign1, J);

J2 = deconvwnr(Ign2,h2,numel(Ign2) * noiseSigma^2 ./ (abs(fft2(Ign2)).^2));
% J = deconvwnr(Ign1,h1, noiseSigma^2 / var(Ign1(:)));
figure();   imshow(J2/255);
ISNR_wnr_2 = insr(I, Ign2, J2);


function BSNR = bsnr(z, noisePower)
    % BSNR = signal variance / noise variance
    % Large BSNR => Light noise effect
    % Small BSNR => Heavy noise effect
    BSNR = 10 * log10(var(z(:))/noisePower);
end

function ISNR = insr(f, y, f_hat)
    % f = original clean image
    % y = blurred and noisy image
    % f_hat = restored image
    numerator = sum(sum( (f-y).^2 ));
    denominator = sum(sum( (f-f_hat).^2 ));
    ISNR = 10 * log10(numerator/denominator);
end

function img_deblurred = wienerFiltering(img_blurred, h, noiseSigma)
% The noise variance should be known.
% This noise variance can be known according to the priori knowledge of the image acquisition process.
% Or, may be estimated from the local variance of a flat region of the observed image.
Snn = numel(img_blurred) * noiseSigma .^ 2;
img_fre = fft2(img_blurred);
Sff = (abs(img_fre)) .^ 2;
H_fre = fft2(h,size(img_blurred,1),size(img_blurred,2));
img_deblurred_fre = (Sff .* conj(H_fre) .*  img_fre) ./ (Sff .* (abs(H_fre)).^2 + Snn);
% nsr = Snn ./ Sff;
% img_deblurred_fre = (conj(H_fre) .*  img_fre) ./ ((abs(H_fre)).^2 + nsr);
img_deblurred = real(ifft2(img_deblurred_fre));
end

function img_deblurred = inverseFiltering(img_blurred, h)
    % img_blurred is the degraded image
    % h is the degradation model
    img_fre = fft2(img_blurred);
    h_fre = fft2(h,size(img_blurred,1),size(img_blurred,2));
    img_deblurred = real(ifft2(img_fre ./ h_fre));
end

function img_deblurred = pseudoinverseFiltering(img_blurred, h, epsilon)
    % img_blurred is the degraded image
    % h is the degradation model
    % epsilon is the threshold for pseudo-inverse filtering
    img_fre = fft2(img_blurred);
    h_fre = fft2(h,size(img_blurred,1),size(img_blurred,2));
    img_deblurred_fre = zeros(size(img_fre));
    nz = find(abs(h_fre) >= epsilon);
    img_deblurred_fre(nz) = img_fre(nz) ./ h_fre(nz);
    img_deblurred = real(ifft2(img_deblurred_fre));
end
