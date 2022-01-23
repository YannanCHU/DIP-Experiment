% Hadamard Transform
clc;
clear;
close all;

img = imread("Cameraman.tif");
figure();
image(img);

N = size(img,1);    % Length of Walsh (Hadamard) matrix
Hn = hadamard(N);
figure();   imagesc(Hn/N);  title(N +" \times " + N + " Hadamard matrix (non-ordered)");
colormap(gray);
% Hn = 1;
% for i = 1:log2(N)
%     Hn = [Hn Hn; Hn -Hn];
% end

img = double(img);
% Two-Dimensional Discrete Hadamard Transform
img_H = Hn * img * Hn' / N;

figure(); imagesc(log(abs(img_H))); title("The logarithmic magnitude of the 2D non-ordered DHT of Cameraman");
colormap(jet(64));  colorbar;   daspect([1,1,1]);

img_H_recover = Hn' * img_H * Hn / N;
figure(); imagesc(img_H_recover); title("The logarithmic magnitude of the 2D non-ordered DHT of Cameraman");
colorbar;   daspect([1,1,1]);


% cumulative energy sequences
Hadamard_unord_ces = zeros(1, size(img_H,1));
for i = 1:1:size(img_H,1)
    Hadamard_unord_ces(i) = sum(sum(abs(img_H(1:i,1:i)).^2));
end
Hadamard_unord_ces = Hadamard_unord_ces / sum(sum(abs(img_H).^2));


%% ordered Hadamard matrix
Hn_ordered = orderedHadamard(N);
figure();   imagesc(Hn_ordered/N);  title(N +" \times " + N + " Hadamard matrix (Ordered)");
colormap(gray); daspect([1,1,1]);

img_H_ordered = Hn_ordered * img * Hn_ordered' / N;

figure(); imagesc(log(abs(img_H_ordered))); title("The logarithmic magnitude of the 2D Ordered DHT of Cameraman");
colormap(jet(64));  colorbar; daspect([1,1,1]);

% cumulative energy sequences
Hadamard_ord_ces = zeros(1, size(img_H_ordered,1));
for i = 1:1:size(img_H_ordered,1)
    Hadamard_ord_ces(i) = sum(sum(abs(img_H_ordered(1:i,1:i)).^2));
end
Hadamard_ord_ces = Hadamard_ord_ces / sum(sum(abs(img_H_ordered).^2));

%% cumulative energy sequences comparison
img_dct = dct2(img);
% cumulative energy sequences
img_dct_ces = zeros(1, size(img_dct,1));
for i = 1:1:size(img_dct,1)
    img_dct_ces(i) = sum(sum(abs(img_dct(1:i,1:i)).^2));
end
img_dct_ces = img_dct_ces / sum(sum(abs(img_dct).^2));


figure();   plot(1:size(img_H,1), Hadamard_unord_ces, 1:size(img_H,1), Hadamard_ord_ces, 1:size(img_H,1), img_dct_ces, 'LineWidth', 2);  
title("The cumulative energy sequences of three transforms", 'fontSize', 24)
legend("non-ordered DHT", "Ordered DHT", "DCT", 'fontSize', 12);
ylabel("Fraction of entire image energy", 'fontSize', 12);
xlabel("Lenth of patch", 'fontSize', 12)
xlim([1, 256]);

function orderedHadamardMatrix = orderedHadamard(N)
    hadamardMatrix = hadamard(N);

    HadIdx = 0:N-1;                          % Hadamard index
    M = log2(N)+1;                           % Number of bits to represent the index

    binHadIdx = fliplr(dec2bin(HadIdx,M))-'0'; % Bit reversing of the binary index
    binSeqIdx = zeros(N,M-1);                  % Pre-allocate memory
    for k = M:-1:2
        % Binary sequency index
        binSeqIdx(:,k) = xor(binHadIdx(:,k),binHadIdx(:,k-1));
    end
    SeqIdx = binSeqIdx*pow2((M-1:-1:0)');    % Binary to integer sequency index
    orderedHadamardMatrix = hadamardMatrix(SeqIdx+1,:); % 1-based indexing
end

