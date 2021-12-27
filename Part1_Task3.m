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
colormap(jet(64));  colorbar;

%% ordered Hadamard matrix
Hn_ordered = orderedHadamard(N);
figure();   imagesc(Hn_ordered/N);  title(N +" \times " + N + " Hadamard matrix (Ordered)");
colormap(gray);

img_H_ordered = Hn_ordered * img * Hn_ordered' / N;

figure(); imagesc(log(abs(img_H_ordered))); title("The logarithmic magnitude of the 2D Ordered DHT of Cameraman");
colormap(jet(64));  colorbar;

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

