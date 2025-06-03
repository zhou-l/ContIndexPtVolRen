%% For image plane idx pt files, use columns 1 and 2
% M = load('ImgPlaneIdxPtDebug_load2.txt');
%X = M(:,1:2);

%% For idxPtDebug compute files, use subdimensions
% M = load('idxPtDebug_compute045.txt');

%  M = load('idxPtDebug_computeTooth06.txt');
% M = load('idxPtDebug_computeToothFull.txt');
%  M = load('idxPtDebug_computeTooth045.txt');
% M = load('idxPtDebug_computeBraTs045.txt');
M = load('idxPtDebug_computeIsabel045.txt');
%% Resolution of the density image
H = 1000;
W = 1000;
% X = M(:,2:3); % idx pts of subdim 0
X = M(:,6:7); % idx pts of subdim 1
% X = M(:,10:11); % idx pts of subdim 2
minX = min(X);
maxX = max(X);
disp(minX);
disp(maxX);
nX = (X - minX)./(maxX - minX);
nX = nX .* [H-1 W-1];
disp('# original samples:')
disp(size(nX,1));
% get the input image
I = zeros(H,W);
orgiPos = double(int32(nX + [0.5 0.5]));
orgiPos = unique(orgiPos, 'rows');
disp('# unique samples');
disp(size(orgiPos,1));
orgiPosNoNAN = orgiPos;
cnt = 1;
%  remove nan values
for i = 1:size(orgiPos,1)
    if sum(isnan(orgiPos(i,:))) > 0
        continue;
    end
    orgiPosNoNAN(cnt,:) = orgiPos(i,:);
    cnt = cnt + 1;
end
orgiPosNoNAN(cnt:end,:) = [];
orgiPos = orgiPosNoNAN;

% create the density image from original samples
for i = 1:size(nX,1)
    if sum(isnan(nX(i,:))) > 0
        continue;
    end
    pos = int32(nX(i,:)+ [0.5 0.5]);
%     pos = int32(orgiPos(i,:));
    I(pos(1),pos(2)) = I(pos(1),pos(2)) + 1;
end

% IX = I';
% IX = flip(IX,1);
% gaussSig = 2;
% IXg = imgaussfilt(IX, gaussSig);
% figure, 
% subplot(1,2,1); imshow(log(IX+1),[]);
% subplot(1,2,2); imshow(log(IXg+1),[]);

% create importance samples
% N = 5 * 10^5;
N = round(0.5 * size(orgiPos, 1));
disp('# important samples:');
disp(N);
target = zeros(N,1);
samplesR = randi(H, N,1); %rand(N,1);%(rand(N,1)-0.5)*5; 
samplesC = randi(W, N,1); %rand(N,1); %rand(N, 1) * 10; 

[G, Gdir] = imgradient(I);
for i = 1:N
    target(i,:) = G(samplesR(i),samplesC(i));
end
%try the gradient magnitude as the target
proposal = 1/N;
w = target ./ proposal;
sample_ind = randsample([1:N],N,true,w);
nImportant = [samplesR(sample_ind) samplesC(sample_ind)];
dispI = flip(I',1);
figure, imshow(log(dispI+1),[]);hold on;
% plot(nImportant(:,1),H - nImportant(:,2), '.');
s = scatter(nImportant(:,1),nImportant(:,2), 'filled');
hold off;
alpha(s,.5);

% The list of kernel types
% Kf = {@gaussKernel, @epanechnikovKernel, @expKernel, @cosKernel, @triangleKernel, @avgKernel};
kernelList = [2];

%% In general, epanechikov and cosine kernels are superior to others.
% test with original unique samples
useDynBandwidth = true;
Iunique = myReconstructKDE(orgiPos, I, H, W, kernelList, useDynBandwidth);
% write out hdr file
Ihdr = cat(3,Iunique,Iunique,Iunique);
% hdrwrite(Ihdr,'BraTs_flat1_subdim0.hdr');
% hdrwrite(Ihdr,'BraTs_flat1_subdim1.hdr');
% hdrwrite(Ihdr,'tooth_flat1_subdim0.hdr');
 hdrwrite(Ihdr,'flat1_subdim1.hdr');
%% test with important sampling

% Iimport = myReconstructKDE(nImportant, I, H, W, kernelList);

