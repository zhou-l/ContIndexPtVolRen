function Id = myReconstructKDE(nSamples, I, H, W, kernelList, useDynBandwidth)
% The list of kernel types
% Kf = {@gaussKernel, @epanechnikovKernel, @expKernel, @cosKernel, @triangleKernel, @avgKernel};

if nargin < 5
     zkernelList = [1,2,3,4];
else
   zkernelList = kernelList;
end

KinSearch = 50;
rmax = 50;
h = 9;
lookupPos = zeros(H*W,2);
cnt = 1;
for r = 1:H
    for c = 1:W
        lookupPos(cnt,:) = [r c];
        cnt = cnt + 1;
        %         I(r,c) = myLocalKDE(nX, [r,c], nmax, rmax, h);
    end
end
[NIdx,D] = knnsearch(nSamples, lookupPos, 'K', KinSearch);
IX  = zeros(H,W);
%% The choice of nmax affects the result!!!
% parameter of the maximum density included in a kernel
nmax = 6000;
% useDynBandwidth = false;
for nk = 1:length(zkernelList)
    for i = 1:size(lookupPos,1)
        pos = lookupPos(i,:);
        nn = nSamples(NIdx(i,:),:);
        %     I(pos(1),pos(2)) = myLocalKDEwNN(pos, nn, D(i,:), rmax, h);
        IX(pos(1),pos(2)) = myLocalKDEwNN_dataVal(pos, nn, D(i,:), nmax, rmax, h, I, zkernelList(nk),useDynBandwidth);
    end
    figure; title('reconstruction result');
    Id= IX';
    Id = flip(Id,1);
    
    subplot(1,2,1); imshow(log(Id+1),[]); title('recontruct');
    Iorg = I';Iorg = flip(Iorg,1);
    subplot(1,2,2); imshow(log(Iorg+1), []);title('original');
    
end

% imwrite(Id2,'BraTs_flat1_subdim2_imp.png');
end