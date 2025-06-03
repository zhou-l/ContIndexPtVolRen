function pcaVols = testOctreeEigReconstruct(zfileList, zvalDiff)
%% load data volumes
% TODO: need to test on downsampled volumes
if nargin < 1
    fileList = {'./Toothuc0.20.nrrd'; './Toothucderiv10.20.nrrd'; './Toothucderiv20.20.nrrd'};
else
    fileList = zfileList;
    if nargin < 2
        valDiff = 0.01;
    else
        valDiff = zvalDiff;
    end
end
% fileList = {'../data/Toothuc.nrrd'; '../data/Toothucderiv1.nrrd'; '../data/Toothucderiv2.nrrd'};
%  fileList = {'../data/Pf25cpUs.nhdr'; '../data/TCf25cpUs.nhdr'; '../data/QVAPORf25cpUs.nhdr'; '../data/PRECIPf25cpUs.nhdr'};
% myEigInfoReconstruct(fileList)
nVols = length(fileList);

% fn = sprintf('eigVolDS%i_%0.2f.nrrd', 0, sampleRatePerAxis);
[Vlist, nVols, volInfoList] = loadVolumes(fileList);
% build an gradient magnitude volume combining all attributes
volSize = size(Vlist{1});
VgradSum = zeros(volSize(1),volSize(2),volSize(3)); 
% VgradList = Vlist;
for i = 1:nVols
    [Gmag,Gaz,Gel] = imgradient3(Vlist{i});
%     VgradList{i} = Gmag;
    gmin = min(min(min(Gmag)));
    gmax = max(max(max(Gmag)));
 
    Gmag = (Gmag - gmin) ./ (gmax-gmin);
    VgradSum = VgradSum + Gmag;
end
clear Gmag Gaz Gel;
% convolve with a filter with a size of 2*R+1
R = 3;
Kernel = zeros(2*R+1,2*R+1,2*R+1);
Kctr = [R+1 R+1 R+1];
cnt = 0;
for z = 1:2*R+1
    for r = 1:2*R+1
        for c = 1:2*R+1
            v = [r c z] - Kctr;
            if norm(v) <= R
                Kernel(r,c,z) = 1;
                cnt = cnt + 1;
            end
        end
    end
end
if cnt == 0
    Kernel = 0;
else
    Kernel = Kernel ./ cnt;
end
% try with convolution
convMImg = imfilter(VgradSum,Kernel);
VgradSum=convMImg;
figure, 
histogram(convMImg);
figure,
montage(reshape(VgradSum,volSize(1),volSize(2),1,volSize(3)),'DisplayRange',[])
title('Filtered Gradient magnitude')
%% test octree
tic
disp('Building the octree.');

% create padded volumes

[odimX, odimY, odimZ] = size(Vlist{1});
nextPow2 = 2^ceil(log2( max(size(Vlist{1}))));
padDimFirstHalf = zeros(3,1);
padDimSecondHalf = zeros(3,1);
VpadList = cell(nVols,1);
V=Vlist{1};
for i = 1:3
    if mod(size(V,i),2) == 0
        padDimFirstHalf(i) = int16((nextPow2 - size(VgradSum,i))/2);
        padDimSecondHalf(i) = int16((nextPow2 - size(VgradSum,i))/2);
    else
        
        padDimFirstHalf(i) = int16(nextPow2/2 - floor(size(VgradSum,i)/2));
        padDimSecondHalf(i) = int16(nextPow2 - size(VgradSum,i) - padDimFirstHalf(i));
    end
end

% pad in 3D
VP = zeros(nextPow2,nextPow2,nextPow2);
VP(1+padDimFirstHalf(1):padDimFirstHalf(1)+odimX,...
    1+padDimFirstHalf(2):padDimFirstHalf(2)+odimY,...
    1+padDimFirstHalf(3):padDimFirstHalf(3)+odimZ) = VgradSum;
% V = VP;
figure,
montage(reshape(VP,nextPow2,nextPow2,1,nextPow2),'DisplayRange',[])
dimX = nextPow2; dimY = nextPow2; dimZ = nextPow2;
% create the octree
% Value diff is the only parameter of octree!
% valDiff = 0.01;
OT = VolOctree(VP, [dimX dimY dimZ], 'binValDiff', valDiff);
toc
% show the octree. NOTE: can be quite slow!
% figure, OT.plot();
%%
tic
disp('creating padded volumes');
for nv = 1:nVols
    V = Vlist{nv}; %VgradSum;
    maxVal = max(max(max(V)));
%     disp(maxVal);
    [dimX, dimY, dimZ] = size(V);
    maxSize = max(size(V));
    
    nextPow2 = 2^ceil(log2(maxSize));
    for i = 1:3
        if mod(size(V,i),2) == 0
            padDimFirstHalf(i) = (nextPow2 - size(V,i))/2;
            padDimSecondHalf(i) = (nextPow2 - size(V,i))/2;
        else
            
            padDimFirstHalf(i) = nextPow2/2 - floor(size(V,i)/2);
            padDimSecondHalf(i) = nextPow2 - size(V,i) - padDimFirstHalf(i);
        end
    end
    % pad in 3D
    VP = zeros(nextPow2,nextPow2,nextPow2);
    VP(1+padDimFirstHalf(1):padDimFirstHalf(1)+dimX,...
        1+padDimFirstHalf(2):padDimFirstHalf(2)+dimY,...
        1+padDimFirstHalf(3):padDimFirstHalf(3)+dimZ) = V;
    VpadList{nv} = VP;
    dimX = nextPow2; dimY = nextPow2; dimZ = nextPow2;
    
    figure,montage(reshape(VP,dimX,dimY,dimZ),'DisplayRange',[]);

%     valDiff = 100;
    % create the octree
%     valDiff = 0.2;
%     OT = VolOctree(V, [dimX dimY dimZ], 'binValDiff', valDiff);
     % plot the octree
%      OT.plot();
end
% figure, OT.plot(1,1,6);
toc
%% Compute PCA with octree
clear Vlist;
orgVolMinPos = 1+ padDimFirstHalf;
orgVolMaxPos = [odimX odimY odimZ] + padDimFirstHalf;
minBlockSize = inf;
maxBlockSize = -inf;
corsLvl = inf;
fineLvl = -inf;
% pca list where each item is an nVols*nVols 
% pcaList = zeros(OT.BinCount,(nVols+1)*nVols + 1); % use the last column to record 
% pca list records only the Principal component!
pcaList = zeros(OT.BinCount, 8*nVols+1);  % pcaList: [cnt, #posCalculated|E1_p0|E1_p1|...|E1_p7|]
otLeafNodeIdList = zeros(OT.BinCount,1);
cnt = 1;
disp('Calculating local PCA for octree nodes.');
disp("# of all nodes");
disp(OT.BinCount)
tic
inNodeValThres = 0.01;
for i = 1:OT.BinCount
    hasChildren = find(OT.BinParents == i,1);
    isEightValsNeeded = false;
    % NOTE: DEBUG!!!
    % debug for nodes with side of 16
    %     if OT.BinDepths(i)~=2
    %         continue;
    %     end
    % end debug
    if ~isempty(hasChildren)
        continue;
    end
    % is leaf node
    bounds = OT.BinBoundaries(i,:);
    %         disp(OT.BinDepths(i));
    boundSize = diff(bounds([1:3;4:6])) + [1 1 1];

    % check the size of the block
    if boundSize(1) >= 2 && OT.BinValMinMax(2)-OT.BinValMinMax(1) >= inNodeValThres
        isEightValsNeeded = true;
    end
    
    % check if outside of the original volume?
    % TODO: We need to do an AABB-AABB test!!!
    if ~aabbIntersect(bounds(1:3), bounds(4:6), orgVolMinPos, orgVolMaxPos)
        continue;
    end
    
    if isEightValsNeeded
        % use eight pos
        %    3-E2-2
        %   /|   /|
        %  0-E1-1 |
        %  | 4-E3-5
        %  |/   |/
        %  7-E4-6
        pcaList(cnt,1) = 8;
        pos = bounds(1:3); % 0
        [eigVecList, eigValList] = localPCAspaceDom(pos, VpadList, R); % spatial domain local PCA
        pcaList(cnt,1 + 1:1+nVols) = eigVecList(:,1);
        
        pos = bounds(1:3)+[boundSize(1),0,0]; % 1
        [eigVecList, eigValList] = localPCAspaceDom(pos, VpadList, R); % spatial domain local PCA
        pcaList(cnt,1+nVols+1:1 + 2*nVols) = eigVecList(:,1);
        
        pos = bounds(1:3)+[boundSize(1),boundSize(2),0];% 2
        [eigVecList, eigValList] = localPCAspaceDom(pos, VpadList, R); % spatial domain local PCA
        pcaList(cnt,1 + 2*nVols+1:1+3*nVols) = eigVecList(:,1);
        
        pos = bounds(1:3)+[0,boundSize(2),0];% 3
        [eigVecList, eigValList] = localPCAspaceDom(pos, VpadList, R); % spatial domain local PCA
        pcaList(cnt,1 + 3*nVols+1:1+4*nVols) = eigVecList(:,1);
        
        pos = bounds(1:3)+[0,boundSize(2),boundSize(3)];% 4
        [eigVecList, eigValList] = localPCAspaceDom(pos, VpadList, R); % spatial domain local PCA
        pcaList(cnt,1 + 4*nVols+1:1+5*nVols) = eigVecList(:,1);
        
        pos = bounds(1:3)+[boundSize(1),boundSize(2),boundSize(3)];% 5
        [eigVecList, eigValList] = localPCAspaceDom(pos, VpadList, R); % spatial domain local PCA
        pcaList(cnt,1 + 5*nVols+1:1+6*nVols) = eigVecList(:,1);
        
        pos = bounds(1:3)+[boundSize(1),0,boundSize(3)];% 6
        [eigVecList, eigValList] = localPCAspaceDom(pos, VpadList, R); % spatial domain local PCA
        pcaList(cnt,1 + 6*nVols+1:1+7*nVols) = eigVecList(:,1);
        
        pos = bounds(1:3)+[0,0,boundSize(3)];% 7
        [eigVecList, eigValList] = localPCAspaceDom(pos, VpadList, R); % spatial domain local PCA
        pcaList(cnt,1 + 7*nVols+1:1+8*nVols) = eigVecList(:,1);
    else
        % Only one pos at the center!
        pos = floor(0.5 .* (bounds([4:6]) + bounds([1:3]))); %#ok<NBRAK>
        [eigVecList, eigValList] = localPCAspaceDom(pos, VpadList, R); % spatial domain local PCA
        %         record PCA result
%         pcaList(cnt,1:nVols*nVols) = reshape(eigVecList,[1,nVols*nVols]);
%         pcaList(cnt,nVols*nVols+1:end) = reshape(eigValList,[1,nVols]);
        allEigVecs = reshape(eigVecList,[1,nVols*nVols]);
        pcaList(cnt,1) = 1;
        pcaList(cnt,1+1:1+nVols) = allEigVecs(1:nVols);
%         pcaList(cnt,nVols*nVols+1:end) = reshape(eigValList,[1,nVols]);
%         pcaList(cnt, nVols+1:end) = []; % remove other slots
    end
%     %         record PCA result
%     pcaList(cnt,1:nVols*nVols) = reshape(eigVecList,[1,nVols*nVols]);
%     pcaList(cnt,nVols*nVols+1:end) = reshape(eigValList,[1,nVols]);
    
    otLeafNodeIdList(cnt,:) = i;
    cnt = cnt + 1;
    %         minBlockSize = min(minBlockSize, boundSize(1));
    %         maxBlockSize = max(maxBlockSize, boundSize(1));
    %
    %         corsLvl = min(corsLvl, OT.BinDepths(i));
    %         fineLvl = max(fineLvl, OT.BinDepths(i));
end
disp("# of leaf nodes");
disp(cnt)
toc
pcaList(cnt+1:end,:) = [];
otLeafNodeIdList(cnt+1:end,:) = [];
%% reconstruction pca volume from octree
pcaVols = cell(nVols,1);
pcaPadVols=cell(nVols,1);
for i = 1:nVols
    pcaVols{i} = zeros(odimX,odimY,odimZ);
    pcaPadVols{i} = zeros(dimX,dimY,dimZ);
end
minLevel=1;
maxLevel=inf;
for i = 1:length(otLeafNodeIdList)
    nodeId=otLeafNodeIdList(i);
    if nodeId<1
        continue;
    end
    %DEBUG: draw nodes that higher than the max level
    if OT.BinDepths(nodeId) < minLevel || OT.BinDepths(nodeId) > maxLevel
        continue;
    end
        
    bounds = OT.BinBoundaries(nodeId,:);
    boundSize = diff(bounds([1:3;4:6])) + [1 1 1];
    for bz=bounds(3):bounds(6)
        for by = bounds(2):bounds(5)
            for bx = bounds(1):bounds(4)
                pos = [bx by bz];
                opos = pos - padDimFirstHalf';
                 if  opos(1) < 1 || opos(2) < 1 || opos(3) < 1 ||...
                     opos(1) > odimX || opos(2) > odimY || opos(3) > odimZ
    % outside of the data 
                    continue;   
                 end
                 % get only the primary Eigenvecctor info
%                  pcaInfo = pcaList(i,1:nVols);
                 pcaInfo = pcaList(i,2:1+nVols);
                 % set pca info for the pcaVols 
                 % Should use 8 corners of the octree and interpolate!
%                  eigInterp3D();

                 for nv = 1:nVols
                     pcaPadVols{nv}(pos(1),pos(2),pos(3)) = pcaInfo(nv);
                     pcaVols{nv}(opos(1),opos(2),opos(3)) = pcaInfo(nv);%pcaInfo((nv-1)*nVols+1:(nv)*nVols);
                 end
            end
        end
    end
end

% write out pca vols
for i = 1:nVols
    [filepath, name, ext] = fileparts(fileList{i});
    fileName = sprintf('%s_Oct_eigVolReconst%i.nrrd', name, i);
    ok = nrrdWriter(fileName, single(pcaVols{i}), [1 1 1], [0 0 0], 'raw');
    fileName = sprintf('%s_PadOct_eigVolReconst%i.nrrd', name, i);
    ok = nrrdWriter(fileName, single(pcaPadVols{i}), [1 1 1], [0 0 0], 'raw');
end
end

