%% For image plane idx pt files, use columns 1 and 2
% M = load('ImgPlaneIdxPtDebug_load2.txt');
%X = M(:,1:2);
inputName = 'idxPtDebug_computeIsabel045.txt';

pat = 'compute';
k = strfind(inputName, pat);
if isempty(k)
    keyName = 'contIdxPtImg';
else
    [filepath, fileName, ext] = fileparts(inputName);
    keyName = fileName(1,k(1)+length(pat):end);
end

%% For idxPtDebug compute files, use subdimensions
% M = load('idxPtDebug_compute045.txt');
% M = load('idxPtDebug_computeBraTs.txt');
M = load(inputName);
%% Resolution of the density image
H = 1000;
W = 1000;
numSubDim = int32((size(M,2)/4));
% process each subdimension
for d = 1:numSubDim
    startCol = (d-1)*4+2;
    endCol = startCol+1;
    % X = M(:,2:3); % idx pts of subdim 1
    % X = M(:,6:7); % idx pts of subdim 2
    X = M(:,startCol:endCol); % idx pts of subdim 2
    [rows, cols] = find(~isnan(X));
    X = X(unique(rows),:);
    [rows, cols] = find(isinf(X));
    rows = unique(rows);
    X(rows,:) = [];

    minX = min(X);
    maxX = max(X);
    disp(minX);
    disp(maxX);
    nX = (X - minX)./(maxX - minX);
    % test kde density
    % figure;ksdensity(nX);
    
    nX = nX .* [H W];
    I = zeros(H,W);
    nmax = 50;
    rmax = 40;
    h = 3.3;
    lookupPos = zeros(H*W,2);
    cnt = 1;
    for r = 1:H
        for c = 1:W
            lookupPos(cnt,:) = [r c];
            cnt = cnt + 1;
            %         I(r,c) = myLocalKDE(nX, [r,c], nmax, rmax, h);
        end
    end
    [NIdx,D] = knnsearch(nX, lookupPos, 'K', nmax);
    for i = 1:size(lookupPos,1)
        pos = lookupPos(i,:);
        nn = nX(NIdx(i,:),:);
        I(pos(1),pos(2)) = myLocalKDEwNN(pos, nn, D(i,:), rmax, h);
    end
    figure;
    Id= I';
    Id = flip(Id,1);
    imshow(Id,[]); hold on;
    % plot(nX(:,2),nX(:,1),'.');
    % imwrite(Id,'flat1_subdim1.png');
    % imwrite(Id,'BraTs_flat1_subdim1.png');
    outFileName = sprintf('%s_subdim%i.png',keyName, d);
    imwrite(Id,outFileName);
    hdrFileName = sprintf('%s_flat1_subdim%d.hdr', keyName, d-1);
    I3 = cat(3, Id,Id,Id);
    hdrwrite(I3, hdrFileName);
end
% write out continuous idx points images

