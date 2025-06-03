%% For image plane idx pt files, use columns 1 and 2
%% The function for converting discrete indexed points to continuous indexed
%% points with dynamic KDE
function buildContIdxPt_dynamicKDE(scrPointInfoInputName)

% M = load('ImgPlaneIdxPtDebug_load2.txt');
%X = M(:,1:2);
%  M = load('idxPtDebug_computeTooth06.txt');
% M = load('idxPtDebug_computeBraTs.txt');
% inputName = 'idxPtDebug_computeTooth06.txt';
% inputName = 'idxPtDebug_computeIsabel045.txt';
% inputName = 'idxPtDebug_computeSynthNew.txt';
% inputName = 'idxPtDebug_computeIsabel045.txt';
% inputName = 'idxPtDebug_computeSynthNew.txt';
%  inputName = 'idxPtDebug_computeTooth.txt';
% inputName = './DTI/idxPtDebug_computeGK2.txt';
% inputName = './DTI/idxPtDebug_compute_gk2clcpcs.txt';
%  inputName = 'idxPtDebug_computeFaParaPeri.txt';
% inputName = 'idxPtDebug_compute_p2FAParaPerp.txt';
% inputName = 'idxPtDebug_compute_p1RAFATr.txt';
%  inputName = 'idxPtDebug_compute_p2_RAFATr.txt';
% inputName = 'idxPtDebug_compute_p2LinPlaSph.txt';
% inputName = 'idxPtDebug_compute_p1LinPlaSph.txt';
% inputName = 'idxPtDebug_compute_p2TrFaLin.txt';
% inputName = 'idxPtDebug_compute_p2PerpFATr.txt';
% inputName = 'idxPtDebug_computeDTI2_FaTrPerp.txt';
inputName = scrPointInfoInputName; %'idxPtDebug_computeBratsT1_gradmag.txt';
%% 2-flats
% inputName = 'twoFlats_idxPtDebug_computeHalfCyl3d.txt';
% inputName = 'twoFlats_idxPtDebug_computeCtbl3d09.txt';
% inputName = 'twoFlats_idxPtDebug_computeTornado10.txt';
% inputName = './DTI/twoFlats_idxPtDebug_computeGK2.txt';
% inputName = './DTI/twoFlats_idxPtDebug_compute_gk2clcpcs.txt';
% inputName = 'idxPtDebug_computeHalfCyl3d.txt';
% inputName = '../data/BraTS/indexedPtVols/idxPtDebug_computeBraTsT1.txt';
% inputName = 'idxPtDebug_computeIsabel045.txt';
% inputName = 'idxPtDebug_computeSynthNew.txt';
% inputName = 'idxPtDebug_computeBraTs045.txt';
% inputName = 'idxPtDebug_computeBraTs.txt';
% inputName = 'twoFlats_idxPtDebug_computeSynthNew.txt';
% inputName = 'twoFlats_idxPtDebug_computeTooth.txt';
% inputName = 'twoFlats_idxPtDebug_computeDTI12D.txt';
% inputName = 'twoFlats_idxPtDebug_computeFaParaPeri.txt';
% inputName = 'twoFlats_idxPtDebug_compute_p1FAParaPerp.txt';
%  inputName = 'twoFlats_idxPtDebug_compute_p1RAFATr.txt';
%  inputName = 'twoFlats_idxPtDebug_compute_p2_RAFATr.txt';
% inputName = 'twoFlats_idxPtDebug_compute_p2LinPlaSph.txt';
% inputName = 'twoFlats_idxPtDebug_compute_p1LinPlaSph.txt';
% inputName = 'twoFlats_idxPtDebug_compute_p2TrFaLin.txt';
% inputName = 'twoFlats_idxPtDebug_compute_p2PerpFATr.txt';
% inputName = 'twoFlats_idxPtDebug_compute_dtiEyeP2.txt';
% inputName = 'twoFlats_idxPtDebug_compute_synthetic3D.txt';
% inputName = 'twoFlats_idxPtDebug_compute_FaTrPerp.txt';
% inputName = 'twoFlats_idxPtDebug_computeIsabel3Vars.txt';

[filepath, fileName, ext] = fileparts(inputName);
pat = 'compute';
k = strfind(fileName, pat);
if isempty(k)
    keyName = 'contIdxPtImg';
else
    keyName = fileName(1,k(1)+length(pat):end);
end
twof = 'twoFlats';
kf = strfind(fileName, twof);
if isempty(kf)
    typeN = 'flat1';
else
    typeN = 'flat2';
end
    
%% For idxPtDebug compute files, use subdimensions
M = load(inputName);
%% Resolution of the density image
H = 1000;
W = 1000;
numSubDim = int32((size(M,2)/4));
startDim = 1;
for d = startDim:numSubDim
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
    % NOTE: We have to match the value ranges with the upstream indexed points compuation!!!
    % Especially, for Y, as X typically uses the full ranges
% _minPCInselbergPt = FLOATVECTOR2(-0.15f, -1.25f);
% _maxPCInselbergPt = FLOATVECTOR2(1.15f, 1.25f);
    minX = min(X); 
    
    minX(2) = -1.25;
    maxX = max(X); 
    maxX(2) = 1.25;

    disp(minX);
    disp(maxX);
    
%     % remove the (0,0) points for d > 1
%     if d > 1
%         idx = find(X(:,1) == 0);
%         X(idx,:) = [];
%     end
%     minX(1) = min(X(:,1)); maxX(1) = max(X(:,1));
%     disp(minX); disp(maxX);
    
    nX = (X - minX)./(maxX - minX);
    nX = nX .* [W-1 H-1];
    disp('# original samples:')
    disp(size(nX,1));
    % get the input image
    I = zeros(W,H);
    orgiPos = double(int32(nX + [0.5 0.5]));
    orgiPos = unique(orgiPos, 'rows');
    orgiPos(:,1) = min(W,orgiPos(:,1));
    orgiPos(:,2) = min(H,orgiPos(:,2));
    disp('# unique samples');
    disp(size(orgiPos,1));    
    % create the density image from original samples
    for i = 1:size(nX,1)
%         if sum(isnan(nX(i,:))) > 0
%             continue;
%         end
        pos = int32(nX(i,:)+ [0.5 0.5]);
        pos(1) = max(1,min(W,pos(1)));
        pos(2) = max(1,min(H,pos(2)));
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
    s = scatter(nImportant(:,1),H-nImportant(:,2), 'filled');
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
    hdrFileName = sprintf('%s_%s_subdim%d.hdr', keyName, typeN, d-1);
    hdrwrite(Ihdr,hdrFileName);
end
end
%% test with important sampling

% Iimport = myReconstructKDE(nImportant, I, H, W, kernelList);

