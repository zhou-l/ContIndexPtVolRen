useCase = 'isabel';

if strcmp(useCase, 'tooth') == 1
    % Uncomment to use the tooth data
    I0 = imread('VolData0-80.ppm');
    I1 = imread('VolData1-80.ppm');
    I2 = imread('VolData2-80.ppm');
    I0 = double(rgb2gray(I0));
    I1 = double(rgb2gray(I1));
    I2 = double(rgb2gray(I2));
elseif strcmp(useCase, 'isabel') == 1
    % Uncomment  to use the isabel data
    I0 = imread('isabel/VolData0-45.ppm');
    I1 = imread('isabel/VolData1-45.ppm');
    I2 = imread('isabel/VolData2-45.ppm');
    I0 = double(rgb2gray(I0));
    I1 = double(rgb2gray(I1));
    I2 = double(rgb2gray(I2));
else
    % Uncomment to use the testdata
    Imv = buildTestData();
    I0 = Imv{1}; I1 = Imv{2}; I2 = Imv{3};
end
% Use gradient as the target pdf
gI0 = imgradient(I0, 'sobel');
gI1 = imgradient(I1, 'sobel');
gI2 = imgradient(I2, 'sobel');
figure;

subplot(2, 3, 1); imshow(I0, []);title('Scalar Data');
subplot(2, 3, 2); imshow(I1, []);
subplot(2, 3, 3); imshow(I2, []);


subplot(2, 3, 4); imshow(gI0, []);title('Gradient Magnitude Data');
subplot(2, 3, 5); imshow(gI1, []);
subplot(2, 3, 6); imshow(gI2, []);

N = 10^5;
% Multiple importance sampling, using the maximum?
heursiticType = 0;

multiImg = max(max(gI0,gI1), gI2);
% Or trying the power heuristic?
% inline Float PowerHeuristic(int nf, Float fPdf, int ng, Float gPdf) {
%     Float f = nf * fPdf, g = ng * gPdf;
%     return (f * f) / (f * f + g * g);
% }
sumImg = gI0 + gI1 + gI2;

maxR = size(multiImg, 1) ;
maxC = size(multiImg, 2) ;


figure,imshow(multiImg, []);
mean_samples = randi(maxR, N,1); %rand(N,1);%(rand(N,1)-0.5)*5; 
sigma_samples = randi(maxC,N,1); %rand(N,1); %rand(N, 1) * 10; 
proposal = 1/N;
% idxR = int32(maxR .* mean_samples) + 1;
% idxC = int32(maxC .* sigma_samples)+ 1;
idxR = mean_samples;
idxC = sigma_samples;

targetPower = zeros(N,1);
targetSum = zeros(N,1);
targetMax = zeros(N,1);

nList = [N; N; N];
% m is the dimension of the input data
m = size(nList, 1);

% test neighborhood computation
scalarImgs = cell(3,1);
scalarImgs{1,:} = I0;
scalarImgs{2,:} = I1;
scalarImgs{3,:} = I2;
 % store the top 3 eigenvectors of m dimensions (3*m layers)
pca_map_impSamp = zeros(maxR, maxC,3 * m);
pca_map_gold = zeros(maxR, maxC, 3*m );
R = 3;

% store the associated pca result
nEV = 3; % number of eigenvectors

% xrange = [240, 269];
% yrange = [220, 249];
xrange = [448 455];
yrange = [412 420];
pca_focusMap = zeros(yrange(2)-yrange(1),xrange(2)-xrange(1),3 * m);
% for r = yrange(1,1):yrange(1,2)%1:maxR
for r = 1:maxR
    %     for c = xrange(1,1):xrange(1,2)%1:maxC
    for c = 1:maxC
        pos = [r,c];
        pca_coeff_gold = findSamplesInBallGold(pos, scalarImgs, R);
        
        offPos = pos - [yrange(1,1),xrange(1,1)] + [1 1];
        for k = 1:nEV
            % store the eigenvector information
            pca_map_gold(pos(1),pos(2),(k-1)*m+1:k * m) = pca_coeff_gold(:,k);
%             pca_map_gold(pos(1),pos(2),m+1:2*m) = pca_coeff_gold(:,2);
%             pca_map_gold(pos(1),pos(2),2*m+1:3*m) = pca_coeff_gold(:,3);

            if offPos(1) > 0 && offPos(2) > 0 &&... 
                offPos(1) <=size(pca_focusMap,1) && offPos(2) <=size(pca_focusMap,2)
                pca_focusMap(offPos(1),offPos(2),(k-1)*m+1:k * m) = pca_coeff_gold(:,k);
%                 pca_focusMap(offPos(1),offPos(2),m+1:2*m) = pca_coeff_gold(:,2);
%                 pca_focusMap(offPos(1),offPos(2),2*m+1:3*m) = pca_coeff_gold(:,3);
            end
        end
        %     diffM = pca_coeff_ImpSample - pca_coeff_gold;
        %     dfn = norm(diffM, 'fro');
        %     disp(pca_coeff_gold - pca_coeff_ImpSample);
        %     disp(dfn);
    end
end
% calculate measurements for the vector maps.
figure,
subplot(2,3,1);imshow(pca_focusMap(:,:,1),[]);
subplot(2,3,2);imshow(pca_focusMap(:,:,2),[]);
subplot(2,3,3);imshow(pca_focusMap(:,:,3),[]);

subplot(2,3,4);imshow(pca_map_gold(:,:,1),[]); title("Gold Sample PCA");
subplot(2,3,5);imshow(pca_map_gold(:,:,2),[]);
subplot(2,3,6);imshow(pca_map_gold(:,:,3),[]);

function w= PowerHeuristic2Pdfs(nf, fPdf, ng, gPdf) 
 f = nf * fPdf;
 g = ng * gPdf;
 w = (f * f) / (f * f + g * g);
end

function wList = PowerHeuristicMPdfs(nList, pdfList)
m = size(nList,1);
fList = zeros(m,1);
wList = zeros(m,1);
denom = 0;
for i =1:m
    fList(i) = nList(i) * pdfList(i);
    denom = denom + fList(i)*fList(i);
end
%  w = (f * f) / (f * f + g * g);
if denom == 0
    return;
end
for i = 1:m
    wList(i) = (fList(i)*fList(i)) / denom;
end
end


function pca_coeff = findSamplesInBallwSampleMap(pos, imgSamples, imgDataList, R)
    % finding all samples within the Radius
    mDim = size(imgDataList,1);
    [maxY, maxX] = size(imgSamples);
    neighborVals = zeros((2*R+1)* (2*R+1), mDim);
    cnt = 1;
    pca_coeff = eye(3); %Eigenvectors are axis alignd and the ellipse is anisotropic
    for y = -int32(R):int32(R)
        for x = -int32(R):int32(R)
               % check if inside the disk and not on the central pixel
            if (x == 0 && y == 0) || (x*x + y*y > R * R)
                continue;
            end
            npos = int32(pos) + [x,y];
            
            if npos(1) < 1 || npos(2) < 1 ||... 
                npos(1) > maxY || npos(2) > maxX
                continue;
            end

            mask = imgSamples(npos(1),npos(2));
            if mask == 0
                continue;
            end
            % we care only positions with samples
            val = zeros(1,mDim);
            for m = 1:mDim
                val(m) = imgDataList{m}(npos(1),npos(2),:);
            end
            neighborVals(cnt,:)= val;
            cnt = cnt + 1;
        end
    end
    neighborVals(cnt:end,:) = []; % remove empty locations
    if(size(neighborVals,1)<=mDim)% if the number of neighbors is smaller than the dimensionality
%         pca_coeff = zeros(3); %Eigen vectors are zeros!
        return;
    end
    % compute statistics / PCA
    X = neighborVals;   
    [pca_coeff,score,latent] = pca(X);
%      m = size(latent,1);
%      for i = 1:m
%          pca_coeff(:,i) = pca_coeff(:,i) .* latent(i);
%      end
end

function pca_coeff = findSamplesInBallGold(pos, imgDataList, R)
    pca_coeff = zeros(3); %Eigen vectors are zeros!
    % finding all samples within the Radius
    mDim = size(imgDataList,1);
    if mDim == 0 
        return ;
    end
    [maxY, maxX] = size(imgDataList{1});
    neighborVals = zeros((2*R+1)* (2*R+1), mDim);
    cnt = 1;

    for y = -int32(R):int32(R)
        for x = -int32(R):int32(R)
               % check if inside the disk and not on the central pixel
            if (x == 0 && y == 0) || (x*x + y*y > R * R)
                continue;
            end
            npos = int32(pos) + [x,y];
            
            if npos(1) < 1 || npos(2) < 1 ||... 
                npos(1) > maxY || npos(2) > maxX
                continue;
            end

            % we care only positions with samples
            val = zeros(1,mDim);
            for m = 1:mDim
                val(m) = imgDataList{m}(npos(1),npos(2),:);
            end
            neighborVals(cnt,:)= val;
            cnt = cnt + 1;
        end
    end
    neighborVals(cnt:end,:) = []; % remove empty locations
    if(isempty(neighborVals))
%         pca_coeff = zeros(3); %Eigen vectors are zeros!
        return;
    end
    % compute statistics / PCA
    X = neighborVals;   
    [pca_coeff,score,latent] = pca(X);
%      m = size(latent,1);
%      for i = 1:m
%          pca_coeff(:,i) = pca_coeff(:,i) .* latent(i);
%      end
end

function M = sortEigVec(D)
M = D;   
for i = 1:size(D,2)
    V = D(:,i);
    [mv, ind] = max(V);   
    M(:,ind) = V;
end
end

function [corrsEvs,corrsbc] = sortAssocEvs(EvList, m, bc)
corrsbc = bc;
corrsEvs = EvList;
numVerts = size(EvList,1);
numVecPerVert = 3;
%The association matrix
matAssoc = zeros(numVerts, numVecPerVert);
for i = 1:numVecPerVert
    for j = 1:numVerts-1
        if matAssoc(j,i) == 0
            maxDot = -realmax;
            assocVecId = 1;
            assocVertId = 1;
            for k = j+1:numVerts
                for ii = 1:numVecPerVert
                    dotV = dot(EvList(j,(i-1)*m+1:i*m),EvList(k,(ii-1)*m+1:ii*m));
                    if dotV > maxDot
                        assocVecId = ii;
                        assocVertId = k;
                        maxDot = max(dotV, maxDot);
                    end
                end
            end
            matAssoc(j,i) = 1;
            matAssoc(assocVertId, assocVecId)  = 1;
            if assocVecId == i
                continue;
            end
            corrsEvs(assocVertId, (i-1)*m+1:i*m) = EvList(assocVertId, (assocVecId - 1)*m+1:assocVecId * m);
%             corrsbc(1,assocVertId)
        end
    end
end
end