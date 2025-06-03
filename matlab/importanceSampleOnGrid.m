% useCase = 'isabel';
useCase = 'tooth';
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

%  close all;
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
% the radius of PCA analysis
R = 3;
multiImg = max(max(gI0,gI1), gI2);
% figure, imshow(multiImg);
Kernel = zeros(2*R+1);
Kctr = [R+1 R+1];
cnt = 0;
for r = 1:2*R+1
    for c = 1:2*R+1
        v = [r c] - Kctr;
        if norm(v) <= R
            Kernel(r,c) = 1;
            cnt = cnt + 1;
        end
    end
end
Kernel = Kernel ./ cnt;
% try with convolution
convMImg = filter2(Kernel,multiImg);
% figure, imshow(convMImg, []);
multiImg = convMImg;
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
scalarImgs = cell(3,1);
scalarImgs{1,:} = I0;
scalarImgs{2,:} = I1;
scalarImgs{3,:} = I2;
% test interpolation on a single triangle
% testInterpTriangle(scalarImgs, m);


for i=1:N
    targetSum(i,:) = sumImg(idxR(i), idxC(i));
%     % with maximum heuristic
    targetMax(i,:) = multiImg(idxR(i), idxC(i));
    % try MIP with power heuristic
    pdfList = [gI0(idxR(i), idxC(i)); gI1(idxR(i), idxC(i)); gI2(idxR(i), idxC(i))];
    wList = PowerHeuristicMPdfs(nList, pdfList);
    targetPower(i,:) = 0;
    for d = 1:m
        targetPower(i,:) = targetPower(i,:) + wList(d) * pdfList(d) ;
    end
end


%target = targetPower;
target = targetMax;
% calculate importance weight
w = target ./ proposal;
% normalization
w = w ./ sum(w);

wPower = targetPower ./ proposal;
% wPower = wPower ./ sum(wPower);

wSum = targetSum ./ proposal;
wSum = wSum ./ sum(wSum);

wMax = targetMax./ proposal;
wMax = wMax ./ sum(wMax);

% resample, with replacement, according to importance weight
sample_ind = randsample([1:N],N,true,w);

sample_ind_Power = randsample([1:N],N,true,wPower);
sample_ind_Sum = randsample([1:N],N,true,wSum);
sample_ind_Max = randsample([1:N],N,true,wMax);


y_max_samples = mean_samples(sample_ind_Max);
y_sum_samples = mean_samples(sample_ind_Sum);
y_power_samples = mean_samples(sample_ind_Power);
x_max_samples = sigma_samples(sample_ind_Max);
x_sum_samples = sigma_samples(sample_ind_Sum);
x_power_samples = sigma_samples(sample_ind_Power);

%plot(true_mean, true_sigma, 'r.','MarkerSize',5^2)
figure;
 plot((y_max_samples), x_max_samples,'.'); 
 figure;
subplot(1,3,1); plot((y_max_samples), x_max_samples,'.'); title("Max MIP");
subplot(1,3,2); plot((y_sum_samples), x_sum_samples, '.'); title("Sum MIP");
subplot(1,3,3); plot((y_power_samples), x_power_samples, '.'); 
title('Power MIP');
%%
% set an image for samples
imgSamples = zeros(maxR, maxC);
for i = 1:size(x_max_samples,1)
    imgSamples(y_max_samples(i), x_max_samples(i)) = 1;
end
% test neighborhood computation

 % store the top 3 eigenvectors of m dimensions (3*m layers)
pca_map_impSamp = zeros(maxR, maxC,3 * m);
pca_map_gold = zeros(maxR, maxC, 3*m );

% concat two arrays
sampleList = cat(2, y_max_samples, x_max_samples);
% remove redundant positions
sampleList = unique(sampleList, 'rows');
disp('sampling Rate;');
disp(length(sampleList)/(maxR*maxC));
% store the associated pca result
nEV = 3; % number of eigenvectors
pcaImpSampList = zeros(size(sampleList,1),nEV * m);
pcaGoldList = pcaImpSampList;
for i = 1:size(sampleList,1)
    pos = sampleList(i,:);
    pca_coeff_ImpSample = findSamplesInBallwSampleMap(pos, imgSamples, scalarImgs, R);
    pca_coeff_gold = findSamplesInBallGold(pos, scalarImgs, R);

%     disp('ImpSample Coeff');
%     disp(pca_coeff_ImpSample);
    % sort the matrix by its maximum component
%     pca_coeff_ImpSample = sortEigVec(pca_coeff_ImpSample);
%     pca_coeff_gold = sortEigVec(pca_coeff_gold);
    
    for k = 1:m
         % store the eigenvector information
        pcaImpSampList(i,1:m) = pca_coeff_ImpSample(:,1);
        pcaImpSampList(i,m+1:2*m) = pca_coeff_ImpSample(:,2);
        pcaImpSampList(i,2*m+1:3*m) = pca_coeff_ImpSample(:,3);
        pcaGoldList(i,1:m) = pca_coeff_gold(:,1);
        pcaGoldList(i,m+1:2*m) = pca_coeff_gold(:,2);
        pcaGoldList(i,2*m+1:3*m) = pca_coeff_gold(:,3);


        pca_map_impSamp(pos(1),pos(2),1:m) = pca_coeff_ImpSample(:,1);
        pca_map_impSamp(pos(1),pos(2),m+1:2*m) = pca_coeff_ImpSample(:,2);
        pca_map_impSamp(pos(1),pos(2),2*m+1:3*m) = pca_coeff_ImpSample(:,3);
        
        pca_map_gold(pos(1),pos(2),1:m) = pca_coeff_gold(:,1);
        pca_map_gold(pos(1),pos(2),m+1:2*m) = pca_coeff_gold(:,2);
        pca_map_gold(pos(1),pos(2),2*m+1:3*m) = pca_coeff_gold(:,3);
    end
%     diffM = pca_coeff_ImpSample - pca_coeff_gold;
%     dfn = norm(diffM, 'fro');
%     disp(pca_coeff_gold - pca_coeff_ImpSample);
%     disp(dfn);
end

%%
% test interpolation
X = sampleList;%cat(2, y_max_samples, x_max_samples);
DT = delaunayTriangulation(X);
figure,triplot(DT);
% % find the triangle for testing
% Pq = [252, 272];
% [ti,bc] = pointLocation(DT,Pq);
% ptIds = DT(ti,:);
% % figure, triplot(DT(ti,:));
% posXY = DT.Points(ptIds,:);
% figure, triplot(DT(ti,:), posXY(:,1), posXY(:,2));
% triVals = pcaImpSampList(ptIds,:);
% triValsGold = pcaGoldList(ptIds,:);

figure,
subplot(1,3,1); title('ImportSample PCA samples'); imshow(pca_map_impSamp(:,:,1),[]);
subplot(1,3,2); imshow(pca_map_impSamp(:,:,2),[]);
subplot(1,3,3); imshow(pca_map_impSamp(:,:,3),[]);
% storing interpolated Eigenvectors
EvInterp = zeros(1, nEV * m);

xrange = [240, 269];
yrange = [220, 249];
pca_focusMap = zeros(yrange(2)-yrange(1),xrange(2)-xrange(1),3 * m);
% for y = yrange(1,1):yrange(1,2)
for y = 1:maxR
    %     for x = xrange(1,1):xrange(1,2)
    for x = 1:maxC
        Pq = [y,x];
        offPos = [y,x] - [yrange(1,1),xrange(1,1)] + [1 1];
        if imgSamples(y,x) == 0
            evList = [];
            % bc is the barycentric coordinate of Pq and ti is the triangle
            % index in DT
            [ti,bc] = pointLocation(DT,Pq);
            if isnan(ti) 
             continue;
            end
            % fetch the data values associated with the triangle
            ptIds = DT(ti,:);
            triVals = pcaImpSampList(ptIds,:);
            triValsGold = pcaGoldList(ptIds,:);
            % find the correct correspondence of vectors
%             corresEvs = triVals;
            corresEvs = sortAssocEvs(triVals, m, bc);
%             corresEvGold = triValsGold;
            corresEvGold = sortAssocEvs(triValsGold, m, bc);
            % and then perform the interpolation
            for nn = 1:nEV
                EvInterp(1,(nn-1)*m + 1:nn*m) = eigInterp2D(corresEvs(:,(nn-1)*m+1:nn*m), bc);
                pca_map_impSamp(y,x,(nn-1)*m + 1:nn*m) = EvInterp(1,(nn-1)*m + 1:nn*m) ;
                % interpolate the gold computation
                EV = eigInterp2D(corresEvGold(:,(nn-1)*m+1:nn*m), bc);
                pca_map_gold(y,x,(nn-1)*m + 1:nn*m) = EV;
                if nn == 1 && offPos(1) > 0 && offPos(2) > 0 &&... 
                offPos(1) <size(pca_focusMap,1) && offPos(2) <size(pca_focusMap,2)
                    pca_focusMap(offPos(1),offPos(2),1) = EV(1,1);
                    pca_focusMap(offPos(1),offPos(2),2) = EV(1,2);
                    pca_focusMap(offPos(1),offPos(2),3) = EV(1,3);
                end
            end
      
        
%             ev1 = eigInterp2D(corresEvs(:,1:m), bc);
%             ev2 = eigInterp2D(corresEvs(:,m+1:2*m), bc);
%             ev3 = eigInterp2D(corresEvs(:,2*m+1:3*m), bc);           
%             Idx = knnsearch(sampleList,Pq,'K',4);
%             knownPos = sampleList(Idx,:);
%             for nn = 1:size(knownPos,1)
%                 evList(end+1,:) = pca_map_impSamp(knownPos(nn,1),knownPos(nn,2),:);
%             end
            % do the interpolation for [x,y]
%             [Ei,evi]=eigInterp(evList, knownPos,);

        end
    end
end
figure,
subplot(2,3,1); imshow(pca_focusMap(:,:,1),[]);
subplot(2,3,2); imshow(pca_focusMap(:,:,2),[]);
subplot(2,3,3); imshow(pca_focusMap(:,:,3),[]);

% calculate measurements for the vector maps.
% figure,
% subplot(2,3,1); title('ImportSample PCA'); imshow(pca_map_impSamp(:,:,1),[]);
% subplot(2,3,2); imshow(pca_map_impSamp(:,:,2),[]);
% subplot(2,3,3); imshow(pca_map_impSamp(:,:,3),[]);
% 
% figure;
subplot(2,3,4); imshow(pca_map_impSamp(:,:,1),[]);title('Importance Sample PCA'); 
% hold on; plot( x_max_samples, y_max_samples, '.'); 
subplot(2,3,5); imshow(pca_map_impSamp(:,:,2),[]);
% hold on; plot(x_max_samples, y_max_samples, '.'); 
subplot(2,3,6); imshow(pca_map_impSamp(:,:,3),[]);
% hold on; plot(x_max_samples, y_max_samples, '.'); 
hold off;

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

function testInterpTriangle(imgDataList, m)
%     pts = [253, 261; 257 261; 255 258];
    pts = [452, 412; 448 420; 455 418; 455 412; 448, 412; 455, 420 ];
    nEV = 3; % the number of eigenvectors 
    Vlist = zeros(size(pts,1),3*m);
    xmin = min(pts(:,1));
    xmax = max(pts(:,1));
    ymin = min(pts(:,2));
    ymax = max(pts(:,2));
    R = 3;
    pca_test_gold = zeros(ymax-ymin+1, xmax-xmin+1, 3*m );
    pca_test_interp = zeros(ymax-ymin+1, xmax-xmin+1, 3*m );

    % find the gold solution
    for x = xmin:xmax
        for y = ymin:ymax
              pca_gold = findSamplesInBallGold([y,x], imgDataList, R);
              offpos = [y,x] - [ymin, xmin] + 1;
              for k = 1:nEV
                pca_test_gold(offpos(1),offpos(2),(k-1)*m+1:k * m) = pca_gold(:,k);
%                 pca_test_gold(offpos(1),offpos(2),m+1:2*m) = pca_gold(:,2);
%                 pca_test_gold(offpos(1),offpos(2),2*m+1:3*m) = pca_gold(:,3);
              end
        end
    end

    % interpolation scheme
    pp = pts - [xmin ymin] + 1;
    for i = 1:size(pp,1)
        Vlist(i,:) = pca_test_gold(pp(i,2),pp(i,1),:);
    end
    meanV = mean(Vlist, 1);
    C = [1 2 3; 1 3 4; 1 2 5; 2 3 6];
    TR = triangulation(C,pts);
    triplot(TR);
    for x = xmin:xmax
        for y = ymin:ymax
            Pq = [x y];
            [ti,bc] = pointLocation(TR,Pq);
            offpos = [y,x] - [ymin,xmin] + 1;
            if isnan(ti)
                % outside of the triangle
                pca_test_interp(offpos(1),offpos(2),:) = meanV;
                continue;
            end
            % find the correct correspondence of vectors
            ptIdx = TR.ConnectivityList(ti,:);
            corresEvs = sortAssocEvs(Vlist(ptIdx,:), m, bc);

            % and then perform the interpolation
            for nn = 1:nEV
                EvInterp(1,(nn-1)*m + 1:nn*m) = eigInterp2D(corresEvs(:,(nn-1)*m+1:nn*m), bc);
                pca_test_interp(offpos(1),offpos(2),(nn-1)*m + 1:nn*m) = EvInterp(1,(nn-1)*m + 1:nn*m) ;
                % interpolate the gold computation          
            end
        end
    end
   figure;
   subplot(2,3,1);imshow(pca_test_interp(:,:,1),[]);
   subplot(2,3,2);imshow(pca_test_interp(:,:,2),[]);
   subplot(2,3,3);imshow(pca_test_interp(:,:,3),[]);

   subplot(2,3,4);imshow(pca_test_gold(:,:,1),[]);
   subplot(2,3,5);imshow(pca_test_gold(:,:,2),[]);
   subplot(2,3,6);imshow(pca_test_gold(:,:,3),[]);
end

function pca_coeff = findSamplesInBallwSampleMap(pos, imgSamples, imgDataList, R) % pos is [y, x]
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
            npos = int32(pos) + [y,x];
            
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

function pca_coeff = findSamplesInBallGold(pos, imgDataList, R) % pos is [y, x]
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
            npos = int32(pos) + [y,x];
            
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