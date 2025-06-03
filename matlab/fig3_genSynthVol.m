close all;
featureSize = 40;
num_featuresamples = featureSize * featureSize * featureSize;
% num_samples = 1000;
X = linspace(0,pi,num_featuresamples);
X = X';
% Y = sin(X);
noise_std = 0.1;
noisePt_perX = 20;
noise = rand(num_featuresamples,1);
Y = sin(X);
Xnoise = X; %repmat(X,noisePt_perX,1);
Ynoise = Y; %repmat(Y,noisePt_perX,1);
Ynoise = noise .* noise_std + Ynoise;
Znoise = 4-Ynoise;
Xnoise = Xnoise ./ pi;
Ynoise = Ynoise ./ 1.2;
Znoise = Znoise ./ 4;
figure, scatter(X,Y);
% writeCSV(X,Y,'sin.csv');
figure, scatter(Ynoise, Znoise);
% feature = [Xnoise Ynoise Znoise];
feature = [X Y Znoise];

% another feature
lineFsize = 20;
num_lineFsamples = lineFsize * lineFsize * lineFsize;
Xl = linspace(0, 1, num_lineFsamples);
Xl = Xl';
Yl = 1 - 0.6 .* Xl;
Zl= 0.2 .* Xl;
Fl = [Xl Yl Zl];
% third feautre
lineF2Size = 12;
Fl2 = [Xl Zl Yl];
figure, scatter(Xl, Yl);
figure, scatter(Yl, Zl);
% writeCSV(Xnoise, Ynoise, 'sinNoise.csv');

volDim = 64;
close all;
numVols = 3;
Vlist = cell(numVols, 1);
numSamples = volDim * volDim *volDim;
%% create features with patterns in the value domain
mu = [0.1 0.2 0.1];
sig1 = 0.01;
sig2 = 0.02;
sig3 = 0.02;
rou1 = 0;%0.0;
rou2 = 0;%1;
rou3 = 0;%0.8;
cov = [ sig1*sig1, rou1*sig1*sig2, rou2*sig1*sig3;...
        rou1*sig2*sig1, sig2*sig2, rou3*sig2*sig3;...
        rou2*sig3*sig1, rou3*sig3*sig2, sig3*sig3;...
     ];
S = mvnrnd(mu,cov, numSamples);
figure, scatter(S(:,1), S(:,2)), title('S distr'); 

%%
rng(0,'twister');

ctr = [volDim/2 volDim/2 volDim/2];
for d  = 1:numVols
    V = zeros(volDim, volDim, volDim);
    Vlist{d,:} = V;
end

 % the background is gaussian 
for z = 1:volDim
    for y = 1:volDim
        for x = 1:volDim
            idx = randi(numSamples);
            idx2 = randi(num_featuresamples);
            idx3 = randi(num_lineFsamples);
            for d = 1:numVols
                Vlist{d}(x,y,z) = S(idx,d);
                % create features
                toCtr = [x y z] - ctr;
                dist = norm(toCtr);
                % the sine feature
                if abs(z - ctr(3)) <= featureSize/2 && abs(y - ctr(2)) <= featureSize/2 &&...
                        abs(x - ctr(1)) <= featureSize/2
                    Vlist{d}(x,y,z) = feature(idx2,d);
                end
                
                % % Line feature 1
                % if z >= volDim-lineFsize && y >= volDim-lineFsize && x >= volDim-lineFsize
                %     Vlist{d}(x,y,z) = Fl(idx3,d);
                % end
                
               % Line feature 2
                % if z < lineF2Size && y < lineF2Size && x < lineF2Size
                %     Vlist{d}(x,y,z) = Fl2(idx3,d);
                % end
      
            end

        end
    end
end

% collect the whole data
 cnt = 1;
 fullDist = zeros(numSamples, numVols);
for z = 1:volDim
    for y = 1:volDim
        for x = 1:volDim
            for d = 1:numVols
                fullDist(cnt,d) = Vlist{d}(x,y,z);
            end
            
            cnt = cnt + 1;
        end
    end
end
figure, hold on;
scatter(fullDist(:,1),fullDist(:,2), '.'), title('Whole data');
% scatter(F1(:,1),F1(:,2), 'filled', 'r');
% scatter(F2(:,1),F2(:,2), 'filled', 'g');
% scatter(F3(:,1),F3(:,2), 'filled', 'b');
hold off;
%%
name = 'fig3_synthVol3D';
for d = 1:int16(numVols)
    fileName = sprintf('%s%d.nrrd', name, d);
    % float
    ok = nrrdWriter(fileName, single(Vlist{d}), [1 1 1], [0 0 0], 'raw');
    % uint8
    % ok = nrrdWriter(fileName, uint8(255.*Vlist{i}), [1 1 1], [0 0 0], 'raw');
    % if i > 1
    %     figure; h = histogram2(Vlist{i-1},Vlist{i});
    % end
end

% Validataion
volInfoList = cell(numVols,1);
for d = 1:int16(numVols)
    % filename = fileList{i};
    fileName = sprintf('%s%d.nrrd', name, d);
    disp(fileName)
    headerInfo = nhdr_nrrd_read(fileName, true);
    volInfoList{d} = headerInfo;
    V = headerInfo.data;
    dimX = size(V,1);
    dimY = size(V,2);
    dimZ = size(V,3);
     cnt = 1;
    for z = 1:dimZ
        for y = 1:dimY
            for x = 1:dimX
               if isnan(V(x,y,z))
                   V(x,y,z) = -6;
               end
               fullDist(cnt,d) = V(x,y,z);
               cnt = cnt+1;
            end
        end
    end
%     volumeViewer(V);
    Vlist{d} = V;
end 
figure, hold on;
scatter(fullDist(:,1),fullDist(:,2), '.'), title('Loaded whole data');
hold off;