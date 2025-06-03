% fileList = {'../data/BraTS/BraTS20_Training_001_t1.nhdr'; '../data/BraTS/BraTS20_Training_001_t1ce.nhdr';...
%     '../data/BraTS/BraTS20_Training_001_t2.nrrd'; '../data/BraTS/BraTS20_Training_001_flair.nhdr'};
fileList = {'./Toothuc0.20.nrrd'; './Toothucderiv10.20.nrrd'; './Toothucderiv20.20.nrrd'};
majEigFileList = fileList;
nameOnlyList = fileList;
for i = 1:length(majEigFileList)
    [filepath, name, ext] = fileparts(fileList{i});
%     fileName = sprintf('%s_Oct_eigVolReconst.nrrd', name);
    fileName = sprintf('%s_eigVol_gold%i_1.00.nrrd', name, i);
    nameOnlyList{i} = name;
    majEigFileList{i} = fileName;
end

commonName = commonString(nameOnlyList);

[Vlist, nVols, volInfoList] = loadVolumes(fileList);
[majEVlist, mVols, majEvolInfoList] = loadVolumes(majEigFileList);
volSize = size(Vlist{1});
numVoxels = volSize(1)*volSize(2)*volSize(3);
rangeMu = zeros(nVols, 2);
for i = 1:nVols
    minMu = min(min(min(Vlist{i})));
    maxMu = max(max(max(Vlist{i})));
    rangeMu(i,:) = [minMu, maxMu];
    disp(rangeMu(i,:));
end
mu = zeros(numVoxels,nVols);
majEig = zeros(numVoxels,mVols);

% maxScalVal = 255;
for i = 1:nVols
    range = rangeMu(i,2) - rangeMu(i,1);
%     if range == 0
%         range = 255;
%     end
%     
%     mu(:,i) = (reshape(Vlist{i},numVoxels,1) - rangeMu(i,1)) ./ range;
    majEig(:,i) = reshape(majEVlist{i},numVoxels,1);
end

% compute the mean value using the neighborhood information!
% mu has to be computed from the input volume data!!!
R = 3;
for z = 1:volSize(3)
    for y = 1:volSize(2)
        for x = 1:volSize(1)
            pos = [x y z];
            id = ((z-1)*volSize(2) + (y-1))*volSize(1)+x;
            localMu = localMu_spaceDom(pos, Vlist, R); % spatial domain local PCA
              mu(id,:) = localMu;
        end
    end
end
 
% %% load mu from C++ output file
% majPCAData = readmatrix('majEigV_mu_tooth02.csv');
% % majPCA = majPCAData(:,5:7);
% mu = majPCAData(:,2:4);
%%
% mu = mu ./ 255;
% convert to indexed points
idxPtData = calcPflats(mu, majEig);
nIdxVols = size(idxPtData,2);
for i = 1:nIdxVols/2
    figure,plot(idxPtData(:,2*(i-1)+1),idxPtData(:,2*i),'.');
end

% write out idx vols
idxVols = cell(nIdxVols,1);
idxVolFileList = cell(nIdxVols,1);
for i = 1:nIdxVols
    idxVol = reshape(idxPtData(:,i), volSize(1),volSize(2),volSize(3));
    if mod(i-1,2) == 0
        coordName = 'x';
    else
        coordName = 'y';
    end
    idxVols{i} = idxVol; 
    name = commonName;
    id = int16(floor((i-1)/2));
    fileName = sprintf('%s_idxVol%i_%s.nrrd', name, id, coordName);
    idxVolFileList{i} = fileName; 
end
[outFnames, okList]= writeVolumes(idxVols, idxVolFileList);