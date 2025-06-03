function [majEigVols, secEigVols, majEigValVols, fullDim, dsDim] = myEigInfoCompute(fileList, sampleRatePerAxis)
% load the nhdr/nrrd file
nVols = length(fileList);
Vlist = cell(nVols, 1);
for i = 1:nVols
    filename = fileList{i};
    headerInfo = nhdr_nrrd_read(filename, true);
    [filepath, name, ext] = fileparts(filename);
    V = headerInfo.data;
    dimX = size(V,1);
    dimY = size(V,2);
    dimZ = size(V,3);
    for z = 1:dimZ
        for y = 1:dimY
            for x = 1:dimX
               if isnan(V(x,y,z))
                   V(x,y,z) = -6;
               end
            end
        end
    end
%     volumeViewer(V);
    Vlist{i} = V;
end 
fullDim = [dimX dimY dimZ];
% downsample the data volumes
sratePerAxis = sampleRatePerAxis;
VsList = myVolDownSample(Vlist,sratePerAxis);
dsDimX = size(VsList{1},1);
dsDimY = size(VsList{1},2);
dsDimZ = size(VsList{1},3);
dsDim = [dsDimX dsDimY dsDimZ];
% compute the eigen information for every voxel in the downsampled volumes
R = 3; % radius of the neighborhood
majEigVols = cell(nVols,1);
secEigVols = cell(nVols,1);
for i = 1:nVols
    majEigVols{i} = zeros(dsDimX,dsDimY,dsDimZ);
    secEigVols{i} = zeros(dsDimX,dsDimY,dsDimZ);
end    
majEigValVols = zeros(dsDimX,dsDimY,dsDimZ);

for z = 1:dsDimZ
    for y = 1:dsDimY
        for x = 1:dsDimX
            pos = [x y z];
            [eigVecList, eigValList] = localPCAspaceDom(pos, VsList, R); % spatial domain local PCA
            for d = 1:nVols
                majEigVols{d}(pos(1),pos(2),pos(3)) = eigVecList(d,1); % the major eigenvector
                secEigVols{d}(pos(1),pos(2),pos(3)) = eigVecList(d,2); % the second eigenvector
                majEigValVols(pos(1),pos(2),pos(3)) = eigValList(1); % the eigenvalue of the major eigenvector
            end
        end
    end
end
end
