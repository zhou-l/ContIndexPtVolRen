majPCAData = load('majEigVDebug_computeTooth02.txt');
majPCA = majPCAData(:,2:4);

fileList = {'./Toothuc0.20.nrrd'; './Toothucderiv10.20.nrrd'; './Toothucderiv20.20.nrrd'};
[Vlist, nVols, volInfoList] = loadVolumes(fileList);
volSize = size(Vlist{1});
mu = zeros(volSize(1)*volSize(2)*volSize(3),nVols);
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
%%
mu = mu ./ 255;
idxPtData = calcPflats(mu, majPCA);
writematrix(idxPtData, 'matlabConvertedIdxPt.csv');

nIdxVols = size(idxPtData,2);
for i = 1:nIdxVols/2
    figure,plot(idxPtData(:,2*(i-1)+1),idxPtData(:,2*i),'.');
end
