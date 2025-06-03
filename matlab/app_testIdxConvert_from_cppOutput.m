majPCAData = readmatrix('majEigV_mu_tooth02.csv');
majPCA = majPCAData(:,5:7);
mu = majPCAData(:,2:4);
% mu = mu ./ 255;
idxPtData = calcPflats(mu, majPCA);
writematrix(idxPtData, 'matlabConvertedIdxPtX.csv');

nIdxVols = size(idxPtData,2);
for i = 1:nIdxVols/2
    figure,plot(idxPtData(:,2*(i-1)+1),idxPtData(:,2*i),'.');
end
