function [fullSzEig, majEigVols, Vlist] = myEigInfoReconstruct(fileList)
nVols = length(fileList);
sampleRatePerAxis = 0.1;

[Vlist, nVols, volInfoList] = loadVolumes(fileList);
% compute the eigen information of 
[majEigVols, majEigValVols, fullDim, dsDim] = myEigInfoCompute(fileList, sampleRatePerAxis);
for i = 1:length(majEigVols)
    fileName = sprintf('eigVolDS%i.nrrd', i);
    ok = nrrdWriter(fileName, single(majEigVols{i}), [1 1 1], [0 0 0], 'raw');
end
% reconstruction with interpolate
reconstDim = floor(fullDim .* sampleRatePerAxis * 3);
fullSzEig = myInterpEigVols(nVols, reconstDim, dsDim, majEigVols);

for i = 1:length(fullSzEig)
    fileName = sprintf('eigVolReconst%i.nrrd', i);
    ok = nrrdWriter(fileName, single(fullSzEig{i}), [1 1 1], [0 0 0], 'raw');
end

end
