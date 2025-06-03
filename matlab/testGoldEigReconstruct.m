function goldMajEigVols = testGoldEigReconstruct(fileList, zsampleRatePerAxis)
if nargin < 2
     goldSampleRatePerAxis = 1;
else
   goldSampleRatePerAxis = zsampleRatePerAxis;
end
%% create the gold version
disp('Computing the gold standard.');
tic

[goldMajEigVols, goldSecEigVols, majEigValVols, fullDim, goldDim] = myEigInfoCompute(fileList, goldSampleRatePerAxis);
toc
for i = 1:length(goldMajEigVols)
    [filepath, name, ext] = fileparts(fileList{i});
    fileName = sprintf('%s_eigMaj_gold%i_%0.2f.nrrd', name, i, goldSampleRatePerAxis);
    ok = nrrdWriter(fileName, single(goldMajEigVols{i}), [1 1 1], [0 0 0], 'raw');
    
    fileName2 = sprintf('%s_eigSec_gold%i_%0.2f.nrrd', name, i, goldSampleRatePerAxis);
    ok = nrrdWriter(fileName2, single(goldSecEigVols{i}), [1 1 1], [0 0 0], 'raw');
end
end