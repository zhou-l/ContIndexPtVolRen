function volEigInfoGold(fileList, goldSampleRatePerAxis)
%% create the gold version
disp('Computing the gold standard.');
tic
% goldSampleRatePerAxis = newSampleRatePerAxis;
[goldMajEigVols, majEigValVols, fullDim, goldDim] = myEigInfoCompute(fileList, goldSampleRatePerAxis);
toc
for i = 1:length(goldMajEigVols)
    [filepath, name, ext] = fileparts(fileList{i});
    fileName = sprintf('%s_eigVolGold%i_%0.2f.nrrd', name, i, goldSampleRatePerAxis);
    ok = nrrdWriter(fileName, single(goldMajEigVols{i}), [1 1 1], [0 0 0], 'raw');
end
end