%% load data volumes
fileList = {'../data/Toothuc.nrrd'; '../data/Toothucderiv1.nrrd'; '../data/Toothucderiv2.nrrd'};
% fileList = {'../data/Pf25cpUs.nhdr'; '../data/TCf25cpUs.nhdr'; '../data/QVAPORf25cpUs.nhdr'; '../data/PRECIPf25cpUs.nhdr'};
% myEigInfoReconstruct(fileList)
nVols = length(fileList);

% fn = sprintf('eigVolDS%i_%0.2f.nrrd', 0, sampleRatePerAxis);
[Vlist, nVols, volInfoList] = loadVolumes(fileList);
% build an gradient magnitude volume combining all attributes
volSize = size(Vlist{1});
VgradSum = zeros(volSize(1),volSize(2),volSize(3));
VgradList = Vlist;
for i = 1:nVols
    [Gmag,Gaz,Gel] = imgradient3(Vlist{i});
    VgradList{i} = Gmag;
    gmin = min(min(min(Gmag)));
    gmax = max(max(max(Gmag)));
    Gmag = (Gmag - gmin) ./ (gmax-gmin);
    VgradSum = VgradSum + Gmag;
end
figure, 
montage(reshape(VgradSum,volSize(1),volSize(2),1,volSize(3)),'DisplayRange',[])
title('Gradient magnitude')
%% compute downsampled eigen info volumes
tic
sampleRatePerAxis = 0.2;
% compute the eigen information of  
disp('Computing the downsampled eigen information volumes.');
[majEigVols, majEigValVols, fullDim, dsDim] = myEigInfoCompute(fileList, sampleRatePerAxis);
for i = 1:length(majEigVols)
    [filepath, name, ext] = fileparts(fileList{i});
    fileName = sprintf('%s_eigVolDS%i_%0.2f.nrrd', name, i, sampleRatePerAxis);
    ok = nrrdWriter(fileName, single(majEigVols{i}), [1 1 1], [0 0 0], 'raw');
end
toc


%% reconstruction with interpolate
tic
disp('Working on the interpolation.');
newSampleRatePerAxis = sampleRatePerAxis * 2;
reconstDim = floor(fullDim .* newSampleRatePerAxis);
fullSzEig = myInterpEigVols(nVols, reconstDim, dsDim, majEigVols);
toc
for i = 1:length(fullSzEig)
    [filepath, name, ext] = fileparts(fileList{i});
    fileName = sprintf('%s_eigVolReconst%i_%0.2f.nrrd', name, i, newSampleRatePerAxis);
    ok = nrrdWriter(fileName, single(fullSzEig{i}), [1 1 1], [0 0 0], 'raw');
end


%% create the gold version
disp('Computing the gold standard.');
tic
goldSampleRatePerAxis = newSampleRatePerAxis;
[goldMajEigVols, majEigValVols, fullDim, goldDim] = myEigInfoCompute(fileList, goldSampleRatePerAxis);
toc
for i = 1:length(goldMajEigVols)
    [filepath, name, ext] = fileparts(fileList{i});
    fileName = sprintf('%s_eigVolGold%i_%0.2f.nrrd', name, i, goldSampleRatePerAxis);
    ok = nrrdWriter(fileName, single(goldMajEigVols{i}), [1 1 1], [0 0 0], 'raw');
end
