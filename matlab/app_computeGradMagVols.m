% fileList = {'../data/Toothuc.nrrd'; '../data/Toothucderiv1.nrrd'; '../data/Toothucderiv2.nrrd'};
fileList = {'../data/BraTS/BraTS20_Training_001_t1.nhdr';...
    '../data/BraTS/BraTS20_Training_001_t1ce.nhdr';... 
    '../data/BraTS/BraTS20_Training_001_t2.nrrd';... 
    '../data/BraTS/BraTS20_Training_001_flair.nhdr'};
% myEigInfoReconstruct(fileList)
% fn = sprintf('eigVolDS%i_%0.2f.nrrd', 0, sampleRatePerAxis);
[Vlist, nVols, volInfoList] = loadVolumes(fileList);
gmVlist = cell(nVols,1);
for i=1:nVols
    [gmVlist{i},Gaz, Gelev] = imgradient3(Vlist{i},'central');
    % normalize the gradient magnitude
    gmVlist{i} = rescale(gmVlist{i});
end
% write out downSampled vols
for i = 1:nVols
    [filepath, name, ext] = fileparts(fileList{i});
    fileName = sprintf('%s/GradMag_%s.nrrd', filepath, name);
    ok = nrrdWriter(fileName, uint8(255.*gmVlist{i}), [1 1 1], [0 0 0], 'raw');
end
