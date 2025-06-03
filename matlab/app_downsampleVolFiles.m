% fileList = {'../data/Toothuc.nrrd'; '../data/Toothucderiv1.nrrd'; '../data/Toothucderiv2.nrrd'};
fileList = {'../data/Pf25cpUs.nhdr'; '../data/TCf25cpUs.nhdr'; '../data/QVAPORf25cpUs.nhdr'; '../data/PRECIPf25cpUs.nhdr'};
% myEigInfoReconstruct(fileList)
% fn = sprintf('eigVolDS%i_%0.2f.nrrd', 0, sampleRatePerAxis);
[Vlist, nVols, volInfoList] = loadVolumes(fileList);
dsVlist = cell(nVols,1);
scale = 0.4;
for i=1:nVols
    dsVlist{i} = imresize3(Vlist{i},scale);
end
% write out downSampled vols
for i = 1:nVols
    [filepath, name, ext] = fileparts(fileList{i});
    fileName = sprintf('%s%0.2f.nrrd', name, scale);
    ok = nrrdWriter(fileName, single(dsVlist{i}), [1 1 1], [0 0 0], 'raw');
end
