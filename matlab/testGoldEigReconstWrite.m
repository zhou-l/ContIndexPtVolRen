
% % toy examples
fileList = {'./Toothuc0.20.nrrd'; './Toothucderiv10.20.nrrd'; './Toothucderiv20.20.nrrd'};
% fileList = {'./Pf25cpUs0.40.nrrd'; './TCf25cpUs0.40.nrrd'; './QVAPORf25cpUs0.40.nrrd'; './PRECIPf25cpUs0.40.nrrd'};

% octree eigenvectors computation
valDiff = 0.02;
%% gold-standard
[majEigVols, majEigValVols, dim1, dim2] = testGoldEigReconstruct(fileList);
for i = 1:length(majEigVols)
end


% also write to a text file
write();