% toy examples
% fileList = {'./Toothuc0.20.nrrd'; './Toothucderiv10.20.nrrd'; './Toothucderiv20.20.nrrd'};
fileList = {'./Pf25cpUs0.20.nrrd'; './TCf25cpUs0.20.nrrd'; './QVAPORf25cpUs0.20.nrrd'; './PRECIPf25cpUs0.20.nrrd'};

% full-sized
% fileList = {'../data/Toothuc.nrrd'; '../data/Toothucderiv1.nrrd'; '../data/Toothucderiv2.nrrd'};
%fileList = {'../data/Pf25cpUs.nhdr'; '../data/TCf25cpUs.nhdr'; '../data/QVAPORf25cpUs.nhdr'; '../data/PRECIPf25cpUs.nhdr'};

% octree eigenvectors computation
valDiff = 0.02;
% pcaOctVols = testOctreeEigReconstruct(fileList, valDiff);
pcaOctVols = testOctreeEigReconstInterp(fileList, valDiff);
%% gold-standard
pcaGoldVols = testGoldEigReconstruct(fileList);