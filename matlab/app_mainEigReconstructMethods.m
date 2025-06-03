% PCA computation for local linear correlations in volumes with octree
% accesleration.

% % toy examples
fileList = {'./Toothuc0.20.nrrd'; './Toothucderiv10.20.nrrd'; './Toothucderiv20.20.nrrd'};
% fileList = {'./Pf25cpUs0.40.nrrd'; './TCf25cpUs0.40.nrrd'; './QVAPORf25cpUs0.40.nrrd'; './PRECIPf25cpUs0.40.nrrd'};

% full-sized
% "BraTS20_Training_001_t1.nhdr" "BraTS20_Training_001_t1ce.nhdr" "BraTS20_Training_001_t2.nrrd" "BraTS20_Training_001_flair.nhdr"
% fileList = {'../data/BraTS/BraTS20_Training_001_t1.nhdr'; '../data/BraTS/BraTS20_Training_001_t1ce.nhdr';...
%     '../data/BraTS/BraTS20_Training_001_t2.nrrd'; '../data/BraTS/BraTS20_Training_001_flair.nhdr'};
% fileList = {'../data/Toothuc.nrrd'; '../data/Toothucderiv1.nrrd'; '../data/Toothucderiv2.nrrd'};
% fileList = {'../data/Pf25cpUs.nhdr'; '../data/TCf25cpUs.nhdr'; '../data/QVAPORf25cpUs.nhdr'; '../data/PRECIPf25cpUs.nhdr'};
%%
% octree eigenvectors computation
valDiff = 0.02;
% pcaOctVols = testOctreeEigReconstruct(fileList, valDiff);
pcaOctVols = testOctreeEigReconstInterp(fileList, valDiff);
%% gold-standard
pcaGoldVols = testGoldEigReconstruct(fileList);