% The entry point of the conversion application

% M = load('ImgPlaneIdxPtDebug_load2.txt');
%X = M(:,1:2);
%  M = load('idxPtDebug_computeTooth06.txt');
% M = load('idxPtDebug_computeBraTs.txt');
% inputName = 'idxPtDebug_computeTooth06.txt';
% inputName = 'idxPtDebug_computeIsabel045.txt';
% inputName = 'idxPtDebug_computeSynthNew.txt';
% inputName = 'idxPtDebug_computeIsabel045.txt';
% inputName = 'idxPtDebug_computeSynthNew.txt';
%  inputName = 'idxPtDebug_computeTooth.txt';
% inputName = './DTI/idxPtDebug_computeGK2.txt';
% inputName = './DTI/idxPtDebug_compute_gk2clcpcs.txt';
%  inputName = 'idxPtDebug_computeFaParaPeri.txt';
% inputName = 'idxPtDebug_compute_p2FAParaPerp.txt';
% inputName = 'idxPtDebug_compute_p1RAFATr.txt';
%  inputName = 'idxPtDebug_compute_p2_RAFATr.txt';
% inputName = 'idxPtDebug_compute_p2LinPlaSph.txt';
% inputName = 'idxPtDebug_compute_p1LinPlaSph.txt';
% inputName = 'idxPtDebug_compute_p2TrFaLin.txt';
% inputName = 'idxPtDebug_compute_p2PerpFATr.txt';
% inputName = 'idxPtDebug_computeDTI2_FaTrPerp.txt';
inputName = 'idxPtDebug_computeBratsT1_gradmag.txt';
%% 2-flats
% inputName = 'twoFlats_idxPtDebug_computeHalfCyl3d.txt';
% inputName = 'twoFlats_idxPtDebug_computeCtbl3d09.txt';
% inputName = 'twoFlats_idxPtDebug_computeTornado10.txt';
% inputName = './DTI/twoFlats_idxPtDebug_computeGK2.txt';
% inputName = './DTI/twoFlats_idxPtDebug_compute_gk2clcpcs.txt';
% inputName = 'idxPtDebug_computeHalfCyl3d.txt';
% inputName = '../data/BraTS/indexedPtVols/idxPtDebug_computeBraTsT1.txt';
% inputName = 'idxPtDebug_computeIsabel045.txt';
% inputName = 'idxPtDebug_computeSynthNew.txt';
% inputName = 'idxPtDebug_computeBraTs045.txt';
% inputName = 'idxPtDebug_computeBraTs.txt';
% inputName = 'twoFlats_idxPtDebug_computeSynthNew.txt';
% inputName = 'twoFlats_idxPtDebug_computeTooth.txt';
% inputName = 'twoFlats_idxPtDebug_computeDTI12D.txt';
% inputName = 'twoFlats_idxPtDebug_computeFaParaPeri.txt';
% inputName = 'twoFlats_idxPtDebug_compute_p1FAParaPerp.txt';
%  inputName = 'twoFlats_idxPtDebug_compute_p1RAFATr.txt';
%  inputName = 'twoFlats_idxPtDebug_compute_p2_RAFATr.txt';
% inputName = 'twoFlats_idxPtDebug_compute_p2LinPlaSph.txt';
% inputName = 'twoFlats_idxPtDebug_compute_p1LinPlaSph.txt';
% inputName = 'twoFlats_idxPtDebug_compute_p2TrFaLin.txt';
% inputName = 'twoFlats_idxPtDebug_compute_p2PerpFATr.txt';
% inputName = 'twoFlats_idxPtDebug_compute_dtiEyeP2.txt';
% inputName = 'twoFlats_idxPtDebug_compute_synthetic3D.txt';
% inputName = 'twoFlats_idxPtDebug_compute_FaTrPerp.txt';
% inputName = 'twoFlats_idxPtDebug_computeIsabel3Vars.txt';
%% The working function
buildContIdxPt_dynamicKDE(inputName);

