function [eigVecs, eigVals] = localPCAspaceDom(pos, vols, R, maxCnt)
if nargin < 4
    zmaxCnt = 100;
else
    zmaxCnt = maxCnt;
end

nVols = length(vols);
if nVols < 1
    eigVecs = [];
    eigVals = [];
    return;
end

D = zeros(zmaxCnt,nVols);
% cnt = 1;
dim = size(vols{1});
% compute with Monte Carlo sampling!!!
s = rng;
for i =1:zmaxCnt
    % get random position
    offset = 2.*R.*rand(1,3) - [R,R,R];
    % with integral position
    %             npos = pos + [x y z];
    % with random position
    npos = pos + offset;
    if npos(1) < 1 || npos(2) < 1 || npos(3) < 1 ||...
            npos(1) > dim(1) || npos(2) > dim(2) || npos(3) > dim(3)
        continue;
    end
    val = zeros(1,nVols);
    for d = 1:nVols
        % with value fetching
        %                val(1,d) = vols{d}(npos(1),npos(2),npos(3));
        % with interpolation
        val(1,d) = volInterp(vols{d}, npos);
    end
    D(i,:)  = val;
end
% the columns of the input data are variables, rows are data samples
DM = D;% D';
% each column is a eigenvector
% we need only top 'nVols' components
eigVecs = pca(DM, 'NumComponents', nVols);
eigVals = vecnorm(eigVecs);
% if abs(eigVals(1)- eigVals(2))<1e-8 && abs(eigVals(2)- eigVals(3))<1e-8
%   eigVecs = sortEigVec(eigVecs);
% end

end