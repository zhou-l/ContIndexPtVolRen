DM = load('neighborTest0.txt');

eigVecs = pca(DM, 'NumComponents', 3);
mu = mean(DM)
eigVals = vecnorm(eigVecs);
if abs(eigVals(1)- eigVals(2))<1e-8 && abs(eigVals(2)- eigVals(3))<1e-8
  eigVecs = sortEigVec(eigVecs);
end
majEigVec = transpose(eigVecs(:,1))
secEigVec = eigVecs(:,2)
idxPtData = calcPflats(mu, majEigVec);
disp(idxPtData);
% gold idx point pos:
% id == 44616: -0.359255 -0.0206616 1.56832 0.579919  