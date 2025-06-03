% sort eigenvectors by its largest coordinate
% NOTE: Use only when the principal components are isotropic!
function M = sortEigVec(D)
M = D;   
for i = 1:size(D,2)
    V = D(:,i);
    [mv, ind] = max(V);   
    M(:,ind) = V;
end
end