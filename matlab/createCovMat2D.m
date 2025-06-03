function Sig = createCovMat2D(lambda1, lambda2, theta)
sig1 = lambda1; %0.01;
sig2 = lambda2; %0.0001;

%L: eigenvalues
%V: eigenvectors
L = [sig1, 0;0, sig2];

lx = cos(theta);
ly = sin(theta);
V = [lx,-ly;ly,lx];
Sig = V * L * inv(V);
end