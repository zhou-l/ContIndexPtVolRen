function Vlist = gen2FlatSynVolumes(zvolDim)
if nargin < 1
    volDim = 64;
else
    volDim= zvolDim;
end
close all;
numVols = 3;
Vlist = cell(numVols, 1);
numSamples = volDim * volDim *volDim;
%% create features with patterns in the value domain
mu = [0.1 0.2 0.1];
sig1 = 0.01;
sig2 = 0.02;
sig3 = 0.02;
rou1 = 0;%0.0;
rou2 = 0;%1;
rou3 = 0;%0.8;
cov = [ sig1*sig1, rou1*sig1*sig2, rou2*sig1*sig3;...
        rou1*sig2*sig1, sig2*sig2, rou3*sig2*sig3;...
        rou2*sig3*sig1, rou3*sig3*sig2, sig3*sig3;...
     ];
S = mvnrnd(mu,cov, numSamples);
figure, scatter(S(:,1), S(:,2)), title('S distr'); 

%% test 3D eigenvectors
mu3D = [0.8 0.4 0.6];
cov3D = createCovMat3D(0.3, 0.05, 0.02, pi/6, 0, pi/4);
F3D = mvnrnd(mu3D, cov3D, numSamples);
figure, scatter(F3D(:,1), F3D(:,2)), title('F3D subD0');
figure, scatter(F3D(:,2), F3D(:,3)), title('F3D subD1');

mu2 = [0.8 0.4];
cov2 = createCovMat2D(0.3, 0.05, pi/4);
F1_0 = mvnrnd(mu2, cov2, numSamples);
figure, scatter(F1_0(:,1), F1_0(:,2)), title('F1s0 distr'); 



mu2_1 = [0.4 0.5];
cov2_1 = createCovMat2D(0.3, 0.1, -pi/4);
F1_1 = mvnrnd(mu2_1, cov2_1, numSamples);
F1 = cat(2,F1_0,F1_1(:,2));

mu3 = [0.4 0.6];
cov3 = createCovMat2D(0.01, 0.2, pi/6);
F2_0 = mvnrnd(mu3, cov3, numSamples);
figure, scatter(F2_0(:,1), F2_0(:,2)), title('F2s0 distr'); 
mu3_1 = [0.6 0.3];
cov3_1 = createCovMat2D(0.2, 0.05, pi/3);
F2_1 = mvnrnd(mu3_1, cov3_1, numSamples);
F2 = [F2_0, F2_1(:,2)];

mu4 = [0.6 0.3];
cov4 = createCovMat2D(0.2,0.1, -pi/3);
F3_0 = mvnrnd(mu4, cov4, numSamples);
figure, scatter(F3_0(:,1), F3_0(:,2)), title('F3s0 distr'); 
mu4_1 = [0.3 0.8];
cov4_1 = createCovMat2D(0.1,0.2, pi/3);
F3_1 = mvnrnd(mu4_1, cov4_1, numSamples);
F3 = [F3_0, F3_1(:,2)];

FX = F1_0 + F2_0 + F3_0;
figure, scatter(FX(:,1),FX(:,2)), title('Without S');
rng(0,'twister');
% P1rnd = rand(numSamples);
p1x = rand(numSamples,1);
c_1 = 1;
c_2 = 1;
c_3 = 1;
% the normal of the first plane
nc = [c_1, c_2, c_3];
nc =  nc ./ norm(nc);
c_0 = 2;
% the normal of the second plane
% c_0 = c_1 * p1x + c_2 * p1y + c_3 * p1z;
nc2 = [-1 -1 1];
nc2 = nc2 ./ norm(nc2);
p1y = rand(numSamples,1);
p1z = (c_0 - nc(1) .* p1x - nc(2) .* p1y)./nc(3);
P1 = [p1x, p1y, p1z];

p2x = rand(numSamples,1);
p2z = 2 - rand(numSamples,1);
p2y = (0 - nc2(1) .* p2x - nc2(3) .* p2z) ./ nc2(2) ;
P2 = [p2x, p2y, p2z];


%% Set three

R0 = 25;
R1 = 10;
ctr = [volDim/2 volDim/2 volDim/2];
for i  = 1:numVols
    V = zeros(volDim, volDim, volDim);
    Vlist{i,:} = V;
end

 % the background is gaussian 
for z = 1:volDim
    for y = 1:volDim
        for x = 1:volDim
            idx = randi(numSamples);
            idx2 = randi(numSamples);
            idx3 = randi(numSamples);
            idx4 = randi(numSamples);
            for d = 1:numVols
                Vlist{d}(x,y,z) = S(idx,d);
                % create features
                toCtr = [x y z] - ctr;
                dist = norm(toCtr);
                if dist >= R0
                    continue;
                end
                if dist <= R1
%                     Vlist{d}(x,y,z) = F3(idx4,d); % plane 3
                    Vlist{d}(x,y,z) = P1(idx4,d);
                    continue;
                end
%                 if x < volDim/2
%                     Vlist{d}(x,y,z) = F1(idx2,d); % feature 1
%                 else
%                     Vlist{d}(x,y,z) = F2(idx3,d); % feature 2
% 
%                 end
                Vlist{d}(x,y, z) = P2(idx3,d);
            end

        end
    end
end

% collect the whole data
 cnt = 1;
 fullDist = zeros(numSamples, numVols);
for z = 1:volDim
    for y = 1:volDim
        for x = 1:volDim
            for d = 1:numVols
                fullDist(cnt,d) = Vlist{d}(x,y,z);
            end
            
            cnt = cnt + 1;
        end
    end
end
figure, hold on;
scatter(fullDist(:,1),fullDist(:,2)), title('Whole data');
% scatter(F1(:,1),F1(:,2), 'filled', 'r');
% scatter(F2(:,1),F2(:,2), 'filled', 'g');
scatter(P1(:,1),P1(:,2), 'filled', 'r');
scatter(P2(:,1),P2(:,2), 'filled', 'g');
scatter(P1(:,2),P1(:,3), 'filled', 'b');
scatter(P2(:,2),P2(:,3), 'filled', 'y');
hold off;
%%
name = 'PlaneSynthVolN';
for i = 1:numVols
    fileName = sprintf('%s%d.nrrd', name, i);
    ok = nrrdWriter(fileName, single(Vlist{i}), [1 1 1], [0 0 0], 'raw');
    if i > 1
        figure; h = histogram2(Vlist{i-1},Vlist{i});
    end
end

end
% 
% function h = histo2D(Vol1, Vol2, numBins)
% 
% end

% psi is for angle in XY
% theta is for YZ
% phi is for XZ

function Sig = createCovMat3D(lambda1, lambda2, lambda3, theta, phi, psi)
sig1 = lambda1; %0.01;
sig2 = lambda2; %0.0001;
sig3 = lambda3;
%L: eigenvalues
%V: eigenvectors
L = [sig1, 0 0;...
     0, sig2, 0;...
     0, 0, sig3];

Vx = [1 0 0;...
    0 cos(theta), -sin(theta);...
    0 sin(theta), cos(theta)];
Vy = [cos(phi) 0 sin(phi);...
      0 1 0;...
      -sin(phi) 0 cos(phi)];
Vz = [cos(psi) -sin(psi) 0;...
      sin(psi) cos(psi) 0;...
      0 0 1];
V= Vz * Vy * Vx;
% V = [lx,-ly;ly,lx];
Sig = V * L * inv(V);
end


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