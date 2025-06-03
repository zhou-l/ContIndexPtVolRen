% gridx1 = -0.25:.05:1.25;
% gridx2 = 0:.1:15;
% [x1,x2] = meshgrid(gridx1, gridx2);
% x1 = x1(:);
% x2 = x2(:);
% xi = [x1 x2];
% rng('default')  % For reproducibility
% x = [0+.5*rand(20,1) 5+2.5*rand(20,1);
%             .75+.25*rand(10,1) 8.75+1.25*rand(10,1)];
% figure
% ksdensity(x,xi);   
M = load('ImgPlaneIdxPtDebug_load2.txt');
X = M(:,1:2);

figure, plot(X(:,1),X(:,2),'.');
% fx = ones(length(X),2);
figure;
ksdensity(X, 'kernel', 'triangle', 'BandWidth', 0.02, 'PlotFcn', 'contour' );
% plot(xi, f);
% disp(h);