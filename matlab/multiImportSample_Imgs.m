% I0 = imread('VolData0-80.ppm');
% I1 = imread('VolData1-80.ppm');
% I2 = imread('VolData2-80.ppm');

I0 = imread('isabel/VolData0-45.ppm');
I1 = imread('isabel/VolData1-45.ppm');
I2 = imread('isabel/VolData2-45.ppm');

I0 = rgb2gray(I0);
I1 = rgb2gray(I1);
I2 = rgb2gray(I2);
close all;
% Use gradient as the target pdf
gI0 = imgradient(I0, 'sobel');
gI1 = imgradient(I1, 'sobel');
gI2 = imgradient(I2, 'sobel');
figure;

subplot(2, 3, 1); imshow(I0, []);title('Scalar Data');
subplot(2, 3, 2); imshow(I1, []);
subplot(2, 3, 3); imshow(I2, []);


subplot(2, 3, 4); imshow(gI0, []);title('Gradient Magnitude Data');
subplot(2, 3, 5); imshow(gI1, []);
subplot(2, 3, 6); imshow(gI2, []);

N = 10^5;
% Multiple importance sampling, using the maximum?

heursiticType = 0;


multiImg = max(max(gI0,gI1), gI2);
% Or trying the power heuristic?
% inline Float PowerHeuristic(int nf, Float fPdf, int ng, Float gPdf) {
%     Float f = nf * fPdf, g = ng * gPdf;
%     return (f * f) / (f * f + g * g);
% }
sumImg = gI0 + gI1 + gI2;

maxR = size(multiImg, 1) - 1;
maxC = size(multiImg, 2) -1;


figure,imshow(multiImg, []);
mean_samples = rand(N,1);%(rand(N,1)-0.5)*5; 
sigma_samples = rand(N,1); %rand(N, 1) * 10; 
proposal = 1/N;
idxR = int32(maxR .* mean_samples) + 1;
idxC = int32(maxC .* sigma_samples)+ 1;

targetPower = zeros(N,1);
targetSum = zeros(N,1);
targetMax = zeros(N,1);

nList = [N; N; N];
m = size(nList, 1);
for i=1:N
    targetSum(i,:) = sumImg(idxR(i), idxC(i));
%     % with maximum heuristic
    targetMax(i,:) = multiImg(idxR(i), idxC(i));
    % try MIP with power heuristic
    pdfList = [gI0(idxR(i), idxC(i)); gI1(idxR(i), idxC(i)); gI2(idxR(i), idxC(i))];
    wList = PowerHeuristicMPdfs(nList, pdfList);
    targetPower(i,:) = 0;
    for d = 1:m
        targetPower(i,:) = targetPower(i,:) + wList(d) * pdfList(d) ./ m;
    end
end


target = targetPower;
% calculate importance weight
w = target ./ proposal;
% normalization
w = w ./ sum(w);

wPower = targetPower ./ proposal;
wPower = wPower ./ sum(wPower);

wSum = targetSum ./ proposal;
wSum = wSum ./ sum(wSum);

wMax = targetMax./ proposal;
wMax = wMax ./ sum(wMax);

% resample, with replacement, according to importance weight
sample_ind = randsample([1:N],N,true,w);

sample_ind_Power = randsample([1:N],N,true,wPower);
sample_ind_Sum = randsample([1:N],N,true,wSum);
sample_ind_Max = randsample([1:N],N,true,wMax);

mean_samples = mean_samples(sample_ind);
sigma_samples = sigma_samples(sample_ind);

y_max_samples = mean_samples(sample_ind_Max);
y_sum_samples = mean_samples(sample_ind_Sum);
y_power_samples = mean_samples(sample_ind_Power);
x_max_samples = sigma_samples(sample_ind_Max);
x_sum_samples = sigma_samples(sample_ind_Sum);
x_power_samples = sigma_samples(sample_ind_Power);

%% plot
hold on;
plot(maxC.*sigma_samples+1, maxR.*(mean_samples)+1,'.')
xlabel('mean')
ylabel('sigma')
hold on

%plot(true_mean, true_sigma, 'r.','MarkerSize',5^2)
axis square

figure;
subplot(1,3,1); plot(maxC.*x_max_samples+1, maxR.*(y_max_samples)+1,'.'); title("Max MIP");
subplot(1,3,2); plot(maxC.*x_sum_samples+1, maxR.*(y_sum_samples)+1,'.'); title("Sum MIP");
subplot(1,3,3); plot(maxC.*x_power_samples+1, maxR.*(y_power_samples)+1,'.'); 
title('Power MIP');


function w= PowerHeuristic2Pdfs(nf, fPdf, ng, gPdf) 
 f = nf * fPdf;
 g = ng * gPdf;
 w = (f * f) / (f * f + g * g);
end

function wList = PowerHeuristicMPdfs(nList, pdfList)
m = size(nList,1);
fList = zeros(m,1);
wList = zeros(m,1);
denom = 0;
for i =1:m
    fList(i) = nList(i) * pdfList(i);
    denom = denom + fList(i)*fList(i);
end
%  w = (f * f) / (f * f + g * g);
if denom == 0
    return;
end
for i = 1:m
    wList(i) = (fList(i)*fList(i)) / denom;
end
end