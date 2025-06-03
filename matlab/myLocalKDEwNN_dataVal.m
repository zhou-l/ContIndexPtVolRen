function f = myLocalKDEwNN_dataVal(pos, nn, D, nmax, rmax, h, I, kernelType, useDynamicBandwidth)

if nargin < 7
    zkernelType = 1;
else
    zkernelType = kernelType;
end

if nargin < 8
    zUseDynBandwidth = true;
else
    zUseDynBandwidth = useDynamicBandwidth;
end
% I: the image of data values with 'pos'
% maxD = -realmax;
neighbors = nn;
zh = h;
% discard neighbors if outside of the max radius
cnt = 1;
Dmax = D(1);
% seems that D(i) are sorted...
valSum = 0;
for i = 1:size(nn,1)
    %     Dmax = max(Dmax, D(i));
    if nn(i,1) <= 0 || nn(i,2) <= 0
        continue;
    end
    valSum = valSum + I(nn(i,1),nn(i,2));
    if valSum > nmax && cnt > 1 % should have at least one sample!
        break;
    end
    neighbors(cnt,:) = nn(i,:);
    cnt = cnt+1;
    
    Dmax = max(Dmax, D(i));
end
neighbors(cnt:end,:) = [];
% rthres = 20;
% if Dmax < rthres

% end

f = 0;
% assuming an normal distribution kernel
n = size(neighbors,1);
if n == 0
    return;
end

% dynamic bandwidth for dense regions
if zUseDynBandwidth
  zh = max(1, min(h, Dmax/2));
%   if Dmax > 6
%       zh = h;% min(h,Dmax/3);
%   end
else  
  % fixed bandwidth scheme  
  zh = h;
end
% if Dmax <= 10
%     zh = max(1,min(Dmax/3, 2));%Dmax /4;
% % zh = min(Dmax,2);
% % elseif Dmax >= rmax
% %     zh = h;%min(h, Dmax/3);
% else
%     zh = min(h, Dmax/3);
% end
Kf = {@gaussKernel, @epanechnikovKernel, @expKernel, @cosKernel, @triangleKernel, @avgKernel};

for i = 1:n
    xi = neighbors(i,:);
    %     if xi(1) <= 0 || xi(2) <= 0
    %         continue;
    %     end
    xx = norm(pos - xi);%/h;
    val = I(xi(1),xi(2));
    % gaussian kernel
    %     fh = factor * val * exp(-(xx * xx)/(2 * h * h));
    % try averaging
    
    %% Warning: the function version seems to be much slower...
    fh = Kf{zkernelType}(xx,zh,val);
    f = f + fh;
end
% denom = 1/min(Dmax,n);
denom = 1/zh;%1/min(rmax,Dmax);
% denom = 1/min(Dmax,rmax);
f = f * denom;% (n * h);
end

function  f = epanechnikovKernel(x, h, v)
if x > h
    f = 0;
    return;
end
f = v * (1 - (x*x)/(h*h)) ;

end

function f = gaussKernel(x, h, v)
factor = 1.0;% 1.0/(h*sqrt(2*pi));
f = factor * v * exp(-(x * x)/(2 * h * h));
end

function f = triangleKernel(x,h,v)
f = v * (1 - x/h);
end

function f = cosKernel(x,h,v)
if x > h
    f = 0;
    return;
end
f = max(0,v * cos(pi * x / (2 * h)));
end

function f = expKernel(x,h,v)
f = v * exp(-x/h);
end

function f = avgKernel(x,h,v)
if x > h
    f = 0;
    return;
end
f = v;
end