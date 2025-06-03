function f = myLocalKDEwNN(pos, nn, D, rmax, h)

% maxD = -realmax;
neighbors = nn;
cnt = 1;
for i = 1:size(nn,1)
    if D(i)< rmax
        neighbors(cnt,:) = nn(i,:);
        cnt = cnt+1;
    end
end

neighbors(cnt:end,:) = [];
f = 0;
% assuming an normal distribution kernel
n = size(neighbors,1);
if n == 0
    return;
end
factor = 1.0;% 1.0/(h*sqrt(2*pi));
for i = 1:n
    xi = neighbors(i,:);
    xx = norm(pos - xi);%/h;
    fh = factor * exp(-(xx * xx)/(2 * h * h));
%     fh = gaussKernel(xx,h);
%     fh = epanechnikovKernel(xx);
    f = f + fh;
end
denom = 1/min(rmax,n);
f = f * denom;% (n * h);
% if n >= 1
%     f = f / (n * h);
% end

end

function  f = epanechnikovKernel(x)
f = 0.75 * (1 - x*x) ;
end

function f = gaussKernel(x, h)
factor = 1.0;% 1.0/(h*sqrt(2*pi));
f = factor * exp(-(x * x)/(2 * h * h));
end

function triangleKernel()
end