function f = myLocalKDE(X, pos, nmax, rmax, h)

% find the K-NN of the position

Idx = knnsearch(X, pos, 'K', nmax);
nn = X(Idx,:);
% maxD = -realmax;
neighbors = nn;
cnt = 1;
for i = 1:size(nn,1)
    v = pos - nn(i,:);
    d = norm(v);
    if d < rmax
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
    xx = norm(pos - xi)/h;
    fh = factor * exp(-(xx * xx)/(2 * h * h));
    f = f + fh;
end
f = f / (n*h);
% if n >= 1
%     f = f / (n * h);
% end

end