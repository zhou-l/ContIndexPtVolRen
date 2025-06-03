function Ev = eigInterp2D(EtriList, bc) % EtriList: Eigenvectors of a triangle (or any polygon?); bc: barycentric coordinates
n = size(EtriList,1);
theta = -realmax;
for i = 1:n

    for j = i+1:n
        u = EtriList(i,:);
        v = EtriList(j,:);
        % a robust evaluation of the angle of two vectors
        % https://www.mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab
        thetai = atan2(norm(cross(u,v)),dot(u,v));
%         thetai = acos(dot(EtriList(i,:),EtriList(j,:)));
        theta = max(theta, thetai);

    end
end
Ev = zeros(1,size(EtriList,2));
% if theta is 0, meaning that all vectors before interpolation are the
% same.
if theta == 0
    Ev = EtriList(1,:);
    return;
end
for i=1:n
    w = sin(bc(i) * theta) / sin(theta);
    Ev = Ev + w .* EtriList(i,:);
end
Ev = normalize(Ev,'norm');
