function Ev = eigInterp3D(EtriList, bc) % EtriList: Eigenvectors of a triangle (or any polygon?); bc: barycentric coordinates
n = size(EtriList,1);
theta = -realmax;
for i = 1:n

    for j = i+1:n
        u = EtriList(i,:);
        v = EtriList(j,:);
        % a robust evaluation of the angle of two vectors (of 3D)
        % https://www.mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab
%         thetai = atan2(norm(cross(u,v)),dot(u,v));
        % N-D vectors
        thetai = acos(dot(u,v));
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
%trilinear interpolation
% for i=1:length(bc)
%     w = sin(bc(i) * theta) / sin(theta);
%     Ev = Ev + w .* EtriList(i,:);
% end

w1 = sin(bc(1) * theta)/sin(theta);
w2 = sin(bc(2) * theta)/sin(theta);
w3 = sin(bc(3) * theta)/sin(theta);
w4 = sin(bc(4) * theta)/sin(theta);
w5 = sin(bc(5) * theta)/sin(theta);
w6 = sin(bc(6) * theta)/sin(theta);

%    3-E2-2
%   /|   /|
%  0-E1-1 |
%  | 4-E3-5
%  |/   |/
%  7-E4-6
% Be careful with the order!!!
Ex1 = w1 .* EtriList(1,:) + w2 .* EtriList(2,:);
Ex2 = w2 .* EtriList(3,:) + w1 .* EtriList(4,:);
Ex3 = w1 .* EtriList(5,:) + w2 .* EtriList(6,:);
Ex4 = w2 .* EtriList(7,:) + w1 .* EtriList(8,:);

Ey1 = w3 .* Ex1 + w4 .* Ex2;
Ey2 = w4 .* Ex3 + w3 .* Ex4;

Ez = w5 .* Ey1 + w6 .* Ey2;
Ev = Ez;

Ev = normalize(Ev,'norm');
end