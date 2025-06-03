function Ev = eigInterpPerElement3D(EtriList, bc) % EtriList: Eigenvectors of a triangle (or any polygon?); bc: barycentric coordinates
n = size(EtriList,1);
w1 = bc(1); %sin(bc(1) * theta)/sin(theta);
w2 = bc(2); %sin(bc(2) * theta)/sin(theta);
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

% Be careful with the order!!!
w3 = bc(3); %sin(bc(3) * theta)/sin(theta);
w4 = bc(4); %sin(bc(4) * theta)/sin(theta);
Ey1 = w3 .* Ex1 + w4 .* Ex2;
Ey2 = w4 .* Ex3 + w3 .* Ex4;

w5 = bc(5); %sin(bc(5) * theta)/sin(theta);
w6 = bc(6); %sin(bc(6) * theta)/sin(theta);
Ez = w5 .* Ey1 + w6 .* Ey2;
Ev = Ez;

Ev = normalize(Ev,'norm');
end