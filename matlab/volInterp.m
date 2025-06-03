function Ev = volInterp(Vol, pos)

p = zeros(8,3);
p(1,:) = floor(pos);
a = pos - p(1,:);
offset =...
[0 0 0;
 1 0 0;
 1 1 0;
 0 1 0;
 0 1 1;
 1 1 1;
 1 0 1;
 0 0 1
];
vals = zeros(8,1);
for i = 1:8
    p(i,:) = p(1,:) + offset(i,:);
    vals(i) = Vol(p(i,1),p(i,2),p(i,3));
end


w1 = 1 - a(1);
w2 = 1 - w1;
w3 = 1 - a(2);
w4 = 1 - w3;
w5 = 1 - a(3);
w6 = 1 - w5;

%    3-E2-2
%   /|   /|
%  0-E1-1 |
%  | 4-E3-5
%  |/   |/
%  7-E4-6
% Be careful with the order!!!
Ex1 = w1 .* vals(1,:) + w2 .* vals(2,:);
Ex2 = w2 .* vals(3,:) + w1 .* vals(4,:);
Ex3 = w1 .* vals(5,:) + w2 .* vals(6,:);
Ex4 = w2 .* vals(7,:) + w1 .* vals(8,:);

Ey1 = w3 .* Ex1 + w4 .* Ex2;
Ey2 = w4 .* Ex3 + w3 .* Ex4;

Ez = w5 .* Ey1 + w6 .* Ey2;
Ev = Ez;


end