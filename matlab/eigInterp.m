function [Ei,evi] = eigInterp(E1, E2, eigVal1, eigVal2, t)
cosTheta = dot(E1, E2);
theta = acos(cosTheta);
Ei = sin((1-t)*theta)/sin(theta) .* E1 + sin(t * theta)/sin(theta) .* E2;
evi = exp((1-t)*log(eigVal1)+t*log(eigVal2));
