
function testEigInterp
E1 = [0, 0, 1, 0];
E2 = [0, 0, 0, 1];
[Ei,ev] = eigInterp(E1, E2, 1, 1, 0.5);
disp(Ei);
disp(ev);