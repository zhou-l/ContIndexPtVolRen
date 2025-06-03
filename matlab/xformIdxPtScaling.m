function pt = xformIdxPtScaling(c1, c2, c3)

x = [-0.5, -0.25, 0, 0.1, 0.25, 0.5, 0.75, 0.9, 1, 1.25, 1.5];
y = [1.306  1.153  1 0.9312 0.8555 0.812 0.8555 0.9312 1 1.153 1.306];
if c2 == 0
    u = 0;
    v = -c3/c1;
    if c1 == 0
        v = 0;
    end
    pt = [u,v];
    return;
end

if c1 == -c2
    if rand() < 0.5
        u = -0.5;
    else
        u = 1.5;
    end
    v = 0.0;
     pt = [u,v];
    return;
end
a = -c1/c2;
b = 2 * c3 / (c1 - c2);
theta = atan(a);
if theta > pi/4
    u = theta*2/pi - 1;
else
    u = theta*2/pi+1;
end
s = spline(x,y,u);
v = (u - 0.5) * b * s;
 pt = [u,v];
end
