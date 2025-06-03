function Imv = buildTestData
H = 20;
W = 20;
I0 = zeros(H,W);
I1 = I0;
I2 = I0;
Ctr = [H/2, W/2];
% draw a single Guassian distribution in I0
sigma0 = 10.0;
for r = 1:H
    for c = 1:W
        x = [r c] - Ctr;
        d = norm(x);
%         f = exp(-0.5 *(d * d)/(sigma0*sigma0)) * 1/sigma0*sqrt(2*pi);
        if d <= 8
            f = 10;
        else
            f = 0;
        end
        I0(r,c) = f;
    end
end
% figure, imshow(I0,[]);

sigma1 = 2.0;
for r = 1:H
    for c = 1:W
        x = [r c] - Ctr;
           d = norm(x);
%         f = exp(-0.5 *(d * d)/(sigma1*sigma1)) * 1/sigma1*sqrt(2*pi);
        if d <= 4
            f = 20;
        else
            f = 0;
        end
        I1(r,c) = f;
    end
end
% figure, imshow(I1,[]);

sigma1 = 2.0;
for r = 1:H
    for c = 1:W
        x = [r c] - Ctr;
        d = norm(x);
        if d >= 5 && d <= 8
         f = exp(-0.5 *((d-6)* (d-6))/(sigma1*sigma1)) * 1/sigma1*sqrt(2*pi);
        else
          f = 0;
        end
    
        I2(r,c) = f;
    end
end
figure, imshow(I0,[]);
Imv = cell(3,1);
Imv{1} = I0;
Imv{2} = I1;
Imv{3} = I2;
end