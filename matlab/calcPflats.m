function idxVol = calcPflats(mu, majEig, zSecEig)
if nargin < 3
    secEig = zeros(size(majEig,1));
else
    secEig = zSecEig;
end
ptMajEig = mu + majEig;
vecLen = vecnorm(ptMajEig,2,2);
ptMajEig = ptMajEig ./ vecLen;

%
ptSecEig = mu + secEig;

dimM = size(mu,2);
N = size(mu,1);
dimIdxPt = 2 * (dimM - 1);
idxVol = zeros(N,dimIdxPt);
% compute 1-flats
pw = 3; % width of homogeneous coordinates of a subdimension 
ww = 2; % width of coordinates of the image plane 
oneFlats = calcOneFlatGeneralForm(mu, ptMajEig);

% % get twoFlats
% twoFlats = calcTwoFlatGeneralForm(mu, ptMatEig, ptSecEig);

%calc for each data item
for i = 1:size(oneFlats,1)
    for j = 1:dimM-1
       col = pw*(j-1)+1;
       coords = oneFlats(i,col:col+pw-1);
%        idPt = xformIdxPtScaling(coords(1),coords(2),coords(3));
        idPt = xformIdxPtNoScaling(coords(1),coords(2),coords(3));
       idxVol(i,ww*(j-1)+1:ww*j) = idPt;
    end
end


end


function oneFlats =  calcOneFlatGeneralForm(cartP1, cartP2)

    pw = 3; % width of homogeneous coordinates of a subdimension 
    oneFlats = zeros(size(cartP1,1), pw * (size(cartP1,2)-1));

    for i = 1:size(cartP1, 2)-1
        p1 = cartP1(:,i:i+1);
        p2 = cartP2(:,i:i+1);
        
        c1 = p2(:,2) - p1(:,2);
        c2 = -(p2(:,1) - p1(:,1));
        c3 = p1(:,2) .* (p2(:,1) - p1(:,1)) - p1(:,1) .* (p2(:,2) - p1(:,2));

        col = pw*(i-1)+1;
%         nzIdx = find(c2);
        
            m = -c1./c2;
            b = -c3./c2;
            c1 = -m;
            c2 = ones(length(c1),1);
            c3 = -b;
       
        Cc = [c1,c2,c3];
        oneFlats(:, col:col+pw-1) = Cc;
    end
end

function twoFlats = calcTwoFlatGeneralForm()
end
